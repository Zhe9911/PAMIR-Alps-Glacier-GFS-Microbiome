# MAG phylogeny visualization
# Manuscript output: Figure 3a circular tree.

rm(list = ls())
gc()
graphics.off()

out_dir <- file.path("results", "03_mags", "phylogeny")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- Packages ---
library(tidyverse)
library(phyloseq)
library(RColorBrewer)
library(ggplot2)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)
library(ape)

# --- Load inputs and set display parameters ---

PAMIR_MAGs_rela <- readRDS(
  file.path("results", "03_mags", "intermediate", "PAMIR_MAGs_rela.rds")
)
print(PAMIR_MAGs_rela)

# Load MAG metadata and define the ecological-trend display order.
MAGs_info <- read.csv(
  file.path("results", "03_mags", "intermediate", "MAGs_info.csv"),
  stringsAsFactors = FALSE
)
stopifnot(all(c("MAGs", "Trend", "lfc_GI", "Phylum") %in% names(MAGs_info)))

MAGs_info$Trend <- factor(
  MAGs_info$Trend,
  levels = c("Glaciophiles", "Non-sig", "Glaciophobes"),
  ordered = TRUE
)

# Use a shared symmetric LFC-GI scale across the phylogenetic ring (Fig. 3a)
# and the LFC-GI histogram (Fig. 3b).
lfc_GI_color_limit <- max(abs(MAGs_info$lfc_GI), na.rm = TRUE)
lfc_GI_color_limits <- c(-lfc_GI_color_limit, lfc_GI_color_limit)

phylo_tree <- read.tree(
  file.path("data", "processed", "metagenomics", "iqtree_tree.treefile")
)

PAMIR_MAGs_tree <- merge_phyloseq(PAMIR_MAGs_rela, phylo_tree)
stopifnot(setequal(taxa_names(PAMIR_MAGs_tree), phylo_tree$tip.label))

tax_df <- as.data.frame(as(tax_table(PAMIR_MAGs_tree), "matrix"))
tax_df$OTUID <- rownames(tax_df)
stopifnot("Phylum" %in% names(tax_df))

# --- Root the tree ---

# Root the tree with Patescibacteriota as the outgroup.

patescibacteria_tips <- tax_df %>%
  filter(
    !is.na(Phylum),
    str_detect(Phylum, fixed("Patescibacteriota"))
  ) %>%
  pull(OTUID) %>%
  intersect(phylo_tree$tip.label) %>%
  unique()
stopifnot(length(patescibacteria_tips) > 0)

# Resolve the root as a bifurcation for clearer visualization.
rooted_tree <- root(phylo_tree, outgroup = patescibacteria_tips, resolve.root = TRUE)

# Store the rooted tree in the phyloseq object.
phy_tree(PAMIR_MAGs_tree) <- rooted_tree

# --- 1. Prepare taxonomy and abundance annotations ---

# Extract taxonomy with an explicit MAG identifier.
tax_df <- as.data.frame(tax_table(PAMIR_MAGs_tree)) %>%
  rownames_to_column("OTUID") # Use the MAG identifier for downstream joins.

# Sum relative abundance across samples for each MAG.
otu_abundance <- data.frame(
  OTUID = taxa_names(PAMIR_MAGs_tree),
  TotalAbundance = taxa_sums(PAMIR_MAGs_tree)
)

# Combine taxonomy and total abundance.
merged_df <- left_join(tax_df, otu_abundance, by = "OTUID")

# Retain the most abundant phyla and collapse the remainder into "Other".
target_rank <- "Phylum"
N_top_phyla <- 17

# Rank phyla by total abundance.
top_phyla_list <- merged_df %>%
  group_by(!!sym(target_rank)) %>%
  summarise(TotalReads = sum(TotalAbundance, na.rm = TRUE)) %>%
  arrange(desc(TotalReads)) %>%
  slice_head(n = N_top_phyla) %>% # Retain the top N phyla.
  pull(!!sym(target_rank))

# Create simplified phylum labels.
tax_df <- tax_df %>%
  mutate(
    Simplified_Phylum = ifelse(
      !!sym(target_rank) %in% top_phyla_list,
      !!sym(target_rank),
      "Other"
    )
  )

# --- 2. Define the phylum color palette ---

# Retain abundance-based phylum order.
phyla_names <- intersect(top_phyla_list, unique(tax_df$Simplified_Phylum))

# Interpolate the palette to the required number of phyla.
main_colors <- colorRampPalette(brewer.pal(12, "Paired"))(length(phyla_names))

# Create a named color vector.
phylum_colors <- setNames(main_colors, phyla_names)

# Assign light gray to collapsed phyla.
phylum_colors["Other"] <- "#D3D3D3" # Light gray.

# Keep "Other" at the end of the legend.
tax_df$Simplified_Phylum <- factor(tax_df$Simplified_Phylum, levels = c(phyla_names, "Other"))

# --- 3. Group tree branches by phylum ---

# Group MAG identifiers by simplified phylum.
otu_list_for_grouping <- split(tax_df$OTUID, tax_df$Simplified_Phylum)

# Map phylum groups onto the phylo object.
tree_obj <- phy_tree(PAMIR_MAGs_tree)
tree_grouped <- groupOTU(tree_obj, otu_list_for_grouping)

# --- 4. Build the circular phylogenetic tree ---

# Build the base tree with phylum-colored branches.
p1 <- ggtree(
  tree_grouped,
  layout = "fan",
  open.angle = 8,       # Leave space for the legend.
  size = 0.5,           # Use thin branch lines.
  aes(color = group)    # Use phylum groups assigned by groupOTU().
) +
  geom_treescale(width = 0.1, x = NULL, y = NULL) +
  scale_color_manual(
    values = phylum_colors,
    name = "Phylum",
    # Enlarge branch symbols in the legend.
    guide = guide_legend(override.aes = list(size = 4, linewidth = 2), order = 1)
  ) +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12, face = "bold")
  )

# Add the phylum annotation ring.
p2 <- p1 +
  ggnewscale::new_scale_fill() + # Start a new fill scale.
  geom_fruit(
    data = tax_df,
    geom = geom_tile,            # Draw a continuous color strip.
    mapping = aes(y = OTUID, fill = Simplified_Phylum),
    width = 0.07,                # Set the tile width.
    offset = 0.05,               # Set the distance from tree tips.
    pwidth = 0.1                 # Set the relative ring width.
  ) +
  scale_fill_manual(
    values = phylum_colors,
    name = "Phylum",
    guide = "none"               # Suppress the duplicate phylum legend.
  )

# --- 5. Mark MAGs with significant directional GI associations ---

# Standardize the MAG identifier.
MAGs_list <- MAGs_info %>% rename(OTUID = MAGs)

# Retain robust positive and negative GI associations.
Trend_MAGs <- MAGs_list %>%
  select(OTUID, Trend) %>%
  filter(Trend %in% c("Glaciophiles", "Glaciophobes"))

print(paste("Number of MAGs to annotate:", nrow(Trend_MAGs)))

# Attach trend metadata to tree tips.
p3 <- p2 %<+% Trend_MAGs +

  # Start a fill scale for ecological trends.
  ggnewscale::new_scale_fill() +

  # Mark directionally associated MAGs.
  geom_tippoint(
    # Plot only tips represented in Trend_MAGs.
    mapping = aes(
      subset = !is.na(Trend),
      fill = Trend
    ),
    size = 3.0,
    shape = 21,       # Use a filled circle with an independent border.
    color = "grey",
    stroke = 0.2,
    alpha = 0.9
  ) +

  # Use colors consistent with the LFC-GI scale.
  scale_fill_manual(
    values = c("Glaciophiles" = "#D53E4F", "Glaciophobes" = "#3288BD"),
    name = "Ecological Trend"
  )

# --- Add the LFC-GI annotation ring ---

tree_tips <- p3$data %>% filter(isTip) %>% pull(label)
MAGs_info_filtered <- MAGs_list %>% filter(OTUID %in% tree_tips)

p4 <- p3 +
  # Start a fill scale for continuous LFC-GI.
  ggnewscale::new_scale_fill() +

  # Add the outer LFC-GI ring.
  ggtreeExtra::geom_fruit(
    data = MAGs_info_filtered,
    geom = geom_tile,
    mapping = aes(
      y = OTUID,
      fill = lfc_GI,    # Map LFC-GI estimates to color.
      x = "LFC"         # Use a fixed position for the annotation ring.
    ),

    # Set the ring layout.
    pwidth = 0.07,
    offset = 0.03,

    # Set the ring axis and grid.
    axis.params = list(
      axis = "x",
      text.size = 3,
      text.angle = 0,
      vjust = 0.5,
      line.color = NA   # Retain the label without drawing an axis line.
    ),
    grid.params = list()
  ) +

  # Apply a zero-centered diverging LFC-GI scale.
  scale_fill_gradient2(
    low = "#3288BD",   # Negative LFC-GI estimate.
    mid = "white",     # Zero LFC-GI estimate.
    high = "#D53E4F",  # Positive LFC-GI estimate.
    midpoint = 0,
    limits = lfc_GI_color_limits,
    name = "LFC GI",
    guide = guide_colorbar(
      order = 4,
      barwidth = 0.5,
      barheight = 5,
      title.position = "top",
      frame.colour = "black",
      frame.linewidth = 0.5
    )
  )

# --- 6. Position phylum labels on the circular tree ---

# Extract tip coordinates from the plotted tree.
plot_data <- p2$data %>%
  dplyr::filter(isTip) %>%
  dplyr::select(OTUID = label, angle)

# Add simplified phylum assignments.
label_data <- left_join(plot_data,
                        tax_df %>% select(OTUID, Simplified_Phylum),
                        by = "OTUID")

# Calculate circular means across the 0°/360° boundary.
circular_mean_degrees <- function(angle) {
  angle_radians <- angle * pi / 180
  mean_angle <- atan2(
    mean(sin(angle_radians), na.rm = TRUE),
    mean(cos(angle_radians), na.rm = TRUE)
  ) * 180 / pi
  (mean_angle + 360) %% 360
}

# Select the tip nearest the circular mean of each phylum.
phylum_labels <- label_data %>%
  filter(Simplified_Phylum != "Other", !is.na(angle)) %>%
  group_by(Simplified_Phylum) %>%
  summarise(
    mean_angle = circular_mean_degrees(angle),
    representative_index = which.min(
      abs(((angle - mean_angle + 180) %% 360) - 180)
    ),
    OTUID_rep = OTUID[representative_index],
    angle_rep = angle[representative_index],
    .groups = "drop"
  ) %>%
  select(-representative_index)

# Rotate labels to remain upright on both sides of the tree.
phylum_labels <- phylum_labels %>%
  mutate(
    text_angle = case_when(
      angle_rep >= 90 & angle_rep <= 270 ~ angle_rep + 180,
      TRUE ~ angle_rep
    ),
    hjust_param = 0.5
  )

# --- 7. Add phylum labels ---

# Place labels over the phylum ring without allocating another track.

p5 <- p4 +
  geom_fruit(
    data = phylum_labels,
    geom = geom_text,
    mapping = aes(
      y = OTUID_rep,
      label = Simplified_Phylum,
      angle = text_angle
    ),
    hjust = 0.5,
    vjust = 0.5,
    size = 4,
    fontface = "bold",
    color = "black",

    # Set label placement.
    offset = -0.25,
    pwidth = 0
  )

# --- Save figure ---
ggsave(
  file.path(out_dir, "Final_MAGs_Phylo_iqtree.pdf"),
  plot = p5,
  width = 16,
  height = 16, # Use a large canvas for dense tip annotations.
  device = cairo_pdf, # Preserve vector graphics and transparency.
  limitsize = FALSE
)
