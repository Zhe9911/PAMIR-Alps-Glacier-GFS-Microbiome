# 16S NMDS analysis
# Manuscript output: Figure 2a.

rm(list = ls())
gc()
graphics.off()

input_path <- file.path(
  "results", "01_16s", "intermediate", "PAMIR_16S_final.rds"
)
result_dir <- file.path("results", "01_16s", "nmds")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

# --- Setup ---
library(phyloseq)
library(ggplot2)
library(vegan)
library(tidyverse)
library(ggrepel)

# --- Load input data ---

ps <- readRDS(input_path)
print(ps)

# --- Merge replicate samples ---

PAMIR_16S_merged <- ps

# Collapse technical replicates to a shared base sample ID.
group_vec <- gsub("_.*", "", sample_names(PAMIR_16S_merged))

sample_data(PAMIR_16S_merged)$Sample_Base <- group_vec

# Merge counts across replicates.
PAMIR_16S_merged_phyloseq <- merge_samples(PAMIR_16S_merged, "Sample_Base")

# Restore metadata after merging samples.
original_meta <- data.frame(sample_data(PAMIR_16S_merged))

# Keep one metadata record per merged sample.
restored_meta2 <- original_meta[!duplicated(original_meta$Sample_Base), ]

# Match metadata row names to the merged phyloseq object.
rownames(restored_meta2) <- restored_meta2$Sample_Base
restored_meta2$Sample_ID <- NULL
restored_meta2 <- rename(restored_meta2, Sample_ID = Sample_Base)

# Reattach metadata to the merged phyloseq object.
sample_data(PAMIR_16S_merged_phyloseq) <- sample_data(restored_meta2)

# Preserve OTU table orientation if merge_samples flips it.
if (!taxa_are_rows(PAMIR_16S_merged_phyloseq) && taxa_are_rows(PAMIR_16S_merged)) {
  otu_table(PAMIR_16S_merged_phyloseq) <- otu_table(
    t(otu_table(PAMIR_16S_merged_phyloseq)),
    taxa_are_rows = TRUE
  )
}

print(PAMIR_16S_merged_phyloseq)

# --- NMDS analysis ---

# Restrict the ordination to ice and sediment samples.
ps_NMDS <- subset_samples(PAMIR_16S_merged_phyloseq, Source != "water")
ps_NMDS <- prune_taxa(taxa_sums(ps_NMDS) > 0, ps_NMDS)
ps_NMDS

sample_data_NMDS <- as.data.frame(sample_data(ps_NMDS))

# Transform counts to relative abundance before Bray-Curtis calculation.
ps_NMDS_rel <- transform_sample_counts(ps_NMDS, function(x) x / sum(x))

# Calculate the Bray-Curtis distance matrix.
set.seed(666)
bray_dist_NMDS <- distance(ps_NMDS_rel, method = "bray")

# Report stress for k = 2:4 as a diagnostic; this does not select dimensionality.
stress_values <- list()
k_values <- 2:4

for (k in k_values) {
  message(paste("Trying NMDS with k =", k))
  nmds_try <- metaMDS(bray_dist_NMDS, k = k, trymax = 100)
  stress_values[[as.character(k)]] <- nmds_try$stress
  message(paste("  Stress =", round(nmds_try$stress, 2)))
}

# Use the manuscript-specified two-dimensional solution for the final plot.
nmds <- metaMDS(bray_dist_NMDS, k = 2, trymax = 100)
message(paste("Final stress value:", round(nmds$stress, 2)))

pdf(
  file.path(result_dir, "figure_2a_nmds_stressplot.pdf"),
  width = 8,
  height = 8
)
stressplot(nmds)
dev.off()

# Extract ordination coordinates and join sample metadata.
nmds_coords <- as.data.frame(nmds$points)
colnames(nmds_coords) <- paste0("NMDS", 1:2)

nmds_coords$sample_id <- rownames(nmds_coords)
nmds_coords <- merge(
  nmds_coords,
  sample_data_NMDS,
  by.x = "sample_id",
  by.y = "row.names"
)
write.csv(
  data.frame(k = as.integer(names(stress_values)), stress = unlist(stress_values)),
  file.path(result_dir, "figure_2a_nmds_dimension_stress.csv"),
  row.names = FALSE
)
write.csv(
  nmds_coords,
  file.path(result_dir, "figure_2a_nmds_coordinates.csv"),
  row.names = FALSE
)

# Order glacier-position groups from ice to lower-glacier sediments.
desired_order <- c("ice", "High_Glacial_sed", "Mid_Glacial_sed", "Low_Glacial_sed")

# Apply the plotting order to the glacier-position groups.
nmds_coords$Group_GI <- factor(nmds_coords$Group_GI, levels = desired_order)

# Order points by glacier-position class, not raw GI.
# Ice samples have missing GI, so sorting by -GI would push them to the end
# and make geom_path connect the last sediment point back to ice.
nmds_coords <- nmds_coords[
  order(nmds_coords$gl_name, nmds_coords$Group_GI, -nmds_coords$GI),
]

# Plot NMDS points and connect samples from the same glacier.
p6 <- ggplot(nmds_coords, aes(x = NMDS1, y = NMDS2)) +
  geom_point(
    aes(color = Group_GI, shape = Region),
    size = 3,
    alpha = 0.9,
    stroke = 1.1
  ) +
  geom_path(aes(group = gl_name, color = Group_GI), alpha = 1) +
  geom_text_repel(
    data = subset(nmds_coords, Group_GI == "ice"),
    aes(label = gl_name),
    size = 3,
    color = "black",
    max.overlaps = Inf
  ) +
  scale_shape_manual(values = c(16, 17)) +
  scale_color_viridis_d(option = "D", begin = 0.15, end = 0.85) +
  theme_minimal() +
  labs(
    title = "NMDS Ordination of 16S Microbial Communities",
    subtitle = paste("Stress =", round(nmds$stress, 2)),
    color = "Group by GI",
    shape = "Region"
  ) +
  theme(legend.position = "right")

ggsave(
  file.path(result_dir, "figure_2a_nmds_connected.pdf"),
  p6,
  width = 8,
  height = 7
)
