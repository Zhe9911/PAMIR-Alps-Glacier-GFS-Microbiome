# Project significant KEGG GSEA signals back onto MAGs and lineage summaries
# Manuscript outputs: Figure 4 and Extended Data Table 4.

rm(list = ls())
gc()
graphics.off()

out_dir <- file.path("results", "04_functions", "phylogenetic_projection")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- Setup ---

library(tidyverse)
library(ape)
library(phyloseq)
library(ggtree)
library(patchwork)
library(RColorBrewer)

# --- 1. Read source tables and significant GSEA outputs ---

ko_mag_dt <- read_tsv(
  file.path("data", "processed", "metagenomics", "ko_mag_dt.tsv"),
  show_col_types = FALSE
)
MAGs_info <- read.csv(
  file.path("results", "03_mags", "intermediate", "MAGs_info.csv"),
  stringsAsFactors = FALSE
)

pathway_gsea_sig_df <- read_csv(
  file.path("results", "04_functions", "ko_gsea", "KO_ranked_pathway_GSEA_sig.csv"),
  show_col_types = FALSE
)

module_gsea_sig_df <- read_csv(
  file.path("results", "04_functions", "ko_gsea", "KO_ranked_module_GSEA_sig.csv"),
  show_col_types = FALSE
)

main_text_feature_metadata <- read_csv(
  file.path("results", "04_functions", "ko_gsea", "Supplementary Data 2.csv"),
  show_col_types = FALSE
) %>%
  # Reuse the original KEGG names, broad functional groups, and manually
  # specified result order from the GSEA script.
  transmute(
    feature_id = as.character(KEGG_ID),
    KEGG_name = as.character(KEGG_name),
    Broad_functional_group = as.character(Broad_functional_group),
    feature_order = row_number()
  ) %>%
  distinct(feature_id, .keep_all = TRUE)

# --- 2. Collect significant KEGG features from pathway and module GSEA ---

# Each feature is represented by its GSEA core-enrichment KO set rather than its
# full KEGG membership. This is a descriptive projection of the KOs driving
# enrichment in 01_ko_gsea.R, not a pathway activity score.
#
# Exclude weakly significant and ecologically less interpretable features from
# all projection outputs. These near-threshold or broad central-metabolism
# signals are better treated as omitted low-priority hits for the glacier
# ice-to-downstream sediment system rather than carried into lineage summaries.
excluded_features <- c(
  "ko03060", # Protein export; weak FDR and only a small core-enrichment KO set.
  "M00620",  # Incomplete reductive citrate cycle; near-threshold module signal.
  "M00932",  # Phylloquinone biosynthesis; near-threshold module signal.
  "ko00620", # Pyruvate metabolism; broad central metabolism and weak FDR.
  "ko00270", # Cysteine and methionine metabolism; broad and weakly significant.
  "ko00670", # Folate-mediated one-carbon metabolism; weakly significant.
  "ko00791", # Atrazine degradation; retained only in the complete GSEA table.
  "ko00627", # Aminobenzoate degradation; retained only in the complete GSEA table.
  "ko00196", # Photosynthetic antenna proteins; lineage-bound signal.
  "M00145"   # Cyanobacterial-type NAD(P)H:quinone oxidoreductase; lineage-bound signal.
)

significant_features <- bind_rows(
  pathway_gsea_sig_df %>%
    transmute(
      feature_type = "pathway",
      feature_id = as.character(ID),
      expected_direction = as.character(Direction),
      core_enrichment = as.character(core_enrichment)
    ),
  module_gsea_sig_df %>%
    transmute(
      feature_type = "module",
      feature_id = as.character(ID),
      expected_direction = as.character(Direction),
      core_enrichment = as.character(core_enrichment)
    )
) %>%
  left_join(main_text_feature_metadata, by = "feature_id") %>%
  filter(
    !is.na(feature_id),
    feature_id != "",
    !is.na(KEGG_name),
    KEGG_name != "",
    !is.na(Broad_functional_group),
    Broad_functional_group != "",
    !feature_id %in% excluded_features,
    expected_direction %in% c("GI-positive", "GI-negative"),
    !is.na(core_enrichment),
    core_enrichment != ""
) %>%
  distinct(feature_type, feature_id, .keep_all = TRUE)

# --- 3. Prepare MAG metadata and KO presence lookup ---

mag_annot_df <- MAGs_info %>%
  transmute(
    MAGs = as.character(MAGs),
    Trend = as.character(Trend),
    lfc_GI = as.numeric(lfc_GI),
    MAG_completeness = as.numeric(Completeness) / 100,
    genome_size = as.numeric(Genome_Size),
    log10_genome_size = log10(as.numeric(Genome_Size)),
    Phylum = as.character(Phylum),
    Class = as.character(Class),
    Order = as.character(Order),
    Family = as.character(Family),
    Genus = as.character(Genus)
  ) %>%
  # Keep only MAGs that can be used both for GI-gradient analyses and for
  # completeness/genome-size adjustment.
  filter(
    !is.na(MAGs),
    !is.na(lfc_GI),
    !is.na(MAG_completeness),
    !is.na(log10_genome_size)
  ) %>%
  distinct(MAGs, .keep_all = TRUE)

ko_presence_long <- ko_mag_dt %>%
  pivot_longer(
    cols = -KO,
    names_to = "MAGs",
    values_to = "KO_present"
  ) %>%
  transmute(
    KO = as.character(KO),
    MAGs = as.character(MAGs),
    KO_present = as.integer(KO_present)
  ) %>%
  # Store only positive KO calls; missing joins below are treated as absences.
  filter(KO_present == 1) %>%
  distinct(KO, MAGs, .keep_all = TRUE)

# --- 4. Build core-enrichment feature-to-KO mappings ---

feature_ko_map <- significant_features %>%
  select(
    feature_type,
    feature_id,
    KEGG_name,
    Broad_functional_group,
    feature_order,
    core_enrichment
  ) %>%
  # fgsea stores core-enrichment KOs as slash-separated strings.
  separate_rows(core_enrichment, sep = "/") %>%
  transmute(
    feature_type,
    feature_id,
    KEGG_name,
    Broad_functional_group,
    feature_order,
    KO = as.character(core_enrichment)
  ) %>%
  filter(!is.na(KO), KO != "") %>%
  distinct(
    feature_type,
    feature_id,
    KEGG_name,
    Broad_functional_group,
    feature_order,
    KO
  )

# --- 5. Project core-enrichment features onto MAGs ---

# core_enrichment_coverage is the fraction of core-enrichment KOs detected in each MAG.

feature_sizes <- feature_ko_map %>%
  count(
    feature_type,
    feature_id,
    KEGG_name,
    Broad_functional_group,
    feature_order,
    name = "n_total_kos"
  )

mag_feature_scores <- expand_grid(
  # Evaluate every retained feature in every MAG, including zero-coverage cases.
  MAGs = mag_annot_df$MAGs,
  feature_sizes
) %>%
  left_join(
    significant_features %>%
      select(feature_type, feature_id, expected_direction),
    by = c("feature_type", "feature_id")
  ) %>%
  left_join(
    feature_ko_map %>%
      select(
        feature_type,
        feature_id,
        KEGG_name,
        Broad_functional_group,
        feature_order,
        KO
      ),
    by = c(
      "feature_type",
      "feature_id",
      "KEGG_name",
      "Broad_functional_group",
      "feature_order"
    ),
    relationship = "many-to-many"
  ) %>%
  left_join(ko_presence_long, by = c("MAGs", "KO")) %>%
  mutate(KO_present = replace_na(KO_present, 0L)) %>%
  group_by(
    MAGs,
    feature_type,
    feature_id,
    KEGG_name,
    Broad_functional_group,
    feature_order,
    expected_direction
  ) %>%
  summarise(
    n_total_kos = {
      u <- unique(n_total_kos)
      stopifnot("n_total_kos not unique within feature" = length(u) == 1L)
      u
    },
    n_present_kos = sum(KO_present),
    core_enrichment_coverage = n_present_kos / n_total_kos,
    .groups = "drop"
  ) %>%
  left_join(mag_annot_df, by = "MAGs")

write_csv(mag_feature_scores, file.path(out_dir, "MAG_feature_scores_long.csv"))

# --- 6. Plot phylogeny and KEGG core-enrichment KO coverage ---

phylo_tree <- read.tree(
  file.path("data", "processed", "metagenomics", "iqtree_tree.treefile")
)
tree_tips <- phylo_tree$tip.label

tree_tax_base_df <- tibble(MAGs = tree_tips) %>%
  left_join(
    mag_annot_df %>%
      select(MAGs, Phylum),
    by = "MAGs"
  )

# Root the tree on Patescibacteriota to keep the phylogenetic display
# biologically interpretable and stable across plotting runs.
outgroup_tips <- tree_tax_base_df %>%
  filter(Phylum == "Patescibacteriota") %>%
  pull(MAGs)

plot_tree <- root(phylo_tree, outgroup = outgroup_tips, resolve.root = TRUE)

PAMIR_MAGs_rela_for_colors <- readRDS(
  file.path("results", "03_mags", "intermediate", "PAMIR_MAGs_rela.rds")
)
PAMIR_MAGs_tree_for_colors <- merge_phyloseq(
  PAMIR_MAGs_rela_for_colors,
  phylo_tree
)

tax_df_for_colors <- as.data.frame(
  as(tax_table(PAMIR_MAGs_tree_for_colors), "matrix")
) %>%
  rownames_to_column("MAGs")

otu_abundance_for_colors <- tibble(
  MAGs = taxa_names(PAMIR_MAGs_tree_for_colors),
  TotalAbundance = taxa_sums(PAMIR_MAGs_tree_for_colors)
)

phyla_names <- tax_df_for_colors %>%
  left_join(otu_abundance_for_colors, by = "MAGs") %>%
  group_by(Phylum) %>%
  summarise(TotalReads = sum(TotalAbundance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(TotalReads)) %>%
  # Color the dominant phyla explicitly; collapse rarer phyla to "Other" to
  # avoid an unreadable legend.
  slice_head(n = 17) %>%
  pull(Phylum)

tree_tax_plot_df <- tree_tax_base_df %>%
  mutate(
    Simplified_Phylum = if_else(Phylum %in% phyla_names, Phylum, "Other"),
    Simplified_Phylum = factor(Simplified_Phylum, levels = c(phyla_names, "Other"))
  )

phylum_colors <- setNames(
  colorRampPalette(brewer.pal(12, "Paired"))(length(phyla_names)),
  phyla_names
)
phylum_colors["Other"] <- "#D3D3D3"

tree_grouped <- groupOTU(
  plot_tree,
  split(tree_tax_plot_df$MAGs, tree_tax_plot_df$Simplified_Phylum)
)

tree_layout_df <- ggtree(tree_grouped, layout = "rectangular")$data

# Convert ggtree coordinates to ordinary ggplot segments so the tree and the
# feature heatmap can share the same MAG order on the x-axis.
tree_segment_df <- tree_layout_df %>%
  filter(!is.na(parent)) %>%
  left_join(
    tree_layout_df %>%
      select(parent = node, parent_x = x, parent_y = y),
    by = "parent"
  ) %>%
  rename(child_x = x, child_y = y) %>%
  mutate(
    group = as.character(group),
    group = if_else(group %in% c(phyla_names, "Other"), group, "Other"),
    group = factor(group, levels = c(phyla_names, "Other"))
  )

tree_vertical_df <- tree_segment_df %>%
  transmute(x = child_y, xend = child_y, y = parent_x, yend = child_x, group)

tree_horizontal_df <- tree_segment_df %>%
  transmute(x = parent_y, xend = child_y, y = parent_x, yend = parent_x, group)

tip_layout_df <- tree_layout_df %>%
  filter(isTip) %>%
  transmute(MAGs = label, tip_x = y) %>%
  left_join(
    tree_tax_plot_df %>%
      select(MAGs, Simplified_Phylum),
    by = "MAGs"
  )

stopifnot(
  "Tree tips must have one unique horizontal coordinate each" =
    nrow(tip_layout_df) == length(tree_tips) &&
      n_distinct(tip_layout_df$tip_x) == length(tree_tips),
  "Tree tip coordinates must span the complete MAG order" =
    isTRUE(all.equal(sort(tip_layout_df$tip_x), seq_along(tree_tips)))
)

phylum_block_df <- tip_layout_df %>%
  filter(Simplified_Phylum != "Other") %>%
  group_by(Simplified_Phylum) %>%
  summarise(
    xmin = min(tip_x) - 0.5,
    xmax = max(tip_x) + 0.5,
    .groups = "drop"
  )

tree_ymax <- max(tree_layout_df$x, na.rm = TRUE)
tree_ymin <- min(tree_layout_df$x, na.rm = TRUE)

feature_label_df <- significant_features %>%
  arrange(feature_order) %>%
  transmute(
    feature_type,
    feature_id,
    GSEA_direction = case_when(
      expected_direction == "GI-positive" ~ "High-GI enriched",
      expected_direction == "GI-negative" ~ "Low-GI enriched",
      TRUE ~ NA_character_
    ),
    feature_display = case_when(
      expected_direction == "GI-positive" ~ paste0(feature_id, " + ", KEGG_name),
      expected_direction == "GI-negative" ~ paste0(feature_id, " - ", KEGG_name),
      TRUE ~ paste0(feature_id, " ", KEGG_name)
    )
  )

feature_display_levels <- rev(feature_label_df$feature_display)

feature_direction_df <- feature_label_df %>%
  transmute(
    feature_display = factor(feature_display, levels = feature_display_levels),
    GSEA_direction = factor(
      GSEA_direction,
      levels = c("High-GI enriched", "Low-GI enriched")
    )
  )

p_feature_labels <- ggplot(
  feature_direction_df,
  aes(x = 1, y = feature_display, label = feature_display)
) +
  geom_text(hjust = 1, size = 2.2, color = "black") +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_discrete(drop = FALSE) +
  coord_cartesian(clip = "off") +
  theme_void() +
  theme(plot.margin = margin(t = 2, r = 2, b = 10, l = 10))

feature_heatmap_df <- mag_feature_scores %>%
  filter(MAGs %in% tree_tips) %>%
  left_join(feature_label_df, by = c("feature_type", "feature_id")) %>%
  filter(!is.na(feature_display)) %>%
  left_join(
    tip_layout_df %>%
      select(MAGs, tip_x),
    by = "MAGs"
  ) %>%
  mutate(
    feature_display = factor(feature_display, levels = feature_display_levels),
    # Color intensity represents core-enrichment KO coverage only. The expected
    # GI direction remains encoded by the + or - sign in each feature label.
    core_enrichment_coverage_percent = 100 * core_enrichment_coverage
  )

stopifnot(
  "All heatmap MAGs must map to a tree-tip coordinate" =
    !anyNA(feature_heatmap_df$tip_x)
)

p_tree_top <- ggplot() +
  geom_rect(
    data = phylum_block_df,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = Simplified_Phylum),
    alpha = 0.50,
    show.legend = FALSE
  ) +
  geom_segment(
    data = tree_horizontal_df,
    aes(x = x, xend = xend, y = y, yend = yend, color = group),
    linewidth = 0.25
  ) +
  geom_segment(
    data = tree_vertical_df,
    aes(x = x, xend = xend, y = y, yend = yend, color = group),
    linewidth = 0.25
  ) +
  scale_color_manual(
    values = phylum_colors,
    name = "Phylum",
    guide = guide_legend(override.aes = list(linewidth = 1.2), order = 1)
  ) +
  scale_fill_manual(
    values = phylum_colors,
    guide = "none"
  ) +
  scale_x_continuous(limits = c(0.5, length(tree_tips) + 0.5), expand = c(0, 0)) +
  scale_y_reverse(limits = c(tree_ymax * 1.03, tree_ymin), expand = c(0, 0)) +
  labs(
    title = "KEGG Core-enrichment KO Coverage Across the MAG Phylogeny",
    subtitle = "Heatmap columns follow the horizontal MAG order"
  ) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.key.height = grid::unit(0.35, "lines"),
    legend.spacing.y = grid::unit(0.05, "lines"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    plot.margin = margin(t = 10, r = 10, b = 2, l = 2)
  )

p_feature_direction <- ggplot(
  feature_direction_df,
  aes(x = 1, y = feature_display, fill = GSEA_direction)
) +
  geom_tile(width = 0.82, height = 0.92) +
  scale_fill_manual(
    values = c(
      "High-GI enriched" = "#D53E4F",
      "Low-GI enriched" = "#3288BD"
    ),
    drop = FALSE,
    name = "GSEA enrichment direction",
    guide = guide_legend(
      title.position = "top",
      order = 2
    )
  ) +
  scale_x_continuous(limits = c(0.5, 1.5), expand = c(0, 0)) +
  scale_y_discrete(drop = FALSE) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 9) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right",
    plot.margin = margin(t = 2, r = 2, b = 10, l = 2)
  )

p_feature_heatmap <- ggplot(
  feature_heatmap_df,
  aes(x = tip_x, y = feature_display, fill = core_enrichment_coverage_percent)
) +
  geom_tile(width = 1, height = 0.92) +
  scale_fill_gradient(
    low = "white",
    high = "grey25",
    limits = c(0, 100),
    name = "Core-enrichment KO coverage (%)",
    guide = guide_colorbar(
      barwidth = 0.6,
      barheight = 6,
      title.position = "top",
      frame.colour = "black",
      frame.linewidth = 0.3,
      order = 3
    )
  ) +
  scale_x_continuous(limits = c(0.5, length(tree_tips) + 0.5), expand = c(0, 0)) +
  scale_y_discrete(drop = FALSE) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 9) +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right",
    plot.margin = margin(t = 2, r = 10, b = 10, l = 2)
  )

feature_tree_layout <- c(
  area(t = 1, l = 1, b = 1, r = 1),
  area(t = 1, l = 2, b = 1, r = 2),
  area(t = 1, l = 3, b = 1, r = 3),
  area(t = 2, l = 1, b = 2, r = 1),
  area(t = 2, l = 2, b = 2, r = 2),
  area(t = 2, l = 3, b = 2, r = 3)
)

p_feature_tree <-
  plot_spacer() + plot_spacer() + p_tree_top +
  p_feature_labels + p_feature_direction + p_feature_heatmap +
  plot_layout(
    design = feature_tree_layout,
    widths = c(5.3, 0.25, 6.95),
    heights = c(0.23, 1.25),
    guides = "collect"
  ) &
  theme(legend.position = "right")

ggsave(
  file.path(out_dir, "KEGG_core_enrichment_KO_coverage_tree.pdf"),
  p_feature_tree,
  width = 12.5,
  height = 8,
  units = "in",
  limitsize = FALSE
)

# --- 7. Summarize feature distribution at phylum level ---

phylum_feature_summary <- mag_feature_scores %>%
  # Summarize how broadly each core-enrichment KO set is distributed within each
  # lineage before any trend modeling.
  group_by(
    feature_type,
    feature_id,
    KEGG_name,
    Broad_functional_group,
    feature_order,
    expected_direction,
    Phylum
  ) %>%
  summarise(
    n_MAGs = n(),
    n_MAGs_with_core_enrichment_coverage = sum(core_enrichment_coverage > 0, na.rm = TRUE),
    fraction_MAGs_with_core_enrichment_coverage = n_MAGs_with_core_enrichment_coverage / n_MAGs,
    mean_core_enrichment_coverage = mean(core_enrichment_coverage, na.rm = TRUE),
    median_core_enrichment_coverage = median(core_enrichment_coverage, na.rm = TRUE),
    mean_lfc_GI = mean(lfc_GI, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(
    feature_order,
    desc(fraction_MAGs_with_core_enrichment_coverage),
    desc(mean_core_enrichment_coverage)
  )

write_csv(
  phylum_feature_summary,
  file.path(out_dir, "phylum_core_enrichment_coverage_summary.csv")
)

# --- 8. Within-phylum trend analysis for all features ---

# Testability criteria are applied separately to each phylum-feature
# combination. The quasi-binomial GLM has four parameters: the intercept,
# lfc_GI, MAG completeness, and log10 genome size. Sample-size support is
# assessed through the fitted model's residual degrees of freedom rather than
# through a direct MAG-count cutoff. Requiring at least four residual degrees
# of freedom avoids saturated or nearly saturated fits and generally requires
# about eight complete MAG observations. Requiring at least four distinct
# core-enrichment coverage values prevents inference from combinations with
# insufficient within-phylum response variation.
min_model_df_residual <- 4L
min_distinct_core_enrichment_coverage <- 4L

feature_key_cols <- c(
  "feature_type",
  "feature_id",
  "KEGG_name",
  "Broad_functional_group",
  "feature_order",
  "expected_direction"
)

# For each KEGG feature within each phylum, test whether core-enrichment KO
# coverage changes with GI preference after controlling for MAG completeness and
# genome size. The response is the number of detected core-enrichment KOs out of
# the total core-enrichment KOs, modeled with a quasi-binomial GLM to retain the
# natural count/proportion structure and allow overdispersion:
#   cbind(n_present_kos, n_total_kos - n_present_kos) ~
#     lfc_GI + MAG_completeness + log10_genome_size
# Raw Spearman correlation is retained as an unadjusted, nonparametric direction
# check.
spearman_or_na <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]

  if (length(x) < 3 || n_distinct(x) <= 1 || n_distinct(y) <= 1) {
    return(NA_real_)
  }

  suppressWarnings(cor(x, y, method = "spearman"))
}

fit_within_phylum_glm <- function(dat) {
  fit <- tryCatch(
    glm(
      cbind(n_present_kos, n_total_kos - n_present_kos) ~
        lfc_GI + MAG_completeness + log10_genome_size,
      family = quasibinomial(link = "logit"),
      data = dat
    ),
    error = function(e) NULL
  )

  fit_summary <- if (is.null(fit)) NULL else summary(fit)

  if (is.null(fit_summary) || !"lfc_GI" %in% rownames(fit_summary$coefficients)) {
    return(tibble(
      lfc_GI_logit_estimate = NA_real_,
      lfc_GI_logit_std_error = NA_real_,
      lfc_GI_logit_t = NA_real_,
      lfc_GI_logit_p = NA_real_,
      model_dispersion = NA_real_,
      model_converged = FALSE,
      model_df_residual = NA_real_
    ))
  }

  coef_tab <- fit_summary$coefficients
  tibble(
    lfc_GI_logit_estimate = coef_tab["lfc_GI", "Estimate"],
    lfc_GI_logit_std_error = coef_tab["lfc_GI", "Std. Error"],
    lfc_GI_logit_t = coef_tab["lfc_GI", "t value"],
    lfc_GI_logit_p = coef_tab["lfc_GI", "Pr(>|t|)"],
    model_dispersion = fit_summary$dispersion,
    model_converged = isTRUE(fit$converged),
    model_df_residual = df.residual(fit)
  )
}

min_finite_or_na <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) {
    return(NA_real_)
  }
  min(x)
}

feature_phylum_trend_df <- mag_feature_scores %>%
  group_by(across(all_of(c(feature_key_cols, "Phylum")))) %>%
  summarise(
    n_MAGs_phylum = n(),
    n_complete_model_cases = sum(
      is.finite(core_enrichment_coverage) &
        is.finite(lfc_GI) &
        is.finite(MAG_completeness) &
        is.finite(log10_genome_size)
    ),
    n_distinct_lfc_GI = n_distinct(lfc_GI[is.finite(lfc_GI)]),
    n_distinct_core_enrichment_coverage = n_distinct(
      core_enrichment_coverage[is.finite(core_enrichment_coverage)]
    ),
    n_MAGs_with_core_enrichment_coverage = sum(
      core_enrichment_coverage > 0,
      na.rm = TRUE
    ),
    n_MAGs_without_core_enrichment_coverage = sum(
      core_enrichment_coverage == 0,
      na.rm = TRUE
    ),
    fraction_MAGs_with_core_enrichment_coverage =
      n_MAGs_with_core_enrichment_coverage / n_MAGs_phylum,
    mean_core_enrichment_coverage = mean(core_enrichment_coverage, na.rm = TRUE),
    raw_coverage_lfc_GI_spearman = spearman_or_na(core_enrichment_coverage, lfc_GI),
    model = list(fit_within_phylum_glm(pick(everything()))),
    .groups = "drop"
  ) %>%
  unnest(model) %>%
  mutate(
    # Although a model is attempted for every phylum-feature combination, a
    # combination is considered testable only if the GLM converges, retains
    # sufficient residual degrees of freedom, contains adequate within-phylum
    # coverage variation, and yields a finite lfc_GI coefficient. Estimates
    # from combinations failing any criterion are retained for diagnostics but
    # are not interpreted or included in multiple-testing correction.
    testable_phylum = case_when(
      !model_converged ~ FALSE,
      !is.finite(model_df_residual) |
        model_df_residual < min_model_df_residual ~ FALSE,
      n_distinct_core_enrichment_coverage <
        min_distinct_core_enrichment_coverage ~ FALSE,
      !is.finite(lfc_GI_logit_estimate) ~ FALSE,
      TRUE ~ TRUE
    ),
    trend_testability_reason = case_when(
      !model_converged ~ "glm_not_converged",
      !is.finite(model_df_residual) |
        model_df_residual < min_model_df_residual ~ "no_model_residual_df",
      n_distinct_core_enrichment_coverage <
        min_distinct_core_enrichment_coverage ~
        "too_few_distinct_coverage_values",
      !is.finite(lfc_GI_logit_estimate) ~ "undefined_glm",
      TRUE ~ "testable"
    ),
    expected_direction_supported = case_when(
      testable_phylum & expected_direction == "GI-positive" ~ lfc_GI_logit_estimate > 0,
      testable_phylum & expected_direction == "GI-negative" ~ lfc_GI_logit_estimate < 0,
      TRUE ~ NA
    )
  ) %>%
  group_by(across(all_of(feature_key_cols))) %>%
  mutate(
    # Within each KEGG feature, control FDR across only the testable phylum-
    # feature combinations. Untestable combinations are assigned NA before
    # Benjamini-Hochberg adjustment and therefore do not enter the correction.
    lfc_GI_logit_fdr = p.adjust(
      if_else(testable_phylum, lfc_GI_logit_p, NA_real_),
      method = "BH"
    )
  ) %>%
  ungroup() %>%
  arrange(feature_order, Phylum)

# Exclude all photosynthesis-related features from the within-phylum trend
# summary. ko00196 and M00145 were already removed from every projection output
# above; the remaining four stay in the projection and phylogenetic
# visualization but are omitted from this cross-lineage summary.
photosynthesis_features_excluded_from_trend_summary <- c(
  "ko00195", # Photosynthesis.
  "ko00196", # Photosynthetic antenna proteins.
  "M00145",  # Cyanobacterial-type NAD(P)H:quinone oxidoreductase.
  "M00161",  # Photosystem II.
  "M00162",  # Cytochrome b6f complex.
  "M00163"   # Photosystem I.
)

min_testable_phyla_for_cross_lineage_summary <- 3L
within_phylum_fdr_threshold <- 0.05

feature_within_phylum_trend_summary <- feature_phylum_trend_df %>%
  # Report only broadly distributed feature trends. Two photosynthesis features
  # were excluded globally; the remaining four are retained upstream but
  # omitted from this table.
  filter(!feature_id %in% photosynthesis_features_excluded_from_trend_summary) %>%
  group_by(across(all_of(feature_key_cols))) %>%
  summarise(
    n_testable_phyla = sum(testable_phylum, na.rm = TRUE),
    n_direction_concordant_phyla = sum(expected_direction_supported, na.rm = TRUE),
    direction_concordant_phyla = paste(
      sort(unique(Phylum[expected_direction_supported %in% TRUE])),
      collapse = "; "
    ),
    n_fdr_significant_concordant_phyla = sum(
      expected_direction_supported %in% TRUE &
        is.finite(lfc_GI_logit_fdr) &
        lfc_GI_logit_fdr < within_phylum_fdr_threshold
    ),
    fdr_significant_concordant_phyla = paste(
      sort(unique(Phylum[
        expected_direction_supported %in% TRUE &
          is.finite(lfc_GI_logit_fdr) &
          lfc_GI_logit_fdr < within_phylum_fdr_threshold
      ])),
      collapse = "; "
    ),
    n_fdr_significant_opposite_phyla = sum(
      expected_direction_supported %in% FALSE &
        is.finite(lfc_GI_logit_fdr) &
        lfc_GI_logit_fdr < within_phylum_fdr_threshold
    ),
    fdr_significant_opposite_phyla = paste(
      sort(unique(Phylum[
        expected_direction_supported %in% FALSE &
          is.finite(lfc_GI_logit_fdr) &
          lfc_GI_logit_fdr < within_phylum_fdr_threshold
      ])),
      collapse = "; "
    ),
    median_lfc_GI_logit_estimate = median(
      lfc_GI_logit_estimate[testable_phylum],
      na.rm = TRUE
    ),
    # This is the strongest direction-concordant phylum-level FDR within the
    # feature, not a feature-level FDR.
    best_direction_concordant_phylum_fdr = min_finite_or_na(
      lfc_GI_logit_fdr[expected_direction_supported %in% TRUE]
    ),
    .groups = "drop"
  ) %>%
  mutate(
    direction_concordance_fraction = if_else(
      n_testable_phyla > 0,
      n_direction_concordant_phyla / n_testable_phyla,
      NA_real_
    ),
    median_lfc_GI_logit_estimate = if_else(
      n_testable_phyla > 0,
      median_lfc_GI_logit_estimate,
      NA_real_
    ),
    best_direction_concordant_phylum_fdr = if_else(
      n_direction_concordant_phyla > 0,
      best_direction_concordant_phylum_fdr,
      NA_real_
    )
  ) %>%
  filter(n_testable_phyla >= min_testable_phyla_for_cross_lineage_summary) %>%
  arrange(feature_order)

feature_within_phylum_trend_extended_table <- feature_within_phylum_trend_summary %>%
  transmute(
    KEGG_ID = feature_id,
    Feature_type = feature_type,
    KEGG_name,
    Broad_functional_group,
    Expected_GI_direction = recode(
      expected_direction,
      "GI-positive" = "High-GI enriched",
      "GI-negative" = "Low-GI enriched"
    ),
    Testable_phyla = n_testable_phyla,
    Direction_concordant_phyla = n_direction_concordant_phyla,
    Direction_concordance_percent = round(100 * direction_concordance_fraction, 1),
    FDR_significant_concordant_phyla = n_fdr_significant_concordant_phyla,
    FDR_significant_opposite_phyla = n_fdr_significant_opposite_phyla,
    Median_lfc_GI_logit_slope_all_testable_phyla =
      signif(median_lfc_GI_logit_estimate, 3),
    Best_direction_concordant_phylum_FDR =
      signif(best_direction_concordant_phylum_fdr, 3)
  )

write_csv(
  feature_phylum_trend_df,
  file.path(out_dir, "feature_within_phylum_trend_details.csv")
)

write_csv(
  feature_within_phylum_trend_summary,
  file.path(out_dir, "feature_within_phylum_trend_summary.csv")
)

write_csv(
  feature_within_phylum_trend_extended_table,
  file.path(out_dir, "feature_within_phylum_trend_extended_data_table.csv")
)
