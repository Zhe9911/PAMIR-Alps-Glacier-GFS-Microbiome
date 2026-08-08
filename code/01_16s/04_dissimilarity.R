# Dissimilarity analysis using Bray-Curtis and weighted UniFrac
# Manuscript outputs: Figure 1c and sediment-ice inputs for Figures 2c-d.

rm(list = ls())
gc()
graphics.off()

input_path <- file.path(
  "results", "01_16s", "intermediate", "PAMIR_16S_final.rds"
)
tree_path <- file.path("data", "processed", "16s", "dna-sequences.tree")
result_dir <- file.path("results", "01_16s", "dissimilarity")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

# --- Setup ---
library(phyloseq)
library(ggplot2)
library(vegan)
library(tidyverse)
library(ape)
library(phytools)

# --- Load input data ---

ps <- readRDS(input_path)
print(ps)

sample_metadata_df <- data.frame(sample_data(ps))
sample_metadata_df$Source <- as.character(sample_metadata_df$Source)
sample_names <- sample_names(ps)

# --- Functions ---

get_aligned_phyloseq_with_rooted_tree <- function(physeq_obj, tree_path) {
  raw_tree <- read.tree(tree_path)
  shared_taxa <- intersect(taxa_names(physeq_obj), raw_tree$tip.label)

  physeq_aligned <- prune_taxa(shared_taxa, physeq_obj)
  rooted_tree <- midpoint.root(keep.tip(raw_tree, shared_taxa))
  phy_tree(physeq_aligned) <- rooted_tree

  physeq_aligned
}

matrix_to_unique_pair_df <- function(dist_matrix, value_col = "Dissimilarity") {
  pair_idx <- which(upper.tri(dist_matrix) & !is.na(dist_matrix), arr.ind = TRUE)

  out <- data.frame(
    Sample1 = rownames(dist_matrix)[pair_idx[, 1]],
    Sample2 = colnames(dist_matrix)[pair_idx[, 2]],
    stringsAsFactors = FALSE
  )
  out[[value_col]] <- dist_matrix[pair_idx]
  out
}

# Transform counts to relative abundance.
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))

# --- Build dissimilarity matrices ---

# Calculate the Bray-Curtis dissimilarity matrix.
bray_dist <- distance(ps_rel, method = "bray")

bray_dist_table <- as.matrix(bray_dist)

# Align taxa to the phylogenetic tree before calculating weighted UniFrac.
ps_rel_phylo <- get_aligned_phyloseq_with_rooted_tree(ps_rel, tree_path)
print(ps_rel_phylo)

wunifrac_dist <- distance(ps_rel_phylo, method = "wunifrac")

wunifrac_dist_table <- as.matrix(wunifrac_dist)

# --- Extract pairwise distances for grouped comparisons ---

extract_dissimilarity <- function(dist_matrix, metadata, group_col) {
  # Convert the distance matrix to long format.
  dist_df <- as.data.frame(as.table(dist_matrix))
  colnames(dist_df) <- c("Sample1", "Sample2", "Distance")

  # Remove comparisons in which both entries have the same Sample_ID.
  dist_df <- dist_df[dist_df$Sample1 != dist_df$Sample2, ]

  # Add metadata for both samples.
  dist_df$Source1 <- metadata[match(dist_df$Sample1, metadata$Sample_ID), group_col]
  dist_df$Source2 <- metadata[match(dist_df$Sample2, metadata$Sample_ID), group_col]
  dist_df$Region1 <- metadata[match(dist_df$Sample1, metadata$Sample_ID), "Region"]
  dist_df$Region2 <- metadata[match(dist_df$Sample2, metadata$Sample_ID), "Region"]

  # Create comparison categories.
  dist_df$Comparison <- paste(
    pmin(dist_df$Source1, dist_df$Source2),
    pmax(dist_df$Source1, dist_df$Source2),
    sep = " vs "
  )

  # Retain one record per unordered sample pair.
  dist_df$pair_id <- paste(
    pmin(as.character(dist_df$Sample1), as.character(dist_df$Sample2)),
    pmax(as.character(dist_df$Sample1), as.character(dist_df$Sample2)),
    sep = "_"
  )
  dist_df <- dist_df[!duplicated(dist_df$pair_id), ]

  return(dist_df)
}

# --- Bray-Curtis distances for Figure 1c ---

# Extract beta-diversity data.
Bray_Curtis_data <- extract_dissimilarity(bray_dist_table, sample_metadata_df, "Source")

# Create within-group and between-group categories.
Bray_Curtis_data$ComparisonType <- ifelse(
  Bray_Curtis_data$Source1 == Bray_Curtis_data$Source2,
  "Within-group",
  "Between-group"
)

# Create detailed comparison labels.
Bray_Curtis_data$DetailedComparison <- ifelse(
  Bray_Curtis_data$Source1 == Bray_Curtis_data$Source2,
  paste("Within", Bray_Curtis_data$Source1),
  Bray_Curtis_data$Comparison
)

# Set the factor-level order.
Bray_Curtis_data$DetailedComparison <- factor(
  Bray_Curtis_data$DetailedComparison,
  levels = c(
    "Within ice", "Within water", "Within sed",
    "ice vs water", "ice vs sed", "sed vs water"
  )
)

# Define the planned pairwise comparisons.
my_comparisons <- list(
  c("Within ice", "Within sed"),
  c("Within ice", "Within water"),
  c("Within sed", "Within water"),
  c("ice vs sed", "ice vs water"),
  c("ice vs sed", "sed vs water"),
  c("ice vs water", "sed vs water")
)

library(ggpubr)

# Run sample-label permutation tests on pairwise distance contrasts.
# This keeps the comparison target as mean Bray-Curtis distance differences,
# while avoiding a Wilcoxon test on non-independent pairwise distances.
run_distance_contrast_permutation <- function(
    dist_matrix,
    metadata,
    group_col,
    comparisons,
    detailed_levels,
    n_perm = 9999,
    seed = 123) {
  pair_idx <- which(upper.tri(dist_matrix) & !is.na(dist_matrix), arr.ind = TRUE)
  pair_data <- data.frame(
    Sample1 = rownames(dist_matrix)[pair_idx[, 1]],
    Sample2 = colnames(dist_matrix)[pair_idx[, 2]],
    Distance = dist_matrix[pair_idx],
    stringsAsFactors = FALSE
  )

  group_by_sample <- metadata[[group_col]]
  names(group_by_sample) <- metadata$Sample_ID

  mean_distance_by_group <- function(group_vector) {
    source1 <- group_vector[pair_data$Sample1]
    source2 <- group_vector[pair_data$Sample2]
    comparison <- paste(pmin(source1, source2), pmax(source1, source2), sep = " vs ")
    detailed_comparison <- ifelse(source1 == source2, paste("Within", source1), comparison)
    detailed_comparison <- factor(detailed_comparison, levels = detailed_levels)
    tapply(pair_data$Distance, detailed_comparison, mean, na.rm = TRUE)
  }

  set.seed(seed)
  observed_means <- mean_distance_by_group(group_by_sample)

  perm_diffs <- replicate(n_perm, {
    group_perm <- group_by_sample
    group_perm[] <- sample(group_perm)
    perm_means <- mean_distance_by_group(group_perm)
    vapply(comparisons, function(comp) {
      perm_means[comp[1]] - perm_means[comp[2]]
    }, numeric(1))
  })

  observed_diffs <- vapply(comparisons, function(comp) {
    observed_means[comp[1]] - observed_means[comp[2]]
  }, numeric(1))

  p_values <- vapply(seq_along(comparisons), function(i) {
    (sum(abs(perm_diffs[i, ]) >= abs(observed_diffs[i]), na.rm = TRUE) + 1) /
      (sum(!is.na(perm_diffs[i, ])) + 1)
  }, numeric(1))

  out <- data.frame(
    group1 = vapply(comparisons, `[`, character(1), 1),
    group2 = vapply(comparisons, `[`, character(1), 2),
    observed_diff = observed_diffs,
    p = p_values,
    stringsAsFactors = FALSE
  )
  out$contrast_family <- ifelse(
    grepl("^Within", out$group1) & grepl("^Within", out$group2),
    "within_group",
    "between_group"
  )
  out$p.adj <- ave(
    out$p,
    out$contrast_family,
    FUN = function(x) p.adjust(x, method = "BH")
  )
  out$p.adj.signif <- cut(
    out$p.adj,
    breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
    labels = c("***", "**", "*", "ns")
  )
  out$y.position <- seq(1.03, by = 0.04, length.out = nrow(out))
  out
}

stat.test <- run_distance_contrast_permutation(
  bray_dist_table,
  sample_metadata_df,
  "Source",
  my_comparisons,
  levels(Bray_Curtis_data$DetailedComparison)
)

print(stat.test)
write.csv(
  stat.test,
  file.path(result_dir, "figure_1c_permutation_contrasts.csv"),
  row.names = FALSE
)

# Plot Bray-Curtis comparisons across habitats.
p3 <- ggplot(
  Bray_Curtis_data,
  aes(x = DetailedComparison, y = Distance, fill = ComparisonType)
) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.2, size = 0.2) +
  # Add adjusted significance labels to the plot.
  stat_pvalue_manual(
    stat.test,
    label = "p.adj.signif",
    tip.length = 0.01,
    hide.ns = TRUE,
    step.increase = 0.04,
    bracket.nudge.y = 0.03
  ) +
  labs(
    title = "Beta Diversity Dissimilarity Comparison",
    subtitle = "Sample-label permutation contrast test with BH correction (FDR)",
    y = "Bray-Curtis Dissimilarity",
    x = NULL,
    fill = "Group Type"
  ) +
  scale_y_continuous(
    breaks = seq(0, 1, 0.25),
    limits = c(0, 1.2),
    expand = c(0, 0)
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 9),
    legend.position = "top"
  ) +
  scale_fill_manual(
    values = c("Within-group" = "#c25e73", "Between-group" = "#6b81a6")
  )

ggsave(
  file.path(result_dir, "figure_1c_bray_curtis_dissimilarity.pdf"),
  p3,
  width = 8,
  height = 8,
  dpi = 300
)

# --- Bray-Curtis between sample types within the same glacier ---

filtered_bray_dist <- matrix(
  NA,
  nrow = nrow(bray_dist_table),
  ncol = ncol(bray_dist_table),
  dimnames = list(
    rownames(bray_dist_table),
    colnames(bray_dist_table)
  )
)

# Retain pairs from the same glacier but different sample types.
for (i in seq_along(sample_names)) {
  for (j in seq_along(sample_names)) {
    if (i == j) next

    sample_i <- sample_names[i]
    sample_j <- sample_names[j]

    gl_name_i <- as.character(sample_metadata_df[sample_i, "gl_name"])
    gl_name_j <- as.character(sample_metadata_df[sample_j, "gl_name"])

    sample_type_i <- as.character(sample_metadata_df[sample_i, "Source"])
    sample_type_j <- as.character(sample_metadata_df[sample_j, "Source"])

    if (gl_name_i == gl_name_j && sample_type_i != sample_type_j) {
      filtered_bray_dist[sample_i, sample_j] <- bray_dist_table[sample_i, sample_j]
    }
  }
}

filtered_dist_df <- matrix_to_unique_pair_df(filtered_bray_dist, "Dissimilarity")
filtered_dist_df$Sample1_type <- sample_metadata_df$Source[
  match(filtered_dist_df$Sample1, rownames(sample_metadata_df))
]
filtered_dist_df$Sample2_type <- sample_metadata_df$Source[
  match(filtered_dist_df$Sample2, rownames(sample_metadata_df))
]
filtered_dist_df$gl_name <- sample_metadata_df$gl_name[
  match(filtered_dist_df$Sample1, rownames(sample_metadata_df))
]
filtered_dist_df$Region <- sample_metadata_df$Region[
  match(filtered_dist_df$Sample1, rownames(sample_metadata_df))
]
filtered_dist_df <- filtered_dist_df %>%
  select(Sample1, Sample2, Sample1_type, Sample2_type, gl_name, Region, Dissimilarity)

# Create an order-independent sample-type comparison label.
filtered_dist_df <- filtered_dist_df %>%
  mutate(comparison = ifelse(
    Sample1_type < Sample2_type,
    paste(Sample1_type, Sample2_type, sep = "_vs_"),
    paste(Sample2_type, Sample1_type, sep = "_vs_")
  )) %>%
  filter(comparison == "ice_vs_sed")

write.csv(
  filtered_dist_df,
  file.path(result_dir, "within_glacier_ice_sed_bray_curtis.csv"),
  row.names = FALSE
)

# --- Weighted UniFrac between sample types within the same glacier ---

filtered_wunifrac_dist <- matrix(
  NA,
  nrow = nrow(wunifrac_dist_table),
  ncol = ncol(wunifrac_dist_table),
  dimnames = list(
    rownames(wunifrac_dist_table),
    colnames(wunifrac_dist_table)
  )
)

# Retain pairs from the same glacier but different sample types.
for (i in seq_along(sample_names)) {
  for (j in seq_along(sample_names)) {
    if (i == j) next

    sample_i <- sample_names[i]
    sample_j <- sample_names[j]

    gl_name_i <- as.character(sample_metadata_df[sample_i, "gl_name"])
    gl_name_j <- as.character(sample_metadata_df[sample_j, "gl_name"])

    sample_type_i <- as.character(sample_metadata_df[sample_i, "Source"])
    sample_type_j <- as.character(sample_metadata_df[sample_j, "Source"])

    if (gl_name_i == gl_name_j && sample_type_i != sample_type_j) {
      filtered_wunifrac_dist[sample_i, sample_j] <- wunifrac_dist_table[sample_i, sample_j]
    }
  }
}

filtered_dist_df <- matrix_to_unique_pair_df(filtered_wunifrac_dist, "Dissimilarity")
filtered_dist_df$Sample1_type <- sample_metadata_df$Source[
  match(filtered_dist_df$Sample1, rownames(sample_metadata_df))
]
filtered_dist_df$Sample2_type <- sample_metadata_df$Source[
  match(filtered_dist_df$Sample2, rownames(sample_metadata_df))
]
filtered_dist_df$gl_name <- sample_metadata_df$gl_name[
  match(filtered_dist_df$Sample1, rownames(sample_metadata_df))
]
filtered_dist_df$Region <- sample_metadata_df$Region[
  match(filtered_dist_df$Sample1, rownames(sample_metadata_df))
]
filtered_dist_df <- filtered_dist_df %>%
  select(Sample1, Sample2, Sample1_type, Sample2_type, gl_name, Region, Dissimilarity)

# Create an order-independent sample-type comparison label.
filtered_dist_df <- filtered_dist_df %>%
  mutate(comparison = ifelse(
    Sample1_type < Sample2_type,
    paste(Sample1_type, Sample2_type, sep = "_vs_"),
    paste(Sample2_type, Sample1_type, sep = "_vs_")
  )) %>%
  filter(comparison == "ice_vs_sed")

write.csv(
  filtered_dist_df,
  file.path(result_dir, "within_glacier_ice_sed_weighted_unifrac.csv"),
  row.names = FALSE
)
