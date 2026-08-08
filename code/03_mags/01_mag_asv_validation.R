# MAG-ASV concordance validation
# Manuscript outputs: Extended Data Figures 3-4.

rm(list = ls())
gc()
graphics.off()

out_dir <- file.path("results", "03_mags", "asv_mag_validation")
intermediate_dir <- file.path("results", "03_mags", "intermediate")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(intermediate_dir, recursive = TRUE, showWarnings = FALSE)

# --- Setup ---

library(tidyverse)
library(phyloseq)
library(vegan)
library(ggplot2)

# --- Load and clean CoverM relative-abundance data ---

mags_rela_ab_raw <- read.table(
  file.path("data", "processed", "metagenomics", "coverm_mags_relative_abundance.tsv"),
  header = TRUE,
  sep = "\t",
  quote = "",
  fill = TRUE,
  stringsAsFactors = FALSE
)

mags_rela_ab <- mags_rela_ab_raw %>%
  filter(Genome != "unmapped")

rownames(mags_rela_ab) <- mags_rela_ab$Genome
mags_rela_ab$Genome <- NULL

# Clean sample names.
colnames(mags_rela_ab) <- gsub(
  "\\_1.fq\\.gz\\.Relative\\.Abundance.*$", "",
  gsub("^cat_mags\\.fa\\.", "", colnames(mags_rela_ab))
)

mags_rela_ab <- as.data.frame(mags_rela_ab)

# --- Load and clean metadata ---

metadata_raw <- read_csv(
  file.path("data", "processed", "metagenomics", "metadata_metagenomics.csv"),
  show_col_types = FALSE
) %>%
  column_to_rownames(var = colnames(.)[1]) %>%
  filter(!is.na(Source)) %>%
  mutate(Source = factor(Source)) %>%
  select(where(~ !all(is.na(.))))

# Calculate sediment GI tertiles.
quantile_breaks <- quantile(
  metadata_raw$GI[metadata_raw$Source == "sed" & !is.na(metadata_raw$GI)],
  probs = seq(0, 1, by = 1 / 3),
  na.rm = TRUE
)

print(quantile_breaks)

# Assign Group_GI from sediment GI tertiles and existing habitat labels.
metadata <- metadata_raw %>%
  mutate(
    Group_GI = case_when(
      Source == "ice" ~ "ice",
      Source %in% c("sed", "water") & !is.na(GI) ~ paste(
        as.character(
          cut(
            GI,
            breaks = quantile_breaks,
            include.lowest = TRUE,
            labels = c("Low_Glacial", "Mid_Glacial", "High_Glacial")
          )
        ),
        Source,
        sep = "_"
      ),
      !is.na(GI) ~ as.character(
        cut(
          GI,
          breaks = quantile_breaks,
          include.lowest = TRUE,
          labels = c("Low_Glacial", "Mid_Glacial", "High_Glacial")
        )
      ),
      TRUE ~ NA_character_
    )
  )

stats_quant <- metadata %>%
  group_by(Group_GI) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(Method = "Quantile")

print(stats_quant)

# --- Load and clean GTDB-Tk taxonomy ---

taxonomy_MAGs_raw <- read.delim(
  file.path("data", "processed", "metagenomics", "gtdbtk.bac120.summary.tsv"),
  header = TRUE,
  quote = ""
) %>%
  select(user_genome, classification) %>%
  separate(
    classification,
    c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
    sep = ";"
  ) %>%
  rename(MAGs = user_genome)

rownames(taxonomy_MAGs_raw) <- taxonomy_MAGs_raw$MAGs
taxonomy_MAGs_raw$MAGs <- NULL

taxonomy_MAGs <- apply(
  taxonomy_MAGs_raw,
  2,
  function(x) {
    # Remove GTDB rank prefixes such as d__, p__, and c__.
    x <- gsub("^[a-z]__", "", x)
    x[x == "" | x == ""] <- "unknown"
    x[x == "s__"] <- "unknown"
    x
  }
)

# Relabel major Pseudomonadota classes at phylum level for plotting.
taxonomy_MAGs <- data.frame(taxonomy_MAGs) %>%
  mutate(
    Phylum = case_when(
      Phylum == "Pseudomonadota" &
        Class %in% c("Alphaproteobacteria", "Gammaproteobacteria") ~ Class,
      TRUE ~ Phylum
    )
  )
table(taxonomy_MAGs$Phylum)

taxonomy_MAGs <- as.matrix(taxonomy_MAGs)

# --- Build the MAG phyloseq object ---

MAGs_table <- otu_table(mags_rela_ab, taxa_are_rows = TRUE)
sample_data_MAGs <- sample_data(metadata)
tax_MAGs <- tax_table(taxonomy_MAGs)

PAMIR_MAGs <- merge_phyloseq(
  MAGs_table,
  sample_data_MAGs,
  tax_MAGs
)

print(PAMIR_MAGs)

# Renormalize to relative abundance (%).
PAMIR_MAGs_rela <- transform_sample_counts(
  PAMIR_MAGs,
  function(x) x / sum(x) * 100
)

# Define plotting order for all sample groups, including water GI classes.
Group_GI_order <- c(
  "ice",
  "High_Glacial_sed",
  "Mid_Glacial_sed",
  "Low_Glacial_sed",
  "High_Glacial_water",
  "Mid_Glacial_water",
  "Low_Glacial_water"
)

sample_data(PAMIR_MAGs_rela)$Group_GI <- factor(
  sample_data(PAMIR_MAGs_rela)$Group_GI,
  levels = Group_GI_order
)

saveRDS(
  PAMIR_MAGs_rela,
  file = file.path(intermediate_dir, "PAMIR_MAGs_rela.rds")
)

# --- Order-level ASV-MAG concordance: functions ---

# Apply only format-level standardization; no taxonomic synonyms are merged.
harmonize_order_names <- function(x) {
  x <- as.character(x)
  x <- gsub("^[a-z]__", "", x)
  x[is.na(x) | x == "" | x == "NA"] <- "unknown"
  x
}

# Construct a stable identifier for matching 16S and metagenomic samples.
make_order_match_id <- function(meta_df) {
  paste(
    meta_df$Source,
    meta_df$Region,
    meta_df$gl_name,
    round(meta_df$distance_to_glacier_snout, 3),
    round(meta_df$GI, 6),
    sep = "|"
  )
}

# Merge ecological 16S replicates before matching them to metagenomic samples.
prepare_16s_order_phyloseq <- function(rds_path) {
  ps_16s <- readRDS(rds_path)
  sample_data(ps_16s)$Sample_Base <- gsub("_.*", "", sample_names(ps_16s))
  ps_16s_merged <- merge_samples(ps_16s, "Sample_Base")

  original_meta <- as.data.frame(sample_data(ps_16s))
  restored_meta <- original_meta[!duplicated(original_meta$Sample_Base), ]
  rownames(restored_meta) <- restored_meta$Sample_Base
  restored_meta$Sample_ID <- restored_meta$Sample_Base
  restored_meta$Sample_Base <- NULL
  sample_data(ps_16s_merged) <- sample_data(restored_meta)

  if (!taxa_are_rows(ps_16s_merged) && taxa_are_rows(ps_16s)) {
    otu_table(ps_16s_merged) <- otu_table(
      t(otu_table(ps_16s_merged)),
      taxa_are_rows = TRUE
    )
  }

  ps_16s_merged %>%
    subset_samples(Source %in% c("ice", "sed")) %>%
    prune_taxa(taxa_sums(.) > 0, .) %>%
    prune_samples(sample_sums(.) > 0, .) %>%
    transform_sample_counts(function(x) x / sum(x))
}

# Record every original-to-standardized order label for reproducibility.
build_order_nomenclature_audit <- function(ps_obj, dataset) {
  raw_order <- as.character(tax_table(ps_obj)[, "Order"])
  harmonized_order <- harmonize_order_names(raw_order)
  action <- case_when(
    is.na(raw_order) | raw_order == "" | raw_order == "NA" ~
      "missing label standardised to unknown",
    grepl("^[a-z]__", raw_order) ~ "rank prefix removed",
    TRUE ~ "unchanged"
  )

  tibble(
    Dataset = dataset,
    Original_order = if_else(
      is.na(raw_order) | raw_order == "",
      "<missing>",
      raw_order
    ),
    Harmonized_order = harmonized_order,
    Standardization = action
  ) %>%
    distinct()
}

# Aggregate relative abundances at order rank for sample-wise comparison.
aggregate_order_matrix <- function(ps_obj) {
  ps_order <- tax_glom(ps_obj, taxrank = "Order", NArm = FALSE)
  order_names <- harmonize_order_names(tax_table(ps_order)[, "Order"])
  mat <- as(otu_table(ps_order), "matrix")
  if (taxa_are_rows(ps_order))
    mat <- t(mat)
  colnames(mat) <- order_names
  as.data.frame(mat, check.names = FALSE)
}

# --- Order-level concordance: inputs and sample matching ---

asv_object_path <- file.path(
  "results", "01_16s", "intermediate", "PAMIR_16S_final.rds"
)
ps_16s_order <- prepare_16s_order_phyloseq(asv_object_path)
ps_mags_order <- PAMIR_MAGs_rela %>%
  subset_samples(Source %in% c("ice", "sed")) %>%
  prune_taxa(taxa_sums(.) > 0, .) %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  transform_sample_counts(function(x) x / sum(x))

nomenclature_audit <- bind_rows(
  build_order_nomenclature_audit(ps_16s_order, "16S ASVs"),
  build_order_nomenclature_audit(ps_mags_order, "MAGs")
)
write_csv(
  nomenclature_audit,
  file.path(out_dir, "Order_nomenclature_crosswalk.csv")
)

nomenclature_summary <- nomenclature_audit %>%
  group_by(Dataset) %>%
  summarise(
    Unique_original_names = n_distinct(Original_order),
    Unique_harmonized_names = n_distinct(Harmonized_order),
    Names_with_rank_prefix_removed = sum(Standardization == "rank prefix removed"),
    Missing_names_standardised_to_unknown = sum(
      Standardization == "missing label standardised to unknown"
    ),
    .groups = "drop"
  )

meta_16s <- data.frame(sample_data(ps_16s_order)) %>%
  rownames_to_column("Sample_16S") %>%
  mutate(Match_ID = make_order_match_id(.))
meta_mags <- data.frame(sample_data(ps_mags_order)) %>%
  rownames_to_column("Sample_MAG") %>%
  mutate(Match_ID = make_order_match_id(.))
order_sample_map <- inner_join(
  meta_16s %>% select(Sample_16S, Match_ID),
  meta_mags %>% select(Sample_MAG, Match_ID),
  by = "Match_ID"
) %>%
  arrange(Match_ID)
write_csv(
  order_sample_map,
  file.path(out_dir, "Order_ASV_MAG_sample_map.csv")
)

ps_16s_order <- prune_samples(order_sample_map$Sample_16S, ps_16s_order)
ps_mags_order <- prune_samples(order_sample_map$Sample_MAG, ps_mags_order)
sample_names(ps_16s_order) <- order_sample_map$Match_ID[
  match(sample_names(ps_16s_order), order_sample_map$Sample_16S)
]
sample_names(ps_mags_order) <- order_sample_map$Match_ID[
  match(sample_names(ps_mags_order), order_sample_map$Sample_MAG)
]

asv_order_mat <- aggregate_order_matrix(ps_16s_order)[
  order_sample_map$Match_ID,
  ,
  drop = FALSE
]
mag_order_mat <- aggregate_order_matrix(ps_mags_order)[
  order_sample_map$Match_ID,
  ,
  drop = FALSE
]
all_orders <- union(colnames(asv_order_mat), colnames(mag_order_mat))
# Add structural zeros so both matrices contain the same order set.
asv_order_mat[setdiff(all_orders, colnames(asv_order_mat))] <- 0
mag_order_mat[setdiff(all_orders, colnames(mag_order_mat))] <- 0
asv_order_mat <- asv_order_mat[, all_orders, drop = FALSE]
mag_order_mat <- mag_order_mat[, all_orders, drop = FALSE]

# --- Order-level concordance: coverage and statistical tests ---

asv_detected_orders <- names(which(colSums(asv_order_mat) > 0))
mag_detected_orders <- names(which(colSums(mag_order_mat) > 0))
shared_orders <- intersect(asv_detected_orders, mag_detected_orders)
asv_mean_order <- colMeans(asv_order_mat)
mag_mean_order <- colMeans(mag_order_mat)

nomenclature_summary <- bind_rows(
  nomenclature_summary,
  tibble(
    Dataset = "Matched abundance matrices",
    Unique_original_names = NA_integer_,
    Unique_harmonized_names = NA_integer_,
    Names_with_rank_prefix_removed = NA_integer_,
    Missing_names_standardised_to_unknown = NA_integer_
  )
) %>%
  mutate(
    ASV_detected_orders = c(
      rep(NA_integer_, n() - 1),
      length(asv_detected_orders)
    ),
    MAG_detected_orders = c(
      rep(NA_integer_, n() - 1),
      length(mag_detected_orders)
    ),
    Shared_orders_after_harmonization = c(
      rep(NA_integer_, n() - 1),
      length(shared_orders)
    )
  )
write_csv(
  nomenclature_summary,
  file.path(out_dir, "Order_nomenclature_summary.csv")
)

order_coverage <- tibble(
  Metric = c(
    "ASV detected orders",
    "MAG detected orders",
    "Shared orders",
    "ASV-only orders",
    "MAG-only orders",
    "ASV abundance covered by shared orders",
    "MAG abundance covered by shared orders"
  ),
  Value = c(
    length(asv_detected_orders),
    length(mag_detected_orders),
    length(shared_orders),
    length(setdiff(asv_detected_orders, mag_detected_orders)),
    length(setdiff(mag_detected_orders, asv_detected_orders)),
    sum(asv_mean_order[shared_orders]) / sum(asv_mean_order[asv_detected_orders]),
    sum(mag_mean_order[shared_orders]) / sum(mag_mean_order[mag_detected_orders])
    )
  )
write_csv(
  order_coverage,
  file.path(out_dir, "Order_ASV_MAG_coverage_summary.csv")
)

order_correlations <- tibble(Order = shared_orders) %>%
  mutate(
    ASV_mean_abundance = asv_mean_order[Order],
    MAG_mean_abundance = mag_mean_order[Order],
    Spearman_rho = map_dbl(
      Order,
      ~ cor(asv_order_mat[[.x]], mag_order_mat[[.x]], method = "spearman")
    ),
    P_value = map_dbl(
      Order,
      ~ cor.test(
        asv_order_mat[[.x]],
        mag_order_mat[[.x]],
        method = "spearman",
        exact = FALSE
      )$p.value
    ),
    # Apply BH correction across shared orders.
    FDR = p.adjust(P_value, method = "BH"),
    Mean_abundance = (ASV_mean_abundance + MAG_mean_abundance) / 2
  ) %>%
  arrange(desc(Mean_abundance))
write_csv(
  order_correlations,
  file.path(out_dir, "Order_ASV_MAG_spearman_correlations.csv")
)

mantel_seed <- 666L
set.seed(mantel_seed)
mantel_result <- mantel(
  vegdist(asv_order_mat, method = "bray"),
  vegdist(mag_order_mat, method = "bray"),
  method = "spearman",
  permutations = 9999
)
capture.output(
  cat("Permutation seed:", mantel_seed, "\n"),
  print(mantel_result),
  file = file.path(out_dir, "Order_ASV_MAG_mantel_results.txt")
)

# --- Order-level concordance: visualization ---

order_correlation_plot <- order_correlations %>%
  mutate(Order = fct_reorder(Order, Spearman_rho)) %>%
  ggplot(aes(x = Spearman_rho, y = Order)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(aes(size = Mean_abundance, color = FDR < 0.05), alpha = 0.9) +
  scale_x_continuous(limits = c(-1, 1)) +
  scale_size_continuous(labels = scales::percent) +
  scale_color_manual(values = c("TRUE" = "#4F8FC0", "FALSE" = "grey60")) +
  labs(
    x = "Spearman rho",
    y = NULL,
    size = "Mean abundance",
    color = "FDR < 0.05"
  ) +
  theme_bw() +
  theme(panel.grid.minor = element_blank())
ggsave(
  file.path(out_dir, "Order_ASV_MAG_spearman_dotplot.pdf"),
  order_correlation_plot,
  width = 7,
  height = 8
)

# --- Procrustes validation: functions ---

make_sample_match_id <- function(meta_df) {
  # Build a stable cross-dataset key because 16S and MAG sample IDs differ.
  paste(
    meta_df$Source,
    meta_df$Region,
    meta_df$gl_name,
    round(meta_df$distance_to_glacier_snout, 3),
    round(meta_df$GI, 6),
    sep = "|"
  )
}

prepare_16s_procrustes_phyloseq <- function(rds_path) {
  ps_16s <- readRDS(rds_path)

  # Match the 16S NMDS script: collapse technical replicates by base sample ID.
  sample_data(ps_16s)$Sample_Base <- gsub("_.*", "", sample_names(ps_16s))
  ps_16s_merged <- merge_samples(ps_16s, "Sample_Base")

  original_meta <- as.data.frame(sample_data(ps_16s))
  restored_meta <- original_meta[!duplicated(original_meta$Sample_Base), ]
  rownames(restored_meta) <- restored_meta$Sample_Base
  restored_meta$Sample_ID <- restored_meta$Sample_Base
  restored_meta$Sample_Base <- NULL

  sample_data(ps_16s_merged) <- sample_data(restored_meta)

  if (!taxa_are_rows(ps_16s_merged) && taxa_are_rows(ps_16s)) {
    otu_table(ps_16s_merged) <- otu_table(
      t(otu_table(ps_16s_merged)),
      taxa_are_rows = TRUE
    )
  }

  ps_16s_merged %>%
    subset_samples(Source != "water") %>%
    prune_taxa(taxa_sums(.) > 0, .) %>%
    prune_samples(sample_sums(.) > 0, .) %>%
    transform_sample_counts(function(x) x / sum(x))
}

rename_samples_by_match_id <- function(ps_obj, sample_map, sample_col) {
  # Replace dataset-specific sample IDs with the shared matching key.
  new_names <- sample_map$Match_ID[
    match(sample_names(ps_obj), sample_map[[sample_col]])
  ]

  sample_names(ps_obj) <- new_names
  ps_obj
}

order_dist_by_samples <- function(dist_obj, sample_order) {
  # Procrustes compares rows directly, so both distance matrices must use
  # exactly the same sample order after matching.
  dist_mat <- as.matrix(dist_obj)

  as.dist(dist_mat[sample_order, sample_order])
}

# --- Procrustes validation: inputs and sample matching ---

ps_16s_procrustes <- prepare_16s_procrustes_phyloseq(asv_object_path)

ps_mags_procrustes <- PAMIR_MAGs_rela %>%
  subset_samples(!(Source %in% c("water", "cryoconite"))) %>%
  prune_taxa(taxa_sums(.) > 0, .) %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  transform_sample_counts(function(x) x / sum(x))

meta_16s_procrustes <- data.frame(sample_data(ps_16s_procrustes)) %>%
  rownames_to_column("Sample_16S") %>%
  mutate(Match_ID = make_sample_match_id(.))

meta_mags_procrustes <- data.frame(sample_data(ps_mags_procrustes)) %>%
  rownames_to_column("Sample_MAG") %>%
  mutate(Match_ID = make_sample_match_id(.))

# Stop early if the metadata key does not identify samples one-to-one.
duplicated_match_ids <- union(
  meta_16s_procrustes$Match_ID[duplicated(meta_16s_procrustes$Match_ID)],
  meta_mags_procrustes$Match_ID[duplicated(meta_mags_procrustes$Match_ID)]
)

if (length(duplicated_match_ids) > 0) {
  stop(
    "Duplicated sample matching keys detected. Check metadata for: ",
    paste(duplicated_match_ids, collapse = "; ")
  )
}

# Keep only samples present in both datasets and preserve the intended order.
procrustes_sample_map <- inner_join(
  meta_16s_procrustes %>% select(Sample_16S, Match_ID),
  meta_mags_procrustes %>%
    select(Sample_MAG, Match_ID, Source, Region, gl_name, GI, Group_GI),
  by = "Match_ID"
) %>%
  arrange(Region, gl_name, Source, desc(GI))

write_csv(
  procrustes_sample_map,
  file.path(out_dir, "Procrustes_16S_MAG_sample_map.csv")
)

ps_16s_procrustes <- prune_samples(
  procrustes_sample_map$Sample_16S,
  ps_16s_procrustes
)
ps_mags_procrustes <- prune_samples(
  procrustes_sample_map$Sample_MAG,
  ps_mags_procrustes
)

ps_16s_procrustes <- rename_samples_by_match_id(
  ps_16s_procrustes,
  procrustes_sample_map,
  "Sample_16S"
)

ps_mags_procrustes <- rename_samples_by_match_id(
  ps_mags_procrustes,
  procrustes_sample_map,
  "Sample_MAG"
)

shared_match_ids <- procrustes_sample_map$Match_ID
ps_16s_procrustes <- prune_samples(shared_match_ids, ps_16s_procrustes)
ps_mags_procrustes <- prune_samples(shared_match_ids, ps_mags_procrustes)

# Remove taxa absent from the matched subset before calculating Bray-Curtis.
ps_16s_procrustes <- prune_taxa(
  taxa_sums(ps_16s_procrustes) > 0,
  ps_16s_procrustes
)
ps_mags_procrustes <- prune_taxa(
  taxa_sums(ps_mags_procrustes) > 0,
  ps_mags_procrustes
)

# --- Procrustes validation: NMDS and permutation test ---

nmds_protest_seed <- 666L
set.seed(nmds_protest_seed)

# Compute Bray-Curtis distances and explicitly reorder labels to prevent
# accidental sample misalignment between 16S and MAG ordinations.
bray_16s_procrustes <- distance(ps_16s_procrustes, method = "bray") %>%
  order_dist_by_samples(shared_match_ids)
bray_mags_procrustes <- distance(ps_mags_procrustes, method = "bray") %>%
  order_dist_by_samples(shared_match_ids)

nmds_16s_procrustes <- metaMDS(bray_16s_procrustes, k = 2, trymax = 100)
nmds_mags_procrustes <- metaMDS(bray_mags_procrustes, k = 2, trymax = 100)

# Final guard: procrustes() assumes corresponding rows are paired samples.
if (!identical(
  rownames(nmds_16s_procrustes$points),
  rownames(nmds_mags_procrustes$points)
)) {
  stop("16S and MAG NMDS point matrices are not aligned by sample order.")
}

# Symmetric Procrustes allows rotation, translation, reflection, and scaling
# when comparing the two NMDS configurations.
procrustes_16s_mags <- procrustes(
  nmds_16s_procrustes,
  nmds_mags_procrustes,
  symmetric = TRUE
)

protest_16s_mags <- protest(
  nmds_16s_procrustes,
  nmds_mags_procrustes,
  permutations = 9999
)

capture.output(
  {
    cat("Procrustes analysis: 16S vs MAG abundance community structure\n")
    cat("Matched samples:", nrow(procrustes_sample_map), "\n\n")
    cat("16S NMDS stress:", nmds_16s_procrustes$stress, "\n")
    cat("MAG NMDS stress:", nmds_mags_procrustes$stress, "\n\n")
    print(procrustes_16s_mags)
    cat("\n")
    print(protest_16s_mags)
  },
  file = file.path(out_dir, "Procrustes_16S_MAG_results.txt")
)

procrustes_residuals <- procrustes_sample_map %>%
  mutate(Procrustes_residual = residuals(procrustes_16s_mags)[Match_ID])

write_csv(
  procrustes_residuals,
  file.path(out_dir, "Procrustes_16S_MAG_residuals.csv")
)

# --- Procrustes validation: visualization ---

# Supplementary Procrustes figure: paired configurations and sample residuals.
# The GI palette and region shapes match the main MAG NMDS figure. 16S and MAG
# configurations are distinguished by open and filled symbols, respectively.
procrustes_display_levels <- c(
  "ice",
  "High_Glacial_sed",
  "Mid_Glacial_sed",
  "Low_Glacial_sed"
)
procrustes_display_labels <- c(
  "Ice",
  "High GI",
  "Mid GI",
  "Low GI"
)

procrustes_coordinates <- bind_rows(
  as.data.frame(procrustes_16s_mags$X) %>%
    setNames(c("Procrustes1", "Procrustes2")) %>%
    rownames_to_column("Match_ID") %>%
    mutate(Method = "16S ASVs"),
  as.data.frame(procrustes_16s_mags$Yrot) %>%
    setNames(c("Procrustes1", "Procrustes2")) %>%
    rownames_to_column("Match_ID") %>%
    mutate(Method = "MAGs")
) %>%
  left_join(
    procrustes_sample_map %>%
      select(Match_ID, Source, Region, Group_GI),
    by = "Match_ID"
  ) %>%
  mutate(
    Group_GI = factor(
      Group_GI,
      levels = procrustes_display_levels,
      labels = procrustes_display_labels
    ),
    Method = factor(Method, levels = c("16S ASVs", "MAGs"))
  )

procrustes_segments <- procrustes_coordinates %>%
  mutate(Method = recode(Method, "16S ASVs" = "ASV", "MAGs" = "MAG")) %>%
  select(Match_ID, Region, Group_GI, Method, Procrustes1, Procrustes2) %>%
  pivot_wider(
    names_from = Method,
    values_from = c(Procrustes1, Procrustes2),
    names_sep = "_"
  )

procrustes_overlay_plot <- ggplot() +
  geom_segment(
    data = procrustes_segments,
    aes(
      x = Procrustes1_ASV,
      y = Procrustes2_ASV,
      xend = Procrustes1_MAG,
      yend = Procrustes2_MAG
    ),
    color = "grey70",
    linewidth = 0.35,
    alpha = 0.7
  ) +
  geom_point(
    data = filter(procrustes_coordinates, Method == "16S ASVs"),
    aes(
      x = Procrustes1,
      y = Procrustes2,
      color = Group_GI,
      shape = Region
    ),
    fill = "white",
    size = 2.8,
    stroke = 0.8
  ) +
  geom_point(
    data = filter(procrustes_coordinates, Method == "MAGs"),
    aes(
      x = Procrustes1,
      y = Procrustes2,
      color = Group_GI,
      fill = Group_GI,
      shape = Region
    ),
    size = 2.8,
    stroke = 0.8
  ) +
  scale_shape_manual(values = c("Alps" = 21, "Kyrgyzstan" = 24)) +
  scale_color_viridis_d(option = "D", begin = 0.15, end = 0.85) +
  scale_fill_viridis_d(option = "D", begin = 0.15, end = 0.85) +
  coord_equal() +
  labs(
    title = "ASV-MAG Procrustes comparison",
    subtitle = paste0(
      "n = ", nrow(procrustes_sample_map),
      "; PROTEST r = ", sprintf("%.3f", protest_16s_mags$t0),
      "; P = ", format(protest_16s_mags$signif, scientific = FALSE),
      "; open = 16S ASVs, filled = MAGs"
    ),
    x = "Procrustes axis 1",
    y = "Procrustes axis 2",
    color = "Glacial Index group",
    shape = "Region"
  ) +
  guides(fill = "none") +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

ggsave(
  file.path(out_dir, "Procrustes_16S_MAG_NMDS.pdf"),
  procrustes_overlay_plot,
  width = 10,
  height = 11,
  dpi = 300
)
