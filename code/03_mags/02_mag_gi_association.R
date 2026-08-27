# MAG GI-association analysis
# Manuscript output: Figure 3b.
# The MAG-level classifications are also used by Figures 3a, 3c, 4, and 5.

rm(list = ls())
gc()
graphics.off()

out_dir <- file.path("results", "03_mags", "gi_association")
intermediate_dir <- file.path("results", "03_mags", "intermediate")
generated_dir <- file.path("results", "generated", "03_mags", "ancombc2")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(intermediate_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(generated_dir, recursive = TRUE, showWarnings = FALSE)

# --- Packages ---
library(tidyverse)
library(phyloseq)
library(ggplot2)

# --- 1. Read and clean the CoverM count table ---

mags_count <- read.table(
  file.path("data", "processed", "metagenomics", "coverm_mags_count.tsv"),
  header = TRUE,
  sep = "\t",
  quote = "",
  fill = TRUE,
  stringsAsFactors = FALSE
)

rownames(mags_count) <- mags_count$Genome
mags_count$Genome <- NULL

# Clean sample names.
colnames(mags_count) <- gsub(
  "\\_1.fq\\.gz\\.Read\\.Count.*$", "",
  gsub("^cat_mags\\.fa\\.", "", colnames(mags_count))
)

mags_count <- as.data.frame(mags_count)

# --- 2. Read and clean metadata ---

metadata <- read_csv(
  file.path("data", "processed", "metagenomics", "metadata_metagenomics.csv"),
  show_col_types = FALSE
) %>%
  column_to_rownames(var = colnames(.)[1]) %>%
  filter(!is.na(Source)) %>% # Remove samples without a habitat assignment.
  mutate(Source = factor(Source)) %>% # Store Source as a factor.
  mutate(gl_name = factor(gl_name)) %>% # Store glacier identity as a factor.
  select(where(~ !all(is.na(.)))) %>%
  mutate(Sample_ID = rownames(.))

# --- 3. Read and clean GTDB-Tk taxonomy ---

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
  taxonomy_MAGs_raw, 2,
  function(x) {
    x <- gsub("^[a-z]__", "", x)   # Remove GTDB prefixes such as d__, p__, and c__.
    x[x == "" | x == ""] <- "unknown"
    x[x == "s__"] <- "unknown"
    x
  }
)

taxonomy_MAGs <- as.matrix(taxonomy_MAGs)

# --- 4. Build the phyloseq object ---

MAGs_table <- otu_table(mags_count, taxa_are_rows = TRUE)
sample_data_MAGs <- sample_data(metadata)
tax_MAGs <- tax_table(taxonomy_MAGs)

PAMIR_MAGs_count <- merge_phyloseq(
  MAGs_table,
  sample_data_MAGs,
  tax_MAGs
)

print(PAMIR_MAGs_count)

# --- ANCOM-BC2 analysis of MAG responses to GI ---

library(ANCOMBC)

# Retain sediment samples only.
PAMIR_MAGs_count_sed <- subset_samples(PAMIR_MAGs_count, Source == "sed")

# Remove MAGs with zero counts.
PAMIR_MAGs_count_sed <- prune_taxa(
  taxa_sums(PAMIR_MAGs_count_sed) > 0,
  PAMIR_MAGs_count_sed
)

PAMIR_MAGs_count_sed

ancombc_out <- ancombc2(
  data = PAMIR_MAGs_count_sed,
  # Core model specification.
  # Region adjusts intercepts; GI estimates the shared slope.
  fix_formula = "GI + Region",
  rand_formula = "(1 | gl_name)", # Account for non-independence within glaciers.

  # Comparison settings.
  group = NULL,          # Do not define a group factor or pairwise comparison.
  struc_zero = FALSE,    # Structural-zero testing is not required.
  # A continuous predictor does not require a negative lower bound.
  neg_lb = FALSE,
  global = FALSE,        # Use coefficients directly without a global multi-group test.
  pairwise = FALSE,      # Do not test pairwise contrasts for a continuous predictor.
  trend = FALSE,         # Ordered-factor trend testing is not applicable.

  # Statistical and filtering settings.
  # Use BH-FDR correction, which is less conservative than Holm.
  p_adj_method = "fdr",
  prv_cut = 0.10,        # Remove MAGs present in fewer than 10% of samples.
  lib_cut = 1000,        # Remove samples with library sizes below 1,000.
  alpha = 0.05,          # Use a significance threshold of 0.05.

  # Sensitivity analysis and computation.
  pseudo_sens = TRUE,    # Enable the recommended pseudo-count sensitivity analysis.
  n_cl = 8,              # Use eight parallel workers.
  verbose = TRUE         # Show progress.
)

saveRDS(ancombc_out, file = file.path(generated_dir, "ancom_res_sed.rds"))

res_df <- ancombc_out$res
log_table <- ancombc_out$bias_correct_log_table

# --- Select robust GI-associated MAGs ---
# Select robust MAGs that are significant after correction (q < 0.05)
# and pass the sensitivity test.

final_mags_trend <- res_df %>%
  rename(MAGs = taxon) %>%
  # Define response categories.
  mutate(Trend = case_when(
    diff_robust_GI != TRUE ~ "Non-sig",  # Explicit FALSE values are non-significant.
    lfc_GI > 0 ~ "Glaciophiles",         # Significant increase with GI.
    lfc_GI < 0 ~ "Glaciophobes"          # Significant decrease with GI.
  )) %>%
  arrange(Trend == "Non-sig", desc(lfc_GI)) # Sort significant MAGs first.

write_csv(
  final_mags_trend,
  file.path(out_dir, "MAG_GI_association_results.csv")
)

# Merge taxonomy and ANCOM-BC2 results.
taxonomy_MAGs <- data.frame(taxonomy_MAGs) %>% rownames_to_column("MAGs")

MAGs_info <- taxonomy_MAGs %>%
  left_join(
    final_mags_trend %>%
      select(MAGs, lfc_GI, se_GI, Trend),
    by = "MAGs"
  )

MAGs_info$Trend <- factor(
  MAGs_info$Trend,
  levels = c("Glaciophiles", "Non-sig", "Glaciophobes"),
  ordered = TRUE
)

count_num <- table(MAGs_info$Trend)
print(count_num)
summary(MAGs_info$lfc_GI)

# For Pseudomonadota, relabel Alpha/Gamma classes at phylum level.
MAGs_info <- MAGs_info %>%
  mutate(Phylum = case_when(
    Phylum == "Pseudomonadota" & Class %in% c("Alphaproteobacteria", "Gammaproteobacteria") ~ Class,
    TRUE ~ Phylum
  ))

table(MAGs_info$Phylum)

# Add MAG quality metrics.
MAGs_quality_report <- read.delim(
  file.path("data", "processed", "metagenomics", "quality_report.tsv"),
  header = TRUE,
  quote = ""
)

MAGs_quality_report <- dplyr::rename(MAGs_quality_report, MAGs = Name)

HQ_MAGs_count <- sum(MAGs_quality_report$Completeness > 90 & MAGs_quality_report$Contamination < 5)
print(HQ_MAGs_count)

MAGs_info <- merge(MAGs_info, MAGs_quality_report, by = "MAGs", all.x = TRUE)

# Save the MAG-level results.
write.csv(
  MAGs_info,
  file.path(intermediate_dir, "MAGs_info.csv"),
  row.names = FALSE
)

# --- Export Polaromonas inputs ---

# Export Polaromonas MAG phenotypes for downstream pangenome analysis.
MAGs_phenotype_polar <- MAGs_info %>%
  select(MAGs, lfc_GI, Genus) %>%
  filter(Genus == "Polaromonas") %>%
  select(MAGs, lfc_GI)

write.table(
  MAGs_phenotype_polar,
  file = file.path(intermediate_dir, "MAGs_phenotype_polar.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# Export Polaromonas MAG completeness values for downstream pangenome analysis.
MAGs_completeness_polar <- MAGs_info %>%
  select(MAGs, Completeness, Genus) %>%
  filter(Genus == "Polaromonas") %>%
  select(MAGs, Completeness)

write.table(
  MAGs_completeness_polar,
  file = file.path(intermediate_dir, "MAGs_completeness_polar.tsv"),
  sep = "\t",
  row.names = FALSE, # Do not write row names.
  quote = FALSE
)

# --- LFC-GI histogram ---

# Color each histogram bin by its LFC-GI midpoint. This continuous diverging
# scale matches the LFC-GI ring in Final_MAGs_Phylo_iqtree (Fig. 3a).
# Use a symmetric range derived from all MAGs so identical LFC-GI values have
# identical colors in Fig. 3a and Fig. 3b.
lfc_GI_color_limit <- max(abs(MAGs_info$lfc_GI), na.rm = TRUE)
lfc_GI_color_limits <- c(-lfc_GI_color_limit, lfc_GI_color_limit)

hist_plot <- ggplot(MAGs_info, aes(x = lfc_GI)) +
  geom_histogram(
    aes(fill = after_stat(x)),
    bins = 100,
    color = "grey35",
    linewidth = 0.15
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey25",
    linewidth = 0.4
  ) +
  scale_fill_gradient2(
    low = "#3288BD",   # Negative LFC-GI values.
    mid = "white",     # Zero.
    high = "#D53E4F",  # Positive LFC-GI values.
    midpoint = 0,
    limits = lfc_GI_color_limits,
    oob = scales::squish,
    name = "LFC GI",
    guide = guide_colorbar(
      order = 4,
      barwidth = 0.5,
      barheight = 5,
      title.position = "top",
      frame.colour = "black",
      frame.linewidth = 0.5
    )
  ) +
  theme_bw() +
  labs(
    title = "Distribution of MAG Responses to GI",
    x = "log fold change (GI)",
    y = "Count"
  )

# Save the LFC-GI distribution as a vector PDF.
ggsave(
  file.path(out_dir, "MAG_LFC_GI_distribution.pdf"),
  plot = hist_plot,
  width = 6,
  height = 4.5,
  device = grDevices::cairo_pdf,
  bg = "white"
)
