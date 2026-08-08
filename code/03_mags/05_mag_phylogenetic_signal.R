# Phylogenetic signal in MAG-level GI-association effect sizes
# Manuscript analysis: primary uncertainty-adjusted phylogenetic signal and
# unadjusted validation reported with Figure 3.

rm(list = ls())
gc()
graphics.off()

out_dir <- file.path("results", "03_mags", "phylogenetic_signal")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

library(tidyverse)
library(phyloseq)
library(ape)
library(phytools)

# --- Load data ---

PAMIR_MAGs_rela <- readRDS(
  file.path("results", "03_mags", "intermediate", "PAMIR_MAGs_rela.rds")
)
MAGs_info <- read.csv(
  file.path("results", "03_mags", "intermediate", "MAGs_info.csv"),
  stringsAsFactors = FALSE
)
phylo_tree <- read.tree(
  file.path("data", "processed", "metagenomics", "iqtree_tree.treefile")
)

tax_df <- as.data.frame(as(tax_table(PAMIR_MAGs_rela), "matrix"))
tax_df$OTUID <- rownames(tax_df)

# --- Root the tree with Patescibacteriota as the outgroup ---

patescibacteria_tips <- tax_df %>%
  filter(
    !is.na(Phylum),
    str_detect(Phylum, fixed("Patescibacteriota"))
  ) %>%
  pull(OTUID) %>%
  intersect(phylo_tree$tip.label) %>%
  unique()

rooted_tree <- root(phylo_tree, outgroup = patescibacteria_tips, resolve.root = TRUE)

# --- Align tree tips with MAG-level environmental preference ---

trait_df <- MAGs_info %>%
  select(MAGs, lfc_GI, se_GI, Phylum, Class, Order, Family, Genus) %>%
  filter(
    !is.na(MAGs),
    nzchar(MAGs),
    is.finite(lfc_GI),
    is.finite(se_GI),
    se_GI >= 0
  )

common_mags <- intersect(rooted_tree$tip.label, trait_df$MAGs)

aligned_tree <- keep.tip(rooted_tree, common_mags)

trait_named <- trait_df %>%
  filter(MAGs %in% common_mags) %>%
  { setNames(.$lfc_GI, .$MAGs) } %>%
  .[aligned_tree$tip.label]

trait_se_named <- trait_df %>%
  filter(MAGs %in% common_mags) %>%
  { setNames(.$se_GI, .$MAGs) } %>%
  .[aligned_tree$tip.label]

# --- Phylogenetic-signal analysis ---

# Primary analysis: use the ANCOM-BC2 standard error associated with each
# MAG-specific LFC-GI estimate. phylosig() internally constructs
# the diagonal error-variance matrix as diag(se_GI^2).
set.seed(1)
lambda_fit_se <- phytools::phylosig(
  aligned_tree,
  trait_named,
  method = "lambda",
  test = TRUE,
  se = trait_se_named
)

if (!is.null(lambda_fit_se$convergence) && lambda_fit_se$convergence != 0) {
  stop(
    "Pagel lambda optimization with se_GI did not converge: ",
    lambda_fit_se$message
  )
}

set.seed(1)
k_fit_se <- phytools::phylosig(
  aligned_tree,
  trait_named,
  method = "K",
  test = TRUE,
  nsim = 999,
  se = trait_se_named
)

# Validation analysis: treat the estimated LFC-GI values as error-free.
set.seed(42)
lambda_fit_unadjusted <- phytools::phylosig(
  aligned_tree,
  trait_named,
  method = "lambda",
  test = TRUE
)

set.seed(42)
k_fit_unadjusted <- phytools::phylosig(
  aligned_tree,
  trait_named,
  method = "K",
  test = TRUE,
  nsim = 999
)

signal_results <- bind_rows(
  tibble(
    analysis = "se_GI_adjusted",
    analysis_role = "primary",
    metric = c("Pagel_lambda", "Blomberg_K"),
    estimate = c(unname(lambda_fit_se$lambda), unname(k_fit_se$K)),
    p_value = c(lambda_fit_se$P, k_fit_se$P),
    test = c(
      "likelihood_ratio_vs_lambda_0",
      "permutation_vs_no_signal"
    ),
    n_tips = length(trait_named),
    uses_se_GI = TRUE,
    median_se_GI = median(trait_se_named),
    estimated_sigma2 = c(lambda_fit_se$sig2, k_fit_se$sig2),
    optimization_convergence = c(
      as.integer(lambda_fit_se$convergence),
      NA_integer_
    )
  ),
  tibble(
    analysis = "unadjusted",
    analysis_role = "validation",
    metric = c("Pagel_lambda", "Blomberg_K"),
    estimate = c(
      unname(lambda_fit_unadjusted$lambda),
      unname(k_fit_unadjusted$K)
    ),
    p_value = c(lambda_fit_unadjusted$P, k_fit_unadjusted$P),
    test = c(
      "likelihood_ratio_vs_lambda_0",
      "permutation_vs_no_signal"
    ),
    n_tips = length(trait_named),
    uses_se_GI = FALSE,
    median_se_GI = NA_real_,
    estimated_sigma2 = NA_real_,
    optimization_convergence = NA_integer_
  )
)

write.csv(
  signal_results,
  file.path(out_dir, "phylogenetic_signal_summary.csv"),
  row.names = FALSE
)

write.csv(
  filter(signal_results, analysis == "unadjusted"),
  file.path(out_dir, "phylogenetic_signal_unadjusted.csv"),
  row.names = FALSE
)

write.csv(
  filter(signal_results, analysis == "se_GI_adjusted"),
  file.path(out_dir, "phylogenetic_signal_se_GI_adjusted.csv"),
  row.names = FALSE
)

matching_summary <- tibble(
  category = c(
    "tree_tips", "valid_trait_MAGs", "matched_MAGs",
    "tree_only", "trait_only"
  ),
  n = c(
    length(rooted_tree$tip.label),
    nrow(trait_df),
    length(common_mags),
    length(setdiff(rooted_tree$tip.label, trait_df$MAGs)),
    length(setdiff(trait_df$MAGs, rooted_tree$tip.label))
  )
)

write.csv(
  matching_summary,
  file.path(out_dir, "tree_trait_matching_summary.csv"),
  row.names = FALSE
)

analysis_data <- tibble(
  MAGs = aligned_tree$tip.label,
  lfc_GI = as.numeric(trait_named),
  se_GI = as.numeric(trait_se_named)
)

write.csv(
  analysis_data,
  file.path(out_dir, "phylogenetic_signal_input_data.csv"),
  row.names = FALSE
)

write.tree(
  aligned_tree,
  file = file.path(out_dir, "aligned_rooted_tree.treefile")
)

# --- Plain-text summary ---

summary_text <- c(
  paste0("Matched MAGs: ", length(trait_named)),
  paste0("Median se_GI: ", signif(median(trait_se_named), 4)),
  "",
  "Primary analysis (se_GI-adjusted):",
  paste0(
    "Pagel lambda: ", round(lambda_fit_se$lambda, 3),
    " (LRT p = ", signif(lambda_fit_se$P, 3),
    ", convergence = ", lambda_fit_se$convergence, ")"
  ),
  paste0(
    "Blomberg K: ", round(k_fit_se$K, 3),
    " (permutation p = ", signif(k_fit_se$P, 3),
    ", nsim = 999)"
  ),
  "",
  "Validation analysis (unadjusted):",
  paste0(
    "Pagel lambda: ", round(lambda_fit_unadjusted$lambda, 3),
    " (LRT p = ", signif(lambda_fit_unadjusted$P, 3), ")"
  ),
  paste0(
    "Blomberg K: ", round(k_fit_unadjusted$K, 3),
    " (permutation p = ", signif(k_fit_unadjusted$P, 3),
    ", nsim = 999)"
  )
)

writeLines(
  summary_text,
  file.path(out_dir, "phylogenetic_conservatism_summary.txt")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(out_dir, "sessionInfo.txt")
)
