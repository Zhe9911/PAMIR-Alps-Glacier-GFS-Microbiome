# Polaromonas clade LFC-GI contrast, permutation test, and ANI assessment
# Manuscript output: Figure 5b with quality, uncertainty, and ANI validation.

rm(list = ls())
gc()
graphics.off()

library(tidyverse)

# --- 1. Paths and analysis parameters ---

genus <- "Polaromonas"
input_dir <- file.path("data", "processed", "polaromonas")
out_dir <- file.path("results", "05_polaromonas", "clade_contrast")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

clade_file <- file.path(input_dir, "MAGs_Group_clade.tsv")
mags_info_file <- file.path(
  "results",
  "03_mags",
  "intermediate",
  "MAGs_info.csv"
)
ani_file <- file.path(input_dir, "ANIb_percentage_identity.tsv")

n_background_perm <- 9999
n_error_draws <- 1000
background_permutation_seed <- 123
measurement_error_seed <- 124
metadata_jitter_seed <- 125
lfc_jitter_seed <- 126
max_exact_partitions <- 1e6
ani_species_threshold <- 95

# --- 2. Read clade assignments and Polaromonas LFC-GI estimates ---

# Labels clade1 and clade2 correspond to Clade-D and Clade-G,
# respectively.

clade_df <- read_tsv(clade_file, show_col_types = FALSE) %>%
  transmute(
    MAGs = as.character(MAGs),
    Group_clade = factor(Group_clade, levels = c("clade1", "clade2"))
  )

MAGs_info <- read.csv(mags_info_file, stringsAsFactors = FALSE)

genus_lfc_df <- MAGs_info %>%
  transmute(
    MAGs = as.character(MAGs),
    Genus = as.character(Genus),
    lfc_GI = as.numeric(lfc_GI),
    se_GI = as.numeric(se_GI)
  ) %>%
  filter(Genus == genus, !is.na(lfc_GI))

clade_lfc_df <- genus_lfc_df %>%
  inner_join(clade_df, by = "MAGs") %>%
  filter(!is.na(Group_clade))

other_lfc_df <- genus_lfc_df %>%
  anti_join(clade_df, by = "MAGs") %>%
  mutate(Group_clade = "other_Polaromonas")

lfc_plot_df <- bind_rows(clade_lfc_df, other_lfc_df) %>%
  mutate(
    Group_clade = factor(
      as.character(Group_clade),
      levels = c("clade1", "clade2", "other_Polaromonas")
    )
  )

clade_stats <- clade_lfc_df %>%
  group_by(Group_clade) %>%
  summarise(
    n = n(),
    median = median(lfc_GI),
    mean = mean(lfc_GI),
    iqr = IQR(lfc_GI),
    .groups = "drop"
  )

# --- 3. Quantify the Clade-G minus Clade-D LFC-GI contrast ---

clade1_vals <- clade_lfc_df$lfc_GI[clade_lfc_df$Group_clade == "clade1"]
clade2_vals <- clade_lfc_df$lfc_GI[clade_lfc_df$Group_clade == "clade2"]

observed_contrast <- median(clade2_vals) - median(clade1_vals)

# A positive median contrast or Cliff's delta indicates higher LFC-GI in
# Clade-G than in Clade-D.
cliffs_delta <- mean(outer(clade2_vals, clade1_vals, ">")) -
  mean(outer(clade2_vals, clade1_vals, "<"))

# --- 4. Test MAG quality and genome-property differences between clades ---

# These comparisons assess potential technical or genomic confounding of the
# ecological LFC-GI contrast; BH correction is applied across metadata metrics.

metadata <- clade_df %>%
  left_join(
    MAGs_info %>%
      transmute(
        MAGs = as.character(MAGs),
        lfc_GI = as.numeric(lfc_GI),
        se_GI = as.numeric(se_GI),
        Trend = as.character(Trend),
        Completeness = as.numeric(Completeness),
        Contamination = as.numeric(Contamination),
        Coding_Density = as.numeric(Coding_Density),
        Contig_N50 = as.numeric(Contig_N50),
        Average_Gene_Length = as.numeric(Average_Gene_Length),
        Genome_Size = as.numeric(Genome_Size),
        log10_genome_size = log10(as.numeric(Genome_Size)),
        GC_Content = as.numeric(GC_Content),
        Total_Coding_Sequences = as.numeric(Total_Coding_Sequences),
        Total_Contigs = as.numeric(Total_Contigs),
        Max_Contig_Length = as.numeric(Max_Contig_Length)
      ),
    by = "MAGs"
  ) %>%
  filter(
    !is.na(lfc_GI),
    !is.na(Completeness),
    !is.na(log10_genome_size)
  )

metadata_metric_specs <- tibble(
  metric = c(
    "Completeness",
    "Contamination",
    "Genome_Size",
    "GC_Content",
    "Contig_N50",
    "Total_Contigs"
  ),
  label = c(
    "Completeness (%)",
    "Contamination (%)",
    "Genome size (bp)",
    "GC content",
    "Contig N50 (bp)",
    "Total contigs"
  )
)

metadata_plot_df <- metadata %>%
  select(Group_clade, all_of(metadata_metric_specs$metric)) %>%
  pivot_longer(
    cols = all_of(metadata_metric_specs$metric),
    names_to = "metric",
    values_to = "value"
  ) %>%
  filter(!is.na(value))

format_p_value <- function(p) {
  ifelse(is.na(p), "NA", format.pval(p, digits = 2, eps = 1e-4))
}

metadata_metric_tests <- metadata_plot_df %>%
  group_by(metric) %>%
  group_modify(~ {
    if (n_distinct(.x$Group_clade) < 2 || sd(.x$value) == 0) {
      return(tibble(median_diff_clade2_minus_clade1 = NA_real_, p.value = NA_real_))
    }

    tibble(
      median_diff_clade2_minus_clade1 = unname(diff(tapply(.x$value, .x$Group_clade, median))),
      p.value = wilcox.test(value ~ Group_clade, data = .x, exact = FALSE)$p.value
    )
  }) %>%
  ungroup() %>%
  left_join(metadata_metric_specs, by = "metric") %>%
  mutate(
    adjusted_q_value = p.adjust(p.value, method = "BH"),
    facet_label = paste0(label, "\nWilcoxon p = ", format_p_value(p.value))
  )

metadata_plot_df <- metadata_plot_df %>%
  left_join(
    metadata_metric_tests %>% select(metric, facet_label),
    by = "metric"
  ) %>%
  mutate(
    facet_label = factor(
      facet_label,
      levels = metadata_metric_tests$facet_label[
        match(metadata_metric_specs$metric, metadata_metric_tests$metric)
      ]
    )
  )

p_metadata_metrics <- ggplot(
  metadata_plot_df,
  aes(x = Group_clade, y = value, fill = Group_clade)
) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.85, color = "grey25") +
  geom_jitter(
    aes(color = Group_clade),
    position = position_jitter(width = 0.12, height = 0, seed = metadata_jitter_seed),
    size = 1.8,
    alpha = 0.85
  ) +
  facet_wrap(~ facet_label, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c("clade1" = "#3288BD", "clade2" = "#D53E4F")) +
  scale_color_manual(values = c("clade1" = "#2166AC", "clade2" = "#B2182B")) +
  scale_x_discrete(labels = c("clade1" = "Clade-D", "clade2" = "Clade-G")) +
  scale_y_continuous(labels = scales::label_number()) +
  labs(
    title = "Polaromonas clade-level metadata comparisons",
    x = NULL,
    y = NULL
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.background = element_rect(fill = "grey92", color = "grey70"),
    strip.text = element_text(size = 8, face = "bold")
  )

# --- 5. Test the Clade-G minus Clade-D median LFC-GI contrast by permutation ---

n_clade1 <- sum(clade_lfc_df$Group_clade == "clade1")
n_clade2 <- sum(clade_lfc_df$Group_clade == "clade2")

# Primary test: exhaustively reassign the observed LFC-GI estimates while
# preserving the two observed clade sizes. This provides the exact null
# distribution for the topology-defined groups.
target_lfc <- clade_lfc_df$lfc_GI
target_se <- clade_lfc_df$se_GI
target_group_space <- choose(length(target_lfc), n_clade1)

if (!is.finite(target_group_space) || target_group_space > max_exact_partitions) {
  stop(
    "Exact target-clade enumeration requires ", target_group_space,
    " partitions, exceeding max_exact_partitions = ", max_exact_partitions
  )
}

if (any(!is.finite(target_se)) || any(target_se < 0)) {
  stop("All target-clade MAGs require finite, non-negative se_GI values")
}

clade1_partition_indices <- combn(seq_along(target_lfc), n_clade1)
clade2_partition_indices <- apply(
  clade1_partition_indices,
  2,
  function(clade1_idx) setdiff(seq_along(target_lfc), clade1_idx)
)

calculate_partition_contrasts <- function(values) {
  clade1_value_matrix <- matrix(
    values[clade1_partition_indices],
    nrow = n_clade1
  )
  clade2_value_matrix <- matrix(
    values[clade2_partition_indices],
    nrow = n_clade2
  )

  matrixStats::colMedians(clade2_value_matrix) -
    matrixStats::colMedians(clade1_value_matrix)
}

target_contrasts <- calculate_partition_contrasts(target_lfc)
target_n_partitions <- length(target_contrasts)
target_extreme_count <- sum(abs(target_contrasts) >= abs(observed_contrast))
target_p_two_sided_absolute <- target_extreme_count / target_n_partitions

# Sensitivity analysis 1: propagate each MAG-specific LFC-GI standard error and
# repeat the exact permutation test for every simulated dataset.
set.seed(measurement_error_seed)
simulated_lfc_matrix <- matrix(
  rnorm(
    length(target_lfc) * n_error_draws,
    mean = rep(target_lfc, times = n_error_draws),
    sd = rep(target_se, times = n_error_draws)
  ),
  nrow = length(target_lfc),
  ncol = n_error_draws
)

observed_clade1_indices <- which(clade_lfc_df$Group_clade == "clade1")
observed_clade2_indices <- which(clade_lfc_df$Group_clade == "clade2")

measurement_error_draws <- map_dfr(seq_len(n_error_draws), function(draw_id) {
  simulated_lfc <- simulated_lfc_matrix[, draw_id]
  simulated_observed_contrast <-
    median(simulated_lfc[observed_clade2_indices]) -
    median(simulated_lfc[observed_clade1_indices])
  simulated_partition_contrasts <- calculate_partition_contrasts(simulated_lfc)
  simulated_extreme_count <- sum(
    abs(simulated_partition_contrasts) >= abs(simulated_observed_contrast)
  )

  tibble(
    draw_id = draw_id,
    simulated_observed_contrast = simulated_observed_contrast,
    exact_extreme_count = simulated_extreme_count,
    exact_n_partitions = target_n_partitions,
    exact_p_two_sided_absolute = simulated_extreme_count / target_n_partitions
  )
})

measurement_error_summary <- measurement_error_draws %>%
  summarise(
    n_draws = n(),
    seed = measurement_error_seed,
    observed_point_estimate_contrast = observed_contrast,
    simulated_contrast_median = median(simulated_observed_contrast),
    simulated_contrast_q025 = quantile(simulated_observed_contrast, 0.025),
    simulated_contrast_q975 = quantile(simulated_observed_contrast, 0.975),
    probability_contrast_greater_than_zero = mean(simulated_observed_contrast > 0),
    exact_p_median = median(exact_p_two_sided_absolute),
    exact_p_q025 = quantile(exact_p_two_sided_absolute, 0.025),
    exact_p_q975 = quantile(exact_p_two_sided_absolute, 0.975),
    proportion_exact_p_below_0.05 = mean(exact_p_two_sided_absolute < 0.05),
    method_note = paste(
      "Independent normal measurement-error propagation using each MAG's lfc_GI and se_GI;",
      "all exact target-clade partitions were re-evaluated in every draw."
    )
  )

# Sensitivity analysis 2: compare the observed contrast with Monte Carlo groups
# of identical sizes sampled without replacement from all Polaromonas MAGs.
set.seed(background_permutation_seed)
genus_contrasts <- replicate(n_background_perm, {
  sampled_lfc <- sample(genus_lfc_df$lfc_GI, size = n_clade1 + n_clade2, replace = FALSE)
  median(sampled_lfc[(n_clade1 + 1):(n_clade1 + n_clade2)]) -
    median(sampled_lfc[1:n_clade1])
})

genus_extreme_count <- sum(abs(genus_contrasts) >= abs(observed_contrast))
genus_p_two_sided_absolute <- (genus_extreme_count + 1) /
  (n_background_perm + 1)

lfc_test_summary <- tibble(
  test = "exhaustive_target_clade_permutation",
  comparison = "clade2_minus_clade1",
  n_clade1 = n_clade1,
  n_clade2 = n_clade2,
  median_clade1 = clade_stats$median[clade_stats$Group_clade == "clade1"],
  median_clade2 = clade_stats$median[clade_stats$Group_clade == "clade2"],
  mean_clade1 = clade_stats$mean[clade_stats$Group_clade == "clade1"],
  mean_clade2 = clade_stats$mean[clade_stats$Group_clade == "clade2"],
  iqr_clade1 = clade_stats$iqr[clade_stats$Group_clade == "clade1"],
  iqr_clade2 = clade_stats$iqr[clade_stats$Group_clade == "clade2"],
  effect_size_type = "median_lfc_GI_difference",
  effect_size = observed_contrast,
  auxiliary_effect_size_type = "cliffs_delta",
  auxiliary_effect_size = cliffs_delta,
  n_exact_partitions = target_n_partitions,
  exact_extreme_count = target_extreme_count,
  p_value_definition = "Pr(|permuted median contrast| >= |observed median contrast|)",
  p.value = target_p_two_sided_absolute,
  method_note = paste(
    "Primary analysis; exhaustive enumeration of all fixed-size partitions.",
    "Positive effect means clade2 has higher LFC-GI."
  )
)

null_distribution <- bind_rows(
  tibble(
    background = "target_clades_only",
    perm_id = seq_along(target_contrasts),
    permuted_contrast = target_contrasts
  ),
  tibble(
    background = "all_Polaromonas_MAGs_with_lfc_GI",
    perm_id = seq_along(genus_contrasts),
    permuted_contrast = genus_contrasts
  )
)

permutation_summary <- bind_rows(
  tibble(
    test = "target_clade_lfcGI_permutation",
    background = "target_clades_only",
    comparison = "clade2_minus_clade1",
    n_clade1 = n_clade1,
    n_clade2 = n_clade2,
    background_n = nrow(clade_lfc_df),
    permutation_type = "exhaustive_exact",
    n_partitions = target_n_partitions,
    seed = NA_integer_,
    observed_contrast = observed_contrast,
    null_mean = mean(target_contrasts),
    null_sd = sd(target_contrasts),
    null_q025 = quantile(target_contrasts, 0.025),
    null_q500 = quantile(target_contrasts, 0.5),
    null_q975 = quantile(target_contrasts, 0.975),
    extreme_count = target_extreme_count,
    p_two_sided_absolute = target_p_two_sided_absolute,
    method_note = "Primary analysis: exhaustive enumeration of all target-clade partitions with observed clade sizes; statistic is median(clade2) - median(clade1), and the two-sided P value uses the absolute contrast."
  ),
  tibble(
    test = "genus_background_random_partition",
    background = "all_Polaromonas_MAGs_with_lfc_GI",
    comparison = "clade2_minus_clade1",
    n_clade1 = n_clade1,
    n_clade2 = n_clade2,
    background_n = nrow(genus_lfc_df),
    permutation_type = "monte_carlo",
    n_partitions = n_background_perm,
    seed = background_permutation_seed,
    observed_contrast = observed_contrast,
    null_mean = mean(genus_contrasts),
    null_sd = sd(genus_contrasts),
    null_q025 = quantile(genus_contrasts, 0.025),
    null_q500 = quantile(genus_contrasts, 0.5),
    null_q975 = quantile(genus_contrasts, 0.975),
    extreme_count = genus_extreme_count,
    p_two_sided_absolute = genus_p_two_sided_absolute,
    method_note = "Sensitivity analysis: Monte Carlo random partitions sampled without replacement from all Polaromonas MAGs while preserving target clade sizes; statistic is median(random clade2-sized group) - median(random clade1-sized group), and the two-sided P value uses the absolute contrast."
  )
)

# --- 6. Assess within- and between-clade ANI structure ---

ani_df <- read.table(
  ani_file,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

ani_mat <- as.matrix(ani_df[, -1, drop = FALSE])
storage.mode(ani_mat) <- "numeric"
rownames(ani_mat) <- ani_df[[1]]

target_mags_for_ani <- intersect(clade_df$MAGs, intersect(rownames(ani_mat), colnames(ani_mat)))

# Convert ANI proportions to percentages when required, then summarize pairwise
# values relative to the predefined species-level threshold.
ani_scale <- if (max(ani_mat, na.rm = TRUE) <= 1.5) 100 else 1
ani_pairs <- combn(target_mags_for_ani, 2, simplify = FALSE)

ani_pairwise <- map_dfr(ani_pairs, function(pair) {
  mag1 <- pair[1]
  mag2 <- pair[2]
  clade_mag1 <- as.character(clade_df$Group_clade[match(mag1, clade_df$MAGs)])
  clade_mag2 <- as.character(clade_df$Group_clade[match(mag2, clade_df$MAGs)])
  ani_forward <- ani_mat[mag1, mag2] * ani_scale
  ani_reverse <- ani_mat[mag2, mag1] * ani_scale

  tibble(
    MAG1 = mag1,
    MAG2 = mag2,
    clade_MAG1 = clade_mag1,
    clade_MAG2 = clade_mag2,
    pair_type = case_when(
      clade_mag1 == "clade1" & clade_mag2 == "clade1" ~ "within_clade1",
      clade_mag1 == "clade2" & clade_mag2 == "clade2" ~ "within_clade2",
      TRUE ~ "between_clades"
    ),
    ANI_forward_percent = ani_forward,
    ANI_reverse_percent = ani_reverse,
    ANI_mean_percent = mean(c(ani_forward, ani_reverse), na.rm = TRUE)
  )
})

ani_species_summary <- ani_pairwise %>%
  group_by(pair_type) %>%
  summarise(
    n_pairs = n(),
    min_ANI_percent = min(ANI_mean_percent, na.rm = TRUE),
    median_ANI_percent = median(ANI_mean_percent, na.rm = TRUE),
    mean_ANI_percent = mean(ANI_mean_percent, na.rm = TRUE),
    max_ANI_percent = max(ANI_mean_percent, na.rm = TRUE),
    n_pairs_at_or_above_species_threshold = sum(
      ANI_mean_percent >= ani_species_threshold,
      na.rm = TRUE
    ),
    species_threshold_percent = ani_species_threshold,
    interpretation = case_when(
      first(pair_type) == "between_clades" & max_ANI_percent < ani_species_threshold ~
        "Between-clade ANI is below the species boundary; the two clades are genomically separated at species level.",
      first(pair_type) != "between_clades" & max_ANI_percent < ani_species_threshold ~
        "Within-clade ANI is below the species boundary; this clade should not be interpreted as one single species by ANI alone.",
      TRUE ~
        "At least one pair reaches the species boundary."
    ),
    .groups = "drop"
  ) %>%
  arrange(match(pair_type, c("within_clade1", "within_clade2", "between_clades")))

# --- 7. Plot the ecological contrast and robustness analyses ---

plot_group_counts <- lfc_plot_df %>%
  count(Group_clade)

plot_group_labels <- c(
  "clade1" = paste0(
    "Clade-D\n(n=",
    plot_group_counts$n[plot_group_counts$Group_clade == "clade1"],
    ")"
  ),
  "clade2" = paste0(
    "Clade-G\n(n=",
    plot_group_counts$n[plot_group_counts$Group_clade == "clade2"],
    ")"
  ),
  "other_Polaromonas" = paste0(
    "other\nPolaromonas\n(n=",
    plot_group_counts$n[
      plot_group_counts$Group_clade == "other_Polaromonas"
    ],
    ")"
  )
)

p_lfc <- ggplot(
  lfc_plot_df,
  aes(x = Group_clade, y = lfc_GI, fill = Group_clade)
) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.85, color = "grey25") +
  geom_jitter(
    aes(color = Group_clade),
    position = position_jitter(width = 0.12, height = 0, seed = lfc_jitter_seed),
    size = 2.6,
    alpha = 0.9
  ) +
  scale_fill_manual(
    values = c(
      "clade1" = "grey75",
      "clade2" = "grey75",
      "other_Polaromonas" = "grey75"
    )
  ) +
  scale_color_manual(
    values = c(
      "clade1" = "grey25",
      "clade2" = "grey25",
      "other_Polaromonas" = "grey25"
    )
  ) +
  scale_x_discrete(labels = plot_group_labels) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3, color = "grey40") +
  labs(
    title = "Polaromonas clade-level glacier-index association",
    subtitle = paste0(
      "Median contrast = ",
      signif(observed_contrast, 3),
      "; target two-sided permutation p = ",
      signif(target_p_two_sided_absolute, 3)
    ),
    x = NULL,
    y = "Polaromonas lfc_GI"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  )

permutation_plot_df <- null_distribution %>%
  mutate(
    background_label = recode(
      background,
      target_clades_only = "target clades only",
      all_Polaromonas_MAGs_with_lfc_GI = "all Polaromonas MAGs with lfc_GI"
    )
  )

p_perm <- ggplot(
  permutation_plot_df,
  aes(x = permuted_contrast)
) +
  geom_histogram(bins = 40, color = "white", fill = "grey70") +
  geom_vline(
    xintercept = observed_contrast,
    linetype = "dashed",
    linewidth = 0.8,
    color = "#D53E4F"
  ) +
  facet_wrap(~ background_label, ncol = 1) +
  labs(
    title = "Permutation test of Polaromonas sister-clade lfc_GI contrast",
    subtitle = paste0(
      "Observed median contrast = ",
      signif(observed_contrast, 3),
      "\ntarget-only two-sided p = ",
      signif(target_p_two_sided_absolute, 3),
      "; all-Polaromonas two-sided p = ",
      signif(genus_p_two_sided_absolute, 3)
    ),
    x = "Permuted median lfc_GI contrast: Clade-G - Clade-D",
    y = "Permutation count"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 9, lineheight = 0.95)
  )

p_measurement_error <- ggplot(
  measurement_error_draws,
  aes(x = simulated_observed_contrast)
) +
  geom_histogram(bins = 40, color = "white", fill = "grey70") +
  geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.6) +
  geom_vline(
    xintercept = observed_contrast,
    linetype = "dashed",
    linewidth = 0.8,
    color = "#D53E4F"
  ) +
  labs(
    title = "LFC-GI measurement-error sensitivity analysis",
    subtitle = paste0(
      "Pr(contrast > 0) = ",
      signif(measurement_error_summary$probability_contrast_greater_than_zero, 3),
      "; proportion exact P < 0.05 = ",
      signif(measurement_error_summary$proportion_exact_p_below_0.05, 3)
    ),
    x = "Simulated median LFC-GI contrast: Clade-G - Clade-D",
    y = "Simulation count"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  )

# --- 8. Export statistical summaries, pairwise ANI results, and figures ---

write.csv(
  metadata,
  file.path(out_dir, "Polaromonas_clade_metadata.csv"),
  row.names = FALSE
)

write.csv(
  metadata_metric_tests,
  file.path(out_dir, "Polaromonas_clade_metadata_metric_wilcox_summary.csv"),
  row.names = FALSE
)

write.csv(
  lfc_test_summary,
  file.path(out_dir, "Polaromonas_clade_lfcGI_test_summary.csv"),
  row.names = FALSE
)

write.csv(
  permutation_summary,
  file.path(out_dir, "Polaromonas_clade_lfcGI_permutation_summary.csv"),
  row.names = FALSE
)

write.csv(
  null_distribution,
  file.path(out_dir, "Polaromonas_clade_lfcGI_permutation_null_distribution.csv"),
  row.names = FALSE
)

write.csv(
  measurement_error_draws,
  file.path(out_dir, "Polaromonas_clade_lfcGI_measurement_error_draws.csv"),
  row.names = FALSE
)

write.csv(
  measurement_error_summary,
  file.path(out_dir, "Polaromonas_clade_lfcGI_measurement_error_summary.csv"),
  row.names = FALSE
)

write.csv(
  ani_pairwise,
  file.path(out_dir, "Polaromonas_clade_ANI_pairwise.csv"),
  row.names = FALSE
)

write.csv(
  ani_species_summary,
  file.path(out_dir, "Polaromonas_clade_ANI_species_summary.csv"),
  row.names = FALSE
)

ggsave(
  file.path(out_dir, "Polaromonas_clade_metadata_metrics_boxplot.pdf"),
  p_metadata_metrics,
  width = 10,
  height = 6.5
)

ggsave(
  file.path(out_dir, "Polaromonas_clade_lfcGI_boxplot.pdf"),
  p_lfc,
  width = 5,
  height = 5
)

ggsave(
  file.path(out_dir, "Polaromonas_clade_lfcGI_permutation_null_distribution.pdf"),
  p_perm,
  width = 6,
  height = 5
)

ggsave(
  file.path(out_dir, "Polaromonas_clade_lfcGI_measurement_error_sensitivity.pdf"),
  p_measurement_error,
  width = 6,
  height = 4.5
)

print(lfc_test_summary)
print(permutation_summary)
print(measurement_error_summary)
print(ani_species_summary)
