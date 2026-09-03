# iCAMP pairwise process analysis for Figure 2f
# Manuscript outputs: Figure 2f and Extended Data Table 3.

# --- Analysis overview ---

# This script tests associations between sediment site-pair process fractions and
# their mean GI (`GI_mean`) and absolute GI difference (`GI_difference`). Both
# terms enter the same model, so each coefficient is adjusted for the other.
# Replicate-pair iCAMP results are averaged to unordered site pairs; iCAMP is not
# rerun.

# Because site pairs share nodes, inference uses Freedman-Lane MRQAP with
# synchronous node permutations within each glacier. This preserves Region,
# glacier membership and dyadic structure. The model also controls same-glacier
# status, and each site pair receives equal weight.

project_root <- normalizePath(
  Sys.getenv("PAMIR_PAPER_ROOT", unset = "."),
  winslash = "/",
  mustWork = TRUE
)

out_dir <- file.path(project_root, "results", "02_community", "icamp")
ps_path <- file.path(
  project_root, "results", "01_16s", "intermediate", "PAMIR_16S_final.rds"
)
generated_replicate_pair_path <- file.path(
  out_dir, "sediment_pairwise_process_fractions.csv"
)
frozen_replicate_pair_path <- file.path(
  project_root, "data", "derived", "frozen", "16s", "icamp",
  "sediment_pairwise_process_fractions.csv"
)
default_site_pair_path <- file.path(
  out_dir, "sediment_site_pair_process_fractions.csv"
)

# --- Setup ---

library(dplyr)

# --- Analysis parameters ---

params <- list(
  n_permutations = 9999L,
  seed = 20260715L,
  write_outputs = TRUE
)

# --- Load input data ---

if (file.exists(generated_replicate_pair_path)) {
  replicate_pair_path <- generated_replicate_pair_path
} else if (file.exists(frozen_replicate_pair_path)) {
  replicate_pair_path <- frozen_replicate_pair_path
} else {
  stop(
    "Missing iCAMP pairwise process fractions. Expected either ",
    generated_replicate_pair_path, " (full iCAMP result) or ",
    frozen_replicate_pair_path, " (Zenodo frozen result).",
    call. = FALSE
  )
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
message("Using iCAMP pairwise process fractions: ", replicate_pair_path)

# --- Build site-pair data ---

processes <- c("HoS", "HeS", "HD", "DL", "DR")

# Recover site, glacier and GI metadata for each biological replicate.
ps <- readRDS(ps_path)
sample_metadata <- data.frame(
  phyloseq::sample_data(ps),
  check.names = FALSE
)
sample_metadata <- sample_metadata %>%
  dplyr::transmute(
    Sample_ID = as.character(Sample_ID),
    Source = as.character(Source),
    Region = as.character(Region),
    glacier = as.character(gl_name),
    GI = as.numeric(as.character(GI)),
    Site = stringr::str_replace(Sample_ID, "_[A-Za-z]+$", "")
  ) %>%
  dplyr::filter(Source == "sed")

# Load and standardize the completed replicate-pair iCAMP fractions.
replicate_pairs <- readr::read_csv(replicate_pair_path, show_col_types = FALSE)
replicate_pairs <- replicate_pairs %>%
  dplyr::mutate(
    Region = as.character(Region),
    Sample1 = as.character(Sample1),
    Sample2 = as.character(Sample2),
    dplyr::across(dplyr::all_of(processes), as.numeric)
  )

# Attach site, glacier and GI metadata to both endpoints of each pair.
endpoint1 <- sample_metadata %>%
  dplyr::transmute(
    Sample1 = Sample_ID,
    Site1 = Site,
    Region1 = Region,
    GI1 = GI
  )
endpoint2 <- sample_metadata %>%
  dplyr::transmute(
    Sample2 = Sample_ID,
    Site2 = Site,
    Region2 = Region,
    GI2 = GI
  )
replicate_pairs <- replicate_pairs %>%
  dplyr::left_join(endpoint1, by = "Sample1") %>%
  dplyr::left_join(endpoint2, by = "Sample2")
if (anyNA(replicate_pairs[c("Site1", "Site2", "GI1", "GI2")]) ||
    any(replicate_pairs$Region != replicate_pairs$Region1) ||
    any(replicate_pairs$Region != replicate_pairs$Region2)) {
  stop("Replicate-pair endpoints do not match the current sediment metadata.", call. = FALSE)
}

# Average replicate comparisons within each unordered pair of distinct sites.
site_pairs <- replicate_pairs %>%
  dplyr::filter(Site1 != Site2) %>%
  dplyr::mutate(
    SiteA = pmin(Site1, Site2),
    SiteB = pmax(Site1, Site2)
  ) %>%
  dplyr::group_by(Region, SiteA, SiteB) %>%
  dplyr::summarise(
    n_replicate_pairs = dplyr::n(),
    dplyr::across(dplyr::all_of(processes), mean),
    .groups = "drop"
  )

# Reduce replicate metadata to one record per sediment site before rejoining it.
site_metadata <- sample_metadata %>%
  dplyr::group_by(Region, Site) %>%
  dplyr::summarise(
    glacier = dplyr::first(glacier),
    GI = mean(GI),
    .groups = "drop"
  )
site_a <- site_metadata %>%
  dplyr::transmute(
    Region,
    SiteA = Site,
    glacier_A = glacier,
    GI_A = GI
  )
site_b <- site_metadata %>%
  dplyr::transmute(
    Region,
    SiteB = Site,
    glacier_B = glacier,
    GI_B = GI
  )

site_pairs <- site_pairs %>%
  dplyr::left_join(site_a, by = c("Region", "SiteA")) %>%
  dplyr::left_join(site_b, by = c("Region", "SiteB")) %>%
  dplyr::mutate(
    same_glacier = glacier_A == glacier_B,
    GI_mean = (GI_A + GI_B) / 2,
    GI_difference = abs(GI_A - GI_B),
    Selection_total = HoS + HeS,
    total_process_fraction = HoS + HeS + HD + DL + DR
  ) %>%
  dplyr::relocate(
    Region, SiteA, SiteB, glacier_A, glacier_B, same_glacier,
    GI_A, GI_B, GI_mean, GI_difference, n_replicate_pairs
  )

# --- Prepare and validate site-pair data ---

# Each row is an unordered distinct-site pair averaged across replicate pairs.
site_pairs <- site_pairs %>%
  dplyr::mutate(
    Region = factor(Region),
    SiteA = as.character(SiteA),
    SiteB = as.character(SiteB),
    same_glacier = factor(as.logical(same_glacier), levels = c(FALSE, TRUE)),
    dplyr::across(
      c(GI_A, GI_B, GI_mean, GI_difference, dplyr::all_of(processes)),
      as.numeric
    )
  )

pair_count_audit <- site_pairs %>%
  dplyr::group_by(Region) %>%
  dplyr::summarise(
    n_sites = dplyr::n_distinct(c(SiteA, SiteB)),
    observed_pairs = dplyr::n(),
    expected_pairs = choose(n_sites, 2),
    .groups = "drop"
  )
if (any(pair_count_audit$observed_pairs != pair_count_audit$expected_pairs)) {
  # Node relabelling requires a complete within-region pair network.
  stop("The input is not a complete site-pair network within each Region.", call. = FALSE)
}

# Standardize GI predictors so coefficients represent a one-SD change.
GI_mean_center <- mean(site_pairs$GI_mean)
GI_mean_scale <- stats::sd(site_pairs$GI_mean)
GI_difference_center <- mean(site_pairs$GI_difference)
GI_difference_scale <- stats::sd(site_pairs$GI_difference)

site_pairs <- site_pairs %>%
  dplyr::mutate(
    GI_mean_z = (GI_mean - GI_mean_center) / GI_mean_scale,
    GI_difference_z = (GI_difference - GI_difference_center) / GI_difference_scale
  )

data_summary <- site_pairs %>%
  tidyr::pivot_longer(
    dplyr::all_of(processes),
    names_to = "process",
    values_to = "fraction"
  ) %>%
  dplyr::group_by(process) %>%
  dplyr::summarise(
    n_site_pairs = dplyr::n(),
    mean_fraction = mean(fraction),
    sd_fraction = stats::sd(fraction),
    min_fraction = min(fraction),
    max_fraction = max(fraction),
    .groups = "drop"
  )

print(pair_count_audit)
print(data_summary)

# The five fractions sum to one, so an increase in one process reflects
# redistribution among processes rather than independent responses.

# --- Joint GI model with MRQAP inference ---

# Use the same additive Region structure as the preceding analysis; each GI term
# is adjusted for the other and for same-glacier status.
X_full <- stats::model.matrix(
  ~ Region + same_glacier + GI_mean_z + GI_difference_z,
  data = site_pairs
)
focal_terms <- c(
  GI_mean = "GI_mean_z",
  GI_difference = "GI_difference_z"
)
# Reuse one permutation set across all processes and both GI tests. Synchronously
# permute node labels within each glacier to preserve Region, glacier membership,
# dyadic dependence and glacier-pair structure.
# Collapse the two endpoint columns into one site-level membership table.
endpoint_strata <- dplyr::bind_rows(
  site_pairs %>%
    dplyr::transmute(
      Region = as.character(Region),
      Site = as.character(SiteA),
      glacier = as.character(glacier_A)
    ),
  site_pairs %>%
    dplyr::transmute(
      Region = as.character(Region),
      Site = as.character(SiteB),
      glacier = as.character(glacier_B)
    )
) %>%
  dplyr::distinct()

regions <- unique(endpoint_strata$Region)
sites_by_region_glacier <- lapply(
  stats::setNames(regions, regions),
  function(region_i) {
    region_strata <- endpoint_strata[
      endpoint_strata$Region == region_i,
      ,
      drop = FALSE
    ]
    split(region_strata$Site, region_strata$glacier)
  }
)

set.seed(as.integer(params$seed))
# Each vector maps sites to permuted labels within the same glacier.
node_permutations <- replicate(
  as.integer(params$n_permutations),
  lapply(sites_by_region_glacier, function(glacier_groups) {
    glacier_maps <- lapply(
      glacier_groups,
      function(nodes) stats::setNames(sample(nodes), nodes)
    )
    do.call(c, unname(glacier_maps))
  }),
  simplify = FALSE
)

# --- Repeated MRQAP helpers ---

pair_key <- function(region, site1, site2) {
  # Site-pair orientation is arbitrary, so use a canonical unordered key.
  paste(region, pmin(site1, site2), pmax(site1, site2), sep = "::")
}

squish_fraction <- function(x, eps = 1e-6) {
  # Keep the logit finite if a future input contains an exact zero or one.
  pmin(pmax(as.numeric(x), eps), 1 - eps)
}

permute_edge_values <- function(values, data, node_permutation) {
  # Relabel the dyadic network without treating pair rows as independent.
  permuted <- numeric(length(values))
  for (region_i in names(node_permutation)) {
    index <- which(data$Region == region_i)
    node_map <- node_permutation[[region_i]]
    lookup <- stats::setNames(
      values[index],
      pair_key(data$Region[index], data$SiteA[index], data$SiteB[index])
    )
    mapped_key <- pair_key(
      data$Region[index],
      unname(node_map[data$SiteA[index]]),
      unname(node_map[data$SiteB[index]])
    )
    permuted[index] <- unname(lookup[mapped_key])
  }
  permuted
}

fit_lm_matrix <- function(X, y) {
  # Fit from the design matrix to avoid repeated formula parsing.
  fit <- stats::lm.fit(x = X, y = y)
  residual_df <- nrow(X) - fit$rank
  sigma2 <- sum(fit$residuals^2) / residual_df
  rank_index <- seq_len(fit$rank)
  R <- qr.R(fit$qr)[rank_index, rank_index, drop = FALSE]
  # Restore QR-pivoted standard errors to the original matrix-column order.
  standard_error <- rep(NA_real_, ncol(X))
  standard_error[fit$qr$pivot[rank_index]] <- sqrt(
    diag(chol2inv(R)) * sigma2
  )
  list(
    coefficients = stats::setNames(fit$coefficients, colnames(X)),
    standard_error = stats::setNames(standard_error, colnames(X)),
    statistic = stats::setNames(fit$coefficients / standard_error, colnames(X)),
    fitted = fit$fitted.values,
    residuals = fit$residuals,
    r_squared = 1 - sum(fit$residuals^2) / sum((y - mean(y))^2),
    residual_df = residual_df
  )
}

freedman_lane_test <- function(
  y,
  X_full,
  X_reduced,
  focal_term,
  data,
  permutations
) {
  # Add permuted reduced-model residuals to the reduced fitted values.
  observed <- fit_lm_matrix(X_full, y)
  reduced <- fit_lm_matrix(X_reduced, y)
  permuted_statistics <- vapply(permutations, function(node_permutation) {
    y_permuted <- reduced$fitted + permute_edge_values(
      reduced$residuals, data, node_permutation
    )
    fit_lm_matrix(X_full, y_permuted)$statistic[[focal_term]]
  }, numeric(1))

  observed_statistic <- observed$statistic[[focal_term]]
  tibble::tibble(
    estimate_logit_per_sd = observed$coefficients[[focal_term]],
    model_based_se = observed$standard_error[[focal_term]],
    t_observed = observed_statistic,
    # Use two-sided absolute t values and a plus-one P-value correction.
    p_permutation = (1 + sum(abs(permuted_statistics) >= abs(observed_statistic))) /
      (length(permuted_statistics) + 1),
    r_squared = observed$r_squared,
    residual_df = observed$residual_df
  )
}

model_results <- dplyr::bind_rows(
  lapply(processes, function(process_i) {
    # Express each fraction relative to all other processes; bound zero and one
    # to keep the logit finite.
    response_logit <- stats::qlogis(squish_fraction(site_pairs[[process_i]]))
    dplyr::bind_rows(
      lapply(names(focal_terms), function(predictor_i) {
        focal_term <- unname(focal_terms[[predictor_i]])
        # Remove only the focal GI term; retain Region, same-glacier status and
        # the other GI predictor in the reduced model.
        X_reduced <- X_full[, colnames(X_full) != focal_term, drop = FALSE]
        freedman_lane_test(
          y = response_logit,
          X_full = X_full,
          X_reduced = X_reduced,
          focal_term = focal_term,
          data = site_pairs,
          permutations = node_permutations
        ) %>%
          dplyr::mutate(
            process = process_i,
            predictor = predictor_i,
            .before = 1
          )
      })
    )
  })
) %>%
  # Apply BH correction across processes within each GI predictor.
  dplyr::group_by(predictor) %>%
  dplyr::mutate(
    q_BH_within_predictor = stats::p.adjust(p_permutation, method = "BH")
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    odds_ratio_per_sd = exp(estimate_logit_per_sd),
    n_sites = dplyr::n_distinct(c(site_pairs$SiteA, site_pairs$SiteB)),
    n_glaciers = dplyr::n_distinct(c(site_pairs$glacier_A, site_pairs$glacier_B)),
    n_site_pairs = nrow(site_pairs),
    n_permutations = as.integer(params$n_permutations),
    permutation_resolution = 1 / (n_permutations + 1),
    permutation_scheme = "Node labels constrained within exact glacier",
    model = "logit(fraction) ~ Region + same_glacier + GI_mean_z + GI_difference_z"
  ) %>%
  dplyr::arrange(predictor, match(process, processes))

# --- Result tables ---

paper_table <- model_results %>%
  dplyr::select(
    process, predictor, estimate_logit_per_sd, odds_ratio_per_sd,
    p_permutation, q_BH_within_predictor, n_sites, n_glaciers,
    n_site_pairs, n_permutations
  )

extended_process_labels <- c(
  HoS = "Homogeneous selection (HoS)",
  HeS = "Heterogeneous selection (HeS)",
  HD = "Homogenizing dispersal (HD)",
  DL = "Dispersal limitation (DL)",
  DR = "Drift and other undominated processes (DR)"
)

extended_predictor_labels <- c(
  GI_mean = "pairwise mean GI",
  GI_difference = "absolute GI difference"
)

# Extended Data table for both prespecified GI predictors. Keep the existing
# schema and identify the predictor within the process label.
# Because inference uses MRQAP, report observed t but omit model-based SE.
extended_data_table <- model_results %>%
  dplyr::filter(predictor %in% names(extended_predictor_labels)) %>%
  dplyr::mutate(
    predictor_order = match(predictor, names(extended_predictor_labels)),
    process_order = match(process, processes)
  ) %>%
  dplyr::arrange(predictor_order, process_order) %>%
  dplyr::transmute(
    `Assembly process` = paste0(
      unname(extended_process_labels[process]), ": ",
      unname(extended_predictor_labels[predictor])
    ),
    `Coefficient (logit scale per predictor SD)` = estimate_logit_per_sd,
    `Odds ratio per predictor SD` = odds_ratio_per_sd,
    `Observed t statistic` = t_observed,
    `MRQAP P value` = p_permutation,
    `BH q value (within predictor)` = q_BH_within_predictor,
    `Model R2` = r_squared,
    `Number of sites` = n_sites,
    `Number of glaciers` = n_glaciers,
    `Number of site pairs` = n_site_pairs,
    `Number of permutations` = n_permutations,
    `Permutation resolution` = permutation_resolution,
    `Random seed` = as.integer(params$seed),
    `Permutation scheme` = permutation_scheme,
    `Model specification` = model
  )

print(paper_table)
print(extended_data_table)

# `GI_mean` represents position along the GI gradient; `GI_difference` represents
# separation between pair endpoints. Coefficients are standardized and mutually
# adjusted, with permutation P values as the primary inference.

# --- Main figure: pairwise HoS along mean GI ---

# Plot observed site pairs and joint-model predictions for both regions.
# Predictions set GI_difference to its mean (z = 0) and same_glacier to FALSE,
# the dominant pair type. Confidence ribbons are omitted because inference uses
# glacier-constrained MRQAP.
region_palette <- c(Alps = "#E69F00", Kyrgyzstan = "#56B4E9")

hos_response_logit <- stats::qlogis(squish_fraction(site_pairs$HoS))
hos_full_fit <- fit_lm_matrix(X_full, hos_response_logit)

hos_prediction_grid <- dplyr::bind_rows(
  lapply(levels(site_pairs$Region), function(region_i) {
    region_data <- site_pairs %>%
      dplyr::filter(Region == region_i)
    gi_sequence <- seq(
      min(region_data$GI_mean),
      max(region_data$GI_mean),
      length.out = 200
    )
    new_data <- tibble::tibble(
      Region = factor(region_i, levels = levels(site_pairs$Region)),
      same_glacier = factor(FALSE, levels = levels(site_pairs$same_glacier)),
      GI_mean = gi_sequence,
      GI_mean_z = (gi_sequence - GI_mean_center) / GI_mean_scale,
      GI_difference_z = 0
    )
    X_prediction <- stats::model.matrix(
      ~ Region + same_glacier + GI_mean_z + GI_difference_z,
      data = new_data
    )
    X_prediction <- X_prediction[, colnames(X_full), drop = FALSE]
    new_data %>%
      dplyr::mutate(
        pred_HoS = stats::plogis(
          drop(X_prediction %*% hos_full_fit$coefficients)
        )
      )
  })
)

hos_gi_result <- paper_table %>%
  dplyr::filter(process == "HoS", predictor == "GI_mean")

format_permutation_p <- function(x) {
  ifelse(
    x < 0.001,
    formatC(x, format = "f", digits = 4),
    formatC(x, format = "f", digits = 3)
  )
}

hos_stat_label <- sprintf(
  "beta = %.3f per SD\nMRQAP P = %s; BH q = %s\n49 sites; 19 glaciers; 588 site pairs",
  hos_gi_result$estimate_logit_per_sd,
  format_permutation_p(hos_gi_result$p_permutation),
  format_permutation_p(hos_gi_result$q_BH_within_predictor)
)

pairwise_hos_main_plot <- ggplot2::ggplot(
  site_pairs,
  ggplot2::aes(x = GI_mean, y = HoS, colour = Region)
) +
  ggplot2::geom_point(size = 1.2, alpha = 0.24) +
  ggplot2::geom_line(
    data = hos_prediction_grid,
    ggplot2::aes(x = GI_mean, y = pred_HoS, colour = Region),
    linewidth = 1.25,
    inherit.aes = FALSE
  ) +
  ggplot2::annotate(
    "label",
    x = min(site_pairs$GI_mean),
    y = max(site_pairs$HoS),
    label = hos_stat_label,
    hjust = 1,
    vjust = 1,
    size = 3.5,
    lineheight = 1.1,
    linewidth = 0,
    fill = "white",
    alpha = 0.9,
    colour = "black"
  ) +
  ggplot2::scale_colour_manual(values = region_palette) +
  ggplot2::scale_x_reverse() +
  ggplot2::labs(
    x = "Pairwise mean Glacial Index (GI)",
    y = "Homogeneous selection fraction",
    colour = "Region"
  ) +
  ggplot2::guides(
    colour = ggplot2::guide_legend(
      override.aes = list(alpha = 1, linewidth = 1.1, size = 2.4)
    )
  ) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(
    legend.position = "top",
    legend.justification = "center",
    panel.grid.minor = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(7, 9, 7, 7)
  )

# --- Save outputs ---

if (isTRUE(params$write_outputs)) {
  # Write tables and the PDF figure without modifying upstream iCAMP results.
  readr::write_csv(site_pairs, default_site_pair_path)
  readr::write_csv(
    model_results,
    file.path(out_dir, "extended_data_table_3_complete_mrqap_results.csv")
  )
  readr::write_csv(
    paper_table,
    file.path(out_dir, "figure_2f_mrqap_paper_table.csv")
  )
  readr::write_csv(
    extended_data_table,
    file.path(out_dir, "extended_data_table_3.csv")
  )
  ggplot2::ggsave(
    file.path(out_dir, "figure_2f_homogeneous_selection_vs_gi.pdf"),
    pairwise_hos_main_plot,
    width = 7,
    height = 5.6,
    units = "in"
  )
  readr::write_csv(
    tibble::tibble(
      parameter = c("n_permutations", "seed", "permutation_scheme"),
      value = c(
        as.character(params$n_permutations),
        as.character(params$seed),
        "Freedman-Lane node permutations constrained within exact glacier"
      )
    ),
    file.path(out_dir, "icamp_mrqap_parameters.csv")
  )
  writeLines(
    capture.output(sessionInfo()),
    file.path(out_dir, "icamp_mrqap_session_info.txt")
  )
}

# Manuscript unit: 49 sediment sites in 19 glaciers and 588 within-Region pairs.
# MRQAP treats pairs as dyads sharing 49 nodes and preserves glacier membership.
