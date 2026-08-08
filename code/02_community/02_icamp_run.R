# Main server-side iCAMP run using sediment samples only
#
# This script computes iCAMP and saves the standardized replicate-pair result:
#   results/02_community/icamp/sediment_pairwise_process_fractions.csv
#
# Run from the repository root:
#   Rscript code/02_community/02_icamp_run.R

# --- 1. Package requirements ---

required_packages <- c(
  "ape",
  "dplyr",
  "tidyr",
  "purrr",
  "readr",
  "tibble",
  "stringr",
  "rlang",
  "phyloseq",
  "iCAMP"
)

# Package versions are recorded in renv.lock; loading them here also fails
# immediately when the active R environment is incomplete.
invisible(lapply(required_packages, library, character.only = TRUE))

# --- 2. Runtime configuration and project paths ---

# Allow runtime settings to be overridden without editing the script.
env_int <- function(name, default) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || value == "")
    return(default)
  as.integer(value)
}

env_num <- function(name, default) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || value == "")
    return(default)
  as.numeric(value)
}

project_root <- normalizePath(
  Sys.getenv("PAMIR_PAPER_ROOT", unset = "."),
  winslash = "/",
  mustWork = TRUE
)

# Keep compact workflow outputs separate from large, regenerable iCAMP objects.
ps_path <- file.path(
  project_root, "results", "01_16s", "intermediate", "PAMIR_16S_final.rds"
)
tree_path <- file.path(project_root, "data", "processed", "16s", "dna-sequences.tree")
out_dir <- file.path(project_root, "results", "02_community", "icamp")
generated_dir <- file.path(
  project_root, "results", "generated", "02_community", "icamp"
)
raw_dir <- file.path(generated_dir, "raw")
pd_dir <- file.path(generated_dir, "phylogenetic_distance")
pairwise_output_file <- file.path(out_dir, "sediment_pairwise_process_fractions.csv")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pd_dir, recursive = TRUE, showWarnings = FALSE)

asv_min_prevalence <- env_int("ASV_MIN_PREVALENCE", 3)
asv_min_total_count <- env_int("ASV_MIN_TOTAL_COUNT", 10)
asv_min_prevalence_fraction <- env_num("ASV_MIN_PREVALENCE_FRACTION", 0.02)
icamp_rand <- env_int("ICAMP_RAND", 1000)
icamp_nworker <- env_int("ICAMP_NWORKER", 48)
icamp_memory_G <- env_num("ICAMP_MEMORY_G", 300)
icamp_input_type <- "unrarefied_counts"
icamp_prefix <- "PAMIR_sediment_iCAMP"
icamp_seed <- 20260416L
# Fix the random seed used by iCAMP's null-model randomizations.
set.seed(icamp_seed)

readr::write_csv(
  tibble::tibble(
    parameter = c(
      "asv_min_prevalence", "asv_min_total_count",
      "asv_min_prevalence_fraction", "randomizations", "seed",
      "input_type", "nworker", "memory_G"
    ),
    value = as.character(c(
      asv_min_prevalence, asv_min_total_count,
      asv_min_prevalence_fraction, icamp_rand, icamp_seed,
      icamp_input_type, icamp_nworker, icamp_memory_G
    ))
  ),
  file.path(out_dir, "icamp_parameters.csv")
)

# --- 3. Validation and iCAMP-output helpers ---

# Small validation helpers keep failures early and explicit, before the
# expensive iCAMP step starts.
stop_if_missing <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0) {
    stop(
      "Missing required input file(s):\n",
      paste(missing, collapse = "\n"),
      call. = FALSE
    )
  }
}

require_columns <- function(df, cols, df_name) {
  missing <- setdiff(cols, names(df))
  if (length(missing) > 0) {
    stop(
      df_name,
      " is missing required column(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
}

to_numeric_safely <- function(x) {
  # Some iCAMP tables store numeric fractions as factors or character strings.
  if (is.numeric(x))
    return(x)
  suppressWarnings(as.numeric(as.character(x)))
}

sanitize_file_label <- function(x) {
  # Region names are reused in directory names and iCAMP file prefixes.
  gsub("[^A-Za-z0-9_]+", "_", x)
}

process_label <- function(x) {
  # Map version-dependent process labels to the five abbreviations used here.
  z <- stringr::str_to_lower(as.character(x))
  z <- stringr::str_replace_all(z, "[^a-z]+", "_")
  dplyr::case_when(
    stringr::str_detect(
      z,
      "homogeneous_selection|homo_selection|^hos$|hselection|homogeneous"
    ) ~ "HoS",
    stringr::str_detect(
      z,
      "heterogeneous_selection|variable_selection|var_selection|^hes$|^vs$|vselection|heterogeneous"
    ) ~ "HeS",
    stringr::str_detect(z, "homogenizing_dispersal|homo_dispersal|^hd$") ~ "HD",
    stringr::str_detect(z, "dispersal_limitation|^dl$|limitation") ~ "DL",
    stringr::str_detect(z, "drift|undominated|dominated|^dr$|^drift$") ~ "DR",
    TRUE ~ NA_character_
  )
}

# iCAMP output tables may omit processes with zero contribution for a pair.
# Completing the five-process grid prevents missing categories from being
# mistaken for missing samples during later averaging.
complete_process_grid <- function(df, id_cols) {
  process_levels <- c("HoS", "HeS", "HD", "DL", "DR")
  df %>%
    tidyr::complete(
      tidyr::nesting(!!!rlang::syms(id_cols)),
      Process = process_levels,
      fill = list(Fraction = 0)
    )
}

normalize_pairwise_fractions <- function(df) {
  # Accept proportions, percentages, or unscaled positive weights and convert
  # every sample pair to proportions summing to one.
  df %>%
    dplyr::group_by(Region, Sample1, Sample2) %>%
    dplyr::mutate(
      pair_total_raw = sum(Fraction, na.rm = TRUE),
      Fraction = dplyr::case_when(
        is.na(pair_total_raw) | pair_total_raw <= 0 ~ NA_real_,
        dplyr::between(pair_total_raw, 95, 105) ~ Fraction / 100,
        abs(pair_total_raw - 1) <= 0.05 ~ Fraction,
        TRUE ~ Fraction / pair_total_raw
      )
    ) %>%
    dplyr::ungroup()
}

audit_pairwise_orientation <- function(df) {
  # Pair order is arbitrary downstream, so A-B and B-A must not coexist.
  pair_orientation <- df %>%
    dplyr::distinct(Region, Sample1, Sample2) %>%
    dplyr::mutate(
      pair_min = pmin(Sample1, Sample2),
      pair_max = pmax(Sample1, Sample2)
    ) %>%
    dplyr::count(Region, pair_min, pair_max, name = "n_orientations")

  duplicated_pairs <- pair_orientation %>%
    dplyr::filter(n_orientations > 1)

  if (nrow(duplicated_pairs) > 0) {
    stop(
      "The extracted iCAMP pairwise table contains both sample-pair orientations ",
      "for at least one region-specific pair (for example A-B and B-A). ",
      "This script assumes one unordered pair per comparison when averaging ",
      "replicate-level process fractions, so please deduplicate the pairwise ",
      "table explicitly before continuing.",
      call. = FALSE
    )
  }

  invisible(pair_orientation)
}

# Different iCAMP versions expose slightly different formal arguments.
# This wrapper keeps the script portable while still passing all supported
# parameters on installations that accept them.
call_iCAMP_function <- function(fun, args, fun_name) {
  formal_names <- names(formals(fun))
  if (!("..." %in% formal_names)) {
    unsupported <- setdiff(names(args), formal_names)
    if (length(unsupported) > 0) {
      message(
        fun_name,
        " does not support argument(s): ",
        paste(unsupported, collapse = ", "),
        ". Skipping them."
      )
    }
    args <- args[intersect(names(args), formal_names)]
  }
  do.call(fun, args)
}

collect_data_frames <- function(x, path = "result") {
  # Recursively flatten tabular components while retaining their object paths.
  out <- list()
  if (is.data.frame(x)) {
    out[[path]] <- x
  } else if (is.matrix(x)) {
    out[[path]] <- as.data.frame(x, check.names = FALSE)
  } else if (is.list(x)) {
    for (nm in names(x)) {
      child_path <- paste(path, nm, sep = "$")
      out <- c(out, collect_data_frames(x[[nm]], child_path))
    }
  }
  out
}

first_existing <- function(nms, candidates) {
  # Candidate order defines which known column alias takes precedence.
  candidates[candidates %in% nms][1]
}

extract_pairwise_process <- function(icamp_result, region_name, valid_samples) {
  # Search all nested tables because their locations differ among iCAMP builds.
  dfs <- collect_data_frames(icamp_result)

  # Search nested iCAMP result tables for a sample-pair process table. The
  # extractor accepts both long and wide outputs because package versions vary.
  for (df_path in names(dfs)) {
    df <- dfs[[df_path]]
    if (nrow(df) == 0 || ncol(df) < 3)
      next

    original_names <- names(df)
    # Normalize names only for matching; preserve the originals for wide tables.
    normalized_names <- stringr::str_to_lower(gsub("[^A-Za-z0-9]+", "_", original_names))
    names(df) <- normalized_names

    sample1_col <- first_existing(
      names(df),
      c(
        "sample1", "sample_1", "name1", "comm1", "site1", "site_1",
        "sampleid1", "sample_id1", "id1", "id_1"
      )
    )
    sample2_col <- first_existing(
      names(df),
      c(
        "sample2", "sample_2", "name2", "comm2", "site2", "site_2",
        "sampleid2", "sample_id2", "id2", "id_2"
      )
    )

    if (is.na(sample1_col) || is.na(sample2_col))
      next

    candidate <- df %>%
      dplyr::mutate(
        Sample1 = as.character(.data[[sample1_col]]),
        Sample2 = as.character(.data[[sample2_col]])
      ) %>%
      dplyr::filter(Sample1 %in% valid_samples, Sample2 %in% valid_samples, Sample1 != Sample2)

    if (nrow(candidate) == 0)
      next

    process_col <- first_existing(
      names(candidate),
      c("process", "processes", "process_category", "processes_category", "process_name")
    )
    fraction_col <- first_existing(
      names(candidate),
      c(
        "fraction", "fractions", "frac", "frc", "weight", "wt",
        "contribution", "percentage", "relative_importance"
      )
    )

    if (!is.na(process_col) && !is.na(fraction_col)) {
      # Preferred long format: one process and one fraction per row.
      extracted <- candidate %>%
        dplyr::transmute(
          Region = region_name,
          Sample1,
          Sample2,
          Process = process_label(.data[[process_col]]),
          Fraction = to_numeric_safely(.data[[fraction_col]])
        ) %>%
        dplyr::filter(!is.na(Process), !is.na(Fraction))

      if (nrow(extracted) > 0)
        return(extracted)
    }

    wide_process_cols <- original_names[!is.na(process_label(original_names))]
    wide_process_cols_norm <- normalized_names[match(wide_process_cols, original_names)]

    if (length(wide_process_cols_norm) >= 3) {
      # Alternative wide format: process fractions occupy separate columns.
      extracted <- candidate %>%
        dplyr::select(Sample1, Sample2, dplyr::all_of(wide_process_cols_norm)) %>%
        tidyr::pivot_longer(
          cols = dplyr::all_of(wide_process_cols_norm),
          names_to = "raw_process",
          values_to = "Fraction"
        ) %>%
        dplyr::mutate(
          Region = region_name,
          Process = process_label(raw_process),
          Fraction = to_numeric_safely(Fraction)
        ) %>%
        dplyr::select(Region, Sample1, Sample2, Process, Fraction) %>%
        dplyr::filter(!is.na(Process), !is.na(Fraction))

      if (nrow(extracted) > 0)
        return(extracted)
    }
  }

  stop(
    "Could not identify a pairwise process-fraction table in the iCAMP result for region ",
    region_name,
    call. = FALSE
  )
}

# --- 4. Load sediment samples and construct the count matrix ---

stop_if_missing(c(ps_path, tree_path))

message("Project root: ", project_root)
message("Pairwise output file: ", pairwise_output_file)
message("iCAMP rand=", icamp_rand, ", nworker=", icamp_nworker, ", memory.G=", icamp_memory_G)

ps <- readRDS(ps_path)
# The cleaned 16S object retains only identifiers required for sediment
# selection and regional iCAMP pools; site, glacier, and GI metadata are
# attached downstream in 03_icamp_analysis.R.
required_sample_vars <- c(
  "Sample_ID", "Source", "Region"
)
sample_data_df <- data.frame(phyloseq::sample_data(ps), check.names = FALSE)
require_columns(sample_data_df, required_sample_vars, "sample_data(PAMIR_16S_final.rds)")

# Phyloseq row names are internal keys; Sample_ID is the identifier written to
# iCAMP inputs and final pairwise outputs.
metadata_all <- sample_data_df %>%
  tibble::rownames_to_column("phyloseq_sample_name") %>%
  dplyr::mutate(
    Sample_ID = as.character(Sample_ID),
    Source = as.character(Source),
    Region = as.factor(Region)
  )

sed_sample_ids <- metadata_all %>%
  dplyr::filter(Source == "sed") %>%
  dplyr::pull(phyloseq_sample_name)

if (length(sed_sample_ids) == 0) {
  stop("No sediment samples found with `Source == 'sed'`.", call. = FALSE)
}

ps_sed <- phyloseq::prune_samples(sed_sample_ids, ps)
otu <- as(phyloseq::otu_table(ps_sed), "matrix")
# Normalize the OTU table orientation to ASVs in rows and samples in columns.
if (!phyloseq::taxa_are_rows(ps_sed)) {
  otu <- t(otu)
}

sample_meta_sed <- metadata_all %>%
  dplyr::filter(phyloseq_sample_name %in% sed_sample_ids) %>%
  dplyr::select(
    phyloseq_sample_name, Sample_ID, Region
  )

if (anyDuplicated(sample_meta_sed$Sample_ID) > 0) {
  stop("Sediment Sample_ID values are not unique.", call. = FALSE)
}

# --- 5. Filter ASVs and align the phylogenetic tree ---

# Filtering is applied once across all sediment samples. Region-specific taxa
# absent from a regional pool are removed later inside run_icamp_region().
tree <- ape::read.tree(tree_path)
prevalence <- rowSums(otu > 0)
total_count <- rowSums(otu)
# Apply the stricter of the absolute and proportional prevalence thresholds.
prevalence_threshold <- max(asv_min_prevalence, ceiling(ncol(otu) * asv_min_prevalence_fraction))
keep_asv <- prevalence >= prevalence_threshold & total_count >= asv_min_total_count
otu_filt <- otu[keep_asv, , drop = FALSE]

# Keep only ASVs retained by abundance/prevalence filtering and represented in
# the phylogenetic tree used by iCAMP.
common_asv <- intersect(rownames(otu_filt), tree$tip.label)
if (length(common_asv) < 100) {
  stop("Fewer than 100 filtered ASVs overlap with tree tip labels.", call. = FALSE)
}

otu_filt <- otu_filt[common_asv, , drop = FALSE]
tree_filt <- ape::keep.tip(tree, common_asv)

# iCAMP expects communities in rows and taxa in columns.
otu_comm_count <- t(otu_filt)
rownames(otu_comm_count) <- sample_meta_sed$Sample_ID[
  match(rownames(otu_comm_count), sample_meta_sed$phyloseq_sample_name)
]

if (any(is.na(rownames(otu_comm_count)))) {
  stop("Could not map all phyloseq sample names to Sample_ID after filtering.", call. = FALSE)
}

otu_comm_count <- otu_comm_count[, tree_filt$tip.label, drop = FALSE]

message("Sediment samples: ", nrow(otu_comm_count))
message("Filtered ASVs after tree alignment: ", ncol(otu_comm_count))
message(
  "Prevalence threshold: ", prevalence_threshold,
  " samples; total count threshold: ", asv_min_total_count
)

# --- 6. Run and preserve regional iCAMP analyses ---

run_icamp_region <- function(region_name, comm_count, meta, tree_obj) {
  region_label <- sanitize_file_label(region_name)
  # Isolate regional distance files and raw results to prevent name collisions.
  pd_region_dir <- file.path(pd_dir, region_label)
  region_output_dir <- file.path(raw_dir, region_label)
  dir.create(pd_region_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(region_output_dir, recursive = TRUE, showWarnings = FALSE)

  region_samples <- meta %>%
    dplyr::filter(as.character(Region) == region_name) %>%
    dplyr::pull(Sample_ID)

  comm_region <- comm_count[region_samples, , drop = FALSE]
  # Remove taxa absent from the current regional species pool.
  comm_region <- comm_region[, colSums(comm_region > 0) > 0, drop = FALSE]
  run_prefix <- paste0(icamp_prefix, "_", icamp_input_type, "_", region_label)

  if (nrow(comm_region) < 6) {
    stop(
      "Region ", region_name,
      " has fewer than 6 sediment samples after filtering.",
      call. = FALSE
    )
  }
  if (ncol(comm_region) < 100) {
    stop(
      "Region ", region_name,
      " has fewer than 100 ASVs after filtering/alignment.",
      call. = FALSE
    )
  }
  # Prune the phylogeny to exactly the taxa retained for this regional run.
  tree_region <- ape::keep.tip(tree_obj, colnames(comm_region))

  message(
    "Running iCAMP for ", region_name,
    ": samples=", nrow(comm_region),
    ", ASVs=", ncol(comm_region)
  )

  icamp_result <- call_iCAMP_function(
    iCAMP::icamp.big,
    list(
      comm = as.data.frame(comm_region, check.names = FALSE),
      # Main run uses unrarefied counts; unit.sum supplies sample library sizes
      # for iCAMP's sample-size normalization.
      unit.sum = rowSums(comm_region),
      tree = tree_region,
      pd.wd = pd_region_dir,
      rand = icamp_rand,
      prefix = run_prefix,
      ds = 0.2,
      pd.cut = NA,
      sp.check = TRUE,
      phylo.rand.scale = "within.bin",
      taxa.rand.scale = "across.all",
      sig.index = "SES.RC",
      bin.size.limit = 24,
      nworker = icamp_nworker,
      memory.G = icamp_memory_G,
      detail.save = FALSE,
      qp.save = FALSE,
      detail.null = FALSE,
      output.wd = region_output_dir,
      pd.backingfile = paste0(run_prefix, "_pd.bin"),
      pd.desc.file = paste0(run_prefix, "_pd.desc"),
      pd.spname.file = paste0(run_prefix, "_pd_taxon_name.csv"),
      treepath.file = paste0(run_prefix, "_tree_path.rda")
    ),
    "icamp.big"
  )

  # Preserve the complete regional result returned by iCAMP so the standardized
  # pairwise CSV can be traced back to an explicit, reusable raw result object.
  icamp_result_path <- file.path(
    region_output_dir,
    paste0(run_prefix, "_icamp_result.rds")
  )
  saveRDS(icamp_result, icamp_result_path)
  message("Saved regional iCAMP result: ", icamp_result_path)

  icamp_result
}

regions <- sort(unique(as.character(sample_meta_sed$Region)))

# Run iCAMP separately within each region so null processes are estimated
# within, not across, the Alps/Kyrgyzstan divide.
icamp_results <- purrr::set_names(
  purrr::map(
    regions,
    run_icamp_region,
    comm_count = otu_comm_count,
    meta = sample_meta_sed,
    tree_obj = tree_filt
  ),
  regions
)

# --- 7. Standardize and save replicate-pair process fractions ---

# Convert version-dependent regional iCAMP result structures into one common
# five-process table before validating and writing the downstream input.
# Rows resolving to the same pair and process are collapsed before omitted
# zero-contribution processes are restored.
pairwise_process <- purrr::imap_dfr(
  icamp_results,
  ~ extract_pairwise_process(
    icamp_result = .x,
    region_name = .y,
    valid_samples = sample_meta_sed %>%
      dplyr::filter(as.character(Region) == .y) %>%
      dplyr::pull(Sample_ID)
  )
) %>%
  dplyr::group_by(Region, Sample1, Sample2, Process) %>%
  dplyr::summarise(Fraction = sum(Fraction, na.rm = TRUE), .groups = "drop") %>%
  complete_process_grid(id_cols = c("Region", "Sample1", "Sample2")) %>%
  normalize_pairwise_fractions()

audit_pairwise_orientation(pairwise_process)

# Save one row per unordered replicate pair. Endpoint- and site-level
# summaries are derived downstream in 03_icamp_analysis.R.
# The wide table also records selection and all-process totals for downstream
# analysis and final integrity checks.
pairwise_process_complete <- pairwise_process %>%
  dplyr::select(Region, Sample1, Sample2, Process, Fraction) %>%
  tidyr::pivot_wider(names_from = Process, values_from = Fraction, values_fill = 0) %>%
  dplyr::mutate(
    Selection_total = HoS + HeS,
    total_process_fraction = HoS + HeS + HD + DL + DR
  )

process_cols <- c("HoS", "HeS", "HD", "DL", "DR")
# Allow only negligible floating-point excursions beyond the valid range.
if (any(as.matrix(pairwise_process_complete[process_cols]) < -1e-8 |
        as.matrix(pairwise_process_complete[process_cols]) > 1 + 1e-8, na.rm = TRUE)) {
  stop("Some pairwise sediment process fractions are outside [0,1].", call. = FALSE)
}
if (max(abs(pairwise_process_complete$total_process_fraction - 1), na.rm = TRUE) > 0.01) {
  stop("Pairwise sediment process fractions do not sum close to 1.", call. = FALSE)
}

readr::write_csv(pairwise_process_complete, pairwise_output_file)
# Record the exact R and package environment alongside the standardized result.
writeLines(
  capture.output(sessionInfo()),
  file.path(out_dir, "icamp_session_info.txt")
)

message("Saved standardized replicate-pair sediment iCAMP result: ", pairwise_output_file)
message("Summary: n_replicate_pairs=", nrow(pairwise_process_complete))
