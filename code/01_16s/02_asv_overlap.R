# ASV overlap and taxonomic composition analysis
# Manuscript outputs: Figure 1a and Extended Data Figures 1-2.

rm(list = ls())
gc()
graphics.off()

input_path <- file.path(
  "results", "01_16s", "intermediate", "PAMIR_16S_final.rds"
)
result_dir <- file.path("results", "01_16s", "asv_overlap")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

# --- Setup ---

library(phyloseq)
library(ggplot2)
library(VennDiagram)
library(tidyverse)
library(grid)
library(gridExtra)
library(patchwork)
# --- Load input data ---

ps <- readRDS(input_path)
print(ps)

categories <- c("ice", "water", "sed")
region_categories <- c("Alps", "Kyrgyzstan")
region_venn_colors <- c("#4F8FC0", "#D98C5F")

# --- Environment-specific region-level ASV overlap ---

# Compare the two regions within each environment so that geographic overlap is
# not confounded by turnover among ice, streamwater, and sediment communities.
region_display_labels <- c("Alps", "Central Asia")

run_environment_region_overlap <- function(physeq_obj, source_categories,
                                           region_categories, output_dir,
                                           venn_colors,
                                           display_region_labels) {
  sample_metadata <- data.frame(
    sample_data(physeq_obj),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  regional_summary_list <- list()
  per_sample_list <- list()
  shared_asv_list <- list()
  venn_grob_list <- list()

  for (source_category in source_categories) {
    source_sample_ids <- rownames(sample_metadata)[
      as.character(sample_metadata$Source) == source_category
    ]

    ps_source <- prune_samples(source_sample_ids, physeq_obj)
    source_metadata <- data.frame(
      sample_data(ps_source),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )

    region_list <- setNames(vector("list", length(region_categories)),
                            region_categories)

    for (region in region_categories) {
      region_sample_ids <- rownames(source_metadata)[
        as.character(source_metadata$Region) == region
      ]

      ps_source_region <- prune_samples(region_sample_ids, ps_source)
      region_taxa_sums <- taxa_sums(ps_source_region)
      region_list[[region]] <- names(region_taxa_sums)[region_taxa_sums > 0]
    }

    shared_asvs <- Reduce(intersect, region_list)

    # Save a Venn diagram for the regional comparison within this environment.
    source_region_venn <- venn.diagram(
      x = region_list,
      category.names = display_region_labels,
      filename = NULL,
      output = TRUE,
      fill = venn_colors,
      alpha = 0.5,
      cex = 1.5,
      cat.cex = 1.2,
      cat.pos = c(-20, 20),
      main = paste("Regional ASV Overlap in", source_category)
    )

    venn_grob_list[[source_category]] <- grobTree(
      children = source_region_venn
    )

    ggsave(
      file.path(
        output_dir,
        paste0("venn_diagram_regions_within_", source_category, ".pdf")
      ),
      source_region_venn,
      width = 5,
      height = 5,
      dpi = 300
    )

    # Calculate the fraction of each environment-specific sample represented by
    # ASVs detected in both regions for that same environment.
    ps_source_rel <- transform_sample_counts(ps_source, function(x) x / sum(x))
    otu_rel_source_matrix <- as(otu_table(ps_source_rel), "matrix")
    if (taxa_are_rows(ps_source_rel)) {
      otu_rel_source_matrix <- t(otu_rel_source_matrix)
    }

    shared_rel_abundance <- rowSums(
      otu_rel_source_matrix[, shared_asvs, drop = FALSE]
    )

    metadata_order <- match(
      rownames(otu_rel_source_matrix),
      rownames(source_metadata)
    )

    per_sample_source <- tibble(
      Sample_ID = rownames(otu_rel_source_matrix),
      Source = as.character(source_metadata$Source[metadata_order]),
      Region = as.character(source_metadata$Region[metadata_order]),
      SharedRegionASVCount = length(shared_asvs),
      SharedRegionASVRelativeAbundance = shared_rel_abundance,
      SharedRegionASVRelativeAbundancePercent = 100 * shared_rel_abundance
    )

    source_summary <- per_sample_source %>%
      group_by(Source, Region) %>%
      summarise(
        SampleCount = n(),
        MeanSharedRegionASVRelativeAbundance =
          mean(SharedRegionASVRelativeAbundance),
        SDSharedRegionASVRelativeAbundance =
          sd(SharedRegionASVRelativeAbundance),
        MedianSharedRegionASVRelativeAbundance =
          median(SharedRegionASVRelativeAbundance),
        MinSharedRegionASVRelativeAbundance =
          min(SharedRegionASVRelativeAbundance),
        MaxSharedRegionASVRelativeAbundance =
          max(SharedRegionASVRelativeAbundance),
        .groups = "drop"
      ) %>%
      mutate(
        RegionalASVPoolCount = vapply(
          Region,
          function(region) length(region_list[[region]]),
          integer(1)
        ),
        SharedRegionASVCount = length(shared_asvs),
        SharedASVPercentOfRegionalPool =
          100 * SharedRegionASVCount / RegionalASVPoolCount,
        OverallMeanSharedRegionASVRelativeAbundance =
          mean(per_sample_source$SharedRegionASVRelativeAbundance),
        OverallMedianSharedRegionASVRelativeAbundance =
          median(per_sample_source$SharedRegionASVRelativeAbundance),
        .after = SampleCount
      ) %>%
      mutate(
        MeanSharedRegionASVRelativeAbundancePercent =
          100 * MeanSharedRegionASVRelativeAbundance,
        OverallMeanSharedRegionASVRelativeAbundancePercent =
          100 * OverallMeanSharedRegionASVRelativeAbundance
      )

    regional_summary_list[[source_category]] <- source_summary
    per_sample_list[[source_category]] <- per_sample_source
    shared_asv_list[[source_category]] <- tibble(
      Source = rep(source_category, length(shared_asvs)),
      ASV = shared_asvs
    )

    cat(
      source_category, ": ",
      length(region_list[[region_categories[1]]]), " ", region_categories[1],
      " ASVs; ",
      length(region_list[[region_categories[2]]]), " ", region_categories[2],
      " ASVs; ", length(shared_asvs), " shared ASVs\n",
      sep = ""
    )
  }

  combined_region_venn <- arrangeGrob(grobs = venn_grob_list, nrow = 1)
  ggsave(
    file.path(output_dir, "extended_data_figure_1_region_overlap.pdf"),
    combined_region_venn,
    width = 15,
    height = 5,
    dpi = 300
  )

  regional_summary <- bind_rows(regional_summary_list)
  per_sample_summary <- bind_rows(per_sample_list)
  shared_asv_table <- bind_rows(shared_asv_list)

  print(regional_summary)

  write.csv(
    regional_summary,
    file.path(output_dir, "environment_region_overlap_summary.csv"),
    row.names = FALSE
  )
  write.csv(
    per_sample_summary,
    file.path(
      output_dir,
      "environment_region_shared_asv_relative_abundance_per_sample.csv"
    ),
    row.names = FALSE
  )
  write.csv(
    shared_asv_table,
    file.path(output_dir, "environment_region_shared_asv_ids.csv"),
    row.names = FALSE
  )

  invisible(
    list(
      summary = regional_summary,
      per_sample = per_sample_summary,
      shared_asvs = shared_asv_table
    )
  )
}

environment_region_overlap_results <- run_environment_region_overlap(
  physeq_obj = ps,
  source_categories = categories,
  region_categories = region_categories,
  output_dir = result_dir,
  venn_colors = region_venn_colors,
  display_region_labels = region_display_labels
)

# --- Habitat-level ASV overlap and taxonomic composition ---

# Set the taxonomic rank used for composition summaries.
target_rank <- "Phylum"
# Set the number of dominant phyla to display.
N_top_phyla <- 13

# Define the palette for phylum-level plots.
color_palette <- c(
  Cyanobacteria = "#3385BB",
  Chloroflexi = "#DE9E83",
  Proteobacteria = "#A6CEE3",
  Actinobacteriota = "#F3E587",
  Acidobacteriota = "#7F9D55",
  Bacteroidota = "#F57C7C",
  Planctomycetota = "#E42622",
  Bdellovibrionota = "#84BF96",
  Myxococcota = "#FE8D19",
  Deinococcota = "#6DBD57",
  Patescibacteria = "#9D7BBA",
  Gemmatimonadota = "#977899",
  Verrucomicrobiota = "#FBB268",
  Other = "grey70"
)

# Define fill colors for habitat-level Venn diagrams.
venn_colors <- c("#6BADFD", "#8AD8B0", "#EDA871")

# Build a presence/absence ASV list for each habitat.
build_category_list <- function(physeq_obj, categories) {
  category_list <- list()

  for (cat in categories) {
    ps_subset <- prune_samples(sample_data(physeq_obj)$Source == cat, physeq_obj)
    present_asvs <- names(taxa_sums(ps_subset))[taxa_sums(ps_subset) > 0]
    category_list[[cat]] <- present_asvs
  }

  category_list
}

# Identify dominant taxa at the selected rank.
get_top_taxa_info <- function(physeq_obj, target_rank, n_top) {
  tax_df_local <- as.data.frame(as(tax_table(physeq_obj), "matrix"))
  tax_df_local$OTUID <- rownames(tax_df_local)

  otu_abundance <- data.frame(
    OTUID = names(taxa_sums(physeq_obj)),
    TotalAbundance = taxa_sums(physeq_obj)
  )

  merged_df <- left_join(tax_df_local, otu_abundance, by = "OTUID")

  top_taxa_list <- merged_df %>%
    group_by(!!sym(target_rank)) %>%
    summarise(TotalReads = sum(TotalAbundance, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(TotalReads)) %>%
    slice(1:n_top) %>%
    pull(!!sym(target_rank)) %>%
    stats::na.omit() %>%
    as.character()

  list(
    top_taxa_list = top_taxa_list,
    display_levels = c(top_taxa_list, "Other")
  )
}

# Plot the taxonomic composition of a selected ASV subset.
plot_taxonomic_composition <- function(
    asv_list,
    physeq_obj,
    top_taxa_list,
    display_levels,
    plot_title = "Taxonomic Composition") {
  if (length(asv_list) == 0) {
    return(
      ggplot() +
        annotate("text", x = 1, y = 1, label = "No ASVs in this subset", size = 3) +
        labs(title = plot_title, x = NULL, y = "Proportion Within Subset") +
        theme_void() +
        theme(plot.title = element_text(hjust = 0.5, size = 9))
    )
  }

  # Restrict the phyloseq object to the selected ASVs.
  ps_sub <- prune_taxa(asv_list, physeq_obj)

  # Aggregate taxa at the selected rank while retaining missing labels.
  ps_glom <- tax_glom(ps_sub, taxrank = target_rank, NArm = FALSE)

  # Extract rank-level taxonomy labels.
  raw_tax_vec <- as.data.frame(tax_table(ps_glom), stringsAsFactors = FALSE)[[target_rank]]

  # Collapse non-dominant and missing taxa into "Other".
  final_group_name <- dplyr::if_else(
    !is.na(raw_tax_vec) & raw_tax_vec %in% top_taxa_list,
    raw_tax_vec,
    "Other"
  )

  # Aggregate abundance by the final taxonomic group.
  df_plot <- data.frame(
    Taxon = final_group_name,
    Abundance = taxa_sums(ps_glom)
  ) %>%
    group_by(Taxon) %>%
    summarise(Abundance = sum(Abundance), .groups = "drop")

  # Convert abundance to within-subset relative proportions.
  df_plot$RelativeAbundance <- df_plot$Abundance / sum(df_plot$Abundance)

  # Set the plotting order.
  df_plot$Taxon <- factor(df_plot$Taxon, levels = display_levels)

  # Draw the stacked bar chart.
  p <- ggplot(df_plot, aes(x = "Overall Composition", y = RelativeAbundance, fill = Taxon)) +
    geom_bar(stat = "identity", position = "stack", width = 0.9) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(title = plot_title, y = "Proportion Within Subset", x = NULL, fill = target_rank) +
    theme_minimal(base_size = 8) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.major.x = element_blank(),
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 9)
    ) +
    scale_fill_manual(values = color_palette, limits = display_levels, drop = FALSE)

  return(p)
}

# Generate habitat-level Venn diagrams and composition panels.
run_habitat_analysis <- function(
    physeq_obj,
    analysis_label,
    venn_file,
    combined_file,
    shared_abundance_file) {
  category_list <- build_category_list(physeq_obj, categories)

  venn_plot <- venn.diagram(
    x = category_list,
    category.names = categories,
    filename = NULL,
    disable.logging = TRUE,
    output = TRUE,
    fill = venn_colors,
    alpha = 0.5,
    cex = 1.5,
    cat.cex = 1.5,
    main = paste("ASV Numbers in Sample Comparisons", analysis_label)
  )

  ggsave(venn_file, venn_plot, width = 5, height = 5, dpi = 300)

  for (cat in categories) {
    cat(analysis_label, "-", cat, ": ", length(category_list[[cat]]), " ASVs\n", sep = "")
  }

  only_ice_asvs <- setdiff(
    category_list[["ice"]],
    union(category_list[["water"]], category_list[["sed"]])
  )
  only_water_asvs <- setdiff(
    category_list[["water"]],
    union(category_list[["ice"]], category_list[["sed"]])
  )
  only_sed_asvs <- setdiff(
    category_list[["sed"]],
    union(category_list[["ice"]], category_list[["water"]])
  )
  shared_all_three_asvs <- Reduce(intersect, category_list)

  physeq_rel <- transform_sample_counts(physeq_obj, function(x) x / sum(x))

  otu_rel_habitat_matrix <- as(otu_table(physeq_rel), "matrix")
  if (taxa_are_rows(physeq_rel)) {
    otu_rel_habitat_matrix <- t(otu_rel_habitat_matrix)
  }

  shared_all_three_rel_abundance <- if (length(shared_all_three_asvs) > 0) {
    rowSums(otu_rel_habitat_matrix[, shared_all_three_asvs, drop = FALSE])
  } else {
    rep(0, nrow(otu_rel_habitat_matrix))
  }

  sample_info_local <- data.frame(sample_data(physeq_obj)) %>%
    rownames_to_column(var = "Sample_ID_from_rownames")

  if (!"Sample_ID" %in% colnames(sample_info_local)) {
    sample_info_local$Sample_ID <- sample_info_local$Sample_ID_from_rownames
  }

  habitat_shared_abundance_summary <- tibble(
    Sample_ID = rownames(otu_rel_habitat_matrix),
    SharedAllThreeHabitatASVRelAbundance = shared_all_three_rel_abundance
  ) %>%
    left_join(
      sample_info_local %>%
        select(Sample_ID, Source),
      by = "Sample_ID"
    ) %>%
    group_by(Source) %>%
    summarise(
      MeanSharedAllThreeHabitatASVRelAbundance =
        mean(SharedAllThreeHabitatASVRelAbundance, na.rm = TRUE),
      .groups = "drop"
    )

  overall_shared_all_three_abundance <-
    mean(shared_all_three_rel_abundance, na.rm = TRUE)

  habitat_shared_abundance_summary <- bind_rows(
    tibble(
      Source = "All samples",
      MeanSharedAllThreeHabitatASVRelAbundance =
        overall_shared_all_three_abundance
    ),
    habitat_shared_abundance_summary
  ) %>%
    mutate(
      Analysis = analysis_label,
      SharedAllThreeHabitatASVCount = length(shared_all_three_asvs),
      .before = 1
    )

  print(habitat_shared_abundance_summary)

  write.csv(
    habitat_shared_abundance_summary,
    shared_abundance_file,
    row.names = FALSE
  )

  ice_water_shared <- setdiff(
    intersect(category_list[["ice"]], category_list[["water"]]),
    category_list[["sed"]]
  )

  ice_sed_shared <- setdiff(
    intersect(category_list[["ice"]], category_list[["sed"]]),
    category_list[["water"]]
  )

  water_sed_shared <- setdiff(
    intersect(category_list[["water"]], category_list[["sed"]]),
    category_list[["ice"]]
  )

  top_taxa_info <- get_top_taxa_info(physeq_rel, target_rank, N_top_phyla)

  plot_only_ice <- plot_taxonomic_composition(
    asv_list = only_ice_asvs,
    physeq_obj = physeq_rel,
    top_taxa_list = top_taxa_info$top_taxa_list,
    display_levels = top_taxa_info$display_levels,
    plot_title = "Ice Only"
  ) +
    xlab(paste("", length(only_ice_asvs), " ASVs", sep = ""))

  plot_only_water <- plot_taxonomic_composition(
    asv_list = only_water_asvs,
    physeq_obj = physeq_rel,
    top_taxa_list = top_taxa_info$top_taxa_list,
    display_levels = top_taxa_info$display_levels,
    plot_title = "Water Only"
  ) +
    xlab(paste("", length(only_water_asvs), " ASVs", sep = ""))

  plot_only_sed <- plot_taxonomic_composition(
    asv_list = only_sed_asvs,
    physeq_obj = physeq_rel,
    top_taxa_list = top_taxa_info$top_taxa_list,
    display_levels = top_taxa_info$display_levels,
    plot_title = "Sediment Only"
  ) +
    xlab(paste("", length(only_sed_asvs), " ASVs", sep = ""))

  plot_shared_all_three <- plot_taxonomic_composition(
    asv_list = shared_all_three_asvs,
    physeq_obj = physeq_rel,
    top_taxa_list = top_taxa_info$top_taxa_list,
    display_levels = top_taxa_info$display_levels,
    plot_title = "Shared: All Three"
  ) +
    xlab(paste("", length(shared_all_three_asvs), " ASVs", sep = ""))

  plot_ice_water_shared <- plot_taxonomic_composition(
    asv_list = ice_water_shared,
    physeq_obj = physeq_rel,
    top_taxa_list = top_taxa_info$top_taxa_list,
    display_levels = top_taxa_info$display_levels,
    plot_title = "ice&water Only"
  ) +
    xlab(paste("", length(ice_water_shared), " ASVs", sep = ""))

  plot_ice_sed_shared <- plot_taxonomic_composition(
    asv_list = ice_sed_shared,
    physeq_obj = physeq_rel,
    top_taxa_list = top_taxa_info$top_taxa_list,
    display_levels = top_taxa_info$display_levels,
    plot_title = "ice&sed Only"
  ) +
    xlab(paste("", length(ice_sed_shared), " ASVs", sep = ""))

  plot_water_sed_shared <- plot_taxonomic_composition(
    asv_list = water_sed_shared,
    physeq_obj = physeq_rel,
    top_taxa_list = top_taxa_info$top_taxa_list,
    display_levels = top_taxa_info$display_levels,
    plot_title = "Water&Sed Only"
  ) +
    xlab(paste("", length(water_sed_shared), " ASVs", sep = ""))

  combined_plot <-
    plot_only_ice + theme(legend.position = "none") +
    plot_only_water + theme(legend.position = "none") +
    plot_only_sed + theme(legend.position = "none") +
    plot_ice_water_shared + theme(legend.position = "none") +
    plot_ice_sed_shared + theme(legend.position = "none") +
    plot_water_sed_shared + theme(legend.position = "none") +
    plot_shared_all_three + theme(legend.position = "right") +
    plot_layout(nrow = 1, guides = "collect") &
    theme(
      plot.title = element_text(size = 9, hjust = 0.5),
      axis.title.y = element_text(size = 8)
    )

  combined_plot <- combined_plot + plot_annotation(title = analysis_label)

  ggsave(combined_file, plot = combined_plot, width = 10, height = 5)
}

# Run the habitat analysis across all samples.
run_habitat_analysis(
  physeq_obj = ps,
  analysis_label = "All Regions",
  venn_file = file.path(result_dir, "figure_1a_habitat_overlap.pdf"),
  combined_file = file.path(result_dir, "extended_data_figure_2_intersection_taxonomy.pdf"),
  shared_abundance_file = file.path(
    result_dir,
    "habitat_shared_all_three_asv_relative_abundance_summary.csv"
  )
)
