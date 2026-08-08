# Polaromonas all-KO prevalence analysis across clades using Fisher's exact tests
# Manuscript outputs: Figures 5c-d and the Polaromonas portion of
# Supplementary Data 2.

rm(list = ls())
gc()
graphics.off()

library(tidyverse)
library(patchwork)

# --- 1. Define input and output paths ---

genus <- "Polaromonas"
input_dir <- file.path("data", "processed", "polaromonas")
result_dir <- file.path("results", "05_polaromonas")
out_dir <- file.path(result_dir, "ko_clade_presence")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ko_copy_number_file <- file.path(
  result_dir,
  paste0(genus, "_KO_copy_number_dt.tsv")
)
pangenome_summary_file <- file.path(
  input_dir,
  paste0(genus, "_Pan_gene_clusters_summary.txt")
)
gene_cluster_classification_file <- file.path(
  result_dir,
  paste0(genus, "_gene_cluster_prevalence_classification.tsv")
)
clade_file <- file.path(input_dir, "MAGs_Group_clade.tsv")
marker_gene_file <- file.path("config", "polaromonas-marker-genes.tsv")
mags_info_file <- file.path(
  "results",
  "03_mags",
  "intermediate",
  "MAGs_info.csv"
)
supplementary_data_2_file <- file.path(out_dir, "Supplementary Data 2.csv")
heatmap_out_dir <- file.path(result_dir, "marker_gene_clade_heatmap")
dir.create(heatmap_out_dir, recursive = TRUE, showWarnings = FALSE)

clade_df <- read_tsv(clade_file, show_col_types = FALSE) %>%
  transmute(
    MAGs = as.character(MAGs),
    Group_clade = factor(Group_clade, levels = c("clade1", "clade2"))
  )

ko_copy_number <- read_tsv(ko_copy_number_file, show_col_types = FALSE) %>%
  mutate(KO = as.character(KO))

# --- 2. Convert KO copy number to presence/absence and summarize prevalence ---

# A KO is considered present when its copy number is greater than zero. All MAGs
# assigned to Clade-D or Clade-G are retained in every KO comparison.

ko_presence_long <- ko_copy_number %>%
  select(KO, all_of(clade_df$MAGs)) %>%
  pivot_longer(
    cols = all_of(clade_df$MAGs),
    names_to = "MAGs",
    values_to = "copy_number"
  ) %>%
  mutate(
    copy_number = as.numeric(copy_number),
    presence = copy_number > 0
  ) %>%
  left_join(clade_df, by = "MAGs")

ko_clade_summary <- ko_presence_long %>%
  group_by(KO, Group_clade) %>%
  summarise(
    n_mags = n(),
    n_present = sum(presence),
    prevalence = n_present / n_mags,
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Group_clade,
    values_from = c(n_mags, n_present, prevalence),
    names_glue = "{.value}_{Group_clade}"
  ) %>%
  mutate(
    prevalence_diff_clade2_minus_clade1 =
      prevalence_clade2 - prevalence_clade1,
    abs_prevalence_diff = abs(prevalence_diff_clade2_minus_clade1)
  )

# --- 3. Test KO prevalence differences between Clade-D and Clade-G ---

# Two-sided Fisher exact tests are applied to every KO, followed by BH-FDR
# correction across the complete tested KO set. Both raw and adjusted P values
# are retained for transparent reporting.

run_fisher <- function(
  n_present_clade1,
  n_mags_clade1,
  n_present_clade2,
  n_mags_clade2
) {
  contingency_table <- matrix(
    c(
      n_mags_clade1 - n_present_clade1,
      n_present_clade1,
      n_mags_clade2 - n_present_clade2,
      n_present_clade2
    ),
    nrow = 2,
    byrow = TRUE
  )

  fisher_result <- fisher.test(contingency_table, alternative = "two.sided")

  tibble(
    fisher_odds_ratio_clade2_vs_clade1 =
      unname(fisher_result$estimate),
    fisher_p_value = fisher_result$p.value
  )
}

fisher_results <- pmap_dfr(
  ko_clade_summary %>%
    select(
      n_present_clade1,
      n_mags_clade1,
      n_present_clade2,
      n_mags_clade2
    ),
  run_fisher
)

all_ko_results <- bind_cols(ko_clade_summary, fisher_results) %>%
  mutate(
    fisher_q_value = p.adjust(fisher_p_value, method = "BH"),
    prevalence_direction = case_when(
      prevalence_diff_clade2_minus_clade1 > 0 ~ "clade2_higher",
      prevalence_diff_clade2_minus_clade1 < 0 ~ "clade1_higher",
      TRUE ~ "equal"
    )
  ) %>%
  arrange(fisher_p_value, desc(abs_prevalence_diff), KO)

results_file <- file.path(
  out_dir,
  paste0(genus, "_all_KO_clade_presence_fisher_summary.csv")
)
summary_file <- file.path(out_dir, paste0(genus, "_analysis_summary.txt"))

write_csv(all_ko_results, results_file, na = "NA")

n_clade1 <- sum(clade_df$Group_clade == "clade1")
n_clade2 <- sum(clade_df$Group_clade == "clade2")

summary_lines <- c(
  paste0("Genus: ", genus),
  paste0("All clade MAGs retained: yes (", nrow(clade_df), ")"),
  paste0("Clade-D MAGs: ", n_clade1),
  paste0("Clade-G MAGs: ", n_clade2),
  paste0("KOs tested: ", nrow(all_ko_results)),
  paste0(
    "KOs with Fisher raw P < 0.05: ",
    sum(all_ko_results$fisher_p_value < 0.05)
  ),
  paste0(
    "KOs with Fisher BH q < 0.05: ",
    sum(all_ko_results$fisher_q_value < 0.05)
  ),
  paste0("Result table: ", results_file)
)

write_lines(summary_lines, summary_file)
cat(paste(summary_lines, collapse = "\n"), "\n")

# --- 4. Define candidate clade-associated functional modules ---

# Candidate clade-associated functional modules were identified by integrating
# pronounced inter-clade prevalence differences, co-occurrence of functionally
# related genes, and established biological knowledge.
# config/polaromonas-marker-genes.tsv records the curated KO set, functional
# categories, canonical focal or supporting roles, and the Figure 5c-d display
# subset and order.

# --- 5. Map focal marker KOs to pangenome gene clusters ---

# Only focal KOs contribute to this summary. Unique associated gene clusters are
# classified as Soft_Core, Shell, or Cloud using the categories defined by
# 01_pangenome_and_tree.R; Shell and Cloud together constitute the accessory
# fraction.

marker_genes <- read_tsv(
  marker_gene_file,
  show_col_types = FALSE,
  col_types = cols(
    KO = col_character(),
    gene_name = col_character(),
    marker_role = col_character(),
    functional_category = col_character(),
    include_in_heatmap = col_logical(),
    heatmap_order = col_integer()
  )
)

heatmap_orders <- marker_genes$heatmap_order[
  marker_genes$include_in_heatmap
]

stopifnot(
  nrow(marker_genes) == 98L,
  n_distinct(marker_genes$KO) == 98L,
  sum(marker_genes$marker_role == "focal") == 43L,
  sum(marker_genes$marker_role == "supporting") == 55L,
  !anyNA(marker_genes$include_in_heatmap),
  sum(marker_genes$include_in_heatmap) == 75L,
  identical(sort(heatmap_orders), seq_len(75L)),
  all(is.na(
    marker_genes$heatmap_order[!marker_genes$include_in_heatmap]
  ))
)

focal_marker_genes <- marker_genes %>%
  filter(marker_role == "focal") %>%
  distinct(KO, gene_name)

pangenome_summary <- read_tsv(
  pangenome_summary_file,
  quote = "",
  show_col_types = FALSE
)

gene_cluster_classification <- read_tsv(
  gene_cluster_classification_file,
  show_col_types = FALSE
)

focal_marker_gene_clusters <- pangenome_summary %>%
  select(gene_cluster_id, KOfam_ACC) %>%
  filter(!is.na(KOfam_ACC), KOfam_ACC != "") %>%
  distinct() %>%
  inner_join(focal_marker_genes, by = c("KOfam_ACC" = "KO")) %>%
  left_join(gene_cluster_classification, by = "gene_cluster_id") %>%
  transmute(
    KO = KOfam_ACC,
    gene_name,
    gene_cluster_id,
    pangenome_category
  ) %>%
  distinct() %>%
  arrange(KO, gene_name, pangenome_category, gene_cluster_id)

focal_gene_cluster_category_counts <- focal_marker_gene_clusters %>%
  count(pangenome_category, name = "n_gene_clusters") %>%
  mutate(
    total_focal_gene_clusters = sum(n_gene_clusters),
    percent = round(100 * n_gene_clusters / total_focal_gene_clusters, 1)
  )

focal_gene_cluster_accessory_summary <- focal_marker_gene_clusters %>%
  summarise(
    n_focal_kos = n_distinct(KO),
    n_focal_gene_clusters = n_distinct(gene_cluster_id),
    n_soft_core_gene_clusters = sum(pangenome_category == "Soft_Core"),
    n_shell_gene_clusters = sum(pangenome_category == "Shell"),
    n_cloud_gene_clusters = sum(pangenome_category == "Cloud"),
    n_accessory_gene_clusters = sum(pangenome_category %in% c("Shell", "Cloud")),
    accessory_percent = round(
      100 * n_accessory_gene_clusters / n_focal_gene_clusters,
      1
    )
  )

write_tsv(
  focal_marker_gene_clusters,
  file.path(
    heatmap_out_dir,
    paste0(genus, "_focal_KO_gene_cluster_prevalence_categories.tsv")
  )
)
write_tsv(
  focal_gene_cluster_category_counts,
  file.path(
    heatmap_out_dir,
    paste0(genus, "_focal_KO_gene_cluster_category_counts.tsv")
  )
)
write_tsv(
  focal_gene_cluster_accessory_summary,
  file.path(
    heatmap_out_dir,
    paste0(genus, "_focal_KO_gene_cluster_accessory_summary.tsv")
  )
)

# --- 6. Export the complete curated marker set as Supplementary Data 2 ---

# The table combines functional interpretation with clade-specific prevalence,
# prevalence contrast, Fisher statistics, and the clade with higher prevalence.
marker_ko_functions <- pangenome_summary %>%
  transmute(
    KO = KOfam_ACC,
    KO_function = str_squish(str_remove(KOfam, "\\s*\\[EC:[^]]+\\]$"))
  ) %>%
  filter(!is.na(KO), KO != "", !is.na(KO_function), KO_function != "") %>%
  distinct(KO, .keep_all = TRUE)

supplementary_data_2 <- marker_genes %>%
  left_join(marker_ko_functions, by = "KO") %>%
  left_join(all_ko_results, by = "KO") %>%
  transmute(
    Functional_category = functional_category,
    Marker_role = marker_role,
    KO,
    Gene_name = gene_name,
    KO_function,
    `Clade-D_prevalence` = round(prevalence_clade1, 3),
    `Clade-G_prevalence` = round(prevalence_clade2, 3),
    `Delta_prevalence_Clade-G_minus_Clade-D` = round(
      prevalence_diff_clade2_minus_clade1,
      3
    ),
    Fisher_raw_P = signif(fisher_p_value, 4),
    Fisher_BH_Q = signif(fisher_q_value, 4),
    Higher_prevalence_clade = recode(
      prevalence_direction,
      clade1_higher = "Clade-D",
      clade2_higher = "Clade-G",
      equal = "Equal"
    )
  )

write_csv(supplementary_data_2, supplementary_data_2_file, na = "NA")

# --- 7. Define the ecologically informed subset for the main-text heatmap ---

# Guided by ecological context, a subset of candidate functional modules was
# selected for visualization. This subset is a presentation choice rather than
# an additional statistical filter. Its membership and display order are stored
# with the complete marker set in config/polaromonas-marker-genes.tsv.

# --- 8. Select heatmap markers and order MAGs for display ---

# Clade-G is displayed above Clade-D. Within each clade, MAGs are ordered from
# highest to lowest LFC-GI, and the same order is used by all annotation tracks.

heatmap_marker_genes <- marker_genes %>%
  filter(include_in_heatmap) %>%
  arrange(heatmap_order) %>%
  transmute(
    marker_order = heatmap_order,
    KO,
    gene_name,
    duplicated_label = duplicated(gene_name) |
      duplicated(gene_name, fromLast = TRUE),
    row_label = if_else(
      duplicated_label,
      paste0(gene_name, " [", KO, "]"),
      gene_name
    ),
    row_label = make.unique(row_label, sep = " #")
  )

mag_metadata <- clade_df %>%
  left_join(
    read_csv(mags_info_file, show_col_types = FALSE) %>%
      transmute(
        MAGs = as.character(MAGs),
        Genus = as.character(Genus),
        lfc_GI = as.numeric(lfc_GI)
      ),
    by = "MAGs"
  ) %>%
  filter(Genus == genus | is.na(Genus)) %>%
  arrange(desc(Group_clade), desc(lfc_GI), MAGs) %>%
  mutate(MAGs = factor(MAGs, levels = MAGs))

# --- 9. Convert marker KOs to heatmap copy-number states ---

# Copy numbers are displayed as absent (0), single-copy (1), or multicopy (2+);
# exact copy-number labels are omitted from individual heatmap cells.

heatmap_df <- heatmap_marker_genes %>%
  select(marker_order, KO, gene_name, row_label) %>%
  left_join(ko_copy_number, by = "KO") %>%
  select(
    marker_order,
    KO,
    gene_name,
    row_label,
    all_of(as.character(mag_metadata$MAGs))
  ) %>%
  pivot_longer(
    cols = all_of(as.character(mag_metadata$MAGs)),
    names_to = "MAGs",
    values_to = "copy_number"
  ) %>%
  mutate(
    copy_number = replace_na(as.numeric(copy_number), 0),
    copy_state = case_when(
      copy_number == 0 ~ "0",
      copy_number == 1 ~ "1",
      copy_number >= 2 ~ "2+"
    ),
    copy_state = factor(copy_state, levels = c("0", "1", "2+")),
    MAGs = factor(MAGs, levels = rev(levels(mag_metadata$MAGs))),
    row_label = factor(row_label, levels = heatmap_marker_genes$row_label)
  )

clade_boundary <- sum(mag_metadata$Group_clade == "clade1") + 0.5

# --- 10. Plot the horizontal marker-gene copy-number heatmap ---

# Marker genes are columns and MAGs are rows. Clade membership and LFC-GI are
# shown as aligned side annotations, with copy-number state encoded by color.

clade_annotation_df <- mag_metadata %>%
  transmute(
    MAGs = factor(MAGs, levels = rev(levels(mag_metadata$MAGs))),
    annotation = "Clade",
    Clade = recode(
      as.character(Group_clade),
      clade1 = "Clade-D",
      clade2 = "Clade-G"
    )
  )

clade_label_df <- clade_annotation_df %>%
  group_by(Clade) %>%
  slice(ceiling(n() / 2)) %>%
  ungroup()

p_clade_annotation <- ggplot(
  clade_annotation_df,
  aes(x = annotation, y = MAGs, fill = Clade)
) +
  geom_tile(color = "white", linewidth = 0.18) +
  geom_text(
    data = clade_label_df,
    aes(label = Clade),
    fontface = "bold",
    size = 2.8,
    angle = 90
  ) +
  geom_hline(yintercept = clade_boundary, linewidth = 0.4, color = "grey20") +
  scale_fill_manual(
    values = c("Clade-D" = "#9ECAE1", "Clade-G" = "#FCAE91"),
    guide = "none"
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 9) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none",
    plot.margin = margin(10, 0, 10, 10)
  )

lfc_annotation_df <- mag_metadata %>%
  transmute(
    MAGs = factor(MAGs, levels = rev(levels(mag_metadata$MAGs))),
    annotation = "lfc_GI",
    lfc_GI
  )

p_lfc_annotation <- ggplot(
  lfc_annotation_df,
  aes(x = annotation, y = MAGs, fill = lfc_GI)
) +
  geom_tile(color = "grey92", linewidth = 0.18) +
  geom_hline(yintercept = clade_boundary, linewidth = 0.4, color = "grey20") +
  scale_fill_gradient2(
    low = "#3288BD",
    mid = "white",
    high = "#D53E4F",
    midpoint = 0,
    name = "lfc_GI"
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 9) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right",
    plot.margin = margin(10, 0, 10, 0)
  )

p_marker_heatmap <- ggplot(
  heatmap_df,
  aes(x = row_label, y = MAGs, fill = copy_state)
) +
  geom_tile(color = "grey92", linewidth = 0.18) +
  geom_hline(yintercept = clade_boundary, linewidth = 0.4, color = "grey20") +
  scale_x_discrete(
    position = "top",
    guide = guide_axis(angle = 60)
  ) +
  scale_y_discrete(position = "right") +
  scale_fill_manual(
    values = c("0" = "grey96", "1" = "#f6c95f", "2+" = "#f7941d"),
    name = "KO copy number"
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 9) +
  theme(
    panel.grid = element_blank(),
    axis.text.x.top = element_text(
      size = 6,
      margin = margin(b = 0)
    ),
    axis.text.y = element_text(size = 6),
    legend.position = "bottom",
    plot.margin = margin(0, 10, 10, 0)
  )

p_marker_heatmap_with_lfc <- wrap_plots(
  p_clade_annotation,
  p_lfc_annotation,
  p_marker_heatmap,
  nrow = 1,
  widths = c(0.06, 0.08, 1),
  guides = "collect"
) &
  theme(legend.position = "right")

ggsave(
  file.path(
    heatmap_out_dir,
    "Polaromonas_marker_gene_clade_copy_number_heatmap.pdf"
  ),
  p_marker_heatmap_with_lfc,
  width = max(14, 0.20 * nrow(heatmap_marker_genes) + 3),
  height = max(5.5, 0.25 * nrow(mag_metadata) + 2)
)
