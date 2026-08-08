# KO-level association and KEGG GSEA using continuous LFC-GI
# Manuscript outputs: Figure 4, Extended Data Figure 6, and Supplementary Data 2.
# Writes complete pathway/module results and leading-edge KO tables.

rm(list = ls())
gc()
graphics.off()

out_dir <- file.path("results", "04_functions", "ko_gsea")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- Setup ---

library(tidyverse)
library(brglm2)
library(clusterProfiler)
library(ggrepel)

# --- Load input data ---

ko_mag_dt <- read_tsv(
  file.path("data", "processed", "metagenomics", "ko_mag_dt.tsv"),
  show_col_types = FALSE
)
MAGs_info <- read.csv(
  file.path("results", "03_mags", "intermediate", "MAGs_info.csv"),
  stringsAsFactors = FALSE
)
kegg_brite <- read.csv(
  file.path("data", "processed", "metagenomics", "kegg_brite.csv"),
  header = TRUE,
  sep = ",",
  stringsAsFactors = FALSE
)

# --- 1. Metadata for KO-by-MAG modeling ---

meta_df <- MAGs_info %>%
  transmute(
    MAGs = as.character(MAGs),
    lfc_GI = as.numeric(lfc_GI),
    MAG_completeness = as.numeric(Completeness) / 100,
    genome_size = as.numeric(Genome_Size),
    log10_genome_size = log10(as.numeric(Genome_Size))
  ) %>%
  filter(
    !is.na(MAGs),
    !is.na(lfc_GI),
    !is.na(MAG_completeness),
    !is.na(log10_genome_size)
  ) %>%
  distinct(MAGs, .keep_all = TRUE)

valid_mags <- intersect(meta_df$MAGs, colnames(ko_mag_dt))
meta_df <- meta_df %>%
  filter(MAGs %in% valid_mags) %>%
  arrange(match(MAGs, valid_mags))

cat("Valid MAG count:", nrow(meta_df), "\n")

# --- 2. KO presence/absence matrix ---

ko_data_subset <- ko_mag_dt[, c("KO", valid_mags)]

total_counts <- rowSums(ko_data_subset[, valid_mags])

# Keep KOs present in at least 5 valid MAGs. Very rare KOs provide too few
# presences for stable logistic-regression estimates and can make enrichment
# results sensitive to individual MAGs. KOs present in all MAGs are also
# uninformative because there is no presence/absence variation to model.

keep_ko_idx <- which(total_counts >= 5 & total_counts < length(valid_mags))

target_kos <- ko_data_subset$KO[keep_ko_idx]
ko_mat <- t(as.matrix(ko_data_subset[keep_ko_idx, valid_mags]))
storage.mode(ko_mat) <- "integer"
colnames(ko_mat) <- target_kos

cat("KO count tested:", length(target_kos), "\n")

# --- 3. KO-level Firth logistic regression ---

run_glm_robust <- function(ko_vec, metadata) {
  curr_dat <- metadata
  curr_dat$KO_presence <- as.integer(ko_vec)

  tryCatch({
    fit <- glm(
      KO_presence ~ lfc_GI + MAG_completeness + log10_genome_size,
      data = curr_dat,
      family = binomial(link = "logit"),
      method = "brglmFit"
    )

    coef_tab <- summary(fit)$coefficients
    target_term <- "lfc_GI"

    if (target_term %in% rownames(coef_tab)) {
      return(c(
        estimate = coef_tab[target_term, "Estimate"],
        std_error = coef_tab[target_term, "Std. Error"],
        statistic = coef_tab[target_term, "z value"],
        p.value = coef_tab[target_term, "Pr(>|z|)"]
      ))
    }

    c(estimate = NA, std_error = NA, statistic = NA, p.value = NA)
  }, error = function(e) {
    c(estimate = NA, std_error = NA, statistic = NA, p.value = NA)
  })
}

cat("Running KO-level Firth logistic models...\n")

glm_res_list <- apply(ko_mat, 2, function(x) run_glm_robust(x, meta_df))
glm_res_df <- as.data.frame(t(glm_res_list))
glm_res_df$KO <- rownames(glm_res_df)

# --- 4. Tidy regression results ---

df_stats <- glm_res_df %>%
  filter(!is.na(p.value)) %>%
  mutate(padj = p.adjust(p.value, method = "BH")) %>%
  left_join(
    kegg_brite %>%
      select(KO, Description) %>%
      distinct(),
    by = "KO"
  ) %>%
  mutate(
    Direction = case_when(
      padj < 0.05 & estimate > 0 ~ "GI-positive",
      padj < 0.05 & estimate < 0 ~ "GI-negative",
      TRUE ~ "NS"
    )
  ) %>%
  arrange(padj, desc(abs(estimate)))

write.csv(
  df_stats,
  file.path(out_dir, "KO_lfcGI_firth_all.csv"),
  row.names = FALSE
)

write.csv(
  df_stats %>%
    filter(Direction != "NS"),
  file.path(out_dir, "KO_lfcGI_firth_sig.csv"),
  row.names = FALSE
)

sig_gi_positive_ko <- df_stats %>%
  filter(Direction == "GI-positive") %>%
  pull(KO)

sig_gi_negative_ko <- df_stats %>%
  filter(Direction == "GI-negative") %>%
  pull(KO)

cat("GI-positive KOs (FDR < 0.05):", length(sig_gi_positive_ko), "\n")
cat("GI-negative KOs (FDR < 0.05):", length(sig_gi_negative_ko), "\n")

# --- 5. Volcano plot ---

volcano_df <- df_stats %>%
  mutate(
    log10p = -log10(pmax(p.value, 1e-300))
  )

fdr_p_threshold <- df_stats %>%
  filter(padj < 0.05) %>%
  pull(p.value) %>%
  max(na.rm = TRUE)

label_data <- volcano_df %>%
  filter(Direction != "NS") %>%
  arrange(padj) %>%
  slice_head(n = 20) %>%
  separate(
    Description,
    into = c("gene", "description"),
    sep = "; ",
    remove = FALSE,
    extra = "merge",
    fill = "right"
  )

p_volcano <- ggplot(volcano_df, aes(x = estimate, y = log10p)) +
  geom_point(
    data = subset(volcano_df, Direction == "NS"),
    color = "grey80",
    size = 1,
    alpha = 0.5
  ) +
  geom_point(
    data = subset(volcano_df, Direction != "NS"),
    aes(color = Direction),
    size = 2,
    alpha = 0.8
  ) +
  scale_color_manual(values = c(
    "GI-positive" = "#D53E4F",
    "GI-negative" = "#3288BD",
    "NS" = "#777777"
  )) +
  geom_hline(
    yintercept = -log10(fdr_p_threshold),
    linetype = "dashed",
    color = "black",
    linewidth = 0.3
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.3) +
  geom_text_repel(
    data = label_data,
    aes(label = gene),
    size = 3,
    box.padding = 0.5,
    point.padding = 0.3,
    max.overlaps = 20,
    show.legend = FALSE,
    color = "black"
  ) +
  theme_classic() +
  labs(
    x = "Log-odds slope for lfc_GI",
    y = "-Log10 (raw p value)",
    title = "KO Associations Along Continuous lfc_GI",
    subtitle = "Dashed horizontal line marks the largest raw p value among KOs passing FDR < 0.05"
  ) +
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  )

ggsave(
  file.path(out_dir, "KO_lfcGI_firth_volcano.pdf"),
  p_volcano,
  width = 8,
  height = 8
)

# --- 6. Stable KO ranking for pre-ranked GSEA ---

ko_rank_df <- df_stats %>%
  transmute(
    KO = as.character(KO),
    estimate = as.numeric(estimate),
    statistic = as.numeric(statistic),
    p.value = as.numeric(p.value),
    padj = as.numeric(padj),
    Description = as.character(Description),
    Direction = as.character(Direction)
  ) %>%
  filter(
    !is.na(KO),
    !is.na(statistic),
    is.finite(statistic)
  ) %>%
  distinct(KO, .keep_all = TRUE) %>%
  arrange(desc(statistic), KO) %>%
  group_by(statistic) %>%
  mutate(
    # Add a deterministic, tiny offset within tied statistics so fgsea does not
    # break ties arbitrarily. The alphabetically earlier KO stays earlier.
    rank_statistic = statistic +
      (n() - row_number() + 1L) * 1e-12 * as.integer(n() > 1L)
  ) %>%
  ungroup()

ko_rank <- ko_rank_df$rank_statistic
names(ko_rank) <- ko_rank_df$KO
ko_rank <- sort(ko_rank, decreasing = TRUE)

write_tsv(
  tibble(KO = names(ko_rank), rank_statistic = unname(ko_rank)),
  file.path(out_dir, "KO_ranked_statistic.tsv")
)

# --- 7. KEGG pathway GSEA from KO ranking ---

pathway_term2gene <- kegg_brite %>%
  transmute(
    term = as.character(PATH),
    gene = as.character(KO),
    name = as.character(C)
  ) %>%
  filter(
    !is.na(term),
    term != "",
    !is.na(gene),
    gene != "",
    gene %in% names(ko_rank)
  ) %>%
  distinct(term, gene, name)

pathway_term2name <- pathway_term2gene %>%
  select(term, name) %>%
  distinct()

pathway_term2gene <- pathway_term2gene %>%
  select(term, gene) %>%
  distinct()

pathway_size_df <- pathway_term2gene %>%
  count(term, name = "n_KOs_tested") %>%
  left_join(pathway_term2name, by = "term") %>%
  arrange(desc(n_KOs_tested), term)

min_ko_per_pathway <- 5
max_ko_per_pathway <- 500

pathway_term2gene_filtered <- pathway_term2gene %>%
  semi_join(
    pathway_size_df %>%
      filter(
        n_KOs_tested >= min_ko_per_pathway,
        n_KOs_tested <= max_ko_per_pathway
      ),
    by = "term"
  )

pathway_term2name_filtered <- pathway_term2name %>%
  semi_join(
    pathway_term2gene_filtered %>%
      distinct(term),
    by = "term"
  )

set.seed(42)
pathway_gsea <- GSEA(
  geneList = ko_rank,
  TERM2GENE = pathway_term2gene_filtered,
  TERM2NAME = pathway_term2name_filtered,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  minGSSize = min_ko_per_pathway,
  maxGSSize = max_ko_per_pathway,
  eps = 0,
  by = "fgsea",
  seed = TRUE,
  verbose = FALSE
)

pathway_gsea_df <- as_tibble(pathway_gsea@result) %>%
  arrange(p.adjust, desc(abs(NES))) %>%
  mutate(
    Direction = case_when(
      NES > 0 ~ "GI-positive",
      NES < 0 ~ "GI-negative",
      TRUE ~ "Neutral"
    )
  )

pathway_gsea_sig_df <- pathway_gsea_df %>%
  filter(!is.na(p.adjust), p.adjust < 0.05)

write.csv(
  pathway_gsea_df,
  file.path(out_dir, "KO_ranked_pathway_GSEA_all.csv"),
  row.names = FALSE
)

write.csv(
  pathway_gsea_sig_df,
  file.path(out_dir, "KO_ranked_pathway_GSEA_sig.csv"),
  row.names = FALSE
)

pathway_leading_edge_df <- pathway_gsea_sig_df %>%
  select(ID, Description, NES, p.adjust, core_enrichment) %>%
  filter(!is.na(core_enrichment), core_enrichment != "") %>%
  separate_longer_delim(core_enrichment, delim = "/") %>%
  rename(KO = core_enrichment) %>%
  left_join(
    ko_rank_df,
    by = "KO",
    suffix = c("_pathway", "_KO")
  ) %>%
  arrange(p.adjust, desc(abs(NES)), desc(statistic))

write.csv(
  pathway_leading_edge_df,
  file.path(out_dir, "KO_ranked_pathway_GSEA_leading_edge_KOs.csv"),
  row.names = FALSE
)

# --- 8. KEGG module GSEA from KO ranking ---

term2gene_file <- file.path(
  "data", "processed", "metagenomics", "kegg_module_TERM2GENE.tsv"
)
term2name_file <- file.path(
  "data", "processed", "metagenomics", "kegg_module_TERM2NAME.tsv"
)

module_term2gene <- read_tsv(term2gene_file, show_col_types = FALSE) %>%
  transmute(
    term = str_remove(term, "^md:"),
    gene = str_remove(KO, "^ko:")
  ) %>%
  distinct() %>%
  filter(gene %in% names(ko_rank))

module_term2name <- read_tsv(term2name_file, show_col_types = FALSE) %>%
  transmute(
    term = str_remove(term, "^md:"),
    name = name
  ) %>%
  distinct()

module_size_df <- module_term2gene %>%
  count(term, name = "n_KOs_tested") %>%
  left_join(module_term2name, by = "term") %>%
  arrange(desc(n_KOs_tested), term)

min_ko_per_module <- 5
max_ko_per_module <- 500

module_term2gene_filtered <- module_term2gene %>%
  semi_join(
    module_size_df %>%
      filter(
        n_KOs_tested >= min_ko_per_module,
        n_KOs_tested <= max_ko_per_module
      ),
    by = "term"
  )

module_term2name_filtered <- module_term2name %>%
  semi_join(
    module_term2gene_filtered %>%
      distinct(term),
    by = "term"
  )

set.seed(42)
module_gsea <- GSEA(
  geneList = ko_rank,
  TERM2GENE = module_term2gene_filtered,
  TERM2NAME = module_term2name_filtered,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  minGSSize = min_ko_per_module,
  maxGSSize = max_ko_per_module,
  eps = 0,
  by = "fgsea",
  seed = TRUE,
  verbose = FALSE
)

module_gsea_df <- as_tibble(module_gsea@result) %>%
  arrange(p.adjust, desc(abs(NES))) %>%
  mutate(
    Direction = case_when(
      NES > 0 ~ "GI-positive",
      NES < 0 ~ "GI-negative",
      TRUE ~ "Neutral"
    )
  )

module_gsea_sig_df <- module_gsea_df %>%
  filter(!is.na(p.adjust), p.adjust < 0.05)

write.csv(
  module_gsea_df,
  file.path(out_dir, "KO_ranked_module_GSEA_all.csv"),
  row.names = FALSE
)

write.csv(
  module_gsea_sig_df,
  file.path(out_dir, "KO_ranked_module_GSEA_sig.csv"),
  row.names = FALSE
)

module_leading_edge_df <- module_gsea_sig_df %>%
  select(ID, Description, NES, p.adjust, core_enrichment) %>%
  filter(!is.na(core_enrichment), core_enrichment != "") %>%
  separate_longer_delim(core_enrichment, delim = "/") %>%
  rename(KO = core_enrichment) %>%
  left_join(
    ko_rank_df,
    by = "KO",
    suffix = c("_module", "_KO")
  ) %>%
  arrange(p.adjust, desc(abs(NES)), desc(statistic))

write.csv(
  module_leading_edge_df,
  file.path(out_dir, "KO_ranked_module_GSEA_leading_edge_KOs.csv"),
  row.names = FALSE
)

# --- 9. Main-text KEGG GSEA result table ---

clean_kegg_description <- function(description) {
  description %>%
    str_remove("^\\d+\\s+") %>%
    str_remove("\\s+\\[PATH:.*?\\]$")
}

assign_main_text_functional_group <- function(kegg_id) {
  case_when(
    kegg_id %in% c(
      "ko00195", "ko00196", "M00145", "M00163", "M00161", "M00162"
    ) ~ "Photosynthetic light reactions and electron transport",
    kegg_id %in% c(
      "ko00500", "ko00630", "ko00650", "ko00720", "ko00620", "M00620"
    ) ~ "Carbon storage, organic-acid metabolism, and anaplerotic carbon pathways",
    kegg_id %in% c(
      "ko00906", "M01045", "M00097"
    ) ~ "Carotenoid pigment biosynthesis and photoprotection",
    kegg_id %in% c(
      "ko00860", "M00924", "M00932"
    ) ~ "Cofactor, tetrapyrrole, and quinone biosynthesis",
    kegg_id %in% c(
      "ko05111", "ko03070", "ko02040", "ko03060"
    ) ~ "Surface colonization, motility, secretion, and protein export",
    kegg_id %in% c(
      "ko04112"
    ) ~ "Cell cycle and cellular regulation",
    kegg_id %in% c(
      "ko00260", "ko00280", "ko00670", "ko00270", "M00975"
    ) ~ "Amino-acid, betaine, and one-carbon metabolism",
    kegg_id %in% c(
      "M00144", "M00151", "M00595"
    ) ~ "Respiratory energy metabolism",
    kegg_id %in% c(
      "ko00791", "ko00627"
    ) ~ "Aromatic compound metabolism and xenobiotic-transformation enzymes",
    TRUE ~ "Other functional enrichment"
  )
}

MAIN_TEXT_KEGG_ORDER <- c(
  "ko00195", "ko00196", "M00145", "M00163", "M00161", "M00162",
  "ko00906", "M01045", "M00097",
  "ko00860", "M00924", "M00932",
  "ko00500",
  "ko00791", "ko00627",
  "ko05111", "ko03070", "ko02040", "ko03060",
  "ko04112",
  "M00975", "ko00260", "ko00280", "ko00670", "ko00270",
  "ko00630", "ko00650", "ko00720", "ko00620", "M00620",
  "M00144", "M00151", "M00595"
)

format_gsea_for_supplement <- function(gsea_df, result_type) {
  gsea_df %>%
    transmute(
      Result_type = result_type,
      KEGG_ID = as.character(ID),
      KEGG_name = clean_kegg_description(Description),
      Broad_functional_group = assign_main_text_functional_group(ID),
      Broad_functional_group_source = "Author-assigned",
      GI_enrichment = as.character(Direction),
      # Number of KOs from this KEGG set that were represented in the ranked
      # KO test universe; this is not the complete known KEGG set size.
      KOs_in_ranked_test_set = as.integer(setSize),
      Core_enrichment_KO_count = str_count(core_enrichment, "/") + 1L,
      Core_enrichment_KOs = as.character(core_enrichment),
      # Retain full numerical precision in the supplementary data table.
      NES = as.numeric(NES),
      P_value = as.numeric(pvalue),
      FDR = as.numeric(p.adjust)
    )
}

main_text_kegg_gsea_table <- bind_rows(
  format_gsea_for_supplement(pathway_gsea_sig_df, "Pathway"),
  format_gsea_for_supplement(module_gsea_sig_df, "Module")
) %>%
  mutate(
    GI_enrichment = case_when(
      GI_enrichment == "GI-positive" ~ "High-GI enriched",
      GI_enrichment == "GI-negative" ~ "Low-GI enriched",
      TRUE ~ GI_enrichment
    ),
    GI_enrichment = factor(
      GI_enrichment,
      levels = c("High-GI enriched", "Low-GI enriched")
    ),
    Result_type = factor(Result_type, levels = c("Pathway", "Module"))
  ) %>%
  mutate(Main_text_order = match(KEGG_ID, MAIN_TEXT_KEGG_ORDER)) %>%
  arrange(
    is.na(Main_text_order),
    Main_text_order,
    GI_enrichment,
    Broad_functional_group,
    Result_type,
    FDR,
    desc(abs(NES))
  ) %>%
  mutate(
    Result_type = as.character(Result_type),
    GI_enrichment = as.character(GI_enrichment)
  ) %>%
  select(
    -Main_text_order
  )

write.csv(
  main_text_kegg_gsea_table,
  file.path(out_dir, "Supplementary Data 2.csv"),
  row.names = FALSE
)

# --- 10. Console summary ---

cat("Ranked KOs used:", length(ko_rank), "\n")
cat("Pathways tested:", n_distinct(pathway_term2gene_filtered$term), "\n")
cat("Significant GI-positive pathways (FDR < 0.05):", sum(pathway_gsea_sig_df$NES > 0), "\n")
cat("Significant GI-negative pathways (FDR < 0.05):", sum(pathway_gsea_sig_df$NES < 0), "\n")
cat("Modules tested:", n_distinct(module_term2gene_filtered$term), "\n")
cat("Significant GI-positive modules (FDR < 0.05):", sum(module_gsea_sig_df$NES > 0), "\n")
cat("Significant GI-negative modules (FDR < 0.05):", sum(module_gsea_sig_df$NES < 0), "\n")
