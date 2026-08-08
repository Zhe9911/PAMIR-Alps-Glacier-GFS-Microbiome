# Polaromonas MAG pangenome analysis
# Manuscript outputs: Figure 5a, pangenome matrices, and prevalence categories.

rm(list = ls())
gc()
graphics.off()

# --- 1. Paths and parameters ---

genus <- "Polaromonas"
input_dir <- file.path("data", "processed", "polaromonas")
genus_dir <- file.path("results", "05_polaromonas")
dir.create(genus_dir, recursive = TRUE, showWarnings = FALSE)

soft_core_prevalence_threshold <- 0.75
cloud_prevalence_threshold <- 0.15

trend_colors <- c(
  "Glaciophiles" = "#D53E4F",
  "Glaciophobes" = "#3288BD",
  "Non-sig" = "grey80"
)

lfc_low_color <- "#3288BD"
lfc_mid_color <- "white"
lfc_high_color <- "#D53E4F"

tree_file <- file.path(input_dir, paste0(genus, "_tree.treefile"))
outgroup_tree_file <- file.path(
  input_dir,
  paste0(genus, "_tree_with_outgroup.treefile")
)
summary_file <- file.path(
  input_dir,
  paste0(genus, "_Pan_gene_clusters_summary.txt")
)
mag_info_file <- file.path(
  "results",
  "03_mags",
  "intermediate",
  "MAGs_info.csv"
)

# --- 2. Load packages ---

library(tidyverse)
library(ggtree)
library(ape)

# --- 3. Read pangenome, phylogeny, and MAG metadata inputs ---

MAGs_info <- read.csv(mag_info_file, stringsAsFactors = FALSE)
df <- read_tsv(summary_file, quote = "", show_col_types = FALSE)

tree_tips <- read.tree(tree_file)$tip.label

# --- 4. Build gene-cluster presence/absence and KO copy-number matrices ---

# Gene-cluster counts are converted to binary presence/absence for prevalence
# analyses, whereas KO counts are retained for downstream copy-number plots.

gc_counts <- df %>%
  transmute(
    genome_name = as.character(genome_name),
    gene_cluster_id = as.character(gene_cluster_id),
    KOfam_ACC = as.character(KOfam_ACC)
  )

gc_copy_counts <- gc_counts %>%
  select(genome_name, gene_cluster_id) %>%
  count(genome_name, gene_cluster_id, name = "copy_number") %>%
  arrange(genome_name, gene_cluster_id)

pa_matrix <- gc_copy_counts %>%
  pivot_wider(
    names_from = gene_cluster_id,
    values_from = copy_number,
    values_fill = list(copy_number = 0)
  ) %>%
  column_to_rownames("genome_name")

pa_binary <- pa_matrix
pa_binary[] <- lapply(pa_binary, function(x) as.integer(x > 0))

genus_GC_binary <- rownames_to_column(
  as.data.frame(t(pa_binary)),
  var = "GC"
)

write_tsv(
  genus_GC_binary,
  file.path(genus_dir, paste0(genus, "_GC_presence_absence_dt.tsv"))
)

ko_copy_counts <- gc_counts %>%
  filter(!is.na(KOfam_ACC), KOfam_ACC != "") %>%
  count(genome_name, KOfam_ACC, name = "copy_number") %>%
  arrange(genome_name, KOfam_ACC)

ko_matrix <- ko_copy_counts %>%
  pivot_wider(
    names_from = KOfam_ACC,
    values_from = copy_number,
    values_fill = list(copy_number = 0)
  ) %>%
  column_to_rownames("genome_name")

genus_KO_copy_number <- rownames_to_column(
  as.data.frame(t(ko_matrix)),
  var = "KO"
)

write_tsv(
  genus_KO_copy_number,
  file.path(genus_dir, paste0(genus, "_KO_copy_number_dt.tsv"))
)

# --- 5. Prepare ecological and genome-quality metadata for tree annotation ---

mag_metadata <- MAGs_info %>%
  transmute(
    MAGs = as.character(MAGs),
    lfc_GI = as.numeric(lfc_GI),
    Trend = as.character(Trend),
    MAG_completeness = as.numeric(Completeness) / 100,
    log10_genome_size = log10(as.numeric(Genome_Size))
  ) %>%
  filter(
    MAGs %in% tree_tips,
    MAGs %in% rownames(pa_matrix),
    !is.na(lfc_GI),
    !is.na(MAG_completeness),
    !is.na(log10_genome_size)
  ) %>%
  distinct(MAGs, .keep_all = TRUE)

tree_metadata <- mag_metadata %>%
  transmute(
    label = MAGs,
    lfc_GI = lfc_GI,
    Trend = Trend
  )

# --- 6. Plot the rectangular rooted Polaromonas phylogeny for Figure 5a ---

outgroup_phylo_tree_full <- read.tree(outgroup_tree_file)

# These Rhodoferax outgroups define the root and are removed for Figure 5a.
outgroup_tips_to_drop <- c(
  "PAMIR_meta11_S19_comebin_16850_sub",
  "PAMIR_meta11_S19_concoct_30"
)

missing_outgroup_tips <- setdiff(
  outgroup_tips_to_drop,
  outgroup_phylo_tree_full$tip.label
)
if (length(missing_outgroup_tips) > 0) {
  stop(
    "Missing expected rooting outgroup tip(s): ",
    paste(missing_outgroup_tips, collapse = ", ")
  )
}

outgroup_phylo_tree <- drop.tip(outgroup_phylo_tree_full, outgroup_tips_to_drop)

p_rectangular <- ggtree(
  outgroup_phylo_tree,
  layout = "rectangular",
  linewidth = 0.55
) %<+% tree_metadata

root_edge_length <- 0.03
root_node <- Ntip(outgroup_phylo_tree) + 1L

p_rectangular_tree <- p_rectangular +
  geom_rootedge(rootedge = root_edge_length, linewidth = 0.8) +
  geom_point2(
    aes(subset = node == root_node),
    shape = 23,
    size = 3,
    stroke = 0.7,
    color = "black",
    fill = "white"
  ) +
  geom_treescale(width = 0.1, fontsize = 3.5) +
  geom_tiplab(
    aes(label = label),
    size = 2.6,
    align = TRUE,
    linetype = 0,
    offset = 0.085
  ) +
  geom_tippoint(
    aes(fill = Trend),
    shape = 21,
    size = 2.5,
    color = "grey30",
    stroke = 0.2,
    alpha = 0.9,
    na.rm = TRUE
  ) +
  scale_fill_manual(
    values = trend_colors,
    name = "Ecological Trend",
    na.value = "black"
  ) +
  theme(legend.position = "right")

rectangular_heatmap_df <- mag_metadata %>%
  transmute(
    OTUID = MAGs,
    lfc_GI_rect = lfc_GI
  )

p_rectangular_heat <- p_rectangular_tree +
  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(
    data = rectangular_heatmap_df,
    geom = geom_tile,
    mapping = aes(
      y = OTUID,
      fill = lfc_GI_rect,
      x = "LFC"
    ),
    pwidth = 0.055,
    offset = 0.04,
    axis.params = list(
      axis = "x",
      text.size = 3,
      text.angle = 0,
      vjust = 0.5,
      line.color = NA
    ),
    grid.params = list()
  ) +
  scale_fill_gradient2(
    low = lfc_low_color,
    mid = lfc_mid_color,
    high = lfc_high_color,
    midpoint = 0,
    name = "LFC GI",
    guide = guide_colorbar(
      barwidth = 0.5,
      barheight = 5,
      title.position = "top",
      frame.colour = "black",
      frame.linewidth = 0.5
    )
  )

rectangular_tree_pdf <- file.path(
  genus_dir,
  paste0(genus, "_phylo_tree_with_outgroup_rectangular_heatmap.pdf")
)
ggsave(
  rectangular_tree_pdf,
  p_rectangular_heat,
  width = 14,
  height = 11,
  device = cairo_pdf
)

# --- 7. Summarize the number of analyzed MAGs and gene clusters ---

analysis_summary <- tibble(
  genus = genus,
  n_mags = nrow(mag_metadata),
  n_gene_clusters = ncol(pa_matrix),
  n_kos = ncol(ko_matrix),
  rooting_outgroups = paste(outgroup_tips_to_drop, collapse = ";"),
  figure_5a_tree_pdf = rectangular_tree_pdf
)

write.csv(
  analysis_summary,
  file.path(genus_dir, "pangenome_genus_analysis_summary.csv"),
  row.names = FALSE
)

print(analysis_summary)

# --- 8. Inspect gene-cluster prevalence and define pangenome categories ---

# The prevalence thresholds were selected as explicit operational cutoffs after
# inspection of the observed gene-cluster frequency distribution. The
# distribution contains a dominant low-prevalence peak at 2-6 MAGs, followed by
# an intermediate long tail, and a high-prevalence enrichment beginning near
# 35 MAGs. Gene clusters occurring in <=15% of MAGs are classified as Cloud,
# those occurring in >=75% as Soft_Core, and all intermediate clusters as Shell.
#
# For the 46 Polaromonas MAGs analyzed here, the integer boundaries are
# floor(46 * 0.15) = 6 MAGs for Cloud and ceiling(46 * 0.75) = 35 MAGs for
# Soft_Core. These are distribution-supported operational categories rather
# than statistically estimated breakpoints. Soft_Core is a relaxed, MAG-aware
# category and is not equivalent to a strict 100% core genome.

total_genomes <- length(unique(df$genome_name))
soft_core_min_genomes <- ceiling(total_genomes * soft_core_prevalence_threshold)
cloud_max_genomes <- floor(total_genomes * cloud_prevalence_threshold)

pangenome_classification_thresholds <- tibble(
  genus = genus,
  total_genomes = total_genomes,
  cloud_prevalence_threshold = cloud_prevalence_threshold,
  cloud_max_genomes = cloud_max_genomes,
  soft_core_prevalence_threshold = soft_core_prevalence_threshold,
  soft_core_min_genomes = soft_core_min_genomes
)

gene_cluster_prevalence <- df %>%
  select(gene_cluster_id, num_genomes_gene_cluster_has_hits) %>%
  distinct()

# Retain zero-count prevalence levels so the exported table spans the full
# range from one to all 46 MAGs.
gene_cluster_frequency_distribution <- gene_cluster_prevalence %>%
  count(num_genomes_gene_cluster_has_hits, name = "n_gene_clusters") %>%
  complete(
    num_genomes_gene_cluster_has_hits = seq_len(total_genomes),
    fill = list(n_gene_clusters = 0)
  ) %>%
  arrange(num_genomes_gene_cluster_has_hits)

gene_cluster_frequency_distribution_file <- file.path(
  genus_dir,
  paste0(genus, "_gene_cluster_prevalence_frequency.tsv")
)

classification_thresholds_file <- file.path(
  genus_dir,
  paste0(genus, "_pangenome_classification_thresholds.tsv")
)

write_tsv(
  gene_cluster_frequency_distribution,
  gene_cluster_frequency_distribution_file
)

write_tsv(
  pangenome_classification_thresholds,
  classification_thresholds_file
)

gene_cluster_classification <- gene_cluster_prevalence %>%
  mutate(
    pangenome_category = case_when(
      num_genomes_gene_cluster_has_hits >= soft_core_min_genomes ~ "Soft_Core",
      num_genomes_gene_cluster_has_hits <= cloud_max_genomes ~ "Cloud",
      TRUE ~ "Shell"
    )
  ) %>%
  select(gene_cluster_id, pangenome_category) %>%
  arrange(pangenome_category, gene_cluster_id)

gene_cluster_classification_file <- file.path(
  genus_dir,
  paste0(genus, "_gene_cluster_prevalence_classification.tsv")
)

write_tsv(gene_cluster_classification, gene_cluster_classification_file)

print(pangenome_classification_thresholds)
print(table(gene_cluster_classification$pangenome_category))
print(gene_cluster_frequency_distribution_file)
print(classification_thresholds_file)
print(gene_cluster_classification_file)
