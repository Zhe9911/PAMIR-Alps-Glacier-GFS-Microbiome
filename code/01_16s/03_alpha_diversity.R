# 16S rRNA alpha-diversity analysis
# Manuscript outputs: Figures 1b and 2b.

rm(list = ls())
gc()
graphics.off()

input_path <- file.path(
  "results", "01_16s", "intermediate", "PAMIR_16S_final.rds"
)
result_dir <- file.path("results", "01_16s", "alpha_diversity")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

# --- Setup ---

library(phyloseq)
library(iNEXT)
library(ggplot2)
library(dplyr)
library(tibble)
library(FSA)
library(broom)
library(tidyr)
library(ggpubr)

# --- Load input data ---

PAMIR_16S <- readRDS(input_path)
PAMIR_16S

# Extract the abundance matrix and orient it as taxa x samples for iNEXT.
otu_mat <- as(otu_table(PAMIR_16S), "matrix")
# iNEXT expects taxa in rows and samples in columns.
if (!taxa_are_rows(PAMIR_16S)) {
  otu_mat <- t(otu_mat)
}

# Preserve sample identifiers as an explicit column for later metadata joins.
meta_df <- as(sample_data(PAMIR_16S), "data.frame") %>%
  rownames_to_column("SampleID")

# --- Standardized Hill numbers by habitat ---

# Standardize samples to 95% coverage so diversity is compared at equal completeness.
# Set a fixed seed immediately before estimateD so its bootstrap intervals and
# downstream jitter are reproducible.
set.seed(666)
div_est_all <- estimateD(
  otu_mat,
  q = c(0, 1),
  datatype = "abundance",
  base = "coverage",
  level = 0.95
)

# Join metadata and define the plotting order.
stats_df_all <- div_est_all %>%
  left_join(meta_df, by = c("Assemblage" = "SampleID")) %>%
  mutate(
    Diversity_Type = case_when(
      Order.q == 0 ~ "Richness (q=0)",
      Order.q == 1 ~ "Shannon (q=1)"
    )
  ) %>%
  mutate(Source = factor(Source, levels = c("ice", "sed", "water"))) %>%
  mutate(
    Group_GI = factor(
      Group_GI,
      levels = c(
        "ice", "High_Glacial_water", "Mid_Glacial_water", "Low_Glacial_water",
        "High_Glacial_sed", "Mid_Glacial_sed", "Low_Glacial_sed"
      )
    )
  )

# Test whether each diversity metric differs globally among ice, sediment, and water.
stats_summary <- stats_df_all %>%
  group_by(Diversity_Type) %>%
  group_modify(~ tidy(kruskal.test(qD ~ Source, data = .x)))

print(stats_summary)
write.csv(
  stats_summary,
  file.path(result_dir, "figure_1b_kruskal_wallis.csv"),
  row.names = FALSE
)

# Follow significant global differences with Dunn pairwise tests and BH correction.
pairwise_results <- stats_df_all %>%
  group_by(Diversity_Type) %>%
  group_modify(~ {
    d <- dunnTest(qD ~ Source, data = .x, method = "bh")
    res <- as.data.frame(d$res)
    # Convert adjusted P values to significance labels.
    res <- res %>%
      mutate(p.signif = case_when(
        P.adj > 0.05 ~ "ns",
        P.adj <= 0.001 ~ "***",
        P.adj <= 0.01 ~ "**",
        P.adj <= 0.05 ~ "*"
      ))
    return(res)
  }) %>%
  ungroup()

# Retain only significant pairwise comparisons for annotation.
pairwise_clean <- pairwise_results %>%
  select(Diversity_Type, Comparison, P.adj, p.signif) %>%
  separate(Comparison, into = c("group1", "group2"), sep = " - ") %>%
  filter(P.adj < 0.05)

write.csv(
  pairwise_results,
  file.path(result_dir, "figure_1b_dunn_bh.csv"),
  row.names = FALSE
)

# Build facet-specific y positions for significance brackets.
max_vals <- stats_df_all %>%
  group_by(Diversity_Type) %>%
  summarise(max_y = max(qD, na.rm = TRUE))

sig_df <- pairwise_clean %>%
  left_join(max_vals, by = "Diversity_Type") %>%
  group_by(Diversity_Type) %>%
  mutate(
    # Stack brackets upward within each facet.
    y.position = max_y * (1.05 + 0.08 * (row_number() - 1))
  ) %>%
  ungroup()

# Plot standardized Hill numbers across habitats.
p_all_metrics <- ggplot(stats_df_all, aes(x = Source, y = qD, fill = Source)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.9) +
  geom_jitter(
    aes(shape = Region),
    width = 0.2,
    alpha = 0.4,
    size = 1.5,
    color = "black"
  ) +
  scale_fill_manual(values = c(
    ice = "#6BADFD",
    sed = "#EDA871",
    water = "#8AD8B0"
  )) +
  scale_shape_manual(values = c(Alps = 16, Kyrgyzstan = 17)) +
  facet_wrap(~ Diversity_Type, scales = "free_y") +
  stat_pvalue_manual(
    sig_df,
    label = "p.signif",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    tip.length = 0.01,
    vjust = 0.5
  ) +
  labs(
    y = "Standardized Diversity (Coverage = 95%)",
    x = NULL,
    title = NULL
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11, color = "black"),
    axis.text.y = element_text(color = "black"),
    strip.text = element_text(size = 12, face = "bold"),
    strip.background = element_rect(fill = "grey95")
  )

# Save the habitat comparison figure and its underlying standardized estimates.
ggsave(
  file.path(result_dir, "figure_1b_hill_diversity.pdf"),
  p_all_metrics,
  width = 11,
  height = 8,
  dpi = 300
)

write.csv(
  stats_df_all,
  file.path(result_dir, "coverage_95_hill_estimates.csv"),
  row.names = FALSE
)

# --- Ice and sediment group comparison ---

# Retain ice and sediment samples and separate richness and Shannon estimates.
stats_df <- subset(stats_df_all, Source != "water")
stats_df_1 <- subset(stats_df, Order.q == 1)

# Define the directional spatial sequence used as an ordered predictor.
trend_levels <- c("ice", "High_Glacial_sed", "Mid_Glacial_sed", "Low_Glacial_sed")

# Encode the four groups as ranks 1-4 for the Shannon-only trend test.
sed_trend_df <- subset(stats_df_1, Group_GI %in% trend_levels)
sed_trend_df$Group_GI <- factor(
  sed_trend_df$Group_GI,
  levels = trend_levels,
  ordered = TRUE
)
sed_trend_df$Group_order <- as.integer(sed_trend_df$Group_GI)

test_one_metric <- function(data) {
  # Use a two-sided Spearman test because the trend direction was not prespecified.
  # exact = FALSE accommodates tied group ranks in the ordered predictor.
  complete_data <- data[complete.cases(data[, c("Group_order", "qD")]), ]
  test <- cor.test(
    complete_data$Group_order,
    complete_data$qD,
    method = "spearman",
    alternative = "two.sided",
    exact = FALSE
  )

  data.frame(
    Order.q = unique(complete_data$Order.q),
    Diversity_Type = unique(complete_data$Diversity_Type),
    n = nrow(complete_data),
    rho = unname(test$estimate),
    p_value = test$p.value
  )
}

sed_trend_results <- do.call(
  rbind,
  lapply(split(sed_trend_df, sed_trend_df$Diversity_Type), test_one_metric)
)

sed_trend_results$p.signif <- ifelse(
  sed_trend_results$p_value <= 0.001, "***",
  ifelse(
    sed_trend_results$p_value <= 0.01, "**",
    ifelse(sed_trend_results$p_value <= 0.05, "*", "ns")
  )
)

sed_trend_results <- sed_trend_results[order(sed_trend_results$Order.q), ]

print(sed_trend_results)

# Export the Shannon effect size and raw two-sided P value for reporting.
write.csv(
  sed_trend_results,
  file.path(result_dir, "figure_2b_ice_to_sediment_spearman.csv"),
  row.names = FALSE
)

# Use concise display labels while preserving the predefined group order.
group_labels <- c("Ice", "High_GI_sed", "Mid_GI_sed", "Low_GI_sed")

shannon_df <- subset(stats_df_1, Group_GI %in% trend_levels)
shannon_df$Group_GI <- factor(
  shannon_df$Group_GI,
  levels = trend_levels,
  labels = group_labels
)
trend_shannon <- subset(sed_trend_results, Order.q == 1)

# Report the Shannon effect size and raw two-sided P value on the plot.
trend_label <- sprintf(
  "Spearman trend: rho = %.3f, P = %.2e",
  trend_shannon$rho,
  trend_shannon$p_value
)

y_max <- max(shannon_df$qD, na.rm = TRUE)
# Round the upper limit upward to leave stable space for the trend annotation.
y_upper <- ceiling(y_max * 1.15 / 100) * 100
y_label <- y_upper * 0.98

shannon_plot <- ggplot(shannon_df, aes(x = Group_GI, y = qD, fill = Group_GI)) +
  geom_boxplot(
    width = 0.62,
    outlier.shape = NA,
    color = "grey20",
    linewidth = 1.5,
    alpha = 0.85
  ) +
  geom_jitter(
    aes(shape = Region),
    width = 0.18,
    size = 3.0,
    alpha = 0.8,
    color = "grey25"
  ) +
  annotate(
    "text",
    x = 2.5,
    y = y_label,
    label = trend_label,
    size = 5.1,
    color = "black"
  ) +
  scale_fill_viridis_d(option = "D", begin = 0.15, end = 0.85) +
  scale_shape_manual(values = c(21, 24)) +
  scale_y_continuous(
    limits = c(0, y_upper),
    breaks = seq(0, y_upper, by = 250),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  labs(x = NULL, y = "Shannon diversity") +
  theme_classic(base_size = 22) +
  theme(
    legend.position = "none",
    axis.line = element_line(color = "black", linewidth = 1.4),
    axis.ticks = element_line(color = "black", linewidth = 1.3),
    axis.ticks.length = grid::unit(0.12, "cm"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title.y = element_text(color = "black", margin = margin(r = 12)),
    plot.margin = margin(t = 10, r = 12, b = 10, l = 10)
  )

ggsave(
  file.path(result_dir, "figure_2b_shannon_diversity.pdf"),
  shannon_plot,
  width = 7,
  height = 7
)
