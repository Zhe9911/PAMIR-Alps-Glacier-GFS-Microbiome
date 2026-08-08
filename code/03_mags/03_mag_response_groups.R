# MAG response-group abundance analysis
# Manuscript output: Figure 3c.

rm(list = ls())
gc()
graphics.off()

out_dir <- file.path("results", "03_mags", "response_groups")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- Setup ---

library(tidyverse)
library(phyloseq)
library(ggplot2)

# --- Load input data ---

PAMIR_MAGs_rela <- readRDS(
  file.path("results", "03_mags", "intermediate", "PAMIR_MAGs_rela.rds")
)
print(PAMIR_MAGs_rela)

MAGs_info <- read.csv(
  file.path("results", "03_mags", "intermediate", "MAGs_info.csv"),
  stringsAsFactors = FALSE
)
count_num <- table(MAGs_info$Trend)
print(count_num)

# --- Shared settings ---

group_levels <- c(
  "ice",
  "High_Glacial_sed",
  "Mid_Glacial_sed",
  "Low_Glacial_sed"
)
trend_levels <- c("Glaciophiles", "Glaciophobes", "Non-sig")
desired_order <- group_levels
boxplot_trend_levels <- c("Glaciophiles", "Non-sig", "Glaciophobes")
trend_colors <- c(
  "Glaciophiles" = "#D53E4F",
  "Glaciophobes" = "#3288BD",
  "Non-sig" = "#777777"
)
percent_axis <- function(x) paste0(x, "%")

# --- Prepare ice and sediment relative-abundance data ---

# Exclude water and cryoconite, then renormalize retained samples to sum to one.
ps_ice_sed <- subset_samples(
  PAMIR_MAGs_rela,
  !(Source %in% c("water", "cryoconite"))
)
ps_ice_sed <- prune_taxa(taxa_sums(ps_ice_sed) > 0, ps_ice_sed)
ps_ice_sed_rel <- transform_sample_counts(
  ps_ice_sed,
  function(x) x / sum(x)
)

# --- Aggregate individual MAGs into GI-response groups ---

# Index response labels by MAG identifier.
MAGs_info_indexed <- MAGs_info %>%
  column_to_rownames(var = "MAGs")

otu_df <- as.data.frame(otu_table(ps_ice_sed_rel))

# Orient the abundance matrix with MAGs in rows.
if (!taxa_are_rows(ps_ice_sed_rel)) {
  otu_df <- t(otu_df)
}

# Attach the inferred GI-response group to each MAG.
otu_df$Trend <- MAGs_info_indexed[rownames(otu_df), "Trend"]

# Sum relative abundances within each response group and sample.
otu_merged_df <- otu_df %>%
  group_by(Trend) %>%
  summarise(across(everything(), sum)) %>%
  tibble::column_to_rownames("Trend")

# Rebuild a phyloseq object with the three response groups as taxa.
OTU_new <- otu_table(as.matrix(otu_merged_df), taxa_are_rows = TRUE)

TAX_new_df <- data.frame(
  Trend = rownames(otu_merged_df),
  row.names = rownames(otu_merged_df)
)
TAX_new <- tax_table(as.matrix(TAX_new_df))

SAM_old <- sample_data(ps_ice_sed_rel)

ps_trend_merged <- phyloseq(OTU_new, TAX_new, SAM_old)

print(ps_trend_merged)

# Convert grouped abundances to long format and set plotting order.
plot_data <- psmelt(ps_trend_merged)

desired_order <- c(
  "ice",
  "High_Glacial_sed",
  "Mid_Glacial_sed",
  "Low_Glacial_sed"
)
plot_data$Group_GI <- factor(plot_data$Group_GI, levels = desired_order)
plot_data$OTU <- factor(plot_data$OTU, levels = boxplot_trend_levels)

# --- Ordered-stage Spearman analysis ---

# Encode the four ecological stages as 1-4 for two-sided Spearman tests, then
# apply BH correction across the three response groups.

spearman_trend_results <- plot_data %>%
  mutate(Stage = match(as.character(Group_GI), desired_order)) %>%
  filter(!is.na(Stage), !is.na(Abundance), !is.na(OTU)) %>%
  group_by(OTU) %>%
  group_modify(~ {
    trend_test <- cor.test(
      .x$Stage,
      .x$Abundance,
      method = "spearman",
      alternative = "two.sided",
      exact = FALSE
    )

    tibble(
      n = nrow(.x),
      rho = unname(trend_test$estimate),
      p_value = trend_test$p.value
    )
  }) %>%
  ungroup() %>%
  mutate(
    p_adjusted_BH = p.adjust(p_value, method = "BH"),
    significance = case_when(
      p_adjusted_BH <= 0.001 ~ "***",
      p_adjusted_BH <= 0.01 ~ "**",
      p_adjusted_BH <= 0.05 ~ "*",
      TRUE ~ "ns"
    )
  ) %>%
  arrange(match(as.character(OTU), boxplot_trend_levels))

write.csv(
  spearman_trend_results,
  file.path(out_dir, "Ordered_Stage_Spearman_Trend.csv"),
  row.names = FALSE
)

# --- Prepare figure annotations ---

format_spearman_fdr <- function(x) {
  ifelse(
    x < 0.001,
    format(x, scientific = TRUE, digits = 3),
    formatC(x, format = "f", digits = 3)
  )
}

spearman_annotation <- spearman_trend_results %>%
  mutate(
    Group_GI = factor("ice", levels = desired_order),
    annotation_y = case_when(
      OTU == "Glaciophiles" ~ 1.00,
      OTU == "Non-sig" ~ 0.94,
      OTU == "Glaciophobes" ~ 0.88
    ),
    label = paste0(
      as.character(OTU),
      ": \u03c1 = ", sprintf("%.3f", rho),
      ", FDR = ", format_spearman_fdr(p_adjusted_BH)
    )
  )

# --- Figure 3c: response-group abundance ---

# Plot all samples together; diamonds denote response-group means.
p_single <- ggplot(plot_data, aes(x = Group_GI, y = Abundance, fill = OTU)) +
  geom_boxplot(
    alpha = 0.8,
    outlier.shape = NA,
    position = position_dodge(0.8),
    width = 0.7
  ) +
  geom_point(
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
    size = 1,
    alpha = 0.4,
    show.legend = FALSE
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    aes(group = OTU),
    position = position_dodge(0.8),
    shape = 23,
    size = 2.8,
    stroke = 0.5,
    color = "black",
    fill = "black",
    show.legend = FALSE
  ) +
  geom_text(
    data = spearman_annotation,
    aes(
      x = Group_GI,
      y = annotation_y,
      label = label,
      color = OTU
    ),
    inherit.aes = FALSE,
    hjust = 0,
    size = 3.4,
    fontface = "bold",
    show.legend = FALSE
  ) +
  scale_fill_manual(values = trend_colors) +
  scale_color_manual(values = trend_colors, guide = "none") +
  scale_x_discrete(
    labels = c("Ice", "High", "Mid", "Low")
  ) +
  scale_y_continuous(
    limits = c(-0.05, 1.05),
    breaks = seq(0, 1, by = 0.25),
    labels = scales::label_percent(accuracy = 1),
    expand = c(0, 0)
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 11,
      color = "black"
    ),
    axis.title = element_text(size = 12, face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.position = "top"
  ) +
  labs(
    x = "",
    y = "Total Relative Abundance",
    fill = ""
  )

# --- Save figure ---

ggsave(
  file.path(out_dir, "Group_Relative_Abundance.pdf"),
  p_single,
  width = 7,
  height = 6,
  device = grDevices::cairo_pdf,
  bg = "white"
)
