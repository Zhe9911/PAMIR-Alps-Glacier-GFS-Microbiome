# GAMM and SourceTracker analysis
# Models community dissimilarity, SourceTracker mixing proportions, and
# bacterial abundance along the Glacial Index gradient.

rm(list = ls())
gc()
graphics.off()

ps_path <- file.path(
  "results", "01_16s", "intermediate", "PAMIR_16S_final.rds"
)
dissimilarity_dir <- file.path("results", "01_16s", "dissimilarity")
input_dir <- file.path("data", "processed", "16s")
sourcetracker_input_dir <- file.path(
  "results", "generated", "02_community", "sourcetracker"
)
sourcetracker_result_path <- file.path(
  input_dir, "sourcetracker", "all_results_depth_10000.csv"
)
result_dir <- file.path("results", "02_community", "gamm")

dir.create(sourcetracker_input_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

library(phyloseq)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(vegan)

library(mgcv)
library(MuMIn)
library(DHARMa)
if (!requireNamespace("mgcViz", quietly = TRUE)) {
  stop("mgcViz is required by DHARMa for mgcv GAM diagnostics")
}

ps <- readRDS(ps_path)
ps
sample_metadata_df <- data.frame(sample_data(ps))

# Assign each sediment sample its dissimilarity-to-ice metrics. The pairwise
# tables are direction-free, so the sediment sample can appear in either column.
extract_sediment_metrics_to_ice <- function(df, value_cols) {
  df %>%
    filter(
      (Sample1_type == "ice" & Sample2_type == "sed") |
        (Sample1_type == "sed" & Sample2_type == "ice")
    ) %>%
    mutate(
      Sample_ID = case_when(
        Sample1_type == "sed" & Sample2_type == "ice" ~ Sample1,
        Sample1_type == "ice" & Sample2_type == "sed" ~ Sample2,
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(Sample_ID)) %>%
    group_by(Sample_ID) %>%
    summarise(
      across(all_of(value_cols), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    )
}

count_complete_cases <- function(data, vars) {
  sum(complete.cases(data[, vars]))
}

read_env_input_data <- function(path) {
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

extract_sediment_env_by_site <- function(df, env_vars, suffix_pattern = "_[A-Za-z]+$") {
  required_cols <- c("Source", "Sample_ID", env_vars)
  missing_cols <- setdiff(required_cols, names(df))

  if (length(missing_cols) > 0) {
    stop(
      paste(
        "Missing required columns in env_GI.csv:",
        paste(missing_cols, collapse = ", ")
      )
    )
  }

  df %>%
    filter(Source == "sed") %>%
    mutate(Sample_ID = sub(suffix_pattern, "", Sample_ID)) %>%
    mutate(across(all_of(env_vars), ~ as.numeric(as.character(.x)))) %>%
    group_by(Sample_ID) %>%
    summarise(
      across(all_of(env_vars), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    mutate(across(all_of(env_vars), ~ if_else(is.nan(.x), NA_real_, .x)))
}

# --- Bray-Curtis dissimilarity to ice samples ---

BC_dissimilarity_original <- read.csv(
  file.path(dissimilarity_dir, "within_glacier_ice_sed_bray_curtis.csv"),
  sep = ","
)

# Keep all ice-sediment pairs regardless of row direction, assign the metric to
# the sediment sample, and average across multiple ice partners if present.
BC_dissimilarity_df <- extract_sediment_metrics_to_ice(
  BC_dissimilarity_original,
  value_cols = "Dissimilarity"
) %>%
  rename(BC_dissimilarity_to_ice = Dissimilarity)

combined_df <- sample_metadata_df %>%
  left_join(BC_dissimilarity_df %>% select(Sample_ID, BC_dissimilarity_to_ice),
            by = "Sample_ID")

# --- Weighted UniFrac dissimilarity to ice samples ---

UF_dissimilarity_original <- read.csv(
  file.path(dissimilarity_dir, "within_glacier_ice_sed_weighted_unifrac.csv"),
  sep = ","
)

# Use the same direction-agnostic handling as for Bray-Curtis.
UF_dissimilarity_df <- extract_sediment_metrics_to_ice(
  UF_dissimilarity_original,
  value_cols = "Dissimilarity"
) %>%
  rename(UF_dissimilarity_to_ice = Dissimilarity)

combined_df <- combined_df %>%
  left_join(UF_dissimilarity_df %>% select(Sample_ID, UF_dissimilarity_to_ice),
            by = "Sample_ID")

# --- SourceTracker results ---

# Convert the count table to BIOM format and export metadata for SourceTracker.
library(biomformat)
biom_path <- file.path(sourcetracker_input_dir, "ST_otu_table.biom")
biom_object <- make_biom(data = as(otu_table(ps), "matrix"))
write_biom(biom_object, biom_path)

meta_ST <- as.data.frame(as.matrix(sample_data(ps)))
write.csv(
  meta_ST,
  file = file.path(sourcetracker_input_dir, "ST_metadata.csv"),
  row.names = TRUE
)

# SourceTracker2 can be rerun with the companion Python script.
# This analysis reads the frozen mixing-proportion table used for Figure 2e.

Proportion <- read.csv(sourcetracker_result_path, sep = ",")
Proportion <- Proportion %>%
  select(-Glacier, -Group_GI) %>%
  pivot_wider(
    id_cols = SampleID,
    names_from = Source,
    values_from = Proportion
  ) %>%
  rename(
    Sample_ID = SampleID,
    Ice_proportion = ice,
    Un_proportion = Unknown
  )

combined_df <- combined_df %>% left_join(Proportion, by = "Sample_ID")

# --- Merge replicate samples ---

combined_df_sed <- subset(combined_df, Source == "sed")

if (!"BA" %in% names(combined_df_sed)) {
  stop("BA was not found in sediment sample metadata.")
}

combined_df_sed$BA <- as.numeric(as.character(combined_df_sed$BA))

combined_df_sed$Sample <- sub("_[A-Za-z]+$", "", combined_df_sed$Sample_ID)

combined_final <- combined_df_sed %>%
  group_by(Sample) %>%
  summarise(
    # Summarize numeric variables across replicate measurements.
    across(
      where(is.numeric),
      list(mean = ~mean(., na.rm = TRUE), sd = ~sd(., na.rm = TRUE)),
      .names = "{.col}_{.fn}"
    ),
    # Carry over the first value for non-numeric metadata fields.
    across(
      !where(is.numeric),
      ~ first(.)
    ),
    .groups = "drop"
  )

combined_final <- combined_final %>%
  select(-Sample_ID) %>%
  # Geographic coordinates remain as metadata but are not used as predictors.
  # Elevation is renamed because it is an explicit candidate predictor.
  rename(
    Sample_ID = Sample,
    GI = GI_mean,
    Elevation = Elevation_mean,
    distance = distance_to_glacier_snout_mean,
    gl_size = gl_size_mean
  ) %>%
  # Convert distance from meters to kilometers.
  mutate(distance = distance / 1000)

combined_final$gl_name <- as.factor(combined_final$gl_name)
combined_final$Region <- as.factor(combined_final$Region)
# Ensure GI is stored as numeric.
combined_final$GI <- as.numeric(as.character(combined_final$GI))

# --- Add key sediment environmental variables from env_GI.csv ---

env_main_vars <- c("temp", "turb", "DOC_mean", "pH", "NO3")

env_main_df <- read_env_input_data(file.path(input_dir, "env_GI.csv")) %>%
  extract_sediment_env_by_site(env_vars = env_main_vars)

combined_final <- combined_final %>%
  left_join(env_main_df, by = "Sample_ID")

env_join_summary <- tibble(
  variable = env_main_vars,
  non_missing_after_join = map_int(env_main_vars, ~ sum(!is.na(combined_final[[.x]])))
)

print(env_join_summary)

metric_sample_sizes <- tibble(
  metric = c("BC_dissimilarity_to_ice", "UF_dissimilarity_to_ice", "Ice_proportion"),
  sediment_replicates_with_values = c(
    sum(!is.na(combined_df_sed$BC_dissimilarity_to_ice)),
    sum(!is.na(combined_df_sed$UF_dissimilarity_to_ice)),
    sum(!is.na(combined_df_sed$Ice_proportion))
  ),
  merged_samples_with_values = c(
    sum(!is.na(combined_final$BC_dissimilarity_to_ice_mean)),
    sum(!is.na(combined_final$UF_dissimilarity_to_ice_mean)),
    sum(!is.na(combined_final$Ice_proportion_mean))
  )
)

model_sample_size_specs <- tribble(
  ~model,                        ~required_vars,
  "mod_BCdis_empty_random",      c("BC_dissimilarity_to_ice_mean", "gl_name"),
  "mod_BCdis_region_baseline",   c("BC_dissimilarity_to_ice_mean", "Region", "gl_name"),
  "mod_BCdis_1",                 c("BC_dissimilarity_to_ice_mean", "Region", "GI", "gl_name"),
  "mod_BCdis_2",                 c("BC_dissimilarity_to_ice_mean", "Region", "distance", "gl_name"),
  "mod_BCdis_3",                 c("BC_dissimilarity_to_ice_mean", "Region", "gl_size", "gl_name"),
  "mod_BCdis_4",                 c(
    "BC_dissimilarity_to_ice_mean", "Region", "Elevation", "gl_name"
  ),
  "mod_BCdis_5",                 c("BC_dissimilarity_to_ice_mean", "Region", "temp", "gl_name"),
  "mod_BCdis_6",                 c("BC_dissimilarity_to_ice_mean", "Region", "pH", "gl_name"),
  "mod_BCdis_7",                 c(
    "BC_dissimilarity_to_ice_mean", "Region", "Chla_mean", "gl_name"
  ),
  "mod_BCdis_8",                 c(
    "BC_dissimilarity_to_ice_mean", "Region", "distance", "gl_size", "gl_name"
  ),
  "mod_BCdis_GI_only",           c("BC_dissimilarity_to_ice_mean", "GI", "gl_name"),
  "mod_BCdis_9",                 c("BC_dissimilarity_to_ice_mean", "Region", "GI", "gl_name"),
  "mod_BCdis_final_GI",          c("BC_dissimilarity_to_ice_mean", "Region", "GI", "gl_name"),
  "mod_BCdis_final_GI_no_random", c("BC_dissimilarity_to_ice_mean", "Region", "GI"),
  "mod_BCdis_baseline_Region",   c("BC_dissimilarity_to_ice_mean", "Region", "gl_name"),
  "mod_UFdis_empty_random",      c("UF_dissimilarity_to_ice_mean", "gl_name"),
  "mod_UFdis_region_baseline",   c("UF_dissimilarity_to_ice_mean", "Region", "gl_name"),
  "mod_UFdis_1",                 c("UF_dissimilarity_to_ice_mean", "Region", "GI", "gl_name"),
  "mod_UFdis_2",                 c("UF_dissimilarity_to_ice_mean", "Region", "distance", "gl_name"),
  "mod_UFdis_3",                 c("UF_dissimilarity_to_ice_mean", "Region", "gl_size", "gl_name"),
  "mod_UFdis_4",                 c(
    "UF_dissimilarity_to_ice_mean", "Region", "Elevation", "gl_name"
  ),
  "mod_UFdis_5",                 c("UF_dissimilarity_to_ice_mean", "Region", "temp", "gl_name"),
  "mod_UFdis_6",                 c("UF_dissimilarity_to_ice_mean", "Region", "pH", "gl_name"),
  "mod_UFdis_7",                 c(
    "UF_dissimilarity_to_ice_mean", "Region", "Chla_mean", "gl_name"
  ),
  "mod_UFdis_8",                 c(
    "UF_dissimilarity_to_ice_mean", "Region", "distance", "gl_size", "gl_name"
  ),
  "mod_UFdis_GI_only",           c("UF_dissimilarity_to_ice_mean", "GI", "gl_name"),
  "mod_UFdis_9",                 c("UF_dissimilarity_to_ice_mean", "Region", "GI", "gl_name"),
  "mod_UFdis_final_GI",          c("UF_dissimilarity_to_ice_mean", "Region", "GI", "gl_name"),
  "mod_UFdis_final_GI_no_random", c("UF_dissimilarity_to_ice_mean", "Region", "GI"),
  "mod_UFdis_baseline_Region",   c("UF_dissimilarity_to_ice_mean", "Region", "gl_name"),
  "mod_prop_empty_random",       c("Ice_proportion_mean", "gl_name"),
  "mod_prop_region_baseline",    c("Ice_proportion_mean", "Region", "gl_name"),
  "mod_prop_1",                  c("Ice_proportion_mean", "Region", "GI", "gl_name"),
  "mod_prop_2",                  c("Ice_proportion_mean", "Region", "distance", "gl_name"),
  "mod_prop_3",                  c("Ice_proportion_mean", "Region", "gl_size", "gl_name"),
  "mod_prop_4",                  c("Ice_proportion_mean", "Region", "Elevation", "gl_name"),
  "mod_prop_5",                  c("Ice_proportion_mean", "Region", "temp", "gl_name"),
  "mod_prop_6",                  c("Ice_proportion_mean", "Region", "pH", "gl_name"),
  "mod_prop_7",                  c("Ice_proportion_mean", "Region", "Chla_mean", "gl_name"),
  "mod_prop_8",                  c(
    "Ice_proportion_mean", "Region", "distance", "gl_size", "gl_name"
  ),
  "mod_prop_GI_only",            c("Ice_proportion_mean", "GI", "gl_name"),
  "mod_prop_9",                  c("Ice_proportion_mean", "Region", "GI", "gl_name"),
  "mod_prop_final_GI",           c("Ice_proportion_mean", "Region", "GI", "gl_name"),
  "mod_prop_final_GI_no_random",  c("Ice_proportion_mean", "Region", "GI"),
  "mod_prop_baseline_Region",    c("Ice_proportion_mean", "Region", "gl_name")
)

model_sample_sizes <- model_sample_size_specs %>%
  mutate(n = map_int(required_vars, ~ count_complete_cases(combined_final, .x))) %>%
  select(model, n)

print(metric_sample_sizes)
print(model_sample_sizes)

# --- GAMM analyses ---

# Build a prediction grid across the observed GI range in each region.
new_data <- expand.grid(
  GI = seq(min(combined_final$GI), max(combined_final$GI), length.out = 200),
  Region = levels(as.factor(combined_final$Region))
)
new_data$gl_name <- combined_final$gl_name[1]

# --- Bray-Curtis GAMM: model selection ---

# Screen linear and nonlinear structures with GAM.

# Fit the empty random-effects model.
mod_BCdis_empty_random <- gam(
  formula = BC_dissimilarity_to_ice_mean ~ 1 + s(gl_name, bs = "re"),
  family = betar(link = "logit"),
  data = combined_final,
  method = "ML"
)

# Fit the Region-only baseline model.
mod_BCdis_region_baseline <- gam(
  formula = BC_dissimilarity_to_ice_mean ~ Region + s(gl_name, bs = "re"),
  family = betar(link = "logit"),
  data = combined_final,
  method = "ML"
)

# Fit the candidate environmental models.
mod_BCdis_1 <- gam(
  formula = BC_dissimilarity_to_ice_mean ~ Region + s(GI, k = 5) + s(gl_name, bs = "re"),
  family = betar(link = "logit"),
  data = combined_final,
  method = "ML"
)

mod_BCdis_2 <- gam(
  formula = BC_dissimilarity_to_ice_mean ~ Region + s(distance, k = 5) + s(gl_name, bs = "re"),
  family = betar(link = "logit"),
  data = combined_final,
  method = "ML"
)

mod_BCdis_3 <- gam(
  formula = BC_dissimilarity_to_ice_mean ~ Region + s(gl_size, k = 5) + s(gl_name, bs = "re"),
  family = betar(link = "logit"),
  data = combined_final,
  method = "ML"
)

mod_BCdis_4 <- gam(
  formula = BC_dissimilarity_to_ice_mean ~ Region + s(Elevation, k = 5) + s(gl_name, bs = "re"),
  family = betar(link = "logit"),
  data = combined_final,
  method = "ML"
)

mod_BCdis_5 <- gam(
  formula = BC_dissimilarity_to_ice_mean ~ Region + s(temp, k = 5) + s(gl_name, bs = "re"),
  family = betar(link = "logit"),
  data = combined_final,
  method = "ML"
)

mod_BCdis_6 <- gam(
  formula = BC_dissimilarity_to_ice_mean ~ Region + s(pH, k = 5) + s(gl_name, bs = "re"),
  family = betar(link = "logit"),
  data = combined_final,
  method = "ML"
)

mod_BCdis_7 <- gam(
  formula = BC_dissimilarity_to_ice_mean ~ Region + s(Chla_mean, k = 5) + s(gl_name, bs = "re"),
  family = betar(link = "logit"),
  data = combined_final,
  method = "ML"
)

mod_BCdis_8 <- gam(
  formula = BC_dissimilarity_to_ice_mean ~ Region +
    s(distance, k = 5) + s(gl_size, k = 5) + s(gl_name, bs = "re"),
  family = betar(link = "logit"),
  data = combined_final,
  method = "ML"
)

AIC(
  mod_BCdis_empty_random, mod_BCdis_region_baseline,
  mod_BCdis_1, mod_BCdis_2, mod_BCdis_3, mod_BCdis_4,
  mod_BCdis_5, mod_BCdis_6, mod_BCdis_7, mod_BCdis_8
)
BIC(
  mod_BCdis_empty_random, mod_BCdis_region_baseline,
  mod_BCdis_1, mod_BCdis_2, mod_BCdis_3, mod_BCdis_4,
  mod_BCdis_5, mod_BCdis_6, mod_BCdis_7, mod_BCdis_8
)

summary(mod_BCdis_empty_random)
summary(mod_BCdis_region_baseline)
summary(mod_BCdis_1)
summary(mod_BCdis_5)

# `mod_BCdis_1` has the strongest AIC/BIC support among the candidate
# environmental models, and the fitted GI smooth is effectively linear
# (`edf = 1.00`). Next, test the Region × GI interaction.

mod_BCdis_GI_only <- gam(
  formula = BC_dissimilarity_to_ice_mean ~ GI + s(gl_name, bs = "re"),
  family = betar(link = "logit"),
  data = combined_final,
  method = "ML"
)

mod_BCdis_9 <- gam(
  formula = BC_dissimilarity_to_ice_mean ~ Region * GI + s(gl_name, bs = "re"),
  family = betar(link = "logit"),
  data = combined_final,
  method = "ML"
)

summary(mod_BCdis_9)

AIC(mod_BCdis_empty_random, mod_BCdis_region_baseline,
    mod_BCdis_GI_only, mod_BCdis_1, mod_BCdis_9)
BIC(mod_BCdis_empty_random, mod_BCdis_region_baseline,
    mod_BCdis_GI_only, mod_BCdis_1, mod_BCdis_9)

# The Region × GI interaction is not significant and is not supported by
# AIC/BIC. Retain the additive model
# `BC dissimilarity ~ Region + GI + (1 | gl_name)` and refit with REML.

# --- Bray-Curtis GAMM: final model ---

mod_BCdis_final_GI <- gam(
  BC_dissimilarity_to_ice_mean ~ Region + GI + s(gl_name, bs = "re"),
  family = betar(link = "logit"),
  data = combined_final,
  method = "REML"
)

summary(mod_BCdis_final_GI)

mod_BCdis_final_GI_no_random <- gam(
  BC_dissimilarity_to_ice_mean ~ Region + GI,
  family = betar(link = "logit"),
  data = combined_final,
  method = "REML"
)

summary(mod_BCdis_final_GI_no_random)

# --- Bray-Curtis GAMM: diagnostics ---

# Keep non-exported DHARMa plots off the default batch graphics device.
grDevices::pdf(file = NULL)
sim_res <- simulateResiduals(mod_BCdis_final_GI, n = 1000)
plotQQunif(sim_res) # Plot the uniform QQ diagnostic.
testDispersion(sim_res)  # Test residual dispersion.
plotResiduals(sim_res, combined_final$GI)

# The DHARMa residual diagnostics show no obvious deviation from the expected
# residual distribution by visual inspection.

# --- Bray-Curtis GAMM: incremental GI contribution ---

# Refit the baseline with REML = TRUE so the R2 comparison matches the final model.
mod_BCdis_baseline_Region <- gam(
  BC_dissimilarity_to_ice_mean ~ Region + s(gl_name, bs = "re"),
  family = betar(link = "logit"),
  data = combined_final,
  method = "REML"
)

r2_BCdis_final <- summary(mod_BCdis_final_GI)$r.sq
r2_BCdis_baseline <- summary(mod_BCdis_baseline_Region)$r.sq

GI_BCdis_increment <- r2_BCdis_final - r2_BCdis_baseline
cat("Final Model Adj. R-squared:", round(r2_BCdis_final, 4), "\n")
cat("Baseline Model Adj. R-squared:", round(r2_BCdis_baseline, 4), "\n")
cat("Incremental contribution of GI:", round(GI_BCdis_increment, 4), "\n")

dev_BCdis_final <- summary(mod_BCdis_final_GI)$dev.expl
dev_BCdis_baseline <- summary(mod_BCdis_baseline_Region)$dev.expl
cat("Incremental Deviance Explained by GI:", round(dev_BCdis_final - dev_BCdis_baseline, 4), "\n")

# --- Bray-Curtis GAMM: visualization ---

preds_BCdis <- predict(
  mod_BCdis_final_GI,
  newdata = new_data,
  type = "link",
  se.fit = TRUE,
  exclude = "s(gl_name)"
)

# Compute confidence intervals on the link scale, then back-transform once.
new_data$pred_BCdis  <- plogis(preds_BCdis$fit)
new_data$lower_BCdis <- plogis(preds_BCdis$fit - 1.96 * preds_BCdis$se.fit)
new_data$upper_BCdis <- plogis(preds_BCdis$fit + 1.96 * preds_BCdis$se.fit)

BCdis_GI_plot <- ggplot() +
  geom_ribbon(data = new_data,
              aes(x = GI,
                  ymin = lower_BCdis,
                  ymax = upper_BCdis,
                  fill = Region),
              alpha = 0.2) +

  geom_errorbar(data = combined_final,
                aes(x = GI,
                    ymin = BC_dissimilarity_to_ice_mean - BC_dissimilarity_to_ice_sd,
                    ymax = BC_dissimilarity_to_ice_mean + BC_dissimilarity_to_ice_sd,
                    color = Region),
                width = 0, alpha = 0.5, size = 0.5) +
  geom_point(data = combined_final,
             aes(x = GI, y = BC_dissimilarity_to_ice_mean, color = Region),
             size = 3, alpha = 0.8) +
  geom_line(data = new_data,
            aes(x = GI, y = pred_BCdis, color = Region),
            linewidth = 1.2) +

  scale_color_manual(values = c("#E69F00", "#56B4E9")) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9")) +
  scale_x_reverse() +
  labs(x = "Glacial Index (GI)",
       y = "Bray-Curtis Dissimilarity (to Ice)",
       title = "Bray-Curtis Dissimilarity ~ GI + Region") +
  theme_bw() +
  theme(legend.position = "top",
        text = element_text(size = 12))

# --- Weighted UniFrac GAMM: model selection ---

# Screen linear and nonlinear structures with GAM.

# Fit the empty random-effects model.
mod_UFdis_empty_random <- gam(
  formula = UF_dissimilarity_to_ice_mean ~ 1 + s(gl_name, bs = "re"),
  family = betar(link = "logit"),
  data = combined_final,
  method = "ML"
)

# Fit the Region-only baseline model.
mod_UFdis_region_baseline <- gam(
  formula = UF_dissimilarity_to_ice_mean ~ Region + s(gl_name, bs = "re"),
  family = betar(link = "logit"),
  data = combined_final,
  method = "ML"
)

# Fit the candidate environmental models.
mod_UFdis_1 <- gam(
  formula = UF_dissimilarity_to_ice_mean ~ Region + s(GI, k = 5) + s(gl_name, bs = "re"),
  family = betar(link = "logit"),
  data = combined_final,
  method = "ML"
)

mod_UFdis_2 <- gam(
  formula = UF_dissimilarity_to_ice_mean ~ Region + s(distance, k = 5) + s(gl_name, bs = "re"),
  family = betar(link = "logit"),
  data = combined_final,
  method = "ML"
)

mod_UFdis_3 <- gam(
  formula = UF_dissimilarity_to_ice_mean ~ Region + s(gl_size, k = 5) + s(gl_name, bs = "re"),
  family = betar(link = "logit"),
  data = combined_final,
  method = "ML"
)

mod_UFdis_4 <- gam(
  formula = UF_dissimilarity_to_ice_mean ~ Region + s(Elevation, k = 5) + s(gl_name, bs = "re"),
  family = betar(link = "logit"),
  data = combined_final,
  method = "ML"
)

mod_UFdis_5 <- gam(
  formula = UF_dissimilarity_to_ice_mean ~ Region + s(temp, k = 5) + s(gl_name, bs = "re"),
  family = betar(link = "logit"),
  data = combined_final,
  method = "ML"
)

mod_UFdis_6 <- gam(
  formula = UF_dissimilarity_to_ice_mean ~ Region + s(pH, k = 5) + s(gl_name, bs = "re"),
  family = betar(link = "logit"),
  data = combined_final,
  method = "ML"
)

mod_UFdis_7 <- gam(
  formula = UF_dissimilarity_to_ice_mean ~ Region + s(Chla_mean, k = 5) + s(gl_name, bs = "re"),
  family = betar(link = "logit"),
  data = combined_final,
  method = "ML"
)

mod_UFdis_8 <- gam(
  formula = UF_dissimilarity_to_ice_mean ~ Region +
    s(distance, k = 5) + s(gl_size, k = 5) + s(gl_name, bs = "re"),
  family = betar(link = "logit"),
  data = combined_final,
  method = "ML"
)

AIC(
  mod_UFdis_empty_random, mod_UFdis_region_baseline,
  mod_UFdis_1, mod_UFdis_2, mod_UFdis_3, mod_UFdis_4,
  mod_UFdis_5, mod_UFdis_6, mod_UFdis_7, mod_UFdis_8
)
BIC(
  mod_UFdis_empty_random, mod_UFdis_region_baseline,
  mod_UFdis_1, mod_UFdis_2, mod_UFdis_3, mod_UFdis_4,
  mod_UFdis_5, mod_UFdis_6, mod_UFdis_7, mod_UFdis_8
)

summary(mod_UFdis_empty_random)
summary(mod_UFdis_region_baseline)
summary(mod_UFdis_1)
summary(mod_UFdis_5)

# `mod_UFdis_1` has the strongest AIC/BIC support among the candidate
# environmental models, and the fitted GI smooth is effectively linear
# (`edf = 1.00`). Next, test the Region × GI interaction.

mod_UFdis_GI_only <- gam(
  formula = UF_dissimilarity_to_ice_mean ~ GI + s(gl_name, bs = "re"),
  family = betar(link = "logit"),
  data = combined_final,
  method = "ML"
)

mod_UFdis_9 <- gam(
  formula = UF_dissimilarity_to_ice_mean ~ Region * GI + s(gl_name, bs = "re"),
  family = betar(link = "logit"),
  data = combined_final,
  method = "ML"
)

summary(mod_UFdis_9)

AIC(mod_UFdis_empty_random, mod_UFdis_region_baseline,
    mod_UFdis_GI_only, mod_UFdis_1, mod_UFdis_9)
BIC(mod_UFdis_empty_random, mod_UFdis_region_baseline,
    mod_UFdis_GI_only, mod_UFdis_1, mod_UFdis_9)

# The Region × GI interaction is not significant and is not supported by
# AIC/BIC. The additive `Region + GI` model is retained for consistency with the
# sediment GAMM framework; it has nearly identical AIC to the GI-only model and
# slightly stronger BIC support. Refit the additive model with REML.

# --- Weighted UniFrac GAMM: final model ---

mod_UFdis_final_GI <- gam(
  UF_dissimilarity_to_ice_mean ~ Region + GI + s(gl_name, bs = "re"),
  family = betar(link = "logit"),
  data = combined_final,
  method = "REML"
)

summary(mod_UFdis_final_GI)

mod_UFdis_final_GI_no_random <- gam(
  UF_dissimilarity_to_ice_mean ~ Region + GI,
  family = betar(link = "logit"),
  data = combined_final,
  method = "REML"
)

summary(mod_UFdis_final_GI_no_random)
# --- Weighted UniFrac GAMM: diagnostics ---

sim_res <- simulateResiduals(mod_UFdis_final_GI, n = 1000)
plotQQunif(sim_res) # Plot the uniform QQ diagnostic.
testDispersion(sim_res)  # Test residual dispersion.
plotResiduals(sim_res, combined_final$GI)

# The DHARMa residual diagnostics show no obvious deviation from the expected
# residual distribution by visual inspection.

# --- Weighted UniFrac GAMM: incremental GI contribution ---

# Refit the baseline with REML = TRUE so the R2 comparison matches the final model.
mod_UFdis_baseline_Region <- gam(
  UF_dissimilarity_to_ice_mean ~ Region + s(gl_name, bs = "re"),
  family = betar(link = "logit"),
  data = combined_final,
  method = "REML"
)

r2_UFdis_final <- summary(mod_UFdis_final_GI)$r.sq
r2_UFdis_baseline <- summary(mod_UFdis_baseline_Region)$r.sq

GI_UFdis_increment <- r2_UFdis_final - r2_UFdis_baseline
cat("Final Model Adj. R-squared:", round(r2_UFdis_final, 4), "\n")
cat("Baseline Model Adj. R-squared:", round(r2_UFdis_baseline, 4), "\n")
cat("Incremental contribution of GI:", round(GI_UFdis_increment, 4), "\n")

dev_UFdis_final <- summary(mod_UFdis_final_GI)$dev.expl
dev_UFdis_baseline <- summary(mod_UFdis_baseline_Region)$dev.expl
cat("Incremental Deviance Explained by GI:", round(dev_UFdis_final - dev_UFdis_baseline, 4), "\n")

# --- Weighted UniFrac GAMM: visualization ---

preds_UFdis <- predict(
  mod_UFdis_final_GI,
  newdata = new_data,
  type = "link",
  se.fit = TRUE,
  exclude = "s(gl_name)"
)

# Compute confidence intervals on the link scale, then back-transform once.
new_data$pred_UFdis  <- plogis(preds_UFdis$fit)
new_data$lower_UFdis <- plogis(preds_UFdis$fit - 1.96 * preds_UFdis$se.fit)
new_data$upper_UFdis <- plogis(preds_UFdis$fit + 1.96 * preds_UFdis$se.fit)

UFdis_GI_plot <- ggplot() +
  geom_ribbon(data = new_data,
              aes(x = GI,
                  ymin = lower_UFdis,
                  ymax = upper_UFdis,
                  fill = Region),
              alpha = 0.2) +

  geom_errorbar(data = combined_final,
                aes(x = GI,
                    ymin = UF_dissimilarity_to_ice_mean - UF_dissimilarity_to_ice_sd,
                    ymax = UF_dissimilarity_to_ice_mean + UF_dissimilarity_to_ice_sd,
                    color = Region),
                width = 0, alpha = 0.5, size = 0.5) +
  geom_point(data = combined_final,
             aes(x = GI, y = UF_dissimilarity_to_ice_mean, color = Region),
             size = 3, alpha = 0.8) +
  geom_line(data = new_data,
            aes(x = GI, y = pred_UFdis, color = Region),
            linewidth = 1.2) +
  coord_cartesian(ylim = c(0, 0.3)) +

  scale_color_manual(values = c("#E69F00", "#56B4E9")) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9")) +
  scale_x_reverse() +
  labs(x = "Glacial Index (GI)",
       y = "Weighted UniFrac Dissimilarity (to Ice)",
       title = "Weighted UniFrac Dissimilarity ~ GI + Region") +
  theme_bw() +
  theme(legend.position = "top",
        text = element_text(size = 12))

# --- SourceTracker ice-source proportion GAMM: model selection ---

# Screen linear and nonlinear structures with GAM.

# Fit the empty random-effects model.
mod_prop_empty_random <- gam(
  formula = Ice_proportion_mean ~ 1 + s(gl_name, bs = "re"),
  family = betar(link = "logit"),
  data = combined_final,
  method = "ML"
)

# Fit the Region-only baseline model.
mod_prop_region_baseline <- gam(
  formula = Ice_proportion_mean ~ Region + s(gl_name, bs = "re"),
  family = betar(link = "logit"),
  data = combined_final,
  method = "ML"
)

# Fit the candidate environmental models.
mod_prop_1 <- gam(
  formula = Ice_proportion_mean ~ Region + s(GI, k = 5) + s(gl_name, bs = "re"),
  family = betar(link = "logit"),
  data = combined_final,
  method = "ML"
)

mod_prop_2 <- gam(
  formula = Ice_proportion_mean ~ Region + s(distance, k = 5) + s(gl_name, bs = "re"),
  family = betar(link = "logit"),
  data = combined_final,
  method = "ML"
)

mod_prop_3 <- gam(
  formula = Ice_proportion_mean ~ Region + s(gl_size, k = 5) + s(gl_name, bs = "re"),
  family = betar(link = "logit"),
  data = combined_final,
  method = "ML"
)

mod_prop_4 <- gam(
  formula = Ice_proportion_mean ~ Region + s(Elevation, k = 5) + s(gl_name, bs = "re"),
  family = betar(link = "logit"),
  data = combined_final,
  method = "ML"
)

mod_prop_5 <- gam(
  formula = Ice_proportion_mean ~ Region + s(temp, k = 5) + s(gl_name, bs = "re"),
  family = betar(link = "logit"),
  data = combined_final,
  method = "ML"
)

mod_prop_6 <- gam(
  formula = Ice_proportion_mean ~ Region + s(pH, k = 5) + s(gl_name, bs = "re"),
  family = betar(link = "logit"),
  data = combined_final,
  method = "ML"
)

mod_prop_7 <- gam(
  formula = Ice_proportion_mean ~ Region + s(Chla_mean, k = 5) + s(gl_name, bs = "re"),
  family = betar(link = "logit"),
  data = combined_final,
  method = "ML"
)

mod_prop_8 <- gam(
  formula = Ice_proportion_mean ~ Region +
    s(distance, k = 5) + s(gl_size, k = 5) + s(gl_name, bs = "re"),
  family = betar(link = "logit"),
  data = combined_final,
  method = "ML"
)

AIC(mod_prop_empty_random, mod_prop_region_baseline, mod_prop_1, mod_prop_2,
    mod_prop_3, mod_prop_4, mod_prop_5, mod_prop_6, mod_prop_7, mod_prop_8)
BIC(mod_prop_empty_random, mod_prop_region_baseline, mod_prop_1, mod_prop_2,
    mod_prop_3, mod_prop_4, mod_prop_5, mod_prop_6, mod_prop_7, mod_prop_8)

summary(mod_prop_empty_random)
summary(mod_prop_region_baseline)
summary(mod_prop_1)
summary(mod_prop_5)

# `mod_prop_1` has the strongest AIC/BIC support among the candidate
# environmental models, and the fitted GI smooth is effectively linear
# (`edf = 1.00`). Next, test the Region × GI interaction.

mod_prop_GI_only <- gam(
  formula = Ice_proportion_mean ~ GI + s(gl_name, bs = "re"),
  family = betar(link = "logit"),
  data = combined_final,
  method = "ML"
)

mod_prop_9 <- gam(
  formula = Ice_proportion_mean ~ Region * GI + s(gl_name, bs = "re"),
  family = betar(link = "logit"),
  data = combined_final,
  method = "ML"
)

summary(mod_prop_9)

AIC(mod_prop_empty_random, mod_prop_region_baseline,
    mod_prop_GI_only, mod_prop_1, mod_prop_9)
BIC(mod_prop_empty_random, mod_prop_region_baseline,
    mod_prop_GI_only, mod_prop_1, mod_prop_9)
# The Region × GI interaction is not significant and is not supported by
# AIC/BIC. Retain the additive model
# `prop ~ Region + GI + (1 | gl_name)` and refit with REML.

# --- SourceTracker ice-source proportion GAMM: final model ---

mod_prop_final_GI <- gam(
  Ice_proportion_mean ~ Region + GI + s(gl_name, bs = "re"),
  family = betar(link = "logit"),
  data = combined_final,
  method = "REML"
)

summary(mod_prop_final_GI)

mod_prop_final_GI_no_random <- gam(
  Ice_proportion_mean ~ Region + GI,
  family = betar(link = "logit"),
  data = combined_final,
  method = "REML"
)

summary(mod_prop_final_GI_no_random)

# --- SourceTracker ice-source proportion GAMM: diagnostics ---

sim_res <- simulateResiduals(mod_prop_final_GI, n = 1000)
plotQQunif(sim_res) # Plot the uniform QQ diagnostic.
testDispersion(sim_res)  # Test residual dispersion.
plotResiduals(sim_res, combined_final$GI)
grDevices::dev.off()

# The DHARMa residual diagnostics show no obvious deviation from the expected
# residual distribution by visual inspection.

# --- SourceTracker ice-source proportion GAMM: incremental GI contribution ---

# Refit the baseline with REML = TRUE so the R2 comparison matches the final model.
mod_prop_baseline_Region <- gam(
  Ice_proportion_mean ~ Region + s(gl_name, bs = "re"),
  family = betar(link = "logit"),
  data = combined_final,
  method = "REML"
)

# Extract the adjusted R-squared values from the two models.
r2_final    <- summary(mod_prop_final_GI)$r.sq
r2_baseline <- summary(mod_prop_baseline_Region)$r.sq

# Calculate the incremental contribution of GI.
GI_contribution <- r2_final - r2_baseline

cat("Final Model Adj. R-squared:", round(r2_final, 4), "\n")
cat("Baseline Model Adj. R-squared:", round(r2_baseline, 4), "\n")
cat("Incremental contribution of GI:", round(GI_contribution, 4), "\n")

# Deviance explained is a core goodness-of-fit metric for GAMs.
dev_final <- summary(mod_prop_final_GI)$dev.expl
dev_baseline <- summary(mod_prop_baseline_Region)$dev.expl
cat("Incremental Deviance Explained by GI:", round(dev_final - dev_baseline, 4), "\n")

# --- SourceTracker ice-source proportion GAMM: visualization ---

preds_prop <- predict(
  mod_prop_final_GI,
  newdata = new_data,
  type = "link",
  se.fit = TRUE,
  # Exclude the random-effect smooth in mgcv predictions.
  exclude = "s(gl_name)"
)

# Compute confidence intervals on the link scale, then back-transform once.
new_data$pred_prop  <- plogis(preds_prop$fit)
new_data$lower_prop <- plogis(preds_prop$fit - 1.96 * preds_prop$se.fit)
new_data$upper_prop <- plogis(preds_prop$fit + 1.96 * preds_prop$se.fit)

prop_GI_plot <- ggplot() +
  geom_ribbon(data = new_data,
              aes(x = GI,
                  ymin = lower_prop,
                  ymax = upper_prop,
                  fill = Region),
              alpha = 0.2) +

  geom_errorbar(data = combined_final,
                aes(x = GI,
                    ymin = Ice_proportion_mean - Ice_proportion_sd,
                    ymax = Ice_proportion_mean + Ice_proportion_sd,
                    color = Region),
                width = 0, alpha = 0.5, size = 0.5) +
  geom_point(data = combined_final,
             aes(x = GI, y = Ice_proportion_mean, color = Region),
             size = 3, alpha = 0.8) +
  geom_line(data = new_data,
            aes(x = GI, y = pred_prop, color = Region),
            linewidth = 1.2) +
  coord_cartesian(ylim = c(0, 0.8)) +

  scale_color_manual(values = c("#E69F00", "#56B4E9")) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9")) +
  scale_x_reverse() +
  labs(x = "Glacial Index (GI)",
       y = "Proportion of Ice",
       title = "Ice Proportion ~ GI + Region") +
  theme_bw() +
  theme(legend.position = "top",
        text = element_text(size = 12))

ggsave(file.path(result_dir, "figure_2c_e_gamm_correlations.pdf"),
       grid.arrange(BCdis_GI_plot, UFdis_GI_plot, prop_GI_plot, ncol = 2),
       width = 9, height = 10)

# --- Sediment bacterial-abundance GAMM ---

# This model supports the bacterial-abundance result reported in Extended Data
# Table 1.

if (!"BA_mean" %in% names(combined_final)) {
  stop("BA_mean was not found in combined_final. Check whether BA is numeric in sample metadata.")
}

ba_sed_df <- combined_final %>%
  transmute(
    Sample_ID,
    Region,
    gl_name,
    GI,
    BA_mean,
    BA_sd
  ) %>%
  filter(!is.na(GI), !is.na(BA_mean))

ba_sed_pseudocount <- min(ba_sed_df$BA_mean[ba_sed_df$BA_mean > 0], na.rm = TRUE) / 2
if (!is.finite(ba_sed_pseudocount)) {
  ba_sed_pseudocount <- 1
}

ba_sed_df <- ba_sed_df %>%
  mutate(log10_BA = log10(BA_mean + ba_sed_pseudocount))

# Fit the empty and baseline random-effects models.
mod_BA_sed_empty_random <- gam(
  log10_BA ~ 1 + s(gl_name, bs = "re"),
  data = ba_sed_df,
  method = "ML"
)

mod_BA_sed_region_baseline <- gam(
  log10_BA ~ Region + s(gl_name, bs = "re"),
  data = ba_sed_df,
  method = "ML"
)

# Fit candidate GI models; all retain the glacier random intercept.
mod_BA_sed_GI_only <- gam(
  log10_BA ~ GI + s(gl_name, bs = "re"),
  data = ba_sed_df,
  method = "ML"
)

mod_BA_sed_region_GI <- gam(
  log10_BA ~ Region + GI + s(gl_name, bs = "re"),
  data = ba_sed_df,
  method = "ML"
)

mod_BA_sed_9 <- gam(
  log10_BA ~ Region * GI + s(gl_name, bs = "re"),
  data = ba_sed_df,
  method = "ML"
)

summary(mod_BA_sed_empty_random)
summary(mod_BA_sed_region_baseline)
summary(mod_BA_sed_GI_only)
summary(mod_BA_sed_region_GI)
summary(mod_BA_sed_9)

sed_BA_model_summary <- tibble(
  model = c(
    "Null + glacier RE",
    "Region + glacier RE",
    "GI + glacier RE",
    "Region + GI + glacier RE",
    "Region x GI + glacier RE"
  ),
  AIC = c(
    AIC(mod_BA_sed_empty_random),
    AIC(mod_BA_sed_region_baseline),
    AIC(mod_BA_sed_GI_only),
    AIC(mod_BA_sed_region_GI),
    AIC(mod_BA_sed_9)
  ),
  BIC = c(
    BIC(mod_BA_sed_empty_random),
    BIC(mod_BA_sed_region_baseline),
    BIC(mod_BA_sed_GI_only),
    BIC(mod_BA_sed_region_GI),
    BIC(mod_BA_sed_9)
  )
)

print(sed_BA_model_summary)

# Refit the final additive model with REML.
mod_BA_sed_final_GI <- gam(
  log10_BA ~ Region + GI + s(gl_name, bs = "re"),
  data = ba_sed_df,
  method = "REML"
)

summary(mod_BA_sed_final_GI)

# Sensitivity analysis: omit the glacier random intercept only after the
# primary model structure has been selected.
mod_BA_sed_final_GI_no_random <- gam(
  log10_BA ~ Region + GI,
  data = ba_sed_df,
  method = "REML"
)

summary(mod_BA_sed_final_GI_no_random)

# --- Extended Data Table 1 ---

format_p_value <- function(p) {
  case_when(
    is.na(p) ~ NA_character_,
    p < 0.001 ~ "<0.001",
    TRUE ~ sprintf("%.3f", p)
  )
}

get_model_r2 <- function(model, model_summary) {
  if (inherits(model, "lm") && !inherits(model, "gam")) {
    return(unname(model_summary$adj.r.squared))
  }

  unname(model_summary$r.sq)
}

get_model_deviance <- function(model_summary) {
  if (is.null(model_summary$dev.expl)) {
    return(NA_real_)
  }

  unname(model_summary$dev.expl)
}

get_coefficient_table <- function(model_summary) {
  if (!is.null(model_summary$coefficients)) {
    return(model_summary$coefficients)
  }

  model_summary$p.table
}

get_p_column <- function(coefficient_table) {
  grep("^Pr\\(", colnames(coefficient_table), value = TRUE)[1]
}

# For the current two-region design, the Region coefficient P value represents
# the fixed regional contrast. If Region is expanded to more than two levels,
# this helper reports the smallest coefficient-level P value rather than an
# omnibus factor test.
extract_coefficient_p_value <- function(coefficient_table, term_pattern) {
  term_rows <- grep(term_pattern, rownames(coefficient_table), value = TRUE)
  if (length(term_rows) == 0) {
    return(NA_character_)
  }

  p_col <- get_p_column(coefficient_table)
  format_p_value(min(unname(coefficient_table[term_rows, p_col]), na.rm = TRUE))
}

# Random effects are represented in mgcv as smooth terms; the glacier random
# effect P value is therefore extracted from the smooth table row for s(gl_name).
extract_random_effect_p_value <- function(model_summary, smooth_name = "s(gl_name)") {
  if (is.null(model_summary$s.table) || !smooth_name %in% rownames(model_summary$s.table)) {
    return(NA_character_)
  }

  smooth_row <- model_summary$s.table[smooth_name, , drop = FALSE]
  p_col <- grep("^p", colnames(smooth_row), value = TRUE, ignore.case = TRUE)[1]
  format_p_value(unname(smooth_row[1, p_col]))
}

extract_model_table_row <- function(
  model,
  response_label,
  model_label,
  gi_term = c("none", "linear")
) {
  gi_term <- match.arg(gi_term)
  model_summary <- summary(model)
  coefficient_table <- get_coefficient_table(model_summary)

  gi_effect <- NA_character_
  gi_ci <- NA_character_
  gi_p_value <- NA_character_
  region_p_value <- extract_coefficient_p_value(coefficient_table, "^Region")
  random_effect_p_value <- extract_random_effect_p_value(model_summary)

  if (gi_term == "linear") {
    gi_row <- coefficient_table["GI", , drop = FALSE]
    p_col <- get_p_column(gi_row)
    gi_estimate <- unname(gi_row[1, "Estimate"])
    gi_se <- unname(gi_row[1, "Std. Error"])
    gi_effect <- sprintf("linear estimate = %.3f", gi_estimate)
    gi_ci <- sprintf("%.3f to %.3f", gi_estimate - 1.96 * gi_se, gi_estimate + 1.96 * gi_se)
    gi_p_value <- format_p_value(unname(gi_row[1, p_col]))
  }

  tibble(
    Response = response_label,
    Model = model_label,
    `GI effect` = gi_effect,
    `GI 95% CI` = gi_ci,
    `GI P value` = gi_p_value,
    `Region P value` = region_p_value,
    `Glacier RE P value` = random_effect_p_value,
    `Adj. R2` = get_model_r2(model, model_summary),
    `Deviance explained` = get_model_deviance(model_summary),
    AIC = AIC(model),
    BIC = BIC(model)
  )
}

extended_gam_results_table <- bind_rows(
  extract_model_table_row(
    mod_BCdis_empty_random,
    "Bray-Curtis dissimilarity to ice",
    "Null + glacier RE",
    "none"
  ),
  extract_model_table_row(
    mod_BCdis_region_baseline,
    "Bray-Curtis dissimilarity to ice",
    "Region + glacier RE",
    "none"
  ),
  extract_model_table_row(
    mod_BCdis_final_GI,
    "Bray-Curtis dissimilarity to ice",
    "Region + GI + glacier RE",
    "linear"
  ),
  extract_model_table_row(
    mod_UFdis_empty_random,
    "Weighted UniFrac dissimilarity to ice",
    "Null + glacier RE",
    "none"
  ),
  extract_model_table_row(
    mod_UFdis_region_baseline,
    "Weighted UniFrac dissimilarity to ice",
    "Region + glacier RE",
    "none"
  ),
  extract_model_table_row(
    mod_UFdis_final_GI,
    "Weighted UniFrac dissimilarity to ice",
    "Region + GI + glacier RE",
    "linear"
  ),
  extract_model_table_row(
    mod_prop_empty_random,
    "Proportion of ice source",
    "Null + glacier RE",
    "none"
  ),
  extract_model_table_row(
    mod_prop_region_baseline,
    "Proportion of ice source",
    "Region + glacier RE",
    "none"
  ),
  extract_model_table_row(
    mod_prop_final_GI,
    "Proportion of ice source",
    "Region + GI + glacier RE",
    "linear"
  ),
  extract_model_table_row(
    mod_BA_sed_empty_random,
    "Bacterial abundance",
    "Null + glacier RE",
    "none"
  ),
  extract_model_table_row(
    mod_BA_sed_region_baseline,
    "Bacterial abundance",
    "Region + glacier RE",
    "none"
  ),
  extract_model_table_row(
    mod_BA_sed_final_GI,
    "Bacterial abundance",
    "Region + GI + glacier RE",
    "linear"
  )
) %>%
  mutate(
    `Adj. R2` = round(`Adj. R2`, 3),
    `Deviance explained` = round(`Deviance explained`, 3),
    AIC = round(AIC, 1),
    BIC = round(BIC, 1)
  )

# Add Region x GI sensitivity models to the Extended Data Table 1 export.
# The GI coefficient is the slope in the reference region; the interaction
# coefficient is the between-region slope difference.
extract_interaction_table_row <- function(model, response_label, model_label) {
  model_summary <- summary(model)
  coefficient_table <- get_coefficient_table(model_summary)
  p_col <- get_p_column(coefficient_table)
  interaction_rows <- grep("^Region.*:GI$|^GI:Region", rownames(coefficient_table), value = TRUE)
  if (length(interaction_rows) != 1) {
    stop("Could not uniquely identify the Region x GI coefficient.")
  }

  interaction_row <- coefficient_table[interaction_rows, , drop = FALSE]
  interaction_estimate <- unname(interaction_row[1, "Estimate"])
  interaction_se <- unname(interaction_row[1, "Std. Error"])
  region_rows <- grep("^Region[^:]*$", rownames(coefficient_table), value = TRUE)
  region_p_value <- if (length(region_rows) == 0) {
    NA_character_
  } else {
    format_p_value(min(unname(coefficient_table[region_rows, p_col]), na.rm = TRUE))
  }

  extract_model_table_row(model, response_label, model_label, "linear") %>%
    mutate(
      `GI effect` = sub("linear estimate", "reference-region slope", `GI effect`, fixed = TRUE),
      `Region P value` = region_p_value,
      `Region x GI effect` = sprintf("slope difference = %.3f", interaction_estimate),
      `Region x GI 95% CI` = sprintf(
        "%.3f to %.3f",
        interaction_estimate - 1.96 * interaction_se,
        interaction_estimate + 1.96 * interaction_se
      ),
      `Region x GI P value` = format_p_value(unname(interaction_row[1, p_col]))
    )
}

extended_gam_results_table <- extended_gam_results_table %>%
  mutate(
    `Region x GI effect` = NA_character_,
    `Region x GI 95% CI` = NA_character_,
    `Region x GI P value` = NA_character_
  ) %>%
  bind_rows(
    extract_interaction_table_row(
      mod_BCdis_9,
      "Bray-Curtis dissimilarity to ice",
      "Region x GI + glacier RE"
    ),
    extract_interaction_table_row(
      mod_UFdis_9,
      "Weighted UniFrac dissimilarity to ice",
      "Region x GI + glacier RE"
    ),
    extract_interaction_table_row(
      mod_prop_9,
      "Proportion of ice source",
      "Region x GI + glacier RE"
    ),
    extract_interaction_table_row(mod_BA_sed_9, "Bacterial abundance", "Region x GI + glacier RE")
  )

tryCatch(
  write.csv(
    extended_gam_results_table,
    file.path(result_dir, "extended_data_table_1_gamm_models.csv"),
    row.names = FALSE
  ),
  error = function(e) {
    warning(
      "Could not write extended_GAM_model_summary_table.csv. ",
      "Close the file if it is open in another program and rerun this chunk. ",
      "The table is still printed below. Original error: ",
      conditionMessage(e)
    )
  }
)

# Extended sediment GAMM model summary. Beta GAMM estimates are on the
# logit-link scale; bacterial abundance is modeled as log10-transformed BA.
# Glacier RE denotes a glacier random-effect smooth; Glacier RE P values are
# approximate smooth-term P values from mgcv.
print(extended_gam_results_table)

# --- Supplementary PERMANOVA across habitat-GI categories ---

# These tests provide a supplementary group-level comparison across composite
# habitat-GI categories, including ice, water GI categories, and sediment GI
# categories. Because `Group_GI` combines habitat type and categorical GI
# classes, this analysis is intended only to test broad compositional
# differences among these predefined groups. It is not used to infer continuous
# GI effects or to replace the sediment GAMM analyses above.

ps_rel_initial <- transform_sample_counts(ps, function(x) x / sum(x))

set.seed(666)
bray_dist_initial <- distance(ps_rel_initial, method = "bray")

env_vars <- c("Group_GI")
env_data <- sample_metadata_df[, env_vars, drop = FALSE]
env_data <- data.frame(env_data, row.names = rownames(env_data))
env_data <- env_data[complete.cases(env_data), , drop = FALSE]
env_data$Group_GI <- factor(env_data$Group_GI)

bray_dist_initial_ordered <- as.dist(
  as.matrix(bray_dist_initial)[rownames(env_data), rownames(env_data)]
)

# --- Overall PERMANOVA and dispersion check ---

set.seed(666)

adonis_group_gi <- adonis2(
  bray_dist_initial_ordered ~ Group_GI,
  data = env_data,
  permutations = 999
)

betadisper_group_gi <- betadisper(
  bray_dist_initial_ordered,
  env_data$Group_GI
)
betadisper_group_gi_anova <- anova(betadisper_group_gi)

adonis_group_gi
betadisper_group_gi_anova

# --- Pairwise PERMANOVA ---

pairwise_adonis_custom <- function(dist_matrix, meta_df, group_col) {
  groups <- factor(meta_df[[group_col]])
  u_groups <- levels(groups)
  comb <- combn(u_groups, 2)

  result <- data.frame(
    comparison = character(),
    F_value = numeric(),
    R2 = numeric(),
    p_value = numeric(),
    adj_p_value = numeric(),
    stringsAsFactors = FALSE
  )

  for (i in 1:ncol(comb)) {
    group_pair <- comb[, i]
    subset_idx <- groups %in% group_pair
    subset_samples <- rownames(meta_df)[subset_idx]

    dist_mat <- as.matrix(dist_matrix)
    subset_dist <- as.dist(dist_mat[subset_samples, subset_samples])
    subset_groups <- droplevels(groups[subset_idx])
    temp_df <- data.frame(group = subset_groups, row.names = subset_samples)

    perm_res <- adonis2(subset_dist ~ group, data = temp_df, permutations = 999)

    result[i, "comparison"] <- paste(group_pair[1], "vs", group_pair[2])
    result[i, "F_value"] <- perm_res[1, "F"]
    result[i, "R2"] <- perm_res[1, "R2"]
    result[i, "p_value"] <- perm_res[1, "Pr(>F)"]
  }

  result$adj_p_value <- p.adjust(result$p_value, method = "BH")
  result
}

pairwise_results <- pairwise_adonis_custom(bray_dist_initial_ordered, env_data, "Group_GI")
print(pairwise_results)

# --- Extended Data Table 2 ---

overall_permanova_row <- tibble(
  Analysis = "Overall PERMANOVA",
  Comparison = "All habitat-GI categories",
  Df = unname(adonis_group_gi[1, "Df"]),
  F = unname(adonis_group_gi[1, "F"]),
  R2 = unname(adonis_group_gi[1, "R2"]),
  `P value` = format_p_value(unname(adonis_group_gi[1, "Pr(>F)"])),
  `Adjusted P value` = NA_character_,
  `Dispersion P value` = format_p_value(unname(betadisper_group_gi_anova[1, "Pr(>F)"]))
)

pairwise_permanova_rows <- pairwise_results %>%
  transmute(
    Analysis = "Pairwise PERMANOVA",
    Comparison = comparison,
    Df = 1,
    F = F_value,
    R2 = R2,
    `P value` = format_p_value(p_value),
    `Adjusted P value` = format_p_value(adj_p_value),
    `Dispersion P value` = NA_character_
  )

extended_permanova_results_table <- bind_rows(
  overall_permanova_row,
  pairwise_permanova_rows
) %>%
  mutate(
    F = round(F, 3),
    R2 = round(R2, 3)
  )

tryCatch(
  write.csv(
    extended_permanova_results_table,
    file.path(result_dir, "extended_data_table_2_permanova.csv"),
    row.names = FALSE
  ),
  error = function(e) {
    warning("Could not write extended PERMANOVA table: ", conditionMessage(e))
  }
)

# Extended Data Table 2 tests differences in microbial community composition
# among habitat and glacier-influence categories using Bray-Curtis
# dissimilarity.
print(extended_permanova_results_table)
