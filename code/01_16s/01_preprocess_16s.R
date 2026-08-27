# 16S rRNA sequence preprocessing
# Manuscript output: Extended Data Figure 6.
# Produces the filtered phyloseq object used by downstream 16S analyses.

gc()
graphics.off()

library(phyloseq)
library(qiime2R)
library(tidyverse)
library(decontam)
library(biomformat)

# --- Import and process 16S data ---

input_dir <- file.path("data", "processed", "16s")
result_dir <- file.path("results", "01_16s")
intermediate_dir <- file.path(result_dir, "intermediate")
qc_dir <- file.path(result_dir, "qc")

dir.create(intermediate_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

ASV16SFull <- read.csv(
  file.path(input_dir, "amplicons_all_16S_table_tsv.csv"),
  sep = ","
)
ASV16S <- ASV16SFull %>%
  select(-taxonomy) %>%
  column_to_rownames(var = "OTUID")

taxonomy16SFull <- ASV16SFull %>%
  select(OTUID, taxonomy) %>%
  separate(
    taxonomy,
    c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
    "; "
  ) %>%
  column_to_rownames(var = "OTUID") %>%
  as.matrix()
taxonomy16SFull <- gsub(".__", "", taxonomy16SFull)

# Load sample metadata.
metadata16S_raw <- read.csv(file.path(input_dir, "metadata_16S.csv"), sep = ",") %>%
  column_to_rownames(var = "Sample")

# Assign sediment GI groups while retaining habitat labels for other samples.

# Calculate sediment GI tertiles.
quantile_breaks <- quantile(
  metadata16S_raw$GI[metadata16S_raw$Source == "sed" & !is.na(metadata16S_raw$GI)],
  probs = seq(0, 1, by = 1 / 3),
  na.rm = TRUE
)

print(quantile_breaks)

# Assign Group_GI from sediment GI tertiles and existing habitat labels.
metadata16S <- metadata16S_raw %>%
  mutate(
    Group_GI = case_when(
      Source == "ice" ~ "ice",
      Source %in% c("sed", "water") & !is.na(GI) ~ paste(
        as.character(cut(
          GI,
          breaks = quantile_breaks,
          include.lowest = TRUE,
          labels = c("Low_Glacial", "Mid_Glacial", "High_Glacial")
        )),
        Source,
        sep = "_"
      ),
      !is.na(GI) ~ as.character(cut(
        GI,
        breaks = quantile_breaks,
        include.lowest = TRUE,
        labels = c("Low_Glacial", "Mid_Glacial", "High_Glacial")
      )),
      TRUE ~ NA_character_
    )
  )

stats_quant <- metadata16S %>%
  group_by(Group_GI) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(Method = "Quantile")

print(stats_quant)

# Create the initial phyloseq object.

ASV_table_16S <- otu_table(data.frame(ASV16S), taxa_are_rows = TRUE)
sample_data_16S <- sample_data(data.frame(metadata16S))
tax_table_16SFull <- tax_table(taxonomy16SFull)

PAMIR_16S <- merge_phyloseq(ASV_table_16S, sample_data_16S, tax_table_16SFull)
PAMIR_16S <- prune_taxa(taxa_sums(PAMIR_16S) > 0, PAMIR_16S)
print(PAMIR_16S)

# Replace missing taxonomy labels with "unknown".

tax_mat <- as.matrix(tax_table(PAMIR_16S))
tax_mat[is.na(tax_mat)] <- "unknown"
tax_table(PAMIR_16S) <- tax_table(tax_mat)

# --- Taxonomic filtering ---

# Retain only ASVs assigned to Bacteria.
PAMIR_16S_filtered <- subset_taxa(
  PAMIR_16S,
  !is.na(Kingdom) & Kingdom == "Bacteria"
)
print(paste("After keeping only Bacteria:", ntaxa(PAMIR_16S_filtered), "taxa remaining"))

# Remove chloroplast ASVs.
PAMIR_16S_filtered <- subset_taxa(
  PAMIR_16S_filtered,
  is.na(Order) | Order != "Chloroplast"
)
print(paste("After chloroplast filtering:", ntaxa(PAMIR_16S_filtered), "taxa remaining"))

# Remove mitochondrial ASVs.
PAMIR_16S_filtered <- subset_taxa(
  PAMIR_16S_filtered,
  is.na(Family) | Family != "Mitochondria"
)
print(paste("After mitochondria filtering:", ntaxa(PAMIR_16S_filtered), "taxa remaining"))

# --- Control processing and contaminant removal ---

# Separate positive controls from all remaining samples.

PAMIR_16S_positive_controls <- subset_samples(
  PAMIR_16S_filtered,
  Source == "Positive_Control"
)
PAMIR_16S_bacteria <- subset_samples(PAMIR_16S_filtered, Source != "Positive_Control")

# Remove empty ASVs from both objects.
PAMIR_16S_bacteria <- prune_taxa(taxa_sums(PAMIR_16S_bacteria) > 0, PAMIR_16S_bacteria)
PAMIR_16S_positive_controls <- prune_taxa(
  taxa_sums(PAMIR_16S_positive_controls) > 0,
  PAMIR_16S_positive_controls
)

# Identify negative controls.
PAMIR_16S_negative_controls <- subset_samples(
  PAMIR_16S_bacteria,
  Source == "Negative_Control"
)
PAMIR_16S_negative_controls <- prune_taxa(
  taxa_sums(PAMIR_16S_negative_controls) > 0,
  PAMIR_16S_negative_controls
)
print(PAMIR_16S_negative_controls)

sample_data(PAMIR_16S_bacteria)$is.ng <-
  sample_data(PAMIR_16S_bacteria)$Source == "Negative_Control"

# Identify contaminants with decontam.
contamdf.prev <- isContaminant(
  PAMIR_16S_bacteria,
  method = "prevalence",
  neg = "is.ng",
  threshold = 0.5
)
print(paste("Contaminants identified:", sum(contamdf.prev$contaminant)))

# Calculate the proportion of reads assigned to contaminants.
contam_reads <- sum(taxa_sums(PAMIR_16S_bacteria)[contamdf.prev$contaminant])
total_reads <- sum(taxa_sums(PAMIR_16S_bacteria))
print(contam_reads / total_reads)

# Remove contaminants and negative controls.
PAMIR_16S_bacteria_clean <- prune_taxa(!contamdf.prev$contaminant, PAMIR_16S_bacteria)
PAMIR_16S_bacteria_clean <- subset_samples(PAMIR_16S_bacteria_clean, !is.ng)
PAMIR_16S_bacteria_clean <- prune_taxa(
  taxa_sums(PAMIR_16S_bacteria_clean) > 0,
  PAMIR_16S_bacteria_clean
)
print(
  paste(
    "After decontamination and removing negative controls:",
    ntaxa(PAMIR_16S_bacteria_clean),
    "taxa remaining"
  )
)

# Retain ASVs with more than one read in at least two samples.
PAMIR_16S_bacteria_clean_filt <- PAMIR_16S_bacteria_clean %>%
  filter_taxa(function(x) sum(x > 1) > 1, TRUE)

print(
  sprintf(
    "Removed: %d rare ASVs (%.2f%%)",
    ntaxa(PAMIR_16S_bacteria_clean) - ntaxa(PAMIR_16S_bacteria_clean_filt),
    100 * (1 - ntaxa(PAMIR_16S_bacteria_clean_filt) / ntaxa(PAMIR_16S_bacteria_clean))
  )
)

# Retain ASVs with at least 10 reads across all samples.
PAMIR_16S_bacteria_final <- filter_taxa(
  PAMIR_16S_bacteria_clean_filt,
  function(x) sum(x) >= 10,
  TRUE
)

print(
  sprintf(
    "Removed: %d rare ASVs",
    ntaxa(PAMIR_16S_bacteria_clean_filt) - ntaxa(PAMIR_16S_bacteria_final)
  )
)

print(PAMIR_16S_bacteria_final)

# Define the sample-type order.
Source_order <- c("ice", "water", "sed")

sample_data(PAMIR_16S_bacteria_final)$Source <- factor(
  sample_data(PAMIR_16S_bacteria_final)$Source,
  levels = Source_order
)

# --- Save cleaned 16S objects ---

cleaned_path <- file.path(intermediate_dir, "PAMIR_16S_bacteria_cleaned.rds")
positive_control_path <- file.path(intermediate_dir, "PAMIR_16S_positive_controls.rds")
final_path <- file.path(intermediate_dir, "PAMIR_16S_final.rds")

saveRDS(PAMIR_16S_bacteria_final, file = cleaned_path)
saveRDS(PAMIR_16S_positive_controls, file = positive_control_path)

asv_ids <- taxa_names(PAMIR_16S_bacteria_final)
writeLines(asv_ids, file.path(intermediate_dir, "ASV_list.txt"))

# --- Positive-control validation ---

gc()
graphics.off()

library(ggplot2)
library(viridis)

# Load the positive-control phyloseq object.
PAMIR_16S_positive_controls <- readRDS(positive_control_path)
print(PAMIR_16S_positive_controls)

# Export the combined OTU and taxonomy table.
otu_table <- as.data.frame(otu_table(PAMIR_16S_positive_controls)) %>%
  rownames_to_column("OTUID")

tax_table <- as.data.frame(tax_table(PAMIR_16S_positive_controls)) %>%
  rownames_to_column("OTUID") %>%
  mutate(across(where(is.factor), as.character))

combined_table <- left_join(otu_table, tax_table, by = "OTUID")
write_csv(combined_table, file.path(qc_dir, "positive_control_otu_tax_table.csv"))

# Convert counts to relative abundance.
PAMIR_16S_positive_controls_rel <- transform_sample_counts(
  PAMIR_16S_positive_controls,
  function(x) x / sum(x, na.rm = TRUE)
)

# Aggregate taxa at genus level.
PAMIR_16S_positive_controls_genus <- tax_glom(
  PAMIR_16S_positive_controls_rel,
  taxrank = "Genus"
)

# Convert the object to long format.
pos_ctrl_genus_df <- psmelt(PAMIR_16S_positive_controls_genus)

# Load the theoretical mock-community composition.
bc_data <- read_csv(
  file.path(input_dir, "mock_theoretical_data.csv"),
  show_col_types = FALSE
)

# Format the theoretical composition and combine it with observed controls.
bc_data_formatted <- bc_data %>%
  rename(Abundance = Theory) %>%
  mutate(Sample = "Theory") %>%
  select(Sample, Abundance, Genus)

combined_genus_df <- bind_rows(pos_ctrl_genus_df, bc_data_formatted) %>%
  mutate(
    Sample = factor(
      Sample,
      levels = c(
        "Theory", "posCON_151_S1", "posCON_152_S2", "posCON_94_S44",
        "PAMIR__45_S333", "PAMIR__72_S360"
      )
    )
  )

# Identify the expected mock-community genera.
top8_genera <- bc_data$Genus
top8_genera <- unique(top8_genera)

# Recode all other genera as "Others".
combined_genus_df_new <- combined_genus_df %>%
  mutate(Genus = ifelse(Genus %in% top8_genera, as.character(Genus), "Others"))

# Generate the genus-level stacked bar plot.
genus_plot <- ggplot(combined_genus_df_new, aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  scale_fill_brewer(palette = "Set3") +
  labs(
    title = "Genus-level Composition of Positive Controls",
    x = "Sample",
    y = "Relative Abundance",
    fill = "Genus"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave(
  file.path(qc_dir, "extended_data_figure_6_mock_composition.pdf"),
  genus_plot,
  width = 6,
  height = 6
)

# --- Merge ice sample replicates ---

# Combine replicate ice samples from the Trift and Tortin glaciers.

# Reload the cleaned study object before creating the merged analysis object.
ps <- readRDS(cleaned_path)

# Replace phyloseq sample names with Sample_ID values.
new_sample_names <- as.character(sample_data(ps)$Sample_ID)
sample_names(ps) <- new_sample_names

# Initialize the merge groups with the original sample names.
sample_data(ps)$Merge_Group <- sample_names(ps)

# Define the replicate sample sets.
s20_replicates <- c("S20_A", "S20_B", "S20_C", "S20_D", "S20_E")
s21_replicates <- c("S21_A", "S21_B", "S21_C", "S21_D", "S21_E")

# Assign only the designated replicates to shared merge groups.
sample_data(ps)$Merge_Group[sample_names(ps) %in% s20_replicates] <- "S20"
sample_data(ps)$Merge_Group[sample_names(ps) %in% s21_replicates] <- "S21"

# Merge samples and sum their ASV counts.
PAMIR_16S_merged_ice <- merge_samples(ps, "Merge_Group")

# Restore taxa-as-rows orientation if merge_samples changes it.
if (!taxa_are_rows(PAMIR_16S_merged_ice) && taxa_are_rows(ps)) {
  otu_table(PAMIR_16S_merged_ice) <- otu_table(
    t(otu_table(PAMIR_16S_merged_ice)),
    taxa_are_rows = TRUE
  )
}

# Restore sample metadata.
original_meta <- data.frame(sample_data(ps))

# Keep one representative metadata row per Merge_Group.
restored_meta <- original_meta[!duplicated(original_meta$Merge_Group), ]

# Set metadata row names to Merge_Group values.
rownames(restored_meta) <- restored_meta$Merge_Group

# Update Sample_ID to match the merged row names.
restored_meta$Sample_ID <- rownames(restored_meta)

# Remove the auxiliary merge column.
restored_meta$Merge_Group <- NULL

# Attach the restored metadata to the merged phyloseq object.
sample_data(PAMIR_16S_merged_ice) <- sample_data(restored_meta)

# Inspect and save the final merged object.
print(PAMIR_16S_merged_ice)
saveRDS(PAMIR_16S_merged_ice, file = final_path)

# Analysis complete.
