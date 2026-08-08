# Genome-resolved workflow

[Code index](../README.md) |
[Setup](../../README.md#quick-start) |
[Workflow manifest](../../manifest/workflow-dependencies.tsv)

## Run

From the repository root:

```text
Rscript code/03_mags/01_mag_asv_validation.R
Rscript code/03_mags/02_mag_gi_association.R
Rscript code/03_mags/03_mag_response_groups.R
Rscript code/03_mags/04_mag_phylogeny_plot.R
Rscript code/03_mags/05_mag_phylogenetic_signal.R
```

Steps 1 and 2 create `PAMIR_MAGs_rela.rds` and `MAGs_info.csv`; steps 3–5
require both. Outputs are written to `results/03_mags/`.

## Inputs

- `data/processed/metagenomics/coverm_mags_relative_abundance.tsv`
- `data/processed/metagenomics/coverm_mags_count.tsv`
- `data/processed/metagenomics/metadata_metagenomics.csv`
- `data/processed/metagenomics/gtdbtk.bac120.summary.tsv`
- `data/processed/metagenomics/quality_report.tsv`
- `data/processed/metagenomics/iqtree_tree.treefile`
- `results/01_16s/intermediate/PAMIR_16S_final.rds`

## Main settings

- ASV–MAG validation: 62 matched samples, Spearman/BH, and 9,999 Mantel
  permutations.
- Community concordance: Bray–Curtis NMDS, symmetric Procrustes, and PROTEST
  with 9,999 permutations.
- ANCOM-BC2 uses eight workers, sediment MAG counts, `GI + Region` fixed
  effects, glacier random intercepts, 10% prevalence, library cutoff 1,000,
  and BH-FDR 0.05.
- Response groups follow ice, high-, mid-, and low-GI; trees use a
  Patescibacteriota root.
- Signal tests use MAG-specific `se_GI`; Blomberg's K uses 999 permutations.
