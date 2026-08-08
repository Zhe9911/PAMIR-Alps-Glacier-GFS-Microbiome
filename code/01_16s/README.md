# 16S workflow

[Code index](../README.md) |
[Setup](../../README.md#quick-start) |
[Workflow manifest](../../manifest/workflow-dependencies.tsv)

## Run

From the repository root:

```text
Rscript code/01_16s/01_preprocess_16s.R
Rscript code/01_16s/02_asv_overlap.R
Rscript code/01_16s/03_alpha_diversity.R
Rscript code/01_16s/04_dissimilarity.R
Rscript code/01_16s/05_nmds.R
```

Outputs are written to `results/01_16s/`. The main downstream object is
`results/01_16s/intermediate/PAMIR_16S_final.rds`.

## Inputs

- `data/processed/16s/amplicons_all_16S_table_tsv.csv`
- `data/processed/16s/metadata_16S.csv`
- `data/processed/16s/mock_theoretical_data.csv`
- `data/processed/16s/dna-sequences.tree`

The 95%-coverage `estimateD()` bootstrap and plotted jitter use seed `666`.
