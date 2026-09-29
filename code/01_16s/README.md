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

Step 1 generates Extended Data Figure 6 (mock validation); step 2 generates
Extended Data Figure 1 (Alps–Central Asia ASV overlap within each habitat).
Step 5 generates Figure 2a and its PERMANOVA/PERMDISP tables.

## Inputs

- `data/processed/16s/amplicons_all_16S_table_tsv.csv`
- `data/processed/16s/metadata_16S.csv`
- `data/processed/16s/mock_theoretical_data.csv`
- `data/processed/16s/dna-sequences.tree`

## Fixed settings

- The 95%-coverage `estimateD()` bootstrap and plotted jitter use seed `666`.
- Step 5 calculates Bray-Curtis dissimilarity for merged ice and sediment
  samples, excluding water. It tests the four `Group_GI` classes with 9,999
  permutations restricted within `gl_name`, then applies the same permutation
  design to PERMDISP. Both tests use seed `666`; the PERMANOVA is an overall
  group test and does not establish pairwise or strictly monotonic differences.
