# Functional workflow

[Code index](../README.md) |
[Setup](../../README.md#quick-start) |
[Workflow manifest](../../manifest/workflow-dependencies.tsv)

## Run

After preparing the three KEGG files below, run from the repository root:

```text
Rscript code/04_functions/01_ko_gsea.R
Rscript code/04_functions/02_kegg_phylogenetic_projection.R
```

Step 1 writes the complete GSEA tables and the source volcano for Extended Data
Figure 5, whose manuscript layout was manually fine-tuned. Step 2 uses step 1's
significant GSEA tables. Outputs are written to `results/04_functions/`.

## Inputs

The companion data package supplies:

- `data/processed/metagenomics/ko_mag_dt.tsv`
- `data/processed/metagenomics/iqtree_tree.treefile`

Module 03 supplies:

- `results/03_mags/intermediate/MAGs_info.csv`
- `results/03_mags/intermediate/PAMIR_MAGs_rela.rds`

### External KEGG inputs

These required files are not redistributed:

- `data/processed/metagenomics/kegg_brite.csv`: KEGG `ko00001` BRITE hierarchy,
  downloaded in htext format on 2025-06-24 and flattened to
  `A,B,C,KO,Description,PATH,BR` (62,998 rows in the manuscript snapshot).
- `data/processed/metagenomics/kegg_module_TERM2GENE.tsv`:
  [`link/ko/module`](https://rest.kegg.jp/link/ko/module), columns `term,KO`
  (4,156 rows in the manuscript snapshot).
- `data/processed/metagenomics/kegg_module_TERM2NAME.tsv`:
  [`list/module`](https://rest.kegg.jp/list/module), columns `term,name`
  (570 rows in the manuscript snapshot).

Both module mappings were queried on 2026-04-02 using `KEGGREST` 1.50.0, with
prefixes removed and rows deduplicated and sorted. This client version does not
identify a KEGG database release. Live content may differ; obtain files under the
[KEGG terms](https://www.kegg.jp/kegg/legal.html).

## Main settings

- KOs occur in at least five valid MAGs and are absent from at least one.
- Bias-reduced binomial models use LFC-GI, completeness, and log10 genome size;
  KO-level Wald P values use BH correction.
- Finite signed Wald z statistics enter GSEA; gene sets contain 5–500 KOs,
  with `clusterProfiler`/`fgsea`, seed `42`, and BH correction.
- Feature scores are detected fractions of core-enrichment KOs.
- Within-phylum quasi-binomial models adjust for completeness and genome size;
  BH correction is applied within each feature.
