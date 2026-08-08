# *Polaromonas* workflow

[Code index](../README.md) |
[Setup](../../README.md#quick-start) |
[Workflow manifest](../../manifest/workflow-dependencies.tsv)

## Run

From the repository root:

```text
Rscript code/05_polaromonas/01_pangenome_and_tree.R
Rscript code/05_polaromonas/02_clade_contrast.R
Rscript code/05_polaromonas/03_ko_marker_analysis.R
```

Step 1 creates `Polaromonas_KO_copy_number_dt.tsv` and
`Polaromonas_gene_cluster_prevalence_classification.tsv`, both required by step
3. Outputs are written to `results/05_polaromonas/`.

## Inputs

- `data/processed/polaromonas/Polaromonas_tree.treefile`
- `data/processed/polaromonas/Polaromonas_tree_with_outgroup.treefile`
- `data/processed/polaromonas/Polaromonas_Pan_gene_clusters_summary.txt`
- `data/processed/polaromonas/MAGs_Group_clade.tsv`
- `data/processed/polaromonas/ANIb_percentage_identity.tsv`
- `results/03_mags/intermediate/MAGs_info.csv`
- `config/polaromonas-marker-genes.tsv`

## Main settings

- Gene clusters present in ≥75% of MAGs are `Soft_Core`, those present in ≤15%
  are `Cloud`, and the remainder are `Shell`.
- Figure 5a is a rectangular tree rooted with two *Rhodoferax* outgroups,
  which are omitted from the displayed 46-MAG tree.
- `clade1` and `clade2` correspond to Clade-D (n = 6) and Clade-G (n = 10).
- The primary two-sided clade test evaluates all 8,008 fixed-size partitions.
- Sensitivity analyses use seeds `123` and `124`; plot jitter uses `125` and
  `126`.
- KO prevalence uses two-sided Fisher exact tests with BH correction across
  2,759 KOs.
