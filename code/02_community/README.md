# Community turnover and assembly workflow

[Code index](../README.md) |
[Setup](../../README.md#quick-start) |
[Workflow manifest](../../manifest/workflow-dependencies.tsv)

Complete the [`01_16s` workflow](../01_16s/README.md) first.

| Route | Script | Manuscript output |
|---|---|---|
| Default | `01_gamm_sourcetracker.R` | Figures 2c-e; sediment GAMMs in Extended Data Table 1; water GAMMs in Extended Data Table 2 |
| Optional Linux/HPC | `02_icamp_run.R` | Full upstream iCAMP calculation for Figure 2f |
| Default | `03_icamp_analysis.R` | Figure 2f and Extended Data Table 3 |

## Default run

This route reads frozen SourceTracker2 and iCAMP results and skips
`02_icamp_run.R`:

```text
Rscript code/02_community/01_gamm_sourcetracker.R
Rscript code/02_community/03_icamp_analysis.R
```

`03_icamp_analysis.R` uses
`results/02_community/icamp/sediment_pairwise_process_fractions.csv` when
present; otherwise it reads
`data/derived/frozen/16s/icamp/sediment_pairwise_process_fractions.csv`.

## Optional recalculations

The optional [SourceTracker2 audit](sourcetracker2/README.md) does not replace
the frozen input used by `01_gamm_sourcetracker.R`.

For a full Linux/HPC iCAMP calculation:

```text
ICAMP_RAND=1000 ICAMP_NWORKER=48 ICAMP_MEMORY_G=300 \
  Rscript code/02_community/02_icamp_run.R
Rscript code/02_community/03_icamp_analysis.R
```

Outputs are written to `results/02_community/`; large intermediates are written
to `results/generated/`.

## Inputs

- `results/01_16s/intermediate/PAMIR_16S_final.rds`
- `results/01_16s/dissimilarity/within_glacier_ice_sed_bray_curtis.csv`
- `results/01_16s/dissimilarity/within_glacier_ice_sed_weighted_unifrac.csv`
- `data/processed/16s/env_GI.csv`
- `data/processed/16s/sourcetracker/all_results_depth_10000.csv`
- `data/derived/frozen/16s/icamp/sediment_pairwise_process_fractions.csv`
- `data/processed/16s/dna-sequences.tree` (water weighted UniFrac and optional full iCAMP run)

The frozen SourceTracker table must contain two rows (`ice` and `Unknown`) for
each of 146 sediment and 37 water sinks (183 sinks total). The workflow validates
this coverage and stops if any required sink is absent.

## Fixed settings

- Unrarefied sediment counts; ASV total >=10; prevalence >=3 samples or 2%.
- Water and sediment responses use the same GAMM framework: additive Region +
  GI fixed effects, glacier identity as a random-effect smooth, ML candidate
  comparisons, and REML final models. Beta responses use a logit-link beta
  family; bacterial abundance uses `log10(BA + half the smallest positive BA)`.
- The secondary water environmental screen also compares Region-adjusted smooth
  terms for distance from the glacier snout, glacier size, and both predictors
  together. These exploratory candidates are printed for review but are not
  included in Extended Data Table 2.
- For each water beta response, the marginal GI contribution is reported as the
  difference in adjusted R2 and deviance explained between the REML additive
  `Region + GI + glacier RE` model and its REML `Region + glacier RE` baseline.
- iCAMP: 1,000 randomizations, seed `20260416`; HoS, HeS, HD, DL and DR.
- Freedman-Lane MRQAP for pairwise mean GI and absolute GI difference: 9,999
  glacier-constrained permutations, seed `20260715`, predictor-wise BH correction.
