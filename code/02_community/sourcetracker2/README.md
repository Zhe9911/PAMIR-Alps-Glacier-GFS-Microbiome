# SourceTracker2 runner

[Community workflow](../README.md) |
[Code index](../../README.md) |
[Setup](../../../README.md#quick-start)

This optional runner recalculates the Figure 2e SourceTracker2 result. The
default workflow reads
`data/processed/16s/sourcetracker/all_results_depth_10000.csv`.

## Run

From the repository root:

```text
Rscript code/02_community/01_gamm_sourcetracker.R
python code/02_community/sourcetracker2/run_sourcetracker2_script.py
```

The R command writes `ST_otu_table.biom` and `ST_metadata.csv` under
`results/generated/02_community/sourcetracker/`; the Python command writes
`Results/all_results.csv` below the same directory.

## Fixed settings

- 19 glacier-wise runs; ice sources and non-ice sinks; depth 10,000.
- Burn-in 100; 10 restarts; 50 draws/restart; 8 workers.
- Alpha1 0.001; alpha2 0.1; beta 10.

Reference environment: Python 3.8.18; SourceTracker2 `2.0.1.dev0`, commit
`a921a0563391203581f48e56f750391898f4a625`; biom-format 2.1.16; NumPy 1.24.4;
pandas 2.0.3; tqdm 4.67.1.

SourceTracker2 is not managed by `renv` and exposes no seed, so exact
regeneration is not guaranteed. Generated results do not replace the frozen
input.
