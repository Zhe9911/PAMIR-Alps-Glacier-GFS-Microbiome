# Computational environment

The numbered R workflows use R 4.5.2, Bioconductor 3.22 and renv 1.1.5.

## Restore

From the repository root in a clean checkout with R 4.5.2:

```text
Rscript -e "renv::restore(prompt = FALSE)"
```

Package versions are locked in [`renv.lock`](../renv.lock).

## System requirements

The reference environment uses Windows with Rtools45. Figure generation
requires Cairo:

```text
Rscript -e "stopifnot(capabilities('cairo'))"
```

Linux and macOS source installations may need the build tools and native
libraries in [`r-system-requirements.tsv`](r-system-requirements.tsv). The table
reproduces upstream package metadata, including optional and vignette-only
requirements.

## Optional workflows

Full iCAMP recalculation requires Linux/HPC; the default uses the frozen result
described in the [`02_community` guide](../code/02_community/README.md).

SourceTracker2 uses a separate Python environment; see the
[optional audit guide](../code/02_community/sourcetracker2/README.md).
