# PAMIR-Alps-Glacier-GFS-Microbiome

Analysis code and reproducibility materials for the manuscript:

> *Glacier-to-downstream continuum shapes microbial connectivity and adaptive genomic traits*

## Quick start

1. Get the code and enter the repository root:

   ```text
   git clone https://github.com/Zhe9911/PAMIR-Alps-Glacier-GFS-Microbiome.git
   cd PAMIR-Alps-Glacier-GFS-Microbiome
   ```

2. Download `PAMIR-Alps-Glacier-GFS-Microbiome-data.zip` from
   [Zenodo](https://doi.org/10.5281/zenodo.21703331). Extract the archive, then
   copy its `data/` directory into the repository root.

   The archive includes the complete 183-sink SourceTracker input required by
   the [`02_community` workflow](code/02_community/README.md#inputs).

3. Install R 4.5.2 and restore the locked R packages:

   ```text
   Rscript -e "renv::restore(prompt = FALSE)"
   ```

4. Run the modules in order:

   | Order | Module | Requirements |
   |---|---|---|
   | 1 | [`01_16s`](code/01_16s/README.md) | Zenodo data |
   | 2 | [`02_community`](code/02_community/README.md) | Module 01 outputs and Zenodo data |
   | 3 | [`03_mags`](code/03_mags/README.md) | Module 01 outputs and Zenodo data |
   | 4 | [`04_functions`](code/04_functions/README.md) | Module 03 outputs, Zenodo data and three external KEGG files |
   | 5 | [`05_polaromonas`](code/05_polaromonas/README.md) | Module 03 outputs and Zenodo data |

## Repository layout

| Path | Contents |
|---|---|
| [`code/`](code/README.md) | Analysis modules |
| [`config/`](config/README.md) | Analysis parameters |
| [`environment/`](environment/README.md) | R and system environment |
| [`manifest/`](manifest/README.md) | Workflow dependencies and manuscript output map |
| [`results/`](results/README.md) | Git-ignored runtime outputs |

## License and citation

Code is released under the [MIT License](LICENSE). The
[Zenodo data package](https://doi.org/10.5281/zenodo.21703331) is released under
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Citation metadata are
provided in [`CITATION.cff`](CITATION.cff).
