# RICE_TDP43-prion-seeding

This repository contains the analysis code associated with the manuscript “ALS-associated mutations in the TDP-43 C-terminal domain alter in vitro aggregation kinetics but not prion-seeding in a cell model.” The figure-input data are distributed separately through an associated Zenodo record.

## Repository contents

- `src/` contains the Python scripts for data processing, statistical analysis, and figure generation.
- `output/` contains generated figures and analysis products.
- `other_analysis_resources/` contains CellProfiler projects and related example resources.

## Data

Download the `input_data/` directory from the associated Zenodo record [`10.5281/zenodo.24271618`] and place it in the repository root before running the analysis scripts. Preserve the directory structure in the data archive so that the scripts can find their inputs. The data archive includes processed figure-input tables and examples for selected analysis steps; it does not include every original experimental measurement or microscopy image.

## Environment
Included in the repository are 4 .yml environment files. 
To run Cellpose via `src/cellpose/run_cellpose.py` please use the `cellpose-environment.yml` which has cross-platform reproducibility, or the frozen version also included (with specific dependency and build details).
To run all other scripts, please use `environment.yml` or refer to the corresponding frozen version. This environment has more standard packages for plotting and statistics. 


## Reproducing analyses

Run scripts from the repository root so that their relative paths resolve as expected. The scripts read prepared inputs from `input_data/` and write generated files under `output/`. Some upstream workflows require additional experimental files or manual preparation before they can be run.

The [script guide](README_scripts.md) describes the purpose and inputs of the analysis scripts. Figure-level scripts can be run independently when their corresponding processed input tables are available. The Figure 1 plate-reader workflow is documented under `src/Figure_1/data_analysis/`; the Cellpose environment specification is provided under `src/cellpose/`.

## Manuscript

ALS-associated mutations in the TDP-43 C-terminal domain alter in vitro aggregation kinetics but not prion-seeding in a cell model.

