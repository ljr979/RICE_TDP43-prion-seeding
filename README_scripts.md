# Analysis and figure scripts

This guide summarizes the analysis scripts in `src/`, their roles, and the input data associated with each workflow. Paths are relative to the repository root. Download the separate Zenodo data archive and place its `input_data/` directory in the repository root before running scripts that use these files. Figure-generation scripts use prepared CSV tables; upstream analysis steps may require raw images, plate-reader exports, or other experiment-specific files.

## Figure 1: plate-reader assays

The Figure 1 analysis workflow processes plate-reader fluorescence measurements and fits aggregation curves.

| Script | Purpose | Input |
|---|---|---|
| [`src/Figure_1/data_analysis/0_process_raw.py`](src/Figure_1/data_analysis/0_process_raw.py) | Labels wells using a plate map and processes fluorescence measurements over time. | Raw data and plate map in `input_data/Figure_1/data_analysis/0_process_raw/`. |
| [`src/Figure_1/data_analysis/1_fit_curves.py`](src/Figure_1/data_analysis/1_fit_curves.py) | Fits curves to processed fluorescence traces and calculates kinetic parameters. | `input_data/Figure_1/data_analysis/1_fit_curves/1_change_in_intensity_plus_info.csv`. |
| [`src/Figure_1/data_analysis/2_visualise_params.py`](src/Figure_1/data_analysis/2_visualise_params.py) | Plots fitted lag-time and elongation-rate parameters. | `input_data/Figure_1/data_analysis/2_visualise_params/fit_parameters.csv`. |
| [`src/Figure_1/B_THT-assays.py`](src/Figure_1/B_THT-assays.py) | Plots ThT fluorescence time courses. | `input_data/Figure_1/1_change_in_intensity_plus_info.csv`. |
| [`src/Figure_1/B-ii_kaplan_meier.py`](src/Figure_1/B-ii_kaplan_meier.py) | Plots lag-time Kaplan–Meier curves and calculates pairwise log-rank statistics. | `input_data/Figure_1/B-ii_tlag_km.csv`. |
| [`src/Figure_1/Fig1C-Suppfig1C_fit_params.py`](src/Figure_1/Fig1C-Suppfig1C_fit_params.py) | Plots fitted kinetic parameters and reports statistical summaries. | `input_data/Figure_1/SuppFig1C_fitparams_elongation-rate.csv`. |

The example plate-reader export and plate map are from replicate 1. The fitting and parameter-visualization directories contain prepared example tables for those steps.

## Figure 2: cell-seeding analyses

| Script | Purpose | Input |
|---|---|---|
| [`src/Figure_2/C_seeded-cells_transuntrans.py`](src/Figure_2/C_seeded-cells_transuntrans.py) | Compares aggregate-positive cell measurements across seeding and transfection conditions. | `input_data/Figure_2/Fig2C_seeding_transfection.csv`. |
| [`src/Figure_2/F_dNLS_vs_parental.py`](src/Figure_2/F_dNLS_vs_parental.py) | Compares aggregate-positive measurements in dNLS and parental cells. | `input_data/Figure_2/Fig2F_dNLSvsparental.csv`. |
| [`src/Figure_2/H_variant_seeding.py`](src/Figure_2/H_variant_seeding.py) | Compares variant-seeding measurements. | `input_data/Figure_2/Fig2H_variant-seeding.csv`. |
| [`src/Figure_2/data_analysis/2_normalise_to_wt_figure2.py`](src/Figure_2/data_analysis/2_normalise_to_wt_figure2.py) | Normalizes aggregate-positive measurements to a WT reference. | Aggregate-count data from the shared widefield workflow. |

The shared widefield analysis scripts in [`src/aggregate_identification/data_analysis/widefield/`](src/aggregate_identification/data_analysis/widefield/) collect CellProfiler exports, count aggregate-positive cells, and calculate aggregate-to-cell intensity measures. Example input exports and a plate map are provided under `input_data/aggregate_identification/data_analysis/widefield/`.

## Figure 3: particle number, variants, and fractions

| Script | Purpose | Input |
|---|---|---|
| [`src/Figure_3/B_sonicate_transfect.py`](src/Figure_3/B_sonicate_transfect.py) | Plots seeding measurements for sonicated and control samples. | `input_data/Figure_3/Fig3B_sonicated_aggcount.csv`. |
| [`src/Figure_3/D_E_particle_titration.py`](src/Figure_3/D_E_particle_titration.py) | Plots particle-titration and dilution-response measurements. | `input_data/Figure_3/Fig3D_particle_titration.csv`. |
| [`src/Figure_3/F_matrix_normed_seeding.py`](src/Figure_3/F_matrix_normed_seeding.py) | Plots the particle-number-normalized comparison matrix. | `input_data/Figure_3/Fig3F_matrix_seeding_normed.csv`. |
| [`src/Figure_3/G_variants_normalised.py`](src/Figure_3/G_variants_normalised.py) | Plots normalized variant-seeding measurements. | `input_data/Figure_3/Fig3G_variant_particlenum_norm.csv`. |
| [`src/Figure_3/H_variants-soluble-insoluble-seeding.py`](src/Figure_3/H_variants-soluble-insoluble-seeding.py) | Plots seeding measurements for whole, soluble, and insoluble fractions. | `input_data/Figure_3/Fig3H_relative_to_whole.csv`. |
| [`src/Figure_3/SuppFig5_inclusion_intensity.py`](src/Figure_3/SuppFig5_inclusion_intensity.py) | Plots inclusion-intensity measurements. | `input_data/Figure_3/SuppFig5_inclusion_intensity.csv`. |
| [`src/Figure_3/SuppFig7A_normed_particles.py`](src/Figure_3/SuppFig7A_normed_particles.py) | Plots particle-number-normalized seeding measurements. | `input_data/Figure_3/SuppFig7A_particle_num_seeding.csv`. |

The Figure 3 data-analysis scripts calculate relative-seeding and soluble/insoluble-fraction summaries from intermediate aggregate-count tables in `input_data/Figure_3/data_analysis/`.

## Figure 4: fibrils

| Script | Purpose | Input |
|---|---|---|
| [`src/Figure_4/B_fibrils_dNLS_coloc.py`](src/Figure_4/B_fibrils_dNLS_coloc.py) | Plots labelled-fibril and inclusion colocalization measurements. | `input_data/Figure_4/Fig4B_collated_coloc_fibs.csv`. |
| [`src/Figure_4/D-E_digested_fibs_seeding.py`](src/Figure_4/D-E_digested_fibs_seeding.py) | Plots digested-fibril seeding measurements and paired comparisons. | `input_data/Figure_4/Fig4D_digested_fibrils.csv` and `input_data/Figure_4/Fig4Di_allreps_paired_dotplot_input.csv`. |

Figure 4 data-analysis scripts include aggregate-count and fold-change processing steps under `src/Figure_4/data_analysis/`.

## Shared tools

| Script | Purpose |
|---|---|
| [`src/paperfig_utils.py`](src/paperfig_utils.py) | Shared plotting, palette, and figure-export helpers. |
| [`src/remove_metadata.py`](src/remove_metadata.py) | Removes selected path and filename metadata columns from CSV files in `other_analysis_resources/cellprofiler_pipelines/`; the script overwrites those files. |
| [`src/cellpose/run_cellpose.py`](src/cellpose/run_cellpose.py) | Processes TIFF image channels with Cellpose to generate cell masks. See the [Cellpose workflow notes](src/cellpose/README.md). |

The script inputs and outputs are configured in each source file. Figure-level tables are available under `input_data/`; the repository does not include every raw experimental measurement or microscopy image needed to regenerate all processed tables.
