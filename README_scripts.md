# Script-by-script reproduction guide

This is an additional guide to the repository's existing `README.md`. It
describes the scripts currently under `src/`, their declared input and output paths, and the work that precedes each script. It does not change any script's analysis settings or behavior.

## Important reproducibility notes

- Run commands from the repository root unless a script explicitly requires a different working directory. Most paths in the scripts are relative paths.
- The final-figure scripts generally read processed CSV tables from
  `input_data/Figure_1/` through `input_data/Figure_4/`; these tables are
  already present. They do not reconstruct the upstream raw experiments.
- Several Figure 1 and Figure 2 analysis scripts point to files or folders
  that do not match the files currently in this repository. Those mismatches
  are listed as **FLAG** items below; they have not been guessed at or silently
  corrected.
- The widefield scripts expect intermediate data under `results/`, but no
  repository-root `results/` data tree is present. `output/` is not a
  substitute for that input tree.
- `src/cellpose/run_cellpose.py` points to a Z: drive outside this repository.
  No TIFF input images were found in the repository. Its example result images
  are outputs, not source images.
- Figure 1 data-analysis scripts have editable experiment identifiers and
  processing options in their source. Confirm those settings against the
  experiment being reproduced before running the pipeline.

## Figure 1

| Script | What it does | Input data | What comes before it |
|---|---|---|---|
| [`src/Figure_1/data_analysis/0_process_raw.py`](src/Figure_1/data_analysis/0_process_raw.py) | Converts a plate-reader CSV to treatment-labelled fluorescence traces, optionally averages technical replicates or subtracts background, and writes processed tables and QC plots. | Configured as `input_data/Figure_1/data_analysis/0_process_raw/rep_1_raw_data.csv` and `rep_1_plate_map.csv`. | Export the raw plate-reader measurements and prepare a plate map with matching well coordinates and sample identifiers. The data examples are from replicate 1.
| [`src/Figure_1/data_analysis/1_fit_curves.py`](src/Figure_1/data_analysis/1_fit_curves.py) | Fits candidate Boltzmann-based curves to normalized traces, records fit metrics, selects the best fit per sample, and saves fit plots and CSV summaries. | `input_data/Figure_1/data_analysis/1_fit_curves/1_change_in_intensity_plus_info.csv` (present). | Process the matching raw plate-reader export with `0_process_raw.py`; verify the resulting traces and the script's subset/truncation/quality-threshold settings. |
| [`src/Figure_1/data_analysis/2_visualise_params.py`](src/Figure_1/data_analysis/2_visualise_params.py) | Plots lag phase and elongation-rate parameters from fit results. | Configured as `input_data/Figure_1/data_analysis/3_visualise_params/fit_parameters.csv`. | Run `1_fit_curves.py` and place its final selected-fit table at the configured input. 
| [`src/Figure_1/B_THT-assays.py`](src/Figure_1/B_THT-assays.py) | Plots ThT fluorescence time courses for mutation and monomer/control groups. | Configured as `input_data/Figure_1/1_change_in_intensity_plus_info.csv`. | Produce the processed fluorescence table with `0_process_raw.py` (or otherwise supply the same processed format). **FLAG:** the repository's similarly named file is `input_data/Figure_1/data_analysis/1_fit_curves/1_change_in_intensity_plus_info.csv`; the configured top-level file is absent. Confirm the intended source and measurement interval before running. |
| [`src/Figure_1/B-ii_kaplan_meier.py`](src/Figure_1/B-ii_kaplan_meier.py) | Plots Kaplan–Meier curves for lag times and writes pairwise log-rank results with FDR-adjusted p-values. | `input_data/Figure_1/B-ii_tlag_km.csv` (present). | Supply or prepare the lag-time/event table. The script describes it as coming from curve fitting; confirm how the event/censoring column was assigned. |
| [`src/Figure_1/Fig1C-Suppfig1C_fit_params.py`](src/Figure_1/Fig1C-Suppfig1C_fit_params.py) | Plots fitted lag-time and elongation-rate parameters and prints ANOVA/Tukey summaries. | `input_data/Figure_1/SuppFig1C_fitparams_elongation-rate.csv` (present). | Run or otherwise provide the curated fit-parameter table with treatment and replicate columns. |

## Figure 2

| Script | What it does | Input data | What comes before it |
|---|---|---|---|
| [`src/Figure_2/C_seeded-cells_transuntrans.py`](src/Figure_2/C_seeded-cells_transuntrans.py) | Plots seeded-cell measurements by treatment and prints ANOVA/Tukey results. | `input_data/Figure_2/Fig2B_seeding_transfection.csv` (present). | Prepare/collate the replicate-level measurements and treatment labels. |
| [`src/Figure_2/F_dNLS_vs_parental.py`](src/Figure_2/F_dNLS_vs_parental.py) | Calculates aggregate-positive percentages per cell type, plots unnormalized and no-seed-normalized comparisons. | `input_data/Figure_2/Fig2F_dNLSvsparental.csv` (present). | Collate CellProfiler-derived cell counts across experiments and provide the required cell-type, seeding, replicate, and aggregate-count columns. The script excludes one hard-coded replicate from normalization. |
| [`src/Figure_2/H_variant_seeding.py`](src/Figure_2/H_variant_seeding.py) | Filters and averages variant-seeding data, plots the results, and prints ANOVA/Tukey output. | `input_data/Figure_2/Fig2H_variant-seeding.csv` (present). | Prepare normalized seeding measurements with `t3`, `plot_treat`, `Experiment`, and `normed_to_WT` columns. |
| [`src/Figure_2/data_analysis/widefield/0_collect_files.py`](src/Figure_2/data_analysis/widefield/0_collect_files.py) | Collects per-image CellProfiler `_Cell.` and `_Aggs.` CSV exports, maps image folders to treatments using a plate map, and writes combined tables. | Configured under `input_data/Figure_2/data_analysis/0_collect_files/replicate/`, including `plate_map.csv` and matching CellProfiler CSV exports. | Run the CellProfiler pipelines on the relevant microscopy images and prepare the plate map. **FLAG:** this configured input tree is not present. The script's experiment label is also a placeholder (`Expnum-replicate`). |
| [`src/Figure_2/data_analysis/widefield/1_aggregate_count.py`](src/Figure_2/data_analysis/widefield/1_aggregate_count.py) | Counts aggregate-positive cells and writes per-timepoint plots and a melted count table using a treatment layout embedded in the script. | `results/0_collated_cellprofiler/Exp12-2/Cell_collated.csv`. | Run the upstream CellProfiler collation and ensure its output is placed at the expected `results/` path with the required columns. **FLAG:** no such root-level `results/` input is present. The selected columns use `Children_Aggs_Count`, while the subsequent count expression refers to `Children_aggs_Count`; verify the intended CellProfiler field before running. |
| [`src/Figure_2/data_analysis/widefield/1_aggregate_count_plate-map.py`](src/Figure_2/data_analysis/widefield/1_aggregate_count_plate-map.py) | Counts aggregate-positive cells by timepoint, treatment, and well, then writes a count table and plots. | `results/0_collated_cellprofiler/Exp30-6/Cell_collated.csv`. | Supply a collated cell table that already contains `treat` and `fold1_fold2`. **FLAG:** no such root-level `results/` input is present. Its `Exp30-6` configuration differs from `1_aggregate_count.py`'s `Exp12-2`. |
| [`src/Figure_2/data_analysis/widefield/2_normalise_to_wt_figure2.py`](src/Figure_2/data_analysis/widefield/2_normalise_to_wt_figure2.py) | Normalizes percentage-positive counts to a WT reference within each timepoint. | `results/1_aggregate_count/Exp30-6/Exp30-6_count_melt.csv`. | Run the matching count script and provide a `t3` value matching the configured WT label (`UD-WT`). **FLAG:** the upstream table and `results/` tree are absent. Confirm the treatment coding is consistent before running. |
| [`src/Figure_2/data_analysis/widefield/3_intensity_in-v-out.py`](src/Figure_2/data_analysis/widefield/3_intensity_in-v-out.py) | Matches aggregate and cell intensity records, calculates ratios, and writes plots and intermediate tables. | `results/0_collated_cellprofiler/Exp30-6/Cell_collated.csv` and `Aggs_collated.csv`. | Run CellProfiler and collate both object tables with matching parent identifiers, experiment metadata, treatment, and intensity columns. **FLAG:** no such root-level `results/` inputs are present. |
| [`src/Figure_2/data_analysis/widefield/all_replicates.py`](src/Figure_2/data_analysis/widefield/all_replicates.py) | Collects count, normalized-count, and intensity tables across replicate folders and makes cross-replicate summaries. | Searches `results/` for the preceding count and intensity outputs. | Complete the replicate-level widefield steps and keep the expected outputs in `results/`. **FLAG:** the required tree is absent. The script creates `normalized_mean_intensity` on a filtered dataframe but later groups the original `aggs_cells` table by that field; verify the intended intermediate table before running. |

## Figure 3

These are figure-level plotting/statistics scripts. Their processed CSV inputs are
already present under `input_data/Figure_3/`; the scripts do not derive those
tables from raw microscopy or assay data.

| Script | What it does | Input data | What comes before it |
|---|---|---|---|
| [`src/Figure_3/B_sonicate_transfect.py`](src/Figure_3/B_sonicate_transfect.py) | Plots aggregate-positive measurements and values normalized to the unsonicated condition; prints ANOVA/Tukey results. | `Fig3B_sonicated_aggcount.csv` (present). | Generate and collate replicate-level sonication measurements. |
| [`src/Figure_3/D_particle_titration.py`](src/Figure_3/D_particle_titration.py) | Plots the measured effect of dilution by seed batch. | `Fig3D_particle_titration.csv` (present). | Prepare the dilution-response summary table with batch, dilution, and response fields. |
| [`src/Figure_3/E_inclusion_intensity.py`](src/Figure_3/E_inclusion_intensity.py) | Plots aggregate-to-cell intensity ratios by batch and concentration. | `Fig3E_inclusion_intensity.csv` (present). | Calculate or collate inclusion and whole-cell intensity measurements. |
| [`src/Figure_3/F_normed_particles.py`](src/Figure_3/F_normed_particles.py) | Plots the effect of particle-number normalization across batches and conditions. | `Fig3F_particle_num_seeding.csv` (present). | Calculate the original and normalized seeding responses and assemble the plotting columns. |
| [`src/Figure_3/G_variants_normalised.py`](src/Figure_3/G_variants_normalised.py) | Averages normalized variant-seeding values per experiment, plots the comparison, and prints ANOVA/Tukey results. | `Fig3G_variant_particlenum_norm.csv` (present). | Prepare per-experiment normalized values and treatment labels. |
| [`src/Figure_3/H_variants-soluble-insoluble-seeding.py`](src/Figure_3/H_variants-soluble-insoluble-seeding.py) | Plots seeding values for whole, insoluble, and soluble fractions. | `Fig3H_variant_sol_insol_seeded.csv` (present). | Prepare the fraction-level summary table with the plotting columns used by the script. |

## Figure 4

These scripts consume processed summary CSVs already present under
`input_data/Figure_4/`.

| Script | What it does | Input data | What comes before it |
|---|---|---|---|
| [`src/Figure_4/B_fibrils_dNLS_coloc.py`](src/Figure_4/B_fibrils_dNLS_coloc.py) | Plots colocalization measurements between labelled fibrils and inclusions. | `Fig4B_collated_coloc_fibs.csv` (present). | Process microscopy in CellProfiler and collate replicate/channel colocalization values. |
| [`src/Figure_4/E_digested_fibs_seeding.py`](src/Figure_4/E_digested_fibs_seeding.py) | Plots fold-change and paired digested/undigested seeding data, and prints ANOVA, paired t-test, and Wilcoxon summaries. | `Fig4E_digested_fibrils.csv` and `Fig4Ei_allreps_paired_dotplot_input.csv` (present). | Prepare the summary and paired-replicate tables. Confirm the paired identifiers and the intended one-sided test direction before interpreting results. |

## Other scripts in `src/`

| Script | What it does | Data needed / prior work |
|---|---|---|
| [`src/paperfig_utils.py`](src/paperfig_utils.py) | Provides shared figure-export, plotting-style, and palette helpers imported by figure scripts. | No data input. It is a support module, not a figure-generation step. |
| [`src/remove_metadata.py`](src/remove_metadata.py) | Recursively removes selected path/filename columns from CSV files under `other_analysis_resources/cellprofiler_pipelines/`. It writes back to the same filenames. | No `input_data/` input; it modifies those CSVs in place. Do not run on the only copy of any source data. It is not required for generating the plotted figures. |
| [`src/cellpose/run_cellpose.py`](src/cellpose/run_cellpose.py) | Adjusts TIFF image intensities, runs Cellpose, and writes mask TIFFs and optional overlays. | A Cellpose environment, raw microscopy TIFFs, and an appropriate input folder. **FLAG:** the source points to a Z: drive and a named SOD1 experiment, not to this repository's `input_data/`; no TIFFs are present. Confirm the dataset, channel layout, and model settings before treating this as part of the TDP-43 figure workflow. |

## Files to add for raw-data reproducibility

The repository currently contains processed figure-input CSVs but no TIFF image
files and no complete raw widefield CellProfiler input tree. Bring the following
data into the repository (or provide a stable, documented download location) if
the goal is to regenerate all processed inputs from original measurements:

1. **Figure 1 plate-reader source data:** the raw plate-reader CSV and matching
   plate map for every replicate used in the paper. The current pipeline
   configuration asks for `rep_3_raw_data.csv` and `rep_3_plate_map.csv`, but
   only `rep_1_raw_data.csv` and `rep_1_plate_map.csv` are currently present.
   The existing rep 1 pair alone cannot establish whether it represents all
   experiments used for the published fit tables.
2. **Figure 2 widefield source exports:** the microscopy images used for
   CellProfiler, each corresponding CellProfiler `_Cell.` and `_Aggs.` CSV
   export, and the matching `plate_map.csv` for each replicate/experiment. The
   scripts refer to these under `input_data/Figure_2/data_analysis/0_collect_files/`
   but that input tree is absent. The CellProfiler project files and example
   outputs under `other_analysis_resources/cellprofiler_pipelines/` are not a
   replacement for the experiment-specific source exports.
3. **Figure 2 replicate-level intermediate tables:** the `Cell_collated.csv`
   and `Aggs_collated.csv` inputs expected beneath `results/0_collated_cellprofiler/`,
   plus the count and intensity intermediate CSVs expected by the later
   widefield scripts. Alternatively, document and implement the precise
   upstream route that creates these files from item 2. **FLAG:** current
   `0_collect_files.py` writes to `output/`, while downstream scripts read from
   `results/`; the intended handoff is unclear.
4. **Cellpose source images, if that analysis belongs to this paper:** the
   original Leica LAS-X TIFF image set for the confirmed experiment, plus any
   required channel metadata. **FLAG:** the script currently identifies an
   Experiment115 SOD1 dataset and uses a local Z: path; I cannot determine which
   TDP-43 paper images, if any, should be substituted.
5. **Any upstream raw data for Figures 2–4:** the processed CSV tables used by
   the figure scripts are checked in, but the original assay measurements or
   image-derived measurements that generated those summaries are not all
   traceable from the current scripts. The required source files and processing
   history need confirmation from the experiment records rather than inference
   from the final table names.

## Existing processed CSV inputs

The direct figure-level tables currently present in `input_data/` are:

- Figure 1: `B-ii_tlag_km.csv`,
  `SuppFig1C_fitparams_elongation-rate.csv`.
- Figure 2: `Fig2B_seeding_transfection.csv`,
  `Fig2F_dNLSvsparental.csv`, `Fig2H_variant-seeding.csv`.
- Figure 3: `Fig3B_sonicated_aggcount.csv`,
  `Fig3D_particle_titration.csv`, `Fig3E_inclusion_intensity.csv`,
  `Fig3F_particle_num_seeding.csv`, `Fig3G_variant_particlenum_norm.csv`,
  `Fig3H_variant_sol_insol_seeded.csv`.
- Additional Figure 3 table: `Fig3DF_particle_num_seeding.csv` is present but is
  not read directly by any current script under `src/`; its intended use is
  unconfirmed.
- Figure 4: `Fig4B_collated_coloc_fibs.csv`,
  `Fig4E_digested_fibrils.csv`,
  `Fig4Ei_allreps_paired_dotplot_input.csv`.

There are also Figure 1 intermediate CSVs under
`input_data/Figure_1/data_analysis/`; see the path mismatches above before using
them with the configured scripts.

## Items requiring author confirmation

- Which Figure 1 raw replicate is the intended example/reproduction run:
  `rep_1` (currently present) or `rep_3` (configured in the script)?
- Should the Figure 1 visualization input folder be `2_visualise_params` or
  `3_visualise_params`?
- Should the Figure 1 ThT plot read the nested processed CSV or a separate
  top-level file?
- What experiment-specific data and file layout should replace the missing
  Figure 2 `results/` inputs, and how should `0_collect_files.py` hand its
  outputs to the count/intensity scripts?
- Does the Cellpose script belong to this TDP-43 manuscript workflow, and which
  TIFF dataset and channel configuration are correct?
- Which original, non-processed source files underlie the already checked-in
  Figure 2–4 summary CSVs?
