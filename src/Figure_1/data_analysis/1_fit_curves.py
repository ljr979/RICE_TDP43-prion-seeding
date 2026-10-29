"""Fit sigmoidal models to fluorescence time-series data from seeding assays.

This script prepares time-resolved fluorescence traces, optionally truncates the
series to an early fitting window, and fits several constrained Boltzmann-based
models to each biological replicate. It then exports per-fit summary metrics,
including lag time, elongation rate, and goodness of fit, and produces figures
for visual inspection of each fitted curve.

The fitting strategy is intentionally constrained to avoid biologically
implausible solutions such as negative baselines, negative maxima, or negative
slopes. After fitting, the script compares model performance using R-squared and
retains the best-performing fits above a user-defined quality threshold.
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger
from scipy.optimize import curve_fit
from src.paperfig_utils import save_figure_svg_png, set_params_for_plotting

logger.info('Import ok')


def truncate_timeseries(df, desired_length_h, time_interval):
    """Return a copy of the input table limited to an early time window.

    Parameters
    ----------
    df : pandas.DataFrame
        Wide-format dataframe containing metadata columns and timepoint columns.
    desired_length_h : float
        Duration, in hours, to retain from the beginning of the time series.
    time_interval : float
        Sampling interval, in minutes, between consecutive timepoints.

    Returns
    -------
    pandas.DataFrame
        Dataframe containing the original metadata columns plus only the
        requested number of timepoint columns.
    """
    keep_timepoints = int((60 / time_interval) * desired_length_h)

    timepoint_cols = [col for col in change_in.columns if col not in info_cols]
    keep_timepoint_cols = timepoint_cols[:keep_timepoints]
    df_truncated = df[info_cols + keep_timepoint_cols]
    return df_truncated


def boltzmann(x, bottom, top, V50, slope):
    """Evaluate a standard four-parameter Boltzmann sigmoid."""
    return bottom + (top - bottom) / (1 + np.exp((V50 - x) / slope))


def boltzmann_with_exp_decay(x, bottom, top, V50, slope, decay):
    """Evaluate a Boltzmann sigmoid with an exponential decay term."""
    return boltzmann(x, bottom, top, V50, slope) * np.exp(-decay * x)


def boltzmann_with_linear(x, bottom, top, V50, slope, m):
    """Evaluate a Boltzmann sigmoid with an added linear component."""
    return boltzmann(x, bottom, top, V50, slope) + m * x


def boltzmann_with_both(x, bottom, top, V50, slope, decay, m):
    """Evaluate a Boltzmann sigmoid with decay and linear correction terms."""
    return (boltzmann(x, bottom, top, V50, slope) * np.exp(-decay * x)) + m * x


def r_squared_calculator(x, y, func, popt):
    """Calculate the coefficient of determination for a fitted model.

    Parameters
    ----------
    x : numpy.ndarray
        Independent variable values used for fitting.
    y : numpy.ndarray
        Observed dependent variable values.
    func : callable
        Model function used to generate fitted values.
    popt : sequence
        Optimized model parameters returned by ``curve_fit``.

    Returns
    -------
    float
        R-squared value for the fitted model.
    """
    residuals = y - func(x, *popt)
    ss_res = np.sum(residuals ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    return 1 - (ss_res / ss_tot)


def fit_curve(x_vals, y_vals, model_func, p0, bounds, inform, fit_type, extra_params=[], plot_failed=False):
    """Fit one candidate model to a single fluorescence trace.

    Parameters
    ----------
    x_vals : numpy.ndarray
        Time values for the trace.
    y_vals : numpy.ndarray
        Measured fluorescence values for the trace.
    model_func : callable
        Candidate model function to fit.
    p0 : list
        Initial parameter guesses supplied to ``curve_fit``.
    bounds : tuple[list, list]
        Lower and upper parameter bounds for constrained optimization.
    inform : str
        Human-readable identifier for the current sample and replicate.
    fit_type : str
        Name of the model being fitted.
    extra_params : list, optional
        Reserved for future extensions. Currently unused.
    plot_failed : bool, optional
        If ``True``, display traces for failed fits to aid debugging.

    Returns
    -------
    list | None
        A list containing the sample identifier, fitted parameters, derived fit
        metrics, fit type, and elongation rate if fitting succeeds; otherwise
        ``None``.
    """
    try:
        popt, pcov = curve_fit(
            model_func,
            x_vals,
            y_vals,
            p0=p0,
            bounds=bounds,
            method='trf',
            loss='soft_l1',
            f_scale=0.1
        )

        rsq = r_squared_calculator(x_vals, y_vals, model_func, popt)
        bottom, top, V50, slope = popt[:4]
        tlag = ((4 * slope * bottom) / (top - bottom)) + (V50 - (2 * slope))
        elongation = (top - bottom) / (4 * slope)

        return [inform, *popt, tlag, round(rsq, 3), fit_type, elongation]

    except RuntimeError as e:
        print(f"[{inform}] Fit failed ({fit_type}): {e}")
        if plot_failed:
            plt.figure(figsize=(6, 4))
            plt.plot(x_vals, y_vals, 'o-')
            plt.title(f"Failed Fit: {inform}\n{fit_type}")
            plt.xlabel("Time")
            plt.ylabel("Intensity")
            plt.show()
        return None


def is_valid_initial_guess(p0, bounds):
    """Return ``True`` if all initial guesses fall within the provided bounds."""
    lower, upper = bounds
    return all(l <= p <= u for p, l, u in zip(p0, lower, upper))


if __name__ == '__main__':

    filter_cols = []

    input_path = 'input_data/Figure_1/data_analysis/1_fit_curves/1_change_in_intensity_plus_info.csv'
    output_folder = 'output/Figure_1/data_analysis/1_fit_curves/'
    subset = True
    save_plot_fits = True
    set_params_for_plotting(fontsize=10)

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Load the processed fluorescence traces and retain the normalized signal
    # used for curve fitting.
    clusters_summary = pd.read_csv(input_path)
    clusters_summary = clusters_summary[clusters_summary['t1'] != 'seeds']
    change_in = clusters_summary[clusters_summary['type'] == 'normalised_to_max_change']
    info_cols = ['ID', 'treat2', 'treat1', 'type', 't1', 't2', 't3', 'replicate']
    change_in.drop([col for col in change_in.columns.tolist() if 'Unnamed: 0' in col], axis=1, inplace=True)

    truncate = True
    if truncate:
        # Optionally limit the fit to the early phase of the trace when late
        # plateaus or declines reduce fit stability.
        change_in = truncate_timeseries(df=change_in, desired_length_h=12, time_interval=5)

    # Convert wide-format timepoint columns into long format so each replicate
    # can be iterated as a standard x/y trace for fitting.
    clean_data = pd.melt(
        change_in,
        id_vars=info_cols,
        value_name='intensity',
        var_name='timepoint'
    )
    clean_data['timepoint'] = (clean_data['timepoint'].astype(int)) * 5 / 60
    clean_data.to_csv(f'{output_folder}clean_for_fitting.csv')

    fit_results = []
    no_fit = []

    # Optional subset mode is useful while testing parameter bounds or fit
    # behavior on a small number of example traces.
    test_list = ['xseeding-WT-1']
    if subset == True:
        clean_data['inform'] = clean_data[['t1', 't2', 't3']].apply(lambda row: '-'.join(row.values.astype(str)), axis=1)
        clean_data = clean_data[clean_data['inform'].isin(test_list)]

    for (sample, t2, t3), df in clean_data.groupby(['t1', 't2', 't3']):
        for replicate, df1 in df.groupby('replicate'):
            df1
            inform = f'{sample}-{t2}-{t3}-{replicate}'

            x_vals = df1['timepoint'].values
            y_vals = df1['intensity'].values
            # Clamp negative intensities to zero so constrained fits are not
            # distorted by small background-subtraction artifacts.
            y_vals[y_vals < 0] = 0
            if len(y_vals) < 4 or np.std(y_vals) < 1e-3 or np.any(np.isnan(y_vals)):
                no_fit.append((inform, 'skipped: too flat or invalid'))
                continue

            y_min, y_max = np.nanmin(y_vals), np.nanmax(y_vals)
            x_median = np.median(x_vals)
            x_max = np.max(x_vals)

            fit_configs = [
                {
                    "model": boltzmann,
                    "p0": [y_min, y_max, x_median, 1],
                    "bounds": ([0, 0, 0, 0], [np.inf, np.inf, x_max, 100]),
                    "type": "boltzmann",
                    "color": "red"
                },
                {
                    "model": boltzmann_with_exp_decay,
                    "p0": [y_min, y_max, x_median, 1, 0.001],
                    "bounds": ([0, 0, 0, 0, 0], [np.inf, np.inf, x_max, 100, 1]),
                    "type": "boltzmann_with_exp_decay",
                    "color": "green"
                },
                {
                    "model": boltzmann_with_linear,
                    "p0": [y_min, y_max, x_median, 1, 0],
                    "bounds": ([0, 0, 0, 0, -10], [np.inf, np.inf, x_max, 100, 10]),
                    "type": "boltzmann_with_linear",
                    "color": "blue"
                },
                {
                    "model": boltzmann_with_both,
                    "p0": [y_min, y_max, x_median, 1, 0.001, 0],
                    "bounds": ([0, 0, 0, 0, 0, -10], [np.inf, np.inf, x_max, 100, 1, 10]),
                    "type": "boltzmann_with_both",
                    "color": "purple"
                },
            ]

            # Plot the raw data once and overlay each candidate fit so model
            # performance can be inspected visually for every replicate.
            fig, ax = plt.subplots(figsize=(8, 5))
            ax.plot(x_vals, y_vals, 'ko', label='Data')
            success = False

            for fit in fit_configs:
                if not is_valid_initial_guess(fit["p0"], fit["bounds"]):
                    no_fit.append((inform, fit["type"] + ' — invalid p0'))
                    continue

                result = fit_curve(
                    x_vals,
                    y_vals,
                    model_func=fit["model"],
                    p0=fit["p0"],
                    bounds=fit["bounds"],
                    inform=inform,
                    fit_type=fit["type"],
                    plot_failed=False
                )

                if result:
                    subset_data = result[0:1] + result[-4:]
                    fit_results.append(subset_data)

                    popt = result[1:1 + len(fit["p0"])]
                    rsq = result[-3]

                    x_fit = np.linspace(min(x_vals), max(x_vals), 500)
                    y_fit = fit["model"](x_fit, *popt)

                    label = f"{fit['type']} (R²={rsq:.3f})"
                    ax.plot(x_fit, y_fit, '-', color=fit["color"], label=label)
                    success = True
                else:
                    no_fit.append((inform, fit["type"]))

            ax.set_title(f'Fits for {inform}')
            ax.set_xlabel('Timepoint')
            ax.set_ylabel('Fluorescence Intensity')

            if success:
                ax.legend()
            else:
                ax.text(
                    0.5, 0.5, 'No successful fits',
                    ha='center', va='center',
                    transform=ax.transAxes
                )

            fig.tight_layout()

            if save_plot_fits:
                save_figure_svg_png(fig, f'{output_folder}fitted_curve_{inform}', width=3.5, dpi=600)

            plt.show()

fit_params = pd.DataFrame(fit_results, columns=['inform','tlag', 'rsq', 'fit_type', 'elongation_rate'])



fit_params.to_csv(f'{output_folder}fit_parameters.csv')


# Compare R-squared values across candidate models to identify the most
# consistently successful fitting strategy.
melted_params = pd.melt(fit_params, id_vars=['inform',
 'tlag', 'fit_type', 'elongation_rate'
 ], value_vars=['rsq'])

sns.barplot(melted_params, x="fit_type", y="value")


for fit_type, df in melted_params.groupby('fit_type'):
    df
    fit_type
    amazing = df[df['value']>0.995]
    not_amazing = df[df['value']<0.995]

    prop_good = (len(amazing)/len(df))*100
    print(fit_type, prop_good)


# Retain the single best-performing model for each sample based on the maximum
# observed R-squared value.
good_fits = []
for sample, df in melted_params.groupby('inform'):
    sample
    df
    maximum = max(df['value'])
    best_fit = df[df['value']==maximum]
    if len(best_fit) > 1:
        best_fit = best_fit.iloc[:1]
    good_fits.append(best_fit)

good_fits = pd.concat(good_fits)
len(good_fits)

# Separate acceptable fits from lower-confidence fits using the fixed quality
# threshold applied throughout this analysis.
not_so_good_fits = good_fits[good_fits['value']<0.995]
best_filtered = good_fits[good_fits['value']>0.995]

perc = len(not_so_good_fits)/len(good_fits)*100
print(f'{perc} are not good enough fits to keep')
best_filtered.rename(columns={'value':'rquare'}, inplace=True)
best_filtered.drop([col for col in best_filtered.columns.tolist() if 'variable' in col], axis=1, inplace=True)

best_filtered.to_csv(f'{output_folder}fit_parameters.csv')
not_so_good_fits.to_csv(f'{output_folder}no_fit.csv')


