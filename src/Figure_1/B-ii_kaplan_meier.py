
"""Plot lag-time Kaplan-Meier curves and compare protein groups statistically.

Reads the prepared lag-time/event table for Figure 1B-ii, exports the survival
curves, and writes pairwise log-rank tests with Benjamini-Hochberg adjustment.
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from src.paperfig_utils import save_figure_svg_png, set_params_for_plotting
from lifelines import KaplanMeierFitter 
from itertools import combinations
from lifelines.statistics import logrank_test
from statsmodels.stats.multitest import multipletests

if __name__ == '__main__':

    input_folder = 'input_data/Figure_1/'

    output_folder = 'output/Figure_1/B-ii_kaplan_meier/'

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    set_params_for_plotting(fontsize=10)

    # Read lag times and event indicators prepared by the curve-fitting workflow.
    combined = pd.read_csv(f'{input_folder}B-ii_tlag_km.csv')

    #define palette for manuscript
    palette_dict = { 'LCD' : 'Purples', 'A382P':'Greens', 'G335D': 'Reds', 'W385G':'Blues'}

    # Get a single color from the palette name — middle color of 6-tone palette
    #dictionary comprehension
    single_color_palette = {
    #this line assigns to the 'protein' variable we have looped over, 6 RGB color tuples to each protein, but then grabs only the middle (3rd) one so it's not too dark or light.
    protein: sns.color_palette(palette_dict[protein], n_colors=6)[3]
    #for each protein in my previously described dictionary
    for protein in palette_dict
    }


    #load in the kaplan meier function from lifelines
    kmf = KaplanMeierFitter()
    #prepare to plot/initialise figure object
    Fig, ax = plt.subplots(figsize=(10, 6))
    # Group by desired condition (protein + treatment, or just treatment)
    summary_rows = []
    for name, group in combined.groupby(['protein']):
        # Use fitted lag time as the survival duration.
        T = group['tlag']
        # Use the supplied event indicator to distinguish observations from censoring.
        E = group['event']
        label = f'{name[0]}'
        #fit the kaplan meier function to data
        kmf.fit(T, E, label=label)
        # Median survival (lag phase) time
        median = kmf.median_survival_time_
        # Quantiles from KM curve - to get information about the steepness of the curve i.e. once it starts aggregating, how long until MOST wells have begun etc.
        t10 = kmf.percentile(0.90)  # 10% events = 90% survival
        t50 = kmf.percentile(0.50)
        t90 = kmf.percentile(0.10)

        summary_rows.append({
            "protein": name,
            "median_lag": median,
            "t10": t10,
            "t50": t50,
            "t90": t90,
            "cascade_width": t90 - t10 #this is asking - do the events cascade at the same rate
        })
        #plot kaplan meier as survival curve
        kmf.plot(ci_show=False, color=single_color_palette.get(label, 'gray'))

    summary_df = pd.DataFrame(summary_rows)
    #save

    #plotting km curves all on same figure
    plt.xlabel("Time (h)")
    plt.ylabel("Fraction Not Aggregated")
    plt.tight_layout()
    plt.legend(bbox_to_anchor=(1.0,1.0))

    save_figure_svg_png(Fig, f'{output_folder}E_kaplan_meier_plot', width=3.5)


    # Compare each pair of protein groups using the log-rank test.
    #Log-rank test → compares entire survival processes
    groups = combined['protein'].unique() #group by variant
    results = []
    for g1, g2 in combinations(groups, 2):
        # Take one pair of proteins at a time so every unique pair is tested once.
        d1 = combined[combined['protein'] == g1]
        # Subset the dataframe to rows belonging to the second protein in the pair.
        d2 = combined[combined['protein'] == g2]

        # Run a two-sample log-rank test comparing time-to-event distributions.
        res = logrank_test(
            d1['tlag'], d2['tlag'],
            # Event indicator for group A (1=event observed, 0=censored).
            event_observed_A=d1['event'],
            # Event indicator for group B (1=event observed, 0=censored).
            event_observed_B=d2['event']
        )

        # Store a compact summary of this pairwise test for later tabulation/output.
        results.append({
            # Name of the first protein in the pairwise comparison.
            "group1": g1,
            # Name of the second protein in the pairwise comparison.
            "group2": g2,
            # p-value from the log-rank test (evidence for survival-curve difference).
            "p": res.p_value,
            # True if this comparison meets the significance threshold p < 0.01.
            "significant_p_lt_0_01": res.p_value < 0.01,
            # Test statistic from the log-rank test (magnitude of separation).
            "stat": res.test_statistic
        })

    pairwise_df = pd.DataFrame(results)
    print(pairwise_df)

    # Adjust the pairwise p-values for multiple comparisons.
    pairwise_df["p_fdr"] = multipletests(pairwise_df["p"], method="fdr_bh")[1]
    #final 'is corrected p value < 0.01???) column
    pairwise_df['p_fdr_corrected_0_01'] = pairwise_df['p_fdr']<0.01
    pairwise_df.to_csv(f'{output_folder}Fig1B-ii_kaplan_stats.csv')
