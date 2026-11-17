"""Plot aggregate-positive cell measurements for Figure 2B.

Reads the prepared replicate-level table, calculates within-experiment ratios
to the monomer control, and exports the treatment comparison and statistical
summary.
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from src.paperfig_utils import save_figure_svg_png, set_params_for_plotting
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from scipy import stats
    
def bardotplot(fold_increase_df, custom_palette, orders, ax):
    """Draw treatment-level bars and replicate points on the supplied axes."""
    sns.barplot(data=fold_increase_df,
        x='plot_label', y='value',
        palette=custom_palette,
        order=orders,
        alpha=0.7,
        dodge=True, edgecolor='black', errcolor='black', linewidth=1,capsize=0.15,
        ax=ax)
    sns.stripplot(data=fold_increase_df, x='plot_label', y='value',order=orders, dodge=True, jitter=False, color='#fff',s=6,edgecolor='black', alpha=0.6, linewidth=1, ax=ax)
    ax.set_ylim(0,10)
    ax.set_ylabel('Cells with seeds (%)')
    ax.set_xlabel('Time (h)')
    plt.tight_layout()


def execute_anova_with_tukey(fold_increase_df):
    """Run one-way ANOVA and Tukey HSD on the plotted treatment groups."""
    stats_df = fold_increase_df[['plot_label', 'value']].dropna()
    group_values = [g['value'].to_numpy()for _, g in stats_df.groupby('plot_label')]
    # Run one-way ANOVA
    f_stat, p_value = stats.f_oneway(*group_values)
    print(f"One-way ANOVA F = {f_stat:.3f}, p = {p_value:.5g}")

    # Run Tukey's HSD test
    tukey = pairwise_tukeyhsd(
    endog=stats_df['value'],
    groups=stats_df['plot_label'],
    alpha=0.05
    )
    
    return tukey

def calculate_foldchange_seeding(input_data):
    """Calculate each experiment's response relative to its monomer control."""
    fold_increase_df = []
    for replicate, df in input_data.groupby('Experiment'):
        df['fold_increase'] = df['value']/df['value'][df['treat']=='cont_monomer_trans'].values[0]
        fold_increase_df.append(df)
    fold_increase_df = pd.concat(fold_increase_df)
    return fold_increase_df

if __name__ == '__main__':

    input_folder = 'input_data/Figure_2/'
    output_folder = 'output/Figure_2/C_transuntrans/'

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    set_params_for_plotting(fontsize=10)

    # Read cleaned measurements collated across experimental replicates.
    input_data = pd.read_csv(f'{input_folder}Fig2C_seeding_transfection.csv')

    #normalise to the monomer control within each replicate to account for variability between days - option to plot but final plot currently plots % positive, not normalised.
    fold_increase_df = calculate_foldchange_seeding(input_data)
    
    #rename treatment names for plotting purposes
    plot_dict = {
        'cont_vehicle_trans':'Vehicle', 'med_10s_notrans':'UT', 'med_10s_trans':'T',
       'cont_monomer_trans':'Monomer', 'cont_BSA_trans':'BSA'
    }
    fold_increase_df['plot_label'] = fold_increase_df['treat'].map(plot_dict)

    # Define the palette used for treatment categories.
    hex_colors = ["#4a006f","#470a77","#45147e","#421e86","#3f288d","#3d3195","#3a3b9c","#3745a4","#354fab","#3259b3"]
    custom_palette = sns.color_palette(hex_colors)
    sns.palplot(custom_palette)


    orders = ['Vehicle', 'BSA', 'Monomer', 'T', 'UT']
    # Plot the measurements and save the figure.
    Fig, ax = plt.subplots()
    bardotplot(fold_increase_df, custom_palette, orders, ax)
    save_figure_svg_png(Fig, f'{output_folder}2C_seeded-cells_transfection', width=3.5)

 
    #statistics 
    # Group values by treatment
    tukey = execute_anova_with_tukey(fold_increase_df)
    print(tukey.summary())
