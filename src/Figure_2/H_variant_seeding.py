"""Plot normalized variant-seeding measurements for Figure 2H.

Filters the prepared table to the non-particle-equalized conditions, summarizes
technical observations within each experiment, and exports the comparison plot.
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from statannotations.Annotator import Annotator
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from scipy import stats 
import matplotlib.patches as mpatches
from src.paperfig_utils import save_figure_svg_png, set_params_for_plotting

def bardotplot(data, xcol, ycol, order, hue=None, hue_order=None, scat_hue=None, scat_hue_order=None, palette=False, xlabel='', ylabel=False, pairs=False, correction=None, xticks=None, groups=None, group_label_y=-0.18, group_line_y=-0.05, ax=None, legend='', dot_size=5, cap_size=0.2, cap_width=2):
    """Draw a bar plot with individual observations and optional annotations."""
    if ax is None:
        fig, ax = plt.subplots()
    sns.barplot(
        data=data,
        x=xcol,
        y=ycol,
        hue=hue,
        palette=palette,
        capsize=cap_size,
        errwidth=cap_width,
        errcolor='black',
        ax=ax,
        dodge=True,
        order=order,
        hue_order=hue_order,
        edgecolor='black'
    )
    sns.stripplot(
        data=data,
        x=xcol,
        y=ycol,
        hue=scat_hue,
        color='white',
        ax=ax,
        edgecolor='black',
        linewidth=1,
        s=dot_size,
        order=order,
        hue_order=scat_hue_order,
        dodge=True,
    )

    if pairs:
        annotator = Annotator(
            ax=ax, pairs=pairs, data=data, x=xcol, y=ycol, order=order, hue=hue, hue_order=hue_order)
        annotator.configure(test='t-test_ind', text_format='star',
                            loc='inside', comparisons_correction=correction, line_width=0.5)
        annotator.apply_and_annotate()

    ax.set(ylabel=ylabel)

    ax.set_xlabel(xlabel)
    if xticks:
        ax.set_xticks(xticks)
        ax.set_xticklabels(hue_order*len(order))
    if groups:
        for group_label, (x0, x1, x2) in groups.items():
            ax.annotate(group_label, xy=(x0, group_label_y),
                        xycoords='data', ha='center', annotation_clip=False)
            trans = ax.get_xaxis_transform()
            ax.plot([x1, x2], [group_line_y, group_line_y],
                    color="black", transform=trans, clip_on=False)

    if legend == '':
        ax.legend('', frameon=False)
    else:    
        handles, labels = plt.gca().get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        ax.legend(by_label.values(), by_label.keys())
   
    return fig, ax

def anova_tukeyshsd(df, value_col, grouping_col):
    """Print one-way ANOVA and Tukey HSD results for the selected columns."""
    groups = [group[value_col].values for name, group in df.groupby(grouping_col)]
    # Run one-way ANOVA
    f_stat, p_value = stats.f_oneway(*groups)

    print(f"F-statistic: {f_stat:.3f}")
    print(f"P-value: {p_value:.5f}")


    # Run Tukey's HSD test
    tukey = pairwise_tukeyhsd(endog=df[value_col], groups=df[grouping_col], alpha=0.05)

    # Print summary
    print(tukey)

if __name__ == '__main__':


    input_folder = 'input_data/Figure_2/'

    output_folder = 'output/Figure_2/H_variants_seeding/'

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    set_params_for_plotting(fontsize=10)

    collated_coloc = pd.read_csv(f'{input_folder}Fig2H_variant-seeding.csv')

    # Keep only treatments that were not particle-equalized.
    not_norm_only = collated_coloc[collated_coloc['t3'].isin(['only', 'vehicle','UD-A382P', 'UD-G335D', 'UD-W385G','UD-WT','UD-mon'])]
    not_norm_only.drop([col for col in not_norm_only.columns.tolist() if 'Unnamed:' in col], axis=1 , inplace=True)

    # Average technical wells within each experiment and treatment.
    normed_aggs_avg = not_norm_only.groupby(['Experiment', 'plot_treat'])['normed_to_WT'].agg(
    mean='mean',
    std_dev='std'
    ).reset_index()
    
    palette_dict = { 

        'LCD WT' : 'Purples', 
        'A382P':'Greens', 
        'G335D': 'Reds', 
        'W385G':'Blues', 
        'vehicle': 'Greys',
        'untreated': 'Greys'

    }


    xaxis_order = ['untreated', 'vehicle', 'LCD WT', 'G335D', 'A382P','W385G']

    # This will hold the full palette
    combined_palette = {}
    # For each protein, assign a base palette
    for protein in xaxis_order:
        base_colors = sns.color_palette(palette_dict[protein], n_colors=3)
        combined_palette[protein] = base_colors[1]


    fig, ax = bardotplot(normed_aggs_avg, 'plot_treat', 'mean', xaxis_order, palette=combined_palette,  xlabel='', scat_hue='Experiment', ylabel='Fold difference in seeding (%) relative to WT', pairs=False, correction=None, xticks=None, groups=None, group_label_y=-0.18, group_line_y=-0.05, ax=None, legend='', dot_size=5, cap_size=0.15, cap_width=2)
    ax.tick_params(axis='x', rotation=90)

    save_figure_svg_png(fig, f'{output_folder}Fig2H_variants_seeding', width=3.5)

    anova_tukeyshsd(df=normed_aggs_avg, grouping_col='plot_treat', value_col='mean')
                                        
                                        
                                        