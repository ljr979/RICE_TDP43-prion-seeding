"""Plot aggregate-to-cell intensity ratios by batch and concentration.

Reads the processed inclusion-intensity table and exports the Figure 3E
violin-and-point plot.
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from statannotations.Annotator import Annotator
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from scipy import stats 
from src.paperfig_utils import save_figure_svg_png, set_params_for_plotting

def violindotplot(data, xcol, ycol, order, hue=None, hue_order=None, scat_hue=None, scat_hue_order=None, palette=False, xlabel='', ylabel=False, pairs=False, correction=None, xticks=None, groups=None, group_label_y=-0.18, group_line_y=-0.05, ax=None, legend='', dot_size=5):
    """Draw a violin plot with individual observations and optional labels."""
    if ax is None:
        fig, ax = plt.subplots()
        sns.violinplot(
    data=data,
    x=xcol,
    y=ycol,
    hue=hue,
    hue_order=hue_order,
    palette=palette,
    ax=ax,
    order=order,
    edgecolor='black'
    )
    sns.stripplot(
    data=data,
    x=xcol,
    y=ycol,
    hue=scat_hue,
    color='white',
    alpha=0.3,
    ax=ax,
    edgecolor='black',
    linewidth=1,
    s=dot_size,
    order=order,
    hue_order=scat_hue_order,
    dodge=True,
)


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


    input_folder ='input_data/Figure_3/'

    output_folder = 'output/Figure_3/E_inclusion_intensity/'

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    set_params_for_plotting(fontsize=10)

    # Read the prepared inclusion-intensity measurements.
    collated_intensity = pd.read_csv(f'{input_folder}Fig3E_inclusion_intensity.csv')

    titrations = collated_intensity[collated_intensity['type'].isin(['base', 'test'])]

    xaxis_order=['batch1', 'batch2', 'batch3']

    original_palette = 'Purples'
    # For each protein, assign a base palette
    base_colors = sns.color_palette(original_palette, n_colors=2)
    concentration_palette = {'0.3': base_colors[1], '0.15': base_colors[0]}

    fig, ax = violindotplot(titrations, 'original_batch', 'ratio_agg_total', xaxis_order, hue='concentration', hue_order=['0.3', '0.15'], palette=concentration_palette,  xlabel='',ylabel='Aggregate intensity:whole cell intensity', scat_hue='concentration', scat_hue_order=['0.3', '0.15'], pairs=False, correction=None, xticks=None, groups=None, group_label_y=-0.18, group_line_y=-0.05, ax=None, legend=True, dot_size=2)
    ax.tick_params(axis='x', rotation=90)
    ax.set_ylim(-5,90)
    plt.tight_layout()
    save_figure_svg_png(fig, f'{output_folder}3E_agg_intensity_cell_ratio', width=3.5)