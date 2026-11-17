"""Plot seeding measurements for whole, soluble, and insoluble fractions.

Reads the prepared Figure 3H summary, applies the manuscript palette, and
exports a fraction-level comparison for the TDP-43 variants.
"""

import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from statannotations.Annotator import Annotator
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from scipy import stats 
import matplotlib.patches as mpatches
from src.paperfig_utils import save_figure_svg_png, set_params_for_plotting

def bardotplot(data, xcol, ycol, order, hue=None, hue_order=None, scat_hue=None, scat_hue_order=None, palette=False, xlabel='', ylabel=False, pairs=False, correction=None, xticks=None, groups=None, group_label_y=-0.18, group_line_y=-0.05, ax=None, legend='', dot_size=5, cap_size=0.2, cap_width=2):
    """Plot bars with individual observations and optional statistical labels.

    Args:
        data: Dataframe containing the measurements and plotting metadata.
        xcol: Column used for x-axis categories.
        ycol: Column containing the plotted response.
        order: Display order for x-axis categories.
        hue: Optional column used to split bars into groups.
        hue_order: Display order for hue groups.
        scat_hue: Optional column used to color individual observations.
        scat_hue_order: Display order for scatter groups.
        palette: Seaborn palette specification.
        xlabel: X-axis label.
        ylabel: Y-axis label.
        pairs: Optional group pairs for statistical annotation.
        correction: Multiple-comparison correction passed to the annotator.
        xticks: Optional x-axis tick positions.
        groups: Optional grouped-axis labels and line coordinates.
        group_label_y: Vertical position for optional group labels.
        group_line_y: Vertical position for optional group lines.
        ax: Existing Matplotlib axes, or ``None`` to create axes.
        legend: Empty string hides the legend; another value displays it.
        dot_size: Size of individual observations.
        cap_size: Width of bar error-bar caps.
        cap_width: Line width of bar error bars.

    Returns:
        tuple: The Matplotlib figure and axes.
    """
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
    """Print one-way ANOVA and Tukey HSD results for the selected columns.

    Args:
        df: Dataframe containing observations and grouping labels.
        value_col: Column containing the numeric response.
        grouping_col: Column defining the groups to compare.
    """
    groups = [group[value_col].values for name, group in df.groupby(grouping_col)]
    # Run one-way ANOVA
    f_stat, p_value = stats.f_oneway(*groups)

    print(f"F-statistic: {f_stat:.3f}")
    print(f"P-value: {p_value:.5f}")


    # Run Tukey's HSD test
    tukey = pairwise_tukeyhsd(endog=df[value_col], groups=df[grouping_col], alpha=0.05)

    # Print summary
    print(tukey)

def set_custom_mutant_palette(xaxis_order, t1_order=['whole', 'insoluble', 'soluble']):
    """Build a protein/fraction palette for the grouped bar plot.

    Args:
        xaxis_order: Protein labels in plotting order.
        t1_order: Fraction labels used to select shades from each palette.

    Returns:
        dict: Palette keyed by ``(protein, fraction)`` pairs.
    """
    palette_dict = { 'WT' : 'Purples', 'A382P':'Greens', 'G335D': 'Reds', 'W385G':'Blues', 'x': 'Oranges'
    }

    t1_order = ['whole', 'insoluble', 'soluble']


    # This will hold the full palette
    combined_palette = {}
    # For each protein, assign a base palette
    for protein in xaxis_order:
        base_colors = sns.color_palette(palette_dict[protein], n_colors=len(t1_order))
        for i, treatment in enumerate(t1_order):
            combined_palette[(protein, treatment)] = base_colors[i]
    return combined_palette

if __name__ == '__main__':


    input_folder = 'input_data/Figure_3/'

    output_folder = 'output/Figure_3/H_variants-solinsol_seeding/'

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    set_params_for_plotting(fontsize=10)
    # Read the prepared whole, soluble, and insoluble fraction summaries.
    collated_solinsol_seeding = pd.read_csv(f'{input_folder}Fig3H_relative_to_whole.csv')
    # Remove any CSV index column from the input table.
    collated_solinsol_seeding.drop([col for col in collated_solinsol_seeding.columns.tolist() if 'Unnamed:' in col], axis=1, inplace=True)

    # Define a purple palette for the plotted conditions.
    hex_colors = ["#4a006f","#460677","#45147e","#421e86","#3f288d","#3d3195","#3a3b9c","#3745a4","#354fab","#3259b3"]
    custom_palette = sns.color_palette(hex_colors)
    sns.palplot(custom_palette)
    # Set the variant display order.
    xaxis_order = ['WT', 'G335D', 'A382P', 'W385G']

    # Plot seeding values normalized to the wild-type particle number.
    fig, ax = bardotplot(collated_solinsol_seeding, 'protein', 'mean', xaxis_order, palette=custom_palette,  hue=
    'plot_treat', hue_order=['whole', 'insoluble', 'soluble'], xlabel='', scat_hue='plot_treat', scat_hue_order=['whole', 'insoluble', 'soluble'],ylabel='Seeding (%) relative to whole fraction', pairs=False, correction=None, xticks=None, groups=None, group_label_y=-0.18, group_line_y=-0.05, ax=None, legend='', dot_size=5, cap_size=0.15, cap_width=2)
    # Build the manuscript palette for variant and fraction combinations.
    combined_palette = set_custom_mutant_palette(xaxis_order,t1_order=['whole', 'insoluble', 'soluble'])
    bar_colors = [combined_palette[(protein, 'insoluble')] for protein in xaxis_order]
    newcol_dict = dict(zip(xaxis_order, bar_colors))
    # Encode fraction type with fill patterns and variant with color.
    for stype, patch in zip([('whole', mut) for mut in xaxis_order]+ [('insoluble', mut) for mut in xaxis_order]+[('soluble', mut) for mut in xaxis_order], ax.patches):
        if stype[0] =='insoluble':
            patch.set_facecolor(newcol_dict[stype[1]])
            patch.set_edgecolor('#fff')
            patch.set_linewidth(1)
            patch.set_hatch('///')

        if stype[0] =='soluble':
            patch.set_facecolor(newcol_dict[stype[1]])
            patch.set_edgecolor('#fff')
            patch.set_linewidth(1)
            patch.set_hatch('o')

        else:
            patch.set_facecolor(newcol_dict[stype[1]])
            patch.set_edgecolor('#fff')
            patch.set_linewidth(1)


    #create legends to match colour/pattern scheme
    mon_patch = mpatches.Patch(facecolor='dimgray', label='whole', edgecolor='#fff', linewidth=1)
    seed_patch = mpatches.Patch(facecolor='dimgray', label='insoluble', edgecolor='#fff', linewidth=0, hatch='///')
    sol_patch = mpatches.Patch(facecolor='dimgray', label='soluble', edgecolor='#fff', linewidth=0, hatch='o')

    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(reversed(labels), reversed(handles)))
    ax.legend([mon_patch, seed_patch, sol_patch], ['whole','insoluble', 'soluble'], frameon=False, loc='upper left', handlelength=2., handletextpad=0.25, bbox_to_anchor=(-0.005, 1.15), ncols=3)
    plt.tight_layout()
    #save figure
    save_figure_svg_png(fig, f'{output_folder}3_sup_seeding_mutants', width=3.5)

    #anova_tukeyshsd(df=collated_solinsol_seeding, grouping_col='plot_treat', value_col='mean')
                                        

                                        
