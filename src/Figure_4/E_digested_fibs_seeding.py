"""Plot Figure 4E digestion comparisons and paired replicate measurements.

Reads processed summary and paired data tables, exports unpaired and paired
plots, and prints the configured statistical tests.
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from statannotations.Annotator import Annotator
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from scipy import stats 
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from src.paperfig_utils import make_palette
from src.paperfig_utils import save_figure_svg_png, set_params_for_plotting
from scipy.stats import wilcoxon

def bardotplot(data, xcol, ycol, order, hue=None, hue_order=None, scat_hue=None, scat_hue_order=None, palette=False, xlabel='', ylabel=False, pairs=False, correction=None, xticks=None, groups=None, group_label_y=-0.18, group_line_y=-0.05, ax=None, legend='', dot_size=5, cap_size=0.2, cap_width=2):
    """Draw bars and individual observations with optional grouping labels."""
    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.figure
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

def plot_paired_points_nolines(paired_wide, protein_order_present, protein_palette, ax=None, figsize=(10, 6)):
    """Plot digestion states as paired categories without connecting lines."""
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    # Color = protein identity, Marker = digestion state.
    for protein, sub in paired_wide.groupby('protein'):
        color = protein_palette.get(protein, 'gray')
        ax.scatter(
            sub['x_undigested'],
            sub['undigested'],
            color=color,
            marker='o',
            edgecolor='black',
            s=70,
            zorder=3,
        )
        ax.scatter(
            sub['x_PKdigested'],
            sub['PKdigested'],
            color=color,
            marker='^',
            edgecolor='black',
            s=70,
            zorder=3,
        )

    # Label each observation with its experiment identifier.
    for _, row in paired_wide.iterrows():
        ax.text(row['x_undigested'] + 0.015, row['undigested'], row['Experiment'], fontsize=7, va='center')
        ax.text(row['x_PKdigested'] + 0.015, row['PKdigested'], row['Experiment'], fontsize=7, va='center')

    # Set category tick labels.
    ax.set_xticks(range(len(protein_order_present)))
    ax.set_xticklabels(protein_order_present, rotation=45, ha='right')
    ax.set_ylabel('Mean % cells with aggregates')
    ax.set_xlabel('Protein')
    ax.set_title('Paired observations without connecting lines')
    protein_handles = [
        mpatches.Patch(facecolor=protein_palette.get(p, 'gray'), edgecolor='black', label=p)
        for p in protein_order_present
    ]
    digestion_handles = [
        Line2D([0], [0], marker='o', color='white', markerfacecolor='white', markeredgecolor='black', linestyle='None', label='undigested'),
        Line2D([0], [0], marker='^', color='white', markerfacecolor='white', markeredgecolor='black', linestyle='None', label='PKdigested')
    ]
    first_legend = ax.legend(handles=protein_handles, title='Protein', loc='upper left', bbox_to_anchor=(1.02, 1.0))
    ax.add_artist(first_legend)
    ax.legend(handles=digestion_handles, title='State', loc='upper left', bbox_to_anchor=(1.02, 0.58))
    ax.grid(axis='y', linestyle='--', alpha=0.4)
    plt.tight_layout()
    return fig, ax

def plot_paired_points_lines(paired_wide, protein_order_present, protein_palette, ax=None, figsize=(10, 6)):
    """Plot paired digestion-state observations connected within replicates."""
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure


    # Connect each experiment's undigested and digested observations.
    for _, row in paired_wide.iterrows():
        ax.plot(
            [row['x_undigested'], row['x_PKdigested']],
            [row['undigested'], row['PKdigested']],
            color='grey',
            alpha=0.7,
            linewidth=1.5,
            zorder=2,
        )

    # Overlay the observations, using color for protein and shape for digestion.
    for protein, sub in paired_wide.groupby('protein'):
        color = protein_palette.get(protein, 'gray')
        ax.scatter(
            sub['x_undigested'],
            sub['undigested'],
            color=color,
            marker='o',
            edgecolor='black',
            s=35,
            zorder=3,
        )
        ax.scatter(
            sub['x_PKdigested'],
            sub['PKdigested'],
            color=color,
            marker='^',
            edgecolor='black',
            s=35,
            zorder=3,
        )

    # Replicate labels are omitted from this plot to reduce visual clutter.
    # for _, row in paired_wide.iterrows():
    #     ax.text(row['x_undigested'] + 0.015, row['undigested'], row['Experiment'], fontsize=7, va='center')
    #     ax.text(row['x_PKdigested'] + 0.015, row['PKdigested'], row['Experiment'], fontsize=7, va='center')

    # Format axes and legends for the paired comparison.
    ax.set_xticks(range(len(protein_order_present)))
    ax.set_xticklabels(protein_order_present, rotation=90)
    ax.set_ylabel('Cells with aggregates (%)')
    ax.set_xlabel('')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    protein_handles = [
        mpatches.Patch(facecolor=protein_palette.get(p, 'gray'), edgecolor='black', label=p)
        for p in protein_order_present
    ]
    digestion_handles = [
        Line2D([0], [0], marker='o', color='white', markerfacecolor='white', markeredgecolor='black', linestyle='None', label='undigested'),
        Line2D([0], [0], marker='^', color='white', markerfacecolor='white', markeredgecolor='black', linestyle='None', label='PKdigested')
    ]
    first_legend = ax.legend(handles=protein_handles, title='Protein', loc='upper left', bbox_to_anchor=(1.02, 1.0))
    ax.add_artist(first_legend)
    ax.legend(handles=digestion_handles, title='State', loc='upper left', bbox_to_anchor=(1.02, 0.1))

    plt.tight_layout()
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

    input_folder = 'input_data/Figure_4/'

    output_folder = 'output/Figure_4/E_digested_fibs_seeding/'

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    set_params_for_plotting(fontsize=10)

    collated_coloc = pd.read_csv(f'{input_folder}Fig4E_digested_fibrils.csv')
    
    xaxis_order = ['WT', 'G335D', 'A382P', 'W385G']

    # Use manuscript-consistent colors: one representative color per protein palette.
    paper_palette = make_palette()

    fig, ax = bardotplot(collated_coloc, 'protein', 'mean', xaxis_order, palette=paper_palette,  xlabel='', ylabel='Fold change in % cells with aggregates', pairs=False, correction=None, xticks=None, groups=None, group_label_y=-0.18, group_line_y=-0.05, ax=None, legend=True, dot_size=5, cap_size=0.15, cap_width=2)

    ax.tick_params(axis='x', rotation=90)
    ax.set_ylim(0,3.5)
    plt.tight_layout()
    save_figure_svg_png(fig, f'{output_folder}4_digested_foldchange_seeding', width=3.5)

    anova_tukeyshsd(df=collated_coloc, grouping_col='protein', value_col='mean')



    # Prepare the paired replicate-level plot and tests.
    paired_df = pd.read_csv(f'{input_folder}Fig4Ei_allreps_paired_dotplot_input.csv')

    #SANITY check just points alone no lines (plus annotation of experiment number)
    fig, ax = plot_paired_points_nolines(
        paired_wide=paired_df,
        protein_order_present=['WT', 'G335D', 'A382P', 'W385G'],
        protein_palette=paper_palette
    )
    save_figure_svg_png(fig, f'{output_folder}4Ei_digested_paired_nolines', width=3.5)

    # Plot the paired observations with connecting lines.
    fig, ax = plot_paired_points_lines(paired_wide=paired_df, protein_order_present=['WT', 'G335D', 'A382P', 'W385G'], protein_palette=paper_palette)
    save_figure_svg_png(fig, f'{output_folder}4Ei_prot_K_paired_lines', width=3.5)

    # Compare paired digestion states within each protein group.
    for protein, df in paired_df.groupby('protein'):
        t_statistic, p_value = stats.ttest_rel(df['undigested'], df['PKdigested'], alternative='less')
        print(f'T-TEST: {protein} - tstat = {t_statistic}, p_val = {p_value}')


        stat, p = wilcoxon(
            df["PKdigested"],
            df["undigested"]
        )
        print(f'wilcoxon: {protein} - stat = {stat}, p_val = {p}')
