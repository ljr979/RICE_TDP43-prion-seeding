import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger
from statannotations.Annotator import Annotator
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from scipy import stats 
import matplotlib.patches as mpatches
from src.paperfig_utils import make_palette
from src.paperfig_utils import save_figure_svg_png, set_params_for_plotting

def bardotplot(data, xcol, ycol, order, hue=None, hue_order=None, scat_hue=None, scat_hue_order=None, palette=False, xlabel='', ylabel=False, pairs=False, correction=None, xticks=None, groups=None, group_label_y=-0.18, group_line_y=-0.05, ax=None, legend='', dot_size=5, cap_size=0.2, cap_width=2):
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
    groups = [group[value_col].values for name, group in df.groupby(grouping_col)]
    # Run one-way ANOVA
    f_stat, p_value = stats.f_oneway(*groups)

    print(f"F-statistic: {f_stat:.3f}")
    print(f"P-value: {p_value:.5f}")


    # Run Tukey's HSD test
    tukey = pairwise_tukeyhsd(endog=df[value_col], groups=df[grouping_col], alpha=0.05)

    # Print summary
    print(tukey)
    
    palette_dict = palette_dict
    t1_order = t1_order

    # This will hold the full palette
    combined_palette = {}
    # For each protein, assign a base palette
    for protein in xaxis_order:
        protein
        print(palette_dict[protein])
        base_colors = sns.color_palette(palette_dict[protein], n_colors=len(t1_order))
        for i, treatment in enumerate(t1_order):
            combined_palette[(protein, treatment)] = base_colors[i]
    return combined_palette

if __name__ == '__main__':


    input_folder = 'input_data/Figure_4/'

    output_folder = 'output/Figure_4/B_fibrils_dNLS_coloc/'

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    set_params_for_plotting(fontsize=10)
    #read in collated data from all replicates (colocalisation in each channel. see data analysis and cellprofiler pipelines for all image processing and data wrangling)
    collated_coloc = pd.read_csv(f'{input_folder}Fig4B_collated_coloc_fibs.csv')

    #plotting order
    xaxis_order = ['Atto647 alone', 'WT', 'G335D', 'A382P', 'W385G']

    fig, ax = bardotplot(collated_coloc, 'protein', 'value', xaxis_order, palette='viridis',  hue='plot_sample', hue_order=['pre-formed fibril(LCD-Atto647)', 'dNLS-mGFP inclusions'], xlabel='', scat_hue='plot_sample',scat_hue_order=['pre-formed fibril(LCD-Atto647)', 'dNLS-mGFP inclusions'], ylabel='Colocalisation of aggregates with each other(%)', pairs=False, correction=None, xticks=None, groups=None, group_label_y=-0.18, group_line_y=-0.05, ax=None, legend=True, dot_size=5, cap_size=0.15, cap_width=2)

    ax.tick_params(axis='x', rotation=90)
    ax.set_ylim(0,100)

    palette = make_palette()
    #add lines to bar faces and make colour correct from paper palette
    for stype, patch in zip([('pre-formed fibril(LCD-Atto647)', mut) for mut in xaxis_order]+ [('dNLS-mGFP inclusions', mut) for mut in xaxis_order], ax.patches):

        print(stype[1])
        if stype[0] =='pre-formed fibril(LCD-Atto647)':
            patch.set_facecolor(palette[stype[1]])
            patch.set_edgecolor('#fff')
            patch.set_linewidth(1)
            patch.set_hatch('///')

        else:
            patch.set_facecolor(palette[stype[1]])
            patch.set_edgecolor('#fff')
            patch.set_linewidth(1)


    #Make legend match plotting
    mon_patch = mpatches.Patch(facecolor='dimgray', label='dNLS-mGFP inclusions', edgecolor='#fff', linewidth=1)
    seed_patch = mpatches.Patch(facecolor='dimgray', label='pre-formed fibril(LCD-Atto647)', edgecolor='#fff', linewidth=0, hatch='///')
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(reversed(labels), reversed(handles)))
    ax.legend([mon_patch, seed_patch], ['dNLS-mGFP inclusions','pre-formed fibril (Atto647)'], frameon=False, loc='upper left', handlelength=2., handletextpad=0.25, bbox_to_anchor=(-0.005, 1.15), ncols=3)

    plt.tight_layout()

    save_figure_svg_png(fig, f'{output_folder}4_coloc_label_fibrils', width=3.5)