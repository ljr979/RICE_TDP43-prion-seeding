"""Script which does the replicate-wise analysis for figure 3H - meaning normalising the amount of seeding from the supernatant or pellet to the WHOLE fraction for that same sample, for a single replicate. This is the data that then gets collated across replicates for panel H.
This script plots based on 'count_melt.csv' file output from 1_aggregate_count
"""
import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger
from statannotations.Annotator import Annotator
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from scipy import stats
from matplotlib.colors import LogNorm

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

def save_figure_svg_png(fig, filename, width=3.5, height=None, dpi=300, tight=True):
    """
    Save a matplotlib figure as both SVG and PNG with standardized formatting.
    
    Parameters:
    - fig: matplotlib.figure.Figure
        The figure object to save.
    - filename: str
        Path (with or without extension). Saves as .svg and .png.
    - width: float
        Width in inches. Default is 3.5" (single column).
    - height: float or None
        Height in inches. If None, uses golden ratio.
    - dpi: int
        Dots per inch for PNG (raster export).
    - tight: bool
        Whether to use tight_layout before saving.
    """
    # Set height using golden ratio if not provided
    if height is None:
        height = width * 0.618  # Golden ratio

    fig.set_size_inches(width, height)

    if tight:
        plt.tight_layout()

    # Ensure filename has no extension so we can add both .svg and .png
    base, _ = os.path.splitext(filename)

    svg_path = f"{base}.svg"
    png_path = f"{base}.png"

    # Save SVG
    fig.savefig(svg_path, format='svg', dpi=dpi)
    # Save PNG
    fig.savefig(png_path, format='png', dpi=dpi)

    print(f"Saved SVG to: {svg_path}")
    print(f"Saved PNG to: {png_path}")

if __name__ == '__main__':
    
    input_folder = 'input_data/Figure_3/data_analysis/'
    output_folder = 'output/Figure_3/H_sol-insol_relative_seeding/data_analysis/'

    combined = pd.read_csv(f'{input_folder}3H_replicate_count_melt.csv')
    
    
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    rel_to_whole = []
    for protein, df in combined.groupby('protein'):
        protein
        df
        norm_val = df[df['types']=='UD']['value'].mean()
        df['relative_to_whole'] = df['value']/norm_val
        rel_to_whole.append(df)
    rel_to_whole = pd.concat(rel_to_whole)

    
    rel_to_whole.to_csv(f'{output_folder}H_relative-to-whole.csv')

    Fig, ax = plt.subplots()
    sns.barplot(data=rel_to_whole, x='protein', y='relative_to_whole', palette='Blues', hue='types', hue_order=['UD', 'sup', 'pel'],order=['WT', 'A382P', 'G335D', 'W385G'],capsize=0.1,linewidth=2, edgecolor='black', ax=ax)
    sns.stripplot(data=rel_to_whole, x='protein', y='relative_to_whole', palette='Blues',hue='types', hue_order=['UD', 'sup', 'pel'],order=['WT', 'A382P', 'G335D', 'W385G'], dodge=True, ax=ax, edgecolor='black', alpha=0.8, linewidth=1)
    #ax.set_ylim(0,1)
    #ax.set_ylim(0,15)
    ax.set_ylabel('Aggregation relative to unseparated fibrils')
    ax.tick_params(axis='x', rotation=90)

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles=handles[0:3], labels=labels[0:3])
    plt.tight_layout()
    plt.savefig(f'{output_folder}H_replicate_relative-seeding.png')

