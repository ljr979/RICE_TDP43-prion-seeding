"""Compare aggregate-positive parental and dNLS cells for Figure 2F.

Calculates per-replicate positive-cell percentages, plots the unnormalized
values, then plots seeded values relative to the no-seed condition.
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from src.paperfig_utils import save_figure_svg_png, set_params_for_plotting

def calculate_cells_seeded(collated_data, agg_count):
    """Append aggregate-positive cell percentages for each experimental group.

    Args:
        collated_data: Collated CellProfiler cell data and treatment metadata.
        agg_count: List that receives one summary row per group.
    """
    for (celltype, exp_type, replicate), df in collated_data.groupby(['celltype', 'seed_noseed', 'expnum-replicate']):
        print(celltype, exp_type, replicate)
        #find the total number of cells for each replicate and cell type
        total_cells = len(df['ObjectNumber'].unique())
        if celltype == 'parental':
            #find the number of cells where aggregates were identified (this uses column names from cellprofiler output)
            cells_w_aggs = len(df[df['Children_aggs_init_Count']>0])
        elif celltype =='dNLS':
            #find number of cells identified via cellprofiler in the dnls cells
            cells_w_aggs = len(df[df['Children_Combinedaggs_Count']>0])
        print(f'total cells: {total_cells}, with aggs: {cells_w_aggs}')
        #calculate the percentage of cells which had aggregates
        perc_pos = ((cells_w_aggs/total_cells)*100)+1
        #add to dataframe
        agg_count.append([celltype, exp_type, replicate, perc_pos])

def nonnorm_barplot(agg_count, orders,  ax, yvar, pal, ylimit, ylabel, hue_var='seed_noseed'):
    """Draw bars and replicate points for the supplied cell-level summary."""
    use_hue = hue_var is not None

    sns.barplot(data=agg_count,
        x='celltype', y=yvar,
        hue=hue_var,
        alpha=0.7,
        dodge=use_hue, edgecolor='black', errcolor='black', linewidth=1,capsize=0.15,palette=pal,
        ax=ax)
    sns.stripplot(data=agg_count, x='celltype', y=yvar,order=orders,hue=hue_var, dodge=use_hue, jitter=False, color='#fff',s=6,edgecolor='black', alpha=0.6, linewidth=1, ax=ax)
    ax.set_ylim(0,ylimit)
    ax.set_ylabel(f"{ylabel}")
    ax.set_xlabel('Cell type')
    sns.despine()

    if use_hue:
        plt.legend(title='')
        handles, labels = ax.get_legend_handles_labels()
        l = plt.legend(handles[2:], labels[2:], bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        # Overwrite the text items dynamically
        new_labels = ["- seeds", "+ seeds"]
        for text, new_text in zip(l.get_texts(), new_labels):
            text.set_text(new_text)
    else:
        # No hue legend needed for single-bar plots.
        legend = ax.get_legend()
        if legend is not None:
            legend.remove()

if __name__ == '__main__':

    input_folder = 'input_data/Figure_2/'
    output_folder = 'output/Figure_2/F_dNLS_vs_parental/'

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    set_params_for_plotting(fontsize=10)

    # Read the collated CellProfiler measurements.
    collated_data = pd.read_csv(f'{input_folder}Fig2F_dNLSvsparental.csv')
    #find out the number of aggregates per experimental replicate, per treatment, per cell type.
    agg_count = []
    calculate_cells_seeded(collated_data, agg_count)
    agg_count = pd.DataFrame(agg_count, columns=['celltype', 'seed_noseed', 'expnum-replicate', 'percentage_positive_cells'])

    # Set the cell-type display order.
    orders = ['dNLS', 'parental']

    # Define the two-color palette used for the cell types.
    hex_colors = ["#4a006f","#3259b3"]
    custom_palette = sns.color_palette(hex_colors)
    sns.palplot(custom_palette)

    # Plot and save the unnormalized aggregate-positive percentages.
    Fig, ax = plt.subplots()
    nonnorm_barplot(agg_count, orders, ax, yvar='percentage_positive_cells', pal=custom_palette,ylimit=20, ylabel='Cells with aggregates (%)')
    save_figure_svg_png(Fig, f'{output_folder}SuppFig2_parental_vs_dNLS_non-normed', width=3.5, height=2.8)


    # Normalize seeded values to the no-seed condition within each replicate.
    # Exclude the replicate that lacks a no-seed measurement.
    agg_count = agg_count[agg_count['expnum-replicate']!='27-2']

    ratios = []
    for celltype, df in agg_count.groupby('celltype'):
        for replicate, dfs in df.groupby('expnum-replicate'):
            #divide all treatments by the no seeds treatment 
            divide_by = dfs[dfs['seed_noseed']=='noseeds']['percentage_positive_cells'].values[0]
            dfs['ratio'] = dfs['percentage_positive_cells']/divide_by
            ratios.append(dfs)
    
    ratios = pd.concat(ratios)

    


    orders = ['dNLS', 'parental']
    ratios_seedonly = ratios[ratios['seed_noseed']=='seeds']
    # Plot and save the seeded-to-no-seed ratios.
    Fig, ax = plt.subplots()
    nonnorm_barplot(ratios_seedonly, orders, ax, yvar='ratio', pal=custom_palette, ylimit=5, ylabel='Fold change in cells with aggregates', hue_var=None)
    save_figure_svg_png(Fig, f'{output_folder}2E_parental_vs_dNLS', width=3.5)
