"""Combine per-replicate widefield outputs into collated summary tables.

This script searches replicate output folders for the expected CSV files from
the widefield analysis pipeline, concatenates matching tables across replicates,
and produces summary plots for aggregate counts, normalized seeding, and
fluorescence intensity ratios.
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns



def gather_replicates(input_folder, rep_type, file_I_want):
    """Collect matching CSV files from replicate subfolders.

    Args:
        input_folder (str): Root directory containing replicate outputs.
        rep_type (str): Substring used to select replicate folders.
        file_I_want (str): Filename fragment used to select the desired CSV.

    Returns:
        pandas.DataFrame: Concatenated dataframe of all matching files.
    """
    filepaths = [[f'{root}/{name}' for name in files] for root, dirs, files in os.walk(f'{input_folder}/')]

    filepaths = [item for sublist in filepaths for item in sublist if rep_type in item and '.csv' in item]
    y = []
    for name in filepaths:
        name = name.replace('\\', '/')
        exp_num = name.split('/')[-2]
        if file_I_want in name:
            x = pd.read_csv(f'{name}')
            if 'Experiment' not in x.columns.tolist():
                x['Experiment'] = exp_num
            
            y.append(x)
    y=pd.concat(y)
    return y

def violins(output_folder, data, col_of_interest, pal, ylabel, save=False, saveas = '', figtype='png',ylim=(0,1)):
    """Plot a violin summary for a chosen response variable."""
    Fig, ax = plt.subplots()
    sns.violinplot(data, x='plot_treat', y=col_of_interest, palette=pal, ax=ax)
    ax.set_ylim(ylim)
    ax.set_ylabel(f'{ylabel}')
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles=handles[0:], labels=labels[0:])
    if save==True:
        plt.savefig(f'{output_folder}{saveas}.{figtype}')

def bar(d, x_var, y_var, hue_var, ylim, ylab, save=False, plt_title='', pal='BuGn',figtype='png'):
    """Plot a bar-and-dot summary with a hue dimension."""
    Fig, ax = plt.subplots()
    ax = sns.stripplot(x=x_var, y=y_var, data=d, dodge=True, jitter=True, size=8, palette=pal, hue=hue_var)
    sns.barplot(x=x_var, y=y_var, data=d, capsize=.1, hue=hue_var, alpha=0.6, linewidth=2, palette=pal, edgecolor='0', errorbar="se", ax=ax)
    handles, labels = ax.get_legend_handles_labels()
    l = plt.legend(handles[0:3], labels[0:3], bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    ax.tick_params(axis='x', rotation=55)
    ax.set_xlabel('')
    ax.set_ylim(0,ylim)
    ax.set_ylabel(ylab)
    ax.set_title(plt_title)
    plt.tight_layout()
    if save==True:
        plt.savefig(f'{output_folder}{plt_title}.{figtype}')
    plt.show()

def bar_2(d, x_var, y_var, hue_var, ylim, ylab, save=False, plt_title='', pal='BuGn',figtype='png'):
    """Plot a simplified bar-and-dot summary without a hue-dense legend."""
    Fig, ax = plt.subplots()
    ax = sns.stripplot(x=x_var, y=y_var, data=d, jitter=True, size=10, palette=pal)
    sns.barplot(x=x_var, y=y_var, data=d, capsize=.1, alpha=0.6, linewidth=2, palette=pal, edgecolor='0', errorbar="sd", ax=ax)
    handles, labels = ax.get_legend_handles_labels()
    l = plt.legend(handles[0:1], labels[0:1], bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    ax.tick_params(axis='x', rotation=55)
    ax.set_xlabel('')
    ax.set_ylim(0,2)
    ax.set_ylabel(ylab)
    ax.set_title(plt_title)
    plt.tight_layout()
    return Fig, ax
   
def lineplot(d, x_var='Children_mRUBYaggs_big_Count', y_var='Children_Maskedsmall_Count', hue_var='plot_treat', ylim=(0,20), ylab='# of TDP43-2KQ-mRUBY puncta per cell', xlab='# of TDP43-2KQ-mRUBY aggregates per cell', save=False, plt_title='lineplot_smallcount-largecount.png', figtype='png'):
    """Plot paired aggregate counts per cell as a line summary."""
    Fig, ax = plt.subplots()
    ax = sns.lineplot(data=d, x=x_var, y=y_var, hue=hue_var)
    ax.set_ylim(ylim)
    ax.set_ylabel(ylab)
    ax.set_xlabel(xlab)
    if save == True:
        plt.savefig(f'{output_folder}{plt_title}.{figtype}')

if __name__ == '__main__':
    
    input_folder = 'results/'
    output_folder = 'results/collated_replicates/'

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Gather the per-replicate outputs produced by the upstream widefield
    # analysis scripts.
    ratios = gather_replicates(input_folder, rep_type='1_aggregate_count', file_I_want='count_melt.csv')

    normed_aggs = gather_replicates(input_folder, rep_type='1_aggregate_count', file_I_want='normed_amount_agg.csv')

    aggs_cells = gather_replicates(input_folder, rep_type='2_intensity_in-v-out', file_I_want='cellprofiler_output_cols_subset.csv')
    
    in_vs_out = gather_replicates(input_folder, rep_type='2_intensity_in-v-out', file_I_want='_Mean_in-v-out.csv')
    #violins(output_folder, ratios, col_of_interest='count_ratio_small_big', pal='BuPu', ylabel='puncta count : aggregate count', save=False, saveas='ratio_count_puncta_aggregate.png', ylim=(-1,10), figtype='png')


    custom_palette = ['#5184A4']  
    sns.set_palette(custom_palette)


    # Map compact treatment labels to reader-friendly names used in the final
    # figures.
    treats_dict = {'only':'untreated', 'vehicle':'vehicle', 'WT':'LCD monomer', 'old':'prior batch (WT)', 'UD-A382P':'A382P', 'UD-G335D':'G335D', 'UD-W385G':'W385G',
       'UD-WT':'LCD WT', 'normed-A382P': 'A382P (normalised)', 'normed-G335D': 'G335D (normalised)', 'normed-W385G':'W385G (normalised)'}
    
    cells_only = ratios[ratios['variable'] == 'perc_pos']
    cells_only['plot_treat'] = cells_only['t3'].map(treats_dict)


    counts_average = cells_only.groupby(['Experiment', 'plot_treat'])['value'].agg(
    mean='mean',
    std_dev='std'
    ).reset_index()

    orders = ['untreated', 'vehicle', 'prior batch (WT)', 'LCD WT', 'A382P', 'A382P (normalised)', 'G335D', 'G335D (normalised)',
        'W385G', 'W385G (normalised)']
    counts_average = counts_average[counts_average['Experiment']!='Exp30-5']
#
    Fig, ax = plt.subplots()
    sns.barplot(data=counts_average, x='plot_treat', y='mean', order=orders, palette=custom_palette, edgecolor='black', linewidth=1, capsize=0.3,ax=ax)
    sns.stripplot(x='plot_treat', y='mean', data=counts_average, order=orders, palette='viridis', hue='Experiment',size=7, edgecolor='black', linewidth=1, alpha=0.8,  ax=ax)
    ax.tick_params(axis='x', rotation=90)
    ax.set_xlabel('')
    ax.set_ylabel('Cells with aggregates (%)')
    ax.set_ylim(0,15)
    plt.tight_layout()
    plt.savefig(f'{output_folder}Agg_count_collated.png', bbox_inches = "tight")
    plt.savefig(f'{output_folder}Agg_count_collated.svg', bbox_inches = "tight")

    # Repeat the aggregate-count summary for the WT-normalized output.

    normed_aggs['plot_treat'] = normed_aggs['t3'].map(treats_dict)

    normed_aggs = normed_aggs[normed_aggs['Experiment']!='Exp30-5']
    counts_average = normed_aggs.groupby(['Experiment', 'plot_treat'])['normed_to_WT'].agg(
    mean='mean',
    std_dev='std'
    ).reset_index()

    Fig, ax = plt.subplots()
    sns.barplot(data=counts_average, x='plot_treat', y='mean', order=orders, palette=custom_palette, edgecolor='black', linewidth=1, capsize=0.3,ax=ax)
    sns.stripplot(x='plot_treat', y='mean', data=counts_average, order=orders, palette='viridis', hue='Experiment',size=7, edgecolor='black', linewidth=1, alpha=0.8,  ax=ax)
    ax.tick_params(axis='x', rotation=90)
    ax.set_xlabel('')
    ax.set_ylabel('Seeding relative to WT')
    ax.set_ylim(0,2.5)
    plt.tight_layout()
    plt.savefig(f'{output_folder}normed_relative_WT.png')
    plt.savefig(f'{output_folder}normed_relative_WT.svg')

    # Compare the intensity captured within aggregates across experiments.
    aggs_only = aggs_cells[aggs_cells['types']=='Aggs_collated']
    # Normalize within each experiment before plotting aggregate intensity
    # distributions.
    aggs_only['normalized_mean_intensity'] = aggs_only.groupby(['Experiment', 'types'])['Intensity_MeanIntensity_GaussCorr_'].transform(lambda x: x / x.max())

    aggs_only[['t1', 't2', 't3']] = aggs_only['treat'].str.split('_', expand=True)
    aggs_only['plot_treat'] = aggs_only['t3'].map(treats_dict)

    # Plot the per-experiment aggregate-intensity distribution.
    for experiment, df in aggs_only.groupby('Experiment'):
        fig,ax = plt.subplots()
        sns.violinplot(data=df, x='plot_treat', y='normalized_mean_intensity')
        ax.tick_params(axis='x', rotation=90)
        ax.set_ylim(0, 1.2)
        plt.tight_layout()
        plt.title(f'{experiment}')
        plt.savefig(f'{output_folder}norm_mean_intensity_box_collated.png', bbox_inches = "tight")
        plt.savefig(f'{output_folder}norm_mean_intensity_box_collated.svg', bbox_inches = "tight")

    # Average the per-experiment distributions for the final summary plot.
    aggs_cells_aggonly_average = aggs_cells[aggs_cells['types']=='Aggs_collated'].groupby(['Experiment', 't3'])['normalized_mean_intensity'].agg(
    mean='mean',
    std_dev='std'
    ).reset_index()
    
    aggs_cells_aggonly_average['plot_treat'] = aggs_cells_aggonly_average['t3'].map(treats_dict)

    Fig, ax = plt.subplots()
    sns.barplot(data=aggs_cells_aggonly_average, x='plot_treat', y='mean',palette=custom_palette, edgecolor='black', linewidth=1, capsize=0.3, order=orders,ax=ax)
    sns.stripplot(x='plot_treat', y='mean', data=aggs_cells_aggonly_average, order=orders, palette='Blues', hue='Experiment', size=7, edgecolor='black', linewidth=1, alpha=0.8,  ax=ax)
    ax.tick_params(axis='x', rotation=90)
    ax.set_xlabel('')
    ax.set_ylabel('Relative aggregate intensity')
    ax.set_ylim(0.4,1)
    plt.tight_layout()
    plt.savefig(f'{output_folder}avg_norm_mean_intensity_bar_collated.png', bbox_inches = "tight")
    plt.savefig(f'{output_folder}avg_norm_mean_intensity_bar_collated.svg', bbox_inches = "tight")


    
    in_vs_out[['t1', 't2', 't3']] = in_vs_out['treat'].str.split('_', expand=True)
    in_vs_out['plot_treat'] = in_vs_out['t3'].map(treats_dict)

    # Plot the aggregate-to-diffuse intensity ratio distribution.

    fig,ax = plt.subplots()
    sns.violinplot(data=in_vs_out, x='plot_treat', y='ratio_agg_total', palette=custom_palette, density_norm='area',order=orders)
    ax.tick_params(axis='x', rotation=90)
    ax.set_ylim(-0.2, 55)

    plt.tight_layout()
    plt.ylabel('aggregate:diffuse fluorescence')
    plt.savefig(f'{output_folder}intensity_in-vs-out_box_collated.png', bbox_inches = "tight")
    plt.savefig(f'{output_folder}intensity_in-vs-out_box_collated.svg', bbox_inches = "tight")


    # Average the ratio values across experiments for the final plot.
    ratio_in_out_av = in_vs_out.groupby(['Experiment', 't3'])['ratio_agg_total'].agg(
    mean='mean',
    std_dev='std'
    ).reset_index()
    
    ratio_in_out_av['plot_treat'] = ratio_in_out_av['t3'].map(treats_dict)

    Fig, ax = plt.subplots()
    sns.barplot(data=ratio_in_out_av, x='plot_treat', y='mean', order=orders, palette=custom_palette, edgecolor='black', linewidth=1, capsize=0.3, ax=ax)
    sns.stripplot(x='plot_treat', y='mean', data=ratio_in_out_av, order=orders, hue='Experiment', palette='Blues', size=7, edgecolor='black', linewidth=1, alpha=0.8, ax=ax)
    ax.tick_params(axis='x', rotation=90)
    ax.set_xlabel('')
    ax.set_ylabel('Ratio of aggregate:diffuse fluorescence')
    ax.set_ylim(-0.2,30)
    plt.tight_layout()
    plt.savefig(f'{output_folder}intensity_in-vs-out_bar_collated.png', bbox_inches = "tight")
    plt.savefig(f'{output_folder}intensity_in-vs-out_bar_collated.svg', bbox_inches = "tight")



    
    from statsmodels.stats.multicomp import pairwise_tukeyhsd
    from scipy import stats

    normed_aggs = normed_aggs[normed_aggs['variable']=='perc_pos']

    normed_aggs_avg = normed_aggs.groupby(['Experiment', 'treat'])['normed_to_WT'].agg(
    mean='mean',
    std_dev='std'
    ).reset_index()
    
    
    groups = [group['value'].values for name, group in normed_aggs.groupby('plot_treat')]

    f_stat, p_value = stats.f_oneway(*groups)

    print(f"F-statistic: {f_stat:.3f}")
    print(f"P-value: {p_value:.5f}")


    tukey = pairwise_tukeyhsd(endog=counts_average['mean'], groups=counts_average['plot_treat'], alpha=0.05)

    print(tukey)

    normed_aggs.to_csv(f'{output_folder}collated_mutants_seeded_normtoWT.csv')