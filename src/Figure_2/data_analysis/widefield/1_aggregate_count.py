"""Count aggregate-positive cells and summarize the fraction per treatment.

This script reads the collated CellProfiler cell table, assigns each well to a
treatment using a hard-coded layout map, calculates aggregate-positive cell
fractions, and exports summary tables and plots for each time point.
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

exp_num = 'Exp12-2'
input_folder = f'results/0_collated_cellprofiler/{exp_num}'
output_folder = f'results/1_aggregate_count/{exp_num}/'


if not os.path.exists(output_folder):
   os.makedirs(output_folder)

cells = pd.read_csv(f'{input_folder}/Cell_collated.csv')

# Keep only the columns needed to identify wells, treatments, and aggregate
# positivity.
subset = cells[['ImageNumber',
 'ObjectNumber','Metadata_Folder1',
 'Metadata_Folder2','Metadata_Region','Children_Aggs_Count',
 'Intensity_IntegratedIntensity_GaussCorr_',
 'Intensity_MeanIntensity_GaussCorr_',
 'Number_Object_Number',
 'exp_num', 'tp',
 'Unique_parent']]

subset.dropna(subset=['Metadata_Folder2'], inplace=True)
subset['Metadata_Folder2'] = subset['Metadata_Folder2'].astype(int)
subset['fold1_fold2'] = subset['Metadata_Folder1'].astype(str) + '_' + subset['Metadata_Folder2'].astype(str)

treatments = {
    'B_10':'med_10s_notrans', 'B_11':'med_10s_notrans', 'B_2':'cont_0s_only', 'B_3':'cont_0s_only', 'B_4':'cont_0s_only', 'B_5':'cont_0s_cells', 'B_6':'cont_0s_only', 'B_7':'med_10s_notrans', 'B_8':'med_10s_notrans',
       'B_9':'med_10s_notrans', 'C_10':'old_10s_trans', 'C_11':'old_10s_trans', 'C_2':'cont_0s_trans', 'C_3':'cont_0s_trans', 'C_4':'cont_0s_trans', 'C_5':'cont_0s_trans', 'C_6':'cont_0s_trans', 'C_7':'old_10s_trans',
       'C_8':'old_10s_trans', 'C_9':'old_10s_trans', 'D_10':'old_10s_notrans', 'D_11':'old_10s_notrans', 'D_2':'small_10s_trans', 'D_3':'small_10s_trans', 'D_4':'small_10s_trans', 'D_5':'small_10s_trans', 'D_6':'small_10s_trans',
       'D_7':'old_10s_notrans', 'D_8':'old_10s_notrans', 'D_9':'old_10s_notrans', 'E_2':'med_10s_trans', 'E_3':'med_10s_trans', 'E_4':'med_10s_trans', 'E_5':'med_10s_trans', 'E_6':'med_10s_trans', 'F_2':'cont_0s_buffer',
       'F_3':'cont_0s_buffer', 'F_4':'cont_0s_buffer', 'F_5':'cont_0s_buffer', 'F_6':'cont_0s_buffer', 'G_2':'small_10s_notrans', 'G_3':'small_10s_notrans', 'G_4':'small_10s_notrans', 'G_5':'small_10s_notrans', 'G_6':'small_10s_notrans'
}
subset['treat'] = subset['fold1_fold2'].map(treatments)

# Count the number of aggregate-positive cells per time point and treatment.
counts = []
for tp, y in subset.groupby('tp'):
    for treatment, df in y.groupby('treat'):
        num_cells_pos = len(df[df['Children_aggs_Count'] > 0])
        total_cells = len(df)


        counts.append([tp, treatment, num_cells_pos, total_cells])

counts = pd.DataFrame(counts,  columns=['tp','treat', 'positive_count', 'total_cells'])
counts['prop_pos'] = counts['positive_count']/counts['total_cells']
counts['perc_pos'] = counts['prop_pos']*100

# Normalize to the maximum observed fraction so the relative scale can be
# plotted alongside the absolute percentage values.
maxo = max(counts['prop_pos'])
counts['rel_prop'] = counts['prop_pos']/maxo

count_melt = pd.melt(counts, id_vars=['treat', 'tp'], value_vars=['perc_pos' , 'rel_prop'])

count_melt.to_csv(f'{output_folder}{exp_num}_count_melt.csv')


count_melt[['t1', 't2', 't3'] ]= count_melt[f'treat'].str.split(f'_',expand=True)
for tp, df in count_melt.groupby('tp'):
    Fig, ax = plt.subplots()
    sns.barplot(data=df[df['variable']=='perc_pos'], x='t3', y='value', order=['only', 'buffer', 'notrans', 'trans'], palette='BuGn',hue='t1', hue_order=['cont', 'old', 'small', 'med'], ax=ax)
    sns.stripplot(data=df[df['variable']=='perc_pos'], x='t3', y='value', dodge=True, order=['only', 'buffer', 'notrans', 'trans'], palette='BuGn',hue='t1', hue_order=['cont', 'old', 'small', 'med'], ax=ax, edgecolor='black', alpha=0.8, linewidth=1)
    ax.set_ylim(0,max(df['value']+1))
    ax.set_ylabel('Percentage of cells with aggregates')
    ax.tick_params(axis='x', rotation=45)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles=handles[0:], labels=labels[0:])
    plt.savefig(f'{output_folder}{tp}_Aggs_count.png')
    plt.savefig(f'{output_folder}{tp}_Aggs_count.svg')

    Fig, ax = plt.subplots()
    sns.barplot(data=df[df['variable']=='rel_prop'], x='t3', y='value', order=['only', 'buffer', 'notrans', 'trans'], palette='BuGn',hue='t1', hue_order=['cont', 'old', 'small', 'med'], ax=ax)
    ax.set_ylim(0,1)
    ax.set_ylabel('relative to each other')
    ax.tick_params(axis='x', rotation=45)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles=handles[0:], labels=labels[0:])