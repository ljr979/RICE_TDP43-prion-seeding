"""Count aggregate-positive cells using a plate-map-based treatment layout.

This script reads a collated CellProfiler cell table, maps wells to treatments
using a plate map encoded in the script, summarizes aggregate-positive cell
fractions, and exports both per-well tables and comparison plots.
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

exp_num = 'Expnum-replicate'

input_folder = f'output/0_collect_files/{exp_num}'
output_folder = f'output/1_aggregate_count/{exp_num}/'
cells = pd.read_csv(f'{input_folder}/Cell_collated.csv')

if not os.path.exists(output_folder):
   os.makedirs(output_folder)

subset = cells[['ImageNumber',
 'ObjectNumber','Metadata_Folder1',
 'Metadata_Folder2','Metadata_Region','Children_Aggs_Count',
 'Intensity_IntegratedIntensity_GaussCorr_',
 'Intensity_MeanIntensity_GaussCorr_',
 'Number_Object_Number',
 'exp_num', 'tp',
 'Unique_parent', 'treat', 'fold1_fold2']]

subset.to_csv(f'{output_folder}cells_subset.csv')

# Count aggregate-positive cells within each time point, treatment, and well.
counts = []
for (tp,treatment), y in subset.groupby(['tp', 'treat']):
    for well, df in y.groupby('fold1_fold2'):
        num_cells_pos = len(df[df['Children_Aggs_Count'] > 0])
        total_cells = len(df)
        
        counts.append([tp, treatment, well, num_cells_pos, total_cells])

counts = pd.DataFrame(counts,  columns=['tp','treat', 'well','positive_count', 'total_cells'])
counts['prop_pos'] = counts['positive_count']/counts['total_cells']
counts['perc_pos'] = counts['prop_pos']*100

# Normalize to the maximum observed fraction so the relative scale can be used
# for comparison across conditions.
maxo = max(counts['prop_pos'])
counts['rel_prop'] = counts['prop_pos']/maxo

count_melt = pd.melt(counts, id_vars=['treat', 'tp', 'well'], value_vars=['perc_pos' , 'rel_prop'])
count_melt[['t1', 't2', 't3'] ]= count_melt[f'treat'].str.split(f'_',expand=True)
count_melt.to_csv(f'{output_folder}{exp_num.split('/')[-1]}_count_melt.csv')

# Plot the percentage of aggregate-positive cells for every treatment.
# Timepoint is shown as hue so the same code works for one or multiple timepoints.
perc_pos_data = count_melt[count_melt['variable'] == 'perc_pos'].copy()

Fig, ax = plt.subplots()
sns.barplot(data=perc_pos_data, x='treat', y='value', color='#c64d4d', edgecolor='black', ax=ax)
sns.stripplot(
    data=perc_pos_data,
    x='treat',
    y='value',
    hue='tp',
    dodge=True,
    palette='dark:gray',
    edgecolor='black',
    alpha=0.8,
    linewidth=1,
    ax=ax,
)
ax.set_ylim(0, max(perc_pos_data['value']) + 1)
ax.set_ylabel('Percentage of cells with aggregates')
ax.set_xlabel('Treatment')
ax.tick_params(axis='x', rotation=90)

handles, labels = ax.get_legend_handles_labels()
if labels:
    ax.legend(handles=handles, labels=labels, title='Timepoint', bbox_to_anchor=(1.02, 1), loc='upper left')

plt.tight_layout()
plt.savefig(f'{output_folder}all_treatments_Aggs_count.png')
plt.savefig(f'{output_folder}all_treatments_Aggs_count.svg')
