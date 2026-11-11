#now calculate the % of each TYPE of aggregate that is colocalised with the other type of aggregate, and plot as a barplot.
#noting that example data in this repository doesn't include the full suite, just example files from a subset. 
import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger

exp_num = 'Expnum-replicate'
input_folder = f'output/Figure_4/B_fibrils_dNLS_coloc/data_analysis/0_collect_files/'
output_folder = f'output/Figure_4/B_fibrils_dNLS_coloc/data_analysis/1_aggregate_count/'


if not os.path.exists(output_folder):
   os.makedirs(output_folder)

data = pd.read_csv(f'{input_folder}/0_collated_LCD_dNLS_cyto-aggs.csv')

#filter out irrelevant cellprofiler columns.
subset = data[['ImageNumber', 'ObjectNumber', 
       'Metadata_folder1', 'Metadata_sample',
       'Metadata_type', 'AreaShape_Area',
       'AreaShape_Compactness', 'AreaShape_Eccentricity',
       'AreaShape_MeanRadius', 'AreaShape_MedianRadius', 'AreaShape_Solidity',
       'Children_LCD_aggs_Count',
       'Intensity_IntegratedIntensity_tdp_correct_smooth',
       'Intensity_IntegratedIntensity_v5_chaps_gausscorr',
       'Intensity_MeanIntensity_tdp_correct_smooth',
       'Intensity_MeanIntensity_v5_chaps_gausscorr',
       'Neighbors_FirstClosestDistance_Nuclei_Expanded',
       'Neighbors_FirstClosestObjectNumber_Nuclei_Expanded',
       'Parent_LCD_aggs', 'Parent_cytoplasm', 'exp_num', 'type', 'protein',
        'Children_Greenaggs_Count', 'Parent_Greenaggs',
       'Intensity_MedianIntensity_tdp_correct_smooth',
       'Intensity_MedianIntensity_v5_chaps_gausscorr']]
#save!
subset.to_csv(f'{output_folder}1_subset_data.csv')

#this gives the proportion of aggregates with one another.
counts = []
#loop over the sample (from the metadata)
for tp, y in subset.groupby('Metadata_type'):
    tp
    y
    print(tp)
    #filter so that it's not including the cytoplasm information for this sample- because we are only looking at aggregates alone right now not any cell by cell data.
    y = y[y['type']!='cytoplasm']
    for aggtype, df in y.groupby('type'):
        aggtype
        print(aggtype)
        total_count = len(df)
        if aggtype == 'Greenaggs':
            other_child = 'LCD_aggs'
        if aggtype == 'LCD_aggs':
            other_child = 'Greenaggs'
        coloc_count = len(df[df[f'Children_{other_child}_Count'] > 0])
        counts.append([tp, aggtype, coloc_count, total_count])

counts = pd.DataFrame(counts,  columns=['sample','Aggregate', 'coloc_count', 'total_count'])
counts['protein'] = counts['sample'].str.split('seeds_', expand=True)[0]
counts['prop_pos'] = counts['coloc_count']/counts['total_count']
counts['perc_pos'] = counts['prop_pos']*100

count_melt = pd.melt(counts, id_vars=['sample', 'protein','Aggregate'], value_vars=['perc_pos'])

count_melt.to_csv(f'{output_folder}{exp_num}_proportion_positive_aggs.csv')

plotting_names = {'Greenaggs':'dNLS-mGFP inclusions', 'LCD_aggs':'pre-formed fibril(LCD-Atto647)'}
count_melt['plot_sample'] = count_melt['Aggregate'].map(plotting_names)

custom_palette = {
    'dNLS-mGFP inclusions': '#14591D',
    'pre-formed fibril(LCD-Atto647)': '#660708'
    }


count_melt['protein'].replace({'dyeonlyno':'Atto647 alone'}, inplace=True)

Fig, ax = plt.subplots()
sns.barplot(data=count_melt, x='protein', y='value', hue='plot_sample', palette=custom_palette,  edgecolor='black', linewidth=2, ax=ax)
sns.stripplot(data=count_melt, x='protein', y='value', hue='plot_sample',  palette=custom_palette,  dodge=True, jitter=False, ax=ax, edgecolor='black', alpha=0.8, linewidth=1)
ax.set_xlabel('')
ax.set_ylim(0,100)
ax.set_ylabel('Colocalisation of aggregates with each other(%)')
ax.tick_params(axis='x', rotation=45)
handles, labels = ax.get_legend_handles_labels()
l = plt.legend(handles[0:2], labels[0:2], bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.title('')
plt.tight_layout()
plt.savefig(f'{output_folder}Aggs_count.png')
plt.savefig(f'{output_folder}Aggs_count.svg')
