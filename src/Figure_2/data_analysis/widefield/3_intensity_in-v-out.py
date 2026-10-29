"""Quantify aggregate-localized fluorescence relative to total cell intensity.

This script combines CellProfiler cell and aggregate tables, matches each
aggregate to its parent cell, calculates the ratio of aggregate fluorescence to
whole-cell fluorescence, and exports summary tables and visualizations for each
condition.
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

exp_num = 'Exp30-6'

input_folder = 'results/0_collated_cellprofiler/'
output_folder = f'results/2_intensity_in-v-out/{exp_num}/'

if not os.path.exists(output_folder):
   os.makedirs(output_folder)


read_in=['Cell_collated', 'Aggs_collated']

alls = []
agg_neg = []
for d in read_in:
   x = pd.read_csv(f'{input_folder}{exp_num}/{d}.csv')
   x['types'] = d
   if d == 'Cell_collated':
   # Keep track of cells without aggregates so the final ratio table can
   # include zero-aggregate observations.
      neg = x[x['Children_Aggs_Count'] == 0]
      agg_neg.append(neg)
      pos = x[x['Children_Aggs_Count'] > 0]
      x = pos
   alls.append(x)

alls = pd.concat(alls)
alls.to_csv(f'{output_folder}collated_cellprofiler_output.csv')

agg_neg = pd.concat(agg_neg)


cell_subset = ['ImageNumber',
 'ObjectNumber',
 #'AreaShape_Area',
 #'AreaShape_Compactness',
 #'AreaShape_MeanRadius',
 'Children_Aggs_Count',
 'Intensity_IntegratedIntensity_GaussCorr_',
 
 'Intensity_MeanIntensity_GaussCorr_',
 'exp_num',
 'Unique_parent',
 'types', 'treat', 'fold1_fold2']


alls = alls[cell_subset]
agg_neg = agg_neg[cell_subset]
alls.to_csv(f'{output_folder}cellprofiler_output_cols_subset.csv')


#now want to find PER matched cell, the ratio of total cell intensity to intensity in small aggregates, big aggregates and combined
melted = pd.melt(alls, id_vars=['Unique_parent', 'ImageNumber', 'ObjectNumber',  'types', 'treat','fold1_fold2'], value_vars=[ 'Intensity_IntegratedIntensity_GaussCorr_','Intensity_MeanIntensity_GaussCorr_', 'Children_Aggs_Count'])
neg_melt = pd.melt(agg_neg, id_vars=['Unique_parent', 'ImageNumber', 'ObjectNumber',  'types', 'treat','fold1_fold2'], value_vars=[ 'Intensity_IntegratedIntensity_GaussCorr_','Intensity_MeanIntensity_GaussCorr_', 'Children_Aggs_Count'])
melted.to_csv(f'{output_folder}melted_aggs_cells.csv')
neg_melt.to_csv(f'{output_folder}cells_no_aggs.csv')


ratio_info = []
whats_up = []
for parent, df in melted.groupby('Unique_parent'):
   treat = df['treat'].unique()[0]
   intensity_df = df[df['variable']=='Intensity_MeanIntensity_GaussCorr_']
   if len(intensity_df) > 2:
        total = intensity_df['value'][intensity_df['types'] == 'Cell_collated'].values[0]
        agg_intensity_all = intensity_df['value'][intensity_df['types'].isin(['Aggs_collated'])].values.sum()
        inform = [parent, treat, total, agg_intensity_all]
        ratio_info.append(inform)

ratios = pd.DataFrame(ratio_info, columns=['Unique_parent', 'treat', 'total_intensity', 'agg_intensity'])



ratios['ratio_agg_total'] = ratios['agg_intensity']/ratios['total_intensity']

ratios[['Experiment', 'Treat', 'folder1','folder2', 'image', 'parent']] = ratios['Unique_parent'].str.split('_', expand=True)

#need to match the aggregate negative data back on
neg_melt[['Experiment', 'Treat', 'folder1', 'folder2', 'image', 'parent']] = neg_melt['Unique_parent'].str.split('_', expand=True)

neg_melt_intensity = neg_melt[neg_melt['variable']=='Intensity_MeanIntensity_GaussCorr_']

neg_melt_intensity.rename(columns={'value':'total_intensity'}, inplace=True)

neg_melt_intensity['agg_intensity'] = 0

neg_melt_intensity['ratio_agg_total'] = neg_melt_intensity['agg_intensity']/neg_melt_intensity['total_intensity']



neg_melt_intensity = neg_melt_intensity[['Unique_parent',
 'total_intensity',
 'agg_intensity',
 'ratio_agg_total',
 'Experiment',
 'treat',
 'Treat',
 'folder1',
 'folder2',
 'image',
 'parent']]


ratios_all = pd.concat([ratios, neg_melt_intensity])

#make a new df to look at the ones that have an impossible ratio (i.e. higher than the total cell intensity)
ratios_whack = ratios[ratios['ratio_agg_total']>1]

#select for only those that HAVE aggregates, and that the aggregate intensity is not higher than the CELL intensity (SHOULDN'T BE POSSIBLE?)
#not possible with MEAN intensity as the mean intensity of an agg is way brighter than the whole cell :) 
#ratios_aggs = ratios[ratios['ratio_agg_total']<=1]
ratios_aggs = ratios[ratios['ratio_agg_total']>0]
ratios_aggs['log_ratio_agg_total'] = np.log10(ratios_aggs['ratio_agg_total'])

# Also keep a combined table that includes zero-aggregate cells so readers can
# inspect the full distribution if needed.
ratios_all = ratios_all[ratios_all['ratio_agg_total']<=1]









# Plotting begins below.


ratios_aggs.to_csv(f'{output_folder}_Mean_in-v-out.csv')

# First plot only cells that contain aggregates, so the ratio reflects the
# aggregate-localized signal within aggregate-positive cells.
Fig, ax = plt.subplots()
sns.boxplot(data=ratios_aggs, x='treat', y='ratio_agg_total', palette='YlOrRd', ax=ax)
ax.set_ylim(0,50)
ax.set_ylabel('aggregate intensity : total cell intensity')
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles=handles[0:], labels=labels[0:])
ax.tick_params(axis='x', rotation=90)
plt.tight_layout()
plt.savefig(f'{output_folder}mean_ratio_in_agg_v_aggonly.png')
plt.savefig(f'{output_folder}mean_ratio_in_agg_v_aggonly.svg')




# Plot the aggregate intensity alone for the same aggregate-positive subset.
Fig, ax = plt.subplots()
sns.boxplot(data=ratios_aggs, x='treat', y='agg_intensity', palette='YlOrRd', ax=ax)
#ax.set_ylim(0,1.5)
ax.set_ylabel('aggregate intensity')
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles=handles[0:], labels=labels[0:])
ax.tick_params(axis='x', rotation=90)
plt.tight_layout()
plt.savefig(f'{output_folder}mean_agg_intensity.png')
plt.savefig(f'{output_folder}mean_agg_intensity.svg')
