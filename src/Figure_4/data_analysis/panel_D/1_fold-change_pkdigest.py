#CALCULATE the fold-change in a single replicate (not collated) in percentage of cells with aggregates, AFTER PK digestion.
#these are collated after all replicates and plotted for the figure in panel D
#prior to this step, the same 'percentage colocalisation' scripts were run as in all other figures (pipeline included for Figure 2 data, but performed the same with varied plate maps/treatments as the first step for all other figures). This script picks up from the 'count_melt' file. 
import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger

exp_num = 'Expnum-replicate'

input_folder = f'input_data/Figure_4/data_analysis/'
output_folder = f'output/Figure_4/D-E_digested_fibs_seeding/data_analysis/'

if not os.path.exists(output_folder):
    os.makedirs(output_folder)
dataframe = pd.read_csv(f'{input_folder}{exp_num}_count_melt.csv')

#fold increase
dataframe[['digestion', 'protein']] = dataframe['t3'].str.split('-', expand=True)

foldchange = []
for protein, df in dataframe[dataframe['variable']=='perc_pos'].groupby('protein'):
    protein
    df
    #find the average value of the undigested sample
    udig = df['value'][df['digestion']=='undigested'].values.mean()
    for well, values in df[df['digestion']=='PKdigested'].groupby('well'):
        well
        values

        values['foldchange'] = values['value']/udig
        foldchange.append(values)

foldchange = pd.concat(foldchange)
foldchange.to_csv(f'{output_folder}fold_change_aggcount.csv')

Fig, ax = plt.subplots()
sns.barplot(data=foldchange, x='protein', y='foldchange', palette='Blues', order=['WT', 'A382P', 'G335D', 'W385G'],capsize=0.1,linewidth=2, edgecolor='black', ax=ax)
sns.stripplot(data=foldchange, x='protein', y='foldchange', palette='Blues',order=['WT', 'A382P', 'G335D', 'W385G'], ax=ax, edgecolor='black', alpha=0.8, linewidth=1)
#ax.set_ylim(0,1)
#ax.set_ylim(0,15)
ax.set_ylabel('Fold change in aggregation after digestion')
ax.tick_params(axis='x', rotation=90)

handles, labels = ax.get_legend_handles_labels()
ax.legend(handles=handles[0:5], labels=labels[0:5])
plt.tight_layout()
plt.savefig(f'{output_folder}{tp}_fold-change_digestion.png')
plt.savefig(f'{output_folder}{tp}_fold-change_digestion.svg')