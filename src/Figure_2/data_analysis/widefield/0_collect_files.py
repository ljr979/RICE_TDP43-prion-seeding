"""Collate CellProfiler exports from widefield experiments into treatment tables.

This script walks a replicate-specific folder of CellProfiler CSV exports,
attaches sample metadata from a plate map, and writes combined Cell and Aggs
tables that can be used by downstream widefield analysis scripts.
"""

"""Collate CellProfiler widefield exports and attach plate-map annotations.

Reads per-image Cell and aggregate CSV exports from a replicate folder and
writes combined tables for subsequent widefield analysis steps.
"""

import os
import pandas as pd
from src.paperfig_utils import create_sample_dictionary

Expnum = 'Expnum-replicate'
input_folder = "input_data/Figure_2/data_analysis/0_collect_files/replicate/"
output_folder = f'output/Figure_2/data_analysis/0_collect_files/{Expnum}/'

if not os.path.exists(output_folder):
   os.makedirs(output_folder)

# CellProfiler exports use consistent suffixes that identify cell and aggregate tables. The script collects each type separately and preserves them as two combined outputs for downstream steps.
namelist = ['_Cell.', '_Aggs.']
filepaths = [[f'{root}/{name}' for name in files if '.csv' in name] for root, dirs, files in os.walk(f'{input_folder}/')]
filepaths = [item for sublist in filepaths for item in sublist]


for name in namelist:
    print(name)
    xlist = [item for item in filepaths if name in item]
    print(xlist)
    name = name.split('.')[0].split('_')[1]

    collated = []
    for files in xlist:
        files = files.replace('\\', '/')
        files = files.replace('//', '/')

        # Read the plate map once per file so CellProfiler coordinates can be mapped back to the original sample names used in the experiment.
        plate_map = f'{input_folder}/plate_map.csv'
        locs = pd.read_csv(plate_map)
        locs.rename(columns={'Unnamed: 0':'row'}, inplace=True)
        params = pd.read_csv(f'{files}')
        params['exp_num'] = Expnum
        params['tp'] = files.split('/')[-1].split('_')[1]

        # Some exports store folder metadata in the image path rather than a dedicated field, so recover the last two path components here.
        params[['Metadata_Folder1', 'Metadata_Folder2']] = params['PathName_dNLS'].str.split('\\', expand=True).iloc[:, (len(params['PathName_dNLS'].str.split('\\', expand=True).columns.tolist())-2):len(params['PathName_dNLS'].str.split('\\', expand=True).columns.tolist())]

        params['fold1_fold2'] = params['Metadata_Folder1'].astype(str) + '_' + params['Metadata_Folder2'].astype(str)
        col_name_dict = create_sample_dictionary(map_df=locs)
        params['treat'] = params['fold1_fold2'].map(col_name_dict)
        
        if name == 'Cell':
            params['Unique_parent'] = params['exp_num'] + '_'+params['tp']+ '_' + params['Metadata_Folder1'].astype(str) + '_' + params['Metadata_Folder2'].astype(str) + '_' + params['ImageNumber'].astype(str)+'_'+params['ObjectNumber'].astype(str)
        if name != 'Cell':
            w_parent = params[params['Parent_Cell']!=0]
            w_parent['Unique_parent'] = params['exp_num'] + '_'+params['tp']+ '_' + params['Metadata_Folder1'].astype(str) + '_' + params['Metadata_Folder2'].astype(str) + '_' + params['ImageNumber'].astype(str)+'_'+params['Parent_Cell'].astype(str)
            params = w_parent
        collated.append(params)     

    collated = pd.concat(collated)

    print(len(collated))
    collated.to_csv(f'{output_folder}{name}_collated.csv')

