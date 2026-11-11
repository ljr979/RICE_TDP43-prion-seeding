#handy script if you need to combine separate treatments whose output is the same, to compare between.
#this script takes input as confocal colocalisation cellprofiler output, combines them and saves as input for the next script in this pipeline.
#note, this would only be for one replicate. The final figure takes the collated replicates from multiple replicates of this pipeline.
#note that input data for this analysis lives in 'other_analysis_resources/cellprofiler_pipelines/colocalisation/cellprofiler_output/'. The actual input for the corresponding figure can be found in input_data/Figure_4/
import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

Expnum = 'Expnum-replicate'
input_folder = f"other_analysis_resources/cellprofiler_pipelines/colocalisation/"

output_folder = f'output/Figure_4/B_fibrils_dNLS_coloc/data_analysis/0_collect_files/'

if not os.path.exists(output_folder):
   os.makedirs(output_folder)

#define the filenames we are interested in
namelist = ['Greenaggs.csv', 'LCD_aggs.csv', 'cytoplasm.csv']
#first want to collate the csv files that cell profiler output 
filepaths = [[f'{root}/{name}' for name in files if '.csv' in name] for root, dirs, files in os.walk(f'{input_folder}/')]

filepaths = [item for sublist in filepaths for item in sublist]

collated = []
for name in namelist:
#flatten list
    xlist = [item for item in filepaths if name in item]
    
    for files in xlist:
        files = files.replace('\\', '/')
        files = files.replace('//', '/')
        
        params = pd.read_csv(f'{files}')
       
        params['exp_num'] = Expnum
        params['type'] = f'{name}'.split('.csv')[0]
        params['protein'] = files.split('//')[-1].split('_')[-1].split(f'{name}')[0]

        collated.append(params)
        
collated = pd.concat(collated)
#save for next step
collated.to_csv(f'{output_folder}0_collated_LCD_dNLS_cyto-aggs.csv')


