#little script to overwrite the raw paths in these dataframes, for publication.
import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger

input_folder = 'other_analysis_resources/cellprofiler_pipelines/'
output_folder = 'other_analysis_resources/cellprofiler_pipelines/'
if not os.path.exists(output_folder):
   os.makedirs(output_folder)


new_rc_params = {'text.usetex': False,
'svg.fonttype': 'none'
}
plt.rcParams.update(new_rc_params)

cp_files = [[f'{root}/{name}' for name in files if '.csv' in name] for root, dirs, files in os.walk(f'{input_folder}/')]
cp_files=[item for sublist in cp_files for item in sublist]

cols_to_remove = ['_FileLocation', 'PathName', 'FileName' ]

for count, filename in enumerate(cp_files):

   df = pd.read_csv(f'{filename}')
   counter = count + 1
   df = df.loc[:, ~df.columns.str.contains('|'.join(cols_to_remove))]
   #df.drop([col for col in df.columns.tolist() if col.str.contains()])
   print(counter, filename, df.columns.tolist() )
   df.to_csv(f'{filename}')
