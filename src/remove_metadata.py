"""Remove path and filename metadata columns from CellProfiler CSV exports.

The script recursively edits CSV files beneath
``other_analysis_resources/cellprofiler_pipelines/`` in place. Keep backups of
any files that must retain their original metadata.
"""

import os
import pandas as pd

input_folder = 'other_analysis_resources/cellprofiler_pipelines/'
output_folder = 'other_analysis_resources/cellprofiler_pipelines/'
if not os.path.exists(output_folder):
   os.makedirs(output_folder)


cp_files = [[f'{root}/{name}' for name in files if '.csv' in name] for root, dirs, files in os.walk(f'{input_folder}/')]
cp_files=[item for sublist in cp_files for item in sublist]

cols_to_remove = ['_FileLocation', 'PathName', 'FileName' ]

for count, filename in enumerate(cp_files):

   df = pd.read_csv(f'{filename}')
   counter = count + 1
   df = df.loc[:, ~df.columns.str.contains('|'.join(cols_to_remove))]
   print(counter, filename, df.columns.tolist() )
   df.to_csv(f'{filename}')
