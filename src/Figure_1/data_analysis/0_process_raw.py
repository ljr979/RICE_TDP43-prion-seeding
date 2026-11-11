"""Process raw plate-reader fluorescence data for downstream seeding analysis.

This script reads manually prepared raw assay exports and their matching plate
maps, converts the wide plate-reader layout into treatment-resolved tables, and
derives normalized fluorescence traces suitable for visualization and kinetic
curve fitting. It can optionally average technical replicates, subtract blank
signals, and generate exploratory time-course plots for quality control.

The workflow assumes that the raw plate-reader CSV and plate map have already
been standardized manually: data should begin in the top-left of the sheet,
timepoints should be stored as column headers, and the plate map should only
contain wells for which measurements are present.
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


def read_add_time(input_folder, filename, interval_mins):
    """Read a raw plate-reader export and convert it into a time-indexed table.

    Args:
        input_folder (str): Directory containing the input CSV file.
        filename (str): Base name of the raw plate-reader file, without suffix.
        interval_mins (int): Number of minutes between successive measurements.

    Returns:
        pandas.DataFrame: Transposed fluorescence data with one row per timepoint
        and one column per measured well.
    """
    df = pd.read_csv(f'{input_folder}{filename}.csv' )
    # Remove instrument-generated identifier columns that are not part of the
    # numeric time-series measurements.
    df.drop([col for col in df.columns.tolist() if 'Sample' in col], axis=1, inplace=True)
    df.replace(r'\N', np.nan, inplace = True)
    df_transpose = df.T.reset_index().rename(columns={'index':'time(h)'})
    data_only = df_transpose[~df_transpose['time(h)'].isin(['row', 'column'])]
    data_only.columns = ['time(h)'] + list(df['row'].astype(str) + '_' + df['column'].astype(str))
    data_only.insert(0, 'time(min)', range(0, 0 + len(data_only)*interval_mins,interval_mins))
    data_only['time(h)'] = data_only['time(min)']/60
    data_only.drop([col for col in data_only.columns.tolist() if 'time(h)' in col], axis=1, inplace=True)

    return data_only

def create_sample_dictionary(map_df):
    """Map plate coordinates to sample identifiers from a plate map.

    Args:
        map_df (pandas.DataFrame): Plate map with one row per plate row and one
        column per plate column.

    Returns:
        dict: Dictionary mapping keys such as ``A_1`` to user-defined sample
        identifiers.
    """
    sample_dict = {}

    for index, row in map_df.iterrows():
        row_label = row['row']

        for col in map_df.columns:
            if col == 'row':
                continue
            col
            if pd.notna(row[col]):
                coordinate = f"{row_label}_{col}"
                sample_dict[coordinate] = row[col]

    return sample_dict

def rename_cols(df, col_name_dict):
    """Rename well-coordinate columns to sample identifiers.

    Args:
        df (pandas.DataFrame): Dataframe returned by ``read_add_time``.
        col_name_dict (dict): Mapping from plate coordinates to sample names.

    Returns:
        pandas.DataFrame: Transposed dataframe with one row per sample and an
        ``ID`` column containing the assigned sample identifier.
    """
    df.rename(columns=col_name_dict, inplace=True)
    df.drop([col for col in df.columns.tolist() if col in ['Well', 'Unnamed: 1', 'time(min)']], axis=1, inplace=True)
    df = df.T.reset_index().rename(columns={'index':'ID'})
    return df

def add_info (df):
    """Split composite sample identifiers into treatment metadata columns.

    Args:
        df (pandas.DataFrame): Sample-level dataframe whose ``ID`` column stores
        composite treatment labels separated by underscores.

    Returns:
        pandas.DataFrame: Dataframe augmented with treatment and replicate
        annotation columns.
    """
    new=[]
    for ids, d in df.groupby('ID'):    
        ids
        d
        if 'time(h)' not in ids:
            d['treat2'] = ids.split('_')[1]
            d['treat1'] = ids.split('_')[0]
            #d['treat3'] = ids.split('-')[2]

            d['replicate'] = d.groupby('ID').cumcount()+ 1
            new.append(d)
    new = pd.concat(new)
    return new

def average_tech_replicates(new):
    """Average numeric measurements across technical replicates.

    Args:
        new (pandas.DataFrame): Dataframe containing replicate-level sample
        measurements and metadata.

    Returns:
        pandas.DataFrame: Mean fluorescence trace for each treatment group.
    """
    averages = []
    for treat, x in new.groupby('treat2'):
        x
        to_av = [col for col in x if col not in ['ID', 'replicate', 'treat1', 'treat2']]
        numeric = x[[c for c in x.columns if c in to_av]].astype(int)
        average = pd.DataFrame(numeric.mean(axis=0)).T
        average['treat'] = treat
        averages.append(average)
    averages = pd.concat(averages)
    return averages

def subtract_background(averages, blank_id):
    """Subtract the mean blank signal from all non-blank traces.

    Args:
        averages (pandas.DataFrame): Processed fluorescence traces.
        blank_id (str): Identifier corresponding to the blank control wells.

    Returns:
        tuple: Background-subtracted dataframe, list of retained treatment IDs,
        and a transposed version of the subtracted data.
    """
    blanks = averages[averages['ID'] == blank_id]
    numeric_cols = [col for col in averages.columns if isinstance(col, (int, float))]
    string_cols = [col for col in averages.columns if col not in numeric_cols]
    blank_av = blanks[numeric_cols].mean()

    treats = averages[averages['ID'] != blank_id].copy()

    treats['replicate'] = treats['replicate'].astype(str)
    treats['ID'] = treats['ID'] + '_' + treats['replicate']

    subtracted = treats[numeric_cols].sub(blank_av, axis=1)

    subtracted = pd.concat([treats[string_cols], subtracted], axis=1)
    subtracted['type'] = 'blank_subtracted'

    subtracted_T = subtracted.set_index('ID').T

    return subtracted, treats['ID'].tolist(), subtracted_T

def plot_fluo_over_time(d, interval_mins, output_folder, palette, Exp_num, length_exp, types=['change_in_intensity'], ylim=False, savefig=True):
    """Plot fluorescence traces over time for selected treatment groups.

    Args:
        d (pandas.DataFrame): Wide-format processed fluorescence table.
        interval_mins (int): Time between consecutive measurements in minutes.
        output_folder (str): Directory used for figure export.
        palette (str | dict): Seaborn palette or palette mapping.
        Exp_num (str): Experiment identifier for output filenames.
        length_exp (float): Upper limit for the x-axis in hours.
        types (list, optional): Data types to include in the plot.
        ylim (tuple | bool, optional): Optional y-axis limits.
        savefig (bool, optional): If ``True``, save the generated figures.

    Returns:
        pandas.DataFrame: Long-format dataframe used to generate the plots.
    """

    test_Df = pd.melt(d, id_vars= ['ID', 'type', 't1',  't2', 't3', 'replicate'], var_name='timepoint', value_name='intensity' )

    test_Df['time (h)'] = (test_Df['timepoint'].astype(int)*interval_mins)/60 
    if ylim == False:
        ylim = (min(test_Df['intensity']), max(test_Df['intensity']))

    for t1, xyz in test_Df[test_Df['type'].isin(types)].groupby(['t1']):
        print(t1)
        t1 = t1[0]
        sub = len(xyz['t2'].unique().tolist())
        if sub == 1:
            row = 1
            col = 2
        if sub > 1: 
            row = 1
            col = sub
        Fig, axes = plt.subplots(row, col, figsize=(16,5))
        if isinstance(palette, dict):
            new_pal = palette[t1]
        else:
            new_pal = palette
        for i, (t2, abc) in enumerate(xyz.groupby('t2')):
            i
            t2
            abc
            sns.lineplot(data=abc, x='time (h)', y='intensity', palette=new_pal, hue='t3', ax=axes[i])
            axes[i].set_xlim(0,length_exp)
            axes[i].set_ylim(ylim)
            axes[i].annotate(f'{t2}', xy=(0.8, 1))
            axes[i].legend()
        Fig.suptitle(f'{t1}')        
        if savefig == True:
            plt.savefig(f'{output_folder}{Exp_num}_{t1}timeseries.png')
            plt.savefig(f'{output_folder}{Exp_num}_{t1}timeseries.svg')
            plt.show()
    return test_Df

def change_over_time_normalise(treats):
    """Calculate baseline-subtracted and normalized fluorescence traces.

    Args:
        treats (pandas.DataFrame): Dataframe containing metadata columns and
        numeric timepoint columns.

    Returns:
        pandas.DataFrame: Concatenated dataframe containing both absolute change
        from baseline and values normalized to each trace maximum.
    """
    numeric_cols = [col for col in treats.columns if isinstance(col, (int, float))]
    string_cols = [col for col in treats.columns if col not in numeric_cols]
    change_in_intensity = treats[numeric_cols].sub(treats[numeric_cols].iloc[:, 0], axis=0)
    max_change = change_in_intensity.max(axis=1)

    prop_of_total = change_in_intensity.div(max_change.where(max_change > 0, 1), axis=0)
    prop_of_total.sub(prop_of_total.iloc[:, 0], axis=0)

    change_in_intensity = pd.concat([treats[string_cols], change_in_intensity], axis=1)

    change_in_intensity['type'] = 'change_in_intensity'
    
    prop_of_total = pd.concat([treats[string_cols], prop_of_total], axis=1)
    prop_of_total['type'] = 'normalised_to_max_change'

    return pd.concat([change_in_intensity, prop_of_total])

def get_more_info(df, col_to_split, what_to_split):
    """Split a composite column into ``t1``, ``t2``, and ``t3`` annotations.

    Args:
        df (pandas.DataFrame): Dataframe containing a composite treatment column.
        col_to_split (str): Name of the column to split.
        what_to_split (str): Delimiter used to separate treatment components.

    Returns:
        pandas.DataFrame: Input dataframe with added treatment annotation
        columns.
    """
    df[['t1', 't2', 't3', ]] = df[f'{col_to_split}'].str.split(f'{what_to_split}',expand=True)
    return df

if __name__ == '__main__':


    replicate = 'rep_1'
    Exp_num = 'Exp31-1'
    input_folder = 'input_data/Figure_1/data_analysis/0_process_raw/'
    output_folder = 'output/Figure_1/data_analysis/0_process_raw/'

    # Define the manually curated plate map that maps well positions to sample
    # identities.
    plate_map = f'{input_folder}{replicate}_plate_map.csv'

    filename = f'{replicate}_raw_data'

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    interval_mins = 5
    bg_subtract = True
    blank_id = 'control_blank-x-x'
    df = read_add_time(input_folder, filename, interval_mins)

    # Read the plate map and generate the well-to-sample lookup used to label
    # the raw fluorescence traces.
    locs = pd.read_csv(plate_map)
    locs.rename(columns={'Unnamed: 0':'row'}, inplace=True)
    col_name_dict = create_sample_dictionary(map_df=locs)
    df = rename_cols(df, col_name_dict)
    df = df.fillna(0)
    new = add_info(df)
    new.drop(new.index[:1], inplace=True)

    # Optional averaging collapses technical replicates before downstream
    # normalization and plotting.
    average = False
    if average:
        averages = average_tech_replicates(new)
    if not average:
        averages = new

    bg_subtract = False
    if bg_subtract:
        subtracted, treats_list, treats_transp = subtract_background(averages, blank_id=blank_id)
        subtracted.to_csv(f'{output_folder}1_bg_subtracted.csv')
        alls = change_over_time_normalise(treats=subtracted)

    if not bg_subtract:
        treats = averages 
        treats['type'] = 'no_subtract'
        treats.to_csv(f'{output_folder}1_processed_raw.csv')

        # Reshape the processed traces into a form suitable for normalization and
        # later plotting.
        selected_cols = ['ID']+[col for col in averages.columns if isinstance(col, (int, float))]
        treats_transp = treats[selected_cols]
        treats_transp = treats_transp.T
        header = treats_transp.iloc[0]
        treats_transp = treats_transp[1:].rename(columns=header)
        treats_transp = treats_transp.apply(pd.to_numeric)
        treats_transp = df.T.reset_index().rename(columns = {'index':'ID'})

        alls = change_over_time_normalise(treats)

    get_more_info(df=alls, col_to_split='treat2', what_to_split='-')
    alls.to_csv(f'{output_folder}1_change_in_intensity_plus_info.csv')

    x = alls
    controls = x[x['treat1'].isin(['control'])]
    tests = x[x['treat1'].isin(['test'])]

    tests = tests[tests['type']=='normalised_to_max_change']


    test_Df = pd.melt(tests, id_vars= ['ID', 'treat1',  'treat2','type','replicate', 't1', 't2', 't3'], var_name='timepoint', value_name='intensity' )
    test_Df['time (h)'] = (test_Df['timepoint'].astype(int)*interval_mins)/60

    palette_dict = {

        'WT' : 'Purples', 
        'A382P':'Greens', 
        'G335D': 'Reds', 
        'W385G':'Blues', 
        'x': 'coolwarm'

    }

    test_Df = get_more_info(test_Df, col_to_split='treat2', what_to_split='-')

    # Generate per-mutation quality-control plots of normalized fluorescence
    # trajectories for visual inspection.
    for mut, df in test_Df.groupby('t2'):
        df
        print(mut)
        
        new_pal = palette_dict[mut]
        Fig, ax = plt.subplots()
        lines = sns.lineplot(data=df, x='time (h)', y='intensity', palette=new_pal, hue='t1', ax=ax)
        ax.set_title(f'{mut}')
        ax.set_ylim(-0.05, 1)
        plt.tight_layout() 
        plt.savefig(f'{output_folder}{mut}_change-in-intensity.png')
