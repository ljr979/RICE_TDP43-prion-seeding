"""Plot ThT fluorescence time courses for Figure 1B and Supplementary Figure 1.

Reads a processed time-course table, summarizes signal by treatment and time,
and exports separate plots for each variant and for monomer/control samples.
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def save_figure_svg_png(fig, filename, width=3.5, height=None, dpi=300, tight=True):
    """Save a Matplotlib figure as SVG and PNG at the requested size."""
    """
    Save a matplotlib figure as both SVG and PNG with standardized formatting.
    
    Parameters:
    - fig: matplotlib.figure.Figure
        The figure object to save.
    - filename: str
        Path (with or without extension). Saves as .svg and .png.
    - width: float
        Width in inches. Default is 3.5" (single column).
    - height: float or None
        Height in inches. If None, uses golden ratio.
    - dpi: int
        Dots per inch for PNG (raster export).
    - tight: bool
        Whether to use tight_layout before saving.
    """
    # Set height using golden ratio if not provided
    if height is None:
        height = width * 0.618  # Golden ratio

    fig.set_size_inches(width, height)

    if tight:
        plt.tight_layout()

    # Ensure filename has no extension so we can add both .svg and .png
    base, _ = os.path.splitext(filename)

    svg_path = f"{base}.svg"
    png_path = f"{base}.png"

    # Save SVG
    fig.savefig(svg_path, format='svg', dpi=dpi)
    # Save PNG
    fig.savefig(png_path, format='png', dpi=dpi)

    print(f"Saved SVG to: {svg_path}")
    print(f"Saved PNG to: {png_path}")

def set_params_for_plotting(fontsize=10):
    """Set the plotting style used by the Figure 1 fluorescence panels."""
    new_rc_params = {'text.usetex': False,
    'font.family': 'Arial',
    'font.size': fontsize,            # Base font size
    'axes.titlesize': fontsize,
    'axes.labelsize': fontsize,
    'xtick.labelsize': fontsize,
    'ytick.labelsize': fontsize,      
    'legend.fontsize': fontsize,
    'figure.titlesize': fontsize,
    'svg.fonttype': 'none',     # Optional: keeps text as text in SVG
    }
    
    plt.rcParams.update(new_rc_params)

if __name__ == '__main__':

    input_folder = 'input_data/Figure_1/'
    output_folder = 'output/Figure_1/B_THT-assays/'

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    set_params_for_plotting(fontsize=10)

    
    # Read the processed fluorescence measurements and treatment metadata.
    data = pd.read_csv(f'{input_folder}1_change_in_intensity_plus_info.csv')
    # Drop the CSV index column if it was included during export.
    data.drop([col for col in data.columns.tolist() if 'Unnamed: 0' in col], axis=1, inplace=True)
    # Restrict the plots to absolute change in intensity.
    data = data[data['type'] == 'change_in_intensity']
    # Set the acquisition interval used to convert measurement indices to hours.
    interval_mins = 6
    # Reshape timepoint columns into one row per sample and measurement.
    test_Df = pd.melt(data, id_vars= ['ID','treat2', 'treat1','t1', 't2', 't3','type','replicate'], var_name='timepoint', value_name='intensity' )
    # Convert the measurement index to elapsed hours.
    test_Df['time (h)'] = (test_Df['timepoint'].astype(int)*interval_mins)/60
    # Define a sequential color palette for each variant.
    palette_dict = {
        'WT' : 'Purples',
        'A382P':'Greens', 
        'G335D': 'Reds', 
        'W385G':'Blues', 'x': 
        'Oranges'

    }
    # Define marker shapes for treatment classes.
    marker_dict = {
        'blank' : 'X', 
        'monomer':'o', 
        'seeding': 'D', 
        'seeds':'s', 
        'xseeding': 'X'
                   }
    
    # Calculate mean signal and standard error for each condition and time.
    summary_df = (
        test_Df
        .groupby(['time (h)', 't1', 't2'])
        .agg(mean_intensity=('intensity', 'mean'),
            sem_intensity=('intensity', 'sem'))
        .reset_index()
    )

    # Retain regularly spaced points for a less crowded plot.
    scatter_df = (
        summary_df
        .groupby(['t1', 't2'])
        .apply(lambda g: g.iloc[::6])
        .reset_index(drop=True)
    )


    # Generate the per-variant Supplementary Figure 1 plots.
    for mut, df in test_Df.groupby('t2'):
        print(mut)

        # Set up figure
        Fig, ax = plt.subplots(figsize=(8, 6))
        # Get a single color from the palette name — middle color of 6-tone palette
        single_color = sns.color_palette(palette_dict[mut], n_colors=6)[3]
        # Use only one color for all lines, but vary marker with `style`
        new_pal = {k: single_color for k in df['t1'].unique()}
        # Lineplot (faint lines)
        sns.lineplot(
            data=df,
            x='time (h)', y='intensity',
            hue='t1',
            style='t1',
            dashes=False,
            palette=new_pal,
            linewidth=0.5,
            alpha=0.8,
            legend=False,
            ax=ax
        )

        # Scatterplot (strong markers)
        sns.scatterplot(
            data=scatter_df[scatter_df['t2'] == f'{mut}'],
            x='time (h)', y='mean_intensity',
            hue='t1',
            style='t1',
            markers=marker_dict,
            palette=new_pal,
            s=15,
            edgecolor='black',
            alpha=0.6, 
            ax=ax
        )
        
        ax.set_title(f'{mut}')
        ax.set_ylabel('Change in ThT intensity (x$10^4$)')
        ax.set_xlabel('Time (h)')
        ax.set_xlim(0,15)
        ax.set_ylim(-300, 30000)
        plt.tight_layout()
        ax.get_legend().set_visible(False) 
        ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0)) 
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        save_figure_svg_png(Fig, f'{output_folder}SuppFig1_change-in-intensity_{mut}', height=1.9, width=2.1)

        # Prepare a second plot containing only blank and monomer controls.
        monomer_only = test_Df[test_Df['t1'].isin(['blank', 'monomer'])].copy()
        mon_only_palette_dict = {}
        for mut, df in monomer_only.groupby('t2'):
            print(mut)
            # Get a single color from the palette name — middle color of 6-tone palette
            single_color = sns.color_palette(palette_dict[mut], n_colors=6)[3]
            print(single_color)
            mon_only_palette_dict[mut] = single_color

        # Set up figure
        Fig, ax = plt.subplots(figsize=(8, 6))
        # Lineplot (faint lines)
        sns.lineplot(
            data=monomer_only,
            x='time (h)', y='intensity',
            hue='t2',
            style='t1',
            dashes=False,
            palette=mon_only_palette_dict,
            linewidth=0.5,
            alpha=0.8,
            legend=False,
            ax=ax
        )

        # Scatterplot (strong markers)
        sns.scatterplot(
            data=scatter_df[scatter_df['t1'].isin(['blank', 'monomer'])].copy(),
            x='time (h)', y='mean_intensity',
            hue='t2',
            style='t1',
            markers=marker_dict,
            palette=mon_only_palette_dict,
            s=15,
            edgecolor='black',
            alpha=0.6, 
            ax=ax
        )
        
        #ax.set_title(f'{mut}')
        ax.set_ylabel('Change in ThT intensity (x$10^4$)')
        ax.set_xlabel('Time (h)')
        ax.set_xlim(0,15)
        ax.set_ylim(-300, 30000)
        plt.tight_layout()
        ax.get_legend().set_visible(False) 
        ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0)) 
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        save_figure_svg_png(Fig, f'{output_folder}Fig1B_change-in-intensity_monomer-only', height=1.9, width=2.1)