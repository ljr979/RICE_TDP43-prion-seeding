import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger
from statannotations.Annotator import Annotator
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from scipy import stats
from matplotlib.colors import LogNorm

def bardotplot(data, xcol, ycol, order, hue=None, hue_order=None, scat_hue=None, scat_hue_order=None, palette=False, xlabel='', ylabel=False, pairs=False, correction=None, xticks=None, groups=None, group_label_y=-0.18, group_line_y=-0.05, ax=None, legend='', dot_size=5, cap_size=0.2, cap_width=2):
    if ax is None:
        fig, ax = plt.subplots()
    sns.barplot(
        data=data,
        x=xcol,
        y=ycol,
        hue=hue,
        palette=palette,
        capsize=cap_size,
        errwidth=cap_width,
        errcolor='black',
        ax=ax,
        dodge=True,
        order=order,
        hue_order=hue_order,
        edgecolor='black'
    )
    sns.stripplot(
        data=data,
        x=xcol,
        y=ycol,
        hue=scat_hue,
        color='white',
        ax=ax,
        edgecolor='black',
        linewidth=1,
        s=dot_size,
        order=order,
        hue_order=scat_hue_order,
        dodge=True,
    )

    if pairs:
        annotator = Annotator(
            ax=ax, pairs=pairs, data=data, x=xcol, y=ycol, order=order, hue=hue, hue_order=hue_order)
        annotator.configure(test='t-test_ind', text_format='star',
                            loc='inside', comparisons_correction=correction, line_width=0.5)
        annotator.apply_and_annotate()

    ax.set(ylabel=ylabel)

    ax.set_xlabel(xlabel)
    if xticks:
        ax.set_xticks(xticks)
        ax.set_xticklabels(hue_order*len(order))
    if groups:
        for group_label, (x0, x1, x2) in groups.items():
            ax.annotate(group_label, xy=(x0, group_label_y),
                        xycoords='data', ha='center', annotation_clip=False)
            trans = ax.get_xaxis_transform()
            ax.plot([x1, x2], [group_line_y, group_line_y],
                    color="black", transform=trans, clip_on=False)

    if legend == '':
        ax.legend('', frameon=False)
    else:    
        handles, labels = plt.gca().get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        ax.legend(by_label.values(), by_label.keys())
   
    return fig, ax

def save_figure_svg_png(fig, filename, width=3.5, height=None, dpi=300, tight=True):
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

if __name__ == '__main__':
    
    input_folder = 'input_data/Figure_3/data_analysis/'
    output_folder = 'output/Figure_3/D_E_F_particles_titration/data_analysis/'

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    input_data = pd.read_csv(f'{input_folder}1_collated_aggregate_counts.csv')
    input_data.drop([col for col in input_data.columns.tolist() if 'Unnamed:' in col], axis=1, inplace=True)
    #find averages of each treatment across wells
    averages = pd.DataFrame(input_data.groupby(['Experiment', 'treat'])['value'].mean()).reset_index()
    #split treatment data up
    averages[['type', 'original_batch', 'concentration']] = averages['treat'].str.split('_', expand=True)

    #first subset for the normalisation data with their original batches too.
    norm_data = averages[averages['type'].isin(['base', 'normed'])]
    titrations = averages[averages['type'].isin(['base', 'test'])]

    #make a sensible naming convention for when plotting
    grouping_dict = {
            'batch1':'Seeds 1', 
            'batch2' : 'Seeds 2',
            'batch3' : 'Seeds 3'
    }


    #group the batchs and their normalisation together.
    norm_data['norm_group'] = 'filler'
    for count,batch in enumerate(['batch1', 'batch2', 'batch3']):
        batch
        count = count+1
        norm_data['norm_group'][(norm_data['type']=='base')&(norm_data['original_batch']==f'{batch}')] = grouping_dict[batch]
        norm_data['norm_group'][(norm_data['type']=='normed')&(norm_data['concentration']==f'to{batch}')] = grouping_dict[batch]


    #for each replicate, need to calculate a) the original ratio between the batches (so, for seeds1 group, the ratio between seeds 1 and seeds 2, and seeds 1 and seeds 3). Set this as 1. Then calculate the normalised ratio between the same treatments (i.e after normalising to batch 1, did this ratio change?). then, report that as the 'correction' i.e I guess reporting the degree of change that normalisation affords? 
    #this part relates each batch to each other batch. this correlates with PANEL F of figure 3, and supplementary figure 7.
    # Collect rows as records
    records = []
    new_basedf = []

    for replicate, df in norm_data.groupby('Experiment'):
        base_df = df[df['type'] == 'base']
        
        for batch in ['batch1', 'batch2', 'batch3']:
            try:
                copy = base_df.copy()
                ref_val = copy.loc[copy['original_batch'] == batch, 'value'].values[0]
                copy['divided'] = copy['value'] / ref_val
                copy['type_difference'] = 'original_difference'
                copy['div_by_batch'] = batch
                new_basedf.append(copy)
                print(batch, ref_val, copy['divided'].values)
            
            except IndexError:
                continue  # skip if the reference batch isn't found
            

        for norm_group, norm_df in df.groupby('norm_group'):
            norm_batch_num = norm_group.split()[-1]  # get number from 'Seeds 2'
            match_batch = f'batch{norm_batch_num}'

            try:
                ref_val = norm_df.loc[norm_df['original_batch'] == match_batch, 'value'].values[0]

                norm_df['divided'] = norm_df['value']/ref_val
                norm_df['type_difference'] = 'normed_difference'
                norm_df['div_by_batch'] = match_batch
                records.append(norm_df)

            except IndexError:
                norm_diffs = [None] * len(norm_df)
            

    # Convert to DataFrame
    df_differences = pd.concat(records)
    original_diffs = pd.concat(new_basedf)

    # Preview
    print(df_differences.head())

    combined_dfs = pd.concat([df_differences, original_diffs])

    #from this point where I'm at, this current dataframe needs to pivot so that I can divide the normalised difference by the original difference. Maybe add in a column again that groups them specifically together rather than the 'source key' column then can group on the 'grouping' and divide by the original difference and plot that. 
    norm_correction = []
    for experiment, data in combined_dfs.groupby('Experiment'):
        data
        if experiment == 'Exp28-1':
            data = data[data['treat']!='normed_batch1_tobatch3']
        notself = data[data['original_batch']!=data['div_by_batch']]
        for grouping, df in notself.groupby('div_by_batch'):
            grouping
            df
            for batch in ['batch1', 'batch2', 'batch3']:
                newdf = df[df['original_batch']==batch]
                if len(newdf) > 1:
                    newdf['norm_correction'] = newdf['divided'][newdf['type_difference']=='normed_difference'].values[0]/newdf['divided'][newdf['type_difference']=='original_difference'].values[0]
                    norm_correction.append(newdf)
                    newdf['plot_group'] = newdf['original_batch']+'_'+newdf['div_by_batch']
                else:
                    print(batch, grouping, experiment, f'length of dataframe is 0')

    norm_correction = pd.concat(norm_correction)
    norm_correction.to_csv(f'{output_folder}collated_normalised_degree_correction.csv')

#------------------------------------plotting below
    #Now the data is all organised, we can plot the titration alone (i.e dilution of particles = less seeding. Figure 3 E)
    #now here we want to be plotting the titration data. 
    titrations['dilution_factor'] = titrations['concentration'].astype(float)/float(titrations['concentration'][titrations['type']=='base'].unique()[0])

    dilution_impact = []
    for experiment, df in titrations.groupby('Experiment'):
        df
        for original_batch, data in df.groupby('original_batch'):
            data
            data['effect_of_dilution'] = data['value']/data['value'][data['type']=='base'].values[0]
            dilution_impact.append(data)
    dilution_impact = pd.concat(dilution_impact)
    #save!
    dilution_impact.to_csv(f'{output_folder}dilution_impact_on_seeding.csv')

    dilution_impact['dilution_factor'] = dilution_impact['dilution_factor'].astype(str)
    #and plot it here. 

    #FIGURE 3, PANEL E
    fig,ax = plt.subplots()
    sns.barplot(data=dilution_impact, x='original_batch', y='effect_of_dilution', hue='dilution_factor', palette='viridis')
    plt.savefig(f'{output_folder}E_basic_diluted_seeds.png')

    #for plotting only the normalised data (not dilutions, i.e. panel F)
    plotting_df = norm_correction[norm_correction['type']!='base']



    #------------------- FIGURE 3F: REARRANGING THE EXACT SAME DATA FROM ABOVE, BUT MAKING AVAILABLE TO PLOT as a matrix. SEE Figure_3/F_normed_particles for final plotting, and input_data/Figure_3/Fig3F_particle_num_seeding.csv for final data for figure.
    #Assuming a tidy DataFrame df with columns condition ("original"/"equalised" ('type_difference')), ref ('div_by_batch'), seed ('original batch'), rep('Experiment'), value('divided'):
    summ = norm_correction.groupby(["type_difference", "div_by_batch", "original_batch"])["divided"].agg(["mean", "std"]).reset_index()
    batch_map = {'batch1':1, 'batch2':2, 'batch3':3}
    summ['div_by_batch'] = summ['div_by_batch'].map(batch_map)
    summ['original_batch'] = summ['original_batch'].map(batch_map)
    #save this data as the exact version we use to plot figure 3F. 
    summ.to_csv(f'input_data/Figure_3/Fig3F_matrix_seeding_normed.csv')

    #only plotting below here.
    batches = [1,2,3]
    norm = LogNorm(vmin=0.5, vmax=2)

    fig, axes = plt.subplots(1, 2, figsize=(7, 3.2), sharey=True, layout="constrained")


    for ax, cond, title in zip(axes, ["original_difference", "normed_difference"],
                            ["Original", "Particle # equalised"]):
        s = summ[summ.type_difference == cond]
        mean = s.pivot(index="div_by_batch", columns="original_batch", values="mean").reindex(index=batches, columns=batches)
        sd = s.pivot(index="div_by_batch", columns="original_batch", values="std").reindex(index=batches, columns=batches)
        annot = (mean.round(2).astype(str) + "\n± " + sd.round(2).astype(str)).where(mean.notna(), "")

        sns.heatmap(mean, ax=ax, annot=False, cmap="PRGn", norm=norm, cbar=False,
                linewidths=1, linecolor="white", square=True)
        ax.set_facecolor("#e6e6e6")   # greys out the undefined diagonal
            
        cmap = ax.collections[0].cmap
        for i, ref in enumerate(batches):
            for j, seed in enumerate(batches):
                m, sdev = mean.loc[ref, seed], sd.loc[ref, seed]
                if pd.notna(m):
                    r, g, b, _ = cmap(norm(m))
                    lum = 0.299 * r + 0.587 * g + 0.114 * b
                    ax.text(j + 0.5, i + 0.5, f"{m:.2f}\n± {sdev:.2f}",
                            ha="center", va="center", fontsize=9,
                            color="white" if lum < 0.45 else "black")

        ax.set_title(title)
        ax.set_xlabel("Seed batch")
        ax.set_ylabel("Batch normalised to" if ax is axes[0] else "")
        ax.set_xticklabels([f"{b}" for b in batches], rotation=0)
        ax.set_yticklabels([f"{b}" for b in batches], rotation=0)

    cb = fig.colorbar(axes[0].collections[0], ax=axes, shrink=0.8,
                    label="Degree of change due to normalisation")
    cb.set_ticks([0.5, 0.75, 1, 1.5, 2])
    cb.set_ticklabels(["0.5", "0.75", "1", "1.5", "2"])


    save_figure_svg_png(fig, f'{output_folder}3F_particle-norm_seeding', width=3.5, height=None, dpi=600, tight=False)


    #-----------------------supp figure correlating to this data (actual individual data points.)

    #SUPPLEMENTARY FIGURE 7 PANEL A: plot the original vs. normalised  differences
    xaxis_order = ['batch2_batch1', 'batch3_batch1', 'batch1_batch2','batch3_batch2', 'batch1_batch3','batch2_batch3' ]
    fig, ax = bardotplot(norm_correction, 'plot_group', 'divided', xaxis_order, palette='viridis',  hue='type_difference', hue_order=['original_difference', 'normed_difference'], xlabel='',scat_hue='type_difference',scat_hue_order=['original_difference', 'normed_difference'], ylabel='Degree of change due to normalisation', pairs=False, correction=None, xticks=None, groups=None, group_label_y=-0.18, group_line_y=-0.05, ax=None, legend=True, dot_size=5, cap_size=0.15, cap_width=2)
    ax.tick_params(axis='x', rotation=90)



