import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from statannotations.Annotator import Annotator
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from scipy import stats 
import matplotlib.patches as mpatches
from src.paperfig_utils import save_figure_svg_png, set_params_for_plotting

def bardotplot(data, xcol, ycol, order, hue=None, hue_order=None, scat_hue=None, scat_hue_order=None, palette=False, xlabel='', ylabel=False, pairs=False, correction=None, xticks=None, groups=None, group_label_y=-0.18, group_line_y=-0.05, ax=None, legend='', dot_size=5, cap_size=0.2, cap_width=2):
    """plots a BARPLOT with a scatterplot overlaid

    Args:
        data (df): dataframe containing information of interest
        xcol (str): column name for xaxis
        ycol (str): column name for y axis
        order (list): _description_
        hue (_type_, optional): if hue is being used to get split bars, pass in a string for the column name here. Defaults to None.
        hue_order (_type_, optional): should pass a list to this to ensure that the hue splitting is plotted in the correct order. Defaults to None.
        scat_hue (_type_, optional): _description_. Defaults to None.
        scat_hue_order (_type_, optional): if using hue, need to also have hue for scatterplot and pass the order of the hue so they match. Defaults to None.
        palette (bool, optional): _description_. Defaults to False.
        xlabel (str, optional): _description_. Defaults to ''.
        ylabel (bool, optional): _description_. Defaults to False.
        pairs (bool, optional): _description_. Defaults to False.
        correction (_type_, optional): _description_. Defaults to None.
        xticks (_type_, optional): _description_. Defaults to None.
        groups (_type_, optional): _description_. Defaults to None.
        group_label_y (float, optional): _description_. Defaults to -0.18.
        group_line_y (float, optional): _description_. Defaults to -0.05.
        ax (_type_, optional): _description_. Defaults to None.
        legend (str, optional): _description_. Defaults to ''.
        dot_size (int, optional): _description_. Defaults to 5.
        cap_size (float, optional): _description_. Defaults to 0.2.
        cap_width (int, optional): _description_. Defaults to 2.

    Returns:
        _type_: a figure (plot)
    """
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

def anova_tukeyshsd(df, value_col, grouping_col):
    """run a one-way anova with Tukey's HSD test afterward. Print out results 

    Args:
        df (_type_): _description_
        value_col (_type_): _description_
        grouping_col (_type_): _description_
    """
    groups = [group[value_col].values for name, group in df.groupby(grouping_col)]
    # Run one-way ANOVA
    f_stat, p_value = stats.f_oneway(*groups)

    print(f"F-statistic: {f_stat:.3f}")
    print(f"P-value: {p_value:.5f}")


    # Run Tukey's HSD test
    tukey = pairwise_tukeyhsd(endog=df[value_col], groups=df[grouping_col], alpha=0.05)

    # Print summary
    print(tukey)

def set_custom_mutant_palette(xaxis_order, t1_order=['whole', 'insoluble', 'soluble']):
    """sets palette according to the PROTEIN column, to be passed to barplot function when the hue variable is being used. sets a different colour for each x axis value, one colour from a different palette each.

    Args:
        xaxis_order (_type_): _description_
        t1_order (list, optional): _description_. Defaults to ['whole', 'insoluble', 'soluble'].

    Returns:
        _type_: _description_
    """
    palette_dict = { 'WT' : 'Purples', 'A382P':'Greens', 'G335D': 'Reds', 'W385G':'Blues', 'x': 'Oranges'
    }

    t1_order = ['whole', 'insoluble', 'soluble']


    # This will hold the full palette
    combined_palette = {}
    # For each protein, assign a base palette
    for protein in xaxis_order:
        base_colors = sns.color_palette(palette_dict[protein], n_colors=len(t1_order))
        for i, treatment in enumerate(t1_order):
            combined_palette[(protein, treatment)] = base_colors[i]
    return combined_palette

if __name__ == '__main__':


    input_folder = 'input_data/Figure_3/'

    output_folder = 'output/Figure_3/H_variants-solinsol_seeding/'

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    set_params_for_plotting(fontsize=10)
    #read in data
    collated_solinsol_seeding = pd.read_csv(f'{input_folder}Fig3H_variant_sol_insol_seeded.csv')
    #drop unnecessary columns
    collated_solinsol_seeding.drop([col for col in collated_solinsol_seeding.columns.tolist() if 'Unnamed:' in col], axis=1, inplace=True)

    #set new purple / blue palette
    hex_colors = ["#4a006f","#460677","#45147e","#421e86","#3f288d","#3d3195","#3a3b9c","#3745a4","#354fab","#3259b3"]
    custom_palette = sns.color_palette(hex_colors)
    sns.palplot(custom_palette)
    #set order of plotting (x ax)
    xaxis_order = ['WT', 'G335D', 'A382P', 'W385G']

    #plot seeding of variants when particle numbers have been normalised to the wild-type
    fig, ax = bardotplot(collated_solinsol_seeding, 'protein', 'mean', xaxis_order, palette=custom_palette,  hue=
    'plot_treat', hue_order=['whole', 'insoluble', 'soluble'], xlabel='', scat_hue='plot_treat', scat_hue_order=['whole', 'insoluble', 'soluble'],ylabel='Seeding (%) relative to whole fraction', pairs=False, correction=None, xticks=None, groups=None, group_label_y=-0.18, group_line_y=-0.05, ax=None, legend='', dot_size=5, cap_size=0.15, cap_width=2)
    #set palette according to paper standard
    combined_palette = set_custom_mutant_palette(xaxis_order,t1_order=['whole', 'insoluble', 'soluble'])
    bar_colors = [combined_palette[(protein, 'insoluble')] for protein in xaxis_order]
    newcol_dict = dict(zip(xaxis_order, bar_colors))
    #set bar faces: patterns according to fraction being represented, colour according to palette.
    for stype, patch in zip([('whole', mut) for mut in xaxis_order]+ [('insoluble', mut) for mut in xaxis_order]+[('soluble', mut) for mut in xaxis_order], ax.patches):
        if stype[0] =='insoluble':
            patch.set_facecolor(newcol_dict[stype[1]])
            patch.set_edgecolor('#fff')
            patch.set_linewidth(1)
            patch.set_hatch('///')

        if stype[0] =='soluble':
            patch.set_facecolor(newcol_dict[stype[1]])
            patch.set_edgecolor('#fff')
            patch.set_linewidth(1)
            patch.set_hatch('o')

        else:
            patch.set_facecolor(newcol_dict[stype[1]])
            patch.set_edgecolor('#fff')
            patch.set_linewidth(1)


    #create legends to match colour/pattern scheme
    mon_patch = mpatches.Patch(facecolor='dimgray', label='whole', edgecolor='#fff', linewidth=1)
    seed_patch = mpatches.Patch(facecolor='dimgray', label='insoluble', edgecolor='#fff', linewidth=0, hatch='///')
    sol_patch = mpatches.Patch(facecolor='dimgray', label='soluble', edgecolor='#fff', linewidth=0, hatch='o')

    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(reversed(labels), reversed(handles)))
    ax.legend([mon_patch, seed_patch, sol_patch], ['whole','insoluble', 'soluble'], frameon=False, loc='upper left', handlelength=2., handletextpad=0.25, bbox_to_anchor=(-0.005, 1.15), ncols=3)
    plt.tight_layout()
    #save figure
    save_figure_svg_png(fig, f'{output_folder}3_sup_seeding_mutants', width=3.5)

    #anova_tukeyshsd(df=collated_solinsol_seeding, grouping_col='plot_treat', value_col='mean')
                                        

                                        

