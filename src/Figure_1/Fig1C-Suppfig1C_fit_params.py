
"""Plot fitted lag-time and elongation-rate parameters for Figure 1C.

Reads the curated parameter table, separates conditions by variant and seeding
type, saves the manuscript plots, and prints per-variant statistical summaries.
"""
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from src.paperfig_utils import save_figure_svg_png, set_params_for_plotting
import matplotlib.patches as mpatches
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from scipy import stats


def bardotplot(data, xcol, ycol, order, hue=None, hue_order=None, scat_hue=None, scat_hue_order=None, palette=False, xlabel='', ylabel=False, xticks=None, group_label_y=-0.18, group_line_y=-0.05, ax=None, legend='', dot_size=5, cap_size=0.2, cap_width=2):
    """Plot a parameter as bars with individual observations overlaid.

    Args:
        data: Table containing fit parameters and treatment metadata.
        xcol: Column used for the x-axis categories.
        ycol: Column containing the parameter to plot.
        order: Display order for x-axis categories.
        hue: Optional column used to split bars into groups.
        hue_order: Display order for hue groups.
        scat_hue: Optional column used to color individual observations.
        scat_hue_order: Display order for scatter groups.
        palette: Seaborn palette specification.
        xlabel: X-axis label.
        ylabel: Y-axis label.
        xticks: Optional x-axis tick positions.
        group_label_y: Vertical position for optional group labels.
        group_line_y: Vertical position for optional group lines.
        ax: Existing Matplotlib axes, or ``None`` to create axes.
        legend: Empty string hides the legend; another value displays it.
        dot_size: Size of the individual data points.
        cap_size: Width of bar error-bar caps.
        cap_width: Line width of bar error bars.

    Returns:
        tuple: The Matplotlib figure and axes.
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


    ax.set(ylabel=ylabel)

    ax.set_xlabel(xlabel)

    if xticks:
        ax.set_xticks(xticks)
        ax.set_xticklabels(hue_order*len(order))
    

    if legend == '':
        ax.legend('', frameon=False)
    else:    
        handles, labels = plt.gca().get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        ax.legend(by_label.values(), by_label.keys())
   
    return fig, ax

if __name__ == '__main__':

    input_folder = 'input_data/Figure_1/'
    output_folder = 'output/Figure_1/'

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)


    set_params_for_plotting(fontsize=10)

    fit_params = pd.read_csv(f'{input_folder}SuppFig1C_fitparams_elongation-rate.csv')

    palette_dict = { 

        'WT' : 'Purples', 
        'A382P':'Greens', 
        'G335D': 'Reds', 
        'W385G':'Blues', 
        'x': 'Oranges'

    }

    t1_order = ['monomer', 'seeding', 'xseeding']
    xaxis_order = ['WT', 'A382P', 'G335D', 'W385G']

    # This will hold the full palette
    combined_palette = {}
    # For each protein, assign a base palette
    for protein in xaxis_order:
        base_colors = sns.color_palette(palette_dict[protein], n_colors=len(t1_order))
        for i, treatment in enumerate(t1_order):
            combined_palette[(protein, treatment)] = base_colors[i]

    #now assign each bar a colour based on this palette and the type of seeding
    bar_colors = [combined_palette[(protein, 'seeding')] for protein in xaxis_order]
    newcol_dict = dict(zip(xaxis_order, bar_colors))

    #change the plotted y label depending on the parameter being plotted
    ylabel_dict = {
        'tlag': 'Lag Phase (h)',
        'elongation_rate': 'Elonation rate (change in ThT/h)'
    }

    for param in ['tlag', 'elongation_rate']:
        #plot elonation rate as a barplot with scatterplot overlaid
        fig, ax = bardotplot(fit_params, 't2', f'{param}', xaxis_order, hue='t1', hue_order=t1_order, scat_hue='t1', scat_hue_order=t1_order, palette='tab10', xlabel='', ylabel=ylabel_dict[f'{param}'], xticks=None, group_label_y=-0.18, group_line_y=-0.05, ax=None, legend='True', dot_size=5, cap_size=0.15, cap_width=2)

        #assign different patches to each bar depending on the type of seeding occurring
        for stype, patch in zip([('monomer', mut) for mut in xaxis_order]+ [('seeding', mut) for mut in xaxis_order]+[('xseeding', mut) for mut in xaxis_order], ax.patches):
            if stype[0] =='seeding':
                patch.set_facecolor(newcol_dict[stype[1]])
                patch.set_edgecolor('#fff')
                patch.set_linewidth(1)
                patch.set_hatch('///')
            elif stype[0] == 'xseeding':
                patch.set_facecolor(newcol_dict[stype[1]])
                patch.set_edgecolor('#fff')
                patch.set_linewidth(1)
                patch.set_hatch('o')
            else:
                patch.set_facecolor(newcol_dict[stype[1]])
                patch.set_edgecolor('#fff')
                patch.set_linewidth(1)
            
            #add legend entries for each of the patches defined above
            mon_patch = mpatches.Patch(facecolor='dimgray', label='Monomer', edgecolor='#fff', linewidth=1)
            seed_patch = mpatches.Patch(facecolor='dimgray', label='Self-seed', edgecolor='#fff', linewidth=0, hatch='///')
            xseed_patch = mpatches.Patch(facecolor='dimgray', label='Cross-seed', edgecolor='#fff', linewidth=0, hatch='o')
            #apply patches to the legend
            handles, labels = ax.get_legend_handles_labels()
            by_label = dict(zip(reversed(labels), reversed(handles)))
            ax.legend([mon_patch, seed_patch, xseed_patch], ['Monomer', 'Self-seed', 'Cross-seed'], frameon=False, loc='upper left', handlelength=2., handletextpad=0.25, bbox_to_anchor=(-0.005, 1.15), ncols=3)

        plt.tight_layout()
        if param=='elongation_rate':
            save_figure_svg_png(fig, f'{output_folder}SuppFig1C_elongation', width=3.5, dpi=600)

        elif param=='tlag':
            save_figure_svg_png(fig, f'{output_folder}Fig1C_lagphase', width=3.5, dpi=600)


    #statistics to add to the plot (comparing for each protein, between the different types of seeding e.g. cross seed vs. monomer only vs. self seed)

    fit_params ['for_stats'] = fit_params[['t1', 't2']].apply('-'.join, axis=1)
    # Group the lag_phase values by treatment
    for protein, df in fit_params.groupby('t2'):
        print(protein)
        groups = [group['tlag'].values for name, group in df.groupby('for_stats')]
        # Run one-way ANOVA
        f_stat, p_value = stats.f_oneway(*groups)

        print(f"F-statistic: {f_stat:.3f}")
        print(f"P-value: {p_value:.5f}")


        # Run Tukey's HSD test
        tukey = pairwise_tukeyhsd(endog=df['tlag'], groups=df['for_stats'], alpha=0.05)

        # Print summary
        print(tukey)
