"""Compare original and particle-normalized seeding responses for Figure 3F.

Reads the prepared paired-response table and exports a patterned
bar-and-point comparison across batches.
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patches as mpatches
from statannotations.Annotator import Annotator
from src.paperfig_utils import save_figure_svg_png, set_params_for_plotting

def bardotplot(data, xcol, ycol, order, hue=None, hue_order=None, scat_hue=None, scat_hue_order=None, palette=False, xlabel='', ylabel=False, pairs=False, correction=None, xticks=None, groups=None, group_label_y=-0.18, group_line_y=-0.05, ax=None, legend='', dot_size=5, cap_size=0.2, cap_width=2):
    """Draw a bar plot with individual observations and optional annotations."""
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

def set_purple_palette():
    """Display the purple palette used for the normalized-particle comparison."""
    hex_colors = ["#4a006f","#470a77","#45147e","#421e86","#3f288d","#3d3195","#3a3b9c","#3745a4","#354fab","#3259b3"]
    custom_palette = sns.color_palette(hex_colors)
    sns.palplot(custom_palette)

if __name__ == '__main__':

    input_folder = 'input_data/Figure_3/'
    output_folder = 'output/Figure_3/F_normed_particles/'

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    set_params_for_plotting(fontsize=10)


    # Read the prepared original and particle-normalized response values.
    norm_correction = pd.read_csv(f'{input_folder}Fig3F_particle_num_seeding.csv')


    # plot the original vs. normalised  differences
    #format being 'normed batch _ original batch' where original batch refers to the batch of fibrils which the normalised batch has had particle numbers equalised to
    xaxis_order = ['batch2_batch1', 'batch3_batch1', 'batch1_batch2','batch3_batch2', 'batch1_batch3','batch2_batch3' ]

    fig, ax = bardotplot(norm_correction, 'plot_group', 'divided', xaxis_order, palette='viridis',  hue='type_difference', hue_order=['original_difference', 'normed_difference'], xlabel='',scat_hue='type_difference',scat_hue_order=['original_difference', 'normed_difference'], ylabel='Degree of change due to normalisation', pairs=False, correction=None, xticks=None, groups=None, group_label_y=-0.18, group_line_y=-0.05, ax=None, legend=True, dot_size=5, cap_size=0.15, cap_width=2)
    ax.tick_params(axis='x', rotation=90)
    # Display the palette used by the manuscript.
    set_purple_palette()

    for stype, patch in zip([('original_difference', mut) for mut in xaxis_order]+ [('normed_difference', mut) for mut in xaxis_order], ax.patches):
        if stype[0] =='normed_difference':
            patch.set_facecolor("#45147e")
            patch.set_edgecolor('#fff')
            patch.set_linewidth(1)
            patch.set_hatch('///')

        else:
            patch.set_facecolor("#45147e")
            patch.set_edgecolor('#fff')
            patch.set_linewidth(1)
        
    # Add a legend identifying the original and normalized conditions.
    mon_patch = mpatches.Patch(facecolor='dimgray', label='original_difference', edgecolor='#fff', linewidth=1)
    seed_patch = mpatches.Patch(facecolor='dimgray', label='normed_difference', edgecolor='#fff', linewidth=0, hatch='///')

    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(reversed(labels), reversed(handles)))
    ax.legend([mon_patch, seed_patch], ['Original', 'Particle # equalised'], frameon=False, loc='upper left', handlelength=2., handletextpad=0.25, bbox_to_anchor=(-0.005, 1.15), ncols=3)
    plt.tight_layout()
    save_figure_svg_png(fig, f'{output_folder}3F_normed_particles', width=3.5)
