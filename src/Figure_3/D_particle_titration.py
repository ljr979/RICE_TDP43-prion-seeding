import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger
import matplotlib.patches as mpatches
from src.paperfig_utils import save_figure_svg_png, set_params_for_plotting

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

def set_purple_palette():
    hex_colors = ["#4a006f","#470a77","#45147e","#421e86","#3f288d","#3d3195","#3a3b9c","#3745a4","#354fab","#3259b3"]
    custom_palette = sns.color_palette(hex_colors)
    sns.palplot(custom_palette)

if __name__ == '__main__':

    input_folder = 'input_data/Figure_3/'
    output_folder = 'output/Figure_3/D_particles_titration/'

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    set_params_for_plotting(fontsize=10)

    titrations = pd.read_csv(f'{input_folder}Fig3D_particle_titration.csv')
    titrations['dilution_factor'] = titrations['dilution_factor'].astype(str)
    xaxis_order = ['batch1', 'batch2', 'batch3']

    fig, ax = bardotplot(titrations, 'original_batch', 'effect_of_dilution', xaxis_order, palette='viridis',  hue='dilution_factor', hue_order=['1.0', '0.5'], xlabel='Seed Batch',scat_hue='dilution_factor',scat_hue_order=['1.0', '0.5'], ylabel='Effect of dilution on seeding', pairs=False, correction=None, xticks=None, groups=None, group_label_y=-0.18, group_line_y=-0.05, ax=None, legend=True, dot_size=5, cap_size=0.15, cap_width=2)
    ax.tick_params(axis='x', rotation=90)

    for stype, patch in zip([('1.0', mut) for mut in xaxis_order]+ [('0.5', mut) for mut in xaxis_order], ax.patches):
        if stype[0] =='0.5':
            patch.set_facecolor("#45147e")
            patch.set_edgecolor('#fff')
            patch.set_linewidth(1)
            patch.set_hatch('///')

        else:
            patch.set_facecolor("#45147e")
            patch.set_edgecolor('#fff')
            patch.set_linewidth(1)
        

    #need to add in STATISTICS and also LEGEND ENTRIES
    mon_patch = mpatches.Patch(facecolor='dimgray', label='1.0', edgecolor='#fff', linewidth=1)
    seed_patch = mpatches.Patch(facecolor='dimgray', label='0.5', edgecolor='#fff', linewidth=0, hatch='///')


    #ax.axhline(1, linestyle='--', color='lightgrey')
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(reversed(labels), reversed(handles)))
    ax.legend([mon_patch, seed_patch], ['1.0x', '0.5x'], frameon=False, loc='upper left', handlelength=2., handletextpad=0.25, bbox_to_anchor=(-0.005, 1.15), ncols=3)

    plt.tight_layout()
    save_figure_svg_png(fig, f'{output_folder}3B_dilution_seeding', width=3.5)