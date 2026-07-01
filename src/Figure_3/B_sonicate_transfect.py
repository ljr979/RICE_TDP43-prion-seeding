import os, re
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger
from src.paperfig_utils import save_figure_svg_png, set_params_for_plotting
import matplotlib.pyplot as plt
from statannotations.Annotator import Annotator
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from scipy import stats
    
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
        dodge=False,
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
    groups = [group[value_col].values for name, group in df.groupby(grouping_col)]
    # Run one-way ANOVA
    f_stat, p_value = stats.f_oneway(*groups)

    print(f"F-statistic: {f_stat:.3f}")
    print(f"P-value: {p_value:.5f}")


    # Run Tukey's HSD test
    tukey = pairwise_tukeyhsd(endog=df[value_col], groups=df[grouping_col], alpha=0.05)

    # Print summary
    print(tukey)


if __name__ == '__main__':

    input_folder = 'input_data/Figure_3/'
    output_folder = 'output/Figure_3/B_sonicate_transfect/'


    if not os.path.exists(output_folder):
        os.makedirs(output_folder)


    set_params_for_plotting(fontsize=10)

    #read in the data collated from three replicates seeding variants. These have been pre-processed (cellprofiler + data analysis scripts) and cleaned before collating into one dataframe for plotting 
    collated_sonication_data = pd.read_csv(f'{input_folder}Fig3B_sonicated_aggcount.csv')

    xaxis_order = ['vehicle', 'unsonicated', '10s', '2min']

    #plot the raw counts  (rather than normalised to unsonicated)
    fig, ax = bardotplot(collated_sonication_data, 'treatment1', 'value', xaxis_order, palette='Purples',  scat_hue='Experiment', xlabel='', ylabel='Cells with aggregates (%)', pairs=False, correction=None, xticks=None, groups=None, group_label_y=-0.18, group_line_y=-0.05, ax=None, legend='', dot_size=5, cap_size=0.15, cap_width=2)


    #now calculate for each biological replicate, the ratio between unsonicated and sonicated (i.e. the fold change in aggregation as a result of the sonication)
    normalised_to_unson = []
    for replicate, df in collated_sonication_data.groupby('Experiment'):
        unson_value = df[df['treatment1'] == 'unsonicated']['value'].values[0]
        df['relative_to_unson_ratio'] = df['value']/unson_value
        normalised_to_unson.append(df)
    normalised_to_unson = pd.concat(normalised_to_unson)

    #set new purple / blue palette
    hex_colors = ["#4a006f","#470a77","#45147e","#421e86","#3f288d","#3d3195","#3a3b9c","#3745a4","#354fab","#3259b3"]
    custom_palette = sns.color_palette(hex_colors)
    sns.palplot(custom_palette)

    fig, ax = bardotplot(normalised_to_unson, 'treatment1', 'relative_to_unson_ratio', xaxis_order, palette=custom_palette,  xlabel='', ylabel='Seeding relative to unsonicated', pairs=False, correction=None, xticks=None, groups=None, group_label_y=-0.18, group_line_y=-0.05, ax=None, legend='', dot_size=5, cap_size=0.15, cap_width=2)

    save_figure_svg_png(fig, f'{output_folder}Fig3B_sonicate_transfect', width=3.5)

    anova_tukeyshsd(df=normalised_to_unson, grouping_col='treatment1', value_col='relative_to_unson_ratio')