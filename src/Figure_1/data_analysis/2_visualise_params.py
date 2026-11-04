"""Visualize fitted kinetic parameters across seeding conditions.

This script reads the curated parameter table produced by the fitting workflow,
parses treatment annotations from sample identifiers, and generates bar-and-dot
plots for lag time and elongation rate across monomer, self-seeding, and
cross-seeding conditions. The resulting figures are intended for exploratory
comparison of kinetic trends across TDP-43 variants and treatment groups.
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def bar(fit_params, output_folder, monomer_only, group='t1', hue_order= None, order=None,hue_var='t2', pal='BuGn', x_var='t3', y_var='tlag', savefig=False, ylabel='Lag Phase (h)', xlabel='Seed added (%)', ylim=None):
    """Generate grouped bar-and-dot plots for one kinetic parameter.

    Args:
        fit_params (pandas.DataFrame): Table of fitted parameters and treatment
            annotations.
        output_folder (str): Directory for optional figure export.
        monomer_only (pandas.DataFrame): Control dataframe to append to each
            grouped subset before plotting.
        group (str, optional): Column used to define subplot groups.
        hue_order (list | None, optional): Hue plotting order.
        order (list | None, optional): X-axis category order.
        hue_var (str, optional): Column used for color grouping.
        pal (str | list | dict, optional): Seaborn palette specification.
        x_var (str, optional): Column plotted on the x-axis.
        y_var (str, optional): Column plotted on the y-axis.
        savefig (bool, optional): If ``True``, save the generated figure.
        ylabel (str, optional): Y-axis label.
        xlabel (str, optional): X-axis label.
        ylim (float | None, optional): Optional upper y-axis limit.
    """
    fig, axes = plt.subplots(len(fit_params[group].unique()), 1, sharex=False, sharey=True, figsize=(5,15))
    for i, (t1, df) in enumerate(fit_params.groupby(group)):
        if i>1:
            i = i-1
        if i==1:
            i = 0
        df = df[df['t3']!='seeds']

        control_df = monomer_only
        control_df['t1'] = t1
        df1 = pd.concat([df, control_df])

        sns.barplot(df1, x=x_var, y=y_var, hue=hue_var, hue_order=hue_order, order=order, palette=pal, dodge=True, edgecolor='darkgray', linewidth=0.5,capsize=0.1, ax=axes[i])
        axes[i].annotate(f'{t1}', xy=(-0.4, 14.8))
        sns.stripplot(data=df1, x=x_var, y=y_var, hue=hue_var, hue_order=hue_order, order=order, dodge=True, jitter=False, palette=pal, edgecolor='black', alpha=0.8, linewidth=1, ax=axes[i])
        if ylim:
            axes[i].set_ylim(0,ylim)
        axes[i].set_ylabel(ylabel)
        axes[i].set_xlabel(xlabel)
        plt.tight_layout()
    if savefig==True:
        plt.savefig(f'{output_folder}{y_var}_plot_all.svg')
        plt.savefig(f'{output_folder}{y_var}_plot_all.png')

def simple_bar(fit_params, output_folder, monomer_only, group='t1', hue_order= None, order=None,hue_var='t2', pal='YlOrRd', x_var='t3', y_var='tlag', savefig=False, ylabel='Lag Phase (h)', xlabel='Seed added (%)', ylim=None):
    """Generate a single-panel bar-and-dot plot for one kinetic parameter.

    Args:
        fit_params (pandas.DataFrame): Table of fitted parameters and treatment
            annotations.
        output_folder (str): Directory for optional figure export.
        monomer_only (pandas.DataFrame): Control dataframe to append to the
            plotted subset.
        group (str, optional): Column used to iterate over grouped subsets.
        hue_order (list | None, optional): Hue plotting order.
        order (list | None, optional): X-axis category order.
        hue_var (str, optional): Column used for color grouping.
        pal (str | list | dict, optional): Seaborn palette specification.
        x_var (str, optional): Column plotted on the x-axis.
        y_var (str, optional): Column plotted on the y-axis.
        savefig (bool, optional): If ``True``, save the generated figure.
        ylabel (str, optional): Y-axis label.
        xlabel (str, optional): X-axis label.
        ylim (float | None, optional): Optional upper y-axis limit.
    """
    fig, ax = plt.subplots(1, 1, sharex=False, sharey=True, figsize=(5,5))
    for i, (t1, df) in enumerate(fit_params.groupby(group)):
        i = i-1
        df = df[df['t3']!='seeds']

        
        control_df = monomer_only
        control_df['t1'] = t1
        df1 = pd.concat([df, control_df])

        ax=sns.barplot(df1, x=x_var, y=y_var, hue=hue_var, hue_order=hue_order, order=order, palette=pal, dodge=True, edgecolor='black', linewidth=0.5,capsize=0.1)
        sns.stripplot(data=df1, x=x_var, y=y_var, hue=hue_var, hue_order=hue_order, order=order, dodge=True, jitter=False, palette=pal, edgecolor='black', alpha=0.8, linewidth=1, ax=ax)
        handles, labels = ax.get_legend_handles_labels()
        l = plt.legend(handles[0:5], labels[0:3], bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        if ylim:
            ax.set_ylim(0,ylim)
        ax.set_ylabel(ylabel)
        ax.set_xlabel(xlabel)
        ax.tick_params(axis='x', rotation=45)
        plt.tight_layout()
    if savefig==True:
        plt.savefig(f'{output_folder}{y_var}_plot_all.svg')
        plt.savefig(f'{output_folder}{y_var}_plot_all.png')

if __name__ == "__main__":    
    input_folder = 'input_data/Figure_1/data_analysis/3_visualise_params/'
    output_folder = 'output/Figure_1/data_analysis/3_visualise_params/'

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    fit_params = pd.read_csv(f'{input_folder}fit_parameters.csv')

    fit_params.drop([col for col in fit_params.columns.tolist() if col in ['Unnamed: 0.1', 'Unnamed: 0']], axis=1, inplace=True)

    # Recover treatment annotations from the composite sample identifier stored
    # in the fitting output.
    fit_params[['t1', 't2', 't3', 'replicate']] = fit_params['inform'].str.split('-', expand=True)

    # Apply the empirical fit-quality threshold used for downstream comparison.
    fit_params = fit_params[fit_params['rquare']>0.99]
    fit_params['tlag'][fit_params['tlag']<=0] = 0
    params_to_plot = ['tlag', 'elongation_rate']

    # Plot each kinetic parameter using a shared monomer reference condition.
    for param in params_to_plot:
        monomer_only = fit_params[fit_params['inform']=='monomer-WT-0-1.0']
        monomer_only['t2'] = 'LCD_monomer'
        simple_bar(fit_params, output_folder, monomer_only, group='t3', hue_var='t1', x_var='t2', order=['WT', 'A382P', 'G335D', 'W385G'] , hue_order=['monomer', 'seeding', 'xseeding'],y_var=param, savefig=False, ylabel=param, xlabel='TDP MUTANT')


