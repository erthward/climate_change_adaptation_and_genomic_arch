import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import re
import os

"""
Plots the progression of the species' population size, mean fitness, or
mean phenotype of the trait adapted to the shifting gradient,
from a length of time prior to the environmental change and equal in length
to the environmental change event until the end of the sim
"""

# plot params
title_fontsize = 17
axislab_fontsize = 13
ticklab_fontsize = 6
annot_fontsize = 8
cbar_fontsize = 9
fig_width = 9
fig_height = 6
dpi = 400
n_ticklabels = 5

# time step when env change event starts
env_change_start = 50
#env_change_start = 2500

# data directory
datadir = '/home/deth/Desktop/CAL/research/projects/sim/geonomics'
#datadir = '/global/scratch/users/drewhart/ch2/output/analysis_dir'
filename_patt = 'PRACTICE_TS_DATA.*csv'
#filename_patt = 'output_PID-.*_TS_DATA.csv'

# list of files to plot
files = [os.path.join(datadir, f) for f in os.listdir(datadir) if re.search(
                                                        filename_patt, f)]

# way to hard-code the vmax values for the scenarios, given that it would
# be kind of a pain to write code to do just this
linkages = ['independent', 'weak', 'strong']
genicities =  [4, 20, 100]

# just read in and concatenate all files at once, rather than writing a loop,
# since even for 1000 its, each file should only be 100K and there should be
# ~1000 files, so reading them all in together should use only ~100MB RAM
df = pd.concat([pd.read_csv(f, na_filter=False) for f in files])

# get list of the ordered time steps
time_steps = np.sort(np.unique(df.time_step))


def plot_ts_for_all_scenarios(df, var, show_plots=False):

    assert var in ('Nt', 'mean_fit', 'mean_z'), ('var must be one of "Nt", '
                                                 '"mean_fit", or "mean_z"')

    # get min and max y values
    ymin = np.min(df[var])
    ymax = np.max(df[var])

    # create a fig for non-null sims and another for null
    fig_nonnull = plt.figure(dpi=dpi, figsize=(fig_width, fig_height))
    fig_null = plt.figure(dpi=dpi, figsize=(fig_width, fig_height))
    gs_nonnull = fig_nonnull.add_gridspec(nrows=3, ncols=3,
                                          width_ratios=[1]*3)
    gs_null = fig_null.add_gridspec(nrows=3, ncols=3,
                                    width_ratios=[1]*3)

    # package fig objects to be looped over
    fig_objs = {'non-null': [fig_nonnull, gs_nonnull],
                'null': [fig_null, gs_null]
               }

    uncertainty_colors = {'non-null':'#c4626e77', # red
                          'null': '#79c2d977', # blue
                         }

    # loop over linkage, then over genicity
    for linkage_i, linkage in enumerate(linkages):
        for genicity_j, genicity in enumerate(genicities):

            # loop over nullnesses
            for nullness, objs in fig_objs.items():
                fig, gs = objs

                # subset the df
                subdf = df[(df.nullness == nullness) &
                           (df.genicity == genicity) &
                           (df.linkage == linkage)]

                # find the mean 5th and 95th percentiles for each time step
                means = []
                pctiles_5 = []
                pctiles_95 = []
                for ts in time_steps:
                    vals = subdf[subdf.time_step == ts][var]
                    means.append(np.mean(vals))
                    pctiles = np.percentile(vals, [5, 95])
                    pctiles_5.append(pctiles[0])
                    pctiles_95.append(pctiles[1])

                # plot this scenario
                ax = fig.add_subplot(gs[linkage_i, genicity_j])
                uncertainty = plt.Polygon([*zip(time_steps, pctiles_5)] +
                                [*zip(time_steps[::-1], pctiles_95[::-1])] +
                                [[time_steps[0], pctiles_5[0]]])
                uncertainty.set_color(uncertainty_colors[nullness])
                ax.add_patch(uncertainty)
                ax.plot([env_change_start, env_change_start],
                        [ymin, ymax], ':r', alpha = 1)
                ax.plot(time_steps, means, 'k')

                # titles and labels, other axes params
                if linkage_i == 0:
                    ax.set_title(genicity,
                                 fontdict={'fontsize': title_fontsize})
                if linkage_i == 2:
                    ax.set_xlabel('time step',
                                  fontdict={'fontsize': axislab_fontsize})
                if genicity_j == 0:
                    ax.set_ylabel(linkage,
                                  fontdict={'fontsize': axislab_fontsize})
                ax.tick_params(labelsize = ticklab_fontsize)
                if linkage_i < 2:
                    ax.set_xticks([])
                    ax.set_xticklabels([])
                if genicity_j > 0:
                    ax.set_yticks([])
                    ax.set_yticklabels([])
                ax.set_ylim([ymin,ymax])
                ax.set_xlim([np.min(time_steps), np.max(time_steps)])


    if show_plots:
        plt.show()
    else:
        # adjust suplot spacing, and add suptitle
        for nullness, fig in zip(['non-null', 'null'],
                                 [fig_nonnull, fig_null]):
            fig.subplots_adjust(left=0.08,
                                bottom=0.1,
                                right=0.99,
                                top=0.88,
                                wspace=0.1,
                                hspace=0.1)
            fig.suptitle('%s: %s' % (var, nullness))
        fig_nonnull.savefig(os.path.join(datadir,
                                         'ch2_%s_over_time.jpg' % var),
                            dpi=dpi, orientation='landscape')
        fig_null.savefig(os.path.join(datadir,
                                      'ch2_%s_over_time_NULL.jpg' % var),
                            dpi=dpi, orientation='landscape')


plot_ts_for_all_scenarios(df, 'Nt')
plot_ts_for_all_scenarios(df, 'mean_fit')
plot_ts_for_all_scenarios(df, 'mean_z')
