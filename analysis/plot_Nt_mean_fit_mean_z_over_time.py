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
title_fontsize = 18
y_axislab_fontsize = 18
x_axislab_fontsize = 15
ticklab_fontsize = 12
annot_fontsize = 8
cbar_fontsize = 9
fig_width = 9
fig_height = 6
dpi = 400
n_ticklabels = 5


# data directory
if os.getcwd().split('/')[1] == 'home':
    datadir = '/home/deth/Desktop/tmp_ch2_stats_tests_dev/'
    analysis_dir = ('/home/deth/Desktop/CAL/research/projects/sim/ch2/'
                    'climate_change_adaptation_and_genomic_arch/analysis')
else:
    with open(('/global/scratch/users/drewhart/ch2/climate_change_adaptation_'
               'and_genomic_arch/analysis/outputdir.txt'), 'r') as f:
        datadir = f.read().strip()
    with open(('/global/scratch/users/drewhart/ch2/climate_change_adaptation_'
               'and_genomic_arch/analysis/analysisdir.txt'), 'r') as f:
        analysis_dir = f.read().strip()
    #datadir = '/global/scratch/users/drewhart/ch2/output/output'
    #analysis_dir = '/global/scratch/users/drewhart/ch2/output/analysis'
filename_patt = 'output_PID-.*_TS_DATA.csv'

# list of files to plot
files = [os.path.join(datadir, f) for f in os.listdir(datadir) if re.search(
                                                        filename_patt, f)]

# set linkages
linkages = ['independent', 'weak', 'strong']

# just read in and concatenate all files at once, rather than writing a loop,
# since even for 1000 its, each file should only be 100K and there should be
# ~1000 files, so reading them all in together should use only ~100MB RAM
df = pd.concat([pd.read_csv(f, na_filter=False) for f in files])

# get list of the ordered time steps
env_change_start = 2000
len_env_change_event = int(np.max(df.time_step)/3)
time_steps = [*range(env_change_start-len_env_change_event,
                     env_change_start+(2*len_env_change_event))]


def plot_ts_for_all_scenarios(df, redundancy, var, show_plots=False):

    assert var in ('Nt', 'mean_fit', 'mean_z'), ('var must be one of "Nt", '
                                                 '"mean_fit", or "mean_z"')

    # set genicities, based on redundancy
    assert redundancy in ['lo', 'hi']
    genicities = [4, 20, 100]
    if redundancy == 'hi':
        genicities = [g*2 for g in genicities]

    # get min and max y values
    ymin = np.min(df[var])
    ymax = np.max(df[var])

    # create a fig for non-null sims and another for null
    fig = plt.figure(dpi=dpi, figsize=(fig_width, fig_height))
    gs = fig.add_gridspec(nrows=3, ncols=3,
                                    width_ratios=[1]*3)

    uncertainty_colors = {'non-null':'#c4626e77', # red
                          'null': '#79c2d977', # blue
                         }

    # loop over linkage, then over genicity
    for linkage_i, linkage in enumerate(linkages):
        for genicity_j, genicity in enumerate(genicities):

            ax = fig.add_subplot(gs[linkage_i, genicity_j])

            # plot a dim gray box over the post-environmental change time period
            after = plt.Polygon([[time_steps[-1], ymin],
                                 [time_steps[-1], ymax],
                                 [time_steps[-len_env_change_event-1], ymax],
                                 [time_steps[-len_env_change_event-1], ymin],
                                 [time_steps[-1], ymin]])
            after.set_color('#aaaaaa66')
            ax.add_patch(after)

            # loop over nullnesses
            for nullness in ['null', 'non-null']:

                # subset the df
                subdf = df[(df.nullness == nullness) &
                           (df.genicity == genicity) &
                           (df.linkage == linkage)]

                # find the mean 5th and 95th percentiles for each time step
                means = []
                pctiles_5 = []
                pctiles_95 = []
                for ts_i, ts in enumerate(time_steps):
                    vals = subdf[subdf.time_step == ts_i][var]
                    means.append(np.mean(vals))
                    pctiles = np.percentile(vals, [5, 95])
                    pctiles_5.append(pctiles[0])
                    pctiles_95.append(pctiles[1])

                # plot this scenario
                uncertainty = plt.Polygon([*zip(time_steps, pctiles_5)] +
                                [*zip(time_steps[::-1], pctiles_95[::-1])] +
                                [[time_steps[0], pctiles_5[0]]])
                uncertainty.set_color(uncertainty_colors[nullness])
                ax.add_patch(uncertainty)
                ax.plot([env_change_start, env_change_start],
                        [ymin, ymax], ':r', alpha = 1)
                ax.plot([env_change_start+len_env_change_event,
                         env_change_start+len_env_change_event],
                        [ymin, ymax], ':r', alpha = 1)
                ax.plot(time_steps, means, 'k', linewidth=0.75)

                # titles and labels, other axes params
                if linkage_i == 0:
                    ax.set_title(genicity,
                                 fontdict={'fontsize': title_fontsize})
                if linkage_i == 2:
                    ax.set_xlabel('time step',
                                  fontdict={'fontsize': x_axislab_fontsize})
                if genicity_j == 0:
                    ax.set_ylabel(linkage,
                                  fontdict={'fontsize': y_axislab_fontsize})
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
        fig.subplots_adjust(left=0.13,
                            bottom=0.1,
                            right=0.99,
                            top=0.88,
                            wspace=0.1,
                            hspace=0.1)
        #fig.suptitle('%s' % var)
        fig.savefig(os.path.join(analysis_dir,
                    'ch2_%s_%sREDUND_over_time.jpg' % (var, redundancy)),
                            dpi=dpi, orientation='landscape')

for redundancy in ['lo', 'hi']:
    plot_ts_for_all_scenarios(df, redundancy, 'Nt')
    plot_ts_for_all_scenarios(df, redundancy, 'mean_fit')
    plot_ts_for_all_scenarios(df, redundancy, 'mean_z')
