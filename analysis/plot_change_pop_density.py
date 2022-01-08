import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as mplPolygon
from collections import Counter as C
from palettable.cartocolors.sequential import PinkYl_7
from palettable.scientific.sequential import Acton_20_r
from shapely.geometry import Polygon as shapelyPolygon
import statsmodels.api as sm
import sys
import re
import os

# TODO:
    # figure out why expectation line and heatmap grid uneven when run on Savio
    # get rid of right-hand black strip in panelled plot output


# plot params
suptitle_fontsize = 50
title_fontsize = 40
axislab_fontsize = 20
ticklab_fontsize = 14
annot_fontsize = 14
cbar_fontsize = 14
fig_width = 11.2
fig_height = 5.6
dpi = 400
n_ticklabels = 5
linewidth = 1
linealpha = 0.8
marker_size = 5

# data directory
if os.getcwd().split('/')[1] == 'home':
    datadir = '/home/deth/Desktop/tmp_ch2_stats_tests_dev/'
    analysis_dir = '/home/deth/Desktop/tmp_ch2_stats_tests_dev/'
else:
    with open(('/global/scratch/users/drewhart/ch2/climate_change_adaptation_'
               'and_genomic_arch/analysis/outputdir.txt'), 'r') as f:
        datadir = f.read().strip()
    with open(('/global/scratch/users/drewhart/ch2/climate_change_adaptation_'
               'and_genomic_arch/analysis/analysisdir.txt'), 'r') as f:
        analysis_dir = f.read().strip()
    #datadir = '/global/scratch/users/drewhart/ch2/output/output'
    #analysis_dir = '/global/scratch/users/drewhart/ch2/output/analysis'

# lists of all possible linkage and genicity values
linkages = ['independent', 'weak', 'strong']
genicities =  [2, 4, 10, 20, 50, 100]

# environmental change time steps
change_start_t = 999
change_half_t = 1124
change_end_t = 1249


def plot_pop_density_shift(linkage, genicity, just_get_max_dens_per_run=False,
                           overall_max_dens_per_run=None):

    assert linkage in ['independent', 'weak', 'strong']
    assert genicity in [2, 4, 10, 20, 50, 100]

    # get candidate filenames for change-start-time-step files
    dirname_patt = 'mod-non-null_L%s_G%i_its0_' % (linkage, genicity)
    filename_patt = ('mod-non-null_L%s_G%i_its0_randID\d{7}PID-'
                     '\d{5,6}_it--1_t-2\d{3}_spp-spp_0.csv') % (linkage, genicity)

    filenames = {}
    for dirname in os.listdir(datadir):
        if (dirname.startswith('GNX_') and
            re.search(dirname_patt, dirname)):
            candidate_filenames = [fn for fn in os.listdir(os.path.join(
                datadir, dirname, 'it--1', 'spp-spp_0')) if re.search(filename_patt,
                                                             fn)]
            # drop the middle-timestep files
            candidate_filenames = [fn for fn in candidate_filenames if not
                                   re.search('%i' % change_half_t, fn)]
            candidate_filenames = [os.path.join(datadir, dirname, 'it--1', 'spp-spp_0',
                                            fn) for fn in candidate_filenames]
            # only add this directory and its files to the analysis if I got all 3 timeteps,
            # otherwise print warning
            if len(candidate_filenames) == 2:
                assert len([fn for fn in candidate_filenames if '-%i_' % change_start_t in fn])== 1
                #assert len([fn for fn in candidate_filenames if '-%i_' % change_half_t in fn])== 1
                assert len([fn for fn in candidate_filenames if '-%i_' % change_end_t in fn])== 1
                # sort files by timestep
                candidate_filenames = [*np.sort(candidate_filenames)]
                filenames[dirname] = candidate_filenames
            else:
                print(('\n\nWARNING: following directory did not contain'
                       '2 valid files, 1 for each of the right timesteps:\n\n\t'
                       '%s\n\n') % dirname)

    # create the figure
    fig = plt.figure(dpi=dpi,
                     figsize=(fig_width, fig_height)
                    )
    #fig.suptitle(genicity, fontdict={'fontsize': suptitle_fontsize})
    gs = fig.add_gridspec(nrows=1, ncols=2, width_ratios=[0.95, 1.15])

    # loop through the three time steps to be analyzed
    time_steps = {'before': change_start_t,
                  #'during': change_half_t,
                  'after': change_end_t}
    for time_step_n, time_step_info in enumerate(time_steps.items()):
        title, time_step = time_step_info

        # array equal in size to landscape grid, to keep count of pop densities
        cts = np.zeros((50,50))

        # var to keep track of file count
        file_ct = 0

        ax = fig.add_subplot(gs[0, time_step_n])
        #ax.axis('off')
        if linkage == 'independent':
            ax.set_title(title, fontdict={'fontsize': title_fontsize})

        # loop through all the sims' directories
        for sim_dirname, sim_filenames in filenames.items():
            print('\nNOW PROCESSING FILES IN: %s\n\n' % sim_dirname)
            df = pd.read_csv(sim_filenames[time_step_n])

            # get location data
            xs = df['x']
            ys = df['y']

            # add to cts
            for i, j in zip(ys.astype(int), xs.astype(int)):
                cts[i,j] += 1

            file_ct += 1

        # just return the max density, if requested
        if just_get_max_dens_per_run:
            return np.max(cts)/file_ct

        # plot the normalized raster
        if overall_max_dens_per_run is None:
            overall_max_dens_per_run = 4*file_ct
        im = ax.imshow(cts/file_ct,
                       #cmap=mpl.cm.bone,
                       cmap=mpl.cm.gray,
                       vmin=0,
                       vmax=overall_max_dens_per_run,
                      )
        # add colorbar to rightmost plot
        if time_step > change_half+2:
            plt.colorbar(im)

        # set ticks and ticklabels and axis labels
        ax.set_xticks([])
        ax.set_xticklabels([])
        ax.set_yticks([])
        ax.set_yticklabels([])
        ax.set_xlim([0,50])
        ax.set_ylim([0,50])
        ax.set_xlabel('')
        ax.set_ylabel('')

    # adjust colorbar label fontsize
    fig.get_axes()[-1].tick_params(labelsize=cbar_fontsize)

    # adjust suplot spacing
    plt.subplots_adjust(left=0.05,
                        bottom=0.13,
                        right=0.96,
                        top=0.85,
                        wspace=0.14,
                        hspace=0.05)

    # return fig
    return fig


# first get max densities, to normalize plots against
all_max_dens_per_run = []
for linkage in ['independent', 'weak', 'strong']:
    for genicity in [2, 4, 10, 20, 50, 100]:
        print('\n\n======================\n\n')
        print('\tLINKAGE: %s' % linkage)
        print('\tGENICITY: %i' % genicity)
        dirname_patt = 'mod-non-null_L%s_G%i_its0_' % (linkage, genicity)
        dirs = os.listdir(datadir)
        candidate_dirs = [d for d in dirs if re.search(dirname_patt, d)]
        if len(candidate_dirs) > 0:
            # make the fig
            max_dens_per_run = plot_pop_density_shift(linkage, genicity,
                                              just_get_max_dens_per_run=True)
            all_max_dens_per_run.append(max_dens_per_run)
overall_max_dens_per_run = np.max(all_max_dens_per_run)


# dict to store density change metrics
dens_change_dict = {'linkage': [], 'genicity': [], 'dens_change': []}
# produce plots for all scenarios
for linkage in ['independent', 'weak', 'strong']:
    for genicity in [2, 4, 10, 20, 50, 100]:
        print('\n\n======================\n\n')
        print('\tLINKAGE: %s' % linkage)
        print('\tGENICITY: %i' % genicity)
        dirname_patt = 'mod-non-null_L%s_G%i_its0_' % (linkage, genicity)
        dirs = os.listdir(datadir)
        candidate_dirs = [d for d in dirs if re.search(dirname_patt, d)]
        if len(candidate_dirs) > 0:
            # make the fig
            fig = plot_pop_density_shift(linkage, genicity,
                            overall_max_dens_per_run=overall_max_dens_per_run)
            # save the fig
            fig.savefig(os.path.join(analysis_dir,
                        'pop_density_shift_L%s_G%s.png' % (linkage,
                                                    str(genicity).zfill(2))),
                        dpi=dpi,
                        orientation='landscape',
                       )
