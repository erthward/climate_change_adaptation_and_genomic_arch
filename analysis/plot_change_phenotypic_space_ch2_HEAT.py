import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from palettable.cartocolors.sequential import PinkYl_7
from palettable.scientific.sequential import Acton_20_r
import seaborn as sns
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
fig_width = 14.5
fig_height = 5.6
dpi = 400
n_ticklabels = 5

# data directory
if os.getcwd().split('/')[1] == 'home':
    datadir = '/home/deth/Desktop/tmp_ch2_stats_tests_dev/'
    analysis_dir = '/home/deth/Desktop/tmp_ch2_stats_tests_dev/'
else:
    datadir = '/global/scratch/users/drewhart/ch2/output/output'
    analysis_dir = '/global/scratch/users/drewhart/ch2/output/analysis'

# way to hard-code the vmax values for the scenarios, given that it would
# be kind of a pain to write code to do just this
linkages = ['independent', 'weak', 'strong']
genicities =  [2, 4, 10, 20, 50, 100]
vmaxes = {l: {g: 0.3 for g in genicities} for l in linkages}

# get the horizontal values used in the non-shifting landscape layer
# and in the before- and after-shifting shifting layer,
# then set as dict of x and y values for plotting horizontal expectation lines
stable = np.linspace(0, 1, 50)[[0,-1]]
b4 = np.linspace(0, 1, 50)[[0,-1]]
af = np.linspace(0.5, 1, 50)[[0,-1]]
expec_lines = {2499: (b4, stable),
               2624: ((b4+af)/2, stable),
               2749: (af, stable),
              }


def get_min_pw_diff(vals):
    """
    get the minimum pairwise diff in a set of values
    """
    diffs = []
    for i, val in enumerate(vals):
        for val2 in vals[i+1:]:
            diffs.append(np.abs(val - val2))
    diffs = [diff for diff in diffs if diff > 0]
    return(min(diffs))


def plot_phenotypic_shift(linkage, genicity):

    assert linkage in ['independent', 'weak', 'strong']
    assert genicity in [2, 4, 10, 20, 50, 100]

    # get candidate filenames for time-step-2499 files
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
            candidate_filenames = [os.path.join(datadir, dirname, 'it--1', 'spp-spp_0',
                                            fn) for fn in candidate_filenames]
            # only add this directory and its files to the analysis if I got all 3 timeteps,
            # otherwise print warning
            if len(candidate_filenames) == 3:
                assert len([fn for fn in candidate_filenames if '-2499_' in fn])== 1
                assert len([fn for fn in candidate_filenames if '-2624_' in fn])== 1
                assert len([fn for fn in candidate_filenames if '-2749_' in fn])== 1
                # sort files by timestep
                candidate_filenames = [*np.sort(candidate_filenames)]
                filenames[dirname] = candidate_filenames
            else:
                print(('\n\nWARNING: following directory did not contain'
                       '3 files, 1 for each of the right timesteps:\n\n\t'
                       '%s\n\n') % dirname)

    # create the figure
    fig = plt.figure(dpi=dpi,
                     figsize=(fig_width, fig_height)
                    )
    #fig.suptitle(genicity, fontdict={'fontsize': suptitle_fontsize})
    gs = fig.add_gridspec(nrows=1, ncols=3, width_ratios=[0.95, 0.95, 1.15])

    # loop through the three time steps to be analyzed
    time_steps = {'before': 2499,
                  'during': 2624,
                  'after': 2749}
    for time_step_n, time_step_info in enumerate(time_steps.items()):
        title, time_step = time_step_info

        ax = fig.add_subplot(gs[0, time_step_n])
        #ax.axis('off')
        if linkage == 'independent':
            ax.set_title(title, fontdict={'fontsize': title_fontsize})

        # set up breaks for the 2D histogram (i.e., heatmap)
        # NOTE: ADDING 2 TO genicity GENERATES A NUMBER OF BINS EQUAL TO
        #       THE NUMBER OF POSSIBLE GENOTYPES AND CREATES A SET OF BINS
        #       THAT CAPTURES EXACTLY 1 POSSIBLE GENOTYPE IN EACH
        brks = np.linspace(0, 1, genicity + 2)
        # increase the last break (1.0) slightly, so that < captures phenotypes
        # of 1.0
        brks[-1] *= 1.000001
        # make the heatmap array
        arr = np.zeros([len(brks)-1]*2)

        # loop through all the sims' directories
        for sim_dirname, sim_filenames in filenames.items():
            print('\nNOW PROCESSING FILES IN: %s\n\n' % sim_dirname)
            df = pd.read_csv(sim_filenames[time_step_n])

            # get phenotype data
            zs = np.stack([np.array([float(n) for n in val.lstrip(
                '[').rstrip(']').split(', ')]) for val in df['z']])
            z_trt0 = zs[:, 0]
            z_trt1 = zs[:, 1]

            # fill the heatmap array (i.e., 2d histogram)
            # NOTE: to make shift appear L-->R (matching spatial shift across
            # landscape), i-->rows-->trait1-->stable trait,
            # whereas j-->cols-->trait0-->shifting trait
            for i, brk_trt1 in enumerate(brks[:-1]):
                for j, brk_trt0 in enumerate(brks[:-1]):
                   ct = np.sum((zs[:,0]>=brk_trt0) *
                                (zs[:,0]<brks[j+1]) *
                                (zs[:,1]>=brk_trt1) *
                                (zs[:,1]<brks[i+1]))
                   arr[i, j] += ct

        arr = arr/np.sum(arr)
        assert np.allclose(np.sum(arr), 1)

        # create labels array
        annot = np.array(['%0.3f' % n for n in arr.ravel()]).reshape(arr.shape)

        # plot the heatmap
        sns.heatmap(arr,
                    #vmin=0,
                    #vmax=vmaxes[linkage][genicity],
                    #annot=annot,
                    #annot_kws = {'fontsize': annot_fontsize},
                    fmt="",
                    cmap=Acton_20_r.mpl_colormap,
                    #cmap=PinkYl_7.mpl_colormap,
                    #cmap=mpl.cm.Wistia,
                    #cbar=False,
                    cbar=time_step_n==2,
                    linewidths=0.3,
                    ax=ax)

        # add diagonal 1:1, and expectation line for comparison at later steps
        ax.plot(max(ax.get_xlim()) * expec_lines[2499][0],
                (max(ax.get_ylim())*expec_lines[2499][1]),
                '-k', alpha=0.4, linewidth=0.5)
        if time_step > 2500:
            ax.plot(max(ax.get_xlim()) * expec_lines[time_step][0],
                    (max(ax.get_ylim())*expec_lines[time_step][1]),
                    '-k', alpha=0.4, linewidth=0.5)

        # set ticks and ticklabels and axis labels
        ticks = np.linspace(0, genicity+1, n_ticklabels)
        ticklabels = np.linspace(0, 1, n_ticklabels)
        ax.set_xticks(ticks)
        ax.set_xticklabels(ticklabels, fontdict={'fontsize':ticklab_fontsize})
        for tick in ax.get_xticklabels():
            tick.set_rotation(0)
        if linkage == 'strong':
            ax.set_xlabel('shifting trait',
                        fontdict={'fontsize': axislab_fontsize})
        else:
            ax.set_xlabel('')
        if time_step_n == 0 and genicity == 4:
            ax.set_yticks(ticks)
            ax.set_yticklabels([tl if (n > 0) else '' for n,
                                        tl in enumerate(ticklabels)],
                               fontdict={'fontsize':ticklab_fontsize})
            for tick in ax.get_yticklabels():
                    tick.set_rotation(90)
            ax.set_ylabel('stable trait',
                          fontdict={'fontsize': axislab_fontsize})
        else:
            ax.set_yticks(())
            ax.set_ylabel('')

        # other plot controls
        ax.invert_yaxis()
        ax.tick_params(axis=u'both', which=u'both',length=0) # make ticks len==0

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
            fig = plot_phenotypic_shift(linkage, genicity)
            # save the fig
            fig.savefig(os.path.join(analysis_dir,
                        'phenotypic_shift_L%s_G%s_HEAT.png' % (linkage,
                                                    str(genicity).zfill(2))),
                        dpi=dpi,
                        orientation='landscape',
                       )
