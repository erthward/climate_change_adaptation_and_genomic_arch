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

"""
Enter name of filename corresponding to time step 2499.
Script will plot for time steps 2499, 2624, and 2749.

"""

# plot params
title_fontsize = 24
axislab_fontsize = 16
ticklab_fontsize = 11
annot_fontsize = 8
cbar_fontsize = 11
fig_width = 14
fig_height = 5
dpi = 400

# way to hard-code the vmax values for the scenarios, given that it would
# be kind of a pain to write code to do just this
linkages = ['independent', 'weak', 'strong']
genicities =  [2, 4, 10, 20, 50, 100]
vmaxes = {l: {g: 0.3 for g in genicities} for l in linkages}


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
    for dirname in os.listdir('.'):
        if (dirname.startswith('GNX_') and
            re.search(dirname_patt, dirname)):
            candidate_filenames = [fn for fn in os.listdir(os.path.join(
                dirname, 'it--1', 'spp-spp_0')) if re.search(filename_patt,
                                                             fn)]
            candidate_filenames = [os.path.join(dirname, 'it--1', 'spp-spp_0',
                                            fn) for fn in candidate_filenames]
            # only add this directory and its files to the analysis if I got all 3 timeteps,
            # otherwise print warning
            if len(candidate_filenames) == 3:
                assert len([fn for fn in candidate_filenames if '-2499_' in fn])== 1
                assert len([fn for fn in candidate_filenames if '-2624_' in fn])== 1
                assert len([fn for fn in candidate_filenames if '-2749_' in fn])== 1
                filenames[dirname] = candidate_filenames
            else:
                print(('\n\nWARNING: following directory did not contain'
                       '3 files, 1 for each of the right timesteps:\n\n\t'
                       '%s\n\n') % dirname)

    # create the figure
    fig = plt.figure(dpi=dpi,
                     figsize=(fig_width, fig_height)
                    )
    gs = fig.add_gridspec(nrows=1, ncols=3, width_ratios=[0.95, 0.95, 1.15])

    # loop through the three time steps to be analyzed
    time_steps = {'before': 2499,
                  'during': 2624,
                  'after': 2749}
    for time_step_n, time_step_info in enumerate(time_steps.items()):
        title, time_step = time_step_info

        ax = fig.add_subplot(gs[0, time_step_n])
        #ax.axis('off')
        ax.set_title(title, fontdict={'fontsize': title_fontsize})

        # read in a sample dataset and use that to get the
        # characteristic break for this dataset
        samp_df = pd.read_csv(filenames[[*filenames.keys()][0]][0])
        # get phenotype data
        zs = np.stack([np.array([float(n) for n in val.lstrip(
            '[').rstrip(']').split(', ')]) for val in samp_df['z']])
        # (i.e., min diff between phenotypes)
        brk = get_min_pw_diff(zs.ravel())

        # make the heatmap array
        brks = np.arange(0, 1+brk, brk)
        arr = np.zeros([len(brks)]*2)

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
            for i, brk_trt0 in enumerate(brks):
                for j, brk_trt1 in enumerate(brks):
                    ct = sum((zs[:, 0] == brk_trt0) * (zs[:, 1] == brk_trt1))
                    arr[i, j] += ct

        arr = arr/np.sum(arr)
        assert np.allclose(np.sum(arr), 1)

        # create labels array
        annot = np.array(['%0.3f' % n for n in arr.ravel()]).reshape(arr.shape)

        # plot the heatmap
        sns.heatmap(arr,
                    #vmin=0,
                    vmax=vmaxes[linkage][genicity],
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

        # set ticks and ticklabels and axis labels
        ticks = np.linspace(0, len(brks)- 1, len(brks), dtype=np.int)
        ticks = ticks + 0.5 # shift ticks to centers of heatmap cells
        ax.set_xticks(ticks)
        ax.set_xticklabels([str(brk) for brk in brks],
                           fontdict={'fontsize':ticklab_fontsize})
        ax.set_xlabel('shifting trait',
                      fontdict={'fontsize': axislab_fontsize})
        if time_step_n == 0:
            ax.set_yticks(ticks)
            ax.set_yticklabels([str(brk) for brk in brks],
                               fontdict={'fontsize':ticklab_fontsize})
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
                        bottom=0.11,
                        right=0.99,
                        top=0.90,
                        wspace=0.05,
                        hspace=None)

    # return fig
    return fig


# produce plots for all scenarios
for linkage in ['independent', 'weak', 'strong']:
    for genicity in [2, 4, 10, 20, 50, 100]:
        dirname_patt = 'mod-non-null_L%s_G%i_its0_' % (linkage, genicity)
        dirs = os.listdir('.')
        candidate_dirs = [d for d in dirs if re.search(dirname_patt, d)]
        if len(candidate_dirs) > 0:
            # make the fig
            fig = plot_phenotypic_shift(linkage, genicity)
            # save the fig
            fig.savefig('phenotypic_shift_L%s_G%i.png' % (linkage, genicity),
                        dpi=dpi,
                        orientation='landscape',
                       )

plt.show()
