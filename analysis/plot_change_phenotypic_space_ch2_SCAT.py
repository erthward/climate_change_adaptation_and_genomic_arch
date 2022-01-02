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
fig_width = 14.5
fig_height = 5.6
dpi = 400
n_ticklabels = 5
linewidth = 1
linealpha = 0.8
marker_sizes = {4: 500,
                20: 40,
                100: 3,
               }

# data directory
if os.getcwd().split('/')[1] == 'home':
    datadir = '/home/deth/Desktop/tmp_ch2_stats_tests_dev/'
else:
    datadir = '/global/scratch/users/drewhart/ch2/output/analysis'

# lists of all possible linkage and genicity values
linkages = ['independent', 'weak', 'strong']
genicities =  [2, 4, 10, 20, 50, 100]

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


def plot_phenotypic_shift(linkage, genicity, fix_ur_corner=True):

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

        xs = []
        ys = []

        polys = []
        areas = []

        # loop through all the sims' directories
        for sim_dirname, sim_filenames in filenames.items():
            print('\nNOW PROCESSING FILES IN: %s\n\n' % sim_dirname)
            df = pd.read_csv(sim_filenames[time_step_n])

            # get phenotype data
            zs = np.stack([np.array([float(n) for n in val.lstrip(
                '[').rstrip(']').split(', ')]) for val in df['z']])
            curr_xs = [*zs[:, 0]]
            xs.extend(curr_xs)
            curr_ys = [*zs[:, 1]]
            ys.extend(curr_ys)

            # TODO: get polygon area for this iteration
            if time_step > 2700:
                if not fix_ur_corner:
                    X = sm.add_constant(np.vstack(curr_xs))
                    y = np.vstack(curr_ys)
                else:
                    # NOTE: subtract 1 from x values, and add no constant term, 
                    #       then later add 1 again to trasnlate back, to be able to
                    #       fix upper-left corner at (1,1)
                    X = np.vstack(curr_xs) - 1
                    y = np.vstack(curr_ys) - 1
                mod = sm.OLS(y, X).fit()
                pred_xs = [0, 1]
                if fix_ur_corner:
                    pred_xs = [pred_x-1 for pred_x in pred_xs]
                    pred_ys = mod.predict(pred_xs)
                    pred_xs = [pred_x +1 for pred_x in pred_xs]
                    pred_ys = [pred_y +1 for pred_y in pred_ys]
                else:
                    # NOTE: need to add in the constant term to get predictions
                    pred_ys = mod.predict(np.array([*zip([1,1], pred_xs)]))
                poly_xs = [*pred_xs, expec_lines[time_step][0][0], pred_xs[0]]
                poly_ys = [*pred_ys, expec_lines[time_step][1][0], pred_ys[0]]
                poly = shapelyPolygon([*zip(poly_xs, poly_ys)])
                area = poly.area
                polys.append(poly)
                areas.append(area)

        # plot the scatter
        uniqs = [*set([*zip(xs, ys)])]
        cts = C([*zip(xs, ys)])
        sizes = [cts[uniq] for uniq in uniqs]
        sizes = [(s-np.min(sizes))/(np.max(sizes)-np.min(sizes)) for s in sizes]
        max_alpha = 0.9
        alphas = [s*max_alpha for s in sizes]
        max_size = 25
        sizes = np.array(sizes) * max_size
        plot_color = '#000000'
        for coords, size, alpha in zip(uniqs, sizes, alphas):
            x, y = coords
            ax.scatter(x, y,
                       s=marker_sizes[genicity],
                       c=plot_color,
                       alpha=alpha,
                       edgecolors='none')

        # fit linear regression, if in second or third time point
        if time_step > 2500:
            if not fix_ur_corner:
                X = sm.add_constant(np.vstack(xs))
                y = np.vstack(ys)
            else:
                # NOTE: subtract 1 from x values, and add no constant term, 
                #       then later add 1 again to trasnlate back, to be able to
                #       fix upper-left corner at (1,1)
                X = np.vstack(xs) - 1
                y = np.vstack(ys) - 1
            mod = sm.OLS(y, X).fit()
            pred_xs = [0, 1]
            if fix_ur_corner:
                pred_xs = [pred_x-1 for pred_x in pred_xs]
                pred_ys = mod.predict(pred_xs)
                pred_xs = [pred_x +1 for pred_x in pred_xs]
                pred_ys = [pred_y +1 for pred_y in pred_ys]
            else:
                # NOTE: need to add in the constant term to get predictions
                pred_ys = mod.predict(np.array([*zip([1,1], pred_xs)]))
            ax.plot(pred_xs, pred_ys, ':', color=plot_color, alpha=linealpha,
                    linewidth=linewidth)

        # add expectation and trend lines, as needed
        if time_step > 2500:
            # calculate area between trend line and expectation line
            # (i.e. 'phenotypic undershoot')
            poly_xs = [*pred_xs, expec_lines[time_step][0][0], pred_xs[0]]
            poly_ys = [*pred_ys, expec_lines[time_step][1][0], pred_ys[0]]
            poly = shapelyPolygon([*zip(poly_xs, poly_ys)])
            area = poly.area
            # save the area for output, if final timestep
            # and plot the area
            poly = mplPolygon([*zip(poly_xs, poly_ys)])
            poly_color = '#f24156'
            poly.set_color(poly_color)
            poly.set_alpha(0.5)
            ax.add_patch(poly)
            # expectation line
            ax.plot(*expec_lines[time_step], '-k', alpha=linealpha,
                    linewidth=linewidth)
        # expectation line
        ax.plot(*expec_lines[2499], '-k', alpha=linealpha, linewidth=linewidth)


        # set ticks and ticklabels and axis labels
        ticks = np.linspace(0, 1, n_ticklabels)
        ticklabels = np.linspace(0, 1, n_ticklabels)
        ax.set_xticks(ticks)
        ax.set_xticklabels(ticklabels, fontdict={'fontsize':ticklab_fontsize})
        ax.set_xlim([0,1])
        ax.set_ylim([0,1])
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
        #ax.invert_yaxis()
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
    return fig, areas


# dict to store pheno undershoot area data
pheno_undershoot_dict = {'linkage': [], 'genicity': [], 'undershoot': []}
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
            fig, undershoot = plot_phenotypic_shift(linkage, genicity)
            # save the fig
            fig.savefig(os.path.join(datadir,
                        'phenotypic_shift_L%s_G%s_SCAT.png' % (linkage,
                                                    str(genicity).zfill(2))),
                        dpi=dpi,
                        orientation='landscape',
                       )
            # save the undershoot data
            pheno_undershoot_dict['linkage'].extend([linkage]*len(undershoot))
            pheno_undershoot_dict['genicity'].extend([genicity]*len(undershoot))
            pheno_undershoot_dict['undershoot'].extend(undershoot)

# analyze undershoot data
linkage_dict = {'independent': 0.5, 'weak': 0.05, 'strong': 0.005}
pheno_undershoot_df = pd.DataFrame.from_dict(pheno_undershoot_dict).replace(
    {'linkage': linkage_dict})
print(pheno_undershoot_df)
pheno_undershoot_df.to_csv(os.path.join(datadir,
                                        'phenotypic_shift_undershoot.csv'),
                           index=False)
#pheno_undershoot_df['intxn'] = (pheno_undershoot_df.linkage *
#                                pheno_undershoot_df.genicity)

y = pheno_undershoot_df['undershoot']
X = sm.add_constant(pheno_undershoot_df[['linkage', 'genicity']])#, 'intxn']])
mod = sm.OLS(y, X).fit()
print(mod.summary())


