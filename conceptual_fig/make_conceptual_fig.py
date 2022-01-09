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
fig_width = 10.5
fig_height = 5.6
dpi = 400
n_ticklabels = 5
gamete_linestyle = '-k'
gamete_linewidth = 1.5
gamete_linealpha = 1
recomb_linestyle = ':k'
recomb_linewidth = 4
recomb_linealpha = 0.5
subplots_adj_left=0.08
subplots_adj_bottom=0.13
subplots_adj_right=0.96
subplots_adj_top=0.85
subplots_adj_wspace=0.14
subplots_adj_hspace=0.05
min_x=0
max_x=1
min_y=0
max_y=1
savefig=True
savefig=True
x_buff=max_x/20
orientation='landscape'
savefig=True

# lists of all possible linkage and genicity values
linkages = {
    'independent': 0.5,
    'weak': 0.05,
    'strong': 0.005
           }
genicities = {
    'low': 4,
    'moderate': 20,
    'high': 100
}

# get the horizontal values used in the non-shifting landscape layer
# and in the before- and after-shifting shifting layer,
# then set as dict of x and y values for plotting horizontal expectation lines
stable = np.linspace(0, 1, 50)[[0,-1]]
b4 = np.linspace(0, 1, 50)[[0,-1]]
af = np.linspace(0.5, 1, 50)[[0,-1]]
expec_lines = {
    'before': (b4, stable),
    'after': (af, stable),
              }

# marker colors
trait_colors = {
    'shifting': ['#29b8ff', '#ff2b56'],
    'stable': ['#000000', '#ffffff'],
              }


# function to calculate phenotypic value of a set of allele colors, and return
# corresponding color on linear color ramp
def get_phenotype(cols):
    effect_size = 1/(len(cols)/2)
    trait_labels = {0: 'shifting', 1: 'stable'}
    zs = {}
    for trt_n in [0, 1]:
        z = 0.5
        for n, col in enumerate(cols[trt_n::2]):
            if n %2:
                effect_size_factor = -1
            else:
                effect_size_factor = +1
            z += effect_size_factor * effect_size * np.where(np.array(
                            trait_colors[trait_labels[trt_n]])==col)[0][0]
        zs[trt_n] = z
    return zs


# marker types
trait_markers = {}
for col in trait_colors['shifting']:
    trait_markers[col] = '*'
for col in trait_colors['stable']:
    trait_markers[col] = 'o'

# marker sizes
trait_sizes = {}
for col in trait_colors['shifting']:
    trait_sizes[col] = {4: 100,
                        20: 20,
                        100: 4,
                       }
for col in trait_colors['stable']:
    trait_sizes[col] = {4: 100,
                        20: 20,
                        100: 4,
                       }

# marker edgewidths
trait_linewidths = {}
for col in trait_colors['shifting']:
    trait_linewidths[col] = {4: 0,
                             20: 0,
                             100: 0,
                            }
for col in trait_colors['stable']:
    trait_linewidths[col] = {4: 0.2,
                             20: 0.2,
                             100: 0.2,
                            }


def get_gamete_params(genicity,
                      min_x=0, max_x=1,
                      inv=False):
    allele_colors = [col for col_list in [*zip(trait_colors['shifting'],
            trait_colors['stable'])] for col in col_list] * (int(genicity/2))
    allele_xs = np.linspace(min_x, max_x, int(genicity*2))
    if inv:
        allele_colors = allele_colors[2:] + allele_colors[:2]
    assert len(allele_colors) == len(allele_xs)
    return allele_colors, allele_xs


def get_recomb_gamete_params(genicity, r,
                             min_x=0, mid_x=0.5, max_x=1,
                             min_y=0, mid_y=0.5, max_y=1,
                             x_buff = 0.5):
    # get x interval
    delta_x = (mid_x-x_buff-min_x)/(genicity*2-1)
    # get regular gamete params
    cols, xs = get_gamete_params(genicity, min_x, mid_x-x_buff)
    # get inverted params
    inv_cols, xs = get_gamete_params(genicity, min_x, mid_x-x_buff,
                                     inv=True)
    # simulate recombination
    recombined_idxs = np.cumsum(np.random.binomial(1, r, len(cols)))%2
    recomb_cols = [{0:cols, 1:inv_cols}[idx][i] for i,
                                        idx in enumerate(recombined_idxs)]
    recomb_xs = xs + mid_x + x_buff*2
    recomb_track_xs = []
    recomb_track_ys = []
    for i, x, idx in zip(range(len(recombined_idxs)), xs, recombined_idxs):
        if i == 0:
            recomb_track_xs.append(x)
            recomb_track_ys.append({0: min_y, 1: max_y}[idx])
        if i > 0:
            recomb_track_xs.append(x - (0.5*delta_x))
            recomb_track_ys.append({0: min_y, 1: max_y}[recombined_idxs[i-1]])
            if idx != recombined_idxs[i-1]:
                recomb_track_xs.append(x - (0.5*delta_x))
                recomb_track_ys.append({0: min_y, 1: max_y}[idx])
        if i == len(recombined_idxs)-1:
            recomb_track_xs.append(x)
            recomb_track_ys.append({0: min_y, 1: max_y}[idx])

    assert len(cols) == len(inv_cols) == len(recomb_cols)
    assert len(xs) == len(recomb_xs)
    assert len(recomb_track_xs) == len(recomb_track_ys)

    return (cols, xs,
            inv_cols,
            recomb_cols, recomb_xs,
            recomb_track_xs, recomb_track_ys)


def plot_one_gamete(genicity, ax=None,
                    xs=None, cols=None,
                    min_x=0, max_x=1, y=0.5,
                    edgecolor='black'):
    if ax is None:
        no_ax = True
        fig,ax = plt.subplots()
    else:
        no_ax = False
    if xs is None and cols is None:
        cols, xs =get_gamete_params(genicity, min_x, max_x)
    else:
        pass
    ax.plot([min_x, max_x], [y, y], gamete_linestyle,
            linewidth=gamete_linewidth,
            alpha=gamete_linealpha)
    for x, col in zip(xs, cols):
        ax.scatter(x, y, c=col, marker=trait_markers[col],
                   s=trait_sizes[col][genicity],
                   edgecolors=edgecolor,
                   linewidths=trait_linewidths[col][genicity],
                  )
    if no_ax:
        plt.show()


def plot_recombining_gametes(genicity, r, ax=None,
                             min_x=0, max_x=1,
                             min_y=0, max_y=1,
                             x_buff=0.05):
    if ax is None:
        no_ax = True
        fig,ax = plt.subplots()
    else:
        no_ax = False
    # get x and y midpoints
    mid_x = min_x + ((max_x - min_x)/2)
    mid_y = min_y + ((max_y - min_y)/2)
    (cols, xs,
     inv_cols,
     recomb_cols, recomb_xs,
     recomb_track_xs, recomb_track_ys) = get_recomb_gamete_params(genicity, r,
                                                    min_x=min_x,
                                                    mid_x=mid_x,
                                                    max_x=max_x,
                                                    min_y=min_y,
                                                    mid_y=mid_y,
                                                    max_y=max_y,
                                                    x_buff=x_buff)
    plot_one_gamete(genicity, ax, xs, cols,
                    min_x=min_x, max_x=mid_x-x_buff, y=min_y)
    plot_one_gamete(genicity, ax, xs, inv_cols,
                    min_x=min_x, max_x=mid_x-x_buff, y=max_y)
    plot_one_gamete(genicity, ax, recomb_xs, recomb_cols,
                    min_x=mid_x+(2*x_buff), max_x=max_x, y=mid_y)
    ax.plot(recomb_track_xs, recomb_track_ys, recomb_linestyle,
            linewidth=recomb_linewidth, alpha=recomb_linealpha)

    print(get_phenotype(recomb_cols))

    if no_ax:
        plt.show()


def plot_concept(min_x=0, max_x=1, min_y=0, max_y=1, x_buff=0.05, savefig=True):

    # create the figure
    fig = plt.figure(dpi=dpi,
                     figsize=(fig_width, fig_height)
                    )
    #fig.suptitle(genicity, fontdict={'fontsize': suptitle_fontsize})
    gs = fig.add_gridspec(nrows=3, ncols=3, width_ratios=[1,1,1])

    i = 0
    j = 0
    for linkage, linkage_n in linkages.items():
        for genicity, genicity_n in genicities.items():
            print(i, j)
            ax = fig.add_subplot(gs[i, j%3])

            # get rid of the ticks and labels
            ax.set_xticks([])
            ax.set_xticklabels([])
            ax.set_yticks([])
            ax.set_yticklabels([])

            # plot the scenario
            plot_recombining_gametes(genicity=genicity_n, r=linkage_n, ax=ax,
                                     min_x=min_x, max_x=max_x,
                                     min_y=min_y, max_y=max_y, x_buff=x_buff)

            # other plot controls
            #ax.invert_yaxis()
            ax.tick_params(axis=u'both', which=u'both',length=0) # delete ticks

            j+=1
        i+=1

    # adjust suplot spacing
    plt.subplots_adjust(left=subplots_adj_left,
                        bottom=subplots_adj_bottom,
                        right=subplots_adj_right,
                        top=subplots_adj_top,
                        wspace=subplots_adj_wspace,
                        hspace=subplots_adj_hspace)
    if savefig:
        fig.savefig(('/home/deth/Desktop/CAL/research/projects/sim/ch2/'
                     'climate_change_adaptation_and_genomic_arch/conceptual_fig'
                     '/conceptual_fig.png'),
                    dpi=dpi, orientation=orientation)
    # return fig
    return fig


if __name__ == '__main__':
    plot_concept(min_x=min_x,
                 max_x=max_x,
                 min_y=min_y,
                 max_y=max_y,
                 x_buff=x_buff,
                 savefig=savefig
                )
    plt.close('all')
