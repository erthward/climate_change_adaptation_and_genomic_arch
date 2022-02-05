import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from matplotlib.patches import Polygon as mplPolygon
from collections import Counter as C
from palettable.cartocolors.sequential import PinkYl_7
from palettable.scientific.sequential import Acton_20_r
from shapely.geometry import Polygon as shapelyPolygon
import statsmodels.api as sm
import sys
import re
import os

"""
TODO:

- just make the before-after change highlighted by gray->red edgecolor on a gray dist?
- decide on axis order, then align with the pheno-shift fig axes

"""


# plot params
suptitle_fontsize = 50
title_fontsize = 40
contour_axislab_fontsize = 10
contour_ticklab_fontsize = 7
annot_fontsize = 14
cbar_fontsize = 14
fig_width = 3
fig_height = 5
dpi = 400
n_ticklabels = 5
contour_alpha = 0.5
contour_linewidth = 0.1
contour_linecolor = 'gray'
contour_lims = [0,1]
vert_connector_linewidth = 0.25
subplots_adj_left=0.01
subplots_adj_bottom=0.01
subplots_adj_right=0.99
subplots_adj_top=0.99
subplots_adj_wspace=0.01
subplots_adj_hspace=0.01
min_x=0
max_x=1
min_y=0
max_y=1
savefig=True
savefig=True
x_buff=max_x/20
orientation='landscape'
savefig=True


fig = plt.figure(dpi=dpi,
                 figsize=(fig_width, fig_height))

gs = fig.add_gridspec(nrows=8, ncols=3)

# functional form of fitness kernel (assumed isotropic in 2d space)
# NOTE: NOT ACTUALLY THE FUNCTIONAL FORM USED IN THE MODEL,
#       BECAUSE gamma=1 IN OUR SIMS, SO FITNESS DECREASES LINEARLY;
#       JUST A CONCEPTUAL FRAMEWORK

def get_fit(max_fit_z, actual_z, σ=0.025, c=3):
    max_fit_z = np.array(max_fit_z).reshape((2,1))
    actual_z = np.array(actual_z).reshape((2,1))
    Σ = np.array([[σ, 0], [0, σ]])
    diff_vec = actual_z - max_fit_z
    fit = (1/(np.sqrt(2*np.pi*Σ)))*np.exp(
        (-1/2)*np.matmul(np.matmul((diff_vec).T, np.linalg.inv(Σ)), diff_vec))
    fit = fit[0,0]
    #fit = c * np.exp((-(((max_fit_z[0]-actual_z[0])**2)/(2*(sigma**2)))) +
    #                   (((max_fit_z[1]-actual_z[1])**2)/2*(sigma**2)))
    #fit = ((1-np.abs(max_fit_z[0]-actual_z[0])**2)*
           #(1-np.abs(max_fit_z[1]-actual_z[1])**2))
    return fit


# figure to get surface values for a fitness peak centered on max_fit
def get_fit_surf(max_fit, sigma=0.2):
    x = np.linspace(max_fit[0]-(2*sigma), max_fit[0]+(2*sigma))
    y = np.linspace(max_fit[1]-(2*sigma), max_fit[1]+(2*sigma))
    X, Y = np.meshgrid(x, y)
    Z = np.ones(X.shape)*np.nan
    for i in range(X.shape[0]):
        for j in range(X.shape[1]):
            fit = get_fit(max_fit, [X[i,j], Y[i,j]])
            Z[i,j] = fit
    # normalize 0 to 1, and 'mask' ~0 values out
    Z = (Z-np.min(Z))/(np.max(Z)-np.min(Z))
    for i in range(Z.shape[0]):
        for j in range(Z.shape[1]):
            if np.allclose(0, Z[i,j], atol=1e-1):
                Z[i,j] = np.nan
    return X, Y, Z


# get rid of ticks and ticklabels
def delete_ticks(ax):
    ax.set_xticks(())
    ax.set_xticklabels(())
    ax.set_yticks(())
    ax.set_yticklabels(())
    try:
        ax.set_zticks(())
        ax.set_zticklabels(())
    except Exception:
        pass
    return


# get the horizontal values used in the non-shifting landscape layer
# and in the before- and after-shifting shifting layer,
# then set as dict of x and y values for plotting horizontal expectation lines
stable = np.linspace(1, 0, 52)
b4 = np.linspace(1, 0, 52)
af = np.linspace(1, 0.5, 52)
# stack trt_1_b4, trt_2_b4, trt_1_af, trt_2_af
time_diff_landscape = np.vstack([b4 for _ in range(13)] +
                                [af for _ in range(13)] +
                                [b4 for _ in range(26)])

# colors for before- and after-climate change surfaces
highlight_col = '#ffe436'
contcolors = {'b4': 'gray',
              'af': 'white',
             }
contalphas = {'b4': 0.8,
              'af': 0.2,
             }
edgecolors = {'b4': 'gray',
              'af': highlight_col,
             }

# PART I:
# plot the before/after-split landscape at top
ax = fig.add_subplot(gs[:3, :])
im = ax.imshow(time_diff_landscape, cmap='gray')
time_diff_landscape[:13, :] = np.nan
time_diff_landscape[26:39, :] = np.nan
im = ax.imshow(time_diff_landscape, cmap='bwr')
ax.plot([0,52], [26,26], '--k', linewidth=2)
ax.plot([0,52], [13,13], '--k', linewidth=1)
ax.plot([0,52], [39,39], '--k', linewidth=1)
ax.plot([0.25,2,2,0.25,0.25], [0.25,0.25,49.75,49.75,0.25], ':r', linewidth=0.5)
ax.plot([24,26,26,24,24], [0.25,0.25,49.75,49.75,0.25], ':r', linewidth=0.5)
ax.plot([48,49.75,49.75,48,48], [0.25,0.25,49.75,49.75,0.25], ':r', linewidth=0.5)
ax.set_xlim([0,50])
ax.set_ylim([0,50])
#plt.colorbar(im)
delete_ticks(ax)


# PART II:
positions_and_fits = {0.5:{'b4':(1,1),
                           'af':(1,1),
                          },
                      24.5:{'b4':(0.5,0.5),
                            'af':(0.75,0.5),
                           },
                      49.5:{'b4':(0,0),
                            'af':(0.5,0),
                           },
                     }
# plot before-after fitness landscapes at ends and middle of landscape
for ax_i, x_pos in enumerate([0.5, 24.5, 49.5]):
    ax = fig.add_subplot(gs[3:5, ax_i], projection='3d')
    #ax.set_aspect('equal')
    ax.invert_yaxis()
    #ax.invert_xaxis()
    ax._axis3don = False
    ax.set_zticks(())
    ax.set_zticklabels(())
    #ax.set_ylabel('shifting', fontdict={'fontsize': contour_axislab_fontsize})
    #ax.set_xlabel('stable', fontdict={'fontsize': contour_axislab_fontsize})
    #ax.set_xticks([0,0.25,0.5,0.75,1])
    #ax.set_yticks([0,0.25,0.5,0.75,1])
    ax.set_ylim([1,0])
    ax.set_xlim(contour_lims)
    ax.set_zlim(contour_lims)
    ax.grid(False)
    ax.view_init(60, 90)

    # get fitness values, based on landscape positions
    fits = positions_and_fits[x_pos]

    # plot manual gridding and axes
    ax.plot([0,1,1,0,0], [0,0,1,1,0], [0,0,0,0,0], 'k', linewidth=0.5)
    ax.plot([fits['b4'][0]]*2, [fits['b4'][1]]*2, [0,1], 'k',
            linewidth=vert_connector_linewidth)
    ax.plot([fits['af'][0]]*2, [fits['af'][1]]*2, [0,1], 'k',
            linewidth=vert_connector_linewidth)
    for grid_val in [0.25, 0.5, 0.75]:
        ax.plot([grid_val, grid_val], [0,1], [0,0], ':', color='gray',
                linewidth=0.1, alpha=0.65)
        ax.plot([0,1], [grid_val, grid_val], [0,0], ':', color='gray',
                linewidth=0.1, alpha=0.65)

    # plot before- and after-climate change fitness surfs
    b4_surf = get_fit_surf(fits['b4'])
    af_surf = get_fit_surf(fits['af'])
    ax.plot_surface(*b4_surf, color=contcolors['b4'],
                    edgecolor=edgecolors['b4'],
                    alpha=contalphas['b4'],
                    linewidth=0.05)
    ax.plot_surface(*af_surf, color=contcolors['af'],
                    edgecolor=edgecolors['af'],
                    alpha=contalphas['af'],
                    linewidth=0.05)
fig.show()


# PART III:
# plot before-after comparison of fitness 'ridge' across real landscape
ax = fig.add_subplot(gs[5:,:], projection='3d')
ax._axis3don = False
ax.invert_yaxis()
#ax.invert_xaxis()
#ax.set_aspect('equal')
ax.view_init(80, 90)
#ax.set_ylabel('shifting', fontdict={'fontsize': contour_axislab_fontsize})
#ax.set_xlabel('stable', fontdict={'fontsize': contour_axislab_fontsize})
ax.set_ylim([1,0])
ax.set_xlim(contour_lims)
ax.set_zlim(contour_lims)
b4_ys = np.stack([b4]*25)
b4_xs = np.linspace(b4-0.05, b4+0.05, 25)
af_xs = np.linspace(af-0.05, af+0.05, 25)
b4_fits = np.ones(af_xs.shape)*np.nan
af_fits = np.ones(af_xs.shape)*np.nan
for i in range(b4_fits.shape[0]):
    for j in range(b4_fits.shape[1]):
        b4_fit = get_fit((b4_xs[13,j], b4_ys[i, j]),
                         (b4_xs[i,j], b4_ys[i,j]),
                         σ=0.01)
        af_fit = get_fit((af_xs[13,j], b4_ys[i, j]),
                         (af_xs[i,j], b4_ys[i,j]),
                         σ=0.01)
        b4_fits[i,j] = b4_fit
        af_fits[i,j] = af_fit
b4_fits = (b4_fits - np.min(b4_fits))/(np.max(b4_fits)-np.min(b4_fits))
af_fits = (af_fits - np.min(af_fits))/(np.max(af_fits)-np.min(af_fits))


# plot manual axes and grid and before and after lines, below the surfaces
# axes
ax.plot([0,1,1,0,0], [0,0,1,1,0], [0,0,0,0,0], 'k', linewidth=0.5)
# before and after lines
ax.plot([0,1], [0,1], [0,0], edgecolors['b4'], linewidth=1)
ax.plot([0.5,1], [0,1], [0,0], edgecolors['af'], linewidth=1)
# vertical connectors
ax.plot([0,0], [0,0], [0,1], 'k', linewidth=vert_connector_linewidth, alpha=0.5)
ax.plot([0.5,0.5], [0.5,0.5], [0,1], 'k', linewidth=vert_connector_linewidth, alpha=0.5)
ax.plot([1,1], [1,1], [0,1], 'k', linewidth=vert_connector_linewidth, alpha=0.5)
ax.plot([0,0], [0,0], [0,1], 'k', linewidth=vert_connector_linewidth, alpha=0.5)
ax.plot([0.5,0.5], [0,0], [0,1], 'k', linewidth=vert_connector_linewidth, alpha=0.5)
ax.plot([0.75,0.75], [0.5,0.5], [0,1], 'k', linewidth=vert_connector_linewidth, alpha=0.5)
# grid lines
for grid_val in [0.25, 0.5, 0.75]:
        ax.plot([grid_val, grid_val], [0,1], [0,0], color='gray',
                linewidth=0.1, alpha=0.65)
        ax.plot([0,1], [grid_val, grid_val], [0,0], color='gray',
                linewidth=0.1, alpha=0.65)

ax.plot_surface(b4_xs, b4_ys, b4_fits, color=contcolors['b4'],
                edgecolor=edgecolors['b4'],
                linewidth=0.05,
                alpha=contalphas['b4'])
ax.plot_surface(af_xs, b4_ys, af_fits, color=contcolors['af'],
                edgecolor=edgecolors['af'],
                linewidth=0.05,
                 alpha=contalphas['af'])


plt.subplots_adjust(left=subplots_adj_left,
                    bottom=subplots_adj_bottom,
                    right=subplots_adj_right,
                    top=subplots_adj_top,
                    wspace=subplots_adj_wspace,
                    hspace=subplots_adj_hspace)


fig.savefig('ch2_conceptual_fig_raw.png',
            dpi=dpi, orientation='portrait')
