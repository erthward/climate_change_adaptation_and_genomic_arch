import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from matplotlib.patches import FancyArrowPatch
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


save_it = True


# plot params
suptitle_fontsize = 50
title_fontsize = 40
contour_axislab_fontsize = 10
contour_ticklab_fontsize = 7
annot_fontsize = 14
cbar_fontsize = 14
fig_width = 5
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


# some 3D arrow code stolen from:
    #https://stackoverflow.com/questions/22867620/putting-arrowheads-on-vectors-in-matplotlibs-3d-plot
class Arrow3D(FancyArrowPatch):

    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = mplot3d.proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)


# functional form of fitness kernel (assumed isotropic in 2d space)
# NOTE: NOT ACTUALLY THE FUNCTIONAL FORM USED IN THE MODEL,
#       BECAUSE gamma=1 IN OUR SIMS, SO FITNESS DECREASES LINEARLY;
#       JUST A CONCEPTUAL FRAMEWORK
def get_fit(max_fit_z, actual_z, σ=0.005, c=3):
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
def get_fit_surf(max_fit, sigma=0.1):
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


fig = plt.figure(dpi=dpi,
                 figsize=(fig_width, fig_height))

gs = fig.add_gridspec(nrows=9, ncols=5,
                      height_ratios=[1,1,1,1,1.5,1,1,1,1])

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
ax = fig.add_subplot(gs[:4, 1:4])
im = ax.imshow(time_diff_landscape, cmap='gray')
time_diff_landscape[:13, :] = np.nan
time_diff_landscape[26:39, :] = np.nan
im = ax.imshow(time_diff_landscape, cmap='bwr')
#ax.plot([0,52], [26,26], '--k', linewidth=2)
ax.plot([0,52], [12.5,12.5], '--k', linewidth=1)
ax.plot([0,52], [38.5,38.5], '--k', linewidth=1)
ax.plot([0.25,2,2,0.25,0.25], [26.25,26.25,49.75,49.75,26.25], '-',
        color='#555555', linewidth=1)
ax.plot([24,26,26,24,24], [26.25,26.25,49.75,49.75,26.25], '-',
        color='#555555', linewidth=1)
ax.plot([48,49.75,49.75,48,48], [26.25,26.25,49.75,49.75,26.25], '-',
        color='#555555', linewidth=1)
ax.plot([0.25,2,2,0.25,0.25], [0.25,0.25,24.75,24.75,0.25], '-',
        color=highlight_col, linewidth=1)
ax.plot([24,26,26,24,24], [0.25,0.25,24.75,24.75,0.25], '-',
        color=highlight_col, linewidth=1)
ax.plot([48,49.75,49.75,48,48], [0.25,0.25,24.75,24.75,0.25], '-',
        color=highlight_col, linewidth=1)

ax.set_xlim([-0.5, 50.5])
ax.set_ylim([-0.5, 50.5])
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
ax = fig.add_subplot(gs[5:, :], projection='3d')
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
ax.view_init(55, 90)

# add manual colormaps along axes
manual_cmap_arr_xs = np.linspace(0, 1+(15*(1/200)), 210)
manual_cmap_arr_ys = np.linspace(0-(20*(1/200)), 1, 210)
manual_cmap_arr_X, manual_cmap_arr_Y = np.meshgrid(manual_cmap_arr_xs,
                                                   manual_cmap_arr_ys)
manual_cmap_arr = np.ones((210,210))*np.nan
manual_cmap_arr[0:10, :200] = np.linspace(0,1,200)
for j in range(200, 210):
    manual_cmap_arr[10:, j] = np.linspace(0,1,200)
cmap_arr_img = np.ones((210,210,4))*np.nan

for i in range(0, 10):
    for j in range(0, 194):
        cmap_arr_img[i,j,:] = plt.cm.bwr(manual_cmap_arr[i,j])
for i in range(21, 210):
    for j in range(200, 210):
        cmap_arr_img[i,j,:] = plt.cm.gray(manual_cmap_arr[i,j])
ax.plot_surface(manual_cmap_arr_X, manual_cmap_arr_Y, np.zeros((210, 210)),
                facecolors=cmap_arr_img, shade=False)

# plot manual gridding and axes
ax.plot([0,1,1,0,0], [0,0,1,1,0], [0,0,0,0,0], 'k', linewidth=0.5)
for grid_val in [0.25, 0.5, 0.75]:
    ax.plot([grid_val, grid_val], [0,1], [0,0], ':', color='gray',
            linewidth=0.1, alpha=0.8)
    ax.plot([0,1], [grid_val, grid_val], [0,0], ':', color='gray',
            linewidth=0.1, alpha=0.8)


for ax_i, x_pos in enumerate([0.5, 24.5, 49.5]):

    # get fitness values, based on landscape positions
    fits = positions_and_fits[x_pos]

    # plot before- and after-climate change fitness surfs
    b4_surf = get_fit_surf(fits['b4'])
    af_surf = get_fit_surf(fits['af'])
    ax.plot([fits['b4'][0]]*2, [fits['b4'][1]]*2, [0,1], 'k',
            linewidth=vert_connector_linewidth)
    ax.plot([fits['af'][0]]*2, [fits['af'][1]]*2, [0,1], 'k',
            linewidth=vert_connector_linewidth)
    ax.plot_surface(*b4_surf, color=contcolors['b4'],
                    edgecolor=edgecolors['b4'],
                    alpha=contalphas['b4'],
                    linewidth=0.05)
    ax.plot_surface(*af_surf, color=contcolors['af'],
                    edgecolor=edgecolors['af'],
                    alpha=contalphas['af'],
                    linewidth=0.05)
    # plot change arrows
    #if ax_i > 0:
    #    arw = Arrow3D([np.mean(b4_surf[0]), np.mean(af_surf[0])],
    #                [np.mean(b4_surf[1]), np.mean(af_surf[1])],
    #                [1,1],
    #                mutation_scale=20, lw=3, arrowstyle="-|>", color="k")
    #    ax.add_artist(arw)

# before and after lines
ax.plot([0,1], [0,1], [0,0], '#555555', linewidth=1)
ax.plot([0.5,1], [0,1], [0,0], edgecolors['af'], linewidth=1)
fig.subplots_adjust(left=subplots_adj_left,
                    bottom=subplots_adj_bottom,
                    right=subplots_adj_right,
                    top=subplots_adj_top,
                    wspace=subplots_adj_wspace,
                    hspace=subplots_adj_hspace)


fig.show()

if save_it:
    fig.savefig('ch2_conceptual_fig_RAW.png',
                dpi=dpi, orientation='portrait')
