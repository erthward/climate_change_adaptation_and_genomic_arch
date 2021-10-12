#!/usr/bin/python
#batch_script.py

# flake8: noqa

'''
TODO:
    - clean up all the now extremely gross data-structure code I've cobbled
    together
    - run on Savio
    - code up more/better statistical tests
'''

SAVE=True

import matplotlib.pyplot as plt
import matplotlib as mpl
import geonomics as gnx
import pandas as pd
import numpy as np
import scipy.stats
import bisect
import tskit
import copy
import time
import os


#-------------------
# params for figures
#-------------------
# flag indicating whether or not to plot the lineage map
map_gene_flow = True
# set colors
colors = {'null': {'nonneut': '#237ade',  # dark blue
                   'neut': '#94b3d6'},    # light blue
          'non-null': {'nonneut': '#d60957',  # dark rose
                   'neut': '#d98daa'}}    # light rose
suptitle_size = 7 #DETH
rowlab_size = 6 #DETH
collab_size = 6 #DETH
ticklabel_size = 6 #DETH

#dir_tick_locs = np.linspace(0, 360, 9)
dir_tick_locs = np.linspace(0, 360, 5)
#dir_tick_labs = ['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW', 'N']
dir_tick_labs = ['N', 'E', 'S', 'W', 'N']


#------------
# hist figure
#------------
fig_hist = plt.figure('hist')
fig_hist.suptitle(('histograms of gene-flow direction, distance,'
                  ' and speed, across linkage-genicity combinations'),
                  fontsize=suptitle_size) #DETH
gs_hist = mpl.gridspec.GridSpec(3, 9+1) #DETH

for genicity_n, genicity in enumerate([5, 20, 100]):
    for linkage_n, linkage in enumerate(['independent', 'weak', 'strong']):
        for stat_n, stat in enumerate(['dir', 'dist', 'speed']):
            row_idx = linkage_n
            col_idx = 3* genicity_n + stat_n
            ax = fig_hist.add_subplot(gs_hist[row_idx, col_idx])
            if col_idx == 0:
                ax.set_ylabel('linkage:\n%s' % str(linkage), size=rowlab_size)
            if row_idx == 0 and col_idx in [0, 3, 6]:
                ax.set_title('|| genicity: %i' % genicity, size=collab_size) #DETH
            # plot a histogram of the unlinked and then linked data
            for nullness in ['non-null', 'null']:
                for neut_idx, neutrality in enumerate(['neut', 'nonneut']):
                    data = np.random.beta(np.random.uniform(0,10),
                                          np.random.uniform(0,10), 10000)
                    if stat == 'dir':
                        data = data * 360
                    if neutrality == 'neut': #DETH
                        kde = scipy.stats.gaussian_kde(data)
                        xx = np.linspace(min(data), max(data), 1000)
                        kde_vals = kde(xx)
                        vals, breaks, bars = ax.hist(data, bins=50, alpha=0)
                        kde_plot_factor = max(vals)/max(kde_vals)
                        kde_plot_vals = [val * kde_plot_factor for val in kde_vals]
                        ax.plot(xx, kde_plot_vals, alpha=0.5,
                        #sns.kdeplot(data, ax=ax, alpha=0.5,
                                    label='%s: %s' % (nullness, neutrality),
                                    color=colors[nullness][neutrality]) #DETH
                    else:
                        ax.hist(data, bins = 50, alpha=0.5,
                            label= '%s: %s' % (nullness, neutrality),
                            color=colors[nullness][neutrality])
                    ax.set_xlabel(stat, size=8)
                    if row_idx==2 and col_idx==8:
                        ax.legend(prop={'size': 6},
                                  fontsize=5, #DETH
                                  bbox_to_anchor=(1.5,0.2)) # DETH
                    if col_idx in [0, 3, 6]:
                        ax.set_xticks(dir_tick_locs)
                        ax.set_xticklabels(dir_tick_labs)
                    ax.tick_params(labelsize=ticklabel_size) # DETH
                    ax.tick_params(axis='y', rotation=60) # DETH


fig_hist.show()

# save all the figures
# set the spacing between subplots # DETH
plt.subplots_adjust(left=0.1,
                    bottom=0.1,
                    right=0.9,
                    top=0.9,
                    wspace=1.0,
                    hspace=0.9)
plt.show()

if SAVE:
    try:
        rcParams['figure.figsize'] = 40, 12
    except Exception as e:
        pass
    fig_hist.savefig('practicing_hist_fig.png', format='png', dpi=1000)
