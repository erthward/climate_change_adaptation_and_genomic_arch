import numpy as np
from numpy import pi as pi
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import vonmises
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import statsmodels.api as sm
from statsmodels.formula.api import ols
from statsmodels.stats.multicomp import pairwise_tukeyhsd
import os



# NOTE:
    # I use N,E,W,S to talk about directions in this script, for var-naming
    # purposes, but the directions simulated in the study do not truly
    # correspond to N,E,W,S (especially in the sense that climate change
    # simulation is meant to depict climate shifting upslope locally,
    # not northward regionally!)


# TODO:

    # write linear model to look at multivariate response of dists to the
    # scenarios?

    # consider if/how to depict uncertainty envelopes on the dists?



# plot params
col_label_fontsize = 18
row_label_fontsize = 18
fig_width = 13.8
fig_height = 4.6
dpi = 400
# factor by which to multiply the max absolute values in the gray reference
# circles in order to determine the plots' axis limits
axlim_graycirc_factor = 1.7

# data directory
if os.getcwd().split('/')[1] == 'home':
    analysis_dir = '/media/deth/SLAB/ch2/analysis/'
else:
    with open(('/global/scratch/users/drewhart/ch2/climate_change_adaptation_'
               'and_genomic_arch/analysis/analysisdir.txt'), 'r') as f:
        analysis_dir = f.read().strip()
    #analysis_dir = '/global/scratch/users/drewhart/ch2/output/analysis'


def viz_vM_mix_dist(row=None, df=None, mu=None, kappa=None, alpha=None,
                    nullness=None,
                    plot_type='circ', to_deg=True, plot_dirlabels=False,
                    on_curr_ax=False, return_axlims=True, n_rand=10000):
    """
    Visualize a von Mises mixture distribution

    NOTE: n_rand is only going to be approximate (but very close)
    to true total number of draws used to visualize the dist
    """
    # need mu, kappa, alpha, and nullness if row is None
    if row is None and df is None:
        assert (mu is not None and
                kappa is not None and
                alpha is not None and
                nullness is not None), ("MU, KAPPA, ALPHA, AND NULLNESS "
                                        "ARGS MISSING")
        # make all numeric params objects into 2d np objects,
        # so that looping works below
        mu = np.atleast_2d(mu)
        kappa = np.atleast_2d(kappa)
        alpha = np.atleast_2d(alpha)

    # otherwise, grab them from the row
    elif row is not None and df is None:
        df = pd.DataFrame(row).T
        mu = df.iloc[:, -12:-8].values
        kappa = df.iloc[:, -8:-4].values
        alpha = df.iloc[:, -4:].values
        if nullness is None:
            nullness = row['nullness']
        else:
            pass

    # or from the df
    elif row is None and df is not None:
        mu = df.iloc[:, -12:-8].values
        kappa = df.iloc[:, -8:-4].values
        alpha = df.iloc[:, -4:].values
        if nullness is None:
            nulless = df.loc[:, 'nullness']
        else:
            pass

    # create empty lists to be filled with the N/S- and E-facing PDF densities,
    # then returned
    E_densities = []
    NS_densities = []

    # get the current axes, for plotting, if requested
    if on_curr_ax:
        ax = plt.gca()
    else:
        ax=None

    # get plot color, based on nullness
    if nullness == 'null':
        col = '#79c2d9' # blue
    elif nullness == 'non-null':
        col = '#c4626e' # red
    plot_vals = []
    if plot_type == 'hist':
        # for each set of fitted params
        for i in range(mu.shape[0]):
            # determine how many random vals to draw from each dist
            samp_sizes = [int(n) for n in np.floor(np.array(alpha[i,:])*n_rand)]
            # draw vals
            for mix_n, samp_size in enumerate(samp_sizes):
                plot_vals.extend(vonmises.rvs(kappa[i, mix_n],
                                         loc=mu[i, mix_n],
                                         scale=1,
                                         size=samp_size))
                # correct to the interval [0, 2pi]
                plot_vals = [val+(2*pi) if val <0 else val for val in plot_vals]
    elif plot_type == 'circ':
        # get ~3000-ish evenly spaced densities
        dirs = np.arange(0, 2*pi, 0.002)
        # get indices of directions closest to E, N, and S
        E_idx = 0
        N_idx = np.abs(dirs-(pi/2)).argmin()
        S_idx = np.abs(dirs-(3*pi/2)).argmin()
        # numpy array to contain the output mixture-dist densities
        mix_dens = np.zeros(dirs.shape)
        # for each set of fitted params
        for i in range(mu.shape[0]):
            sub_mix_dens = np.zeros(dirs.shape)
            for mix_n in range(mu.shape[1]):
                # calculate densities
                dens = vonmises.pdf(dirs,
                                    kappa=kappa[i, mix_n],
                                    loc=mu[i, mix_n],
                                    scale=1)
                # weight densities by this distribution's alpha
                weight_dens = alpha[i, mix_n] * dens
                # add the vals to the weighted sum
                sub_mix_dens = sub_mix_dens + weight_dens
            # save this mixture distribution's E-facing density
            # and the mean of its N- and S-facing densities
            E_dens = sub_mix_dens[E_idx]
            N_dens = sub_mix_dens[N_idx]
            S_dens = sub_mix_dens[S_idx]
            NS_mean_dens = np.mean((N_dens, S_dens))
            E_densities.append(E_dens)
            NS_densities.append(NS_mean_dens)
            # store the weighted sum densities
            mix_dens = mix_dens + sub_mix_dens
        # get the mean of the multiple mixture-dist densities
        mix_dens = mix_dens/mu.shape[0]
        # store the mix-dist densities
        plot_vals.extend(mix_dens)

    # convert to degrees, if necessary
    if to_deg:
        plot_vals_deg = [dir/2/pi*360 for dir in plot_vals]
        plot_vals_deg = [dir + 360 if dir <0 else dir for dir in plot_vals_deg]

    # plot it
    if plot_type == 'hist':
        if to_deg:
            plt.hist(plot_vals_deg, bins=200, color=col, alpha=0.6)
        else:
            plt.hist(plot_vals, bins=200, color=col, alpha=0.6)

    elif plot_type == 'circ':
        # NOTE: need to rotate the angles artificially in order to line up
        # 'true north' (derived from gnx as 0 degrees at compass N, increasing
        # clockwise) with what would be the default output of plotting
        # cos(dirs) as xs and sin(dirs) as ys ([0-->2pi] would start at compass
        # east and increase counterclockwise);
        # to do this, just use cos(dirs) as ys and sin(dirs) as xs instead!
        # convert dirs and densities to coords on the unit circle
        ys = np.cos(dirs) * plot_vals
        ref_ys = np.cos(dirs) * np.mean(plot_vals)
        xs = np.sin(dirs) * plot_vals
        ref_xs = np.sin(dirs) * np.mean(plot_vals)
        # make polygon patchcollection to display the dist
        patches = []
        poly = Polygon(np.array([*zip(xs, ys)]), True, color=col)
        patches.append(poly)
        pc = PatchCollection(patches, alpha=0.6, match_original=True)
        # add axes, for reference
        if ax is None:
            fig = plt.figure(figsize=(fig_width, fig_height))
            ax = fig.add_subplot(111)
            plt.plot([0,0], [-1,1], '-', color='gray')
            plt.plot([-1,1], [0,0], '-', color='gray')
            # plot the resulting density on top of mean-scaled reference circle
            plt.plot(ref_xs, ref_ys, ':', color='gray', alpha=0.5)
            ax.add_collection(pc)
            #plt.plot(xs, ys, color=col)
        else:
            # plot the resulting density on top of mean-scaled reference circle,
            # with an arrow indicating the direction of climate shift
            ax.plot(ref_xs, ref_ys, ':', color='gray', alpha=0.5)
            ax.plot([0,0], [-1,1], ':', color='gray', linewidth=0.5, alpha=0.5)
            ax.plot([-1,1], [0,0], ':', color='gray', linewidth=0.5, alpha=0.5)
            ax.add_collection(pc)
            #ax.plot(xs, ys, color=col)
        lim_val = axlim_graycirc_factor*np.max(np.abs(np.concatenate((ref_xs, ref_ys))))
        axlims = [-1*lim_val, lim_val]
        plt.xlim(axlims)
        plt.ylim(axlims)
        # get rid of the ticks and label cardinal directions instead
        ax.set_xticks([])
        ax.set_xticklabels([])
        ax.set_yticks([])
        ax.set_yticklabels([])
        xlims = ax.get_xlim()
        ylims = ax.get_ylim()

    # make sure the E-facing and N/S-facing density data are equal in length
    assert len(E_densities) == len(NS_densities)

    if return_axlims and plot_type == 'circ':
        return axlims, E_densities, NS_densities
    else:
        return E_densities, NS_densities


def make_vM_mix_dist_comparison_grid(df, genicities,
                                     it=None,
                                     col_labelsize=col_label_fontsize,
                                     row_labelsize=row_label_fontsize,
                                     plot_type='circ',
                                     return_E_NS_dens_df=True):
    """
    Visualize the null/non-null overlain von Mises mixture distributions
    across all 9 simulation scenarios

    NOTE: just defaults to only non-neutral loci for now
    """
    # dict to store all E-facing and mean N/S-facing densities and their covars
    E_NS_dens_dict = {'nullness': [],
                      'genicity': [],
                      'linkage': [],
                      'E_dens': [],
                      'NS_mean_dens': [],
                     }
    # list to store axlims
    axlims = []
    # make the gridded fig
    fig = plt.figure()
    # counter to keep track of plot number
    plt_num = 1
    linkages = ['independent', 'weak', 'strong']
    for linkage in linkages:
        for genicity in df.genicity.unique():
            # make the next axis instance
            ax = fig.add_subplot(3, 3, plt_num)

            # manage row and column labels
            ax.set_xlabel('')
            ax.set_ylabel('')
            ax.set_title('')
            col_labs = dict(zip([1,2,3],
                                [g for g in genicities]))
            if plt_num in [1, 2, 3]:
                ax.set_title(col_labs[plt_num], size=col_labelsize)
            if plt_num in [1, 4, 7]:
                ax.set_ylabel(linkage, fontsize=row_labelsize)

            # subset the df
            subdf = df[(df.linkage == linkage) &
                       (df.genicity == genicity)]
            assert len(subdf.linkage.unique()) == 1
            if len(ax.get_ylabel()) > 0:
                assert subdf.linkage.unique()[0] == ax.get_ylabel()
            # subset further, by nullness
            for nullness in ['null', 'non-null']:
                nullness_subdf = subdf[subdf.nullness == nullness]
                # subset for only a certain it, if requested
                if it is not None:
                    plotdf = nullness_subdf[nullness_subdf.it == it]
                # otherwise, plot mean densities (calculated within
                # viz_vM_mix_dist) across all its
                else:
                    plotdf = nullness_subdf
                curr_axlims, E_densities, NS_densities = viz_vM_mix_dist(
                                              df=plotdf,
                                              nullness=nullness,
                                              on_curr_ax=True,
                                              plot_type=plot_type)
                axlims.append(curr_axlims)

                # store the E-facing and mean N/S-facing density data
                for E_dens_val, NS_dens_val in zip(E_densities, NS_densities):
                    E_NS_dens_dict['nullness'].append(nullness)
                    E_NS_dens_dict['genicity'].append(genicity)
                    E_NS_dens_dict['linkage'].append(linkage)
                    E_NS_dens_dict['E_dens'].append(E_dens_val)
                    E_NS_dens_dict['NS_mean_dens'].append(NS_dens_val)


            # after overplotting null & non-null, increment to next plot number
            plt_num += 1

    # get max axlims & set for all axes, add direction labels, fix aspect ratios
    if plot_type == 'circ':
        abs_axlim_vals = np.array([lims[1] for lims in axlims])
        max_axlims_idx = np.where(abs_axlim_vals == np.max(abs_axlim_vals))[0][0]
        univ_axlims = axlims[max_axlims_idx]
        for ax in fig.get_axes():
            ax.set_xlim(univ_axlims)
            ax.set_ylim(univ_axlims)
            # dir label
            smidge = univ_axlims[1]*0.05
            edge = univ_axlims[1] - smidge
            arrow_xs = [0.9*edge, edge, 0.9*edge]
            arrow_ys = [0.1*univ_axlims[0],
                      np.mean(univ_axlims),
                      0.1*univ_axlims[1]]
            ax.plot(arrow_xs, arrow_ys, color='#a83644', alpha=0.7, linewidth=2)
            # 1/1 aspect ratio
            ax.set_aspect('equal')

    # adjust spacing of subplots
    plt.subplots_adjust(left=0.01,
                        bottom=0.04,
                        right=0.99,
                        top=0.93,
                        wspace=0.02,
                        hspace=None
                       )
    plt.show()
    if return_E_NS_dens_df:
        E_NS_dens_df = pd.DataFrame.from_dict(E_NS_dens_dict)
        linkage_numeric_dict = {'independent':0.5, 'weak':0.05, 'strong':0.005}
        E_NS_dens_df['linkage'] = [linkage_numeric_dict[l] for
                                                    l in E_NS_dens_df.linkage]
        return (fig, E_NS_dens_df)
    return (fig, None)


# run for both high and low redundancy
for redundancy in ['hi', 'lo']:

    print(('\n\n\n' +
           '~'*80 +
           'MAKING VON MISES VIZ FOR %s-REDUND SCENARIOS...\n\n' % redundancy))

    genicities = [4, 20, 100]
    if redundancy == 'hi':
        genicities = [g*2 for g in genicities]

    # load data
    df = pd.read_csv(os.path.join(analysis_dir,
                    'ch2_all_fitted_vM_params_%sREDUND.csv' % redundancy),
                     na_filter=False)

    # plot the full grid
    grid_fig, E_NS_dens_df = make_vM_mix_dist_comparison_grid(df,
                                                              genicities,
                                                              plot_type='circ',
                                                              it=None,
                                                           return_E_NS_dens_df=True)
    E_NS_dens_df.to_csv(os.path.join(analysis_dir,
                    'ch2_E_NS_gene_flow_densities_%sREDUND.csv' % redundancy),
                        index=False)

    grid_fig.savefig(os.path.join(analysis_dir,
                    'ch2_gene_flow_dir_analysis_%sREDUND.png' % redundancy),
                     dpi=400)
