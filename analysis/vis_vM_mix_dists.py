import numpy as np
from numpy import pi as pi
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import vonmises
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

# TODO:
    # get R script running across whole dir output files (~3.4 GB) on
    # savio

    # then loop that over all dir files as a job

    # then combine all results into a master df

    # then write linear model to look at multivariate response of dists to the
    # scenarios

    # then consider if/how to depict uncertainty envelopes on the dists

    # then produce a single, 3x3 plot matrix depicting the null and non-null
    # characteristic dists across all scenarios


def vis_vM_mix_dist(row=None, df=None, mu=None, kappa=None, alpha=None, nullness=None,
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
            fig = plt.figure()
            ax = fig.add_subplot(111)
            plt.plot([0,0], [-1,1], '-', color='gray')
            plt.plot([-1,1], [0,0], '-', color='gray')
            # plot the resulting density on top of mean-scaled reference circle
            plt.plot(ref_xs, ref_ys, ':', color='gray', alpha=0.5)
            ax.add_collection(pc)
            #plt.plot(xs, ys, color=col)
        else:
            # plot the resulting density on top of mean-scaled reference circle
            ax.plot(ref_xs, ref_ys, ':', color='gray', alpha=0.5)
            ax.add_collection(pc)
            #ax.plot(xs, ys, color=col)
        lim_val = 1.1*np.max(np.abs(np.concatenate((xs, ys, ref_xs, ref_ys))))
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
        if not on_curr_ax or plot_dirlabels:
            plt.text(0+0.03*xlims[0], 1.05*ylims[1], 'N', size=18, color='gray')
            plt.text(1.05*xlims[1], 0+0.05*ylims[0], 'E', size=18, color='gray')
            plt.text(0+0.03*xlims[0], 1.13*ylims[0], 'S', size=18, color='gray')
            plt.text(1.15*xlims[0], 0+0.05*ylims[0], 'W', size=18, color='gray')

    if return_axlims and plot_type == 'circ':
        return axlims
    else:
        return


def make_vM_mix_dist_comparison_grid(df, neutrality='nonneut', it=None,
                                     labelsize=16, plot_type='circ'):
    """
    Visualize the null/non-null overlain von Mises mixture distributions
    across all 9 simulation scenarios

    NOTE: just defaults to only non-neutral loci for now
    """
    # list to store axlims
    axlims = []
    # make the gridded fig
    fig = plt.figure()
    # counter to keep track of plot number
    plt_num = 1
    for linkage in df.linkage.unique():
        for genicity in df.genicity.unique():
            # make the next axis instance
            ax = fig.add_subplot(3, 3, plt_num)

            # manage row and column labels
            ax.set_xlabel('')
            ax.set_ylabel('')
            ax.set_title('')
            col_labs = dict(zip([1,2,3],
                                ['genicity: %i' % g for g in [4, 20, 100]]))
            row_labs = dict(zip([1, 4, 7],
                                 ['linkage: %s' % l for l in ['independent',
                                                             'weak',
                                                             'strong']]))
            if plt_num in [1, 2, 3]:
                ax.set_title(col_labs[plt_num], size=labelsize)
            if plt_num in [1, 4, 7]:
                ax.set_ylabel(row_labs[plt_num], fontsize=labelsize)

            # subset the df
            subdf = df[(df.linkage == linkage) &
                       (df.genicity == genicity) &
                       (df.neutrality == neutrality)]
            # subset further, by nullness
            for nullness in ['null', 'non-null']:
                nullness_subdf = subdf[subdf.nullness == nullness]
                # subset for only a certain it, if requested
                if it is not None:
                    plotdf = nullness_subdf[nullness_subdf.it == it]
                # otherwise, plot mean densities (calculated within
                # vis_vM_mix_dist) across all its
                else:
                    plotdf = nullness_subdf
                curr_axlims = vis_vM_mix_dist(df=plotdf,
                                              nullness=nullness,
                                              on_curr_ax=True,
                                              plot_type=plot_type)
                axlims.append(curr_axlims)


            # after overplotting null & non-null, increment to next plot number
            plt_num += 1

    # get max axlims and set for all axes, and add direction labels, and fix
    # aspect ratios
    if plot_type == 'circ':
        abs_axlim_vals = np.array([lims[1] for lims in axlims])
        max_axlims_idx = np.where(abs_axlim_vals == np.max(abs_axlim_vals))[0][0]
        univ_axlims = axlims[max_axlims_idx]
        for ax in fig.get_axes():
            ax.set_xlim(univ_axlims)
            ax.set_ylim(univ_axlims)
            # dir labels 
            # NOTE: for now, just adjusting labels a 'smidge' to approx. center
            #       them with the ref circle, but this will easily break
            smidge = univ_axlims[1]*0.1
            ax.text(0-smidge, univ_axlims[1]*0.5-smidge, 'N', size=14, color='gray')
            ax.text(univ_axlims[1]*0.5-smidge, 0-smidge, 'E', size=14, color='gray')
            ax.text(0-smidge, univ_axlims[0]*0.5-smidge, 'S', size=14, color='gray')
            ax.text(univ_axlims[0]*0.5-smidge, 0-smidge, 'W', size=14, color='gray')
            # 1/1 aspect ratio
            ax.set_aspect('equal')

    # show and return fig
    plt.show()
    return fig


# load data
df = pd.read_csv('./ch2_fitted_vM_params.csv', na_filter=False)

# plot the full grid
grid_fig = make_vM_mix_dist_comparison_grid(df, plot_type='circ', it=None)