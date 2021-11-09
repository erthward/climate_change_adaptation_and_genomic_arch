import numpy as np
from numpy import pi as pi
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import vonmises
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

# TODO:
    # now wrap the making of a single null/non-null comparison in a fn
    # then make a 3x3 plot matrix, loop over the scenarios, and plot comparison
    # example of each (for non-neutral loci only)

    # then get R script running across whole dir output files (~3.4 GB) on
    # savio

    # then loop that over all dir files as a job

    # then combine all results into a master df

    # then write linear model to look at multivariate response of dists to the
    # scenarios

    # then consider if/how to depict uncertainty envelopes on the dists

    # then produce a single, 3x3 plot matrix depicting the null and non-null
    # characteristic dists across all scenarios!

# load data
dat = pd.read_csv('./TEST_output.csv', na_filter=False)


def vis_vM_mix_dist(row=None, mu=None, kappa=None, alpha=None, nullness=None,
                    plot_type='circ', to_deg=True,
                    on_curr_ax=False, n_rand=10000):
    """
    Visualize a von Mises mixture distribution

    NOTE: n_rand is only going to be approximate (but very close)
    to true total number of draws used to visualize the dist
    """
    # need mu, kappa, alpha, and nullness if row is None
    if row is None:
        assert (mu is not None and
                kappa is not None and
                alpha is not None and
                nullness is not None), ("MU, KAPPA, ALPHA, AND NULLNESS "
                                        "ARGS MISSING")

    # otherwise, grab them from the row
    else:
        mu = [*row[-12:-8]]
        kappa = [*row[-8:-4]]
        alpha = [*row[-4:]]
        nullness = row['nullness']

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
        # determine how many random vals to draw from each dist
        samp_sizes = [int(n) for n in np.floor(np.array(alpha)*n_rand)]
        # draw vals
        for mix_n, samp_size in enumerate(samp_sizes):
            plot_vals.extend(vonmises.rvs(kappa[mix_n],
                                     loc=mu[mix_n],
                                     scale=1,
                                     size=samp_size))
            # correct to the interval [0, 2pi]
            plot_vals = [val+(2*pi) if val <0 else val for val in plot_vals]
    elif plot_type == 'circ':
        # get ~3000-ish evenly spaced densities
        dirs = np.arange(0, 2*pi, 0.002)
        # numpy array to contain the output mixture-dist densities
        mix_dens = np.zeros(dirs.shape)
        weighted_densities = []
        for mix_n in range(len(mu)):
            # calculate densities
            dens = vonmises.pdf(dirs,
                                kappa=kappa[mix_n],
                                loc=mu[mix_n],
                                scale=1)
            # weight densities by this distribution's alpha
            weight_dens = alpha[mix_n] * dens
            # store the vals
            mix_dens = mix_dens + weight_dens
        # store the mix-dist densities
        plot_vals.extend(mix_dens)

    # convert to degrees, if necessary
    if to_deg:
        plot_vals_deg = [dir/2/pi*360 for dir in plot_vals]
        plot_vals_deg = [dir + 360 if dir <0 else dir for dir in plot_vals_deg]

    # plot it
    if plot_type == 'hist':
        if to_deg:
            plt.hist(plot_vals_deg, bins=200, color=col)
        else:
            plt.hist(plot_vals, bins=200, color=col)

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
        minmax_xy = [-1*lim_val, lim_val]
        plt.xlim(minmax_xy)
        plt.ylim(minmax_xy)
        # get rid of the ticks and label cardinal directions instead
        ax.set_xticks([])
        ax.set_xticklabels([])
        ax.set_yticks([])
        ax.set_yticklabels([])
        xlims = ax.get_xlim()
        ylims = ax.get_ylim()
        if not on_curr_ax:
            plt.text(0+0.03*xlims[0], 1.05*ylims[1], 'N', size=18)
            plt.text(1.05*xlims[1], 0+0.05*ylims[0], 'E', size=18)
            plt.text(0+0.03*xlims[0], 1.13*ylims[0], 'S', size=18)
            plt.text(1.15*xlims[0], 0+0.05*ylims[0], 'W', size=18)

    return ax


# plot an example!
row = dat.iloc[50,:] # null
vis_vM_mix_dist(row)
row = dat.iloc[49,:] # non-null
vis_vM_mix_dist(row, on_curr_ax=True)
plt.show()
