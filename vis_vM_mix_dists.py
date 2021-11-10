import numpy as np
from numpy import pi as pi
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import vonmises
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

# TODO:
    # debug the grid plots:
        # just bc little data???
        # something fundamentally screwy???
        # do I need to keep my sort_multicol_groups fn?

    # then get R script running across whole dir output files (~3.4 GB) on
    # savio

    # then loop that over all dir files as a job

    # then combine all results into a master df

    # then write linear model to look at multivariate response of dists to the
    # scenarios

    # then consider if/how to depict uncertainty envelopes on the dists

    # then produce a single, 3x3 plot matrix depicting the null and non-null
    # characteristic dists across all scenarios!


def sort_multicol_groups(arr, sortcols, n_groups):
    """
    sort rows of an array in column-groups across, based on the order that
    sorts the numbers in the sortcols

    NOTE: n_groups indicates the TOTAL number of col groups, including the
          sortcols' group!

    NOTE: THE SECOND LINE OF CODE MAKES THE HARD ASSUMPTION THAT THE SORTCOLS
          ARE THE LEFTMOST COL GROUP IN THE ARRAY!
    """
    # get list of the indices needed to sort the sortcols' vals in each row
    sortcols_idxs = np.stack([np.argsort(arr[i,
                                    sortcols]) for i in range(arr.shape[0])])
    # get the list of the indices needed to sort all cols' vals in each row,
    # based on the sortcols
    sort_idxs = np.hstack([sortcols_idxs] + [
                        sortcols_idxs+len(sortcols) for _ in range(n_groups-1)])

    # sort the array using the sort_idxs
    sorted = np.array([arr[np.repeat(i, len(sortcols)*n_groups),
                           sort_idxs[i, :]] for i in range(arr.shape[0])])

    return sorted


def vis_vM_mix_dist(row=None, mu=None, kappa=None, alpha=None, nullness=None,
                    plot_type='circ', to_deg=True, plot_dirlabels=False,
                    on_curr_ax=False, return_axlims=True, n_rand=10000):
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
        if nullness is None:
            nullness = row['nullness']
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

    if return_axlims:
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
                # otherwise, get means across all its
                else:
                    # TODO!
                    plotdf = pd.DataFrame(nullness_subdf.mean(axis=0)).T
                # plot this df's data
                assert len(plotdf) == 1, (">1 row found for this plot!\n%s" %(
                                                                   str(plotdf)))
                print(linkage, genicity, nullness)
                print(nullness_subdf)
                print(plotdf)
                print('\n' + '='*80 + '\n')
                curr_axlims = vis_vM_mix_dist(plotdf.iloc[0,:],
                                              nullness=nullness,
                                              on_curr_ax=True,
                                              plot_type=plot_type)
                axlims.append(curr_axlims)


            # after overplotting null & non-null, increment to next plot number
            plt_num += 1

    # get max axlims and set for all axes, and add direction labels
    abs_axlim_vals = np.array([lims[1] for lims in axlims])
    max_axlims_idx = np.where(abs_axlim_vals == np.max(abs_axlim_vals))[0][0]
    univ_axlims = axlims[max_axlims_idx]
    for ax in fig.get_axes():
        ax.set_xlim(univ_axlims)
        ax.set_ylim(univ_axlims)
        # dir labels 
        ax.text(0, 1, 'N', size=18, color='gray')
        ax.text(1, 0, 'E', size=18, color='gray')
        ax.text(0, -1, 'S', size=18, color='gray')
        ax.text(-1, 0, 'W', size=18, color='gray')

    # show and return fig
    plt.show()
    return fig


# load data
df = pd.read_csv('./TEST_output.csv', na_filter=False)

# resort the params' columns so that vM mix dists' params are listed across
# the rows in order of increasing mu/loc param, so that 
# the params' values can then be averaged down the cols
# NOTE: sortcols are [0,1,2,3] because they're expressed with reference to the
#       input array, not the DataFrame
#df.iloc[:,5:] = sort_multicol_groups(df.iloc[:, 5:].values, [0,1,2,3], 3)

# plot the full grid
grid_fig = make_vM_mix_dist_comparison_grid(df, plot_type='circ', it=None)
#grid_fig = make_vM_mix_dist_comparison_grid(df, plot_type='hist', it=1)

# plot an example!
#row = df.iloc[50,:] # null
#vis_vM_mix_dist(row)
#row = df.iloc[49,:] # non-null

