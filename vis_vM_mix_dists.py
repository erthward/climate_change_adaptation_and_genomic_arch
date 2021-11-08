import numpy as np
from numpy import pi as pi
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import vonmises

# TODO:
    # figure out why doesn't work on 'row mode'

    # almost certain that somehow N-S and E-W axes are swapped in my results;
    # figure out why

    # consider making polygon patches and plotting those instead...


# load data
dat = pd.read_csv('./TEST_output.csv')


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
        print(mu)
        print(kappa)
        print(alpha)
        print(nullness)

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
    elif plot_type == 'circ':
        # get ~3000-ish evenly spaced densities
        dirs = np.arange(-pi, pi, 0.002)
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
            plt.hist(plot_vals_deg, bins=200, col=col, ax=ax)
        else:
            plt.hist(plot_vals, bins=200, col=col, ax=ax)

    elif plot_type == 'circ':
        # convert dirs and densities to coords on the unit circle
        xs = np.cos(dirs) * plot_vals
        ref_xs = np.cos(dirs) * np.mean(plot_vals)
        ys = np.sin(dirs) * plot_vals
        ref_ys = np.sin(dirs) * np.mean(plot_vals)
        # add axes, for reference
        if ax is None:
            plt.plot([0,0], [-1,1], '-k')
            plt.plot([-1,1], [0,0], '-k')
            # plot the resulting density on top of mean-scaled reference circle
            plt.plot(ref_xs, ref_ys, ':', color='gray', alpha=0.5)
            plt.plot(xs, ys, color=col)
        else:
            # plot the resulting density on top of mean-scaled reference circle
            ax.plot(ref_xs, ref_ys, ':', color='gray', alpha=0.5)
            ax.plot(xs, ys, color=col)
        min_xy = 1.05*min([*xs]+[*ref_xs]+[*xs]+[*ref_xs])
        max_xy = 1.05*max([*xs]+[*ref_xs]+[*xs]+[*ref_xs])
        minmax_xy = [-1*max(np.abs(np.array([min_xy, max_xy]))),
                     max(np.abs(np.array([min_xy, max_xy])))]
        plt.xlim(minmax_xy)
        plt.ylim(minmax_xy)

    return ax



