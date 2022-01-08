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


# plot params
col_label_fontsize = 15
row_label_fontsize = 13
cbar_fontsize = 9
fig_width = 13.8
fig_height = 4.6
dpi = 400
n_ticklabels = 5

# data directory
if os.getcwd().split('/')[1] == 'home':
    analysis_dir = ('/home/deth/Desktop/CAL/research/projects/sim/'
               'ch2/climate_change_adaptation_and_genomic_arch/analysis')
else:
    with open(('/global/scratch/users/drewhart/ch2/climate_change_adaptation_'
               'and_genomic_arch/analysis/analysisdir.txt'), 'f') as f:
        analysis_dir = f.read().strip()
    #analysis_dir = '/global/scratch/users/drewhart/ch2/output/analysis'

# TODO:

    # write linear model to look at multivariate response of dists to the
    # scenarios?

    # consider if/how to depict uncertainty envelopes on the dists?


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

    # make sure the E-facing and N/S-facing density data are equal in length
    assert len(E_densities) == len(NS_densities)

    if return_axlims and plot_type == 'circ':
        return axlims, E_densities, NS_densities
    else:
        return E_densities, NS_densities


def make_vM_mix_dist_comparison_grid(df, it=None,
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
    for linkage in df.linkage.unique():
        for genicity in df.genicity.unique():
            # make the next axis instance
            ax = fig.add_subplot(3, 3, plt_num)

            # manage row and column labels
            ax.set_xlabel('')
            ax.set_ylabel('')
            ax.set_title('')
            col_labs = dict(zip([1,2,3],
                                #['genicity: %i' % g for g in [4, 20, 100]]))
                                [g for g in [4, 20, 100]]))
            row_labs = dict(zip([1, 4, 7],
                                 #['linkage: %s' % l for l in ['independent',
                                 [l for l in ['independent',
                                                            'weak', 'strong']]))
            if plt_num in [1, 2, 3]:
                ax.set_title(col_labs[plt_num], size=col_labelsize)
            if plt_num in [1, 4, 7]:
                ax.set_ylabel(row_labs[plt_num], fontsize=row_labelsize)

            # subset the df
            subdf = df[(df.linkage == linkage) &
                       (df.genicity == genicity)]
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


# load data
df = pd.read_csv(os.path.join(analysis_dir, 'ch2_all_fitted_vM_params.csv'),
                 na_filter=False)

# plot the full grid
grid_fig, E_NS_dens_df = make_vM_mix_dist_comparison_grid(df,
                                                          plot_type='circ',
                                                          it=None,
                                                       return_E_NS_dens_df=True)
E_NS_dens_df.to_csv(os.path.join(analysis_dir, 'ch2_E_NS_gene_flow_densities.csv'),
                    index=False)

grid_fig.savefig(os.path.join(analysis_dir, 'ch2_gene_flow_dir_analysis.png'),
                 dpi=400)


# run null and non-null statistical models, and gather results into a small df
results = {'nullness': [],
           'response_var': [],
           'covar': [],
           'coeff': [],
           'p': [],
           'R2': [],
           }
    # NOTES:
        # response vars aren't exactly normal (a bit skewed), but not awful
        # either, so just using a basic linear model for now
all_tuk_data = []
for nullness in ['null', 'non-null']:
    print('\n%s MODELS\n:' % nullness.upper())
    for response_var in ['E_dens', 'NS_mean_dens']:
        print('\n\tresponse var: %s\n' % response_var)
        subdf = E_NS_dens_df[E_NS_dens_df['nullness'] == nullness]

        # code for two-way ANOVA with post-hoc Tukey HSD test stolen from:
        # https://towardsdatascience.com/anova-tukey-test-in-python-b3082b6e6bda

        # R-style formula specification using a string
        # NOTE: C() operator coerces a variable to categorical
        formula = ('%s ~ C(linkage) + C(genicity) + '
                   'C(linkage):C(genicity)') % response_var
        lm = ols(formula, subdf).fit()
        table = sm.stats.anova_lm(lm, typ=2)
        print(table)
        # perform multiple pairwise comparison (Tukey HSD)
        subdf['combination'] = (subdf.linkage.astype(str) + " / " +
                                subdf.genicity.astype(str))
        tuk = pairwise_tukeyhsd(endog=subdf[response_var],
                                groups=subdf['combination'], alpha=0.05)
        # coerce the tukeyhsd table to a DataFrame, then filter and sort
        tukey_data = pd.DataFrame(data=tuk._results_table.data[1:],
                                  columns = tuk._results_table.data[0])
        group1_comp =tukey_data.loc[tukey_data.reject ==
                                    True].groupby('group1').reject.count()
        group2_comp = tukey_data.loc[tukey_data.reject ==
                                     True].groupby('group2').reject.count()
        tukey_data_summ = pd.concat([group1_comp, group2_comp], axis=1)
        tukey_data_summ = tukey_data_summ.fillna(0)
        tukey_data_summ.columns = ['reject1', 'reject2']
        tukey_data_summ['total_sum'] = (tukey_data_summ.reject1 +
                                        tukey_data_summ.reject2)
        print('\n\n\tTukey HSD pairwise comparison test results:\n\n')
        print(tukey_data_summ.sort_values('total_sum',ascending=False))
        # save these Tukey results to our list
        tukey_data['nullness'] = nullness
        tukey_data['dir'] = response_var
        all_tuk_data.append(tukey_data)


        # build identical model, but with numeric instead of categorical vars
        y = subdf[response_var]
        X = sm.add_constant(subdf[['linkage', 'genicity']])
        mod = sm.OLS(y, X).fit()
        coeffs = mod.params
        print(mod.summary())
        print('\n'*4 + '+'*80 + '\n'*4)
        p_vals = mod.pvalues
        R2 = mod.rsquared
        for coeff in coeffs.index:
            results['nullness'].append(nullness)
            results['response_var'].append(response_var)
            results['covar'].append(coeff)
            results['coeff'].append(coeffs[coeff])
            results['p'].append(coeffs[coeff])
            results['R2'].append(R2)

results_df = pd.DataFrame.from_dict(results)
results_df.to_csv(os.path.join(analysis_dir, 'ch2_gene_flow_dir_analysis_stats.csv'),
                  index=False)
tuk_df = pd.concat(all_tuk_data)
tuk_df.to_csv(os.path.join(analysis_dir, 'ch2_gene_flow_dir_tukey_results.csv'),
                  index=False)
