#!/usr/bin/python
#prelim.py

# flake8: noqa

import geonomics as gnx
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import copy
import tskit
import bisect
import time


#################
# set main params
#################

# set time when environmental change begins
change_T = 10 #500
# set total time
T = 20 #800
# set number of iterations for each sim
n_its = 1
# set the different numbers of loci to use
n_loci = [5, 20, 100]
# reset the K_factor (if desired)
K_factor=0.25

# flag indicating whether or not to plot the lineage map
map_gene_flow = True


# create ParamsDict objects
filepath=('/home/drew/Desktop/stuff/berk/research/projects/sim/'
          'ch2_adapt_clim_chng_genarch/prelim_params.py')
unlinked_params = gnx.read_parameters_file(filepath=filepath)

# tweak the carrying capacity, if requested
if K_factor is not None:
    unlinked_params['comm']['species']['spp_0']['init']['K_factor'] = K_factor

linked_params = copy.deepcopy(unlinked_params)

# set the linkage values for the linked scenarios' params
linked_params['comm']['species']['spp_0']['gen_arch']['r_distr_alpha'] = 5000
linked_params['comm']['species']['spp_0']['gen_arch']['r_distr_beta'] = 1e7





##################
# define functions
##################

def save_data(l, mod, data_dict):
    neut_loci = np.random.choice(mod.comm[0].gen_arch.neut_loci, l,
                                 replace=False)
    nonneut_loci = mod.comm[0].gen_arch.traits[0].loci
    neut_stats = mod.comm[0]._calc_lineage_stats(
        stats=['dir', 'dist', 'speed'], max_time_ago=300, loci=neut_loci)
    nonneut_stats = mod.comm[0]._calc_lineage_stats(
        stats=['dir', 'dist', 'speed'], max_time_ago=300, loci=nonneut_loci)
    for stat, dataset in nonneut_stats.items():
        # add the non-neutral stats to their lists
        for data in dataset.values():
            data_dict[stat][0].extend(data)
        # add the neutral data too
        for data in neut_stats[stat].values():
            data_dict[stat][1].extend(data)


# function to run both the linked and unlinked sim for a given number of loci
def run_sims_l_loci(l, u_params, l_params, n_its=1,
                    with_env_change=True):
    #create deep copies of the params
    unlinked_params = copy.deepcopy(u_params)
    linked_params = copy.deepcopy(l_params)

    # set the model names (so that they save data correctly in separate
    # directories)
    unlinked_params.model['name'] = 'unlinked_%iloc' % l
    linked_params.model['name'] = 'linked_%iloc' % l

    # get the start-time of the environmental change event (for plotting later)
    change_start_t = unlinked_params.landscape.layers.shift.change[0].start_t

    # edit the parameters so that no environmental change occurs,
    # and there is no gradient, if with_env_change is False
    if not with_env_change:
        unlinked_params['landscape']['layers']['shift'].pop('change');
        linked_params['landscape']['layers']['shift'].pop('change');
        unlinked_params['landscape']['layers']['shift']['init'][
            'defined']['rast'] = np.ones((50, 50))*0.5
        unlinked_params['landscape']['layers']['stable']['init'][
            'defined']['rast'] = np.ones((50, 50))*0.5
        linked_params['landscape']['layers']['shift']['init'][
            'defined']['rast'] = np.ones((50, 50))*0.5
        linked_params['landscape']['layers']['stable']['init'][
            'defined']['rast'] = np.ones((50, 50))*0.5


    # set the numbers of loci for both scenarios, the
    # effect-size distribution to match the number of loci,
    # and the total genome lengths
    for p in [unlinked_params, linked_params]:
        p['comm']['species']['spp_0']['gen_arch']['L'] = 3 * l
        for trt in p['comm']['species']['spp_0']['gen_arch']['traits'].values():
            trt['n_loci'] = l
            trt['alpha_distr_mu'] = 1/l

    # set the plotting colors
    if with_env_change:
        line_color = '#f59f2f' # orange
    else:
        line_color = '#f2d43a' # yellow


    # create the main and lineage figures
    fig = plt.figure('main %i' % l)
    fig.suptitle('%i loci per trait' % l)

    if map_gene_flow:
        gf_fig = plt.figure('gene flow %i' % l)
        gf_fig.suptitle(('%i loci per trait\nneutral loci in light '
                         'gray, nonneutral in dark gray') % l)
    else:
        gf_fig = None


    # create lineage-stat data structures
    # for each stat, the first list will hold non-neutral locus data, the
    # second will hold paired neutral locus data
    unlinked_data_dict = {'dir': [[], []],
                        'dist': [[], []],
                        'speed': [[], []],
                       }

    linked_data_dict = {'dir': [[], []],
                        'dist': [[], []],
                        'speed': [[], []],
                       }

    #################
    #### UNLINKED SIM
    #################

    for it in range(n_its):

        # create the model
        mod = gnx.make_model(gnx.make_params_dict(unlinked_params))

        # check that the number of trait 0 loci is correct
        assert len(mod.comm[0].gen_arch.traits[0].loci) == l, (
                'EXPECETED %i LOCI BUT GOT %i!') % (
                l, len(mod.comm[0].gen_arch.traits[0].loci))

        # create data struct for the fitness values
        unlinked_fit = []

        # burn the model in
        mod.walk(10000, mode='burn')


        # run the model
        for t in range(change_T):

            # keep printing the number of loci,
            # to help me track things while it's running
            n_loci = len(mod.comm[0].gen_arch.traits[0].loci)
            print('\n[unlinked, %i loci, it %i]' % (n_loci, it))

            mod.walk(1, mode='main')

            # save the fitness value
            unlinked_fit.append(np.mean(mod.comm[0]._get_fit()))

        # add the pre-change population to the fig
        if it == 0:
            plt.figure('main %i' % l)
            ax1 = fig.add_subplot(251)
            ax1.set_title('Pre-change, Unlinked, trait 0', size=16)
            mod.plot_phenotype(0, 0, alpha=0.5)

        # run the model
        for t in range(change_T, T):

            # keep printing the number of loci,
            # to help me track things while it's running
            n_loci = len(mod.comm[0].gen_arch.traits[0].loci)
            print('\n[unlinked, %i loci, it %i]' % (n_loci, it))

            mod.walk(1, mode='main')

            # save the fitness value
            unlinked_fit.append(np.mean(mod.comm[0]._get_fit()))

        # save the data
        save_data(l, mod, unlinked_data_dict)

        if it == 0:

            plt.figure('main %i' % l)
            # make the fitness and map figure
            ax2 = fig.add_subplot(252)
            ax2.set_title('Unlinked, trait 0', size=16)
            mod.plot_phenotype(0, 0, alpha=0.5)
            ax3 = fig.add_subplot(253)
            ax3.set_title('Unlinked, trait 1', size=16)
            mod.plot_phenotype(0, 1, alpha=0.5)
            ax4 = fig.add_subplot(254)
            ax4.set_title('Unlinked, pop size over time', size=16)
            mod.plot_pop_growth(0, actual=False)
            ax4.plot(range(len(mod.comm[0].Nt)), mod.comm[0].Nt,
                     color=line_color)
            ax4.plot([change_start_t] * 2, [0, max(mod.comm[0].Nt)], ':k',
                     alpha=0.5)
            ax5 = fig.add_subplot(255)
            ax5.set_title('Unlinked, fitness over time', size=16)
            ax5.plot(range(len(unlinked_fit)), unlinked_fit, color=line_color)
            ax5.plot([change_start_t] * 2, [0, 1], ':k', alpha=0.5)
            ax5.set_xlabel('t')
            ax5.set_ylabel('fitness')

            # make the gene-flow figure
            plt.figure('gene flow %i' % l)
            gf_ax1 = gf_fig.add_subplot(121)
            gf_ax1.set_title('Unlinked, trait 0', size=16)
            # define number of neutral and non-neutral loci to randomly draw
            n_locs=5
            nonneut_loci_to_plot = np.random.choice(
                        mod.comm[0].gen_arch.traits[0].loci, n_locs,
                        replace=False)

            neut_loci_to_plot = np.random.choice(
                        mod.comm[0].gen_arch.neut_loci, n_locs,
                        replace=False)
            nonneut_individs = np.random.choice([*mod.comm[0]], 10*n_locs,
                                                replace=False)
            neut_individs = np.random.choice([*mod.comm[0]], 10*n_locs,
                                             replace=False)
            for n, loc in enumerate(nonneut_loci_to_plot):
                mod.comm[0]._plot_gene_flow(loc, 'vector', mod.land,
                        individs=nonneut_individs[n*n_locs:n*n_locs+n_locs],
                        color='#616161', phenotype=0, size=35)
            for n, loc in enumerate(neut_loci_to_plot):
                mod.comm[0]._plot_gene_flow(loc, 'vector', mod.land,
                        individs=nonneut_individs[n*n_locs:n*n_locs+n_locs],
                        color='#f2f2f2', phenotype=0, size=35)

    ###############
    #### LINKED SIM
    ###############

        # make the model
        mod = gnx.make_model(gnx.make_params_dict(linked_params))

        # check that the number of trait 0 loci is correct
        assert len(mod.comm[0].gen_arch.traits[0].loci) == l, (
                    'EXPECETED %i LOCI BUT GOT %i!') % (
                    l, len(mod.comm[0].gen_arch.traits[0].loci))

        # create data struct for the linkage values
        linked_fit = []

        # burn the model in 
        mod.walk(10000, mode='burn')

        # run the model
        for t in range(change_T):

            # print out the number of trait 0 loci, to help me track things
            # while it's running
            n_loci = len(mod.comm[0].gen_arch.traits[0].loci)
            print('\n[linked, %i loci, it %i]' % (n_loci, it))

            mod.walk(1, mode='main')

            # store the fitness value
            linked_fit.append(np.mean(mod.comm[0]._get_fit()))

        # add the pre-change population to the fig
        if it == 0:
            plt.figure('main %i' % l)
            ax1 = fig.add_subplot(256)
            ax1.set_title('Pre-change, Linked, trait 0', size=16)
            mod.plot_phenotype(0, 0, alpha=0.5)

        # run the model
        for t in range(change_T, T):

            # print out the number of trait 0 loci, to help me track things
            # while it's running
            n_loci = len(mod.comm[0].gen_arch.traits[0].loci)
            print('\n[linked, %i loci, it %i]' % (n_loci, it))

            mod.walk(1, mode='main')

            # store the fitness value
            linked_fit.append(np.mean(mod.comm[0]._get_fit()))


        # save the data
        save_data(l, mod, linked_data_dict)

        if it == 0:
            plt.figure('main %i' % l)
            print('NOW MAKING PLOTS')
            # make the fitness and map figure
            ax7 = fig.add_subplot(257)
            ax7.set_title('Linked, trait 0', size=16)
            mod.plot_phenotype(0, 0, alpha=0.5)
            ax8 = fig.add_subplot(258)
            ax8.set_title('Linked, trait 1', size=16)
            mod.plot_phenotype(0, 1, alpha=0.5)
            ax9 = fig.add_subplot(259)
            ax9.set_title('Linked, pop size over time', size=16)
            mod.plot_pop_growth(0, actual=False)
            ax9.plot(range(len(mod.comm[0].Nt)), mod.comm[0].Nt,
                     color=line_color)
            ax9.plot([change_start_t] * 2, [0, max(mod.comm[0].Nt)], ':k',
                     alpha=0.5)
            ax10 = fig.add_subplot(2, 5, 10)
            ax10.set_title('Linked, fitness over time', size=16)
            ax10.plot(range(len(linked_fit)), linked_fit, color=line_color)
            ax10.plot([change_start_t] * 2, [0, 1], ':k', alpha=0.5)
            ax10.set_xlabel('t')
            ax10.set_ylabel('fitness')

            # make the gene-flow figure
            plt.figure('gene flow %i' % l)
            gf_ax2 = gf_fig.add_subplot(122)
            gf_ax2.set_title('Linked, trait 0', size=16)
            # define number of neutral and non-neutral loci to randomly draw
            n_locs=5
            nonneut_loci_to_plot = np.random.choice(
                        mod.comm[0].gen_arch.traits[0].loci, n_locs,
                        replace=False)
            neut_loci_to_plot = np.random.choice(
                        mod.comm[0].gen_arch.neut_loci, n_locs,
                        replace=False)
            nonneut_individs = np.random.choice([*mod.comm[0]], 10*n_locs,
                                                replace=False)
            neut_individs = np.random.choice([*mod.comm[0]], 10*n_locs,
                                             replace=False)
            for n, loc in enumerate(nonneut_loci_to_plot):
                mod.comm[0]._plot_gene_flow(loc, 'vector', mod.land,
                        individs=nonneut_individs[n*n_locs:n*n_locs+n_locs],
                        color='#616161', phenotype=0, size=35)
            for n, loc in enumerate(neut_loci_to_plot):
                mod.comm[0]._plot_gene_flow(loc, 'vector', mod.land,
                        individs=nonneut_individs[n*n_locs:n*n_locs+n_locs],
                        color='#f2f2f2', phenotype=0, size=35)


    return(fig, gf_fig, unlinked_data_dict, linked_data_dict)



###################################
# create containers for the figures
###################################

fig_dict = {}
gene_flow_fig_dict = {}
null_fig_dict = {}
null_gene_flow_fig_dict = {}

# create containers for the spatial-pedigree data that will be calculated
unlinked_data_dicts = {}
linked_data_dicts = {}
null_unlinked_data_dicts = {}
null_linked_data_dicts = {}



##############
# run the sims
##############

for l in n_loci:
    # run a linked/unlinked pair of scenarios for each number of loci
    print('#'*80 + '\nNOW RUNNING FOR %i LOCI\n' % l)
    main_fig, gene_flow_fig, unlinked_data, linked_data = run_sims_l_loci(l,
                                                        unlinked_params,
                                                        linked_params,
                                                        n_its=n_its)
    fig_dict[l] = main_fig
    gene_flow_fig_dict[l] = gene_flow_fig

    unlinked_data_dicts[l] = unlinked_data
    linked_data_dicts[l] = linked_data

    # then run a linked/unlinked pair without environmental change,
    # for comparison
    main_fig, gene_flow_fig, unlinked_data, linked_data = run_sims_l_loci(l,
                                                        unlinked_params,
                                                        linked_params,
                                                        n_its=n_its,
                                                        with_env_change=False)
    null_fig_dict[l] = main_fig
    null_gene_flow_fig_dict[l] = gene_flow_fig

    null_unlinked_data_dicts[l] = unlinked_data
    null_linked_data_dicts[l] = linked_data



##############################
# make gene-flow analysis figs
##############################

data_dict_dict = {'linked': linked_data_dicts,
                  'unlinked': unlinked_data_dicts
                 }

null_data_dict_dict = {'linked': null_linked_data_dicts,
                       'unlinked': null_unlinked_data_dicts
                      }

metric_text_dict = {'dir': 'direction\n(degrees, clockwise from N)',
                    'dist': 'distance\n(cell widths)',
                    'speed': 'gene flow rate\n(in cell widths per time step)'
                   }

dir_tick_locs = np.linspace(0, 360, 9)
dir_tick_labs = ['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW', 'N']

gf_hist_fig = plt.figure('hists')
ax_num = 1
for metric in ['dir', 'dist', 'speed']:
    for l in n_loci:
        for linkage, data_dict in data_dict_dict.items():
            # set the 4 histogram colors
            hist_cols = ['#a83252', '#99547d', '#25f5a2', '#2d7d52']
            ax = gf_hist_fig.add_subplot(3, 6, ax_num)
            ax.set_title('%i loci, %s' % (l, linkage))
            # plot a histogram of the unlinked and then linked data
            for n, neutrality in enumerate(['non-neutral', 'neutral']):
                ax.hist(data_dict[l][metric][n], bins = 50, alpha=0.7,
                        label='shift: ' + neutrality, color=hist_cols.pop())
            for n, neutrality in enumerate(['non-neutral', 'neutral']):
                ax.hist(null_data_dict_dict[linkage][l][metric][n], bins=50,
                        alpha=0.7, label='null: ' +
                        neutrality,color=hist_cols.pop())
            if ax_num == 18:
                ax.legend(prop={'size': 10})
            if ax_num in [1, 7, 13]:
                ax.set_ylabel(metric_text_dict[metric])
            if 1 <= ax_num <= 6:
                ax.set_xticks(dir_tick_locs)
                ax.set_xticklabels(dir_tick_labs)
            if 7 <= ax_num <= 12:
                ax.set_xlim([0, 50])
            elif 13 <= ax_num:
                ax.set_xlim([0,3])
            ax_num += 1

# show all the figures
for fig in fig_dict.values():
    fig.show()

if map_gene_flow:
    for fig in gene_flow_fig_dict.values():
        fig.show()

gf_hist_fig.show()
