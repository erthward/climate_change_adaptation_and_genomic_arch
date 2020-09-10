#!/usr/bin/python
#prelim.py

# flake8: noqa

'''
TODO:
    - start test-running and debugging
'''

import geonomics as gnx
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import copy
import tskit
import bisect
import time


#/\/\/\/\/\/\/\/\
# set main params
#\/\/\/\/\/\/\/\/

#------------
# main params
#------------
# set number of iterations for each sim
n_its = 2
# set the different numbers of loci to use
genicities = [5, 20, 100]
# set the different linkage levels to use
linkages = ['none', 'weak', 'strong']
linkages_dict = {'none': {'r_distr_alpha': 0.5,
                     'r_distr_beta': None},
            'weak': {'r_distr_alpha': 0.05,
                     'r_distr_beta': None},
            'strong': {'r_distr_alpha': 0.005,
                       'r_distr_beta': None}
           }

#---------------------------------------------------
# params to reduce runtime for debugging/development
#---------------------------------------------------
# set time when environmental change begins
change_T = 10 #500
# set total time over which environmental change takes place
T = 20 #800
# calc length of environmental change period
deltaT_env_change = T - change_T
# reset the K_factor (if desired)
K_factor=0.25

#--------------------------
# params for output control
#--------------------------
verbose = True
script_println_header = '=====>'

#----------------------------------
# params for stats to be calculated
#----------------------------------
gene_flow_stats = ['dir', 'dist', 'speed']
other_stats = ['Nt', 'fit']
stats = gene_flow_stats + other_stats
use_individs_curr_pos=True

#-------------------
# params for figures
#-------------------
# flag indicating whether or not to plot the lineage map
map_gene_flow = True
# set colors
colors = {'null': {'nonneut': '#237ade',  # dark blue
                   'neut': '#94b3d6'},    # light blue
          'main': {'nonneut': '#d60957',  # dark rose
                   'neut': '#d98daa'}}    # light rose
rowlab_size = 16
collab_size = 16


#----------
# filepaths
#----------
# create ParamsDict objects
filepath=('/home/drew/Desktop/stuff/berk/research/projects/sim/'
          'ch2_adapt_clim_chng_genarch/prelim_params.py')

#---------------------------------
# read, tweak, and copy the params
#---------------------------------
params = gnx.read_parameters_file(filepath=filepath)
# tweak the carrying capacity, if requested
if K_factor is not None:
    params['comm']['species']['spp_0']['init']['K_factor'] = K_factor


#/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
# create output data structure
#\/\/\/\/\/\/\/\/\/\/\/\/\/\/\

# output data structure for the gene flow stats
output = {nullness: {linkage: {genicity: {stat: {it: {neutrality: []
                                                      for neutrality in ['neut',
                                                                     'nonneut']}
                                                 for it in range(n_its)}
                                          for stat in gene_flow_stats}
                               for genicity in genicities}
                     for linkage in linkages}
          for nullness in ['main', 'null']}
# add sections for the non-gene flow stats
for nullness in ['main', 'null']:
    for linkage in linkages:
        for genicity in genicities:
            for stat in other_stats:
                output[nullness][linkage][genicity][stat] = {}
                for n_it in range(n_its):
                    output[nullness][linkage][genicity][stat][n_it] = []




#/\/\/\/\/\/\/\/
# create figures
#\/\/\/\/\/\/\/\

#------------------
# timecourse figure
#------------------
fig_time = plt.figure('time')
fig_time.suptitle(('example simulation time-courses, across linkage-genicity'
                   ' combinations'))
gs_time = mpl.gridspec.GridSpec(2*len(linkages), 2*len(genicities))
# use dicts to store the row and col indices in the gridspec that
# correspond with each of the linkage-genicity combo scenarios
gs_time_row_idxs = {linkage:i*2 for i, linkage in enumerate(linkages)}
gs_time_col_idxs = {genicity:i*2 for i, genicity in enumerate(genicities)}
# and use another dict to store the (row,col) offsets for each of the
# 4 plots to be created for each linkage-genicity combo scenario
gs_time_rowcol_idx_offsets = {'before': (0,0),
                              'after': (0,1),
                              'Nt': (1,0),
                              'fit': (1,1)}

def get_fig_time_row_col_idxs(linkage, genicity, plot_type):
    row_idx = gs_time_row_idxs[linkage]
    col_idx = gs_time_col_idxs[genicity]
    offsets = gs_time_rowcol_idx_offsets[plot_type]
    row_idx += offsets[0]
    col_idx += offsets[1]
    return (row_idx, col_idx)

#------------
# hist figure
#------------
fig_hist = plt.figure('hist')
fig_hist.suptitle(('histograms of gene-flow direction, distance,'
                  ' and speed, across linkage-genicity combinations'))
gs_hist = mpl.gridspec.GridSpec(len(linkages), 3*len(genicities))


#/\/\/\/\/\/\/\/\/
# define functions
#\/\/\/\/\/\/\/\/\

def save_data(nullness, genicity, linkage, n_it, mod, output, max_time_ago,
              fit_data):
    '''
    Save the current model's neutral and non-neutral locus data to its
    appropriate spot in the output dict
    '''
    # grab the non-neutral loci
    nonneut_loci = mod.comm[0].gen_arch.traits[0].loci
    # grab an equal amt of random neutral loci to the num of non-neutral loci
    neut_loci = np.random.choice(mod.comm[0].gen_arch.neut_loci, genicity,
                                 replace=False)
    # calculate gene-flow stats for both the non-neutral and neutral loci
    nonneut_stats = mod.comm[0]._calc_lineage_stats(stats=['dir', 'dist',
                                                           'speed'],
                                    use_individs_curr_pos=use_individs_curr_pos,
                                                    max_time_ago=max_time_ago,
                                                    loci=nonneut_loci)
    neut_stats = mod.comm[0]._calc_lineage_stats(stats=['dir', 'dist', 'speed'],
                                    use_individs_curr_pos=use_individs_curr_pos,
                                                 max_time_ago=max_time_ago,
                                                 loci=neut_loci)
    # save the gene-flow data in the right places
    for stat, dataset in nonneut_stats.items():
        # add the non-neutral stats to their lists
        for data in dataset.values():
            stat_list=output[nullness][linkage][genicity][stat][n_it]['nonneut']
            stat_list.extend(data)
        # add the neutral data too
        for data in neut_stats[stat].values():
            stat_list = output[nullness][linkage][genicity][stat][n_it]['neut']
            stat_list.extend(data)

    # save the non-gene-flow data in the right places
    output[nullness][linkage][genicity]['Nt'][n_it].extend(mod.comm[0].Nt)
    output[nullness][linkage][genicity]['fit'][n_it].extend(fit_data)


def set_params(params, linkage, genicity, nullness):
    copy_params = copy.deepcopy(params)
    # set the model name (so that it saves data correctly in separate dirs)
    model_name = '_%s'.join([nullness, linkage, str(genicity), str(n_it)])
    model_name = model_name % ('L', 'G', 'I')
    if verbose:
        print('\n%sNOW RUNNING MODEL %s...\n' % (script_println_header,
                                               model_name))
    copy_params.model['name'] = model_name

    # set the linkage params correctly
    r_distr_params = linkages_dict[linkage]
    for k,v in r_distr_params.items():
        copy_params['comm']['species']['spp_0']['gen_arch'][k] = v

    # set the num of loci, and set the effect-size distribution
    # and the total genome length to match the num of loci
    copy_params['comm']['species']['spp_0']['gen_arch']['L'] = 3 * genicity
    for trt in copy_params['comm']['species']['spp_0']['gen_arch']['traits'].values():
        trt['n_loci'] = genicity
        trt['alpha_distr_mu'] = 1/genicity

    # if this is a null sim, edit the params so that
    # no env change occurs and there is no gradient
    if nullness == 'null':
        copy_params['landscape']['layers']['shift'].pop('change');
        copy_params['landscape']['layers']['shift']['init'][
            'defined']['rast'] = np.ones((50, 50))*0.5
        copy_params['landscape']['layers']['stable']['init'][
            'defined']['rast'] = np.ones((50, 50))*0.5

    return copy_params


def run_sim(nullness, linkage, genicity, n_its, params, output):
    '''Run the simulations for a given set of parameters and hyperparameters
    '''
    # get the correct params for this sim
    copy_params = set_params(params, linkage, genicity, nullness)

    # loop through the iterations
    for n_it in range(n_its):
        print('\n%sITERATION NUMBER %i\n' % (script_println_header, n_it))

        # empty list for the fitness data for this model
        fit_data = []

        # create the model
        mod = gnx.make_model(gnx.make_params_dict(copy_params))

        # check that the number of trait 0 loci is correct
        assert len(mod.comm[0].gen_arch.traits[0].loci) == genicity, (
                'EXPECETED %i LOCI BUT GOT %i!') % (
                genicity, len(mod.comm[0].gen_arch.traits[0].loci))

        # burn the model in
        mod.walk(10000, mode='burn')

        # run the model up to the env change event
        for t in range(change_T):

            # keep printing the number of loci,
            # to help me track things while it's running
            n_loci = len(mod.comm[0].gen_arch.traits[0].loci)
            print('\n%s[%s, %i loci, it %i]' % (script_println_header,
                                                linkage, genicity, n_it))
            # walk 1 step
            mod.walk(1, mode='main')

            # save the fitness data for this timestep
            fit_data.append(np.mean(mod.comm[0]._get_fit()))

        # add the pre-change population to the fig,
        # if this is the first iteration
        if n_it == 0:
            # before-change phenotypic map
            row_idx, col_idx = get_fig_time_row_col_idxs(linkage, genicity,
                                                         'before')
            ax = fig_time.add_subplot(gs_time[row_idx, col_idx])
            #mod.plot_phenotype(0, 0, alpha=0.5)
            ax.imshow(mod.land[0].rast, cmap='coolwarm')
            ax.scatter(mod.comm[0]._get_x()-0.5,
                       mod.comm[0]._get_y()-0.5,
                       c=mod.comm[0]._get_z(0),
                       cmap='coolwarm',
                       s=25)
            ax.set_xticks([])
            ax.set_yticks([])

            #set the row and/or col labels, if needed
            if col_idx == 0 and row_idx in [0, 2, 4]:
                ax.set_ylabel('linkage: %s' % str(linkage), size=rowlab_size)
            if row_idx == 0 and col_idx in [0, 2, 4]:
                ax.set_title('genicity: %s' % str(genicity), size=collab_size)

        # run the model through the env change event
        for t in range(change_T, T):

            # keep printing the number of loci,
            # to help me track things while it's running
            n_loci = len(mod.comm[0].gen_arch.traits[0].loci)
            print('\n%s[%s, %i loci, it %i]' % (script_println_header,
                                                linkage, genicity, n_it))
            # walk 1 step
            mod.walk(1, mode='main')

            # save the fitness value
            fit_data.append(np.mean(mod.comm[0]._get_fit()))

        # save the data
        save_data(nullness, genicity, linkage, n_it, mod, output,
                  deltaT_env_change, fit_data)

        # add the post-change population and other plots to the fig,
        # if this is the first iteration
        if n_it == 0:
            # after-change phenotypic map
            row_idx, col_idx = get_fig_time_row_col_idxs(linkage, genicity,
                                                         'after')
            ax = fig_time.add_subplot(gs_time[row_idx, col_idx])
            #mod.plot_phenotype(0, 0, alpha=0.5)
            ax.imshow(mod.land[0].rast, cmap='coolwarm')
            ax.scatter(mod.comm[0]._get_x()-0.5,
                       mod.comm[0]._get_y()-0.5,
                       c=mod.comm[0]._get_z(0),
                       cmap='coolwarm',
                       s=25)
            ax.set_xticks([])
            ax.set_yticks([])

            # pop-growth plot
            row_idx, col_idx = get_fig_time_row_col_idxs(linkage, genicity,
                                                         'Nt')
            ax = fig_time.add_subplot(gs_time[row_idx, col_idx])
            ts = [*np.linspace(-len(mod.comm[0].Nt), T-1, len(mod.comm[0].Nt))]
            x0 = mod.comm[0].Nt[0]/mod.comm[0].K.sum() 
            logistic = [gnx.structs.species._calc_logistic_soln(x0,
                                    mod.comm[0].R, t) for t in range(len(ts))]
            ax.plot(ts, logistic, color=colors[nullness]['neut'], alpha=0.5)
            ax.plot(ts, mod.comm[0].Nt, color=colors[nullness]['nonneut'])
            #mod.plot_pop_growth(0, expected_color='gray',
            #                    actual_color=colors[nullness]['nonneut'])
            ax.plot([change_T]*2, [0, max(mod.comm[0].Nt)], ':k')
            ax.plot([-1]*2, [0, max(mod.comm[0].Nt)], '-k')
            ax.set_xlabel('t', size=8)
            ax.set_ylabel('N(t)', size=8)

            # mean fitness change plot
            row_idx, col_idx = get_fig_time_row_col_idxs(linkage, genicity,
                                                         'fit')
            ax = fig_time.add_subplot(gs_time[row_idx, col_idx])
            ax.plot(range(len(fit_data)), fit_data,
                    color=colors[nullness]['nonneut'])
            ax.plot([change_T]*2, [0, 1], ':k')
            # TODO: perhaps change this to calculate the min y val in some way?
            ax.set_ylim([0.5, 1])
            ax.set_xlabel('t', size=8)
            ax.set_ylabel('mean fitness', size=8)

            # TODO: calc and store neut vs. non-neut tests??

            # TODO: still do something like this figure?
            # make the gene-flow figure
            #plt.figure('gene flow %i' % l)
            #gf_ax1 = gf_fig.add_subplot(121)
            #gf_ax1.set_title('Unlinked, trait 0', size=16)
            ## define number of neutral and non-neutral loci to randomly draw
            #n_locs=5
            #nonneut_loci_to_plot = np.random.choice(
            #            mod.comm[0].gen_arch.traits[0].loci, n_locs,
            #            replace=False)

            #neut_loci_to_plot = np.random.choice(
            #            mod.comm[0].gen_arch.neut_loci, n_locs,
            #            replace=False)
            #nonneut_individs = np.random.choice([*mod.comm[0]], 10*n_locs,
            #                                    replace=False)
            #neut_individs = np.random.choice([*mod.comm[0]], 10*n_locs,
            #                                 replace=False)
            #for n, loc in enumerate(nonneut_loci_to_plot):
            #    mod.comm[0]._plot_gene_flow(loc, 'vector', mod.land,
            #            individs=nonneut_individs[n*n_locs:n*n_locs+n_locs],
            #            color='#616161', phenotype=0, size=35)
            #for n, loc in enumerate(neut_loci_to_plot):
            #    mod.comm[0]._plot_gene_flow(loc, 'vector', mod.land,
            #            individs=nonneut_individs[n*n_locs:n*n_locs+n_locs],
            #            color='#f2f2f2', phenotype=0, size=35)



#/\/\/\/\/\/\/\/\/\/\/\/\
# loop through iterations
#\/\/\/\/\/\/\/\/\/\/\/\/

for linkage in linkages:
    for genicity in genicities:

        #-----------------------------------------
        # run simulation with environmental change
        #-----------------------------------------
        run_sim('main', linkage, genicity, n_its, params, output)
        
        #--------------------
        # run null simulation
        #--------------------
        run_sim('null', linkage, genicity, n_its, params, output)

        #-------------------------
        # TODO null/non-null comparison
        #-------------------------
        
        # calc and store null vs non-null comparison


# TODO: run overall tests/summaries of stored test results





#/\/\/\/\/\/\/\
# finalize figs
#\/\/\/\/\/\/\/

#dir_tick_locs = np.linspace(0, 360, 9)
dir_tick_locs = np.linspace(0, 360, 5)
#dir_tick_labs = ['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW', 'N']
dir_tick_labs = ['N', 'E', 'S', 'W', 'N']

for genicity_n, genicity in enumerate(genicities):
    for linkage_n, linkage in enumerate(linkages):
        for stat_n, stat in enumerate(gene_flow_stats):
            row_idx = linkage_n
            col_idx = (genicity_n*3)+stat_n
            ax = fig_hist.add_subplot(gs_hist[row_idx, col_idx])
            if col_idx == 0:
                ax.set_ylabel('linkage: %s' % str(linkage), size=rowlab_size)
            if row_idx == 0 and col_idx in [0, 3, 6]:
                ax.set_title('genicity: %i' % genicity, size=collab_size)
            # plot a histogram of the unlinked and then linked data
            for nullness in ['main', 'null']:
                for neut_idx, neutrality in enumerate(['neut', 'nonneut']):
                    data_dict = output[nullness][linkage][genicity][stat]
                    data = [v[neutrality] for v in data_dict.values()]
                    data = [val for sublist in data for val in sublist]
                    data = [np.nan if val is None else val for val in data]
                    ax.hist(data, bins = 50, alpha=0.5,
                            label= '%s: %s' % (nullness, neutrality),
                            color=colors[nullness][neutrality])
                    ax.set_xlabel(stat, size=8)
                    if row_idx==2 and col_idx==8:
                        ax.legend(prop={'size': 10})
                    if col_idx in [0, 3, 6]:
                        ax.set_xticks(dir_tick_locs)
                        ax.set_xticklabels(dir_tick_labs)
                    # TODO standardize axes

# show all the figures
fig_time.show()
fig_hist.show()
