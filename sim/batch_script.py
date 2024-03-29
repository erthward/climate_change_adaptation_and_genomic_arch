#!/usr/bin/python
#batch_script.py

# flake8: noqa

from collections import Counter as C
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

#/\/\//\/\/\
# set params
#\/\/\/\/\/\

#------------
# main params
#------------
# set number of iterations for each sim
n_its = 1

# reset the K_factor
K_factor = 2.5

# set the different numbers of loci to use,
# and the factor to multiply effect size by
# to either implement high genetic redundancy or not
#redundancy = 'hi'
redundancy = 'lo'
assert redundancy in ['lo', 'hi']
if redundancy == 'lo':
    genicities = [4, 20, 100]
    alpha_factor = 1
elif redundancy == 'hi':
    genicities = [8, 40, 200]
    alpha_factor = 2
print('\n\n%s REDUND:\n\tGENICITIES: %s\nα FACTOR: %i\n\n' % (redundancy,
                                                              str(genicities),
                                                              alpha_factor))

# set the different linkage levels to use
linkages = ['strong', 'weak', 'independent']
linkages_dict = {'independent': {'r_distr_alpha': 0.5,
                                 'r_distr_beta': None},
                 'weak': {'r_distr_alpha': 0.05,
                          'r_distr_beta': None},
                 'strong': {'r_distr_alpha': 0.005,
                            'r_distr_beta': None}
           }


#-------------------------------------
# params and settings for data storage
#-------------------------------------

# get process ID as a string
pid = str(os.getpid())

# prep the counts table table (for regular printout)
border_patt = ">~~<"
max_len_gen_str = max([len(str(gen)) for gen in genicities])
gen_row_headers = ["%s ||" % (str(gen) + (" " * (max_len_gen_str -
                                        len(str(gen))))) for gen in genicities]
gen_rows = [head + (" nul: %s  non: %s |" * len(
                                        linkages)) for head in gen_row_headers]
max_len_it_str = len(str(n_its))
linkage_header = (" " * (max_len_gen_str+1)) + "||"
linkage_header += "".join(([" " + link.upper() + ":" + (" " *
                (13 - len(link) +max_len_it_str)) + "|" for link in linkages]))
pid_header = "FROM PID: %s\n" % pid
pid_header = pid_header + "="*(len(pid_header)-1) + '\n'
cts_table = "\n".join([linkage_header] + gen_rows) + "\n"
cts_table = (("\n"*3) + (border_patt*20) + ("\n"*3) + pid_header +
             cts_table + ("\n"*2))
cts_table_list = [0] * (2 * len(linkages) * len(genicities))


#-----------------------------------------------------
# params to reduce runtime/for debugging & development
#-----------------------------------------------------

# get timesteps
if os.getcwd().split('/')[1] == 'home':
    steps = pd.read_csv(('/home/deth/Desktop/CAL/research/projects/sim/'
                         'ch2/climate_change_adaptation_and_genomic_arch/sim/'
                         'time_steps.CSV'))
# or else get filepaths on Savio
else:
    steps = pd.read_csv(('/global/scratch/users/drewhart/'
                         'ch2/climate_change_adaptation_and_genomic_arch/sim/'
                         'time_steps.CSV'))

# set time when environmental change begins
change_T = int(steps[steps['name']=='start']['num'].values[0])

# set time when environmental change ends
T = int(steps[steps['name']=='end']['num'].values[0])
print('climate change will start at time step %i' % change_T)
print('simulation will end at time step %i' % T)

# calc length of environmental change period
deltaT_env_change = T - change_T

# calc time prior to climate change to start storing Nt, mean fitness, and mean
# phenotype
t_record_data_b4_env_change = change_T - deltaT_env_change

# print out debugging info?
debug_verbose = False

#--------------------------
# params for stdout control
#--------------------------
mod_verbose = False
script_println_header = '=====>'

#----------------------------------
# params for stats to be calculated
#----------------------------------
gene_flow_stats = ['dir']
other_stats = ['Nt', 'fit']
stats = gene_flow_stats + other_stats
use_individs_curr_pos=True

#-------------------
# params for figures
#-------------------

make_plots = False
# set colors
colors = {'null': {'nonneut': '#237ade',  # dark blue
                   'neut': '#94b3d6'},    # light blue
          'non-null': {'nonneut': '#d60957',  # dark rose
                   'neut': '#d98daa'}}    # light rose
suptitle_size = 7
rowlab_size = 6
collab_size = 6
ticklabel_size = 6


#--------------------------
# filepaths and file-saving
#--------------------------
# get filepaths on my laptop
if os.getcwd().split('/')[1] == 'home':
    params_filepath=('/home/deth/Desktop/CAL/research/projects/sim/'
              'ch2/climate_change_adaptation_and_genomic_arch/sim/template_params.py')
    # path to dir for output CSVs
    output_path = ('/home/drew/Desktop/stuff/berk/research/projects/sim/'
               'ch2_adapt_clim_chng_genarch')
# or else get filepaths on Savio
else:
    params_filepath=(('/global/scratch/users/drewhart/ch2/'
                      'climate_change_adaptation_and_genomic_arch/sim/template_params.py'))
    # path to dir for output CSVs
    with open(('/global/scratch/users/drewhart/ch2/climate_change_adaptation_'
               'and_genomic_arch/analysis/outputdir.txt'), 'r') as f:
        output_path = f.read().strip()
    #output_path = '/global/scratch/users/drewhart/ch2/output/output'

# make sure job is running in the right directory
#os.chdir(output_path)
#print("NOW IN THIS DIR! ", os.getcwd())

#---------------------------------
# read, tweak, and copy the params
#---------------------------------
# create ParamsDict object
params = gnx.read_parameters_file(filepath=params_filepath)
data_times = params['model']['data']['sampling']['when']
print('data will be collected at timesteps %i, %i, and %i' % (data_times[0],
                                                           data_times[1],
                                                           data_times[2],
                                                          ))
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
          for nullness in ['non-null', 'null']}
# add sections for the non-gene flow stats
for nullness in ['non-null', 'null']:
    for linkage in linkages:
        for genicity in genicities:
            for stat in other_stats:
                output[nullness][linkage][genicity][stat] = {}
                for n_it in range(n_its):
                    output[nullness][linkage][genicity][stat][n_it] = []

#/\/\/\/\/\/\/\/
# create figures
#\/\/\/\/\/\/\/\

if make_plots:
    #------------------
    # timecourse figure
    #------------------

    fig_time = {'non-null': plt.figure('time'),
                'null': plt.figure('time, null')}
    [fig.suptitle(('example simulation time-courses, across linkage-genicity'
                   ' combinations; NULLNESS:'
                   ' %s' % nullness)) for nullness, fig in zip(['non-null', 'null'],
                                                               fig_time.values())]

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
                      ' and speed, across linkage-genicity combinations'),
                      fontsize=suptitle_size)
    # NOTE: extra column gives space for legend to live in; hacky, but works
    gs_hist = mpl.gridspec.GridSpec(len(linkages), 3*len(genicities)+1)


#/\/\/\/\/\/\/\/\/
# define functions
#\/\/\/\/\/\/\/\/\

def convert_compass_ang_2_quadrant_ang(ang, in_rads=True):
    out_ang = -ang+90
    if out_ang <0:
        out_ang = out_ang+360
    if in_rads:
        out_ang = 2*np.pi*(out_ang/360)
    return out_ang


def calc_NSness_Eness(dirs):
    # calculate and store 'north/southness' and 'eastness'
    sines = [np.sin(convert_compass_ang_2_quadrant_ang(d)) for d in dirs]
    cosines = [np.cos(convert_compass_ang_2_quadrant_ang(d)) for d in dirs]
    NSness = np.mean(np.abs(sines))
    Eness = np.mean([d for d in cosines if d >= 0])
    return (NSness, Eness)


def store_data(nullness, genicity, linkage, n_it, mod, output, max_time_ago,
              fit_data):
    '''
    Store the current model's non-neutral locus data to its
    appropriate spot in the output dict
    '''
    # output data structure for summary stats
    stats_output = {'dir': {k: np.nan for k in ['NSness', 'Eness']}}

    # grab the non-neutral loci
    nonneut_loci = mod.comm[0].gen_arch.traits[0].loci
    # and reduce them down to a number of loci equal to half of the minimum
    # genicity (to ensure equal sample sizes in my ultimate fitted von Mises
    # distributions),
    # using only positive-effect loci (since the typical expectation is that
    # they would flow upslope (i.e., east)
    sample_loci = mod.comm[0].gen_arch.traits[0].loci[
        np.random.choice(np.where(mod.comm[0].gen_arch.traits[0].alpha>0)[0],
                         size=2, replace=False)]

    # calculate gene-flow stats for the non-neutral loci
    # NOTE: each stat is structured as:
    # {'dir': {loc_0: [dir_ind0_chrom0, dir_ind0_chrom1,
    #                  dir_ind1_chrom0, ..., dir_indN_chrom1],
    #          loc_1:
    #          ...
    #          loc_n},
    # }
    # in other words, for each stat, for each locus, gene flow stats for all
    # chromosomes in current pop
    nonneut_stats = mod.comm[0]._calc_lineage_stats(stats=['dir'],
                                    use_individs_curr_pos=use_individs_curr_pos,
                                                    max_time_ago=max_time_ago,
                                                    loci=sample_loci)
    # store the gene-flow data in the right places
    for stat, dataset in nonneut_stats.items():
        nonneut_data = []
        # add the non-neutral stats to their lists
        for data in dataset.values():
            stat_list=output[nullness][linkage][genicity][stat][n_it]['nonneut']
            stat_list.extend(data)
            nonneut_data.extend(data)

        nonneut_data = [val if val is not None else np.nan for val in nonneut_data]

        # calculate and store the 'north/southness' and 'eastness'
        NSness, Eness = calc_NSness_Eness(nonneut_data)
        stats_output[stat]['NSness'] = NSness
        stats_output[stat]['Eness'] = Eness

    # store the non-gene-flow data in the right places
    output[nullness][linkage][genicity]['Nt'][n_it].extend(mod.comm[0].Nt)
    output[nullness][linkage][genicity]['fit'][n_it].extend(fit_data)

    return(stats_output)


def make_custom_genarch_file(genicity, recomb_rate, pid, write=True):
    assert(isinstance(genicity, int))
    L = int(genicity*2)
    assert(L%4==0)
    cols = ['locus', 'p', 'dom', 'r', 'trait', 'alpha']
    locus = [*range(L)]
    p = [0.5]*L # all loci start at MAF==0.5
    dom = [0]*L # all codominant
    r = [0]+[recomb_rate]*(L-1) # first locus must be 0
    trait = ['trait_0', 'trait_1']*int((L/2)) # alternate traits
    alpha = [1/genicity,
             1/genicity,
             -1/genicity,
             -1/genicity]*int(L/4) # alpha=1/genicity reaches all phenotypes [0,1] precisely
                                   # (greater than = redundancy @ extremes)
    alpha = [alpha_factor*n for n in alpha]
    df = pd.DataFrame.from_dict(dict(zip(cols, [locus, p, dom, r, trait, alpha])))
    if write:
        tmp_filename = 'tmp_genarch_%s.csv' % pid
        df.to_csv(tmp_filename, index=False)
    return tmp_filename


def set_params(params, linkage, genicity, nullness):
    copy_params = copy.deepcopy(params)
    #create random number ID, so that directories don't get overwritten
    rand_id = int(np.random.uniform()*1e7)

    # set the model name (so that it saves data correctly in separate dirs)
    model_name = ('_%s'.join([nullness, linkage, str(genicity), str(n_it), str(rand_id)]) +
                   "PID-%s" % pid)
    model_name = model_name % ('L', 'G', 'its', 'randID')
    if debug_verbose:
        print('\n%sNOW RUNNING MODEL %s...\n' % (script_println_header,
                                               model_name))
    copy_params['model']['name'] = model_name

    # create the custom genomic architecture file and connect it to the params
    recomb_rate = linkages_dict[linkage]['r_distr_alpha']
    genarch_filename = make_custom_genarch_file(genicity, recomb_rate, pid)
    copy_params['comm']['species']['spp_0']['gen_arch'][
                                        'gen_arch_file'] = genarch_filename

    # set the linkage params correctly
    r_distr_params = linkages_dict[linkage]
    for k,v in r_distr_params.items():
        copy_params['comm']['species']['spp_0']['gen_arch'][k] = v

    # set the num of loci, and set the effect-size distribution
    # and the total genome length to match the num of loci
    copy_params['comm']['species']['spp_0']['gen_arch']['L'] = 2 * genicity
    for trt in copy_params['comm']['species']['spp_0']['gen_arch']['traits'].values():
        trt['n_loci'] = genicity
        trt['alpha_distr_mu'] = alpha_factor*(1/genicity)

    # if this is a null sim, edit the params so that
    # no env change occurs 
    if nullness == 'null':
        copy_params['landscape']['layers']['shift'].pop('change');

    return copy_params, genarch_filename


def run_sim(nullness, linkage, genicity, n_its, params, output,
            cts_table_list, cts_table_list_idx):
    '''
    Run the simulations for a given set of parameters and hyperparameters
    '''
    # get the correct params for this sim
    copy_params, genarch_filename = set_params(params, linkage,
                                               genicity, nullness)
    assert str(pid) in copy_params['model']['name'], ('PID MISSING')

    # set up stats lists to be returned
    delta_Nt = []
    delta_fit = []
    NSness = []
    Eness = []

    # Nt, mean fitness, and mean phenotype
    # (to be recorded for each time step during climate change, and for an 
    # equal-length petiod beforehand, for viz and stats analysis)
    Nt_list = []
    mean_fit_list = []
    mean_z_list = []

    # loop through the iterations
    for n_it in range(n_its):

        # empty list for the fitness data for this model
        fit_data = []

        # create the model
        mod = gnx.make_model(gnx.make_params_dict(copy_params),
                             name=copy_params['model']['name'])
        assert 'unnamed' not in mod.name, 'STILL UNNAMED!'

        # delete the temporary genarch file
        os.remove(genarch_filename)

        # check that loci alternate between traits 0 and 1
        assert np.all(mod.comm[0].gen_arch.traits[1].loci ==
                      mod.comm[0].gen_arch.traits[0].loci+1), ('Traits\' loci '
                                                               'do not '
                                                               'appear to be '
                                                               'intercalated!')

        # save the original carrying capacity raster (to use in plotting later)
        if make_plots:
            orig_K = np.copy(mod.comm[0].K)

        # check that the number of trait 0 loci is correct
        assert len(mod.comm[0].gen_arch.traits[0].loci) == genicity, (
                'EXPECETED %i LOCI BUT GOT %i!') % (
                genicity, len(mod.comm[0].gen_arch.traits[0].loci))

        # print info to verify parameterization is correct, if requested
        # and on the first iteration
        if debug_verbose and n_it == 0:
            gen_arch_params = copy_params['comm']['species']['spp_0']['gen_arch']
            genicities = [gen_arch_params['traits']['trait_0']['n_loci'],
                          gen_arch_params['traits']['trait_1']['n_loci']]
            print('\t\t\tGENICITIES: %s' % str(genicities))
            print(('\t\t\tNON-NEUT GENES:\n'
                   '\t\t\t\tt0:%s\n'
                   '\t\t\t\tt1:%s') % (str(mod.comm[0].gen_arch.traits[0].loci),
                                      str(mod.comm[0].gen_arch.traits[1].loci)))

            linkage_nums = [gen_arch_params['r_distr_alpha'],
                            gen_arch_params['r_distr_beta']]

        #print iteration number
        if debug_verbose:
            print('\n%sITERATION NUMBER %i\n' % (script_println_header, n_it))

        # burn the model in
        mod.walk(10000, mode='burn', verbose=mod_verbose)

        # run the model up to the env change event
        for t in range(change_T):

            # record the Nt, mean fit, and mean z, if <= deltaT_env_change
            # prior to climate change
            if t >= t_record_data_b4_env_change:
                Nt_list.append(mod.comm[0].Nt[-1])
                mean_fit_list.append(np.mean(mod.comm[0]._get_fit()))
                mean_z_list.append(np.mean(mod.comm[0]._get_z()[:,0]))

            # keep printing the number of loci,
            # to help me track things while it's running
            n_loci = len(mod.comm[0].gen_arch.traits[0].loci)
            if debug_verbose:
                print('\n%s[%s, %i loci, it %i]' % (script_println_header,
                                                    linkage, genicity, n_it))
            # walk 1 step
            mod.walk(1, mode='main', verbose=mod_verbose)

            # store the fitness data for this timestep
            fit_data.append(np.mean(mod.comm[0]._get_fit()))

        # calculate the mean vars before the shift
        mean_Nt_b4 = np.mean(mod.comm[0].Nt[-deltaT_env_change:])
        mean_fit_b4 = np.mean(fit_data[-deltaT_env_change:])

        # add the pre-change population to the fig,
        # if this is the first iteration and plots are requested
        if make_plots and n_it == 0:


            # before-change phenotypic map
            row_idx, col_idx = get_fig_time_row_col_idxs(linkage, genicity,
                                                         'before')
            ax = fig_time[nullness].add_subplot(gs_time[row_idx, col_idx])
            #mod.plot_phenotype(0, 0, alpha=0.5)
            ax.imshow(mod.land[0].rast, cmap='coolwarm', vmin=0, vmax=1)
            ax.scatter(mod.comm[0]._get_x()-0.5,
                       mod.comm[0]._get_y()-0.5,
                       c=mod.comm[0]._get_z(0),
                       cmap='coolwarm',
                       vmin=0, vmax=1, s=5, alpha=0.2)
            ax.set_xticks([])
            ax.set_yticks([])

            #set the row and/or col labels, if needed
            if col_idx == 0 and row_idx in [0, 2, 4]:
                ax.set_ylabel('linkage: %s' % (str(linkage) + ': ' + str(mod.comm[0].gen_arch.recombinations._rates[3])), size=rowlab_size)
                #ax.set_ylabel('linkage: %s' % str(linkage), size=rowlab_size)
            if row_idx == 0 and col_idx in [0, 2, 4]:
                ax.set_title('genicity: %s' % str(genicity), size=collab_size)

        # run the model through the env change event
        for t in range(change_T, T):

            # record the Nt, mean fit, and mean z
            Nt_list.append(mod.comm[0].Nt[-1])
            mean_fit_list.append(np.mean(mod.comm[0]._get_fit()))
            mean_z_list.append(np.mean(mod.comm[0]._get_z()[:,0]))

            # keep printing the number of loci,
            # to help me track things while it's running
            n_loci = len(mod.comm[0].gen_arch.traits[0].loci)
            if debug_verbose:
                print('\n%s[%s, %i loci, it %i]' % (script_println_header,
                                                    linkage, genicity, n_it))
            # walk 1 step
            mod.walk(1, mode='main', verbose=mod_verbose)

            # store the fitness value
            fit_data.append(np.mean(mod.comm[0]._get_fit()))

        # record the Nt, mean fit, and mean z
        Nt_list.append(mod.comm[0].Nt[-1])
        mean_fit_list.append(np.mean(mod.comm[0]._get_fit()))
        mean_z_list.append(np.mean(mod.comm[0]._get_z()[:,0]))

        # calculate and store delta_Nt and delta_fit
        mean_Nt_af = np.mean(mod.comm[0].Nt[-deltaT_env_change:])
        mean_fit_af = np.mean(fit_data[-deltaT_env_change:])
        delta_Nt.append(mean_Nt_af - mean_Nt_b4)
        delta_fit.append(mean_fit_af - mean_fit_b4)

        # store the data
        stats_output = store_data(nullness, genicity, linkage, n_it,
                                  mod, output, deltaT_env_change, fit_data)

        # store the summary-stats output
        dir_stats = stats_output['dir']
        NSness.append(dir_stats['NSness'])
        Eness.append(dir_stats['Eness'])

        # add the post-change population and other plots to the fig,
        # if this is the first iteration and plots are requested
        if make_plots and n_it == 0:


            # after-change phenotypic map
            row_idx, col_idx = get_fig_time_row_col_idxs(linkage, genicity,
                                                         'after')
            ax = fig_time[nullness].add_subplot(gs_time[row_idx, col_idx])
            #mod.plot_phenotype(0, 0, alpha=0.5)
            ax.imshow(mod.land[0].rast, cmap='coolwarm', vmin=0, vmax=1)
            ax.scatter(mod.comm[0]._get_x()-0.5,
                       mod.comm[0]._get_y()-0.5,
                       c=mod.comm[0]._get_z(0),
                       cmap='coolwarm',
                       vmin=0, vmax=1, s=5, alpha=0.2)
            ax.set_xticks([])
            ax.set_yticks([])


            # pop-growth plot
            row_idx, col_idx = get_fig_time_row_col_idxs(linkage, genicity,
                                                         'Nt')
            ax = fig_time[nullness].add_subplot(gs_time[row_idx, col_idx])
            ts = [*range(T)]
            x0 = mod.comm[0].Nt[-T]/orig_K.sum()
            logistic = [gnx.structs.species._calc_logistic_soln(x0,
                                                mod.comm[0].R, t) for t in ts]
            logistic = [n * orig_K.sum() for n in logistic]
            ax.plot(ts, logistic, color=colors[nullness]['neut'], alpha=0.5)
            ax.plot(ts, mod.comm[0].Nt[-T:], color=colors[nullness]['nonneut'])
            #mod.plot_pop_growth(0, expected_color='gray',
            #                    actual_color=colors[nullness]['nonneut'])
            ax.plot([change_T]*2, [0, max(mod.comm[0].Nt)], ':k')
            ax.set_xlabel('t', size=8)
            ax.set_ylabel('N(t)', size=8)


            # mean fitness change plot
            row_idx, col_idx = get_fig_time_row_col_idxs(linkage, genicity,
                                                         'fit')
            ax = fig_time[nullness].add_subplot(gs_time[row_idx, col_idx])
            ax.plot(range(len(fit_data)), fit_data,
                    color=colors[nullness]['nonneut'])
            ax.plot([change_T]*2, [0, 1], ':k')
            ax.set_ylim([0.5, 1])
            ax.set_xlabel('t', size=8)
            ax.set_ylabel('mean fitness', size=8)

        cts_table_list[cts_table_list_idx] += 1
        print(cts_table % tuple([str(n) + " " * (max_len_it_str -
                                        len(str(n))) for n in cts_table_list]))

        # run the model for another 250 time steps after the climate change
        # period, just to see if/where impacts on Nt, mean_fit, and mean_z
        # are lagged and/or persistent
        for t in range(T, T+deltaT_env_change):

            # record the Nt, mean fit, and mean z
            Nt_list.append(mod.comm[0].Nt[-1])
            mean_fit_list.append(np.mean(mod.comm[0]._get_fit()))
            mean_z_list.append(np.mean(mod.comm[0]._get_z()[:,0]))

            # walk 1 step
            mod.walk(1, mode='main', verbose=mod_verbose)
    return (delta_Nt, delta_fit,
            NSness, Eness,
            Nt_list, mean_fit_list, mean_z_list)


# gather gene flow data into DataFrames
def make_stat_df(stat, output):
    # create empty columns for gene flow df
    genicity_col = []
    linkage_col = []
    nullness_col = []
    neutrality_col = []
    it_col = []
    stat_col = []
    # loop through and get data
    for nullness, nullness_dict in output.items():
        for linkage, linkage_dict in nullness_dict.items():
            for genicity, genicity_dict in linkage_dict.items():
                for genicity, genicity_dict in linkage_dict.items():
                    for it, it_dict in genicity_dict[stat].items():
                        for neutrality, data in it_dict.items():
                            nrows = len(data)
                            genicity_col.extend([genicity]*nrows)
                            linkage_col.extend([linkage]*nrows)
                            nullness_col.extend([nullness]*nrows)
                            neutrality_col.extend([neutrality]*nrows)
                            it_col.extend([it]*nrows)
                            # replace Nones with NAs
                            data = [np.nan if v is None else v for v in data]
                            stat_col.extend(data)
    # make into a dataframe
    df = pd.DataFrame({'genicity': genicity_col,
                       'linkage': linkage_col,
                       'nullness': nullness_col,
                       'neutrality': neutrality_col,
                       'it': it_col,
                       stat: stat_col})
    return df



#/\/\/\/\/\/\/\/\/\/\/\/\
# loop through iterations
#\/\/\/\/\/\/\/\/\/\/\/\/

# create empty columns for Nt and fit dataframe
linkage_col = []
genicity_col = []
nullness_col = []
delta_Nt_col = []
delta_fit_col = []
NSness_col = []
Eness_col = []

# dict to become dataframe for mean phenotype data
ts_data_dict = {k:[] for k in ['time_step',
                               'nullness',
                               'genicity',
                               'linkage',
                               'Nt',
                               'mean_fit',
                               'mean_z']}

cts_table_list_idx = 0
cts_table_list = [0] * (2 * len(linkages) * len(genicities))
for genicity in genicities:
    for linkage in linkages:

        #-----------------------------------------
        # run simulation with environmental change
        #-----------------------------------------
        (delta_Nt, delta_fit, NSness, Eness,
         Nt_list, mean_fit_list, mean_z_list) = run_sim( 'non-null',
                                                        linkage,
                                                        genicity,
                                                        n_its,
                                                        params,
                                                        output,
                                                        cts_table_list,
                                                        cts_table_list_idx)
        cts_table_list_idx += 1

        #--------------------
        # run null simulation
        #--------------------
        (delta_Nt_null, delta_fit_null, NSness_null, Eness_null,
         Nt_list_null, mean_fit_list_null, mean_z_list_null) = run_sim(
                                                        'null',
                                                         linkage,
                                                         genicity,
                                                         n_its,
                                                         params,
                                                         output,
                                                         cts_table_list,
                                                         cts_table_list_idx)
        cts_table_list_idx += 1

        #-------------------------
        # null/non-null comparison
        #-------------------------

        assert (len(delta_Nt) == len(delta_fit) ==
                len(delta_Nt_null) == len(delta_fit_null) == n_its)

        linkage_col.extend([linkage]*n_its*2)
        genicity_col.extend([genicity]*n_its*2)
        nullness_col.extend(['non_null']*n_its)
        nullness_col.extend(['null']*n_its)
        delta_Nt_col.extend(delta_Nt)
        delta_Nt_col.extend(delta_Nt_null)
        delta_fit_col.extend(delta_fit)
        delta_fit_col.extend(delta_fit_null)
        NSness_col.extend(NSness)
        NSness_col.extend(NSness_null)
        Eness_col.extend(Eness)
        Eness_col.extend(Eness_null)

        # store mean phenotype data
        ts_data_dict['time_step'].extend([*range(len(mean_z_list))]*2)
        ts_data_dict['nullness'].extend((['non-null']*len(mean_z_list)) + (
                                        ['null']*len(mean_z_list_null)))
        ts_data_dict['genicity'].extend([genicity]*len(mean_z_list)*2)
        ts_data_dict['linkage'].extend([linkage]*len(mean_z_list)*2)
        ts_data_dict['Nt'].extend(Nt_list)
        ts_data_dict['Nt'].extend(Nt_list_null)
        ts_data_dict['mean_fit'].extend(mean_fit_list)
        ts_data_dict['mean_fit'].extend(mean_fit_list_null)
        ts_data_dict['mean_z'].extend(mean_z_list)
        ts_data_dict['mean_z'].extend(mean_z_list_null)
        assert len(np.unique([len(v) for v in ts_data_dict.values()])) == 1, (
            'all columns in ts_data_dict are not equal length!')

        assert (len(linkage_col) == len(genicity_col) == len(nullness_col) ==
                len(delta_Nt_col) == len(delta_fit_col))


print('\n\n\n%s\nNOW PREPPING DATA FOR PID %s\n\n%s\n' % ('<>'*40, str(pid), '<>'*40))
# gather delta Nt and delta fit data into a DataFrame
df = pd.DataFrame({'linkage': linkage_col,
                   'genicity': genicity_col,
                   'nullness': nullness_col,
                   'delta_Nt': delta_Nt_col,
                   'delta_fit': delta_fit_col,
                   'NSness': NSness_col,
                   'Eness': Eness_col,
                  })

df_dir = make_stat_df('dir', output)
df_ts_data = pd.DataFrame(ts_data_dict)

#/\/\/\/\/\/\/\
# finalize figs
#\/\/\/\/\/\/\/
if make_plots:

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
                    ax.set_title('|| genicity: %i' % genicity, size=collab_size)
                # plot a histogram of the unlinked and then linked data
                for nullness in ['non-null', 'null']:
                    #for neut_idx, neutrality in enumerate(['neut', 'nonneut']):
                    for neutrality in ['nonneut']:
                        data_dict = output[nullness][linkage][genicity][stat]
                        data = [v[neutrality] for v in data_dict.values()]
                        data = [val for sublist in data for val in sublist]
                        data = [np.nan if val is None else val for val in data]
                        if neutrality == 'neut':
                            try:
                                kde = scipy.stats.gaussian_kde(data)
                                xx = np.linspace(min(data), max(data), 1000)
                                kde_vals = kde(xx)
                                # NOTE: make a completely transparent hist,
                                #       to steal the bar heights
                                #       from it and use them to scale the kde!
                                vals, breaks, bars = ax.hist(data, bins=50,
                                                             alpha=0)
                                kde_plot_factor = max(vals)/max(kde_vals)
                                kde_plot_vals = [val * kde_plot_factor for val
                                                                in kde_vals]
                                ax.plot(xx, kde_plot_vals, alpha=0.5,
                                        label='%s: %s' % (nullness, neutrality),
                                        color=colors[nullness][neutrality])
                            except Exception as e:
                                print(('\n\nCOULD NOT PLOT KDE\n\nERROR '
                                      'THROWN:\n\t%s') % e)
                        else    :
                            ax.hist(data, bins = 50, alpha=0.5,
                                    label= '%s: %s' % (nullness, neutrality),
                                    color=colors[nullness][neutrality])
                        ax.set_xlabel(stat, size=8)
                        if row_idx==2 and col_idx==8:
                            ax.legend(prop={'size': 10},
                                      fontsize=5,
                                      bbox_to_anchor=(1.5, 0.2))
                        if col_idx in [0, 3, 6]:
                            ax.set_xticks(dir_tick_locs)
                            ax.set_xticklabels(dir_tick_labs)
                        ax.tick_params(labelsize=ticklabel_size)
                        ax.tick_params(axis='y', rotation=60)
                        # TODO standardize axes

    # show all the figures
    [fig.show() for fig in fig_time.values()]
    fig_hist.show()

    # set the subplot spacing
    plt.subplots_adjust(left=0.1,
                        bottom=0.1,
                        right=0.9,
                        top=0.9,
                        wspace=1.0,
                        hspace=0.9)
    try:
        rcParams['figure.figsize'] = 40, 12
    except Exception as e:
        pass
    for name, fig in fig_time.items():
        fig.savefig((output_path + '/fig_time_' + name + '_PID-%s' % pid +
                     '.png'), format='png', dpi=1000)
    fig_hist.savefig(output_path + '/fig_hist' + '_PID-%s' % pid + '.png',
                     format='png', dpi=1000)


#/\/\/\/\/\/\/\/\/\
# output dataframes
#\/\/\/\/\/\/\/\/\/

# save dfs to disk, including:

# high-level stats output
print('\n\n\n%s\nNOW SAVING MAIN DF FOR PID %s\n\n%s\n' % ('<>'*40, str(pid), '<>'*40))
df.to_csv(os.path.join(output_path, 'output_PID-%s.csv' % pid), index=False)

# df containing raw gene flow dir data from individual locus-chrom combos
print('\n\n\n%s\nNOW SAVING DIR DF FOR PID %s\n\n%s\n' % ('<>'*40, str(pid), '<>'*40))
df_dir.to_csv(os.path.join(output_path, 'output_PID-%s_DIR.csv' % pid),
              index=False)

# time-series df
print('\n\n\n%s\nNOW SAVING TS DF FOR PID %s\n\n%s\n' % ('<>'*40, str(pid), '<>'*40))
df_ts_data.to_csv(os.path.join(output_path, 'output_PID-%s_TS_DATA.csv' % pid),
                 index=False)
