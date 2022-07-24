import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import statsmodels.api as sm
from copy import deepcopy
import vcf
import sys
import re
import os


# get arg determining whether to plot null or non-null sims
nullness = sys.argv[1].lower()
assert nullness in ['null', 'non-null']

# get arg determining whether or not sims are redundant
redundant = sys.argv[2].lower()
assert redundant in ['t', 'f']
if redundant == 't':
    redundant = True
else:
    redundant = False

# plot params
suptitle_fontsize = 50
title_fontsize = 40
axislab_fontsize = 18
ticklab_fontsize = 12
cbarlab_fontsize = 14
cbarticklab_fontsize=14
fig_width = 10
fig_height = 12
dpi = 400
n_ticklabels = 5
linewidth = 1
linealpha = 0.8
marker_sizes = {4: 500,
                8: 500,
                20: 40,
                40: 40,
                100: 3,
                200: 3,
               }
subplots_adj_left=0.12
subplots_adj_bottom=0.08
subplots_adj_right=0.90
subplots_adj_top=0.95
subplots_adj_wspace=0.14
subplots_adj_hspace=0.05

# data directory
if os.getcwd().split('/')[1] == 'home':
    datadir = '/media/deth/SLAB/ch2/dense_ts'
else:
    datadir = 'global/scratch/users/drewhart/ch2/output/output_dense_ts'

# save results?
save_res = True

# x-threshold of right side of landscape tracked for the analysis
l_edge_of_r_side = 40

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
change_start_t = change_T-1
change_end_t = T-1
change_half_t = int((change_start_t + change_end_t)/2)

# lists of all possible linkage and genicity values
linkages = ['independent', 'weak', 'strong']
# NOTE: will be 8, 40, 200 if redundant scenarios, otherwise 4, 20, 100
genicities =  [*np.array([4, 20, 100]) * (1 + 1*redundant)]

# get the horizontal values used in the non-shifting landscape layer
# and in the before- and after-shifting shifting layer,
# then set as dict of x and y values for plotting horizontal expectation lines
stable = np.linspace(0, 1, 50)[[0,-1]]
b4 = np.linspace(0, 1, 50)[[0,-1]]
af = np.linspace(0.5, 1, 50)[[0,-1]]

# get candidate filenames for change-start-time-step files
dirname_patt_template = 'mod-%s_L%s_G%i_its0_'
filename_patt_template = ('mod-%s_L%s_G%i_its0_randID\d{7}PID-'
                          '\d{5,6}_it--1_t-\d{3,4}_spp-spp_0.csv')

# number of time steps for which we should have files
timesteps = [*range(change_start_t, change_end_change_end_t, 2)]
n_timesteps = len(timesteps)

# alpha factor multiplied by 1/genicity to derive effect sizes in sims
alpha_factor=2

# whether to set vmin and vmax
set_vmin_vmax = False


def calc_allele_freqs_at_landscape_ends(vcf_filepath, csv_filepath,
                                        l_edge_of_r_side=40):
    reader = vcf.Reader(filename=vcf_filepath)
    df = pd.read_csv(csv_filepath)
    # NOTE: setting individs' idx as df index and then querying if x >=
    # l_edge_of_r_side is >2x faster than getting valid individuals
    # and then finding if each sample's individual is in that list
    df = df.set_index('idx')
    # get freqs of 1 alleles within that side, for each locus
    freqs = []

    for rec in reader:
        ct_1_allele = 0
        ct_all_alleles = 0
        for samp in rec.samples:
            if df.loc[int(samp.sample), 'x'] >= l_edge_of_r_side:
                ct = len([n for n in samp.gt_nums.split('|') if n == '1']) 
                ct_1_allele += ct
                ct_all_alleles += 2
        freq = ct_1_allele / ct_all_alleles
        freqs.append(freq)
    # keep only every the 1-allele freqs of loci 0, 4, ... and the 0-allele
    # freqs of loci 2, 6, ..., because the even loci were set to influence
    # trait 0 (the trait responding to climate change), and the effect sizes
    # were set to +, -, +, -, ... (such that 1 alleles could contribute
    # to climate adaptation for loci 0, 4, ... and 0 alleles could contribute
    # to climate adaptation for loci 2, 6, ...)
    freqs = freqs[::2]
    freqs = [freq if n%2==0 else 1-freq for n, freq in enumerate(freqs)]
    return freqs


def calc_contributions_to_divergence(freqs, start_freqs, α):
    ds = [2*α*(freq - start_freq) for freq, start_freq in zip(freqs,
                                                              start_freqs)]
    return ds



def plot_allele_freq_changes(linkage, genicity, ax):

    assert linkage in ['independent', 'weak', 'strong']
    assert genicity in [2, 4, 8, 10, 20, 40, 50, 100, 200]

    dirname_patt = dirname_patt_template % (nullness, linkage, genicity)
    filename_patt = filename_patt_template % (nullness, linkage, genicity)

    filenames = {}
    for dirname in os.listdir(datadir):
        if (dirname.startswith('GNX_') and
            re.search(dirname_patt, dirname)):
            candidate_filenames = [fn for fn in os.listdir(os.path.join(
                datadir, dirname, 'it--1', 'spp-spp_0')) if re.search(filename_patt, fn)]
            # only add this directory and its files to the analysis if I got
            # all time steps, otherwise print warning
            if len(candidate_filenames) == n_timesteps:
                # sort files by timestep
                candidate_filenames = [*np.sort(candidate_filenames)]
                filenames[dirname] = candidate_filenames
            else:
                print(('\n\nWARNING: following directory did not contain'
                       '%i valid files, 1 for each of the right timesteps:\n\n\t'
                       '%s\n\n') % (n_timesteps, dirname))

    # keep just a single complete set of files
    dirname = [*filenames.keys()][0]
    filenames = filenames[dirname]

    # get the full dirname for these files
    filesdir = os.path.join(datadir,
                           dirname,
                           'it--1',
                           'spp-spp_0')

    # calculate the α value
    α = 1/genicity*alpha_factor

    # loop through the time steps to be analyzed
    ds_arr = np.ones((genicity, len(timesteps)))*np.nan
    for j, ts in enumerate(timesteps):

        print('\n\ttime step %i\n\n' % ts)

        # get vcf and csv filenames
        csv_filename = [f for f in filenames if re.search(
                                                't-%i_spp-spp_0.csv' % ts,f)]
        assert len(csv_filename) == 1
        csv_filename = csv_filename[0]
        vcf_filename = os.path.splitext(csv_filename)[0]+'.vcf'

        # concat with paths
        csv_filepath = os.path.join(filesdir, csv_filename)
        vcf_filepath = os.path.join(filesdir, vcf_filename)

        # get the allele freqs on the right side of the landscape
        freqs = calc_allele_freqs_at_landscape_ends(vcf_filepath, csv_filepath,
                                                    l_edge_of_r_side=40)

        # if this is the first timestep, store these for comparison with all
        # later timesteps
        if ts == change_start_t:
            start_freqs = deepcopy(freqs)
            ds_arr[:,j] = 0
        # else, compare to first timestep and store 'divergence' contribution
        # values
        else:
            ds = calc_contributions_to_divergence(freqs, start_freqs, α)
            ds_arr[:,j] = ds

    # display the ds array on the given axes
    if set_vmin_vmax:
        vmin, vmax = (-2*α, 2*α)
    else:
        vmin, vmax = (None, None)
    im = ax.imshow(ds_arr,
              cmap='coolwarm',
              vmin=vmin,
              vmax=vmax)
    ax.set_aspect('auto')
    ax.invert_yaxis()

    # save array to file
    arr_filename = 'allele_freqs_L%s_G%s_%s%s.txt' % (linkage,
                                                      genicity,
                                                      nullness,
                                                      ('_redundant' * redundant))
    np.savetxt(os.path.join(datadir, arr_filename), ds_arr)

    # return image, to be used for colorbar if needed
    return im


fig = plt.figure(figsize=(fig_height,fig_width))
gs = fig.add_gridspec(nrows=3, ncols=3, width_ratios=[1,1,1.2])
for i, linkage in enumerate(linkages):
    for j, genicity in enumerate(genicities):
        ax = fig.add_subplot(gs[i,j])
        print('\n\nProcessing file for linkage=%s, genicity=%i...\n\n' % (
                                                            linkage, genicity))
        im = plot_allele_freq_changes(linkage, genicity, ax)
        # add colorbar on last column
        if j == 2:
            α = 1/genicity*alpha_factor
            cbar_tickpos = [-2*α, 0, 2*α]
            cbar_ticklabs = ['$-2\\alpha$', '$0$', '$2\\alpha$']
            cbar = plt.colorbar(im, ax=ax, ticks=cbar_tickpos)
            cbar.set_ticklabels(cbar_ticklabs)
            cbar.set_label('contribution to\nadaptive divergence',
                           size=cbarlab_fontsize)
        #if i == 0:
        #    ax.set_title('%i' % genicity, fontdict={'fontsize': title_fontsize})
        if i < 2:
            ax.set_xticks(())
            ax.set_xticklabels(())
        else:
            ax.set_xticks(np.linspace(-0.5, len(timesteps)-0.5, 5),
                          [str(int(n)) for n in np.linspace(timesteps[0]+1,
                                                       timesteps[-1]+1, 5)])
            ax.set_xlabel('time step', fontdict={'fontsize': axislab_fontsize})
        if j > 0:
            ax.set_yticks(())
            ax.set_yticklabels(())
        else:
            #ax.set_ylabel('%s' % linkage, fontdict={'fontsize': title_fontsize})
            ax.set_ylabel('locus position', fontdict={'fontsize':
                                                     axislab_fontsize})
            ax.set_yticks(np.linspace(-0.5, genicity-0.5, 5),
                          [str(int(n)) for n in np.linspace(0, genicity*2, 5)])
        ax.set_xlim(-0.5, len(timesteps)-0.5)
        ax.set_ylim(-0.5, genicity-0.5)
        ax.tick_params(labelsize=ticklab_fontsize)
        ax.tick_params(axis=u'both', which=u'both',length=0)

fig.subplots_adjust(left=subplots_adj_left,
                    bottom=subplots_adj_bottom,
                    right=subplots_adj_right,
                    top=subplots_adj_top,
                    wspace=subplots_adj_wspace,
                    hspace=subplots_adj_hspace)
fig.show()

if save_res:
    fig.savefig('allele_freq_changes_%s%s.png' % (nullness,
                                                  '_redundant' * redundant),
                dpi=500)
