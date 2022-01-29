import numpy as np
import xarray as xr
import os, re, sys

# CLI arg determine whether or not to run asserts
# (run without asserts to determine PIDs to drop;
#  run with asserts to make sure clean prior to running analysis)
run_asserts = sys.argv[1]
assert run_asserts.lower() in ['t', 'f']
if run_asserts.lower() == 't':
    run_asserts = True
else:
    run_asserts = False

# regex patts
dir_patt = 'GNX_mod-%s_L%s_G%i_its0_randID\d+PID-\d{6}'
timestep_patt = '(?<=t-)\d+(?=_spp)'
summary_csv_patt = 'output_PID-%s\w{0,8}\.csv'
pid_patt = '(?<=PID-)\d{6}'

# ouput directory
if os.getcwd().split('/')[1] == 'home':
    datadir = '/home/deth/Desktop/tmp_ch2_stats_tests_dev/'
else:
    with open(('/global/scratch/users/drewhart/ch2/climate_change_adaptation_'
               'and_genomic_arch/analysis/outputdir.txt'), 'r') as f:
        datadir = f.read().strip()

# store iteration counts
cts = xr.DataArray(np.zeros((2, 3, 3)),
                   {'nullness': ['null', 'non-null'],
                   'linkage': ['independent', 'weak', 'strong'],
                    'genicity': [4, 20, 100]
                   })

# store all the distinct timesteps for which data has been generated
all_timesteps = set()

# get all unique PIDs in dir
all_pids_in_dir = set([re.search(pid_patt,
            f).group() for f in os.listdir(datadir) if re.search(pid_patt, f)])
# track which PIDs have complete datasets
pids_complete = {pid: 0 for pid in all_pids_in_dir}

# for each PID
for pid in all_pids_in_dir:
    # for each nullness and scenario
    for nullness in cts.indexes['nullness']:
        for linkage in cts.indexes['linkage']:
            for genicity in cts.indexes['genicity']:
                # get all model-generated dirs
                curr_dir_patt = dir_patt % (nullness, linkage, genicity)
                dirs = [f for f in os.listdir(datadir) if re.search(curr_dir_patt,
                                                                    f)]

                # check that the dir has a complete dataset
                for d in dirs:
                    d = os.path.join(datadir, d, 'it--1', 'spp-spp_0')
                    files = os.listdir(d)
                    # get the distinct timesteps
                    timesteps_this_dir = set([int(re.search(timestep_patt,
                                                       f).group()) for f in files])
                    if run_asserts:
                        assert len(timesteps_this_dir) == 3
                    # add to overall set
                    for ts in timesteps_this_dir:
                        all_timesteps.add(ts)

                    # add 1 to this scenario-nullness combination's count
                    if len(timesteps_this_dir) == 3:
                        cts.loc[dict(nullness=nullness,
                                 linkage=linkage,
                                 genicity=genicity)] += 1

    # check there's exactly one set of 4 CSV files for this PID
    summary_csvs = [f for f in os.listdir(datadir) if
                    re.search(summary_csv_patt % pid, f)]
    if run_asserts:
        assert len(summary_csvs) == 4, ('Did not find exactly 4 summary '
                                    'CSVs for PID %s\n\nCSVs: '
                                   '%s') % (pid, ', '.join(summary_csvs))
    if len(summary_csvs) == 4:
        pids_complete[pid] += 1


# make sure only 3 timesteps covered across all files
if run_asserts:
    assert len(all_timesteps)==3, ('data files found for more than 3 timesteps!\n\n'
                                'timesteps: %s\n\n') % (', '.join([*all_timesteps]))

    # make sure no summary CSV files found for PIDS not covered by data dirs
    assert (len(np.unique([*pids_complete.values()])) == 1
        and np.unique([*pids_complete.values()]) == np.array(18)), ("Some PIDs do "
        "not have a full complement of 18 complete data dirs. Identify, cull "
        " then rerun!\n\nPID complete dir counts\n%s\n\n") % (str(pids_complete))

    # make sure that the same number of iterations is available for all scenarios
    assert np.all(cts.values == np.max(cts.values)), ("Not equal coverage across "
                                                      "all combinations of nullness"
                                                      " and scenario!\n\ncounts:\n"
                                "%s\n\n") % str(cts.values)

    print(('\n\n\tAll clean! %i complete datasets '
           'available.\n\n') % np.unique(cts.values)[0])
else:
    print('\n\ntimesteps with data files:/n%s/n/n' % str(all_timesteps))
    print('\n\ndataset counts:\n%s\n\n' % str(cts))
    print('\n\nPID completeness:\n%s\n\n' % str(pids_complete))
