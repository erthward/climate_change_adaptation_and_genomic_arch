import re, os
import pandas as pd
import multiprocessing as mp
from estimate_von_mises_mixture_params import estimate_von_mises_mixture_params

if os.getcwd().split('/')[1] == 'home':
    datadir = '/home/deth/Desktop/CAL/research/projects/sim/ch2/climate_change_adaptation_and_genomic_arch/analysis'
else:
    with open(('/global/scratch/users/drewhart/ch2/climate_change_adaptation_'
               'and_genomic_arch/analysis/analysisdir.txt'), 'r') as f:
        datadir = f.read().strip()
    #datadir = '/global/scratch/users/drewhart/ch2/output/analysis'

# concatenate all fitted params files into a single file
fitted_params_files = [f for f in os.listdir(datadir) if re.search(
                                                    'DIR_FITTED_PARAMS', f)]
dfs = []
for f in fitted_params_files:
    df = pd.read_csv(os.path.join(datadir, f),
                     index_col=None,
                     header=0,
                     na_filter=False)
    dfs.append(df)
outdf = pd.concat(dfs, axis=0, ignore_index=True)
# write out concatenated params
outdf.sort_values(['linkage', 'genicity', 'nullness', 'it'])
outdf.to_csv(os.path.join(datadir, 'ch2_all_fitted_vM_params.csv'),
             index=False)
