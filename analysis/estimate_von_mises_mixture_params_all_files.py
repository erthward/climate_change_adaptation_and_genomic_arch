import re, os
import pandas as pd
import multiprocessing as mp
from estimate_von_mises_mixture_params import estimate_von_mises_mixture_params

datadir = '../../output/output/'
file_patt = 'DIR_short\.csv'

if __name__ == '__main__':
    # get all valid files to be processed
    files = [os.path.join(datadir, f) for f in os.listdir(datadir) if re.search(
                                                        'DIR_short\.csv', f)]
    # how many CPUs?
    ncpu = mp.cpu_count()
    print("\n\nUSING %i CPUs\n" % ncpu)
    print('-'*80, '\n\n')
    # set method to 'spawn' instead of 'fork', to avoid deadlock
    # (see: https://pythonspeed.com/articles/python-multiprocessing/) 
    mp.set_start_method('spawn')
    # make the pool
    pool = mp.Pool(int(ncpu/2))
    #pool = mp.Pool(ncpu)
    #pool = mp.Pool(1)
    # map a serial list of batch numbers into the pool
    pool.map_async(estimate_von_mises_mixture_params, files)
    # close the pool and join 'results' (although nothing to join)
    pool.close()
    pool.join()



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
