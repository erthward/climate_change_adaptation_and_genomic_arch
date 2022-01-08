import re, os
import pandas as pd
import multiprocessing as mp
from estimate_von_mises_mixture_params import estimate_von_mises_mixture_params

if os.getcwd().split('/')[1] == 'home':
    datadir = '/home/deth/Desktop/CAL/research/projects/sim/ch2/climate_change_adaptation_and_genomic_arch/analysis'
else:
    with open(('/global/scratch/users/drewhart/ch2/climate_change_adaptation_'
               'and_genomic_arch/analysis/outputdir.txt'), 'r') as f:
        datadir = f.read().strip()
    #datadir = '/global/scratch/users/drewhart/ch2/output/output'

if __name__ == '__main__':
    # get all valid files to be processed
    files = [os.path.join(datadir, f) for f in os.listdir(datadir) if re.search(
                                                        'DIR\.csv', f)]
    # how many CPUs?
    # NOTE: divide in half to avoid annoying OOM errors (for now...)
    ncpu = mp.cpu_count()/2
    # reduce to just the number of files, if < ncpu
    if len(files) < ncpu:
        ncpu = len(files)
    print("\n\nUSING %i CPUs\n" % ncpu)
    print('-'*80, '\n\n')
    # set method to 'spawn' instead of 'fork', to avoid deadlock
    # (see: https://pythonspeed.com/articles/python-multiprocessing/) 
    mp.set_start_method('spawn')
    # make the pool
    pool = mp.Pool(int(ncpu))
    #pool = mp.Pool(ncpu)
    #pool = mp.Pool(1)
    # map a serial list of batch numbers into the pool
    pool.map_async(estimate_von_mises_mixture_params, files)
    # close the pool and join 'results' (although nothing to join)
    pool.close()
    pool.join()
