import multiprocessing as mp
from run_batch_script import run_batch_script

if __name__ == '__main__':
    # how many CPUs?
    ncpu = mp.cpu_count()
    print("\n\nUSING %i CPUs\n" % ncpu)
    print('-'*80, '\n\n')

    # set method to 'spawn' instead of 'fork', to avoid deadlock
    # (see: https://pythonspeed.com/articles/python-multiprocessing/) 
    mp.set_start_method('spawn')

    # make the pool
    # NOTE: USING 5 JOBS WITH A POOL OF 10 WORKERS FOR HI-REDUNDANCY SCENARIOS (100 DATASETS TOTAL)
    #pool = mp.Pool(10)
    # NOTE: USING 3 JOBS WITH A POOL OF 17 WORKERS FOR LO-REDUNDANCY SCENARIOS (102 DATASETS TOTAL)
    pool = mp.Pool(4)

    # map a serial list of batch numbers into the pool
    pool.map_async(run_batch_script, [*range(ncpu)])

    # close the pool and join 'results' (although nothing to join)
    pool.close()
    pool.join()
