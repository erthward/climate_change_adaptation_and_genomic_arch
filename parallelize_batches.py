import multiprocessing as mp
from run_one_batch import run_one_batch

if __name__ == '__main__':
    # how many CPUs?
    ncpu = mp.cpu_count()
    print("USING %i CPUs\n" % ncpu)
    print('-'*80, '\n\n')

    # set method to 'spawn' instead of 'fork', to avoid deadlock
    # (see: https://pythonspeed.com/articles/python-multiprocessing/) 
    mp.set_start_method('spawn')

    # make the pool
    pool = mp.Pool(ncpu)

    # map a serial list of batch numbers into the pool
    pool.map_async(run_one_batch, [*range(ncpu)])

    # close the pool and join 'results' (although nothing to join)
    pool.close()
    pool.join()
