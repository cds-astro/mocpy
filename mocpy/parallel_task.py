from multiprocessing import Pool

def send_to_multiple_processes(func, args, workers):
    pool = Pool(workers)
    res = pool.map(func, args)
    pool.close()
    pool.join()
    return res