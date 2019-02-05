from multiprocessing import Pool

def send_to_multiple_processes(func, args, workers):
    pool = Pool(processes=workers)
    res = pool.map(func, args)
    return res