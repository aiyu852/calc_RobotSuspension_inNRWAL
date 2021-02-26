import multiprocessing as mp
import itertools
from functools import partial

def job(z, r, item):
    (x, y) = item
    return x * y + z + r
    
def multicore(z):
    x_y = list(itertools.product(range(10), range(10)))
    r = 2
    func = partial(job, z, r)
    pool = mp.Pool() 
    res = pool.map(func, x_y)
    return res

if __name__ == '__main__':
    res = multicore(1)
    print(res)