import multiprocessing
import time


def f(x):
    return x**12


if __name__ == '__main__':

    start = time.time()
    cores = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(processes=cores)
    xs = range(200)

    # method 2: imap
    for y in pool.imap(f, xs):
        print(y)            # 0, 1, 4, 9, 16, respectively

    end = time.time()
    print(end-start)
