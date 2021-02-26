import os
from multiprocessing import Process, Manager
import multiprocessing


def f(d, l):
    d[os.getpid()] = os.getpid()
    l.append(os.getpid())
    print(d)
    print(l)


if __name__ == '__main__':
    with Manager() as manager:
        d = manager.dict()
        l = manager.list(range(5))

        # pool = multiprocessing.Pool()
        # pool.map(f, d, l)
        process_list=[]
        for process_ in range(3):
            process_ = Process(target=f, args=(d, l))
            process_.start()
            process_list.append(process_)
        for p in process_list:
            p.join()
