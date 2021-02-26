from multiprocessing import Process, Value, Lock


def f(s, x, lock):
    for x_ in x:
        lock.acquire()
        if x_ > s.value:
            s.value = x_
        lock.release()

if __name__ == '__main__':
    xs = range(201)

    sum = Value('d', 0)
    lock = Lock()

    num = 12
    processes = [Process(target=f, args=(sum, xs[i*len(xs)//num:(i+1)*len(xs)//num+1],  lock))
                 for i in range(num)]
    for p in processes:
        p.start()
    for p in processes:
        p.join()
    print(sum.value)
