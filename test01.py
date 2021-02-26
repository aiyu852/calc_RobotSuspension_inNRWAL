from multiprocessing import Process, Value, Lock
from time import sleep


def count(x, lock):
    for i in range(20):
        sleep(0.01)
        with lock:
            x.value += 1


def main():
    counter = Value('i', 0)
    lock = Lock()
    processes = [Process(target=count, args=(counter, lock))
                 for i in range(10)]
    for p in processes:
        p.start()
    for p in processes:
        p.join()
    print(counter.value)


if __name__ == '__main__':
    for i in range(10):
        main()
