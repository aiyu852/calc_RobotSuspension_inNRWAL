import numpy as np
from multiprocessing import Process, Value, Lock
import time
import itertools
from functools import partial

# 计算的坐标遍历细分值
sub_num = 10
# 遍历的坐标空间
Lowlim_x = -40.7
Upplim_x = 52
Lowlim_y = 0
Upplim_y = 45.8
# 转轴中心的坐标值
xy0 = (94.104, 11.439)
# 轮子中心与转轴的距离s
s = np.sqrt(xy0[0]**2+xy0[1]**2)
# 轮中心与转轴连线和水平线的夹角θ
theta = np.arctan(xy0[1]/xy0[0])
# 弹簧刚度
G = 70000
'''
find_avaliable_point
设置遍历的空间

Interfere_WheelandSpring
判断弹簧是否和轮胎干涉

calc_counterforce
计算地面支撑力
'''


def find_avaliable_point(Lx, Ux, Ly, Uy, sn, D, d):
    for x in np.arange(Lx, Ux, sn):
        for y in np.arange(Ly, Uy, sn):
            # 坐标点不在车轮圆内，不在电机所在的圆内，不在两个圆切线下
            if (x**2+y**2 <= (39.5+(D+d)/2)**2) or ((x-48.1057)**2+(y-26.4438)**2 <= (17.5+(d+D)/2)**2) or (x > 0 and y <= 0.2627*x+31.9017+(D+d)/2/np.cos(np.arctan(0.2627))):
                continue
            yield x, y


def Interfere_WheelandSpring(xy1, xy2, length_spring):
    l1 = np.sqrt(xy1[0]**2+xy1[1]**2)
    l2 = np.sqrt(xy2[0]**2+xy2[1]**2)
    p = (l1+l2+length_spring)/2
    if np.sqrt(p*(p-l1)*(p-l2)*(p-length_spring))*2/length_spring <= 42:
        return 1
    else:
        return 0


def calc_counterforce(xy0, k, r, H, alpha, L0_spring):
    F_N = []
    F_spring = []
    for index, height_WheeltoChassis in enumerate(range(10, 36)):
        alpha_delta = np.arcsin(
            (height_WheeltoChassis-10+xy0[1])/s)-np.arcsin(xy0[1]/s)
        alpha_now = alpha - alpha_delta
        L_spring_now = np.sqrt(r**2+H**2-2*r*H*np.cos(alpha_now))
        # 判断弹簧是否失去弹力
        if L_spring_now < L0_spring:
            return None
        else:
            # 计算弹簧力F_spring，弹簧力臂arm_F_spring,判断弹簧力是否单调
            F_spring.append(k*(L_spring_now-L0_spring))
            if index > 0 and F_spring[index-1] < F_spring[index]:
                return None
            beta = np.arccos(
                (L_spring_now**2+H**2-r**2)/(2*L_spring_now*H))
            arm_F_spring = H*np.sin(beta)

            # 得到平衡扭矩T
            T = arm_F_spring * F_spring[index]

            # 计算地面支反力力臂arm_F_N， 地面支反力F_N
            arm_F_N = s * np.cos(theta)
            F_N.append(T/arm_F_N)
            # print(index, F_spring[index], F_N[index])
    return F_N


def find_spring(k, xy0, xy1, xy2, L0_spring):
    # 计算弹力最大状态下的弹簧长度（认为遍历的点是轮胎收上去的状态）
    L_spring = np.sqrt((xy1[0]-xy2[0])**2+(xy1[1]-xy2[1])**2)
    # 弹簧太短就提前跳出
    if L_spring < L0_spring:
        return None, 'The spring is too short'
    # 计算弹簧固定点与转轴的间距H和弹簧铰接点与转轴的间距r和H与r的夹角alpha
    r = np.sqrt((xy0[0]-xy2[0])**2+(xy0[1]-xy2[1])**2)
    H = np.sqrt((xy0[0]-xy1[0])**2+(xy0[1]-xy1[1])**2)
    alpha = np.arccos((r**2+H**2-L_spring**2)/(2*r*H))
    # 计算变化的地面反力
    return calc_counterforce(xy0, k, r, H, alpha, L0_spring)


def process(find_F_N_max, find_F_N_var, items, lock):
    for item in items:
        C = item[1]/item[2]
        k = G*item[2]/(8*(C**3)*item[0])
        # 弹簧原长
        L0_spring = (item[0]+1)*item[2]
        xy1_source = find_avaliable_point(
            Lowlim_x, Upplim_x, Lowlim_y, Upplim_y, sub_num, item[1], item[2])
        while 1:
            try:
                # xy1 是弹簧的固定点
                xy1 = next(xy1_source)
                xy2_source = find_avaliable_point(
                    Lowlim_x, Upplim_x, Lowlim_y, Upplim_y, sub_num, item[1], item[2])
                while 1:
                    try:
                        # xy2 是弹簧的铰接点
                        xy2 = next(xy2_source)

                        # print("计算第{}种弹簧,第{}个点".format(
                        #     num_springs, num_points))
                        if find_spring(k, xy0, xy1, xy2, L0_spring) is None:
                            continue
                        F_N = np.array(find_spring(
                            k, xy0, xy1, xy2, L0_spring))
                        lock.acquire()
                        if F_N.size < 26:
                            continue
                        else:
                            if np.min(F_N) == 0:
                                continue
                            if np.max(F_N) > 15:
                                continue
                            with lock:
                                if (np.max(F_N) > find_F_N_max.value and find_F_N_max.value < 12) or (np.var(F_N) < find_F_N_var.value and find_F_N_max.value > 12):
                                    # find_F_N = F_N
                                    find_F_N_max.value = np.max(F_N)
                                    find_F_N_var.value = np.var(F_N)
                                    # find_n = item[0]
                                    # find_D = item[1]
                                    # find_d = item[2]
                                    # find_xy1 = xy1
                                    # find_xy2 = xy2
                    except:
                        break
            except:
                break


if __name__ == '__main__':
    # 最终数值

    find_F_N_max = Value('d', 0)
    find_F_N_var = Value('d', 0)
    # find_F_N_max = 0
    # find_F_N_var = 0
    # find_n = 0
    # find_D = 0
    # find_d = 0
    # find_xy1 = (0, 0)
    # find_xy2 = (0, 0)
    lock = Lock()

    items = list(itertools.product(range(30, 55),
                                   range(3, 7), [0.3, 0.5, 0.6, 0.8, 1, 1.2]))

    num_cpus = 20

    start = time.time()
    processes = [Process(target=process, args=(find_F_N_max, find_F_N_var, items[i*len(items)//num_cpus:(i+1)*len(items)//num_cpus+1], lock))
                 for i in range(num_cpus)]
    for p in processes:
        p.start()
    for p in processes:
        p.join()
    end = time.time()
    # print("D,d,n={},{},{},xy1={},xy2={}, F_max={},F_var={}".format(
    #     find_D, find_d, find_n, find_xy1, find_xy1, find_F_N_max, find_F_N_var))
    # print(find_F_N)

    print(end-start)

    print(find_F_N_max.value, find_F_N_var.value)
