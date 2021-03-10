import numpy as np
from multiprocessing import Process, Manager
import multiprocessing
import time
import itertools
from functools import partial

# 计算的坐标遍历细分值
sub_num = 3
# 遍历的坐标空间
Lowlim_x = -40.5
Upplim_x = 33.5
Lowlim_y = 39
Upplim_y = 45
# 轮胎中心的坐标值
xy0 = (94.104, -11.439)

# 弹簧刚度
G = 70000

'''
find_avaliable_point
设置遍历的空间

Interfere_WheelandSpring
判断弹簧是否和轮胎干涉

calc_counterforce
计算地面支撑力

find_spring
根据传入弹簧参数调用calc_counterforce

process
定义了单进程的循环计算过程
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
        # 判断弹簧拉伸的长度是否过长
        if (L_spring_now-L0_spring)/L0_spring > 1:
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
    if Interfere_WheelandSpring(xy1, xy2, L_spring):
        return None, 'The Wheel and spring is interfered!'
    # 计算弹簧固定点与转轴的间距H和弹簧铰接点与转轴的间距r和H与r的夹角alpha
    r = np.sqrt((xy0[0]-xy2[0])**2+(xy0[1]-xy2[1])**2)
    H = np.sqrt((xy0[0]-xy1[0])**2+(xy0[1]-xy1[1])**2)
    if r > H:
        return None, 'The spring is opposite direction'
    alpha = np.arccos((r**2+H**2-L_spring**2)/(2*r*H))
    # 计算变化的地面反力
    return calc_counterforce(xy0, k, r, H, alpha, L0_spring)


def process(items, find_values):
    find_F_N = []
    find_F_N_max = 1
    find_F_N_min = 0
    find_n = 0
    find_D = 0
    find_d = 0
    find_xy1 = (0, 0)
    find_xy2 = (0, 0)
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
                        if F_N.size < 26:
                            continue
                        else:
                            if np.min(F_N) == 0:
                                continue
                            if np.max(F_N) > 15:
                                continue
                            if (((np.max(F_N)-np.min(F_N))/np.max(F_N) < (find_F_N_max-find_F_N_min)/find_F_N_max) and np.max(F_N) > 11 and np.max(F_N) < 13):
                                find_F_N = F_N
                                find_F_N_max = np.max(F_N)
                                find_F_N_min = np.min(F_N)
                                find_n = item[0]
                                find_D = item[1]
                                find_d = item[2]
                                find_xy1 = xy1
                                find_xy2 = xy2
                    except:
                        break
            except:
                break
    f = [find_F_N, find_F_N_max, find_F_N_min,
         find_n, find_D, find_d, find_xy1, find_xy2]
    find_values.append(f)


if __name__ == '__main__':
    # 最终数值

    find_values = Manager().list()

    items = list(itertools.product(range(30, 55),
                                   range(3, 7), [0.3, 0.5, 0.6, 0.8, 1, 1.2]))
    # items = list(itertools.product(range(47, 48),
    #                                range(4, 5), [0.6]))

    num_cpus = multiprocessing.cpu_count()
    print(num_cpus)

    start = time.time()
    processes = [Process(target=process, args=(items[i*len(items)//num_cpus:(i+1)*len(items)//num_cpus+1], find_values))
                 for i in range(num_cpus//2*3)]
    for p in processes:
        p.start()
    for p in processes:
        p.join()
    end = time.time()

    print(end-start)
    # N = []
    # N_max = 1
    # N_min = 0
    # N_n = 0
    # N_D = 0
    # N_d = 0
    # N_xy1 = (0, 0)
    # N_xy2 = (0, 0)
    # for F_N, F_N_max, F_N_min, n, D, d, xy1, xy2 in find_values:
    #     if ((F_N_max-F_N_min)/F_N_max < (N_max-N_min)/N_max):
    #         N = F_N
    #         N_max = F_N_max
    #         N_min = F_N_min
    #         N_n = n
    #         N_D = D
    #         N_d = d
    #         N_xy1 = xy1
    #         N_xy2 = xy2

    # print("D,d,n={},{},{},xy1={},xy2={}, F_max={},对地压力损失率={}".format(
    #     N_D, N_d, N_n, N_xy1, N_xy2, N_max, (N_max-N_min)/N_max))
    # print(N)
