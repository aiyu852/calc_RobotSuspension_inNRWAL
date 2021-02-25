import numpy as np
# import matplotlib.pyplot as plt

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

# 最终数值
find_F_N_max = 0
find_F_N_var = 0
find_n = 0
find_D = 0
find_d = 0
find_xy1 = (0, 0)
find_xy2 = (0, 0)

'''
find_avaliable_point
设置遍历的空间

Interfere_WheelandSpring
判断弹簧是否和轮胎干涉

calc_counterforce
计算地面支撑力
'''


def find_avaliable_point(Lx, Ux, Ly, Uy, sn):
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


def calc_counterforce(xy0, r, H, alpha, theta):
    F_N = []
    F_spring = []
    for i, height_WheeltoChassis in enumerate(np.arange(10, 36)):
        alpha_delta = np.arcsin(
            (height_WheeltoChassis-10+xy0[1])/s)-np.arcsin(xy0[1]/s)
        xy2_now = (xy2[0]*np.cos(alpha_delta)-xy2[1]*np.sin(alpha_delta),
                   xy2[0]*np.sin(alpha_delta)+xy2[1]*np.cos(alpha_delta))
        alpha_now = alpha - alpha_delta
        L_spring_now = np.sqrt(r**2+H**2-2*r*H*np.cos(alpha_now))
        if L_spring_now < L0_spring:
            print("弹簧失去弹力")
            break
        elif Interfere_WheelandSpring(xy1, xy2_now, L_spring_now):
            print("弹簧和轮胎干涉")
            break
        else:
            # 计算弹簧力F_spring，弹簧力臂arm_F_spring,判断弹簧力是否单调
            F_spring.append(k*(L_spring_now-L0_spring))
            if i > 0 and F_spring[i-1] < F_spring[i]:
                print("弹簧力不单调")
                break
            beta = np.arccos(
                (L_spring_now**2+H**2-r**2)/(2*L_spring_now*H))
            arm_F_spring = H*np.sin(beta)

            # 得到平衡扭矩T
            T = arm_F_spring * F_spring

            # 计算地面支反力力臂arm_F_N， 地面支反力F_N
            arm_F_N = s * np.cos(theta)
            F_N.append(T/arm_F_N)
            return F_N


def find_spring(xy0, xy1, xy2, L0_spring):
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
    return calc_counterforce(xy0, r, H, alpha, theta)


def fanded_evaluation(n, D, d, xy1, xy2, F_N, find_xy1, find_xy2,
                      find_F_N_max, find_F_N_var, find_n, find_D, find_d):
    find_F_N_max = np.nax(F_N)
    find_F_N_var = np.var(F_N)
    find_xy1 = xy1
    find_xy2 = xy2
    find_n = n
    find_D = D
    find_d = d


num_points = 0
num_springs = 0
while 1:
    try:
        # xy1 是弹簧的固定点
        xy1 = next(find_avaliable_point(
            Lowlim_x, Upplim_x, Lowlim_y, Upplim_y, sub_num))
        while 1:
            try:
                # xy2 是弹簧的铰接点
                xy2 = next(find_avaliable_point(
                    Lowlim_x, Upplim_x, Lowlim_y, Upplim_y, sub_num))
                num_points += 1
                num_springs = 0
                for n in np.arange(40, 51, 2):
                    for D in range(4, 7):
                        for d in [0.3, 0.5, 0.6, 0.8, 1]:
                            num_springs += 1
                            C = D/d
                            k = G*d/(8*(C**3)*n)
                            # 弹簧原长
                            L0_spring = (n+1)*d

                            print("在计算{}个点,第{}种弹簧".format(
                                num_points, num_springs))

                            if find_spring(xy0, xy1, xy2, L0_spring) is None:
                                continue
                            else:
                                F_N = np.array(find_spring(
                                    xy0, xy1, xy2, L0_spring))
                                if np.min(F_N) == 0:
                                    continue
                                elif np.max(F_N) > 15:
                                    continue
                                elif np.max(F_N) > find_F_N_max or find_F_N_var > np.var(F_N):
                                    fanded_evaluation(n, D, d, xy1, xy2, F_N, find_xy1, find_xy2,
                                                      find_F_N_max, find_F_N_var, find_n, find_D, find_d)
            except:
                break
    except:
        break
