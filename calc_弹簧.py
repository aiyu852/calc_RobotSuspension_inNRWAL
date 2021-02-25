import numpy as np
import matplotlib.pyplot as plt

# 计算的坐标遍历细分值
sub_num = 10


# 遍历的坐标空间
Lowlim_x = -40.7
Upplim_x = 52
Lowlim_y = 0
Upplim_y = 45.8


# 转轴中心的坐标值
x0 = 94.104
y0 = 11.439

# 轮子中心与转轴的距离s
s = np.sqrt(x0**2+y0**2)
# 轮中心与转轴连线和水平线的夹角θ
theta = -6.9/180*np.pi

# 得到和轮外形不干涉的遍历点


def find_avaliable_point(Lx, Ux, Ly, Uy, sn):
    for x in np.arange(Lx, Ux, sn):
        for y in np.arange(Ly, Uy, sn):
            if (x**2+y**2 <= (39.5+(D+d)/2)**2) or ((x-48.1057)**2+(y-26.4438)**2 <= (17.5+(d+D)/2)**2) or (x > 0 and y <= 0.2627*x+31.9017+(D+d)/2/np.cos(np.arctan(0.2627))):
                continue
            yield x, y


# 最终数值
find_F_N_mean = 0
find_F_N_var = 0
find_n = 0
find_D = 0
find_d = 0
find_x1, find_y1 = (0, 0)
find_x2, find_y2 = (0, 0)

# 弹簧刚度
G = 70000
# plt.figure()
for n in np.arange(40, 51, 2):
    for D in range(4, 7):
        for d in [0.3, 0.5, 0.6, 0.8, 1]:

            # for d in [1]:
            if d == 1:
                R_m = 1450
            elif d == 0.8:
                R_m = 1550
            elif d == 0.3:
                R_m = 1650
            else:
                R_m = 1600
            C = D/d
            k = G*d/(8*(C**3)*n)
            # 弹簧原长
            L0_spring = (n+1)*d
            F_N_max = 0

            x1_y1 = find_avaliable_point(
                Lowlim_x, Upplim_x, Lowlim_y, Upplim_y, sub_num)
            num = 0
            while 1:
                try:
                    # x1，y1 是弹簧的固定点
                    x1, y1 = next(x1_y1)
                    x2_y2 = find_avaliable_point(
                        Lowlim_x, Upplim_x, Lowlim_y, Upplim_y, sub_num)
                    while 1:
                        num += 1
                        try:
                            # x2，y2 是弹簧的铰接点
                            x2, y2 = next(x2_y2)
                            print("n={},D={},d={},已经计算{}个点,({},{}),({},{})".format(
                                n, D, d, num, x1, y1, x2, y2))

                            # 计算当前状态下的弹簧长度
                            L_spring = np.sqrt((x1-x2)**2+(y1-y2)**2)

                            #
                            #
                            # 提前跳出当前循环
                            if L_spring < 35:
                                continue
                            # 弹簧位置不能和轮胎干涉
                            l1 = np.sqrt(x1**2+y1**2)
                            l2 = np.sqrt(x2**2+y2**2)
                            p = (l1+l2+L_spring)/2
                            if np.sqrt(p*(p-l1)*(p-l2)*(p-L_spring))*2/L_spring <= 42:
                                continue
                            #
                            #

                            # 计算弹簧铰接点与转轴的距离r
                            r = np.sqrt((x0-x2)**2+(y0-y2)**2)
                            # 计算弹簧固定点与转轴的间距H
                            H = np.sqrt((x0-x1)**2+(y0-y1)**2)

                            # 计算角度变化量
                            # 以地面距离车底壳距离为循环量
                            F_N = []
                            list_F_spring = []
                            i = 0
                            for height_wheeltochassis in np.arange(10, 36):

                                # 计算r和H的初始角度值
                                alpha = np.arccos(
                                    (r**2+H**2-L_spring**2)/(2*r*H))
                                alpha_delta = np.arcsin(
                                    (height_wheeltochassis-10+11.4)/s)-np.arcsin(11.4/s)
                                alpha -= alpha_delta
                                # 计算弹簧弹力
                                L_spring_now = np.sqrt(
                                    r**2+H**2-2*r*H*np.cos(alpha))
                                F_spring = k*(L_spring_now-L0_spring)
                                if F_spring < 0:
                                    F_spring = 0.0
                                    print("弹簧弹力为0")
                                    break
                                list_F_spring.append(F_spring)

                                # 计算弹力力臂
                                beta = np.arccos(
                                    (L_spring_now**2+H**2-r**2)/(2*L_spring_now*H))
                                arm_F_spring = H*np.sin(beta)
                                # 弹簧对转轴施加的扭矩
                                T = arm_F_spring * F_spring
                                # 对地压力反力力臂
                                arm_F_N = s * np.cos(theta)
                                # 对地压力
                                F_N.append(T/arm_F_N)
                                if i > 0:
                                    if F_N[i-1] < F_N[i]:
                                        print("不单调")
                                        break
                                i += 1
                                # print(height_wheeltochassis, alpha*180/np.pi,
                                #       L_spring_now, F_spring, T, T/arm_F_N)
                            if i < 26:
                                continue
                            if np.min(list_F_spring) == 0:
                                continue
                            # 计算F_N列表中的F_N的最大最小数值，平均值，方差
                            arr_F_N = np.array(F_N)
                            #  设置优选条件
                            if np.min(arr_F_N) == 0:
                                print("最小力为0")
                                continue
                            # if np.max(arr_F_N) > 12:
                            #     print("最大力太大")
                            #     continue

                            # 疲劳条件设计
                            # arr_F_spring = np.array(list_F_spring)
                            # K_quduxishu = (4*C-1)/(4*C-4)+0.615/C
                            # F_S = 0.5*R_m*np.pi*(d**3)/(K_quduxishu*8*D)
                            # if np.max(arr_F_spring) > F_S:
                            #     continue

                            if np.mean(arr_F_N) > find_F_N_mean and find_F_N_mean < 10:
                                find_F_N_mean = np.mean(arr_F_N)
                                find_F_N_var = np.var(arr_F_N)
                                find_x1 = x1
                                find_y1 = y1
                                find_x2 = x2
                                find_y2 = y2
                                find_n = n
                                find_D = D
                                find_d = d
                            if np.var(arr_F_N) < find_F_N_var and np.mean(arr_F_N) > 10:
                                find_F_N_mean = np.mean(arr_F_N)
                                find_F_N_var = np.var(arr_F_N)
                                find_x1 = x1
                                find_y1 = y1
                                find_x2 = x2
                                find_y2 = y2
                                find_n = n
                                find_D = D
                                find_d = d

                        except:
                            break
                except:
                    break

# plt.show()
if D > 0:
    # 输出优选参数
    print("弹簧参数为D={},d={},n={}".format(find_D, find_d, find_n))
    print("此弹簧固定点（x1,y1)=({},{})".format(find_x1, find_y1))
    print("此弹簧铰接点（x2,y2)=({},{})".format(find_x2, find_y2))
    print("此时对地压力平均值为{}".format(find_F_N_mean))
    print("此时对地压力方差值为{}".format(find_F_N_var))

    # 计算弹簧铰接点与转轴的距离r
    r = np.sqrt((x0-find_x2)**2+(y0-find_y2)**2)
    # 计算弹簧固定点与转轴的间距H
    H = np.sqrt((x0-find_x1)**2+(y0-find_y1)**2)

    print("r={},H={}".format(r, H))

    if find_d == 0:
        print('未找到')
    else:
        find_C = find_D/find_d
        k = G*find_d/(8*(find_C**3)*find_n)

        # 计算当前状态下的弹簧长度和弹簧初始长度
        L_spring = np.sqrt((find_x1-find_x2)**2+(find_y1-find_y2)**2)
        L0_spring = (find_n+1)*find_d
        plt.figure()
        for height_wheeltochassis in np.arange(10, 36):
            # 计算r和H的初始角度值
            alpha = np.arccos((r**2+H**2-L_spring**2)/(2*r*H))
            print("alpha={},{}".format(alpha, alpha/np.pi*180))
            alpha_delta = np.arcsin(
                (height_wheeltochassis-10+11.4)/s)-np.arcsin(11.4/s)
            alpha -= alpha_delta
            # 计算弹簧弹力
            L_spring_now = np.sqrt(r**2+H**2-2*r*H*np.cos(alpha))
            F_spring = k*(L_spring_now-L0_spring)

            # 计算弹力力臂
            beta = np.arccos((L_spring_now**2+H**2-r**2)/(2*L_spring_now*H))
            arm_F_spring = H*np.sin(beta)
            # 弹簧对转轴施加的扭矩
            T = arm_F_spring * F_spring
            # 对地压力反力力臂
            arm_F_N = s * np.cos(theta)
            # 对地压力
            print("弹簧长度:{}，对地压力为:{}".format(L_spring_now, T/arm_F_N))
            plt.scatter(height_wheeltochassis, T/arm_F_N)

# plt.show()
