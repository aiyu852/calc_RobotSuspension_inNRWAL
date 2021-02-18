import numpy as np
from matplotlib import pyplot as plt
plt.rcParams['font.sans-serif'] = ['SimHei']  # 用来正常显示中文标签
plt.rcParams['axes.unicode_minus'] = False  # 用来正常显示负号


'''
    ===========================================
    |             规定初始参数值               |
    |    p 后缀代表针轮   c 后缀代表摆线轮      |
    |  a 为偏心距  b 为摆线轮厚度 k1 为短副系数 |
    |         T=减速比 乘 输入扭矩             |
    |  F_max 存储每一步迭代计算出的最大应力值   |
    ===========================================
'''
r_p = 25.5
z_p = 16
z_c = z_p-1
i_h = z_p/z_c
r_rp = 2.5
a = 0.9
k1 = a*z_p/r_p
b = 4.5

T = 25785  # N·mm
F_a = 2.2*T/k1/z_c/r_p  # 理想受力
F_b = 0.55*T/a/z_c  # 单齿受力
F_max = [0]*100  # 初始化最大啮合力，使其成为一个列表对象
F_max[0] = (F_a + F_b) / 2  # 得到初始迭代值F_max[0]

miu = 0.3  # 材料的泊松比
E = 206000  # 弹性模量
d_z = 3  # 针齿销轴的直径
l_z = 16  # 针齿支点的间距
# 以上两个参数建议参考李力行《摆线针轮行星传动的齿形修正及受力分析》

'''
    ===========================================
    |             最佳修形量公式               |
    |           D_j 是给定的径向间隙           |
    |      修形量的计算公式出自关天民论文       |
    ===========================================
'''
D_j = 0.01
D_Rrp = D_j/(1-(1-k1**2)**(1/2))
D_Rp = -D_j*(1-k1**2)**(1/2)/(1-(1-k1**2)**(1/2))
# 输出修形量
print("正等距修形量Δr_rp：", D_Rrp)
print("负移距修形量Δr_p：", D_Rp)

'''
    #######公式中参数我将以latex格式描述########
    ===========================================
    calc_delta_max
    -------------------------------------------
    输入：(r_p, k1, z_p, r_rp, b, F_max)
    返回: int or float
    作用：计算最大变形量 \delta_{max}
    变量:   rou： 计算 \rho 齿廓曲率半径
            miu： 定义 \mu 材料泊松比
            E: 定义 弹性模量
            c: 计算中间量
            w_max: 计算接触点公法线方向总接触变形
            J: 计算中间量
            f_max: 计算针齿弯曲变形
            delta_max: 输出最大变形量
    ===========================================
    ===========================================
    clac_delta_i
    -------------------------------------------
    输入：(delta_max, k1)
    返回: list
    作用：计算每个啮合点的总变形量 \delta_{i}
    变量：   delta: 初始化数列用于存储各个啮合处的接触总变形 \delta_{i}
                j: 角度弧度制
    ===========================================
    ===========================================
    find_phi_mn
    -------------------------------------------
    输入：(delta, Q)
    返回: list
    作用：计算总变形量大于初始间隙的角度区间 (m,n)
    变量：   phi_mn: 初始化数列用于存储各个啮合处的接触总变形 \delta_{i}
            Q: 初始间隙值的列表
    ===========================================
    ===========================================
    find_mANDn
    -------------------------------------------
    输入：(phi_mn)
    返回: int, int 
    作用：计算总变形量大于初始间隙的角度区间内真正首个啮合的针齿位置，区间内发生啮合的针齿位置数量
    变量：   phi_m: 根据每个针齿的相距角度遍历 得到处于啮合区间的第一个啮合针齿的角度位置
            m_n: 得到在啮合区间内的针齿数量
    ===========================================
    ===========================================
    calc_F_max
    -------------------------------------------
    输入：(m_n, a, z_c, phi_m, z_p, k1, D_Rp, D_Rrp, delta_max)
    返回: int or float 
    作用：根据每个啮合针齿位置总变形计算得到最大的啮合力
    变量：   r_c: 计算摆线轮节圆半径 r^{\prime}_{c}
            j: 区间内第i个针齿的角度位置
            k: 换算弧度制
            l_i: 第i个针齿啮合时到摆线轮转动中心的距离
            delta_phi_i: 计算第i个针齿啮合位置摆线轮的初始间隙
            U: 中间量
            T: 减速比乘以输入扭矩
            F_max：根据每个啮合针齿位置总变形计算得到的最大啮合力
    ===========================================
    ===========================================
    calc_Fi
    -------------------------------------------
    输入：(j)
    返回: int or float, int or float
    作用：根据最大的啮合力计算给定位置的变形量与初始间隙
    变量：   k: 根据输入的j转换弧度制
            delta_i: 根据最大总变形量计算得到k位置的总变形量
            Q_i: 计算k位置摆线轮的初始间隙
    ===========================================
    ===========================================
    calc_Q
    -------------------------------------------
    输入：(k)
    返回: int or float
    作用：计算给定k位置的初始间隙
    变量：  S_i: 中间量
            Q_i: 计算初始间隙 \Delta(\phi_{i})
    ===========================================
'''


def calc_delta_max(r_p, k1, z_p, r_rp, b, F_max):
    rou = r_p*(1+k1**2-2*k1**2)**(3/2) / \
        np.abs(k1*(z_p+1)*k1-1-z_p*k1**2)+r_rp  # +D_Rrp
    # miu = 0.3
    # E = 206000
    c = 4.99/1000*np.sqrt(2*(1-miu**2)/E*F_max/b*2*rou*r_rp/(r_rp+rou))
    w_max = 2*(1-miu**2)/E*F_max/np.pi/b*(2/3+np.log(16*r_rp*rou/c**2))

    # f_max
    J = np.pi * d_z**4/64
    # 其中3是针齿销（不是针齿销外面的轴套）的直径
    # 具体参考李力行的论文《摆线针轮行星传动的齿形修正及受力分析》

    f_max = F_max*31*l_z**3/(48*E*J*64)
    delta_max = f_max + w_max
    return delta_max


def clac_delta_i(delta_max, k1):
    delta = [0]*180
    for i in range(180):
        j = np.pi/180*i
        delta[i] = delta_max*np.sin(j)/np.sqrt(1+k1**2-2*k1*np.cos(j))
    return delta


def find_phi_mn(delta, Q):
    phi_mn = []
    for i in range(180):
        if delta[i] > Q[i]:
            phi_mn.append(i)
    return phi_mn


def find_mANDn(phi_mn):
    for i in range(48):
        phi_m = i * 360/z_p
        if phi_m >= phi_mn[0]:
            for j in range(len(phi_mn)):
                if (phi_m+j*360/z_p) > phi_mn[-1]:
                    m_n = j
                    break
            break
    return phi_m, m_n


def calc_F_max(m_n, a, z_c, phi_m, z_p, k1, D_Rp, D_Rrp, delta_max):
    U = 0
    for i in range(m_n):
        r_c = a*z_c
        j = phi_m+i*360/z_p
        k = j*np.pi/180
        S_i = 1/np.sqrt(1+k1**2-2*k1*np.cos(k))
        l_i = r_c*np.sin(k)*S_i
        delta_phi_i = calc_Q(k)
        U += l_i*(l_i/r_c-delta_phi_i/delta_max)
    F_max = 0.55 * T / U
    return F_max


def calc_Fi(j):
    k = j*np.pi/180
    delta_i = delta_max*np.sin(k)/np.sqrt(1+k1**2-2*k1*np.cos(k))
    Q_i = calc_Q(k)
    return delta_i, Q_i


def calc_Q(k):
    S_i = 1/np.sqrt(1+k1**2-2*k1*np.cos(k))
    Q_i = D_Rrp*(1-np.sin(k)*S_i)+D_Rp*S_i * \
        (1-k1*np.cos(k)-np.sqrt(1-k1**2)*np.sin(k))
    return Q_i


# 计算初始间隙
Q = [0]*180
for i in range(180):
    j = i*np.pi/180
    Q[i] = calc_Q(j)

plt.plot(Q, label=r'$\Delta\phi_i$')
plt.grid(True)

print("迭代次数\t\tF_max\t\tDelta_Fmax\t\tF_max/1000")

for p in range(1, 100):
    delta_max = calc_delta_max(r_p, k1, z_p, r_rp, b, F_max[p-1])
    delta = clac_delta_i(delta_max, k1)
    phi_mn = find_phi_mn(delta, Q)
    phi_m, m_n = find_mANDn(phi_mn)
    F_max[p] = calc_F_max(m_n, a, z_c, phi_m, z_p,
                          k1, D_Rp, D_Rrp, delta_max)
    # 计算迭代结果与输入量的差值
    delta_F_max = np.abs(F_max[p]-F_max[p-1])
    print("\t{0}\t{1}\t{2}\t{3}".format(
        p, F_max[p], delta_F_max, F_max[p]/1000))
    if delta_F_max < F_max[p]/1000:
        F_max[-1] = (F_max[p] + F_max[p-1]) / 2
        break

# 根据最后得到的真实最大啮合力计算真正的啮合区间
delta_max = calc_delta_max(r_p, k1, z_p, r_rp, b, F_max[-1])
delta = clac_delta_i(delta_max, k1)
phi_mn = find_phi_mn(delta, Q)
phi_m, m_n = find_mANDn(phi_mn)

# 将真实啮合的针齿位置用彩色点描出
F_i = []
print("受力针齿编号\t针齿上的啮合力")
for i in range(m_n):
    j = phi_m+i*360/z_p
    delta_i, Q_i = calc_Fi(j)
    Fi = (delta_i-Q_i)/delta_max*F_max[-1]
    print("\t{0}\t{1}".format(i+1, Fi))
    F_i.append(Fi)
    # 绘图描点
    plt.scatter(j, delta_i)
    plt.scatter(j, Q_i)

# 绘图
plt.plot(delta, label=r'$\Delta\delta_i$')
plt.xlim([0, 180])
plt.ylim([0, 0.01])
plt.xlabel("角参量"+r'$\phi$/rad')
plt.ylabel("初始间隙"+r'$\Delta\phi_i$/mm'+'\n'+"变形量"+r'$\Delta\delta_i$/mm')
plt.legend()

rou = r_p*(1+k1**2-2*k1**2)**(3/2) / \
    np.abs(k1*(z_p+1)*k1-1-z_p*k1**2)+r_rp
rou_ei = np.abs(rou*r_rp/(rou-r_rp))
sigma = 0.418*(206000*F_max[-1]/b/rou_ei)**(1/2)

print("最大齿面接触强度：", sigma)

plt.show()