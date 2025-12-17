# -*- coding: utf-8 -*-
# 此为第二题的第二问，使用四阶龙格库塔法
import numpy as np

# 浮子与振子在初始时刻的状态
x_fu = 0.0              # 浮子的位移
x_zh = 0.0              # 振子的位移
v_fu = 0.0              # 浮子的速度
v_zh = 0.0              # 振子的速度

wave_fn = 2.2143                                  # 波浪的圆频率
m_fu = 4866                                       # 浮子的质量
m_fu_add = 1165.992                               # 浮子的附加质量
m_zhen = 2433                                     # 振子的质量
k = 80000                                         # pto弹簧刚度
c_wave = 167.8395                                 # 波浪的兴波阻尼系数
# c_pto                                 --------- # pto的阻尼系数,这个是代求值
# kkk                                   --------- # 幂指数取值

dt = 0.01                                        # 时间步长
end = 500                                         # 结束时间
power_start = 400.0
# 假设400s之后是稳定的，在计算400s到500s之间的平均功率

# 波浪激励力的计算
def wave_force(t):
    return 4890*np.cos(wave_fn*t)

# 静水恢复力的计算,disp为相对于平衡位置的位移，
def water_force(disp):
    if disp <= 2:
        v_water = np.pi*1*1 * disp
    else:
        # 大圆锥的体积
        v1 = np.pi*1*1*0.8/3
        # 小圆锥的体积
        v2 = pow((2.8-disp)/0.8, 2)*np.pi*(2.8-disp)/3
        v_water = 2*np.pi*1*1 + (v1 - v2)
    return v_water * 1025 * 9.8


# 计算浮子的加速度
def acc_fu(time, x_fu, v_fu, x_zh, v_zh, c_pto, kkk):
    return (wave_force(time) - c_wave*v_fu - c_pto*pow(abs(v_fu-v_zh), kkk)*(v_fu-v_zh) - water_force(x_fu) - k*(x_fu-x_zh))/(m_fu + m_fu_add)


# 计算振子的加速度，振子与浮子，关于阻尼、弹簧产生的力互为相互作用力
def acc_zh(time, x_fu, v_fu, x_zh, v_zh, c_pto, kkk):
    return (c_pto*pow(abs(v_fu-v_zh), kkk)*(v_fu-v_zh)+k*(x_fu-x_zh))/m_zhen


# 等价的方法
def deriv(run_time, y, c_pto, kkk):
    x_fu = y[0]
    v_fu = y[1]
    x_zh = y[2]
    v_zh = y[3]
    a_fu = acc_fu(run_time, x_fu, v_fu, x_zh, v_zh, c_pto, kkk)
    a_zh = acc_zh(run_time, x_fu, v_fu, x_zh, v_zh, c_pto, kkk)

    return np.array([v_fu, a_fu, v_zh, a_zh], dtype=float)

# 四阶龙格库塔法主要做的事情是，给定当前时刻的t以及当前时刻的系统状态，使用四阶龙格库塔法将此状态推进一个步长，得到下一时刻的状态
def rk4_step(t, y, h, c_pto, kkk):
    k1 = deriv(t, y, c_pto, kkk)
    k2 = deriv(t+0.5*h, y+0.5*h*k1, c_pto, kkk)
    k3 = deriv(t+0.5*h, y+0.5*h*k2, c_pto, kkk)
    k4 = deriv(t+h, y+h*k3, c_pto, kkk)
    return y + h*(k1 + 2*k2 + 2*k3 + k4)/6

def cal_power(c_pto, kkk):
    count = 0       # 时间步
    power_sum = 0.0 # 功率
    time = 0.0      # 起始时间
    y = np.array([0.0, 0.0, 0.0, 0.0], dtype=float)
    steps = int(np.floor(end / dt))

    for zz in range(steps):
        y = rk4_step(time, y, dt, c_pto, kkk)
        time = time+dt

        if time >= power_start:
            v_rel = y[1] - y[3]
            pp = c_pto * (abs(v_rel) ** kkk) * (v_rel ** 2)
            power_sum = power_sum + pp
            count = count + 1

    return power_sum / count



if __name__ == '__main__':

    gol_best_c = 0.0        # 最佳的阻尼
    gol_best_k = 0.0        # 最佳的指数
    gol_best_p = -1.0       # 最佳的功率

    for kkk in np.arange(0.0, 1.0001, 0.1):
        best_c = 0.0
        best_p = -1.0
        evaluated = {}

        for c in np.arange(0, 100001, 1000):
            P = cal_power(c, kkk)
            evaluated[c] = P
            if P > best_p:
                best_p = P
                best_c = c


        for c in range(int(max(0, best_c - 1000)), int(min(100000, best_c + 1000)) + 1, 100):
            if c in evaluated:
                continue
            P = cal_power(c, kkk)
            evaluated[c] = P
            if P > best_p:
                best_p = P
                best_c = c


        for c in range(int(max(0, best_c - 100)), int(min(100000, best_c + 100)) + 1, 10):
            if c in evaluated:
                continue
            P = cal_power(c, kkk)
            evaluated[c] = P
            if P > best_p:
                best_p = P
                best_c = c

        for c in range(int(max(0, best_c - 10)), int(min(100000, best_c + 10)) + 1, 1):
            if c in evaluated:
                continue
            P = cal_power(c, kkk)
            evaluated[c] = P
            if P > best_p:
                best_p = P
                best_c = c

        if best_p > gol_best_p:
            gol_best_p = best_p
            gol_best_c = best_c
            gol_best_k = kkk

    print("最优阻尼系数:{}".format(gol_best_c))
    print("最优幂指数:{}".format(gol_best_k))
    print("对平均输出功率:{}".format(gol_best_p))
