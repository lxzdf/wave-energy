# -*- coding: utf-8 -*-
# 此为第一题的第一问，使用四阶龙格库塔法
from math import trunc

import numpy as np
import matplotlib.pyplot as plt

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

time = 0                                          # 开始时间
dt = 0.01                                        # 时间步长
end = 500                                         # 结束时间

# 假设400s之后是稳定的，在计算400s到500s之间的平均功率

# 波浪激励力的计算
def wave_force(t):
    return 4890*np.cos(wave_fn*t)

# 静水恢复力的计算
def water_force(x_fu):
    # x_fu 就是 u1
    if x_fu <= 2:
        v_water = np.pi * x_fu
    else:
        v_water = 2 * np.pi + (1 - ((2.8 - x_fu)/0.8)**3) * (0.8 * np.pi / 3)
    return v_water * 1025 * 9.8



# 计算浮子的加速度
def acc_fu(time, x_fu, v_fu, x_zh, v_zh, c_pto, kkk):
    return (wave_force(time) - c_wave*v_fu - c_pto*pow(abs(v_fu-v_zh), kkk)*(v_fu-v_zh) - water_force(x_fu)-k*(x_fu-x_zh))/(m_fu + m_fu_add)

# 计算振子的加速度
def acc_zh(time, x_fu, v_fu, x_zh, v_zh, c_pto, kkk):
    return (c_pto*pow(abs(v_fu-v_zh), kkk)*(v_fu-v_zh)+k*(x_fu-x_zh))/m_zhen


def rk4(time, x_fu, v_fu, x_zh, v_zh, c_pto, kkk):
    # 浮子的速度(L)、加速度(K)----[t+Δt时刻]
    fu_K1 = acc_fu(time, x_fu, v_fu, x_zh, v_zh, c_pto, kkk)
    fu_L1 = v_fu                                                                                                # 使用t时刻的速度
    # 振子的速度(L)、加速度(K)----[t+Δt时刻]
    zh_K1 = acc_zh(time, x_fu, v_fu, x_zh, v_zh, c_pto, kkk)
    zh_L1 = v_zh


    fu_K2 = acc_fu(time+0.5*dt, x_fu+0.5*dt*fu_L1, v_fu+0.5*dt*fu_K1, x_zh+0.5*dt*zh_L1, v_zh+0.5*dt*zh_K1, c_pto, kkk)
    fu_L2 = v_fu + 0.5*dt*fu_K1
    zh_K2 = acc_zh(time+0.5*dt, x_fu+0.5*dt*fu_L1, v_fu+0.5*dt*fu_K1, x_zh+0.5*dt*zh_L1, v_zh+0.5*dt*zh_K1, c_pto, kkk)
    zh_L2 = v_zh + 0.5*dt*zh_K1


    fu_K3 = acc_fu(time+0.5*dt, x_fu+0.5*dt*fu_L2, v_fu+0.5*dt*fu_K2, x_zh+0.5*dt*zh_L2, v_zh+0.5*dt*zh_K2, c_pto, kkk)
    fu_L3 = v_fu + 0.5*dt*fu_K2
    zh_K3 = acc_zh(time+0.5*dt, x_fu+0.5*dt*fu_L2, v_fu+0.5*dt*fu_K2, x_zh+0.5*dt*zh_L2, v_zh+0.5*dt*zh_K2,c_pto, kkk)
    zh_L3 = v_zh + 0.5*dt*zh_K2

    fu_K4 = acc_fu(time+dt, x_fu+dt*fu_L3, v_fu+dt*fu_K3, x_zh+dt*zh_L3, v_zh+dt*zh_K3, c_pto, kkk)
    fu_L4 = v_fu + fu_K3*dt
    zh_K4 = acc_zh(time+dt, x_fu+dt*fu_L3, v_fu+dt*fu_K3, x_zh+dt*zh_L3, v_zh+dt*zh_K3, c_pto, kkk)
    zh_L4 = v_zh + zh_K3*dt

    # 计算得到t+Δt时刻的浮子、振子的速度与位移
    v_fu_new = v_fu + dt*(fu_K1 + 2*fu_K2 + 2*fu_K3 + fu_K4)/6
    x_fu_new = x_fu + dt*(fu_L1 + 2*fu_L2 + 2*fu_L3 + fu_L4)/6

    v_zh_new = v_zh + dt*(zh_K1 + 2*zh_K2 + 2*zh_K3 + zh_K4)/6
    x_zh_new = x_zh + dt*(zh_L1 + 2*zh_L2 + 2*zh_L3 + zh_L4)/6

    return v_fu_new, x_fu_new, v_zh_new, x_zh_new


def cal_power(c_pto, kkk):
    count = 0   # 时间步
    power_sum = 0.0 # 功率
    time = 0.0
    x_fu = 0.0
    v_fu = 0.0
    x_zh = 0.0
    v_zh = 0.0

    while time <= end:
        v_fu, x_fu, v_zh, x_zh = rk4(time, x_fu, v_fu, x_zh, v_zh, c_pto, kkk)
        if time >= 400:
            power_sum = power_sum + c_pto*pow(abs(v_fu-v_zh), kkk)*(v_fu-v_zh)*(v_fu-v_zh)
            count = count + 1
        time = time + dt

    return power_sum/count

if __name__ == '__main__':

    gol_best_c = 0.0
    gol_best_k = 0.0
    gol_best_p = -1.0

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
