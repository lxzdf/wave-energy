# -*- coding: utf-8 -*-
# 此为第一题的第一问，使用四阶龙格库塔法
import os
import numpy as np
import matplotlib.pyplot as plt

# 浮子与振子在初始时刻的状态
x_fu = 0.0              # 浮子的位移
x_zh = 0.0              # 振子的位移
v_fu = 0.0              # 浮子的速度
v_zh = 0.0              # 振子的速度

wave_fn = 1.4005                                  # 波浪的圆频率
m_fu = 4866                                       # 浮子的质量
m_fu_add = 1335.535                               # 浮子的附加质量
m_zhen = 2433                                     # 振子的质量
k = 80000                                         # pto弹簧刚度
c_wave = 656.3616                                 # 波浪的兴波阻尼系数
c_pto = 10000                                     # pto的阻尼系数,不是定值,第一问这是个定制

time = 0                                          # 开始时间
dt = 0.001                                        # 时间步长
end = 200                                         # 结束时间

# 波浪激励力的计算
def wave_force(t):
    return 6250*np.cos(wave_fn*t)

# 静水恢复力的计算
def water_force(t):
    return 1025*9.8*v_fu*t

# 计算浮子的加速度
def acc_fu(time, x_fu, v_fu, x_zh, v_zh):
    return (wave_force(time) - c_wave*v_fu - c_pto*(v_fu-v_zh) - water_force()-k*(x_fu-x_zh))/(m_fu + m_fu_add)

# 计算振子的加速度
def acc_zh(time, x_fu, v_fu, x_zh, v_zh):
    return 0-c_pto*(v_fu-v_zh)-k*(x_fu-x_zh)/m_zhen


# 浮子的速度(L)、加速度(K)----[t+Δt时刻]
fu_K1 = acc_fu(time, x_fu, v_fu, x_zh, v_zh)
fu_L1 = v_fu                                                                                                # 使用t时刻的速度
# 振子的速度(L)、加速度(K)----[t+Δt时刻]
zh_K1 = acc_zh(time, x_fu, v_fu, x_zh, v_zh)
zh_L1 = v_zh


fu_K2 = acc_fu(time+0.5*dt, x_fu+0.5*dt*fu_L1, v_fu+0.5*dt*fu_K1, x_zh+0.5*dt*zh_L1, v_zh+0.5*dt*zh_K1)
fu_L2 = v_fu + 0.5*dt*fu_K1
zh_K2 = acc_zh(time+0.5*dt, x_fu+0.5*dt*fu_L1, v_fu+0.5*dt*fu_K1, x_zh+0.5*dt*zh_L1, v_zh+0.5*dt*zh_K1)
zh_L2 = v_zh + 0.5*dt*zh_K1


fu_K3 = acc_fu(time+0.5*dt, x_fu+0.5*dt*fu_L2, v_fu+0.5*dt*fu_K2, x_zh+0.5*dt*zh_L2, v_zh+0.5*dt*zh_K2)
fu_L3 = v_fu + 0.5*dt*fu_K2
zh_K3 = acc_zh(time+0.5*dt, x_fu+0.5*dt*fu_L2, v_fu+0.5*dt*fu_K2, x_zh+0.5*dt*zh_L2, v_zh+0.5*dt*zh_K2)
zh_L3 = v_zh + 0.5*dt*zh_K2


fu_K4 = acc_fu(time+dt, x_fu+dt*fu_L3, v_fu+dt*fu_K3, x_zh+dt*zh_L3, v_zh+dt*zh_K3)
fu_L4 = v_fu + fu_K3*dt
zh_K4 = acc_zh(time+dt, x_fu+dt*fu_L3, v_fu+dt*fu_K3, x_zh+dt*zh_L3, v_zh+dt*zh_K3)
zh_L4 = v_zh + zh_K3*dt


# 计算得到t+Δt时刻的浮子、振子的速度与位移
v_fu = v_fu + dt*(fu_K1 + 2*fu_K2 + 2*fu_K3 + fu_K4)/6
x_fu = x_fu + dt*(fu_L1 + 2*fu_L2 + 2*fu_L3 + fu_L4)/6

v_zh = v_zh + dt*(zh_K1 + 2*zh_K2 + 2*zh_K3 + zh_K4)/6
x_zh = x_zh + dt*(zh_L1 + 2*zh_L2 + 2*zh_L3 + zh_L4)/6