# -*- coding: utf-8 -*-
# 此为第一题的第一问，使用newmark_beta法
import numpy as np
import pandas as pd

# 浮子与振子在初始时刻的状态
x_fu_old = 0.0              # 浮子的位移，这里简单说一下，即相对于平衡位置的位移，其>0则表示浮子吃水加深
v_fu_old = 0.0              # 浮子的速度
a_fu_old = 0.0
x_zh_old = 0.0              # 振子的位移
v_zh_old = 0.0              # 振子的速度
a_zh_old = 0.0

x_fu_new = 0.0              # 浮子的位移
v_fu_new = 0.0              # 浮子的速度
a_fu_new = 0.0
x_zh_new = 0.0              # 振子的位移
v_zh_new = 0.0              # 振子的速度
a_zh_new = 0.0

wave_fn = 1.4005                                  # 波浪的圆频率
m_fu = 4866                                       # 浮子的质量
m_fu_add = 1335.535                               # 浮子的附加质量
m_zhen = 2433                                     # 振子的质量
k = 80000                                         # pto弹簧刚度
c_wave = 656.3616                                 # 波浪的兴波阻尼系数
c_pto = 10000                                     # pto的阻尼系数,不是定值,第一问这是个定制

beta = 0.25
gma = 0.5

time = 0                                          # 开始时间
dt = 0.001                                        # 时间步长
end = 200                                         # 结束时间

a0 = 1 / (beta * dt * dt)
a1 = gma / (beta * dt)
a2 = 1 / (beta * dt)
a3 = 1 / (2 * beta) - 1
a4 = gma / beta - 1
a5 = (gma / (2 * beta) - 1) * dt
a6 = dt * (1 - gma)
a7 = dt * gma


# 波浪激励力的计算，振子的力与浮子的力互为作用力与反作用力，振子的力为负的
def wave_force(t):
    return 6250*np.cos(wave_fn*t)

# 静水恢复力的计算,可以把海水想象成一个弹簧,这里的静水恢复力是变化的量，关于圆柱的变化
def water_force(x_fu):
    if 0 < x_fu <= 2:
        v_water = np.pi*1*1 * x_fu
    else:
        # 大圆锥的体积
        v1 = np.pi*1*1*0.8/3
        # 小圆锥的体积
        v2 = pow((2.8-x_fu)/0.8, 2)*np.pi*(2.8-x_fu)/3
        v_water = 2*np.pi*1*1 + v1 - v2
    return v_water * 1025 * 9.8


# 计算浮子的加速度
def acc_fu(time, x_fu, v_fu, x_zh, v_zh):
    return (wave_force(time) - c_wave*v_fu - c_pto*(v_fu-v_zh) - water_force(x_fu)-k*(x_fu-x_zh))/(m_fu + m_fu_add)

# 计算振子的加速度
def acc_zh(time, x_fu, v_fu, x_zh, v_zh):
    return (c_pto*(v_fu-v_zh)+k*(x_fu-x_zh))/m_zhen


def rk4(time, x_fu, v_fu, x_zh, v_zh):
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
    v_fu_new = v_fu + dt*(fu_K1 + 2*fu_K2 + 2*fu_K3 + fu_K4)/6
    x_fu_new = x_fu + dt*(fu_L1 + 2*fu_L2 + 2*fu_L3 + fu_L4)/6

    v_zh_new = v_zh + dt*(zh_K1 + 2*zh_K2 + 2*zh_K3 + zh_K4)/6
    x_zh_new = x_zh + dt*(zh_L1 + 2*zh_L2 + 2*zh_L3 + zh_L4)/6

    return v_fu_new, x_fu_new, v_zh_new, x_zh_new


t_list = []
x_fu_list, v_fu_list, x_zh_list, v_zh_list = [], [], [], []
while time <= end:
    v_fu, x_fu, v_zh, x_zh = rk4(time, x_fu, v_fu, x_zh, v_zh)
    t_list.append(time)
    x_fu_list.append(x_fu)
    v_fu_list.append(v_fu)
    x_zh_list.append(x_zh)
    v_zh_list.append(v_zh)
    time = time + dt


# 写入数据
sample_step = int(0.2/dt)   # 0.2 / 0.001 = 200

rows = []
for i, t in enumerate(t_list):
    if i % sample_step == 0:
        rows.append([t,x_fu_list[i], v_fu_list[i],x_zh_list[i], v_zh_list[i]])

df = pd.DataFrame(rows,columns=['t(s)','x_fu(m)', 'v_fu(m/s)','x_zh(m)', 'v_zh(m/s)'])

df.to_excel("./result1-1.xlsx", index=False)

print("计算完成！")