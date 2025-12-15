# -*- coding: utf-8 -*-
# 此为第一题的第一问，使用newmark_beta法
import numpy as np
import pandas as pd

# 浮子与振子在初始时刻的状态
x_fu = 0.0              # 浮子的位移，这里简单说一下，即相对于平衡位置的位移，其>0则表示浮子吃水加深
v_fu = 0.0              # 浮子的速度
a_fu = 0.0

x_zh = 0.0              # 振子的位移
v_zh = 0.0              # 振子的速度
a_zh = 0.0              # 振子的加速度

wave_fn = 1.4005                                  # 波浪的圆频率
m_fu = 4866                                       # 浮子的质量
m_fu_add = 1335.535                               # 浮子的附加质量
m_zhen = 2433                                     # 振子的质量
pto_k = 80000                                         # pto弹簧刚度
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

# 想象成弹簧，求一下静水的等效刚度
def water_k(x_fu):
    return np.pi*1*1 * 1025*9.8


def newmark(t, x_fu, v_fu, a_fu, x_zh, v_zh, a_zh):
    x_fu_old = x_fu
    v_fu_old = v_fu
    a_fu_old = a_fu

    x_zh_old = x_zh
    v_zh_old = v_zh
    a_zh_old = a_zh

    # 对于浮子
    k_fu = (water_k(x_fu_old) + pto_k) + a0*(m_fu + m_fu_add) + a1*(c_pto + c_wave)
    p_fu = wave_force(t) + (a0*x_fu_old + a2*v_fu_old + a3*a_fu_old)*(m_fu + m_fu_add) + (a1*x_fu_old + a4*v_fu_old + a5*a_fu_old)*(c_wave + c_pto)
    x_fu_new = p_fu/k_fu
    a_fu_new = a0*(x_fu_new - x_fu_old) - a2*v_fu_old - a3*a_fu_old
    v_fu_new = v_fu_old + a6*a_fu_old + a7*a_fu_new

    # 对于振子，其阻尼与弹簧刚度与浮子的互为作用力与反作用力
    k_zh= -pto_k + a0*m_zhen - a1*c_pto


    return x_fu_new, v_fu_new, x_zh_new, v_zh_new


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