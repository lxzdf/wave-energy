# -*- coding: utf-8 -*-
# 此为第一题的第一问，使用四阶龙格库塔法，但是这是最新的方法
import numpy as np
import pandas as pd

# 浮子与振子在初始时刻的状态,即在平衡时的状态的状态
x_fu = 0.0              # 浮子的位移
x_zh = 0.0              # 振子的位移
v_fu = 0.0              # 浮子的速度
v_zh = 0.0              # 振子的速度

wave_fn = 1.4005                                  # 波浪的圆频率
m_fu = 4866                                       # 浮子的质量
m_add = 1335.535                                  # 浮子的附加质量
m_zhen = 2433                                     # 振子的质量
c_wave = 656.3616                                 # 波浪的兴波阻尼系数
c_pto = 10000                                     # pto的阻尼系数,不是定值,第一问这是个定值

time = 0                                          # 开始时间
dt = 0.01                                         # 时间步长
end = 200                                         # 结束时间

# 波浪激励力的计算
def wave_force(run_time):
    return 6250*np.cos(wave_fn*run_time)

# 静水恢复力的计算,可以把海水想象成一个弹簧,这里的静水恢复力是变化的量，关于圆柱的变化，此外可以验算一下，实际上这个变化只在圆柱段，不会到圆锥段的
# 静水恢复力的计算,disp为相对于平衡位置的位移，设向下为正
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
def acc_fu(time, x_fu, v_fu, x_zh, v_zh):
    return (wave_force(time) - c_wave*v_fu - c_pto*(v_fu-v_zh) - water_force(x_fu) - 80000*(x_fu-x_zh))/(m_fu + m_add)

# 计算振子的加速度，振子与浮子，关于阻尼、弹簧产生的力互为相互作用力
def acc_zh(time, x_fu, v_fu, x_zh, v_zh):
    return (c_pto*(v_fu-v_zh)+80000*(x_fu-x_zh))/m_zhen

# 等价的方法
def deriv(run_time, y):
    x_fu = y[0]
    v_fu = y[1]
    x_zh = y[2]
    v_zh = y[3]
    a_fu = acc_fu(run_time, x_fu, v_fu, x_zh, v_zh)
    a_zh = acc_zh(run_time, x_fu, v_fu, x_zh, v_zh)

    return np.array([v_fu, a_fu, v_zh, a_zh], dtype=float)

# 四阶龙格库塔法主要做的事情是，给定当前时刻的t以及当前时刻的系统状态，使用四阶龙格库塔法将此状态推进一个步长，得到下一时刻的状态
def rk4_step(t, y, h):
    k1 = deriv(t, y)
    k2 = deriv(t+0.5*h, y+0.5*h*k1)
    k3 = deriv(t+0.5*h, y+0.5*h*k2)
    k4 = deriv(t+h, y+h*k3)
    return y + h*(k1 + 2*k2 + 2*k3 + k4)/6


y = np.array([x_fu, v_fu, x_zh, v_zh], dtype=float)
t_list = [0.0]
x_fu_list, v_fu_list, x_zh_list, v_zh_list = [y[0]], [y[1]], [y[2]], [y[3]]

while time <= end:
    y = rk4_step(time, y, dt)  # 推到 time+dt
    time += dt

    t_list.append(time)
    x_fu_list.append(y[0])
    v_fu_list.append(y[1])
    x_zh_list.append(y[2])
    v_zh_list.append(y[3])


# 写入数据
sample_step = int(0.2/dt)   # 0.2 / 0.01 = 20

rows = []
for i, t in enumerate(t_list):
    if i % sample_step == 0:
        rows.append([t,x_fu_list[i], v_fu_list[i],x_zh_list[i], v_zh_list[i]])


df = pd.DataFrame(rows,columns=['t(s)','x_fu(m)', 'v_fu(m/s)','x_zh(m)', 'v_zh(m/s)'])
df.to_excel("./result1-1.xlsx", index=False)
print("计算完成！")