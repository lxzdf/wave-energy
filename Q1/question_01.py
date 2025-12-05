# -*- coding: utf-8 -*-
# 此为第一题的第一问，使用四阶龙格库塔法
import os
import numpy as np
import matplotlib.pyplot as plt

x_fu = 0.0              # 浮子的位移
x_zh = 0.0              # 振子的位移
v_fu = 0.0              # 浮子的速度
v_zh = 0.0              # 振子的速度

wave_fn = 1.4005                                  # 波浪的圆频率
f_wave = 6250*np.cos(0.2*wave_fn)                 # 波浪激励力
f_re =  0                                         # 静水恢复力
m_fu = 4866                                       # 浮子的质量
m_fu_add = 1335.535                               # 浮子的附加质量
m_zhen = 2433                                     # 振子的质量
k = 80000                                         # pto弹簧刚度
c_pto = 10000                                     # pto的阻尼系数
c_wave = 656.3616                                 # 波浪的兴波阻尼系数

time = 0                                             # 开始时间
dt = 0.001                                        # 时间步长
end = 200                                         # 结束时间

# 计算浮子的加速度
def acc_fu(time, x_fu, v_fu, x_zh, v_zh):
    pass

# 计算振子的加速度
def acc_zh(time, x_fu, v_fu, x_zh, v_zh):
    pass


# 浮子的速度(L)、加速度(K)----[t+Δt时刻]
fu_K1 = acc_fu(time, x_fu, v_fu, x_zh, v_zh)
fu_L1 = v_fu                                                                                                # 使用t时刻的速度
# 振子的速度(L)、加速度(K)----[t+Δt时刻]
zh_K1 = acc_zh(time, x_fu, v_fu, x_zh, v_zh)
zh_L1 = v_zh


fu_k2 = acc_fu(time+0.5*dt, x_fu+0.5*dt*fu_L1, v_fu+0.5*dt*fu_K1, x_zh+0.5*dt*zh_L1, v_zh+0.5*dt*zh_K1)
fu_L2 = fu_L1 + 0.5*dt*fu_K1;
zh_K2 = acc_fu(time+0.5*dt, x_fu+0.5*dt*fu_L1, v_fu+0.5*dt*fu_K1, x_zh+0.5*dt*zh_L1, v_zh+0.5*dt*zh_K1)




