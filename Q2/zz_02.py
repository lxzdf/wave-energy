# -*- coding: utf-8 -*-
# 此为第二题的第一问，使用四阶龙格库塔法
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar

wave_fn = 2.2143                                  # 波浪的圆频率
m_fu = 4866                                       # 浮子的质量
m_add = 1165.992                                  # 浮子的附加质量
m_zhen = 2433                                     # 振子的质量
c_wave = 167.8395                                 # 波浪的兴波阻尼系数
c_pto = None                                      # pto的阻尼系数,这个是代求值
kkk = None                                        # 幂指数,代求值

dt = 0.01                                         # 时间步长
end_time = 500                                         # 结束时间
power_start = 400.0
# 假设400s之后是稳定的，在计算400s到500s之间的平均功率
# 波浪激励力的计算
def wave_force(run_time):
    return 4890*np.cos(wave_fn*run_time)

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
    return (wave_force(time) - c_wave*v_fu - c_pto*pow(abs((v_fu-v_zh)), kkk)*(v_fu-v_zh) - water_force(x_fu)-80000*(x_fu-x_zh))/(m_fu + m_add)

# 计算振子的加速度
def acc_zh(time, x_fu, v_fu, x_zh, v_zh, c_pto, kkk):
    return (c_pto*pow(abs((v_fu-v_zh)), kkk)*(v_fu-v_zh)+80000*(x_fu-x_zh))/m_zhen


# y = [x_fu, v_fu, x_zh, v_zh]
def deriv(run_time, y, c_pto, kkk):
    x_fu = y[0]
    v_fu = y[1]
    x_zh = y[2]
    v_zh = y[3]
    a_fu = acc_fu(run_time, x_fu, v_fu, x_zh, v_zh, c_pto, kkk)
    a_zh = acc_zh(run_time, x_fu, v_fu, x_zh, v_zh, c_pto, kkk)

    return [v_fu, a_fu, v_zh, a_zh]


def cal_power(c_pto, kkk):
    t_list = np.arange(0.0, end_time+0.5*dt, dt)
    y0 = [0.0, 0.0, 0.0, 0.0]

    sol = solve_ivp(lambda t, y:deriv(t, y, c_pto, kkk),
                    t_span=(0.0, end_time),
                    y0=y0,  # 初始值
                    t_eval=t_list,
                    method="RK45",
                    rtol=1e-6,
                    atol=1e-6
                    )

    v_fu = sol.y[1]
    v_zh = sol.y[3]
    cal_time = sol.t

    power_sum = 0.0
    count = 0.0
    for i in range(len(cal_time)):
        if cal_time[i] > power_start:
            vel = v_fu[i] - v_zh[i]
            power_sum = power_sum + c_pto*pow(abs(vel), kkk)*vel*vel
            count = count + 1

    return power_sum/count


def best_for_k(kkk):
    res = minimize_scalar(lambda c: -1*cal_power(c, kkk),
                   bounds=(0.0, 100000.0),
                   method="bounded")
    best_c = res.x
    best_p = -res.fun

    return best_c, best_p



if __name__ == '__main__':
    best_k = None
    best_c = None
    best_p = -1.0

    k_grid = np.arange(0.0, 1.0001, 0.05)
    for k_exp in k_grid:
        c_opt, p_opt = best_for_k(k_exp)
        if p_opt > best_p:
            best_p = p_opt
            best_k = k_exp
            best_c = c_opt

    # 2) 细化 k：在最优附近 ±0.05，用步长 0.01 再搜一遍
    k_left = max(0.0, best_k - 0.05)
    k_right = min(1.0, best_k + 0.05)
    k_fine = np.arange(k_left, k_right + 1e-12, 0.01)
    for k_exp in k_fine:
        c_opt, p_opt = best_for_k(k_exp)
        if p_opt > best_p:
            best_p = p_opt
            best_k = k_exp
            best_c = c_opt

    print("最优阻尼系数:{}".format(best_c))
    print("最优幂指数:{}".format(best_k))
    print("对平均输出功率:{}".format(best_p))

