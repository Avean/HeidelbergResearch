import numpy as np

def extended_brusselator_model():
    # j_22 and J_12 (type-3) unstable:
    p = 1.5 #5 # 1.5
    q = 3.5 #15  # 3.5
    # J_12 (type-2) unstable
    # p = 1.5
    # q = 1.5
    # J_12 (type-2) unstable
    # p = 1.5
    # q = 1.7
    # j_22, J_12 (type-3), J_23 (type-1) unstable
    # p = 1.2
    # q = 3.5
    # j_22, J_12 (type-2) unstable
    # p = 1.5
    # q = 2.7
    # wierd stuff happening
    # p = 5
    # q = 15
    func_f = f"(y - x)"  # "(y**2 - x**2)"
    func_g = f"{p} - {2+q}*y + y**2*z + x"
    func_h = f"{q} * y - y**2*z"
    return [func_f, func_g, func_h], [p, p, q/p]

def extended_brusselator_model_v2():
    p = 1.5
    q = 1.5
    func_f = f"{p} - {2+q}*x + x**2*z + y"
    func_g = f"x - y"
    func_h = f"{q} * x - x**2*z"
    return [func_f, func_g, func_h], [p, p, q/p]

def receptor_based_model():
    # muf = 0.716 # 0.95  # 0.87
    # mub = 0.731  # 0.7  # 0.68
    # mul = 0.945
    # mue = 0.172
    # m1 = 6.766
    # m2 = 16.764
    # m3 = 12.219
    muf = 0.54677082 # 0.95  # 0.87
    mub = 0.63416923  # 0.7  # 0.68
    mul = 0.20643638
    mue = 0.28009361
    m1 = 8.71467558
    m2 = 6.57493195
    m3 = 17.03622515
    ubar_guess = (m1*m2*mue - m3*muf + np.sqrt((-m1*m2*mue + m3*muf)**2 - 4*m1*m2*mue**2*muf*mul))/(2*m2*mue*muf)
    vbar_guess =  muf / (m1 - muf*ubar_guess)
    wbar_guess = m3*muf / (mue*m1) * ubar_guess
    func_f = f"-{muf}*x + {m1}*(x*y)/(1+x*y) - {mub}*x*y"
    func_g = f"-{mul}*y + {m2}*(x*y)/(1+x*y) - {mub}*x*y - y*z"
    func_h = f"-{mue}*z + {m3}*(x*y)/(1+x*y)"
    return [func_f, func_g, func_h], [ubar_guess, vbar_guess, wbar_guess]

def receptor_based_model_simpl():
    # mu1, mu2, mu3, m1, m2, m3 = 0.06605506,  0.97269018,  0.16120774,  7.07918872,  2.26488268, 10.66501855
    mu1, mu2, mu3, m1, m2, m3 = 1.0, 1.0, 0.6, 2.50, 9.68, 7.00
    ubar_guess = (m1*m2*mu3 - m3*mu1 + np.sqrt((-m1*m2*mu3 + m3*mu1)**2 - 4*m1*m2*mu3**2*mu1*mu2))/(2*m2*mu3*mu1)
    vbar_guess =  mu1 / (m1 - mu1*ubar_guess)
    wbar_guess = m3*mu1 / (mu3*m1) * ubar_guess
    func_f = f"-{mu1}*x + {m1}*(x*y)/(1+x*y)"
    func_g = f"-{mu2}*y + {m2}*(x*y)/(1+x*y) - y*z"
    func_h = f"-{mu3}*z + {m3}*(x*y)/(1+x*y)"
    return [func_f, func_g, func_h], [ubar_guess, vbar_guess, wbar_guess]

def some_model():
    k2 = 0.5
    k3 = 0.5
    k4 = 7
    k5 = 0.1
    k7 = 0.45
    k9 = 0.1
    gamma = 500
    d = 10
    h0 = 10
    p0 = 10
    q0 = 10
    fgf8 = 1
    hoxd13 = 1
    x0 = 1
    y0 = 1
    z0 = 1
    func_f = f"{k2}*(y-{y0}) - {k3}*(z-{z0}) - (x-{x0})**3"
    func_g = f"-({k4} - 2.5*{fgf8}*{hoxd13}*(x-{x0}) - {k5}*(y-{y0}))"
    func_h = f"-({k7} + 1.25*{fgf8}*{hoxd13}*(x-{x0}) - {k9}*(z-{z0}))"
    return [func_f, func_g, func_h], [1, 1, 1]

def linear_model():
    func_f = "-1*x - 2*y - 1*z"
    func_g = " 1*x + 1.5*y - 2*z"
    func_h = " 1*x + 2*y - 2*z"
    return [func_f, func_g, func_h], [0, 0, 0]

# initial_guess_steady_state = (1,1,1) # (ubar_guess, vbar_guess, wbar_guess) # (p, p, q/p) # (5, 5, 15/5)  # (1.7, 0.3, 9.6) #(0, 0, 0)
# guessed_steady_state = [1.6527232, 0.29941888, 9.5284198]
# print(evaluate_symbolic_at_value_three_vars(func_f, *steady_state))
# print(evaluate_symbolic_at_value_three_vars(func_g, *steady_state))
# print(evaluate_symbolic_at_value_three_vars(func_h, *steady_state))