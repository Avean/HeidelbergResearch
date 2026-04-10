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


def DKK_model_three_equs():
    a = [1.5, 1.2, 1.5, 4.0]
    func_f = f"{a[3]}*y / ((1 + z)*(1 + {a[1]}*x)) - x"
    func_g = f"{a[0]}*x**2 - y"
    func_h = f"y / (1 + {a[2]}*x) - z"
    initial_guess_steady_state = [3.62960198320098, 19.76101583468475, 3.066384258084593] #[1, a[0]*1.0**2, a[0]*1.0**2 / (1 + a[2]*1.0)]
    return [func_f, func_g, func_h], initial_guess_steady_state

def vegetation_model_subcrit_PF():
    beta = [0.9090909090909092, 0.4722550177095633, 8.586454867446605, 0.9090909090909092]
    # For Diff_2 = 1.0, max bif point is subcritical for wave number 1, rest bif points supercritical.
    # For Diff_2 = 2.0  max bif point is supercritical for wave number 1, all superctritical.
    # For Diff_2 = 0.5, max bif point is supercritical for wave number 1, all superctritical.
    # For Diff_2 = 0.5, Diff_2 = 2.0 TP, for Diff_2 = 1.0 not? No, also TP but after fold bif, so shifted peak for x and y.
    # old:[0.7936507936507936, 0.9358240650757659, 38.621310622174455, 0.7936507936507936]
    alp = beta[0]
    bet = beta[1]
    gam = beta[2]
    lam = beta[3]
    func_f = f"{alp}*y - (1 - {bet}*y)*x"  # f"-{alp}*x + {bet}*(x*y)/(1+x*y)"
    func_g = f"{gam}*y**2*z - ({alp} + {bet}*x)*y"
    func_h = f"1 - {gam}*y**2*z - {lam}*z"
    Delta = np.sqrt(gam*(gam - 4*(alp+bet)*alp*lam))
    initial_guess_steady_state = [(-2*alp*bet*lam + gam + Delta) / (2*(bet**2*lam + gam)), (gam + Delta) / (2*(alp + bet)*gam), (2*(alp+bet)*bet*lam + gam - Delta) / (2*lam*(bet**2*lam+gam))]
    return [func_f, func_g, func_h], initial_guess_steady_state

def mod_vegetation_model_SDDI():
    # Modify vegetation model such that zero branch in f exists.
    beta = [0.5, 0.2, 0.5, 5.0, 1.0]
    #S-DDI (D_2 = 1.0): [0.5, 0.2, 0.5, 5.0, 1.0]:
    # but zero branch has autocatalysis condition at y_bar = 2.0797,
    # so for stable jump patterns: y has to be sufficiently far away from y_bar (y < alp/bet = 0.5/0.5 = 1.0).
    # Also for Jacobian to be stable at zero branch, we need y<=0.5.
    alp = beta[0]
    alp_2 = beta[1]
    bet = beta[2]
    gam = beta[3]
    lam = beta[4]
    func_f = f"-{alp}*x + {bet}*(x*y)/(1+x*y)"
    func_g = f"{gam}*y**2*z - ({alp_2} + {bet}*x)*y"
    func_h = f"1 - {gam}*y**2*z - {lam}*z"
    Delta = np.sqrt(gam*(gam - 4*(alp+bet)*alp*lam))
    initial_guess_steady_state = [(-2*alp*bet*lam + gam + Delta) / (2*(bet**2*lam + gam)), (gam + Delta) / (2*(alp + bet)*gam), (2*(alp+bet)*bet*lam + gam - Delta) / (2*lam*(bet**2*lam+gam))]
    return [func_f, func_g, func_h], initial_guess_steady_state

def mod_vegetation_model_WDDI():
    # Modify vegetation model such that zero branch in f exists.
    beta = [0.5, 0.5, 0.5, 5.0, 1.0]
    #Wave-DDI (D_2 = 1.0):
    # Max bif point of p1p2-p3=0: 0.00625 for j=2. d_j^+ > 0 for j>=2.
    # Max bif point of p3=0: 0.00182 for j=3. d_j > 0 for j>=2. B(3)=-0.2255, B(j)<0 for j>=2.
    alp = beta[0]
    alp_2 = beta[1]
    bet = beta[2]
    gam = beta[3]
    lam = beta[4]
    func_f = f"-{alp}*x + {bet}*(x*y)/(1+x*y)"
    func_g = f"{gam}*y**2*z - ({alp_2} + {bet}*x)*y"
    func_h = f"1 - {gam}*y**2*z - {lam}*z"
    Delta = np.sqrt(gam*(gam - 4*(alp+bet)*alp*lam))
    initial_guess_steady_state = [(-2*alp*bet*lam + gam + Delta) / (2*(bet**2*lam + gam)), (gam + Delta) / (2*(alp + bet)*gam), (2*(alp+bet)*bet*lam + gam - Delta) / (2*lam*(bet**2*lam+gam))]
    return [func_f, func_g, func_h], initial_guess_steady_state

def mod_vegetation_model_WDDI_old():
    # Modify vegetation model such that zero branch in f exists.
    beta = [0.5, 0.5, 0.47, 8.6, 0.9]
    #Wave-DDI (D_2 = 1.0): [0.5, 0.47, 8.6, 0.9]
    # Max bif point of p1p2-p3=0: 0.002875 for j=3. d_j^+ > 0 for j>=2.
    # Max bif point of p3=0: 0.00115 for j=4. d_j > 0 for j>=3. B(4)=-0.22, B(j) < 0 for j>=3.
    alp = beta[0]
    alp_2 = beta[1]
    bet = beta[2]
    gam = beta[3]
    lam = beta[4]
    func_f = f"-{alp}*x + {bet}*(x*y)/(1+x*y)"
    func_g = f"{gam}*y**2*z - ({alp_2} + {bet}*x)*y"
    func_h = f"1 - {gam}*y**2*z - {lam}*z"
    Delta = np.sqrt(gam*(gam - 4*(alp+bet)*alp*lam))
    initial_guess_steady_state = [(-2*alp*bet*lam + gam + Delta) / (2*(bet**2*lam + gam)), (gam + Delta) / (2*(alp + bet)*gam), (2*(alp+bet)*bet*lam + gam - Delta) / (2*lam*(bet**2*lam+gam))]
    return [func_f, func_g, func_h], initial_guess_steady_state

def mod_vegetation_model_autocat():
    # Modify vegetation model such that zero branch in f exists.
    beta = [0.05, 0.2, 0.5, 8.0, 1.0]
    # Steady state: (0.9633217955403929, 1.8119101848200547, 0.03667820445960722)
    # Jacobian: [[0.07019289, 0.06390186, 0.],[-0.90595509,0.3816609,26.26414814],[0., -1.0633218, -27.26414814]]
    # Unstable subsystems: J_1, J_2, J_12, J_13.
    alp = beta[0]
    alp_2 = beta[1]
    bet = beta[2]
    gam = beta[3]
    lam = beta[4]
    func_f = f"-{alp}*x + {bet}*(x*y)/(1+x*y)"
    func_g = f"{gam}*y**2*z - ({alp_2} + {bet}*x)*y"
    func_h = f"1 - {gam}*y**2*z - {lam}*z"
    Delta = np.sqrt(gam*(gam - 4*(alp+bet)*alp*lam))
    initial_guess_steady_state = [(-2*alp*bet*lam + gam + Delta) / (2*(bet**2*lam + gam)), (gam + Delta) / (2*(alp + bet)*gam), (2*(alp+bet)*bet*lam + gam - Delta) / (2*lam*(bet**2*lam+gam))]
    return [func_f, func_g, func_h], initial_guess_steady_state