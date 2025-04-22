import numpy as np
# from scipy import integrate
from scipy.optimize import minimize
from scipy.optimize import minimize_scalar

global EPS
EPS = 1e-10

def calc_argmax_N_v2(kap):
    def objective(n):
        return -(kap**2 * n * (np.sqrt(2) - 1) - kap * np.log(4*n)) / (kap * 4*np.pi**2*n**2 * np.log(4*n))

    res = minimize_scalar(objective) # , bounds=(1.0, float('inf')), method='bounded'
    argmax_N = int(res.x)
    return argmax_N

def calc_argmax_N(kap):
    argmax_N = np.argmax([(kap**2* n * (np.sqrt(2) - 1) - kap * np.log(4*n)) \
                          / (kap * 4*np.pi**2*n**2 * np.log(4*n)) for n in range(1, 1000)])
    return argmax_N + 1

def calc_D_1(kap):
    D_1 = np.max([(kap * n * (np.sqrt(2) - 1) / np.log(4*n) - 1) / (4*np.pi**2*n**2) for n in range(1, 10000)])
    return D_1

def calc_D_1_v2(kap):
    argmax_N = calc_argmax_N_v2(kap)
    if argmax_N > 1:
        D_1_N_small = (kap * (argmax_N-1) * (np.sqrt(2) - 1) / np.log(4*(argmax_N-1)) - 1) / (4*np.pi**2*(argmax_N-1)**2)
    else:
        D_1_N_small = float('-inf')
    D_1_N_middle = (kap * argmax_N * (np.sqrt(2) - 1) / np.log(4*argmax_N) - 1) / (4*np.pi**2*argmax_N**2)
    D_1_N_large = (kap * (argmax_N+1) * (np.sqrt(2) - 1) / np.log(4*(argmax_N+1)) - 1) / (4*np.pi**2*(argmax_N+1)**2)
    return max(D_1_N_small, D_1_N_middle, D_1_N_large)

def calc_D_2(kap):
    D_2 = (kap - 1 ) / (4*np.pi**2)
    return D_2

def calc_D_min(kappa):
    D_min = max(calc_D_1(kappa), calc_D_2(kappa))
    return D_min

def calc_D_min_v2(kappa):
    D_min = max(calc_D_1_v2(kappa), calc_D_2(kappa))
    return D_min

def C(kap, d):
    summand = 0
    k = 0
    while summand > EPS:
        summand += 1/(1 + d*4*k**2*np.pi**2 - 2*kap)

def calc_D_max_lin(kappa):
    # TODO: Implement this function
    D_max = 15*kappa / (4*np.pi**2)
    # if C(kappa, D_max) > (8*kappa)^(-1)+EPS:
    #     print('Error in calculations!')
    # while C(kappa, D_max) < (8*kappa)^(-1):
    #     if C(kappa, D_max/2) < (8*kappa)^(-1):
    #         D_max = D_max/2
    #     elif C(kappa, 2*D_max/3) > (8*kappa)^(-1):
    #         D_max = 2*D_max/3
    #     else:
    #         D_max = D_max - 1e-2*D_max
    return D_max

def calc_D_max_v2(kappa):
    def C(kap, d):
        summand = 0
        k = 1
        while True:
            term = 1 / (1 + d * 4 * np.pi**2 * k**2 - 2 * kap)
            if term < EPS:
                break
            summand += term
            k += 1
        return summand

    D_upper_estim = 15 * kappa / (4 * np.pi**2)
    while C(kappa, D_upper_estim/2) < 1 / (8 * kappa):
        D_upper_estim /= 2

    D_lower_estim = D_upper_estim / 2
    while D_upper_estim - D_lower_estim > EPS:
        D_mid = (D_upper_estim + D_lower_estim) / 2
        if C(kappa, D_mid) < 1 / (8 * kappa):
            D_upper_estim = D_mid
        else:
            D_lower_estim = D_mid

    return D_upper_estim

# kappa_list = [k for k in np.linspace(1.0, 2, 100)]
# Diff_list = [calc_D_2(k)-calc_D_1(k) for k in kappa_list]
# print(kappa_list[np.argmax(Diff_list)])
# print(np.max(Diff_list))


# Define functional to minimize

def functional(R, kap, diff):
    func_1 = 1/2*R[0]**2 - kap * R[0]
    func_2 = sum([1/2 * R[i]**2 * (1+diff*4*np.pi**2*i**2) for i in range(1, len(R))])
    integral_value = np.trapezoid([np.exp(np.dot(R[1:], np.sqrt(2) * np.cos(2 * np.pi * np.arange(1, len(R)) * x))) for x in np.linspace(0, 1, 1000)], dx=1/1000)
    func_3 = - kap * np.log(integral_value)
    # print(f'func_3 with trapezoid: {func_3}')
    # integral_value, _ = integrate.quad(lambda x: np.exp(np.dot(R[1:], np.sqrt(2) * np.cos(2 * np.pi * np.arange(1, len(R)) * x))), 0, 1)
    # func_3 = - kap * np.log(integral_value)
    # print(f'func_3 with quad: {func_3}')
    return func_1 + func_2 + func_3

def minimize_functional(kappa, D):
    # Initial guess
    N = calc_argmax_N(kappa) + 15
    x0 = np.array([kappa] + [np.sqrt(2) * kappa / (1 + D * 4 * np.pi**2 * i**2) for i in range(1, N + 1)])

    # Perform the minimization
    def functional_dep_only_on_R(R):
        return functional(R, kappa, D)
    
    options = {'maxiter': 1000, 'disp': False}
    result = minimize(functional_dep_only_on_R, x0, method='L-BFGS-B', options=options)
    print(f'Minimization complete for kappa={kappa}, D={D}')

    return result.fun, result.x, result.fun < functional_dep_only_on_R(np.array([kappa]))

def main():
    # # Define constraints
    # def constraint1(x):
    #     # Example constraint: x[0] + x[1] = 3
    #     return x[0] + x[1] - 3

    # def constraint2(x):
    #     # Example constraint: x[0] - x[1] >= 0
    #     return x[0] - x[1]

    # # Initial guess
    # N = calc_argmax_N_v2(kappa) + 20
    # x0 = [kappa]
    # x0.extend([np.sqrt(2)*kappa / (1+ D*4*np.pi**2*i**2) for i in range(1, N+1)])

    # # Define constraints in the form of dictionaries
    # cons = [{'type': 'eq', 'fun': constraint1},
    #         {'type': 'ineq', 'fun': constraint2}]

    # # Perform the minimization
    # result = minimize(functional, x0) #, constraints=cons)

    kappa_test = 1
    D_test = max(calc_D_1(kappa_test), calc_D_2(kappa_test)) + 0.5e-5
    # If kappa = 0.5, then D = max(D_1, D_2) + 0.5e-5 is still working --> max(D_1, D_2) = D_2 not sharp
    # If kappa = 1, then D = max(D_1, D_2) + 0.5e-5 is still working --> max(D_1, D_2) = D_2 not sharp
    # If kappa = 1.1, then D = max(D_1, D_2) + 0.5e-5 = D_1 + 0.5e-5 is still working --> max(D_1, D_2) not sharp
    # If kappa = 1.5, then D = max(D_1, D_2) + 0.5e-6 is not working --> max(D_1, D_2)=D_1 sharp for kappa large enough?

    # If kappa = 1.002, then D_1 < D_2 and J(Phi) < J(kappa) for D_2 - 0.5e-5
    # If kappa = 1.003, then D_1 > D_2 and J(Phi) < J(kappa) for D_1 - 0.5e-5 > D_2

    print(f'D < D_1: {D_test} < {calc_D_1(kappa_test)}: {D_test < calc_D_1(kappa_test)}')
    print(f'D < D_2: {D_test} < {calc_D_2(kappa_test)}: {D_test < calc_D_2(kappa_test)}')

    result_fun, result_x, nonconst_minima_found = minimize_functional(kappa_test, D_test)
    print("Optimal value:", result_fun)
    print("Optimal solution:", result_x)

    print(f'J(min) < J(kappa) = {functional([kappa_test], kappa_test, D_test)}: {nonconst_minima_found}')

if __name__ == "__main__":
    main()
