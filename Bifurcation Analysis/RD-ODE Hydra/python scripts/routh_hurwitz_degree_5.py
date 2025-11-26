import numpy as np
import sympy as sp
from sympy import Matrix
from sympy.abc import x, y
import matplotlib.pyplot as plt
# from scipy.optimize import fsolve


def ruth_hurwitz_degree_5(p):
    """Check stability of a degree 5 polynomial using the Routh-Hurwitz criterion.

    Args:
        p (list): Coefficients of the polynomial in decreasing order of degree.

    Returns:
        list: First column of the Routh array.
    """
    if len(p) != 6:
        raise ValueError("Polynomial must be of degree 5 (6 coefficients).")

    # Construct the Routh array
    routh_array = [[0 for _ in range(3)] for _ in range(6)]
    routh_array[0] = [p[0], p[2], p[4]]
    routh_array[1] = [p[1], p[3], p[5]]

    for i in range(2, 6):
        for j in range(2):
            numerator = (routh_array[i-1][0] * routh_array[i-2][j+1] -
                         routh_array[i-2][0] * routh_array[i-1][j+1])
            denominator = routh_array[i-1][0]
            if denominator == 0:
                routh_array[i][j] = 0
            else:
                routh_array[i][j] = numerator / denominator

    # print("Routh Array:")
    # for row in routh_array:
    #     print(row)
    # The first column of the Routh array gives the necessary conditions for stability
    stab_conditions = []
    for i in range(6):
        stab_conditions.append(routh_array[i][0])
    # print(f"Ruth-Hurwitz conditions {stab_conditions}")
    return stab_conditions

def charac_poly(Jac, D, wavenumber):
    M = Jac - wavenumber**2 * sp.pi**2 * D
    # poly = M.charpoly(x).as_expr()
    # poly_coef = sp.Poly(poly, x).all_coeffs()
    poly_coef = M.charpoly(x).all_coeffs()
    poly_coef = [sp.simplify(coef) for coef in poly_coef]
    return poly_coef

def ruth_hurwitz_conds_matrix(Jac, D, wavenumber):
    poly_coef = charac_poly(Jac, D, wavenumber)
    stab_conds = ruth_hurwitz_degree_5(poly_coef)
    stab_conds = [sp.simplify(cond) for cond in stab_conds]
    return stab_conds

def get_point_zero_eigv(Jac: Matrix, D: Matrix, wavenumber: int, d: sp.Symbol):
    M = Jac - wavenumber**2 * sp.pi**2 * D
    det = M.det()
    root = sp.poly(det, d, domain='RR').nroots()
    return root[0].evalf()

def get_bifurcation_point(Jac: Matrix, D: Matrix, wavenumber: int, d: sp.Symbol, fast: bool = True) -> float:
    stab_conds = ruth_hurwitz_conds_matrix(Jac, D, wavenumber)

    candidates = []
    for i, cond in enumerate(stab_conds):
        # print(f"Condition {i}: {cond}")
        if i == 0:
            continue  # skip the first condition (always positive for this system)
        # print("")
        # print(f"Condition {i}: {cond} > 0")

        # def func(d_val):
        #     return cond.subs(d, d_val)
        # func_lam = sp.lambdify(d, func(d))
        # cand_fsolve = fsolve(func_lam, 0.008)
        # candidates.extend([c for c in cand_fsolve if np.imag(c) == 0])

        # cond = cond.evalf()
        num, den = sp.fraction(sp.cancel(sp.simplify(cond)))
        try:
            num_roots = sp.poly(num, d, domain='RR').nroots()
            num_roots = [np.complex128(root.evalf()) for root in num_roots]
            den_roots = sp.poly(den, d, domain='RR').nroots()
            den_roots = [np.complex128(root.evalf()) for root in den_roots]
        except:
            num_coef = sp.Poly(num, d).all_coeffs()
            den_coef = sp.Poly(den, d).all_coeffs()
            num_roots = np.roots([float(coef.evalf()) for coef in num_coef])
            den_roots = np.roots([float(coef.evalf()) for coef in den_coef])
            # num_roots = np.roots(sp.poly(num, d).all_coeffs())
            # den_roots = np.roots(sp.poly(den, d).all_coeffs())
        
        cond_roots = [root for root in num_roots if all(abs(root - s) > 1e-10 for s in den_roots)]
        # print(f"Zeros: {cond_roots}")
        candidates.extend([np.real(c) for c in cond_roots if np.imag(c) == 0]) # c.is_real

        # print(f"Candidates for bifurcation point from this condition: {cond_roots}")
    # print("")
    # print(f"All candidates for bifurcation point: {candidates}")
    if fast:
        if [c for c in candidates if c > 0 and c < 1e5]:
            bif_point = min([c for c in candidates if c > 0 and c < 1e5])
        else:
            bif_point = np.nan
    else:
        bif_point = max([c for c in candidates if c < 1e5])
    print(f"{"Minimum" if fast else "Maximum"} positive bifurcation point for wavenumber {wavenumber}: {bif_point}")

    return bif_point

def get_bif_points_fast(J, D, max_wn = 200):
    # d = sp.Symbol("d")
    wn = sp.Symbol("wn")
    Mw = J - wn**2 * sp.pi**2 * D
    detw = Mw.det()
    dw = sp.solve(detw, d)[0]
    # print(dw)
    bif_points = []
    wave_numbers = []
    for k in range(1, max_wn + 1):
        dk = dw.subs(wn, k)
        if dk > 0:
            bif_points.append(dk)
            wave_numbers.append(k)
        else:
            break
    return wave_numbers, bif_points


if __name__ == '__main__':
    d = sp.Symbol("d")
    U0 = [0.08167547924306925, 0.54587711524379, 3.6049476660047937, 1.8514050545898464, 0.08167547924306925]
    beta = [1.0629, 540.4003, 1.1596, 11.5964, 11.5964, 4.8254]
    Wl, A, Wd, C, S = U0
    Jac_allg = [
        [- beta[2]*beta[5]*S / ((1+A)*(1+C)*(1+beta[2]*Wl)**2) - 1.0, -beta[5]*S / ((1+A)**2*(1+C)*(1+beta[2]*Wl)), 0.0, -beta[5]*S / ((1+A)*(1+C)**2*(1+beta[2]*Wl)), beta[5] / ((1+A)*(1+C)*(1+beta[2]*Wl))],
        [-beta[0]*beta[3] / (1+beta[3]*Wl)**2, -1.0, 0.0, 0.0, 0.0],
        [beta[1]*S, 0.0, -1.0, 0.0, beta[1]*Wl],
        [-beta[4]*Wd / (1+beta[4]*Wl)**2, 0.0, 1 / (1+beta[4]*Wl), -1.0, 0.0],
        [1.0, 0.0, 0.0, 0.0, -1.0]
        ]
    Jac = Matrix([
        [-1.0865168027146117, -0.05283439313362805, 0.0, -0.028643941383072864, 0.9999999999999997],
        [-3.2510268461888066, -1.0, 0.0, 0.0, 0.0],
        [44.1374534855984, 0.0, -1.0, 0.0, 44.1374534855984],
        [-11.026231669288988, 0.0, 0.5135733514383242, -1.0, 0.0],
        [1.0, 0.0, 0.0, 0.0, -1.0]
        ])
    
    # u0_variant = [0.00260781026235803, 1.03170013849885, 1.4092614481213583, 1.36789465726587, 0.00260781026235803]
    # Jac_variant = Matrix([[-1.0030148996730284, -0.0012835606066773457, 0.0, -0.0011013202189362518, 1.000000000000001],
    #            [-11.61282169573653, -1.0, 0.0, 0.0, 0.0],
    #            [270.20015, 0.0, -1.0, 0.0, 270.20015],
    #            [-15.397028807704197, 0.0, 0.9706464752082521, -1.0, 0.0],
    #            [1.0, 0.0, 0.0, 0.0, -1.0]])
    
    # Jac_variant = Matrix([[-1.54835418090850, -0.968673468738539, 0, -0.519279522931973, 1.00000000000000], [-0.0713698848523722, -1, 0, 0, 0], [1.15960000000000, 0, -1, 0, 11.5964000000000], [-0.896792686154216, 0, 0.0760938750398222, -1, 0], [1, 0, 0, 0, -1]])
    # Jac_variant_2 = Matrix([[-1.00301489967303, -0.00128356060667735, 0, -0.00110132021893625, 1.00000000000000], [-11.6128216957365, -1, 0, 0, 0], [540.400300000000, 0, -1, 0, 0], [-15.3970288077042, 0, 0.970646475208252, -1, 0], [1, 0, 0, 0, -1]])

    # Jac = Matrix([
    #     [-1.08652, -0.0528344,  0.0,      -0.0286439,  1.0],
    #     [-3.25103, -1.0,        0.0,       0.0,        0.0],
    #     [44.1375,   0.0,       -1.0,       0.0,        44.1375],
    #     [-11.0262,  0.0,        0.513573, -1.0,        0.0],
    #     [1.0,       0.0,        0.0,       0.0,       -1.0]
    # ])
    # nu = [0.0, 3.8154e-05, 0.4433, 6.0713e-08, 0.0004]
    nu = [0.0, 3.8154e-05, 0.4433, 6.0713e-08, 0.0004]
    # D = Matrix.diag(0.0, d, 1e-5, 10.0, 1e-6)
    D = Matrix.diag(0.0, 3.8154e-05, d, 6.0713e-08, 0.0004)   

    # list_of_wavenumbers = range(1, 13)
    # plt.plot(list_of_wavenumbers, [get_bifurcation_point(Jac, D, k, d, fast=True) for k in list_of_wavenumbers], 'bo', markersize=4)
    # plt.xlabel("Wavenumber")
    # plt.ylabel("Bifurcation point")
    # plt.title("Bifurcation points vs Wavenumber")
    # plt.grid()
    # plt.show()
     
    for k in range(1, 13):
        bif_point = get_bifurcation_point(Jac, D, k, d, fast=True)
        print(f"{k}: {bif_point}")
        # print(f"{dw.subs(wn, k)}")
        