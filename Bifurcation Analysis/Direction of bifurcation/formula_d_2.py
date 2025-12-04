import numpy as np
import scipy.integrate as integrate


###################################################################################
# Calculate d_2 in d(s) = d_hat + d_2*s**2 + ...
def calculate_d_2(fu, fv, fuu, fuv, fvv, fuuu, fuuv, fuvv, fvvv, gu, gv, guu, guv, gvv, guuu, guuv, guvv, gvvv, d_hat,
                  D_2, L, l_j):
    # We have branch of non-constant steady states for diffusion coefficient d(s) = d_hat + d_2*s^2 + d_3*s^3 + ...
    # Calculate d_2 in order to approximate behavior of branch.
    const_f = _constant_f(gu, gv, D_2, l_j)
    const_f_bar = _constant_f_bar(fv, gv, D_2, l_j)

    A_1, A_2, A_3 = _calc_all_A_i(fuu, fuv, fvv, fuuu, fuuv, fuvv, fvvv, guu, guv, gvv, guuu, guuv, guvv, gvvv, const_f,
                                  const_f_bar)

    eta, xi = _calc_eta_xi(fu, fv, fuu, fuv, fvv, gu, gv, guu, guv, gvv, const_f, d_hat, D_2, l_j)

    integral_1 = integrate.quad(
        lambda x: _P_2(eta, xi, l_j, x)[0] * np.cos(np.sqrt(l_j) * x) ** 2, 0, L)[0]
    integral_2 = integrate.quad(
        lambda x: _P_2(eta, xi, l_j, x)[1] * np.cos(np.sqrt(l_j) * x) ** 2, 0, L)[0]

    d_2 = 2 / (l_j * L) * (A_1 * integral_1 + A_2 * integral_2 + (3 * L) / 8 * A_3)
    return d_2

# Calculate d_2 in d(s) = d_hat + d_2*s**2 + ..., hopefully faster version
def calculate_d_2_faster(fu, fv, fuu, fuv, fvv, fuuu, fuuv, fuvv, fvvv, gu, gv, guu, guv, gvv, guuu, guuv, guvv, gvvv, d_hat,
                         D_2, L, l_j):
    # We have branch of non-constant steady states for diffusion coefficient d(s) = d_hat + d_2*s^2 + d_3*s^3 + ...
    # Calculate d_2 in order to approximate behavior of branch.
    const_f = _constant_f(gu, gv, D_2, l_j)
    const_f_bar = _constant_f_bar(fv, gv, D_2, l_j)

    A_1, A_2, A_3 = _calc_all_A_i(fuu, fuv, fvv, fuuu, fuuv, fuvv, fvvv, guu, guv, gvv, guuu, guuv, guvv, gvvv, const_f,
                                  const_f_bar)

    eta, xi = _calc_eta_xi(fu, fv, fuu, fuv, fvv, gu, gv, guu, guv, gvv, const_f, d_hat, D_2, l_j)

    l_j_root = np.emath.sqrt(l_j)
    # It holds j*pi = sqrt{l_j}*L. If j natural number, then np.sin(2*l_j_root*L) = np.sin(4*l_j_root*L) = 0, i.e. part_2 = 0
    part_1 = 3/4 * A_3 + A_1/2 * (2*eta[0] + xi[0]) + A_2/2 * (2*eta[1] + xi[1])
    part_2 = (4*A_1*(eta[0]+xi[0]) + 4*A_2*(eta[1]+xi[1])) * np.sin(2*l_j_root*L) + (A_1*xi[0] + A_2*xi[1]) * np.sin(4*l_j_root*L)

    d_2 = 1/l_j * part_1 + 1 / (8*L*l_j*l_j_root) * part_2
    return d_2

def calculate_limit_d_2(fu, fv, fuu, fuv, fvv, fuuu, fuuv, fuvv, fvvv, gu, gv, guu, guv, gvv, guuu, guuv, guvv, gvvv, d_max, D_2, l_max):
    c_f = _constant_f(gu, gv, D_2, l_max)
    c_f_bar = _constant_f_bar(fv, gv, D_2, l_max)
    A_1, A_2, A_3 = _calc_all_A_i(fuu, fuv, fvv, fuuu, fuuv, fuvv, fvvv, guu, guv, gvv, guuu, guuv, guvv, gvvv, c_f, c_f_bar)
    eta, xi = _calc_eta_xi(fu, fv, fuu, fuv, fvv, gu, gv, guu, guv, gvv, c_f, d_max, D_2, l_max)
    
    d_2 = A_1/(2*l_max) * (2*eta[0] + xi[0]) + A_2/(2*l_max) * (2*eta[1] + xi[1]) + (3*A_3)/(4*l_max)

    return d_2

def _constant_f(gu, gv, D_2, l_j):
    # Help function to calculate d_2 in d(s) = d_hat + d_2*s^2 + d_3*s^3 + ...
    return gu / (D_2 * l_j - gv)


def _constant_f_bar(fv, gv, D_2, l_j):
    # Help function to calculate d_2 in d(s) = d_hat + d_2*s^2 + d_3*s^3 + ...
    return fv / (D_2 * l_j - gv)


def _calc_A_1(fuu, fuv, guu, guv, f, f_bar):
    # Help function to calculate d_2 in d(s) = d_hat + d_2*s^2 + d_3*s^3 + ...
    return fuu + fuv * f + guu * f_bar + guv * f * f_bar


def _calc_A_2(fuv, fvv, guv, gvv, f, f_bar):
    # Help function to calculate d_2 in d(s) = d_hat + d_2*s^2 + d_3*s^3 + ...
    return fuv + fvv * f + guv * f_bar + gvv * f_bar * f


def _calc_A_3(fuuu, fuuv, fuvv, fvvv, guuu, guuv, guvv, gvvv, f, f_bar):
    # Help function to calculate d_2 in d(s) = d_hat + d_2*s^2 + d_3*s^3 + ...
    A_3_1 = 1 / 6 * fuuu + 1 / 2 * fuuv * f + 1 / 2 * fuvv * f ** 2 + 1 / 6 * fvvv * f ** 3
    A_3_2 = (1 / 6 * guuu * f_bar + 1 / 2 * guuv * f * f_bar + 1 / 2 * guvv * f ** 2 * f_bar
             + 1 / 6 * gvvv * f ** 3 * f_bar)
    return A_3_1 + A_3_2


def _calc_all_A_i(fuu, fuv, fvv, fuuu, fuuv, fuvv, fvvv, guu, guv, gvv, guuu, guuv, guvv, gvvv, f, f_bar):
    # Help function to calculate d_2 in d(s) = d_hat + d_2*s^2 + d_3*s^3 + ...
    A_1 = _calc_A_1(fuu, fuv, guu, guv, f, f_bar)
    A_2 = _calc_A_2(fuv, fvv, guv, gvv, f, f_bar)
    A_3 = _calc_A_3(fuuu, fuuv, fuvv, fvvv, guuu, guuv, guvv, gvvv, f, f_bar)
    return A_1, A_2, A_3


def _calc_eta_xi(fu, fv, fuu, fuv, fvv, gu, gv, guu, guv, gvv, const_f, d_hat, D_2, l_j):
    # Help function to calculate d_2 in d(s) = d_hat + d_2*s^2 + d_3*s^3 + ...
    h = np.array([1 / 2 * (1 / 2 * fuu + fuv * const_f + 1 / 2 * fvv * const_f ** 2),
                  1 / 2 * (1 / 2 * guu + guv * const_f + 1 / 2 * gvv * const_f ** 2)])
    L_0 = np.array([[fu, fv], [gu, gv]])
    L_2 = np.array([[fu - 4 * d_hat * l_j, fv], [gu, gv - 4 * D_2 * l_j]])

    # def eta_first_equ(y, z):
    #     return L_0[0][0] * y + L_0[0][1] * z + h[0]
    #
    # def eta_sec_equ(y, z):
    #     return L_0[1][0] * y + L_0[1][1] * z + h[1]
    #
    # eta = solve_system_of_equations_numerically([eta_first_equ, eta_sec_equ], [0, 0])
    eta = np.linalg.solve(L_0, -h)

    # def xi_first_equ(y, z):
    #     return L_2[0][0] * y + L_2[0][1] * z + h[0]
    #
    # def xi_sec_equ(y, z):
    #     return L_2[1][0] * y + L_2[1][1] * z + h[1]
    #
    # xi = solve_system_of_equations_numerically([xi_first_equ, xi_sec_equ], [0, 0])
    xi = np.linalg.solve(L_2, -h)
    return eta, xi


def _P_2(eta, xi, l_j, x):
    u_2 = eta[0] + xi[0] * np.cos(2 * np.sqrt(l_j) * x)
    v_2 = eta[1] + xi[1] * np.cos(2 * np.sqrt(l_j) * x)
    return u_2, v_2

