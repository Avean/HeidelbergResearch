import numpy as np
from numpy import linalg as LA

from surpress_print import suppress_print
from properties_bifurcating_branch import TuringBifurcationThreeEqu
from derivatives import get_first_derivs_at_steady_state, get_all_hat_derivatives_at_steady_state
from formula_d_2 import calculate_d_2, calculate_d_2_faster

# These functions help if you want to calculate max bifurcation point, corr. d_2 and other stuff for changing domain size
# Also some functions only relay on Jacobian instead of nonlinearities

def calc_d_2_at_max_bif_var_domain_size(function_f, function_g, function_h, D_2, domain_length, initial_guess_steady_state):
    with suppress_print():
        turing_bif = TuringBifurcationThreeEqu(function_f, function_g, function_h, D_2, domain_length, initial_guess_steady_state)
    return turing_bif.calc_d_2_at_max_bif()

def calc_det_L_2_at_max_bif(function_f, function_g, function_h, D_2, domain_length, initial_guess_steady_state):
    if type(domain_length) is np.ndarray:
        domain_length = domain_length[0]
    with suppress_print():
        turing_bif = TuringBifurcationThreeEqu(function_f, function_g, function_h, D_2, domain_length, initial_guess_steady_state)
        max_wave_number, max_bif_point = turing_bif.search_max_bifurcation_point()
    
    # p function such that f(p(v,w), v, w) = 0 around constant steady state.
    # Such a function exists due to fu \ne 0 and implicit function theorem.
    # Calculate derivative of hat_g(v,w) = g(p(v,w),v,w) and hat_h(v,w) = h(p(v,w),v,w).
    fu, fv, fw, gu, gv, gw, hu, hv, hw = get_first_derivs_at_steady_state(function_f, function_g, function_h, initial_guess_steady_state)
    return _calc_det_L_2(fu, fv, fw, gu, gv, gw, hu, hv, hw, D_2, domain_length, max_wave_number, max_bif_point)

def approx_max_wave_number_from_Jac_p3(fu, fv, fw, gu, gv, gw, hu, hv, hw, D_2, domain_length, rounded: bool, both: bool=False):
    func_f = f"{fu}*x + {fv}*y - {fw}*z"
    func_g = f"{gu}*x + {gv}*y - {gw}*z"
    func_h = f"{hu}*x + {hv}*y - {hw}*z"
    return TuringBifurcationThreeEqu(func_f, func_g, func_h, D_2, domain_length, [0,0,0]).approx_max_wave_number_p3(rounded=rounded, both=both)

def det_L_2_direct_from_Jacobian(fu, fv, fw, gu, gv, gw, hu, hv, hw, D_2, domain_length, rounded_wave_number=True):
    j_max, d_max = approx_max_wave_number_from_Jac_p3(fu, fv, fw, gu, gv, gw, hu, hv, hw, D_2, domain_length, rounded_wave_number)
    return _calc_det_L_2(fu, fv, fw, gu, gv, gw, hu, hv, hw, D_2, domain_length, j_max, d_max)

def get_det_L_2_limit(fu, fv, fw, gu, gv, gw, hu, hv, hw, D_2):
    domain_length = 1000

    det_L_2 = det_L_2_direct_from_Jacobian(fu, fv, fw, gu, gv, gw, hu, hv, hw, D_2, domain_length, rounded_wave_number=False)
    
    return det_L_2

def d_2_direct_from_nonlinearities_approx(model_func_f, model_func_g, model_func_h, D_2, domain_length, initial_guess_steady_state, faster=True, rdd_wave_number=True):
    # Not really limit, only for a subsequence. For det L_2 it worked much better, so here some new dependencies on domain_size but one subsequence more or less uneffected.
    turing_bif = TuringBifurcationThreeEqu(model_func_f, model_func_g, model_func_h, D_2, domain_length, initial_guess_steady_state)
    j_max, d_max = turing_bif.approx_max_wave_number_p3(rounded=rdd_wave_number, both=False)
    l_max = j_max**2 * np.pi ** 2 / domain_length ** 2

    derivs_hat_g, derivs_hat_h = get_all_hat_derivatives_at_steady_state(model_func_f, model_func_g, model_func_h, initial_guess_steady_state)
    if faster:
        d_2 = calculate_d_2_faster(*derivs_hat_g, *derivs_hat_h, d_max, D_2, domain_length, l_max)
    else:
        d_2 = calculate_d_2(*derivs_hat_g, *derivs_hat_h, d_max, D_2, domain_length, l_max)

    return d_2

def _calc_det_L_2(fu, fv, fw, gu, gv, gw, hu, hv, hw, D_2, domain_length, max_wave_number, max_bif_point):
    hat_gv = gv - gu * fv / fu
    hat_gw = gw - gu * fw / fu
    hat_hv = hv - hu * fv / fu
    hat_hw = hw - hu * fw / fu

    l_j = max_wave_number ** 2 * np.pi ** 2 / domain_length ** 2
    detL2 = (hat_gv - 4 * max_bif_point * l_j) * (hat_hw - 4 * D_2 * l_j) - hat_gw * hat_hv
    return detL2
