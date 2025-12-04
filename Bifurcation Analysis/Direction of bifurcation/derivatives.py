import sympy as sp

from symbolic_functions import evaluate_symbolic_at_value_three_vars, get_constant_steady_state

##################################################################################
# Derivative functions
def compute_derivatives_three_vars(func):
    """
    Computes the first-order partial derivatives of a given function in three variables.

    Parameters:
    func : str
        The function as a string.

    Returns:
    tuple
        A tuple containing the partial derivatives with respect to x, y and z.
    """
    # Define the symbols
    x, y, z = sp.symbols('x y z')
    # Parse the function string into a SymPy expression
    f = sp.sympify(func)
    # Compute partial derivatives
    return sp.diff(f, x), sp.diff(f, y), sp.diff(f, z)


def calc_derivatives_of_nonlinearity_to_order_three_symb(func):
    fu, fv, fw = compute_derivatives_three_vars(func)

    fuu, fuv, fuw = compute_derivatives_three_vars(fu)
    _, fvv, fvw = compute_derivatives_three_vars(fv)
    _, _, fww = compute_derivatives_three_vars(fw)

    fuuu, fuuv, fuuw = compute_derivatives_three_vars(fuu)
    _, fuvv, fuvw = compute_derivatives_three_vars(fuv)
    _, _, fuww = compute_derivatives_three_vars(fuw)
    _, fvvv, fvvw = compute_derivatives_three_vars(fvv)
    _, _, fvww = compute_derivatives_three_vars(fvw)
    _, _, fwww = compute_derivatives_three_vars(fww)

    return fu, fv, fw, fuu, fuv, fuw, fvv, fvw, fww, fuuu, fuuv, fuuw, fuvv, fuvw, fuww, fvvv, fvvw, fvww, fwww

def get_first_derivs_eval(model_func_f, model_func_g, model_func_h, u_bar, v_bar, w_bar):
    # u_bar, v_bar, w_bar = get_constant_steady_state(model_func_f, model_func_g, model_func_h, initial_guess)
    f_derv_symb = compute_derivatives_three_vars(model_func_f)
    g_derv_symb = compute_derivatives_three_vars(model_func_g)
    h_derv_symb = compute_derivatives_three_vars(model_func_h)
    fu, fv, fw = [evaluate_symbolic_at_value_three_vars(derv, u_bar, v_bar, w_bar) for derv in f_derv_symb]
    gu, gv, gw = [evaluate_symbolic_at_value_three_vars(derv, u_bar, v_bar, w_bar) for derv in g_derv_symb]
    hu, hv, hw = [evaluate_symbolic_at_value_three_vars(derv, u_bar, v_bar, w_bar) for derv in h_derv_symb]
    return fu, fv, fw, gu, gv, gw, hu, hv, hw

def get_first_derivs_at_steady_state(model_func_f, model_func_g, model_func_h, initial_guess):
    u_bar, v_bar, w_bar = get_constant_steady_state(model_func_f, model_func_g, model_func_h, initial_guess)
    fu, fv, fw, gu, gv, gw, hu, hv, hw = get_first_derivs_eval(model_func_f, model_func_g, model_func_h, u_bar, v_bar, w_bar)
    return fu, fv, fw, gu, gv, gw, hu, hv, hw


def calc_derivatives_of_p(fu, fv, fw, fuu, fuv, fuw, fvv, fvw, fww, fuuu, fuuv, fuuw, fuvv, fuvw, fuww, fvvv,
                          fvvw, fvww, fwww):
    pv = - fv / fu
    pw = - fw / fu
    pvv = - (fuu*pv**2 + 2*fuv*pv + fvv) / fu
    pvw = - (fuu*pv*pw + fuw*pv + fuv*pw + fvw) / fu
    pww = - (fuu*pw**2 + 2*fuw*pw + fww) / fu
    pvvv = - (fuuu*pv**3 + 3*fuuv*pv**2 + 3*fuu*pv*pvv + 3*fuvv*pv + 3*fuv*pvv + fvvv) / fu
    pvvw = - (fuuu*pv**2*pw + fuuw*pv**2 + fuu*(2*pv*pvw+pvv*pw) + 2*fuuv*pv*pw + 2*fuvw*pv + 2*fuv*pvw + fuw*pvv +
                fuvv*pw + fvvw) / fu
    pvww = - (fuuu*pv*pw**2 + 2*fuuw*pv*pw + fuu*(2*pvw*pw+pv*pww) + fuww*pv + 2*fuw*pvw + fuuv*pw**2 + 2*fuvw*pw +
                fuv*pww + fvww) / fu
    pwww = - (fuuu*pw**3 + 3*fuuw*pw**2 + 3*fuu*pw*pww + 3*fuww*pw + 3*fuw*pww + fwww) / fu

    derivatives_p = [pv, pw, pvv, pvw, pww, pvvv, pvvw, pvww, pwww]
    return derivatives_p


def calc_hat_g_first_derv(pv, pw, gu, gv, gw):
    hat_gv = gu*pv + gv   # gv - gu * fv / fu
    hat_gw = gu*pw + gw   # gw - gu * fw / fu
    return hat_gv, hat_gw

def calc_hat_g_sec_derv(pv, pw, pvv, pvw, pww, gu, guu, guv, guw, gvv, gvw, gww):
    hat_gvv = guu*pv**2 + 2*guv*pv + gu*pvv + gvv
    hat_gvw = guu*pv*pw + guw*pv + gu*pvw + guv*pw + gvw
    hat_gww = guu*pw**2 + 2*guw*pw + gu*pww + gww
    return hat_gvv, hat_gvw, hat_gww

def calc_hat_g_third_derv(pv, pw, pvv, pvw, pww, pvvv, pvvw, pvww, pwww, gu, guu, guv, guw, guuu, guuv, guuw, guvv, guvw, guww, gvvv, gvvw, gvww, gwww):
    hat_gvvv = guuu*pv**3 + 3*guuv*pv**2 + 3*guu*pv*pvv + 3*guvv*pv + 3*guv*pvv + gu*pvvv + gvvv
    hat_gvvw = (guuu*pv**2*pw + guuw*pv**2 + guu*(2*pv*pvw+pvv*pw) + 2*guuv*pv*pw + 2*guvw*pv + 2*guv*pvw + guw*pvv
                + gu*pvvw + guvv*pw + gvvw)
    hat_gvww = (guuu*pv*pw**2 + 2*guuw*pv*pw + guu*(2*pvw*pw+pv*pww) + guww*pv + 2*guw*pvw + gu*pvww + guuv*pw**2
                + 2*guvw*pw + guv*pww + gvww)
    hat_gwww = guuu*pw**3 + 3*guuw*pw**2 + 3*guu*pw*pww + 3*guww*pw + 3*guw*pww + gu*pwww + gwww
    return hat_gvvv, hat_gvvw, hat_gvww, hat_gwww

def calc_all_hat_g_dervs(pv, pw, pvv, pvw, pww, pvvv, pvvw, pvww, pwww, gu, gv, gw, guu, guv, guw, gvv, gvw, gww, guuu, guuv, guuw, guvv, guvw, guww, gvvv, gvvw, gvww, gwww):
    derivs = []
    derivs.extend(calc_hat_g_first_derv(pv, pw, gu, gv, gw))
    derivs.extend(calc_hat_g_sec_derv(pv, pw, pvv, pvw, pww, gu, guu, guv, guw, gvv, gvw, gww))
    derivs.extend(calc_hat_g_third_derv(pv, pw, pvv, pvw, pww, pvvv, pvvw, pvww, pwww, gu, guu, guv, guw, guuu, guuv, guuw, guvv, guvw, guww, gvvv, gvvw, gvww, gwww))
    return derivs


def get_all_hat_derivatives(nonlin_f, nonlin_g, nonlin_h, u_bar, v_bar, w_bar):
    derivs_f_symb = calc_derivatives_of_nonlinearity_to_order_three_symb(nonlin_f)
    derivs_f = [evaluate_symbolic_at_value_three_vars(derv, u_bar, v_bar, w_bar) for derv in derivs_f_symb]
    derivs_g_symb = calc_derivatives_of_nonlinearity_to_order_three_symb(nonlin_g)
    derivs_g = [evaluate_symbolic_at_value_three_vars(derv, u_bar, v_bar, w_bar) for derv in derivs_g_symb]
    derivs_h_symb = calc_derivatives_of_nonlinearity_to_order_three_symb(nonlin_h)
    derivs_h = [evaluate_symbolic_at_value_three_vars(derv, u_bar, v_bar, w_bar) for derv in derivs_h_symb]

    derivatives_p = calc_derivatives_of_p(*derivs_f)

    derivs_hat_g = calc_all_hat_g_dervs(*derivatives_p, *derivs_g)
    derivs_hat_h = calc_all_hat_g_dervs(*derivatives_p, *derivs_h)

    return derivs_hat_g, derivs_hat_h

def get_all_hat_derivatives_at_steady_state(nonlin_f, nonlin_g, nonlin_h, initial_guess):
    u_bar, v_bar, w_bar = get_constant_steady_state(nonlin_f, nonlin_g, nonlin_h, initial_guess)
    derivs_hat_g, derivs_hat_h = get_all_hat_derivatives(nonlin_f,nonlin_g, nonlin_h, u_bar, v_bar, w_bar)
    return derivs_hat_g, derivs_hat_h



# def get_derivs_at_constant_steady_state(self, func):
#     """
#     Evaluates all derivatives of function func at the constant steady state.

#     Returns:
#     list
#         A list containing the evaluated derivatives.
#     """

#     derivatives = calc_derivatives_of_nonlinearity_to_order_three_symb(func)
#     return [evaluate_symbolic_at_value_three_vars(derv, self.u, self.v, self.w) for derv in derivatives]