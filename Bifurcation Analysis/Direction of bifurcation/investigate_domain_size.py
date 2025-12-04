import numpy as np
import matplotlib.pyplot as plt

from surpress_print import suppress_print
from properties_bifurcating_branch import TuringBifurcationThreeEqu

# Try to find domain size such that maximal bifurcation point satisfies suff conditions for Turing patterns
def find_domain_size_max_bif_satisfies_suff_conds(function_f, function_g, function_h, D_2, initial_guess_steady_state): 
    with suppress_print():
        possible_turing_patterns = dict()
        if TuringBifurcationThreeEqu(function_f, function_g, function_h, D_2, 1, initial_guess_steady_state).fu < 0:
            for length in np.linspace(100, 1000, 451):
                turing_bif = TuringBifurcationThreeEqu(function_f, function_g, function_h, D_2, length, initial_guess_steady_state)
                d_max_ind, d_max = turing_bif.search_max_bifurcation_point()
                if d_max_ind != 0:
                    d_max, d_2_branch = turing_bif.turing_bifurcation(d_max_ind)
                else:
                    d_2_branch = 0.0

                if d_max > 0 > d_2_branch:
                    # In example model, we have direction DDI = -1 as J_12 unstable and D_1/D_2 sufficiently small.
                    # Hence, maximal bifurcation point d_max positive and d_2 from d(s) = d_max + d_2 * s^2 + ... negative is sufficient for Turing patterns if fu<0.
                    candid_max_bif = (d_max_ind, d_max, d_2_branch)
                    possible_turing_patterns[str(length)] = candid_max_bif
    print('\n\n')
    print('Possible (safe?) candidates of domain size for Turing patterns:')
    print(possible_turing_patterns)



# Consider model with Turing patterns, i.e. for max. bif point the branch is in right direction for some domain size.
# Try to find domain size such that branch is in wrong direction for the maximal bif. point.
# Or more general try to find sign change of d_2 corresponding to maximal bifurcation point by changing domain size.
# For any bif point possible, see L = 13 and Ext. Brusselator but also for maximal bif point? 
# Wrong: Found it, e.g. domain length >= 1034.05 in comparison to <=1034.04. Nope!!! Just didn't far enough to calculate max. bif. point!!!

# def plot_domain_size_lambda_max(domain_size_list):
#     max_corr_eigenvalue_rdd = [approx_max_wave_number(fu, fv, fw, gu, gv, gw, hu, hv, hw, Diff_2, length, True)[0]**2 * np.pi**2 / length**2 for length in domain_size_list]
#     max_corr_eigenvalue = [approx_max_wave_number(fu, fv, fw, gu, gv, gw, hu, hv, hw, Diff_2, length, False)[0]**2 * np.pi**2 / length**2 for length in domain_size_list]
#     with suppress_print():
#         max_wave_number_exact = [TuringBifurcationThreeEqu(func_f, func_g, func_h, Diff_2, length, initial_guess_steady_state).calc_max_bifurcation_point()[0][0]**2 * np.pi**2 / length**2 for length in domain_size_list]
#     plt.plot(domain_size_list, max_corr_eigenvalue, 'g*', ms=10, label='Lapl. eigv to approx max wave number')
#     plt.plot(domain_size_list, max_corr_eigenvalue_rdd, 'r*', ms=8, label='Lapl. eigv to rounded approx max wave number')
#     plt.plot(domain_size_list, max_wave_number_exact, 'b*', label='Lapl. eigv to exact max wave number')
#     plt.xlabel('domain size')
#     plt.ylabel('lambda_max')
#     plt.legend(loc="upper right")
#     plt.figure()

# def plot_domain_size_d_max(domain_size_list):
#     max_bif_point_rdd = [approx_max_wave_number(fu, fv, fw, gu, gv, gw, hu, hv, hw, Diff_2, length, True)[1] for length in domain_size_list]
#     max_bif_point = [approx_max_wave_number(fu, fv, fw, gu, gv, gw, hu, hv, hw, Diff_2, length, False)[1] for length in domain_size_list]
#     with suppress_print():
#         max_bif_point_exact = [TuringBifurcationThreeEqu(func_f, func_g, func_h, Diff_2, length, initial_guess_steady_state).calc_max_bifurcation_point()[0][1] for length in domain_size_list]
#     plt.plot(domain_size_list, max_bif_point, 'g*', ms=10, label='max bif point approx.')
#     plt.plot(domain_size_list, max_bif_point_rdd, 'r*', ms=8, label='max bif point approx with rounded wave number')
#     plt.plot(domain_size_list, max_bif_point_exact, 'b*', label='max bif point exact')
#     plt.xlabel('domain size')
#     plt.ylabel('d_max')
#     plt.legend()
#     plt.figure()

# def plot_domain_size_d_2_corr_dmax(domain_size_list):
#     plt.plot(domain_size_list, [d_2_direct_from_nonlin_without_rdd(func_f, func_g, func_h, Diff_2, domain_length=length, initial_guess_steady_state=initial_guess_steady_state, faster=False) for length in domain_size_list], 'g*', ms=10, label='d_2 at max bif point approx.')
#     plt.plot(domain_size_list, [d_2_direct_from_nonlin_without_rdd(func_f, func_g, func_h, Diff_2, domain_length=length, initial_guess_steady_state=initial_guess_steady_state, faster=True) for length in domain_size_list], 'g*', label='d_2 at max bif point approx. faster')
#     plt.plot(domain_size_list, [d_2_direct_from_nonlinearities_approx(func_f, func_g, func_h, Diff_2, length, initial_guess_steady_state) for length in domain_size_list], 'r*', ms=8, label='d_2 at max bif point approx with rdd wave number')
#     plt.plot(domain_size_list, [calc_d_2_at_max_bif(func_f, func_g, func_h, Diff_2, length, initial_guess_steady_state) for length in domain_size_list], 'b*', label='d_2 at max bif point exact')
#     plt.xlabel('domain size')
#     plt.ylabel('d_2 corr. to d_max')
#     plt.legend()
#     plt.figure()

# def plot_domain_size_A_i_rdd_wave_number(domain_size_list):
#     plt.plot(domain_size_list, [_calc_all_A_i(hat_guu, hat_guv, hat_gvv, hat_guuu, hat_guuv, hat_guvv, hat_gvvv, hat_huu, hat_huv, hat_hvv, hat_huuu, hat_huuv, hat_huvv, hat_hvvv, _constant_f(hat_hu, hat_hv, Diff_2, approx_max_wave_number(fu, fv, fw, gu, gv, gw, hu, hv, hw, Diff_2, length, True)[0]**2 * np.pi**2 / length**2),
#                               _constant_f_bar(hat_gv, hat_hv, Diff_2, approx_max_wave_number(fu, fv, fw, gu, gv, gw, hu, hv, hw, Diff_2, length, True)[0]**2 * np.pi**2 / length**2))[0] for length in domain_size_list], 'g*', label='A_1 at max bif point approx. with rdd wave number')
#     plt.plot(domain_size_list, [_calc_all_A_i(hat_guu, hat_guv, hat_gvv, hat_guuu, hat_guuv, hat_guvv, hat_gvvv, hat_huu, hat_huv, hat_hvv, hat_huuu, hat_huuv, hat_huvv, hat_hvvv, _constant_f(hat_hu, hat_hv, Diff_2, approx_max_wave_number(fu, fv, fw, gu, gv, gw, hu, hv, hw, Diff_2, length, True)[0]**2 * np.pi**2 / length**2),
#                                 _constant_f_bar(hat_gv, hat_hv, Diff_2, approx_max_wave_number(fu, fv, fw, gu, gv, gw, hu, hv, hw, Diff_2, length, True)[0]**2 * np.pi**2 / length**2))[1] for length in domain_size_list], 'r*', label='A_2 at max bif point approx. with rdd wave number')
#     plt.plot(domain_size_list, [_calc_all_A_i(hat_guu, hat_guv, hat_gvv, hat_guuu, hat_guuv, hat_guvv, hat_gvvv, hat_huu, hat_huv, hat_hvv, hat_huuu, hat_huuv, hat_huvv, hat_hvvv, _constant_f(hat_hu, hat_hv, Diff_2, approx_max_wave_number(fu, fv, fw, gu, gv, gw, hu, hv, hw, Diff_2, length, True)[0]**2 * np.pi**2 / length**2),
#                                 _constant_f_bar(hat_gv, hat_hv, Diff_2, approx_max_wave_number(fu, fv, fw, gu, gv, gw, hu, hv, hw, Diff_2, length, True)[0]**2 * np.pi**2 / length**2))[2] for length in domain_size_list], 'b*', label='A_3 at max bif point approx. with rdd wave number')
#     plt.legend()
#     plt.figure()

# def plot_domain_size_A_i_not_rdd_wave_number(domain_size_list):
#     plt.plot(domain_size_list, [_calc_all_A_i(hat_guu, hat_guv, hat_gvv, hat_guuu, hat_guuv, hat_guvv, hat_gvvv, hat_huu, hat_huv, hat_hvv, hat_huuu, hat_huuv, hat_huvv, hat_hvvv, _constant_f(hat_hu, hat_hv, Diff_2, approx_max_wave_number(fu, fv, fw, gu, gv, gw, hu, hv, hw, Diff_2, length, False)[0]**2 * np.pi**2 / length**2),
#                                   _constant_f_bar(hat_gv, hat_hv, Diff_2, approx_max_wave_number(fu, fv, fw, gu, gv, gw, hu, hv, hw, Diff_2, length, False)[0]**2 * np.pi**2 / length**2))[0] for length in domain_size_list], 'g*', label='A_1 at max bif point approx., not rdd wave number')
#     plt.plot(domain_size_list, [_calc_all_A_i(hat_guu, hat_guv, hat_gvv, hat_guuu, hat_guuv, hat_guvv, hat_gvvv, hat_huu, hat_huv, hat_hvv, hat_huuu, hat_huuv, hat_huvv, hat_hvvv, _constant_f(hat_hu, hat_hv, Diff_2, approx_max_wave_number(fu, fv, fw, gu, gv, gw, hu, hv, hw, Diff_2, length, False)[0]**2 * np.pi**2 / length**2),
#                                   _constant_f_bar(hat_gv, hat_hv, Diff_2, approx_max_wave_number(fu, fv, fw, gu, gv, gw, hu, hv, hw, Diff_2, length, False)[0]**2 * np.pi**2 / length**2))[1] for length in domain_size_list], 'r*', label='A_2 at max bif point approx., not rdd wave number')
#     plt.plot(domain_size_list, [_calc_all_A_i(hat_guu, hat_guv, hat_gvv, hat_guuu, hat_guuv, hat_guvv, hat_gvvv, hat_huu, hat_huv, hat_hvv, hat_huuu, hat_huuv, hat_huvv, hat_hvvv, _constant_f(hat_hu, hat_hv, Diff_2, approx_max_wave_number(fu, fv, fw, gu, gv, gw, hu, hv, hw, Diff_2, length, False)[0]**2 * np.pi**2 / length**2),
#                                   _constant_f_bar(hat_gv, hat_hv, Diff_2, approx_max_wave_number(fu, fv, fw, gu, gv, gw, hu, hv, hw, Diff_2, length, False)[0]**2 * np.pi**2 / length**2))[2] for length in domain_size_list], 'b*', label='A_3 at max bif point approx., not rdd wave number')
#     plt.legend()
#     plt.figure()

# def plot_domain_size_integral_part_d_2(domain_size_list):
#     plt.plot(domain_size_list, [1/length*integrate.quad(
#         lambda x: _P_2(
#             _calc_eta_xi(hat_gu, hat_gv, hat_guu, hat_guv, hat_gvv, hat_hu, hat_hv, hat_huu, hat_huv, hat_hvv, _constant_f(hat_hu, hat_hv, Diff_2, approx_max_wave_number(fu, fv, fw, gu, gv, gw, hu, hv, hw, Diff_2, length, False)[0]**2 * np.pi**2 / length**2), approx_max_wave_number(fu, fv, fw, gu, gv, gw, hu, hv, hw, Diff_2, length, False)[1], Diff_2, approx_max_wave_number(fu, fv, fw, gu, gv, gw, hu, hv, hw, Diff_2, length, False)[0]**2 * np.pi**2 / length**2)[0],
#             _calc_eta_xi(hat_gu, hat_gv, hat_guu, hat_guv, hat_gvv, hat_hu, hat_hv, hat_huu, hat_huv, hat_hvv, _constant_f(hat_hu, hat_hv, Diff_2, approx_max_wave_number(fu, fv, fw, gu, gv, gw, hu, hv, hw, Diff_2, length, False)[0]**2 * np.pi**2 / length**2), approx_max_wave_number(fu, fv, fw, gu, gv, gw, hu, hv, hw, Diff_2, length, False)[1], Diff_2, approx_max_wave_number(fu, fv, fw, gu, gv, gw, hu, hv, hw, Diff_2, length, False)[0]**2 * np.pi**2 / length**2)[1],
#             approx_max_wave_number(fu, fv, fw, gu, gv, gw, hu, hv, hw, Diff_2, length, False)[0]**2 * np.pi**2 / length**2, 
#             x)[0] * np.cos(np.sqrt(approx_max_wave_number(fu, fv, fw, gu, gv, gw, hu, hv, hw, Diff_2, length, False)[0]**2 * np.pi**2 / length**2) * x) ** 2, 0, length)[0] for length in domain_size_list], 'g*', label='Integral_1 approx with NOT rdd wave number')
#     plt.plot(domain_size_list, [1/length*integrate.quad(
#         lambda x: _P_2(
#             _calc_eta_xi(hat_gu, hat_gv, hat_guu, hat_guv, hat_gvv, hat_hu, hat_hv, hat_huu, hat_huv, hat_hvv, _constant_f(hat_hu, hat_hv, Diff_2, approx_max_wave_number(fu, fv, fw, gu, gv, gw, hu, hv, hw, Diff_2, length, False)[0]**2 * np.pi**2 / length**2), approx_max_wave_number(fu, fv, fw, gu, gv, gw, hu, hv, hw, Diff_2, length, False)[1], Diff_2, approx_max_wave_number(fu, fv, fw, gu, gv, gw, hu, hv, hw, Diff_2, length, False)[0]**2 * np.pi**2 / length**2)[0],
#             _calc_eta_xi(hat_gu, hat_gv, hat_guu, hat_guv, hat_gvv, hat_hu, hat_hv, hat_huu, hat_huv, hat_hvv, _constant_f(hat_hu, hat_hv, Diff_2, approx_max_wave_number(fu, fv, fw, gu, gv, gw, hu, hv, hw, Diff_2, length, False)[0]**2 * np.pi**2 / length**2), approx_max_wave_number(fu, fv, fw, gu, gv, gw, hu, hv, hw, Diff_2, length, False)[1], Diff_2, approx_max_wave_number(fu, fv, fw, gu, gv, gw, hu, hv, hw, Diff_2, length, False)[0]**2 * np.pi**2 / length**2)[1],
#             approx_max_wave_number(fu, fv, fw, gu, gv, gw, hu, hv, hw, Diff_2, length, False)[0]**2 * np.pi**2 / length**2, 
#             x)[1] * np.cos(np.sqrt(approx_max_wave_number(fu, fv, fw, gu, gv, gw, hu, hv, hw, Diff_2, length, False)[0]**2 * np.pi**2 / length**2) * x) ** 2, 0, length)[0] for length in domain_size_list], 'r*', label='Integral_2 approx with NOT rdd wave number')
#     plt.legend()
#     # Without rounding j_max = sqrt(z_bar) to an integer, l_{j_max} and d_max are constant w.r.t. domain size.
#     # --> c_f's and thus A_i's are constant when j_max not rounded
#     # --> But d_2 not constant even when rounding since integrals 1/L*I_1, 1/L*I_2 not constant

# def plot_domain_size_detL2_max_bif(domain_size_list):
#     detL2_at_max_bif_list = [calc_det_L_2_at_max_bif(func_f, func_g, func_h, Diff_2, length, initial_guess_steady_state) for length in domain_size_list]
#     detL2_at_max_bif_approx = [det_L_2_direct_from_Jacobian(fu, fv, fw, gu, gv, gw, hu, hv, hw, Diff_2, length) for length in domain_size_list]
#     detL2_at_max_bif_limit = [get_det_L_2_limit(fu, fv, fw, gu, gv, gw, hu, hv, hw, Diff_2) for _ in domain_size_list]
#     plt.plot(domain_size_list, detL2_at_max_bif_list, 'r*')
#     plt.plot(domain_size_list, detL2_at_max_bif_approx, 'g*')
#     plt.plot(domain_size_list, detL2_at_max_bif_limit, 'b*')
#     plt.figure()

# def plot_fixed_diff_2(diff_2_list):
#     plt.plot(diff_2_list, [calc_det_L_2_at_max_bif(func_f, func_g, func_h, D_2, length_domain, initial_guess_steady_state) for D_2 in diff_2_list], 'r*')
#     plt.plot(diff_2_list, [det_L_2_direct_from_Jacobian(fu, fv, fw, gu, gv, gw, hu, hv, hw, D_2, length_domain) for D_2 in diff_2_list], 'g*')
#     plt.plot(diff_2_list, [get_det_L_2_limit(fu, fv, fw, gu, gv, gw, hu, hv, hw, D_2) for D_2 in diff_2_list], 'b*')
#     plt.figure()
#     plt.plot(diff_2_list, [calc_d_2_at_max_bif(func_f, func_g, func_h, D_2, length_domain, initial_guess_steady_state=initial_guess_steady_state) for D_2 in diff_2_list], 'g+')
#     plt.plot(diff_2_list, [d_2_direct_from_nonlinearities_approx(func_f, func_g, func_h, D_2, length_domain, initial_guess_steady_state=initial_guess_steady_state) for D_2 in diff_2_list], 'r+')
#     plt.plot(diff_2_list, [d_2_direct_from_nonlin_without_rdd(func_f, func_g, func_h, D_2, domain_length=100, initial_guess_steady_state=initial_guess_steady_state) for D_2 in diff_2_list], 'b+')
#     plt.plot(diff_2_list, [d_2_limit(func_f, func_g, func_h, D_2, initial_guess_steady_state=initial_guess_steady_state) for D_2 in diff_2_list], 'b+')
#     plt.figure()
