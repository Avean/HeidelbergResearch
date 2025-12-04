import matplotlib.pyplot as plt
import numpy as np
# import scipy.integrate as integrate
import sympy as sp

# from scipy.optimize import fsolve
from numpy import linalg as LA

from surpress_print import suppress_print
from derivatives import get_constant_steady_state, get_first_derivs_eval, get_all_hat_derivatives
from submatrices import check_for_unstable_submatrix
from formula_d_2 import calculate_d_2, calculate_d_2_faster, calculate_limit_d_2

class TuringBifurcationThreeEqu:
    def __init__(self, nonlinear_f, nonlinear_g, nonlinear_h, fixed_diffusion_coeff, domain_size, initial_guess):
        self.nonlinear_f = nonlinear_f
        self.nonlinear_g = nonlinear_g
        self.nonlinear_h = nonlinear_h
        self.D_2 = fixed_diffusion_coeff
        self.domain_size = domain_size
        # self.l_j = wavenumber ** 2 * np.pi ** 2 / L ** 2
        # self.wavenumber = wavenumber

        print(f'Domain size: {domain_size}, Diffusion D_2: {fixed_diffusion_coeff}')

        self.u, self.v, self.w = get_constant_steady_state(nonlinear_f, nonlinear_g, nonlinear_h, initial_guess)
        print(f'Constant steady state: ({self.u}, {self.v}, {self.w}).')
        self.fu, self.fv, self.fw, self.gu, self.gv, self.gw, self.hu, self.hv, self.hw = get_first_derivs_eval(nonlinear_f, nonlinear_g, nonlinear_h, self.u, self.v, self.w)
        self.Jacobian = np.array([[self.fu, self.fv, self.fw], [self.gu, self.gv, self.gw],
                                  [self.hu, self.hv, self.hw]])
        print('Jacobian:')
        print(self.Jacobian)

    def turing_bifurcation(self, wavenumber):
        l_j = wavenumber**2 * np.pi**2 / self.domain_size**2
        # At bifurcation point d_hat, we have a branch of non-constant steady states for diffusion coefficient
        # d(s) = d_hat + d_2*s^2 + d_3*s^3 + ...
        # Save d_hat and d_2 as attributes of class object to avoid unnecessary calculation steps.
        d_hat = self.get_bifurcation_point(l_j)
        d_2 = self.calculate_d_2(d_hat, l_j)
        print(f'Bifurcation point {d_hat} with branch d(s) = {d_hat} + {d_2}*s^2 + ... '
              f'for wavenumber {wavenumber}')
        print('')

        self.check_Turing_patterns_possible(d_hat, d_2, l_j, wavenumber)
        print(f'It should hold: sign derivative mu_bar(d_hat) = direction DDI: '
              f'{np.sign(self._calc_derv_mu_j(l_j))} = {self._get_direction_of_DDI(d_hat, l_j)}: '
              f'{np.sign(self._calc_derv_mu_j(l_j)) == self._get_direction_of_DDI(d_hat, l_j)}')
        print(f'If DDI, it holds sign critical eigv = -1 * sign d_2 * sign derivative mu_bar(d_hat) = -1 * sign d_2 * direction DDI. '
              f'Hence, if direction of branch in direction of DDI, we automatically have negative critical eigenvalue.')

        return d_hat, d_2

    def check_stability_of_Jacobian(self):
        det_J12 = self.fu * self.gv - self.fv * self.gu
        det_J13 = self.fu * self.hw - self.fw * self.hu
        det_J23 = self.gv * self.hw - self.gw * self.hv
        det_J = LA.det(self.Jacobian)
        p1 = - (self.fu + self.gv + self.hw)
        p2 = det_J12 + det_J23 + det_J13
        p3 = - det_J

        if p1 > 0 and p3 > 0 and p1 * p2 - p3 > 0:
            print(f'Jacobian at steady state ({self.u}, {self.v}, {self.w}) stable')
            return True
        else:
            print(f'Jacobian at steady state ({self.u}, {self.v}, {self.w}) unstable')
            return False

    def check_DDI_possible(self):
        Jac_stable = self.check_stability_of_Jacobian()
        if Jac_stable and check_for_unstable_submatrix(self.Jacobian):
            print('Jacobian stable and at least one submatrix unstable. Hence, DDI possible.')
            return True
        elif not Jac_stable:
            print('Jacobian not stable, thus, no DDI possible.')
            return False
        else:
            print('Jacobian stable but also all submatrices, thus, no DDI since system of three equations.')
            return False

    def get_bifurcation_point(self, l_j, info=False):
        # Three possible bifurcation points of two different types.
        # According to Ruth-Hurwitz, two conditions can break under assumption that Jacobian J is stable.
        d_j = self._get_bifurcation_p3(l_j)
        print(f'd_j = {d_j}')
        d_j_plus, d_j_minus = self._get_bifurcation_p1p2_p3(l_j)
        print(f'd_j_minus = {d_j_minus} and d_j_plus = {d_j_plus}')
        if np.isreal(d_j_plus) and max(d_j_plus, d_j_minus)> d_j:
            if info:
                return max(d_j_plus, d_j_minus), 'complex eigvs'
            else:
                return max(d_j_plus, d_j_minus)
        else:
            if info:
                return d_j, 'simple eigv'
            else:
                return d_j

    def calculate_d_2(self, d_hat, l_j, faster=True):
        derivs_hat_g, derivs_hat_h = get_all_hat_derivatives(self.nonlinear_f, self.nonlinear_g, self.nonlinear_h, self.u, self.v, self.w)

        # print('Jacobian of g(p(v,w),v,w) and h(p(v,w),v,w) at steady state, where f(p(v,w),v,w)=0 around steady state:')
        # print(np.array([derivs_hat_g[0:2], derivs_hat_h[0:2]]))

        if faster:
            d_2 = calculate_d_2_faster(*derivs_hat_g, *derivs_hat_h, d_hat, self.D_2, self.domain_size, l_j)
        else:
            d_2 = calculate_d_2(*derivs_hat_g, *derivs_hat_h, d_hat, self.D_2, self.domain_size, l_j)
        return d_2
    
    def d_2_limit_wrt_domain_size(self):
        # We know that l_max and d_max are constant w.r.t. domain length if we do not round the approximated j_max to an integer. (j_max is of course not constant but O(L^2))
        j_max, d_max = self.approx_max_wave_number_p3(rounded=False, both=False)
        l_max = j_max**2 * np.pi ** 2 / self.domain_size ** 2
        derivs_hat_g, derivs_hat_h = get_all_hat_derivatives(self.nonlinear_f, self.nonlinear_g, self.nonlinear_h, self.u, self.v, self.w)
        d_2 = calculate_limit_d_2(*derivs_hat_g, *derivs_hat_h, d_max=d_max, D_2=self.D_2, l_max=l_max)
        return d_2
    
    def B_k_limit_wrt_domain_size(self):
        # We know that l_max and d_max are constant w.r.t. domain length if we do not round the approximated j_max to an integer. (j_max is of course not constant but O(L^2))
        j_max, d_max = self.approx_max_wave_number_p3(rounded=False, both=False)
        l_max = j_max**2 * np.pi ** 2 / self.domain_size ** 2
        derivs_hat_g, derivs_hat_h = get_all_hat_derivatives(self.nonlinear_f, self.nonlinear_g, self.nonlinear_h, self.u, self.v, self.w)
        d_2 = calculate_limit_d_2(*derivs_hat_g, *derivs_hat_h, d_max=d_max, D_2=self.D_2, l_max=l_max)
        return 2*l_max*d_2

    def calc_d_2_at_max_bif(self):
        with suppress_print():
            max_wave_number, max_bif_point = self.search_max_bifurcation_point()
            # d_2_at_max_bif = turing_bif.turing_bifurcation(max_wave_number)
            d_2_at_max_bif = self.calculate_d_2(max_bif_point, max_wave_number**2*np.pi**2 / self.domain_size**2)
        return d_2_at_max_bif

    def check_Turing_patterns_possible(self, d_hat, d_2, l_j, wavenumber):
        correct_direction = self._direction_of_bifurcation(d_hat, d_2, l_j)
        stability_crit_eigv = self._calculate_sign_critical_eigv(d_hat, d_2, l_j)
        if self.fu < 0 and correct_direction and stability_crit_eigv == -1.0:
            print('Turing patterns may emerge (but not for sure) under assumption of DDI!')
        elif self.fu >= 0:
            print('Autocatalysis condition, thus all regular patterns crossing constant steady state unstable.')
        elif not correct_direction or stability_crit_eigv == 1.0:
            print(f'No Turing patterns of wavenumber {wavenumber} possible! '
                  f'Either bifurcation in wrong direction or pattern unstable.')
        else:
            print(f'Direction of DDI and branch coincide.'
                  f'However we cannot calculate sign of critical eigenvalue as we are in the case d_hat = d_j_plus')

    def check_Turing_patterns_sure(self, max_range=801):
        with suppress_print():
            ind_max, d_max = self.search_max_bifurcation_point(max_range)
        l_max = np.pi**2 * ind_max**2 / self.domain_size**2
        d_2 = self.calculate_d_2(d_max, l_max)
        print(f'Maximal bifurcation point is {d_max} with corresponding wavenumber {ind_max}.')
        correct_direction = self._direction_of_bifurcation(d_max, d_2, l_max)
        stability_crit_eigv = self._calculate_sign_critical_eigv(d_max, d_2, l_max)
        if self.fu < 0 and correct_direction and stability_crit_eigv == -1.0:
            print('Turing patterns emerge for max. wave number!')
        elif self.fu >= 0:
            print('Autocatalysis condition, thus all regular patterns crossing constant steady state unstable.')
        elif not correct_direction or stability_crit_eigv == 1.0:
            print(f'No Turing patterns of max. wavenumber {ind_max} possible! '
                  f'Either bifurcation in wrong direction or pattern unstable.')
        else:
            print(f'Direction of DDI and branch coincide.'
                  f'However we cannot calculate sign of critical eigenvalue as we are in the case d_hat = d_j_plus')


    def calculate_B_j_limit_wrt_j(self):
        derivs_hat_g, derivs_hat_h = get_all_hat_derivatives(self.nonlinear_f, self.nonlinear_g, self.nonlinear_h, self.u, self.v, self.w)
        hat_gv, hat_gw, hat_gvv, hat_gvw, hat_gww, hat_gvvv, hat_gvvw, hat_gvww, hat_gwww = derivs_hat_g
        hat_hv, hat_hw, hat_hvv, hat_hvw, hat_hww, hat_hvvv, hat_hvvw, hat_hvww, hat_hwww = derivs_hat_h

        det = hat_gv * hat_hw - hat_gw * hat_hv
        limit_B_j = - 1 / (2 * det) * (hat_hw * hat_gvv**2 - hat_gw * hat_hvv * hat_gvv - hat_hv * hat_gvv * hat_gvw + hat_gv * hat_hvv * hat_gvw) + 1/4 * hat_gvvv + hat_gvv**2 / (12*hat_gv)
        return limit_B_j

    def search_max_bifurcation_point(self, max_range=801):
        [d_p3_ind_max, d_p3_max], [d_p1p2_p3_ind_max, d_p1p2_p3_max] = self._search_max_bifurcation_point_both(max_range)
        if d_p3_max >= d_p1p2_p3_max:
            return d_p3_ind_max, d_p3_max
        else:
            return d_p1p2_p3_ind_max, d_p1p2_p3_max

    def approx_max_wave_number_p3(self, rounded: bool, both: bool=False):
        detJ12 = self.fu*self.gv - self.gu*self.fv
        detJ13 = self.fu*self.hw - self.hu*self.fw
        detJ = LA.det(self.Jacobian)

        c1 = detJ12 * self.D_2 * np.pi**2 * self.domain_size**2
        c2 = detJ * self.domain_size**4
        c3 = self.fu * self.D_2 * np.pi**4
        c4 = detJ13 * self.domain_size**2 * np.pi**2

        z_bar = c2/c1 + np.emath.sqrt((c2/c1)**2 -(c2*c4)/(c1*c3))
        if z_bar <= 0:
            print(f'z_bar is not positive: {z_bar} for domain size: {self.domain_size}')
            return np.nan, np.nan
        if rounded:
            j_max_rdd = int(np.round(np.emath.sqrt(z_bar)))
            if j_max_rdd <= 0:
                j_max_rdd = 1
            if both:
                if j_max_rdd**2 - z_bar < 0:
                    j_max = [j_max_rdd, j_max_rdd+1]
                else:
                    j_max = [max(j_max_rdd-1,0), j_max_rdd]
            else:
                poss_j_max = [max(j_max_rdd - 1,1), j_max_rdd, j_max_rdd + 1]
                poss_wave_numbs = [np.pi**2 * j**2 / self.domain_size**2 for j in poss_j_max]
                poss_bif_points = [self._get_bifurcation_p3(wave_numb) for wave_numb in poss_wave_numbs]
                index_max_wave = np.argmax(poss_bif_points)
                j_max = poss_j_max[index_max_wave]
                j_max_rdd = poss_j_max[index_max_wave]
        else:
            j_max = np.emath.sqrt(z_bar)
        # sqrt_for_dmax = c1 * np.emath.sqrt((c2/c1)**2 - (c2*c4)/(c1*c3))
        # d_max_2 = (c1**2 * sqrt_for_dmax) / (2*c2*(c2*c3-c1*c4) + (2*c2*c3-c1*c4)*sqrt_for_dmax)
        if rounded:
            d_max = (c1*j_max_rdd**2 - c2) / (c3*j_max_rdd**4 - c4*j_max_rdd**2)
        else:
            d_max = (c1*j_max**2 - c2) / (c3*j_max**4 - c4*j_max**2)
        # l_k = j_max**2 * np.pi ** 2 / domain_length ** 2
        return j_max, d_max

    def bifurcation_at_max_wave_num_p3(self):
        wave_j = self.approx_max_wave_number_p3(rounded=True, both=True)[0]
        print("Investigation of branch for wavenumber j=" + str(wave_j[0]))
        self.turing_bifurcation(max(wave_j[0],1))
        print('\n')
        print("Investigation of branch for wavenumber j=" + str(wave_j[1]))
        self.turing_bifurcation(max(wave_j[1],1))
    
    def plot_B_k_for_multiple_domain_sizes(self, domain_sizes=np.linspace(0.01, 15, 100)):
        saved_domain_size = self.domain_size
        B_k_list_pos = []
        B_k_list_neg = []
        with suppress_print():
            for domain_size in domain_sizes:
                self.domain_size = domain_size
                k_max, d_max = self.approx_max_wave_number_p3(rounded=True, both=False)
                l_k = k_max**2 * np.pi**2 / self.domain_size**2
                ind_max_searched, d_max_searched = self.search_max_bifurcation_point(max_range=10+2*k_max)
                if d_max > 0 and d_max >= d_max_searched - 1e-3:
                    d_2 = self.calculate_d_2(d_max, l_k)
                    B_k = 2 * l_k * d_2
                    if B_k >= 0:
                        B_k_list_pos.append(B_k)
                        B_k_list_neg.append(np.nan)
                    else:
                        B_k_list_pos.append(np.nan)
                        B_k_list_neg.append(B_k)
                else:
                    B_k_list_neg.append(np.nan)
                    B_k_list_pos.append(np.nan)
        self.domain_size = saved_domain_size

        limit_B_k = self.B_k_limit_wrt_domain_size()

        plt.figure(figsize=(6, 4))
        colors = ['#1b9e77', '#d95f02', '#7570b3']
        markers = ['o', 's', '^', 'D']
        plt.plot(domain_sizes, B_k_list_neg, color=colors[0], linestyle='-', linewidth=5, markeredgecolor='black', label='$\hat{d} = d_{k(L)}^0 > 0$ and $\mathcal{B}(k(L))<0$')
        plt.plot(domain_sizes, B_k_list_pos, color=colors[1], linestyle='-', linewidth=5, markeredgecolor='black', label='$\hat{d} = d_{k(L)}^0 > 0$ and $\mathcal{B}(k(L))\geq0$')
        plt.plot(domain_sizes, [limit_B_k for _ in domain_sizes], color=colors[2], linestyle='--', label='Limit $\mathcal{B}(k(L))$ w.r.t. $L$')
        plt.xlabel('Domain size $L$', fontsize=16)
        plt.ylabel('$\mathcal{B}(k(L))$', fontsize=16)
        # ax = plt.gca()
        # plt.xticks(domain_sizes, fontsize=10)
        # ticks = ax.get_xticks()
        # labels = [str(int(tick)) if i % 4 == 0 else '' for i, tick in enumerate(ticks)]
        # ax.set_xticklabels(labels)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.ylim(-4, 0)
        plt.legend(fontsize=14)
        plt.tight_layout()
        # plt.savefig("plot_B_k_multiple_domain_sizes.pdf", dpi=300)
        plt.show()

    def plot_B_j_for_multiple_wave_numbers(self, j_int = np.linspace(1, 21, 21)):
        B_j_neg = []
        B_j_pos = []
        B_j_neg_noDDI = []
        B_j_pos_noDDI = []
        with suppress_print():
            for j in j_int:
                lam_j = j**2*np.pi**2 / self.domain_size**2
                d_hat_j = self._get_bifurcation_p3(lam_j)
                d_2_j = self.calculate_d_2(d_hat_j, lam_j, faster=True)
                B_j = 2 * lam_j * d_2_j
                if d_hat_j > 0:
                    if d_2_j < 0:
                        B_j_neg.append(B_j)
                        B_j_pos.append(np.nan)
                        B_j_neg_noDDI.append(np.nan)
                        B_j_pos_noDDI.append(np.nan)
                    else:
                        B_j_neg.append(np.nan)
                        B_j_pos.append(B_j)
                        B_j_neg_noDDI.append(np.nan)
                        B_j_pos_noDDI.append(np.nan)
                else:
                    if d_2_j < 0:
                        B_j_neg.append(np.nan)
                        B_j_pos.append(np.nan)
                        B_j_neg_noDDI.append(B_j)
                        B_j_pos_noDDI.append(np.nan)
                    else:
                        B_j_neg.append(np.nan)
                        B_j_pos.append(np.nan)
                        B_j_neg_noDDI.append(np.nan)
                        B_j_pos_noDDI.append(B_j)

        max_bif_point_j, bif_point_max = self.search_max_bifurcation_point(max_range=301)
        max_bif_l_j = max_bif_point_j**2 * np.pi**2 / self.domain_size**2
        
        B_j_limit = self.calculate_B_j_limit_wrt_j()

        plt.figure(figsize=(6, 4))
        colors = ['#1b9e77', '#d95f02', '#7570b3']
        # colors2 = ['#66c2a5', '#fc8d62', '#8da0cb']
        # colors3 = ['#009E73', '#D55E00', '#0072B2']
        markers = ['o', 's', '^', 'D']
        plt.plot(j_int, B_j_neg, marker=markers[0], color=colors[0], linestyle='None', ms=8, markeredgecolor='black', label='$d_j^0>0$ and $\mathcal{B}(j)<0$')
        plt.plot(j_int, B_j_pos, marker=markers[1], color=colors[1], linestyle='None', ms=8, markeredgecolor='black', label='$d_j^0>0$ and $\mathcal{B}(j) \geq 0$')
        plt.plot(j_int, [B_j_limit for _ in j_int], color=colors[2], linestyle='--', label='Limit $\mathcal{B}(j)$ w.r.t. j')
        # plt.plot(j_int, d_2_neg_noDDI, marker=markers[2], color=colors[0], linestyle='None', ms=8, markeredgecolor='black', label='$d_j^0 \leq 0$ and $\mathcal{B}(j)<0$')
        # plt.plot(j_int, d_2_pos_noDDI, marker=markers[3], color=colors[1], linestyle='None', ms=8, markeredgecolor='black', label='$d_j^0 \leq 0$ and $\mathcal{B}(j) \geq 0$')
        plt.plot(max_bif_point_j, 2*max_bif_l_j*self.calculate_d_2(bif_point_max, max_bif_l_j, faster=True), marker=markers[0], color=colors[2], linestyle='None', ms=14, fillstyle='none', label='Maximal bifurcation point')
        plt.xlabel('wave number j', fontsize=12)
        plt.ylabel('$\mathcal{B}(j)$', fontsize=12)
        ax = plt.gca()
        plt.xticks(j_int, fontsize=10)
        ticks = ax.get_xticks()
        labels = [str(int(tick)) if i % 4 == 0 else '' for i, tick in enumerate(ticks)]
        ax.set_xticklabels(labels)
        plt.yticks(fontsize=10)
        plt.legend(fontsize=10)
        plt.tight_layout()
        # plt.savefig("plot_B_j_multiple_wavenumbers.pdf", dpi=300)
        plt.show()

    def plot_d_2_for_multiple_wave_numbers(self, j_int = np.linspace(1, 21, 21)):  # np.linspace(350, 450, 101)
        # d_2_vals = [self.calculate_d_2_only_depend_l_j(j**2*np.pi**2 /
        # self.domain_size**2) for j in j_int]
        d_2_neg = []
        d_2_pos = []
        d_2_neg_noDDI = []
        d_2_pos_noDDI = []
        with suppress_print():
            for j in j_int:
                lam_j = j**2*np.pi**2 / self.domain_size**2
                d_hat_j = self._get_bifurcation_p3(lam_j)
                d_2_j = self.calculate_d_2(d_hat_j, lam_j, faster=True)
                if d_hat_j > 0:
                    if d_2_j < 0:
                        d_2_neg.append(d_2_j)
                        d_2_pos.append(np.nan)
                        d_2_neg_noDDI.append(np.nan)
                        d_2_pos_noDDI.append(np.nan)
                    else:
                        d_2_neg.append(np.nan)
                        d_2_pos.append(d_2_j)
                        d_2_neg_noDDI.append(np.nan)
                        d_2_pos_noDDI.append(np.nan)
                else:
                    if d_2_j < 0:
                        d_2_neg.append(np.nan)
                        d_2_pos.append(np.nan)
                        d_2_neg_noDDI.append(d_2_j)
                        d_2_pos_noDDI.append(np.nan)
                    else:
                        d_2_neg.append(np.nan)
                        d_2_pos.append(np.nan)
                        d_2_neg_noDDI.append(np.nan)
                        d_2_pos_noDDI.append(d_2_j)

        max_bif_point_j, bif_point_max = self.search_max_bifurcation_point(max_range=301)

        limit_B_j = self.calculate_B_j_limit_wrt_j()

        plt.figure(figsize=(6, 4))
        colors = ['#1b9e77', '#d95f02', '#7570b3']
        # colors2 = ['#66c2a5', '#fc8d62', '#8da0cb']
        # colors3 = ['#009E73', '#D55E00', '#0072B2']
        markers = ['o', 's', '^', 'D']
        plt.plot(j_int, d_2_neg, marker=markers[0], color=colors[0], linestyle='None', ms=8, markeredgecolor='black', label='$d_j^0>0$ and $d_2<0$')
        plt.plot(j_int, d_2_pos, marker=markers[1], color=colors[1], linestyle='None', ms=8, markeredgecolor='black', label='$d_j^0>0$ and $d_2 \geq 0$')
        plt.plot(j_int, [limit_B_j/(2*j**2*np.pi**2 / self.domain_size**2) for j in j_int], color=colors[2], linestyle='--', label='Limit $\mathcal{B}(j)$ w.r.t. j / (2l_j)')
        # plt.plot(j_int, d_2_neg_noDDI, marker=markers[2], color=colors[0], linestyle='None', ms=8, markeredgecolor='black', label='$d_j^0 \leq 0$ and $d_2<0$')
        # plt.plot(j_int, d_2_pos_noDDI, marker=markers[3], color=colors[1], linestyle='None', ms=8, markeredgecolor='black', label='$d_j^0 \leq 0$ and $d_2 \geq 0$')
        plt.plot(max_bif_point_j, self.calculate_d_2(bif_point_max, max_bif_point_j**2*np.pi**2/self.domain_size**2, faster=True), marker=markers[0], color=colors[2], linestyle='None', ms=14, fillstyle='none', label='Maximal bifurcation point')
        plt.xlabel('wave number j', fontsize=12)
        plt.ylabel('$d_2$ in $d(s) = d_j^0 + d_2 s^2 + \dots$', fontsize=12)
        ax = plt.gca()
        plt.xticks(j_int, fontsize=10)
        ticks = ax.get_xticks()
        labels = [str(int(tick)) if i % 2 == 0 else '' for i, tick in enumerate(ticks)]
        ax.set_xticklabels(labels)
        plt.yticks(fontsize=10)
        plt.legend(fontsize=10)
        plt.tight_layout()
        # plt.savefig("plot_d_2_multiple_wavenumbers.pdf", dpi=300)
        plt.show()

    def plot_all_possible_bifurcation_points(self, range_wavenumber=range(1, 21)):
        bifurcation_points_p3 = []
        bifurcation_points_hopf = []
        for j in range_wavenumber:
            lam_j = j**2 * np.pi**2 / self.domain_size**2
            d_j_p3 = self._get_bifurcation_p3(lam_j)
            d_j_Hopf = self._get_bifurcation_p1p2_p3(lam_j)[0]
            if np.isreal(d_j_p3) and d_j_p3 > 0:
                bifurcation_points_p3.append(d_j_p3)
            else:
                bifurcation_points_p3.append(np.nan)
            if np.isreal(d_j_Hopf) and d_j_Hopf > 0:
                bifurcation_points_hopf.append(d_j_Hopf)
            else:
                bifurcation_points_hopf.append(np.nan)
        
        max_bif_point_j, bif_point_max = self.search_max_bifurcation_point(max_range=301)
        
        plt.figure(figsize=(6, 4))
        colors = ['#1b9e77', '#d95f02', '#7570b3']
        # colors2 = ['#66c2a5', '#fc8d62', '#8da0cb']
        # colors3 = ['#009E73', '#D55E00', '#0072B2']
        markers = ['o', 's', '^', 'D']
        plt.plot(range_wavenumber, bifurcation_points_p3, marker=markers[0], color=colors[0], linestyle='None', ms=6, markeredgecolor='black', label='Bifurcation candidate $d_j^0$')
        plt.plot(range_wavenumber, bifurcation_points_hopf, marker=markers[3], color=colors[1], linestyle='None', ms=6, markeredgecolor='black', label='Bifurcation candidate $d_j^+$')
        plt.plot(max_bif_point_j, bif_point_max, marker=markers[0], color=colors[2], linestyle='None', ms=11, fillstyle='none', label='Maximal bifurcation point')
        plt.xlabel('wavenumber j', fontsize=12)
        plt.ylabel('Bifurcation candidates', fontsize=12)
        # plt.xticks([min(range_wavenumber), *range(min(range_wavenumber)+1, max(range_wavenumber)+1, 2)], fontsize=10)
        ax = plt.gca()
        plt.xticks(range_wavenumber, fontsize=10)
        ticks = ax.get_xticks()
        labels = [str(int(tick)) if i % 2 == 0 else '' for i, tick in enumerate(ticks)]
        ax.set_xticklabels(labels)
        plt.yticks(fontsize=10)
        plt.legend(fontsize=10)
        plt.tight_layout()
        # plt.savefig("plot_all_bifurcation_candidates.pdf", dpi=300)
        plt.show()

    def plot_bif_point_d_2_for_mult_wave_numbers(self, range_wavenumber=range(1, 21)):
        with suppress_print():
            max_bif_point_j = self.search_max_bifurcation_point(max_range=301)[0]
            bifurcation_points_simple = []
            bifurcation_points_complex = []
            direc_branch_d_2_pos = []
            direc_branch_d_2_neg = []
            for j in range_wavenumber:
                dhat, type_point = self.get_bifurcation_point(j**2 * np.pi**2 / self.domain_size**2, info=True)
                if np.isreal(dhat) and dhat > 0:
                    if type_point == 'simple eigv':
                        bifurcation_points_simple.append(dhat)
                        bifurcation_points_complex.append(np.nan)
                    else:
                        bifurcation_points_simple.append(np.nan)
                        bifurcation_points_complex.append(dhat)
                    d_2 = self.calculate_d_2(dhat, j**2 * np.pi**2 / self.domain_size**2)
                    if d_2 >= 0:
                        direc_branch_d_2_pos.append(d_2)
                        direc_branch_d_2_neg.append(np.nan)
                    else:
                        direc_branch_d_2_pos.append(np.nan)
                        direc_branch_d_2_neg.append(d_2)
                else:
                    bifurcation_points_simple.append(np.nan)
                    bifurcation_points_complex.append(np.nan)
                    direc_branch_d_2_pos.append(np.nan)
                    direc_branch_d_2_neg.append(np.nan)
            plt.plot(range_wavenumber, bifurcation_points_simple, 'b*', label='Simple eigenvalue')
            plt.plot(range_wavenumber, bifurcation_points_complex, 'r*', label='Complex pair of eigenvalues')
            if not np.isnan(bifurcation_points_complex[max_bif_point_j-1]):
                bif_point_max = bifurcation_points_complex[max_bif_point_j-1]
            else:
                bif_point_max = bifurcation_points_simple[max_bif_point_j-1]
            plt.plot(max_bif_point_j, bif_point_max, 'go', ms=9, fillstyle='none', label='Maximal bifurcation point')
            plt.xlabel('wavenumber j')
            plt.xticks(list(range_wavenumber))
            plt.ylabel('Bifurcation candidate b_j')
            plt.legend()
            plt.show()

            d_2_limit = self.d_2_limit_wrt_domain_size()
            plt.plot(range_wavenumber, direc_branch_d_2_pos, 'ro', label='Positive values of d_2')
            plt.plot(range_wavenumber, direc_branch_d_2_neg, 'bo', label='Negative values of d_2')
            #plt.plot(range_wavenumber, direc_branch_d_2, 'bo', label='Value of d_2')
            if not np.isnan(direc_branch_d_2_pos[max_bif_point_j-1]):
                d_2_at_max_bif = direc_branch_d_2_pos[max_bif_point_j-1]
            else:
                d_2_at_max_bif = direc_branch_d_2_neg[max_bif_point_j-1]
            plt.plot(max_bif_point_j, d_2_at_max_bif, 'go', ms=11, fillstyle='none', label='d_2 at maximal bifurcation point')
            plt.plot(range_wavenumber, [d_2_limit for _ in range_wavenumber], 'k--', ms=1, label='Limit of d_2 at max. point w.r.t. domain size')
            plt.xlabel('wavenumber j')
            plt.xticks(list(range_wavenumber))
            plt.ylabel('d_2')
            plt.legend()
            plt.show()

    # def investigate_change_d_2_wrt_wavenumber(self, j_int = np.linspace(1, 50, 50)):
    #     derivs_hat_g, derivs_hat_h = get_all_hat_derivatives(self.nonlinear_f, self.nonlinear_g, self.nonlinear_h, self.u, self.v, self.w)

    #     fu, fv, fuu, fuv, fvv, fuuu, fuuv, fuvv, fvvv = derivs_hat_g
    #     gu, gv, guu, guv, gvv, guuu, guuv, guvv, gvvv = derivs_hat_h

    #     l_j_int = [j**2 * np.pi**2 / self.domain_size**2 for j in j_int]

    #     ## const_f_vals = [_constant_f(gu, gv, self.D_2, lam) for lam in l_j_int]
    #     ## const_f_bar_vals = [_constant_f_bar(fv, gv, self.D_2, lam) for lam in l_j_int]
    #     ##
    #     ## def all_A_i(lam):
    #     ##     const_f = _constant_f(gu, gv, self.D_2, lam)
    #     ##     const_f_bar = _constant_f_bar(fv, gv, self.D_2, lam)
    #     ##     A_1 = _calc_A_1(fuu, fuv, guu, guv, const_f, const_f_bar)
    #     ##     A_2 = _calc_A_2(fuv, fvv, guv ,gvv, const_f, const_f_bar)
    #     ##     A_3 = _calc_A_3(fuuu, fuuv, fuvv, fvvv, guuu, guuv, guvv, gvvv, const_f, const_f_bar)
    #     ##     return A_1, A_2, A_3
    #     ## A_1_vals = [all_A_i(lam)[0] for lam in l_j_int]
    #     ## A_2_vals = [all_A_i(lam)[1] for lam in l_j_int]
    #     ## A_3_vals = [all_A_i(lam)[2] for lam in l_j_int]

    #     def p_2(lam):
    #         with suppress_print():
    #             d_hat = self.get_bifurcation_point(lam)
    #         const_f = _constant_f(gu, gv, self.D_2, lam)
    #         eta, xi = _calc_eta_xi(fu, fv, fuu, fuv, fvv, gu, gv, guu, guv, gvv, const_f, d_hat, self.D_2, lam)
    #         u_2_eta = eta[0]
    #         v_2_eta = eta[1]
    #         u_2_xi = xi[0]
    #         v_2_xi = xi[1]
    #         return u_2_eta, u_2_xi, v_2_eta, v_2_xi

    #     eta_xi = [p_2(lam) for lam in l_j_int]
    #     u2_eta = [eta_xi[k][0] for k in range(0, len(eta_xi))]
    #     u2_xi = [eta_xi[k][1] for k in range(0, len(eta_xi))]
    #     v2_eta = [eta_xi[k][2] for k in range(0, len(eta_xi))]
    #     v2_xi = [eta_xi[k][3] for k in range(0, len(eta_xi))]
    #     p_2_vals = [u2_eta, u2_xi, v2_eta, v2_xi]

    #     def p_2_int(lam):
    #         with suppress_print():
    #             d_hat = self.get_bifurcation_point(lam)
    #         const_f = _constant_f(gu, gv, self.D_2, lam)
    #         eta, xi = _calc_eta_xi(fu, fv, fuu, fuv, fvv, gu, gv, guu, guv, gvv, const_f, d_hat, self.D_2, lam)
    #         integral_1 = integrate.quad(
    #             lambda x: _P_2(eta, xi, lam, x)[0] * np.cos(
    #                 np.sqrt(lam) * x) ** 2, 0, self.domain_size)[0]
    #         integral_2 = integrate.quad(
    #             lambda x: _P_2(eta, xi, lam, x)[1] * np.cos(
    #                 np.sqrt(lam) * x) ** 2, 0, self.domain_size)[0]
    #         return integral_1, integral_2

    #     integral_1_vals = [p_2_int(lam)[0] for lam in l_j_int]
    #     integral_2_vals = [p_2_int(lam)[1] for lam in l_j_int]

    #     ## plt.plot(j_int, const_f_vals, '.', label='Constant f')
    #     ## plt.plot(j_int, const_f_bar_vals, '*', label='Constant f_bar')
    #     ## plt.legend()
    #     ##
    #     ## plt.figure()
    #     ## plt.plot(j_int, A_1_vals, '.', label='A_1')
    #     ## plt.plot(j_int, A_2_vals, '+', label='A_2')
    #     ## plt.plot(j_int, A_3_vals, '*', label='A_3')
    #     ## plt.legend()

    #     plt.figure()
    #     plt.plot(j_int, integral_1_vals, '.', label='Integral 1')
    #     plt.plot(j_int, integral_2_vals, '+', label='Integral 2')
    #     plt.legend()

    #     plt.figure()
    #     labels = ['u_2 eta', 'u_2 xi', 'v_2 eta', 'v_2 xi']
    #     styles = ['.', 'o', '+', '*']
    #     for y_arr, label, style in zip(p_2_vals, labels, styles):
    #         plt.plot(j_int, y_arr, style, label=label)
    #     plt.legend()

    #     def L_2_h(lam):
    #         with suppress_print():
    #             d_hat = self.get_bifurcation_point(lam)
    #         const_f = _constant_f(gu, gv, self.D_2, lam)
    #         L_2 = np.array([[fu - 4 * d_hat * lam, fv], [gu, gv - 4 * self.D_2 * lam]])
    #         h = np.array([1 / 2 * (1 / 2 * fuu + fuv * const_f + 1 / 2 * fvv * const_f ** 2),
    #                       1 / 2 * (1 / 2 * guu + guv * const_f + 1 / 2 * gvv * const_f ** 2)])
    #         return L_2, h

    #     def det_L_2(j):
    #         l_j = j ** 2 * np.pi ** 2 / self.domain_size ** 2
    #         with suppress_print():
    #             d_j = self.get_bifurcation_point(l_j)
    #         return (fu - 4 * d_j * l_j) * (gv - 4 * self.D_2 * l_j) - fv * gu
    #     # print(fsolve(det_L_2, 36))

    #     # J = np.linspace(1, 50, 50)
    #     L_2_1, L_2_2, L_2_3, L_2_4, h_1, h_2, L_2_det, xi_1, xi_2 = [], [], [], [], [], [], [], [], []
    #     for l_j in l_j_int: # j in J:
    #         L_2, h = L_2_h(l_j)
    #         L_2_1.append(L_2[0][0])
    #         L_2_2.append(L_2[0][1])
    #         L_2_3.append(L_2[1][0])
    #         L_2_4.append(L_2[1][1])
    #         h_1.append(h[0])
    #         h_2.append(h[1])
    #         L_2_det.append(LA.det(L_2))
    #         xi = np.linalg.solve(L_2, -h)
    #         xi_1.append(xi[0])
    #         xi_2.append(xi[1])

    #     plt.figure()
    #     labels = ['L_2 det', 'xi_1', 'xi_2']  # 'L_2 1', 'L_2 2', 'L_2 3', 'L_2 4', 'h 1', 'h 2',
    #     L_vals = [L_2_det, xi_1, xi_2]  # L_2_1, L_2_2, L_2_3, L_2_4, h_1, h_2,
    #     for y_arr, label in zip(L_vals, labels):
    #         plt.plot(j_int, y_arr, label=label)
    #     plt.legend()

    #     plt.show()

    def _get_bifurcation_p1p2_p3(self, l_j):
        det_J12 = self.fu * self.gv - self.fv * self.gu
        tr_J12 = self.fu + self.gv
        det_J13 = self.fu * self.hw - self.fw * self.hu
        tr_J13 = self.fu + self.hw
        det_J23 = self.gv * self.hw - self.gw * self.hv
        tr_J = self.fu + self.gv + self.hw
        p1 = - (self.fu + self.gv + self.hw)
        p2 = det_J12 + det_J23 + det_J13
        p3 = - LA.det(self.Jacobian)
        a0 = p1*p2 - p3 + (det_J23 + det_J13 + tr_J * tr_J12) * l_j * self.D_2 - tr_J12 * l_j**2 * self.D_2**2
        a1 = (det_J12 + det_J23 + tr_J * tr_J13) * l_j - 2 * tr_J * self.D_2 * l_j**2 + self.D_2**2 * l_j**3
        a2 = self.D_2 * l_j**3 - tr_J13 * l_j**2
        d_j_plus = (-a1 + np.emath.sqrt(a1**2 - 4*a0*a2)) / (2 * a2)
        d_j_minus = (-a1 - np.emath.sqrt(a1**2 - 4*a0*a2)) / (2 * a2)
        if np.isreal(d_j_plus) and np.isreal(d_j_minus):
            return np.real(d_j_plus), np.real(d_j_minus)
        else:
            return d_j_plus, d_j_minus

    def _get_bifurcation_p3(self, l_j):
        det_J = LA.det(self.Jacobian)
        det_J12 = self.fu * self.gv - self.fv * self.gu
        det_J13 = self.fu * self.hw - self.fw * self.hu
        d_j = (det_J12 * l_j * self.D_2 - det_J) / (l_j * (self.fu * l_j * self.D_2 - det_J13))
        return d_j

    def _direction_of_bifurcation(self, d_hat, d_2, l_j) -> bool:
        direction_branch = self._get_direction_of_branch(d_2)
        direction_DDI = self._get_direction_of_DDI(d_hat, l_j)
        # Check if direction of branch and direction of DDI in same direction
        if direction_branch * direction_DDI > 0:
            branch_DDI_same_direction = True
        elif direction_branch * direction_DDI < 0:
            branch_DDI_same_direction = False
        else:
            print('d_2 = 0 or no change in sign of det(J-l_j*D)!\n' +
                  # f'd_2 = {d_2} and derivative of det(J-l_j*D) in d = {derv_det_J_d_in_d}.\n'
                  'Either way this result not true, in first case need more expansion!')
            branch_DDI_same_direction = False

        print(f'DDI and branch in same direction: {branch_DDI_same_direction}')
        return branch_DDI_same_direction
    
    @staticmethod
    def _get_direction_of_branch(d_2):
        # We have branch of non-constant steady states for diffusion coefficient d(s) = d_hat + d_2*s^2 + ....
        # Hence, d(s) > d_hat or d(s) < d_hat for s small in case of d_2 \ne 0 depending on the sign of d_2.
        direction_branch = np.sign(d_2)
        if direction_branch == 0.0:
            print('d_2 = 0, i.e. we need to expand to higher order to get some results!')
        return direction_branch

    def _get_direction_of_DDI(self, d_hat, l_j):
        # The matrix J - l_j D changes stability at diffusion coefficient d_hat.
        # Either stable for d < d_hat and unstable for d > d_hat or the other way around
        # (or no change of sign in the determinant, i.e. zero eigenvalue in neighborhood of d_hat, but then no DDI).
        # For stable non-constant patterns, we need instability of J - l_j D.
        # The calculation depends on the case of bif. param, i.e. if d_hat = d_j or d_hat = d_j_plus
        if d_hat == self._get_bifurcation_p3(l_j):
            derv_p_3_in_d = (self.fu * self.hw - self.fw * self.hu) * l_j - self.fu * l_j**2 * self.D_2
            if np.sign(derv_p_3_in_d) > 0:
                # In this case det J - l_j D < 0 for d < d_hat and det J - l_j D > 0 for d > d_hat
                # This means, Turing pattern can only emerge for d < d_hat
                direction_DDI = -1.0
            elif np.sign(derv_p_3_in_d) < 0:
                # In this case det J - l_j D > 0 for d < d_hat and det J - l_j D < 0 for d > d_hat
                # This means, Turing pattern can only emerge for d > d_hat
                direction_DDI = 1.0
            else:
                direction_DDI = 0.0
                print('No change in sign of det(J-l_j*D)!')
        else:
            det_J12 = self.fu * self.gv - self.fv * self.gu
            tr_J13 = self.fu + self.hw
            det_J23 = self.gv * self.hw - self.gw * self.hv
            tr_J = self.fu + self.gv + self.hw
            lead_coeff_p1p2_p3 = l_j**3 - tr_J13 * l_j**2
            if lead_coeff_p1p2_p3 > 0.0:
                # In this case p_1*p_2-p_3 < 0 for d < d_hat and p_1*p_2-p_3 > 0 for d > d_hat
                # This means, Turing pattern can only emerge for d < d_hat
                direction_DDI = -1.0
            elif lead_coeff_p1p2_p3 < 0.0:
                # In this case p_1*p_2-p_3 > 0 for d < d_hat and p_1*p_2-p_3 < 0 for d > d_hat
                # This means, Turing pattern can only emerge for d > d_hat
                direction_DDI = +1.0
            else:
                second_coeff_p1p2_p3 = ((det_J12 + det_J23 + tr_J * tr_J13) * l_j
                                        - 2 * tr_J * self.D_2 * l_j**2 + self.D_2**2 * l_j**3)
                if second_coeff_p1p2_p3 > 0.0:
                    # In this case p_1*p_2-p_3 < 0 for d < d_hat and p_1*p_2-p_3 > 0 for d > d_hat
                    # This means, Turing pattern can only emerge for d < d_hat
                    direction_DDI = -1.0
                elif second_coeff_p1p2_p3 < 0.0:
                    # In this case p_1*p_2-p_3 > 0 for d < d_hat and p_1*p_2-p_3 < 0 for d > d_hat
                    # This means, Turing pattern can only emerge for d > d_hat
                    direction_DDI = +1.0
                else:
                    # p_1*p_2 - p_3 constant, so no change in sign!
                    direction_DDI = 0.0
                    print('No change in sign of p_1*p_2-p_3!')

        return direction_DDI

    def _calculate_sign_critical_eigv(self, d_hat, d_2, l_j):
        # If d_hat = d_j:
        # Theorem of Crandall-Rabinowitz: (sd'(s) * derv_mu_j(d_hat)) / crit_eigv(s) to -1 as s to 0,
        # where sd'(s) = 2*d_2*s^2 + 3*d_3*s^3 + ..., mu_j(d) eigenvalue of J-l_j*diag(d, D_2) and crit_eigv eigenvalue
        # of bifurcating non-constant steady states.
        if d_hat == self._get_bifurcation_p3(l_j):
            sign_s_derv_d_close_zero = np.sign(d_2)

            sign_derv_mu_j = np.sign(self._calc_derv_mu_j(l_j))

            sign_nenner = sign_s_derv_d_close_zero * sign_derv_mu_j
            sign_crit_eigv = -1.0 * sign_nenner
        else:
            sign_crit_eigv = None

        print(f'Sign critical eigenvalue: {sign_crit_eigv}')
        return sign_crit_eigv

    def _calc_derv_mu_j(self, l_j):
        det_J12 = self.fu * self.gv - self.fv * self.gu
        tr_J12 = self.fu + self.gv
        det_J13 = self.fu * self.hw - self.fw * self.hu
        tr_J13 = self.fu + self.hw
        det_J23 = self.gv * self.hw - self.gw * self.hv
        d_j = self._get_bifurcation_p3(l_j)
        q = det_J12 + det_J23 + det_J13 + (-tr_J13 * d_j - tr_J12 * self.D_2) * l_j + d_j * self.D_2 * l_j**2
        mu_derv = l_j * (self.fu * l_j * self.D_2 - det_J13) / q
        return mu_derv

    def _search_max_bifurcation_point_both(self, max_range):
        """
        Searches for the maximal bifurcation points within a given range.

        This method calculates the bifurcation points for two different cases:
        1. p_3 = 0
        2. p_1p_2 - p_3 = 0

        It iterates through a range of values, computes the bifurcation points,
        and identifies the maximum bifurcation points for both cases.

        Parameters:
        max_range (int): The upper limit of the range to search for bifurcation points.

        Returns:
        list: A list containing two sublists:
            - The first sublist contains the index and value of the maximal bifurcation point for p_3 = 0.
            - The second sublist contains the index and value of the maximal bifurcation point for p_1p_2 - p_3 = 0.

        Notes:
        - If the maximal bifurcation point is at the boundary of the search range, a warning message is printed.
        - If no real bifurcation points are found, the method returns NaN values.
        """
        bif_points_p3 = []
        bif_points_p1p2_p3 = []
        for j in range(1, max_range):
            l_j = j**2 * np.pi**2 / self.domain_size**2
            d_j = self._get_bifurcation_p3(l_j)
            d_j_plus, d_j_minus = self._get_bifurcation_p1p2_p3(l_j)

            bif_points_p3.append(d_j)
            if np.isreal(d_j_plus):
                bif_points_p1p2_p3.append(d_j_plus)
            else:
                bif_points_p1p2_p3.append(0)
        if bif_points_p3:
            d_p3_ind = np.argmax(bif_points_p3)
            d_p3 = bif_points_p3[d_p3_ind]
        else:
            d_p3_ind = -1
            d_p3 = 0
        if bif_points_p1p2_p3:
            d_p1p2_p3_ind = np.argmax(bif_points_p1p2_p3)
            d_p1p2_p3 = bif_points_p1p2_p3[d_p1p2_p3_ind]
        else:
            d_p1p2_p3_ind = -1
            d_p1p2_p3 = 0
        if d_p3_ind+2 == max_range or d_p1p2_p3_ind+2 == max_range:
            print('We don\'t look far enough for maximal bifurcation point!')
            if d_p3_ind+2 == max_range:
                d_p3_ind, d_p3 = np.nan, np.nan
            if d_p1p2_p3_ind+2 == max_range:
                d_p1p2_p3_ind, d_p1p2_p3 = np.nan, np.nan
            # d_p3_ind, d_p3, d_p1p2_p3_ind, d_p1p2_p3 = np.nan, np.nan, np.nan, np.nan
        print(f'Maximal bifurcation point of p_3=0: j = {d_p3_ind+1}, d_j = {d_p3}')
        print(f'Maximal bifurcation point of p_1p_2-p_3=0: j = {d_p1p2_p3_ind+1}, d_j_plus = {d_p1p2_p3}')
        return [d_p3_ind+1, d_p3], [d_p1p2_p3_ind+1, d_p1p2_p3]

    def _d_2_limit_wrong_version(self):
        # We know that l_max and d_max are constant w.r.t. domain length if we do not round the approximated j_max to an integer. (j_max is of course not constant but O(L^2))
        j_max, d_max = self.approx_max_wave_number_p3(rounded=False, both=False)
        l_max = j_max**2 * np.pi ** 2 / self.domain_size ** 2
        return self.calculate_d_2(d_max, l_max, faster=True)

    def _solve_f_for_u(self):
        x, y, z = sp.symbols('x y z')
        func = sp.sympify(self.nonlinear_f)
        p_symp = sp.solvers.solve(func, x)
        p_right = None
        if type(p_symp) == list:
            for p in p_symp:
                p_exec = sp.lambdify((y, z), p, 'numpy')
                if np.isclose(p_exec(self.v, self.w), self.u):
                    p_right = p
                    break
            if not p_right:
                print('Could not solve f for u, try numerically.')
        else:
            p_right = p_symp
        p_exec = sp.lambdify((y, z), p_right, 'numpy')
        return p_exec
