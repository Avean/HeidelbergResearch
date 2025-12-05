import numpy as np

import models

from properties_bifurcating_branch import TuringBifurcationThreeEqu
from surpress_print import suppress_print


if __name__ == '__main__':

    [func_f, func_g, func_h], initial_guess_steady_state = models.DKK_model_three_equs()

    length_domain = 1.0
    Diff_2 = 1.0 #0.017 #1.14056387

    turing_bif = TuringBifurcationThreeEqu(func_f, func_g, func_h, Diff_2, length_domain,
                                           initial_guess_steady_state)

    
    print('\n')
    if not turing_bif.check_DDI_possible():
        raise SystemExit(0)
    
    print('\n')
    turing_bif.search_max_bifurcation_point()
    #print('Result of approx max wave number' + str(turing_bif.approx_max_wave_number_p3(rounded=True, both=False)))

    with suppress_print():
        j_max, max_bif_point_p3 = turing_bif.approx_max_wave_number_p3(rounded=True, both=False)
        _, d_max_searched = turing_bif.search_max_bifurcation_point(max_range=20+2*j_max)
    print(f'Max bif point p3: {max_bif_point_p3}')
    print(f'D_max searched: {d_max_searched}')
    print(max_bif_point_p3 < d_max_searched - 1e-3)


    # print(f'Approximation of d_2 and limit w.r.t. domain size: {turing_bif.d_2_limit_wrt_domain_size()}')
    with suppress_print():
        B_k_max_bif_point = 2*(j_max**2 * np.pi**2 / length_domain**2)*turing_bif.calculate_d_2(d_hat=max_bif_point_p3, l_j=j_max**2 * np.pi**2 / length_domain**2)
    print(f'B(k) at max bif point: {B_k_max_bif_point}')


    print('\n')
    # max_range = How far to look for maximal bifurcation point:
    # Search for max bif point for wave numbers j in range(1, max_range).
    turing_bif.check_Turing_patterns_sure(max_range=801)

    # Calculate 2*l_j / C^2 d_2 = A_1(2*eta_1+xi_1) + A_2(2*eta_2+xi_2) + 3/2 A_3 for bif point: 
    # print('\n')
    # print('Calculation of 2*l_j / C^2 d_2 for max bifurcation point: ')
    # l_j = j_max**2 * np.pi**2 / length_domain**2
    # print(2*l_j*turing_bif.calculate_d_2(d_hat=max_bif_point_p3, l_j=l_j))
    # print(2*l_j*turing_bif.calc_d_2_at_max_bif())

    # print('\n')
    # turing_bif.bifurcation_at_max_wave_num_p3()
    # print('\n')
    # turing_bif.turing_bifurcation(4)
    # print('\n')
    # turing_bif.turing_bifurcation(6)

    # turing_bif.plot_d_2_for_multiple_wave_numbers(range(1,31))
    turing_bif.plot_B_j_for_multiple_wave_numbers(range(1, 31))
    # turing_bif.plot_B_k_for_multiple_domain_sizes(np.linspace(0.01, 10, 1000))

    turing_bif.plot_all_possible_bifurcation_points(range_wavenumber=range(1, 31))

    # turing_bif.plot_bif_point_d_2_for_mult_wave_numbers(range_wavenumber=range(1, 21))
