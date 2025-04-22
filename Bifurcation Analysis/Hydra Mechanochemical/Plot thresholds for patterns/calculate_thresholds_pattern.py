import numpy as np
# import plotly.express as px
import plotly.graph_objects as go
# from plotly.subplots import make_subplots
# import pandas as pd

from minimize_functional import minimize_functional, calc_D_min_v2, calc_D_max_lin, calc_D_max_v2, EPS


def obtain_threshold_numerically(kappa):
    # Calculate D_min and D_max:
    D_min = calc_D_min_v2(kappa)
    D_max = calc_D_max_v2(kappa) # 0.1 # calc_D_max(kappa)

    # Calculate threshold:
    while D_max - D_min > EPS:
        d_mid = (D_min + D_max) / 2
        _, _, nonconstant_minima_found = minimize_functional(kappa, d_mid)
        if nonconstant_minima_found:
            D_min = d_mid
        else:
            D_max = d_mid
    d_mid = (D_min + D_max) / 2
    print(f'Kappa: {kappa}, Threshold: {d_mid}')
    return d_mid

def plot_numerical_thresholds():
    # kappa_list = np.linspace(0.001, 1.5, 30)
    kappa_list = np.logspace(-1, 0, 20, base=10)
    threshold_list = [obtain_threshold_numerically(k) for k in kappa_list]

    f = open('numerical_calc_thresholds_2.txt', 'x') # Save thresholds to file
    f.write(f'{kappa_list},\n{threshold_list}\n')
    f.write('\n\n')
    f.close()

    D_min_list = [calc_D_min_v2(k) for k in kappa_list]
    D_max_linear = [calc_D_max_lin(k) for k in kappa_list]
    D_max_list = [calc_D_max_v2(k) for k in kappa_list]

    g = open('D_min_2.txt', 'x')
    g.write(f'{kappa_list},\n{D_min_list}\n')
    g.write('\n\n')
    g.close()

    h = open('D_max_2.txt', 'x')
    h.write(f'{kappa_list},\n{D_max_list}\n')
    h.write('\n\n')
    h.close()

    k = open('D_max_lin_2.txt', 'x')
    k.write(f'{kappa_list},\n{D_max_linear}\n')
    k.write('\n\n')
    k.close()

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=kappa_list, y=threshold_list, mode='lines', name='Threshold (numerical)'))
    fig.add_trace(go.Scatter(x=kappa_list, y=D_min_list, mode='lines', name='D_min'))
    fig.add_trace(go.Scatter(x=kappa_list, y=D_max_list, mode='lines', name='D_max'))
    fig.add_trace(go.Scatter(x=kappa_list, y=D_max_linear, mode='lines', name='D_max_linear'))
    fig.update_layout(title='Numerical Thresholds',
                      xaxis_title='Kappa',
                      yaxis_title='Diffusion')
    fig.show()

if __name__ == '__main__':
    # Plot numerical thresholds and save them to text files
    # Text files can be used for plotting with plotting.py
    plot_numerical_thresholds()