import plotly.graph_objects as go
import numpy as np


def read_data(name_file):
    with open(name_file, 'r') as f:
        content = f.read()
    
    blocks = content.strip().split('\n\n[')
    kappa_list = []
    data_list = []

    for block in blocks:
        kappa_string, d_string = block.split(',\n', 1)
        kappa_list.extend([float(x) for x in kappa_string.strip('[]').split()])
        data_list.extend([float(x.split('(')[1].split(')')[0]) for x in d_string.strip('[]').split(', ')])

    return kappa_list, data_list

# Plot of D_min and D_max and numerical threshold for kappa in [0, 1.5], no log scale for any axis
def plot_numerical_threshold_0_1():
    kappa_thresh_list, num_thresh_list = read_data('numerical_calc_thresholds.txt')
    kappa_D_min_list, D_min_list = read_data('D_min.txt')
    fig2 = go.Figure()
    fig2.add_trace(go.Scatter(x=kappa_thresh_list[0:49], y=num_thresh_list[0:49], mode='lines', line=dict(color='blue', width=4), name='Threshold for pattern')) #r'$\large \text{Threshold (numeric.)}$'))
    fig2.add_trace(go.Scatter(x=kappa_D_min_list[0:49], y=D_min_list[0:49], mode='lines', line=dict(color='green', width=4), name=r'$\large D_{min}$'))

    kappa_stab_const = np.linspace(1, 1.5, 100)
    stab_const = [(kappa-1)/(4*np.pi**2) for kappa in kappa_stab_const]
    fig2.add_trace(go.Scatter(x=kappa_stab_const, y=stab_const, mode='lines', line=dict(color='black', width=6, dash='dot'), name='Stability of constant steady state'))

    #fig2.add_trace(go.Scatter(x=kappa_D_max_list, y=D_max_list, mode='lines', name='D_max'))
    #fig2.add_trace(go.Scatter(x=kappa_D_max_linear, y=D_max_linear, mode='lines', name='D_max_linear'))

    fig2.update_layout(title='Difference of numerical threshold and D_min for kappa in [0, 1.5]',
                        xaxis_title=r'$\Large\kappa$',
                        yaxis_title=r'$\Large\text{Diffusion } D$')

    fig2.update_layout(
        font_size=24,
        legend=dict(
            # font=dict(size = 20),
            # itemsizing="constant",
            yanchor="top", yref='paper', y=0.99,
            xanchor="left", xref='paper', x=0.005,
            bordercolor="Black",
            borderwidth=3,
            bgcolor="whitesmoke",
            # tracegroupgap=10,
            #entrywidthmode="pixels",
            # itemwidth=50,
            ),
        plot_bgcolor='rgba(0, 0, 0, 0)',
        paper_bgcolor='rgba(0, 0, 0, 0)',
        )
    fig2.update_xaxes(range=[0,1.5], 
                      showline=True, linewidth=2, linecolor='black', mirror=True, showgrid=True, gridwidth=1, gridcolor='LightGray', zeroline=True, zerolinewidth=1, zerolinecolor='LightGray')
    fig2.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True, showgrid=True, gridwidth=1, gridcolor='LightGray', zeroline=True, zerolinewidth=1, zerolinecolor='LightGray')
    
    fig2.show()

# Plot of D_min and D_max and linear estimate of D_max for kappa in [0, 10]
def plot_analytic_thresholds():
    #kappa_thresh_list, num_thresh_list = read_data('numerical_calc_thresholds.txt')
    kappa_D_min_list, D_min_list = read_data('D_min.txt')
    kappa_D_max_list, D_max_list = read_data('D_max.txt')
    kappa_D_max_linear, D_max_linear = read_data('D_max_lin.txt')
    fig = go.Figure()
    #fig.add_trace(go.Scatter(x=kappa_thresh_list[0:63], y=num_thresh_list[0:63], mode='lines', line=dict(color='blue', width=4), name='Threshold (numerical)'))
    fig.add_trace(go.Scatter(x=kappa_D_min_list[0:101], y=D_min_list[0:101], mode='lines', line=dict(width=2, color='red'), name=r'$\large D_{min}$'))
    fig.add_trace(go.Scatter(x=kappa_D_max_list[0:101], y=D_max_list[0:101], mode='lines', line=dict(width=2, color='green'), name=r'$\large D_{max}$'))
    fig.add_trace(go.Scatter(x=kappa_D_max_linear[0:101], y=D_max_linear[0:101], mode='lines', line=dict(width=2, color='purple'), name=r'$\large \text{Linear estimate of} D_{max}$'))

    fig.update_layout(title='Analytic Thresholds',
                        xaxis_title=r'$\Large\kappa$',
                        yaxis_title=r'$\Large\text{Diffusion } D$')

    # fig.add_trace(go.Scatter(x=[0, 0], y=[0, 0.018], mode='lines', line=dict(color='gray', width=2, dash='dash'), showlegend=False))
    # fig.add_trace(go.Scatter(x=[1.6, 1.6], y=[0, 0.018], mode='lines', line=dict(color='gray', width=2, dash='dash'), showlegend=False))
    # fig.add_trace(go.Scatter(x=[0, 1.6], y=[0, 0], mode='lines', line=dict(color='gray', width=2, dash='dash'), showlegend=False))
    # fig.add_trace(go.Scatter(x=[0, 1.6], y=[0.018, 0.018], mode='lines', line=dict(color='gray', width=2, dash='dash'), showlegend=False))

    fig.update_layout(
        font_size=24,
        legend=dict(
            # font=dict(size = 20),
            # itemsizing="constant",
            # yanchor="top", yref='paper', y=1,
            # xanchor="left", xref='paper', x=1.01,
            bordercolor="Black",
            borderwidth=2,
            # tracegroupgap=10,
            #entrywidthmode="pixels",
            # itemwidth=50,
            ),
        plot_bgcolor='rgba(0, 0, 0, 0)',
        paper_bgcolor='rgba(0, 0, 0, 0)',
        )
    fig.update_xaxes(showline=True, linewidth=2, linecolor='black', mirror=True, showgrid=True, gridwidth=1, gridcolor='LightGray', zeroline=True, zerolinewidth=1, zerolinecolor='LightGray')
    fig.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True, showgrid=True, gridwidth=1, gridcolor='LightGray', zeroline=True, zerolinewidth=1, zerolinecolor='LightGray')
    
    fig.show()

# Plot of D_min and D_max and numerical threshold for kappa in [0, 10], where y_axis = diffusion and x_axis = kappa both log scale
def plot_numerical_threshold_comparison():
    kappa_thresh_list, num_thresh_list = read_data('numerical_calc_thresholds_log.txt')
    kappa_D_min_list, D_min_list = read_data('D_min_log.txt')
    kappa_D_max_list, D_max_list = read_data('D_max_log.txt')
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=kappa_thresh_list[19:94], y=num_thresh_list[19:94], mode='lines', line=dict(color='blue', width=2), name='Threshold (numeric.)'))
    fig.add_trace(go.Scatter(x=kappa_D_min_list[19:117], y=D_min_list[19:117], mode='lines', line=dict(color='red', width=2), name=r'$\large D_{min}$'))
    fig.add_trace(go.Scatter(x=kappa_D_max_list[19:100], y=D_max_list[19:100], mode='lines', line=dict(color='green', width=2), name=r'$\large D_{max}$'))
    # fig.add_trace(go.Scatter(x=kappa_thresh_list[0:49], y=num_thresh_list[0:49], mode='lines', line=dict(color='blue', width=2), name='Threshold (numeric.)'))
    # fig.add_trace(go.Scatter(x=kappa_D_min_list[0:50], y=D_min_list[0:50], mode='lines', line=dict(color='red', width=2), name=r'$\large D_{min}$'))
    # fig.add_trace(go.Scatter(x=kappa_D_max_list[0:50], y=D_max_list[0:50], mode='lines', line=dict(color='green', width=2), name=r'$\large D_{max}$'))
    fig.update_layout(title='Comparison numerical and analytical Thresholds',
                        xaxis_title=r'$\Large\kappa$',
                        yaxis_title=r'$\Large\text{Diffusion } D$')

    fig.add_trace(go.Scatter(x=[0, 0], y=[0, 0.018], mode='lines', line=dict(color='gray', width=2, dash='dash'), showlegend=False))
    fig.add_trace(go.Scatter(x=[1.6, 1.6], y=[0, 0.018], mode='lines', line=dict(color='gray', width=2, dash='dash'), showlegend=False))
    fig.add_trace(go.Scatter(x=[0, 1.6], y=[0, 0], mode='lines', line=dict(color='gray', width=2, dash='dash'), showlegend=False))
    fig.add_trace(go.Scatter(x=[0, 1.6], y=[0.018, 0.018], mode='lines', line=dict(color='gray', width=2, dash='dash'), showlegend=False))

    fig.update_layout(
        font_size=24,
        legend=dict(
            # font=dict(size = 20),
            # itemsizing="constant",
            # yanchor="top", yref='paper', y=1,
            # xanchor="left", xref='paper', x=1.01,
            bordercolor="Black",
            borderwidth=2,
            # tracegroupgap=10,
            #entrywidthmode="pixels",
            # itemwidth=50,
            ),
        plot_bgcolor='rgba(0, 0, 0, 0)',
        paper_bgcolor='rgba(0, 0, 0, 0)',
        )
    fig.update_xaxes(showline=True, linewidth=2, linecolor='black', mirror=True, showgrid=True, gridwidth=1, gridcolor='LightGray', zeroline=True, zerolinewidth=1, zerolinecolor='LightGray', type="log")
    fig.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True, showgrid=True, gridwidth=1, gridcolor='LightGray', zeroline=True, zerolinewidth=1, zerolinecolor='LightGray', type="log")
    
    fig.show()

# Plot of D_min and D_max and numerical threshold for kappa in [0, 10], where y_axis = diffusion and x_axis = kappa both log scale
# Additionally, areas in plot have different colors
def plot_with_areas():
    kappa_thresh_list, num_thresh_list = read_data('numerical_calc_thresholds_log.txt')
    kappa_D_min_list, D_min_list = read_data('D_min_log.txt')
    kappa_D_max_list, D_max_list = read_data('D_max_log.txt')

    fig = go.Figure()
    
    fig.add_trace(go.Scatter(x=kappa_thresh_list[19:94], y=num_thresh_list[19:94], mode='lines', line=dict(color='white', width=0), name='Nonconstant steady states', fill='tozeroy', fillcolor='rgba(176, 196, 222, 0.85)', legendgroup=0))  #lightgray
    fig.add_trace(go.Scatter(x=[0.1, 10], y=[10,10], mode='lines', line=dict(color='white', width=0), showlegend=False))
    fig.add_trace(go.Scatter(x=kappa_thresh_list[19:94], y=num_thresh_list[19:94], mode='lines', line=dict(color='white', width=0), name='Only constant steady state', fill='tonexty', fillcolor='rgba(255, 222, 173, 0.85)', legendgroup=1)) # peachpuff, navajowhite, whitesmoke, lightgoldenrodyellow, moccasin
    fig.add_trace(go.Scatter(x=kappa_D_min_list[19:117], y=D_min_list[19:117], mode='lines', line=dict(color='black', width=1), name=r'$\large D_{min}$', fill='tozeroy', fillpattern=dict(fgcolor='black', fillmode='replace', shape="/", size=10), showlegend=False)) #, fill='tozeroy'
    fig.add_trace(go.Scatter(x=[0.1, 10], y=[10,10], mode='lines', line=dict(color='white', width=0), showlegend=False))
    fig.add_trace(go.Scatter(x=kappa_D_max_list[19:100], y=D_max_list[19:100], mode='lines', line=dict(color='black', width=1), name='proven analytically', fill='tonexty', fillpattern=dict(fgcolor='black', fillmode='replace', shape="/", size=10), legendgroup=2)) # , fill='tonexty'
    
    kappa_stab_const = np.logspace(0, 1, 10000)
    stab_const = [(kappa-1)/(4*np.pi**2) for kappa in kappa_stab_const]
    fig.add_trace(go.Scatter(x=kappa_stab_const, y=stab_const, mode='lines', line=dict(color='black', width=6, dash='dashdot'), name='Stability constant steady state', legendgroup=3))


    fig.update_layout(#title='Comparison numerical and analytical thresholds for pattern formation',
                        xaxis_title=r'$\huge\kappa$',
                        yaxis_title=r'$\huge\text{Diffusion}$')
    
    fig.update_layout(
        font_size=40,
        legend=dict(
            # font=dict(size = 20),
            # itemsizing="constant",
            yanchor="bottom", yref='paper', y=0.01,
            xanchor="right", xref='paper', x=0.995,
            bordercolor="Black",
            borderwidth=5,
            bgcolor="whitesmoke",
            # tracegroupgap=10,
            # entrywidthmode="pixels",
            # itemwidth=50,
            # traceorder="reversed",
            ),
        # plot_bgcolor='rgba(0, 0, 0, 0)',
        # paper_bgcolor='rgba(0, 0, 0, 0)',
        )
    fig.update_xaxes(showline=True, linewidth=2, linecolor='black', mirror=True, type="log",
                     showgrid=True, gridcolor='gray', gridwidth=2,
                     tickmode = 'array', ticks="outside", tickwidth=3, ticklen=6,
                     tickfont=dict(size=48), tickvals=np.logspace(np.log(0.1), np.log(10), 19, base=np.e), ticktext=[0.1, '', '', '', '', '', '', '', '', 1, '', '', '', '', '', '', '', '', 10], nticks=10)
    fig.update_yaxes(range=[-7, 1], tickformat=".0e",
                     showline=True, linewidth=2, linecolor='black', mirror=True, type="log",
                     showgrid=True, gridcolor='gray', gridwidth=2,
                     tickmode = 'array', ticks="outside", tickwidth=3, ticklen=6,
                     tickfont=dict(size=48), tickvals=[1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e+1], ticktext=['1e-7', '', '1e-5', '', '1e-3', '', '1e-1', '', '1e+1'])
    
    fig.show()

# plot_analytic_thresholds()
# plot_numerical_threshold_comparison()
# plot_numerical_threshold_0_1()
plot_with_areas()
