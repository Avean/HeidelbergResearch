import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np

u_bistable = np.linspace(4.1, 4.4, 100)
u_f_bistable = np.linspace(0, 1, 100)
v_g_nullcline_bistable = 1 * (u_bistable)**2
v_f_nullcline_bistable = - 2*(u_f_bistable)**2 + 13



fig = go.Figure()




fig.add_trace(go.Scatter(x=[1, 1], y=[0, 5.5], mode="lines",
                         line=dict(color='#D95319', width=8), showlegend=False))
fig.add_trace(go.Scatter(x=[1, 13], y=[5.5, 5.5], mode="lines",
                         line=dict(color='#D95319', width=8), showlegend=False))
fig.add_trace(go.Scatter(x=[1, 13], y=[0, 0], mode="lines",
                         line=dict(color='#D95319', width=8), showlegend=False))
fig.add_trace(go.Scatter(x=[13, 13], y=[0, 5.5], mode="lines",
                         line=dict(color='#D3D3D3', width=5, dash='dash'), showlegend=False))


fig.add_trace(go.Scatter(x=v_g_nullcline_bistable-5, y=u_bistable+0.7, line=dict(color='#0072BD', width=6), legendgroup='1', showlegend=False))

fig.add_trace(go.Scatter(x=v_f_nullcline_bistable, y=u_f_bistable, line=dict(color='#0072BD', width=6), legendgroup='2', showlegend=False))


fig.update_layout(
        font_family="Times New Roman",
        # legend_font_family="Times New Roman",
        font_size=42,
        plot_bgcolor="rgba(0, 0, 0, 0)",
        paper_bgcolor="rgba(0, 0, 0, 0)",
        
        )

fig.update_annotations(font_size=45)


fig.update_xaxes(range=[0,18], showline=True, linewidth=2, linecolor='black', mirror=False,
                # showgrid=True, gridcolor='rgba(0, 0, 0, 0.13)', gridwidth=1,
                zeroline=True, zerolinecolor='black', zerolinewidth=2,
                title_text=r'$\Huge \kappa$',
                # tickvals=np.linspace(0,18, 37), ticktext=[0, '' ,r'$\huge\kappa_\min$']+['' for _ in range(23)]+[r'$\huge\kappa_1$']+['' for _ in range(10)],
                tickvals=[0, 1, 13], ticktext=[0, r'$\huge\kappa_\min$', r'$\huge\kappa_1$'],
                ticks='outside', ticklen=5, tickwidth=2
                )
fig.update_yaxes(range=[0,6.1], showline=True, linewidth=2, linecolor='black', mirror=False,
                # showgrid=True, gridcolor='rgba(0, 0, 0, 0.13)', gridwidth=1,
                zeroline=True, zerolinecolor='black', zerolinewidth=2,
                title_text=r'$\Huge \|\cdot\|$',
                # tickvals=np.linspace(0,6.5,14), ticktext=[0] + ['' for _ in range(10)] + [r'$\huge \hat C$', '', ''],
                tickvals=[0, 5.5], ticktext=[0, r'$\huge \hat C$'],
                ticks='outside', ticklen=5, tickwidth=2
                )


fig.show()