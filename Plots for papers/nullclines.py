import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np

u_bistable = np.linspace(0, 6.2, 100)
u_f_bistable = np.linspace(0, 5.8, 100)
v_g_nullcline_bistable = 4 * u_bistable
v_f_nullcline_bistable = 0.5*(u_f_bistable - 2)**3 - (u_f_bistable - 2)**2 + (u_f_bistable - 2) + 10

# fig = make_subplots(rows=1, cols=2, subplot_titles=("Bistability Case", "Hysteresis Case"))

# fig.add_trace(go.Scatter(x=v_g_nullcline_bistable, y=u_bistable, line=dict(color='blue', width=3)), row=1, col=1)
# fig.add_trace(go.Scatter(x=v_f_nullcline_bistable, y=u_f_bistable, line=dict(color='red', width=3)), row=1, col=1)
# fig.add_trace(go.Scatter(x=[0], y=[0], mode="markers+text", textposition="top right",
#                          text=[r'$\huge S_0$'],
#                          line=dict(color='black', width=3), showlegend=False),
#                          row=1, col=1)
# fig.add_trace(go.Scatter(x=[10.3431], y=[2.58579], mode="markers+text", textposition="middle right",
#                          text=[r'$\huge S_1$'],
#                          line=dict(color='black', width=3), showlegend=False),
#                          row=1, col=1)
# fig.add_trace(go.Scatter(x=[21.6568], y=[5.41421], mode="markers+text", textposition="top left",
#                          text=[r'$\huge S_2$'],
#                          line=dict(color='black', width=3), showlegend=False),
#                          row=1, col=1)

# fig.update_xaxes(range=[0,25], showline=True, linewidth=2, linecolor='black', mirror=False,
#                 # showgrid=True, gridcolor='gray', gridwidth=2,
#                 zeroline=True, zerolinecolor='black', zerolinewidth=2,
#                 title_text=r'$\huge v$',
#                 # tickvals=[0],
#                 row=1, col=1
#                 )
# fig.update_yaxes(range=[0,6], showline=True, linewidth=2, linecolor='black', mirror=False,
#                 # showgrid=True, gridcolor='gray', gridwidth=2,
#                 zeroline=True, zerolinecolor='black', zerolinewidth=2,
#                 title_text=r'$\huge u$',
#                 # tickvals=[0],
#                 row=1, col=1
#                 )


u_hysteresis = np.linspace(0, 10.32, 100)
u_f_hysteresis = np.linspace(0, 9, 100)
v_g_nullcline_hysteresis = 2 * u_hysteresis
v_f_nullcline_hysteresis = 0.3*(u_f_hysteresis-2)**3 - 2*(u_f_hysteresis - 2)**2 + 0.6*(u_f_hysteresis - 2) + 11.6

# fig.add_trace(go.Scatter(x=v_g_nullcline_hysteresis, y=u_hysteresis, line=dict(color='blue', width=3)), row=1, col=2)
# fig.add_trace(go.Scatter(x=v_f_nullcline_hysteresis, y=u_f_hysteresis, line=dict(color='red', width=3)), row=1, col=2)
# fig.add_trace(go.Scatter(x=[0], y=[0], mode="markers+text", textposition="top right",
#                          text=[r'$\huge S_0$'],
#                          line=dict(color='black', width=3), showlegend=False),
#                          row=1, col=2)
# fig.add_trace(go.Scatter(x=[7.72252], y=[3.86127], mode="markers+text", textposition="top right",
#                          text=[r'$\huge S_1$'],
#                          line=dict(color='black', width=3), showlegend=False),
#                          row=1, col=2)
# fig.add_trace(go.Scatter(x=[17.6108], y=[8.8054], mode="markers+text", textposition="top right",
#                          text=[r'$\huge S_2$'],
#                          line=dict(color='black', width=3), showlegend=False),
#                          row=1, col=2)
# fig.add_trace(go.Scatter(x=[11.6461, 1.05187], y=[2.155, 6.289],
#                         mode="markers+text", textposition="middle right",
#                         text=[r'$\huge(v_H, u_H)$', r'$\huge(v_T, u_T)$'],
#                         line=dict(color='black', width=3), showlegend=False),
#                         row=1, col=2)

# fig.update_xaxes(range=[0,20], showline=True, linewidth=2, linecolor='black', mirror=False,
#                 # showgrid=True, gridcolor='gray', gridwidth=2,
#                 zeroline=True, zerolinecolor='black', zerolinewidth=2,
#                 title_text=r'$\huge v$',
#                 # tickvals=[0],
#                 row=1, col=2
#                 )
# fig.update_yaxes(range=[0,10], showline=True, linewidth=2, linecolor='black', mirror=False,
#                 # showgrid=True, gridcolor='gray', gridwidth=2,
#                 zeroline=True, zerolinecolor='black', zerolinewidth=2,
#                 title_text=r'$\huge u$',
#                 # tickvals=[0],
#                 row=1, col=2
#                 )

# fig.show()


fig = make_subplots(rows=1, cols=2, subplot_titles=("Bistable Case", "Hysteresis Case"))



fig.add_trace(go.Scatter(x=v_g_nullcline_bistable, y=u_bistable, line=dict(color='blue', width=3), name='g(u,v) = 0', legendgroup='1'), #name=r'$\huge g(u,v) = 0$'
              row=1, col=1
              )

fig.add_trace(go.Scatter(x=v_f_nullcline_bistable, y=u_f_bistable, line=dict(color='red', width=3), name='f(u,v) = 0', legendgroup='2'),
              row=1, col=1
              )


fig.add_trace(go.Scatter(x=[0], y=[0], mode="markers", textposition="top right",
                         text=[r'$\Huge S_0$'],
                         marker=dict(color='black', size=13), showlegend=False),
                         row=1, col=1)
fig.add_trace(go.Scatter(x=[10.3431], y=[2.58579], mode="markers", textposition="middle right",
                         text=[r'$\Huge S_1$'],
                         marker=dict(color='black', size=10), showlegend=False),
                         row=1, col=1)
fig.add_trace(go.Scatter(x=[21.6568], y=[5.41421], mode="markers", textposition="top left",
                         text=[r'$\Huge S_2$'],
                         marker=dict(color='black', size=10), showlegend=False),
                         row=1, col=1)


x_data=[0, 10.3431, 21.6568]
y_data=[0, 2.58579, 5.41421]
text_data=[r'$\Huge S_0$', r'$\Huge S_1$', r'$\Huge S_2$',]

data_length = (len(x_data))

xshift = [25, 33, -35] #Position of the annotation (you can also create this for yshift) 
yshift = [50, -3, 31] * data_length #Position of the annotation (you can also create this for yshift) 
showarrow = [False] * data_length #No arrow
font_size = [25] * data_length #Annotation Fontsize

annotations = [dict(x=x, y=y, text=t, showarrow=i, xshift=j, yshift=l, font_size=k)
               for x, y, t, i, j, l, k in zip(x_data, y_data, text_data, showarrow, xshift, yshift, font_size)]

for anno in annotations:
    fig.add_annotation(anno, row=1, col=1)



fig.add_trace(go.Scatter(x=v_g_nullcline_hysteresis, y=u_hysteresis, line=dict(color='blue', width=3), showlegend=False, legendgroup='3'),
              row=1, col=2
              )

fig.add_trace(go.Scatter(x=v_f_nullcline_hysteresis, y=u_f_hysteresis, line=dict(color='red', width=3), showlegend=False, legendgroup='3'),
              row=1, col=2
              )


fig.add_trace(go.Scatter(x=[0], y=[0], mode="markers", textposition="top right",
                         text=[r'$\Huge S_0$'],
                         marker=dict(color='black', size=13), showlegend=False),
                         row=1, col=2)
fig.add_trace(go.Scatter(x=[7.72252], y=[3.86127], mode="markers", textposition="bottom center",
                         text=[r'$\Huge S_1$'],
                         marker=dict(color='black', size=10), showlegend=False),
                         row=1, col=2)
fig.add_trace(go.Scatter(x=[17.6108], y=[8.8054], mode="markers", textposition="top left",
                         text=[r'$\Huge S_2$'],
                         marker=dict(color='black', size=10), showlegend=False),
                         row=1, col=2)
fig.add_trace(go.Scatter(x=[11.6461, 1.05187], y=[2.155, 6.289],
                        mode="markers", textposition="middle right",
                        text=[r'$\Huge(v_H, u_H)$', r'$\Huge(v_T, u_T)$'],
                        marker=dict(color='black', size=10), showlegend=False),
                        row=1, col=2)


x_data=[11.6461, 1.05187, 0, 7.72252, 17.6108]
y_data=[2.155, 6.289, 0, 3.86127, 8.8054]
text_data=[r'$\Huge(v_H, u_H)$', r'$\Huge(v_T, u_T)$', r'$\Huge S_0$', r'$\Huge S_1$', r'$\Huge S_2$',]

data_length = (len(x_data))

xshift = [98, 89, 25, 5, -35] #Position of the annotation (you can also create this for yshift) 
yshift = [2, 2, 50, -35, 31] * data_length #Position of the annotation (you can also create this for yshift) 
showarrow = [False] * data_length #No arrow
font_size = [25] * data_length #Annotation Fontsize

annotations = [dict(x=x, y=y, text=t, showarrow=i, xshift=j, yshift=l, font_size=k)
               for x, y, t, i, j, l, k in zip(x_data, y_data, text_data, showarrow, xshift, yshift, font_size)]

for anno in annotations:
    fig.add_annotation(anno, row=1, col=2)

fig.update_layout(
        font_family="Times New Roman",
        # legend_font_family="Times New Roman",
        font_size=25,
        plot_bgcolor="rgba(0, 0, 0, 0)",
        paper_bgcolor="rgba(0, 0, 0, 0)",
        legend=dict(
            bordercolor="Black",
            borderwidth=2,
            bgcolor="rgba(0, 0, 0, 0)",
            font_size=40,
            ),
        
        )

fig.update_annotations(font_size=45)


fig.update_xaxes(range=[0,25], showline=True, linewidth=2, linecolor='black', mirror=False,
                # showgrid=True, gridcolor='gray', gridwidth=2,
                zeroline=True, zerolinecolor='black', zerolinewidth=2,
                title_text=r'$\Huge v$',
                tickvals=[0],
                row=1, col=1
                )
fig.update_yaxes(range=[0,6], showline=True, linewidth=2, linecolor='black', mirror=False,
                # showgrid=True, gridcolor='gray', gridwidth=2,
                zeroline=True, zerolinecolor='black', zerolinewidth=2,
                title_text=r'$\Huge u$',
                tickvals=[0],
                row=1, col=1
                )

fig.update_xaxes(range=[0,20], showline=True, linewidth=2, linecolor='black', mirror=False,
                # showgrid=True, gridcolor='gray', gridwidth=2,
                zeroline=True, zerolinecolor='black', zerolinewidth=2,
                title_text=r'$\Huge v$',
                tickvals=[0],
                row=1, col=2
                )
fig.update_yaxes(range=[0,10], showline=True, linewidth=2, linecolor='black', mirror=False,
                # showgrid=True, gridcolor='gray', gridwidth=2,
                zeroline=True, zerolinecolor='black', zerolinewidth=2,
                title_text=r'$\Huge u$',
                tickvals=[0],
                row=1, col=2
                )


fig.show()