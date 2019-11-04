from __future__ import print_function
from __future__ import division
from builtins import zip
from builtins import map
from builtins import range
import matplotlib.pyplot as plt
import matplotlib.path as path
import numpy as np

from gpkit.small_scripts import mag

import chart_studio.plotly as py
import plotly.graph_objs as go

def forceAspect(ax, aspect=1):
    im = ax.get_images()
    extent = im[0].get_extent()
    ax.set_aspect(abs((extent[1] - extent[0]) / (extent[3] - extent[2])) / aspect)

def draw_KT_network(sol, coordinates, topology_list):
    coordinate_list = [[coordinates[i[0]], coordinates[i[1]]] for i in topology_list]
    orig, dest = zip(*coordinate_list)
    xmin, ymin, xmax, ymax = 1e10 * np.array([1, 1, -1, -1])
    # Arrow parameters for flow plotting
    hl = 400  # Arrow head length
    try:
        hwf = hl * mag(sol('q')) / max(mag(sol('q')))  # Arrow width for flow solution
        hwd = hl * mag(sol('D')) / max(mag(sol('D')))  # Arrow width for diameter solution
        pwSr = 15*mag(sol('\dot{V}_+')) / max(mag(sol('\dot{V}_+')))  # Point width for sources
        pwSk = 15*mag(sol('\dot{V}_-')) / max(mag(sol('\dot{V}_-')))  # Point width for sinks
    except:
        hwf = hl * mag(sol['q']) / max(mag(sol['q']))  # Arrow width for flow solution
        hwd = hl * mag(sol['D']) / max(mag(sol['D']))  # Arrow width for diameter solution
        pwSr = 15*mag(sol['\dot{V}_+']) / max(mag(sol['\dot{V}_+']))  # Point width for sources
        pwSk = 15*mag(sol['\dot{V}_-']) / max(mag(sol['\dot{V}_-']))  # Point width for sinks
    fig = plt.figure()
    ax1 = fig.add_subplot(1, 2, 1, aspect=1)
    ax2 = fig.add_subplot(1, 2, 2, aspect=1)
    for i in range(len(topology_list)):
        xmin = np.min([xmin, orig[i][0], dest[i][0]])
        xmax = np.max([xmax, orig[i][0], dest[i][0]])
        ymin = np.min([ymin, orig[i][1], dest[i][1]])
        ymax = np.max([ymax, orig[i][1], dest[i][1]])
        # Plotting flow arrows
        ax1.arrow(orig[i][0], orig[i][1], dest[i][0] - orig[i][0], dest[i][1] - orig[i][1], label=i, color='b',
                  width=1. / 3. * hwf[i], head_width=hwf[i], head_length=hl,
                  length_includes_head=True)
        # Plotting diameter arrows
        ax2.arrow(orig[i][0], orig[i][1], dest[i][0] - orig[i][0], dest[i][1] - orig[i][1], label=i, color='r',
                  width=1. / 3. * hwd[i], head_width=hwd[i], head_length=hl,
                  length_includes_head=True)

    for i in coordinates.keys():
        # Plotting sources and sinks
        ax1.plot(coordinates[i][0], coordinates[i][1], 'o', color='b', mfc='none', markersize=pwSr[i])
        ax1.plot(coordinates[i][0], coordinates[i][1], 'o', color='r', mfc='none', markersize=pwSk[i])
        ax1.annotate(i,(coordinates[i][0], coordinates[i][1]))
    ax1.set_xlim([xmin - hl, xmax + hl])
    ax2.set_xlim([xmin - hl, xmax + hl])
    ax1.set_ylim([ymin - hl, ymax + hl])
    ax2.set_ylim([ymin - hl, ymax + hl])
    ax1.set_title('Normalized flows')
    ax2.set_title('Normalized pipe diameters')
    plt.show()


def draw_network(sol, coordinates):
    # Draws a general flow network (GI)
    N = len(coordinates)
    topology_list = []
    n_edges = sum(sum(sol('x') > 1e-10))
    prunedsol = {'q':[], 'D':[], '\dot{V}_+':[], '\dot{V}_-': []}
    for i in range(N):
        for j in range(N):
            if sol('x')[i][j] >= 1e-10:
                topology_list.append([i,j])
                prunedsol['q'] = prunedsol['q'] + [mag(sol('q')[i][j])]
                prunedsol['D'] = prunedsol['D'] + [mag(sol('D')[i][j])]
    prunedsol['\dot{V}_+'] = sol('\dot{V}_+')
    prunedsol['\dot{V}_-'] = sol('\dot{V}_-')
    print(prunedsol)
    draw_KT_network(prunedsol, coordinates, topology_list)

def draw_tree(pathnodes):
    nr_vertices = len(pathnodes)
    v_label = list(map(str, [pathnode.id for pathnode in pathnodes]))
    G = jgraph.Graph.Tree(nr_vertices, [len(i.children) for i in pathnodes]) # 2 stands for children number
    lay = G.layout('rt')

    position = {k: lay[k] for k in range(nr_vertices)}
    Y = [lay[k][1] for k in range(nr_vertices)]
    M = max(Y)

    es = EdgeSeq(G) # sequence of edges
    E = [e.tuple for e in G.es] # list of edges

    L = len(position)
    Xn = [position[k][0] for k in range(L)]
    Yn = [2*M-position[k][1] for k in range(L)]
    Xe = []
    Ye = []
    for edge in E:
        Xe+=[position[edge[0]][0],position[edge[1]][0], None]
        Ye+=[2*M-position[edge[0]][1],2*M-position[edge[1]][1], None]

    labels = v_label

    lines = go.Scatter(x=Xe,
                   y=Ye,
                   mode='lines',
                   line=dict(color='rgb(210,210,210)', width=1),
                   hoverinfo='none'
                   )
    dots = go.Scatter(x=Xn,
                  y=Yn,
                  mode='markers',
                  name='',
                  marker=dict(symbol='dot',
                                size=18,
                                color='#6175c1',    #'#DB4551',
                                line=dict(color='rgb(50,50,50)', width=1)
                                ),
                  text=labels,
                  hoverinfo='text',
                  opacity=0.8
                  )
    def make_annotations(pos, text, font_size=10, font_color='rgb(250,250,250)'):
        L=len(pos)
        if len(text)!=L:
            raise ValueError('The lists pos and text must have the same len')
        annotations = go.Annotations()
        for k in range(L):
            annotations.append(
                go.Annotation(
                text=labels[k], # or replace labels with a different list for the text within the circle
                x=pos[k][0], y=2*M-position[k][1],
                xref='x1', yref='y1',
                font=dict(color=font_color, size=font_size),
                showarrow=False)
            )
        return annotations

    axis = dict(showline=False, # hide axis line, grid, ticklabels and  title
            zeroline=False,
            showgrid=False,
            showticklabels=False,
            )

    layout = dict(title= 'Tree with Reingold-Tilford Layout',
              annotations=make_annotations(position, v_label),
              font=dict(size=12),
              showlegend=False,
              xaxis=go.XAxis(axis),
              yaxis=go.YAxis(axis),
              margin=dict(l=40, r=40, b=85, t=100),
              hovermode='closest',
              plot_bgcolor='rgb(248,248,248)'
              )
    data=go.Data([lines, dots])
    fig=dict(data=data, layout=layout)
    fig['layout'].update(annotations=make_annotations(position, v_label))
    py.iplot(fig, filename='Tree-Reingold-Tilf')


if __name__ == '__main__':
    pass
