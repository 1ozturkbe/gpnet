import matplotlib.pyplot as plt
import matplotlib.path as path
import numpy as np

from gpkit.small_scripts import mag


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
        hwf = hl * mag(sol('F')) / max(mag(sol('F')))  # Arrow width for flow solution
        hwd = hl * mag(sol('D')) / max(mag(sol('D')))  # Arrow width for diameter solution
        pwSr = 15*mag(sol('Sr')) / max(mag(sol('Sr')))  # Point width for sources
        pwSk = 15*mag(sol('Sk')) / max(mag(sol('Sk')))  # Point width for sinks
    except:
        hwf = hl * mag(sol['F']) / max(mag(sol['F']))  # Arrow width for flow solution
        hwd = hl * mag(sol['D']) / max(mag(sol['D']))  # Arrow width for diameter solution
        pwSr = 15*mag(sol['Sr']) / max(mag(sol['Sr']))  # Point width for sources
        pwSk = 15*mag(sol['Sk']) / max(mag(sol['Sk']))  # Point width for sinks
    fig = plt.figure()
    ax1 = fig.add_subplot(1, 2, 1, aspect=1)
    ax2 = fig.add_subplot(1, 2, 2, aspect=1)
    for i in range(len(topology_list)):
        xmin = np.min([xmin, orig[i][0]])
        xmax = np.max([xmax, orig[i][0]])
        ymin = np.min([ymin, orig[i][1]])
        ymax = np.max([ymax, orig[i][1]])
        # Plotting flow arrows
        ax1.arrow(orig[i][0], orig[i][1], dest[i][0] - orig[i][0], dest[i][1] - orig[i][1], label=i, color='b',
                  width=1. / 3. * hwf[i], head_width=hwf[i], head_length=hl,
                  length_includes_head=True)
        # Plotting diameter arrows
        ax2.arrow(orig[i][0], orig[i][1], dest[i][0] - orig[i][0], dest[i][1] - orig[i][1], label=i, color='r',
                  width=1. / 3. * hwd[i], head_width=hwd[i], head_length=hl,
                  length_includes_head=True)

    for i in range(len(coordinates)):
        # Plotting sources and sinks
        ax1.plot(coordinates[i][0], coordinates[i][1], 'o', color='b', mfc='none', markersize=pwSr[i])
        ax1.plot(coordinates[i][0], coordinates[i][1], 'o', color='r', mfc='none', markersize=pwSk[i])
    ax1.set_xlim([xmin - hl, xmax + hl])
    ax2.set_xlim([xmin - hl, xmax + hl])
    ax1.set_ylim([ymin - hl, ymax + hl])
    ax2.set_ylim([ymin - hl, ymax + hl])
    ax1.set_title('Normalized flows')
    ax2.set_title('Normalized pipe diameters')
    plt.show()


def draw_network(sol, coordinates):
    # Draws a general flow network
    N = len(coordinates)
    topology_list = []
    n_edges = sum(sol('x') > 1e-20)
    prunedsol = {'F':[], 'D':[], 'Sr':[], 'Sk': []}
    for i in range(N):
        for j in range(N):
            print sol('x')[i,j]
            if sol('x')[i,j] >= 0:
                topology_list.append([i,j])
                prunedsol['F'] += mag(sol('F')[i,j])
                prunedsol['D'] += mag(sol('D')[i,j])
    prunedsol['Sr'] = sol('Sr')
    prunedsol['Sk'] = sol('Sk')
    print prunedsol
    draw_KT_network(prunedsol, coordinates, topology_list)

if __name__ == '__main__':
    pass
