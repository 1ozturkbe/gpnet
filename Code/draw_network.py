import matplotlib.pyplot as plt
import matplotlib.path as path
import numpy as np

def draw_KT_network(sol, coordinates,topology_list):
    coordinate_list = [[coordinates[i[0]], coordinates[i[1]]] for i in topology_list]
    orig, dest = zip(*coordinate_list)
    xmin,ymin, xmax, ymax = 1e10*np.array([1,1,-1,-1])
    ax = plt.axes()
    hl = 200 # Arrow head length
    hw = 0.5*hl*sol('F').magnitude
    for i in range(len(topology_list)):
        print i
        xmin = np.min([xmin,orig[i][0]])
        xmax = np.max([xmax,orig[i][0]])
        ymin = np.min([ymin,orig[i][1]])
        ymax = np.max([ymax,orig[i][1]])
        plt.arrow(orig[i][0],orig[i][1], dest[i][0]-orig[i][0],dest[i][1]-orig[i][1],label = i,
                  width = 1./3.*hw[i], head_width = hw[i], head_length = hl,
                  length_includes_head = True)
    ax.set_xlim([xmin-hl, xmax+hl])
    ax.set_ylim([ymin-hl, ymax+hl])
    plt.show()



if __name__ == '__main__':
    draw_KT_network(sol,coordinates,topology_list)

