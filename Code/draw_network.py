import matplotlib.pyplot as plt
import matplotlib.path as path
import numpy as np

def forceAspect(ax,aspect=1):
    im = ax.get_images()
    extent =  im[0].get_extent()
    ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)

def draw_KT_network(sol, coordinates,topology_list):
    coordinate_list = [[coordinates[i[0]], coordinates[i[1]]] for i in topology_list]
    orig, dest = zip(*coordinate_list)
    xmin,ymin, xmax, ymax = 1e10*np.array([1,1,-1,-1])
    # Arrow parameters for flow plotting
    hl = 400 # Arrow head length
    hwf = hl*sol('F').magnitude/max(sol('F').magnitude)
    hwd = hl*sol('D').magnitude/max(sol('D').magnitude)
    fig = plt.figure()
    ax1 = fig.add_subplot(1,2,1,aspect=1)
    ax2 = fig.add_subplot(1,2,2,aspect=1)
    for i in range(len(topology_list)):
        xmin = np.min([xmin,orig[i][0]])
        xmax = np.max([xmax,orig[i][0]])
        ymin = np.min([ymin,orig[i][1]])
        ymax = np.max([ymax,orig[i][1]])
        ax1.arrow(orig[i][0],orig[i][1], dest[i][0]-orig[i][0],dest[i][1]-orig[i][1],label = i, color='b',
                  width = 1./3.*hwf[i], head_width = hwf[i], head_length = hl,
                  length_includes_head = True)
        ax2.arrow(orig[i][0],orig[i][1], dest[i][0]-orig[i][0],dest[i][1]-orig[i][1],label = i, color='r',
                  width = 1./3.*hwd[i], head_width = hwd[i], head_length = hl,
                  length_includes_head = True)
    # Formatting
    ax1.set_xlim([xmin-hl, xmax+hl])
    ax2.set_xlim([xmin-hl, xmax+hl])
    ax1.set_ylim([ymin-hl, ymax+hl])
    ax2.set_ylim([ymin-hl, ymax+hl])
    ax1.set_title('Normalized flows')
    ax2.set_title('Normalized pipe diameters')
    plt.show()



if __name__ == '__main__':
    draw_KT_network(sol,coordinates,topology_list)

