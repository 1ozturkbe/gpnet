import numpy as np


def define_topology(topology_list, N):
    topology = [[0 for _ in xrange(N)] for _ in xrange(N)]
    for connection in topology_list:
        topology[connection[0]][connection[1]] = 1
    return topology

def return_undirected(topology_list):
    for i in topology_list:
        if [i[1], i[0]] not in topology_list:
            topology_list.append([i[1], i[0]])
    return sorted(topology_list)


def define_length(coordinates):
    L = {}
    for p1 in coordinates:
        for p2 in coordinates:
            L[p1,p2] = np.sqrt((coordinates[p1][0] - coordinates[p2][0]) ** 2 +
                                        (coordinates[p1][1] - coordinates[p2][1]) ** 2)
    return L

