import numpy as np


def define_topology(topology_list, N):
    topology = [[0 for _ in xrange(N)] for _ in xrange(N)]
    for connection in topology_list:
        topology[connection[0]][connection[1]] = 1
    return topology


def define_length(coordinates):
    L = [[0 for _ in xrange(len(coordinates))] for _ in xrange(len(coordinates))]
    for p1 in coordinates:
        for p2 in coordinates:
            L[p1 - 1][p2 - 1] = np.sqrt((coordinates[p1][0] - coordinates[p2][0]) ** 2 +
                                        (coordinates[p1][1] - coordinates[p2][1]) ** 2)
    return L

