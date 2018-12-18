import numpy as np
import operator
from collections import OrderedDict
from scipy.spatial import ConvexHull

from generate_data import return_undirected, define_length

class Node(object):
    "Generic tree node."
    def __init__(self, id, parents=None, children=None):
        self.id = id
        self.parents = []
        self.children = []
        if parents is not None:
            for parent in parents:
                self.add_parent(parent)
        if children is not None:
            for child in children:
                self.add_child(child)
    def __repr__(self):
        return str(self.id)
    def add_child(self, node):
        assert isinstance(node, Node)
        self.children.append(node)
    def add_parent(self, node):
        assert isinstance(node, Node)
        self.parents.append(node)

def find_convex_hull(coordinate_dict):
    """
    Find the convex hull of a set of points in Euclidian space,
    so we can generate trees from the leaves!
    :return:
    """
    hull_pts = {}
    point_list = [value for key, value in coordinate_dict.iteritems()]
    cvx_hull = ConvexHull(point_list)
    for i in cvx_hull.vertices:
        hull_pts[i] = tuple(cvx_hull.points[i])
    return hull_pts

def calc_total_dist(L_all, coord_inds):
    dists = []
    for i in range(len(coord_inds)-1):
        dists.append(L_all[coord_inds[i],coord_inds[i+1]])
    return sum(dists)


def find_apsp(L_all, topology_list, coordinate_dict):
    """
    WORK IN PROGRESS
    Solves the all-pair-shortest-path problem using Floyd-Warshall algorithm
    :param coordinate_dict:
    :return: apsp_list is a dict of nodes that connect the two indexed nodes with minimal distance
             d_dict is a dict of min distances between the two nodes
    """
    apsp_dict = {}
    d_dict = {}
    n = len(coordinate_dict)
    for i in range(n):
        for j in range(n):
                if i == j:
                    d_dict[i,j] = 0
                    apsp_dict[i,j] = [i,j]
                elif [i,j] in topology_list:
                    d_dict[i,j] = L_all[i,j]
                    apsp_dict[i,j] = [i,j]
                else:
                    d_dict[i,j] = sum(v for i,v in L_all)
                    apsp_dict[i,j] = None
    for k in range(n):
        for i in range(n):
            for j in range(n):
                if i == j:
                    continue
                elif apsp_dict[i,k] and apsp_dict[k,j]:
                    # Algorithm works here
                    newpath = apsp_dict[i,k][0:-1] + apsp_dict[k,j]
                    newdist = calc_total_dist(L_all,newpath)
                    if d_dict[i,j] > newdist:
                        d_dict[i,j] = newdist
                        apsp_dict[i,j] = newpath
    return apsp_dict, d_dict

def nodes_from_topology_list(topology_list, coordinates):
    nodes = [Node(i) for i in coordinates.keys()]
    for i in coordinates.keys():
        for j in topology_list:
            if j[0] == i:
                nodes[i].add_parent(nodes[j[1]])
            if j[1] == i:
                nodes[i].add_child(nodes[j[0]])
    return nodes

def dfs_tree(nodes):
    unexplored = [node.id for node in nodes]
    pathnodes = [nodes[0]]
    stack = pathnodes
    while unexplored != []:
        i = stack.pop()
        while i.id not in unexplored:
            i = stack.pop()
        unexplored.remove(i.id)
        pathnodes = pathnodes + [i]
        for j in sorted(i.children):
            if j.id in unexplored:
                stack = stack + [Node(j.id, [i], j.children)]
    return pathnodes

def topology_list_from_nodes(nodes):
    topology_list = []
    for i in nodes:
        for j in i.parents:
            topology_list.append([j.id, i.id])
        for j in i.children:
            topology_list.append([i.id, j.id])
    return topology_list

def find_single_path_edges(treepath):
    edges = []
    edge = []
    previous = treepath[0]
    visitedlist = []
    forks = []
    for i in range(len(treepath)):
        current = treepath[i]
        visitedlist.append(current.id)
        for j in current.children:
            if previous.id == j.id:
                current.children.remove(j)
            if j.id in visitedlist and edge == []:
                edge = [j]
                current.children.remove(j)
        if len(current.children) > 1:
            forks.append(current)
        for j in current.children:
            if j.id in visitedlist:
                current.children.remove(j)
        edge.append(current)
        previous = current
        if len(current.children) == 1:
            print str(current) +' has single child'
        elif len(current.children) == 0:
            print str(current) +' has no children'
            edges.append(edge)
            edge = []
        else:
            # forks.append(current)
            edges.append(edge)
            edge = [previous]
    return edges, forks

if __name__ == '__main__':
      # Somewhat large problem
    N = 32
    topology_list = [[0, 1], [1, 2], [2, 3], [2, 19], [2, 18], [3, 4], [4, 5], [5, 6], [6, 7], [7, 8], [8, 9], [9, 10],
                     [10, 11], [11, 12], [9, 13], [13, 14], [14, 15], [16, 15], [17, 16], [18, 17], [15, 26], [26, 25],
                     [25, 24], [23, 24], [22, 23], [22, 27], [27, 28], [28, 29], [29, 30], [24, 31], [31, 30], [19, 20],
                     [20, 21], [19, 22]]
    coordinates = {1: (5021.20, 1582.17), 2: (5025.20, 2585.42), 3: (5874.22, 2588.30), 4: (6873.11, 2588.30),
                   5: (8103.51, 2585.42), 6: (8103.51, 3234.67), 7: (8106.66, 4179.28), 8: (8106.66, 5133.78),
                   9: (7318.64, 5133.78), 10: (7319.94, 5831.65), 11: (7319.94, 6671.19), 12: (5636.76, 6676.24),
                   13: (6530.63, 5133.78), 14: (5676.02, 5133.78), 15: (5021.20, 5133.78), 16: (5021.20, 4412.36),
                   17: (5021.20, 3868.52), 18: (5021.20, 3191.49), 19: (3587.87, 2588.30), 20: (3587.87, 1300.84),
                   21: (3587.87, 901.29), 22: (1978.55, 2588.30), 23: (1975.58, 4084.35), 24: (1980.46, 5137.63),
                   25: (3077.46, 5137.63), 26: (3933.52, 5133.78), 27: (846.04, 2588.20), 28: (-552.41, 2588.20),
                   29: (-552.38, 4369.06), 30: (-549.36, 5137.63), 31: (536.45, 5137.63), 0: (5360.71, 1354.05)}

    L_all = define_length(coordinates)
    topology_list = return_undirected(topology_list)
    apsp_dict, d_dict = find_apsp(L_all, topology_list, coordinates)

    nodes = nodes_from_topology_list(topology_list, coordinates)

    pathnodes = dfs_tree(nodes)

    edges, forks = find_single_path_edges(pathnodes)
