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
        # Initialize next node
        current = treepath[i]
        visitedlist.append(current.id)
        # Remove backward looking children
        for j in current.children:
            if previous.id == j.id:
                current.children.remove(j)
        # If starting new edge, find its parent
        if edge == []:
            for j in current.parents:
                if j.id in visitedlist:
                    edge = [j]
                    break
        # Add to forks
        if len(current.children) > 1:
            forks.append(current)
        for j in current.children:
            if j.id in visitedlist:
                current.children.remove(j) # Remove visited children
        edge.append(current)
        previous = current
        if len(current.children) == 1:
            continue # continue on a directly connected edge
        elif len(current.children) == 0: # if in a terminal node
            for fork in forks:
                if previous.id in [u.id for u in fork.children]:
                    edge.append(fork) # check if we have looped back to a fork,
                    break             # or if it is a leaf node
            edges.append(edge)
            edge = []
        else: # if there are many children
            edges.append(edge) # register edge
            edge = [previous] # create new edge
    return edges, forks

if __name__ == '__main__':
    # Stress test case

    #     0 --- 1 --- 2 --- 3
    #     |     |     |     |
    #     4---- 5 --- 6 --- 7
    #     |           |
    #     8 --- 9 --- 10
    #     |     |     |
    #    11---- 12    |
    #     |           |
    #    13---- 14 ---15

    N = 16
    topology_list = [[0,1], [1,2], [2,3],
                     [0,4], [1,5], [2,6], [3,7],
                     [4,5], [5,6], [6,7],
                     [4,8], [6,10],
                     [8,9], [9,10],
                     [8,11], [9,12],
                     [11,12],
                     [11,13], [10,15],
                     [13,14], [14,15]]
    coordinates = {0: [0,0], 1: [1,0], 2: [2,0], 3: [3,0],
                   4: [0,-1], 5: [1,-1], 6: [2,-1], 7: [3,-1],
                   8: [0,-2], 9: [1, -2], 10: [2, -2],
                   11: [0,-3], 12: [1, -3],
                   13: [0, -4], 14: [1, -4], 15: [2,-4]
                  }

    L_all = define_length(coordinates)
    topology_list = return_undirected(topology_list)
    apsp_dict, d_dict = find_apsp(L_all, topology_list, coordinates)

    nodes = nodes_from_topology_list(topology_list, coordinates)

    pathnodes = dfs_tree(nodes)

    edges, forks = find_single_path_edges(pathnodes)
