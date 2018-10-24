import numpy as np
import operator
from collections import OrderedDict
from scipy.spatial import ConvexHull

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


def find_apsp(L_all, coordinate_dict):
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
                apsp_dict[i,j] = [i,j]
                if i == j:
                    d_dict[i,j] = 0
                else:
                    d_dict[i,j] = calc_total_dist(L_all, [i,j])
    for k in range(n):
        for i in range(n):
            for j in range(n):
                if i == j:
                    continue
                elif j >= i:
                    # Algorithm works here
                    newpath = list(set(apsp_dict[i,k] + apsp_dict[k,j]))
                    newdist = calc_total_dist(L_all,newpath)
                    if d_dict[i,j] >= newdist:
                        d_dict[i,j] = newdist
                        apsp_dict[i,j] = newpath
                else:
                    d_dict[i,j] = d_dict[j,i]
                    apsp_dict[i,j] = apsp_dict[j,i].reverse()
    return apsp_dict, d_dict

if __name__ == '__main__':
