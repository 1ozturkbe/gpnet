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
