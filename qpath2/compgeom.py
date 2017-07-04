#
# QPATH2 - a quantitative pathology toolkit
#
# (c) 2017 Vlad Popovici
#

# QPATH2.COMPGEOM: a few computational geometry functions, used mainly for
# annotation processing (e.g. find the intersection between an polygon and
# a ROI in the image).

__all__ = ['simple_polygon_intersection',
           'point_wrt_polygon', 'polygon_equality', 'polygon_is_convex',
           'polygon_is_collinear', 'polygon_is_counterclockwise',
           'rect_inside_polygon', 'polygon_inside_polygon']


from qpath2.compgeom_ import simple_polygon_intersection_, \
    point_wrt_polygon_, \
    polygon_equality_

from CGAL.CGAL_Kernel import Polygon_2, Point_2

import numpy as np
import qpath2.core as core


##-
def polygon_equality(P, Q):
    """Test whether P == Q.

    Args:
        P (numpy.array): (n x 2) The vertex coordinates for the first polygon
            ((x, y) by rows).
        Q (numpy.array): (m x 2) The vertex coordinates for the second polygon
            ((x, y) by rows).

    Returns:
        bool
    """

    n = polygon_equality_(P[:,0].tolist(), P[:,1].tolist(),
                          Q[:,0].tolist(), Q[:,1].tolist())

    if n >= 0:
        return n != 0

    if n == -1 or n == -2:
        raise core.Error("Size mismatch in P or Q")
    else:
        raise core.Error("Unknown error")
##-


##-
def polygon_is_convex(P):
    """Test whether a polygon is convex.
    
    Args:
        P (numpy.array): (n x 2) The vertex coordinates for the polygon
            ((x, y) by rows).
            
    Returns:
        bool
    """
    # use the CGAL official wrapper
    Q = Polygon_2()
    for v in P:
        Q.push_back(Point_2(v[0], v[1]))

    return Q.is_convex()
##-


##-
def polygon_is_collinear(P):
    """Test whether a polygon is degenerated.

    Args:
        P (numpy.array): (n x 2) The vertex coordinates for the polygon
            ((x, y) by rows).

    Returns:
        bool
    """
    # use the CGAL official wrapper
    Q = Polygon_2()
    for v in P:
        Q.push_back(Point_2(v[0], v[1]))

    return Q.is_collinear_oriented()
##-


##-
def polygon_is_counterclockwise(P):
    """Test whether the vertices of a polygon are ordered counterclockwise.

    Args:
        P (numpy.array): (n x 2) The vertex coordinates for the polygon
            ((x, y) by rows).

    Returns:
        bool
    """
    # use the CGAL official wrapper
    Q = Polygon_2()
    for v in P:
        Q.push_back(Point_2(v[0], v[1]))

    return Q.is_counterclockwise_oriented()
##-


##-
def simple_polygon_intersection(P, Q):
    """Compute the intersection of two simple polygons (i.e. without holes and
    a single component).

    Args:
        P (numpy.array): (n x 2) The vertex coordinates for the first polygon
            ((x, y) by rows).
        Q (numpy.array): (m x 2) The vertex coordinates for the second polygon
            ((x, y) by rows).

    Returns:
        a list of numpy.arrays: The intersection may result in several simple
        polygons. Each of them is returned as a two-column numpy.array of
        (x,y) coordinates.
    """
    rx, ry, rn = [], [], []

    n = simple_polygon_intersection_(P[:,0].tolist(), P[:,1].tolist(),
                                     Q[:,0].tolist(), Q[:,1].tolist(),
                                     rx, ry, rn)

    if n == 0:
        return np.array([])
    elif n == -1 or n == -2:
        raise core.Error("Size mismatch in P or Q")
    elif n == -3:
        raise core.Error("The polygons are not simple")
    elif n == -4:
        raise core.Error("The intersection is not bounded")

    if n > 1:
        # several polygons in the intersection
        e = np.cumsum(rn).tolist()
        b = [0] + e[:-1]
        res = [np.array(zip(rx[i:j], ry[i:j])) for i, j in zip(b, e)]
        return res

    return [np.array(zip(rx, ry))]
##-


##-
def point_wrt_polygon(points, Q):
    """Test the position of a list of points with respect to a polygon.

    Args:
        points (list): a list of (x,y) coordinates
        Q (numpy.array): (m x 2) the vertices of a polygon

    Returns:
        a list of codes, one code per point in <points> list:
             1 for a point INSIDE the polygon
             0 for a point on the BOUNDARY of the polygon
            -1 for a point OUTSIDE the polygon
    """

    r = []

    n = point_wrt_polygon_([_p[0] for _p in points],
                           [_p[1] for _p in points],
                           Q[:,0].tolist(), Q[:,1].tolist(), r)

    if n == -1 or n == -2:
        raise core.Error("Size mismatch in p or Q")
    elif n < -2:
        raise core.Error("Unspecific error")

    return r
##-


##-
def rect_inside_polygon(r0, r1, P):
    """Test whether a rectangle is completely inside a polygon.

    Args:
        r0, r1 (pairs): top-left and bottom-right corners, given as
            (x,y) coordinates
        P (numpy.array): (m x 2) the vertices of a polygon

    Returns:
        bool
    """

    if polygon_is_convex(P):
        # a simpler case, it's enough if r0 and r1 are inside P
        # to decide that the whole rectangle lies within the polygon
        r = point_wrt_polygon([r0, r1], P)
        return np.all(np.array(r) == 1)

    # the polygon is not convex: compute the intersection with the
    # rectangle. If this intersection is empty, the rectangle is not
    # inside the polygon. Otherwise, if the intersection equals the
    # rectangle, then the rectangle is inside the polygon.
    r = simple_polygon_intersection(P, np.array([r0, (r0[0],r1[1]), r1, (r1[0],r0[1])]))
    if len(r) == 0 or len(r) > 1:
        # if the intersection is either empty or contains too many components
        # (i.e. the polygons cross) then:
        return False

    return polygon_equality(P, r[0])
##-


##-
def polygon_inside_polygon(P, Q):
    """Test whether a polygon (P) is completely inside a second polygon (Q).

    Args:
        P (numpy.array): (n x 2) The vertex coordinates for the first polygon
            ((x, y) by rows).
        Q (numpy.array): (m x 2) The vertex coordinates for the second polygon
            ((x, y) by rows).

    Returns:
        bool
    """

    # P is contained in Q iif P \cap Q == P
    r = simple_polygon_intersection(P, Q)
    if len(r) == 0 or len(r) > 1:
        # if the intersection is either empty or contains too many components
        # (i.e. the polygons cross) then:
        return False

    return polygon_equality(P, r[0])
##-
