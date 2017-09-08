#
# QPATH2 - a quantitative pathology toolkit
#
# (c) 2017 Vlad Popovici
#

#
# QPATH2.ANNOT.CORE: core classes and functions for slide annotation handling.
#

__all__ = ['Dot', 'Polygon', 'PointSet']

import abc
import numpy as np
import collections

from qpath2.core import Error, TinyVector

##-
class AnnotationObject(object):
    __metaclass__ = abc.ABCMeta
    """Defne the AnnotationObject minimal interface. This class is made
    abstract to force more meaningful names (e.g. Dot, Polygon, etc.) in
    subclasses."""

    xy = np.ndarray((0,0), dtype=np.float64)  # empty array
    _name = None
    _annotation_type = None

    @abc.abstractmethod
    def duplicate(self):
        pass

    def __str__(self):
        """Return a string representation of the object."""
        return str(self.type) + " <" + str(self.name) + ">: \n" + str(self.xy)

    def bounding_box(self):
        """Compute the bounding box of the object."""
        return TinyVector(self.xy.min(axis=0)), TinyVector(self.xy.max(axis=0))

    def translate(self, x_off, y_off=None):
        """Translate the object by a vector [x_off, y_off], i.e.
        the new coordinates will be x' = x + x_off, y' = y + y_off.
        If y_off is None, then the same value as in x_off will be
        used."""
        if y_off is None:
            y_off = x_off
        self.xy += np.array([x_off, y_off])

    def scale(self, x_scale, y_scale=None):
        """Scale the object by a specified factor. The new coordinates
        will be x' = x * x_scale, y' = y * y_scale. If y_scale is
        None, x_scale will be used instead (i.e. isotrop scaling)."""
        if y_scale is None:
            y_scale = x_scale
        self.xy *= np.array([x_scale, y_scale])

    def affine(self, M):
        """Apply an affine transformation to all points of the annotation.

        If M is the affine transformation matrix, the new coordinates
        (x', y') of a point (x, y) will be computed as

        x' = M[1,1] x + M[1,2] y + M[1,3]
        y' = M[2,1] x + M[2,2] y + M[2,3]

        In other words, if P is the 3 x n matrix of n points,
        P = [x; y; 1]
        then the new matrix Q is given by
        Q = M * P

        Parmeters:
        ---------
            M: numpy array [2 x 3]

        Returns:
            nothing
        """

        # simpler than using the matrix operations
        qx = self.xy[:, 0] * M[0, 0] + self.xy[:, 1] * M[0, 1] + M[0, 2]
        qy = self.xy[:, 0] * M[1, 0] + self.xy[:, 1] * M[1, 1] + M[1, 2]

        self.xy[:, 0] = qx
        self.xy[:, 1] = qy

        return

    @property
    def x(self):
        """Return the x coordinate(s) of the object. This is always a
        vector, even for a single point (when it has one element)."""
        return self.xy[:,0]

    @property
    def y(self):
        """Return the y coordinate(s) of the object. This is always a
        vector, even for a single point (when it has one element)."""
        return self.xy[:,1]

    @property
    def size(self):
        """Return the number of points defining the object."""
        return self.xy.shape[0]

    @property
    def name(self):
        """Return the name of the annotation object."""
        return self._name

    @property
    def type(self):
        """Return the annotation type as a string."""
        return self._annotation_type
##-


##-
class Dot(AnnotationObject):
    """Dot: a single position in the image."""

    def __init__(self, x, y=None, name=None):
        self._annotation_type = "DOT"
        self._name = "DOT"

        if name is not None:
            self._name = name

        if y is not None:
            if isinstance(y, str):
                self._name = y
            else:
                self.xy = np.array([[x, y]], dtype=np.float64)
                return

        # check whether x is iterable and build the coords from it:
        if isinstance(x, collections.Iterable):
            self.xy = np.array(x, dtype=np.float64)[:2]
            return
        raise Error('x parameter cannot be interpretated as a 2D vector')

    def duplicate(self):
        return Dot(self.xy, name=self.name)
##-


##-
class PointSet(AnnotationObject):
    """PointSet: an ordered collection of points."""

    def __init__(self, x, y=None, name=None):
        self._annotation_type = "POINTSET"
        self._name = "POINTS"

        if name is not None:
            self._name = name

        if y is not None:
            if isinstance(y, str):
                self._name = y
            else:
                self.xy = np.array(zip(x, y), dtype=np.float64)
                return

        # check whether x is iterable and build the coords from it:
        if isinstance(x, collections.Iterable):
            self.xy = np.array(x, dtype=np.float64)[:,:2]
            return
        raise Error('x parameter cannot be interpretated as a 2D array')


    def duplicate(self):
        return PointSet(self.xy, name=self.name)
##-


##-
class Polygon(PointSet):
    """PointSet: an ordered collection of points."""

    def __init__(self, x, y=None, name=None):
        if name is None:
            name = "POLYGON"
        super(Polygon, self).__init__(x, y=y, name=name)
        self._annotation_type = "POLYGON"

        # ensure a closed contour:
        if not np.all(self.xy[0,] == self.xy[-1,]):
            self.xy = np.concatenate((self.xy, [self.xy[0,]]), axis=0)


    def duplicate(self):
        return Polygon(self.xy, name=self.name)
##-
