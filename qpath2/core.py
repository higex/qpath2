#
# QPATH2 - a quantitative pathology toolkit
#
# (c) 2017 Vlad Popovici
#

"""QPATH2.CORE: Core classes and functions.

Defines exception classes and other basic classes.
"""

__all__ = ['Error', 'WSIInfo', 'MRIBase', 'MRI',
           'MRIExplorer', 'MRISlidingWindow',
           'expandSlicing']


import openslide as osl
import abc
import numpy as np
from math import ceil, floor, pow, log10
import collections
from functools import partial
from itertools import combinations

class Error(Exception):
    """Basic error exception for QPATH2.

    Args:
        msg (str): Human-readable string describing the exception.
        code (:obj:`int`, optional): Error code.

    Attributes:
        msg (str): Human-readable string describing the exception.
        code (int): Error code.
    """

    def __init__(self, msg, code=1, *args):
        self.message = "QPATH2: " + msg
        self.code = code
        super(Error, self).__init__(msg, code, *args)


##-
class WSIInfo(object):
    """Hold some basic info about a WSI.

    Args:
        path (str): full path to WSI file

    Attributes:
        path (str): full path to WSI file
        info (dict): a dictionary containing WSI properties
    """

    path = None
    info = {}

    def __init__(self, path):
        self.path = path
        with osl.OpenSlide(path) as wsi:
            self.info = {'vendor': wsi.properties['openslide.vendor'],
                         'x_mpp': float(wsi.properties['openslide.mpp-x']),
                         'y_mpp': float(wsi.properties['openslide.mpp-y']),
                         'objective': int(wsi.properties['openslide.objective-power']),
                         'x_offset': 0,
                         'y_offset': 0}
            # fill in level data:
            self.info['level_count'] = wsi.level_count
            lv = dict()
            for k in range(wsi.level_count):
                lv[k] = {'x_size': long(wsi.level_dimensions[k][0]),
                         'y_size': long(wsi.level_dimensions[k][1]),
                         'downsample_factor': float(wsi.level_downsamples[k])}
                if 'openslide.level['+str(k)+'].tile-width' in wsi.properties:
                    lv[k]['tile_x_size'] = int(wsi.properties['openslide.level['+str(k)+'].tile-width']),
                    lv[k]['tile_y_size'] = int(wsi.properties['openslide.level[' + str(k) + '].tile-height'])
            self.info['levels'] = lv

            if wsi.properties['openslide.vendor'] == 'hamamatsu':
                self.info['x_offset'] = long(wsi.properties['hamamatsu.XOffsetFromSlideCentre'])
                self.info['y_offset'] = long(wsi.properties['hamamatsu.YOffsetFromSlideCentre'])
##-


##-
class MRIBase(object):
    """Base class for MultiResolutionImages presenting a uniform interface for retrieving
    pixels and regions. Note that changes in image data are not propagated back to the original
    image. This is an abstract class.

    Args:
        wsi_info (WSIInfo): an info object for a whole slide image

    Attributes:
        info (WSIInfo)
    """
    __metaclass__ = abc.ABCMeta

    def __init__(self, wsi_info):
        if not isinstance(wsi_info, WSIInfo):
            raise Error("Only WSIInfo instances are accepted")

        self._wsi_info = wsi_info

    @property
    def info(self):
        return self._wsi_info.info

    @property
    def path(self):
        return self._wsi_info.path

    @property
    def widths(self):
        return [self.info['levels'][l]['x_size'] for l in range(self.nlevels)]

    @property
    def heights(self):
        return [self.info['levels'][l]['y_size'] for l in range(self.nlevels)]

    @property
    def nlevels(self):
        return self.info['level_count']

    @abc.abstractmethod
    def get_region_px(self, x0, y0, width, height, level, as_type=np.uint8):
        """Read a region from the image source. The region is specified in
        pixel coordinates.

        Args:
            x0, y0 (long): top left corner of the region (in pixels, at the specified
            level)
            width, height (long): width and height (in pixels) of the region
            level (int): the magnification level to read from
            as_type: type of the pixels (default numpy.uint8)

        Returns:
            a numpy.ndarray (OpenCV channel ordering: (A)BGR)
        """
        pass

    @abc.abstractmethod
    def get_region(self, x0, y0, width, height, level, as_type=np.uint8):
        """Read a region from the image source. The region is specified in
        slide coordinates.

        Args:
            x0, y0 (long): top left corner of the region (in slide units)
            width, height (long): width and height (in slide units) of the region
            level (int): the magnification level to read from
            as_type: type of the pixels (default numpy.uint8)

        Returns:
            a numpy.ndarray (OpenCV channel ordering: (A)BGR)
        """
        pass

##-


##-
class MRI(MRIBase):
    """A multi-resolution image backed by OpenSlide.

    Args:
        wsi_info (WSIInfo): info about the slide

    Attributes:
        see MRIBase
    """
    _reader = None

    def __init__(self, wsi_info):
        if not isinstance(wsi_info, WSIInfo):
            raise Error("Only WSIInfo instances are accepted")

        self._wsi_info = wsi_info
        self._reader = osl.OpenSlide(self.path)

    def get_region_px(self, x0, y0, width, height, level, as_type=np.uint8):
        """Read a region from the image source. The region is specified in
            pixel coordinates.

            Args:
                x0, y0 (long): top left corner of the region (in pixels, at the specified
                level)
                width, height (long): width and height (in pixels) of the region
                level (int): the magnification level to read from
                as_type: type of the pixels (default numpy.uint8)

            Returns:
                a numpy.ndarray
        """

        # OpenSlide requires specification of (x0,y0) in level-0 coordinates, so
        # convert (x0,y0) if level != 0
        if level > 0:
            # get scaling factors between current level and level 0:
            sx = self.info['levels'][0]['x_size'] / self.info['levels'][level]['x_size']
            sy = self.info['levels'][0]['y_size'] / self.info['levels'][level]['y_size']
            x0 *= sx
            y0 *= sy

        x0, y0, width, height = [long(_x) for _x in [x0, y0, width, height]]
        pil_obj = self._reader.read_region((x0, y0), level, (width, height))
        img_data = np.asarray(pil_obj, dtype=as_type)

        return img_data


    def get_region(self, x0, y0, width, height, level, as_type=np.uint8):
        raise Error("Not yet implemented")
##-


##-
class MRIExplorer(object):
    """Defines an interface for multi-resolution image explorers. An image
    explorer simply returns positions in an image rather than parts of the
    image itself. Hence, it only needs to know about the extent of the image.
    """
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def reset(self):
        """Reset the explore, next call to next() will start from the
        initial conditions.
        """
        pass

    @abc.abstractmethod
    def last(self):
        """Go to last position and return it."""
        pass

    @abc.abstractmethod
    def next(self):
        """Go to next position."""
        pass

    @abc.abstractmethod
    def prev(self):
        """Go to previous position."""
        pass

    @abc.abstractmethod
    def here(self):
        """Returns current position, does not change it."""
        pass

    @abc.abstractmethod
    def total_steps(self):
        """Returns the total number of steps to iterate over all positions
        in the image, according to the specific schedule.
        """
        pass

    def __iter__(self):
        return self

    def __next__(self):
        return self.next()

    def __prev__(self):
        return self.prev()
##-


##-
class MRISlidingWindow(MRIExplorer):
    """A sliding window image explorer. It returns successively the coordinates
    of the sliding window as a tuple (x0, y0, x1, y1).

    Args:
        image_shape : tuple (nrows, ncols)
            Image shape (img.shape).
        w_size : tuple (width, height)
            Window size as a pair of width and height values.
        start : tuple (x0, y0)
            Top left corner of the first window. Defaults to (0,0).
        step : tuple (x_step, y_step)
            Step size for the sliding window, as a pair of horizontal
            and vertical steps. Defaults to (1,1).
    """
    def __init__(self, image_shape, w_size, start=(0,0), step=(1,1)):
        self._image_shape = image_shape
        self._w_size = w_size
        self._start = start
        self._step = step
        self._k = 0

        img_h, img_w = image_shape

        if w_size[0] < 2 or w_size[1] < 2:
            raise ValueError('Window size too small.')

        if img_w < start[0] + w_size[0] or img_h < start[1] + w_size[1]:
            raise ValueError('Start position and/or window size out of image.')

        x, y = np.meshgrid(np.arange(start[0], img_w - w_size[0] + 1, step[0]),
                           np.arange(start[1], img_h - w_size[1] + 1, step[1]))

        self._top_left_corners = zip(x.reshape((-1,)).tolist(),
                                     y.reshape((-1,)).tolist())

    def total_steps(self):
        return len(self._top_left_corners)

    def reset(self):
        self._k = 0

    def here(self):
        if 0 <= self._k < self.total_steps():
            x0, y0 = self._top_left_corners[self._k]
            x1 = min(x0 + self._w_size[0], self._image_shape[1])
            y1 = min(y0 + self._w_size[1], self._image_shape[0])

            return x0, y0, x1, y1
        raise Error("Position outside bounds")

    def last(self):
        if self.total_steps() > 0:
            self._k = self.total_steps() - 1
            x0, y0, x1, y1 = self.here()
            return x0, y0, x1, y1
        else:
            raise Error("Empty iterator")

    def next(self):
        if self._k < self.total_steps():
            x0, y0, x1, y1 = self.here()
            self._k += 1
            return x0, y0, x1, y1
        else:
            raise StopIteration()

    def prev(self):
        if self._k >= 1:
            self._k -= 1
            x0, y0, x1, y1 = self.here()
            return x0, y0, x1, y1
        else:
            raise StopIteration()
##-


## Here are some very handy functions for array slicing taken from ILASTIK/lazyflow.
## Original copyright follows.

###############################################################################
#   lazyflow: data flow based lazy parallel computation framework
#
#       Copyright (C) 2011-2014, the ilastik developers
#                                <team@ilastik.org>
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the Lesser GNU General Public License
# as published by the Free Software Foundation; either version 2.1
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# See the files LICENSE.lgpl2 and LICENSE.lgpl3 for full text of the
# GNU Lesser General Public License version 2.1 and 3 respectively.
# This information is also available on the ilastik web site at:
#		   http://ilastik.org/license/
###############################################################################
class TinyVector(list):
    __slots__ = []

    def copy(self):
        return TinyVector(self)

    def __add__(self, other):
        if isinstance(other, collections.Iterable):
            return TinyVector(x + y for x, y in zip(self, other))
        else:
            return TinyVector(x + other for x in self)

    __radd__ = __add__

    def __iadd__(self, other):
        # Must explicitly override list.__iadd__
        # Others (e.g. isub, imul) can use default implementation.
        if isinstance(other, collections.Iterable):
            self = TinyVector(x + y for x, y in zip(self, other))
            return self
        else:
            self = TinyVector(x + other for x in self)
            return self

    def __sub__(self, other):
        if isinstance(other, collections.Iterable):
            return TinyVector(x - y for x, y in zip(self, other))
        else:
            return TinyVector(x - other for x in self)

    def __rsub__(self, other):
        if isinstance(other, collections.Iterable):
            return TinyVector(y - x for x, y in zip(self, other))
        else:
            return TinyVector(other - x for x in self)

    def __mul__(self, other):
        if isinstance(other, collections.Iterable):
            return TinyVector(x * y for x, y in zip(self, other))
        else:
            return TinyVector(x * other for x in self)

    __rmul__ = __mul__

    def __div__(self, other):
        if isinstance(other, collections.Iterable):
            return TinyVector(x / y for x, y in zip(self, other))
        else:
            return TinyVector(x / other for x in self)

    def __rdiv__(self, other):
        if isinstance(other, collections.Iterable):
            return TinyVector(y / x for x, y in zip(self, other))
        else:
            return TinyVector(other / x for x in self)

    def __mod__(self, other):
        if isinstance(other, collections.Iterable):
            return TinyVector(x % y for x, y in zip(self, other))
        else:
            return TinyVector(x % other for x in self)

    def __rmod__(self, other):
        if isinstance(other, collections.Iterable):
            return TinyVector(y % x for x, y in zip(self, other))
        else:
            return TinyVector(other % x for x in self)

    def __floordiv__(self, other):
        if isinstance(other, collections.Iterable):
            return TinyVector(x // y for x, y in zip(self, other))
        else:
            return TinyVector(x // other for x in self)

    def __rfloordiv__(self, other):
        if isinstance(other, collections.Iterable):
            return TinyVector(y // x for x, y in zip(self, other))
        else:
            return TinyVector(other // x for x in self)

    def __eq__(self, other):
        if isinstance(other, collections.Iterable):
            return TinyVector(x == y for x, y in zip(self, other))
        else:
            return TinyVector(x == other for x in self)

    def __ne__(self, other):
        if isinstance(other, collections.Iterable):
            return TinyVector(x != y for x, y in zip(self, other))
        else:
            return TinyVector(x != other for x in self)

    def __ge__(self, other):
        if isinstance(other, collections.Iterable):
            return TinyVector(x >= y for x, y in zip(self, other))
        else:
            return TinyVector(x >= other for x in self)

    def __le__(self, other):
        if isinstance(other, collections.Iterable):
            return TinyVector(x <= y for x, y in zip(self, other))
        else:
            return TinyVector(x <= other for x in self)

    def __gt__(self, other):
        if isinstance(other, collections.Iterable):
            return TinyVector(x > y for x, y in zip(self, other))
        else:
            return TinyVector(x > other for x in self)

    def __lt__(self, other):
        if isinstance(other, collections.Iterable):
            return TinyVector(x < y for x, y in zip(self, other))
        else:
            return TinyVector(x < other for x in self)

    def __and__(self, other):
        if isinstance(other, collections.Iterable):
            return TinyVector(x & y for x, y in zip(self, other))
        else:
            return TinyVector(x & other for x in self)

    __rand__ = __and__

    def __or__(self, other):
        if isinstance(other, collections.Iterable):
            return TinyVector(x | y for x, y in zip(self, other))
        else:
            return TinyVector(x | other for x in self)

    __ror__ = __or__

    def __xor__(self, other):
        if isinstance(other, collections.Iterable):
            return TinyVector(x ^ y for x, y in zip(self, other))
        else:
            return TinyVector(x ^ other for x in self)

    __rxor__ = __xor__

    def __neg__(self):
        return TinyVector(-x for x in self)

    def __abs__(self):
        return TinyVector(map(abs, self))

    def __pos__(self):
        return TinyVector(self)

    def __invert__(self):
        return TinyVector(~x for x in self)

    def ceil(self):
        return TinyVector(map(ceil, self))
        # return numpy.ceil(numpy.array(self))

    def floor(self):
        return TinyVector(map(floor, self))
        # return numpy.floor(numpy.array(self))

    def _asint(self):
        return TinyVector(map(lambda x: x.__int__(), self))

    def insert(self, index, value):
        l = list(self)
        l.insert(index, value)
        return TinyVector(l)

    def all(self):
        answer = True
        for e in self:
            if not e:
                answer = False
                break
        return answer

    def any(self):
        answer = False
        for e in self:
            if e:
                answer = True
                break
        return answer


def expandSlicing(s, shape):
    """
    Args:
        s: Anything that can be used as a numpy array index:
           - int
           - slice
           - Ellipsis (i.e. ...)
           - Some combo of the above as a tuple or list

        shape: The shape of the array that will be accessed

    Returns:
        A tuple of length N where N=len(shape)
        slice(None) is inserted in missing positions so as not to change the meaning of the slicing.
        e.g. if shape=(1,2,3,4,5):
            0 --> (0,:,:,:,:)
            (0:1) --> (0:1,:,:,:,:)
            : --> (:,:,:,:,:)
            ... --> (:,:,:,:,:)
            (0,0,...,4) --> (0,0,:,:,4)
    """
    if type(s) == list:
        s = tuple(s)
    if type(s) != tuple:
        # Convert : to (:,), or 5 to (5,)
        s = (s,)

    # Compute number of axes missing from the slicing
    if len(shape) - len(s) < 0:
        assert s == (Ellipsis,) or s == (slice(None),), \
            "Slicing must not have more elements than the shape, except for [:] and [...] slices.\n" \
            "Your slicing: {}, your shape: {}".format(s, shape)

        # Replace Ellipsis with (:,:,:)
    if Ellipsis in s:
        ei = s.index(Ellipsis)  # Ellipsis Index
        s = s[0:ei] + (len(shape) - len(s) + 1) * (slice(None),) + s[ei + 1:]

    # Append (:,) until we get the right length
    s += (len(shape) - len(s)) * (slice(None),)

    # Special case: we allow [:] and [...] for empty shapes ()
    if shape == ():
        s = ()

    return s


sTrl1 = lambda x: x if type(x) != slice else x.start if x.start != None else 0
sTrl2 = lambda x, y: y if type(y) != slice else y.stop if y.stop != None else x
sTrl3 = lambda x, y: y + 1 if x == y else y


def sliceToRoi(s, shape, extendSingleton=True):
    """Args:
            slice: slice object (1D) or list of slice objects (N-D)
            shape: the shape of the array to be sliced
            extendSingleton: if True, convert int indexes to slices so the dimension of the slicing matches the dimension of the shape.
       Returns:
            ROI instance corresponding to slice
    """
    s = expandSlicing(s, shape)
    start = map(sTrl1, s)
    stop = map(sTrl2, shape, s)
    if extendSingleton:
        stop = map(sTrl3, start, stop)
    return TinyVector(start), TinyVector(stop)


def roiFromShape(shape):
    start = TinyVector([0] * len(shape))
    stop = TinyVector(shape)
    return (start, stop)


def fullSlicing(shape):
    return roiToSlice(*roiFromShape(shape))


def getIntersection(roiA, roiB, assertIntersect=True):
    start = np.maximum(roiA[0], roiB[0])
    stop = np.minimum(roiA[1], roiB[1])

    if ((stop - start) <= 0).any():
        if assertIntersect:
            assert ((stop - start) > 0).all(), "Rois do not intersect: {} and {}".format(roiA, roiB)
        else:
            return None
    return (start, stop)


rTsl1 = lambda x, y: slice(x.__int__(), y.__int__())


def roiToSlice(start, stop, hardBind=False):
    """Args:
            start (N-D coordinate): inclusive start
            stop  (N-D coordinate): exclusive stop
       Returns:
            list of slice objects describing the [start,stop) range,
            so that numpy.ndarray.__getitem__(roiToSlice(start,stop)) can be used.
    """
    if hardBind:
        res = []
        for sta, stp in zip(start, stop):
            if stp == sta + 1 or stp == sta:
                res.append(int(sta))
            else:
                res.append(slice(int(sta), int(stp)))
        return tuple(res)
    else:
        return tuple(map(rTsl1, start, stop))


def nonzero_bounding_box(data):
    """
    For an array with sparsely distributed non-zero values,
      find the bounding box (a ROI) of the non-zero values.

    Example:
        >>> data = numpy.zeros( (10,100,100) )
        >>> data[4, 30:40, 50:60] = 1
        >>> data[7, 45:55, 30:35] = 255
        >>> nonzero_bounding_box(data)
        array([[ 4, 30, 30],
               [ 8, 55, 60]])

    """
    nonzero_coords = np.nonzero(data)
    if len(nonzero_coords[0]) == 0:
        block_bounding_box_roi = np.array(([0] * data.ndim, [0] * data.ndim))
    else:
        block_bounding_box_roi = np.array([map(np.min, nonzero_coords),
                                           map(np.max, nonzero_coords)])
        block_bounding_box_roi[1, :] += 1
    return block_bounding_box_roi


def containing_rois(rois, inner_roi):
    """
    Given a list of rois and an "inner roi" which may or may not be fully
    contained within some of the rois from the list,
    return the subset of rois that entirely envelop the given inner roi.

    Example:
        >>> rois = [([0,0,0], [10,10,10]),
        ...         ([5,3,2], [11,12,13]),
        ...         ([4,6,4], [5,9,9])]
        >>> containing_rois( rois, ( [4,7,6], [5,8,8] ) )
        array([[[ 0,  0,  0],
                [10, 10, 10]],
        <BLANKLINE>
               [[ 4,  6,  4],
                [ 5,  9,  9]]])
    """
    if not rois:
        return np.array([])
    rois = np.asarray(rois)
    left_matches = (rois[:, 0] <= inner_roi[0])
    right_matches = (rois[:, 1] >= inner_roi[1])
    both_matches = np.logical_and(left_matches, right_matches)
    matching_rows = np.logical_and.reduce(both_matches, axis=1).nonzero()
    return rois[matching_rows]


def enlargeRoiForHalo(start, stop, shape, sigma, window=3.5, enlarge_axes=None, return_result_roi=False):
    """
    Enlarge the given roi (start,stop) with a halo according to the given
    sigma and window size, without exceeding the given total image shape given.

    Except for clipping near the image borders, the halo on all sides of the
    image will have width = sigma*window

    start: ROI start coordinate
    stop: ROI stop coordinate
    shape: Total shape of the image (not to be exceeded)
    sigma: The sigma of the filter.
    window: The window size, expressed in units of sigma.
    enlarge_axes: If provided, indicates which axes to expand with the halo.
                  Should be a list of bools (or 1/0 values).
                  For example, halo_axes=(0,1,1,1,0) means: "enlarge roi for axes 1,2,3 but not axes 0,4"
    return_result_roi: If True, also return the "result roi".
                       That is, the roi which you can use to extract the inner data
                       from an array retrieved using the enlarged roi.
                       For example:
                           roi_with_halo, result_roi = enlargeRoiForHalo(start, stop, sigma, return_result_roi=True)
                           outer_data = myfilter(roi_with_halo)
                           data_without_halo = outer_data[roiToSlice(result_roi)]
                           assert data_without_halo.shape == stop - start
    """
    assert len(start) == len(stop) == len(shape)
    shape = TinyVector(shape)
    if enlarge_axes is None:
        enlarge_axes = TinyVector((1,) * len(start))
    else:
        enlarge_axes = TinyVector(enlarge_axes) * 1

    # non-enlarged axes are zero'd out while we enlarge the rest.
    assert len(enlarge_axes) == len(shape)
    max_spatial_shape = enlarge_axes * shape
    spatial_start = enlarge_axes * start
    spatial_stop = enlarge_axes * stop

    if isinstance(sigma, collections.Iterable):
        sigma = TinyVector(sigma)
    if isinstance(start, collections.Iterable):
        ret_type = type(start[0])
    else:
        ret_type = type(start)

    zeros = TinyVector(start) - start

    enlarged_start = np.maximum(spatial_start - np.ceil(window * sigma), zeros).astype(ret_type)
    enlarged_stop = np.minimum(spatial_stop + np.ceil(window * sigma), max_spatial_shape).astype(ret_type)

    # Restore non-halo elements exactly as they were
    enlarged_start += (enlarge_axes == 0) * start
    enlarged_stop += (enlarge_axes == 0) * stop

    enlarged_roi = np.array((enlarged_start, enlarged_stop))
    if return_result_roi:
        inner_roi = np.asarray((start, stop))
        result_roi = inner_roi - enlarged_roi[0]
        return enlarged_roi, result_roi
    else:
        return enlarged_roi


def getIntersectingBlocks(blockshape, roi, asarray=False):
    """
    Returns the start coordinate of each block that the given roi intersects.
    By default, returned as an array of shape (N,M) (N indexes with M coordinates each).
    If asarray=True, then the blocks are returned as an array of shape (D1,D2,D3,...DN,M)
    such that coordinates of spatially adjacent blocks are returned in adjacent entries of the array.

    (SEE ALSO: ``lazyflow.utility.blockwise_view``)

    For example:

    >>> block_starts = getIntersectingBlocks( (10, 20), [(15, 25),(23, 40)] )
    >>> block_starts.shape
    (2, 2)
    >>> print block_starts
    [[10 20]
     [20 20]]

    >>> block_starts = getIntersectingBlocks( (10, 20), [(15, 25),(23, 41)] )
    >>> block_starts.shape
    (4, 2)
    >>> print block_starts
    [[10 20]
     [10 40]
     [20 20]
     [20 40]]

    Now the same two examples, with asarray=True.  Note the shape of the result.

    >>> block_start_matrix = getIntersectingBlocks( (10, 20), [(15, 25),(23, 40)], asarray=True )
    >>> block_start_matrix.shape
    (2, 1, 2)
    >>> print block_start_matrix
    [[[10 20]]
    <BLANKLINE>
     [[20 20]]]

    >>> block_start_matrix = getIntersectingBlocks( (10, 20), [(15, 25),(23, 41)], asarray=True )
    >>> block_start_matrix.shape
    (2, 2, 2)
    >>> print block_start_matrix
    [[[10 20]
      [10 40]]
    <BLANKLINE>
     [[20 20]
      [20 40]]]


    This function works for negative rois, too.

    >>> block_starts = getIntersectingBlocks( (10, 20), [(-10, -5),(5, 5)] )
    >>> print block_starts
    [[-10 -20]
     [-10   0]
     [  0 -20]
     [  0   0]]
    """
    assert len(blockshape) == len(roi[0]) == len(roi[1]), "blockshape and roi are mismatched: {} vs {}".format(
        blockshape, roi)
    roistart = TinyVector(roi[0])
    roistop = TinyVector(roi[1])
    blockshape = TinyVector(blockshape)

    block_index_map_start = roistart / blockshape
    block_index_map_stop = (roistop + (
    blockshape - 1)) / blockshape  # Add (blockshape-1) first as a faster alternative to ceil()
    block_index_map_shape = block_index_map_stop - block_index_map_start

    num_axes = len(blockshape)
    block_indices = np.indices(block_index_map_shape)
    block_indices = np.rollaxis(block_indices, 0, num_axes + 1)
    block_indices += block_index_map_start

    # Multiply by blockshape to get the list of start coordinates
    block_indices *= blockshape

    if asarray:
        return block_indices
    else:
        # Reshape into N*M matrix for easy iteration
        num_indexes = np.prod(block_indices.shape[0:-1])
        axiscount = block_indices.shape[-1]
        return np.reshape(block_indices, (num_indexes, axiscount))


def getIntersectingRois(dataset_shape, blockshape, roi, clip_blocks_to_roi=True):
    block_starts = getIntersectingBlocks(blockshape, roi)
    block_rois = map(partial(getBlockBounds, dataset_shape, blockshape), block_starts)
    if clip_blocks_to_roi:
        block_rois = map(lambda block_roi: getIntersection(block_roi, roi), block_rois)
    return block_rois


def is_fully_contained(inner_roi, outer_roi):
    inner_roi = np.asarray(inner_roi)
    return (inner_roi[0] >= outer_roi[0]).all() and (inner_roi[1] <= outer_roi[1]).all()


def getBlockBounds(dataset_shape, block_shape, block_start):
    """
    Given a block start coordinate and block shape, return a roi for
    the whole block, clipped to fit within the given dataset shape.

    >>> getBlockBounds( [35,35,35], [10,10,10], [10,20,30] )
    (array([10, 20, 30]), array([20, 30, 35]))
    """
    assert (np.mod(block_start, block_shape) == 0).all(), \
        "Invalid block_start: {}.  Must be a multiple of the block shape: {}" \
            .format(block_start, block_shape)

    entire_dataset_roi = roiFromShape(dataset_shape)
    block_shape = TinyVector(block_shape)
    block_bounds = (block_start, block_start + block_shape)

    # Clip to dataset bounds
    block_bounds = getIntersection(block_bounds, entire_dataset_roi)
    return block_bounds


def determineBlockShape(max_shape, target_block_volume):
    """
    Choose a blockshape that is close to the target_block_volume (in pixels),
    without exceeding max_shape in any dimension.

    The resulting blockshape will be as isometric as possible, except
    for those dimensions must be restricted due to small max_shape.  If possible, the
    block's other dimensions are expanded to achieve a total volume of target_block_volume.

    >>> determineBlockShape( (1000,2000,3000,1), 1e6 )
    (100, 100, 100, 1)

    >>> determineBlockShape( (1,100,5,200,3), 1000 )
    (1, 8, 5, 8, 3)
    """
    assert (TinyVector(max_shape) > 0).all(), "Invalid max_shape: {}".format(max_shape)
    ndims = len(max_shape)

    # Attach indexes to remember where each max_shape element came from
    max_with_index = zip(max_shape, range(len(max_shape)))

    # Sort from smallest to largest, to ensure that larger dimensions
    #   will pick up the slack that smaller dims can't accommodate.
    sorted_max = sorted(max_with_index)

    volume_so_far = 1
    block_shape = []

    for (m, i), num_remaining_axes in zip(sorted_max, range(ndims, 0, -1)):
        # Make a block_shape that is isometric in the remaining dimensions
        remaining_factor = target_block_volume // volume_so_far
        block_side = int(pow(remaining_factor, 1.0 / num_remaining_axes) + 0.5)
        block_side = min(block_side, m)
        block_shape.append(block_side)
        volume_so_far *= block_side

        # Sort block_shape dimensions back to the original axis order
    index_order = zip(*sorted_max)[1]
    indexed_block_shape = zip(index_order, block_shape)
    block_shape = zip(*sorted(indexed_block_shape))[1]
    return tuple(block_shape)


def determine_optimal_request_blockshape(max_blockshape, ideal_blockshape, ram_usage_per_requested_pixel, num_threads,
                                         available_ram):
    """
    Choose a blockshape for requests subject to the following constraints:
    - not larger than max_blockshape in any dimension
    - not too large to run in parallel without exceeding available ram (according to num_threads and available_ram)

    Within those constraints, choose the largest blockshape possible.
    The blockshape will be chosen according to the following heuristics:
    - If any dimensions in ideal_blockshape are 0, prefer to expand those first until max_blockshape is reached.
    (The result is known as atomic_blockshape.)
    - After that, attempt to expand the blockshape by incrementing a dimension according to its width in atomic_blockshape.

    Note: For most use-cases, the ``ram_usage_per_requested_pixel`` parameter refers to the ram consumed when requesting ALL channels of an image.
          Therefore, you probably want to omit the channel dimension from your max_blockshape and ideal_blockshape parameters.

    >>> determine_optimal_request_blockshape( (1000,1000,100), (0,0,1), 4, 10, 1e6 )
    (158, 158, 1)

    >>> determine_optimal_request_blockshape( (1000,1000,100), (0,0,1), 4, 10, 1e9 )
    (1000, 1000, 24)

    """
    assert len(max_blockshape) == len(ideal_blockshape)

    # Convert to numpy for convenience.
    max_blockshape = np.asarray(max_blockshape)
    ideal_blockshape = np.asarray(ideal_blockshape)

    target_block_volume_bytes = available_ram // num_threads
    target_block_volume_pixels = target_block_volume_bytes / ram_usage_per_requested_pixel

    # Replace 0's in the ideal_blockshape with the corresponding piece of max_blockshape
    complete_ideal_blockshape = np.where(ideal_blockshape == 0, max_blockshape, ideal_blockshape)

    # Clip to max
    clipped_ideal_blockshape = np.minimum(max_blockshape, complete_ideal_blockshape)

    atomic_blockshape = determineBlockShape(clipped_ideal_blockshape, target_block_volume_pixels)
    atomic_blockshape = np.asarray(atomic_blockshape)

    if np.prod(clipped_ideal_blockshape) >= target_block_volume_pixels:
        # Target volume is too small for us to stack the atomic blockshape, anyway
        return tuple(atomic_blockshape)

    # Need to stack the ideal_blockshape to come up with something larger.
    # Start with an isotropic block, clipped to the nearest multiple of the atomic_blockshape
    blockshape = np.array(determineBlockShape(max_blockshape, target_block_volume_pixels))
    blockshape -= blockshape % atomic_blockshape

    while True:
        # Find a dimension of atomic_blockshape that isn't already maxed out,
        # And see if we have enough RAM to
        candidate_blockshapes = []
        for index in range(len(blockshape)):
            # If we were to expand the blockshape in this dimension, would the block still fit in RAM?
            candidate_blockshape = blockshape.copy()
            candidate_blockshape[index] += clipped_ideal_blockshape[index]
            if (candidate_blockshape <= max_blockshape).all() and \
                    (np.prod(candidate_blockshape) < target_block_volume_pixels):
                candidate_blockshapes.append(candidate_blockshape)

        if len(candidate_blockshapes) == 0:
            break

        def normalized_surface_area(shape):
            pairs = np.array(list(combinations(shape, 2)))
            surface_area = 2 * (pairs[:, 0] * pairs[:, 1]).sum()
            volume = np.prod(shape)
            return surface_area / volume

        # Choose the best among the canidates
        scores = map(normalized_surface_area, candidate_blockshapes)
        (best_shape, best_score) = min(zip(candidate_blockshapes, scores), key=lambda (shape, score): score)
        blockshape = best_shape

    return tuple(blockshape)


def slicing_to_string(slicing, max_shape=None):
    """
    Returns a string representation of the given slicing, which has been
    formatted with spaces so that multiple such slicings could be printed
    in rows with their columns lined up.

    slicing: A tuple or slice objects, e.g. as returned by numpy.s_
    max_shape: If provided, used to determine how many columns to use
               for each slicing field.
    """
    if max_shape:
        max_digits = map(lambda s: int(log10(s)) + 1, max_shape)
    else:
        max_digits = [0] * len(slicing)
    slice_strings = []
    for i, sl in enumerate(slicing):
        s = ""
        s += ("{:" + str(max_digits[i]) + "}").format(sl.start)
        s += " : "
        s += ("{:" + str(max_digits[i]) + "}").format(sl.stop)
        assert sl.step is None
        slice_strings.append(s)
    return "(" + ",  ".join(slice_strings) + ")"
