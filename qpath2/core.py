#
# QPATH2 - a quantitative pathology toolkit
#
# (c) 2017 Vlad Popovici
#

"""QPATH2.CORE: Core classes and functions.

Defines exception classes and other basic classes.
"""

__all__ = ['Error', 'WSIInfo', 'MRIBase', 'MRI',
           'MRIExplorer', 'MRISlidingWindow']


import openslide as osl
import abc
import numpy as np


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

