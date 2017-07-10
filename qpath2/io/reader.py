#
# QPATH2 - a quantitative pathology toolkit
#
# (c) 2017 Vlad Popovici
#

# QPATH2.IO.READER: various functions for reading whole slide images.

__all__ = ["openslide_read_region_px"]

import numpy as np

from qpath2.core import WSIInfo, Error
from qpath2.io.io_ import osl_read_region_


##-
def openslide_read_region_px(wsi, x0, y0, width, height, level):
    """Read a region of a WSI calling OpenSlide's corresponding C function.

    Args:
        wsi (WSIInfo): meta-data about the slide
        x0, y0 (long): top left corner of the region (in pixels, at level 0)
        width, height (long): width and height (in pixels) of the region
        level (int): the magnification level to read from

    Returns:
        numpy.ndarray (w x h x 4) with dtype=numpy.uint8
    """

    if level >= wsi.info['level_count']:
        raise Error("requested level does not exist")

    if level > 0:
        # top-left corner in level-coordinates:
        x1 = x0 / wsi.info['levels'][level]['downsample_factor']
        y1 = y0 / wsi.info['levels'][level]['downsample_factor']
    else:
        x1, y1 = x0, y0

    # check bounds:
    if x1 >= wsi.info['levels'][level]['x_size'] or \
        y1 >= wsi.info['levels'][level]['y_size'] or \
        x1 + width > wsi.info['levels'][level]['x_size'] or \
        y1 + height > wsi.info['levels'][level]['y_size']:
        raise Error("region out of layer's extent")

    x0, y0, width, height = [long(_x) for _x in [x0, y0, width, height]]
    img = np.zeros((height, width, 4), dtype=np.uint8)
    r = osl_read_region_(wsi.path, img, x0, y0, width, height, level)

    if r != 0:
        raise Error("low-level error in osl_read_region", code=r)

    return img[...,(2,1,0,3)]   # change BGRA into RGBA
##-



