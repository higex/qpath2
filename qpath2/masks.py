#
# QPATH2 - a quantitative pathology toolkit
#
# (c) 2017 Vlad Popovici
#

#
# QPATH2.MASKS - various functions for creating and manipulating image
# masks (i.e. binary images of 0s and 1s).
#

__all__ = ['add_region', 'masked_points', 'apply_mask']

import numpy as np
from skimage.draw import polygon
import vigra

##-
def add_region(mask, poly_line):
    """Add a new masking region by setting to 1 all the
    pixels within the boundaries of a polygon. The changes are
    operated directly in the array.

    Args:
        mask (numpy.array): an array possibly already containing
            some masked regions, to be updated
        poly_line (numpy.array): an N x 2 array with the (x,y)
            coordinates of the polygon vertices as rows

    Returns:
        a numpy.array - the updated mask
    """

    c, r = masked_points(poly_line, mask.shape)
    mask[r, c] = 1

    return mask
##-


##-
def masked_points(poly_line, shape):
    """Compute the coordinates of the points that are inside the polygonal
    region defined by the vertices of the polygon.

    Args:
        poly_line (numpy.array): an N x 2 array with the (x,y)
            coordinates of the polygon vertices as rows
        shape (pair): (width, height) of the rectangular region
            within which the polygon lies (typically image.shape[:2])
    Returns:
        a pair of lists (X, Y) where X[i], Y[i] are the coordinates of a
        point within the mask (polygonal region)
    """

    # check the last point to match the first one
    if (poly_line[0, 0] != poly_line[-1, 0]) or (poly_line[0, 1] != poly_line[-1, 1]):
        np.vstack((poly_line, poly_line[0,]))

    # remeber: row, col in polygon()
    r, c = polygon(poly_line[:,1], poly_line[:,0], shape)

    return c, r
##-


##-
def apply_mask(img, mask):
    """Apply a mask to each channel of an image. Pixels corresponding to 0s in
    the mask will be set to 0. Changes are made in situ.

    Args:
        img (numpy.array or vigra.VigraArray): an image as an N-dim array
            (height x width x no_of_channels)
        mask (numpy.array): a mask as a 2-dim array (height x width)

    Return:
        numpy.array: the modified image
    """
    if mask.dtype is np.bool:
        mask = mask.astype(np.uint8)
    mask[mask > 0] = 1

    if isinstance(img, vigra.VigraArray):
        ch_iter = img.channelIter()
        for c in ch_iter:
            c *= mask
    else:
        if img.ndim == 2:
            img *= mask
        else:
            for k in np.arange(img.shape[2]):
                img[:,:,k] *= mask

    return img
##-
