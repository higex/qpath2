#
# QPATH2 - a quantitative pathology toolkit
#
# (c) 2017 Vlad Popovici
#

"""IO.TILED: tiled storage for large images.

 Instead of storing a multiresolution image in a BigTIFF file, store it in
 the file system as a hierarchy of folders. The basic hierarchy looks like:

.../image name/
              +---- meta.json    <- meta data about the file
              +---- first downsampling level/
                        +---- meta.json
                        | tile_i_j.ppm...
              +---- second downsampling level/
                        +---- meta.json
                        | tile_i_j.ppm...
              ...etc...

 The image type (here .ppm) can be changed.
"""

from __future__ import (absolute_import, division, print_function, unicode_literals)

__all__ = ['save_tiled_image', 'load_tiled_image']

from math import floor
import os
import os.path
import shutil
import simplejson as json
import numpy as np
from skimage.io import imread, imsave

from qpath2.core import MRIBase


##-
def save_tiled_image(img, root, level, tile_geom, img_type="jpeg"):
    """Save an image as a collection of tiles.

    The image is split into a set of fixed-sized (with the exception of right-most and
    bottom-most) tiles.

    *WARNING*: any existing tiles in the path root/level will be deleted!

    Args:
        img (numpy array): an image in OpenCV ordering (BGR). Alpha channel is not
            supported
        root (string): root folder of the image storing hierarchy. The tiles will be
            stored into root/level folder
        level (int): the magnification level
        tile_geom (tuple): (width, height) of the tile
        img_type (string, optional): file type for the tiles

    Returns:
        dict: a dictionary with meta-data about the tiles and original image
    """
    assert(img.ndim == 2 or (img.ndim == 3 and img.shape[2] <= 3))

    n_channels = 1 if img.ndim == 2 else img.shape[2]
    dst_path = root + os.path.sep + 'level_{:d}'.format(level)

    tg = (min(tile_geom[0], img.shape[1]), min(tile_geom[1], img.shape[0]))
    nh = int(floor(img.shape[1] / tg[0])) + (1 if img.shape[1] % tg[0] != 0 else 0)
    nv = int(floor(img.shape[0] / tg[1])) + (1 if img.shape[0] % tg[1] != 0 else 0)

    tile_meta = dict({'level': level,
                      'level_image_width': img.shape[1],
                      'level_image_height': img.shape[0],
                      'level_image_nchannels': 1 if img.ndim == 2 else img.shape[2],
                      'n_tiles_horiz': nh,
                      'n_tiles_vert': nv,
                      'tile_width': tg[0],
                      'tile_height': tg[1]})

    if os.path.exists(dst_path):
        shutil.rmtree(dst_path)
    os.mkdir(dst_path)

    for i in range(nv):
        for j in range(nh):
            i0, j0 = i * tg[1], j * tg[0]
            i1, j1 = min((i + 1) * tg[1], img.shape[0]), min((j + 1) * tg[0], img.shape[1])
            if n_channels == 1:
                im_sub = img[i0:i1, j0:j1]
            else:
                im_sub = img[i0:i1, j0:j1, :]
            tile_meta['tile_' + str(i) + '_' + str(j)] = dict(
                {'name': dst_path + '/tile_' + str(i) + '_' + str(j) + '.' + img_type,
                 'i': i, 'j': j,
                 'x': j0, 'y': i0})
            imsave(dst_path + os.path.sep + 'tile_' + str(i) + '_' + str(j) + '.' + img_type, im_sub)

    with open(dst_path + os.path.sep + 'meta.json', 'w') as fp:
        json.dump(tile_meta, fp, separators=(',', ':'), indent='  ', sort_keys=True)

    return tile_meta
##-end


##-
def load_tiled_image(img_meta):
    """Load a tiled image. All the information about the tile geometry and tile paths is
     taken from img_meta.

    Args:
        img_meta (dict): a descriptor for the tiled image with at least the following
            keys (see save_tiled_image)
                level_image_width
                level_image_height
                level_image_nchannels
                n_tiles_horiz
                n_tiles_vert
            and for each tile, an entry as
                'tile_i_j' which is a dict with keys:
                i
                j
                name
                x
                y

    Returns:
        a numpy.ndarray
    """
    img_w, img_h = long(img_meta['level_image_width']), long(img_meta['level_image_height'])
    nh, nv = long(img_meta['n_tiles_horiz']), long(img_meta['n_tiles_vert'])

    img = np.zeros((img_h, img_w, 3), dtype=np.uint8)

    for i in range(nv):
        for j in range(nh):
            tile_id = 'tile_'+str(i)+'_'+str(j)
            tile = imread(img_meta[tile_id]['name']).astype(np.uint8)
            # the tile might not have the regular default shape, so it's better to use the
            # tile's shape than 'tile_width' and 'tile_height'
            x, y = long(img_meta[tile_id]['x']), long(img_meta[tile_id]['y'])
            img[x:x+tile.width, y:y+tile.height, :] = tile

    return img
##-


##-
class TiledImage(object):
    """A tiled image, loading regions on demand.

    """
    def __init__(self, path):
        pass
##-