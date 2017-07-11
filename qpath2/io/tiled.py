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

__all__ = ['save_tiled_image', 'load_tiled_image', 'TiledImage']

from math import floor
import os
import os.path
import shutil
import simplejson as json
import numpy as np
from skimage.io import imread, imsave

from qpath2.core import MRIBase, Error


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
    _path = None
    _meta = None

    def __init__(self, path):
        self._path = path
        with open(self._path, 'r') as fp:
            self._meta = json.load(fp)


    @property
    def path(self):
        return self._path

    @property
    def height(self):
        return long(self._meta['level_image_height'])

    @property
    def width(self):
        return long(self._meta['level_image_width'])

    @property
    def level(self):
        return int(self._meta['level'])

    @property
    def tile_count_horizontal(self):
        return int(self._meta['n_tiles_horiz'])

    @property
    def tile_count_vertical(self):
        return int(self._meta['n_tiles_vert'])

    @property
    def tile_count(self):
        return self.tile_count_horizontal * self.tile_count_vertical

    @property
    def tile_width(self):
        return int(self._meta['tile_width'])

    @property
    def tile_height(self):
        return int(self._meta['tile_height'])

    def get_tile(self, i, j):
        """Return the (i,j)-th tile.
        Args:
            i, j (int): tile coordinates

        Returns:
            numpy.ndarray
        """
        img = imread(self._meta['tile_' + str(i) + '_' + str(j)]['name'])
        return img


    def get_image(self):
        """Return the whole image, by loading all the tils."""
        return load_tiled_image(self._meta)


    def get_tile_coverage(self, x, y, width, height):
        """Return the indices (i,j) of the tiles covering a given
        rectangular region.

        Args:
            x, y (long): top-left corner coordinates (column, row)
            width, height (long): region extent

        Returns:
            list of pairs: [(i,j), ...] corresponding to tiles_i_j covering
             the region
        """
        x, y, width, height = [long(_z) for _z in [x, y, width, height]]
        if not (0 <= x < self.width):
            raise Error('x out of bounds')
        if not (0 <= y < self.height):
            raise Error('y out of bounds')
        if x + width > self.width or y + height > self.height:
            raise Error('region too large for the image')

        # Find the tiles covering the requested reqion
        start_i = np.int(np.floor(y / self.tile_height))
        start_j = np.int(np.floor(x / self.tile_width))
        end_i = np.int(np.floor((y + height) / self.tile_height) + \
                (1 if (y + height) % self.tile_height != 0 else 0))
        end_j = np.int(np.floor((x + width) / self.tile_width) + \
                (1 if (x + width) % self.tile_width != 0 else 0))

        ij = [(i, j) for i in np.arange(start_i, end_i) for j in np.arange(start_j, end_j)]

        return ij


    def get_region(self, x, y, width, height):
        """Return an arbitrary region within a tiled image.
        Args:
            x, y (long): top-left corner coordinates (column, row)
            width, height (long): region extent

        Returns:
            numpy.ndarray
        """
        x, y, width, height = [long(_z) for _z in [x, y, width, height]]
        if not (0 <= x < self.width):
            raise Error('x out of bounds')
        if not (0 <= y < self.height):
            raise Error('y out of bounds')
        if x + width > self.width or y + height > self.height:
            raise Error('region too large for the image')

        # Algo:
        # -find the tiles to load
        # -load all the tiles
        # -adjust, if needed, the starting and ending points of the
        #  region
        # This is not optimal from a memory usage perspective, but
        # it's simpler.

        # Find the tiles covering the requested reqion
        start_i = np.int(np.floor(y / self.tile_height))
        start_j = np.int(np.floor(x / self.tile_width))
        end_i = np.int(np.floor((y + height) / self.tile_height) + \
                (1 if (y + height) % self.tile_height != 0 else 0))
        end_j = np.int(np.floor((x + width) / self.tile_width) + \
                (1 if (x + width) % self.tile_width != 0 else 0))

        # Load the tiles start_i:end_i, start_j:end_j
        tile = self.get_tile(start_i, start_j)
        nchannels = 1 if tile.ndim == 2 else 3
        if nchannels == 1:
            img = np.zeros((self.tile_height * (end_i - start_i),
                            self.tile_width * (end_j - start_j)), dtype=np.uint8)
        else:
            img = np.zeros((self.tile_height * (end_i - start_i),
                            self.tile_width * (end_j - start_j),
                            tile.shape[2]), dtype=np.uint8)

        if nchannels == 1:
            for i in range(start_i, end_i):
                for j in range(start_i, end_i):
                    tile = self.get_tile(i, j)

                    # last tile in row and last row of tiles might have non-standard
                    # dimensions, so better use the actual tile shape in computing the
                    # end point:
                    img[(i-start_i)*self.tile_height:(i-start_i)*self.tile_height + tile.shape[0],
                        (j-start_j)*self.tile_width:(j-start_j)*self.tile_width + tile.shape[1]] = tile
        else:
            for i in range(start_i, end_i):
                for j in range(start_j, end_j):
                    tile = self.get_tile(i, j)

                    # last tile in row and last row of tiles might have non-standard
                    # dimensions, so better use the actual tile shape in computing the
                    # end point:
                    img[(i-start_i)*self.tile_height:(i-start_i)*self.tile_height + tile.shape[0],
                        (j-start_j)*self.tile_width:(j-start_j)*self.tile_width + tile.shape[1], :] = tile

        # Adjust image to the requested region:
        if nchannels == 1:
            res = img[y - start_i*self.tile_height : y + height - start_i*self.tile_height,
                  x - start_j * self.tile_width : x + width - start_j * self.tile_width].copy()
        else:
            res = img[y - start_i*self.tile_height : y + height - start_i*self.tile_height,
                  x - start_j * self.tile_width : x + width - start_j * self.tile_width, :].copy()

        return res
##-
