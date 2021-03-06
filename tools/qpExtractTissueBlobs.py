#!/usr/bin/env python2

#
# QPATH2 - a quantitative pathology toolkit
#
# (c) 2017 Vlad Popovici
#

"""
Extract images of tissue blobs from a whole slide image. In many cases, on a single slide,
there are several tissue sections (e.g. stained with various antibodies) which need to be
extracted and processed.
"""
from __future__ import print_function, division, with_statement

import openslide as osl
import argparse as opt
import numpy as np
import simplejson as json
import subprocess as sp
import re
import os, os.path

from skimage.transform import resize
from skimage.morphology import remove_small_objects, convex_hull_object
from skimage.filters import threshold_otsu
from skimage.color import rgb2gray
from skimage.measure import label, regionprops
from skimage.external import tifffile
from skimage.io import imsave, imread

from qpath2.core import WSIInfo, MRI
from qpath2.io.tiled import save_tiled_image

import warnings

warnings.simplefilter("ignore", UserWarning)


def main():
    global img_file, res_prefix, s_factors, tile_geom, res_format, n_levels, res_xy

    p = opt.ArgumentParser(description="Extracts tissue blobs from a WSI and stores them as tile images.")
    p.add_argument('img_file', action='store', help='WSI file name')
    p.add_argument('--prefix', action='store', help='path where to store the results', default='./')
    p.add_argument('--min_area', action='store', help='area of the smallest object to keep (in px)',
                   default=4096)
    p.add_argument('--level', action='store', type=int,
                   help='magnification level (0: maximum resolution, default: lowest)',
                   default=-1)
    p.add_argument('-k', '--keep_whole_image', action='store_true',
                   help='keep the whole image as well?')
    p.add_argument('--tile', action='store', help='tile geometry: (w,h)', default='(128,128)')
    p.add_argument('--format', action='store', help='output image format',
                   choices=['ppm', 'tiff', 'jpeg'], default='ppm')
    p.add_argument('-v', '--verbose', action='store_true', help='verbose')
    p.add_argument('-m', '--mask', action='store_true', help='save mask at image resolution?')
    p.add_argument('-n', '--names', action='store', help='a list of tissue blob names (e.g. stains)', nargs='+')

    args = p.parse_args()
    rx = re.compile(r'(\d+,\d+)')
    h, v = rx.findall(args.tile)[0].split(',')
    tile_geom = (int(abs(float(h))), int(abs(float(v))))


    print('Working on ' + args.img_file)
    img = WSIInfo(args.img_file)
    mri = MRI(img)
    #img = osl.OpenSlide(args.img_file)
    #meta = {'objective': img.properties[osl.PROPERTY_NAME_OBJECTIVE_POWER],
    #        'mpp_x': float(img.properties[osl.PROPERTY_NAME_MPP_X]),
    #        'mpp_y': float(img.properties[osl.PROPERTY_NAME_MPP_Y])}
    meta = {'objective': img.info['objective'],
            'mpp_x': img.info['x_mpp'],
            'mpp_y': img.info['y_mpp']}

    if os.path.exists(args.prefix + '/meta.json'):
        with open(args.prefix + '/meta.json', 'r') as fd:
            meta = json.load(fd)

    if args.level == -1 or args.level >= img.info['level_count']:
        args.level = img.info['level_count'] - 1

    # read the lowest resolution and try to detect the tissue pieces
    lowest_res_level = img.info['level_count'] - 1
    #img_full = img.read_region((0, 0), lowest_res_level,
    #                           img.level_dimensions[lowest_res_level])  # mask automatically applied, background is 0
    #img_data = np.asarray(img_full)
    img_data = mri.get_region_px(0, 0,
                                 img.info['levels'][lowest_res_level]['x_size'],
                                 img.info['levels'][lowest_res_level]['y_size'],
                                 lowest_res_level, as_type=np.uint8)

    if img_data.ndim == 3 and img_data.shape[2] > 3:
        img_data = img_data[..., :3]  # drop alpha channel

    if img_data.ndim ==3:
        img_gray = rgb2gray(img_data)
    else:
        img_gray = img_data

    th = threshold_otsu(img_gray)
    img_bin = img_gray >= th  # black background
    # img_bin = dilation(img_bin, selem=disk(3))
    img_bin = remove_small_objects(img_bin, args.min_area)
    img_bin = convex_hull_object(img_bin)

    labels = label(img_bin, neighbors=8, background=0)
    props = regionprops(labels)

    if args.verbose:
        print("Number of tissue blobs: {:d}".format(len(props)))
        imsave(args.prefix + os.path.sep + 'whole_slide.jpeg', img_data)
        imsave(args.prefix + os.path.sep + 'whole_slide_blobs.jpeg', 255*labels)

    # extract tissue regions:
    s = img.info['levels'][lowest_res_level]['downsample_factor'] / \
        img.info['levels'][args.level]['downsample_factor']  # downscale factor for the lowest resolution wrt. desired level
    s0 = img.info['levels'][lowest_res_level]['downsample_factor']  # downscale factor for the lowest resolution wrt. level 0

    mri = None

    # order the regions, from the top-most (smaller y coordinate of the bounding box) to the
    # bottom-most:
    props = sorted(props, key=lambda _pr: _pr.bbox[0])

    k = 0
    for pr in props:
        # get mask:
        msk = img_bin[pr.bbox[0]:pr.bbox[2], pr.bbox[1]:pr.bbox[3]].astype(np.uint8)
        # get image at highest resolution:
        start_x = np.int64(max(s0 * pr.bbox[1], 0))
        start_y = np.int64(max(s0 * pr.bbox[0], 0))
        width = np.int64(s * (pr.bbox[3] - pr.bbox[1]))
        height = np.int64(s * (pr.bbox[2] - pr.bbox[0]))

        # store meta:
        if args.names is not None and k < len(args.names):
            tname = args.names[k]
        else:
            tname = "tissue_{:d}".format(k)

        dst_path = args.prefix + os.path.sep + tname

        if not os.path.exists(dst_path):
            os.mkdir(dst_path)

        meta[tname] = dict({"name": dst_path + os.path.sep + tname + '_level_{:d}.png'.format(args.level),
                            "mask": dst_path + os.path.sep + tname + '_mask_level_{:d}.tiff'.format(args.level),
                            "from_original_level": args.level,
                            "from_original_x": s * pr.bbox[1],
                            "from_original_y": s * pr.bbox[0],
                            "from_original_width": width,
                            "from_original_height": height})

        # Currently (May 2017), a bug in PIL library leads to errors in the case of very
        # large regions. Hence, we rely on the OpenSlide's utility to export the .PNG
        # file. Then re-read it and continue with processing.

        print("Extract tissue blob {:d}".format(k+1))

        # check whether the image already exists - hopefully the correct one, from a
        # previous run
        if not os.path.exists(meta[tname]['name']):
            sp.check_call(["openslide-write-png", args.img_file,
                           '{:d}'.format(start_x),
                           '{:d}'.format(start_y),
                           '{:d}'.format(args.level),
                           '{:d}'.format(width),
                           '{:d}'.format(height),
                           meta[tname]['name']])

        print('reading '+ meta[tname]['name'])
        img_data = imread(meta[tname]['name'])

        # img_tissue = img.read_region((start_x, start_y), args.level,
        #                              (width, height))
        # img_data = mri.get_region_px(s * pr.bbox[1], s * pr.bbox[0], width, height, args.level, as_type=np.uint8)

        msk_from_scanner = None
        if img_data.ndim == 3 and img_data.shape[2] == 4:
            # has alpha-channel
            msk_from_scanner = img_data[:,:,3].astype(np.uint8)  # alpha-channel: 255 for the tissue region, 0 elsewhere
            msk_from_scanner[msk_from_scanner == 255] = 1  # img 0/1
            small_mask = resize(msk_from_scanner, msk.shape, order=0, mode='constant', cval=0, preserve_range=True)
            small_mask[small_mask < 1] = 0
            small_mask *= msk  # 'AND' the masks
            small_mask = small_mask.astype(np.uint8)

            img_data = img_data[:,:,:3]  # drop alpha-channel from image data

        # up-scale the mask and set "True" for foreground
        msk = resize(msk, img_data.shape[:2], mode='constant', order=0, cval=0, preserve_range=True)
        msk[msk < 1] = 0  # drop values due to aliasing

        if msk_from_scanner is not None:
            msk *= msk_from_scanner   # 'AND' the two masks

        msk = msk.astype(np.uint8)

        if img_data.ndim == 2:
            img_data *= msk
        else:
            for ch in np.arange(img_data.shape[2]):
                img_data[:, :, ch] *= msk

        # save
        save_tiled_image(img_data, dst_path, args.level, tile_geom, img_type=args.format)

        # save, anyway, the small mask:
        imsave(dst_path + os.path.sep + tname + '_level_{:d}.ppm'.format(lowest_res_level), 255*small_mask)

        # the large mask is saved only if asked:
        if args.mask:
            with tifffile.TiffWriter(meta[tname]['mask'], bigtiff=True) as tif:
                tif.save(255 * msk, compress=9, tile=(512, 512))

        if not args.keep_whole_image:
           os.remove(meta[tname]['name'])
           meta[tname]['name'] = ''
        # if args.keep_whole_image:
        #     with tifffile.TiffWriter(meta[tname]['name'], bigtiff=True) as tif:
        #         tif.save(img_data, compress=9, tile=(512, 512))

        k += 1
    # end for

    with open(args.prefix + '/meta.json', 'w') as fp:
        json.dump(meta, fp, separators=(',', ':'), indent='  ', sort_keys=True)

    return


## MAIN ##
if __name__ == '__main__':
    main()