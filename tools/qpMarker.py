#!/usr/bin/env python2

#
# QPATH2 - a quantitative pathology toolkit
#
# (c) 2017 Vlad Popovici
#

"""
Mark points or regions in an image. The program needs two images: one RGB image on which
the points/regions will be marked and another one either B/W or gray scale containing
information about points/region of interest. This information depends, of course, on the
kind of marks to be created.
"""

from __future__ import print_function, division, with_statement

import argparse as opt
import numpy as np
import re

from skimage.io import imsave, imread
from skimage.draw import circle, set_color

def main():
    p = opt.ArgumentParser(description="Marks points and regions in a color image.")
    p.add_argument('img_file', action='store', help='target image (RGB)')
    p.add_argument('mask', action='store', help='a B/W or gray scale mask image for points/regions of interest')
    p.add_argument('-o', '--out', action='store', help='name of the result image', default='out.png')
    p.add_argument('-m', '--mark', action='store', choices=['disc', 'cross'],
                   help='type of marking to be used')
    p.add_argument('-c', '--color', action='store', help='color to use for drawing (R,G,B)',
                   default='(255,0,0)')
    p.add_argument('-s', '--size', action='store', default=5, type=int,
                   help='size (in pixels) of the marking')
    p.add_argument('-t', '--threshold', action='store', type=int, default=0,
                   help='all the values >threshold in the mask will be marked in the image')

    args = p.parse_args()

    img = imread(args.img_file)
    msk = imread(args.mask, as_grey=True)

    rx = re.compile(r'(\d+,\d+,\d+)')
    clr = rx.findall(args.color)[0].split(',')
    if len(clr) != 3:
        raise RuntimeError('Incorrect color specification: ' + args.color)

    limit = lambda _x_, mn=0, mx=255: mn if _x_ < mn else mx if _x_ > mx else _x_
    rgb = lambda _x_: limit(int(abs(float(_x_))))
    mk_color = np.array([rgb(clr[0]), rgb(clr[1]), rgb(clr[2])])

    if img.shape[:2] != msk.shape[:2]:
        raise RuntimeError("The image and the mask must have the same extent.")

    [r, c] = np.where(msk > args.threshold)

    if args.mark == 'disc':
        for y, x in zip(r, c):
            rr, cc = circle(y, x, args.size, shape=img.shape)
            img = set_color(img, (rr, cc), mk_color)
    elif args.mark == 'cross':
        pass # to be done
    else:
        raise RuntimeError('Unkown mark shape requested')

    imsave(args.out, img)

    return


if __name__ == "__main__":
    main()