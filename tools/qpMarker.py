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
    p.add_argument('mask', action='store', help='a binary mask as a 2D numpy array (in .npz file)')
    p.add_argument('outfile', action='store', help='name of the result image')

    p.add_argument('-n', '--name', action='store', default='mask',
                   help='name of the mask object in the .npz file')
    p.add_argument('-m', '--mark', action='store', choices=['disk', 'cross'],
                   help='type of marking to be used')
    p.add_argument('-c', '--color', action='store', help='color to use for drawing (R,G,B)',
                   default='(255,0,0)')
    p.add_argument('-s', '--size', action='store', default=5, type=int,
                   help='size (in pixels) of the marking')
    p.add_argument('-a', '--alpha', action='store', type=float, default=0.5,
                   help='alpha level (<1 transparent marking)')

    args = p.parse_args()

    img = imread(args.img_file)
    msk = np.load(args.mask)[args.name]

    rx = re.compile(r'(\d+,\d+,\d+)')
    clr = rx.findall(args.color)[0].split(',')
    if len(clr) != 3:
        raise RuntimeError('Incorrect color specification: ' + args.color)

    limit = lambda _x_, mn=0, mx=255: mn if _x_ < mn else mx if _x_ > mx else _x_
    rgb = lambda _x_: limit(int(abs(float(_x_))))
    mk_color = np.array([rgb(clr[0]), rgb(clr[1]), rgb(clr[2])])

    if img.shape[:2] != msk.shape[:2]:
        raise RuntimeError("The image and the mask must have the same extent.")

    [r, c] = np.where(msk > 0)

    if args.mark == 'disk':
        for y, x in zip(r, c):
            rr, cc = circle(y, x, args.size, shape=img.shape)
            set_color(img, (rr, cc), mk_color, alpha=args.alpha)
    elif args.mark == 'cross':
        pass # to be done
    else:
        raise RuntimeError('Unkown mark shape requested')

    imsave(args.outfile, img)

    return


if __name__ == "__main__":
    main()