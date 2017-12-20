#!/usr/bin/env python2

#
# QPATH2 - a quantitative pathology toolkit
#
# (c) 2017 Vlad Popovici
#

"""
Creates a peak map from a density/intensity map.
"""

from __future__ import print_function, division, with_statement

import argparse as opt
import numpy as np
from skimage.morphology import extrema, disk

def main():
    p = opt.ArgumentParser(description="Detects local peaks in a density/intensity map.")
    p.add_argument('dmap', action='store', help='density/intensity map file (NumPy compressed .npz)')
    p.add_argument('map_name', action='store', help='name of the map object in the file')
    p.add_argument('outfile', action='store', help='name of the result file (.npz)')
    p.add_argument('-r', '--radius', action='store', type=int, default=5,
                   help='radius of the local neighborhood')
    p.add_argument('-z', '--zero', action='store', default=0.0, type=np.float,
                   help='all values not greater than this will be set to 0 (background)')
    p.add_argument('-n', '--name', action='store', default=None,
                   help='name of the resulting object (to be stored in the .npz file)')

    args = p.parse_args()


    dmap = np.load(args.dmap)[args.map_name]
    if dmap.ndim == 3:
        print('Density map with more than 2 dimensions, using first 2.')
        dmap = dmap[:,:,0]
    elif dmap.ndim > 3:
        raise RuntimeError('Data structure unknown.')

    zero = np.cast[dmap.dtype](args.zero)
    dmap[dmap <= zero] = 0

    peaks = extrema.local_maxima(dmap, selem=disk(args.radius))

    if args.name is None:
        res = {args.map_name: peaks}
    else:
        res = {args.name: peaks}

    np.savez_compressed(args.outfile, **res)

    return


if __name__ == '__main__':
    main()