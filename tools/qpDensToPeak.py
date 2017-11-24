#!/usr/bin/env python2

#
# QPATH2 - a quantitative pathology toolkit
#
# (c) 2017 Vlad Popovici
#

"""
Creates a peak map from a density map.
"""

from __future__ import print_function, division, with_statement

import argparse as opt
import numpy as np

from skimage.io import imsave, imread
from skimage.morphology import extrema

def main():
    p = opt.ArgumentParser(description="Detects local peaks in a density/intensity map.")
    p.add_argument('dmap', action='store', help='density/intensity map as a numpy array (.npz) file')
    p.add_argument('-r', '--radius', action='store', type=int, default=5,
                   help='radius of the local neighborhood')

    args = p.parse_args()

    dmap = np.load(args.dmap)
    if dmap.ndim > 2:
        print('Density map with more than 2 dimensions, using first 2.')
        dmap = dmap[:,:,...].squeeze()

    peaks = np.zeros(dmap.shape)



if __name__ == '__main__':
    main()