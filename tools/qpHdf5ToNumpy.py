#!/usr/bin/env python2

#
# QPATH2 - a quantitative pathology toolkit
#
# (c) 2017 Vlad Popovici
#

"""
Extracts an ND array from a HDF5 file and saves it as a NumPy object.
"""
from __future__ import print_function, division, with_statement

import argparse as opt
import numpy as np
import h5py

def main():
    p = opt.ArgumentParser(description="Extracts an ND-array into a NumPy object.")
    p.add_argument('infile', action='store', help='input HDF5 file')
    p.add_argument('object', action='store', help='the name of the object in the input (and output) file(s)')
    p.add_argument('outfile', action='store', help='output NPZ file')

    args = p.parse_args()

    with h5py.File(args.infile, mode='r') as f:
        x = f[args.object][()].squeeze()

    kwds = {args.object: x}
    np.savez_compressed(args.outfile, **kwds)

    return


if __name__ == '__main__':
    main()