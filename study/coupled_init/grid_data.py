#!/usr/bin/env python3

import os
import sys
import argparse

sys.path.append('../../src/thermal')
from open import dataset as open_dataset


def strip_Temp(fp):
    return float(fp.split('Tma_')[-1].split('_')[0])

def strip_δb(fp):
    return float(fp.split('MB_')[-1].split('_')[0])

def main(argv):
    #---------------------------------------------------------------------------
    # Specify command line arguments
    #---------------------------------------------------------------------------
    parser = argparse.ArgumentParser()
    parser.add_argument("in_fn", metavar="path", type=str,
                        help = "Path to .nc file to be gridded enclose in quotes")
    parser.add_argument('-out_fn','--output_filename', type=str,
                        help = "full path to the output figure")

    args, _ = parser.parse_known_args(argv)

    # Raise error input path does not exists
    if not os.path.exists(args.in_fn):
        raise OSError('\n "in_fn" is invalid \n')
    else:
        in_fn = args.in_fn

    # optinal output filename
    if args.output_filename is None:
        # strip file extension from input file
        base_fn = in_fn.split('.nc')[0]
        # added "gridded" extension to in fn if no out_fn is provided
        out_fn  = base_fn + "_gridded.nc"
    else:
        out_fn = args.output_filename

    # find ensemble parameter values
    T_ma = strip_Temp(in_fn)
    δ_b  = strip_δb(in_fn)

    # open and preprocess (i.e grid) the dataset
    ds = open_dataset(in_fn)

    # expand the dimmension to allow concatetnation along
    ds = ds.expand_dims("T_ma").assign_coords(T_ma=('T_ma', [T_ma]))
    ds = ds.expand_dims("Delta_MB").assign_coords(Delta_MB=('Delta_MB', [δ_b]))

    # write the gridded dataset to disk for future use
    ds.to_netcdf(out_fn, "w")

if __name__ == '__main__':
    main(sys.argv[1:])
