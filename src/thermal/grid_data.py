#!/usr/bin/env python3

import os
import sys
import json
import argparse

# functions from local src code
from open import dataset as open_dataset

def main(argv):
    #---------------------------------------------------------------------------
    # Specify command line arguments
    #---------------------------------------------------------------------------
    parser = argparse.ArgumentParser()
    parser.add_argument("in_fn", metavar="path", type=str,
                        help = "Path to .nc file to be gridded enclose in quotes")
    parser.add_argument('-out_fn','--output_filename', type=str,
                        help = "full path to the output figure")
    parser.add_argument('-params','--parameter_dict', type=str,
                        help = "Dictionary of parameters to add to NetCDF object"\
                               "as dimmension. e.g. \'{\"T_ma\":-9.00, \"Delta_MB\":-1.5}\'")


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

    # Check is parameter dictionary was passed
    if args.parameter_dict is not None:
        # should probably add some better error handling here
        param_dict = json.loads(args.parameter_dict)
    else:
        # set to empty dict, loop below will not itterate
        param_dict = dict()

    # open and preprocess (i.e grid) the dataset
    ds = open_dataset(in_fn)

    # loop over ensemble parameter values
    for key in param_dict:
        val = param_dict[key]
        # expand the dimmension along parameter values to concatentate along
        ds = ds.expand_dims(key).assign_coords({key: (key, [val])})

    # write the gridded dataset to disk for future use
    ds.to_netcdf(out_fn, "a")

if __name__ == '__main__':
    main(sys.argv[1:])
