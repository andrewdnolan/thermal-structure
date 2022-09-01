#!/usr/bin/env python3

"""
add_attr.py:

add an the contents of a file as a global attribute a NetCDF file.
allows us to archive the sif and params used to generate a file for future reference

there's probably a way to do this from the command line, using `nco`.
I tried but it was takign to long so went brute force with python
"""

import os
import sys
import argparse
import xarray as xr

def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument("nc_file", metavar="path", type=str,
                        help = "Path to .nc file to add attribute to. Can only append for now")
    parser.add_argument('-f','--file', type=str,
                        help = "path file to be added as attribute")
    parser.add_argument('-a','--attr', type=str,
                        help = "name of the attribute to be created")

    args, _ = parser.parse_known_args(argv)

    # Raise error if input path does not exists
    if not os.path.exists(args.nc_file):
        raise OSError('\n "nc_file" does not exist\n')
    else:
        nc_file = args.nc_file

    # Raise error if input path does not exists
    if not os.path.exists(args.file):
        raise OSError('\n "--file" does not exist\n')
    else:
        file = args.file

    # name of the attr
    attr = args.attr

    # open the NetCDF
    src = xr.open_dataset(nc_file)

    # open the attribute file
    with open(file) as f:
        # dump the attr
        src.attrs[attr] = f.read()

    # write the updatee netcdf
    src.to_netcdf(nc_file, 'a')

if __name__ == '__main__':
    main(sys.argv[1:])
