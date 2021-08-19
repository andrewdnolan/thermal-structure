#!/usr/bin/env python3

"""
read_result.py

    author: andrewdnolan
    date  : 06/11/21

    https://westgrid.github.io/trainingMaterials/materials/xarray20200930.pdfs
    http://cfconventions.org/Data/cf-conventions/cf-conventions-1.8/cf-conventions.html#parametric-v-coord
    http://cfconventions.org/Data/cf-conventions/cf-conventions-1.8/cf-conventions.html#parametric-vertical-coordinate


    https://medium.com/pangeo/vertical-coordinate-transformation-with-pangeo-have-some-pancakes-and-let-xgcm-do-the-work-b0056604d346
"""

import os
import numpy as np
import xarray as xr

from result import *

if __name__ == "__main__":

    # Input files
    # nodes_fp = "mesh_DB/mesh.nodes"
    # result_fp = "mesh_DB/Uncoupled.result"
    result_fp = "./results/mesh/Accumulation_Flux.result"
    nodes_fp  = "./results/mesh/mesh.nodes"
    uncoupled = result_file(result_fp=result_fp, nodes_fp=nodes_fp)

    # # This should be the __repr__ of the result_file class, but haven't gotten there yet.
    # print("")
    # print("Finished reading {}".format(result_fp))
    # print("")
    # print("Coords from mesh.nodes file:")
    # print("----------------------------")
    # for var in uncoupled.coords.dtype.names:
    #     print("  {:<20} ({})".format(var, len(uncoupled.coords[var])))
    # print("Variables from .result file:")
    # print("----------------------------")
    # for var in uncoupled.data.dtype.names:
    #     print("  {:<20} ({})".format(var, len(uncoupled.data[var])))
    # print("")
    # print("")
    #
    # ############################################################################
    # # x-array creation method one
    # ############################################################################
    # dim_x = np.unique(uncoupled.coords["x"])
    # # y is really the "second" coordinate, not the cartesian coordinate $y$
    # dim_z = np.unique(uncoupled.coords["y"])
    #
    # # x-array coords not the attribute coords
    # coords = {"x": dim_x, "z": dim_z}
    #
    # ds = xr.Dataset(
    #     {
    #         var: xr.DataArray(
    #             uncoupled.data[var].reshape(
    #                 uncoupled.mesh_params["Ny"], uncoupled.mesh_params["Nx"]
    #             ),
    #             dims=("z", "x"),
    #             coords=coords,
    #         )
    #         for var in uncoupled.data.dtype.names
    #     }
    # )
    #
    # ds.to_netcdf("method_1.nc")
    # ds.close()
    # # del dim_x, dim_z, coords, ds
    #
    # ############################################################################
    # # x-array creation method two
    # ############################################################################

    dim_x = uncoupled.coords["x"].reshape(
        uncoupled.mesh_params["Ny"], uncoupled.mesh_params["Nx"]
    )
    # y is really the "second" coordinate, not the cartesian coordinate $y$
    dim_z = uncoupled.coords["y"].reshape(
        uncoupled.mesh_params["Ny"], uncoupled.mesh_params["Nx"]
    )

    # x-array coords not the attribute coords
    coords = {
        "t": np.arange(uncoupled.mesh_params["NT"]),
        "coord_1": np.unique(uncoupled.coords["x"]),
        "coord_2": np.unique(uncoupled.coords["y"]),
        "x": (("coord_2", "coord_1"), dim_x),
        "z": (("coord_2", "coord_1"), dim_z),
    }

    ## This is probably the most immediate solution to the problem:
    #       https://stackoverflow.com/questions/57202976/how-to-deal-with-time-dependent-coordinates-in-xarray

    ds = xr.Dataset(
        {
            var: xr.DataArray(
                uncoupled.data[var].reshape(
                    uncoupled.mesh_params["Ny"],
                    uncoupled.mesh_params["Nx"],
                    uncoupled.mesh_params["NT"]
                ),
                dims=("coord_2", "coord_1", "t"),
                coords=coords,
            )
            for var in uncoupled.data.dtype.names
        }
    )

    # Do the vertical coordinate transform
    #ds.coords["z"] = (ds.z * ds.depth) + (ds.freesurf - ds.depth)
    ds.coords["z"] = (ds.z * ds.height) + (ds.zbed)

    # # ## proof that this works!
    # import matplotlib.pyplot as plt
    #
    # plt.pcolormesh(ds.x, ds.z.isel(t=0), np.sqrt(ds["velocity 1"].isel(t=0)**2 + ds["velocity 2"].isel(t=0)**2), edgecolor='k', shading="gouraud")
    # plt.colorbar()
    # plt.tight_layout()
    # plt.show()
    #
    # ds.to_netcdf("method_2.nc")
    # ds.close()
