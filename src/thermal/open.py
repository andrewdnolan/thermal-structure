import numpy as np
import xarray as xr
from glob import glob
from tqdm import tqdm
from derived_fields import calc_percent_temperate, calc_volume, calc_magnitude


"""
# open multiple zarr
    https://discourse.pangeo.io/t/how-to-read-multiple-zarr-archives-at-once-from-s3/2564
"""

# dictionary of 1D variables and vertical index they are define along
vars_1D = { "q_lat" : -2,
            "mass balance" : -1,
            "surface_enthalpy" : -1}

# mapping of variables which we seek to rename.
rename_dict = {'velocity 1'   : 'vel_x',
               'velocity 2'   : 'vel_z',
               'strainrate 1' : 'SR_xx',
               'strainrate 2' : 'SR_zz',
               'strainrate 4' : 'SR_xz',
               'strainrate 5' : 'SR_ii',
               'enthalpy heat diffusivity' : 'kappa_cold'}

def __quads_to_tris(quads):
    """converts quad elements into tri elements
    from: https://stackoverflow.com/a/59971611/10221482
    """

    tris = [[None for j in range(3)] for i in range(2*len(quads))]
    for i in range(len(quads)):
        j = 2*i
        n0 = quads[i][0]
        n1 = quads[i][1]
        n2 = quads[i][2]
        n3 = quads[i][3]
        tris[j][0] = n0
        tris[j][1] = n1
        tris[j][2] = n2
        tris[j + 1][0] = n2
        tris[j + 1][1] = n3
        tris[j + 1][2] = n0
    return tris

def __preprocess_UGRID(ds, out_fn=None):
    """Reshape the NetCDF from Elmer onto a structured grid.

    Parameters
    ----------
    ds : `xarray.Dataset`
        Dataset written by the "NetcdfUGRIDOutputSolver.f90" file

    Returns
    -------
    gridded : `xarray.Dataset`
        Dataset which has been gridded and has common coordinate and dimensions
        as the "result2nc" created NetCDF files

    """

    # Get grid info
    NT = ds.time.size                   # number of timesteps
    NN = ds.nMesh_node.size             # number of nodes
    NX = np.unique(ds.Mesh_node_x).size # number nodes in the x direction
    NZ = NN // NX                       # number nodes in the z direction

    # reshape the array (i.e. unflatten), with fortran convetion
    gridded = ds.coarsen(nMesh_node=NX).construct(nMesh_node=("coord_2", "coord_1"))

    # drop variables no longer of use on structured mesh
    # NOTE:
    # - "BulkElement_Area" is dropped since is on the element faces, whereas
    #    all the data is gridded on nodes.
    # - "BulkElement_Area" is only written once, so not updated as mesh is extruded
    gridded = gridded.drop_vars(['Mesh_face_nodes', 'Mesh', 'BulkElement_Area'])

    # rename the coordinate variable
    gridded = gridded.rename({'time' : 't',
                              'Mesh_node_x' : 'X',
                              'Mesh_node_y' : 'Z'})
    return gridded

def __preprocess(ds, h_min=10.0):
    """Wrapper around various preprocessing functions.
    """
    if 'nMesh_node' in ds.dims:
        ds = __preprocess_UGRID(ds)
    else:
        print('no other support implemented yet. To be Done')

    # loop over the 1D variables, and only store the pertinent info
    for key in vars_1D:
        # double check the variable is in the source file
        if key in ds:
            # get the vertical index the variables is defined a
            index = vars_1D[key]
            ds[key] = ds[key].isel(coord_2=index)

    # rename variables to be more commpact
    for var in rename_dict:
        # double check the variable is in the source file
        if var in src:
            # if so, rename the variable to something more compact
            src = src.rename({var : rename_dict[var]})


    # filter valid ice thickness, add small amount to h_min to deal with roundoff
    ds["height"] = xr.where(ds.height <= h_min + 0.1, 0, ds.height)
    # calculate velocity magnitude from velocity components
    ds['vel_m']  = calc_magnitude(ds['vel_x'], ds['vel_z'])

    # calcute the percent temperate
    ds['percent_temperate'] = calc_percent_temperate(ds)
    # calcute the glacier volume as a function of time
    ds['relative_volume']   = calc_volume(ds)

    return ds

# public functions
def dataset(filename, **kwargs):
    ds = xr.open_dataset(filename, **kwargs)
    ds = __preprocess(ds)
    return ds
