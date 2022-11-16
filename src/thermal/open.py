import numpy as np
import xarray as xr
from glob import glob
from tqdm import tqdm
from derived_fields import calc_percent_temperate, calc_relative_volume


"""
add out of memoery method for gridding the Elmer/Ice NetCDF
    Need for the 2GB+ files written which take up to much memoery when gridded in memoery

    https://github.com/pydata/xarray/issues/1215

# open multiple zarr
    https://discourse.pangeo.io/t/how-to-read-multiple-zarr-archives-at-once-from-s3/2564
"""


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

# def __create_gridded():

def __dump_var(new_ds, old_ds, key):
    """ add gridded variable to the new dataset
    """
    NT = new_ds.t.size         # number of timesteps
    NX = new_ds.coord_1.size   # number nodes in the x direction
    NZ = new_ds.coord_2.size   # number nodes in the z direction

    var  = old_ds[key] # get the coresponding variable
    dims = var.dims    # and it's dimensions

    if ("nMesh_node" in dims) and ("time" in dims):
        new_ds[key] = xr.DataArray(
            var.values.reshape(NT, NZ, NX), dims=["t", "coord_2", "coord_1"]
        )

    # non timedependent variables?
    elif "nMesh_node" in dims:
        new_ds[key] = xr.DataArray(
            var.values.reshape(NZ, NX), dims=["coord_2", "coord_1"]
        )

    return new_ds


def __preprocess_UGRID(ds, out_fn=None):
    """Reshape the UGRID NetCDF results onto a structured grid.

    Parameters
    ----------
    ds : `xarray.Dataset`
        Parsed Dataset written by the "NetcdfUGRIDOutputSolver.f90" file

    Returns
    -------
    new_ds : `xarray.Dataset`
        Processed Dataset which has been gridded and has common coordinate and
        dimensions as the "result2nc" created NetCDF files

    """
    # Get grid info
    NT = ds.time.size                   # number of timesteps
    NN = ds.nMesh_node.size             # number of nodes
    NX = np.unique(ds.Mesh_node_x).size # number nodes in the x direction
    NZ = NN // NX                       # number nodes in the z direction

    return new_ds

def __preprocess(ds, h_min=10.0):
    """Wrapper around various preprocessing functions.
    """
    if 'nMesh_node' in ds.dims:
        ds = __preprocess_UGRID(ds)
    else:
        print('no other support implemented yet. To be Done')

    # filter valid ice thickness, add small amount to h_min to deal with roundoff
    ds["height"] = xr.where(ds.height <= h_min + 0.1, 0, ds.height)
    # calculate velocity magnitude from velocity components
    ds['vel_m']  = np.sqrt(ds['velocity 1']**2 + ds['velocity 2']**2)
    # rename velocities to be more commpact, and match notation of "vel_m"
    ds = ds.rename({"velocity 1" : "vel_x", "velocity 2" : "vel_z"})
    # calcute the percent temperate
    ds['percent_temperate'] = calc_percent_temperate(ds)
    # calcute the percent temperate
    ds['relative_volume']   = calc_relative_volume(ds)

    return ds

# public functions
def dataset(filename, **kwargs):
    ds = xr.open_dataset(filename, **kwargs)
    ds = __preprocess(ds)
    return ds

def mf_dataset(files, preprocess=None, parallel=False, concat_dim=None, **open_kwargs):

    # check if glob was passed to function
    if (type(files) == str) and ('*' in files):
        paths = sorted(glob(files))
    elif (type(files) == str) and ('*' not in files):
        raise ValueError('If `type(files)==str` then must be a glob (i.e. include a *)')

    if parallel:
        import dask

        # wrap the open_dataset, getattr, and preprocess with delayed
        open_ = dask.delayed(dataset)
        getattr_ = dask.delayed(getattr)
        if preprocess is not None:
            preprocess = dask.delayed(preprocess)
    else:
        open_ = dataset
        getattr_ = getattr


    datasets = [open_(p, **open_kwargs) for p in paths]
    closers  = [getattr_(ds, "_close") for ds in datasets]
    if preprocess is not None:
        datasets = [preprocess(ds) for ds in datasets]

    if parallel:
        # calling compute here will return the datasets/file_objs lists,
        # the underlying datasets will still be stored as dask arrays
        datasets, closers = dask.compute(datasets, closers)


    combined = xr.concat(datasets, dim=concat_dim)

    def multi_file_closer():
        for closer in closers:
            closer()

    combined.set_close(multi_file_closer)


    return combined
