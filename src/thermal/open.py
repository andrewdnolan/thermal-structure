
import numpy as np
import xarray as xr

def _quads_to_tris(quads):
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

def _preprocess_UGRID(ds):
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

    # create a new dataset to populate
    new_ds = xr.Dataset()

    # add time dimensions with values from UGRID source
    new_ds = new_ds.expand_dims({"t": NT}).assign_coords({"t": ds.time.values})

    # Grid the nodal X values and add as coordinate
    new_ds.coords["X"] = xr.DataArray(
        ds.Mesh_node_x.values.reshape(NZ, NX), dims=["coord_2", "coord_1"]
    )

    # check if Z coordinate is time dependent
    if "time" in ds.Mesh_node_y.dims:
        # Grid the time dependent nodal Z values and add as coordinate
        new_ds.coords["Z"] = xr.DataArray(
            ds.Mesh_node_y.values.reshape(NT, NZ, NX), dims = ["t", "coord_2", "coord_1"]
        )
    else:
        # Grid the nodal Z values and add as coordinate
        new_ds.coords["Z"] = xr.DataArray(
            ds.Mesh_node_y.values.reshape(NZ, NX), dims = dims
        )

    # Grid the node number and add as a coordinate
    new_ds.coords["NN"] = xr.DataArray(
        ds.nMesh_node.values.reshape(NZ, NX), dims=["coord_2", "coord_1"]
    )
    # extract nodes of quad elements from UGRID source
    new_ds["quad_elements"] = xr.DataArray(
        ds.Mesh_face_nodes.values - 1, dims=["quad_element_number", "quad_element_node"]
    )
    # extract the element area from UGRID source
    new_ds["quad_elements_area"] = xr.DataArray(
        ds.BulkElement_Area.values, dims=["quad_element_number"]
    )
    # split the quad elements into triangles for plotting with matplotlib
    new_ds["tri_elements"] = xr.DataArray(
        _quads_to_tris(new_ds.quad_elements.values),
        dims=["tri_element_number", "tri_element_node"],
    )

    # loop over the elmer vairables
    for key in ds.keys():

        # these fields have already been processed as coordinates.
        if key in ("Mesh_node_x", "Mesh_node_y", "Mesh_face_nodes"):
            continue

        var  = ds[key]   # get the coresponding variable
        dims = var.dims  # and it's dimensions

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

def _preprocess_elmer2nc(ds):
    pass

def _preprocess(ds):
    """Wrapper around various preprocessing functions.
    """
    if 'nMesh_node' in ds.dims:
        ds =  _preprocess_UGRID(ds)
    else:
        print('no other support implemented yet. To be Done')

    ds["height"] = xr.where(ds.height <= 10, 0, ds.height)
    ds['vel_m']  = np.sqrt(ds['velocity 1']**2 + ds['velocity 2']**2)

    return ds

def dataset(filename, **kwargs):
    ds = xr.open_dataset(filename, **kwargs)
    ds = _preprocess(ds)
    return ds

def parameter_ensemble():
    return None

# def mf_dataset(filename, **kwargs):
#     # if
#     return None
