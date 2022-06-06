
import numpy as np
import xarray as xr

def preprocess_UGRID(ds):
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
        dims = ["t", "coord_2", "coord_1"]
    else:
        dims = ["coord_2", "coord_1"]

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
        quads_to_tris(new_ds.quad_elements.values),
        dims=["tri_element_number", "tri_element_node"],
    )

    # loop over the elmer vairables
    for key in ds.keys():

        # these fields have already been processed as coordinates.
        if key in ("Mesh_node_x", "Mesh_node_y"):
            continue

        # get the coresponding variable and it's dimensions
        var  = ds[key]
        dims = var.dims

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
