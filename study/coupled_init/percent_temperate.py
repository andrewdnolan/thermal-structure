#!/usr/bin/env python3

import sys
import numpy as np
import xarray as xr
from glob import glob
from scipy import interpolate
import matplotlib.pyplot as plt


src_fp = '/Volumes/thermal_ssd/DataDir/glc1-a_dx_50_NT_2000_dt_1.0_MB_-1.075_OFF_Tma_-7.00_prog.nc'

src = xr.open_dataset(src_fp)

subset  = src.isel(time=slice(0,None))
subset['Mesh_face_nodes'] = subset.Mesh_face_nodes - 1

# hard code a solution to compare too
element_i_idx  = subset.isel(time=-1).Mesh_face_nodes[0]
element_i_enth = subset.isel(time=-1).enthalpy_h[element_i_idx.values.astype(int)].mean('nMesh_node')

def perc_temp(ds):
    element_enth = xr.where(ds.nMesh_node == ds.Mesh_face_nodes, ds.enthalpy_h, np.nan).mean(dim=('nMesh_node','nMaxMesh_face_nodes'))
    element_H_f  = xr.where(ds.nMesh_node == ds.Mesh_face_nodes, ds['phase change enthalpy'], np.nan).mean(dim=('nMesh_node','nMaxMesh_face_nodes'))

    total_area = ds.BulkElement_Area.sum('nMesh_face')
    temp_area  = xr.where(element_enth>=element_H_f, ds.BulkElement_Area, 0.0).sum('nMesh_face')

    perc_temp = (temp_area/total_area) * 100

    return perc_temp

test = subset.groupby("time").map(perc_temp)

print(test)

# xr.where( subset.nMesh_node == subset.Mesh_face_nodes, 1, 0)

# # Node indexes do not change as a function time or ensemble parameter value
# drop_dims = {dim : 0  for dim in SS_2kya.quad_elements.dims
#                 if dim not in ['quad_element_number', 'quad_element_node']}
# # Select one value per auxiliry dimensions and drop vars for array
# quad_elements = SS_2kya.quad_elements.isel(drop_dims).drop_vars(drop_dims.keys())
#
# stacked = subset.stack(gridcell=('coord_1','coord_2'))
#
# test = xr.where(stacked.NN == quad_elements, stacked.enthalpy_h, 0.0).mean('quad_element_node')
#
# # test = subset.enthalpy_h.where(subset.NN == quad_elements).mean('quad_element_node')
#
# # stacked = subset.stack(gridcell=('coord_1','coord_2'))
# # print(subset)
