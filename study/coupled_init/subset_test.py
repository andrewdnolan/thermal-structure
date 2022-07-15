import xarray as xr

def preprocess(ds):
    return ds.isel(t=[0,-1])

src_fp = 'glc1-a/nc/glc1-a_dx_50_NT_2000_dt_1.0_MB_*_OFF_Tma_*_prog_gridded.nc'
gs_ds  = xr.open_mfdataset(src_fp, parallel=False, preprocess=preprocess)

out_fn = 'glc1-a_dx_50_SS@_2000_dt_1.0_MB_-1.5__-1.0_OFF_Tma_-9.0__-7.0_prog_gridded.nc'

gs_ds.to_netcdf(out_fn, 'w')
