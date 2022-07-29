
# `elmer2nc`

This is program parses `*.result` files written by `Elmer`, and re-writes them `NetCDF` files.

Basic information about the progam can be found by running:
```bash
./elmer2nc.sh --help
```
which returns:
```shell
Usage ./elmer2nc.sh -r <result_fp> [further options...]
     Convert Elmer's .result file to a NetCDF using Fortran code found in the src dir.
     By deafault outfile file will have same name as input but with `.nc` file extension

     [-r <fp>]  => file path to .result file we  want to converts
     [-m <dir>] => mesh directory. If not supplied infered from <result_fp>
     [-t <NT>]  => Number of timesteps within <result_fp> {integer}
     [-o <dir>] => output directory to write .nc file. If not supplied will write to -m <dir>
```

### Notes:
 - currently only works for structured 2-D grids with quadrilateral (404) nodes  
 - currently only works for `ASCII` output formats, no support for binary output yet  
 - The resulting `NetCDF` files have an unstructured $z$ coordinate in. Too apply the `$\zeta$` coordinate transform, with some other basic post-processing, in python use:
  ```python
  import xarray as xr

  with xr.open_dataset(src_fp) as src:
      # correct for minimum ice thickness
      src["height"] = xr.where(src.height <= 10, 0, src.height)
      # apply sigma coordinate transform for vertical coordinate
      src["Z"]     = src.zbed + src.Z * src.height
      # Calculate the magnitude of the velocity vectors
      src['vel_m'] = np.sqrt(src['velocity 1']**2 + src['velocity 2']**2)    
```
 - This program is no longer really needed with the addition of the [`NetcdfUGRIDOutputSolver.f90`](https://github.com/andrewdnolan/thermal-structure/blob/main/src/elmer_src/NetcdfUGRIDOutputSolver.f90), which can write `NetCDF` files during the model exection.  
