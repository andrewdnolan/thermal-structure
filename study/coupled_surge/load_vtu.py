from glob import glob
from paraview.simple import *


dir_fp = '/Users/andrewnolan/Thesis/thermal-structure/study/coupled_surge/'
pseudo_fps = 'glc1-a/VTU/glc1-a_dx_50_NT_60_dt_0.5_MB_-1.225_OFF_Tma_-9.00_B_0.003_pseudo_t*.vtu'
recovery_fps = 'glc1-a/VTU/glc1-a_dx_50_NT_100_dt_1.0_MB_-1.225_OFF_Tma_-9.00_B_0.003_recovery_t*.vtu'

test = OpenDataFile(sorted(glob(dir_fp+pseudo_fps)))
test = OpenDataFile(sorted(glob(dir_fp+recovery_fps)))
