#!/usr/bin/env bash

dx=50         # horizontal gridcell spacing [m]
SP=2          # surge period
QP=2000       # recovery (i.e. quiescent) period
S_dt=0.05     # surging timestep  [a]
Q_dt=1.0      # recovery timestep [a]
key="crmpt12" # glacier identifier
T_ma=0.0      # Dummy value for the air temperature
offset=-0.444    # mass balance anomoly [m a-1]

./surge2steady.py -k $key -dx $dx -SP $SP -QP 100 -beta 1e-4 -SD_dt $S_dt -ST_dt $S_dt \
                    -QT_dt $Q_dt -QD_dt $Q_dt -z_lim 2400 -off $offset -T_ma $T_ma \
                    -RESTART "crmpt12_dx_50_NT_3000_dt_1.0_MB_-0.444_OFF_ISOTHERM_prog.result"

# #corresponds to np.logspace(-3,-4, 9)
# betas=(1.000e-03 7.499e-04 5.623e-04 4.217e-04 3.162e-04 
#        2.371e-04 1.778e-04 1.334e-04 1.000e-04)

# for beta in ${betas[@]}; do
#     ./surge2steady.py -k $key -dx $dx -SP $SP -QP $QP -beta $beta -SD_dt $S_dt -ST_dt $S_dt \
#                       -QT_dt $Q_dt -QD_dt $Q_dt -z_lim 2400 -offset -0.444 -T_ma $T_ma
#                       -RESTART "crmpt12_dx_50_NT_3000_dt_1.0_MB_-0.444_OFF_ISOTHERM_prog.result"
# done