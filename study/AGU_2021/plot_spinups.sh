plot_gridsearch() {
  ##############################################################################
  # Make volume and final Z_s plots from the gridsearch results
  #
  # Input Parameters:
  # -----------------
  # dx      --->   mesh resolution                         [m]
  # dt      --->   time step                               [y]
  # TT      --->   length of the simulation                [y]
  # MB_0    --->   start  of MB gridsearch                 [m i.e.q yr^-1]
  # MB_f    --->   end    of MB gridsearch                 [m i.e.q yr^-1]
  # MB_s    --->   stride of MB gridsearch                 [m i.e.q yr^-1]
  # KEY     --->   glacier id key
  # FIT     --->   fit type to Young et al 2020. MB data
  ##############################################################################


  #---------------------------------------------------------------------------
  # Make the volume plots
  #---------------------------------------------------------------------------
  python3 ../../src/plotting/plot_spinup.py \
           -fp "${FP}/result/${KEY}/nc/${KEY}_${TT}a_dt_${dt}_dx_${dx}_MB_*_OFF.nc" \
           -mb $MB_0 $MB_s $MB_f \
           --plot_volume      \
           --title "Little Kluane $ (dx=${dx}) $" \
           -out_fn "./figs/${KEY}_Vol_${MB_0}--${MB_f}_dx_${dx}m.png"

   #--------------------------------------------------------------------------
   # Make the final z_s plot
   #--------------------------------------------------------------------------
   python3 ../../src/plotting/plot_spinup.py \
           -fp "${FP}/result/${KEY}/nc/${KEY}_${TT}a_dt_${dt}_dx_${dx}_MB_*_OFF.nc" \
           -mb $MB_0 $MB_s $MB_f \
           --plot_Z_s         \
           --title "Little Kluane $ (dx=${dx}) $" \
           -out_fn "./figs/${KEY}_Zs_${MB_0}--${MB_f}_dx_${dx}m.png"
  }

#-------------------------------------------------------------------------------
# grid search parameters
#-------------------------------------------------------------------------------
dx=50                                  # mesh resolution
dt=1                                    # time step size
TT=1000                                 # number of time step
MB_0=0.0                                # MB offset start
MB_f=1.0                                # MB offset finish
MB_s=0.1                                # MB offset stride
SIF='./sifs/simple_spinup.sif'          # template SIF file
KEY='lilk-a'                            # glacier key for input data
FIT='cubic_spline'                      # type of fit to Young et al. 2020 data
PLOT=1                                  # plot gridsearch results (0 false, 1 true)
PROFILE=1                               # profile individual run  (0 false, 1 true)
FP="../MB_tune/"

# # Run the mass balance gridsearch
# gridsearch $dx $dt $TT $MB_0 $MB_f $MB_s $SIF $KEY $FIT
#
# plot the gridsearch data
if [ $PLOT -eq 1 ]; then
  plot_gridsearch dx $dt $TT $MB_0 $MB_f $MB_s $KEY $FIT $fp
fi
