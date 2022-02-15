#!/usr/bin/env bash

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# sample_run.sh:
#   - Mass balance grid-search to find steady state positions for 3
#     test glaciers
#
# To Do:
#   - Add the mass balance offset as a single scalar field, so that we can concat
#     along it later.
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

set +x

gridsearch(){
  ##############################################################################
  # Run a netbalance gridsearch to find a steady state
  #
  # Input Parameters:
  # -----------------
  # dx      --->   mesh resolution                         [m]
  # dt      --->   time step                               [y]
  # TT      --->   length of the simulation                [y]
  # MB_0    --->   start  of MB gridsearch                 [m i.e.q yr^-1]
  # MB_f    --->   end    of MB gridsearch                 [m i.e.q yr^-1]
  # MB_s    --->   stride of MB gridsearch                 [m i.e.q yr^-1]
  # SIF     --->   path to template sif file
  # KEY     --->   glacier id key
  # FIT     --->   fit type to Young et al 2020. MB data
  # PROFILE --->   boolean whether to profile runtime
  ##############################################################################

  NT=$( awk -v TT=$TT -v dt=$dt "BEGIN { print TT/dt }" )

  for OFFSET in $(seq -w $MB_0 $MB_s $MB_f);do
    # Model RUN identifier
    RUN="${KEY}_${TT}a_dt_${dt}_dx_${dx}_MB_${OFFSET}_OFF_${FIT}"

    # File paths to input data
    Zb_fp="../../input_data/${KEY}_bed.dat"
    Zs_fp="../../input_data/${KEY}_surf.dat"

    # Update the .SIF FILE with the model run specifc params
    sed "s#<dt>#"$dt"#g;
         s#<NT>#"$NT"#g;
         s#<DX>#"$dx"#g;
         s#<RUN>#"$RUN"#g;
         s#<KEY>#"$KEY"#g;
         s#<FIT>#"$FIT"#g;
         s#<Zs_fp>#"$Zs_fp"#g;
         s#<Zb_fp>#"$Zb_fp"#g;
         s#<OFFSET>#"$OFFSET"#g" "$SIF" > "./sifs/${RUN}.sif"

   # Start the timer
   start=$(date +%s.%N)

   # Run the model
   ElmerSolver "./sifs/${RUN}.sif"

   # End the timer
   end=$(date +%s.%N)

   # Execution time of the solver
   runtime=$(awk -v start=$start -v end=$end 'BEGIN {print end - start}')

   if [ $PROFILE -eq 1 ]; then
     echo "${dx} ${OFFSET} ${runtime}" |
     awk -v OFS='\t' '{print $1 "\t" $2 "\t" $3}' >> \
                     "result/${KEY}/${KEY}.spinup.time_profile"
   fi

   # Convert result files into NetCDFs
   ../../src/elmer2nc/elmer2nc.sh  -r "./result/${KEY}/mesh_dx${dx}/${RUN}.result" \
                                   -m "./result/${KEY}/mesh_dx${dx}/" \
                                   -t $NT \
                                   -o "./result/${KEY}/nc/"

    # __TO DO__:  Add the mass balance offset variable to the NetCDF file

    # Remove the sif file
    rm "./sifs/${RUN}.sif"
  done

  }


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
           -fp "./result/${KEY}/nc/${KEY}_${TT}a_dt_${dt}_dx_${dx}_MB_*_OFF_${FIT}.nc" \
           -mb $MB_0 $MB_s $MB_f \
           --plot_volume      \
           --title "$ Dx=${dx} $" \
           -out_fn "./figs/${KEY}/Vol_${MB_0}--${MB_f}_dx_${dx}m.png"

   #--------------------------------------------------------------------------
   # Make the final z_s plot
   #--------------------------------------------------------------------------
   python3 ../../src/plotting/plot_spinup.py \
           -fp "./result/${KEY}/nc/${KEY}_${TT}a_dt_${dt}_dx_${dx}_MB_*_OFF_${FIT}.nc" \
           -mb $MB_0 $MB_s $MB_f \
           --plot_Z_s         \
           --title "$ Dx=${dx} $" \
           -out_fn "./figs/${KEY}/Zs_${MB_0}--${MB_f}_dx_${dx}m.png"
  }

#-------------------------------------------------------------------------------
# grid search parameters
#-------------------------------------------------------------------------------
dx=500                                  # mesh resolution
dt=0.5                                  # time step size
TT=1000                                 # number of time step
MB_0=1.0                                # MB offset start
MB_f=1.5                                # MB offset finish
MB_s=0.1                                # MB offset stride
SIF='./sifs/simple_spinup.sif'          # template SIF file
KEY='klut-b'                            # glacier key for input data
FIT='cubic_spline'                      # type of fit to Young et al. 2020 data
PLOT=1                                  # plot gridsearch results (0 false, 1 true)
PROFILE=1                               # profile individual run  (0 false, 1 true)

# # Run the mass balance gridsearch
# gridsearch $dx $dt $TT $MB_0 $MB_f $MB_s $SIF $KEY $FIT
#
# plot the gridsearch data
if [ $PLOT -eq 1 ]; then
  plot_gridsearch dx $dt $TT $MB_0 $MB_f $MB_s $KEY $FIT
fi
