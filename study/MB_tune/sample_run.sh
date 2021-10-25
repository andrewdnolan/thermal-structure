#!/usr/bin/env bash

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# sample_run.sh:
#   - Mass balance grid-search to find steady state positions for 3
#     test glaciers
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

set +x

dt=1                                    # time step size
NT=1000                                 # number of time step
TT=$((NT*dt))                           # total time of simulation

SIF='./simple_spinup.sif'               # template SIF file
BED='farinotti_corrected'               # Mesh DB for the given bed config
KEY='glc1-a'

for dx in 50 100; do
  for OFFSET in $(seq -w -0.1 -0.1 -2.0);do
    # Model RUN identifier
    RUN="${KEY}_${TT}a_dt_${dt}_dx_${dx}_MB_${OFFSET}_OFF"

    # File paths to input data
    Zb_fp="../../input_data/${KEY}_bed.dat"
    Zs_fp="../../input_data/${KEY}_surf.dat"

    # Update the .SIF FILE with the model run specifc params
    sed "s#<NT>#"$NT"#g;
         s#<DX>#"$dx"#g;
         s#<RUN>#"$RUN"#g;
         s#<KEY>#"$KEY"#g;
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

   echo "${dx} ${OFFSET} ${runtime}" |
   awk -v OFS='\t' '{print $1 "\t" $2 "\t" $3}' >> out.time_profile


   # Convert result files into NetCDFs
   ../../src/elmer2nc/elmer2nc.sh  -r "./result/${KEY}/mesh_dx${dx}/${RUN}.result" \
                                   -m "./result/${KEY}/mesh_dx${dx}/" \
                                   -t $NT \
                                   -o "./result/${KEY}/nc/"

    # Remove the sif file
    rm "./sifs/${RUN}.sif"
  done

  #-----------------------------------------------------------------------------
  # Make the volume plots
  #-----------------------------------------------------------------------------
  python3 ../../src/plotting/plot_spinup.py \
           -fp "./result/${KEY}/nc/${KEY}_${TT}a_dt_${dt}_dx_${dx}_MB_*_OFF.nc" \
           -mb -0.1 -0.1 -2.0 \
           --plot_volume      \
           --title "$ Dx=${dx} $" \
           -out_fn "./figs/${KEY}/Vol_-2.0--0.0_dx_${dx}m.png"

   #-----------------------------------------------------------------------------
   # Make the final z_s plot
   #-----------------------------------------------------------------------------
   python3 ../../src/plotting/plot_spinup.py \
           -fp "./result/${KEY}/nc/${KEY}_${TT}a_dt_${dt}_dx_${dx}_MB_*_OFF.nc" \
           -mb -0.1 -0.1 -2.0 \
           --plot_Z_s         \
           --title "$ Dx=${dx} $" \
           -out_fn "./figs/${KEY}/Zs_-2.0--0.0_dx_${dx}m.png"
done
