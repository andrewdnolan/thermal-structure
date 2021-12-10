#!/usr/bin/env bash

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# run.sh:
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

set +x

#-------------------------------------------------------------------------------
# Numerical parameters
#-------------------------------------------------------------------------------
dx=500                                  # mesh resolution
dt=0.1                                  # time step size
NT=251                                  # number of time step
TT=$(awk -v NT=$NT -v dt=$dt "BEGIN { print NT*dt - dt }" )

#-------------------------------------------------------------------------------
# input data parameters
#-------------------------------------------------------------------------------
KEY='twds-b'                            # glacier key for input data
BETA=0.0001                             # Slip coef
OFFSET=10.8                             # MB offset start
RESTART="twds-b_1000a_dt_0.5_dx_500_MB_10.8_OFF_cubic_spline.result"

# File paths to input data
Zb_fp="../../input_data/${KEY}_bed.dat"
Zs_fp="../../input_data/${KEY}_surf.dat"
#-------------------------------------------------------------------------------
# Run the surge
#-------------------------------------------------------------------------------
RUN="${KEY}_${TT}a_dt_${dt}_dx_${dx}_MB_${OFFSET}_OFF_${BETA}_B_pseudo"

# File paths to input data
Zb_fp="../../input_data/${KEY}_bed.dat"
Zs_fp="../../input_data/${KEY}_surf.dat"

# Update the .SIF FILE with the model run specifc params
sed "s#<dt>#"$dt"#g;
     s#<NT>#"$NT"#g;
     s#<DX>#"$dx"#g;
     s#<RUN>#"$RUN"#g;
     s#<KEY>#"$KEY"#g;
     s#<Zs_fp>#"$Zs_fp"#g;
     s#<Zb_fp>#"$Zb_fp"#g;
     s#<OFFSET>#"$OFFSET"#g;
     s#<RESTART>#"$RESTART"#g" "./sifs/pseudo_surge.sif" > "./sifs/${RUN}.sif"

# Start the timer
start=$(date +%s.%N)

# Run the model
ElmerSolver "./sifs/${RUN}.sif"

# End the timer
end=$(date +%s.%N)

# Execution time of the solver
runtime=$(awk -v start=$start -v end=$end 'BEGIN {print end - start}')

echo "${dx} ${OFFSET} ${runtime}" |
awk -v OFS='\t' '{print $1 "\t" $2 "\t" $3}' >> "${KEY}.pseudo.time_profile"


# Convert result files into NetCDFs
../../src/elmer2nc/elmer2nc.sh -r "./mesh_dx${dx}/${RUN}.result" \
                               -m "./mesh_dx${dx}/" \
                               -t $NT \
                               -o "./nc/"

# Remove the sif file
rm "./sifs/${RUN}.sif"


#-------------------------------------------------------------------------------
# Run the recovery from the surge
#-------------------------------------------------------------------------------

STT=$TT                                 # length of surge simulation
dt=1                                    # time step size
NT=1000                                 # number of time step
TT=$((NT*dt))                           # total time of simulation
RESTART="${RUN}.result"

RUN="${KEY}_${TT}a_dt_${dt}_dx_${dx}_MB_${OFFSET}_OFF_${STT}a_${BETA}_B_recovery"

# Update the .SIF FILE with the model run specifc params
sed "s#<dt>#"$dt"#g;
     s#<NT>#"$NT"#g;
     s#<DX>#"$dx"#g;
     s#<RUN>#"$RUN"#g;
     s#<KEY>#"$KEY"#g;
     s#<Zs_fp>#"$Zs_fp"#g;
     s#<Zb_fp>#"$Zb_fp"#g;
     s#<OFFSET>#"$OFFSET"#g;
     s#<RESTART>#"$RESTART"#g" "./sifs/recovery.sif" > "./sifs/${RUN}.sif"

# Start the timer
start=$(date +%s.%N)

# Run the model
ElmerSolver "./sifs/${RUN}.sif"

# End the timer
end=$(date +%s.%N)

# Execution time of the solver
runtime=$(awk -v start=$start -v end=$end 'BEGIN {print end - start}')

echo "${dx} ${OFFSET} ${runtime}" |
awk -v OFS='\t' '{print $1 "\t" $2 "\t" $3}' >> "${KEY}.recovery.time_profile"


# Convert result files into NetCDFs
../../src/elmer2nc/elmer2nc.sh -r "./mesh_dx${dx}/${RUN}.result" \
                               -m "./mesh_dx${dx}/" \
                               -t $NT \
                               -o "./nc/"

# Remove the sif file
rm "./sifs/${RUN}.sif"
