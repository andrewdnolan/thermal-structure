#!/bin/bash
#SBATCH --array=1-442%20                           # 442 jobs that run
#SBATCH --job-name=glc1-a_coupled_init             # base job name for the array
#SBATCH --mem-per-cpu=1500M                        # maximum 2250MMB per job
#SBATCH --time=9:00:00                             # maximum walltime per job
#SBATCH --nodes=1                                  # Only one node is needed
#SBATCH --ntasks=1                                 # These are serial jobs
#SBATCH --mail-type=ALL                            # send all mail (way to much)
#SBATCH --mail-user=andrew.d.nolan@maine.edu       # email to spend updates too
#SBATCH --output=logs/glc1-a_coupled_init_%A_%a.out  # standard output
#SBATCH --error=logs/glc1-a_coupled_init_%A_%a.err   # standard error
# in the previous two lines %A" is replaced by jobID and "%a" with the array index

#Load cedar module file
source ../../config/modulefile.cc.cedar

parse_json()
{
  #-----------------------------------------------------------------------------
  # Parse model parameters from input .json file needed for initialization
  # perl code from: https://stackoverflow.com/a/27127626/10221482
  #-----------------------------------------------------------------------------
  # $1 (json_fp) --->   file path to json parameter file
  #-----------------------------------------------------------------------------

  # parse MB grid search start
  MB_0=$( perl -MJSON -0lnE '
    $params = decode_json $_;
    say $params->{MB}->{range}->[0]
    ' $1 )
  # parse MB grid search stride
  MB_s=$( perl -MJSON -0lnE '
    $params = decode_json $_;
    say $params->{MB}->{range}->[1]
    ' $1 )
  # parse MB grid search end
  MB_f=$( perl -MJSON -0lnE '
    $params = decode_json $_;
    say $params->{MB}->{range}->[2]
    ' $1 )
  # parse horizontal gird spacing
  dx=$(   perl -MJSON -0lnE '
    $params = decode_json $_;
    say $params->{dx}
    ' $1 )
  # parse MB curve fit type
  FIT=$(   perl -MJSON -0lnE '
    $params = decode_json $_;
    say $params->{fit}
    ' $1 )
  # parse MB curve fit type
  k=$(   perl -MJSON -0lnE '
    $params = decode_json $_;
    say $params->{k}
    ' $1 )
}

diagnostic_run()
{
  # Update the .SIF FILE with the model run specifc params
  sed "s#<DX>#"$dx"#g;
       s#<FIT>#"$FIT"#g;
       s#<KEY>#"$KEY"#g;
       s#<T_mean>#"$T_ma"#g;
       s#<offset>#"$offset"#g;
       s#<run_name>#"$run_name"#g;
       s#<SS_itters>#"$SS_itters"#g;" "./sifs/diagnostic.sif" > "./sifs/${run_name}.sif"

  # filepath to log file
  # log_file="${KEY}/logs/${run_name}.log"

  # Run the model
  ElmerSolver "./sifs/${run_name}.sif" #| tee $log_file

  # number of steady state itteration required
  n_itter=$(tac "./${KEY}/mesh_dx${dx}/${run_name}.result" | \
             grep -m1 "Time" | tr -s ' ' | cut -d " " -f 2)

  # Convert result files into NetCDFs
  ../../src/elmer2nc/elmer2nc.sh -r "./${KEY}/mesh_dx${dx}/${run_name}.result" \
                                -m "./${KEY}/mesh_dx${dx}/" \
                                -t $n_itter                 \
                                -o "./${KEY}/nc/"

  # # Remove the sif file
  rm "./sifs/${run_name}.sif"

}

prognostic_run()
{
  # Update the .SIF FILE with the model run specifc params
  sed "s#<DX>#"$dx"#g;
       s#<dt>#"$dt"#g;
       s#<NT>#"$NT"#g;
       s#<KEY>#"$KEY"#g;
       s#<FIT>#"$FIT"#g;
       s#<T_mean>#"$T_ma"#g;
       s#<offset>#"$offset"#g;
       s#<RESTART>#"$RESTART"#g
       s#<run_name>#"$run_name"#g;
       s#<SS_itters>#"$SS_itters"#g;" "./sifs/prognostic.sif" > "./sifs/${run_name}.sif"

  # filepath to log file
  log_file="${KEY}/logs/${run_name}.log"

  # Run the model
  ElmerSolver "./sifs/${run_name}.sif" #| tee $log_file

  # # Convert result files into NetCDFs
  # ../../src/elmer2nc/elmer2nc.sh -r "./${KEY}/mesh_dx${dx}/${run_name}.result" \
  #                               -m "./${KEY}/mesh_dx${dx}/" \
  #                               -t $NT              \
  #                               -o "./${KEY}/nc/"

  # # Remove the sif file
  rm "./sifs/${run_name}.sif"

}

log_runtime()
{

  local OUT_fp="${KEY}/${KEY}.coupled_spinup.time_profile"

  if [ ! -f "$OUT_fp" ]; then
      echo "#dx dt NT offset runtime" |
      awk -v OFS='\t' '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5}' >> \
                      $OUT_fp
  fi

  echo "${dx} ${dt} ${NT} ${offset} ${runtime}" |
  awk -v OFS='\t'  '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5}' >> \
                  $OUT_fp
}

# parse the parameters from the json files
parse_json "params/glc1-a.json"

# # make a list of mass balance offsets
# offsets=($(seq -w $MB_0 $MB_s $MB_f))
# # bash arrays zero indexed
# N=$((SLURM_ARRAY_TASK_ID-1))
# # get ith offset corresponing to SLURM_ARRAY_TASK_ID
# offset=${offsets[$N]}

# get ith offset and T_ma corresponing to SLURM_ARRAY_TASK_ID, nested for loops
# each script is junky. The cartesian product of the two arrays would be a more
# elegant solution, but I was having trouble gettig that to work.
# Refs:
#   - https://unix.stackexchange.com/questions/97814/array-cartesian-product-in-bash
#   - https://stackoverflow.com/questions/23363003/how-to-produce-cartesian-product-in-bash
#   - https://rosettacode.org/wiki/Cartesian_product_of_two_or_more_lists#UNIX_Shell

count=1
# loop over the mass balance offsets
for offset in $(seq -w $MB_0 $MB_s $MB_f); do
  # loop over mean annual air temps
  for T_ma in $(seq -w -9.00 0.1 -7.00); do
    # check if counter equals SLURM_ARRAY_TASK_ID
    if [[ $count -eq $SLURM_ARRAY_TASK_ID ]]; then
        break 2
    fi
    count=$((count+1))
  done
done


# set the glacier key
KEY='glc1-a'
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# steady-state run
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# let's allow for X number of steady state itterations
SS_itters=25

# make the run name based on model params
run_name="${KEY}_dx_${dx}_MB_${offset}_OFF_Tma_${T_ma}_diag"

# run the model for a given offset
diagnostic_run $dx $KEY $offset $run_name $SS_itters


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# transient run
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
dt=1.0

NT=2000
TT=2000
# limit to 10 S.S. itters for transient runs
SS_itters=10
# diagnostic run is now restart variable
RESTART="${run_name}.result"
# prognostic run name
run_name="${KEY}_dx_${dx}_NT_${NT}_dt_${dt}_MB_${offset}_OFF_Tma_${T_ma}_prog"

# Start the timer
start=$(date +%s.%N)

# run the transient model with diagnostic solution as restart fiedl
prognostic_run $dx $KEY $offset $run_name $SS_itters $restart $NT $dt

# End the timer
end=$(date +%s.%N)

# Execution time of the solver
runtime=$(awk -v start=$start -v end=$end 'BEGIN {print end - start}')

# record the prognostic runtimes for future info
log_runtime
