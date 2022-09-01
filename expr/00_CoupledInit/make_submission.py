
"""
# parse the parameters from the json files
parse_json "params/glc1-a.json"

# get ith offset and T_ma corresponing to SLURM_ARRAY_TASK_ID.
# Nested for loops each job in the job-array is junky.
# The cartesian product of the two arrays would be a more elegant solution,
# but I was having trouble gettig that to work.
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
"""
