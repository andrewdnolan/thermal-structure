#!/bin/bash

# load the souce functions
source initialize.sh

# parse the parameters from the json files
parse_json "params/twds-b.json"


dbdt=-0.76428249
T_ma_ref=-9.0


for offset in $(seq -w 0.0 1.0 9.0); do
  ./initialize.py -dx 200 --key "twds-b" -t_f 1000 -off $offset -T_ma $T_ma_ref
done


for offset in $(seq -w 1.0 1.0 9.0); do

  T_ma=$(awk -v T_ma_ref=$T_ma_ref -v dbdt=$dbdt -v offset=$offset \
                'BEGIN {print T_ma_ref + (1.0/dbdt)*offset}')

  ./initialize.py -dx 200 --key "twds-b" -t_f 1000 -off $offset -T_ma $T_ma
done
