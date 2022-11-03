#!/usr/bin/env bash


offset=-0.41
T_ma=-8.5

for key in 'crmpt12_2500_clip' 'crmpt12_2600_clip'; do

  ./initialize.py -k $key -dx 50 -off $offset -T_ma $T_ma -dt 1.0 -t_f 2000

done
