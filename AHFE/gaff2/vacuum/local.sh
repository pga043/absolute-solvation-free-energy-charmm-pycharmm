#!/bin/bash

#export CHARMM_LIB_DIR=/net/orinoco/pga043/charmm/charmm-latest/build_4090/

#END=8
#for((i=0;i<END;i++))
for i in 1 4 5
do
export i=$i
export lig=lig1

echo win"$i" started
python charmm_dyna.py $i $lig > output_"$i" 2> error_"$i"
python charmm_ener.py $i $lig 

done


