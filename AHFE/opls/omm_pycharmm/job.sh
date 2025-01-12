#!/bin/bash

#export CHARMM_LIB_DIR=/net/orinoco/pga043/charmm/charmm-latest/build_4090/
export CHARMM_LIB_DIR=/net/orinoco/pga043/charmm/charmm-latest/build_charmm/

#END=7
#for((i=0;i<=END;i++))
for i in 13 14 
do
export i=$i
rm -rf win"$i"

python omm_dyna.py $i > output_"$i" 2> error_"$i"
python omm_ener.py $i

done


