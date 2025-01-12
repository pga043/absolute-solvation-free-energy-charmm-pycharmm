#!/bin/bash

export CHARMM_LIB_DIR=/net/orinoco/pga043/charmm/charmm-latest/build_charmm/

END=15
for((i=1;i<END;i++))
do
export i=$i
export lig=lig1

python blade_dyna.py $i $lig > output_"$i" 2> error_"$i"
#python blade_ener.py $i $lig

python domdec_ener.py $i $lig &

done


