#!/bin/bash

#export CHARMM_LIB_DIR=/net/orinoco/pga043/charmm/charmm-latest/build_4090/

END=7
for((i=0;i<=END;i++))
do
export i=$i
export lig=lig1

python charmm_dyna.py $i $lig > output_"$i" 2> error_"$i"
python charmm_ener.py $i $lig &

done


