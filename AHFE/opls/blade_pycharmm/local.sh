#!/bin/bash

export CHARMM_LIB_DIR=/net/orinoco/pga043/charmm/charmm-latest/build_charmm/

END=15
for((i=0;i<END;i++))
do
export i=$i

#python blade_dyna.py $i > output_"$i" 2> error_"$i"
#python blade_ener.py $i

python domdec_ener.py $i &

done


