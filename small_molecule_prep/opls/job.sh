#!/bin/bash

export charmm=/net/orinoco/pga043/charmm/49a1_blade/charmm/build_charmm/charmm

export dir=.
export lig=lig7
export seg=MOL

$charmm dir=$dir lig=$lig seg=$seg -i build_ligand.inp > $dir/build_ligand.out
$charmm dir=$dir lig=$lig -i  solvate.inp > $dir/solvate.out

export box=`grep 'GREATERVALUE' $dir/solvate.out | head -n 1 | awk '{print $4}'i | sed 's/^"\(.*\)"$/\1/'`
export fft=`python fft.py $box`

echo "$box" "$fft"

