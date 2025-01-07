#!/bin/bash

export charmm=/net/orinoco/pga043/charmm/49a1_blade/charmm/build_charmm/charmm
export cgenff=/net/orinoco/apps/silcsbio.2024.1/cgenff/cgenff

export dir=.
export lig=lig2
export seg=MOL

obabel -imol2 ../../../rdkit/"$lig".mol2 -opdb > "$lig".pdb
$cgenff -v ../../../rdkit/"$lig".mol2 > "$lig".str
./regroup.awk "$lig".str > "$lig"_rg.str
sed -i "s/UNL/MOL/g" "$lig".pdb
sed -i "s/UNK/MOL/g" "$lig"_rg.str

$charmm dir=$dir lig=$lig seg=$seg -i build_ligand.inp > $dir/build_ligand.out
$charmm dir=$dir lig=$lig -i  solvate.inp > $dir/solvate.out

export box=`grep 'GREATERVALUE' $dir/solvate.out | head -n 1 | awk '{print $4}'i | sed 's/^"\(.*\)"$/\1/'`
export fft=`python fft.py $box`

echo "$box" "$fft"

