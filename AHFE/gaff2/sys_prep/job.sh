#!/bin/bash

export charmm=/net/orinoco/pga043/charmm/49a1_blade/charmm/build_charmm/charmm

export dir=.
export lig=lig1
export seg=MOL

#obabel -imol2 ../../../rdkit/"$lig".mol2 -opdb > "$lig".pdb
#cp ../../../rdkit/"$lig".sdf .
#cp ../../../rdkit/"$lig".mol2 .

#echo "Now processing mol" $lig

#python ligand2charmm.py -Lname "$lig" -nc 0 -index 1
#rm ANTECHAMBER* sqm.* ATOMTYPE.INF *.ac lig*.inp lig*.old *.frcmod *.original lig*.crd

#python gaff2charmm.py "$lig" 
#sed -i "s/UNL/MOL/g" "$lig".pdb
#sed -i "s/IMPH/IMPR/g" "$lig"_charmm.rtf

$charmm dir=$dir lig=$lig seg=$seg -i build_ligand.inp > $dir/build_ligand.out
$charmm dir=$dir lig=$lig -i  solvate.inp > $dir/solvate.out

#export box=`grep 'GREATERVALUE' $dir/solvate.out | head -n 1 | awk '{print $4}'i | sed 's/^"\(.*\)"$/\1/'`
#export fft=`python fft.py $box`

#echo "$box" "$fft"

