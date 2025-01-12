#!/bin/bash

## LigParGen run on one of the PC-lab systems
## pga043@129.177.190.61 
## conda activate boss

export BOSSdir=/export/softwares/BOSS/boss/

for i in lig7
do

echo "$i"
obabel -imol2 ../../../rdkit/"$i".mol2 -opdb > "$i".pdb
sed -i "s/UNL/MOL/g" "$i".pdb
## cgen = CM1A-LBCC working with BOSS5.0

ligpargen -i "$i".pdb -cgen CM1A-LBCC -c 0 -cgenb CM1A-LBCC -cb 0

rm *.q.* *gmx* *tinker* *openmm* *desmond* *xplor* *lammps* *.log *.pqr

#python replaceAtomTypes.py "$i"

#cp "$i".tmp.rtf ../"$i".rtf
#cp "$i".tmp.prm ../"$i".prm

python opls2charmm.py "$i"

mv "$i"_charmm.rtf ..
mv "$i"_charmm.prm ..
mv "$i".charmm.pdb ..

done 

## ligpargen: http://zarbi.chem.yale.edu/ligpargen/index.html

## charge model:  1.14*CM1A-LBCC


