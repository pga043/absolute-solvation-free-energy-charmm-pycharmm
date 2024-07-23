#!/bin/bash

## LigParGen run on one of the PC-lab systems
## pga043@129.177.190.61 
## conda activate boss

export BOSSdir=/export/softwares/BOSS/boss/

while read i
do

echo "$i"

#obabel -imol2 ../"$i".mol2 -opdb > "$i".pdb

## cgen = CM1A-LBCC working with BOSS5.0

#ligpargen -i "$i".pdb -cgen CM1A-LBCC -c 0 -cgenb CM1A-LBCC -cb 0

#rm *.q.* *gmx* *tinker* *openmm* *desmond* *xplor* *lammps*

python replaceAtomTypes.py "$i"

cp "$i".tmp.rtf ../"$i".rtf
cp "$i".tmp.prm ../"$i".prm

done < <(cat ../mol_list.txt)

## ligpargen: http://zarbi.chem.yale.edu/ligpargen/index.html

## charge model:  1.14*CM1A-LBCC


