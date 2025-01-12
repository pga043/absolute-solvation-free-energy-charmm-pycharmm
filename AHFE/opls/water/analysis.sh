#!/bin/bash
#Analysis part requires python2 so execute conda activate py2

#CHARMMEXEC=/net/orinoco/pga043/CHARMM_47a2/charmm_47a2/build_openmm_blade/charmm
CHARMMEXEC=/net/orinoco/pga043/charmm/49a1_blade/charmm/build_charmm/charmm

DIR=`pwd`

export anadir=post

RUNDIR=$DIR/win$i

### Run the simulation
mkdir $RUNDIR/$anadir
cp -r prep $RUNDIR/
cd $RUNDIR

### timeout -s SIGINT 8h
echo $RUNDIR
mpirun -np 1 -x OMP_NUM_THREADS=1 --map-by node $CHARMMEXEC analysis=$anadir -i ../new_analysis.inp #> output_forflip.out 2> error

rm -rf prep
