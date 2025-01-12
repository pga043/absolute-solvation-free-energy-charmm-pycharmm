#!/bin/bash

CHARMMEXEC=/net/orinoco/pga043/charmm/49a1_blade/charmm/build_charmm/charmm

DIR=`pwd`

RUNDIR=$DIR/run1

### Run the simulation
mkdir $RUNDIR
cp -r variables.inp prep $RUNDIR/
cd $RUNDIR

### timeout -s SIGINT 8h
echo "window$i started"
mpirun -np 1 -x OMP_NUM_THREADS=8 --bind-to none --map-by node $CHARMMEXEC seed=$RANDOM -i ../whampert.inp > output 2> error

$CHARMMEXEC -i ../analysis.inp > output_1



