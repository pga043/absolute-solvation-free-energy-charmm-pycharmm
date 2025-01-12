#!/bin/bash

#export run=2

read -p "Enter the run/window number: " run

## 
grep 'PERTRES' run"$run"/output_0

## EPRTOT is the total energy for this window and all previous (since a PERT RESET).
## EFORWARD is the energy for this current window.  
## EPREF is the initial energy difference (ef(0)-ei(0)) 
grep 'PERTURBATION> TI Windowing result' run"$run"/output_0

## EXPAVE is the time average of exp((ef(t)-ei(t)-ef(0)+ei(0))/kT)
## EXPFLC is the fluctuation of this value about its average
## DIFAVE is the time average of (ef(t)-ei(t)-ef(0)+ei(0))
## DIFFLC is the fluctuation of this value about its average
## Note: this value should not be much larger than kT for a
##       good window schedule.  If this value is too large,
##       then smaller window lambda steps should be used.
grep 'DIFFLC' run"$run"/output_0

#grep 'PAVE>' run"$run"/output_0 # look at 5th column, this value should be lower than it's corresponding fluctuation
#grep 'PFLC>' run"$run"/output_0 # look at 5th column, this is the fluctuation corresponding to "PAVE"



