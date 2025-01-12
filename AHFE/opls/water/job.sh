#!/bin/bash

END=15
for((i=9;i<=END;i++))
do
export i=$i

./run.sh

./analysis.sh 

done


