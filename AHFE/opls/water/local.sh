#!/bin/bash

END=8
for((i=1;i<=END;i++))
do
export i=$i

./run.sh

./analysis.sh 

done


