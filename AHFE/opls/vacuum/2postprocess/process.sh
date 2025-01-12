#!/bin/bash
rm result.txt forward_pmf.pdf reverse_pmf.pdf
touch result.txt

END=12 #`cat ../nwindows`
for((i=1;i<=END-1;i++))
do

# directories
export j=`echo "$i + 1" | bc`
#echo $i $j

#echo window"$i" - window"$j"

# same trajectory (same window), different lambda
# paste ../win"$i"/post/win"$i".dat ../win"$i"/post/win"$j".dat | sed '1d' | awk '{print $4 - $2}' > win_"$i"-"$j".dat
# paste ../win"$j"/post/win"$i".dat ../win"$j"/post/win"$j".dat | sed '1d' | awk '{print $2 - $4}' > win_"$j"-"$i".dat

# python plot.py win_"$i"-"$j".dat win_"$j"-"$i".dat

python mbar.py win_"$i"-"$j".dat win_"$j"-"$i".dat >> 00test.dat

#python zwanzig.py win_"$i"-"$j".dat win_"$j"-"$i".dat >> result.txt

# python distribution.py win_"$i"-"$j".dat


done

grep -i forward result.txt | awk '{print $6}' > forward.txt
grep -i reverse result.txt | awk '{print $6}' > reverse.txt

python get_free_energy.py
#awk -F',' '{sum+=$1;} END{print sum;}' 00test.dat 
