#!/bin/bash

rm *.dat *.pdf *.png *.txt

END=20
for((p=1;p<=END;p++))
do

 echo 100 >> tmp.dat


END=20
for((i=1;i<=END;i++))
do

export j=`sed '1d' ../win"$p"/post/win"$i".dat | awk '{print $2}' | tr "\n" " "`

echo $j >> state"$p".dat

done

done

 tr "\n" " " < tmp.dat > vector.dat

paste -d " " state1.dat state2.dat state3.dat state4.dat state5.dat state6.dat state7.dat state8.dat state9.dat state10.dat state11.dat state12.dat state13.dat state14.dat state15.dat state16.dat state17.dat state18.dat state19.dat state20.dat > matrix.dat


