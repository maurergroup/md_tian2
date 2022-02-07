#!/bin/bash

ntrajs=2000
script='md_tian'

for i in `seq 1 10`
do
	j=$((((i-1)*ntrajs)+1))
	cp ${script}.inp ${script}_${i}.inp
        sed -i -e "s/start.*/start $j/g" ${script}_${i}.inp
        sed -i -e "s/ntrajs.*/ntrajs ${ntrajs}/g" ${script}_${i}.inp
	#nohup ./md_tian2.serial.x ${script}_${i} > output_${i}.out &
done
