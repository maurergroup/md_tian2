#!/bin/bash

ntrajs=10000
script='md_tian.inp'

for i in `seq 1 10`
do
	j=$((1+((i-1)*${ntrajs})))
	cp ${script} ${script}_${i}
        sed -i -e "s/start.*/start $j/g" ${script}
        sed -i -e "s/ntrajs.*/ntrajs ${ntrajs}/g" ${script}
	nohup ./md_tian2.serial.x ${script}_${i} > output_${i}.out &
done
