#!/bin/bash

ntrajs=1000
script='md_tian'
nsteps=2000

for i in `seq 1 1000`
do
    j=$((((i-1)*${ntrajs})+1))
    echo $j
    cp ${script}.inp ${script}_${i}.inp
    sed -i -e "s/start.*/start ${j}/g" ${script}_${i}.inp
    sed -i -e "s/ntrajs.*/ntrajs ${ntrajs}/g" ${script}_${i}.inp
    sed -i -e "s/nsteps.*/nsteps ${nsteps}/g" ${script}_${i}.inp
    #nohup ./md_tian2.serial.x ${script}_${i}.inp > output_${i}.out &
done
