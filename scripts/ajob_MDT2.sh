#!/bin/bash

script='md_tian.inp'

for i in `seq 1 7`
do
	j=$((1+((i-1)*5000)))
        sed -i -e "s/start.*/start $j/g" ${script}
	nohup ~/Kai+Victoria/new_calc/md_tian2.serial.x ${script} > output_${i}.out &
	sleep 0.5
done
