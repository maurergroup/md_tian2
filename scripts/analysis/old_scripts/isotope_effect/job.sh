#!/bin/bash

# intention: old submit script for MD simulations

for i in `seq 11 24`
do
    cd traj_$i
    nohup ./md_tian2.x md_tian.inp > output.dat &
    cd ..
done
