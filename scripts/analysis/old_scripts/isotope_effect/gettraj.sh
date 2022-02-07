#!/bin/bash

# intention: old analyze script to see the progress of MD simulations

echo "" >> jobs.dat
for i in `seq 11 24`
do
    ls -lrta traj_$i/traj/ | tail -f -n 2 | grep mxt >> jobs.dat
done
