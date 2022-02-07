#!/bin/bash

# intention: used to get all the trajectories with asterisks in the mxt_fin file

for i in `seq 1 10`
do
touch recalc_star_traj.dat
grep -rl '[*]' traj_${i}/traj/. >> recalc_star_traj.dat #use dot instead of asterisk to avoid "argument list too long" error
done
