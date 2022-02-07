#!/bin/bash

# intention: get the latest written traj file to see the progress of the MD simulations

for i in `seq 1 10`
do
  cd traj_${i}
  ls conf/
  ls -lrt traj | tail -n 1
  ls traj | wc -l
  cd ..
done
