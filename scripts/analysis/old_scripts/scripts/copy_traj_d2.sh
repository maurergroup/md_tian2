#!/bin/bash

folder=d2_traj_gen8

mkdir d2_traj
for i in `seq 7 10`
do
	cp ${folder}_${i}/traj/* d2_traj/
done
