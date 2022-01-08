#!/bin/bash

array=( "D_2eV_43_D1" "D_2eV_43_D2" "D_2eV_51_D1" "D_2eV_51_D2" "D_2eV_59_D1" "D_2eV_59_D2" "H_2eV_40_D1" "H_2eV_40_D2" "H_2eV_50_D1" "H_2eV_50_D2" "H_2eV_59_D1" "H_2eV_59_D2" )

folder=traj_gen8

for condition in ${array[*]}
do
	cd ${condition}
	mkdir -p traj
	for i in `seq 1 100`
	do
	        cp ${folder}_${i}/traj/* traj/
	done
	cd ..
done
