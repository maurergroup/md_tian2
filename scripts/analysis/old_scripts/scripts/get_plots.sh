#!/bin/bash

array=( "D_2eV_43_D1" "D_2eV_43_D2" "D_2eV_51_D1" "D_2eV_51_D2" "D_2eV_59_D1" "D_2eV_59_D2" "H_2eV_40_D1" "H_2eV_40_D2" "H_2eV_50_D1" "H_2eV_50_D2" "H_2eV_59_D1" "H_2eV_59_D2" )

folder=traj_gen8

for condition in ${array[*]}
do
	cd ${condition}
	mv all_traj traj
	python ../scripts/1_CreateMXTSummary2.00.py
	python ../scripts/2_AnalyzePESTrajectory2.09.py
	gnuplot ../scripts/3_plotAnalysis2.10.plt
	gnuplot ../scripts/4_trypolar2.03.plt
	gnuplot ../scripts/5_trypolar2.04_with_points.plt
	cd ..
done
