#!/bin/bash

array=( "D_2eV_43_D" "D_2eV_51_D" "D_2eV_59_D" "H_2eV_40_D" "H_2eV_50_D" "H_2eV_59_D" )

#folder=traj_gen8

for condition in ${array[*]}
do
	#cd ${condition}
	mkdir ${condition}_all
	cat ${condition}1/MXT2Summary.txt ${condition}2/MXT2Summary.txt | grep -v '#' > ${condition}_all/MXT2Summary.txt
	cd ${condition}_all
	#python ../scripts/1_CreateMXTSummary2.00.py
	python ../scripts/2_AnalyzePESTrajectory2.09.py
	gnuplot ../scripts/3_plotAnalysis2.10.plt
	gnuplot ../scripts/4_trypolar2.03.plt
	gnuplot ../scripts/5_trypolar2.04_with_points.plt
	cd ..
done
