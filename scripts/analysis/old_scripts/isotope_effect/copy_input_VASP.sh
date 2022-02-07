#!/bin/bash

# intention: create input folder for 100 calculations with VASP

for i in `seq 1 100`
do
    cp -r dummy/ struc_$i
done

for j in `seq 1 9`
do
    cp struc_vasp/poscar_0000${j}000.dat struc_$j/POSCAR
done

for k in `seq 10 99`
do
    cp struc_vasp/poscar_000${k}000.dat struc_$k/POSCAR
done

cp struc_vasp/poscar_00100000.dat struc_100/POSCAR
