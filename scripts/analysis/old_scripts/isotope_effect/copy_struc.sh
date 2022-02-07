#!/bin/bash

# intention: copy needed POSCAR files to input folder for VASP calculations

for i in `seq 0 1`
do
    for j in `seq 0 9`
    do
        for k in `seq 0 9`
        do
            cp conf/poscar_00${i}${j}${k}000.dat struc_vasp/
        done
    done
done
