#!/bin/bash

# intention: remove all beforehand copied POSCAR files

for i in `seq 101 999`
do
	rm poscar_00000${i}.POSCAR 
done

rm poscar_00001000.POSCAR
