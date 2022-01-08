#!/bin/bash
datafile=${1-input.data}

bohr2ang=0.529177249

i=0

for file in `grep atom $datafile | awk '{print $5}' | sort | uniq`
do
	i=`echo $i +1 | bc -ql`

	elem[$i]=$file
	elem_count[$i]=`grep atom $datafile | grep $file | wc -l`
done

n_elems=${#elem[@]}

echo ${elem[@]}
printf "   1.00\n"
grep lattice $datafile | awk -v bohr2ang=$bohr2ang '{printf "%18.8f %18.8f %18.8f\n", $2*bohr2ang, $3*bohr2ang, $4*bohr2ang}'
echo ${elem_count[@]}
printf "Cartesian\n"

for file in `echo ${elem[@]}`
do
	grep atom $datafile | grep $file | awk -v bohr2ang=$bohr2ang '{printf "%18.8f %18.8f %18.8f\n", $2*bohr2ang, $3*bohr2ang, $4*bohr2ang}'
done


