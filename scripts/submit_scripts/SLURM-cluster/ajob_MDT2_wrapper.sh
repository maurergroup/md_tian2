#!/bin/bash

isotope_array=(H D)
energy_array=(2eV)
polar_array=(50)
domain_array=(d1 d2)

script="ajob_MDT2"

OWD=`pwd`

for isotope in ${isotope_array[*]}
do
    for energy in ${energy_array[*]}
    do
        for polar in ${polar_array[*]}
        do
            for domain in ${domain_array[*]}
            do
                curname=${isotope}_${energy}_${polar}_${domain} # get full folder name
                cd ${OWD}/${curname} 
                bash scripts/copy_input.sh # prepare input folder
                cp ${OWD}/scripts/${script}.sh ${script}_${curname}.sh # prepare array job script
                ascript=${script}_${curname}.sh
                sed -i -e "s/#SBATCH --job-name=.*/#SBATCH --job-name=${curname}/g" ${ascript}
                sed -i -e "s/#SBATCH --output=.*/#SBATCH --output=output_${curname}.out/g" ${ascript}
                sed -i -e "s/FOLDERNAME=.*/FOLDERNAME='${curname}'/g" ${ascript}
                sbatch ${ascript}
            done
        done
    done
done
