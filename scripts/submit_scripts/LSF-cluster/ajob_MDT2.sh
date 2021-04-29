#!/bin/sh

#BSUB -q mpi
#BSUB -n 1
#BSUB -R scratch
#BSUB -W 01:00
#BSUB -o mode1.out

mkdir /scratch/swille
export TMPDIR="/scratch/swille/"
export TMPDIR4="/scratch/swille/"
export PATH=$HOME/bin:${PATH}
cd $HOME/md_tian2/

# change lines below
ntrajs=1000
script='md_tian.inp'
folder='300K'
offset=0

### SCRIPT ###
JOBID=$((LSB_JOBINDEX+offset))
START=$((((LSB_JOBINDEX-1)*ntrajs)+1))
OWD=`pwd`

# create folder on scratch
mkdir -p /tmp/swille3/${folder}
# copy files to scratch folder
rsync -a ${folder} /tmp/swille3/
# change to scratch folder
cd /tmp/swille3/${folder}
# modify input file
sed -i -e "s/start.*/start ${START}/g" ${script}
sed -i -e "s/ntrajs.*/ntrajs ${ntrajs}/g" ${script}

# execute MDT2 and write to output file
~/MDT2/md_tian2.serial.x ${script} > output_${JOBID}.out

# copy results from scratch to home
rsync -a /tmp/swille3/${folder}/traj ${OWD}/${folder}
rsync -a /tmp/swille3/${folder}/output_${JOBID}.out ${OWD}/${folder}

# change to home (not really necessary but clean)
cd ${OWD}
