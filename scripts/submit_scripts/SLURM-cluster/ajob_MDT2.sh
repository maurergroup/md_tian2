#!/bin/bash
#SBATCH -p all
#SBATCH -N 1
#SBATCH --ntasks-per-node=1 # should match with requested mpirun cores
#SBATCH -t 2:00:00
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=swille2@gwdg.de
#SBATCH --job-name=D_2eV_50_d1
#SBATCH --output=output_D_2eV_50_d1.out
#SBATCH --array=1-10%10

module load intel/19.1.2.254
export I_MPI_FABRICS=shm
export FI_PROVIDER=tcp

# change lines below
FOLDERNAME="D_2eV_50_d1"
FOLDERPRE="gen8-d3_traj"


### SCRIPT ###

FOLDERID=${SLURM_ARRAY_TASK_ID}
OWD=`pwd`

# create scratch folder
mkdir -p /scratch2/swille/${FOLDERNAME}/${FOLDERPRE}_${FOLDERID}
# copy files to scratch
rsync -a ${FOLDERPRE}_${FOLDERID} /scratch2/swille/${FOLDERNAME}/
# change to scratch
cd /scratch2/swille/${FOLDERNAME}/${FOLDERPRE}_${FOLDERID}

# execute command
~/MDT2/md_tian2.serial.x md_tian.inp > output_${FOLDERID}.out

# copy results from scratch to home
rsync -a /scratch2/swille/${FOLDERNAME}/${FOLDERPRE}_${FOLDERID}/traj ${OWD}/${FOLDERPRE}_${FOLDERID}
rsync -a /scratch2/swille/${FOLDERNAME}/${FOLDERPRE}_${FOLDERID}/output_${FOLDERID}.out ${OWD}/${FOLDERPRE}_${FOLDERID}

# change to home (not necessary but clean)
cd ${OWD}
