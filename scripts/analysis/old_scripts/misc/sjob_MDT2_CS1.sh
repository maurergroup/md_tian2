#!/bin/bash
#SBATCH -p all
#SBATCH -N 1
#SBATCH --ntasks-per-node=1 # should match with requested mpirun cores
#SBATCH -t 2:00:00
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=swille2@gwdg.de
#SBATCH --job-name=D_2eV_43_D1
#SBATCH --output=output_D_2eV_43_D1.out
#SBATCH --array=1-50%25

module load intel/19.0.4.243
export I_MPI_FABRICS=shm
export FI_PROVIDER=tcp

OWD=`pwd`


FOLDERID=$SLURM_ARRAY_TASK_ID

mkdir -p /scratch2/swille/D_2eV_43_D1/traj_gen8_$FOLDERID
rsync -a traj_gen8_$FOLDERID /scratch2/swille/D_2eV_43_D1/
cd /scratch2/swille/D_2eV_43_D1/traj_gen8_$FOLDERID

./md_tian2.serial.x md_tian.inp > output.out

rsync -a /scratch2/swille/D_2eV_43_D1/traj_gen8_$FOLDERID/traj $OWD/traj_gen8_$FOLDERID
rsync -a /scratch2/swille/D_2eV_43_D1/traj_gen8_$FOLDERID/output.out $OWD/traj_gen8_$FOLDERID
