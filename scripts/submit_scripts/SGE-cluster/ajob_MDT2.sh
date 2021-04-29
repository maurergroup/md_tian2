#$ -S /bin/bash
#$ -l h_rt=4000:00:00
#$ -l mem_free=5G
#$ -m n
#$ -N CH_300K
#$ -cwd
#$ -R y
#$ -t 1-10
#$ -q all.q
#$ -tc 10
#$ -pe openmp 1

source /etc/profile.d/modules.sh
source /cm/shared/apps/intel/composer_xe/current/bin/compilervars.sh intel64
source /cm/shared/apps/intel/mpi/4.1.3.049/intel64/bin/mpivars.sh


#module load gcc
#module load sge
#module load intel/compiler/64/15.0/2015.3.187
#module load intel/mkl/64/11.2/2015.3.187
#module load shared

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export MKL_DYNAMIC=FALSE

# change lines below
ntrajs=1000
script='md_tian.inp'
folder='300K'
offset=0

### SCRIPT ###
JOBID=$((SGE_TASK_ID+offset))
START=$((((SGE_TASK_ID-1)*ntrajs)+1))
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
