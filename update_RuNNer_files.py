#!/usr/bin/env python

# intention: update RuNNer related files in the src folder of md_tian2

# use: python update_RuNNer_files.py <path/to/RuNNer/src/dir> <path/to/md_tian2/src/dir>

import sys
import os

runner_path=sys.argv[1]
md_tian_path=sys.argv[2]

# change to list? -> give entry number+file list?
file_list = 'abstime.f90 atomsymfunction1.f90 atomsymfunction2.f90 atomsymfunction3Andi.f90 atomsymfunction3.f90 atomsymfunction4.f90 \
atomsymfunction5.f90 atomsymfunction6.f90 atomsymfunction8.f90 atomsymfunction9.f90 calconenn.f90 fileunits.f90 getatomsymfunctions.f90 \
getdnodes_values.f90 getneighboridxatomic_para.f90 getneighborsatomic_para.f90 getshortatomic.f90 globaloptions.f90 \
mpi_dummy.f90 mpi_dummy_routines.f90 mpifitdistribution.f90 neighbor_para_short.f90 nnconstants.f90 nnflags.f90 nnshort_atomic.f90 \
predictionoptions.f90 predictionshortatomic.f90 saturation.f90 symfunctions.f90' # RuNNer source files in alphabetical order, add additional files if needed

runner_dir=os.chdir(runner_path) # change to RuNNer directory

copy_files_command = "cp "+str(file_list)+" "+md_tian_path # set up copy command

os.system(copy_files_command) # copy files
