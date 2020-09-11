#!/usr/bin/python

# intention: update RuNNer related files in a separate folder in the src folder of md_tian2

# use: python update_RuNNer_files.py <path/to/RuNNer/src/dir> <path/to/md_tian2/RuNNer/src/dir>

# 2do: add options and show help like in RuNNerUC.py

import sys
import os

runner_path=sys.argv[1]
md_tian_path=sys.argv[2]

# change to list? -> give entry number+file list?
file_list = \
       'abstime.f90 \
        addatoms.f90 \
        atomsymfunction1.f90 \
        atomsymfunction2.f90 \
        atomsymfunction3Andi.f90 \
        atomsymfunction3.f90 \
        atomsymfunction4.f90 \
        atomsymfunction5.f90 \
        atomsymfunction6.f90 \
        atomsymfunction8.f90 \
        atomsymfunction9.f90 \
        calconenn.f90 \
        checkelement.f90 \
        fileunits.f90 \
        fittingoptions.f90 \
        getatomsymfunctions.f90 \
        getcutoff.f90 \
        getdnodes_values.f90 \
        getlistdim.f90 \
        getneighboridxatomic_para.f90 \
        getneighborsatomic_para.f90 \
        getshortatomic.f90 \
        getvolume.f90 \
        globaloptions.f90 \
        mode1options.f90 \
        mpi_dummy.f90 \
        mpi_dummy_routines.f90 \
        mpifitdistribution.f90 \
        neighbor_para_short.f90 \
        nnconstants.f90 \
        nnewald.f90 \
        nnflags.f90 \
        nnshort_atomic.f90 \
        nuccharge.f90 \
        predictionoptions.f90 \
        predictionshortatomic.f90 \
        saturation.f90 \
        sortelements.f90 \
        sortsymfunctions.f90 \
        structures.f90 \
        symfunctions.f90 \
        timings.f90 \
        translate.f90'
# RuNNer source files in alphabetical order, add additional files if needed

runner_dir=os.chdir(runner_path) # change to RuNNer directory

print("Copying the following files:")

for i,name in enumerate(file_list.split()):
    #print(name)
    #print("Copying the following files")
    #copy_files_command = "cp "+str(name)+" "+md_tian_path # set up copy command
    copy_files_command = "cp "+str(name)+" "+md_tian_path
    print(name)
    #print(copy_files_command)
    os.system(copy_files_command) # copy files

print("done!")
