# Makefile for md_tian2

##############################################
# Compiler:
# Intel ifort compiler
#
# Library:
# Intel mkl library
##############################################

#############
### Settings:
#############

# ifort and mkl library
FC     = ifort 
FFLAGS = -O3 -mkl -xcore-avx2 -fp-model precise -132
LIB    =

###############################################
# Changes below this line developers only
###############################################

DEBUG  = -fpe0 -traceback -O0 -g -check all -check bounds -check uninit -ftrapuv -debug all -gen-interfaces -warn interfaces -mkl

RUNNER_SRC_DIR = RuNNer_src

SUBROUTINES = md_tian.o \
	      abstime.o \
	      write_info.o \


MPI_DUMMY_SUBROUTINES = mpi_dummy_routines.o \

MODULES      =  constants.o \
		open_file.o \
		useful_things.o \
               	universe_mod.o \
               	run_config.o \
               	pes_non_mod.o \
              	pes_ho_mod.o \
               	pes_lj_mod.o \
               	pes_emt_mod.o \
               	pes_rebo_mod.o \
               	fileunits.o \
               	predictionoptions.o \
	       	nnflags.o \
               	nnconstants.o \
               	globaloptions.o \
               	inputnncounters.o \
	       	mode1options.o \
	       	nnshort_atomic.o \
	       	nnewald.o \
	       	nnshort_pair.o \
		saturation.o \
		structures.o \
		symfunctions.o \
		timings.o \
		fittingoptions.o \
		pes_nene_mod_supply.o \
	        pes_nene_mod.o \
	        force.o \
	        rpmd.o \
	        md_init.o \
	        md_algo.o \
	        trajectory_info.o \
	        output_mod.o \
	        fit.o \
	        geometry_opt.o \

MODULES_MPI_DUMMY = mpi_dummy.o


# begin make

all:
	cat compile.msg

# compile modules
constants.o: constants.f90
	rm constants.mod ; \
	$(FC)  -c constants.f90

open_file.o: open_file.f90
	rm open_file.mod ; \
        $(FC)  -c open_file.f90

useful_things.o: useful_things.f90
	rm useful_things.mod ; \
        $(FC)  -c useful_things.f90

universe_mod.o: universe_mod.f90
	rm universe_mod.mod ; \
        $(FC)  -c universe_mod.f90

run_config.o: run_config.f90
	rm run_config.mod ; \
        $(FC)  -c run_config.f90

pes_non_mod.o: pes_non_mod.f90
	rm pes_non_mod.mod ; \
        $(FC)  -c pes_non_mod.f90

pes_ho_mod.o: pes_ho_mod.f90
	rm pes_ho_mod.mod ; \
        $(FC)  -c pes_ho_mod.f90

pes_lj_mod.o: pes_lj_mod.f90
	rm pes_lj_mod.mod ; \
        $(FC)  -c pes_lj_mod.f90

pes_emt_mod.o: pes_emt_mod.f90
	rm pes_emt_mod.mod ; \
        $(FC)  -c pes_emt_mod.f90

pes_rebo_mod.o: pes_rebo_mod.f90
	rm pes_rebo_mod.mod ; \
        $(FC)  -c pes_rebo_mod.f90

fileunits.o: $(RUNNER_SRC_DIR)/fileunits.f90
	rm fileunits.mod ; \
        $(FC) $(FFLAGS_MPI) -c fileunits.f90

fittingoptions.o: fittingoptions.f90
	rm fittingoptions.mod ; \
        $(FC) $(FFLAGS_MPI) -c fittingoptions.f90

nnflags.o: nnflags.f90
	rm nnflags.mod ; \
        $(FC) $(FFLAGS_MPI) -c nnflags.f90

nnconstants.o: nnconstants.f90
        rm nnconstants.mod ; \
        $(FC) $(FFLAGS_MPI) -c nnconstants.f90

globaloptions.o: globaloptions.f90
        rm globaloptions.mod ; \
        $(FC) $(FFLAGS_MPI) -c globaloptions.f90

inputnncounters.o: inputnncounters.f90
        rm inputnncounters.mod ; \
        $(FC) $(FFLAGS_MPI) -c inputnncounters.f90

mode1options.o: mode1options.f90
        rm mode1options.mod ; \
        $(FC) $(FFLAGS_MPI) -c mode1options.f90

nnshort_atomic.o: nnshort_atomic.f90
        rm nnshort_atomic.mod ; \
        $(FC) $(FFLAGS_MPI) -c nnshort_atomic.f90

nnewald.o: nnewald.f90
        rm nnewald.mod ; \
        $(FC) $(FFLAGS_MPI) -c nnewald.f90

nnshort_pair.o: nnshort_pair.f90
        rm nnshort_pair.mod ; \
        $(FC) $(FFLAGS_MPI) -c nnshort_pair.f90

saturation.o: saturation.f90
        rm saturation.mod ; \
        $(FC) $(FFLAGS_MPI) -c saturation.f90

structures.o: structures.f90
        rm structures.mod ; \
        $(FC) $(FFLAGS_MPI) -c structures.f90

mpi_module.o: mpi_module.f90
        rm mpi_module.mod ; \
        $(MPIFC) $(FFLAGS_MPI) -c mpi_module.f90

symfunctions.o: symfunctions.f90
        rm symfunctions.mod ; \
        $(FC) $(FFLAGS_MPI) -c symfunctions.f90

timings.o: timings.f90
	rm timings.mod ; \
        $(FC) $(FFLAGS_MPI) -c timings.f90

predictionoptions.o: predictionoptions.f90
	rm predicitonoptions.mod ; \
        $(FC) $(FFLAGS_MPI) -c predictionoptions.f90

mpi_dummy.o: mpi_dummy.f90
	rm mpi_dummy.mod ; \
        $(FC) $(FFLAGS) -c mpi_dummy.f90

pes_nene_mod.o: pes_nene_mod.f90
	rm pes_nene_mod.mod ; \
        $(FC)  -c pes_nene_mod.f90

force.o: force.f90
	rm force.mod ; \
	$(FC)  -c force.f90

rpmd.o: rpmd.f90
	rm rpmd.mod ; \
        $(FC)  -c rpmd.f90

md_init.o: md_init.f90
	rm md_init.mod ; \
        $(FC)  -c md_init.f90

md_algo.o: md_algo.f90
	rm md_algo.mod ; \
        $(FC)  -c md_algo.f90

trajectory_info.o: trajectory_info.f90
	rm trajectory_info.mod ; \
        $(FC)  -c trajectory_info.f90

output_mod.o: output_mod.f90
	rm output_mod.mod ; \
        $(FC)  -c output_mod.f90

fit.o: fit.f90
	rm fit.mod ; \
	$(FC)  -c fit.f90

geometry_opt.o: geometry_opt.f90
	rm geometry_opt.mod ; \
        $(FC)  -c geometry_opt.f90



# compile subroutines

md_tian.o: md_tian.f90
	$(FC) $(FFLAGS) -c md_tian.f90




serial: ./md_tian2.x
./md_tian2.x: $(MODULES_MPI_DUMMY) $(MODULES) $(SUBROUTINES) $(MPI_DUMMY_SUBROUTINES)
	$(FC) $(FFLAGS) -o $@ \
               $(SUBROUTINES) $(MODULES_MPI_DUMMY) $(MODULES) $(MPI_DUMMY_SUBROUTINES) $(LIB)

debug: ./md_tian2_debug.x
./md_tian2_debug.x: $(MODULES) $(SUBROUTINES)
	$(FC) $(DEBUG) -o $@ \
		$(SUBROUTINES) $(MODULES) $(LIB)

clean:
	rm -f *.o *.mod *.x

help:
	@echo ''
	@echo 'make serial		serial version, default'
	@echo 'make clean		delete all object, module and binary files'
	@echo ''
