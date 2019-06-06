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

FC     = ifort 
FFLAGS = -O3 -mkl -xHost -fp-model precise -132
LIB    =

###############################################
# Changes below this line developers only
###############################################

SUBROUTINES = md_tian.o \


MODULES      = constants.o \
        open_file.o \
        useful_things.o \
        universe_mod.o \
        run_config.o \
        pes_non_mod.o \
        pes_ho_mod.o \
        pes_lj_mod.o \
        pes_emt_mod.o \
        pes_rebo_mod.o \
        pes_nene_mod.o \
        force.o \
        rpmd.o \
        md_init.o \
        md_algo.o \
        trajectory_info.o \
        output_mod.o \
        fit.o \
        geometry_opt.o \


# begin make

all:
	cat compile.msg

# compile modules
constants.o: constants.f90
	rm constants.mod ; \
	$(FC) $(FFLAGS_MPI) -c constants.f90

open_file.o: open_file.f90
	rm open_file.mod ; \
        $(FC) $(FFLAGS_MPI) -c open_file.f90

useful_things.o: useful_things.f90
	rm useful_things.mod ; \
        $(FC) $(FFLAGS_MPI) -c useful_things.f90

universe_mod.o: universe_mod.f90
	rm universe_mod.mod ; \
        $(FC) $(FFLAGS_MPI) -c universe_mod.f90

run_config.o: run_config.f90
	rm run_config.mod ; \
        $(FC) $(FFLAGS_MPI) -c run_config.f90

pes_non_mod.o: pes_non_mod.f90
	rm pes_non_mod.mod ; \
        $(FC) $(FFLAGS_MPI) -c pes_non_mod.f90

pes_ho_mod.o: pes_ho_mod.f90
	rm pes_ho_mod.mod ; \
        $(FC) $(FFLAGS_MPI) -c pes_ho_mod.f90

pes_lj_mod.o: pes_lj_mod.f90
	rm pes_lj_mod.mod ; \
        $(FC) $(FFLAGS_MPI) -c pes_lj_mod.f90

pes_emt_mod.o: pes_emt_mod.f90
	rm pes_emt_mod.mod ; \
        $(FC) $(FFLAGS_MPI) -c pes_emt_mod.f90

pes_rebo_mod.o: pes_rebo_mod.f90
	rm pes_rebo_mod.mod ; \
        $(FC) $(FFLAGS_MPI) -c pes_rebo_mod.f90

pes_nene_mod.o: pes_nene_mod.f90
	rm pes_nene_mod.mod ; \
        $(FC) $(FFLAGS_MPI) -c pes_nene_mod.f90

force.o: force.f90
	rm force.mod ; \
	$(FC) $(FFLAGS_MPI) -c force.f90

rpmd.o: rpmd.f90
	rm rpmd.mod ; \
        $(FC) $(FFLAGS_MPI) -c rpmd.f90

md_init.o: md_init.f90
	rm md_init.mod ; \
        $(FC) $(FFLAGS_MPI) -c md_init.f90

md_algo.o: md_algo.f90
	rm md_algo.mod ; \
        $(FC) $(FFLAGS_MPI) -c md_algo.f90

trajectory_info.o: trajectory_info.f90
	rm trajectory_info.mod ; \
        $(FC) $(FFLAGS_MPI) -c trajectory_info.f90

output_mod.o: output_mod.f90
	rm output_mod.mod ; \
        $(FC) $(FFLAGS_MPI) -c output_mod.f90

fit.o: fit.f90
	rm fit.mod ; \
	$(FC) $(FFLAGS_MPI) -c fit.f90

geometry_opt.o: geometry_opt.f90
	rm geometry_opt.mod ; \
        $(FC) $(FFLAGS_MPI) -c geometry_opt.f90



# compile subroutines

md_tian.o: md_tian.f90
	$(FC) $(FFLAGS) -c md_tian.f90

#serial: ./md_tian2.x
#./md_tian2.x: $(MODULES) $(SUBROUTINES) 

serial: ./md_tian2.x
./md_tian2.x: $(MODULES) $(SUBROUTINES)	
	$(FC) $(FFLAGS) -o $@ \
               $(SUBROUTINES) $(MODULES) $(LIB)

clean:
	rm -f *.o *.mod *.x
