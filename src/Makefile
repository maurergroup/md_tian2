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
FFLAGS = -mkl -132 -check all
LIB    = 

###############################################
# Changes below this line developers only
###############################################

SUBROUTINES = md_tian.o \

PES_MODULES =         pes_non_mod.o \
        pes_ho_mod.o \
        pes_lj_mod.o \
        pes_emt_mod.o \
        pes_rebo_mod.o \


MODULES      = constants.o \
        open_file.o \
        useful_things.o \
        universe_mod.o \
        run_config.o \
        force.o \
        rpmd.o \
        md_init.o \
        md_algo.o \
        trajectory_info.o \
        output_mod.o \
        fit.o \
        geometry_opt.o \
        ODFriction_mod.o\
$(PES_MODULES)



# begin make


# compile modules
constants.o: constants.f90
	rm constants.mod ; \
	$(FC) -c constants.f90

open_file.o: open_file.f90
	rm open_file.mod ; \
        $(FC)  -c open_file.f90

useful_things.o: useful_things.f90 constants.o
	rm useful_things.mod ; \
        $(FC)  -c useful_things.f90

universe_mod.o: universe_mod.f90 useful_things.o
	rm universe_mod.mod ; \
        $(FC)  -c universe_mod.f90

run_config.o: run_config.f90 universe_mod.o useful_things.o 
	rm run_config.mod ; \
        $(FC)  -c run_config.f90

pes_non_mod.o: pes_non_mod.f90 universe_mod.o useful_things.o run_config.o
	rm pes_non_mod.mod ; \
        $(FC)  -c pes_non_mod.f90

pes_ho_mod.o: pes_ho_mod.f90 constants.o universe_mod.o useful_things.o run_config.o
	rm pes_ho_mod.mod ; \
        $(FC)  -c pes_ho_mod.f90

pes_lj_mod.o: pes_lj_mod.f90 constants.o universe_mod.o useful_things.o run_config.o
	rm pes_lj_mod.mod ; \
        $(FC)  -c pes_lj_mod.f90

pes_emt_mod.o: pes_emt_mod.f90 universe_mod.o useful_things.o run_config.o
	rm pes_emt_mod.mod ; \
        $(FC)  -c pes_emt_mod.f90

pes_rebo_mod.o: pes_rebo_mod.f90 universe_mod.o
	rm pes_rebo_mod.mod ; \
        $(FC)  -c pes_rebo_mod.f90


force.o: force.f90 rpmd.o $(PES_MODULES)
	rm force.mod ; \
	$(FC)  -c force.f90

rpmd.o: rpmd.f90 universe_mod.o constants.o useful_things.o run_config.o
	rm rpmd.mod ; \
        $(FC)  -c rpmd.f90

md_init.o: md_init.f90 ODFriction_mod.o  rpmd.o run_config.o $(PES_MODULES)
	rm md_init.mod ; \
        $(FC)  -c md_init.f90

md_algo.o: md_algo.f90 ODFriction_mod.o run_config.o rpmd.o $(PES_MODULES)
	rm md_algo.mod ; \
        $(FC)  -c md_algo.f90

trajectory_info.o: trajectory_info.f90 universe_mod.o constants.o rpmd.o useful_things.o open_file.o run_config.o
	rm trajectory_info.mod ; \
        $(FC)  -c trajectory_info.f90

output_mod.o: output_mod.f90 trajectory_info.o pes_rebo_mod.o 
	rm output_mod.mod ; \
        $(FC)  -c output_mod.f90

fit.o: fit.f90 force.o
	rm fit.mod ; \
	$(FC)  -c fit.f90

geometry_opt.o: geometry_opt.f90 force.o rpmd.o output_mod.o md_algo.o
	rm geometry_opt.mod ; \
        $(FC)  -c geometry_opt.f90

ODFriction_mod.o: ODFriction_mod.f90
	rm ODFriction_mod.mod; \
	$(FC) -c ODFriction_mod.f90

# compile subroutines

md_tian.o: md_tian.f90 force.o md_algo.o fit.o trajectory_info.o geometry_opt.o
	$(FC) $(FFLAGS) -c md_tian.f90


serial: ./md_tian2.x
./md_tian2.x: $(MODULES) $(SUBROUTINES)	
	$(FC) $(FFLAGS) -o $@ \
               $(SUBROUTINES) $(MODULES) $(LIB)

clean:
	rm -f *.o *.mod *.x

help:
	@echo ''
	@echo 'make serial		serial version, default'
	@echo 'make clean		delete all object, module and binary files'
	@echo ''
