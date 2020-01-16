md_tian2
========

md_tian2 (Molecular Dynamics Tian Xia 2) is a program for simulating
the scattering of atoms (and molecules) from a surface.

Purpose: -> maybe also here, but put in manual!
Do molecular dynamics, Langevin dynamics, Ring Polymer dynamics

Source code is in Fortran.

List of modules: -> maybe also here, but put in manual!
2do: add RuNNer related files in the following list

constants.f90           : contains the global constants
fit.f90                 : fitting routine for EMT
force.f90               : get the energies/forces in molecular dynamics simulation
geometry_opt.f90        : relax structure
md_algo.f90             : contains propagation algorithms for molecular dynamics simulation
md_init.f90             : initialize molecular dynamics simulation
md_tian2.f90            : main program
open_file.f90           : input/output routines
output_mod.f90          : output format routines
pes_emt_mod.f90         : contains the EMT potential
pes_ho_mod.f90          : contains the
pes_lj_mod.f90          : contains the Lennard-Jones potential
pes_nene_mod.f90        : contains the high-dimensional neural network potential
pes_nene_mod_supply.f90 : contains subroutines concerning the high-dimensional neural network potential
pes_non_mod.f90         : contains the non-interaction potential
pes_rebo_mod.f90        : contains the empirical reactive bond order potential
rpmd.f90                : contains the ring-polymer molecular dynamics simulation routine
run_config.f90          : initialize simulation parameters, read in input files
trajectory_info.f90     : collect information from the calculated trajectories
universe_mod.f90        : contains definitions of user types and all constants
useful_things.f90       : useful math routines




Input files: -> maybe also here, but put in manual!

md_tian.inp	        :control parameters defining the simulation conditions
<potential>.pes         :control parameters for the specific potential used
structure_file          :starting structure, poscar or mxt format possible
input.nn	        :control parameters for the high-dimensional neural network potential
scaling.data	        :contains the bias weights for the high-dimensional neural network potential
weights.XXX.data        :contains the weights for this element (XXX is the element number in the periodic table; e.g. 001 for H)


Compilation and linking (Intel fortran compiler and mkl library): -> maybe also here, but put in manual!

Change the settings in the makefile to your need. Options:

    make help 		print possible arguments for make
	make serial		serial version of md_tian2
	make clean		remove compilation related files


Program basic units -> maybe also here, but put in manual!

          Length   : Ang
          Time     : fs
          Energy   : eV

Program derived units

          Mass     : eV fs^2 / A^2 = 1/103.6382 amu
          Angle    : radian = 180 deg
          Distance : bohr = 0.5291772 Angstroem


The first working and tested version is put together February 18, 2014
on a Fassberg Hill in Dynamics at Surfaces Dep. of MPIbpc
to the flaming storm of applause muffled by thick institute building walls.

Credits:

Daniel J. Auerbach
Svenja Maria Janke
Marvin Kammler
Sascha Kandratsenka
Sebastian Wille


Dynamics at Surfaces Dep.
MPI for Biophysical Chemistry
Am Fassberg 11
37077 Goettingen
Germany

Dynamics at Surfaces Dep.
Institute for Physical Chemistry
Tammannstr. 6
37077 Goettingen
Germany

Md xia4 tian1 is a very important program. It helps to better the world.
jqrw sxrw=n! wr wj nA r sxrw nTrw!



Annotations from Sebastian Wille: -> put in manual!

md_tian.inp file:

  pip: projectile initial position
  pul: projectile upper limit (in /AA)
  T in POSCAR file stands for "True", means atom is NOT fixed and can move (F means atom is fixed and cannot move)
  Repetition of the slab in init file via: conf merge <path/file> <x_rep> <y_rep> (or conf poscar <path/file> <x_rep> <y_rep>)
