md_tian
=======

md_tian 2 (Molecular Dynamics Xia Tian 2) is a program for simulating 
the scattering of atoms (and molecules) from a surface. 



Source code is in Fortran.
List of modules:

2do: update source file list, add Manual to the main folder, add data to LICENSE.txt 

atom_class.f90		:defines of user types and all constants
*force.f90		:calculates energy and forces
*mdalgo.f90		:contains propagation algorithms
md_init.f90		:sets up everything
md_tian.f90		:the main file governing simulations
open_file.f90		:routines to open files smoothly
*output.f90		:output routines
run_config.f90	        :read in the input files
useful_things.f90	:useful math routines

Input files:

md_tian.inp	:control parameters defining the simulation conditions
*.nml		:contain emt-parameters for a species

Compilation and linking (Intel Fortran):

ifort -O3 -ipo -o md_tian atom_class.f90 open_file.f90 useful_things.f90 run_config.f90 md_init.f90 md_tian.f90

Update:

Change the settings in the Makefile to your need and type for the serial version

	make serial


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



Annotations from Sebastian Wille:

md_tian.inp file:

  pip: projectile initial position
  pul: projectile upper limit (in /AA)
  T in POSCAR file stands for "True", means atom is NOT fixed and can move (F means atom is fixed and cannot move)
  Repetition of the slab in init file via: conf merge <path/file> <x_rep> <y_rep> (or conf poscar <path/file> <x_rep> <y_rep>)

Neural network potentials:

In case you want to use the neural network potential implementation as used in the RuNNer program, please contact us for the corresponding .f90 file with all modules, subroutines etc.

Files needed for using the neural network potentials for molecular dynamics simulations:

1) input.nn containing the keywords, architecture and symmetry functions

2) scaling.data containing the bias values

3) All weights.XXX.data containing the weights for this element (XXX is the element number in the periodic table; e.g. 001 for H)
