md_tian2
========

md_tian2 (Molecular Dynamics Tian Xia 2) is a program for simulating
the scattering of atoms (and molecules) from a surface.

Purpose: -> maybe also here, but put in manual!
Do molecular dynamics, Langevin dynamics, Ring Polymer dynamics

Source code is in Fortran.

List of modules: -> maybe also here, but put in manual!
2do: add RuNNer related files in the following list

	constants.f90		contains the global constants
	fit.f90			fitting routine for EMT
	force.f90		get the energies/forces in molecular dynamics simulation
	geometry_opt.f90	relax structure
	md_algo.f90		contains propagation algorithms for molecular dynamics simulation
	md_init.f90		initialize molecular dynamics simulation
	md_tian2.f90		main program
	open_file.f90		input/output routines
	output_mod.f90		output format routines
	pes_emt_mod.f90		contains the effective medium theory potential
	pes_ho_mod.f90		contains the harmonic oscillator potential
	pes_lj_mod.f90		contains the Lennard-Jones potential
	pes_nene_mod.f90	contains the high-dimensional neural network potential (external call of RuNNer)
	pes_non_mod.f90		contains the non-interaction potential
	pes_rebo_mod.f90	contains the reactive empirical bond order potential
	rpmd.f90		contains the ring-polymer molecular dynamics simulation routine
	run_config.f90		initialize simulation parameters, read in input files
	trajectory_info.f90	collect information from the calculated trajectories
	universe_mod.f90	contains definitions of user types and all constants
	useful_things.f90	useful math routines
        write_info.f90          write header and information concerning program routine




Input files: -> maybe also here, but put in manual!

md_tian.inp	        :control parameters defining the simulation conditions
<potential>.pes         :control parameters for the specific potential used
structure_file          :starting structure, poscar or mxt format possible
input.nn	        :control parameters for the high-dimensional neural network potential
scaling.data	        :contains the bias weights for the high-dimensional neural network potential
weights.XXX.data        :contains the weights for this element (XXX is the element number in the periodic table; e.g. 001 for H)


Compilation and linking (Intel fortran compiler and mkl library): -> maybe also here, but put in manual!

In case the gfortran compiler is used you have to insert in the Makefile for the LIB entry:

-Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_gf_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl

for static or

-L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl

for dynamic linking. Check also the intel link advisor webpage under:

https://software.intel.com/content/www/us/en/develop/articles/intel-mkl-link-line-advisor.html

Although it is possible to use GNU compiler, we strongly recommend to use the intel compiler due to performance reasons!

Example calculation:
To check for the technical and numerical stability of the program, one can run the simulation in the example/md/ folder 
and compare the results one gets with the output archive. For this purpose, a NN PES was fitted to a reduced data set. 
The corresponding input and output can be found in the example/fit/ folder.

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

md_tian.inp file: -> add a list with the complete set of keywords with descriptions

  T in POSCAR file stands for "True", means atom is NOT fixed and can move (F means atom is fixed and cannot move)
  Repetition of the slab in init file via: conf merge <path/file> <x_rep> <y_rep> (or conf poscar <path/file> <x_rep> <y_rep>)

1) Following the complete list of all input files:

    i) md_tian.inp
    ii) pes/nene-HC.pes

    i) Following a complete list of keywords for md_tian.inp:

		a) run
			defines whether to 
				1) minimize a structure using fire algorithm (min)
			 	2) in case of the REBO potential, one can fit the parameters (fit) (check md_tian2.f90)
				3) perform molecular dynamics simulations (md) 
		b) start: gives the number of trajectory to start with (we use this number to feed the RNG, therefore if one wants to recalculate a trajectory, the number can be given here; also useful to split large jobs into small chunks)
		c) ntrajs: total number of trajectories to launch
		d) nsteps: maximum number of steps for the simulation (this prevents infinite simulation time in case of adsorption, because there is no stopping criterion compared to scattering)
		e) step: timestep for the propagator (unit is in fs (femtoseconds))
		f) projectile: the first number gives the different types of projectiles (all elements) followed by the element symbol, the mass and the type of propagator to use for this species
		g) lattice: same as for projectile, but here the lattice types are given
		h) force_beads: set the number of beads using in RPMD
		i) pile_tau: when using the path integral Langevin equation (PILE) this set the thermostat time constant
		j) andersen_time: when using the andersen thermostat the number set the collision time (default is 30)
		k) Einc: incidence kinetic energy of the projectile
		l) polar: incidence polar angle of the projectile
		m) azimuth: incidence azimuth angle of the projectile
		n) pip: projectile initial position. One can either set the position or use 'r' to use a RNG to set the position (x,y,z in Ang)
		o) pul: projectile upper limit. Here one can set the z coordinate of the projectile after scattering to stop the simulation (saves some time when simulating trajectories)
		p) Tsurf: set the temperature of the surface
		q) Tproj: set the temperature of the projectile (in case of a single atom this will do nothing)
		r) annealing: set parameters when using simulated annealing
		s) conf: first entry is the format, second the path to the file as string and in case of merge the number of total structures to choose from (RNG decides which structure to take, check hints below)
		t) pes: you have to give the path to the PES file as string
		u) output: give the output type and how often you want it (1 means every time step, 10 every 10th timestep and so on), check the list below which types are available
		v) fit_training_data: path to training data
		w) fit_validation_data: path to validation data
		x) evasp: reference energy for fit
		y) maxit: maximum number of iterations during fit
		z) nthreads: number of threads used for fitting
		aa) adsorption_distance: define distance of projectile to surface (start and end) to get interaction time 
		ab) rng_type: we have two types of RNG, default is 1 if you give no keyword
		ac) debug: debugging, only interesting for developers

        Annotations to points from list:

                f,g) possible propagators:
			1) ver: standard velocity verlet algorithm
			2) lan: langevin dynamics algorithm
			3) and: andersen thermostat
			4) pil: path integral Langevin equation (PILE) thermostat
		r) annealing: 
	             example: annealing <Tmax> <steps> <interval> (max surface T in K, number of steps per cycle, number of steps per T interval)
		s) possible ways to read in structure(s):
			1) poscar: you have to give the path to the file as string 
			2) mxt: you have to give the path to the file as string
			3) merge: this also uses the mxt format, but the structure is split into projectile and lattice configuration files.
                             example: conf merge '<path_to_projectile>' 1 '<path_to_lattice>' 1000 (string with path to projectile folder, number of projectile files, string with path to lattice folder, number of lattice files)
		u) possible output formats (the number after each keyword tell the frequency in time steps to write this format):
			 1) scatter: (only 1 is possible as option; this has to be always present in MD simulations)
			 2) energy: (write all energies, this will create a lot of data, so be careful)
			 3) xyz: (simple xyz format, one structure file is generated)
			 4) poscar: (for each step a separate POSCAR file is written, contains the beads)
			 5) vasp: (like the poscar format, but without the beads. This can be directly read in any GUI.)
			 6) mxt: (for each step a separate file is written in mxt format)
			 7) adsorption_status: (will write a file with the information if the projectile is adsorbed or not for each step during the simulation)
			 8) nene: (prints every timestep the current step on screen, debugging purpose for the NN PES)
			 9) aims: (for each step a separate geometry.in file is written)
			10) runner: (one input.data file is written in the RuNNer format)
                        11) beads: in RPMD with this keyword every single bead will be written to a separate file (depending on output format keyword(s) given); without this keyword, only the center of mass structure will be written

    ii) Following a complete list of keywords for pes/<potential_name>.pes:

                a) Here it will be defined, which interactions and what potential is used to describe them. Furthermore, we have to set if the element is of type projectile (proj) or part of the surface/bulk (latt).
                   For all possible element permutations we have to define interactions. Depending on the potential used, specific parameters have to be given. 
                   In the following an example for H on graphene is shown using a high-dimensional neural-network potential. In this case, interactions between H and C, H and H in case of multiple projectiles as well as C and C can occur. 
                   For all cases, the neural-networks potential should be used (keyword nene). For the parameters, we have to give the folder to look up for the potential specific files
                   The parameters are given with an indentation of a single white space.
                   Optional, the keyword maxnum_extrapolation_warnings together with an integer will set the threshold when to stop the program if the number for warnings is exceeded.
                   The parameters are given with an indentation of a single white space.

                pes nene
                H       C       proj    latt

                pes nene
                C       C       latt    latt

                pes nene
                H       H       proj    proj
                 inp_dir 'RuNNer_input/'
                 maxnum_extrapolation_warnings 100
