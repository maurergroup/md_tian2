module atom_class
    !
    ! Purpose:
    !           This module contains definitions of user types and all constants
    !
    ! Date          	Author          	    History of Revison
    ! ====          	======          	    ==================
    ! 30.03.2017    	Marvin Kammler		    new data structure
    !                   Sascha Kandratsenka
    ! 18.02.2014    	Svenja M. Janke		    Original
    !			        Sascha Kandratsenka
    !			        Dan J. Auerbach

    use constants

    implicit none


    !  Type atoms
    !   structure to hold the position, velocity, force etc. for multiple atoms
    !       use rank2 array so positions, velocities, forces etc. are stored
    !       in sequentional memory locations for efficient access
    !       each array should be allocated (3, n_beads, n_atom)
    !       mass array has length of n_atoms
    type atoms

        integer                       :: natoms          ! number of atoms
        integer                       :: nbeads          ! number of beads per atom
        real(dp),         allocatable :: m(:)            ! mass
        integer,          allocatable :: atn(:)          ! atomic number
        character(len=3), allocatable :: name(:)         ! atomic name
        real(dp), allocatable         :: r(:,:,:)        ! positions
        real(dp), allocatable         :: v(:,:,:)        ! velocities
        real(dp), allocatable         :: f(:,:,:)        ! forces
        logical,  allocatable         :: is_fixed(:,:,:) ! mask array defining frozen atoms (T is frozen)
        logical,  allocatable         :: is_proj(:)      ! to distinguish between proj and latt
        integer,  allocatable         :: idx             ! program-wide index of atom type


    end type atoms

    type simulation_parameters

        integer :: start                                            ! a trajectory to start with
        integer :: ntrajs                                           ! number of trajectories
        integer :: nsteps                                           ! number of steps
        real(dp):: step                                             ! time step in fs
        integer :: nlattices                                        ! number of lattice species
        integer :: nprojectiles                                     ! number of incident species
        character(len=3), allocatable :: name_l(:), name_p(:)       ! atomic names
        real(dp),         allocatable :: mass_l(:), mass_p(:)       ! atomic masses
        integer,          allocatable :: md_algo_l(:), md_algo_p(:) ! and respective key
        real(dp),         allocatable :: einc(:)                    ! incidence energy (eV)
        real(dp),         allocatable :: inclination(:)             ! incidence polar angle (degree)
        real(dp),         allocatable :: azimuth(:)                 ! incidence azimuthal angle (degree)
        real(dp):: Tsurf                                            ! surface temperature in K
        real(dp):: sa_Tmax                                          ! max. Tsurf for simulated annealing in K
        integer :: sa_nsteps                                        ! number of steps per simulated annealing cycle
        integer :: sa_interval                                      ! number of steps per temperature interval
        character(len=7)    :: confname                             ! configuration key
        character(len=max_string_length) :: confname_file           ! name of the system configuration file or folder
        integer :: rep(2)                                           ! defines in-plane repetitions
        integer :: nconfs                                           ! number of configurations to read in
        character(len=max_string_length) :: pes_file                ! name of the file that stores the potential parameters
        character(len=3) :: run                                     ! what to do
        integer, dimension(2) :: output                             ! what to save
        character(len=15) :: pip(3)                                 ! determine initial projectile position

    end type

    interface atoms
        module procedure new_atoms
    end interface

    interface simulation_parameters
        module procedure new_simulation_parameters
    end interface

contains

    ! Constructor for type atoms
    !    input:  n_beads, n_atoms
    !    allocates arrays as (3,n_beads,n_atom)

    function new_atoms(nbeads, natoms)
        integer, intent(in) :: nbeads, natoms
        type(atoms) new_atoms

        allocate(new_atoms%m(natoms))
        allocate(new_atoms%atn(natoms))
        allocate(new_atoms%name(natoms))
        allocate(new_atoms%r(3,nbeads,natoms))
        allocate(new_atoms%v(3,nbeads,natoms))
        allocate(new_atoms%f(3,nbeads,natoms))
        allocate(new_atoms%is_fixed(3,nbeads,natoms))

        new_atoms%nbeads = nbeads
        new_atoms%natoms = natoms
        new_atoms%m     = -1.0_dp
        new_atoms%atn   = 0
        new_atoms%name  = ""
        new_atoms%r     = 0.0_dp
        new_atoms%v     = 0.0_dp
        new_atoms%f     = 0.0_dp
        new_atoms%is_fixed = .FALSE.


    end function

    ! Constructor for type simulation_parameters

    function new_simulation_parameters()

        type(simulation_parameters) new_simulation_parameters

        new_simulation_parameters%start  = default_int
        new_simulation_parameters%ntrajs = default_int
        new_simulation_parameters%nsteps = default_int
        new_simulation_parameters%step   = default_real
        new_simulation_parameters%nlattices = default_int
        new_simulation_parameters%nprojectiles = default_int
        new_simulation_parameters%Tsurf   = default_real
        new_simulation_parameters%sa_Tmax   = default_real
        new_simulation_parameters%sa_nsteps = default_int
        new_simulation_parameters%sa_interval = default_int
        new_simulation_parameters%confname = default_string
        new_simulation_parameters%rep = [0,0]
        new_simulation_parameters%nconfs  = default_int
        new_simulation_parameters%pes_file = default_string
        new_simulation_parameters%run = default_string
        new_simulation_parameters%output = [default_int,default_int]
        new_simulation_parameters%pip = default_string

    end function


end module
