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

    use, intrinsic :: iso_fortran_env

    implicit none
    public         ! public for performance in accessing components
    save

    integer, parameter :: dp = REAL64

    ! Various useful constants
    real(dp), parameter          :: sqrt2    = sqrt(2._dp)
    real(dp), parameter          :: isqrt2   = 1/sqrt2
    real(dp), parameter          :: sqrt3    = sqrt(3._dp)
    real(dp), parameter          :: pi       = acos(-1._dp)
    real(dp), parameter          :: kB       = 8.61733238496e-5_dp
    real(dp), parameter          :: hbar     = 0.6582119280967_dp
    real(dp), parameter          :: twelfth  = 1.0_dp/12.0_dp
    integer, parameter           :: randseed(13)   = [7,5,3,11,9,1,17,2,9,6,4,5,8]

    integer, parameter :: verlet_prop = 1001
    integer, parameter :: beeman_prop = 1002
    integer, parameter :: langevin_prop = 1003
    integer, parameter :: langevin_series_prop = 1004

    integer, parameter :: max_string_length = 1000

    ! Conversion constants to program units
    !
    ! Program basic units
    !           Length : Ang
    !           Time   : fs
    !           Energy : eV
    ! Program derived units
    !           Mass   : eV fs^2 / A^2 = 1/103.6382 amu
    !           Angle  : radian = 180 deg
    !           bohr   : bohr = 0.5291772 Angstroem
    real(dp), parameter          :: amu2mass = 103.638239276_dp
    real(dp), parameter          :: deg2rad  = pi/180.0_dp
    real(dp), parameter          :: bohr2ang = 0.529177211_dp
    real(dp), parameter          :: p2GPa    = 160.2176565_dp

    !  Type atoms
    !   structure to hold the position, velocity, force etc. for multiple atoms
    !       use rank2 array so positions, velocities, forces etc. are stored
    !       in sequentional memory locations for efficient access
    !       each array should be allocated (3, n_beads, n_atom)
    !       mass array has length of n_atoms
    type atoms

        integer                                     :: natoms   ! number of atoms
        integer                                     :: nbeads   ! number of beads per atom
        real(dp), dimension(:),         allocatable :: m        ! mass
        integer,  dimension(:),         allocatable :: atn      ! atomic number
        character(len=3), dimension(:), allocatable :: name     ! atomic name
        real(dp), dimension(:,:,:),     allocatable :: r        ! positions
        real(dp), dimension(:,:,:),     allocatable :: v        ! velocities
        real(dp), dimension(:,:,:),     allocatable :: f        ! forces
        logical,  dimension(:,:,:),     allocatable :: fixed    ! mask array defining frozen atoms
                                                                !  T is frozen

    end type atoms

    type simulation_parameters

        integer :: start                                                    ! a trajectory to start with
        integer :: ntrajs                                                   ! number of trajectories
        integer :: nsteps                                                   ! number of steps
        real(dp):: step                                                     ! time step in fs
        integer :: nlattices                                                ! number of lattice species
        integer :: nprojectiles                                             ! number of incident species
        character(len=3), dimension(:), allocatable :: name_l, name_p       ! atomic names
        real(dp),         dimension(:), allocatable :: mass_l, mass_p       ! atomic masses
        integer,          dimension(:), allocatable :: md_algo_l, md_algo_p ! and respective key
        real(dp),         dimension(:), allocatable :: einc                 ! incidence energy (eV)
        real(dp),         dimension(:), allocatable :: inclination          ! incidence polar angle (degree)
        real(dp),         dimension(:), allocatable :: azimuth              ! incidence azimuthal angle (degree)
        real(dp):: Tsurf                                                    ! surface temperature in K
        real(dp):: sa_Tmax                                                  ! max. Tsurf for simulated annealing in K
        integer :: sa_nsteps                                                ! number of steps per simulated annealing cycle
        integer :: sa_interval                                              ! number of steps per temperature interval
        character(len=7)    :: confname                                     ! configuration key
        character(len=max_string_length) :: confname_file                   ! name of the system configuration file
        integer, dimension(2) :: rep                                        ! defines in-plane repetitions
        integer :: nconfs                                                   ! Number of configurations to read in


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
        allocate(new_atoms%fixed(3,nbeads,natoms))

        new_atoms%nbeads = nbeads
        new_atoms%natoms = natoms
        new_atoms%m     = -1.0_dp
        new_atoms%atn   = 0
        new_atoms%name  = ""
        new_atoms%r     = 0.0_dp
        new_atoms%v     = 0.0_dp
        new_atoms%f     = 0.0_dp
        new_atoms%fixed = .FALSE.


    end function

    ! Constructor for type simulation_parameters

    function new_simulation_parameters()

        type(simulation_parameters) new_simulation_parameters

        new_simulation_parameters%start  = -1
        new_simulation_parameters%ntrajs = -1
        new_simulation_parameters%nsteps = -1
        new_simulation_parameters%step   = -1_dp
        new_simulation_parameters%nlattices = -1
        new_simulation_parameters%nprojectiles = -1
        new_simulation_parameters%Tsurf   = -1_dp
        new_simulation_parameters%sa_Tmax   = -1_dp
        new_simulation_parameters%sa_nsteps = -1
        new_simulation_parameters%sa_interval = -1
        new_simulation_parameters%confname = ""
        new_simulation_parameters%rep = [0,0]
        new_simulation_parameters%nconfs  = -1

    end function


end module
