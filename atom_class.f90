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

        integer                                     :: n_atoms  ! number of atoms
        integer                                     :: n_beads  ! number of beads per atom
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

        integer :: start
        integer :: ntrajs
        integer :: nsteps
        real(dp):: step
        integer :: nlattices       ! number of lattice species
        integer :: nprojectiles    ! number of incident species
        character(len=3), dimension(:), allocatable :: name_l, name_p ! atomic names
        real(dp),         dimension(:), allocatable :: mass_l, mass_p ! atomic masses
        integer,          dimension(:), allocatable :: md_algo_l, md_algo_p     ! and respective key
        real(dp) :: einc          ! incidence energy (eV)
        real(dp) :: inclination   ! incidence polar angle (degree)
        real(dp) :: azimuth       ! incidence azimuthal angle (degree)


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

    function new_atoms(n_beads, n_atoms)
        integer, intent(in) :: n_beads, n_atoms
        type(atoms) new_atoms

        allocate(new_atoms%m(n_atoms))
        allocate(new_atoms%atn(n_atoms))
        allocate(new_atoms%name(n_atoms))
        allocate(new_atoms%r(3,n_beads,n_atoms))
        allocate(new_atoms%v(3,n_beads,n_atoms))
        allocate(new_atoms%f(3,n_beads,n_atoms))
        allocate(new_atoms%fixed(3,n_beads,n_atoms))

        new_atoms%n_beads = n_beads       !   initialize
        new_atoms%n_atoms = n_atoms
        new_atoms%m     = -1.0_dp
        new_atoms%atn   = 0
        new_atoms%name  = ""
        new_atoms%r     = 0.0_dp
        new_atoms%v     = 0.0_dp
        new_atoms%f     = 0.0_dp
        new_atoms%fixed = .FALSE.


    end function

    ! Constructor for type simulation_parameters

    function new_simulation_parameters

        type(simulation_parameters) new_simulation_parameters

        new_simulation_parameters%start  = -1       ! a trajectory to start with
        new_simulation_parameters%ntrajs = -1       ! number of trajectories
        new_simulation_parameters%nsteps = -1       ! number of steps
        new_simulation_parameters%step   = -1_dp    ! time step in fs
        new_simulation_parameters%nlattices = -1    ! number of lattice species
        new_simulation_parameters%nprojectiles = -1 ! number of incident species
        new_simulation_parameters%einc   = -1_dp    ! incidence energy (eV)
        new_simulation_parameters%inclination =-1_dp! incidence polar angle (degree)
        new_simulation_parameters%azimuth   = -1_dp ! incidence azimuthal angle (degree)

    end function


end module
