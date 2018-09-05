module constants

    use, intrinsic :: iso_fortran_env

    implicit none

    integer, parameter :: dp = REAL64

    ! Physical and mathematical constants
    real(dp), parameter :: sqrt2    = sqrt(2.0_dp)
    real(dp), parameter :: isqrt2   = 1.0_dp/sqrt2
    real(dp), parameter :: sqrt3    = sqrt(3.0_dp)
    real(dp), parameter :: pi       = acos(-1.0_dp)
    real(dp), parameter :: kB       = 8.61733238496e-5_dp       ! eV / K
    real(dp), parameter :: hbar     = 0.6582119514467406_dp     ! eV * fs
    real(dp), parameter :: twelfth  = 1.0_dp/12.0_dp
    integer,  parameter :: dimensionality = 3


    ! PES-related constants
    integer, parameter :: nparams_lj    = 2
    integer, parameter :: nparams_morse = 3
    integer, parameter :: nparams_emt   = 7
    integer, parameter :: nparams_rebo  = 13
    integer, parameter :: nparams_nene  = 1

    integer, parameter :: pes_id_lj             = 2001
    integer, parameter :: pes_id_morse          = 2002
    integer, parameter :: pes_id_emt            = 2003
    integer, parameter :: pes_id_rebo           = 2004
    integer, parameter :: pes_id_ho             = 2005
    integer, parameter :: pes_id_simple_lj      = 2006
    integer, parameter :: pes_id_no_interaction = 2007
    integer, parameter :: pes_id_nene           = 2008

    ! PES names
    character(len=*), parameter :: pes_name_lj             = "lj"
    character(len=*), parameter :: pes_name_simple_lj      = "slj"
    character(len=*), parameter :: pes_name_emt            = "emt"
    character(len=*), parameter :: pes_name_rebo           = "rebo"
    character(len=*), parameter :: pes_name_no_interaction = "non"
    character(len=*), parameter :: pes_name_ho             = "ho"
    character(len=*), parameter :: pes_name_morse          = "morse"
    character(len=*), parameter :: pes_name_nene           = "nene"

     ! Internal program constants
    integer, parameter :: randseed(13) = [7,5,3,11,9,1,17,2,9,6,4,5,8]
    integer, parameter :: max_string_length = 1000
    real(dp), parameter :: tolerance = 1.0e-9_dp


    ! Defaults
    integer, parameter   :: default_int    = huge(0_4)
    real(dp), parameter  :: default_real   = huge(0_dp)
    character, parameter :: default_string = ""


    ! Propagation
    integer, parameter :: prop_id_verlet          = 1001
    integer, parameter :: prop_id_beeman          = 1002
    integer, parameter :: prop_id_langevin        = 1003
    integer, parameter :: prop_id_langevin_series = 1004
    integer, parameter :: prop_id_andersen        = 1005
    integer, parameter :: prop_id_pile            = 1006

    integer, parameter :: energy_and_force = 3001
    integer, parameter :: energy_only      = 3002


    ! Geometry optimization
    integer, parameter :: geometry_opt_fire = 4001


    ! Output
    character(len=*), parameter :: output_key_xyz     = "xyz"
    character(len=*), parameter :: output_key_energy  = "energy"
    character(len=*), parameter :: output_key_poscar  = "poscar"
    character(len=*), parameter :: output_key_mxt     = "mxt"
    character(len=*), parameter :: output_key_scatter = "scatter"

    integer, parameter :: output_id_xyz     = 1
    integer, parameter :: output_id_energy  = 2
    integer, parameter :: output_id_poscar  = 3
    integer, parameter :: output_id_mxt     = 4
    integer, parameter :: output_id_scatter = 5

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
    real(dp), parameter :: amu2mass  = 103.638239276_dp
    real(dp), parameter :: deg2rad   = pi/180.0_dp
    real(dp), parameter :: rad2deg   = 180.0_dp/pi
    real(dp), parameter :: bohr2ang  = 0.529177211_dp
    real(dp), parameter :: ang2bohr  = 1/bohr2ang
    real(dp), parameter :: p2GPa     = 160.2176565_dp
    real(dp), parameter :: joule2ev  = 6.2415093433e+18_dp
    real(dp), parameter :: kelvin2ev = kB
    real(dp), parameter :: ha2ev     = 27.21138602_dp
    real(dp), parameter :: forceconv = ha2ev/bohr2ang !51.4220670398_dp; convert Ha/bohr to eV/ang

    !arrays used in nene
    character(len=93)   :: atom_format = '(A4, X, F12.9, X, F12.9, X, F12.9, X, A3, X, F11.8, X, F11.8, X, F11.8, X, F11.8, X, F11.8)'
    character(len=20)   :: ce_format = '(A6 ,X , F11.8)'
    character(len=36)   :: lattice_format = '(A7, X, F11.8, X, F11.8, X, F11.8)'
    character(len=29)   :: forces_string = '(F11.8, X, F11.8, X, F11.8)'

end module constants
