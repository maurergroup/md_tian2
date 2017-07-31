module constants

    use, intrinsic :: iso_fortran_env

    implicit none

    public         ! public for performance in accessing components

    integer, parameter :: dp = REAL64

    ! Physical and mathematical constants
    real(dp), parameter :: sqrt2    = sqrt(2._dp)
    real(dp), parameter :: isqrt2   = 1/sqrt2
    real(dp), parameter :: sqrt3    = sqrt(3._dp)
    real(dp), parameter :: pi       = acos(-1._dp)
    real(dp), parameter :: kB       = 8.61733238496e-5_dp
    real(dp), parameter :: hbar     = 0.6582119280967_dp
    real(dp), parameter :: twelfth  = 1.0_dp/12.0_dp


    ! Internal program constants
    integer, parameter :: randseed(13)   = [7,5,3,11,9,1,17,2,9,6,4,5,8]
    integer, parameter :: max_string_length = 1000


    ! Defaults
    integer, parameter :: default_int = huge(0_4)
    real(dp), parameter :: default_real = huge(0_dp)
    character, parameter :: default_string = ""


    ! Propagation
    integer, parameter :: prop_verlet          = 1001
    integer, parameter :: prop_beeman          = 1002
    integer, parameter :: prop_langevin        = 1003
    integer, parameter :: prop_langevin_series = 1004


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
    real(dp), parameter :: amu2mass = 103.638239276_dp
    real(dp), parameter :: deg2rad  = pi/180.0_dp
    real(dp), parameter :: bohr2ang = 0.529177211_dp
    real(dp), parameter :: p2GPa    = 160.2176565_dp



end module constants
