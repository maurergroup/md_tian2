!############################################################################
! This routine is part of
! md_tian2 (Molecular Dynamics Tian Xia 2)
! (c) 2014-2020 Dan J. Auerbach, Svenja M. Janke, Marvin Kammler,
!               Sascha Kandratsenka, Sebastian Wille
! Dynamics at Surfaces Department
! MPI for Biophysical Chemistry Goettingen, Germany
! Georg-August-Universitaet Goettingen, Germany
!
! This program is free software: you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by the
! Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
! or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
! for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program. If not, see http://www.gnu.org/licenses.
!############################################################################

module constants

    use, intrinsic :: iso_fortran_env

    implicit none

    integer, parameter :: dp = REAL64

    ! Physical and mathematical constants
    real(dp), parameter :: sqrt2          = sqrt(2.0_dp)
    real(dp), parameter :: isqrt2         = 1.0_dp/sqrt2
    real(dp), parameter :: sqrt3          = sqrt(3.0_dp)
    real(dp), parameter :: pi             = acos(-1.0_dp)
    real(dp), parameter :: kB             = 8.61733238496e-5_dp       ! eV / K
    real(dp), parameter :: hbar           = 0.6582119514467406_dp     ! eV * fs
    real(dp), parameter :: twelfth        = 1.0_dp/12.0_dp
    integer,  parameter :: dimensionality = 3

    ! file units
    !integer, parameter  :: inp_unit = 67
    !integer, parameter  :: inp_unit = 38
    !integer, parameter  :: pes_unit = 38
    !integer, parameter  :: geo_unit = 55
    !integer, parameter  :: out_unit = 86
    !integer, parameter  :: inpnn_unit       = 61
    !integer, parameter  :: scaling_unit     = 62
    !integer, parameter  :: scalinge_unit     = 63
    !integer, parameter  :: weight_unit      = 64
    !integer, parameter  :: weighte_unit      = 65


    ! PES-related constants
    integer, parameter :: nparams_lj            = 2
    integer, parameter :: nparams_morse         = 3
    integer, parameter :: nparams_emt           = 7
    integer, parameter :: nparams_rebo          = 13
    integer, parameter :: nparams_nene          = 1

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
    integer, parameter :: randseed(13)      = [7,5,3,11,9,1,17,2,9,6,4,5,8]
    integer, parameter :: max_string_length = 1000
    real(dp), parameter :: tolerance        = 1.0e-9_dp


    ! Defaults
    integer,   parameter  :: default_int           = huge(0_4)
    real(dp),  parameter  :: default_real          = huge(0_dp)
    character, parameter  :: default_string        = ""
    logical,   parameter  :: default_bool          = .false.


    ! Propagation
    integer, parameter :: prop_id_verlet          = 1001
    integer, parameter :: prop_id_beeman          = 1002
    integer, parameter :: prop_id_langevin        = 1003
    integer, parameter :: prop_id_langevin_series = 1004
    integer, parameter :: prop_id_andersen        = 1005
    integer, parameter :: prop_id_pile            = 1006

    integer, parameter :: energy_and_force        = 3001
    integer, parameter :: energy_only             = 3002


    ! Geometry optimization
    integer, parameter :: geometry_opt_fire       = 4001


    ! Output
    character(len=*), parameter :: output_key_xyz       = "xyz"
    character(len=*), parameter :: output_key_energy    = "energy"
    character(len=*), parameter :: output_key_poscar    = "poscar"
    character(len=*), parameter :: output_key_vasp      = "vasp"
    character(len=*), parameter :: output_key_mxt       = "mxt"
    character(len=*), parameter :: output_key_scatter   = "scatter"

    integer, parameter :: output_id_xyz     = 1
    integer, parameter :: output_id_energy  = 2
    integer, parameter :: output_id_poscar  = 3
    integer, parameter :: output_id_vasp    = 4
    integer, parameter :: output_id_mxt     = 5
    integer, parameter :: output_id_scatter = 6

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
    real(dp), parameter :: amu2mass             = 103.638239276_dp
    real(dp), parameter :: deg2rad              = pi/180.0_dp
    real(dp), parameter :: rad2deg              = 180.0_dp/pi
    real(dp), parameter :: bohr2ang             = 0.529177211_dp
    real(dp), parameter :: ang2bohr             = 1/bohr2ang
    real(dp), parameter :: p2GPa                = 160.2176565_dp
    real(dp), parameter :: joule2ev             = 6.2415093433e+18_dp
    real(dp), parameter :: kelvin2ev            = kB
    real(dp), parameter :: ha2ev                = 27.21138602_dp ! convert Ha to eV
    real(dp), parameter :: habohr2evang         = ha2ev*ang2bohr !51.4220670398_dp; convert Ha/bohr to eV/ang
    real(dp), parameter :: habohrcub2evangcub   = ha2ev*(ang2bohr**3)
    real(dp), parameter :: au2gpa               = 29419.844d0 ! 1 Ha/Bohr3 = 29419.844 GPa

end module constants
