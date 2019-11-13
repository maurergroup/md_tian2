!############################################################################
! This routine is part of
! md_tian2 (Molecular Dynamics Tian Xia 2)
! (c) 2014-2019 Dan J. Auerbach, Svenja M. Janke, Marvin Kammler,
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

subroutine set_defaults()

    use constants, only : default_int, default_real, default_string, default_bool ! to get default values

    ! RuNNer related modules (predictionshortatomic.f90)
    use fileunits
    use globaloptions
    use mpi_mod
    use nnflags
    use nnshort_atomic
    use predictionoptions
    use saturation
    use symfunctions
    use timings

    ! RuNNer related modules (initnn.f90)
    use fittingoptions
    use mode1options
    use nnewald
    use nnconstants

    implicit none

    nn_type_short                   = default_int
    mode                            = default_int
    lshort                          = default_bool
    lelec                           = default_bool
    nn_type_elec                    = default_int
    lfounddebug                     = default_bool
    ldebug                          = default_bool
    lfound_num_layersshort          = default_bool
    maxnum_layers_short_atomic      = default_int
    lfound_num_layersewald          = default_bool
    maxnum_layers_elec              = default_int
    lfound_luseatomenergies         = default_bool
    luseatomenergies                = default_bool
    lfound_luseatomcharges          = default_bool
    luseatomcharges                 = default_bool
    lfound_nelem                    = default_bool
    nelem                           = default_int
    npairs                          = default_int
    max_num_pairs                   = default_
    element                         = default_string
    nucelem                         = default_int
    dmin_element                    = default_real
    nodes_short_local               = default_int
    nodes_ewald_local               = default_int
    num_funcvalues_local            = default_int
    num_funcvaluese_local           = default_int
    num_funcvaluesp_local           = default_int
    elementtemp                     = default_string
    ztemp                           = default_int
    maxnum_funcvalues_short_atomic  = default_int
    maxnum_funcvalues_elec          = default_int
    function_type_local             = default_int
    function_type_temp              = default_int
    funccutoff_local                = default_real
    maxcutoff_local                 = default_real
    elementtemp1                    = default_string
    elementtemp2                    = default_string
    elementtemp3                    = default_string





































    ldebug = default_bool
    !maxnum_layers_short_atomic = default_int
    luseatomenergies = default_bool
    luseatomcharges = default_bool




    ! dummy default
    default_int
    default_string
    default_bool
    default_real

end subroutine set_defaults()
