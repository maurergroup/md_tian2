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

! maybe use this module for more stuff related to RuNNer which has to be present in the pes_nene_mod module
! move the subroutines back into pes_nene_mod.f90, so that only one file is needed!
module pes_nene_mod_supply

    use constants, only : default_int, default_real, default_string, default_bool, max_string_length, dp

    ! RuNNer related modules (predict.f90)
    use fileunits
    use fittingoptions
    use globaloptions
    use inputnncounters
    use mode1options
    use mpi_mod
    use nnconstants
    use nnewald
    use nnflags
    use nnshort_atomic
    !use nnshort_pair ! this module is not needed, Pair NN not implemented!!
    use predictionoptions
    !use saturation
    !use structures
    use symfunctions
    use timings

    implicit none

    integer :: ielem
    integer :: iseed

    ! following all needed variable declarations not listed in any RuNNer related module
    logical :: lelement(102)
    integer :: ztemp

    integer, dimension(:)  , allocatable :: nodes_short_local
    integer, dimension(:)  , allocatable :: nodes_ewald_local
    integer, dimension(:)  , allocatable :: num_funcvalues_local
    integer, dimension(:)  , allocatable :: num_funcvaluese_local

    character*2 :: elementtemp
    character*2 :: elementtemp1
    character*2 :: elementtemp2
    character*2 :: elementtemp3

    logical :: lfounddebug
    logical :: lfound_num_layersshort
    logical :: lfound_num_layersewald
    logical :: lfound_num_layerspair
    logical :: lfound_luseatomenergies
    logical :: lfound_luseatomcharges
    logical :: lfound_nelem
    integer :: function_type_temp
    real(dp) :: funccutoff_local
    real(dp) :: maxcutoff_local

    integer num_atoms
    real*8  lattice(3,3)
    real*8  xyzstruct(3,max_num_atoms)
    integer zelem(max_num_atoms)
    integer num_pairs


    contains

    subroutine set_defaults()

        iseed                               = default_int
        ielem                               = default_int
        lelement                            = default_bool

        nn_type_short                       = default_int
        mode                                = default_int
        lshort                              = default_bool
        lelec                               = default_bool
        nn_type_elec                        = default_int
        lfounddebug                         = default_bool
        ldebug                              = default_bool
        lfound_num_layersshort              = default_bool
        maxnum_layers_short_atomic          = default_int
        lfound_num_layersewald              = default_bool
        maxnum_layers_elec                  = default_int
        lfound_luseatomenergies             = default_bool
        luseatomenergies                    = default_bool
        lfound_luseatomcharges              = default_bool
        luseatomcharges                     = default_bool
        lfound_nelem                        = default_bool
        nelem                               = default_int
        npairs                              = default_int
        max_num_pairs                       = default_int
        element                             = default_string
        nucelem                             = default_int
        dmin_element                        = default_real
        !nodes_short_local                   = default_int
        !nodes_ewald_local                   = default_int
        !num_funcvalues_local                = 0 ! needed so that the max function will work
        !num_funcvaluese_local               = 0 ! needed so that the max function will work
        elementtemp                         = default_string
        ztemp                               = default_int
        maxnum_funcvalues_short_atomic      = 0 ! needed so that the max function will work
        maxnum_funcvalues_elec              = 0 ! needed so that the max function will work
        function_type_local                 = default_int
        function_type_temp                  = default_int
        funccutoff_local                    = 0.0d0 ! needed so that the max function will work
        maxcutoff_local                     = 0.0d0 ! needed so that the max function will work
        elementtemp1                        = default_string
        elementtemp2                        = default_string
        elementtemp3                        = default_string

        maxnodes_short_atomic               = 0 ! needed so that the max function will work
        maxnodes_elec                       = 0 ! needed so that the max function will work



































        pstring = default_string

        ldebug = default_bool
        !maxnum_layers_short_atomic = default_int
        luseatomenergies = default_bool
        luseatomcharges = default_bool


    end subroutine set_defaults

    subroutine inputnndefaults()

        implicit none

        if(lshort.and.(nn_type_short.eq.1))then
            nodes_short_atomic(:,:)=0
        endif

        if(lelec.and.(nn_type_elec.eq.1))then
            nodes_elec(:,:)=0
        endif

        if (analyze_error_energy_step == default_real) then
            analyze_error_energy_step = 0.01d0
        end if
        if (analyze_error_force_step == default_real) then
            analyze_error_force_step = 0.01d0
        end if
        if (analyze_error_charge_step == default_real) then
            analyze_error_charge_step = 0.001d0
        end if
        if (paramode == default_int) then
            paramode = 1
        end if
        if (ewaldalpha == default_real) then
            ewaldalpha = 0.0d0
        end if
        if (ewaldcutoff == default_real) then
            ewaldcutoff = 0.0d0
        end if
        if (ewaldkmax == default_int) then
            ewaldkmax = 0
        end if
        if (nenergygroup == default_int) then
            nenergygroup = 1
        end if
        if (nforcegroup == default_int) then
            nforcegroup = 1
        end if
        if (nchargegroup == default_int) then
            nchargegroup = 1
        end if
        if (energyrnd == default_real) then
            energyrnd = 1.0d0
        end if
        if (forcernd == default_real) then
            forcernd = 1.0d0
        end if
        if (chargernd == default_real) then
            chargernd = 1.0d0
        end if
        if (fitethres == default_real) then
            fitethres = 0.0d0
        end if
        if (fitfthres == default_real) then
            fitfthres = 0.0d0
        end if
        if (rmin == default_real) then
            rmin = 0.5d0
        end if
        if (optmodee == default_int) then
            optmodee = 1
        end if
        if (optmodef == default_int) then
            optmodef = 1
        end if
        if (optmodeq == default_int) then
            optmodeq = 1
        end if
        if (nblock == default_int) then
            nblock = 200
        end if
        if (nepochs == default_int) then
            nepochs = 0
        end if
        if (iwriteweight == default_int) then
            iwriteweight = 1
        end if
        if (kalmanthreshold == default_real) then
            kalmanthreshold = 0.0d0
        end if
        if (kalmanthresholdf == default_real) then
            kalmanthresholdf = 0.0d0
        end if
        if (kalmanthresholde == default_real) then
            kalmanthresholde = 0.0d0
        end if
        if (kalmanthresholdc == default_real) then
            kalmanthresholdc = 0.0d0
        end if
        if (kalman_dampe == default_real) then
            kalman_dampe = 1.0d0
        end if
        if (kalman_dampf == default_real) then
            kalman_dampf = 1.0d0
        end if
        if (kalman_dampq == default_real) then
            kalman_dampq = 1.0d0
        end if
        if (steepeststepe == default_real) then
            steepeststepe = 0.01d0
        end if
        if (steepeststepf == default_real) then
            steepeststepf = 0.01d0
        end if
        if (steepeststepq == default_real) then
            steepeststepq = 0.01d0
        end if
        if (scalefactorf == default_real) then
            scalefactorf = 1.d0
        end if
        if (ngrowth == default_int) then
            ngrowth = 0
        end if
        if (growthstep == default_int) then
            growthstep = 1
        end if
        if (dampw == default_real) then
            dampw = 0.0d0
        end if
        if (all(atomrefenergies == default_real)) then
            atomrefenergies(:) = 0.0d0
        end if
        if (weights_min == default_real) then
            weights_min = -1.d0
        end if
        if (weights_max == default_real) then
            weights_max = 1.d0
        end if
        if (biasweights_min == default_real) then
            biasweights_min = -1.d0
        end if
        if (biasweights_max == default_real) then
            biasweights_max = 1.d0
        end if
        if (weightse_min == default_real) then
            weightse_min = -1.d0
        end if
        if (weightse_max == default_real) then
            weightse_max = 1.d0
        end if
        if (fitting_unit == default_int) then
            fitting_unit = 1
        end if
        if (pstring == default_string) then
            pstring = '00000000000000000000'
        end if
        if (nran == default_int) then
            nran = 5
        end if
        if (all(fixedcharge == default_real)) then
            fixedcharge(:) = 99.0d0
        end if
        if (maxforce == default_real) then
            maxforce = 10000.d0
        end if
        if (maxenergy == default_real) then
            maxenergy = 10000.d0
        end if
        if (restrictw == default_real) then
            restrictw = -100000.d0
        end if
        if (fitmode == default_int) then
            fitmode = 1
        end if
        if (scmin_short_atomic == default_real) then
            scmin_short_atomic = 0.0d0
        end if
        if (scmax_short_atomic == default_real) then
            scmax_short_atomic = 1.0d0
        end if
        if (scmin_elec == default_real) then
            scmin_elec = 0.0d0
        end if
        if (scmax_elec == default_real) then
            scmax_elec = 1.0d0
        end if
        if (noisee == default_real) then
            noisee = 0.0d0
        end if
        if (noisef == default_real) then
            noisef = 0.0d0
        end if
        if (noiseq == default_real) then
            noiseq = 0.0d0
        end if
        if (cutoff_type == default_int) then
            cutoff_type = 1
        end if
        if (cutoff_alpha == default_real) then
            cutoff_alpha = 0.0d0
        end if
        if (rscreen_cut == default_real) then
            rscreen_cut = 0.0d0
        end if
        if (rscreen_onset == default_real) then
            rscreen_onset = 0.0d0
        end if
        if (dynforcegroup_start == default_int) then
            dynforcegroup_start = 20
        end if
        if (dynforcegroup_step == default_int) then
            dynforcegroup_step = 2
        end if
        if (nshuffle_weights_short_atomic == default_int) then
            nshuffle_weights_short_atomic = 10
        end if
        if (shuffle_weights_short_atomic == default_real) then
            shuffle_weights_short_atomic = 0.1d0
        end if
        if (saturation_threshold == default_real) then
            saturation_threshold = 0.99d0
        end if
        if (dataclusteringthreshold1 == default_real) then
            dataclusteringthreshold1 = 1.0d0
        end if
        if (dataclusteringthreshold2 == default_real) then
            dataclusteringthreshold2 = 0.0d0
        end if
        if (inputforcethreshold == default_real) then
            inputforcethreshold = 0.001d0
        end if
        if (kalman_epsilon == default_real) then
            kalman_epsilon = 1.0d0
        end if
        if (kalman_q0 == default_real) then
            kalman_q0 = 0.0d0
        end if
        if (kalman_qtau == default_real) then
            kalman_qtau = 0.0d0
        end if
        if (kalman_qmin == default_real) then
            kalman_qmin = 0.0d0
        end if

    end subroutine inputnndefaults

    subroutine checkinputnn(err_main,err_file)

        implicit none

        integer counter

        character(len=*), parameter, intent(in) :: err_main, err_file
        character(len=*), parameter,            :: err_check = "Error in checkinputnn: "


        if (nran .neq. 5) then
            print *, err_main, err_file, err_check, "random_number_type not implemented, only 5 available"
        end if

      if(lfindcontradictions)then
        if(deltagthres.gt.0.0d0)then
          write(*,*)'ERROR: find_contradictions requires positive deltagthres ',&
            deltagthres
          stop
        endif
        if(deltafthres.gt.0.0d0)then
          write(*,*)'ERROR: find_contradictions requires positive deltafthres ',&
            deltafthres
          stop
        endif
      endif

      if((cutoff_alpha.gt.1.00000001d0).or.(cutoff_alpha.lt.0.0d0))then
        write(*,*)'ERROR: please use cutoff_alpha within 0 and 1 ',cutoff_alpha
        stop
      endif

      if(lusenoisematrix) then
        if((kalman_q0 .le. 0.0d0 ).or.(kalman_qmin.le.0.0d0).or.(kalman_qtau.le.0.0d0)) then
          write(*,*)'ERROR: please define the q0,qmin,qtau for noise matrix ', &
          'and use them larger than zero ',kalman_q0,kalman_qmin,kalman_qtau
          stop
        endif
      endif

      if(lfixweights.and.lshuffle_weights_short_atomic)then
        write(*,*)'ERROR: shuffle_weights_short_atomic cannot be combined with fixed weights'
        stop
      endif

      if(lscreen)then
        if(rscreen_onset.gt.rscreen_cut)then
          write(*,*)'ERROR: rscreen_onset .gt. rscreen_cut in screen_electrostatics'
          stop
        endif
        if(rscreen_onset.lt.0.0d0)then
          write(*,*)'ERROR: rscreen_onset .lt. 0 in screen_electrostatics'
          stop
        endif
        if(rscreen_cut.lt.0.0d0)then
          write(*,*)'ERROR: rscreen_cut .lt. 0 in screen_electrostatics'
          stop
        endif
      endif

      if(noisee.lt.0.0d0)then
        write(*,*)'ERROR: noise_energy must not be negative ',noisee
        stop
      endif

      if(noisef.lt.0.0d0)then
        write(*,*)'ERROR: noise_force must not be negative ',noisef
        stop
      endif

      if(noiseq.lt.0.0d0)then
        write(*,*)'ERROR: noise_charge must not be negative ',noiseq
        stop
      endif

      if(lsysweights.and.lnwweights)then
        write(*,'(a)')'Error: Cannot use systematic_weights_short and nguyen_widrow_weights_short together!'
        stop
      endif

      if(lsysweightse.and.lnwweightse)then
        write(*,'(a)')'Error: Cannot use systematic_weights_ewald and nguyen_widrow_weights_ewald together!'
        stop
      endif

      if(lnormnodes.and.lnwweights)then
        write(*,'(a)')'Error: Cannot use normalize_nodes and nguyen_widrow_weights_short together!'
        stop
      endif

      if(lnormnodes.and.lnwweightse)then
        write(*,'(a)')'Error: Cannot use normalize_nodes and nguyen_widrow_weights_ewald together!'
        stop
      endif

      if((count_kalmanthreshold.eq.1).and.(count_lfixederrore.eq.1))then
        write(*,'(2a)')'Error: short_energy_error_threshold cannot be used ',&
          'in combination with fixed_short_energy_error_threshold'
        stop
      endif

      if((count_kalmanthresholdf.eq.1).and.(count_lfixederrorf.eq.1))then
        write(*,'(a)')'Error: short_force_error_thresholdf cannot be used in combination with fixed_short_force_error_threshold'
        stop
      endif

      if(count_mode.eq.0)then
        write(*,*)'Error: runner_mode is not specified'
        stop
      endif

      if((.not.lshort).and.(.not.lelec))then
        write(*,*)'Error: short range and electrostatic NNs are switched off'
        stop
      endif

      if(lshort.and.(maxnum_layers_short_atomic.eq.0).and.(nn_type_short.eq.1))then
        write(*,*)'Error: global_hidden_layers_short is not specified'
        stop
      endif

      if(lshort.and.(maxnum_layers_short_pair.eq.0).and.(nn_type_short.eq.2))then
        write(*,*)'Error: global_hidden_layers_pair is not specified'
        stop
      endif

      if(lelec.and.(nn_type_elec.eq.0))then
        write(*,*)'Error: electrostatic_type is not specified'
        stop
      endif

      if(lelec.and.(nn_type_elec.eq.1).and.(maxnum_layers_elec.eq.0))then
        write(*,*)'Error: global_hidden_layers_electrostatic is not specified'
        stop
      endif

      if(lshort.and.(count_nodes_short_atomic.eq.0).and.(nn_type_short.eq.1))then
        write(*,*)'Error: global_nodes_short is not specified'
        stop
      endif

      if(lelec.and.(nn_type_elec.eq.1).and.(count_nodes_elec.eq.0))then
        write(*,*)'Error: global_nodes_electrostatic is not specified'
        stop
      endif

      if(lshort.and.(count_nodes_short_pair.eq.0).and.(nn_type_short.eq.2))then
        write(*,*)'Error: global_nodes_pair is not specified'
        stop
      endif

      do i1=1,nelem
        if(lshort.and.(nn_type_short.eq.1))then
          if(nodes_short_atomic(maxnum_layers_short_atomic,i1).gt.1)then
            write(*,*)'Error: More than 1 output node currently does '
            write(*,*)'make sense in short range NN'
            stop
          endif
        endif
      enddo ! i1

      do i1=1,nelem
        if(lshort.and.(nn_type_short.eq.1))then
          if(nodes_short_atomic(maxnum_layers_short_atomic,i1).eq.0)then
            write(*,*)'Error: output_nodes_short is 0'
            stop
          endif
        endif
      enddo ! i1

      do i1=1,nelem
        if(lelec.and.(nn_type_elec.eq.1))then
          if(nodes_elec(maxnum_layers_elec,i1).gt.1)then
            write(*,*)'Error: More than 1 output node currently does '
            write(*,*)'make sense in electrostatic NN'
            stop
          endif
        endif
      enddo ! i1

      do i1=1,nelem
        if(lelec.and.(nn_type_elec.eq.1))then
          if(nodes_elec(maxnum_layers_elec,i1).eq.0)then
            write(*,*)'Error: output_nodes_electrostatic is 0'
            stop
          endif
        endif
      enddo ! i1

      do i1=1,nelem
        if(lshort.and.(nn_type_short.eq.1))then
          if(nodes_short_atomic(0,i1).eq.0)then
            write(*,*)'Error: input_nodes_short is 0'
            stop
          endif
        endif
      enddo ! i1

      do i1=1,nelem
        if(lshort.and.(nn_type_short.eq.1))then
          if(nodes_short_atomic(0,i1).ne.num_funcvalues_short_atomic(i1))then
            write(*,*)'Error: num_funcvalues_short_atomic .ne. nodes_short_atomic(0)',&
              num_funcvalues_short_atomic(i1),nodes_short_atomic(0,i1)
            write(*,*)'Did you set the right number of input nodes?'
            stop
          endif
        endif
      enddo ! i1

      do i1=1,nelem
        if(lelec.and.(nn_type_elec.eq.1))then
          if(nodes_elec(0,i1).eq.0)then
            write(*,*)'Error: input_nodes_electrostatic is 0'
            stop
          endif
        endif
      enddo ! i1

      do i1=1,nelem
        if(lelec.and.(nn_type_elec.eq.1))then
          if(nodes_elec(0,i1).ne.num_funcvalues_elec(i1))then
            write(*,*)'Error: num_funcvalues_elec .ne. nodes_elec(0)',&
              num_funcvalues_elec(i1),nodes_elec(0,i1)
            write(*,*)'Did you set the right number of input nodes?'
            stop
          endif
        endif
      enddo ! i1

      if(lshort.and.(nn_type_short.eq.1))then
        if(count_global_activation_short_atomic.eq.0)then
          write(*,*)'Error: global_activation_short is not specified'
          stop
        endif
      endif

      if(lelec.and.(nn_type_elec.eq.1))then
        if(count_global_activation_elec.eq.0)then
          write(*,*)'Error: global_activation_ewald is not specified'
          stop
        endif
      endif

      if(lelec.and.(count_ewaldalpha.eq.0))then
        write(*,*)'Error: ewald_alpha must be specified for electrostatic NN'
        stop
      endif

      if(lelec.and.(ewaldalpha.le.0))then
        write(*,*)'Error: ewald_alpha must be positive ',ewaldalpha
        stop
      endif
!!
      if(lelec.and.(count_ewaldcutoff.eq.0))then
        write(*,*)'Error: ewald_cutoff must be specified for electrostatic NN'
        stop
      endif
!!
      if(lelec.and.(ewaldcutoff.le.0))then
        write(*,*)'Error: ewald_cutoff must be positive ',ewaldcutoff
        stop
      endif
!!
      if(lelec.and.(count_ewaldkmax.eq.0))then
        write(*,*)'Error: ewald_kmax must be specified for electrostatic NN'
        stop
      endif
!!
      if((.not.lshort).and.(luseforces))then
        write(*,*)'### WARNING ### switching off use_short_forces because no short range NN is used'
        luseforces=.false.
      endif
!!
      if(lelec.and.(.not.luseatomcharges))then
        write(*,*)'### WARNING ### use_atom_charges is switched on for electrostatic NN'
        luseatomcharges=.true.
      endif
!!
      if(lshort.and.(luseatomenergies))then
        write(*,*)'### WARNING ### use_atom_energies is switched off (not implemented)'
        luseatomenergies=.false.
      endif
!!
      if((.not.lshort).and.(lremoveatomenergies))then
        write(*,*)'### WARNING ### remove_atom_energies is switched on without short range NN'
      endif
!!
      if(lelec.and.(lchargeconstraint))then
        write(*,'(a)')' ### WARNING ### use_charge_constraint is not maintained at the moment and might fail'
      endif
!!
      if(count_iseed.eq.0)then
        write(*,*)'### WARNING ### no random_seed specified, using default '
      endif
!!
      if(nenergygroup.gt.nblock)then
        nenergygroup=nblock
        write(*,*)'### WARNING ### reducing nenergygroup to nblock'
      endif

      if(count_nelem.eq.0)then
        write(*,*)'Error: number_of_elements not specified'
        stop
      endif
!!
      if(count_element.eq.0)then
        write(*,*)'Error: elements not specified'
        stop
      endif
!!
      if((mode.eq.1).and.(count_splitthres.eq.0))then
        write(*,*)'Error: test_fraction not specified'
        stop
      endif
!!
      if(lcentersym.and.(count_scmin_short_atomic.gt.0))then
        write(*,'(a)')'Error: center_symmetry_functions cannot be combined with scale_min_short_atomic keyword'
        stop
      endif
!!
      if(lcentersym.and.(count_scmax_short_atomic.gt.0))then
        write(*,'(a)')'Error: center_symmetry_functions cannot be combined with scale_max_short_atomic keyword'
        stop
      endif
!!
      if(lcentersym.and.(count_scmin_short_pair.gt.0))then
        write(*,'(a)')'Error: center_symmetry_functions cannot be combined with scale_min_short_pair keyword'
        stop
      endif
!!
      if(lcentersym.and.(count_scmax_short_pair.gt.0))then
        write(*,'(a)')'Error: center_symmetry_functions cannot be combined with scale_max_short_pair keyword'
        stop
      endif
!!
      if(lcentersym.and.(count_scmin_elec.gt.0))then
        write(*,'(a)')'Error: center_symmetry_functions cannot be combined with scale_min_elec keyword'
        stop
      endif
!!
      if(lcentersym.and.(count_scmax_elec.gt.0))then
        write(*,'(a)')'Error: center_symmetry_functions cannot be combined with scale_max_elec keyword'
        stop
      endif
!!
      if((count_scmin_short_atomic.gt.0).and.(.not.lscalesym))then
        write(*,*)'Error: scale_min_short requires keyword scale_symmetry_functions'
        stop
      endif
!!
      if((count_scmax_short_atomic.gt.0).and.(.not.lscalesym))then
        write(*,*)'Error: scale_max_short requires keyword scale_symmetry_functions'
        stop
      endif
!!
      if((count_scmin_short_pair.gt.0).and.(.not.lscalesym))then
        write(*,*)'Error: scale_min_short_pair requires keyword scale_symmetry_functions'
        stop
      endif
!!
      if((count_scmax_short_pair.gt.0).and.(.not.lscalesym))then
        write(*,*)'Error: scale_max_short_pair requires keyword scale_symmetry_functions'
        stop
      endif
!!
      if((count_scmin_elec.gt.0).and.(.not.lscalesym))then
        write(*,*)'Error: scale_min_elec requires keyword scale_symmetry_functions'
        stop
      endif
!!
      if((count_scmax_elec.gt.0).and.(.not.lscalesym))then
        write(*,*)'Error: scale_max_elec requires keyword scale_symmetry_functions'
        stop
      endif
!!
      if(scmin_short_atomic.ge.scmax_short_atomic)then
        write(*,'(a)')'Error: scale_min_short .ge. scale_max_short'
        stop
      endif
!!
      if(scmin_short_pair.ge.scmax_short_pair)then
        write(*,'(a)')'Error: scale_min_short_pair .ge. scale_max_short_pair'
        stop
      endif
!!
      if(scmin_elec.ge.scmax_elec)then
        write(*,'(a)')'Error: scale_min_elec .ge. scale_max_elec'
        stop
      endif
!!
      if(lupdatebyelement.and.lchargeconstraint)then
        lchargeconstraint=.false.
        if(mode.eq.2)then
          write(*,*)'### WARNING ### lchargeconstraint is switched off because of lupdatebyelement'
        endif
      endif
!!
      if(lshort.and.lupdatebyelement.and.(mode.eq.2))then
        write(*,*)'### WARNING ### lupdatebyelement works only for charges and forces'
      endif
!!
      if(luseworste.and.lshort.and.(energyrnd.lt.1.0d0))then
        energyrnd=1.0d0
        write(*,*)'### WARNING ### luseworste overrides energyrnd: ',energyrnd
      endif
!!
      if(luseworstf.and.lshort.and.luseforces.and.(forcernd.lt.1.0d0))then
        forcernd=1.0d0
        write(*,*)'### WARNING ### luseworstf overrides forcernd: ',forcernd
      endif
!!
      if(luseworstq.and.lelec.and.(nn_type_elec.eq.1).and.(chargernd.lt.1.0d0))then
        chargernd=1.0d0
        write(*,*)'### WARNING ### luseworstq overrides chargernd: ',chargernd
      endif
!!
      if(dampw.gt.1.0d0)then
        write(*,*)'Error: dampw must not be larger than 1.0d0 ',dampw
        stop
      endif
!!
      if(ldostress.and.(.not.ldoforces))then
        write(*,*)'### WARNING ### Analytic stress is requested without forces'
        write(*,*)'Switching on calculation of analytic forces'
        ldoforces=.true.
      endif
!!
      if(ldohessian.and.(.not.ldoforces))then
        write(*,*)'### WARNING ### Analytic Hessian is requested without forces'
        write(*,*)'Switching on calculation of analytic forces'
        ldoforces=.true.
      endif
!!
      if(ldostress.and.(mode.eq.1))then
        write(*,*)'### WARNING ### switching off stress calculation in mode 1 for increased performance'
        ldostress=.false.
      endif
!!
      if((count_wconstraint.gt.0).and.(.not.lfixweights))then
        write(*,*)'Error: weight constraints are specified without fix_weights keyword'
        stop
      endif

      if((count_wconstraint.eq.0).and.(lfixweights))then
        write(*,*)'Error: no weights constrained but keyword fix_weights has been selected'
        stop
      endif
!!
      if(weights_min.ge.weights_max)then
        write(*,*)'Error: weights_min > weights_max'
        stop
      endif
!!
      if(biasweights_min.ge.biasweights_max)then
        write(*,*)'Error: biasweights_min > biasweights_max'
        stop
      endif
!!
      if(weightse_min.ge.weightse_max)then
        write(*,*)'Error: weightse_min > weightse_max'
        stop
      endif
!!
      if(kalman_dampe.lt.0.0d0)then
        write(*,*)'ERROR: kalman_damp_short must be non-negative ',kalman_dampe
        stop
      endif
!!
      if(kalman_dampf.lt.0.0d0)then
        write(*,*)'ERROR: kalman_damp_force must be non-negative ',kalman_dampf
        stop
      endif
!!
      if(kalman_dampq.lt.0.0d0)then
        write(*,*)'ERROR: kalman_damp_charge must be non-negative ',kalman_dampq
        stop
      endif
!!
      if(ljointefupdate.and.lelec)then
        write(*,*)'ERROR: joint_energy_force_update is not implemented for lelec and nn_type_elec 2'
        stop
      endif
!!
      if(ljointefupdate)then
        if(optmodee.ne.optmodef)then
          write(*,*)'Error: joint_energy_force_update requires to use the'
          write(*,*)'same optimization algorithm for energy and forces'
          stop
        endif
        if(.not.luseforces)then
          write(*,*)'Error: switch on use_short_forces for joint_energy_force_update'
          stop
        endif
        if(lrepeate)then
          write(*,*)'ERROR: repeated energy update cannot be combined with joint energy and force update'
          stop
        endif
        if(forcernd.lt.1.0d0)then
          write(*,*)'ERROR: joint energy and force update requires force_fraction = 1.0d0'
          stop
        endif
        if(luseworste)then
          write(*,*)'ERROR: joint energy and force update cannot be combined with update_worst_short_energies'
          stop
        endif
        if(luseworstf)then
          write(*,*)'ERROR: joint energy and force update cannot be combined with update_worst_short_forces'
          stop
        endif
        if(nenergygroup.gt.1)then
          write(*,*)'ERROR: joint energy and force update cannot be combined with short_energy_group > 1'
          stop
        endif
        if(nforcegroup.gt.1)then
          write(*,*)'ERROR: joint energy and force update cannot be combined with short_force_group > 1'
          stop
        endif
        if(kalmanthresholdf.gt.0.0d0)then
          write(*,*)'ERROR: joint energy and force update cannot be combined with short_force_error_threshold > 0.0'
          stop
        endif
      endif
!!
      if(maxforce.le.0.0d0)then
        write(*,*)'Error: max_force must not be negative ',maxforce
        stop
      endif
!!
      if(lshort.and.(nn_type_short.eq.1))then
        do i1=1,nelem
          if(num_funcvalues_short_atomic(i1).ne.nodes_short_atomic(0,i1))then
            write(*,*)'Error: num_funcvalues_short_atomic .ne. nodes_short_atomic(0)'
            write(*,*)i1,num_funcvalues_short_atomic(i1),nodes_short_atomic(0,i1)
            stop
          endif
        enddo! i1
      endif
!!
      if(lelec.and.(nn_type_elec.eq.1))then
        do i1=1,nelem
          if(num_funcvalues_elec(i1).ne.nodes_elec(0,i1))then
            write(*,*)'Error: num_funcvalues_elec .ne. nodes_elec(0)'
            write(*,*)i1,num_funcvalues_elec(i1),nodes_elec(0,i1)
            stop
          endif
        enddo ! i1
      endif
!!
      if((nn_type_elec.eq.4).and.(mode.ne.3))then
        write(*,*)'ERROR: electrostatic_type 4 is only valid for prediction mode'
        stop
      endif
!!
      if((mode.eq.3).and.(max_num_atoms.lt.nblock).and.(nn_type_short.eq.1).and.lshort)then
        write(*,*) 'WARNING: reducing points_in_memory to max_num_atoms ',max_num_atoms
        nblock=max_num_atoms
      endif

    end subroutine checkinputnn

    subroutine printinputnn(iseed,ielem,nodes_short_atomic_temp,nodes_elec_temp,kalmanlambda_local,kalmanlambdae_local,actfunc_short_atomic_dummy,actfunc_elec_dummy)

        use fileunits
        use fittingoptions
        use nnflags
        use globaloptions
        use mode1options
        use predictionoptions
        use nnshort_atomic
        use nnewald
        use inputnncounters

        implicit none

        integer iseed
        integer ielem
        integer i
        integer cnt_1, cnt_2, cnt_3
        integer nodes_short_atomic_temp(0:maxnum_layers_short_atomic)
        integer nodes_elec_temp(0:maxnum_layers_elec)

        real*8 kalmanlambda_local
        real*8 kalmanlambdae_local

        character*1 actfunc_short_atomic_dummy(maxnum_layers_short_atomic)
        character*1 actfunc_elec_dummy(maxnum_layers_elec)

        write(*,*)'General input parameters:'
        write(*,*)'-------------------------------------------------------------'

        if(lshort)then
            write(*,*)'Short range NN is on'
        else
            write(*,*)'Short range NN is off'
        endif

        if(lelec.and.(nn_type_elec.eq.1))then
            write(*,*)'Electrostatic NN is on'
        else
            write(*,*)'Electrostatic NN is off'
        endif

        if(lvdw.and.(nn_type_vdw.eq.1))then
            write(*,*)'vdW corrections switched on'
        else
            write(*,*)'vdW corrections switched off'
        endif

        if((mode.eq.1).and.lcheckinputforces)then
            write(*,'(a,f10.6,a)')&
            ' checking input forces, threshold for force vector is  '&
            ,inputforcethreshold,' Ha/Bohr'
        endif

        write(*,*)'-------------------------------------------------------------'

        if(lshort)then
            if(nn_type_short.le.2)then
                write(*,'(a55,i2)')&
                    ' RuNNer nn_type_short                                              ',nn_type_short
            else
                write(*,*)'ERROR: unknown nn_type_short: ',nn_type_short
                stop
            endif
        endif

        if(mode.eq.1)then
            write(*,*)'RuNNer is started in mode for symmetry function calculation (1)'
        elseif(mode.eq.2)then !'
            write(*,*)'RuNNer is started in mode for fitting (2)'
        elseif(mode.eq.3)then
            write(*,*)'RuNNer is started in mode for prediction (3)'
        else
            write(*,*)'Error: Unknown runner_mode: ',mode
            stop
        endif

        write(*,'(a,l)')' debugging mode is                                       ',ldebug

        write(*,'(a,i4)')' parallelization mode                                  ',paramode

        write(*,'(a,l)')' enable detailed time measurement                        ',lfinetime

        if(mode.eq.2)then
            write(*,'(a,l)')' enable detailed time measurement at epoch level         ',lfinetimeepoch
        endif

        write(*,'(a,l)')' silent mode                                             ',lsilent

        if((mode.eq.2).or.(mode.eq.3))then
            write(*,'(a,l)')' NN force check                                          ',lcheckf
        endif

      if(nelem.lt.ielem)then
        write(*,*)'Error: number of elements in structure(s) is larger than '
        write(*,*)'number of elements in input.nn ',ielem,nelem
        stop
      else
        write(*,'(a,i4)')' number of elements                                    ',nelem
      endif

      write(*,*)'elements (sorted):'
      do i1=1,nelem
        write(*,'(i3,x,a2)')nucelem(i1),element(i1)
      enddo


      write(*,'(a,i10)')' seed for random number generator                ',iseed

      if((nran.lt.0).or.(nran.gt.5))then
        write(*,*)'ERROR: Unknown random number generator ',nran
        stop
      endif
      write(*,'(a,i10)')' random number generator type                    ',nran

      write(*,'(a,l)')' remove free atom reference energies                     ',lremoveatomenergies

      if(lfitethres.and.(mode.eq.1))then
        write(*,'(a,f7.3)')' upper energy threshold per atom (Ha)               ',fitethres
      endif

      if(lfitfthres.and.(mode.eq.1))then
        write(*,'(a,f7.3)')' max force component threshold (Ha/Bohr)            ',fitfthres
      endif

      write(*,'(a,f8.3)')' shortest allowed bond in structure                ',rmin

      if(lnormnodes)then
        write(*,*)'Linear combinations at nodes are normalized'
      endif

      write(*,'(a,i3)')' Cutoff_type for symmetry function is                   ',cutoff_type
      write(*,'(a,f7.3)')' Cutoff_alpha for inner cutoff radius is            ',cutoff_alpha

      if(lenforcemaxnumneighborsatomic)then
        write(*,'(a,i3)')&
        ' Enforcing global max_num_neighors_atomic               ',max_num_neighbors_atomic_input
      endif

      if(lshort.and.(mode.ne.1))then
        write(*,*)'-------------------------------------------------------------'
        write(*,*)'Short range NN specifications:'
        write(*,*)'-------------------------------------------------------------'
      endif

      if(lshort.and.(nn_type_short.eq.1).and.(mode.ne.1))then
        write(*,'(a,10i5)')' global hidden layers short range NN                  ',maxnum_layers_short_atomic-1
        write(*,'(a,10i5)')' global nodes hidden layers short NN             ',&
          (nodes_short_atomic_temp(i1),i1=1,maxnum_layers_short_atomic-1)
      endif

      if(lshort.and.(nn_type_short.eq.1).and.(mode.ne.1))then
        write(*,'(a,x,10a)')' global activation functions short                     ',&
          (actfunc_short_atomic_dummy(i),i=1,maxnum_layers_short_atomic)
      endif

      if(lshort.and.(nn_type_short.eq.2).and.(mode.ne.1))then
        write(*,'(a,10i5)')' global hidden layers short range NN pair             ',maxnum_layers_short_pair-1
        write(*,'(a,10i5)')' global nodes hidden layers short NN pair        ',&
          (nodes_short_pair_temp(i1),i1=1,maxnum_layers_short_pair-1)
      endif

      if(lshort.and.(nn_type_short.eq.2).and.(mode.ne.1))then
        write(*,'(a,x,10a)')' global activation functions short pair                ',&
          (actfunc_short_pair_dummy(i),i=1,maxnum_layers_short_pair)
      endif

      if(lelec)then
        write(*,*)'-------------------------------------------------------------'
        write(*,*)'Electrostatic specifications:'
        write(*,*)'-------------------------------------------------------------'
      endif

      if(lelec)then
        write(*,'(a,i5)')' electrostatic_type (nn_type_elec)                    ',nn_type_elec
        if(nn_type_elec.eq.1)then
          write(*,'(a)')' Using separate set of atomic NNs for atomic charges'
        elseif(nn_type_elec.eq.2)then
          write(*,'(a)')' Constructing atomic charges from short range NN'
        elseif(nn_type_elec.eq.3)then
          write(*,'(a)')' Fixed atomic charges are used:'
          do i1=1,nelem
            write(*,'(a1,a2,x,f14.3)')' ',element(i1),fixedcharge(i1)
          enddo
        elseif(nn_type_elec.eq.4)then
          write(*,'(a)')' Using atomic charges from charges.in file'
        else
          write(*,*)'ERROR: Unknown electrostatic_type ',nn_type_elec
          stop
        endif
      endif

      if(lelec.and.(nn_type_elec.eq.1))then
        write(*,'(a,10i5)')' global hidden layers electrostatic NN                ',maxnum_layers_elec-1
        write(*,'(a,10i5)')' global nodes hidden layers electrostatic NN     ',&
          (nodes_elec_temp(i1),i1=1,maxnum_layers_elec-1)
      endif

      if(lelec.and.(nn_type_elec.eq.1))then
        write(*,'(a,x,10a)')' global activation functions electrostatic             ',&
          (actfunc_elec_dummy(i),i=1,maxnum_layers_elec)
      endif

      if(lelec)then
        write(*,'(a,f8.3)')' Ewald alpha                                       ',ewaldalpha
        write(*,'(a,f8.3)')' Ewald cutoff                                      ',ewaldcutoff
        write(*,'(a,i6)')' Ewald kmax                                          ',ewaldkmax
      endif

      if(lelec.and.(mode.eq.0))then
        write(*,'(a,i4)')' Enforce total charge                                ',enforcetotcharge
      endif

      if(lelec)then
        if(lscreen)then
          write(*,'(a,2f14.6)')' Screening electrostatics                            ',rscreen_onset,rscreen_cut
        else
          write(*,'(a)')' No screening of electrostatics requested                '
        endif
      endif

      if(lelec)then
        if(mode.eq.3)then
          if(nn_type_elec.eq.4)then
            write(*,'(a)')' Using atomic charges from file charges.in!'
          endif
        endif
      endif


      if(lvdw)then
        write(*,*)'-------------------------------------------------------------'
        write(*,*)'vdW specifications:'
        write(*,*)'-------------------------------------------------------------'
      endif

      if(lvdw)then
        write(*,'(a,i4)')' vdW type                                            ',nn_type_vdw
      endif
      if(lvdw.and.(nn_type_vdw.eq.1))then
        write(*,'(a,2f14.6)')' vdw screening                                       ',vdw_screening(1),vdw_screening(2)
      endif







      if(lshort.and.(mode.eq.1))then
        write(*,*)'-------------------------------------------------------------'
        write(*,*)'Parameters for symmetry function generation: short range part:'
        write(*,*)'-------------------------------------------------------------'
      endif

      if(lshort.and.(mode.eq.1)) write(*,'(a,l)')&
        ' using forces for fitting                                ',luseforces

      if(lshort.and.(mode.eq.1)) write(*,'(a,l)')&
        ' using atomic energies for fitting                       ',luseatomenergies

      if(lelec.and.(mode.eq.1))write(*,'(a,l)')&
        ' using atomic charges for fitting                        ',luseatomcharges

      if(mode.eq.1)then
        write(*,*)'-------------------------------------------------------------'
        write(*,'(a,f8.4)') ' percentage of data for testing (%)                ',&
          100.d0*splitthres
      endif

      if(mode.eq.2)then
        write(*,*)'-------------------------------------------------------------'
        write(*,*)'General fitting parameters:'
        write(*,*)'-------------------------------------------------------------'
      endif

      if(mode.eq.2)then
        write(*,'(a,i8)')' number of fitting epochs                          ',nepochs
      endif ! mode.eq.2

      if(mode.eq.2)then
        write(*,'(a,l)')' print date and time for each epoch                      ',lprintdateandtime
      endif ! mode.eq.2

      if((mode.eq.2).and.lenableontheflyinput)then
        write(*,'(a,i8)')' on-the-fly input enabled          '
      endif ! mode.eq.2

      if(mode.eq.2)then
        write(*,'(a,i8)')' number of data sets in memory                     ',nblock
      endif ! mode.eq.2

      if(mode.eq.2)then
        if(fitmode.eq.1)then
          write(*,'(a,i8)')' Fitting mode 1 (online learning) selected         '
        elseif(fitmode.eq.2)then
          write(*,'(a,i8)')' Fitting mode 2 (offline learning) selected        '
        endif
      endif ! mode.eq.2

      if(mode.eq.2)then
        write(*,'(a,l)')' Randomly mixing all points in training set              ',lmixpoints
      endif

      if(mode.eq.2)write(*,'(a,l)')' save Kalman filter data                                 ',lsavekalman

      if(mode.eq.2)write(*,'(a,l)')' restart from old Kalman filter data                     ',lrestkalman

      if(mode.eq.2)write(*,'(a,l)')' rescale symmetry functions                              ',lscalesym

      if((mode.eq.2).and.lscalesym.and.lshort.and.(nn_type_short.eq.1))then
        write(*,'(a,f10.3)')' min value of scaled short range symmetry functions ',scmin_short_atomic
        write(*,'(a,f10.3)')' max value of scaled short range symmetry functions ',scmax_short_atomic
      endif

      if((mode.eq.2).and.lscalesym.and.lshort.and.(nn_type_short.eq.2))then
        write(*,'(a,f10.3)')' min value of scaled pair symmetry functions       ',scmin_short_pair
        write(*,'(a,f10.3)')' max value of scaled pair symmetry functions       ',scmax_short_pair
      endif

      if((mode.eq.2).and.lscalesym.and.lelec.and.(nn_type_elec.eq.1))then
        write(*,'(a,f10.3)')' min value of scaled electrostatic symmetry functions ',scmin_elec
        write(*,'(a,f10.3)')' max value of scaled electrostatic symmetry functions ',scmax_elec
      endif

      if(mode.eq.2)write(*,'(a,l)')&
        ' remove CMS from symmetry functions                      ',lcentersym

      if(mode.eq.2)write(*,'(a,l)')&
        ' calculate symmetry function correlation                 ',lpearson_correlation

      if(mode.eq.2)write(*,'(a,l)')&
        ' weight analysis                                         ',lweightanalysis

      if(mode.eq.2)write(*,'(a,l)')&
        ' environment analysis                                    ',lenvironmentanalysis

      if(mode.eq.2)write(*,'(a,l)')&
        ' find contradictions                                     ',lfindcontradictions
      if((mode.eq.2).and.lfindcontradictions)then
        write(*,'(a,f10.3)')' threshold for |deltaG|                            ',deltagthres
        write(*,'(a,f10.3)')' threshold for delta|F|                            ',deltafthres
      endif

      if(mode.eq.2)write(*,'(a,l)')' fix some weights                                        ',lfixweights

      if(mode.eq.2)write(*,'(a,l)')' using growth mode for fitting                           ',lgrowth
      if((mode.eq.2).and.lgrowth)then
        write(*,'(a,i8)')' number of training structures in each growth step ',ngrowth
      endif
      if(lgrowth.and.(mode.eq.2))then
        write(*,'(a,i4)')' epochs with constant training set size in growth mode ',growthstep
      endif

      if((mode.eq.2).and.ldampw)then
        write(*,'(a,l)')' using weight decay                                      ',ldampw
        write(*,'(a,f18.12)')' balance between error and weight decay  ',dampw
      endif

      if((mode.eq.2).and.lupdatebyelement)then
        write(*,'(a,i3)')' do weight update just for one element                  ',elemupdate
        write(*,*)'### WARNING ### RMSEs will refer only to this element'
      endif

      if(mode.eq.2)write(*,'(a,l)')&
        ' global fit of short and charge NN (not implemented)     ',lglobalfit

      if(mode.eq.2)then
        if(fitting_unit.eq.1)then
          write(*,'(a,a2)')' error unit for fitting                                  ','eV'
        elseif(fitting_unit.eq.2)then
          write(*,'(a,a2)')' error unit for fitting                                  ','Ha'
        else
          write(*,*)'Error: add new energy unit in output of readinput.f90!!!'
          stop
        endif
      endif

      if(mode.eq.2)then
        if(lreadunformatted)then
          write(*,'(a)')' Reading unformatted files '
        else
          write(*,'(a)')' Reading formatted files '
        endif
      endif

      if(mode.eq.2)then
        if(lwriteunformatted)then
          write(*,'(a)')' Writing unformatted files '
        else
          write(*,'(a)')' Writing formatted files '
        endif
      endif

      if(mode.eq.2)then
        if((optmodee.eq.1).or.(optmodef.eq.1).or.(optmodeq.eq.1))then
          write(*,'(a,l)')' Resetting Kalman filter matrices each epoch             ',lresetkalman
        endif
      endif

      if(mode.eq.2)then
        if(nn_type_short.eq.1)then
          if(lshuffle_weights_short_atomic)then
            write(*,'(a,i5,f14.6)')' shuffle_weights_short_atomic                             ',&
              nshuffle_weights_short_atomic,shuffle_weights_short_atomic
          endif
        endif
      endif

      if((mode.eq.2).and.lompmkl)then
        write(*,'(a)')' Using omp mkl for Kalman filter in parallel case'
      endif

      if((mode.eq.2).and.lionforcesonly)then
        write(*,'(a)')' Using only forces for fitting in case of ionic structures'
      endif

      if((mode.eq.2).and.lfitstats)then
        write(*,'(a)')' Writing fitting statistics '
      endif

      if((mode.eq.2).and.(restrictw.gt.0.0d0))then
        write(*,'(a,f14.6)')' Restricting absolute value of weights       ',restrictw
        if((restrictw.gt.0.0d0).and.(restrictw.lt.2.0d0))then
          write(*,*)'Currently restrictw must be larger than 2.0'
          stop
        endif
      endif

      if((mode.eq.2).and.lanalyzeerror)then
        write(*,'(a)')' Error analysis requested for final epoch '
        if(lshort.and.(.not.lwritetrainpoints))then
          write(*,*)'WARNING: trainpoints file is required for short range energy error analysis'
          write(*,*)'=> This analysis will not be done'
        endif
        if(lshort.and.luseforces.and.(.not.lwritetrainforces))then
          write(*,*)'WARNING: trainforces file is required for short range force error analysis'
          write(*,*)'=> This analysis will not be done'
        endif
        if(lelec.and.(.not.lwritetraincharges))then
          write(*,*)'WARNING: traincharges file is required for charge error analysis'
          write(*,*)'=> This analysis will not be done'
        endif
      endif
!!
      if(mode.eq.2.and.((luseoldweightsshort).or.(luseoldweightscharge)))then
        write(*,'(a,l)')' Using old scaling data for restart          ',luseoldscaling
      endif
!!
      if((mode.eq.2).and.(lprecond))then
        write(*,*)'Preconditioning of weights is switched on'
      endif
!!
      if((mode.eq.2).and.(linionly))then
        write(*,*)'Termination of mode 2 after initialization requested'
      endif
!!
      if((mode.eq.2).and.(ldataclustering))then
        write(*,'(a,2f14.10)')'data clustering requested with distance thresholds ',&
          dataclusteringthreshold1,dataclusteringthreshold2
      endif
!!
      if((mode.eq.2).and.(lprintconv))then
        write(*,*)'printing of convergence vector requested'
      endif
!!
      if((mode.eq.2).and.(lanalyzecomposition))then
        write(*,*)'analysis of chemical composition requested'
      endif
!!
      if((mode.eq.2).and.lshort)then
        write(*,*)'-------------------------------------------------------------'
        write(*,*)'Fitting parameters short range part:'
        write(*,*)'-------------------------------------------------------------'
      endif
!!
      if(lshort.and.(mode.eq.2)) write(*,'(a,l)')' using forces for fitting                                ',luseforces
!!
      if((mode.eq.2).and.lshort)then
        if(optmodee.eq.1)then
          write(*,'(a)')' using Kalman filter optimization (1) for short range energy'
          if(luseedkalman)then
            write(*,'(a)')' using element decoupled Kalman filter'
            if(ledforcesv2)then
              write(*,'(a)')' using second variant of ED force fitting'
            endif
          endif
        elseif(optmodee.eq.2)then
          write(*,'(a)')' using conjugate gradient optimization (2) for short range energy'
        elseif(optmodee.eq.3)then
          write(*,'(a)')' using steepest descent optimization (3) for short range energy'
        else
          write(*,*)'Error: Unknown optimization mode ',optmodee
          stop
        endif
      endif ! mode.eq.2
!!
      if((mode.eq.2).and.lshort.and.luseforces)then
        if(optmodef.eq.1)then
          write(*,'(a)')' using Kalman filter optimization (1) for short range forces'
        elseif(optmodef.eq.2)then
          write(*,'(a)')' using conjugate gradient optimization (2) for short range forces'
        elseif(optmodef.eq.3)then
          write(*,'(a)')' using steepest descent optimization (3) for short range forces'
        else
          write(*,*)'Error: Unknown optimization mode ',optmodef
          stop
        endif
      endif ! mode.eq.2
!!
      if((mode.eq.2).and.lshort.and.(.not.lfixederrore))&
        write(*,'(a,f14.8)')' short energy error threshold                ',kalmanthreshold
!!
      if((mode.eq.2).and.lshort.and.(.not.lfixederrorf))&
        write(*,'(a,f14.8)')' short force error threshold                 ',kalmanthresholdf
!!
      if((mode.eq.2).and.lshort.and.lfixederrore)write(*,'(a,f14.8)')&
        ' fixed short energy error threshold          ',fixederrore
!!
      if((mode.eq.2).and.lshort.and.lfixederrorf)&
        write(*,'(a,f14.8)')' fixed short force error threshold           ',fixederrorf
!!
      if(mode.eq.2)then
        if(lshort.and.(nn_type_short.eq.1))kalmanlambda(:)=kalmanlambda_local
        if(lshort.and.(nn_type_short.eq.2))kalmanlambdap(:)=kalmanlambda_local
      endif
!!
      if((mode.eq.2).and.lshort.and.(optmodee.eq.1))&
        write(*,'(a,f14.8)')' Kalman lambda (short)                       ',kalmanlambda_local
!!
      if((mode.eq.2).and.lshort.and.(optmodee.eq.1))&
        write(*,'(a,f14.8)')' Kalman nue (short)                          ',kalmannue
!!
      if((mode.eq.2).and.lshort.and.(optmodee.eq.1))&
        write(*,'(a,l)')' use_noisematrix                                         ',lusenoisematrix  !! modifed by kenko
!!
      if((mode.eq.2).and.lshort.and.(optmodee.eq.1).and.lusenoisematrix)&
        write(*,'(a,f14.8)')' kalman_q0                                   ',kalman_q0  !! modifed by kenko
!!
      if((mode.eq.2).and.lshort.and.(optmodee.eq.1).and.lusenoisematrix)&
        write(*,'(a,f14.8)')' kalman_qmin                                 ',kalman_qmin  !! modifed by kenko
!!
      if((mode.eq.2).and.lshort.and.(optmodee.eq.1).and.lusenoisematrix)&
        write(*,'(a,f14.8)')' kalman_qtau                                 ',kalman_qtau  !! modifed by kenko
!!
      if((mode.eq.2).and.lshort.and.(optmodee.eq.1).and.lusenoisematrix)&
        write(*,'(a,f14.8)')' kalman_epsilon                              ',kalman_epsilon  !! modifed by kenko
!!
      if((mode.eq.2).and.lshort.and.(optmodee.eq.1))&
        write(*,'(a,f14.8)')' Kalman damp (short energy)                  ',kalman_dampe
!!
      if((mode.eq.2).and.lshort.and.(optmodee.eq.1))&
        write(*,'(a,f14.8)')' Kalman damp (short force)                   ',kalman_dampf
!!
      if((mode.eq.2).and.lshort.and.(optmodee.eq.3))then
        write(*,'(a,f14.8)')' steepest descent step size short energy     ',steepeststepe
      endif
!!
      if((mode.eq.2).and.lshort.and.(optmodef.eq.3))then
        write(*,'(a,f14.8)')' steepest descent step size short forces     ',steepeststepf
      endif
!!
      if(mode.eq.2)write(*,'(a,l)')' restart fit with old weights (short)                    ',luseoldweightsshort
!!
      if((mode.eq.2).and.lshort.and.luseworste)&
        write(*,'(a,f8.4)')' fraction of worst short range energies            ',worste
!!
      if((mode.eq.2).and.lshort.and.luseforces.and.luseworstf)&
        write(*,'(a,f8.4)')' fraction of worst short range forces              ',worstf
!!
      if((mode.eq.2).and.luseforces.and.lshort)then
        if(scalefactorf.lt.0.0d0)then
          write(*,'(a)')' automatic scaling factor for force update selected'
        else
          write(*,'(a,f11.8)')' scaling factor for force update (scalefactorf) ',scalefactorf
        endif
      endif
!!
      if(lshort.and.(mode.eq.2))&
        write(*,'(a,i8)')' grouping energies in blocks of                    ',nenergygroup
!!
      if(lshort.and.(mode.eq.2)) &
        write(*,'(a,f8.3)')' fraction of energies used for update              ',energyrnd
!!
      if(lshort.and.(mode.eq.2).and.(.not.lfgroupbystruct).and.(luseforces))then
        write(*,'(a,i8)')' grouping forces in blocks of                      ',nforcegroup
      endif
!!
      if(lshort.and.(mode.eq.2).and.(lfgroupbystruct).and.(luseforces))then
        write(*,'(a,i8)')' automatic grouping forces for update by structure'
      endif
!!
      if(lshort.and.(mode.eq.2))write(*,'(a,f8.3)')' fraction of forces used for update                ',forcernd
!!
      if((mode.eq.2).and.lshort.and.(.not.luseoldweightsshort))&
        write(*,'(a,f14.3)')' weights_min                                 ',weights_min
!!
      if((mode.eq.2).and.lshort.and.(.not.luseoldweightsshort))&
        write(*,'(a,f14.3)')' weights_max                                 ',weights_max
!!
      if((mode.eq.2).and.lshort.and.lseparatebiasini.and.(.not.luseoldweightsshort))&
        write(*,'(a,f14.3)')' biasweights_min                             ',biasweights_min
!!
      if((mode.eq.2).and.lshort.and.lseparatebiasini.and.(.not.luseoldweightsshort))&
        write(*,'(a,f14.3)')' biasweights_max                             ',biasweights_max
!!
      if((mode.eq.2).and.lshort.and.lnwweights)write(*,'(a)')' Using Nguyen Widrow weights for short range NN'
!!
      if((mode.eq.2).and.lshort.and.lsysweights)write(*,'(a)')' Using systematic weights for short range NN'
!!
      if((mode.eq.2).and.lelec.and.(nn_type_elec.eq.1).and.lnwweightse)&
        write(*,'(a)')' Using Nguyen Widrow weights for electrostatic NN'
!!
      if((mode.eq.2).and.lelec.and.(nn_type_elec.eq.1).and.lsysweightse)&
        write(*,'(a)')' Using systematic weights for electrostatic NN'
!!
      if((mode.eq.2).and.lshort.and.luseforces.and.lsepkalman.and.(optmodee.eq.1).and.(optmodef.eq.1))then
        write(*,'(a)')' Using separate Kalman filter matrices for short range energies and forces'
      endif
!!
      if((mode.eq.2).and.lshort.and.luseforces.and.lrepeate)then
        write(*,'(a)')' Using repeated energy updates after each force update'
      endif
!!
      if((mode.eq.2).and.lshort.and.(.not.luseforces).and.lfinalforce)then
        write(*,'(a)')' Calculating force error in final epoch only'
      endif
!!
      if((mode.eq.2).and.(lshort))then
        write(*,'(a,f14.3)')' max_energy                                  ',maxenergy
      endif
!!
      if((mode.eq.2).and.(lshort).and.(luseforces))then
        write(*,'(a,f14.3,a)')' max force component used for fitting        ',&
          maxforce,' Ha/Bohr'
      endif
!!
      if((mode.eq.2).and.lshort)then
        write(*,'(a,f14.8,x,a7)')' noise energy threshold                      ',noisee,'Ha/atom'
      endif
!!
      if((mode.eq.2).and.luseforces.and.lshort)then
        write(*,'(a,f14.8,x,a7)')' noise force threshold                       ',noisef,'Ha/Bohr'
      endif
!!
      if((mode.eq.2).and.ldynforcegroup)then
        write(*,'(a,2i8)')' dynamic force grouping                      ',&
          dynforcegroup_start,dynforcegroup_step
      endif
!!
      if((mode.eq.2).and.ldetect_saturation.and.(lshort).and.(nn_type_short.eq.1))then
        write(*,'(a,f14.6)')' detect saturation of nodes is on            ',&
          saturation_threshold
      endif
!!
      if((mode.eq.2).and.lelec)then
        write(*,*)'-------------------------------------------------------------'
        write(*,*)'Fitting parameters electrostatic part:'
        write(*,*)'-------------------------------------------------------------'
      endif
!!
      if((mode.eq.2).and.lelec.and.(nn_type_elec.eq.1))then
        if(optmodeq.eq.1)then
          write(*,*)'Using Kalman filter optimization (1) for atomic charges'
        elseif(optmodeq.eq.2)then
          write(*,*)'Using conjugate gradient optimization (2) for atomic charges'
        elseif(optmodeq.eq.3)then
          write(*,*)'Using steepest descent optimization (3) for atomic charges'
        elseif(optmodeq.eq.4)then
          write(*,*)'Using Levenberg Marquardt optimization (4) for atomic charges'
        else
          write(*,*)'Error: Unknown optimization mode ',optmodeq
          stop
        endif
      endif ! mode.eq.'2
!!
      if((mode.eq.2).and.lelec)write(*,'(a,f14.8)')' charge error threshold                      ',kalmanthresholde
!!
      if(lelec.and.(nn_type_elec.eq.1).and.(mode.eq.2))kalmanlambdae(:)=kalmanlambdae_local
      if((mode.eq.2).and.lelec.and.(nn_type_elec.eq.1).and.(optmodeq.eq.1))&
        write(*,'(a,f14.8)')' Kalman lambda (charge)                      ',kalmanlambdae_local
!!
      if((mode.eq.2).and.lelec.and.(nn_type_elec.eq.1).and.(optmodeq.eq.1))&
        write(*,'(a,f14.8)')' Kalman nue (charge)                         ',kalmannuee
!!
      if((mode.eq.2).and.lelec.and.(nn_type_elec.eq.1).and.(optmodeq.eq.1))&
        write(*,'(a,f14.8)')' Kalman damp (charge)                        ',kalman_dampq
!!
      if((mode.eq.2).and.lelec.and.(nn_type_elec.eq.1).and.(optmodeq.eq.3))then
        write(*,'(a,f14.8)')' steepest descent step size charges          ',steepeststepq
      endif
!!
      if(mode.eq.2)write(*,'(a,l)')' restart fit with old weights (charge)                   ',luseoldweightscharge
!!
      if((mode.eq.2).and.lelec.and.luseworstq)&
        write(*,'(a,f8.4)')' fraction of worst charges                         ',worstq
!!
      if(lelec.and.(nn_type_elec.eq.1).and.(mode.eq.2).and.(.not.lqgroupbystruct))then
        write(*,'(a,i8)')' grouping charges in blocks of                      ',nchargegroup
      endif
!!
      if(lelec.and.(nn_type_elec.eq.1).and.(mode.eq.2).and.(lqgroupbystruct))then
        write(*,'(a,i8)')' automatic grouping charges for update by structure'
      endif

      if(lelec.and.(mode.eq.2))write(*,'(a,f8.3)')&
        ' fraction of charges used for update               ',chargernd

      if((mode.eq.2).and.lelec.and.(nn_type_elec.eq.1).and.(.not.luseoldweightscharge))&
        write(*,'(a,f14.3)')' weightse_min                                ',weightse_min

      if((mode.eq.2).and.lelec.and.(nn_type_elec.eq.1).and.(.not.luseoldweightscharge))&
        write(*,'(a,f14.3)')' weightse_max                                ',weightse_max

      if((mode.eq.2).and.lelec.and.lchargeconstraint)then
        write(*,*)'-------------------------------------------------------------'
        write(*,*)'Fitting parameters charge constraint part:'
        write(*,*)'-------------------------------------------------------------'
      endif

      if(lelec.and.(mode.eq.2))write(*,'(a,l)')&
        ' using total charge constraint                           ',lchargeconstraint

      if((mode.eq.2).and.lelec.and.(lchargeconstraint))&
        write(*,'(a,f14.8)')' total charge error threshold                ',kalmanthresholdc

      if((mode.eq.2).and.lelec.and.lchargeconstraint.and.(optmodeq.eq.1))then
        write(*,'(a,f14.8)')' Kalman lambda (charge constraint)           ',kalmanlambdac
      endif

      if((mode.eq.2).and.lelec.and.lchargeconstraint.and.(optmodeq.eq.1))then
        write(*,'(a,f14.8)')' Kalman nue (charge constraint)              ',kalmannuec
      endif

      if((mode.eq.2).and.lelec)then
        write(*,'(a,f14.8,x,a2)')' noise charge threshold                      ',noiseq,' e'
      endif

      if(mode.eq.2)then
        write(*,*)'-------------------------------------------------------------'
        write(*,*)'Fitting output options:'
        write(*,*)'-------------------------------------------------------------'
      endif

      if(mode.eq.2)write(*,'(a,i6)')' write weights in every epoch                        ',iwriteweight

      if(mode.eq.2)write(*,'(a,l)')' write temporary weights each epoch                      ',lwritetmpweights

      if(mode.eq.2)write(*,'(a,l)')' write trainpoints.out and testpoints.out                ',lwritetrainpoints

      if((mode.eq.2).and.lelec)write(*,'(a,l)')&
        ' write traincharges.out and testcharges.out              ',lwritetraincharges

      if(mode.eq.2)write(*,'(a,l)')&
        ' write trainforces.out and testforces.out                ',lwritetrainforces

      if(mode.eq.3)then
        write(*,*)'-------------------------------------------------------------'
        write(*,*)'Options for prediction mode:'
        write(*,*)'-------------------------------------------------------------'
      endif

      if(mode.eq.3)write(*,'(a,l)')' rescale symmetry functions                              ',lscalesym

      if(mode.eq.3)write(*,'(a,l)')' remove CMS from symmetry functions                      ',lcentersym

      if(mode.eq.3)then
        if(lreadunformatted)then
          write(*,'(a)')' Reading unformatted files '
        else
          write(*,'(a)')' Reading formatted files '
        endif
      endif

      if(mode.eq.3)&
        write(*,'(a,l)')' calculation of analytic forces                          ',ldoforces

      if(mode.eq.3)&
        write(*,'(a,l)')' calculation of analytic Hessian                         ',ldohessian

      if(mode.eq.3)&
        write(*,'(a,l)')' calculation of analytic stress                          ',ldostress

      if((mode.eq.3).and.ldohessian)&
        write(*,'(a,l)')' calculation of analytic hessian                         ',ldohessian

      if(mode.eq.3)&
        write(*,'(a,l)')' write symmetry functions                                ',lwritesymfunctions

      if((mode.eq.3).and.lsens    )write(*,'(a,l)')' calculation of NN sensitivity                           ',lsens

      if(mode.eq.3)write(*,'(a,l)')' prepare md                                              ',lpreparemd

      write(*,*)'============================================================='





    end subroutine printinputnn

    subroutine readscale(filename,filename_error,ndim,iswitch,maxnum_funcvalues_local,num_funcvalues_local,minvalue_local,maxvalue_local,avvalue_local,eshortmin,eshortmax,chargemin,chargemax)

        implicit none

        integer             :: scaling_unit
        integer             :: ndim
        integer             :: maxnum_funcvalues_local
        integer             :: num_funcvalues_local(ndim)
        integer             :: counter_1, counter_2, counter_3
        integer             :: iswitch

        real(dp)            :: avvalue_local(ndim,maxnum_funcvalues_local)
        real(dp)            :: maxvalue_local(ndim,maxnum_funcvalues_local)
        real(dp)            :: minvalue_local(ndim,maxnum_funcvalues_local)
        real(dp)            :: thres
        real(dp)            :: eshortmin
        real(dp)            :: eshortmax
        real(dp)            :: chargemin(nelem)
        real(dp)            :: chargemax(nelem)

        logical             :: lexist

        character(len=*)                    :: err = "Error in readscale: "

        character(len=max_string_length)    :: filename
        character(len=*)                    :: filename_error

        integer, parameter  :: scaling_unit = 62 ! for scling.data and scalinge.data



        thres=0.00001d0

        call open_for_read(scaling_unit, filename); ios = 0

        do counter_1 = 1,ndim
            do counter_2 = 1,num_funcvalues_local(counter_1)
                read(scaling_unit, '(A)', iostat=ios) buffer
                line = line + 1

                if (ios == 0) then
                    call split_string(buffer, words, nwords)

                    if (nwords == 5) then
                        read(words(1),'(i1000)', iostat=ios) counter_3
                        if (ios /= 0) stop err // err_inpnn // "Error in line ", line, ", first argument value must be integer"
                        read(words(2),'(i1000)', iostat=ios) counter_3
                        if (ios /= 0) stop err // err_inpnn // "Error in line ", line, ", second argument value must be integer"
                        read(words(3),*, iostat=ios) minvalue_local(counter_1,counter_2)
                        if (ios /= 0) stop err // err_inpnn // "Error in line ", line, ", third argument value must be a number"
                        read(words(4),*, iostat=ios) maxvalue_local(counter_1,counter_2)
                        if (ios /= 0) stop err // err_inpnn // "Error in line ", line, ", fourth argument value must be a number"
                        read(words(5),*, iostat=ios) avvalue_local(counter_1,counter_2)
                        if (ios /= 0) stop err // err_inpnn // "Error in line ", line, ", fifth argument value must be a number"
                    else
                        print *, err, filename_error, "Error in line: ", line, "; need exactly 5 arguments"
                        stop
                    end if
                else
                    print *, err // filename_error // 'iostat = ', ios
                    stop
                end if

            end do
        end do

        read(scaling_unit, '(A)', iostat=ios) buffer

        if (ios == 0) then
            call split_string(buffer, words, nwords)

            if (iswitch == 1) then
                if (nwords == 2) then
                    read(words(1),*, iostat=ios) eshortmin
                    if (ios /= 0) stop err // err_inpnn // "Error in last line: " // "first argument value must be a number"
                    read(words(2),*, iostat=ios) eshortmax
                    if (ios /= 0) stop err // err_inpnn // "Error in last line: " // "second argument value must be a number"
                else
                    print *, err, filename_error, "Error in last line: need exactly 2 arguments"
                    stop
                end if
            else if (iswitch == 3) ! 3 to be comparable to RuNNer
                do counter_2 = 1,nelem
                    line = line + 1
                    if (nwords == 2) then
                        read(words(1),*, iostat=ios) chargemin(counter_2)
                        if (ios /= 0) stop err // err_inpnn // "Error in line: " // line // ", first argument value must be a number"
                        read(words(2),*, iostat=ios) chargemax(counter_2)
                        if (ios /= 0) stop err // err_inpnn // "Error in line: " // line // ", second argument value must be a number"
                    else
                        print *, err, filename_error, "Error in line: ", line, " need exactly 2 arguments"
                        stop
                    end if
                end do
            end if
        else
            write(*,*) err // filename_error // 'iostat = ', ios
            stop
        end if

        close(scaling_unit)

        do counter_1 = 1,ndim
            do counter_3 = 1,num_funcvalues_local(counter_1)
                if (minvalue_local(counter_1,counter_3) .gt. maxvalue_local(counter_1,counter_3)) then
                    print *, err // filename_error // 'No pairs of this type have been present in training set'
                else
                    if (abs(minvalue_local(counter_1,counter_3) - maxvalue_local(counter_1,counter_3)) .lt. thres) then
                        if (iswitch == 1) then
                            print *, err // filename_error // '### WARNING ###: minvalue=maxvalue ',counter_1,counter_3,nucelem(counter_1)
                        else if (iswitch == 3) then
                            print *, err // filename_error // '### WARNING ###: minvalue_elec=maxvalue_elec ',counter_1,counter_3,nucelem(counter_1)
                        end if
                        if (lscalesym) then
                            if (iswitch == 1) then
                                print *, err // filename_error // 'scaling symmetry functions cannot be used with minvalue=maxvalue'
                                stop
                            else if (iswitch == 3) then
                                print *, err // filename_error // 'scaling symmetry functions cannot be used with minvalue_elec=maxvalue_elec'
                                stop
                            end if
                        end if
                    end if
                end if
            end do
        end do

    end subroutine readscale

    subroutine readweights(directory,iswitch,ndim,maxnum_weights_local,num_weights_local,weights_local)

        implicit none

        integer             :: ndim
        integer             :: iswitch
        integer             :: icount
        integer             :: maxnum_weights_local
        integer             :: num_weights_local(ndim)
        integer             :: counter_1, counter_2, counter_3, counter_4

        real(dp)            :: weights_local(maxnum_weights_local,ndim)

        logical             :: lexist

        character(len=*)                    :: err = "Error in readweights: "
        character(len=*)                    :: err_weight = "Error when reading the following weight file: "
        character(len=*)                    :: err_weighte = "Error when reading the following weighte file: "

        character(len=max_string_length)    :: directory, filename_weight, filename_weighte
        character*40                        :: filename

        integer, parameter  :: weight_unit  = 64
        integer, parameter  :: weighte_unit = 65

        if (iswitch == 0) then
            do counter_1 = 1,ndim
                filename = 'weights.000.data'
                if (nucelem(counter_1) .gt. 99) then
                    write(filename(9:11),'(i3)') nucelem(counter_1)
                else if (nucelem(counter_1) .gt. 9) then
                    write(filename(10:11),'(i2)') nucelem(counter_1)
                else
                    write(filename(11:11),'(i1)') nucelem(counter_1)
                end if
                filename_weight = trim(directory) // trim(filename)
                if (.not. file_exists(filename_weight)) stop err // err_weight // trim(filename) // 'file does not exist'

                call open_for_read(weight_unit, filename_weight); ios = 0

                do counter_2 = 1,num_weights_local(counter_1)
                    read(weight_unit, '(A)', iostat=ios) buffer
                    line = line + 1

                    if (ios == 0) then
                        call split_string(buffer, words, nwords)

                        if (nwords == 1) then
                            read(words(1),*, iostat=ios) weights_local(counter_2,counter_1)
                            if (ios /= 0) stop err // err_weight // trim(filename) // "Error in line ", line, ", first argument value must be a number"
                        else
                            print *, err, err_weight, trim(filename), "Error in line ", line, "need exactly 1 argument"
                            stop
                        end if
                    else
                         write(*,*) err // err_weight // trim(filename) // 'iostat = ', ios
                         stop
                    end if
                end do

                close(weight_unit)

            end do
        else if (iswitch == 1) then
            do counter_1 = 1,ndim
                filename = 'weightse.000.data'
                if (nucelem(counter_1) .gt. 99) then
                    write(filename(9:11),'(i3)') nucelem(counter_1)
                else if (nucelem(counter_1) .gt. 9) then
                    write(filename(10:11),'(i2)') nucelem(counter_1)
                else
                    write(filename(11:11),'(i1)') nucelem(counter_1)
                end if
                filename_weighte = trim(directory) // trim(filename)
                if (.not. file_exists(filename_weighte)) stop err // err_weighte // trim(filename) // 'file does not exist'

                call open_for_read(weighte_unit, filename_weighte); ios = 0

                do counter_2 = 1,num_weights_local(counter_1)
                    read(weighte_unit, '(A)', iostat=ios) buffer
                    line = line + 1

                    if (ios == 0) then
                        call split_string(buffer, words, nwords)

                        if (nwords == 1) then
                            read(words(1),*, iostat=ios) weights_local(counter_2,counter_1)
                            if (ios /= 0) stop err // err_weighte // trim(filename) // "Error in line ", line, ", first argument value must be a number"
                        else
                            print *, err, err_weighte, trim(filename), "Error in line ", line, "need exactly 1 argument"
                            stop
                        end if
                    else
                         write(*,*) err // err_weighte // trim(filename) // 'iostat = ', ios
                         stop
                    end if
                end do

                close(weighte_unit)

            end do
        else
            write(*,*) err // "Error: unknown iswitch value ", iswitch
            stop
        end if

    end subroutine readweights

end module pes_nene_mod_supply
