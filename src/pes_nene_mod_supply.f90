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

! maybe use this module for more stuff related to RuNNer which has to be present in the pes_nen_mod module
! maybe move the subroutines back into pes_nene_mod.f90?
module pes_nene_mod_supply

    use constants, only : default_int, default_real, default_string, default_bool, max_string_length, dp

    ! RuNNer related modules (predict.f90)
    use fileunits
    use globaloptions
    use mpi_mod
    use nnewald
    use nnflags
    use nnshort_atomic
    !use nnshort_pair ! this module is not needed, Pair NN not implemented!!
    use predictionoptions
    !use saturation
    use symfunctions
    use timings

    implicit none

    integer :: ielem
    integer :: iseed

    logical :: lelement(102)

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
        nodes_short_local                   = default_int
        nodes_ewald_local                   = default_int
        num_funcvalues_local                = default_int
        num_funcvaluese_local               = default_int
        num_funcvaluesp_local               = default_int
        elementtemp                         = default_string
        ztemp                               = default_int
        maxnum_funcvalues_short_atomic      = default_int
        maxnum_funcvalues_elec              = default_int
        function_type_local                 = default_int
        function_type_temp                  = default_int
        funccutoff_local                    = default_real
        maxcutoff_local                     = default_real
        elementtemp1                        = default_string
        elementtemp2                        = default_string
        elementtemp3                        = default_string





































        ldebug = default_bool
        !maxnum_layers_short_atomic = default_int
        luseatomenergies = default_bool
        luseatomcharges = default_bool


    end subroutine set_defaults

    subroutine inputnndefaults()

        use nnflags
        use globaloptions
        use mode1options
        use predictionoptions
        use fittingoptions
        use nnshort_atomic
        use nnewald

        implicit none

        if(lshort.and.(nn_type_short.eq.1))then
            nodes_short_atomic(:,:)=0
        endif

        if(lelec.and.(nn_type_elec.eq.1))then
            nodes_elec(:,:)=0
        endif

        analyze_error_energy_step = 0.01d0
        analyze_error_force_step = 0.01d0
        analyze_error_charge_step = 0.001d0
        paramode      =1
        ewaldalpha=0.0d0
        ewaldcutoff=0.0d0
        ewaldkmax=0
        nenergygroup=1
        nforcegroup=1
        nchargegroup=1
        energyrnd=1.0d0
        forcernd=1.0d0
        chargernd=1.0d0
        !fitethres=0.0d0
        !fitfthres=0.0d0
        rmin=0.5d0
        optmodee=1
        optmodef=1
        optmodeq=1
        nblock=200
        !nepochs=0
        iwriteweight=1
        kalmanthreshold=0.0d0
        kalmanthresholdf=0.0d0
        kalmanthresholde=0.0d0
        kalmanthresholdc=0.0d0
        kalman_dampe=1.0d0
        kalman_dampf=1.0d0
        kalman_dampq=1.0d0
        steepeststepe=0.01d0
        steepeststepf=0.01d0
        steepeststepq=0.01d0
        scalefactorf=1.d0
        !ngrowth=0
        !growthstep=1
        dampw=0.0d0
        atomrefenergies(:)=0.0d0
        weights_min=-1.d0
        weights_max=1.d0
        biasweights_min=-1.d0
        biasweights_max=1.d0
        weightse_min=-1.d0
        weightse_max=1.d0
        fitting_unit=1
        pstring='00000000000000000000'
        nran    =5
        fixedcharge(:)=99.0d0
        maxforce=10000.d0
        maxenergy=10000.d0
        restrictw=-100000.d0
        !fitmode=1
        scmin_short_atomic=0.0d0
        scmax_short_atomic=1.0d0
        !scmin_short_pair=0.0d0
        !scmax_short_pair=1.0d0
        scmin_elec=0.0d0
        scmax_elec=1.0d0
        noisee=0.0d0
        noisef=0.0d0
        noiseq=0.0d0
        cutoff_type=1
        cutoff_alpha = 0.0d0
        rscreen_cut=0.0d0
        rscreen_onset=0.0d0
        dynforcegroup_start=20
        dynforcegroup_step=2
        nshuffle_weights_short_atomic=10
        shuffle_weights_short_atomic=0.1d0
        saturation_threshold=0.99d0
        dataclusteringthreshold1=1.0d0
        dataclusteringthreshold2=0.0d0
        inputforcethreshold=0.001d0
        kalman_epsilon = 1.0d0
        kalman_q0 = 0.0d0
        kalman_qtau = 0.0d0
        kalman_qmin = 0.0d0

    end subroutine inputnndefaults

    subroutine checkinputnn(err_main,err_file)

        use nnflags
        use globaloptions
        use inputnncounters
        use fittingoptions
        use predictionoptions
        use mode1options
        use nnshort_atomic
        use nnewald
        use fileunits

        implicit none

        integer counter

        character(len=*), parameter, intent(in) :: err_main, err_file
        character(len=*), parameter,            :: err_check = "Error in checkinputnn: "


        if (nran .neq. 5) then
            print *, err_main, err_file, err_check, "random_number_type not implemented, only 5 available"
        end if

      if(lfindcontradictions)then
        if(deltagthres.gt.0.0d0)then
          write(ounit,*)'ERROR: find_contradictions requires positive deltagthres ',&
            deltagthres
          stop
        endif
        if(deltafthres.gt.0.0d0)then
          write(ounit,*)'ERROR: find_contradictions requires positive deltafthres ',&
            deltafthres
          stop
        endif
      endif

      if((cutoff_alpha.gt.1.00000001d0).or.(cutoff_alpha.lt.0.0d0))then
        write(ounit,*)'ERROR: please use cutoff_alpha within 0 and 1 ',cutoff_alpha
        stop
      endif

      if(lusenoisematrix) then
        if((kalman_q0 .le. 0.0d0 ).or.(kalman_qmin.le.0.0d0).or.(kalman_qtau.le.0.0d0)) then
          write(ounit,*)'ERROR: please define the q0,qmin,qtau for noise matrix ', &
          'and use them larger than zero ',kalman_q0,kalman_qmin,kalman_qtau
          stop
        endif
      endif

      if(lfixweights.and.lshuffle_weights_short_atomic)then
        write(ounit,*)'ERROR: shuffle_weights_short_atomic cannot be combined with fixed weights'
        stop
      endif

      if(lscreen)then
        if(rscreen_onset.gt.rscreen_cut)then
          write(ounit,*)'ERROR: rscreen_onset .gt. rscreen_cut in screen_electrostatics'
          stop
        endif
        if(rscreen_onset.lt.0.0d0)then
          write(ounit,*)'ERROR: rscreen_onset .lt. 0 in screen_electrostatics'
          stop
        endif
        if(rscreen_cut.lt.0.0d0)then
          write(ounit,*)'ERROR: rscreen_cut .lt. 0 in screen_electrostatics'
          stop
        endif
      endif

      if(noisee.lt.0.0d0)then
        write(ounit,*)'ERROR: noise_energy must not be negative ',noisee
        stop
      endif

      if(noisef.lt.0.0d0)then
        write(ounit,*)'ERROR: noise_force must not be negative ',noisef
        stop
      endif

      if(noiseq.lt.0.0d0)then
        write(ounit,*)'ERROR: noise_charge must not be negative ',noiseq
        stop
      endif

      if(lsysweights.and.lnwweights)then
        write(ounit,'(a)')'Error: Cannot use systematic_weights_short and nguyen_widrow_weights_short together!'
        stop
      endif

      if(lsysweightse.and.lnwweightse)then
        write(ounit,'(a)')'Error: Cannot use systematic_weights_ewald and nguyen_widrow_weights_ewald together!'
        stop
      endif

      if(lnormnodes.and.lnwweights)then
        write(ounit,'(a)')'Error: Cannot use normalize_nodes and nguyen_widrow_weights_short together!'
        stop
      endif

      if(lnormnodes.and.lnwweightse)then
        write(ounit,'(a)')'Error: Cannot use normalize_nodes and nguyen_widrow_weights_ewald together!'
        stop
      endif

      if((count_kalmanthreshold.eq.1).and.(count_lfixederrore.eq.1))then
        write(ounit,'(2a)')'Error: short_energy_error_threshold cannot be used ',&
          'in combination with fixed_short_energy_error_threshold'
        stop
      endif

      if((count_kalmanthresholdf.eq.1).and.(count_lfixederrorf.eq.1))then
        write(ounit,'(a)')'Error: short_force_error_thresholdf cannot be used in combination with fixed_short_force_error_threshold'
        stop
      endif

      if(count_mode.eq.0)then
        write(ounit,*)'Error: runner_mode is not specified'
        stop
      endif

      if((.not.lshort).and.(.not.lelec))then
        write(ounit,*)'Error: short range and electrostatic NNs are switched off'
        stop
      endif

      if(lshort.and.(maxnum_layers_short_atomic.eq.0).and.(nn_type_short.eq.1))then
        write(ounit,*)'Error: global_hidden_layers_short is not specified'
        stop
      endif

      if(lshort.and.(maxnum_layers_short_pair.eq.0).and.(nn_type_short.eq.2))then
        write(ounit,*)'Error: global_hidden_layers_pair is not specified'
        stop
      endif

      if(lelec.and.(nn_type_elec.eq.0))then
        write(ounit,*)'Error: electrostatic_type is not specified'
        stop
      endif

      if(lelec.and.(nn_type_elec.eq.1).and.(maxnum_layers_elec.eq.0))then
        write(ounit,*)'Error: global_hidden_layers_electrostatic is not specified'
        stop
      endif

      if(lshort.and.(count_nodes_short_atomic.eq.0).and.(nn_type_short.eq.1))then
        write(ounit,*)'Error: global_nodes_short is not specified'
        stop
      endif

      if(lelec.and.(nn_type_elec.eq.1).and.(count_nodes_elec.eq.0))then
        write(ounit,*)'Error: global_nodes_electrostatic is not specified'
        stop
      endif

      if(lshort.and.(count_nodes_short_pair.eq.0).and.(nn_type_short.eq.2))then
        write(ounit,*)'Error: global_nodes_pair is not specified'
        stop
      endif

      do i1=1,nelem
        if(lshort.and.(nn_type_short.eq.1))then
          if(nodes_short_atomic(maxnum_layers_short_atomic,i1).gt.1)then
            write(ounit,*)'Error: More than 1 output node currently does '
            write(ounit,*)'make sense in short range NN'
            stop
          endif
        endif
      enddo ! i1

      do i1=1,nelem
        if(lshort.and.(nn_type_short.eq.1))then
          if(nodes_short_atomic(maxnum_layers_short_atomic,i1).eq.0)then
            write(ounit,*)'Error: output_nodes_short is 0'
            stop
          endif
        endif
      enddo ! i1

      do i1=1,nelem
        if(lelec.and.(nn_type_elec.eq.1))then
          if(nodes_elec(maxnum_layers_elec,i1).gt.1)then
            write(ounit,*)'Error: More than 1 output node currently does '
            write(ounit,*)'make sense in electrostatic NN'
            stop
          endif
        endif
      enddo ! i1

      do i1=1,nelem
        if(lelec.and.(nn_type_elec.eq.1))then
          if(nodes_elec(maxnum_layers_elec,i1).eq.0)then
            write(ounit,*)'Error: output_nodes_electrostatic is 0'
            stop
          endif
        endif
      enddo ! i1

      do i1=1,nelem
        if(lshort.and.(nn_type_short.eq.1))then
          if(nodes_short_atomic(0,i1).eq.0)then
            write(ounit,*)'Error: input_nodes_short is 0'
            stop
          endif
        endif
      enddo ! i1

      do i1=1,nelem
        if(lshort.and.(nn_type_short.eq.1))then
          if(nodes_short_atomic(0,i1).ne.num_funcvalues_short_atomic(i1))then
            write(ounit,*)'Error: num_funcvalues_short_atomic .ne. nodes_short_atomic(0)',&
              num_funcvalues_short_atomic(i1),nodes_short_atomic(0,i1)
            write(ounit,*)'Did you set the right number of input nodes?'
            stop
          endif
        endif
      enddo ! i1

      do i1=1,nelem
        if(lelec.and.(nn_type_elec.eq.1))then
          if(nodes_elec(0,i1).eq.0)then
            write(ounit,*)'Error: input_nodes_electrostatic is 0'
            stop
          endif
        endif
      enddo ! i1

      do i1=1,nelem
        if(lelec.and.(nn_type_elec.eq.1))then
          if(nodes_elec(0,i1).ne.num_funcvalues_elec(i1))then
            write(ounit,*)'Error: num_funcvalues_elec .ne. nodes_elec(0)',&
              num_funcvalues_elec(i1),nodes_elec(0,i1)
            write(ounit,*)'Did you set the right number of input nodes?'
            stop
          endif
        endif
      enddo ! i1

      if(lshort.and.(nn_type_short.eq.1))then
        if(count_global_activation_short_atomic.eq.0)then
          write(ounit,*)'Error: global_activation_short is not specified'
          stop
        endif
      endif

      if(lelec.and.(nn_type_elec.eq.1))then
        if(count_global_activation_elec.eq.0)then
          write(ounit,*)'Error: global_activation_ewald is not specified'
          stop
        endif
      endif

      if(lelec.and.(count_ewaldalpha.eq.0))then
        write(ounit,*)'Error: ewald_alpha must be specified for electrostatic NN'
        stop
      endif

      if(lelec.and.(ewaldalpha.le.0))then
        write(ounit,*)'Error: ewald_alpha must be positive ',ewaldalpha
        stop
      endif
!!
      if(lelec.and.(count_ewaldcutoff.eq.0))then
        write(ounit,*)'Error: ewald_cutoff must be specified for electrostatic NN'
        stop
      endif
!!
      if(lelec.and.(ewaldcutoff.le.0))then
        write(ounit,*)'Error: ewald_cutoff must be positive ',ewaldcutoff
        stop
      endif
!!
      if(lelec.and.(count_ewaldkmax.eq.0))then
        write(ounit,*)'Error: ewald_kmax must be specified for electrostatic NN'
        stop
      endif
!!
      if((.not.lshort).and.(luseforces))then
        write(ounit,*)'### WARNING ### switching off use_short_forces because no short range NN is used'
        luseforces=.false.
      endif
!!
      if(lelec.and.(.not.luseatomcharges))then
        write(ounit,*)'### WARNING ### use_atom_charges is switched on for electrostatic NN'
        luseatomcharges=.true.
      endif
!!
      if(lshort.and.(luseatomenergies))then
        write(ounit,*)'### WARNING ### use_atom_energies is switched off (not implemented)'
        luseatomenergies=.false.
      endif
!!
      if((.not.lshort).and.(lremoveatomenergies))then
        write(ounit,*)'### WARNING ### remove_atom_energies is switched on without short range NN'
      endif
!!
      if(lelec.and.(lchargeconstraint))then
        write(ounit,'(a)')' ### WARNING ### use_charge_constraint is not maintained at the moment and might fail'
      endif
!!
      if(count_iseed.eq.0)then
        write(ounit,*)'### WARNING ### no random_seed specified, using default '
      endif
!!
      if(nenergygroup.gt.nblock)then
        nenergygroup=nblock
        write(ounit,*)'### WARNING ### reducing nenergygroup to nblock'
      endif

      if(count_nelem.eq.0)then
        write(ounit,*)'Error: number_of_elements not specified'
        stop
      endif
!!
      if(count_element.eq.0)then
        write(ounit,*)'Error: elements not specified'
        stop
      endif
!!
      if((mode.eq.1).and.(count_splitthres.eq.0))then
        write(ounit,*)'Error: test_fraction not specified'
        stop
      endif
!!
      if(lcentersym.and.(count_scmin_short_atomic.gt.0))then
        write(ounit,'(a)')'Error: center_symmetry_functions cannot be combined with scale_min_short_atomic keyword'
        stop
      endif
!!
      if(lcentersym.and.(count_scmax_short_atomic.gt.0))then
        write(ounit,'(a)')'Error: center_symmetry_functions cannot be combined with scale_max_short_atomic keyword'
        stop
      endif
!!
      if(lcentersym.and.(count_scmin_short_pair.gt.0))then
        write(ounit,'(a)')'Error: center_symmetry_functions cannot be combined with scale_min_short_pair keyword'
        stop
      endif
!!
      if(lcentersym.and.(count_scmax_short_pair.gt.0))then
        write(ounit,'(a)')'Error: center_symmetry_functions cannot be combined with scale_max_short_pair keyword'
        stop
      endif
!!
      if(lcentersym.and.(count_scmin_elec.gt.0))then
        write(ounit,'(a)')'Error: center_symmetry_functions cannot be combined with scale_min_elec keyword'
        stop
      endif
!!
      if(lcentersym.and.(count_scmax_elec.gt.0))then
        write(ounit,'(a)')'Error: center_symmetry_functions cannot be combined with scale_max_elec keyword'
        stop
      endif
!!
      if((count_scmin_short_atomic.gt.0).and.(.not.lscalesym))then
        write(ounit,*)'Error: scale_min_short requires keyword scale_symmetry_functions'
        stop
      endif
!!
      if((count_scmax_short_atomic.gt.0).and.(.not.lscalesym))then
        write(ounit,*)'Error: scale_max_short requires keyword scale_symmetry_functions'
        stop
      endif
!!
      if((count_scmin_short_pair.gt.0).and.(.not.lscalesym))then
        write(ounit,*)'Error: scale_min_short_pair requires keyword scale_symmetry_functions'
        stop
      endif
!!
      if((count_scmax_short_pair.gt.0).and.(.not.lscalesym))then
        write(ounit,*)'Error: scale_max_short_pair requires keyword scale_symmetry_functions'
        stop
      endif
!!
      if((count_scmin_elec.gt.0).and.(.not.lscalesym))then
        write(ounit,*)'Error: scale_min_elec requires keyword scale_symmetry_functions'
        stop
      endif
!!
      if((count_scmax_elec.gt.0).and.(.not.lscalesym))then
        write(ounit,*)'Error: scale_max_elec requires keyword scale_symmetry_functions'
        stop
      endif
!!
      if(scmin_short_atomic.ge.scmax_short_atomic)then
        write(ounit,'(a)')'Error: scale_min_short .ge. scale_max_short'
        stop
      endif
!!
      if(scmin_short_pair.ge.scmax_short_pair)then
        write(ounit,'(a)')'Error: scale_min_short_pair .ge. scale_max_short_pair'
        stop
      endif
!!
      if(scmin_elec.ge.scmax_elec)then
        write(ounit,'(a)')'Error: scale_min_elec .ge. scale_max_elec'
        stop
      endif
!!
      if(lupdatebyelement.and.lchargeconstraint)then
        lchargeconstraint=.false.
        if(mode.eq.2)then
          write(ounit,*)'### WARNING ### lchargeconstraint is switched off because of lupdatebyelement'
        endif
      endif
!!
      if(lshort.and.lupdatebyelement.and.(mode.eq.2))then
        write(ounit,*)'### WARNING ### lupdatebyelement works only for charges and forces'
      endif
!!
      if(luseworste.and.lshort.and.(energyrnd.lt.1.0d0))then
        energyrnd=1.0d0
        write(ounit,*)'### WARNING ### luseworste overrides energyrnd: ',energyrnd
      endif
!!
      if(luseworstf.and.lshort.and.luseforces.and.(forcernd.lt.1.0d0))then
        forcernd=1.0d0
        write(ounit,*)'### WARNING ### luseworstf overrides forcernd: ',forcernd
      endif
!!
      if(luseworstq.and.lelec.and.(nn_type_elec.eq.1).and.(chargernd.lt.1.0d0))then
        chargernd=1.0d0
        write(ounit,*)'### WARNING ### luseworstq overrides chargernd: ',chargernd
      endif
!!
      if(dampw.gt.1.0d0)then
        write(ounit,*)'Error: dampw must not be larger than 1.0d0 ',dampw
        stop
      endif
!!
      if(ldostress.and.(.not.ldoforces))then
        write(ounit,*)'### WARNING ### Analytic stress is requested without forces'
        write(ounit,*)'Switching on calculation of analytic forces'
        ldoforces=.true.
      endif
!!
      if(ldohessian.and.(.not.ldoforces))then
        write(ounit,*)'### WARNING ### Analytic Hessian is requested without forces'
        write(ounit,*)'Switching on calculation of analytic forces'
        ldoforces=.true.
      endif
!!
      if(ldostress.and.(mode.eq.1))then
        write(ounit,*)'### WARNING ### switching off stress calculation in mode 1 for increased performance'
        ldostress=.false.
      endif
!!
      if((count_wconstraint.gt.0).and.(.not.lfixweights))then
        write(ounit,*)'Error: weight constraints are specified without fix_weights keyword'
        stop
      endif

      if((count_wconstraint.eq.0).and.(lfixweights))then
        write(ounit,*)'Error: no weights constrained but keyword fix_weights has been selected'
        stop
      endif
!!
      if(weights_min.ge.weights_max)then
        write(ounit,*)'Error: weights_min > weights_max'
        stop
      endif
!!
      if(biasweights_min.ge.biasweights_max)then
        write(ounit,*)'Error: biasweights_min > biasweights_max'
        stop
      endif
!!
      if(weightse_min.ge.weightse_max)then
        write(ounit,*)'Error: weightse_min > weightse_max'
        stop
      endif
!!
      if(kalman_dampe.lt.0.0d0)then
        write(ounit,*)'ERROR: kalman_damp_short must be non-negative ',kalman_dampe
        stop
      endif
!!
      if(kalman_dampf.lt.0.0d0)then
        write(ounit,*)'ERROR: kalman_damp_force must be non-negative ',kalman_dampf
        stop
      endif
!!
      if(kalman_dampq.lt.0.0d0)then
        write(ounit,*)'ERROR: kalman_damp_charge must be non-negative ',kalman_dampq
        stop
      endif
!!
      if(ljointefupdate.and.lelec)then
        write(ounit,*)'ERROR: joint_energy_force_update is not implemented for lelec and nn_type_elec 2'
        stop
      endif
!!
      if(ljointefupdate)then
        if(optmodee.ne.optmodef)then
          write(ounit,*)'Error: joint_energy_force_update requires to use the'
          write(ounit,*)'same optimization algorithm for energy and forces'
          stop
        endif
        if(.not.luseforces)then
          write(ounit,*)'Error: switch on use_short_forces for joint_energy_force_update'
          stop
        endif
        if(lrepeate)then
          write(ounit,*)'ERROR: repeated energy update cannot be combined with joint energy and force update'
          stop
        endif
        if(forcernd.lt.1.0d0)then
          write(ounit,*)'ERROR: joint energy and force update requires force_fraction = 1.0d0'
          stop
        endif
        if(luseworste)then
          write(ounit,*)'ERROR: joint energy and force update cannot be combined with update_worst_short_energies'
          stop
        endif
        if(luseworstf)then
          write(ounit,*)'ERROR: joint energy and force update cannot be combined with update_worst_short_forces'
          stop
        endif
        if(nenergygroup.gt.1)then
          write(ounit,*)'ERROR: joint energy and force update cannot be combined with short_energy_group > 1'
          stop
        endif
        if(nforcegroup.gt.1)then
          write(ounit,*)'ERROR: joint energy and force update cannot be combined with short_force_group > 1'
          stop
        endif
        if(kalmanthresholdf.gt.0.0d0)then
          write(ounit,*)'ERROR: joint energy and force update cannot be combined with short_force_error_threshold > 0.0'
          stop
        endif
      endif
!!
      if(maxforce.le.0.0d0)then
        write(ounit,*)'Error: max_force must not be negative ',maxforce
        stop
      endif
!!
      if(lshort.and.(nn_type_short.eq.1))then
        do i1=1,nelem
          if(num_funcvalues_short_atomic(i1).ne.nodes_short_atomic(0,i1))then
            write(ounit,*)'Error: num_funcvalues_short_atomic .ne. nodes_short_atomic(0)'
            write(ounit,*)i1,num_funcvalues_short_atomic(i1),nodes_short_atomic(0,i1)
            stop
          endif
        enddo! i1
      endif
!!
      if(lelec.and.(nn_type_elec.eq.1))then
        do i1=1,nelem
          if(num_funcvalues_elec(i1).ne.nodes_elec(0,i1))then
            write(ounit,*)'Error: num_funcvalues_elec .ne. nodes_elec(0)'
            write(ounit,*)i1,num_funcvalues_elec(i1),nodes_elec(0,i1)
            stop
          endif
        enddo ! i1
      endif
!!
      if((nn_type_elec.eq.4).and.(mode.ne.3))then
        write(ounit,*)'ERROR: electrostatic_type 4 is only valid for prediction mode'
        stop
      endif
!!
      if((mode.eq.3).and.(max_num_atoms.lt.nblock).and.(nn_type_short.eq.1).and.lshort)then
        write(ounit,*) 'WARNING: reducing points_in_memory to max_num_atoms ',max_num_atoms
        nblock=max_num_atoms
      endif

    end subroutine checkinputnn

    subroutine readscale(filename,filename_unit,filename_error,ndim,iswitch,maxnum_funcvalues_local,num_funcvalues_local,minvalue_local,maxvalue_local,avvalue_local,eshortmin,eshortmax,chargemin,chargemax)

        use fileunits
        use globaloptions

        implicit none

        integer             :: filename_unit
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



        thres=0.00001d0

        call open_for_read(filename_unit, filename); ios = 0

        do counter_1 = 1,ndim
            do counter_2 = 1,num_funcvalues_local(counter_1)
                read(filename_unit, '(A)', iostat=ios) buffer
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

        read(filename_unit, '(A)', iostat=ios) buffer

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
            else if (iswitch == 3)
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

        close(filename_unit)

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

        use fileunits
        use globaloptions
        use nnflags

        implicit none

        integer             :: filename_unit
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

        integer, parameter  :: weight_unit      = 64
        integer, parameter  :: weighte_unit     = 65

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
