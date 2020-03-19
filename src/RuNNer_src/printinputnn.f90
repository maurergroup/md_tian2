!######################################################################
! This routine is part of
! RuNNer - RuNNer Neural Network Energy Representation
! (c) 2008-2019 Prof. Dr. Joerg Behler 
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
!######################################################################
!! called by:
!! - readinput.f90
!!
      subroutine printinputnn(iseed,ielem,&
       nodes_short_atomic_temp,nodes_elec_temp,nodes_short_pair_temp,&
       kalmanlambda_local,kalmanlambdae_local,&
       actfunc_short_atomic_dummy,actfunc_elec_dummy,actfunc_short_pair_dummy)

!!
      use fileunits
      use fittingoptions
      use nnflags
      use globaloptions
      use mode1options
      use predictionoptions
      use nnshort_atomic
      use nnewald
      use nnshort_pair
      use inputnncounters
!!    
      implicit none
!!
      integer iseed                                     ! in 
      integer ielem                                     ! in 
      integer i                                         ! internal
      integer i1,i2,i3                                  ! internal
      integer nodes_short_atomic_temp(0:maxnum_layers_short_atomic)    ! in
      integer nodes_elec_temp(0:maxnum_layers_elec)    ! in
      integer nodes_short_pair_temp(0:maxnum_layers_short_pair)      ! in

!!
      real*8 kalmanlambda_local                         ! internal
      real*8 kalmanlambdae_local                        ! internal
!!
      character*1 actfunc_short_atomic_dummy(maxnum_layers_short_atomic) ! in
      character*1 actfunc_elec_dummy(maxnum_layers_elec) ! in
      character*1 actfunc_short_pair_dummy(maxnum_layers_short_pair)   ! in

!!
!!
      write(ounit,*)'General input parameters:'
      write(ounit,*)'-------------------------------------------------------------'
!!
      if(lshort)then
        write(ounit,*)'Short range NN is on'        
      else
        write(ounit,*)'Short range NN is off'        
      endif
!!
      if(lelec.and.(nn_type_elec.eq.1))then
        write(ounit,*)'Electrostatic NN is on'        
      else
        write(ounit,*)'Electrostatic NN is off'        
      endif
!!
      if(lvdw.and.(nn_type_vdw.eq.1))then
        write(ounit,*)'vdW corrections switched on'        
      else
        write(ounit,*)'vdW corrections switched off'        
      endif
!!
      if((mode.eq.1).and.lcheckinputforces)then
      write(ounit,'(a,f10.6,a)')&
        ' checking input forces, threshold for force vector is  '&
        ,inputforcethreshold,' Ha/Bohr'
      endif
!!
      write(ounit,*)'-------------------------------------------------------------'
!!
      if(lshort)then
        if(nn_type_short.le.2)then
          write(ounit,'(a55,i2)')&
            ' RuNNer nn_type_short                                              ',nn_type_short 
        else
          write(ounit,*)'ERROR: unknown nn_type_short: ',nn_type_short
          stop
        endif
      endif
!!
      if(mode.eq.1)then
        write(ounit,*)'RuNNer is started in mode for symmetry function calculation (1)' 
      elseif(mode.eq.2)then !'
        write(ounit,*)'RuNNer is started in mode for fitting (2)'
      elseif(mode.eq.3)then
        write(ounit,*)'RuNNer is started in mode for prediction (3)'
      else
        write(ounit,*)'Error: Unknown runner_mode: ',mode
        stop
      endif
!!
      write(ounit,'(a,l)')' debugging mode is                                       ',ldebug
!!
      write(ounit,'(a,i4)')' parallelization mode                                  ',paramode
!!
      write(ounit,'(a,l)')' enable detailed time measurement                        ',lfinetime
!!
      if(mode.eq.2)then
        write(ounit,'(a,l)')' enable detailed time measurement at epoch level         ',lfinetimeepoch
      endif
!!
      write(ounit,'(a,l)')' silent mode                                             ',lsilent
!!
      if((mode.eq.2).or.(mode.eq.3))then
        write(ounit,'(a,l)')' NN force check                                          ',lcheckf
      endif
!!
      if(nelem.lt.ielem)then
        write(ounit,*)'Error: number of elements in structure(s) is larger than '
        write(ounit,*)'number of elements in input.nn ',ielem,nelem
        stop
      else
        write(ounit,'(a,i4)')' number of elements                                    ',nelem
      endif
!!
      write(ounit,*)'elements (sorted):'
      do i1=1,nelem
        write(ounit,'(i3,x,a2)')nucelem(i1),element(i1)
      enddo
!!
!!
      write(ounit,'(a,i10)')' seed for random number generator                ',iseed
!!
      if((nran.lt.0).or.(nran.gt.6))then
        write(ounit,*)'ERROR: Unknown random number generator ',nran
        stop
      endif
      write(ounit,'(a,i10)')' random number generator type                    ',nran
!!
      write(ounit,'(a,l)')' remove free atom reference energies                     ',lremoveatomenergies
!!
      if(lfitethres.and.(mode.eq.1))then
        write(ounit,'(a,f7.3)')' upper energy threshold per atom (Ha)               ',fitethres
      endif
!!
      if(lfitfthres.and.(mode.eq.1))then
        write(ounit,'(a,f7.3)')' max force component threshold (Ha/Bohr)            ',fitfthres
      endif
!!
      write(ounit,'(a,f8.3)')' shortest allowed bond in structure                ',rmin
!!
      if(lnormnodes)then
        write(ounit,*)'Linear combinations at nodes are normalized'
      endif   
!!
      write(ounit,'(a,i3)')' Cutoff_type for symmetry function is                   ',cutoff_type
      write(ounit,'(a,f7.3)')' Cutoff_alpha for inner cutoff radius is            ',cutoff_alpha
!!
      if(lenforcemaxnumneighborsatomic)then
        write(ounit,'(a,i3)')&
        ' Enforcing global max_num_neighors_atomic               ',max_num_neighbors_atomic_input
      endif
!!
      if(lshort.and.(mode.ne.1))then
        write(ounit,*)'-------------------------------------------------------------'
        write(ounit,*)'Short range NN specifications:'
        write(ounit,*)'-------------------------------------------------------------'
      endif
!!
      if(lshort.and.(nn_type_short.eq.1).and.(mode.ne.1))then
        write(ounit,'(a,10i5)')' global hidden layers short range NN                  ',maxnum_layers_short_atomic-1
        write(ounit,'(a,10i5)')' global nodes hidden layers short NN             ',&
          (nodes_short_atomic_temp(i1),i1=1,maxnum_layers_short_atomic-1)
      endif
!!
      if(lshort.and.(nn_type_short.eq.1).and.(mode.ne.1))then
        write(ounit,'(a,x,10a)')' global activation functions short                     ',&
          (actfunc_short_atomic_dummy(i),i=1,maxnum_layers_short_atomic)
      endif 
!!
      if(lshort.and.(nn_type_short.eq.2).and.(mode.ne.1))then
        write(ounit,'(a,10i5)')' global hidden layers short range NN pair             ',maxnum_layers_short_pair-1
        write(ounit,'(a,10i5)')' global nodes hidden layers short NN pair        ',&
          (nodes_short_pair_temp(i1),i1=1,maxnum_layers_short_pair-1)
      endif
!!
!!
      if(lshort.and.(nn_type_short.eq.2).and.(mode.ne.1))then
        write(ounit,'(a,x,10a)')' global activation functions short pair                ',&
          (actfunc_short_pair_dummy(i),i=1,maxnum_layers_short_pair)
      endif 
!!
!!
      if(lelec)then
        write(ounit,*)'-------------------------------------------------------------'
        write(ounit,*)'Electrostatic specifications:'
        write(ounit,*)'-------------------------------------------------------------'
      endif
!!
      if(lelec)then
        write(ounit,'(a,i5)')' electrostatic_type (nn_type_elec)                    ',nn_type_elec
        if(nn_type_elec.eq.1)then
          write(ounit,'(a)')' Using separate set of atomic NNs for atomic charges'
        elseif(nn_type_elec.eq.2)then
          write(ounit,'(a)')' Constructing atomic charges from short range NN'
        elseif(nn_type_elec.eq.3)then
          write(ounit,'(a)')' Fixed atomic charges are used:' 
          do i1=1,nelem
            write(ounit,'(a1,a2,x,f14.3)')' ',element(i1),fixedcharge(i1)
          enddo
        elseif(nn_type_elec.eq.4)then
          write(ounit,'(a)')' Using atomic charges from charges.in file' 
        else
          write(ounit,*)'ERROR: Unknown electrostatic_type ',nn_type_elec
          stop
        endif
      endif
!!
      if(lelec.and.(nn_type_elec.eq.1))then
        write(ounit,'(a,10i5)')' global hidden layers electrostatic NN                ',maxnum_layers_elec-1
        write(ounit,'(a,10i5)')' global nodes hidden layers electrostatic NN     ',&
          (nodes_elec_temp(i1),i1=1,maxnum_layers_elec-1)
      endif
!!
      if(lelec.and.(nn_type_elec.eq.1))then
        write(ounit,'(a,x,10a)')' global activation functions electrostatic             ',&
          (actfunc_elec_dummy(i),i=1,maxnum_layers_elec)
      endif 
!!
      if(lelec)then
        write(ounit,'(a,f8.3)')' Ewald alpha                                       ',ewaldalpha
        write(ounit,'(a,f8.3)')' Ewald cutoff                                      ',ewaldcutoff
        write(ounit,'(a,i6)')' Ewald kmax                                          ',ewaldkmax
      endif
!!
      if(lelec.and.(mode.eq.0))then
        write(ounit,'(a,i4)')' Enforce total charge                                ',enforcetotcharge
      endif
!!
      if(lelec)then
        if(lscreen)then
          write(ounit,'(a,2f14.6)')' Screening electrostatics                            ',rscreen_onset,rscreen_cut
        else
          write(ounit,'(a)')' No screening of electrostatics requested                '
        endif
      endif
!!
      if(lelec)then
        if(mode.eq.3)then
          if(nn_type_elec.eq.4)then
            write(ounit,'(a)')' Using atomic charges from file charges.in!' 
          endif
        endif
      endif
!!
!!
      if(lvdw)then
        write(ounit,*)'-------------------------------------------------------------'
        write(ounit,*)'vdW specifications:'
        write(ounit,*)'-------------------------------------------------------------'
      endif
!!
      if(lvdw)then
        write(ounit,'(a,i4)')' vdW type                                            ',nn_type_vdw
      endif
      if(lvdw.and.(nn_type_vdw.eq.1))then
        write(ounit,'(a,2f14.6)')' vdw screening                                       ',vdw_screening(1),vdw_screening(2)
      endif







!!===============================================================
      if(lshort.and.(mode.eq.1))then
        write(ounit,*)'-------------------------------------------------------------'
        write(ounit,*)'Parameters for symmetry function generation: short range part:'
        write(ounit,*)'-------------------------------------------------------------'
      endif
!!
      if(lshort.and.(mode.eq.1)) write(ounit,'(a,l)')&
        ' using forces for fitting                                ',luseforces
!!
      if(lshort.and.(mode.eq.1)) write(ounit,'(a,l)')&
        ' using atomic energies for fitting                       ',luseatomenergies
!!
      if(lelec.and.(mode.eq.1))write(ounit,'(a,l)')&
        ' using atomic charges for fitting                        ',luseatomcharges
!!
      if(mode.eq.1)then
        write(ounit,*)'-------------------------------------------------------------'
        write(ounit,'(a,f8.4)') ' percentage of data for testing (%)                ',&
          100.d0*splitthres
      endif
!!
      if(mode.eq.2)then
        write(ounit,*)'-------------------------------------------------------------'
        write(ounit,*)'General fitting parameters:'
        write(ounit,*)'-------------------------------------------------------------'
      endif
!!
      if(mode.eq.2)then
        write(ounit,'(a,i8)')' number of fitting epochs                          ',nepochs
      endif ! mode.eq.2
!!
      if(mode.eq.2)then
        write(ounit,'(a,l)')' print date and time for each epoch                      ',lprintdateandtime
      endif ! mode.eq.2
!!
      if((mode.eq.2).and.lenableontheflyinput)then
        write(ounit,'(a,i8)')' on-the-fly input enabled          '
      endif ! mode.eq.2
!!
      if(mode.eq.2)then
        write(ounit,'(a,i8)')' number of data sets in memory                     ',nblock
      endif ! mode.eq.2
!!
      if(mode.eq.2)then
        if(fitmode.eq.1)then
          write(ounit,'(a,i8)')' Fitting mode 1 (online learning) selected         '
        elseif(fitmode.eq.2)then
          write(ounit,'(a,i8)')' Fitting mode 2 (offline learning) selected        '
        endif
      endif ! mode.eq.2
!!
      if(mode.eq.2)then
        write(ounit,'(a,l)')' Randomly mixing all points in training set              ',lmixpoints 
      endif
!!
      if(mode.eq.2)write(ounit,'(a,l)')' save Kalman filter data                                 ',lsavekalman
!!
      if(mode.eq.2)write(ounit,'(a,l)')' restart from old Kalman filter data                     ',lrestkalman
!!
      if(mode.eq.2)write(ounit,'(a,l)')' rescale symmetry functions                              ',lscalesym
!!
      if((mode.eq.2).and.lscalesym.and.lshort.and.(nn_type_short.eq.1))then
        write(ounit,'(a,f10.3)')' min value of scaled short range symmetry functions ',scmin_short_atomic
        write(ounit,'(a,f10.3)')' max value of scaled short range symmetry functions ',scmax_short_atomic
      endif
!!
      if((mode.eq.2).and.lscalesym.and.lshort.and.(nn_type_short.eq.2))then
        write(ounit,'(a,f10.3)')' min value of scaled pair symmetry functions       ',scmin_short_pair
        write(ounit,'(a,f10.3)')' max value of scaled pair symmetry functions       ',scmax_short_pair
      endif
!!
      if((mode.eq.2).and.lscalesym.and.lelec.and.(nn_type_elec.eq.1))then
        write(ounit,'(a,f10.3)')' min value of scaled electrostatic symmetry functions ',scmin_elec
        write(ounit,'(a,f10.3)')' max value of scaled electrostatic symmetry functions ',scmax_elec
      endif
!!
      if(mode.eq.2)write(ounit,'(a,l)')&
        ' remove CMS from symmetry functions                      ',lcentersym
!!
      if(mode.eq.2)write(ounit,'(a,l)')&
        ' calculate symmetry function correlation                 ',lpearson_correlation
!!
      if(mode.eq.2)write(ounit,'(a,l)')&
        ' weight analysis                                         ',lweightanalysis
!!
      if(mode.eq.2)write(ounit,'(a,l)')&
        ' environment analysis                                    ',lenvironmentanalysis
!!
      if(mode.eq.2)write(ounit,'(a,l)')&
        ' find contradictions                                     ',lfindcontradictions
      if((mode.eq.2).and.lfindcontradictions)then
        write(ounit,'(a,f10.3)')' threshold for |deltaG|                            ',deltagthres
        write(ounit,'(a,f10.3)')' threshold for delta|F|                            ',deltafthres
      endif
!!
      if(mode.eq.2)write(ounit,'(a,l)')' fix some weights                                        ',lfixweights
!!
      if(mode.eq.2)write(ounit,'(a,l)')' using growth mode for fitting                           ',lgrowth
      if((mode.eq.2).and.lgrowth)then
        write(ounit,'(a,i8)')' number of training structures in each growth step ',ngrowth
      endif
      if(lgrowth.and.(mode.eq.2))then
        write(ounit,'(a,i4)')' epochs with constant training set size in growth mode ',growthstep
      endif
!!
      if((mode.eq.2).and.ldampw)then
        write(ounit,'(a,l)')' using weight decay                                      ',ldampw
        write(ounit,'(a,f18.12)')' balance between error and weight decay  ',dampw
      endif
!!
      if((mode.eq.2).and.lupdatebyelement)then
        write(ounit,'(a,i3)')' do weight update just for one element                  ',elemupdate
        write(ounit,*)'### WARNING ### RMSEs will refer only to this element'
      endif
!!
      if(mode.eq.2)write(ounit,'(a,l)')&
        ' global fit of short and charge NN (not implemented)     ',lglobalfit
!!
      if(mode.eq.2)then
        if(fitting_unit.eq.1)then
          write(ounit,'(a,a2)')' error unit for fitting                                  ','eV'
        elseif(fitting_unit.eq.2)then
          write(ounit,'(a,a2)')' error unit for fitting                                  ','Ha'
        else
          write(ounit,*)'Error: add new energy unit in output of readinput.f90!!!'
          stop
        endif
      endif
!!
      if(mode.eq.2)then
        if(lreadunformatted)then
          write(ounit,'(a)')' Reading unformatted files ' 
        else
          write(ounit,'(a)')' Reading formatted files ' 
        endif
      endif
!!
      if(mode.eq.2)then
        if(lwriteunformatted)then
          write(ounit,'(a)')' Writing unformatted files ' 
        else
          write(ounit,'(a)')' Writing formatted files ' 
        endif
      endif
!!
      if(mode.eq.2)then
        if((optmodee.eq.1).or.(optmodef.eq.1).or.(optmodeq.eq.1))then
          write(ounit,'(a,l)')' Resetting Kalman filter matrices each epoch             ',lresetkalman
        endif
      endif
!!
      if(mode.eq.2)then
        if(nn_type_short.eq.1)then
          if(lshuffle_weights_short_atomic)then
            write(ounit,'(a,i5,f14.6)')' shuffle_weights_short_atomic                             ',&
              nshuffle_weights_short_atomic,shuffle_weights_short_atomic
          endif
        endif
      endif
!!
      if((mode.eq.2).and.lompmkl)then
        write(ounit,'(a)')' Using omp mkl for Kalman filter in parallel case' 
      endif
!!
      if((mode.eq.2).and.lionforcesonly)then
        write(ounit,'(a)')' Using only forces for fitting in case of ionic structures' 
      endif
!!
      if((mode.eq.2).and.lfitstats)then
        write(ounit,'(a)')' Writing fitting statistics ' 
      endif
!!
      if((mode.eq.2).and.(restrictw.gt.0.0d0))then
        write(ounit,'(a,f14.6)')' Restricting absolute value of weights       ',restrictw 
        if((restrictw.gt.0.0d0).and.(restrictw.lt.2.0d0))then
          write(ounit,*)'Currently restrictw must be larger than 2.0'
          stop
        endif
      endif
!!
      if((mode.eq.2).and.lbindingenergyonly.and.(.not.lwritetrainpoints))then
        write(ounit,*)'WARNING: write_binding_energy_only specified without write_trainpoints'
      endif
!!
      if((mode.eq.2).and.lanalyzeerror)then
        write(ounit,'(a)')' Error analysis requested for final epoch ' 
        if(lshort.and.(.not.lwritetrainpoints))then
          write(ounit,*)'WARNING: trainpoints file is required for short range energy error analysis'
          write(ounit,*)'=> This analysis will not be done'
        endif
        if(lshort.and.luseforces.and.(.not.lwritetrainforces))then
          write(ounit,*)'WARNING: trainforces file is required for short range force error analysis'
          write(ounit,*)'=> This analysis will not be done'
        endif
        if(lelec.and.(.not.lwritetraincharges))then
          write(ounit,*)'WARNING: traincharges file is required for charge error analysis'
          write(ounit,*)'=> This analysis will not be done'
        endif
      endif
!!
      if(mode.eq.2.and.((luseoldweightsshort).or.(luseoldweightscharge)))then
        write(ounit,'(a,l)')' Using old scaling data for restart          ',luseoldscaling 
      endif
!!
      if((mode.eq.2).and.(lprecond))then
        write(ounit,*)'Preconditioning of weights is switched on'
      endif
!!
      if((mode.eq.2).and.(linionly))then
        write(ounit,*)'Termination of mode 2 after initialization requested'
      endif
!!
      if((mode.eq.2).and.(ldataclustering))then
        write(ounit,'(a,2f14.10)')'data clustering requested with distance thresholds ',&
          dataclusteringthreshold1,dataclusteringthreshold2
      endif
!!
      if((mode.eq.2).and.(lprintconv))then
        write(ounit,*)'printing of convergence vector requested'
      endif
!!
      if((mode.eq.2).and.(lanalyzecomposition))then
        write(ounit,*)'analysis of chemical composition requested'
      endif
!!
      if((mode.eq.2).and.lshort)then
        write(ounit,*)'-------------------------------------------------------------'
        write(ounit,*)'Fitting parameters short range part:'
        write(ounit,*)'-------------------------------------------------------------'
      endif
!!
      if(lshort.and.(mode.eq.2)) write(ounit,'(a,l)')' using forces for fitting                                ',luseforces
!!
      if((mode.eq.2).and.lshort)then
        if(optmodee.eq.1)then
          write(ounit,'(a)')' using Kalman filter optimization (1) for short range energy'
          if(luseedkalman)then
            write(ounit,'(a)')' using element decoupled Kalman filter' 
            if(ledforcesv2)then
              write(ounit,'(a)')' using second variant of ED force fitting'
            endif
          endif
        elseif(optmodee.eq.2)then
          write(ounit,'(a)')' using conjugate gradient optimization (2) for short range energy'
        elseif(optmodee.eq.3)then
          write(ounit,'(a)')' using steepest descent optimization (3) for short range energy'
        else
          write(ounit,*)'Error: Unknown optimization mode ',optmodee
          stop
        endif
      endif ! mode.eq.2
!!
      if((mode.eq.2).and.lshort.and.luseforces)then
        if(optmodef.eq.1)then
          write(ounit,'(a)')' using Kalman filter optimization (1) for short range forces'
        elseif(optmodef.eq.2)then
          write(ounit,'(a)')' using conjugate gradient optimization (2) for short range forces'
        elseif(optmodef.eq.3)then
          write(ounit,'(a)')' using steepest descent optimization (3) for short range forces'
        else
          write(ounit,*)'Error: Unknown optimization mode ',optmodef
          stop
        endif
      endif ! mode.eq.2
!!
      if((mode.eq.2).and.lshort.and.(.not.lfixederrore))&
        write(ounit,'(a,f14.8)')' short energy error threshold                ',kalmanthreshold
!!
      if((mode.eq.2).and.lshort.and.(.not.lfixederrorf))&
        write(ounit,'(a,f14.8)')' short force error threshold                 ',kalmanthresholdf
!!
      if((mode.eq.2).and.lshort.and.lfixederrore)write(ounit,'(a,f14.8)')&
        ' fixed short energy error threshold          ',fixederrore
!!
      if((mode.eq.2).and.lshort.and.lfixederrorf)&
        write(ounit,'(a,f14.8)')' fixed short force error threshold           ',fixederrorf
!!
      if(mode.eq.2)then
        if(lshort.and.(nn_type_short.eq.1))kalmanlambda(:)=kalmanlambda_local
        if(lshort.and.(nn_type_short.eq.2))kalmanlambdap(:)=kalmanlambda_local
      endif
!!
      if((mode.eq.2).and.lshort.and.(optmodee.eq.1))&
        write(ounit,'(a,f14.8)')' Kalman lambda (short)                       ',kalmanlambda_local
!!
      if((mode.eq.2).and.lshort.and.(optmodee.eq.1))&
        write(ounit,'(a,f14.8)')' Kalman nue (short)                          ',kalmannue
!!
      if((mode.eq.2).and.lshort.and.(optmodee.eq.1))&
        write(ounit,'(a,l)')' use_noisematrix                                         ',lusenoisematrix  !! modifed by kenko
!!
      if((mode.eq.2).and.lshort.and.(optmodee.eq.1).and.lusenoisematrix)&
        write(ounit,'(a,f14.8)')' kalman_q0                                   ',kalman_q0  !! modifed by kenko
!!
      if((mode.eq.2).and.lshort.and.(optmodee.eq.1).and.lusenoisematrix)&
        write(ounit,'(a,f14.8)')' kalman_qmin                                 ',kalman_qmin  !! modifed by kenko
!!
      if((mode.eq.2).and.lshort.and.(optmodee.eq.1).and.lusenoisematrix)&
        write(ounit,'(a,f14.8)')' kalman_qtau                                 ',kalman_qtau  !! modifed by kenko
!!
      if((mode.eq.2).and.lshort.and.(optmodee.eq.1).and.lusenoisematrix)&
        write(ounit,'(a,f14.8)')' kalman_epsilon                              ',kalman_epsilon  !! modifed by kenko
!!
      if((mode.eq.2).and.lshort.and.(optmodee.eq.1))&
        write(ounit,'(a,f14.8)')' Kalman damp (short energy)                  ',kalman_dampe
!!
      if((mode.eq.2).and.lshort.and.(optmodee.eq.1))&
        write(ounit,'(a,f14.8)')' Kalman damp (short force)                   ',kalman_dampf
!!
      if((mode.eq.2).and.lshort.and.(optmodee.eq.3))then
        write(ounit,'(a,f14.8)')' steepest descent step size short energy     ',steepeststepe
      endif
!!
      if((mode.eq.2).and.lshort.and.(optmodef.eq.3))then
        write(ounit,'(a,f14.8)')' steepest descent step size short forces     ',steepeststepf
      endif
!!
      if(mode.eq.2)write(ounit,'(a,l)')' restart fit with old weights (short)                    ',luseoldweightsshort
!!
      if((mode.eq.2).and.lshort.and.luseworste)&
        write(ounit,'(a,f8.4)')' fraction of worst short range energies            ',worste
!!
      if((mode.eq.2).and.lshort.and.luseforces.and.luseworstf)&
        write(ounit,'(a,f8.4)')' fraction of worst short range forces              ',worstf
!!
      if((mode.eq.2).and.luseforces.and.lshort)then
        if(scalefactorf.lt.0.0d0)then
          write(ounit,'(a)')' automatic scaling factor for force update selected'
        else
          write(ounit,'(a,f11.8)')' scaling factor for force update (scalefactorf) ',scalefactorf
        endif
      endif
!!
      if(lshort.and.(mode.eq.2))&
        write(ounit,'(a,i8)')' grouping energies in blocks of                    ',nenergygroup
!!
      if(lshort.and.(mode.eq.2)) &
        write(ounit,'(a,f8.3)')' fraction of energies used for update              ',energyrnd
!!
      if(lshort.and.(mode.eq.2).and.(.not.lfgroupbystruct).and.(luseforces))then
        write(ounit,'(a,i8)')' grouping forces in blocks of                      ',nforcegroup
      endif
!!
      if(lshort.and.(mode.eq.2).and.(lfgroupbystruct).and.(luseforces))then
        write(ounit,'(a,i8)')' automatic grouping forces for update by structure'                     
      endif
!!
      if(lshort.and.(mode.eq.2))write(ounit,'(a,f8.3)')' fraction of forces used for update                ',forcernd
!!
      if((mode.eq.2).and.lshort.and.(.not.luseoldweightsshort))&
        write(ounit,'(a,f14.3)')' weights_min                                 ',weights_min
!!
      if((mode.eq.2).and.lshort.and.(.not.luseoldweightsshort))&
        write(ounit,'(a,f14.3)')' weights_max                                 ',weights_max
!!
      if((mode.eq.2).and.lshort.and.lseparatebiasini.and.(.not.luseoldweightsshort))&
        write(ounit,'(a,f14.3)')' biasweights_min                             ',biasweights_min
!!
      if((mode.eq.2).and.lshort.and.lseparatebiasini.and.(.not.luseoldweightsshort))&
        write(ounit,'(a,f14.3)')' biasweights_max                             ',biasweights_max
!!
      if((mode.eq.2).and.lshort.and.lnwweights)write(ounit,'(a)')' Using Nguyen Widrow weights for short range NN' 
!!
      if((mode.eq.2).and.lshort.and.lsysweights)write(ounit,'(a)')' Using systematic weights for short range NN' 
!!
      if((mode.eq.2).and.lelec.and.(nn_type_elec.eq.1).and.lnwweightse)&
        write(ounit,'(a)')' Using Nguyen Widrow weights for electrostatic NN' 
!!
      if((mode.eq.2).and.lelec.and.(nn_type_elec.eq.1).and.lsysweightse)&
        write(ounit,'(a)')' Using systematic weights for electrostatic NN' 
!!
      if((mode.eq.2).and.lshort.and.luseforces.and.lsepkalman.and.(optmodee.eq.1).and.(optmodef.eq.1))then
        write(ounit,'(a)')' Using separate Kalman filter matrices for short range energies and forces' 
      endif
!!
      if((mode.eq.2).and.lshort.and.luseforces.and.lrepeate)then
        write(ounit,'(a)')' Using repeated energy updates after each force update' 
      endif
!!
      if((mode.eq.2).and.lshort.and.(.not.luseforces).and.lfinalforce)then
        write(ounit,'(a)')' Calculating force error in final epoch only' 
      endif
!!
      if((mode.eq.2).and.(lshort))then
        write(ounit,'(a,f14.3)')' max_energy                                  ',maxenergy
      endif
!!
      if((mode.eq.2).and.(lshort).and.(luseforces))then
        write(ounit,'(a,f14.3,a)')' max force component used for fitting        ',&
          maxforce,' Ha/Bohr'
      endif
!!
      if((mode.eq.2).and.lshort)then
        write(ounit,'(a,f14.8,x,a7)')' noise energy threshold                      ',noisee,'Ha/atom'
      endif
!!
      if((mode.eq.2).and.luseforces.and.lshort)then
        write(ounit,'(a,f14.8,x,a7)')' noise force threshold                       ',noisef,'Ha/Bohr'
      endif
!!
      if((mode.eq.2).and.ldynforcegroup)then
        write(ounit,'(a,2i8)')' dynamic force grouping                      ',&
          dynforcegroup_start,dynforcegroup_step
      endif
!!
      if((mode.eq.2).and.ldetect_saturation.and.(lshort).and.(nn_type_short.eq.1))then
        write(ounit,'(a,f14.6)')' detect saturation of nodes is on            ',&
          saturation_threshold 
      endif
!!
      if((mode.eq.2).and.lelec)then
        write(ounit,*)'-------------------------------------------------------------'
        write(ounit,*)'Fitting parameters electrostatic part:'
        write(ounit,*)'-------------------------------------------------------------'
      endif
!!
      if((mode.eq.2).and.lelec.and.(nn_type_elec.eq.1))then
        if(optmodeq.eq.1)then
          write(ounit,*)'Using Kalman filter optimization (1) for atomic charges'
        elseif(optmodeq.eq.2)then
          write(ounit,*)'Using conjugate gradient optimization (2) for atomic charges'
        elseif(optmodeq.eq.3)then
          write(ounit,*)'Using steepest descent optimization (3) for atomic charges'
        elseif(optmodeq.eq.4)then
          write(ounit,*)'Using Levenberg Marquardt optimization (4) for atomic charges'
        else
          write(ounit,*)'Error: Unknown optimization mode ',optmodeq
          stop
        endif
      endif ! mode.eq.'2
!!
      if((mode.eq.2).and.lelec)write(ounit,'(a,f14.8)')' charge error threshold                      ',kalmanthresholde
!!
      if(lelec.and.(nn_type_elec.eq.1).and.(mode.eq.2))kalmanlambdae(:)=kalmanlambdae_local
      if((mode.eq.2).and.lelec.and.(nn_type_elec.eq.1).and.(optmodeq.eq.1))&
        write(ounit,'(a,f14.8)')' Kalman lambda (charge)                      ',kalmanlambdae_local
!!
      if((mode.eq.2).and.lelec.and.(nn_type_elec.eq.1).and.(optmodeq.eq.1))&
        write(ounit,'(a,f14.8)')' Kalman nue (charge)                         ',kalmannuee
!!
      if((mode.eq.2).and.lelec.and.(nn_type_elec.eq.1).and.(optmodeq.eq.1))&
        write(ounit,'(a,f14.8)')' Kalman damp (charge)                        ',kalman_dampq
!!
      if((mode.eq.2).and.lelec.and.(nn_type_elec.eq.1).and.(optmodeq.eq.3))then
        write(ounit,'(a,f14.8)')' steepest descent step size charges          ',steepeststepq
      endif
!!
      if(mode.eq.2)write(ounit,'(a,l)')' restart fit with old weights (charge)                   ',luseoldweightscharge
!!
      if((mode.eq.2).and.lelec.and.luseworstq)&
        write(ounit,'(a,f8.4)')' fraction of worst charges                         ',worstq
!!
      if(lelec.and.(nn_type_elec.eq.1).and.(mode.eq.2).and.(.not.lqgroupbystruct))then
        write(ounit,'(a,i8)')' grouping charges in blocks of                      ',nchargegroup
      endif
!!
      if(lelec.and.(nn_type_elec.eq.1).and.(mode.eq.2).and.(lqgroupbystruct))then
        write(ounit,'(a,i8)')' automatic grouping charges for update by structure'                     
      endif
!!
      if(lelec.and.(mode.eq.2))write(ounit,'(a,f8.3)')&
        ' fraction of charges used for update               ',chargernd
!!
      if((mode.eq.2).and.lelec.and.(nn_type_elec.eq.1).and.(.not.luseoldweightscharge))&
        write(ounit,'(a,f14.3)')' weightse_min                                ',weightse_min
!!
      if((mode.eq.2).and.lelec.and.(nn_type_elec.eq.1).and.(.not.luseoldweightscharge))&
        write(ounit,'(a,f14.3)')' weightse_max                                ',weightse_max
!!
      if((mode.eq.2).and.lelec.and.lchargeconstraint)then
        write(ounit,*)'-------------------------------------------------------------'
        write(ounit,*)'Fitting parameters charge constraint part:'
        write(ounit,*)'-------------------------------------------------------------'
      endif
!!
      if(lelec.and.(mode.eq.2))write(ounit,'(a,l)')&
        ' using total charge constraint                           ',lchargeconstraint
!!
      if((mode.eq.2).and.lelec.and.(lchargeconstraint))&
        write(ounit,'(a,f14.8)')' total charge error threshold                ',kalmanthresholdc
!!
      if((mode.eq.2).and.lelec.and.lchargeconstraint.and.(optmodeq.eq.1))then
        write(ounit,'(a,f14.8)')' Kalman lambda (charge constraint)           ',kalmanlambdac
      endif
!!
      if((mode.eq.2).and.lelec.and.lchargeconstraint.and.(optmodeq.eq.1))then
        write(ounit,'(a,f14.8)')' Kalman nue (charge constraint)              ',kalmannuec
      endif
!!
      if((mode.eq.2).and.lelec)then
        write(ounit,'(a,f14.8,x,a2)')' noise charge threshold                      ',noiseq,' e'
      endif
!!
      if(mode.eq.2)then
        write(ounit,*)'-------------------------------------------------------------'
        write(ounit,*)'Fitting output options:'
        write(ounit,*)'-------------------------------------------------------------'
      endif
!!
      if(mode.eq.2)write(ounit,'(a,i6)')' write weights in every epoch                        ',iwriteweight
!!
      if(mode.eq.2)write(ounit,'(a,l)')' write temporary weights each epoch                      ',lwritetmpweights
!!
      if(mode.eq.2)write(ounit,'(a,l)')' write trainpoints.out and testpoints.out                ',lwritetrainpoints
      if(mode.eq.2)write(ounit,'(a,l)')' write binding energies only                             ',lbindingenergyonly
!!
      if((mode.eq.2).and.lelec)write(ounit,'(a,l)')&
        ' write traincharges.out and testcharges.out              ',lwritetraincharges
!!
      if(mode.eq.2)write(ounit,'(a,l)')&
        ' write trainforces.out and testforces.out                ',lwritetrainforces
!!
      if(mode.eq.3)then
        write(ounit,*)'-------------------------------------------------------------'
        write(ounit,*)'Options for prediction mode:'
        write(ounit,*)'-------------------------------------------------------------'
      endif
!!
      if(mode.eq.3)write(ounit,'(a,l)')' rescale symmetry functions                              ',lscalesym
!!
      if(mode.eq.3)write(ounit,'(a,l)')' remove CMS from symmetry functions                      ',lcentersym
!!
      if(mode.eq.3)then
        if(lreadunformatted)then
          write(ounit,'(a)')' Reading unformatted files ' 
        else
          write(ounit,'(a)')' Reading formatted files ' 
        endif
      endif
!!
      if(mode.eq.3)&
        write(ounit,'(a,l)')' calculation of analytic forces                          ',ldoforces
!!
      if(mode.eq.3)&
        write(ounit,'(a,l)')' calculation of analytic Hessian                         ',ldohessian
!!
      if(mode.eq.3)&
        write(ounit,'(a,l)')' calculation of analytic stress                          ',ldostress
!!
      if((mode.eq.3).and.ldohessian)&
        write(ounit,'(a,l)')' calculation of analytic hessian                         ',ldohessian
!!
      if(mode.eq.3)&
        write(ounit,'(a,l)')' write symmetry functions                                ',lwritesymfunctions
!!
      if((mode.eq.3).and.lsens    )write(ounit,'(a,l)')' calculation of NN sensitivity                           ',lsens
!!
      if(mode.eq.3)write(ounit,'(a,l)')' prepare md                                              ',lpreparemd
!!
!!
      write(ounit,*)'============================================================='
!!


!!
      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
