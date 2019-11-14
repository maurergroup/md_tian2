!######################################################################
! This routine is part of
! RuNNer - RuNNer Neural Network Energy Representation
! (c) 2008-2018 Prof. Dr. Joerg Behler 
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
      subroutine inputnndefaults()
!!
      use nnflags 
      use globaloptions 
      use mode1options
      use predictionoptions
      use fittingoptions
      use nnshort_atomic
      use nnewald
      use nnshort_pair
!!
      implicit none
!!
      if(lshort.and.(nn_type_short.eq.1))then
        nodes_short_atomic(:,:)=0
      endif
      if(lshort.and.(nn_type_short.eq.2))then
        nodes_short_pair(:,:)=0
      endif
      if(lelec.and.(nn_type_elec.eq.1))then
        nodes_elec(:,:)=0
      endif
      analyze_error_energy_step = 0.01d0
      analyze_error_force_step = 0.01d0
      analyze_error_charge_step = 0.001d0
      ldebug        =.false.
      paramode      =1
      ewaldalpha=0.0d0
      ewaldcutoff=0.0d0
      ewaldkmax=0
      nenergygroup=1
      nforcegroup=1
      nchargegroup=1
      luseforces=.false.
      energyrnd=1.0d0
      forcernd=1.0d0
      chargernd=1.0d0
      luseatomcharges=.false.
      luseatomenergies=.false.
      lremoveatomenergies=.false.
      lchargeconstraint=.false.
      lfitethres=.false.
      fitethres=0.0d0
      lfitfthres=.false.
      fitfthres=0.0d0
      rmin=0.5d0
      optmodee=1
      optmodef=1
      optmodeq=1
      nblock=200
      nepochs=0
      iwriteweight=1
      lwritetmpweights=.false.
      lwritesymfunctions=.false.
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
      lscalesym=.false.
      lcentersym=.false.
      luseoldweightsshort=.false.
      luseoldweightscharge=.false.
      lglobalfit=.false.
      lsavekalman=.false.
      lrestkalman=.false.
      lupdatebyelement=.false.
      luseworste=.false.
      luseworstf=.false.
      luseworstq=.false.
      lgrowth=.false.
      ngrowth=0
      growthstep=1
      ldampw=.false.
      dampw=0.0d0
      lfixweights=.false.
      ldoforces=.false.
      ldohessian=.false.
      ldostress=.false.
      lfinetime=.false.
      lfinetimeepoch=.false.
      lwritetrainpoints=.false.
      lwritetrainforces=.false.
      lwritetraincharges=.false.
      atomrefenergies(:)=0.0d0
      weights_min=-1.d0
      weights_max=1.d0
      biasweights_min=-1.d0
      biasweights_max=1.d0
      weightse_min=-1.d0
      weightse_max=1.d0
      fitting_unit=1
      ljointefupdate=.false.
      pstring='00000000000000000000'
!!      nn_type_short =0
!!      nn_type_elec=0             ! no electrostatics
      nran    =5
      enforcetotcharge=0
      fixedcharge(:)=99.0d0
      lsysweights=.false.
      lsysweightse=.false.
      lsens=.false.
      lreadunformatted=.false.
      lwriteunformatted=.false.
      lresetkalman=.false.
      lsepkalman=.false.
      lrepeate=.false.
      maxforce=10000.d0
      maxenergy=10000.d0
      lfinalforce=.false.
      lcheckf=.false.
      lfitstats=.false.
      lfixederrore=.false.
      lfixederrorf=.false.
      lompmkl=.false.
      lnormnodes=.false.
      restrictw=-100000.d0
      fitmode=1            ! default is online learning
      lanalyzeerror=.false.
      lnwweights=.false.
      lnwweightse=.false.
      scmin_short_atomic=0.0d0  
      scmax_short_atomic=1.0d0
      scmin_short_pair=0.0d0  
      scmax_short_pair=1.0d0

      scmin_elec=0.0d0     
      scmax_elec=1.0d0    
      luseoldscaling=.false.
      lprecond=.false.
      linionly=.false.
      noisee=0.0d0
      noisef=0.0d0
      noiseq=0.0d0
      lprintconv=.false.
      lprintmad=.false.
      lfgroupbystruct=.false.
      lqgroupbystruct=.false.
      cutoff_type=1
      cutoff_alpha = 0.0d0
      lmixpoints=.false.
      lscreen=.false.
      rscreen_cut=0.0d0
      rscreen_onset=0.0d0
      lsilent=.false.
      lpreparemd=.false.
      lseparatebiasini=.false.
      lpearson_correlation=.false.
      lweightanalysis=.false.
      lenvironmentanalysis=.false.
      lfindcontradictions=.false.
      lmd=.false.
      ldynforcegroup=.false.
      dynforcegroup_start=20
      dynforcegroup_step=2
      lshuffle_weights_short_atomic=.false.
      nshuffle_weights_short_atomic=10
      shuffle_weights_short_atomic=0.1d0
      ldetect_saturation=.false.
      saturation_threshold=0.99d0
      dataclusteringthreshold1=1.0d0
      dataclusteringthreshold2=0.0d0
      ldataclustering=.false.
      lprintdateandtime=.false.
      lenableontheflyinput=.false.
      lcheckinputforces=.false.
      lionforcesonly=.false.
      inputforcethreshold=0.001d0
      luseedkalman=.false.
      ledforcesv2=.false.
      lprintforcecomponents=.false.
      kalman_epsilon = 1.0d0
!! KK:  Noise matrix default setting 
      lusenoisematrix = .false.
      kalman_q0 = 0.0d0
      kalman_qtau = 0.0d0
      kalman_qmin = 0.0d0
!!
      return
      end 
