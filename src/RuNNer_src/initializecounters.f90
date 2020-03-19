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
      subroutine initializecounters()
!!
      use inputnncounters
!!
      implicit none
!!
      count_analyze_error_energy_step   =0
      count_analyze_error_force_step   =0
      count_analyze_error_charge_step   =0
      count_mode    =0
      count_vdwscreening   =0
      count_nn_type_elec   =0
      count_nn_type_vdw    =0
      count_ldebug  =0
      count_paramode=0
      count_lshort  =0
      count_lelec   =0
      count_lvdw    =0
      count_num_layers_short_atomic=0
      count_num_layers_elec=0
      count_num_layers_short_pair=0
      count_nodes_short_atomic=0
      count_nodes_elec=0
      count_nodes_short_pair=0
      count_global_activation_short_atomic=0
      count_global_activation_elec=0
      count_global_activation_short_pair=0
      count_ewaldalpha=0
      count_ewaldcutoff=0
      count_ewaldkmax=0
      count_short_energy_group=0
      count_short_force_group=0
      count_charge_group=0
      count_luseforces=0
      count_energyrnd=0
      count_forcernd=0
      count_chargernd=0
      count_luseatomcharges=0
      count_luseatomenergies=0
      count_lremoveatomenergies=0
      count_lchargeconstraint=0
      count_fitethres=0
      count_fitfthres=0
      count_rmin=0
      count_optmodee=0
      count_optmodef=0
      count_optmodeq=0
      count_iseed=0
      count_nblock=0
      count_epochs=0
      count_nelem=0
      count_element=0
      count_iwriteweight=0
      count_lwritetmpweights=0
      count_lwritesymfunctions=0
      count_splitthres=0
      count_kalmanthreshold=0
      count_kalmanthresholdf=0
      count_kalmanthresholde=0
      count_kalmanthresholdc=0
      count_kalmanlambda=0
      count_kalmanlambdae=0
      count_kalmanlambdac=0
      count_kalmannue=0
      count_kalmannuee=0
      count_kalmannuec=0
      count_kalman_dampe=0
      count_kalman_dampf=0
      count_kalman_dampq=0
      count_steepeststepe=0
      count_steepeststepf=0
      count_steepeststepq=0
      count_scalefactorf=0
      count_scalefactorq=0
      count_lscalesym=0
      count_lcentersym=0
      count_luseoldweightsshort=0
      count_luseoldweightscharge=0

      count_lglobalfit=0
      count_lsavekalman=0
      count_lrestkalman=0
      count_elemupdate=0
      count_luseworste=0
      count_luseworstf=0
      count_luseworstq=0
      count_growth=0
      count_dampw=0
      count_lfixweights=0
      count_ldoforces=0
      count_ldostress=0
      count_ldohessian=0
      count_lfinetime=0
      count_lfinetimeepoch=0
      count_lwritetrainpoints=0
      count_lbindingenergyonly=0
      count_lwritetrainforces=0
      count_lwritetraincharges=0
      count_wconstraint=0
      count_weights_min=0
      count_weights_max=0
      count_lseparatebiasini=0
      count_biasweights_min=0
      count_biasweights_max=0
      count_weightse_min=0
      count_weightse_max=0
      count_fittingunit=0
      count_ljointefupdate=0
      count_nn_type_short=0
      count_nran=0
      count_enforcetotcharge=0  
      count_lsysweights=0 
      count_lsysweightse=0 
      count_lsens=0 
      count_lreadunformatted=0 
      count_lwriteunformatted=0 
      count_lresetkalman=0 
      count_lsepkalman=0 
      count_lrepeate=0 
      count_lfinalforce=0 
      count_maxforce=0
      count_maxenergy=0
      count_lcheckf=0
      count_lfitstats=0
      count_lfixederrore=0
      count_lfixederrorf=0
      count_lompmkl=0
      count_lnormnodes=0
      count_restrictw=0
      count_fitmode=0
      count_lanalyzeerror=0
      count_lnwweights=0
      count_lnwweightse=0
      count_scmin_short_atomic=0  
      count_scmax_short_atomic=0
      count_scmin_short_pair=0  
      count_scmax_short_pair=0
      count_scmin_elec=0     
      count_scmax_elec=0    
      count_shuffle_weights_short_atomic=0    
      count_luseoldscaling=0    
      count_lprecond=0    
      count_linionly=0    
      count_noisee=0    
      count_noisef=0    
      count_noiseq=0    
      count_lprintconv=0
      count_lprintmad=0
      count_lfgroupbystruct=0
      count_lqgroupbystruct=0
      count_cutoff_type=0
      count_cutoff_alpha =0
      count_lmixpoints=0
      count_lscreen=0
      count_lsilent=0
      count_lpreparemd=0
      count_lpearson_correlation=0
      count_lweightanalysis=0
      count_lenvironmentanalysis=0
      count_lfindcontradictions=0
      count_lenforcemaxnumneighborsatomic=0
      count_lmd=0
      count_ldynforcegroup=0
      count_ldetect_saturation=0
      count_ldataclustering=0
      count_lanalyzecomposition=0
      count_lprintdateandtime=0
      count_lenableontheflyinput=0
      count_luseedkalman=0 !! MG: ED-Kalman counter
      count_ledforcesv2=0 !!  MG: Second variant of ED-Kalman force fitting
      count_lcheckinputforces=0
      count_lionforcesonly=0
      count_lprintforcecomponents=0
      count_useipi=0
!!
      return
      end 
