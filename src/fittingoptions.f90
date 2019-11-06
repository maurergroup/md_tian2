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
      module fittingoptions 

      implicit none

      integer ngrowth
      integer growthstep
      integer nenergygroup
      integer nforcegroup
      integer nchargegroup
      integer optmodee
      integer optmodef
      integer optmodeq
      integer nepochs
      integer elemupdate
      integer fitmode
      integer fitting_unit
      integer iwriteweight
      integer dynforcegroup_start
      integer dynforcegroup_step

!! number of point with largest training error in final epoch for short range energy
      integer imaxerror_eshort_train 

!! number of point with largest test error in final epoch for short range energy
      integer imaxerror_eshort_test

!! number of point with largest training error in final epoch for electrostatic energy
      integer imaxerror_elec_train 

!! number of point with largest test error in final epoch for electrostatic energy
      integer imaxerror_elec_test  

!! number of point with largest training error in final epoch for total energy
      integer imaxerror_etot_train 

!! number of point with largest test error in final epoch for total energy
      integer imaxerror_etot_test  

      integer nshuffle_weights_short_atomic  
      real*8 shuffle_weights_short_atomic  
      logical lshuffle_weights_short_atomic  

      real*8 analyze_error_energy_step
      real*8 analyze_error_force_step
      real*8 analyze_error_charge_step
      real*8 forcernd
      real*8 chargernd
      real*8 energyrnd
      real*8 worste
      real*8 worstf
      real*8 worstq
      real*8 fixederrore
      real*8 fixederrorf
      real*8 noisee
      real*8 noisef
      real*8 noiseq
      real*8 steepeststepe
      real*8 steepeststepf
      real*8 steepeststepq
      real*8 kalman_dampe
      real*8 kalman_dampf
      real*8 kalman_dampq
      real*8 scalefactorf
      real*8 scalefactorq
      real*8 kalmanthreshold
      real*8 kalmanthresholdf
      real*8 kalmanthresholde
      real*8 kalmanthresholdc
      real*8 kalmannue
      real*8 kalmannuee
      real*8 kalmannuec
      real*8 dampw
      real*8 weights_min
      real*8 weights_max
      real*8 biasweights_min
      real*8 biasweights_max
      real*8 weightse_min
      real*8 weightse_max
      real*8 restrictw
      real*8 maxforce
      real*8 maxenergy
      real*8 kalmanlambdac
      real*8, dimension(:)  , allocatable :: kalmanlambda
      real*8, dimension(:)  , allocatable :: kalmanlambdap
      real*8, dimension(:)  , allocatable :: kalmanlambdae
      real*8 deltagthres 
      real*8 deltafthres 

!! distance thresholds for data clustering
      real*8 dataclusteringthreshold1 
      real*8 dataclusteringthreshold2 

!! largest error of short range energy of a single point in the training set in the last epoch
      real*8 maxerror_eshort_train                       

!! largest error of short range energy of a single point in the test set in the last epoch
      real*8 maxerror_eshort_test                             

!! largest error of electrostatic energy of a single point in the training set in the last epoch
      real*8 maxerror_elec_train                       

!! largest error of electrostatic energy of a single point in the test set in the last epoch
      real*8 maxerror_elec_test                             

!! largest error of total energy of a single point in the training set in the last epoch
      real*8 maxerror_etot_train                       

!! largest error of total energy of a single point in the test set in the last epoch
      real*8 maxerror_etot_test                             



      logical lwritetmpweights
      logical lsavekalman
      logical lrestkalman
      logical lgrowth
      logical luseworste
      logical luseworstf
      logical luseworstq
      logical lfixederrore
      logical lfixederrorf 
      logical lfgroupbystruct
      logical lqgroupbystruct
      logical lmixpoints
      logical lupdatebyelement
      logical ldampw
      logical ljointefupdate
      logical lsysweights
      logical lsysweightse
      logical luseoldweightsshort
      logical luseoldweightscharge
      logical luseoldscaling
      logical linionly
      logical lprecond
      logical lprintconv
      logical lprintmad
      logical lfinalforce
      logical lnwweights
      logical lnwweightse
      logical lglobalfit
      logical lanalyzeerror
      logical lwritetrainpoints
      logical lwritetraincharges
      logical lwritetrainforces
      logical lfixweights
      logical lrepeate
      logical lsepkalman
      logical lfitstats
      logical lchargeconstraint
      logical lresetkalman 
      logical lseparatebiasini
      logical lfindcontradictions
      logical ldynforcegroup
      logical ldataclustering
      logical lanalyzecomposition
      logical lprintdateandtime
      logical lenableontheflyinput
      logical lionforcesonly
      logical lnoisematrix !! modifed by kenko

!! MG: select element-decoupled Kalman filter in mode 2
      logical luseedkalman
!! MG: select 2nd variant of ED-GEKF force fitting
      logical ledforcesv2 

!! print analysis of weights in mode 2 
      logical lweightanalysis 

!! filename array for short range weights in mode 2
      character*20, dimension (:), allocatable :: filenamews

!! filename array for electrostatic weights in mode 2
      character*20, dimension (:), allocatable :: filenamewe

!! filename array for pair weights in mode 2
      character*23, dimension (:), allocatable :: filenamewp


      end module fittingoptions 

