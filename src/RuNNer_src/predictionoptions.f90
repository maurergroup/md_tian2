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
      module predictionoptions 

      implicit none

      integer enforcetotcharge
      integer nn_type_vdw 

      real*8 vdw_screening(2)
      real*8, dimension(:,:) , allocatable :: vdw_coefficients
      real*8 melem(108) ! atomic masses of elements (required in normal mode calculation)

      logical ldoforces
      logical ldohessian
      logical lwritehessian ! Emir
      logical lcalculatefrequencies !Emir
      logical lcalculatenormalmodes !Emir
      logical lpreparemd
      logical lprintforcecomponents
      logical lwritesymfunctions
      logical lvdw

      end module predictionoptions 

