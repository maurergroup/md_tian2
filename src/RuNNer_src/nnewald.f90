!######################################################################
! This routine is part of
! RuNNer - RuNNer Neural Network Energy Representation
! (c) 2008-2020 Prof. Dr. Joerg Behler 
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
      module nnewald 
!!
      implicit none
!!
      integer, dimension(:)  , allocatable :: num_layers_elec
      integer, dimension(:,:), allocatable :: nodes_elec
      integer, dimension(:,:), allocatable :: windex_elec
      integer, dimension(:)  , allocatable :: num_weights_elec
      integer, dimension(:)  , allocatable :: num_funcvalues_elec
      integer maxnodes_elec

      real*8, dimension(:,:)   , allocatable :: weights_elec
      real*8, dimension(:,:,:) , allocatable :: symfunction_elec_list
      real*8 scmin_elec
      real*8 scmax_elec

      character*1, dimension(:,:,:), allocatable :: actfunc_elec

      end module nnewald 

