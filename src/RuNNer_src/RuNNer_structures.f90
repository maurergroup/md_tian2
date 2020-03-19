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
      module structures 
!!
      use globaloptions
!!
      implicit none

      integer, dimension(:)      , allocatable :: num_atoms_list 
      integer, dimension(:)      , allocatable :: num_pairs_list
!!
      integer, dimension(:,:)    , allocatable :: zelem_list 
      integer, dimension(:,:,:)  , allocatable :: zelemp_list 
      integer, dimension(:,:,:)  , allocatable :: zelemtrip_list

      real*8, dimension(:,:,:)  , allocatable :: lattice_list 
      real*8, dimension(:,:,:)  , allocatable :: xyzstruct_list  
      real*8, dimension(:)      , allocatable :: totalcharge_list  
      real*8, dimension(:)      , allocatable :: totalenergy_list  
      real*8, dimension(:)      , allocatable :: shortenergy_list  
      real*8, dimension(:)      , allocatable :: elecenergy_list  
      real*8, dimension(:,:,:)  , allocatable :: totalforce_list  
      real*8, dimension(:,:,:)  , allocatable :: totforce_list  
      real*8, dimension(:,:,:)  , allocatable :: elecforce_list  
      real*8, dimension(:,:,:)  , allocatable :: shortforce_list  
      real*8, dimension(:,:)    , allocatable :: atomcharge_list  
      real*8, dimension(:,:)    , allocatable :: atomenergy_list  

      logical, dimension(:)     , allocatable :: lperiodic_list  

      character*2, dimension(:,:) , allocatable :: elementsymbol_list  

!!      integer num_atoms_list(nblock)
!!      integer num_pairs_list(nblock)
!!      integer zelem_list(nblock,max_num_atoms)
!!      integer zelemp_list(2,nblock,max_num_pairs)

!!      real*8 lattice_list(3,3,nblock)
!!      real*8 xyzstruct_list(3,max_num_atoms,nblock)
!!      real*8 totalcharge_list(nblock)
!!      real*8 totalenergy_list(nblock)
!!      real*8 shortenergy_list(nblock)
!!      real*8 elecenergy_list(nblock)
!!      real*8 totalforce_list(3,max_num_atoms,nblock)
!!      real*8 elecforce_list(3,max_num_atoms,nblock)
!!      real*8 totforce_list(3,max_num_atoms,nblock)
!!      real*8 atomcharge_list(nblock,max_num_atoms)
!!      real*8 atomenergy_list(nblock,max_num_atoms)

!!      logical lperiodic_list(nblock)

!!      character*2 elementsymbol_list(nblock,max_num_atoms)

      end module structures 

