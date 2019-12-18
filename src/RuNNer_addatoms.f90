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
!!
!! called by:
!! - prediction.f90
!! - predictionpair.f90
!!
      subroutine addatoms(num_atoms,&
             zelem,num_atoms_element,&
             totalenergy,atomenergy)
!!
      use globaloptions
!!
      implicit none
!!
      integer num_atoms                            ! in
      integer num_atoms_element(nelem)             ! in  
      integer zelem(max_num_atoms)                 ! in
      integer i1                                   ! internal
!!
      real*8 totalenergy                           ! in/out
      real*8 atomenergy(max_num_atoms)             ! in/out
!!
!!
!! add atomic energies to total energy
        do i1=1,nelem
          totalenergy=totalenergy&
          +dble(num_atoms_element(i1))*atomrefenergies(i1)
        enddo
!!
!!
!! add atomic energies to atomic energy contributions
        do i1=1,num_atoms ! loop over all atoms 
          atomenergy(i1)=atomenergy(i1)&
            +atomrefenergies(elementindex(zelem(i1)))
        enddo
!!
      return
      end
