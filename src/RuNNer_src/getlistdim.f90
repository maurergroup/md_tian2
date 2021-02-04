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
!! called by: 
!! - main.f90
!!
      subroutine getlistdim()
!!
      use nnflags
      use globaloptions
      use fileunits
!!
      implicit none
!!
      integer listdimtemp1
      integer listdimtemp2
      integer listdimtemp3
      integer listdimtemp4
!!
!! set a more clever value for listdim here
      listdim=0
      if(lshort.and.(nn_type_short.eq.1))then
        listdimtemp1=800*max_num_atoms
!! make sure here that in parallel cases neighbor arrays are not larger than necessary
        if(mode.eq.3)then
          listdimtemp1=200*nblock
        endif
        listdim=max(listdim,listdimtemp1)
      endif
      if(lshort.and.(nn_type_short.eq.2))then
        listdimtemp2=400*max_num_pairs
!! make sure here that in parallel cases neighbor arrays are not larger than necessary
        if(mode.eq.3)then
          listdimtemp2=200*nblock
        endif
        listdim=max(listdim,listdimtemp2)
      endif
!! separate electrostatic NN:
!!      if(lelec.and.(nn_type_elec.eq.1))then ! this does not work because also for fixed charges the Ewald sum needs neighbor lists
      if(lelec)then
        listdimtemp3=400*max_num_atoms
!! make sure here that in parallel cases neighbor arrays are not larger than necessary
        if(mode.eq.3)then
          listdimtemp3=200*nblock
        endif
        listdim=max(listdim,listdimtemp3)
      endif
!! Hamiltonian NN:
!! 
!! TODO: lsta can be reduced in size for prediction mode (at least)
!!    if(mode.eq.3)then
!!      lstadim=ceiling(dble(nblock) / dble(mpisize))
!!    else
!!      lstadim=max_num_atoms
!!    endif

      return
      end
