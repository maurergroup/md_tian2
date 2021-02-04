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
!! - prediction.f90
!! - predictionpair.f90
!! - checkonestructure.f90
!!
      subroutine getvolume(lattice,volume)
!!
      use fileunits
!!
      implicit none
!!
      real*8 lattice(3,3)              ! in
      real*8 volume                    ! out
      real*8 tempvec(3)                ! internal
!!
      volume=0.0d0
!!
      tempvec(1)=lattice(1,2)*lattice(2,3)-lattice(1,3)*lattice(2,2)
      tempvec(2)=lattice(1,3)*lattice(2,1)-lattice(1,1)*lattice(2,3)
      tempvec(3)=lattice(1,1)*lattice(2,2)-lattice(1,2)*lattice(2,1)
      volume=tempvec(1)*lattice(3,1)&
            +tempvec(2)*lattice(3,2)&
            +tempvec(3)*lattice(3,3)
      volume=abs(volume)
!!
      if(volume.lt.0.0d0)then
        write(ounit,*)'Error: volume<0 ',volume
        stop
      endif
!!
      return
      end
