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
!! determine elementindex and pairindex arrays

!! called by
!! - readinput.f90
!! - getdimensions.f90
!!
      subroutine sortelements()
!!
      use fileunits
      use globaloptions
      use nnflags
!!
      implicit none
!!
      integer counter                          ! internal
      integer ztemp                            ! internal
      integer i1,i2,i3                         ! internal
!!
      character*2 elementtemp                  ! internal
!!
!!
      elementindex(:)=0
!!
      if(nelem.gt.1)then
 10     continue
        do i1=1,nelem-1
          if(nucelem(i1).gt.nucelem(i1+1))then
            ztemp=nucelem(i1)
            elementtemp=element(i1)
            nucelem(i1)=nucelem(i1+1)
            element(i1)=element(i1+1)
            nucelem(i1+1)=ztemp
            element(i1+1)=elementtemp
            goto 10
          endif
        enddo
      endif
!!
      do i1=1,102
        do i2=1,nelem
          if(nucelem(i2).eq.i1)then
            elementindex(i1)=i2
          endif 
        enddo
      enddo
!!
      pairindex(:,:)=0
      counter       =0
      do i1=1,nelem
        do i2=1,nelem
         if(nucelem(i2).ge.nucelem(i1))then
            counter = counter +1
            pairindex(nucelem(i1),nucelem(i2)) = counter
            pairindex(nucelem(i2),nucelem(i1)) = counter
         endif
        enddo
      enddo
!!
!!
      return
      end
