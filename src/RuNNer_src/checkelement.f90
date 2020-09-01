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
!! called by:
!! - readinput.f90
!!
      subroutine checkelement(elementtemp)
!!
      use fileunits
      use globaloptions
!!
      implicit none
!!
      integer i1            ! internal
!!
      character*2 elementtemp    ! in
!!
      logical lfound             ! internal
!!
      lfound=.false.
!!
      do i1=1,nelem
        if(elementtemp.eq.element(i1))then
          lfound=.true.
        endif
      enddo
!!
      if(.not.lfound) then
        write(ounit,*)'Error: illegal element specified in some keyword in input.nn: ',elementtemp
        stop !'
      endif
!!
      return
      end
