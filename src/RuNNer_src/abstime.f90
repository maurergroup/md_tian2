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

!! called by many subroutines 

        subroutine abstime(timetot,dayold)

        implicit none

        real*8 timetot             ! out
        real*8 seconds             ! internal
        integer day                ! internal
        integer dayold             ! out
        integer nhours             ! internal
        integer nminutes           ! internal
        integer month              ! internal
!!
        character*10 date,time
!!
!! obviously we cannot get a better time resolution than 0.001 s this way
        call date_and_time(date,time)
!!        write(*,*)date,time
!!
        read(time,'(i2,i2,f6.3)')nhours,nminutes,seconds
!!
        timetot=3600.d0*dble(nhours)+60.d0*dble(nminutes)+seconds
!!
        read(date,'(6x,i2)') day
        read(date,'(4x,i2)') month 
!!
        if(dayold.ne.0)then
          if(day.gt.dayold)then
            timetot=timetot+3600.d0*24.d0*dble(day-dayold)
          elseif(day.lt.dayold)then
            if(month.eq.1)then
              timetot=timetot+3600.d0*24.d0*dble(day+31-dayold)
            elseif(month.eq.2)then
              timetot=timetot+3600.d0*24.d0*dble(day+31-dayold)
            elseif(month.eq.3)then
              timetot=timetot+3600.d0*24.d0*dble(day+28-dayold)
            elseif(month.eq.4)then
              timetot=timetot+3600.d0*24.d0*dble(day+31-dayold)
            elseif(month.eq.5)then
              timetot=timetot+3600.d0*24.d0*dble(day+30-dayold)
            elseif(month.eq.6)then
              timetot=timetot+3600.d0*24.d0*dble(day+31-dayold)
            elseif(month.eq.7)then
              timetot=timetot+3600.d0*24.d0*dble(day+30-dayold)
            elseif(month.eq.8)then
              timetot=timetot+3600.d0*24.d0*dble(day+31-dayold)
            elseif(month.eq.9)then
              timetot=timetot+3600.d0*24.d0*dble(day+31-dayold)
            elseif(month.eq.10)then
              timetot=timetot+3600.d0*24.d0*dble(day+30-dayold)
            elseif(month.eq.11)then
              timetot=timetot+3600.d0*24.d0*dble(day+31-dayold)
            elseif(month.eq.12)then
              timetot=timetot+3600.d0*24.d0*dble(day+30-dayold)
            else
              write(*,*)'Error: unknown month??? ',month
            endif
          else ! no overnight job
          endif
        endif
!!
        dayold=day
!!
        return
        end
