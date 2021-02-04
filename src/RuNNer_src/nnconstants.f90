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
      module nnconstants 
      implicit none

      real*8 pi 
      real*8 rad2deg 
      real*8 ha2ev

      contains
      subroutine get_nnconstants()
      use mpi_mod
      implicit none
!!
      pi        = 4.d0*datan(1.d0)
      rad2deg   = 180.d0/pi
      ha2ev     = 27.211d0

      end subroutine get_nnconstants 
      end module nnconstants 
