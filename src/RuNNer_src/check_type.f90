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

! This checks whether a given SF is radial or angular,
! and returns 0 (radial) or 1 (angular) - Emir

subroutine check_type(ftype,symtype)

    implicit none

    integer ftype ! in
    integer symtype ! out

    if(ftype.eq.1.or.ftype.eq.2.or.ftype.eq.4)then
        symtype = 0 ! radial
    else if(ftype.eq.3.or.ftype.eq.7.or.ftype.eq.8.or.ftype.eq.9)then
        symtype = 1 ! angular
    end if

end subroutine check_type
