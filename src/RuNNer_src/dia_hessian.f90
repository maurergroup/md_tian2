!######################################################################
! This routine is part of
! RuNNer - RuNNer Neural Network Energy Representation
! (c) 2008-2018 Prof. Dr. Joerg Behler
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
!! - getvibrationalfrequencies.f90
!! calls the LAPACK diagonalization subroutine DSYEV
!! input:  hessian(n,n) to be diagonalized
!! output: orthonormal eigenvectors and
!! eigenvalues of a in ascending order

subroutine dia_hessian(hessian,eig,natoms,lmodes)

    implicit none

    integer n,l,inf
    integer natoms
    real*8  work((3*natoms)*(3+(3*natoms)/2)) ! ???
    real*8  hessian(3*natoms,3*natoms) !in
    real*8  eig(3*natoms) ! output
    logical lmodes !in

    character dsyev_mode

    if(lmodes)then
        dsyev_mode = 'V'
    else
        dsyev_mode = 'N'
    end if

    l = (3*natoms) * (3 + (3*natoms)/2)
    call dsyev(dsyev_mode,'U',3*natoms,hessian,3*natoms,eig,work,l,inf)

end subroutine

