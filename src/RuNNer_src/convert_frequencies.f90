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

!! Added by Emir

subroutine convert_frequencies(eig,natoms,frequencies)

    implicit none

    integer natoms ! in
    integer i

    real*8  pi,c,ha2j,bohr2m,u2kg ! internal
    real*8  fac1,fac2 ! internal
    real*8  eig(3*natoms) ! in
    real*8  frequencies(3*natoms) ! out

    pi = 3.14159265
    c = 29979245800.d0
    ha2j = 4.3597482e-18
    bohr2m = 5.29177249e-11
    u2kg = 1.660538782e-27

    fac1 = 1/(2*pi*c) ! check me
    fac2 = ha2j/(bohr2m*bohr2m*u2kg)

    !! convert eigenvalues into frequencies
    do i = 1,3*natoms
        if(eig(i).lt.0)then
            frequencies(i) = -1.d0*(fac1 * dsqrt(-1.d0*eig(i)*fac2))
        else
            frequencies(i) = fac1 * dsqrt(eig(i)*fac2)
        end if
    end do



end subroutine convert_frequencies

