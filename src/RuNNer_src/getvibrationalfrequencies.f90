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
!! - getshortatomic.f90
!! Added by Emir

subroutine getvibrationalfrequencies(hessian_local,natoms,frequencies_local,lmodes)

      implicit none

      integer natoms ! in
      integer i,ii

      real*8  eig(3*natoms)
      real*8  frequencies_local(3*natoms) ! in/out
      real*8  hessian_local(3*natoms,3*natoms) ! in/out

      logical lmodes ! in

      !! calculate eigenvalues and eigenvectors of Hessian matrix for a given structure
      !! if lmodes is true, hessian_local is converted into the eigenvector matrix
      call dia_hessian(hessian_local,eig,natoms,lmodes)

      !! convert these eigenvalues to frequencies
      call convert_frequencies(eig,natoms,frequencies_local)


end subroutine getvibrationalfrequencies
