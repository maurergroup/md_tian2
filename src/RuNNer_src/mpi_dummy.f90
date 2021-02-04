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
      module mpi_mod

      implicit none

      integer mpierror
      integer mpirank
      integer mpisize
      integer mpi_comm_world
      integer mpi_sum
      integer mpi_double_precision
      integer mpi_lor
      integer mpi_integer
      integer mpi_real8
      integer mpi_character
      integer mpi_logical
      integer mpi_in_place

      end module mpi_mod
