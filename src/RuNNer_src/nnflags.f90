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
      module nnflags 
      implicit none

!! mode = runner mode (1=calculate symfunctions, 2=fitting, 3=prediction)
      integer mode

!! nn_type_short: type of neural network: nn_type_short=1 => atom-based NN, nn_type_short=2 => pair-based NN
      integer nn_type_short
!! determine if separate charge NN is used (nn_type_elec=1) 
!! or if charges are additional NN output (nn_type_elec=2)
!! of if fixed charges are used instead of NN output (nn_type_elec=3)
      integer nn_type_elec

      logical lshort
!! calculate electrostatic energy (and forces, stress etc.) for given charges (according to keyword nn_type_elec) 
      logical lelec
      integer originatom_id
      integer zatom_id

      contains
      subroutine distribute_nnflags()
      use mpi_mod
      implicit none
!!
      call mpi_bcast(mode,1,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(nn_type_short,1,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(nn_type_elec,1,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(mode,1,mpi_integer,0,mpi_comm_world,mpierror)
      call mpi_bcast(lelec,1,mpi_logical,0,mpi_comm_world,mpierror)
      call mpi_bcast(lshort,1,mpi_logical,0,mpi_comm_world,mpierror)

      end subroutine distribute_nnflags
      end module nnflags 
