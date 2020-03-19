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
      subroutine mpistart()

      implicit none


      return
      end

!#####################################################3

      subroutine mpi_init(mpierror)
      implicit none

      integer mpierror

      mpierror = 0

      return
      end

!#####################################################3
 
      subroutine mpi_finalize(mpierror)

      implicit none

      integer mpierror

      mpierror =  0

      return
      end

!#####################################################3

      subroutine mpi_comm_size(mpi_comm_world,mpisize,mpierror)

      implicit none

      integer mpi_comm_world
      integer mpisize
      integer mpierror

      mpisize = 1
      mpierror= 0

      return
      end

!#####################################################3

      subroutine mpi_comm_rank(mpi_comm_world,mpirank,mpierror)

      implicit none

      integer mpi_comm_world
      integer mpirank
      integer mpierror

      mpirank = 0
      mpierror= 0

      return
      end

!#####################################################3

      subroutine mpi_barrier(mpi_comm_world,mpierror)

      implicit none

      integer mpi_comm_world
      integer mpierror

      mpierror= 0

      return
      end

!#####################################################3

      subroutine mpi_allreduce(rin,rout,ndim,mpi_type,&
        method,mpi_comm_world,mpierror)

      implicit none

      integer mpi_comm_world
      integer mpierror
      integer ndim
      integer mpi_type
      integer method

      real*8 rin(ndim)
      real*8 rout(ndim)

      mpierror=0

      return
      end

!#####################################################3

      subroutine mpi_bcast(array,ndim,mpi_type,&
        mpirank,mpi_comm_world,mpierror)

      implicit none

      integer mpierror
      integer ndim
      integer mpi_type
      integer mpirank
      integer mpi_comm_world

      real*8 array(ndim)

      mpierror=0


      return
      end

!#####################################################3

      subroutine mpi_scatterv(rin,ndim,ndummy,mpi_type,&
        rout,ndummy2,mpi_type2,ndummy3,mpi_comm_world,mpierror)

      implicit none

      integer mpi_comm_world
      integer mpierror
      integer ndim
      integer ndummy
      integer ndummy2
      integer ndummy3
      integer mpi_type
      integer mpi_type2

      real*8 rin(ndim)
      real*8 rout(ndim)

      mpierror=0

      return
      end

!#####################################################3


      subroutine mpi_gatherv(rin,ndim,mpi_type,&
        rout,ndummy,ndummy2,mpi_type2,ndummy3,mpi_comm_world,mpierror)

      implicit none

      integer mpi_comm_world
      integer mpierror
      integer ndim
      integer ndummy
      integer ndummy2
      integer ndummy3
      integer mpi_type
      integer mpi_type2

      real*8 rin(ndim)
      real*8 rout(ndim)

      mpierror=0

      return
      end

!#####################################################3





