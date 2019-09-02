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
      module symfunctions 
!!
      implicit none
!!
      integer, dimension(:,:), allocatable :: function_type_short_atomic
      integer, dimension(:,:), allocatable :: function_type_elec
      integer, dimension(:,:), allocatable :: function_type_short_pair

      integer, dimension(:,:,:), allocatable :: symelement_short_atomic
      integer, dimension(:,:,:), allocatable :: symelement_elec
      integer, dimension(:,:,:), allocatable :: symelement_short_pair

      real*8, dimension(:,:) , allocatable :: funccutoff_short_atomic
      real*8, dimension(:,:) , allocatable :: funccutoff_elec
      real*8, dimension(:,:) , allocatable :: funccutoff_short_pair

      real*8, dimension(:,:) , allocatable :: eta_short_atomic
      real*8, dimension(:,:) , allocatable :: eta_elec
      real*8, dimension(:,:) , allocatable :: eta_short_pair

      real*8, dimension(:,:) , allocatable :: zeta_short_atomic
      real*8, dimension(:,:) , allocatable :: zeta_elec
      real*8, dimension(:,:) , allocatable :: zeta_short_pair

      real*8, dimension(:,:) , allocatable :: lambda_short_atomic
      real*8, dimension(:,:) , allocatable :: lambda_elec
      real*8, dimension(:,:) , allocatable :: lambda_short_pair

      real*8, dimension(:,:) , allocatable :: rshift_short_atomic
      real*8, dimension(:,:) , allocatable :: rshift_elec
      real*8, dimension(:,:) , allocatable :: rshift_short_pair

      real*8 maxcutoff_short_atomic 
      real*8 maxcutoff_elec   
      real*8 maxcutoff_short_pair

      contains 

      subroutine allocatesymfunctions()
!!
      use nnflags 
      use globaloptions 
!!
      implicit none
!!
      if(lshort.and.(nn_type_short.eq.1))then
        allocate(function_type_short_atomic(maxnum_funcvalues_short_atomic,nelem))
        function_type_short_atomic(:,:)=0
        allocate(symelement_short_atomic(maxnum_funcvalues_short_atomic,2,nelem))
        symelement_short_atomic(:,:,:)=0
        allocate(funccutoff_short_atomic(maxnum_funcvalues_short_atomic,nelem))
        funccutoff_short_atomic(:,:)=0.0d0
        allocate(eta_short_atomic(maxnum_funcvalues_short_atomic,nelem))
        eta_short_atomic(:,:)=0.0d0
        allocate(zeta_short_atomic(maxnum_funcvalues_short_atomic,nelem))
        zeta_short_atomic(:,:)=0.0d0
        allocate(lambda_short_atomic(maxnum_funcvalues_short_atomic,nelem))
        lambda_short_atomic(:,:)=0.0d0
        allocate(rshift_short_atomic(maxnum_funcvalues_short_atomic,nelem))
        rshift_short_atomic(:,:)=0.0d0
      endif
!!
      if(lshort.and.(nn_type_short.eq.2))then
        allocate(function_type_short_pair(maxnum_funcvalues_short_pair,npairs))
        function_type_short_pair(:,:)=0
        allocate(symelement_short_pair(maxnum_funcvalues_short_pair,2,npairs))
        symelement_short_pair(:,:,:)=0
        allocate(funccutoff_short_pair(maxnum_funcvalues_short_pair,npairs))
        funccutoff_short_pair(:,:)=0.0d0
        allocate(eta_short_pair(maxnum_funcvalues_short_pair,npairs))
        eta_short_pair(:,:)=0.0d0
        allocate(zeta_short_pair(maxnum_funcvalues_short_pair,npairs))
        zeta_short_pair(:,:)=0.0d0
        allocate(lambda_short_pair(maxnum_funcvalues_short_pair,npairs))
        lambda_short_pair(:,:)=0.0d0
        allocate(rshift_short_pair(maxnum_funcvalues_short_pair,npairs))
        rshift_short_pair(:,:)=0.0d0
      endif
!!
      if(lelec.and.(nn_type_elec.eq.1))then
        allocate(function_type_elec(maxnum_funcvalues_elec,nelem))
        function_type_elec(:,:)=0
        allocate(symelement_elec(maxnum_funcvalues_elec,2,nelem))
        symelement_elec(:,:,:)=0
        allocate(funccutoff_elec(maxnum_funcvalues_elec,nelem))
        funccutoff_elec(:,:)=0.0d0
        allocate(eta_elec(maxnum_funcvalues_elec,nelem))
        eta_elec(:,:)=0.0d0
        allocate(zeta_elec(maxnum_funcvalues_elec,nelem))
        zeta_elec(:,:)=0.0d0
        allocate(lambda_elec(maxnum_funcvalues_elec,nelem))
        lambda_elec(:,:)=0.0d0
        allocate(rshift_elec(maxnum_funcvalues_elec,nelem))
        rshift_elec(:,:)=0.0d0
      endif
!!
      end subroutine allocatesymfunctions
 
      subroutine deallocatesymfunctions()
!!
      use nnflags 
      use globaloptions
!!
      implicit none

      if(lshort.and.(nn_type_short.eq.1))then
        deallocate(function_type_short_atomic)
        deallocate(symelement_short_atomic)
        deallocate(funccutoff_short_atomic)
        deallocate(eta_short_atomic)
        deallocate(zeta_short_atomic)
        deallocate(lambda_short_atomic)
        deallocate(rshift_short_atomic)
      endif
!!
      if(lshort.and.(nn_type_short.eq.2))then
        deallocate(function_type_short_pair)
        deallocate(symelement_short_pair)
        deallocate(funccutoff_short_pair)
        deallocate(eta_short_pair)
        deallocate(zeta_short_pair)
        deallocate(lambda_short_pair)
        deallocate(rshift_short_pair)
      endif

      if(lelec.and.(nn_type_elec.eq.1))then
        deallocate(function_type_elec)
        deallocate(symelement_elec)
        deallocate(funccutoff_elec)
        deallocate(eta_elec)
        deallocate(zeta_elec)
        deallocate(lambda_elec)
        deallocate(rshift_elec)
      endif
!!
      end subroutine deallocatesymfunctions

      subroutine distribute_symfunctions()
!!
      use mpi_mod
      use nnflags 
      use globaloptions
!!
      implicit none
      integer ndim
!!
      if(lshort.and.(nn_type_short.eq.1))then
        call mpi_bcast(function_type_short_atomic,&
          maxnum_funcvalues_short_atomic*nelem,mpi_integer,0,mpi_comm_world,mpierror)
        call mpi_bcast(symelement_short_atomic,&
          maxnum_funcvalues_short_atomic*nelem*2,mpi_integer,0,mpi_comm_world,mpierror)
        call mpi_bcast(funccutoff_short_atomic,&
          maxnum_funcvalues_short_atomic*nelem,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(eta_short_atomic,&
          maxnum_funcvalues_short_atomic*nelem,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(zeta_short_atomic,&
          maxnum_funcvalues_short_atomic*nelem,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(lambda_short_atomic,&
          maxnum_funcvalues_short_atomic*nelem,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(rshift_short_atomic,&
          maxnum_funcvalues_short_atomic*nelem,mpi_real8,0,mpi_comm_world,mpierror)
      endif
!!
      if(lshort.and.(nn_type_short.eq.2))then
        call mpi_bcast(function_type_short_pair,&
          maxnum_funcvalues_short_pair*npairs,mpi_integer,0,mpi_comm_world,mpierror)
        call mpi_bcast(symelement_short_pair,&
          maxnum_funcvalues_short_pair*npairs*2,mpi_integer,0,mpi_comm_world,mpierror)
        call mpi_bcast(funccutoff_short_pair,&
          maxnum_funcvalues_short_pair*npairs,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(eta_short_pair,&
          maxnum_funcvalues_short_pair*npairs,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(zeta_short_pair,&
          maxnum_funcvalues_short_pair*npairs,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(lambda_short_pair,&
          maxnum_funcvalues_short_pair*npairs,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(rshift_short_pair,&
          maxnum_funcvalues_short_pair*npairs,mpi_real8,0,mpi_comm_world,mpierror)
      endif
!!
      if(lelec.and.(nn_type_elec.eq.1))then
        call mpi_bcast(function_type_elec,&
          maxnum_funcvalues_elec*nelem,mpi_integer,0,mpi_comm_world,mpierror)
        call mpi_bcast(symelement_elec,&
          maxnum_funcvalues_elec*nelem*2,mpi_integer,0,mpi_comm_world,mpierror)
        call mpi_bcast(funccutoff_elec,&
          maxnum_funcvalues_elec*nelem,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(eta_elec,&
          maxnum_funcvalues_elec*nelem,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(zeta_elec,&
          maxnum_funcvalues_elec*nelem,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(lambda_elec,&
          maxnum_funcvalues_elec*nelem,mpi_real8,0,mpi_comm_world,mpierror)
        call mpi_bcast(rshift_elec,&
          maxnum_funcvalues_elec*nelem,mpi_real8,0,mpi_comm_world,mpierror)
      endif
!!
      end subroutine distribute_symfunctions

      end module symfunctions 
