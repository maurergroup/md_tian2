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
!! called by: 
!! - predict.f90
!!
      subroutine predictionshortatomic(&
        num_atoms,num_atoms_element,zelem,&
        lattice,xyzstruct,&
        minvalue_short_atomic,maxvalue_short_atomic,avvalue_short_atomic,&
        eshortmin,eshortmax,&
        nntotalenergy,nnshortforce,&
        nnatomenergy,nnshortenergy,nnstress_short,&
         atomenergysum,sens,lperiodic)
!!
      use mpi_mod
      use fileunits
      use predictionoptions 
      use nnflags
      use globaloptions
      use symfunctions 
      use timings
      use nnshort_atomic
!!
      implicit none
!!
      integer zelem(max_num_atoms)                    ! in
      integer num_atoms                               ! in
      integer num_atoms_element(nelem)                ! in 
      integer npoints                                 ! internal
      integer ncount                                  ! internal
      integer ndone                                   ! internal
      integer i1                                      ! internal
      integer natoms                                  ! internal
      integer, dimension(:), allocatable :: atomindex ! internal
      integer, dimension(:,:), allocatable :: neighboridx_short_atomic    ! internal
      integer, dimension(:,:), allocatable :: invneighboridx_short_atomic ! internal
      integer n_start,n_end                           ! internal
      integer, allocatable :: lsta(:,:)                     ! numbers of neighbors
      integer, allocatable :: lstc(:)                       ! identification of atom
      integer, allocatable :: lste(:)                       ! nuclear charge of atom
      integer, allocatable :: num_neighbors_short_atomic(:) ! internal    
      integer max_num_neighbors_short_atomic                ! internal
!!
      real*8 lattice(3,3)                                ! in
      real*8 xyzstruct(3,max_num_atoms)                  ! in
      real*8 minvalue_short_atomic(nelem,maxnum_funcvalues_short_atomic)           ! in
      real*8 maxvalue_short_atomic(nelem,maxnum_funcvalues_short_atomic)           ! in
      real*8 avvalue_short_atomic(nelem,maxnum_funcvalues_short_atomic)            ! in
!! DFT data (not necessarily provided in predicition mode)
!! data predicted the NN
      real*8 nntotalenergy                                              ! out 
      real*8 nnshortforce(3,max_num_atoms)                              ! out 
      real*8 nnatomenergy(max_num_atoms)                                ! out 
!! symmetry function parameters
      real*8 nnshortenergy                                              ! out 
      real*8 nnstress_short(3,3)                                        ! internal
      real*8 sens(nelem,maxnum_funcvalues_short_atomic)                              ! out  
      real*8 eshortmin                                                  ! in
      real*8 eshortmax                                                  ! in
      real*8 atomenergysum                                              ! out 
      real*8, allocatable :: lstb(:,:)                                  ! xyz and r_ij 
!!
      logical lperiodic                                                 ! in
      logical lextrapolation                                            ! internal
!!
!!======================================================================
!! Initialzations
!!======================================================================
!! initializations for NN data
      nntotalenergy          = 0.0d0
      nnatomenergy(:)        = 0.0d0
      nnstress_short(:,:)    = 0.0d0
      nnshortforce(:,:)      = 0.0d0
      nnshortenergy          = 0.0d0
      atomenergysum          = 0.0d0
!! initialization of timings
      timeshort              = 0.0d0 
      timeeshort             = 0.0d0 
      timefshort             = 0.0d0 
      timesshort             = 0.0d0 
      timeallocshort         = 0.0d0 
      timesymshort           = 0.0d0 
      timeextrapolationshort = 0.0d0 
      timescalesymshort      = 0.0d0 
      timescaledsfuncshort   = 0.0d0 
      timecomm1              = 0.0d0 
!! initialization of sensitivities
      sens(:,:)              = 0.0d0
!!
!!======================================================================
!! Start short range part
!!======================================================================
      if(lshort)then
        if(lfinetime)then
          dayshort=0
          call abstime(timeshortstart,dayshort)
        endif ! lfinetime
!!======================================================================
!! memory management strategy: keep only 'nblock' atoms in memory at once
!! => loop step by step over these blocks of atoms 
!! In parallel runs the block of atoms is further split among the processes and each process does natoms atoms.
!!======================================================================
!! initialize auxiliary counters for splitting of atoms
        ncount  = num_atoms ! total number of atoms in structure
        npoints = 0 ! number of atoms to be calculated in this loop step by all processes together
        ndone   = 0 ! number of atoms calculated in previous loops
!!
!! next block of npoints atoms 
 11     continue
        if(ncount.gt.nblock)then
          npoints=nblock
          ncount=ncount-nblock
        else
          npoints=ncount
          ncount=ncount-npoints
        endif
!!
!!======================================================================
!! Start short range energy part for this block of atoms 
!!======================================================================
!!
!!======================================================================
!! preparations for parallel runs 
!!======================================================================
        if((mpirank.eq.0).and.(.not.lmd))then
          write(ounit,*)'-------------------------------------------------------------'
          write(ounit,'(a,i8,a)')'This cycle calculates total ',npoints,' atoms'
        endif !'
!! each process now picks 'natoms' atoms from the total number of 'npoints' atoms of this block
        call mpifitdistribution(npoints,natoms,n_start,n_end)
!! adjust position by atoms already done:
        n_start=n_start+ndone
        n_end  =n_end  +ndone
        if(.not.lmd)then
          write(ounit,'(a,i6,a,i8,a,i8,a,i8)')&
            ' process ',mpirank,' calculates atoms ',n_start,' to ',n_end,' natoms: ',natoms
        endif
        call mpi_barrier(mpi_comm_world,mpierror)
!! determine the atomindex array '
        allocate(atomindex(natoms))
        do i1=1,natoms
          atomindex(i1)=n_start+i1-1
        enddo
!!======================================================================
!! end preparations for parallel runs 
!!======================================================================
!!
        if(lfinetime)then
          dayallocshort=0 
          call abstime(timeallocshortstart,dayallocshort)
        endif ! lfinetime
        allocate(lsta(2,max_num_atoms))
        allocate(lstc(listdim))
        allocate(lste(listdim))
        allocate(lstb(listdim,4))
        allocate(num_neighbors_short_atomic(num_atoms))
        if(lfinetime)then
          call abstime(timeallocshortend,dayallocshort)
          timeallocshort=timeallocshort+timeallocshortend-timeallocshortstart
        endif ! lfinetime
!!
!!======================================================================
!! get num_neighbors_short_atomic and max_num_neighbors_short_atomic for short range part
!!======================================================================
        call getneighborsatomic_para(n_start,n_end,natoms,&
          num_atoms,num_neighbors_short_atomic,zelem,&
          max_num_neighbors_short_atomic,&
          lsta,lstc,lste,&
          maxcutoff_short_atomic,lattice,xyzstruct,lstb,&
          lperiodic)
!!
!!======================================================================
!! get neighboridx_short_atomic and invneighboridx_short_atomic for short range part 
!!======================================================================
        allocate(neighboridx_short_atomic(natoms,0:max_num_neighbors_short_atomic))  
        allocate(invneighboridx_short_atomic(natoms,max_num_atoms))  
        call getneighboridxatomic_para(n_start,natoms,listdim,&
          max_num_atoms,max_num_neighbors_short_atomic,&
          lsta,lstc,neighboridx_short_atomic,invneighboridx_short_atomic)
!!
!!======================================================================
!! calculation of short range atomic energies (nnatomenergy), 
!! forces (nnshortforce), and stress (nnstress_short)
!!======================================================================
        call getshortatomic(n_start,natoms,atomindex,&
          max_num_neighbors_short_atomic,num_atoms,&
          invneighboridx_short_atomic,num_neighbors_short_atomic,&
          neighboridx_short_atomic,zelem,&
          lsta,lstc,lste,lstb,xyzstruct,&
          sens,nnshortforce,nnstress_short,minvalue_short_atomic,&
          maxvalue_short_atomic,avvalue_short_atomic,&
          scmin_short_atomic,scmax_short_atomic,nnatomenergy,lextrapolation,lperiodic)
!!
!!======================================================================
!! add atomic energy constributions of this block of atoms to nntotalenergy
!!======================================================================
        do i1=1,natoms
          nntotalenergy=nntotalenergy+nnatomenergy(atomindex(i1))
        enddo
!!
!!======================================================================
!! deallocate short range arrays depending on natoms of this block of atoms 
!!======================================================================
        deallocate(atomindex)
        deallocate(lsta)
        deallocate(lstc)
        deallocate(lste)
        deallocate(lstb)
        deallocate(neighboridx_short_atomic)  
        deallocate(invneighboridx_short_atomic)  
        deallocate(num_neighbors_short_atomic)
!!
!! update the number of finished atoms:
        ndone=ndone+npoints
!!======================================================================
!! End short range energy part for this block of atoms 
!!======================================================================
!! if there are atoms left to be done go to back and do next block of atoms
        if(ncount.gt.0) goto 11
!!
!! calculate short range time
        if(lfinetime)then
          call abstime(timeshortend,dayshort)
          timeshort=timeshort+timeshortend-timeshortstart
        endif ! lfinetime
!!
!!======================================================================
!! combine the short range arrays of all processes
!!======================================================================
        if(lfinetime)then
          daycomm1=0
          call abstime(timecomm1start,daycomm1)
        endif ! lfinetime
!!
        call mpi_allreduce(mpi_in_place,nnatomenergy,max_num_atoms,&
          mpi_real8,mpi_sum,mpi_comm_world,mpierror)
        call mpi_allreduce(mpi_in_place,nntotalenergy,1,&
          mpi_real8,mpi_sum,mpi_comm_world,mpierror)
        if(ldoforces)then
          call mpi_allreduce(mpi_in_place,nnshortforce,max_num_atoms*3,&
            mpi_real8,mpi_sum,mpi_comm_world,mpierror)
        endif
        if(lsens)then
          call mpi_allreduce(mpi_in_place,sens,nelem*maxnum_funcvalues_short_atomic,&
            mpi_real8,mpi_sum,mpi_comm_world,mpierror)
        endif
        if(ldostress.and.lperiodic)then
          call mpi_allreduce(mpi_in_place,nnstress_short,9,&
            mpi_real8,mpi_sum,mpi_comm_world,mpierror)
        endif 
!!
        if(lfinetime)then
          call abstime(timecomm1end,daycomm1)
          timecomm1=timecomm1+timecomm1end-timecomm1start
        endif ! lfinetime
!!
      endif ! lshort
!!
!!======================================================================
!!======================================================================
!! General output part 
!! (this is not ideal, all output should be done in predict.f90)
!!======================================================================
!!======================================================================
!!
      nnshortenergy=nntotalenergy 
!!======================================================================
!! check for short range energy extrapolation
!!======================================================================
      if((mpirank.eq.0).and.lshort.and.(.not.lmd))then
        if((nntotalenergy/dble(max_num_atoms)).gt.eshortmax)then
          write(ounit,*)'-------------------------------------------------------------'
          write(ounit,'(a)')' WARNING: eshort is .gt. eshortmax '
          write(ounit,'(a,2f20.10)')' average eshort per atom, eshortmax ',&
            nntotalenergy/dble(max_num_atoms),eshortmax
          write(ounit,*)'-------------------------------------------------------------'
        endif
        if((nntotalenergy/dble(max_num_atoms)).lt.eshortmin)then
          write(ounit,*)'-------------------------------------------------------------'
          write(ounit,'(a)')' WARNING: eshort is .lt. eshortmin '
          write(ounit,'(a,2f20.10)')' average eshort per atom, eshortmin ',&
            nntotalenergy/dble(max_num_atoms),eshortmin
          write(ounit,*)'-------------------------------------------------------------'
        endif
      endif ! mpirank
!!'
      return
      end
