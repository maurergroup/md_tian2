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
!! Purpose: calculate symfunction, dsfuncdxyz and strs

!! called by:
!! - getchargesatomic.f90
!! - getshortatomic.f90 
!! - calconefunction_atomic.f90
!!
      subroutine getatomsymfunctions(i1,i2,iindex,natoms,atomindex,natomsdim,&
        max_num_atoms,max_num_neighbors_local,&
        invneighboridx_local,jcount,listdim,lsta,lstc,lste,&
        symelement_local,maxnum_funcvalues_local,&
        cutoff_type,cutoff_alpha,nelem,function_type_local,&
        lstb,funccutoff_local,xyzstruct,symfunction_temp,dsfuncdxyz_temp,strs_temp,&
        eta_local,zeta_local,lambda_local,rshift_local,rmin,&
        ldoforces,ldostress)
!!
      use fileunits
!! don't use globaloptions here
!! don't use symfunctions here
!!
      implicit none
!!
      integer listdim                                            ! in
      integer maxnum_funcvalues_local                            ! in
      integer function_type_local(maxnum_funcvalues_local,nelem) ! in
      integer symelement_local(maxnum_funcvalues_local,2,nelem)  ! in
      integer nelem                                              ! in
      integer lsta(2,max_num_atoms)                              ! in, numbers of neighbors
      integer lstc(listdim)                                      ! in, identification of atom
      integer lste(listdim)                                      ! in, nuclear charge of atom
      integer iindex                                             ! in
      integer i1,i2                                              ! in
      integer i3                                                 ! internal
      integer natoms                                             ! in
      integer atomindex(natoms)                                  ! in
      integer natomsdim                                          ! in
      integer jcount                                             ! in
      integer cutoff_type                                        ! in 
      real*8  cutoff_alpha
      integer max_num_atoms                                      ! in
      integer max_num_neighbors_local                            ! in
      integer invneighboridx_local(natoms,max_num_atoms)         ! in
!!
      real*8 xyzstruct(3,max_num_atoms)                          ! in
      real*8 symfunction_temp(maxnum_funcvalues_local)           ! out
      real*8 dsfuncdxyz_temp(0:max_num_neighbors_local,3)        ! out
      real*8 strs_temp(3,3,maxnum_funcvalues_local)              ! out
      real*8 lstb(listdim,4)                                     ! in, xyz and r_ij 
      real*8 funccutoff_local(maxnum_funcvalues_local,nelem)     ! in
      real*8 eta_local(maxnum_funcvalues_local,nelem)            ! in
      real*8 rshift_local(maxnum_funcvalues_local,nelem)         ! in
      real*8 lambda_local(maxnum_funcvalues_local,nelem)         ! in
      real*8 zeta_local(maxnum_funcvalues_local,nelem)           ! in
      real*8 rmin                                                ! in

      logical ldoforces                                          ! in
      logical ldostress                                          ! in
      logical lrmin                                              ! internal, to be removed?
!!
      dsfuncdxyz_temp(:,:) = 0.0d0
      lrmin                = .true.
!!
      if(function_type_local(i2,iindex).eq.1)then ! radial function
!!
        call atomsymfunction1(i1,i2,iindex,natoms,nelem,&
          max_num_atoms,max_num_neighbors_local,&
          invneighboridx_local,&
          jcount,listdim,lsta,lstc,lste,symelement_local,maxnum_funcvalues_local,&
          lstb,funccutoff_local,xyzstruct,symfunction_temp,dsfuncdxyz_temp,strs_temp,&
          ldoforces,ldostress)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      elseif(function_type_local(i2,iindex).eq.2)then ! radial function
!!
        call atomsymfunction2(i1,i2,iindex,natoms,atomindex,nelem,&
          max_num_atoms,max_num_neighbors_local,&
          invneighboridx_local,&
          jcount,listdim,lsta,lstc,lste,symelement_local,maxnum_funcvalues_local,&
          cutoff_type,cutoff_alpha,lstb,funccutoff_local,xyzstruct,eta_local,rshift_local,&
          symfunction_temp,dsfuncdxyz_temp,strs_temp,&
          ldoforces,ldostress)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      elseif(function_type_local(i2,iindex).eq.3)then ! angular function
!!
        if(ldostress)then
          call atomsymfunction3(i1,i2,iindex,natoms,atomindex,natomsdim,nelem,&
            max_num_atoms,max_num_neighbors_local,&
            invneighboridx_local,&
            jcount,listdim,lsta,lstc,lste,symelement_local,maxnum_funcvalues_local,&
            cutoff_type,cutoff_alpha,lstb,funccutoff_local,xyzstruct,eta_local,zeta_local,&
            lambda_local,rmin,symfunction_temp,dsfuncdxyz_temp,strs_temp,&
            ldoforces,ldostress,lrmin)
!! old
!!          call atomsymfunction3(i1,i2,iindex,natoms,atomindex,natomsdim,nelem,&
!!            max_num_atoms,max_num_neighbors_local,&
!!            invneighboridx_local,&
!!            jcount,listdim,lsta,lstc,lste,symelement_local,maxnum_funcvalues_local,&
!!            cutoff_type,lstb,funccutoff_local,xyzstruct,eta_local,zeta_local,&
!!            lambda_local,rmin,symfunction_temp,dsfuncdxyz_temp,strs_temp,&
!!            ldoforces,ldostress,lrmin)
        else
          call atomsymfunction3Andi(i1,i2,iindex,natoms,atomindex,natomsdim,nelem,&
            max_num_atoms,max_num_neighbors_local,&
            invneighboridx_local,&
            jcount,listdim,lsta,lstc,lste,symelement_local,maxnum_funcvalues_local,&
            cutoff_type,cutoff_alpha,lstb,funccutoff_local,xyzstruct,eta_local,zeta_local,&
            lambda_local,rmin,symfunction_temp,dsfuncdxyz_temp,strs_temp,&
            ldoforces,ldostress,lrmin)
        endif

!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      elseif(function_type_local(i2,iindex).eq.4)then ! radial function
!!
        call atomsymfunction4(i1,i2,iindex,natoms,atomindex,natomsdim,nelem,&
          max_num_atoms,max_num_neighbors_local,&
          invneighboridx_local,&
          jcount,listdim,lsta,lstc,lste,symelement_local,maxnum_funcvalues_local,&
          cutoff_type,cutoff_alpha,lstb,funccutoff_local,xyzstruct,eta_local,&
          symfunction_temp,dsfuncdxyz_temp,strs_temp,&
          ldoforces,ldostress)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      elseif(function_type_local(i2,iindex).eq.5)then ! just Cartesian coordinate
!!
        call atomsymfunction5(i1,i2,iindex,natoms,atomindex,natomsdim,nelem,&
          max_num_atoms,max_num_neighbors_local,&
          invneighboridx_local,&
          jcount,maxnum_funcvalues_local,&
          xyzstruct,eta_local,symfunction_temp,dsfuncdxyz_temp,strs_temp,&
          ldoforces,ldostress)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      elseif(function_type_local(i2,iindex).eq.6)then ! radial function
!!
        call atomsymfunction6(i1,i2,iindex,natoms,atomindex,natomsdim,nelem,&
          max_num_atoms,max_num_neighbors_local,&
          invneighboridx_local,&
          jcount,listdim,lsta,lstc,lste,symelement_local,maxnum_funcvalues_local,&
          lstb,funccutoff_local,xyzstruct,&
          symfunction_temp,dsfuncdxyz_temp,strs_temp,&
          ldoforces,ldostress)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      elseif(function_type_local(i2,iindex).eq.7)then ! angular function
!!
        write(ounit,*)'Error: function type not implemented in getatomsymfunction ',function_type_local(i1,iindex)
        stop !'
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      elseif(function_type_local(i2,iindex).eq.8)then ! angular function
!!
        call atomsymfunction8(i1,i2,iindex,natoms,atomindex,natomsdim,nelem,&
          max_num_atoms,max_num_neighbors_local,&
          invneighboridx_local,&
          jcount,listdim,lsta,lstc,lste,symelement_local,maxnum_funcvalues_local,&
          cutoff_type,cutoff_alpha,lstb,funccutoff_local,xyzstruct,eta_local,rshift_local,&
          lambda_local,rmin,symfunction_temp,dsfuncdxyz_temp,strs_temp,&
          ldoforces,ldostress,lrmin)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      elseif(function_type_local(i2,iindex).eq.9)then ! angular function
!!
        call atomsymfunction9(i1,i2,iindex,natoms,atomindex,natomsdim,nelem,&
          max_num_atoms,max_num_neighbors_local,&
          invneighboridx_local,&
          jcount,listdim,lsta,lstc,lste,symelement_local,maxnum_funcvalues_local,&
          cutoff_type,cutoff_alpha,lstb,funccutoff_local,xyzstruct,eta_local,zeta_local,&
          lambda_local,rmin,symfunction_temp,dsfuncdxyz_temp,strs_temp,&
          ldoforces,ldostress,lrmin)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      else
        write(*,*)'Error: function_type_local not implemented ',function_type_local(i2,iindex)
      endif
!!
!! check for reasonable interatomic distances
      if(.not.lrmin)then
        write(ounit,*)'ERROR: too short bond in getatomsymfunctions ',i1,i2
        do i3=lsta(1,atomindex(i1)),lsta(2,atomindex(i1))
          write(ounit,'(i5,3f14.8)')i3,lstb(i3,1),lstb(i3,2),lstb(i3,3)
        enddo
        stop
      endif
!!
      return
      end
