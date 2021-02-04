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

!! called by:
!! - getatomsymfunctions.f90
!!
      subroutine atomsymfunction5(i1,i2,iindex,natoms,atomindex,natomsdim,nelem,&
        max_num_atoms,max_num_neighbors_local,invneighboridx,&
        jcount,maxnum_funcvalues_local,&
        xyzstruct,eta_local,symfunction_temp,dsfuncdxyz_temp,strs_temp,&
        ldoforces,ldostress)
!!
      use fileunits
!!
      implicit none
!!
      integer i1,i2
      integer nelem 
      integer natomsdim 
      integer natoms 
      integer iindex
      integer max_num_atoms 
      integer max_num_neighbors_local 
      integer maxnum_funcvalues_local 
      integer jcount
      integer invneighboridx(natoms,max_num_atoms) 
      integer atomindex(natoms)
!!
      real*8 strs_temp(3,3,maxnum_funcvalues_local)
      real*8 symfunction_temp(maxnum_funcvalues_local)
      real*8 xyzstruct(3,max_num_atoms)
      real*8 dsfuncdxyz_temp(0:max_num_neighbors_local,3) 
      real*8 eta_local(maxnum_funcvalues_local,nelem)                            ! in
!!
      logical ldoforces
      logical ldostress
!!

      if(eta_local(i2,iindex).eq.1.0d0)then
        symfunction_temp(i2)=xyzstruct(1,i1)
      elseif(eta_local(i2,iindex).eq.2.0d0)then
        symfunction_temp(i2)=xyzstruct(2,i1)
      elseif(eta_local(i2,iindex).eq.3.0d0)then
        symfunction_temp(i2)=xyzstruct(3,i1)
      else
        write(ounit,*)'Error: undefined value for symfunction type 5 ',eta_local(i2,iindex)
        stop !'
      endif

      if(ldoforces)then
!! Calculation of derivatives for forces
!! dsfunc/dx
        if(eta_local(i2,iindex).eq.1.0d0)then
          dsfuncdxyz_temp(invneighboridx(i1,jcount),1)=1.0d0
        else
          dsfuncdxyz_temp(invneighboridx(i1,jcount),1)=0.0d0
        endif
!! dsfunc/dy
        if(eta_local(i2,iindex).eq.2.0d0)then
          dsfuncdxyz_temp(invneighboridx(i1,jcount),2)=1.0d0
        else
          dsfuncdxyz_temp(invneighboridx(i1,jcount),2)=0.0d0
        endif
!! dsfunc/dz
        if(eta_local(i2,iindex).eq.3.0d0)then
          dsfuncdxyz_temp(invneighboridx(i1,jcount),3)=1.0d0
        else
          dsfuncdxyz_temp(invneighboridx(i1,jcount),3)=0.0d0
        endif
!!
        if(ldostress)then
!! Calculation of derivatives for stress
          write(ounit,*)'Error: stress is not implemented for symfunction type 5'
          stop !'
!!
        endif ! ldostress
      endif ! ldoforces

      return
      end
