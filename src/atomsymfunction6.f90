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
!! - getatomsymfunctions.f90
!!
      subroutine atomsymfunction6(i1,i2,iindex,natoms,atomindex,natomsdim,nelem,&
        max_num_atoms,max_num_neighbors_local,invneighboridx,&
        jcount,listdim,lsta,lstc,lste,symelement_local,maxnum_funcvalues_local,&
        lstb,funccutoff_local,xyzstruct,&
        symfunction_temp,dsfuncdxyz_temp,strs_temp,&
        ldoforces,ldostress)
!!
      use fileunits
!!
      implicit none
!!
      integer i1,i2,i3
      integer listdim
      integer nelem 
      integer natomsdim 
      integer natoms 
      integer iindex
      integer max_num_atoms 
      integer max_num_neighbors_local 
      integer maxnum_funcvalues_local 
      integer symelement_local(maxnum_funcvalues_local,2,nelem)          ! in
      integer lsta(2,max_num_atoms)
      integer lstc(listdim)
      integer lste(listdim)
      integer n
      integer jcount
      integer invneighboridx(natoms,max_num_atoms) 
      integer atomindex(natoms) 
!!
      real*8 strs_temp(3,3,maxnum_funcvalues_local)
      real*8 lstb(listdim,4)
      real*8 funccutoff_local(maxnum_funcvalues_local,nelem)
      real*8 rij
      real*8 symfunction_temp(maxnum_funcvalues_local)
      real*8 xyzstruct(3,max_num_atoms)
      real*8 deltaxj,deltayj,deltazj
      real*8 drijdxi, drijdyi, drijdzi
      real*8 drijdxj, drijdyj, drijdzj
      real*8 dsfuncdxyz_temp(0:max_num_neighbors_local,3) 
!!
      logical ldoforces
      logical ldostress
!!

        do i3=lsta(1,atomindex(i1)),lsta(2,atomindex(i1)) ! loop over neighbors
          if(symelement_local(i2,1,iindex).eq.lste(i3))then
          if(lsta(1,atomindex(i1)).ne.lsta(2,atomindex(i1)))then
            write(ounit,*)'ERROR: symmetry function type 6 is only for dimers'
            stop !'
          endif
          rij=lstb(i3,4)
          if(rij.le.funccutoff_local(i2,iindex))then
            n=lstc(i3)
!!
            deltaxj=-1.d0*(xyzstruct(1,jcount)-lstb(i3,1))
            deltayj=-1.d0*(xyzstruct(2,jcount)-lstb(i3,2))
            deltazj=-1.d0*(xyzstruct(3,jcount)-lstb(i3,3))
            drijdxi=-deltaxj/rij
            drijdyi=-deltayj/rij
            drijdzi=-deltazj/rij
            drijdxj=-1.d0*drijdxi
            drijdyj=-1.d0*drijdyi
            drijdzj=-1.d0*drijdzi
!!
            symfunction_temp(i2)=rij
!!
            if(ldoforces)then
!! Calculation of derivatives for forces
!! dsfunc/dx
              dsfuncdxyz_temp(invneighboridx(i1,jcount),1)=dsfuncdxyz_temp(invneighboridx(i1,jcount),1)+&
                 drijdxi
              dsfuncdxyz_temp(invneighboridx(i1,n),1) =dsfuncdxyz_temp(invneighboridx(i1,n),1)+&
                 drijdxj
!! dsfunc/dy
              dsfuncdxyz_temp(invneighboridx(i1,jcount),2)=dsfuncdxyz_temp(invneighboridx(i1,jcount),2)+&
                 drijdyi
              dsfuncdxyz_temp(invneighboridx(i1,n),2) =dsfuncdxyz_temp(invneighboridx(i1,n),2)+&
                 drijdyj
!! dsfunc/dz
              dsfuncdxyz_temp(invneighboridx(i1,jcount),3)=dsfuncdxyz_temp(invneighboridx(i1,jcount),3)+&
                 drijdzi
              dsfuncdxyz_temp(invneighboridx(i1,n),3) =dsfuncdxyz_temp(invneighboridx(i1,n),3)+&
                 drijdzj
!!
              if(ldostress)then
!! 
                write(ounit,*)'ERROR: no stress implemented for symmetry function type 6'
                stop !'
!!
              endif ! ldostress
            endif ! ldoforces
          endif ! rij.le.funccutoff_local(i2,iindex)
          endif !(symelement_local(i2,1,iindex).eq.lste(i3))then
        enddo ! i3

      return
      end
