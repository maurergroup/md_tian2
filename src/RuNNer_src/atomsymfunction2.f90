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
      subroutine atomsymfunction2(i1,i2,iindex,natoms,atomindex,nelem,&
        max_num_atoms,max_num_neighbors_local,invneighboridx,&
        jcount,listdim,lsta,lstc,lste,symelement_local,maxnum_funcvalues_local,&
        cutoff_type,cutoff_alpha,lstb,funccutoff_local,xyzstruct,eta_local,rshift_local,&
        symfunction_temp,dsfuncdxyz_temp,strs_temp,&
        ldoforces,ldostress,ldohessian_local)
!!
      use fileunits
      use nnconstants
!!
      implicit none
!!
      integer i1                                                         ! in, relative number of current central atom
      integer i2                                                         ! in, current symmetry function number for this atom
      integer i3                                                         ! internal
      integer listdim
      integer nelem 
      integer natoms 
      integer atomindex(natoms) 
      integer iindex
      integer max_num_atoms 
      integer max_num_neighbors_local 
      integer maxnum_funcvalues_local 
      integer symelement_local(maxnum_funcvalues_local,2,nelem)          ! in
      integer lsta(2,max_num_atoms)
      integer lstc(listdim)
      integer lste(listdim)
      integer n
      integer jcount                                                     ! in, absolute number of current central atom 
      integer cutoff_type
      real*8  cutoff_alpha
      integer invneighboridx(natoms,max_num_atoms) 
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
      real*8 temp1
      real*8 temp2
      real*8 dsfuncdxyz_temp(0:max_num_neighbors_local,3) 
      real*8 eta_local(maxnum_funcvalues_local,nelem)                            ! in
      real*8 rshift_local(maxnum_funcvalues_local,nelem)                ! in
      real*8 fcutij
      real*8 dfcutijdxi,dfcutijdyi,dfcutijdzi
      real*8 dfcutijdxj,dfcutijdyj,dfcutijdzj
      real*8 dfcutijdxk,dfcutijdyk,dfcutijdzk
!!
      logical ldoforces
      logical ldostress
      logical ldohessian_local
!!
        do i3=lsta(1,atomindex(i1)),lsta(2,atomindex(i1)) ! loop over neighbors
          if(symelement_local(i2,1,iindex).eq.lste(i3))then
!!
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
              call getcutoff(&
                cutoff_type,cutoff_alpha,maxnum_funcvalues_local,nelem,&
                i2,iindex,&
                funccutoff_local,rij,fcutij,temp1)

              dfcutijdxi=temp1*drijdxi
              dfcutijdyi=temp1*drijdyi
              dfcutijdzi=temp1*drijdzi
              dfcutijdxj=-1.d0*dfcutijdxi
              dfcutijdyj=-1.d0*dfcutijdyi
              dfcutijdzj=-1.d0*dfcutijdzi
!!
          symfunction_temp(i2)=symfunction_temp(i2)&
            +dexp(-1.d0*eta_local(i2,iindex)*(rij-rshift_local(i2,iindex))**2)*fcutij 
!!
          if(ldoforces)then
!! Calculation of derivatives for forces 
            temp1=-2.d0*eta_local(i2,iindex)*(rij-rshift_local(i2,iindex))*dexp(-1.d0*eta_local(i2,iindex)&
                  *(rij-rshift_local(i2,iindex))**2)*fcutij
            temp2= dexp(-1.d0*eta_local(i2,iindex)*(rij-rshift_local(i2,iindex))**2)
!! dsfunc/dx
            dsfuncdxyz_temp(invneighboridx(i1,jcount),1)=dsfuncdxyz_temp(invneighboridx(i1,jcount),1)+&
               (drijdxi*temp1 + temp2*dfcutijdxi)
            dsfuncdxyz_temp(invneighboridx(i1,n),1) =dsfuncdxyz_temp(invneighboridx(i1,n),1)+&
               (drijdxj*temp1 + temp2*dfcutijdxj)
!! dsfunc/dy
            dsfuncdxyz_temp(invneighboridx(i1,jcount),2)=dsfuncdxyz_temp(invneighboridx(i1,jcount),2)+&
               (drijdyi*temp1 + temp2*dfcutijdyi)
            dsfuncdxyz_temp(invneighboridx(i1,n),2) =dsfuncdxyz_temp(invneighboridx(i1,n),2)+&
               (drijdyj*temp1 + temp2*dfcutijdyj)
!! dsfunc/dz
            dsfuncdxyz_temp(invneighboridx(i1,jcount),3)=dsfuncdxyz_temp(invneighboridx(i1,jcount),3)+&
               (drijdzi*temp1 + temp2*dfcutijdzi)
            dsfuncdxyz_temp(invneighboridx(i1,n),3) =dsfuncdxyz_temp(invneighboridx(i1,n),3)+&
               (drijdzj*temp1 + temp2*dfcutijdzj)
!!
          if(ldostress)then
!! Calculation of derivatives for stress
!! dsfunc/dx
            strs_temp(1,1,i2)=strs_temp(1,1,i2)+deltaxj*&
               (drijdxj*temp1&
              + temp2*dfcutijdxj)
            strs_temp(2,1,i2)=strs_temp(2,1,i2)+deltayj*&
               (drijdxj*temp1&
              + temp2*dfcutijdxj)
            strs_temp(3,1,i2)=strs_temp(3,1,i2)+deltazj*&
               (drijdxj*temp1&
              + temp2*dfcutijdxj)
!! dsfunc/dy
            strs_temp(1,2,i2)=strs_temp(1,2,i2)+deltaxj*&
               (drijdyj*temp1&
              + temp2*dfcutijdyj)
            strs_temp(2,2,i2)=strs_temp(2,2,i2)+deltayj*&
               (drijdyj*temp1&
              + temp2*dfcutijdyj)
            strs_temp(3,2,i2)=strs_temp(3,2,i2)+deltazj*&
               (drijdyj*temp1&
              + temp2*dfcutijdyj)
!! dsfunc/dz
            strs_temp(1,3,i2)=strs_temp(1,3,i2)+deltaxj*&
               (drijdzj*temp1&
              + temp2*dfcutijdzj)
            strs_temp(2,3,i2)=strs_temp(2,3,i2)+deltayj*&
               (drijdzj*temp1&
              + temp2*dfcutijdzj)
            strs_temp(3,3,i2)=strs_temp(3,3,i2)+deltazj*&
               (drijdzj*temp1&
              + temp2*dfcutijdzj)
!!
          endif ! ldostress 
          endif ! ldoforces
          endif ! rij.le.funccutoff_local(i2,iindex)
          endif ! (symelement_local(i2,1,iindex).eq.lste(i3))then
        enddo ! i3

      return
      end
