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
      subroutine atomsymfunction9(i1,i2,iindex,natoms,atomindex,natomsdim,nelem,&
        max_num_atoms,max_num_neighbors_local,invneighboridx,&
        jcount,listdim,lsta,lstc,lste,symelement_local,maxnum_funcvalues_local,&
        cutoff_type,cutoff_alpha,lstb,funccutoff_local,xyzstruct,eta_local,zeta_local,&
        lambda_local,rmin,symfunction_temp,dsfuncdxyz_temp,strs_temp,&
        ldoforces,ldostress,ldohessian_local,lrmin)
!!
      use fileunits
      use nnconstants 
!!
      implicit none
!!
      integer i1,i2,i3,i4
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
      integer m
      integer jcount 
      integer cutoff_type
      real*8  cutoff_alpha
      integer invneighboridx(natoms,max_num_atoms) 
      integer atomindex(natoms) 
!!
      real*8 strs_temp(3,3,maxnum_funcvalues_local)
      real*8 lstb(listdim,4)
      real*8 funccutoff_local(maxnum_funcvalues_local,nelem)
      real*8 lambda_local(maxnum_funcvalues_local,nelem)                ! in
      real*8 rij
      real*8 rik
      real*8 rjk
      real*8 symfunction_temp(maxnum_funcvalues_local)
      real*8 xyzstruct(3,max_num_atoms)
      real*8 deltaxj,deltayj,deltazj
      real*8 deltaxk,deltayk,deltazk
      real*8 drijdxi, drijdyi, drijdzi
      real*8 drijdxj, drijdyj, drijdzj
      real*8 drijdxk, drijdyk, drijdzk
      real*8 drikdxi, drikdyi, drikdzi
      real*8 drikdxj, drikdyj, drikdzj
      real*8 drikdxk, drikdyk, drikdzk
      real*8 drjkdxi, drjkdyi, drjkdzi
      real*8 drjkdxj, drjkdyj, drjkdzj
      real*8 drjkdxk, drjkdyk, drjkdzk
      real*8 temp1
      real*8 temp2
      real*8 temp3
      real*8 temp4
      real*8 dsfuncdxyz_temp(0:max_num_neighbors_local,3) 
      real*8 eta_local(maxnum_funcvalues_local,nelem)                            ! in
      real*8 zeta_local(maxnum_funcvalues_local,nelem)                      ! in
      real*8 fcutij
      real*8 fcutik
      real*8 dfcutijdxi,dfcutijdyi,dfcutijdzi
      real*8 dfcutijdxj,dfcutijdyj,dfcutijdzj
      real*8 dfcutijdxk,dfcutijdyk,dfcutijdzk
      real*8 dfcutikdxi,dfcutikdyi,dfcutikdzi
      real*8 dfcutikdxj,dfcutikdyj,dfcutikdzj
      real*8 dfcutikdxk,dfcutikdyk,dfcutikdzk
      real*8 dgdxi,dgdyi,dgdzi
      real*8 dgdxj,dgdyj,dgdzj
      real*8 dgdxk,dgdyk,dgdzk
      real*8 dfdxi,dfdyi,dfdzi
      real*8 dfdxj,dfdyj,dfdzj
      real*8 dfdxk,dfdyk,dfdzk
      real*8 rmin
      real*8 f 
      real*8 g 
      real*8 costheta
      real*8 dcosthetadxi,dcosthetadyi,dcosthetadzi
      real*8 dcosthetadxj,dcosthetadyj,dcosthetadzj
      real*8 dcosthetadxk,dcosthetadyk,dcosthetadzk
      real*8 expxyz
      real*8 dexpxyzdxi
      real*8 dexpxyzdyi
      real*8 dexpxyzdzi
      real*8 dexpxyzdxj
      real*8 dexpxyzdyj
      real*8 dexpxyzdzj
      real*8 dexpxyzdxk
      real*8 dexpxyzdyk
      real*8 dexpxyzdzk
!!
      logical ldoforces
      logical ldostress
      logical ldohessian_local
      logical lrmin
!!
        do i3=lsta(1,atomindex(i1)),lsta(2,atomindex(i1)) ! loop over neighbors
          if(symelement_local(i2,1,iindex).eq.lste(i3))then
          rij=lstb(i3,4)
          n=lstc(i3)
          if(rij.le.funccutoff_local(i2,iindex))then
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
          drijdxk=0.0d0
          drijdyk=0.0d0
          drijdzk=0.0d0
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
          dfcutijdxk=0.0d0
          dfcutijdyk=0.0d0
          dfcutijdzk=0.0d0
!!
!!
          do i4=lsta(1,atomindex(i1)),lsta(2,atomindex(i1)) ! loop over all neighbors of this atom
            if(symelement_local(i2,2,iindex).eq.lste(i4))then
!! new double counting criterion JB and TM 08.03.2010:
            if(((i4.gt.i3).and.(lste(i3).eq.lste(i4)))&
             .or.(lste(i3).ne.lste(i4)))then ! avoid double counting
!! This was a bug, because the order of neighbor atoms is not arbitrary
!!            if(i4.gt.i3)then ! avoid double counting
!!
            m=lstc(i4)
            rik=lstb(i4,4)
!!
            if(rik.le.funccutoff_local(i2,iindex))then
!!
              deltaxk=-1.d0*(xyzstruct(1,jcount)-lstb(i4,1))
              deltayk=-1.d0*(xyzstruct(2,jcount)-lstb(i4,2))
              deltazk=-1.d0*(xyzstruct(3,jcount)-lstb(i4,3))
              drikdxi=-deltaxk/rik
              drikdyi=-deltayk/rik
              drikdzi=-deltazk/rik
              drikdxk=-1.d0*drikdxi
              drikdyk=-1.d0*drikdyi
              drikdzk=-1.d0*drikdzi
              drikdxj=0.0d0
              drikdyj=0.0d0
              drikdzj=0.0d0
!!
              call getcutoff(&
                cutoff_type,cutoff_alpha,maxnum_funcvalues_local,nelem,&
                i2,iindex,&
                funccutoff_local,rik,fcutik,temp1)

              dfcutikdxi=temp1*drikdxi
              dfcutikdyi=temp1*drikdyi
              dfcutikdzi=temp1*drikdzi
              dfcutikdxj=0.0d0
              dfcutikdyj=0.0d0
              dfcutikdzj=0.0d0
              dfcutikdxk=-1.d0*dfcutikdxi
              dfcutikdyk=-1.d0*dfcutikdyi
              dfcutikdzk=-1.d0*dfcutikdzi
!!
!! calculate the third bond length between atom i3 and i4
!!
              rjk=(lstb(i3,1)-lstb(i4,1))**2 +&
                  (lstb(i3,2)-lstb(i4,2))**2 +&
                  (lstb(i3,3)-lstb(i4,3))**2
              rjk=dsqrt(rjk)
              if(rjk.le.rmin) then
                lrmin=.false.
              endif
!!
                drjkdxj=(lstb(i3,1)-lstb(i4,1))/rjk
                drjkdyj=(lstb(i3,2)-lstb(i4,2))/rjk
                drjkdzj=(lstb(i3,3)-lstb(i4,3))/rjk
                drjkdxk=-1.d0*drjkdxj
                drjkdyk=-1.d0*drjkdyj
                drjkdzk=-1.d0*drjkdzj
                drjkdxi=0.0d0
                drjkdyi=0.0d0
                drjkdzi=0.0d0
!!
!! costheta=(rjk**2 - rij**2 -rik**2)/(-2.d0*rij*rik)
                f=lambda_local(i2,iindex)*(rjk**2 - rij**2 -rik**2)
                g=-2.d0*rij*rik
                costheta=f/g
                costheta=1.d0+costheta ! avoid negative values
                temp4   =2.d0**(1.d0-zeta_local(i2,iindex))
                costheta=temp4*(costheta**zeta_local(i2,iindex))
!!
!! calculate the derivatives of costheta already here
!! (f/g)' = (f'g - fg')/g^2
                if(ldoforces) then
                 dfdxi=lambda_local(i2,iindex)*(-2.d0*rij*drijdxi - 2.d0*rik*drikdxi)
                 dfdyi=lambda_local(i2,iindex)*(-2.d0*rij*drijdyi - 2.d0*rik*drikdyi)
                 dfdzi=lambda_local(i2,iindex)*(-2.d0*rij*drijdzi - 2.d0*rik*drikdzi)
!!
                 dfdxj=lambda_local(i2,iindex)*(2.d0*rjk*drjkdxj - 2.d0*rij*drijdxj)
                 dfdyj=lambda_local(i2,iindex)*(2.d0*rjk*drjkdyj - 2.d0*rij*drijdyj)
                 dfdzj=lambda_local(i2,iindex)*(2.d0*rjk*drjkdzj - 2.d0*rij*drijdzj)
!!
                 dfdxk=lambda_local(i2,iindex)*(2.d0*rjk*drjkdxk - 2.d0*rik*drikdxk)
                 dfdyk=lambda_local(i2,iindex)*(2.d0*rjk*drjkdyk - 2.d0*rik*drikdyk)
                 dfdzk=lambda_local(i2,iindex)*(2.d0*rjk*drjkdzk - 2.d0*rik*drikdzk)
!!
                 dgdxi=-2.d0*(drijdxi*rik + rij*drikdxi)
                 dgdyi=-2.d0*(drijdyi*rik + rij*drikdyi)
                 dgdzi=-2.d0*(drijdzi*rik + rij*drikdzi)
!!
                 dgdxj=-2.d0*drijdxj*rik
                 dgdyj=-2.d0*drijdyj*rik
                 dgdzj=-2.d0*drijdzj*rik
!!
                 dgdxk=-2.d0*rij*drikdxk
                 dgdyk=-2.d0*rij*drikdyk
                 dgdzk=-2.d0*rij*drikdzk
!!
                 temp4=temp4*zeta_local(i2,iindex)*(1.d0+f/g)**(zeta_local(i2,iindex)-1)/g**2
                 dcosthetadxi=temp4*(dfdxi*g - f*dgdxi) !/g**2
                 dcosthetadyi=temp4*(dfdyi*g - f*dgdyi) !/g**2
                 dcosthetadzi=temp4*(dfdzi*g - f*dgdzi) !/g**2
                 dcosthetadxj=temp4*(dfdxj*g - f*dgdxj) !/g**2
                 dcosthetadyj=temp4*(dfdyj*g - f*dgdyj) !/g**2
                 dcosthetadzj=temp4*(dfdzj*g - f*dgdzj) !/g**2
                 dcosthetadxk=temp4*(dfdxk*g - f*dgdxk) !/g**2
                 dcosthetadyk=temp4*(dfdyk*g - f*dgdyk) !/g**2
                 dcosthetadzk=temp4*(dfdzk*g - f*dgdzk) !/g**2
                endif ! ldoforces for derivative of costheta
!!
!! calculation of exponential term
                expxyz=dexp(-eta_local(i2,iindex)*(rij**2+rik**2))
                temp1=-eta_local(i2,iindex)*2.0d0*expxyz
                dexpxyzdxi=(rij*drijdxi+rik*drikdxi)*temp1
                dexpxyzdyi=(rij*drijdyi+rik*drikdyi)*temp1
                dexpxyzdzi=(rij*drijdzi+rik*drikdzi)*temp1
                dexpxyzdxj=(rij*drijdxj+rik*drikdxj)*temp1
                dexpxyzdyj=(rij*drijdyj+rik*drikdyj)*temp1
                dexpxyzdzj=(rij*drijdzj+rik*drikdzj)*temp1
                dexpxyzdxk=(rij*drijdxk+rik*drikdxk)*temp1
                dexpxyzdyk=(rij*drijdyk+rik*drikdyk)*temp1
                dexpxyzdzk=(rij*drijdzk+rik*drikdzk)*temp1
!!
!! calculation of the symmetry function
          symfunction_temp(i2)=symfunction_temp(i2)&
                                 +costheta*expxyz*fcutij*fcutik
!!
              if(ldoforces)then
!! Calculation of derivatives for forces 
!! (fgh)' = f'gh + fg'h + fgh'
!! (fghi)'= f'ghi + fg'hi + fgh'i + fghi'
!!
!! dsfunc/dx
               temp1=(+dcosthetadxi* expxyz   * fcutij   * fcutik   &
                       +costheta   *dexpxyzdxi* fcutij   * fcutik   &
                       +costheta   * expxyz   *dfcutijdxi* fcutik   &
                       +costheta   * expxyz   * fcutij   *dfcutikdxi)
               temp2=(+dcosthetadxj* expxyz   * fcutij   * fcutik   &
                       +costheta   *dexpxyzdxj* fcutij   * fcutik   &
                       +costheta   * expxyz   *dfcutijdxj* fcutik   &
                       +costheta   * expxyz   * fcutij   *dfcutikdxj)
               temp3=(+dcosthetadxk* expxyz   * fcutij   * fcutik   &
                       +costheta   *dexpxyzdxk* fcutij   * fcutik   &
                       +costheta   * expxyz   *dfcutijdxk* fcutik   &
                       +costheta   * expxyz   * fcutij   *dfcutikdxk)
            dsfuncdxyz_temp(invneighboridx(i1,jcount),1)=dsfuncdxyz_temp(invneighboridx(i1,jcount),1)+&
                                       temp1
            dsfuncdxyz_temp(invneighboridx(i1,n),1) =dsfuncdxyz_temp(invneighboridx(i1,n),1)+&
                                       temp2
            dsfuncdxyz_temp(invneighboridx(i1,m),1) =dsfuncdxyz_temp(invneighboridx(i1,m),1)+&
                                       temp3
              if(ldostress)then
!! Calculation of derivatives for stress
                strs_temp(1,1,i2)=strs_temp(1,1,i2)+deltaxj*temp2&
                  +deltaxk*temp3
                strs_temp(2,1,i2)=strs_temp(2,1,i2)+deltayj*temp2&
                  +deltayk*temp3
                strs_temp(3,1,i2)=strs_temp(3,1,i2)+deltazj*temp2&
                  +deltazk*temp3
              endif ! ldostress 
!!
!! dsfunc/dy
               temp1=(+dcosthetadyi* expxyz   * fcutij   * fcutik   &
                       +costheta   *dexpxyzdyi* fcutij   * fcutik   &
                       +costheta   * expxyz   *dfcutijdyi* fcutik   &
                       +costheta   * expxyz   * fcutij   *dfcutikdyi)
               temp2=(+dcosthetadyj* expxyz   * fcutij   * fcutik   &
                       +costheta   *dexpxyzdyj* fcutij   * fcutik   &
                       +costheta   * expxyz   *dfcutijdyj* fcutik   &
                       +costheta   * expxyz   * fcutij   *dfcutikdyj)
               temp3=(+dcosthetadyk* expxyz   * fcutij   * fcutik   &
                       +costheta   *dexpxyzdyk* fcutij   * fcutik   &
                       +costheta   * expxyz   *dfcutijdyk* fcutik   &
                       +costheta   * expxyz   * fcutij   *dfcutikdyk)
            dsfuncdxyz_temp(invneighboridx(i1,jcount),2)=dsfuncdxyz_temp(invneighboridx(i1,jcount),2)+&
                                       temp1
            dsfuncdxyz_temp(invneighboridx(i1,n),2) =dsfuncdxyz_temp(invneighboridx(i1,n),2)+&
                                       temp2
            dsfuncdxyz_temp(invneighboridx(i1,m),2) =dsfuncdxyz_temp(invneighboridx(i1,m),2)+&
                                       temp3
              if(ldostress)then
!! Calculation of derivatives for stress
                strs_temp(1,2,i2)=strs_temp(1,2,i2)+deltaxj*temp2&
                  +deltaxk*temp3
                strs_temp(2,2,i2)=strs_temp(2,2,i2)+deltayj*temp2&
                  +deltayk*temp3
                strs_temp(3,2,i2)=strs_temp(3,2,i2)+deltazj*temp2&
                  +deltazk*temp3
              endif ! ldostress 
!!
!! dsfunc/dz
               temp1=(+dcosthetadzi* expxyz   * fcutij   * fcutik   &
                      + costheta   *dexpxyzdzi* fcutij   * fcutik   &
                      + costheta   * expxyz   *dfcutijdzi* fcutik   &
                      + costheta   * expxyz   * fcutij   *dfcutikdzi)
               temp2=(+dcosthetadzj* expxyz   * fcutij   * fcutik   &
                      + costheta   *dexpxyzdzj* fcutij   * fcutik   &
                      + costheta   * expxyz   *dfcutijdzj* fcutik   &
                      + costheta   * expxyz   * fcutij   *dfcutikdzj)
               temp3=(+dcosthetadzk* expxyz   * fcutij   * fcutik   &
                      + costheta   *dexpxyzdzk* fcutij   * fcutik   &
                      + costheta   * expxyz   *dfcutijdzk* fcutik   &
                      + costheta   * expxyz   * fcutij   *dfcutikdzk)
            dsfuncdxyz_temp(invneighboridx(i1,jcount),3)=dsfuncdxyz_temp(invneighboridx(i1,jcount),3)+&
                                       temp1
            dsfuncdxyz_temp(invneighboridx(i1,n),3) =dsfuncdxyz_temp(invneighboridx(i1,n),3)+&
                                       temp2
            dsfuncdxyz_temp(invneighboridx(i1,m),3) =dsfuncdxyz_temp(invneighboridx(i1,m),3)+&
                                       temp3
              if(ldostress)then
!! Calculation of derivatives for stress
                strs_temp(1,3,i2)=strs_temp(1,3,i2)+deltaxj*temp2&
                  +deltaxk*temp3
                strs_temp(2,3,i2)=strs_temp(2,3,i2)+deltayj*temp2&
                  +deltayk*temp3
                strs_temp(3,3,i2)=strs_temp(3,3,i2)+deltazj*temp2&
                  +deltazk*temp3
              endif ! ldostress 
              endif ! ldoforces
!!
            endif ! rik .le. cutoff
!!
          endif ! i4.gt.i3
          endif ! (symelement_local(i2,2,iindex).eq.lste(i4))then
          enddo ! i4
          endif ! rij.le.funccutoff_local(i2,iindex)
          endif ! (symelement_local(i2,1,iindex).eq.lste(i3))then
        enddo ! i3
!!
      return
      end
