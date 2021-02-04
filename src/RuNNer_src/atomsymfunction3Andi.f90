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
      subroutine atomsymfunction3Andi(i1,i2,iindex,natoms,atomindex,natomsdim,nelem,&
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
      integer m
      integer jcount 
      integer cutoff_type
      integer invneighboridx(natoms,max_num_atoms) 
!!
      real*8 strs_temp(3,3,maxnum_funcvalues_local)
      real*8 lstb(listdim,4)
      real*8 funccutoff_local(maxnum_funcvalues_local,nelem)
      real*8 cutoff_alpha                                               ! in
      real*8 lambda_local(maxnum_funcvalues_local,nelem)                ! in
      real*8 rij
      real*8 rik
      real*8 rjk
      real*8 symfunction_temp(maxnum_funcvalues_local)
      real*8 xyzstruct(3,max_num_atoms)
      real*8 dsfuncdxyz_temp(0:max_num_neighbors_local,3) 
      real*8 eta_local(maxnum_funcvalues_local,nelem)                  ! in
      real*8 zeta_local(maxnum_funcvalues_local,nelem)                 ! in
      real*8 rmin
      real*8 costheta
!!!!! EDITED by Andi !!!!!! 
      real*8 rc,rc2
      real*8 rinvijik
      real*8 r2ij,r2ik
      real*8 dxij,dyij,dzij
      real*8 dxik,dyik,dzik
      real*8 dxjk,dyjk,dzjk
      real*8 pexp,plambda_local,pfc,pnorm,pzl,p2etaplambda
      real*8 pfcij,pfcik,pfcjk
      real*8 pdfcij,pdfcik,pdfcjk
      real*8 p1,p2,p3
      real*8 fg
!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
      logical ldoforces
      logical ldostress
      logical ldohessian_local
      logical lrmin                                                     ! in/out
!!
      rc=funccutoff_local(i2,iindex)
      rc2=rc*rc
      if(ldostress) then
          write(*,*) "ERROR: stress calculation not implemented in symmetry function type 3"
          stop
      end if
      pnorm = 2.d0**(1.d0-zeta_local(i2,iindex))
      pzl = zeta_local(i2,iindex)*lambda_local(i2,iindex)

!! this loop is correct (compare with second loop over neighbors below)
      do i3=lsta(1,atomindex(i1)),lsta(2,atomindex(i1))-1
          rij=lstb(i3,4)
          if( ( rij .LE. rc ) .AND. &
              ( (symelement_local(i2,1,iindex).EQ.lste(i3)) .OR. &
                (symelement_local(i2,2,iindex).EQ.lste(i3)) ) ) then
              r2ij=rij*rij

!! calculation of cutoff function
!! original code by Andi:              
!!              if(cutoff_type.eq.1)then
!!                  pfcij =0.5d0*(dcos(pi*rij/rc)+1.d0)
!!                  pdfcij=0.5d0*(-dsin(pi*rij/rc))*(pi/rc)
!!              elseif(cutoff_type.eq.2)then
!!                  pfcij =(tanh(1.d0-rij/rc))**3
!!                  pdfcij=(-3.d0/rc)* &
!!                     ((tanh(1.d0-rij/rc))**2 &
!!                    - (tanh(1.d0-rij/rc))**4)
!!              else
!!                  write(ounit,*)'ERROR: Unknown cutoff_type in calconefunction_para'
!!                  stop  !'
!!              endif
!! we want to use this:
              call getcutoff(&
                cutoff_type,cutoff_alpha,maxnum_funcvalues_local,nelem,&
                i2,iindex,&
                funccutoff_local,rij,pfcij,pdfcij)



!! we only loop over unique pairs of neighbors              
              do i4=i3+1,lsta(2,atomindex(i1))
                  if( ( (symelement_local(i2,1,iindex).EQ.lste(i3)) .AND. &
                        (symelement_local(i2,2,iindex).EQ.lste(i4)) ) .OR. &
                      ( (symelement_local(i2,2,iindex).EQ.lste(i3)) .AND. &
                        (symelement_local(i2,1,iindex).EQ.lste(i4)) ) ) then
                      rik=lstb(i4,4)
                      if( rik .LE. rc ) then
                          dxjk=lstb(i3,1)-lstb(i4,1)
                          dyjk=lstb(i3,2)-lstb(i4,2)
                          dzjk=lstb(i3,3)-lstb(i4,3)
                          rjk=dxjk**2 + dyjk**2 + dzjk**2
              if(dsqrt(rjk).le.rmin) then
                write(ounit,*)'atomsymfunction3Andi atoms too close: rjk=',rjk,rmin
                lrmin=.false.
              endif
                          if( rjk .LE. rc2 ) then
                              
                              ! energy calculation
                              rjk=dsqrt(rjk)
                              rinvijik=1.d0/rij/rik
                              r2ik=rik*rik
                              dxij=xyzstruct(1,jcount)-lstb(i3,1)
                              dyij=xyzstruct(2,jcount)-lstb(i3,2)
                              dzij=xyzstruct(3,jcount)-lstb(i3,3)
                              dxik=xyzstruct(1,jcount)-lstb(i4,1)
                              dyik=xyzstruct(2,jcount)-lstb(i4,2)
                              dzik=xyzstruct(3,jcount)-lstb(i4,3)
!! calculation of cutoff function
!! original code by Andi:              
!!                              if(cutoff_type.eq.1)then
!!                                  pfcik =0.5d0*(dcos(pi*rik/rc)+1.d0)
!!                                  pfcjk =0.5d0*(dcos(pi*rjk/rc)+1.d0)
!!                                  pdfcik=0.5d0*(-dsin(pi*rik/rc))*(pi/rc)
!!                                  pdfcjk=0.5d0*(-dsin(pi*rjk/rc))*(pi/rc)
!!                              elseif(cutoff_type.eq.2)then
!!                                  pfcik =(tanh(1.d0-rik/rc))**3
!!                                  pfcjk =(tanh(1.d0-rjk/rc))**3
!!                                  pdfcik=(-3.d0/rc)* &
!!                                     ((tanh(1.d0-rik/rc))**2 &
!!                                    - (tanh(1.d0-rik/rc))**4)
!!                                  pdfcjk=(-3.d0/rc)* &
!!                                     ((tanh(1.d0-rjk/rc))**2 &
!!                                    - (tanh(1.d0-rjk/rc))**4)
!!                              else
!!                                  write(ounit,*)'ERROR: Unknown cutoff_type in calconefunction_para'
!!                                  stop  !'
!!                              endif
!! we want to use this:
              call getcutoff(&
                cutoff_type,cutoff_alpha,maxnum_funcvalues_local,nelem,&
                i2,iindex,&
                funccutoff_local,rik,pfcik,pdfcik)
              call getcutoff(&
                cutoff_type,cutoff_alpha,maxnum_funcvalues_local,nelem,&
                i2,iindex,&
                funccutoff_local,rjk,pfcjk,pdfcjk)

                              costheta = dxij*dxik + dyij*dyik + dzij*dzik
                              costheta = costheta * rinvijik
                              pfc = pfcij * pfcik * pfcjk
                              pexp = exp(-eta_local(i2,iindex) * (r2ij+r2ik+rjk*rjk))
                              plambda_local = 1.d0 + lambda_local(i2,iindex)*costheta
!!
!! CHANGE ANDI: prevent complex "fg" for non-integer values of "zeta_local"
!!
                             !fg = plambda_local**(zeta_local(i2,iindex)-1.d0) * pexp 
                              if(plambda_local .LE. 0.d0) then
                                  fg = 0.d0
                              else
                                  fg = plambda_local**(zeta_local(i2,iindex)-1.d0) * pexp 
                              endif
!! END CHANGE

                              symfunction_temp(i2) = symfunction_temp(i2) + fg*plambda_local*pfc

                              ! force calculation
                              if(ldoforces) then
                                  fg=fg*pnorm
                                  rinvijik=rinvijik*pzl
                                  costheta=costheta*pzl
                                  p2etaplambda=2.d0*eta_local(i2,iindex)*plambda_local
                                  p1 = pfc*(rinvijik - costheta/r2ij - p2etaplambda) + pfcik*pfcjk*pdfcij*plambda_local/rij
                                  p2 = pfc*(rinvijik - costheta/r2ik - p2etaplambda) + pfcij*pfcjk*pdfcik*plambda_local/rik
                                  p3 = pfc*(rinvijik                 + p2etaplambda) - pfcij*pfcik*pdfcjk*plambda_local/rjk
                                  p1 = p1*fg
                                  p2 = p2*fg
                                  p3 = p3*fg
                                  dxij = dxij*p1
                                  dyij = dyij*p1
                                  dzij = dzij*p1
                                  dxik = dxik*p2
                                  dyik = dyik*p2
                                  dzik = dzik*p2
                                  dxjk = dxjk*p3
                                  dyjk = dyjk*p3
                                  dzjk = dzjk*p3

                                  dsfuncdxyz_temp(invneighboridx(i1,jcount),1) &
                                    = dsfuncdxyz_temp(invneighboridx(i1,jcount),1) + dxij + dxik
                                  dsfuncdxyz_temp(invneighboridx(i1,jcount),2) &
                                    = dsfuncdxyz_temp(invneighboridx(i1,jcount),2) + dyij + dyik
                                  dsfuncdxyz_temp(invneighboridx(i1,jcount),3) &
                                    = dsfuncdxyz_temp(invneighboridx(i1,jcount),3) + dzij + dzik

                                  n=lstc(i3)
                                  dsfuncdxyz_temp(invneighboridx(i1,n),1) &
                                    = dsfuncdxyz_temp(invneighboridx(i1,n),1) - dxij - dxjk 
                                  dsfuncdxyz_temp(invneighboridx(i1,n),2) &
                                    = dsfuncdxyz_temp(invneighboridx(i1,n),2) - dyij - dyjk
                                  dsfuncdxyz_temp(invneighboridx(i1,n),3) &
                                    = dsfuncdxyz_temp(invneighboridx(i1,n),3) - dzij - dzjk

                                  m=lstc(i4)
                                  dsfuncdxyz_temp(invneighboridx(i1,m),1) &
                                    = dsfuncdxyz_temp(invneighboridx(i1,m),1) - dxik + dxjk
                                  dsfuncdxyz_temp(invneighboridx(i1,m),2) &
                                    = dsfuncdxyz_temp(invneighboridx(i1,m),2) - dyik + dyjk
                                  dsfuncdxyz_temp(invneighboridx(i1,m),3) &
                                    = dsfuncdxyz_temp(invneighboridx(i1,m),3) - dzik + dzjk
                              end if ! ldoforces

                          end if ! rjk <= funccutoff_local
                      end if ! rik <= funccutoff_local
                   end if ! symelement_local i3 i4                        
              end do ! i4
          end if ! rij <= funccutoff_local
      end do ! i3
      symfunction_temp(i2) = symfunction_temp(i2) * pnorm
!!
      return
      end
