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
!! Purpose : calculating SF derivatives and stress components based on the type
!! called by:
!! - calconefunction_atomic_sfg.f90
!! - Emir

subroutine getsymfunctionderivatives(iatom,isfg,isf,eindex,jcount,xyzstruct_local,&
        invneighboridx,listdim,natoms,maxnum_funcvalues_local,maxnum_sfgroups_local,max_num_neigh_local,&
        max_num_atoms,group_neigh_local,neigh_idx_local,groupdim,maxdim,grouprij,&
        fcut_local,dfcut_local,grouptype,nelem,eta_local,rshift_local,lambda_local,&
        zeta_local,ftype,cutoff,lstb,dsfuncdxyz_temp,ldostress,strs_temp)

    use fileunits ! debug
    use nnconstants

    implicit none

    integer maxdim ! in
    integer groupdim ! in
    integer maxnum_funcvalues_local ! in
    integer maxnum_sfgroups_local ! in
    integer nelem ! in
    integer natoms ! in
    integer listdim ! in
    integer jcount ! in
    integer iatom ! in
    integer isf ! in
    integer isfg ! in
    integer eindex ! in
    integer ftype ! in
    integer grouptype ! in
    integer group_neigh_local(maxdim,grouptype+1)
    integer neigh_idx_local(maxdim,grouptype+1)
    integer max_group_listdim_local ! in
    integer max_num_atoms ! in
    integer max_num_neigh_local !in
    integer ineigh1,ineigh2 ! internal
    integer i1,i2,i3 ! internal
    integer idx1,idx2 ! internal
    integer invneighboridx(natoms,max_num_atoms) ! in

    real*8  temp1,temp2,temp3,temp4,temp5 ! internal
    real*8  cutoff ! in
    real*8  grouprij(maxdim,3)
    real*8  fcut_local(maxdim,3) ! in
    real*8  dfcut_local(maxdim,3) ! in
    real*8  lstb(listdim,4) ! in
    real*8  eta_local(maxnum_funcvalues_local,maxnum_sfgroups_local,nelem) ! in
    real*8  zeta_local(maxnum_funcvalues_local,maxnum_sfgroups_local,nelem) ! in
    real*8  rshift_local(maxnum_funcvalues_local,maxnum_sfgroups_local,nelem) ! in
    real*8  lambda_local(maxnum_funcvalues_local,maxnum_sfgroups_local,nelem) ! in
    real*8  drijdxi, drijdyi, drijdzi  ! internal ...
    real*8  drijdxj, drijdyj, drijdzj
    real*8  drijdxk, drijdyk, drijdzk
    real*8  drikdxi, drikdyi, drikdzi
    real*8  drikdxj, drikdyj, drikdzj
    real*8  drikdxk, drikdyk, drikdzk
    real*8  drjkdxi, drjkdyi, drjkdzi
    real*8  drjkdxj, drjkdyj, drjkdzj
    real*8  drjkdxk, drjkdyk, drjkdzk
    real*8  deltaxj,deltayj,deltazj
    real*8  deltaxk,deltayk,deltazk
    real*8  dfcutijdxi,dfcutijdyi,dfcutijdzi
    real*8  dfcutijdxj,dfcutijdyj,dfcutijdzj
    real*8  dfcutijdxk,dfcutijdyk,dfcutijdzk
    real*8  dfcutikdxi,dfcutikdyi,dfcutikdzi
    real*8  dfcutikdxj,dfcutikdyj,dfcutikdzj
    real*8  dfcutikdxk,dfcutikdyk,dfcutikdzk
    real*8  dfcutjkdxi,dfcutjkdyi,dfcutjkdzi
    real*8  dfcutjkdxj,dfcutjkdyj,dfcutjkdzj
    real*8  dfcutjkdxk,dfcutjkdyk,dfcutjkdzk ! ...internal
    real*8  xyzstruct_local(3,max_num_atoms) ! in
    real*8  dsfuncdxyz_temp(0:max_num_neigh_local,3) ! OUT
    real*8  strs_temp(3,3,maxnum_funcvalues_local)              ! OUT
    real*8  costheta ! internal..
    real*8  theta
    real*8  expxyz
    real*8  f
    real*8  g
    real*8  dgdxi,dgdyi,dgdzi
    real*8  dgdxj,dgdyj,dgdzj
    real*8  dgdxk,dgdyk,dgdzk
    real*8  dfdxi,dfdyi,dfdzi
    real*8  dfdxj,dfdyj,dfdzj
    real*8  dfdxk,dfdyk,dfdzk
    real*8  dcosthetadxi,dcosthetadyi,dcosthetadzi
    real*8  dcosthetadxj,dcosthetadyj,dcosthetadzj
    real*8  dcosthetadxk,dcosthetadyk,dcosthetadzk
    real*8  dthetadxi,dthetadyi,dthetadzi
    real*8  dthetadxj,dthetadyj,dthetadzj
    real*8  dthetadxk,dthetadyk,dthetadzk
    real*8  dexpxyzdxi,dexpxyzdyi,dexpxyzdzi
    real*8  dexpxyzdxj,dexpxyzdyj,dexpxyzdzj
    real*8  dexpxyzdxk,dexpxyzdyk,dexpxyzdzk
    real*8  dxij,dyij,dzij,dxik,dyik,dzik,dxjk,dyjk,dzjk
    real*8  r2ij,r2ik,r2jk,rinvijik
    real*8  pzl,pnorm,plambda_local,pfc,pexp,p1,p2,p3,p2etaplambda,fg

    logical ldostress !in
    if(ftype.eq.1)then ! radial SF (only cutoff function)

        do i1 = 1,groupdim
            ineigh1 = group_neigh_local(i1,1)
            idx1 = neigh_idx_local(i1,1)
            if(grouprij(i1,1).le.cutoff)then

                deltaxj = -1.d0 * (xyzstruct_local(1,jcount)-lstb(ineigh1,1))
                deltayj = -1.d0 * (xyzstruct_local(2,jcount)-lstb(ineigh1,2))
                deltazj = -1.d0 * (xyzstruct_local(3,jcount)-lstb(ineigh1,3))
                drijdxi = -deltaxj / grouprij(i1,1)
                drijdyi = -deltayj / grouprij(i1,1)
                drijdzi = -deltazj / grouprij(i1,1)
                drijdxj=-1.d0 * drijdxi
                drijdyj=-1.d0 * drijdyi
                drijdzj=-1.d0 * drijdzi

                temp1 = -0.5d0*dsin(pi*grouprij(i1,1)/cutoff)*pi/cutoff

                !! dsfunc/dx
                dsfuncdxyz_temp(invneighboridx(iatom,jcount),1)=dsfuncdxyz_temp(invneighboridx(iatom,jcount),1)&
                        +(temp1*drijdxi)
                dsfuncdxyz_temp(invneighboridx(iatom,idx1),1) =dsfuncdxyz_temp(invneighboridx(iatom,idx1),1)&
                        +(temp1*drijdxj)
                !! dsfunc/dy
                dsfuncdxyz_temp(invneighboridx(iatom,jcount),2)=dsfuncdxyz_temp(invneighboridx(iatom,jcount),2)&
                        +(temp1*drijdyi)
                dsfuncdxyz_temp(invneighboridx(iatom,idx1),2) =dsfuncdxyz_temp(invneighboridx(iatom,idx1),2)&
                        +(temp1*drijdyj)
                !! dsfunc/dz
                dsfuncdxyz_temp(invneighboridx(iatom,jcount),3)=dsfuncdxyz_temp(invneighboridx(iatom,jcount),3)&
                        +(temp1*drijdzi)
                dsfuncdxyz_temp(invneighboridx(iatom,idx1),3) =dsfuncdxyz_temp(invneighboridx(iatom,idx1),3)&
                        +(temp1*drijdzj)

                if(ldostress)then
                    strs_temp(1,1,isf)=strs_temp(1,1,isf)&
                            +deltaxj*(drijdxj*temp1)
                    strs_temp(2,1,isf)=strs_temp(2,1,isf)&
                            +deltayj*(drijdxj*temp1)
                    strs_temp(3,1,isf)=strs_temp(3,1,isf)&
                            +deltazj*(drijdxj*temp1)
                    !! dsfunc/dy
                    strs_temp(1,2,isf)=strs_temp(1,2,isf)&
                            +deltaxj*(drijdyj*temp1)
                    strs_temp(2,2,isf)=strs_temp(2,2,isf)&
                            +deltayj*(drijdyj*temp1)
                    strs_temp(3,2,isf)=strs_temp(3,2,isf)&
                            +deltazj*(drijdyj*temp1)
                    !! dsfunc/dz
                    strs_temp(1,3,isf)=strs_temp(1,3,isf)&
                            +deltaxj*(drijdzj*temp1)
                    strs_temp(2,3,isf)=strs_temp(2,3,isf)&
                            +deltayj*(drijdzj*temp1)
                    strs_temp(3,3,isf)=strs_temp(3,3,isf)&
                            +deltazj*(drijdzj*temp1)
                end if

            end if
        end do

    elseif(ftype.eq.2)then ! radial SF (standard)

        do i1 = 1,groupdim
            ineigh1 = group_neigh_local(i1,1)
            idx1 = neigh_idx_local(i1,1)
            if(grouprij(i1,1).le.cutoff)then

                deltaxj = -1.d0 * (xyzstruct_local(1,jcount)-lstb(ineigh1,1))
                deltayj = -1.d0 * (xyzstruct_local(2,jcount)-lstb(ineigh1,2))
                deltazj = -1.d0 * (xyzstruct_local(3,jcount)-lstb(ineigh1,3))
                drijdxi = -deltaxj / grouprij(i1,1)
                drijdyi = -deltayj / grouprij(i1,1)
                drijdzi = -deltazj / grouprij(i1,1)
                drijdxj=-1.d0 * drijdxi
                drijdyj=-1.d0 * drijdyi
                drijdzj=-1.d0 * drijdzi

                dfcutijdxi = dfcut_local(i1,1) * drijdxi
                dfcutijdyi = dfcut_local(i1,1) * drijdyi
                dfcutijdzi = dfcut_local(i1,1) * drijdzi
                dfcutijdxj = -1.d0 * dfcutijdxi
                dfcutijdyj = -1.d0 * dfcutijdyi
                dfcutijdzj = -1.d0 * dfcutijdzi

                !! Calculation of derivatives for forces
                temp1 = -2.d0 * eta_local(isf,isfg,eindex) * (grouprij(i1,1)-rshift_local(isf,isfg,eindex)) * &
                        dexp(-1.d0*eta_local(isf,isfg,eindex) * (grouprij(i1,1)-rshift_local(isf,isfg,eindex))**2) * &
                        fcut_local(i1,1)
                temp2 = dexp(-1.d0*eta_local(isf,isfg,eindex)*(grouprij(i1,1)-rshift_local(isf,isfg,eindex))**2)

                !! dsfunc/dx
                dsfuncdxyz_temp(invneighboridx(iatom,jcount),1) = dsfuncdxyz_temp(invneighboridx(iatom,jcount),1)+&
                        (drijdxi*temp1 + temp2*dfcutijdxi)
                dsfuncdxyz_temp(invneighboridx(iatom,idx1),1) = dsfuncdxyz_temp(invneighboridx(iatom,idx1),1)+&
                        (drijdxj*temp1 + temp2*dfcutijdxj)
                !! dsfunc/dy
                dsfuncdxyz_temp(invneighboridx(iatom,jcount),2)= dsfuncdxyz_temp(invneighboridx(iatom,jcount),2)+&
                        (drijdyi*temp1 + temp2*dfcutijdyi)
                dsfuncdxyz_temp(invneighboridx(iatom,idx1),2) = dsfuncdxyz_temp(invneighboridx(iatom,idx1),2)+&
                        (drijdyj*temp1 + temp2*dfcutijdyj)
                !! dsfunc/dz
                dsfuncdxyz_temp(invneighboridx(iatom,jcount),3) = dsfuncdxyz_temp(invneighboridx(iatom,jcount),3)+&
                        (drijdzi*temp1 + temp2*dfcutijdzi)
                dsfuncdxyz_temp(invneighboridx(iatom,idx1),3) = dsfuncdxyz_temp(invneighboridx(iatom,idx1),3)+&
                        (drijdzj*temp1 + temp2*dfcutijdzj)

                if(ldostress)then
                    !! dsfunc/dx
                    strs_temp(1,1,isf)=strs_temp(1,1,isf)+deltaxj*&
                            (drijdxj*temp1&
                                    + temp2*dfcutijdxj)
                    strs_temp(2,1,isf)=strs_temp(2,1,isf)+deltayj*&
                            (drijdxj*temp1&
                                    + temp2*dfcutijdxj)
                    strs_temp(3,1,isf)=strs_temp(3,1,isf)+deltazj*&
                            (drijdxj*temp1&
                                    + temp2*dfcutijdxj)
                    !! dsfunc/dy
                    strs_temp(1,2,isf)=strs_temp(1,2,isf)+deltaxj*&
                            (drijdyj*temp1&
                                    + temp2*dfcutijdyj)
                    strs_temp(2,2,isf)=strs_temp(2,2,isf)+deltayj*&
                            (drijdyj*temp1&
                                    + temp2*dfcutijdyj)
                    strs_temp(3,2,isf)=strs_temp(3,2,isf)+deltazj*&
                            (drijdyj*temp1&
                                    + temp2*dfcutijdyj)
                    !! dsfunc/dz
                    strs_temp(1,3,isf)=strs_temp(1,3,isf)+deltaxj*&
                            (drijdzj*temp1&
                                    + temp2*dfcutijdzj)
                    strs_temp(2,3,isf)=strs_temp(2,3,isf)+deltayj*&
                            (drijdzj*temp1&
                                    + temp2*dfcutijdzj)
                    strs_temp(3,3,isf)=strs_temp(3,3,isf)+deltazj*&
                            (drijdzj*temp1&
                                    + temp2*dfcutijdzj)
                end if

            end if
        end do

    elseif(ftype.eq.4)then ! radial SF

        do i1 = 1,groupdim
            ineigh1 = group_neigh_local(i1,1)
            idx1 = neigh_idx_local(i1,1)
            if(grouprij(i1,1).le.cutoff)then

                deltaxj = -1.d0 * (xyzstruct_local(1,jcount)-lstb(ineigh1,1))
                deltayj = -1.d0 * (xyzstruct_local(2,jcount)-lstb(ineigh1,2))
                deltazj = -1.d0 * (xyzstruct_local(3,jcount)-lstb(ineigh1,3))
                drijdxi = -deltaxj / grouprij(i1,1)
                drijdyi = -deltayj / grouprij(i1,1)
                drijdzi = -deltazj / grouprij(i1,1)
                drijdxj=-1.d0 * drijdxi
                drijdyj=-1.d0 * drijdyi
                drijdzj=-1.d0 * drijdzi

                dfcutijdxi = dfcut_local(i1,1) * drijdxi
                dfcutijdyi = dfcut_local(i1,1) * drijdyi
                dfcutijdzi = dfcut_local(i1,1) * drijdzi
                dfcutijdxj = -1.d0 * dfcutijdxi
                dfcutijdyj = -1.d0 * dfcutijdyi
                dfcutijdzj = -1.d0 * dfcutijdzi

                !! Calculation of derivatives for forces
                temp1 = -1.d0 * eta_local(isf,isfg,eindex) * dsin(eta_local(isf,isfg,eindex)*grouprij(i1,1))*&
                        fcut_local(i1,1)
                temp2 = dcos(eta_local(isf,isfg,eindex)*grouprij(i1,1))

                !! dsfunc/dx
                dsfuncdxyz_temp(invneighboridx(iatom,jcount),1) = dsfuncdxyz_temp(invneighboridx(iatom,jcount),1)+&
                        (drijdxi*temp1 + temp2*dfcutijdxi)
                dsfuncdxyz_temp(invneighboridx(iatom,idx1),1) = dsfuncdxyz_temp(invneighboridx(iatom,idx1),1)+&
                        (drijdxj*temp1 + temp2*dfcutijdxj)
                !! dsfunc/dy
                dsfuncdxyz_temp(invneighboridx(iatom,jcount),2)= dsfuncdxyz_temp(invneighboridx(iatom,jcount),2)+&
                        (drijdyi*temp1 + temp2*dfcutijdyi)
                dsfuncdxyz_temp(invneighboridx(iatom,idx1),2) = dsfuncdxyz_temp(invneighboridx(iatom,idx1),2)+&
                        (drijdyj*temp1 + temp2*dfcutijdyj)
                !! dsfunc/dz
                dsfuncdxyz_temp(invneighboridx(iatom,jcount),3) = dsfuncdxyz_temp(invneighboridx(iatom,jcount),3)+&
                        (drijdzi*temp1 + temp2*dfcutijdzi)
                dsfuncdxyz_temp(invneighboridx(iatom,idx1),3) = dsfuncdxyz_temp(invneighboridx(iatom,idx1),3)+&
                        (drijdzj*temp1 + temp2*dfcutijdzj)

                if(ldostress)then
                    !! dsfunc/dx
                    strs_temp(1,1,isf)=strs_temp(1,1,isf)+deltaxj*&
                            (drijdxj*temp1&
                                    + temp2*dfcutijdxj)
                    strs_temp(2,1,isf)=strs_temp(2,1,isf)+deltayj*&
                            (drijdxj*temp1&
                                    + temp2*dfcutijdxj)
                    strs_temp(3,1,isf)=strs_temp(3,1,isf)+deltazj*&
                            (drijdxj*temp1&
                                    + temp2*dfcutijdxj)
                    !! dsfunc/dy
                    strs_temp(1,2,isf)=strs_temp(1,2,isf)+deltaxj*&
                            (drijdyj*temp1&
                                    + temp2*dfcutijdyj)
                    strs_temp(2,2,isf)=strs_temp(2,2,isf)+deltayj*&
                            (drijdyj*temp1&
                                    + temp2*dfcutijdyj)
                    strs_temp(3,2,isf)=strs_temp(3,2,isf)+deltazj*&
                            (drijdyj*temp1&
                                    + temp2*dfcutijdyj)
                    !! dsfunc/dz
                    strs_temp(1,3,isf)=strs_temp(1,3,isf)+deltaxj*&
                            (drijdzj*temp1&
                                    + temp2*dfcutijdzj)
                    strs_temp(2,3,isf)=strs_temp(2,3,isf)+deltayj*&
                            (drijdzj*temp1&
                                    + temp2*dfcutijdzj)
                    strs_temp(3,3,isf)=strs_temp(3,3,isf)+deltazj*&
                            (drijdzj*temp1&
                                    + temp2*dfcutijdzj)
                end if

            end if
        end do

    elseif(ftype.eq.3)then ! angular SF (narrow)

        pnorm = 2.d0**(1.d0-zeta_local(isf,isfg,eindex))
        pzl = zeta_local(isf,isfg,eindex)*lambda_local(isf,isfg,eindex)

        do i1 = 1,groupdim
            ineigh1 = group_neigh_local(i1,1)
            ineigh2 = group_neigh_local(i1,2)
            idx1 = neigh_idx_local(i1,1)
            idx2 = neigh_idx_local(i1,2)
            if(grouprij(i1,1).le.cutoff.and.grouprij(i1,2).le.cutoff.and.grouprij(i1,3).le.cutoff)then

                rinvijik = 1.d0/grouprij(i1,1)/grouprij(i1,2)
                r2ij = grouprij(i1,1)*grouprij(i1,1)
                r2ik = grouprij(i1,2)*grouprij(i1,2)
                r2jk = grouprij(i1,3)*grouprij(i1,3)
                deltaxj = xyzstruct_local(1,jcount) - lstb(ineigh1,1)
                deltayj = xyzstruct_local(2,jcount) - lstb(ineigh1,2)
                deltazj = xyzstruct_local(3,jcount) - lstb(ineigh1,3)
                deltaxk = xyzstruct_local(1,jcount) - lstb(ineigh2,1)
                deltayk = xyzstruct_local(2,jcount) - lstb(ineigh2,2)
                deltazk = xyzstruct_local(3,jcount) - lstb(ineigh2,3)
                dxjk = lstb(ineigh1,1) - lstb(ineigh2,1)
                dyjk = lstb(ineigh1,2) - lstb(ineigh2,2)
                dzjk = lstb(ineigh1,3) - lstb(ineigh2,3)

                costheta = deltaxj*deltaxk + deltayj*deltayk + deltazj*deltazk
                costheta = costheta*rinvijik
                pfc = fcut_local(i1,1)* fcut_local(i1,2) * fcut_local(i1,3)
                pexp = exp(-eta_local(isf,isfg,eindex)* (r2ij+r2ik+r2jk))
                plambda_local = 1.d0 + lambda_local(isf,isfg,eindex)*costheta

                if(plambda_local.le.0.d0)then
                    fg = 0.d0
                else
                    fg = plambda_local**(zeta_local(isf,isfg,eindex)-1.d0)*pexp
                end if

                fg = fg*pnorm
                rinvijik = rinvijik*pzl
                costheta = costheta*pzl
                p2etaplambda =2.d0*eta_local(isf,isfg,eindex)*plambda_local
                p1 = pfc*(rinvijik - costheta/r2ij - p2etaplambda) + &
                        fcut_local(i1,2)*fcut_local(i1,3)*dfcut_local(i1,1)*plambda_local/grouprij(i1,1)
                p2 = pfc*(rinvijik - costheta/r2ik - p2etaplambda) + &
                        fcut_local(i1,1)*fcut_local(i1,3)*dfcut_local(i1,2)*plambda_local/grouprij(i1,2)
                p3 = pfc*(rinvijik                 + p2etaplambda) - &
                        fcut_local(i1,1)*fcut_local(i1,2)*dfcut_local(i1,3)*plambda_local/grouprij(i1,3)
                p1 = p1*fg
                p2 = p2*fg
                p3 = p3*fg
                dxij = deltaxj*p1
                dyij = deltayj*p1
                dzij = deltazj*p1
                dxik = deltaxk*p2
                dyik = deltayk*p2
                dzik = deltazk*p2
                dxjk = dxjk*p3
                dyjk = dyjk*p3
                dzjk = dzjk*p3

                dsfuncdxyz_temp(invneighboridx(iatom,jcount),1) &
                         = dsfuncdxyz_temp(invneighboridx(iatom,jcount),1) + dxij + dxik
                dsfuncdxyz_temp(invneighboridx(iatom,jcount),2) &
                        = dsfuncdxyz_temp(invneighboridx(iatom,jcount),2) + dyij + dyik
                dsfuncdxyz_temp(invneighboridx(iatom,jcount),3) &
                        = dsfuncdxyz_temp(invneighboridx(iatom,jcount),3) + dzij + dzik

                dsfuncdxyz_temp(invneighboridx(iatom,idx1),1) &
                        = dsfuncdxyz_temp(invneighboridx(iatom,idx1),1) - dxij - dxjk !temp2x
                dsfuncdxyz_temp(invneighboridx(iatom,idx1),2) &
                        = dsfuncdxyz_temp(invneighboridx(iatom,idx1),2) - dyij - dyjk !temp2y
                dsfuncdxyz_temp(invneighboridx(iatom,idx1),3) &
                        = dsfuncdxyz_temp(invneighboridx(iatom,idx1),3) - dzij - dzjk !temp2z

                dsfuncdxyz_temp(invneighboridx(iatom,idx2),1) &
                        = dsfuncdxyz_temp(invneighboridx(iatom,idx2),1) - dxik + dxjk !temp3x
                dsfuncdxyz_temp(invneighboridx(iatom,idx2),2) &
                        = dsfuncdxyz_temp(invneighboridx(iatom,idx2),2) - dyik + dyjk !temp3y
                dsfuncdxyz_temp(invneighboridx(iatom,idx2),3) &
                        = dsfuncdxyz_temp(invneighboridx(iatom,idx2),3) - dzik + dzjk !temp3z

                if(ldostress)then
                    strs_temp(1,1,isf)=strs_temp(1,1,isf) + deltaxj*(dxij+dxjk)&
                            - deltaxk*(dxjk - dxik)
                    strs_temp(2,1,isf)=strs_temp(2,1,isf) + deltayj*(dxij+dxjk)&
                            - deltayk*(dxjk - dxik)
                    strs_temp(3,1,isf)=strs_temp(3,1,isf) + deltazj*(dxij+dxjk)&
                            - deltazk*(dxjk - dxik)

                    strs_temp(1,2,isf)=strs_temp(1,2,isf) + deltaxj*(dyij+dyjk)&
                            - deltaxk*(dyjk - dyik)
                    strs_temp(2,2,isf)=strs_temp(2,2,isf) + deltayj*(dyij+dyjk)&
                            - deltayk*(dyjk - dyik)
                    strs_temp(3,2,isf)=strs_temp(3,2,isf) + deltazj*(dyij+dyjk)&
                            - deltazk*(dyjk - dyik)

                    strs_temp(1,3,isf)=strs_temp(1,3,isf) + deltaxj*(dzij+dzjk)&
                             -deltaxk*(dzjk - dzik)
                    strs_temp(2,3,isf)=strs_temp(2,3,isf) + deltayj*(dzij+dzjk)&
                            -deltayk*(dzjk - dzik)
                    strs_temp(3,3,isf)=strs_temp(3,3,isf) + deltazj*(dzij+dzjk)&
                             -deltazk*(dzjk - dzik)
                end if


            end if
        end do

    elseif(ftype.eq.8)then ! angular SF (shifted)

        do i1 = 1,groupdim
            ineigh1 = group_neigh_local(i1,1)
            ineigh2 = group_neigh_local(i1,2)
            idx1 = neigh_idx_local(i1,1)
            idx2 = neigh_idx_local(i1,2)
            if(grouprij(i1,1).le.cutoff.and.grouprij(i1,2).le.cutoff.and.grouprij(i1,3).le.cutoff)then

                ! r_ij
                deltaxj = -1.d0 * (xyzstruct_local(1,jcount)-lstb(ineigh1,1))
                deltayj = -1.d0 * (xyzstruct_local(2,jcount)-lstb(ineigh1,2))
                deltazj = -1.d0 * (xyzstruct_local(3,jcount)-lstb(ineigh1,3))
                drijdxi = -deltaxj/grouprij(i1,1)
                drijdyi = -deltayj/grouprij(i1,1)
                drijdzi = -deltazj/grouprij(i1,1)
                drijdxj =- 1.d0*drijdxi
                drijdyj = -1.d0*drijdyi
                drijdzj = -1.d0*drijdzi
                drijdxk=0.0d0
                drijdyk=0.0d0
                drijdzk=0.0d0

                dfcutijdxi = dfcut_local(i1,1)*drijdxi
                dfcutijdyi = dfcut_local(i1,1)*drijdyi
                dfcutijdzi = dfcut_local(i1,1)*drijdzi
                dfcutijdxj = -1.d0*dfcutijdxi
                dfcutijdyj = -1.d0*dfcutijdyi
                dfcutijdzj = -1.d0*dfcutijdzi
                dfcutijdxk = 0.0d0
                dfcutijdyk = 0.0d0
                dfcutijdzk = 0.0d0

                ! r_ik
                deltaxk = -1.d0 * (xyzstruct_local(1,jcount)-lstb(ineigh2,1))
                deltayk = -1.d0 * (xyzstruct_local(2,jcount)-lstb(ineigh2,2))
                deltazk = -1.d0 * (xyzstruct_local(3,jcount)-lstb(ineigh2,3))
                drikdxi = -deltaxk/grouprij(i1,2)
                drikdyi = -deltayk/grouprij(i1,2)
                drikdzi = -deltazk/grouprij(i1,2)
                drikdxk =- 1.d0*drikdxi
                drikdyk = -1.d0*drikdyi
                drikdzk = -1.d0*drikdzi
                drikdxj = 0.0d0
                drikdyj = 0.0d0
                drikdzj = 0.0d0

                dfcutikdxi = dfcut_local(i1,2)*drikdxi
                dfcutikdyi = dfcut_local(i1,2)*drikdyi
                dfcutikdzi = dfcut_local(i1,2)*drikdzi
                dfcutikdxk = -1.d0*dfcutikdxi
                dfcutikdyk = -1.d0*dfcutikdyi
                dfcutikdzk = -1.d0*dfcutikdzi
                dfcutikdxj = 0.0d0
                dfcutikdyj = 0.0d0
                dfcutikdzj = 0.0d0

                ! r_jk
                drjkdxj = 1.d0*(lstb(ineigh1,1)-lstb(ineigh2,1))/grouprij(i1,3)
                drjkdyj = 1.d0*(lstb(ineigh1,2)-lstb(ineigh2,2))/grouprij(i1,3)
                drjkdzj = 1.d0*(lstb(ineigh1,3)-lstb(ineigh2,3))/grouprij(i1,3)
                drjkdxk =- 1.d0*drjkdxj
                drjkdyk = -1.d0*drjkdyj
                drjkdzk = -1.d0*drjkdzj
                drjkdxi = 0.0d0
                drjkdyi = 0.0d0
                drjkdzi = 0.0d0

                dfcutjkdxj = dfcut_local(i1,3)*drjkdxj
                dfcutjkdyj = dfcut_local(i1,3)*drjkdyj
                dfcutjkdzj = dfcut_local(i1,3)*drjkdzj
                dfcutjkdxk = -1.d0*dfcutjkdxj
                dfcutjkdyk = -1.d0*dfcutjkdyj
                dfcutjkdzk = -1.d0*dfcutjkdzj
                dfcutjkdxi = 0.0d0
                dfcutjkdyi = 0.0d0
                dfcutjkdzi = 0.0d0

                !! costheta=(rjk**2 - rij**2 -rik**2)/(-2.d0*rij*rik)
                f = grouprij(i1,3)**2 - grouprij(i1,1)**2 - grouprij(i1,2)**2
                g = -2.d0 * grouprij(i1,1) * grouprij(i1,2)
                !!                costheta=f/g
                !! filter out numerical noise that results in NaN for acos(x):
                temp1=f/g
                if(temp1.gt.1.0d0)then
                    if(temp1.gt.1.1d0)then
                        write(ounit,*)'Error in calconefunction type 8 ',temp1
                        stop
                    endif
                    temp1=1.0d0
                endif
                if(temp1.lt.-1.0d0)then
                    if(temp1.lt.-1.1d0)then
                        write(ounit,*)'Error in calconefunction type 8 ',temp1
                        stop
                    endif
                    temp1=-1.0d0
                endif
                theta = acos(temp1)
                !! convert to degree
                theta = theta*rad2deg

                ! Cosine derivative
                dfdxi = (-2.d0*grouprij(i1,1)*drijdxi - 2.d0*grouprij(i1,2)*drikdxi)
                dfdyi = (-2.d0*grouprij(i1,1)*drijdyi - 2.d0*grouprij(i1,2)*drikdyi)
                dfdzi = (-2.d0*grouprij(i1,1)*drijdzi - 2.d0*grouprij(i1,2)*drikdzi)
                !!
                dfdxj = (2.d0*grouprij(i1,3)*drjkdxj - 2.d0*grouprij(i1,1)*drijdxj)
                dfdyj = (2.d0*grouprij(i1,3)*drjkdyj - 2.d0*grouprij(i1,1)*drijdyj)
                dfdzj = (2.d0*grouprij(i1,3)*drjkdzj - 2.d0*grouprij(i1,1)*drijdzj)
                !!
                dfdxk = (2.d0*grouprij(i1,3)*drjkdxk - 2.d0*grouprij(i1,2)*drikdxk)
                dfdyk = (2.d0*grouprij(i1,3)*drjkdyk - 2.d0*grouprij(i1,2)*drikdyk)
                dfdzk = (2.d0*grouprij(i1,3)*drjkdzk - 2.d0*grouprij(i1,2)*drikdzk)

                !!
                dgdxi = -2.d0*(drijdxi*grouprij(i1,2) + grouprij(i1,1)*drikdxi)
                dgdyi = -2.d0*(drijdyi*grouprij(i1,2) + grouprij(i1,1)*drikdyi)
                dgdzi = -2.d0*(drijdzi*grouprij(i1,2) + grouprij(i1,1)*drikdzi)
                !!
                dgdxj = -2.d0*drijdxj*grouprij(i1,2)
                dgdyj = -2.d0*drijdyj*grouprij(i1,2)
                dgdzj = -2.d0*drijdzj*grouprij(i1,2)
                !!
                dgdxk = -2.d0*grouprij(i1,1)*drikdxk
                dgdyk = -2.d0*grouprij(i1,1)*drikdyk
                dgdzk = -2.d0*grouprij(i1,1)*drikdzk

                !! derivative of acos:
                !! f(x)=acos(x)
                !! f'(x)=-1/sqrt(1-x**2)
                !!
                !! filter out numerical noise that results in Infty for temp4:
                !!                 temp4=-1.d0/dsqrt(1.d0-f**2/g**2)  ! this line has problems with numerical noise
                temp4=f**2/g**2
                !!                 if(temp4.gt.1.0d0)temp4=1.0d0
                temp4=-1.d0/dsqrt(1.d0-temp4)
                dthetadxi=rad2deg*temp4*(dfdxi*g - f*dgdxi)/g**2
                dthetadyi=rad2deg*temp4*(dfdyi*g - f*dgdyi)/g**2
                dthetadzi=rad2deg*temp4*(dfdzi*g - f*dgdzi)/g**2
                dthetadxj=rad2deg*temp4*(dfdxj*g - f*dgdxj)/g**2
                dthetadyj=rad2deg*temp4*(dfdyj*g - f*dgdyj)/g**2
                dthetadzj=rad2deg*temp4*(dfdzj*g - f*dgdzj)/g**2
                dthetadxk=rad2deg*temp4*(dfdxk*g - f*dgdxk)/g**2
                dthetadyk=rad2deg*temp4*(dfdyk*g - f*dgdyk)/g**2
                dthetadzk=rad2deg*temp4*(dfdzk*g - f*dgdzk)/g**2

                !! Calculation of derivatives for forces
                !! (fghi)'= f'ghi + fg'hi + fgh'i + fghi'
                temp1=theta-rshift_local(isf,isfg,eindex)
                temp2=theta+rshift_local(isf,isfg,eindex)
                expxyz=(dexp(-eta_local(isf,isfg,eindex)*(temp1       )**2)&
                        +dexp(-eta_local(isf,isfg,eindex)*(temp2-360.d0)**2)&
                        +dexp(-eta_local(isf,isfg,eindex)*(temp2       )**2)&
                        +dexp(-eta_local(isf,isfg,eindex)*(temp1-360.d0)**2))
                temp3=-2.d0*eta_local(isf,isfg,eindex)*temp1*dexp(-eta_local(isf,isfg,eindex)*(temp1)**2)&
                        -2.d0*eta_local(isf,isfg,eindex)*temp2*dexp(-eta_local(isf,isfg,eindex)*(temp2-360.d0)**2)&
                        -2.d0*eta_local(isf,isfg,eindex)*temp2*dexp(-eta_local(isf,isfg,eindex)*(temp2       )**2)&
                        -2.d0*eta_local(isf,isfg,eindex)*temp1*dexp(-eta_local(isf,isfg,eindex)*(temp1-360.d0)**2)
                !!
                dexpxyzdxi=temp3*dthetadxi
                dexpxyzdyi=temp3*dthetadyi
                dexpxyzdzi=temp3*dthetadzi
                dexpxyzdxj=temp3*dthetadxj
                dexpxyzdyj=temp3*dthetadyj
                dexpxyzdzj=temp3*dthetadzj
                dexpxyzdxk=temp3*dthetadxk
                dexpxyzdyk=temp3*dthetadyk
                dexpxyzdzk=temp3*dthetadzk

                ! dsfunc/dx
                temp1=(dexpxyzdxi* fcut_local(i1,1)   * fcut_local(i1,2)   * fcut_local(i1,3)&
                        +expxyz   * dfcutijdxi* fcut_local(i1,2)   * fcut_local(i1,3)&
                        +expxyz   * fcut_local(i1,1)   *dfcutikdxi* fcut_local(i1,3)&
                        +expxyz   * fcut_local(i1,1)   * fcut_local(i1,2)   *dfcutjkdxi)
                temp2=(dexpxyzdxj* fcut_local(i1,1)   * fcut_local(i1,2)   * fcut_local(i1,3)&
                        +expxyz   *dfcutijdxj* fcut_local(i1,2)   * fcut_local(i1,3)&
                        +expxyz   * fcut_local(i1,1)   *dfcutikdxj* fcut_local(i1,3)&
                        +expxyz   * fcut_local(i1,1)   * fcut_local(i1,2)   *dfcutjkdxj)
                temp3=(dexpxyzdxk* fcut_local(i1,1)   * fcut_local(i1,2)   * fcut_local(i1,3)&
                        +expxyz   *dfcutijdxk* fcut_local(i1,2)   * fcut_local(i1,3)&
                        +expxyz   * fcut_local(i1,1)   *dfcutikdxk* fcut_local(i1,3)&
                        +expxyz   * fcut_local(i1,1)   * fcut_local(i1,2)   *dfcutjkdxk)
                dsfuncdxyz_temp(invneighboridx(iatom,jcount),1)=dsfuncdxyz_temp(invneighboridx(iatom,jcount),1)+&
                        temp1
                dsfuncdxyz_temp(invneighboridx(iatom,idx1),1) =dsfuncdxyz_temp(invneighboridx(iatom,idx1),1)+&
                        temp2
                dsfuncdxyz_temp(invneighboridx(iatom,idx2),1) =dsfuncdxyz_temp(invneighboridx(iatom,idx2),1)+&
                        temp3

                if(ldostress)then
                    strs_temp(1,1,isf)=strs_temp(1,1,isf)+deltaxj*temp2&
                            +deltaxk*temp3
                    strs_temp(2,1,isf)=strs_temp(2,1,isf)+deltayj*temp2&
                            +deltayk*temp3
                    strs_temp(3,1,isf)=strs_temp(3,1,isf)+deltazj*temp2&
                            +deltazk*temp3
                end if


                ! dsfunc/dy
                temp1=(dexpxyzdyi* fcut_local(i1,1)   * fcut_local(i1,2)   * fcut_local(i1,3)&
                        +expxyz   *dfcutijdyi* fcut_local(i1,2)   * fcut_local(i1,3)&
                        +expxyz   * fcut_local(i1,1)   *dfcutikdyi* fcut_local(i1,3)&
                        +expxyz   * fcut_local(i1,1)   * fcut_local(i1,2)   *dfcutjkdyi)
                temp2=(dexpxyzdyj* fcut_local(i1,1)   * fcut_local(i1,2)   * fcut_local(i1,3)&
                        +expxyz   *dfcutijdyj* fcut_local(i1,2)   * fcut_local(i1,3)&
                        +expxyz   * fcut_local(i1,1)   *dfcutikdyj* fcut_local(i1,3)&
                        +expxyz   * fcut_local(i1,1)   * fcut_local(i1,2)   *dfcutjkdyj)
                temp3=(dexpxyzdyk* fcut_local(i1,1)   * fcut_local(i1,2)   * fcut_local(i1,3)&
                        +expxyz   *dfcutijdyk* fcut_local(i1,2)   * fcut_local(i1,3)&
                        +expxyz   * fcut_local(i1,1)   *dfcutikdyk* fcut_local(i1,3)&
                        +expxyz   * fcut_local(i1,1)   * fcut_local(i1,2)   *dfcutjkdyk)
                dsfuncdxyz_temp(invneighboridx(iatom,jcount),2)=dsfuncdxyz_temp(invneighboridx(iatom,jcount),2)+&
                        temp1
                dsfuncdxyz_temp(invneighboridx(iatom,idx1),2) =dsfuncdxyz_temp(invneighboridx(iatom,idx1),2)+&
                        temp2
                dsfuncdxyz_temp(invneighboridx(iatom,idx2),2) =dsfuncdxyz_temp(invneighboridx(iatom,idx2),2)+&
                        temp3

                if(ldostress)then
                    strs_temp(1,2,isf)=strs_temp(1,2,isf)+deltaxj*temp2&
                            +deltaxk*temp3
                    strs_temp(2,2,isf)=strs_temp(2,2,isf)+deltayj*temp2&
                            +deltayk*temp3
                    strs_temp(3,2,isf)=strs_temp(3,2,isf)+deltazj*temp2&
                            +deltazk*temp3
                end if

                ! dsfunc/dz
                temp1=(dexpxyzdzi* fcut_local(i1,1)   * fcut_local(i1,2)   * fcut_local(i1,3)&
                        + expxyz   *dfcutijdzi* fcut_local(i1,2)   * fcut_local(i1,3)&
                        + expxyz   * fcut_local(i1,1)   *dfcutikdzi* fcut_local(i1,3)&
                        + expxyz   * fcut_local(i1,1)   * fcut_local(i1,2)   *dfcutjkdzi)
                temp2=(dexpxyzdzj* fcut_local(i1,1)   * fcut_local(i1,2)   * fcut_local(i1,3)&
                        + expxyz   *dfcutijdzj* fcut_local(i1,2)   * fcut_local(i1,3)&
                        + expxyz   * fcut_local(i1,1)   *dfcutikdzj* fcut_local(i1,3)&
                        + expxyz   * fcut_local(i1,1)   * fcut_local(i1,2)   *dfcutjkdzj)
                temp3=(dexpxyzdzk* fcut_local(i1,1)   * fcut_local(i1,2)   * fcut_local(i1,3)&
                        + expxyz   *dfcutijdzk* fcut_local(i1,2)   * fcut_local(i1,3)&
                        + expxyz   * fcut_local(i1,1)   *dfcutikdzk* fcut_local(i1,3)&
                        + expxyz   * fcut_local(i1,1)   * fcut_local(i1,2)   *dfcutjkdzk)
                dsfuncdxyz_temp(invneighboridx(iatom,jcount),3)=dsfuncdxyz_temp(invneighboridx(iatom,jcount),3)+&
                        temp1
                dsfuncdxyz_temp(invneighboridx(iatom,idx1),3) =dsfuncdxyz_temp(invneighboridx(iatom,idx1),3)+&
                        temp2
                dsfuncdxyz_temp(invneighboridx(iatom,idx2),3) =dsfuncdxyz_temp(invneighboridx(iatom,idx2),3)+&
                        temp3

                if(ldostress)then
                    strs_temp(1,3,isf)=strs_temp(1,3,isf)+deltaxj*temp2&
                            +deltaxk*temp3
                    strs_temp(2,3,isf)=strs_temp(2,3,isf)+deltayj*temp2&
                            +deltayk*temp3
                    strs_temp(3,3,isf)=strs_temp(3,3,isf)+deltazj*temp2&
                            +deltazk*temp3
                end if

            end if
        end do

    elseif(ftype.eq.9)then ! angular SF (wide)
        do i1 = 1,groupdim
            ineigh1 = group_neigh_local(i1,1)
            ineigh2 = group_neigh_local(i1,2)
            idx1 = neigh_idx_local(i1,1)
            idx2 = neigh_idx_local(i1,2)
            if(grouprij(i1,1).le.cutoff.and.grouprij(i1,2).le.cutoff)then

                ! r_ij
                deltaxj = -1.d0 * (xyzstruct_local(1,jcount)-lstb(ineigh1,1))
                deltayj = -1.d0 * (xyzstruct_local(2,jcount)-lstb(ineigh1,2))
                deltazj = -1.d0 * (xyzstruct_local(3,jcount)-lstb(ineigh1,3))
                drijdxi = -deltaxj/grouprij(i1,1)
                drijdyi = -deltayj/grouprij(i1,1)
                drijdzi = -deltazj/grouprij(i1,1)
                drijdxj =- 1.d0*drijdxi
                drijdyj = -1.d0*drijdyi
                drijdzj = -1.d0*drijdzi
                drijdxk=0.0d0
                drijdyk=0.0d0
                drijdzk=0.0d0

                dfcutijdxi = dfcut_local(i1,1)*drijdxi
                dfcutijdyi = dfcut_local(i1,1)*drijdyi
                dfcutijdzi = dfcut_local(i1,1)*drijdzi
                dfcutijdxj = -1.d0*dfcutijdxi
                dfcutijdyj = -1.d0*dfcutijdyi
                dfcutijdzj = -1.d0*dfcutijdzi
                dfcutijdxk = 0.0d0
                dfcutijdyk = 0.0d0
                dfcutijdzk = 0.0d0

                ! r_ik
                deltaxk = -1.d0 * (xyzstruct_local(1,jcount)-lstb(ineigh2,1))
                deltayk = -1.d0 * (xyzstruct_local(2,jcount)-lstb(ineigh2,2))
                deltazk = -1.d0 * (xyzstruct_local(3,jcount)-lstb(ineigh2,3))
                drikdxi = -deltaxk/grouprij(i1,2)
                drikdyi = -deltayk/grouprij(i1,2)
                drikdzi = -deltazk/grouprij(i1,2)
                drikdxk =- 1.d0*drikdxi
                drikdyk = -1.d0*drikdyi
                drikdzk = -1.d0*drikdzi
                drikdxj = 0.0d0
                drikdyj = 0.0d0
                drikdzj = 0.0d0

                dfcutikdxi = dfcut_local(i1,2)*drikdxi
                dfcutikdyi = dfcut_local(i1,2)*drikdyi
                dfcutikdzi = dfcut_local(i1,2)*drikdzi
                dfcutikdxk = -1.d0*dfcutikdxi
                dfcutikdyk = -1.d0*dfcutikdyi
                dfcutikdzk = -1.d0*dfcutikdzi
                dfcutikdxj = 0.0d0
                dfcutikdyj = 0.0d0
                dfcutikdzj = 0.0d0

                ! r_jk
                drjkdxj = 1.d0*(lstb(ineigh1,1)-lstb(ineigh2,1))/grouprij(i1,3)
                drjkdyj = 1.d0*(lstb(ineigh1,2)-lstb(ineigh2,2))/grouprij(i1,3)
                drjkdzj = 1.d0*(lstb(ineigh1,3)-lstb(ineigh2,3))/grouprij(i1,3)
                drjkdxk =- 1.d0*drjkdxj
                drjkdyk = -1.d0*drjkdyj
                drjkdzk = -1.d0*drjkdzj
                drjkdxi = 0.0d0
                drjkdyi = 0.0d0
                drjkdzi = 0.0d0

                ! Cosine
                f = (lambda_local(isf,isfg,eindex)*(grouprij(i1,3)**2 - grouprij(i1,1)**2 &
                        - grouprij(i1,2)**2))
                g = (-2.d0*grouprij(i1,1)*grouprij(i1,2))
                costheta = f/g
                costheta = 1.d0+costheta ! avoid negative values
                temp4   = 2.d0**(1.d0-zeta_local(isf,isfg,eindex))
                !! avoid negative values due to numerical noise
                if(costheta.le.0.d0) then
                    costheta = 0.d0
                else
                    costheta = temp4*(costheta**zeta_local(isf,isfg,eindex))
                endif

                ! Cosine derivative
                dfdxi = lambda_local(isf,isfg,eindex)*(-2.d0*grouprij(i1,1)*drijdxi - &
                        2.d0*grouprij(i1,2)*drikdxi)
                dfdyi = lambda_local(isf,isfg,eindex)*(-2.d0*grouprij(i1,1)*drijdyi - &
                        2.d0*grouprij(i1,2)*drikdyi)
                dfdzi = lambda_local(isf,isfg,eindex)*(-2.d0*grouprij(i1,1)*drijdzi - &
                        2.d0*grouprij(i1,2)*drikdzi)
                !!
                dfdxj = lambda_local(isf,isfg,eindex)*(2.d0*grouprij(i1,3)*drjkdxj - &
                        2.d0*grouprij(i1,1)*drijdxj)
                dfdyj = lambda_local(isf,isfg,eindex)*(2.d0*grouprij(i1,3)*drjkdyj - &
                        2.d0*grouprij(i1,1)*drijdyj)
                dfdzj = lambda_local(isf,isfg,eindex)*(2.d0*grouprij(i1,3)*drjkdzj - &
                        2.d0*grouprij(i1,1)*drijdzj)
                !!
                dfdxk = lambda_local(isf,isfg,eindex)*(2.d0*grouprij(i1,3)*drjkdxk - &
                        2.d0*grouprij(i1,2)*drikdxk)
                dfdyk = lambda_local(isf,isfg,eindex)*(2.d0*grouprij(i1,3)*drjkdyk - &
                        2.d0*grouprij(i1,2)*drikdyk)
                dfdzk = lambda_local(isf,isfg,eindex)*(2.d0*grouprij(i1,3)*drjkdzk - &
                        2.d0*grouprij(i1,2)*drikdzk)

                !!
                dgdxi = -2.d0*(drijdxi*grouprij(i1,2) + grouprij(i1,1)*drikdxi)
                dgdyi = -2.d0*(drijdyi*grouprij(i1,2) + grouprij(i1,1)*drikdyi)
                dgdzi = -2.d0*(drijdzi*grouprij(i1,2) + grouprij(i1,1)*drikdzi)
                !!
                dgdxj = -2.d0*drijdxj*grouprij(i1,2)
                dgdyj = -2.d0*drijdyj*grouprij(i1,2)
                dgdzj = -2.d0*drijdzj*grouprij(i1,2)
                !!
                dgdxk = -2.d0*grouprij(i1,1)*drikdxk
                dgdyk = -2.d0*grouprij(i1,1)*drikdyk
                dgdzk = -2.d0*grouprij(i1,1)*drikdzk

                temp5 = 1 + f/g
                if(temp5.lt.0.0d0)costheta=0.0d0
                temp4 = temp4*zeta_local(isf,isfg,eindex)*temp5**(zeta_local(isf,isfg,eindex)-1)/g**2
                dcosthetadxi=temp4*(dfdxi*g - f*dgdxi) !/g**2
                dcosthetadyi=temp4*(dfdyi*g - f*dgdyi) !/g**2
                dcosthetadzi=temp4*(dfdzi*g - f*dgdzi) !/g**2
                dcosthetadxj=temp4*(dfdxj*g - f*dgdxj) !/g**2
                dcosthetadyj=temp4*(dfdyj*g - f*dgdyj) !/g**2
                dcosthetadzj=temp4*(dfdzj*g - f*dgdzj) !/g**2
                dcosthetadxk=temp4*(dfdxk*g - f*dgdxk) !/g**2
                dcosthetadyk=temp4*(dfdyk*g - f*dgdyk) !/g**2
                dcosthetadzk=temp4*(dfdzk*g - f*dgdzk) !/g**2

                ! Exponentials
                expxyz = dexp(-eta_local(isf,isfg,eindex)*(grouprij(i1,1)**2+grouprij(i1,2)**2))
                temp3 = -eta_local(isf,isfg,eindex)*2.0d0*expxyz
                dexpxyzdxi=(grouprij(i1,1)*drijdxi+grouprij(i1,2)*drikdxi)*temp3
                dexpxyzdyi=(grouprij(i1,1)*drijdyi+grouprij(i1,2)*drikdyi)*temp3
                dexpxyzdzi=(grouprij(i1,1)*drijdzi+grouprij(i1,2)*drikdzi)*temp3
                dexpxyzdxj=(grouprij(i1,1)*drijdxj+grouprij(i1,2)*drikdxj)*temp3
                dexpxyzdyj=(grouprij(i1,1)*drijdyj+grouprij(i1,2)*drikdyj)*temp3
                dexpxyzdzj=(grouprij(i1,1)*drijdzj+grouprij(i1,2)*drikdzj)*temp3
                dexpxyzdxk=(grouprij(i1,1)*drijdxk+grouprij(i1,2)*drikdxk)*temp3
                dexpxyzdyk=(grouprij(i1,1)*drijdyk+grouprij(i1,2)*drikdyk)*temp3
                dexpxyzdzk=(grouprij(i1,1)*drijdzk+grouprij(i1,2)*drikdzk)*temp3

                ! dsfunc/dx
                temp1=(+dcosthetadxi* expxyz   * fcut_local(i1,1)   * fcut_local(i1,2) &
                        +costheta   *dexpxyzdxi* fcut_local(i1,1)   * fcut_local(i1,2)  &
                        +costheta   * expxyz   *    dfcutijdxi          * fcut_local(i1,2)   &
                        +costheta   * expxyz   * fcut_local(i1,1)   *dfcutikdxi)
                temp2=(+dcosthetadxj* expxyz   * fcut_local(i1,1)   * fcut_local(i1,2)  &
                        +costheta   *dexpxyzdxj* fcut_local(i1,1)   * fcut_local(i1,2)  &
                        +costheta   * expxyz   *dfcutijdxj* fcut_local(i1,2)   &
                        +costheta   * expxyz   * fcut_local(i1,1)   *dfcutikdxj)
                temp3=(+dcosthetadxk* expxyz   * fcut_local(i1,1)   * fcut_local(i1,2) &
                        +costheta   *dexpxyzdxk* fcut_local(i1,1)   * fcut_local(i1,2) &
                        +costheta   * expxyz   *dfcutijdxk* fcut_local(i1,2)   &
                        +costheta   * expxyz   * fcut_local(i1,1)   *dfcutikdxk)
                dsfuncdxyz_temp(invneighboridx(iatom,jcount),1)=dsfuncdxyz_temp(invneighboridx(iatom,jcount),1)+&
                        temp1
                dsfuncdxyz_temp(invneighboridx(iatom,idx1),1) =dsfuncdxyz_temp(invneighboridx(iatom,idx1),1)+&
                        temp2
                dsfuncdxyz_temp(invneighboridx(iatom,idx2),1) =dsfuncdxyz_temp(invneighboridx(iatom,idx2),1)+&
                        temp3

                if(ldostress)then
                    strs_temp(1,1,isf)=strs_temp(1,1,isf)+deltaxj*temp2&
                            +deltaxk*temp3
                    strs_temp(2,1,isf)=strs_temp(2,1,isf)+deltayj*temp2&
                            +deltayk*temp3
                    strs_temp(3,1,isf)=strs_temp(3,1,isf)+deltazj*temp2&
                            +deltazk*temp3
                end if

                ! dsfunc/dy
                temp1=(+dcosthetadyi* expxyz   * fcut_local(i1,1)  * fcut_local(i1,2)  &
                        +costheta   *dexpxyzdyi* fcut_local(i1,1)   * fcut_local(i1,2) &
                        +costheta   * expxyz   *dfcutijdyi* fcut_local(i1,2)   &
                        +costheta   * expxyz   * fcut_local(i1,1)   *dfcutikdyi)
                temp2=(+dcosthetadyj* expxyz   * fcut_local(i1,1)   * fcut_local(i1,2) &
                        +costheta   *dexpxyzdyj* fcut_local(i1,1)   * fcut_local(i1,2) &
                        +costheta   * expxyz   *dfcutijdyj* fcut_local(i1,2)&
                        +costheta   * expxyz   * fcut_local(i1,1)   *dfcutikdyj)
                temp3=(+dcosthetadyk* expxyz   * fcut_local(i1,1)   * fcut_local(i1,2)  &
                        +costheta   *dexpxyzdyk* fcut_local(i1,1)   * fcut_local(i1,2)  &
                        +costheta   * expxyz   *dfcutijdyk* fcut_local(i1,2)  &
                        +costheta   * expxyz   * fcut_local(i1,1)   *dfcutikdyk)
                dsfuncdxyz_temp(invneighboridx(iatom,jcount),2)=dsfuncdxyz_temp(invneighboridx(iatom,jcount),2)+&
                        temp1
                dsfuncdxyz_temp(invneighboridx(iatom,idx1),2) =dsfuncdxyz_temp(invneighboridx(iatom,idx1),2)+&
                        temp2
                dsfuncdxyz_temp(invneighboridx(iatom,idx2),2) =dsfuncdxyz_temp(invneighboridx(iatom,idx2),2)+&
                        temp3

                if(ldostress)then
                    strs_temp(1,2,isf)=strs_temp(1,2,isf)+deltaxj*temp2&
                            +deltaxk*temp3
                    strs_temp(2,2,isf)=strs_temp(2,2,isf)+deltayj*temp2&
                            +deltayk*temp3
                    strs_temp(3,2,isf)=strs_temp(3,2,isf)+deltazj*temp2&
                            +deltazk*temp3
                end if

                ! dsfunc/dz
                temp1=(+dcosthetadzi* expxyz   * fcut_local(i1,1)   * fcut_local(i1,2) &
                        + costheta   *dexpxyzdzi* fcut_local(i1,1)   * fcut_local(i1,2) &
                        + costheta   * expxyz   *dfcutijdzi* fcut_local(i1,2)  &
                        + costheta   * expxyz   * fcut_local(i1,1)   *dfcutikdzi)
                temp2=(+dcosthetadzj* expxyz   * fcut_local(i1,1)   * fcut_local(i1,2) &
                        + costheta   *dexpxyzdzj* fcut_local(i1,1)   * fcut_local(i1,2)&
                        + costheta   * expxyz   *dfcutijdzj* fcut_local(i1,2) &
                        + costheta   * expxyz   * fcut_local(i1,1)   *dfcutikdzj)
                temp3=(+dcosthetadzk* expxyz   * fcut_local(i1,1)   * fcut_local(i1,2)  &
                        + costheta   *dexpxyzdzk* fcut_local(i1,1)   * fcut_local(i1,2)  &
                        + costheta   * expxyz   *dfcutijdzk* fcut_local(i1,2)  &
                        + costheta   * expxyz   * fcut_local(i1,1)   *dfcutikdzk)
                dsfuncdxyz_temp(invneighboridx(iatom,jcount),3)=dsfuncdxyz_temp(invneighboridx(iatom,jcount),3)+&
                        temp1
                dsfuncdxyz_temp(invneighboridx(iatom,idx1),3) =dsfuncdxyz_temp(invneighboridx(iatom,idx1),3)+&
                        temp2
                dsfuncdxyz_temp(invneighboridx(iatom,idx2),3) =dsfuncdxyz_temp(invneighboridx(iatom,idx2),3)+&
                        temp3

                if(ldostress)then
                    strs_temp(1,3,isf)=strs_temp(1,3,isf)+deltaxj*temp2&
                            +deltaxk*temp3
                    strs_temp(2,3,isf)=strs_temp(2,3,isf)+deltayj*temp2&
                            +deltayk*temp3
                    strs_temp(3,3,isf)=strs_temp(3,3,isf)+deltazj*temp2&
                            +deltazk*temp3
                end if

                end if

        end do

    end if

end subroutine getsymfunctionderivatives
