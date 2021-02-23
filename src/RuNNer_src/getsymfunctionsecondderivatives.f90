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
!! Purpose : calculating second derivatives of the SFs for Hessian matrix
!! called by:
!! - getshortatomic.f90
!! - added by Emir

subroutine getsymfunctionsecondderivatives(iatom,isfg,isf,eindex,jcount,xyzstruct_local,&
        invneighboridx,listdim,natoms,maxnum_funcvalues_local,maxnum_sfgroups_local,max_num_neigh_local,&
        max_num_atoms,group_neigh_local,neigh_idx_local,groupdim,maxdim,grouprij,fcut_local,dfcut_local,ddfcut_local,&
        grouptype,nelem,eta_local,rshift_local,lambda_local,zeta_local,ftype,cutoff,lstb,ddsfuncddxyz_temp)

    use nnconstants
    use fileunits

    implicit none

    integer i1,i2,i3,i4,i5 ! internal
    integer ineigh1,ineigh2,idx1,idx2 ! internal
    integer groupdim
    integer maxdim
    integer iatom,isf,isfg,eindex,jcount
    integer ftype
    integer nelem
    integer maxnum_funcvalues_local,maxnum_sfgroups_local
    integer grouptype
    integer max_num_atoms
    integer max_num_neigh_local
    integer listdim
    integer natoms
    integer pairidx(3)
    integer group_neigh_local(maxdim,grouptype+1)
    integer neigh_idx_local(maxdim,grouptype+1)
    integer invneighboridx(natoms,max_num_atoms) ! in

    real*8  grouprij(maxdim,3)
    real*8  fcut_local(maxdim,3)
    real*8  dfcut_local(maxdim,3)
    real*8  ddfcut_local(maxdim,3)
    real*8  xyzstruct_local(3,max_num_atoms) ! in
    real*8  deltaij(9),deltaik(9),deltajk(9)
    real*8  ddeltaij(9,9),ddeltaik(9,9),ddeltajk(9,9)
    real*8  drijd(9),drikd(9),drjkd(9)
    real*8  ddrijdd(9,9),ddrikdd(9,9),ddrjkdd(9,9)
    real*8  dfcutijd(9),dfcutikd(9),dfcutjkd(9)
    real*8  f,g
    real*8  dfdxyz(9),dgdxyz(9),ddfdxyz(9,9),ddgdxyz(9,9)
    real*8  ddfcutijdd(9,9),ddfcutikdd(9,9),ddfcutjkdd(9,9)
    real*8  costheta,theta,expxyz
    real*8  dcostheta(9),dtheta(9),dexpxyz(9)
    real*8  ddcostheta(9,9),ddtheta(9,9),ddexpxyz(9,9)
    real*8  lstb(listdim,4) ! in
    real*8  cutoff
    real*8  temp1,temp2,temp3,temp4,temp5
    real*8  eta_local(maxnum_funcvalues_local,maxnum_sfgroups_local,nelem) ! in
    real*8  zeta_local(maxnum_funcvalues_local,maxnum_sfgroups_local,nelem) ! in
    real*8  rshift_local(maxnum_funcvalues_local,maxnum_sfgroups_local,nelem) ! in
    real*8  lambda_local(maxnum_funcvalues_local,maxnum_sfgroups_local,nelem) ! in
    real*8  ddsfuncddxyz_temp(0:max_num_neigh_local,3,0:max_num_neigh_local,3) ! OUT

    if(ftype.eq.1)then ! only type 1 cutoff function

        do i1 = 1,groupdim
            ineigh1 = group_neigh_local(i1,1)
            idx1 = neigh_idx_local(i1,1)
            pairidx(1) = invneighboridx(iatom,jcount)
            pairidx(2) = invneighboridx(iatom,idx1)
            deltaij(:) = 0.d0

            if(grouprij(i1,1).le.cutoff)then

                deltaij(1) = (xyzstruct_local(1,jcount)-lstb(ineigh1,1))
                deltaij(2) = (xyzstruct_local(2,jcount)-lstb(ineigh1,2))
                deltaij(3) = (xyzstruct_local(3,jcount)-lstb(ineigh1,3))
                deltaij(4) = -1.d0*deltaij(1)
                deltaij(5) = -1.d0*deltaij(2)
                deltaij(6) = -1.d0*deltaij(3)

                do i2 = 1,6
                    drijd(i2) = deltaij(i2) / grouprij(i1,1)
                    dfcutijd(i2) = dfcut_local(i1,1) * drijd(i2)
                end do

                !!! d^2R/d^2 and d^2fc/d^2 terms
                ddeltaij(:,:)=0.d0 !CHECK ME

                do i2 = 1,6
                    ddeltaij(i2,i2) = 1.d0
                    if(i2.le.3)then
                        ddeltaij(i2,i2+3) = -1.d0
                    else
                        ddeltaij(i2,i2-3) = -1.d0
                    end if
                end do

                ddrijdd(:,:) = 0.d0
                ddfcutijdd(:,:) = 0.d0

                do i2 = 1,6
                    do i3 = 1,6
                        ddrijdd(i2,i3) = (ddeltaij(i2,i3)*grouprij(i1,1)-deltaij(i2)*drijd(i3))/&
                                (grouprij(i1,1)**2.d0)
                        ddfcutijdd(i2,i3) = ddfcut_local(i1,1)*drijd(i2)*drijd(i3) + dfcut_local(i1,1)*ddrijdd(i2,i3)
                    end do
                end do

                !! Second derivatives
                do i2 = 1,2
                    do i3 = 1,3
                        do i4 = 1,2
                            do i5 = 1,3
                                ddsfuncddxyz_temp(pairidx(i2),i3,pairidx(i4),i5) = &
                                ddsfuncddxyz_temp(pairidx(i2),i3,pairidx(i4),i5) + ddfcutijdd((i2-1)*3+i3,(i4-1)*3+i5)
                            end do
                        end do
                    end do
                end do

            end if
        end do

    elseif(ftype.eq.2)then

        do i1 = 1,groupdim
            ineigh1 = group_neigh_local(i1,1)
            idx1 = neigh_idx_local(i1,1)
            pairidx(1) = invneighboridx(iatom,jcount)
            pairidx(2) = invneighboridx(iatom,idx1)
            deltaij(:) = 0.d0

            if(grouprij(i1,1).le.cutoff)then
                deltaij(1) = (xyzstruct_local(1,jcount)-lstb(ineigh1,1))
                deltaij(2) = (xyzstruct_local(2,jcount)-lstb(ineigh1,2))
                deltaij(3) = (xyzstruct_local(3,jcount)-lstb(ineigh1,3))
                deltaij(4) = -1.d0*deltaij(1)
                deltaij(5) = -1.d0*deltaij(2)
                deltaij(6) = -1.d0*deltaij(3)

                do i2 = 1,6
                    drijd(i2) = deltaij(i2) / grouprij(i1,1)
                    dfcutijd(i2) = dfcut_local(i1,1) * drijd(i2)
                end do

                !!! d^2R/d^2 and d^2fc/d^2 terms
                ddeltaij(:,:)=0.d0 !CHECK ME

                do i2 = 1,6
                    ddeltaij(i2,i2) = 1.d0
                    if(i2.le.3)then
                        ddeltaij(i2,i2+3) = -1.d0
                    else
                        ddeltaij(i2,i2-3) = -1.d0
                    end if
                end do

                ddrijdd(:,:) = 0.d0
                ddfcutijdd(:,:) = 0.d0

                do i2 = 1,6
                    do i3 = 1,6
                        ddrijdd(i2,i3) = (ddeltaij(i2,i3)*grouprij(i1,1)-deltaij(i2)*drijd(i3))/&
                                (grouprij(i1,1)**2.d0)
                        ddfcutijdd(i2,i3) = ddfcut_local(i1,1)*drijd(i2)*drijd(i3) + dfcut_local(i1,1)*ddrijdd(i2,i3)
                    end do
                end do

                temp1 = -2.d0 * eta_local(isf,isfg,eindex) * (grouprij(i1,1)-rshift_local(isf,isfg,eindex))
                temp2 = dexp(-1.d0*eta_local(isf,isfg,eindex)*((grouprij(i1,1)-rshift_local(isf,isfg,eindex))**2))

                !! Second derivatives
                do i2 = 1,2
                    do i3 = 1,3
                        do i4 = 1,2
                            do i5 = 1,3
                                ddsfuncddxyz_temp(pairidx(i2),i3,pairidx(i4),i5) = &
                                    ddsfuncddxyz_temp(pairidx(i2),i3,pairidx(i4),i5) + &
                                    (-2.d0*eta_local(isf,isfg,eindex)*drijd((i2-1)*3+i3)*drijd((i4-1)*3+i5)*temp2*&
                                    fcut_local(i1,1)) + (dfcutijd((i4-1)*3+i5) * drijd((i2-1)*3+i3) * temp1 * temp2) + &
                                    (fcut_local(i1,1) * drijd((i2-1)*3+i3) * drijd((i4-1)*3+i5) * temp2 * temp1 * temp1) + &
                                    (fcut_local(i1,1) * ddrijdd((i2-1)*3+i3,(i4-1)*3+i5) * temp2 * temp1) + &
                                    (drijd((i4-1)*3+i5) * dfcutijd((i2-1)*3+i3) * temp2 * temp1) + &
                                    (ddfcutijdd((i2-1)*3+i3,(i4-1)*3+i5) * temp2)
                            end do
                        end do
                    end do
                end do


            end if
        end do


    elseif(ftype.eq.3)then

        do i1 = 1,groupdim
            ineigh1 = group_neigh_local(i1,1)
            ineigh2 = group_neigh_local(i1,2)
            idx1 = neigh_idx_local(i1,1)
            idx2 = neigh_idx_local(i1,2)
            pairidx(1) = invneighboridx(iatom,jcount)
            pairidx(2) = invneighboridx(iatom,idx1)
            pairidx(3) = invneighboridx(iatom,idx2)
            deltaij(:) = 0.d0
            deltaik(:) = 0.d0
            deltajk(:) = 0.d0

            if(grouprij(i1,1).le.cutoff.and.grouprij(i1,2).le.cutoff.and.grouprij(i1,3).le.cutoff)then
                ! r_ij
                deltaij(1) = (xyzstruct_local(1,jcount)-lstb(ineigh1,1)) !(xi-xj)
                deltaij(2) = (xyzstruct_local(2,jcount)-lstb(ineigh1,2))
                deltaij(3) = (xyzstruct_local(3,jcount)-lstb(ineigh1,3))
                deltaij(4) = -1.d0*deltaij(1) !(xj-xi)
                deltaij(5) = -1.d0*deltaij(2)
                deltaij(6) = -1.d0*deltaij(3)
                do i2 = 1,6
                    drijd(i2) = deltaij(i2) / grouprij(i1,1)
                    dfcutijd(i2) = dfcut_local(i1,1) * drijd(i2)
                end do

                ! r_ik
                deltaik(1) = (xyzstruct_local(1,jcount)-lstb(ineigh2,1)) !(xi-xk)
                deltaik(2) = (xyzstruct_local(2,jcount)-lstb(ineigh2,2))
                deltaik(3) = (xyzstruct_local(3,jcount)-lstb(ineigh2,3))
                deltaik(7) = -1.d0*deltaik(1) !(xk-xi)
                deltaik(8) = -1.d0*deltaik(2)
                deltaik(9) = -1.d0*deltaik(3)
                drikd(:) = 0.d0
                dfcutikd(:) = 0.d0
                do i2 = 1,9
                    drikd(i2) = deltaik(i2) / grouprij(i1,2)
                    dfcutikd(i2) = dfcut_local(i1,2) * drikd(i2)
                end do

                ! r_jk
                deltajk(4) = (lstb(ineigh1,1)-lstb(ineigh2,1)) !(xj-xk)
                deltajk(5) = (lstb(ineigh1,2)-lstb(ineigh2,2))
                deltajk(6) = (lstb(ineigh1,3)-lstb(ineigh2,3))
                deltajk(7) = -1.d0*deltajk(4) !(xk-xj)
                deltajk(8) = -1.d0*deltajk(5)
                deltajk(9) = -1.d0*deltajk(6)
                drjkd(:) = 0.d0
                dfcutjkd(:) = 0.d0
                do i2 = 4,9
                    drjkd(i2) = deltajk(i2) / grouprij(i1,3)
                    dfcutjkd(i2) = dfcut_local(i1,3) * drjkd(i2)
                end do

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

                !!! d^2R/d^2 and d^2fc/d^2 terms
                ddeltaij(:,:)=0.d0 !CHECK ME
                ddeltaik(:,:)=0.d0
                ddeltajk(:,:)=0.d0

                do i2 = 1,9
                    if(i2.le.3)then
                        ddeltaij(i2,i2+3) = -1.d0
                        ddeltaij(i2,i2) = 1.d0
                        ddeltaik(i2,i2+6) = -1.d0
                        ddeltaik(i2,i2) = 1.d0
                    else if (i2.gt.3.and.i2.le.6)then
                        ddeltaij(i2,i2-3) = -1.d0
                        ddeltajk(i2,i2+3) = -1.d0
                        ddeltaij(i2,i2) = 1.d0
                        ddeltajk(i2,i2) = 1.d0
                    else
                        ddeltaik(i2,i2-6) = -1.d0
                        ddeltaik(i2,i2) = 1.d0
                        ddeltajk(i2,i2-3) = -1.d0
                        ddeltajk(i2,i2) = 1.d0
                    end if
                end do

                do i2 = 1,9
                    do i3 = 1,9
                        ddrijdd(i2,i3) = (ddeltaij(i2,i3)*grouprij(i1,1)-deltaij(i2)*drijd(i3))/&
                                (grouprij(i1,1)**2.d0)
                        ddrikdd(i2,i3) = (ddeltaik(i2,i3)*grouprij(i1,2)-deltaik(i2)*drikd(i3))/&
                                (grouprij(i1,2)**2.d0)
                        ddrjkdd(i2,i3) = (ddeltajk(i2,i3)*grouprij(i1,3)-deltajk(i2)*drjkd(i3))/&
                                (grouprij(i1,3)**2.d0)
                        ddfcutijdd(i2,i3) = ddfcut_local(i1,1)*drijd(i2)*drijd(i3) + dfcut_local(i1,1)*ddrijdd(i2,i3)
                        ddfcutikdd(i2,i3) = ddfcut_local(i1,2)*drikd(i2)*drikd(i3) + dfcut_local(i1,2)*ddrikdd(i2,i3)
                        ddfcutjkdd(i2,i3) = ddfcut_local(i1,3)*drjkd(i2)*drjkd(i3) + dfcut_local(i1,3)*ddrjkdd(i2,i3)
                    end do
                end do

                do i2 = 1,3
                    dfdxyz(i2) = lambda_local(isf,isfg,eindex)*(-2.d0*grouprij(i1,1)*drijd(i2) - &
                            2.d0*grouprij(i1,2)*drikd(i2))
                    dfdxyz(i2+3) = lambda_local(isf,isfg,eindex)*(2.d0*grouprij(i1,3)*drjkd(i2+3) - &
                            2.d0*grouprij(i1,1)*drijd(i2+3))
                    dfdxyz(i2+6) = lambda_local(isf,isfg,eindex)*(2.d0*grouprij(i1,3)*drjkd(i2+6) - &
                            2.d0*grouprij(i1,2)*drikd(i2+6))
                    do i3 = 1,9 !check me
                        ddfdxyz(i2,i3) = -2.d0*lambda_local(isf,isfg,eindex)*(drijd(i2)*drijd(i3) + grouprij(i1,1)*&
                                ddrijdd(i2,i3) + drikd(i2)*drikd(i3) + grouprij(i1,2)*ddrikdd(i2,i3))
                        ddfdxyz(i2+3,i3) = 2.d0*lambda_local(isf,isfg,eindex)*(drjkd(i2+3)*drjkd(i3) + grouprij(i1,3)*&
                                ddrjkdd(i2+3,i3) - drijd(i2+3)*drijd(i3) - grouprij(i1,1)*ddrijdd(i2+3,i3))
                        ddfdxyz(i2+6,i3) = 2.d0*lambda_local(isf,isfg,eindex)*(drjkd(i2+6)*drjkd(i3) + grouprij(i1,3)*&
                                ddrjkdd(i2+6,i3) - drikd(i2+6)*drikd(i3) - grouprij(i1,2)*ddrikdd(i2+6,i3))
                    end do
                end do

                do i2 = 1,9
                    dgdxyz(i2) = -2.d0*(drijd(i2)*grouprij(i1,2) + grouprij(i1,1)*drikd(i2))
                    do i3 = 1,9
                        ddgdxyz(i2,i3) = -2.d0*(drikd(i3)*drijd(i2)+grouprij(i1,2)*ddrijdd(i2,i3)+drijd(i3)*drikd(i2)+&
                                grouprij(i1,1)*ddrikdd(i2,i3))
                    end do
                end do

                temp5 = 1 + f/g
                if(temp5.lt.0.0d0)costheta=0.0d0
                temp4 = temp4*zeta_local(isf,isfg,eindex)*temp5**(zeta_local(isf,isfg,eindex)-1)/g**2
                do i2 = 1,9
                    dcostheta(i2) = temp4*(dfdxyz(i2)*g - f*dgdxyz(i2))
                end do

                ! Exponentials
                expxyz = dexp(-eta_local(isf,isfg,eindex)*(grouprij(i1,1)**2+grouprij(i1,2)**2&
                        +grouprij(i1,3)**2))
                temp3 = -eta_local(isf,isfg,eindex)*2.0d0*expxyz
                do i2 = 1,9
                    dexpxyz(i2) = (grouprij(i1,1)*drijd(i2) + grouprij(i1,2)*drikd(i2) + grouprij(i1,3)*drjkd(i2))*temp3
                end do

                do i2 = 1,9
                    do i3 = 1,9
                        ddcostheta(i2,i3) = (2.d0**(1.d0-zeta_local(isf,isfg,eindex)))*zeta_local(isf,isfg,eindex) * &
                                ((zeta_local(isf,isfg,eindex)-1) * (1+f/g)**(zeta_local(isf,isfg,eindex)-2) * &
                                (dfdxyz(i2)*g - f*dgdxyz(i2)) * (dfdxyz(i3)*g - f*dgdxyz(i3)) + &
                                (1+f/g)**(zeta_local(isf,isfg,eindex)-1) * ((ddfdxyz(i2,i3)*g + dfdxyz(i2)*dgdxyz(i3) -&
                                dgdxyz(i2)*dfdxyz(i3) - ddgdxyz(i2,i3)*f)*g*g - (dfdxyz(i2)*g - dgdxyz(i2)*f)*2.d0*g*&
                                dgdxyz(i3))) / g**4.d0

                        ddexpxyz(i2,i3) = -2.d0*eta_local(isf,isfg,eindex)*(grouprij(i1,1)*drijd(i3) + grouprij(i1,2)*&
                                drikd(i3) + grouprij(i1,3)*drjkd(i3))*dexpxyz(i2) + (-2.d0*expxyz*eta_local(isf,isfg,&
                                eindex)*(drijd(i3)*drijd(i2)+grouprij(i1,1)*ddrijdd(i2,i3) + drikd(i3)*drikd(i2) + &
                                grouprij(i1,2)*ddrikdd(i2,i3) + drjkd(i3)*drjkd(i2) + grouprij(i1,3)*ddrjkdd(i2,i3)))
                    end do
                end do


                !! Second derivatives
                do i2 = 1,3 ! atom,n1,n2
                    do i3 = 1,3 ! xyz
                        do i4 = 1,3 ! atom,n1,n2
                            do i5 = 1,3 ! xyz
                            ddsfuncddxyz_temp(pairidx(i2),i3,pairidx(i4),i5) = &
                            ddsfuncddxyz_temp(pairidx(i2),i3,pairidx(i4),i5) + ddcostheta((i2-1)*3+i3,&
                            (i4-1)*3+i5)*expxyz*fcut_local(i1,1)*fcut_local(i1,2)*fcut_local(i1,3) + &
                            dcostheta((i2-1)*3+i3)*dexpxyz((i4-1)*3+i5)*fcut_local(i1,1)*fcut_local(i1,2)*&
                            fcut_local(i1,3) + dcostheta((i2-1)*3+i3)*expxyz*dfcutijd((i4-1)*3+i5)*fcut_local(i1,2)*&
                            fcut_local(i1,3) + dcostheta((i2-1)*3+i3)*expxyz*fcut_local(i1,1)*dfcutikd((i4-1)*3+i5)*&
                            fcut_local(i1,3) + dcostheta((i2-1)*3+i3)*expxyz*fcut_local(i1,1)*fcut_local(i1,2)*&
                            dfcutjkd((i4-1)*3+i5) + &

                            dcostheta((i4-1)*3+i5)*dexpxyz((i2-1)*3+i3)*fcut_local(i1,1)*fcut_local(i1,2)*&
                            fcut_local(i1,3) + ddexpxyz((i2-1)*3+i3,(i4-1)*3+i5)*costheta*fcut_local(i1,1)*&
                            fcut_local(i1,2)*fcut_local(i1,3)+ dexpxyz((i2-1)*3+i3)*costheta*dfcutijd((i4-1)*3+i5)*&
                            fcut_local(i1,2)*fcut_local(i1,3)+ dexpxyz((i2-1)*3+i3)*costheta*fcut_local(i1,1)*&
                            dfcutikd((i4-1)*3+i5)*fcut_local(i1,3)+ dexpxyz((i2-1)*3+i3)*costheta*fcut_local(i1,1)*&
                            fcut_local(i1,2)*dfcutjkd((i4-1)*3+i5) + &

                            dcostheta((i4-1)*3+i5)*expxyz*dfcutijd((i2-1)*3+i3)*fcut_local(i1,2)*&
                            fcut_local(i1,3) + dexpxyz((i4-1)*3+i5)*costheta*dfcutijd((i2-1)*3+i3)*&
                            fcut_local(i1,2)*fcut_local(i1,3)+ expxyz*costheta*ddfcutijdd((i2-1)*3+i3,(i4-1)*3+i5)*&
                            fcut_local(i1,2)*fcut_local(i1,3)+ expxyz*costheta*dfcutijd((i2-1)*3+i3)*&
                            dfcutikd((i4-1)*3+i5)*fcut_local(i1,3)+ expxyz*costheta*dfcutijd((i2-1)*3+i3)*&
                            fcut_local(i1,2)*dfcutjkd((i4-1)*3+i5) + &

                            dcostheta((i4-1)*3+i5)*expxyz*dfcutikd((i2-1)*3+i3)*fcut_local(i1,1)*&
                            fcut_local(i1,3) + dexpxyz((i4-1)*3+i5)*costheta*dfcutikd((i2-1)*3+i3)*&
                            fcut_local(i1,1)*fcut_local(i1,3)+ expxyz*costheta*dfcutijd((i4-1)*3+i5)*&
                            dfcutikd((i2-1)*3+i3)*fcut_local(i1,3)+ expxyz*costheta*fcut_local(i1,1)*&
                            ddfcutikdd((i2-1)*3+i3,(i4-1)*3+i5)*fcut_local(i1,3)+ expxyz*costheta*fcut_local(i1,1)*&
                            dfcutikd((i2-1)*3+i3)*dfcutjkd((i4-1)*3+i5) + &

                            dcostheta((i4-1)*3+i5)*expxyz*dfcutjkd((i2-1)*3+i3)*fcut_local(i1,1)*&
                            fcut_local(i1,2) + dexpxyz((i4-1)*3+i5)*costheta*dfcutjkd((i2-1)*3+i3)*&
                            fcut_local(i1,1)*fcut_local(i1,2)+ expxyz*costheta*dfcutijd((i4-1)*3+i5)*&
                            dfcutjkd((i2-1)*3+i3)*fcut_local(i1,2)+ expxyz*costheta*fcut_local(i1,1)*&
                            ddfcutjkdd((i2-1)*3+i3,(i4-1)*3+i5)*fcut_local(i1,2)+ expxyz*costheta*fcut_local(i1,1)*&
                            dfcutjkd((i2-1)*3+i3)*dfcutikd((i4-1)*3+i5)
                            end do
                        end do
                    end do
                end do
            end if
        end do

    elseif(ftype.eq.4)then

        ddeltaij(:,:)=0.d0
        do i2 = 1,6
            ddeltaij(i2,i2) = 1.d0
            if(i2.le.3)then
                ddeltaij(i2,i2+3) = -1.d0
            else
                ddeltaij(i2,i2-3) = -1.d0
            end if
        end do

        do i1 = 1,groupdim
            ineigh1 = group_neigh_local(i1,1)
            idx1 = neigh_idx_local(i1,1)
            pairidx(1) = invneighboridx(iatom,jcount)
            pairidx(2) = invneighboridx(iatom,idx1)
            deltaij(:) = 0.d0

            if(grouprij(i1,1).le.cutoff)then
                deltaij(1) = (xyzstruct_local(1,jcount)-lstb(ineigh1,1))
                deltaij(2) = (xyzstruct_local(2,jcount)-lstb(ineigh1,2))
                deltaij(3) = (xyzstruct_local(3,jcount)-lstb(ineigh1,3))
                deltaij(4) = -1.d0*deltaij(1)
                deltaij(5) = -1.d0*deltaij(2)
                deltaij(6) = -1.d0*deltaij(3)

                do i2 = 1,6
                    drijd(i2) = deltaij(i2) / grouprij(i1,1)
                    dfcutijd(i2) = dfcut_local(i1,1) * drijd(i2)
                end do

                ddrijdd(:,:) = 0.d0
                ddfcutijdd(:,:) = 0.d0

                do i2 = 1,6
                    do i3 = 1,6
                        ddrijdd(i2,i3) = (ddeltaij(i2,i3)*grouprij(i1,1)-deltaij(i2)*drijd(i3))/&
                                (grouprij(i1,1)**2.d0)
                        ddfcutijdd(i2,i3) = ddfcut_local(i1,1)*drijd(i2)*drijd(i3) + dfcut_local(i1,1)*ddrijdd(i2,i3)
                    end do
                end do

                temp1 = -1.d0 * eta_local(isf,isfg,eindex) * dsin(eta_local(isf,isfg,eindex)*grouprij(i1,1))
                temp2 = dcos(eta_local(isf,isfg,eindex)*grouprij(i1,1))

                !! Second derivatives
                do i2 = 1,2
                    do i3 = 1,3
                        do i4 = 1,2
                            do i5 = 1,3
                                ddsfuncddxyz_temp(pairidx(i2),i3,pairidx(i4),i5) = &
                                ddsfuncddxyz_temp(pairidx(i2),i3,pairidx(i4),i5) + &
                                (-1.d0*eta_local(isf,isfg,eindex)*eta_local(isf,isfg,eindex)*temp2*fcut_local(i1,1)*&
                                drijd((i2-1)*3+i3)*drijd((i4-1)*3+i5) + temp1*ddrijdd((i2-1)*3+i3,(i4-1)*3+i5)*&
                                fcut_local(i1,1) + temp1*drijd((i2-1)*3+i3)*dfcutijd((i4-1)*3+i5) + &
                                temp2*ddfcutijdd((i2-1)*3+i3,(i4-1)*3+i5) + &
                                temp1*drijd((i4-1)*3+i5)*dfcutijd((i2-1)*3+i3))
                            end do
                        end do
                    end do
                end do

            end if
        end do

    elseif(ftype.eq.8)then

        ddeltaij(:,:)=0.d0
        ddeltaik(:,:)=0.d0
        ddeltajk(:,:)=0.d0

        do i2 = 1,9
            if(i2.le.3)then
                ddeltaij(i2,i2+3) = -1.d0
                ddeltaij(i2,i2) = 1.d0
                ddeltaik(i2,i2+6) = -1.d0
                ddeltaik(i2,i2) = 1.d0
            else if (i2.gt.3.and.i2.le.6)then
                ddeltaij(i2,i2-3) = -1.d0
                ddeltajk(i2,i2+3) = -1.d0
                ddeltaij(i2,i2) = 1.d0
                ddeltajk(i2,i2) = 1.d0
            else
                ddeltaik(i2,i2-6) = -1.d0
                ddeltaik(i2,i2) = 1.d0
                ddeltajk(i2,i2-3) = -1.d0
                ddeltajk(i2,i2) = 1.d0
            end if
        end do

        do i1 = 1,groupdim
            ineigh1 = group_neigh_local(i1,1)
            ineigh2 = group_neigh_local(i1,2)
            idx1 = neigh_idx_local(i1,1)
            idx2 = neigh_idx_local(i1,2)
            pairidx(1) = invneighboridx(iatom,jcount)
            pairidx(2) = invneighboridx(iatom,idx1)
            pairidx(3) = invneighboridx(iatom,idx2)
            deltaij(:) = 0.d0
            deltaik(:) = 0.d0
            deltajk(:) = 0.d0

            if(grouprij(i1,1).le.cutoff.and.grouprij(i1,2).le.cutoff.and.grouprij(i1,3).le.cutoff)then
                ! r_ij
                deltaij(1) = (xyzstruct_local(1,jcount)-lstb(ineigh1,1)) !(xi-xj)
                deltaij(2) = (xyzstruct_local(2,jcount)-lstb(ineigh1,2))
                deltaij(3) = (xyzstruct_local(3,jcount)-lstb(ineigh1,3))
                deltaij(4) = -1.d0*deltaij(1) !(xj-xi)
                deltaij(5) = -1.d0*deltaij(2)
                deltaij(6) = -1.d0*deltaij(3)
                do i2 = 1,6
                    drijd(i2) = deltaij(i2) / grouprij(i1,1)
                    dfcutijd(i2) = dfcut_local(i1,1) * drijd(i2)
                end do

                ! r_ik
                deltaik(1) = (xyzstruct_local(1,jcount)-lstb(ineigh2,1)) !(xi-xk)
                deltaik(2) = (xyzstruct_local(2,jcount)-lstb(ineigh2,2))
                deltaik(3) = (xyzstruct_local(3,jcount)-lstb(ineigh2,3))
                deltaik(7) = -1.d0*deltaik(1) !(xk-xi)
                deltaik(8) = -1.d0*deltaik(2)
                deltaik(9) = -1.d0*deltaik(3)
                drikd(:) = 0.d0
                dfcutikd(:) = 0.d0
                do i2 = 1,9
                    drikd(i2) = deltaik(i2) / grouprij(i1,2)
                    dfcutikd(i2) = dfcut_local(i1,2) * drikd(i2)
                end do

                ! r_jk
                deltajk(4) = (lstb(ineigh1,1)-lstb(ineigh2,1)) !(xj-xk)
                deltajk(5) = (lstb(ineigh1,2)-lstb(ineigh2,2))
                deltajk(6) = (lstb(ineigh1,3)-lstb(ineigh2,3))
                deltajk(7) = -1.d0*deltajk(4) !(xk-xj)
                deltajk(8) = -1.d0*deltajk(5)
                deltajk(9) = -1.d0*deltajk(6)
                drjkd(:) = 0.d0
                dfcutjkd(:) = 0.d0
                do i2 = 4,9
                    drjkd(i2) = deltajk(i2) / grouprij(i1,3)
                    dfcutjkd(i2) = dfcut_local(i1,3) * drjkd(i2)
                end do

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

                do i2 = 1,9
                    do i3 = 1,9
                        ddrijdd(i2,i3) = (ddeltaij(i2,i3)*grouprij(i1,1)-deltaij(i2)*drijd(i3))/&
                                (grouprij(i1,1)**2.d0)
                        ddrikdd(i2,i3) = (ddeltaik(i2,i3)*grouprij(i1,2)-deltaik(i2)*drikd(i3))/&
                                (grouprij(i1,2)**2.d0)
                        ddrjkdd(i2,i3) = (ddeltajk(i2,i3)*grouprij(i1,3)-deltajk(i2)*drjkd(i3))/&
                                (grouprij(i1,3)**2.d0)
                        ddfcutijdd(i2,i3) = ddfcut_local(i1,1)*drijd(i2)*drijd(i3) + dfcut_local(i1,1)*ddrijdd(i2,i3)
                        ddfcutikdd(i2,i3) = ddfcut_local(i1,2)*drikd(i2)*drikd(i3) + dfcut_local(i1,2)*ddrikdd(i2,i3)
                        ddfcutjkdd(i2,i3) = ddfcut_local(i1,3)*drjkd(i2)*drjkd(i3) + dfcut_local(i1,3)*ddrjkdd(i2,i3)
                    end do
                end do

                do i2 = 1,9
                    dfdxyz(i2) = (2.d0*grouprij(i1,3)*drjkd(i2) - 2.d0*grouprij(i1,1)*drijd(i2) - &
                                  2.d0*grouprij(i1,2)*drikd(i2))
                    dgdxyz(i2) = -2.d0*(drijd(i2)*grouprij(i1,2) + grouprij(i1,1)*drikd(i2))
                    do i3 = 1,9
                        ddfdxyz(i2,i3) = 2.d0*(drjkd(i3)*drjkd(i2) + grouprij(i1,3)*ddrjkdd(i2,i3) - &
                                drijd(i2)*drijd(i3) - grouprij(i1,1)*ddrijdd(i2,i3) - drikd(i2)*drikd(i3) - &
                                grouprij(i1,2)*ddrikdd(i2,i3))
                        ddgdxyz(i2,i3) = -2.d0*(drikd(i3)*drijd(i2)+drijd(i3)*drikd(i2)+grouprij(i1,2)*ddrijdd(i2,i3)+&
                                grouprij(i1,1)*ddrikdd(i2,i3))
                    end do
                end do

                temp5 = dsqrt(1.d0 - f**2/g**2)
                do i2 = 1,9
                    dtheta(i2) = -1.d0*rad2deg*(dfdxyz(i2)*g - f*dgdxyz(i2))/((g**2)*temp5)
                end do

                temp1 = theta - rshift_local(isf,isfg,eindex)
                temp2 = theta + rshift_local(isf,isfg,eindex)
                temp3 = -2.d0 * eta_local(isf,isfg,eindex)
                expxyz=(dexp(-eta_local(isf,isfg,eindex)*(temp1       )**2)&
                        +dexp(-eta_local(isf,isfg,eindex)*(temp2-360.d0)**2)&
                        +dexp(-eta_local(isf,isfg,eindex)*(temp2       )**2)&
                        +dexp(-eta_local(isf,isfg,eindex)*(temp1-360.d0)**2))

                do i2 = 1,9
                    dexpxyz(i2) =(temp1 * dexp(-eta_local(isf,isfg,eindex)*(temp1**2)) +&
                         temp2 * dexp(-eta_local(isf,isfg,eindex)*(temp2**2)) +&
                         (temp1-360.d0) * dexp(-eta_local(isf,isfg,eindex)*((temp1-360.d0)**2)) +&
                         (temp2-360.d0) * dexp(-eta_local(isf,isfg,eindex)*((temp2-360.d0)**2)))
                end do

                !! WRITE THIS AT THE END !!
                do i2 = 1,9
                    do i3 = 1,9
                        ddtheta(i2,i3) = -1.d0*rad2deg*(((ddfdxyz(i2,i3)*g + dfdxyz(i2)*dgdxyz(i3) - &
                                dgdxyz(i2)*dfdxyz(i3) - ddgdxyz(i2,i3)*f)*g*g - 2.d0*g*dgdxyz(i3)*&
                                (dfdxyz(i2)*g - dgdxyz(i2)*f))*temp5/(g**4) - (dfdxyz(i2)*g - f*dgdxyz(i2))*&
                                (-1.d0*(temp5**(-1))*(dfdxyz(i3)*g - f*dgdxyz(i3))*(f/(g**3)))/(g**2))/(temp5**2)

                        ddexpxyz(i2,i3) = (temp3*ddtheta(i2,i3)*dexpxyz(i2)) + &
                                      (temp3*dtheta(i2)*dtheta(i3)*expxyz) + &
                                      (temp3*temp3*dtheta(i2)*dtheta(i3)*((temp1**2)*dexp(-eta_local(isf,isfg,eindex)*(temp1**2)) +&
                                      (temp2**2)*dexp(-eta_local(isf,isfg,eindex)*(temp2**2)) +&
                                      (temp1-360.d0)**2 * dexp(-eta_local(isf,isfg,eindex)*((temp1-360.d0)**2)) +&
                                      (temp2-360.d0)**2 * dexp(-eta_local(isf,isfg,eindex)*((temp2-360.d0)**2))))
                    end do
                end do

                do i2 = 1,9
                    dexpxyz(i2) = temp3*dexpxyz(i2)*dtheta(i2)
                end do


                !! Second derivatives
                do i2 = 1,3 ! atom,n1,n2
                    do i3 = 1,3 ! xyz
                        do i4 = 1,3 ! atom,n1,n2
                            do i5 = 1,3 ! xyz
                                ddsfuncddxyz_temp(pairidx(i2),i3,pairidx(i4),i5) = &
                                ddsfuncddxyz_temp(pairidx(i2),i3,pairidx(i4),i5) + &

                                ddexpxyz((i2-1)*3+i3,(i4-1)*3+i5)*fcut_local(i1,1)*fcut_local(i1,2)*fcut_local(i1,3)+ &
                                dexpxyz((i2-1)*3+i3)*dfcutijd((i4-1)*3+i5)*fcut_local(i1,2)*fcut_local(i1,3)+ &
                                dexpxyz((i2-1)*3+i3)*dfcutikd((i4-1)*3+i5)*fcut_local(i1,1)*fcut_local(i1,3)+ &
                                dexpxyz((i2-1)*3+i3)*dfcutjkd((i4-1)*3+i5)*fcut_local(i1,2)*fcut_local(i1,1)+ &

                                dexpxyz((i4-1)*3+i5)*dfcutijd((i2-1)*3+i3)*fcut_local(i1,2)*fcut_local(i1,3)+ &
                                expxyz*dfcutijd((i2-1)*3+i3)*fcut_local(i1,2)*dfcutjkd((i4-1)*3+i5)+ &
                                expxyz*dfcutijd((i2-1)*3+i3)*dfcutikd((i4-1)*3+i5)*fcut_local(i1,3)+ &
                                expxyz*ddfcutijdd((i2-1)*3+i3,(i4-1)*3+i5)*fcut_local(i1,2)*fcut_local(i1,3)+ &

                                dexpxyz((i4-1)*3+i5)*dfcutikd((i2-1)*3+i3)*fcut_local(i1,1)*fcut_local(i1,3)+ &
                                expxyz*dfcutikd((i2-1)*3+i3)*fcut_local(i1,1)*dfcutjkd((i4-1)*3+i5)+ &
                                expxyz*dfcutikd((i2-1)*3+i3)*dfcutijd((i4-1)*3+i5)*fcut_local(i1,3)+ &
                                expxyz*ddfcutikdd((i2-1)*3+i3,(i4-1)*3+i5)*fcut_local(i1,1)*fcut_local(i1,3)+ &

                                dexpxyz((i4-1)*3+i5)*dfcutjkd((i2-1)*3+i3)*fcut_local(i1,2)*fcut_local(i1,1)+ &
                                expxyz*dfcutijd((i4-1)*3+i5)*fcut_local(i1,2)*dfcutjkd((i2-1)*3+i3)+ &
                                expxyz*dfcutikd((i4-1)*3+i5)*dfcutjkd((i2-1)*3+i3)*fcut_local(i1,1)+ &
                                expxyz*ddfcutjkdd((i2-1)*3+i3,(i4-1)*3+i5)*fcut_local(i1,2)*fcut_local(i1,1)
                            end do
                        end do
                    end do
                end do
            end if

        end do

    elseif(ftype.eq.9)then

        ddeltaij(:,:)=0.d0
        ddeltaik(:,:)=0.d0
        ddeltajk(:,:)=0.d0

        do i2 = 1,9
            if(i2.le.3)then
                ddeltaij(i2,i2+3) = -1.d0
                ddeltaij(i2,i2) = 1.d0
                ddeltaik(i2,i2+6) = -1.d0
                ddeltaik(i2,i2) = 1.d0
            else if (i2.gt.3.and.i2.le.6)then
                ddeltaij(i2,i2-3) = -1.d0
                ddeltajk(i2,i2+3) = -1.d0
                ddeltaij(i2,i2) = 1.d0
                ddeltajk(i2,i2) = 1.d0
            else
                ddeltaik(i2,i2-6) = -1.d0
                ddeltaik(i2,i2) = 1.d0
                ddeltajk(i2,i2-3) = -1.d0
                ddeltajk(i2,i2) = 1.d0
            end if
        end do

        do i1 = 1,groupdim
            ineigh1 = group_neigh_local(i1,1)
            ineigh2 = group_neigh_local(i1,2)
            idx1 = neigh_idx_local(i1,1)
            idx2 = neigh_idx_local(i1,2)
            pairidx(1) = invneighboridx(iatom,jcount)
            pairidx(2) = invneighboridx(iatom,idx1)
            pairidx(3) = invneighboridx(iatom,idx2)
            deltaij(:) = 0.d0
            deltaik(:) = 0.d0
            deltajk(:) = 0.d0

            if(grouprij(i1,1).le.cutoff.and.grouprij(i1,2).le.cutoff)then

                ! r_ij
                deltaij(1) = (xyzstruct_local(1,jcount)-lstb(ineigh1,1)) !(xi-xj)
                deltaij(2) = (xyzstruct_local(2,jcount)-lstb(ineigh1,2))
                deltaij(3) = (xyzstruct_local(3,jcount)-lstb(ineigh1,3))
                deltaij(4) = -1.d0*deltaij(1) !(xj-xi)
                deltaij(5) = -1.d0*deltaij(2)
                deltaij(6) = -1.d0*deltaij(3)
                drijd(:) = 0.d0
                dfcutijd(:) = 0.d0
                do i2 = 1,6
                    drijd(i2) = deltaij(i2) / grouprij(i1,1)
                    dfcutijd(i2) = dfcut_local(i1,1) * drijd(i2)
                end do

                ! r_ik
                deltaik(1) = (xyzstruct_local(1,jcount)-lstb(ineigh2,1)) !(xi-xk)
                deltaik(2) = (xyzstruct_local(2,jcount)-lstb(ineigh2,2))
                deltaik(3) = (xyzstruct_local(3,jcount)-lstb(ineigh2,3))
                deltaik(7) = -1.d0*deltaik(1) !(xk-xi)
                deltaik(8) = -1.d0*deltaik(2)
                deltaik(9) = -1.d0*deltaik(3)
                drikd(:) = 0.d0
                dfcutikd(:) = 0.d0
                do i2 = 1,9
                    drikd(i2) = deltaik(i2) / grouprij(i1,2)
                    dfcutikd(i2) = dfcut_local(i1,2) * drikd(i2)
                end do

                ! r_jk
                deltajk(4) = (lstb(ineigh1,1)-lstb(ineigh2,1)) !(xj-xk)
                deltajk(5) = (lstb(ineigh1,2)-lstb(ineigh2,2))
                deltajk(6) = (lstb(ineigh1,3)-lstb(ineigh2,3))
                deltajk(7) = -1.d0*deltajk(4) !(xk-xj)
                deltajk(8) = -1.d0*deltajk(5)
                deltajk(9) = -1.d0*deltajk(6)
                drjkd(:) = 0.d0
                dfcutjkd(:) = 0.d0
                do i2 = 4,9
                    drjkd(i2) = deltajk(i2) / grouprij(i1,3)
                    !dfcutjkd(i2) = dfcut_local(i1,3) * drjkd(i2)
                end do


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

                do i2 = 1,9
                    do i3 = 1,9
                        ddrijdd(i2,i3) = (ddeltaij(i2,i3)*grouprij(i1,1)-deltaij(i2)*drijd(i3))/&
                                (grouprij(i1,1)**2.d0)
                        ddrikdd(i2,i3) = (ddeltaik(i2,i3)*grouprij(i1,2)-deltaik(i2)*drikd(i3))/&
                                (grouprij(i1,2)**2.d0)
                        ddrjkdd(i2,i3) = (ddeltajk(i2,i3)*grouprij(i1,3)-deltajk(i2)*drjkd(i3))/&
                                (grouprij(i1,3)**2.d0)
                        ddfcutijdd(i2,i3) = ddfcut_local(i1,1)*drijd(i2)*drijd(i3) + dfcut_local(i1,1)*ddrijdd(i2,i3)
                        ddfcutikdd(i2,i3) = ddfcut_local(i1,2)*drikd(i2)*drikd(i3) + dfcut_local(i1,2)*ddrikdd(i2,i3)
                        !ddfcutjkdd(i2,i3) = ddfcut_local(i1,3)*drjkd(i2)*drjkd(i3) + dfcut_local(i1,3)*ddrjkdd(i2,i3)
                    end do
                end do

                do i2 = 1,9
                    dfdxyz(i2) = lambda_local(isf,isfg,eindex)*(2.d0*grouprij(i1,3)*drjkd(i2) - &
                            2.d0*grouprij(i1,1)*drijd(i2) - 2.d0*grouprij(i1,2)*drikd(i2))

                    dgdxyz(i2) = -2.d0*(drijd(i2)*grouprij(i1,2) + grouprij(i1,1)*drikd(i2))
                    do i3 = 1,9
                        ddfdxyz(i2,i3) = 2.d0*lambda_local(isf,isfg,eindex)*(drjkd(i3)*drjkd(i2) + grouprij(i1,3)*&
                                ddrjkdd(i2,i3) - drijd(i2)*drijd(i3) - grouprij(i1,1)*ddrijdd(i2,i3) - &
                                drikd(i2)*drikd(i3) - grouprij(i1,2)*ddrikdd(i2,i3))

                        ddgdxyz(i2,i3) = -2.d0*(drikd(i3)*drijd(i2)+drijd(i3)*drikd(i2)+grouprij(i1,2)*ddrijdd(i2,i3)+&
                                grouprij(i1,1)*ddrikdd(i2,i3))
                        !ddgdxyz(i2,i3) = -2.d0*(drikd(i3)*drijd(i2) + drijd(i3)*drikd(i2))
                    end do
                end do

                temp5 = 1.d0 + f/g
                if(temp5.lt.0.0d0)costheta=0.0d0
                temp4 = temp4*zeta_local(isf,isfg,eindex)*temp5**(zeta_local(isf,isfg,eindex)-1)/g**2
                do i2 = 1,9
                    dcostheta(i2) = temp4*(dfdxyz(i2)*g - f*dgdxyz(i2))
                end do

                ! Exponentials
                expxyz = dexp(-eta_local(isf,isfg,eindex)*(grouprij(i1,1)**2+grouprij(i1,2)**2))
                temp3 = -eta_local(isf,isfg,eindex)*2.0d0*expxyz
                do i2 = 1,9
                    dexpxyz(i2) = (grouprij(i1,1)*drijd(i2) + grouprij(i1,2)*drikd(i2))*temp3
                end do

                temp1 = (2.d0**(1.d0-zeta_local(isf,isfg,eindex)))*zeta_local(isf,isfg,eindex)
                temp2 = 1.d0 + f/g

                do i2 = 1,9
                    do i3 = 1,9
                        ddcostheta(i2,i3)= temp1 * ((zeta_local(isf,isfg,eindex)-1) * &
                                (1+f/g)**(zeta_local(isf,isfg,eindex)-2) * (dfdxyz(i2)*g - f*dgdxyz(i2)) * &
                                (dfdxyz(i3)*g - f*dgdxyz(i3)) + (1+f/g)**(zeta_local(isf,isfg,eindex)-1) * &
                                ((ddfdxyz(i2,i3)*g + dfdxyz(i2)*dgdxyz(i3) - dgdxyz(i2)*dfdxyz(i3) - &
                                ddgdxyz(i2,i3)*f)*g*g - (dfdxyz(i2)*g - dgdxyz(i2)*f)*2.d0*g*dgdxyz(i3))) / g**4.d0

                        ddexpxyz(i2,i3) = -2.d0*eta_local(isf,isfg,eindex)*(grouprij(i1,1)*drijd(i3) + grouprij(i1,2)*&
                                drikd(i3))*dexpxyz(i2) + (-2.d0*expxyz*eta_local(isf,isfg,&
                                eindex)*(drijd(i3)*drijd(i2)+grouprij(i1,1)*ddrijdd(i2,i3) + drikd(i3)*drikd(i2) + &
                                grouprij(i1,2)*ddrikdd(i2,i3)))
                    end do
                end do


                !! Second derivatives
                do i2 = 1,3 ! atom,n1,n2
                    do i3 = 1,3 ! xyz
                        do i4 = 1,3 ! atom,n1,n2
                            do i5 = 1,3 ! xyz
                                ddsfuncddxyz_temp(pairidx(i2),i3,pairidx(i4),i5) = &
                                ddsfuncddxyz_temp(pairidx(i2),i3,pairidx(i4),i5) + &

                                ddcostheta((i2-1)*3+i3,(i4-1)*3+i5) * expxyz * fcut_local(i1,1) * fcut_local(i1,2) + &
                                dcostheta((i2-1)*3+i3) * dexpxyz((i4-1)*3+i5) * fcut_local(i1,1) * fcut_local(i1,2) + &
                                dcostheta((i2-1)*3+i3) * expxyz * dfcutijd((i4-1)*3+i5) * fcut_local(i1,2) + &
                                dcostheta((i2-1)*3+i3) * expxyz * fcut_local(i1,1) * dfcutikd((i4-1)*3+i5) + &

                                dcostheta((i4-1)*3+i5) * dexpxyz((i2-1)*3+i3) * fcut_local(i1,1) * fcut_local(i1,2) + &
                                ddexpxyz((i2-1)*3+i3,(i4-1)*3+i5) * costheta * fcut_local(i1,1) * fcut_local(i1,2) + &
                                dexpxyz((i2-1)*3+i3) * costheta * dfcutijd((i4-1)*3+i5) * fcut_local(i1,2) + &
                                dexpxyz((i2-1)*3+i3) * costheta * fcut_local(i1,1) * dfcutikd((i4-1)*3+i5) + &

                                dcostheta((i4-1)*3+i5) * expxyz * dfcutijd((i2-1)*3+i3) * fcut_local(i1,2) + &
                                dexpxyz((i4-1)*3+i5) * costheta * dfcutijd((i2-1)*3+i3) * fcut_local(i1,2) + &
                                expxyz * costheta * ddfcutijdd((i2-1)*3+i3,(i4-1)*3+i5) * fcut_local(i1,2) + &
                                expxyz * costheta * dfcutijd((i2-1)*3+i3) * dfcutikd((i4-1)*3+i5) + &

                                dcostheta((i4-1)*3+i5) * expxyz * dfcutikd((i2-1)*3+i3) * fcut_local(i1,1) + &
                                dexpxyz((i4-1)*3+i5) * costheta * dfcutikd((i2-1)*3+i3) * fcut_local(i1,1) + &
                                expxyz * costheta * dfcutijd((i4-1)*3+i5) * dfcutikd((i2-1)*3+i3) + &
                                expxyz * costheta * fcut_local(i1,1) * ddfcutikdd((i2-1)*3+i3,(i4-1)*3+i5)
                            end do
                        end do
                    end do
                end do

            end if

        end do
    end if

end subroutine getsymfunctionsecondderivatives
