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
!! Purpose: calculate symmetry functions for a given atom (and derivatives if requested)
!! USING SYMMETRY FUNCTION GROUPS (Emir)
!!
subroutine getatomsymfunctions_sfg(iatom,isfg,isf,eindex,natoms,&
        rmin,lrmin,maxnum_funcvalues_local,maxnum_sfgroups_local,nelem,ftype,cutoff,&
        symfunction_temp,eta_local,rshift_local,lambda_local,zeta_local,groupdim,&
        maxgroupdim,grouprij,groupcutoff)

    use fileunits
    use nnconstants

    implicit none

    integer i1 !internal
    integer iatom ! in
    integer isfg ! in
    integer isf ! in
    integer eindex ! in
    integer nelem ! in
    integer natoms ! in
    integer ftype ! in
    integer groupdim ! in
    integer maxgroupdim ! in
    integer maxnum_funcvalues_local ! in
    integer maxnum_sfgroups_local ! in

    real*8  symfunction_temp(maxnum_funcvalues_local) !in/out
    real*8  groupcutoff(maxgroupdim,3) ! in
    real*8  grouprij(maxgroupdim,3) ! in
    real*8  eta_local(maxnum_funcvalues_local,maxnum_sfgroups_local,nelem) ! in
    real*8  zeta_local(maxnum_funcvalues_local,maxnum_sfgroups_local,nelem) ! in
    real*8  rshift_local(maxnum_funcvalues_local,maxnum_sfgroups_local,nelem) ! in
    real*8  lambda_local(maxnum_funcvalues_local,maxnum_sfgroups_local,nelem) ! in
    real*8  cutoff ! in
    real*8  rmin ! in
    real*8  costheta ! internal
    real*8  theta ! internal
    real*8  f,g,temp1,temp2 ! internal
    real*8  expxyz ! internal

    logical lrmin ! in


    if(ftype.eq.1)then ! radial SF (only cutoff function)
        do i1 = 1,groupdim
            if(grouprij(i1,1).le.cutoff)then
                symfunction_temp(isf) = symfunction_temp(isf) + groupcutoff(i1,1)
            end if
        end do

    elseif(ftype.eq.2)then ! radial SF (standard)

        do i1 = 1,groupdim
            if(grouprij(i1,1).le.cutoff)then
                symfunction_temp(isf) = symfunction_temp(isf) + dexp(-1.d0*eta_local(isf,isfg,eindex)*&
                        (grouprij(i1,1) - rshift_local(isf,isfg,eindex))**2) * groupcutoff(i1,1)
            end if
        end do

    elseif(ftype.eq.3)then ! angular SF (narrow)

        do i1 = 1,groupdim
            if(grouprij(i1,1).le.cutoff.and.grouprij(i1,2).le.cutoff.and.grouprij(i1,3).le.cutoff)then
                if(grouprij(i1,3).le.rmin) then
                    write(ounit,*)'Atoms too close: rjk=',grouprij(i1,3)
                    lrmin=.false.
                endif
                ! cosine term
                costheta = (lambda_local(isf,isfg,eindex)*(grouprij(i1,3)**2 - grouprij(i1,1)**2 - &
                        grouprij(i1,2)**2))/(-2.d0 * grouprij(i1,1) * grouprij(i1,2))
                costheta = 1.d0 + costheta ! avoid negative values
                !! avoid negative values due to numerical noise (ASK)
                if(costheta.le.0.d0) then
                    costheta = 0.d0
                else
                    costheta = 2.d0**(1.d0-zeta_local(isf,isfg,eindex))* &
                            (costheta**zeta_local(isf,isfg,eindex))
                endif
                ! exponential term
                expxyz = dexp(-eta_local(isf,isfg,eindex) * (grouprij(i1,1)**2+&
                        grouprij(i1,2)**2 + grouprij(i1,3)**2))

                symfunction_temp(isf) = symfunction_temp(isf) + costheta *&
                        expxyz * groupcutoff(i1,1) * groupcutoff(i1,2) * groupcutoff(i1,3)
            end if
        end do


    elseif(ftype.eq.4)then  ! radial SF

        do i1 = 1,groupdim
            if(grouprij(i1,1).le.cutoff)then
                symfunction_temp(isf) = symfunction_temp(isf) + dcos(eta_local(isf,isfg,eindex)*grouprij(i1,1))*&
                        groupcutoff(i1,1)
            end if
        end do

    elseif(ftype.eq.8)then ! angular SF (shifted)

        do i1 = 1,groupdim
            if(grouprij(i1,1).le.cutoff.and.grouprij(i1,2).le.cutoff.and.grouprij(i1,3).le.cutoff)then
                if(grouprij(i1,3).le.rmin) then
                    write(ounit,*)'Atoms too close: rjk=',grouprij(i1,3)
                    lrmin=.false.
                endif

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

                temp1=theta-rshift_local(isf,isfg,eindex)
                temp2=theta+rshift_local(isf,isfg,eindex)
                expxyz=(dexp(-eta_local(isf,isfg,eindex)*(temp1       )**2)&
                        +dexp(-eta_local(isf,isfg,eindex)*(temp2-360.d0)**2)&
                        +dexp(-eta_local(isf,isfg,eindex)*(temp2       )**2)&
                        +dexp(-eta_local(isf,isfg,eindex)*(temp1-360.d0)**2))
                symfunction_temp(isf) = symfunction_temp(isf)&
                        +expxyz * groupcutoff(i1,1) * groupcutoff(i1,2) * groupcutoff(i1,3)


            end if
        end do

    elseif(ftype.eq.9)then ! angular SF (wide)

        do i1 = 1,groupdim
            if(grouprij(i1,1).le.cutoff.and.grouprij(i1,2).le.cutoff)then
                if(grouprij(i1,3).le.rmin) then
                    write(ounit,*)'Atoms too close: rjk=',grouprij(i1,3)
                    lrmin=.false.
                endif
                ! cosine term
                costheta = (lambda_local(isf,isfg,eindex)*(grouprij(i1,3)**2 - grouprij(i1,1)**2 - &
                        grouprij(i1,2)**2))/(-2.d0 * grouprij(i1,1) * grouprij(i1,2))
                costheta = 1.d0 + costheta ! avoid negative values
                !! avoid negative values due to numerical noise (ASK)
                if(costheta.le.0.d0) then
                    costheta = 0.d0
                else
                    costheta = 2.d0**(1.d0-zeta_local(isf,isfg,eindex))* &
                            (costheta**zeta_local(isf,isfg,eindex))
                endif
                ! exponential term
                expxyz = dexp(-eta_local(isf,isfg,eindex) * (grouprij(i1,1)**2+&
                        grouprij(i1,2)**2))

                symfunction_temp(isf) = symfunction_temp(isf) + costheta *&
                        expxyz * groupcutoff(i1,1) * groupcutoff(i1,2)
            end if
        end do

    end if

end subroutine getatomsymfunctions_sfg
