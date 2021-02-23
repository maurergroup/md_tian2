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
!! Emir
!! called by:
!! - calconefunction_atomic_sfg.f90
!! This calculates cutoff function values for each symmetry function group

subroutine getgroupcutoffvalues(cutoff_type,cutoff_alpha,ielem,isfg,&
           ftype,maxnum_sfg,nelem,fc_sfg,groupdim,maxdim,group_rij,&
           group_fcut,group_fcut_del,group_fcut_del2)


    implicit none

    integer i1 ! internal
    integer ielem ! in
    integer isfg ! in
    integer cutoff_type ! in
    integer ftype ! in
    integer groupdim ! in
    integer maxdim ! in
    integer maxnum_sfg
    integer nelem
    integer fc_sfg(maxnum_sfg,nelem)

    real*8  cutoff_alpha ! in
    real*8  fcut_temp1 ! internal
    real*8  fcut_temp2 ! internal
    real*8  fcut_temp3 ! internal
    real*8  dfcut_temp1 ! internal
    real*8  dfcut_temp2 ! internal
    real*8  dfcut_temp3 ! internal
    real*8  ddfcut_temp1 ! internal
    real*8  ddfcut_temp2 ! internal
    real*8  ddfcut_temp3 ! internal
    real*8  dummy ! this is the derivative FIX ME (this is required in mode 2 but not now)
    real*8  group_rij(maxdim,3) ! in
    real*8  group_fcut(maxdim,3) ! out
    real*8  group_fcut_del(maxdim,3) ! out
    real*8  group_fcut_del2(maxdim,3) ! out (for Hessian)

    fcut_temp1 = 0.d0
    fcut_temp2 = 0.d0
    fcut_temp3 = 0.d0
    dfcut_temp1 = 0.d0
    dfcut_temp2 = 0.d0
    dfcut_temp3 = 0.d0
    ddfcut_temp1 = 0.d0
    ddfcut_temp2 = 0.d0
    ddfcut_temp3 = 0.d0

    do i1 = 1,groupdim
        if(ftype.eq.1.or.ftype.eq.2.or.ftype.eq.4)then ! radial SF (only cutoff function)
            call getcutoff(cutoff_type,cutoff_alpha,maxnum_sfg,nelem,&
                    isfg,ielem,fc_sfg,group_rij(i1,1),fcut_temp1,dfcut_temp1,ddfcut_temp1)
        elseif(ftype.eq.3.or.ftype.eq.8)then ! narrow angular SF
            call getcutoff(cutoff_type,cutoff_alpha,maxnum_sfg,nelem,&
                    isfg,ielem,fc_sfg,group_rij(i1,1),fcut_temp1,dfcut_temp1,ddfcut_temp1)
            call getcutoff(cutoff_type,cutoff_alpha,maxnum_sfg,nelem,&
                    isfg,ielem,fc_sfg,group_rij(i1,2),fcut_temp2,dfcut_temp2,ddfcut_temp2)
            call getcutoff(cutoff_type,cutoff_alpha,maxnum_sfg,nelem,&
                    isfg,ielem,fc_sfg,group_rij(i1,3),fcut_temp3,dfcut_temp3,ddfcut_temp3)
        elseif(ftype.eq.9)then ! wide angular SF
            call getcutoff(cutoff_type,cutoff_alpha,maxnum_sfg,nelem,&
                    isfg,ielem,fc_sfg,group_rij(i1,1),fcut_temp1,dfcut_temp1,ddfcut_temp1)
            call getcutoff(cutoff_type,cutoff_alpha,maxnum_sfg,nelem,&
                    isfg,ielem,fc_sfg,group_rij(i1,2),fcut_temp2,dfcut_temp2,ddfcut_temp2)

        end if
        group_fcut(i1,1) = fcut_temp1
        group_fcut(i1,2) = fcut_temp2
        group_fcut(i1,3) = fcut_temp3
        group_fcut_del(i1,1) = dfcut_temp1
        group_fcut_del(i1,2) = dfcut_temp2
        group_fcut_del(i1,3) = dfcut_temp3
        group_fcut_del2(i1,1) = ddfcut_temp1
        group_fcut_del2(i1,2) = ddfcut_temp2
        group_fcut_del2(i1,3) = ddfcut_temp3
    end do

end subroutine getgroupcutoffvalues
