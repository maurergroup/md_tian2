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

! This subroutine extracts neighbors or neighbor pairs associated with
! a given symmetry function group - Emir

subroutine getgroupneighlist(low_neigh_bound,up_neigh_bound,lstb_local,lstc_local,&
        lste_local,type,sf_type,listdim,zneigh1,zneigh2,group_listdim,maxdim,&
        group_neigh_list_local,group_neigh_idx_local,group_radial_distances_local)

    use fileunits

    implicit none

    integer type ! in (O:radial, 1:angular)
    integer sf_type ! in
    integer group_neigh_list_local(maxdim,type+1)
    integer group_neigh_idx_local(maxdim,type+1)
    integer zneigh1 ! in
    integer zneigh2 ! in
    integer listdim ! in
    integer group_listdim !in/out
    integer maxdim ! in
    integer lstc_local(listdim) ! in
    integer lste_local(listdim) ! in
    integer i1,i2 ! in
    integer low_neigh_bound,up_neigh_bound

    real*8  lstb_local(listdim,4) ! in
    real*8  group_radial_distances_local(maxdim,3)



    if(type.eq.0)then ! radial SFG
        do i1 = low_neigh_bound,up_neigh_bound
            if(lste_local(i1).eq.zneigh1)then
                group_listdim = group_listdim + 1
                group_neigh_list_local(group_listdim,1) = i1 ! atom ID
                group_neigh_idx_local(group_listdim,1) = lstc_local(i1) ! neighbor index
                group_radial_distances_local(group_listdim,1) = lstb_local(i1,4) ! r_ij
            end if
        end do
    else ! angular SFG
        do i1 = low_neigh_bound,up_neigh_bound
            if(lste_local(i1).eq.zneigh1.and.zneigh2.ne.zneigh1)then
                do i2 = low_neigh_bound,up_neigh_bound
                    if(lste_local(i2).eq.zneigh2)then
                        group_listdim = group_listdim + 1
                        group_neigh_list_local(group_listdim,1) = i1
                        group_neigh_list_local(group_listdim,2) = i2
                        group_neigh_idx_local(group_listdim,1) = lstc_local(i1)
                        group_neigh_idx_local(group_listdim,2) = lstc_local(i2)
                        group_radial_distances_local(group_listdim,1) = lstb_local(i1,4) ! r_ij
                        group_radial_distances_local(group_listdim,2) = lstb_local(i2,4) ! r_ik
                        group_radial_distances_local(group_listdim,3) = get_rjk(lstb_local(i1,1),lstb_local(i2,1),&
                                                            lstb_local(i1,2),lstb_local(i2,2),lstb_local(i1,3),&
                                                            lstb_local(i2,3))
                    end if
                end do
             else if(lste_local(i1).eq.zneigh1.and.zneigh2.eq.zneigh1)then
                do i2 = (i1+1),up_neigh_bound
                    if(lste_local(i2).eq.zneigh2)then
                        group_listdim = group_listdim + 1
                        group_neigh_list_local(group_listdim,1) = i1
                        group_neigh_list_local(group_listdim,2) = i2
                        group_neigh_idx_local(group_listdim,1) = lstc_local(i1)
                        group_neigh_idx_local(group_listdim,2) = lstc_local(i2)
                        group_radial_distances_local(group_listdim,1) = lstb_local(i1,4) ! r_ij
                        group_radial_distances_local(group_listdim,2) = lstb_local(i2,4) ! r_ik
                        group_radial_distances_local(group_listdim,3) = get_rjk(lstb_local(i1,1),lstb_local(i2,1),&
                                    lstb_local(i1,2),lstb_local(i2,2),lstb_local(i1,3),&
                                    lstb_local(i2,3))
                    end if
                end do
            end if
        end do

    end if

    contains
        !! For calculating r_jk
        real*8 function get_rjk(x1,x2,y1,y2,z1,z2)
        real*8 rjksq,x1,x2,y1,y2,z1,z2
        rjksq = (x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2
        get_rjk=dsqrt(rjksq)
        end


end subroutine getgroupneighlist


