!######################################################################
! This module is part of
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
module symfunctiongroups
    !! initially just for short_atomic (Emir)
    implicit none


    integer, dimension(:,:), allocatable :: function_type_short_atomic_sfg
    integer, dimension(:,:,:), allocatable :: neighbors_sfg
    integer, dimension(:,:), allocatable :: sf_count_sfg
    integer, dimension(:), allocatable :: num_groups_short_atomic
    integer, dimension(:,:,:), allocatable :: sf_index_sfg ! for re-mapping


    real*8, dimension(:,:) , allocatable :: funccutoff_short_atomic_sfg
    real*8, dimension(:,:,:) , allocatable :: eta_short_atomic_sfg
    real*8, dimension(:,:,:) , allocatable :: zeta_short_atomic_sfg
    real*8, dimension(:,:,:) , allocatable :: lambda_short_atomic_sfg
    real*8, dimension(:,:,:) , allocatable :: rshift_short_atomic_sfg
    real*8 maxcutoff_short_atomic_sfg

contains

    ! check this allocation later
    subroutine allocatesymfunctiongroups()
        !!
        use nnflags
        use globaloptions
        !!
        implicit none

        maxnum_sfgroups_short_atomic = maxnum_funcvalues_short_atomic ! CHECK ME

        !!
        if(lshort.and.(nn_type_short.eq.1).and.lusesfgroups)then
            allocate(sf_count_sfg(maxnum_sfgroups_short_atomic,nelem))
            sf_count_sfg(:,:)=0
            allocate(sf_index_sfg(maxnum_funcvalues_short_atomic,maxnum_sfgroups_short_atomic,nelem))
            sf_index_sfg(:,:,:)=0
            allocate(num_groups_short_atomic(nelem))
            num_groups_short_atomic(:)=0
            allocate(function_type_short_atomic_sfg(maxnum_sfgroups_short_atomic,nelem))
            function_type_short_atomic_sfg(:,:)=0
            allocate(neighbors_sfg(maxnum_sfgroups_short_atomic,2,nelem))
            neighbors_sfg(:,:,:)=0
            allocate(funccutoff_short_atomic_sfg(maxnum_sfgroups_short_atomic,nelem))
            funccutoff_short_atomic_sfg(:,:)=0.0d0
            allocate(eta_short_atomic_sfg(maxnum_funcvalues_short_atomic,maxnum_sfgroups_short_atomic,nelem))
            eta_short_atomic_sfg(:,:,:)=0.0d0
            allocate(zeta_short_atomic_sfg(maxnum_funcvalues_short_atomic,maxnum_sfgroups_short_atomic,nelem))
            zeta_short_atomic_sfg(:,:,:)=0.0d0
            allocate(lambda_short_atomic_sfg(maxnum_funcvalues_short_atomic,maxnum_sfgroups_short_atomic,nelem))
            lambda_short_atomic_sfg(:,:,:)=0.0d0
            allocate(rshift_short_atomic_sfg(maxnum_funcvalues_short_atomic,maxnum_sfgroups_short_atomic,nelem))
            rshift_short_atomic_sfg(:,:,:)=0.0d0
        endif


    end subroutine allocatesymfunctiongroups


    subroutine setupsymfunctiongroups()

        !! Loops over already sorted SFs and assign them one by one to
        !! their associated SF group. Max num of SF groups was already
        !! calculated in the first read of input.nn to allocate the arrays
        !! beforehand.

        use symfunctions
        use globaloptions

        implicit none

        integer i1 ! internal
        integer i2 ! internal
        integer dum_type ! internal
        integer type_local ! internal
        integer num_sfg_local ! internal
        integer id_sfg_local ! internal
        integer num_sf_local ! internal
        integer dum_create ! internal
        integer group_type ! internal (radial:0, angular:1)
        logical lcreategroup ! internal

        num_sfg_local = 0
        id_sfg_local = 1
        num_sf_local = 0
        dum_type = 0
        maxcutoff_short_atomic_sfg = 0
        lcreategroup = .false.

        do i1 = 1,nelem
            do i2 = 1,maxnum_funcvalues_short_atomic
                if(dum_type.ne.function_type_short_atomic(i2,i1))then
                    if(dum_type.ne.0)then ! for the first angular SF
                        id_sfg_local = num_sfg_local + 1
                    end if
                    num_sfg_local = num_sfg_local + 1 !this should be input to next line
                    call create_sfg(i1,num_sfg_local,function_type_short_atomic(i2,i1)&
                        ,funccutoff_short_atomic(i2,i1),symelement_short_atomic(i2,1,i1)&
                            ,symelement_short_atomic(i2,2,i1))
                    dum_type = function_type_short_atomic(i2,i1)
                    num_groups_short_atomic(i1)=num_groups_short_atomic(i1)+1
                    sf_count_sfg(id_sfg_local,i1)=sf_count_sfg(id_sfg_local,i1)+1
                    sf_index_sfg(sf_count_sfg(id_sfg_local,i1),id_sfg_local,i1) = i2
                    call check_type(function_type_short_atomic(i2,i1),group_type)
                    call add_member(function_type_short_atomic(i2,i1),i1,i2,&
                            sf_count_sfg(id_sfg_local,i1),id_sfg_local)
                else
                    call check_type(function_type_short_atomic(i2,i1),group_type)
                    call check_neighborsandcutoff(group_type,symelement_short_atomic(i2,1,i1),&
                            symelement_short_atomic(i2,2,i1),i1,funccutoff_short_atomic(i2,i1),&
                            function_type_short_atomic(i2,i1),id_sfg_local,num_sfg_local,lcreategroup)
                    if(lcreategroup)then
                        num_sfg_local = num_sfg_local + 1 !this should be input to next line
                        call create_sfg(i1,num_sfg_local,function_type_short_atomic(i2,i1)&
                                ,funccutoff_short_atomic(i2,i1),symelement_short_atomic(i2,1,i1)&
                                ,symelement_short_atomic(i2,2,i1))
                        num_groups_short_atomic(i1)=num_groups_short_atomic(i1)+1
                    endif
                    dum_type = function_type_short_atomic(i2,i1)
                    sf_count_sfg(id_sfg_local,i1) = sf_count_sfg(id_sfg_local,i1) + 1
                    sf_index_sfg(sf_count_sfg(id_sfg_local,i1),id_sfg_local,i1) = i2
                    call add_member(function_type_short_atomic(i2,i1),i1,i2,&
                            sf_count_sfg(id_sfg_local,i1),id_sfg_local)
                endif
                lcreategroup=.false. ! set to its default value
            end do
            num_sfg_local = 0
            id_sfg_local = 1
            dum_type = 0
        end do

        call writesfgroups()

    end subroutine setupsymfunctiongroups

    subroutine create_sfg(ielem,num_sfg,ftype,fcutoff,neigh1,neigh2)

        use symfunctions
        use fileunits

        implicit none

        integer num_sfg !in
        integer ielem !in
        integer ftype !in
        real*8 fcutoff !in
        integer neigh1 !in
        integer neigh2 !in

        function_type_short_atomic_sfg(num_sfg,ielem) = ftype
        funccutoff_short_atomic_sfg(num_sfg,ielem) = fcutoff
        neighbors_sfg(num_sfg,1,ielem) = neigh1
        neighbors_sfg(num_sfg,2,ielem) = neigh2

        ! record the maximum cutoff radius
        if(fcutoff.gt.maxcutoff_short_atomic_sfg)then
            maxcutoff_short_atomic_sfg = fcutoff
        end if

    end subroutine create_sfg

    subroutine add_member(ftype,celem,csf,isf,idsfg)

        use symfunctions
        use fileunits

        implicit none

        integer ftype ! in (SF type)
        integer celem ! in (Element counter)
        integer csf ! in (SF counter)
        integer isf ! in (SF index)
        integer idsfg ! in (SFG id)
        !write(ounit,*)type,celem,csf,isf,idsfg

        if (ftype.eq.1)then ! radial
            !this type contains only the cutoff function
        elseif(ftype.eq.2)then ! radial (standard)
            eta_short_atomic_sfg(isf,idsfg,celem) = eta_short_atomic(csf,celem)
            rshift_short_atomic_sfg(isf,idsfg,celem) = rshift_short_atomic(csf,celem)
        elseif(ftype.eq.4)then ! radial
            eta_short_atomic_sfg(isf,idsfg,celem) = eta_short_atomic(csf,celem)
        elseif(ftype.eq.8)then ! angular(shifted)
            eta_short_atomic_sfg(isf,idsfg,celem) = eta_short_atomic(csf,celem)
            rshift_short_atomic_sfg(isf,idsfg,celem) = rshift_short_atomic(csf,celem) ! this is actually theta_shift
        elseif(ftype.eq.3.or.ftype.eq.9)then ! angulars (wide and narrow)
            zeta_short_atomic_sfg(isf,idsfg,celem) = zeta_short_atomic(csf,celem)
            lambda_short_atomic_sfg(isf,idsfg,celem) = lambda_short_atomic(csf,celem)
            eta_short_atomic_sfg(isf,idsfg,celem) = eta_short_atomic(csf,celem)
        elseif(ftype.eq.5.or.ftype.eq.6)then ! assert error as these are not implemented
            write(ounit,*)'This symmetry function type is not available for group-based calculation. Please &
                                switch of use_sf_groups!'
            stop
        end if

    end subroutine add_member


    subroutine check_neighborsandcutoff(type_local,neigh1,neigh2,ielement,&
            fc_type,sf_type,group_ID,num_group,lcreate)

        use fileunits

        implicit none

        integer type_local ! in
        integer num_group ! in
        integer group_ID ! in/out
        integer sf_type ! in
        integer ielement ! in
        integer i ! internal
        integer neigh1 ! in
        integer neigh2 ! in
        logical lcreate ! in/out
        real*8  fc_type ! in


        if (type_local.eq.0)then ! radial SF
            do i = 1,num_group
                if (neigh1.ne.neighbors_sfg(i,1,ielement).or.&
                        fc_type.ne.funccutoff_short_atomic_sfg(i,ielement).or.&
                        sf_type.ne.function_type_short_atomic_sfg(i,ielement))then
                    cycle
                else
                    group_ID = i
                    exit
                end if
            end do
            if (i.eq.(num_group+1)) then
                lcreate = .true.
                group_ID = num_group + 1
            end if
        else ! angular SF
            do i = 1,num_group
                if (neigh1.ne.neighbors_sfg(i,1,ielement).or.&
                        neigh2.ne.neighbors_sfg(i,2,ielement).or.&
                        fc_type.ne.funccutoff_short_atomic_sfg(i,ielement).or.&
                        sf_type.ne.function_type_short_atomic_sfg(i,ielement))then
                    cycle
                else
                    group_ID = i
                    exit
                end if
            end do
            if (i.eq.(num_group+1)) then
                lcreate = .true.
                group_ID = num_group + 1
            end if
        end if

    end subroutine check_neighborsandcutoff

    subroutine writesfgroups()

        use fileunits
        use globaloptions
        use nnflags

        implicit none

        integer i1 ! in
        integer i2 ! in
        integer i3 ! in

    !! Write symmetry function groups to output

    if(lshort.and.(nn_type_short.eq.1))then
        do i1=1,nelem
            write(ounit,*)'-------------------------------------------------------------'
            write(ounit,*)' short range atomic symmetry &
                    &function groups element ',element(i1),' :'
            write(ounit,*)'-------------------------------------------------------------'
            write(ounit,*)'  ID   Type   Cutoff  Neigh1  Neigh2'
            write(ounit,*)'-------------------------------------------------------------'
            do i2=1,num_groups_short_atomic(i1)
                if(function_type_short_atomic_sfg(i2,i1).eq.2)then
                    write(ounit,'(i5,3x,i3,3x,f8.3,3x,a3)')i2,function_type_short_atomic_sfg(i2,i1),&
                            funccutoff_short_atomic_sfg(i2,i1),element(elementindex(&
                            neighbors_sfg(i2,1,i1)))
                    write(ounit,*)'-------------------------------------------------------------'
                    write(ounit,*)'Associated SFs:  ID      eta     rshift'
                    do i3=1,sf_count_sfg(i2,i1)
                        write(ounit,'(17x,i3,3x,f8.3,3f8.3)')i3,eta_short_atomic_sfg(i3,i2,i1)&
                        ,rshift_short_atomic_sfg(i3,i2,i1)
                    end do
                    write(ounit,*)'-------------------------------------------------------------'

                else if(function_type_short_atomic_sfg(i2,i1).eq.1)then
                    write(ounit,'(i5,3x,i3,3x,f8.3,3x,a3)')i2,function_type_short_atomic_sfg(i2,i1),&
                            funccutoff_short_atomic_sfg(i2,i1),element(elementindex(&
                            neighbors_sfg(i2,1,i1)))
                    write(ounit,*)'-------------------------------------------------------------'
                    write(ounit,*)'Associated SFs:  ID  '
                    do i3=1,sf_count_sfg(i2,i1)
                        write(ounit,'(17x,i3,3x)')i3
                    end do
                    write(ounit,*)'-------------------------------------------------------------'

                elseif(function_type_short_atomic_sfg(i2,i1).eq.4)then
                    write(ounit,'(i5,3x,i3,3x,f8.3,3x,a3)')i2,function_type_short_atomic_sfg(i2,i1),&
                            funccutoff_short_atomic_sfg(i2,i1),element(elementindex(&
                            neighbors_sfg(i2,1,i1)))
                    write(ounit,*)'-------------------------------------------------------------'
                    write(ounit,*)'Associated SFs:  ID      eta     '
                    do i3=1,sf_count_sfg(i2,i1)
                        write(ounit,'(17x,i3,3x,f8.3)')i3,eta_short_atomic_sfg(i3,i2,i1)
                    end do
                    write(ounit,*)'-------------------------------------------------------------'

                elseif(function_type_short_atomic_sfg(i2,i1).eq.8)then
                    write(ounit,'(i5,3x,i3,3x,f8.3,3x,a3,6x,a3)')i2,function_type_short_atomic_sfg(i2,i1),&
                            funccutoff_short_atomic_sfg(i2,i1),element(elementindex(&
                            neighbors_sfg(i2,1,i1))),element(elementindex(neighbors_sfg(i2,2,i1)))
                    write(ounit,*)'-------------------------------------------------------------'
                    write(ounit,*)'Associated SFs:  ID      eta     thetashift'
                    do i3=1,sf_count_sfg(i2,i1)
                        write(ounit,'(17x,i3,3x,f8.3,3f8.3)')i3,eta_short_atomic_sfg(i3,i2,i1)&
                                ,rshift_short_atomic_sfg(i3,i2,i1)
                    end do
                    write(ounit,*)'-------------------------------------------------------------'


                else if(function_type_short_atomic_sfg(i2,i1).eq.3.or.&
                    function_type_short_atomic_sfg(i2,i1).eq.9)then
                    write(ounit,'(i5,3x,i3,3x,f8.3,3x,a3,6x,a3)')i2,function_type_short_atomic_sfg(i2,i1),&
                            funccutoff_short_atomic_sfg(i2,i1),element(elementindex(&
                            neighbors_sfg(i2,1,i1))),element(elementindex(neighbors_sfg(i2,2,i1)))
                    write(ounit,*)'-------------------------------------------------------------'
                    write(ounit,*)'Associated SFs:  ID      eta     zeta    lambda'
                    do i3=1,sf_count_sfg(i2,i1)
                        write(ounit,'(17x,i3,3x,f8.3,3f8.3)')i3,eta_short_atomic_sfg(i3,i2,i1)&
                                ,zeta_short_atomic_sfg(i3,i2,i1),lambda_short_atomic_sfg(i3,i2,i1)
                    end do
                    write(ounit,*)'-------------------------------------------------------------'
                end if

            enddo ! i2
        enddo ! i1=1,nelem
    endif ! lshort


    end subroutine writesfgroups

end module symfunctiongroups
