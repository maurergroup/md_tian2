
module pes_lj_module

    use constants
    use useful_things, only : split_string
    use atom_class

    implicit none


    type lj_pes

        private
        real(dp), allocatable :: eps(:,:)
        real(dp), allocatable :: sigma(:,:)

    end type lj_pes

    type(lj_pes) :: pes_lj



contains

    subroutine read_lj(teil, slab, inp_unit)

        use run_config, only : simparams

        type(atoms), intent(in) :: teil, slab
        integer, intent(in) :: inp_unit

        integer :: nwords = 0
        integer :: ios = 0
        character(len=max_string_length) :: buffer
        character(len=max_string_length) :: words(100)
        integer :: idx1, idx2, ntypes
        character(len=*), parameter :: err_read_lj = "Error in read_lj: "

        ntypes = simparams%nprojectiles+simparams%nlattices

        if (.not. allocated(pes_lj%eps)) then
            allocate(pes_lj%eps(ntypes, ntypes), pes_lj%sigma(ntypes, ntypes))
            pes_lj%eps   = default_real
            pes_lj%sigma = default_real
        end if

        ! line should read something like "H   H   proj    proj"
        read(inp_unit, '(A)', iostat=ios) buffer
        call split_string(buffer, words, nwords)

        if (nwords /= 4) stop err_read_lj // "need four entries in interaction-defining lines"

        if (words(3) == "proj" .and. words(4) == "proj") then
            idx1 = get_idx_from_name(teil, words(1))
            idx2 = get_idx_from_name(teil, words(2))

        else if (words(3) == "proj" .and. words(4) == "latt") then
            idx1 = get_idx_from_name(teil, words(1))
            idx2 = get_idx_from_name(slab, words(2))

        else if (words(3) == "latt" .and. words(4) == "proj") then
            idx1 = get_idx_from_name(slab, words(1))
            idx2 = get_idx_from_name(teil, words(2))

        else if (words(3) == "latt" .and. words(4) == "latt") then
            idx1 = get_idx_from_name(slab, words(1))
            idx2 = get_idx_from_name(slab, words(2))

        else
            stop err_read_lj // "interaction must be defined via 'proj' and 'latt' keywords"

        end if

        do
            read(inp_unit, '(A)', iostat=ios) buffer
            call split_string(buffer, words, nwords)

            ! pes block terminated
            if (nwords == 0 .or. ios /= 0) then
                exit

            ! something went wrong
            else if (nwords /= 2) then
                stop "Error in the PES file: PES parameters must consist of key value pairs. A parameter block must be terminated by a blank line."
            end if

            select case (words(1))
                case ('sigma')
                    read(words(2), *) pes_lj%sigma(idx1, idx2)
                    read(words(2), *) pes_lj%sigma(idx2, idx1)

                case ('epsilon')
                    read(words(2), *) pes_lj%eps(idx1, idx2)
                    read(words(2), *) pes_lj%eps(idx2, idx1)

                case default
                    print *, "Error in the PES file: unknown LJ parameter", words(1)
                    stop

            end select
        end do

    end subroutine read_lj


!    subroutine compute_lj(teil, slab, flag)
!
!        type(atoms), intent(inout) :: teil, slab
!        integer, intent(in)        :: flag
!
!        integer :: i, j, b
!        real(dp) :: xr(3),
!
!        ! teil-teil interaction
!        do i = 1, teil%natoms
!            do j = i+1, teil%natoms
!                do b = 1, teil%nbeads
!
!
!                    xr = xr-box*nint (xr/box)
!                    r2 = xr**2
!                    if (r2.1t.rc2) then
!                        r2i=i/r2
!                        r6i=r2i**3
!                        ff=48*r2i*r6i* (r6i-0.5)
!                        f (i)=f (i) +ff*xr
!                        f (j)=f (j) -ff*xr
!                        en=en+4*r6i* (r6i-l) -ecut
!                    endif
!
!
!                end do
!            end do
!        end do
!
!
!
!    end subroutine compute_lj



end module pes_lj_module




!pes lj
!H   H   proj    proj
!sigma   2.4
!epsilon 3.0
!
!pes lj
!H   D   proj    proj
!sigma   2.4
!epsilon 3.0
!sigma   2.4
!epsilon 3.0
!
!pes lj
!D   D   proj    proj
!sigma   2.4
!epsilon 3.0
