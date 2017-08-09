
module pes_lj_module

    use constants
    use useful_things, only : split_string, pbc_distdir
    use atom_class, only : atoms, get_idx_from_name

    implicit none


    type lj_pes

        private
        real(dp), allocatable :: eps(:,:)
        real(dp), allocatable :: sigma(:,:)
        real(dp), allocatable :: shift(:,:)
        real(dp) ::  cutoff


    end type lj_pes

    type(lj_pes) :: pes_lj



contains

    subroutine read_lj(teil, slab, inp_unit)

        use run_config, only : simparams
        use atom_class, only : simbox

        type(atoms), intent(in) :: teil, slab
        integer, intent(in) :: inp_unit

        integer :: nwords = 0, ios = 0, i
        character(len=max_string_length) :: buffer
        character(len=max_string_length) :: words(100)
        integer  :: idx1, idx2, ntypes
        real(dp) :: temp3(3)
        character(len=*), parameter :: err_read_lj = "Error in read_lj: "

        ntypes = simparams%nprojectiles+simparams%nlattices

        if (.not. allocated(pes_lj%eps)) then
            allocate(pes_lj%eps(ntypes, ntypes),   &
                     pes_lj%sigma(ntypes, ntypes), &
                     pes_lj%shift(ntypes, ntypes))
            pes_lj%eps   = default_real
            pes_lj%sigma = default_real
            pes_lj%shift = default_real
            do i = 1, 3
                temp3(i) = sqrt(sum(simbox(:,i)*simbox(:,i)))
            end do
            pes_lj%cutoff = minval(temp3(:))
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

            ! pes block terminated, set shift
            if (nwords == 0 .or. ios /= 0) then
                pes_lj%shift(idx1,idx2) = 4.0_dp * pes_lj%eps(idx1,idx2) * &
                                           ( (pes_lj%sigma(idx1,idx2)/pes_lj%cutoff)**12 &
                                           - (pes_lj%sigma(idx1,idx2)/pes_lj%cutoff)**6 )
                pes_lj%shift(idx2,idx1) = pes_lj%shift(idx1,idx2)
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


    subroutine compute_lj(teil, slab, flag)

        type(atoms), intent(inout) :: teil, slab
        integer, intent(in)        :: flag

        integer :: i, j, b, idx_i, idx_j
        real(dp) :: vec(3)
        real(dp) :: r, r2, fr2, fr6, fpr
        real(dp) :: f(3)

        ! teil-teil interaction
        do i = 1, teil%natoms
            do j = i+1, teil%natoms
                do b = 1, teil%nbeads

                    idx_i = teil%idx(i)
                    idx_j = teil%idx(j)

                    call pbc_distdir(teil%r(:,b,i), teil%r(:,b,j), r, vec)

                    r2 = r*r
                    fr2 = pes_lj%sigma( idx_i , idx_j ) / r2
                    fr6 = fr2 * fr2 * fr2
                    fpr = 48._dp * pes_lj%eps( idx_i , idx_j ) * fr6 * (fr6 - 0.5_dp) / r2   ! f/r

                    f(:) = fpr * vec(:)

                    teil%f(:,b,i) = teil%f(:,b,i) + f(:)
                    teil%f(:,b,j) = teil%f(:,b,j) - f(:)



                end do
            end do
        end do



    end subroutine compute_lj



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
