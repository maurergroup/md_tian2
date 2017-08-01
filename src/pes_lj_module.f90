
module pes_lj_module

    use constants
    use useful_things, only : split_string

    implicit none


    type lj_pes

        private
        real(dp), allocatable :: eps(:,:)
        real(dp), allocatable :: sigma(:,:)

    end type lj_pes

    type(lj_pes) :: pes_lj



contains

    subroutine read_lj(inp_unit)

        use run_config, only : get_atomic_indices, simparams

        integer, intent(in) :: inp_unit

        integer :: nwords = 0
        integer :: ios = 0
        character(len=max_string_length) :: buffer
        character(len=max_string_length) :: words(100)
        integer :: indices(2), ntypes

        ntypes = simparams%nprojectiles+simparams%nlattices

        if (.not. allocated(pes_lj%eps)) then
            allocate(pes_lj%eps(ntypes, ntypes), pes_lj%sigma(ntypes, ntypes))
            pes_lj%eps   = default_real
            pes_lj%sigma = default_real
        end if

        ! line should read something like "H   H   proj    proj"
        read(inp_unit, '(A)', iostat=ios) buffer
        call split_string(buffer, words, nwords)

        if (nwords == 4) then
            indices = get_atomic_indices(words(1:4))
        else
            stop "Error in PES file. Need four entries in interaction-defining lines"
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
                read(words(2), *) pes_lj%sigma(indices(1),indices(2))
                read(words(2), *) pes_lj%sigma(indices(2),indices(1))

            case ('epsilon')
                read(words(2), *) pes_lj%eps(indices(1),indices(2))
                read(words(2), *) pes_lj%eps(indices(2),indices(1))

            case default
                print *, "Error in the PES file: unknown LJ parameter", words(1)
                stop

            end select
        end do

    end subroutine read_lj



    subroutine check_lj(interacting)

        use run_config, only : simparams

        logical :: interacting(:,:)

        interacting = .false.

        where (pes_lj%eps /= default_real) interacting = .true.


    end subroutine check_lj

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
