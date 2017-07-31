
module pes_lj

    use constants
    use useful_things, only : split_string

    implicit none

    type lj_pes

        private
        real(dp), allocatable :: eps(:,:)
        real(dp), allocatable :: sigma(:,:)

    end type lj_pes

    type(lj_pes) :: pes



contains

    subroutine initialize(inp_unit)

        use run_config, only : get_atomic_indices, simparams

        integer, intent(in) :: inp_unit

        integer :: nwords = 0
        integer :: ios = 0
        character(len=max_string_length) :: buffer
        character(len=max_string_length) :: words(100)
        integer :: indices(2)

        if (.not. allocated(pes%eps)) allocate(pes%eps(simparams%nprojectiles+simparams%nlattices,    &
                                                       simparams%nprojectiles+simparams%nlattices),   &
                                               pes%sigma(simparams%nprojectiles+simparams%nlattices,  &
                                                         simparams%nprojectiles+simparams%nlattices))

        ! line should real something like "H   H   proj    proj"
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
            if (nwords == 0) then
                exit

            ! something went wrong
            else if (nwords /= 2) then
                stop "Error in the PES file: PES parameters must consist of key value pairs."
            end if

            select case (words(1))
            case ('sigma')
                read(words(2), *) pes%sigma(indices(1),indices(2))

            case ('epsilon')
                read(words(2), *) pes%eps(indices(1),indices(2))

            case default
                print *, "Error in the PES file: unknown LJ parameter", words(1)
                stop

            end select

        end do



    end subroutine initialize

end module pes_lj




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
