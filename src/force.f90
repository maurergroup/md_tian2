module force

    use constants
    use universe_mod, only : universe
    use pes_lj_mod, only : compute_lj

    implicit none

contains

    subroutine calc_force(atoms)

        type(universe), intent(inout) :: atoms
        integer :: i, j
        character(len=*), parameter :: err = "Error in calc_force(): "

        do i = 1, atoms%natoms-1
            do j = i+1, atoms%natoms

                select case (atoms%pes( atoms%idx(i),atoms%idx(j) ))
                    case (pes_id_lj)
                        call compute_lj(atoms, i, j, energy_and_force)

                    case (default_int)
                        print *, err // "pes not specified for atoms", i, j
                        stop

                    case default
                        print *, err // "unknown force field:", atoms%pes(atoms%idx(i),atoms%idx(j))
                        stop

                end select

            end do
        end do

    end subroutine calc_force


end module force
