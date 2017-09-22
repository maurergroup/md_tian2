module force

    use constants
    use universe_mod, only : universe

    use pes_lj_mod,  only : compute_lj, compute_simple_lj
    use pes_emt_mod, only : compute_emt
    use pes_ho_mod,  only : compute_ho

    implicit none

contains

    subroutine calc_force(atoms, flag)

        type(universe), intent(inout) :: atoms
        integer, intent(in)           :: flag
        character(len=*), parameter   :: err = "Error in calc_force(): "

        atoms%f    = 0.0_dp
        atoms%a    = 0.0_dp
        atoms%epot = 0.0_dp

        if (any(atoms%pes == pes_id_lj))        call compute_lj       (atoms, flag)
        if (any(atoms%pes == pes_id_simple_lj)) call compute_simple_lj(atoms, flag)
        if (any(atoms%pes == pes_id_emt))       call compute_emt      (atoms, flag)
        if (any(atoms%pes == pes_id_ho))        call compute_ho       (atoms, flag)

        call set_acceleration(atoms)

    end subroutine calc_force


    subroutine set_acceleration(atoms)

        type(universe), intent(inout) :: atoms
        integer :: i

        do i = 1, atoms%natoms
            where(.not. atoms%is_fixed(:,:,i)) atoms%a(:,:,i) = atoms%f(:,:,i)/atoms%m(i)
        end do

    end subroutine set_acceleration


end module force
