module md_algo

    use universe_mod
    use constants
    use run_config

    implicit none



contains

    subroutine propagate_1(atoms)

        type(universe), intent(inout) :: atoms
        integer :: i

        do i = 1, atoms%natoms
            select case(atoms%algo(atoms%idx(i)))
                case (prop_id_verlet)
                    call verlet_1(atoms, i)

                case default
                    stop "Error in propagate_1(): Unknown propagation algorithm"

            end select
        end do

    end subroutine propagate_1



    subroutine propagate_2(atoms)

        type(universe), intent(inout) :: atoms
        integer :: i

        do i = 1, atoms%natoms
            select case(atoms%algo(atoms%idx(i)))
                case (prop_id_verlet)
                    call set_acceleration(atoms)
                    call verlet_2(atoms, i)

                case default
                    stop "Error in propagate_2(): Unknown propagation algorithm"

            end select
        end do

    end subroutine propagate_2




    subroutine set_acceleration(atoms)

        type(universe), intent(inout) :: atoms
        integer :: i

        forall (i = 1 : atoms%natoms)
            where(.not. atoms%is_fixed(:,:,i)) atoms%a(:,:,i) = atoms%f(:,:,i)/atoms%m(atoms%idx(i))
        end forall

    end subroutine set_acceleration



    subroutine verlet_1(atoms, i)

        type(universe), intent(inout) :: atoms
        integer, intent(in) :: i

        where(.not. atoms%is_fixed(:,:,i))
            atoms%v(:,:,i) = atoms%v(:,:,i) + 0.5_dp * simparams%step * atoms%a(:,:,i)
            atoms%r(:,:,i) = atoms%r(:,:,i) +          simparams%step * atoms%v(:,:,i)
        end where

    end subroutine verlet_1



    subroutine verlet_2(atoms, i)
        type(universe), intent(inout) :: atoms
        integer, intent(in) :: i

        where(.not. atoms%is_fixed(:,:,i)) atoms%v(:,:,i) = atoms%v(:,:,i) + 0.5_dp * simparams%step * atoms%a(:,:,i)

    end subroutine verlet_2

end module md_algo
