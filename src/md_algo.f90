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

                case (prop_id_andersen)
                    call andersen(atoms, i)

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
                    call verlet_2(atoms, i)

                case default
                    stop "Error in propagate_2(): Unknown propagation algorithm"

            end select
        end do

    end subroutine propagate_2



    subroutine verlet_1(atoms, i)

        type(universe), intent(inout) :: atoms
        integer, intent(in) :: i
!        print *, "in verlet_1"
!        print *, atoms%r
!        print *, "end verlet_1"

        where(.not. atoms%is_fixed(:,:,i))
            atoms%v(:,:,i) = atoms%v(:,:,i) + 0.5_dp * simparams%step * atoms%a(:,:,i)
            atoms%r(:,:,i) = atoms%r(:,:,i) +          simparams%step * atoms%v(:,:,i)
        end where

    end subroutine verlet_1



    subroutine verlet_2(atoms, i)
        type(universe), intent(inout) :: atoms
        integer, intent(in) :: i

        where(.not. atoms%is_fixed(:,:,i))
            atoms%v(:,:,i) = atoms%v(:,:,i) + 0.5_dp * simparams%step * atoms%a(:,:,i)
        end where

    end subroutine verlet_2


    subroutine andersen(atoms, i)

        type(universe), intent(inout) :: atoms
        integer,        intent(in)    :: i

        real(dp), dimension(3, atoms%nbeads) :: rnd1, rnd2, rnd3, choose
        integer :: k, b

        call random_number(rnd1)
        call random_number(rnd2)
        rnd3 = sqrt(-2.0_dp*log(rnd1)) * cos(2.0_dp*pi*rnd2)

        call random_number(choose)

        do b = 1, atoms%nbeads
        do k = 1, 3
            if (choose(k,b) < 1.0_dp/50.0_dp .and. .not. atoms%is_fixed(k,b,i))
                atoms%v(k,b,i) = rnd3(k,b) * sqrt(atoms%m(atoms%idx(i))/Tsurf)
            end if
        end do
        end do



    end subroutine andersen

end module md_algo
