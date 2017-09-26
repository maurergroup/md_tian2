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
                    call verlet_1(atoms, i)
                    call andersen(atoms, i)

                case (prop_id_pile)
                    call verlet_1(atoms, i)
                    call pile_thermo(atoms, i)

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

                case (prop_id_andersen)
                    call verlet_2(atoms, i)

                case (prop_id_pile)
                    call verlet_2(atoms, i)

                case default
                    stop "Error in propagate_2(): Unknown propagation algorithm"

            end select
        end do

    end subroutine propagate_2



    subroutine verlet_1(atoms, i)

        type(universe), intent(inout) :: atoms
        integer, intent(in) :: i

        where(.not. atoms%is_fixed(:,:,i))
            atoms%v(:,:,i) = atoms%v(:,:,i) + 0.5_dp * simparams%step * atoms%a(:,:,i)
        end where

        ! if rpmd, the positions are being updated in the do_ring_polymer_step subroutine
        if (atoms%nbeads == 1) then
            where(.not. atoms%is_fixed(:,:,i))
                atoms%r(:,:,i) = atoms%r(:,:,i) + simparams%step * atoms%v(:,:,i)
            end where
        end if

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
        real(dp) :: ibetaN, mass
        integer  :: k, b

        ibetaN = atoms%nbeads * kB * simparams%Tsurf

        call random_number(rnd1)
        call random_number(rnd2)
        rnd3 = sqrt(-2.0_dp*log(rnd1)) * cos(2.0_dp*pi*rnd2)

        call random_number(choose)
        mass = atoms%m(i)

        do b = 1, atoms%nbeads
            do k = 1, 3
                if (choose(k,b) < simparams%andersen_freq .and. .not. atoms%is_fixed(k,b,i)) then
                    atoms%v(k,b,i) = rnd3(k,b) * sqrt(ibetaN/mass)
                end if
            end do
        end do

    end subroutine andersen


    subroutine pile_thermo(atoms, i)

        use rpmd, only : cjk, build_cjk

        type(universe), intent(inout) :: atoms
        integer       , intent(in)    :: i

        integer :: b, k
        real(8) :: wk, wn, betaN
        real(8), dimension(atoms%nbeads)    :: c1, c2, gammak
        real(8), dimension(3, atoms%nbeads) :: rnd1, rnd2, zeta, newP, atomP

        if (.not. allocated(cjk)) call build_cjk(atoms%nbeads)

        betaN = 1.0_dp / (kB * simparams%Tsurf * atoms%nbeads)

        ! Transform to normal mode space
        newP = 0.0_dp
        atomP = calc_momentum_one(atoms, i)

        do b = 1, atoms%nbeads
            do k = 1, atoms%nBeads
                newP(:,b) = newP(:,b) + atomP(:,k)*cjk(k,b)
            end do
        end do

        ! generate gamma coefficients for all beads
        do k = 0, atoms%nbeads-1
            if (k .eq. 0) then  ! centroid mode
                gammak(k+1) = 1.0_dp/simparams%pile_tau
            else
                wn = 1.0_dp / betaN / hbar
                wk = 2.0_dp * wn * sin(k*pi/atoms%nbeads)
                gammak(k+1) = 2.0_dp * wk
            end if
        end do

        c1 = exp(-0.5_dp * simparams%step*gammak)
        c2 = sqrt(1.0_dp - c1*c1)

        ! generate random number with zero mean and unit stddev
        call random_number(rnd1)
        call random_number(rnd2)
        zeta = sqrt(-2.0_dp*log(rnd1)) * cos(2.0_dp*pi*rnd2)

        do b = 1, atoms%nbeads
            newP(:,b) = c1(b)*newP(:,b) + sqrt(atoms%m(i)/betaN)*c2(b)*zeta(:,b)
        end do

        ! Transform back to Cartesian space
        atomP = 0.0_dp
        do b = 1, atoms%nbeads
            do k = 1, atoms%nbeads
                where (.not. atoms%is_fixed(:,b,i))
                     atomP(:,b) = atomP(:,b) + newP(:,k)*cjk(b,k)
                end where
            end do
        end do
        atoms%v(:,:,i) = atomP/atoms%m(i)

    end subroutine pile_thermo

end module md_algo
