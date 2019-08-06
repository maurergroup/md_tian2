!############################################################################
! This routine is part of
! md_tian2 (Molecular Dynamics Xia Tian 2)
! (c) 2014-2019 Dan J. Auerbach, Sascha Kandratsenka, Svenja M. Janke, Marvin
! Kammler, Sebastian Wille
! Dynamics at Surfaces Department
! MPI for Biophysical Chemistry Goettingen, Germany
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
!############################################################################

module md_algo

    use universe_mod
    use constants
    use run_config
    use useful_things

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

                case (prop_id_langevin)
                    call langevin_1(atoms, i)

                case default
                    stop "Error in propagate_1(): Unknown propagation algorithm"

            end select
        end do

        if (any(atoms%algo == prop_id_andersen) .or. any(atoms%algo == prop_id_pile)) then
            call remove_com_velocity(atoms)
        end if

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

                case (prop_id_langevin)
                    call ldfa(atoms, i)
                    call langevin_2(atoms, i)

                case default
                    stop "Error in propagate_2(): Unknown propagation algorithm"

            end select
        end do

!        if (any(atoms%algo == prop_id_andersen) .or. any(atoms%algo == prop_id_pile)) then
!            call remove_com_velocity(atoms)
!        end if

    end subroutine propagate_2



    subroutine verlet_1(atoms, i)

        type(universe), intent(inout) :: atoms
        integer, intent(in) :: i

        where(.not. atoms%is_fixed(:,:,i))
            atoms%v(:,:,i) = atoms%v(:,:,i) + 0.5 * simparams%step * atoms%a(:,:,i)
        elsewhere
            atoms%v(:,:,i) = 0.0_dp
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
            atoms%v(:,:,i) = atoms%v(:,:,i) + 0.5 * simparams%step * atoms%a(:,:,i)
        elsewhere
            atoms%v(:,:,i) = 0.0_dp
        end where

    end subroutine verlet_2



    subroutine langevin_1(atoms, i)
        !
        ! Purpose:
        !           1st step of Langevin Dynamics algorithm,
        !           Allen & Tildesley, Computer Simulation of Liquids (1987), page 261.
        !           Li & Wahnström, Phys. Rev. B (1992).
        !

        use pes_emt_mod, only : dens

        type(universe), intent(inout) :: atoms
        integer, intent(in)           :: i

        real(dp)                             :: temp
        real(dp), dimension(atoms%nbeads)    :: c0, c1, c2, xidt, xidt2, ixidt, &
            sigma_r, sigma_v, c_rv
        real(dp), dimension(3, atoms%nbeads) :: randy
        integer :: b

        if (atoms%is_proj(atoms%idx(i))) then
            temp = kB * atoms%nbeads * simparams%Tproj / atoms%m(i)
        else
            temp = kB * atoms%nbeads * simparams%Tsurf / atoms%m(i)
        end if

        xidt = dens(:,i) * simparams%step
        xidt2 = xidt * xidt
        ixidt = simparams%step/xidt   ! 1/dens

        ! Preventing problems due to precision issues
        if (all(xidt > 1e-2) .and. temp > tolerance) then

            c0 = exp(-xidt)
            c1 = (1 - c0) * ixidt
            c2 = (1 - c1/simparams%step)*ixidt

            sigma_r = ixidt * sqrt(temp*(2*xidt - 3 + 4*c0 - c0*c0))
            sigma_v = sqrt(temp * (1 - c0*c0))

            c_rv    = ixidt*temp*(1 - c0)**2/(sigma_r*sigma_v)

        else ! use series up to 2nd order in xi*dt

            c0 = 1 - xidt + 0.5*xidt2
            c1 = (1 - 0.5*xidt + 2*twelfth*xidt2)*simparams%step
            c2 = (0.5 - 2*twelfth*xidt + 0.5*twelfth*xidt2)*simparams%step

            sigma_r = simparams%step * sqrt(temp*(8*twelfth*xidt - 0.5*xidt2))
            sigma_v = sqrt(2 * temp * xidt * (1 - xidt))

            c_rv    = 0.5 * sqrt3 * (1 - 0.125*xidt)

        end if

        call normal_deviate(0.0_dp, 1.0_dp, randy)

        ! no rpmd: propagate positions and partially propagate velocities
        if (atoms%nbeads == 1) then

            where (.not. atoms%is_fixed(:,1,i))
                atoms%r(:,1,i) = atoms%r(:,1,i) + c1(1)*atoms%v(:,1,i) + &
                    c2(1)*simparams%step*atoms%a(:,1,i) + sigma_r(1)*randy(:,1)
                atoms%v(:,1,i) = c0(1)*atoms%v(:,1,i) + &
                    (c1(1)-c2(1))*atoms%a(:,1,i) + sigma_v(1)*c_rv(1)*randy(:,1)
            elsewhere
                atoms%v(:,1,i) = 0.0_dp
            end where


        ! rpmd: partially propagate velocities
        else

            do b = 1, atoms%nbeads

                where (.not. atoms%is_fixed(:,b,i))
                    atoms%v(:,b,i) = c0(b)*atoms%v(:,b,i) + &
                    (c1(b)-c2(b))*atoms%a(:,b,i) + sigma_v(b)*c_rv(b)*randy(:,b)
                elsewhere
                    atoms%v(:,b,i) = 0.0_dp
                end where

            end do

        end if




    end subroutine langevin_1



    subroutine langevin_2(atoms, i)
            !
        ! Purpose:
        !           1st step of Langevin Dynamics algorithm,
        !           Dellago et al. JCP 108 (1998) 1964
        !

        use pes_emt_mod, only : dens

        type(universe), intent(inout) :: atoms
        integer, intent(in)           :: i

        real(dp)                             :: temp
        real(dp), dimension(atoms%nbeads)    :: c0, c1, c2, xidt, xidt2, ixidt, &
                                                    sigma_r, sigma_v, c_rv
        real(dp), dimension(3, atoms%nbeads) :: randy
        integer :: b

        if (atoms%is_proj(atoms%idx(i))) then
            temp = kB * atoms%nbeads * simparams%Tproj / atoms%m(i)
        else
            temp = kB * atoms%nbeads * simparams%Tsurf / atoms%m(i)
        end if

        xidt = dens(:,i) * simparams%step
        xidt2 = xidt * xidt
        ixidt = simparams%step / xidt   ! 1/dens

        ! Preventing problems due to precision issues
        if (all(xidt > 1e-2) .and. temp > tolerance) then

            c0 = exp(-xidt)
            c1 = (1 - c0) * ixidt
            c2 = (1 - c1/simparams%step)*ixidt

            sigma_r = ixidt * sqrt(temp*(2*xidt - 3 + 4*c0 - c0*c0))
            sigma_v = sqrt(temp * (1 - c0*c0))

            c_rv    = ixidt*temp*(1 - c0)**2/(sigma_r*sigma_v)

        else ! use series up to 2nd order in xi*dt

            c0 = 1 - xidt + 0.5*xidt2
            c1 = (1 - 0.5*xidt + 2*twelfth*xidt2)*simparams%step
            c2 = (0.5 - 2*twelfth*xidt + 0.5*twelfth*xidt2)*simparams%step

            sigma_v = sqrt(2 * temp * xidt * (1 - xidt))

            c_rv    = 0.5 * sqrt3 * (1 - 0.125*xidt)

        end if

        call normal_deviate(0.0_dp, 1.0_dp, randy)

        ! partially propagate velocities

        do b = 1, atoms%nbeads

            where (.not. atoms%is_fixed(:,b,i))
                atoms%v(:,b,i) = atoms%v(:,b,i) + c2(b)*atoms%a(:,b,i) + &
                    sigma_v(b)*sqrt(1-c_rv(b)*c_rv(b))*randy(:,b)
            elsewhere
                atoms%v(:,b,i) = 0.0_dp
            end where

        end do

    end subroutine langevin_2




    subroutine andersen(atoms, i)

        type(universe), intent(inout) :: atoms
        integer,        intent(in)    :: i

        real(dp), dimension(3, atoms%nbeads) :: new_v
        real(dp), dimension(atoms%nbeads)    :: choose
        real(dp) :: ibetaN, mass, andersen_threshold
        integer  :: b


        if (atoms%is_proj(atoms%idx(i))) then
            ibetaN = atoms%nbeads * kB * simparams%Tproj
        else
            ibetaN = atoms%nbeads * kB * simparams%Tsurf
        end if

        mass = atoms%m(i)

        call normal_deviate(0.0_dp, sqrt(ibetaN/mass), new_v)

        call random_number(choose)
        andersen_threshold = simparams%step / simparams%andersen_time

        do b = 1, atoms%nbeads
            if (choose(b) < andersen_threshold .and. .not. atoms%is_fixed(1,b,i)) then
                atoms%v(:,b,i) = new_v(:,b)
            end if
        end do


    end subroutine andersen


    subroutine pile_thermo(atoms, i)

        use rpmd, only : cjk, build_cjk

        type(universe), intent(inout) :: atoms
        integer       , intent(in)    :: i

        integer :: b, k
        real(8) :: wk, wn, betaN
        real(8), dimension(atoms%nbeads)    :: c1, c2, gammak
        real(8), dimension(3, atoms%nbeads) :: zeta, newP, atomP

        if (.not. allocated(cjk)) call build_cjk(atoms%nbeads)

        if (atoms%is_proj(atoms%idx(i))) then
            betaN = 1.0_dp / (kB * simparams%Tproj * atoms%nbeads)
        else
            betaN = 1.0_dp / (kB * simparams%Tsurf * atoms%nbeads)
        end if

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
                gammak(k+1) = 1 / simparams%pile_tau
            else
                wn = 1 / betaN / hbar
                wk = 2 * wn * sin(k*pi/atoms%nbeads)
                gammak(k+1) = 2 * wk
            end if
        end do

        c1 = exp(-0.5 * simparams%step*gammak)
        c2 = sqrt(1 - c1*c1)

        ! generate random number with zero mean and unit stddev
        call normal_deviate(0.0_dp, 1.0_dp, zeta)

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



    subroutine ldfa(atoms, i)
        !
        ! Purpose:
        !           Calculate the friction coefficient
        !

        use pes_emt_mod, only : dens

        type(universe), intent(in) :: atoms
        integer, intent(in) :: i
        integer :: j, b
        real(dp) :: fric(atoms%nbeads)
        real(dp) :: temp
        ! As implemented here, the friction coefficient is only applicable for the
        ! H-atom as calculated by Li and Wahnstrom(PRB46(1992)14528)
        ! according to Puska and Nieminen (PRB, 27, 1983, 6121), the mass still
        ! needs to be applied
        ! hbar*eta = hbar**2/mass Q(kf) (conversion between PN and LW)
        character(len=*), parameter :: err = "Error in ldfa(): "
        real(dp), parameter :: convert   = 1.00794_dp * amu2mass
        real(dp), parameter :: coefs(12) = [0.0802484_dp, -1.12851_dp, 9.28508_dp,   &
            2.10064_dp, -843.419_dp, 8.85354e3_dp, -4.89023e4_dp, 1.6741e5_dp, &
            -3.67098e5_dp, 5.03476e5_dp, -3.9426e5_dp, 1.34763e5_dp]

        !   12th order cubic spline fit interpolated from DFT data points of friction
        !   coefficient vs. electron density (calculated from DFT with VASP)

        if (.not. allocated(dens)) then
            print *, err, "EMT density array not allocated. Cannot compute friction."
            stop
        end if

        !print *, "pre", dens(:,i)
        fric = dens(:,i)
        ! hbar*xi in eV
        do b = 1, atoms%nbeads
            if (-1e-12 <= fric(b) .and. fric(b) <= 0.36) then  ! removed offset parameter
                fric(b) = 0.0_dp
                temp = dens(b,i)
                do j = 1, size(coefs)
                    fric(b) = fric(b) + coefs(j)*temp
                    temp = temp*dens(b,i)
                end do
            else if (fric(b) > 0.36) then
                fric(b) = 0.001_dp * (4.7131_dp - exp(-4.41305_dp*fric(b)))
            else
                print *, fric(b), dens(b,i), i, shape(dens)
                print *, err, "EMT density array contains negative values."
                stop
            end if
            dens(b,i) = fric(b)
        end do
        dens(:,i) = dens(:,i) * convert / hbar / atoms%m(i) / atoms%nbeads
        !print *, "post", dens(:,i)
        ! xi in 1/fs

        ! For simulated annealing, the Langevin dynamics are used as a heat bath
        ! But using the Au-atomic densities is too inefficient, so in this case
        ! a friction coefficient is set that is of the order of magnitude of
        ! the friction an H-atom experiences when running through Au.
        ! The friction coefficient is now 0.003 1/fs.

        ! if (sasteps > 0) s%dens = 0.003d0 !0.000015231d0/(imass*convert)

    end subroutine ldfa

end module md_algo
