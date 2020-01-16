!############################################################################
! This routine is part of
! md_tian2 (Molecular Dynamics Tian Xia 2)
! (c) 2014-2020 Dan J. Auerbach, Svenja M. Janke, Marvin Kammler,
!               Sascha Kandratsenka, Sebastian Wille
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

module rpmd !ring polymer MD

    use universe_mod
    use run_config


    implicit none

    real(dp), allocatable :: cjk(:,:)

contains

    subroutine build_cjk(nbeads)
        integer, intent(in) :: nbeads

        integer :: j, k

        allocate(cjk(nbeads, nbeads))

        do j = 1, nbeads
            do k = 0, nbeads-1
                if (k == 0) then
                    cjk(j,k+1) = sqrt(1.0_dp/nbeads)
                else if (1 <= k .and. k <= nbeads/2 - 1) then
                    cjk(j,k+1) = sqrt(2.0_dp/nbeads) * cos(2.0_dp*pi*j*k/nbeads)
                else if (k == nbeads/2) then
                    cjk(j,k+1) = sqrt(1.0_dp/nbeads)*(-1.0_dp)**j
                else if (nbeads/2+1 <= k .and. k <= nbeads-1) then
                    cjk(j,k+1) = sqrt(2.0_dp/nbeads) * sin(2.0_dp*pi*j*k/nbeads)
                else
                    stop "Error in build_cjk()"
                end if
            end do
        end do

    end subroutine build_cjk




    real(dp) function calc_bead_temperature(atoms) result(temperature)

        type(universe), intent(in) :: atoms
        real(dp)                   :: temp
        integer                    :: i

        temp = 0.0_dp
        do i = 1, atoms%natoms
            temp = temp + atoms%m(i) * sum(atoms%v(:,:,i)*atoms%v(:,:,i))
        end do

        temperature = temp / kB / atoms%dof

    end function calc_bead_temperature




    function calc_centroid_positions(atoms) result(cents)

        type(universe), intent(in) :: atoms
        real(dp)                   :: cents(3, atoms%natoms)

        integer  :: i

        if (.not. atoms%is_cart) print *, "Warning: centroid positions are being calculated from direct coordinates!"

        ! active rpmd
        if (atoms%nbeads > 1) then

            do i = 1, atoms%natoms
                cents(:,i) = sum(atoms%r(:,:,i), dim=2)/atoms%nbeads
            end do

        ! no rpmd
        else
            cents = atoms%r(:,1,:)
        end if


    end function calc_centroid_positions




    function calc_centroid_velocities(atoms) result(cents)

        type(universe), intent(in) :: atoms
        real(dp)                   :: cents(3, atoms%natoms)

        integer  :: i

        if (.not. atoms%is_cart) print *, "Warning: centroid velocities are being calculated from direct coordinates!"

        ! active rpmd
        if (atoms%nbeads > 1) then

            do i = 1, atoms%natoms
                cents(:,i) = sum(atoms%v(:,:,i), dim=2)/atoms%nbeads
            end do

        ! no rpmd
        else
            cents = atoms%v(:,1,:)
        end if

    end function calc_centroid_velocities



    function calc_centroid_velocity_one(atoms, the_atom) result(cent_v)

        type(universe), intent(in) :: atoms
        integer       , intent(in) :: the_atom

        real(dp)                   :: cent_v(3)
        integer  :: i

        if (.not. atoms%is_cart) print *, "Warning: centroid velocities are being calculated from direct coordinates!"


        cent_v = sum(atoms%v(:,:,the_atom), dim=2)/atoms%nbeads

    end function calc_centroid_velocity_one




    function calc_inter_bead_distances(atoms) result (dists)

        type(universe), intent(in) :: atoms

        real(dp) :: dists(atoms%nbeads, atoms%natoms)
        real(dp) :: vec(3, atoms%nbeads, atoms%natoms)
        integer :: b, k

        ! no rpmd
        if (atoms%nbeads == 1) then
            dists = 0.0_dp

        ! active rpmd
        else
            do b = 1, atoms%nbeads
                k = modulo(b, atoms%nbeads)+1                   ! bead b+1
                vec(:,b,:) = atoms%r(:,b,:) - atoms%r(:,k,:)    ! distance vector
            end do
            dists = sqrt(sum(vec*vec, dim=1))    ! distance
        end if

    end function calc_inter_bead_distances




    real(dp) function calc_bead_epot(atoms) result(epot)

        type(universe), intent(in) :: atoms

        real(dp) :: dx(atoms%nbeads, atoms%natoms), wn2
        integer  :: i

        dx  = calc_inter_bead_distances(atoms)

        epot = 0.0_dp
        do i = 1, atoms%natoms

            if (atoms%is_proj(atoms%idx(i))) then
                wn2 = (kb * simparams%Tproj * atoms%nbeads / hbar)**2
            else
                wn2 = (kb * simparams%Tsurf * atoms%nbeads / hbar)**2
            end if

            epot = epot + atoms%m(i) * sum(dx(:,i) * dx(:,i)) * wn2
        end do

        epot = 0.5_dp * epot / atoms%nbeads

    end function calc_bead_epot




    subroutine virial_quantum_ekin(atoms, ekin_p, ekin_l)

        type(universe), intent(in)  :: atoms
        real(dp)      , intent(out) :: ekin_p, ekin_l

        real(dp) :: cents(3, atoms%natoms)
        integer  :: b, i

        ekin_p = 0.0_dp
        ekin_l = 0.0_dp

        if (atoms%nbeads > 1) then

            cents = calc_centroid_positions(atoms)

            do i = 1, atoms%natoms
                do b = 1, atoms%nbeads
                    if (atoms%is_proj(atoms%idx(i))) then
                        ekin_p = ekin_p + sum((cents(:,i)-atoms%r(:,b,i)) * atoms%f(:,b,i))
                    else
                        ekin_l = ekin_l + sum((cents(:,i)-atoms%r(:,b,i)) * atoms%f(:,b,i))
                    end if
                end do
            end do

            ekin_p = 0.5_dp * ekin_p / atoms%nbeads
            ekin_l = 0.5_dp * ekin_l / atoms%nbeads

        end if

    end subroutine virial_quantum_ekin




    subroutine do_ring_polymer_step(atoms)

        type(universe), intent(inout) :: atoms

        real(dp), dimension(3, atoms%nbeads, atoms%natoms)  :: p, q, newP, newQ

        real(dp) :: poly(4, atoms%nbeads), p_new(3)
        real(dp) :: twown, wk, wt, wm, cos_wt, sin_wt
        real(dp) :: betaN, piN, mass

        integer :: i, b, k

        if (.not. allocated(cjk)) call build_cjk(atoms%nbeads)

        ! Transform to normal mode space
        call calc_momentum_all(atoms, p)
        q = atoms%r

        !        do k = 1, 3
        !            do i = 1, atoms%natoms
        !                call rfft(p(k,:,i), atoms%nbeads)
        !                call rfft(q(k,:,i), atoms%nbeads)
        !            end do
        !        end do

        ! TODO: Intel MKL 3D-FFT should be faster

        newP = 0.0_dp
        newQ = 0.0_dp
        do b = 1, atoms%nbeads
            do k = 1, atoms%nBeads
                newP(:,b,:) = newP(:,b,:) + p(:,k,:)*cjk(k,b)
                newQ(:,b,:) = newQ(:,b,:) + q(:,k,:)*cjk(k,b)
            end do
        end do
        !        print *, "firstp", newP(:,1,1)
        p = newP
        q = newQ

        piN = pi / atoms%nbeads
        do i = 1, atoms%natoms
            mass = atoms%m(i)
            poly(1, 1) = 1.0_dp
            poly(2, 1) = 0.0_dp
            poly(3, 1) = simparams%step / mass
            poly(4, 1) = 1.0_dp

            if (atoms%is_proj(atoms%idx(i))) then
                betaN = 1.0_dp / (kB * simparams%Tproj * atoms%nbeads)
            else
                betaN = 1.0_dp / (kB * simparams%Tsurf * atoms%nbeads)
            end if

            if (atoms%nbeads > 1) then
                twown = 2.0_dp / betaN / hbar
                do b = 1, atoms%nbeads / 2
                    wk = twown * sin(b * piN)
                    wt = wk * simparams%step
                    wm = wk * mass
                    cos_wt = cos(wt)
                    sin_wt = sin(wt)
                    poly(1, b+1) =       cos_wt
                    poly(2, b+1) = -wm * sin_wt
                    poly(3, b+1) =       sin_wt / wm
                    poly(4, b+1) =       cos_wt
                end do
                do b = 1, (atoms%nbeads - 1) / 2
                    poly(1, atoms%nbeads-b+1) = poly(1, b+1)
                    poly(2, atoms%nbeads-b+1) = poly(2, b+1)
                    poly(3, atoms%nbeads-b+1) = poly(3, b+1)
                    poly(4, atoms%nbeads-b+1) = poly(4, b+1)
                end do
            end if

            do b = 1, atoms%nbeads
                p_new = p(:,b,i) * poly(1,b) + q(:,b,i) * poly(2,b)
                !                if (b == 1 .and. i == 1) then
                !                    print *, "secondp", p(:,b,i)
                !                    print *, "secondq", q(:,b,i)
                !                end if
                q(:,b,i) = p(:,b,i) * poly(3,b) + q(:,b,i) * poly(4,b)
                p(:,b,i) = p_new
            end do
        end do

        ! Transform back to Cartesian space
        where(.not. atoms%is_fixed)
            atoms%r = 0.0_dp
            atoms%v = 0.0_dp
        end where

        !print *, "lastp", p(:,1,1)

        do i = 1, atoms%natoms
            do b = 1, atoms%nbeads
                do k = 1, atoms%nbeads
                    where (.not. atoms%is_fixed(:,b,i))
                        atoms%r(:,b,i) = atoms%r(:,b,i) + q(:,k,i)*cjk(b,k)
                        atoms%v(:,b,i) = atoms%v(:,b,i) + p(:,k,i)*cjk(b,k)/atoms%m(i)
                    end where
                end do
            end do
           ! atoms%v(:,:,i) = atoms%v(:,:,i)
        end do

        !        do k = 1, 3
        !            do i = 1, atoms%natoms
        !                call irfft(p(k,:,i), atoms%nbeads)
        !                call irfft(q(k,:,i), atoms%nbeads)
        !            end do
        !        end do
        !
        !        do i = 1, atoms%natoms
        !            where(.not. atoms%is_fixed(:,:,i))
        !                atoms%r(:,:,i) = q(:,:,i)
        !                atoms%v(:,:,i) = p(:,:,i)/atoms%m(i)
        !            end where
        !        end do

    end subroutine do_ring_polymer_step




    subroutine do_ring_polymer_step_with_forces(atoms, forces)

        type(universe), intent(inout) :: atoms
        real(dp), intent(out)         :: forces(3, atoms%nbeads, atoms%natoms)

        real(dp), dimension(3,atoms%nbeads,atoms%natoms) :: momenta_init, momenta_final

        call calc_momentum_all(atoms, momenta_init)
        call do_ring_polymer_step(atoms)
        call calc_momentum_all(atoms, momenta_final)
        forces = (momenta_final - momenta_init) / simparams%step

    end subroutine do_ring_polymer_step_with_forces




    subroutine radius_of_gyration(atoms, rgyr_p, rgyr_l)

        type(universe), intent(in) :: atoms
        real(dp), intent(out)      :: rgyr_p, rgyr_l

        real(dp) :: cents(3, atoms%natoms)
        integer  :: i, b, nprojs

        rgyr_p = 0.0_dp
        rgyr_l = 0.0_dp

        if (atoms%nbeads == 1) return

        cents  = calc_centroid_positions(atoms)
        nprojs = 0

        do i = 1, atoms%natoms
            if (atoms%is_proj(atoms%idx(i))) then
                nprojs = nprojs + 1
                do b = 1, atoms%nbeads
                    rgyr_p = rgyr_p + sum( (atoms%r(:,b,i)-cents(:,i))**2 )
                end do
            else
                do b = 1, atoms%nbeads
                    rgyr_l = rgyr_l + sum( (atoms%r(:,b,i)-cents(:,i))**2 )
                end do
            end if
        end do

        rgyr_p = rgyr_p/atoms%nbeads/max(nprojs, 1)
        rgyr_l = rgyr_l/atoms%nbeads/max((atoms%natoms-nprojs), 1)

    end subroutine radius_of_gyration




    subroutine bead_ekin(atoms, ekin_p, ekin_l)

        type(universe), intent(in) :: atoms
        real(dp), intent(out)      :: ekin_p, ekin_l
        real(dp)                   :: nrg
        integer :: i

        ekin_p = 0.0_dp
        ekin_l = 0.0_dp

        do i = 1, atoms%natoms

            nrg = atoms%m(i)*sum(atoms%v(:,:,i)*atoms%v(:,:,i))

            if (atoms%is_proj(atoms%idx(i))) then
                ekin_p = ekin_p + nrg
            else
                ekin_l = ekin_l + nrg
            end if

        end do

        ekin_p = 0.5_dp * ekin_p / atoms%nbeads
        ekin_l = 0.5_dp * ekin_l / atoms%nbeads

    end subroutine bead_ekin




    ! Compute the real fast Fourier transform of the given array of data.
    ! Parameters:
    !   x - The array of data to transform
    !   N - The length of the array of data
    ! Returns:
    !   x - The transformed array of data, in half-complex form
!    subroutine rfft(x,N)
!        implicit none
!
!        integer, intent(in) :: N
!        double precision, intent(inout) :: x(N)
!
!        integer, parameter :: Nmax = 1024
!        integer :: Np
!        double precision :: copy(Nmax), factor
!        integer(8) :: plan
!
!
!        data Np /0/
!        save copy, factor, plan, Np
!
!        if (N .ne. Np) then
!            if (Np .ne. 0) call dfftw_destroy_plan(plan)
!            call dfftw_plan_r2r_1d(plan,N,copy,copy,0,64)
!            factor = sqrt(1.0d0/N)
!            Np = N
!        end if
!
!        copy(1:N) = x
!        call dfftw_execute(plan)
!        x = factor * copy(1:N)
!
!    end subroutine rfft

    ! Compute the inverse real fast Fourier transform of the given array of data.
    ! Parameters:
    !   x - The array of data to transform, in half-complex form
    !   N - The length of the array of data
    ! Returns:
    !   x - The transformed array of data
!    subroutine irfft(x,N)
!
!        implicit none
!        integer, intent(in) :: N
!        double precision, intent(inout) :: x(N)
!
!        integer, parameter :: Nmax = 1024
!        integer :: Np
!        double precision :: copy(Nmax), factor
!        integer(8) :: plan
!
!        data Np /0/
!        save copy, factor, plan, Np
!
!        if (N .ne. Np) then
!            ! The input array is a different length than the last array, so we
!            ! must generate a new FFTW plan for the transform
!            ! First delete the previous plan
!            if (Np .ne. 0) call dfftw_destroy_plan(plan)
!            call dfftw_plan_r2r_1d(plan,N,copy,copy,1,64)
!            factor = sqrt(1.0d0/N)
!            Np = N
!        end if
!
!        copy(1:N) = x
!        call dfftw_execute(plan)
!        x = factor * copy(1:N)
!
!    end subroutine irfft


end module rpmd
