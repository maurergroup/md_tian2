module rpmd

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




    !    subroutine calc_centroid_virial_ekin(atoms, ekin_p, ekin_l)
    !
    !        type(universe), intent(in)  :: atoms
    !        real(dp),       intent(out) :: ekin_p, ekin_l
    !
    !        real(dp) :: centroids(3, atoms%natoms), vec(3)
    !        integer  :: i, b
    !
    !        call simple_ekin(atoms, ekin_p, ekin_l)
    !
    !        if (atoms%nbeads > 1) then
    !
    !            centroids = calc_centroid_positions(atoms)
    !
    !            do i = 1, atoms%natoms
    !                do b = 1, atoms%nbeads
    !
    !                    vec = atoms%r(:,b,i) - centroids(:,i)
    !
    !                    if (atoms%is_proj(atoms%idx(i))) then
    !                        ekin_p = ekin_p + sum(vec*atoms%f(:,b,i)) / (2.0_dp * atoms%nbeads)
    !                    else
    !                        ekin_l = ekin_l + sum(vec*atoms%f(:,b,i)) / (2.0_dp * atoms%nbeads)
    !                    end if
    !
    !                end do
    !            end do
    !
    !        end if
    !
    !    end subroutine calc_centroid_virial_ekin




    function calc_centroid_positions(atoms) result(cents)

        type(universe), intent(in) :: atoms
        real(dp)                   :: cents(3, atoms%natoms)

        integer  :: i

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




    function calc_inter_bead_distances(atoms) result (dists)

        type(universe), intent(in) :: atoms

        real(dp) :: dists(atoms%nbeads, atoms%natoms)
        real(dp) :: vec(3, atoms%nbeads, atoms%natoms)
        integer :: b, k

        do b = 1, atoms%nbeads
            k = modulo(b, atoms%nbeads)+1                   ! bead b+1
            vec(:,b,:) = atoms%r(:,b,:) - atoms%r(:,k,:)    ! distance vector
        end do

        dists = sqrt(sum(vec*vec, dim=1))    ! distance

    end function calc_inter_bead_distances



    subroutine calc_primitive_quantum_ekin(atoms, ekin_p, ekin_l)

        type(universe), intent(in)  :: atoms
        real(dp)      , intent(out) :: ekin_p, ekin_l

        real(dp) :: dx(atoms%nbeads, atoms%natoms), pref
        integer  :: i

        ekin_p = 0.0_dp
        ekin_l = 0.0_dp

        if (atoms%nbeads > 1) then

            pref = 0.5_dp * kB**2 * calc_instant_temperature(atoms)**2 * atoms%nbeads / hbar**2
            dx = calc_inter_bead_distances(atoms)

            do i = 1, atoms%natoms
               if (atoms%is_proj(atoms%idx(i))) then
                    ekin_p = ekin_p + sum(dx(:,i)*dx(:,i))*atoms%m(atoms%idx(i))
                else
                    ekin_l = ekin_l + sum(dx(:,i)*dx(:,i))*atoms%m(atoms%idx(i))
                end if
            end do

            ekin_p = ekin_p * pref
            ekin_l = ekin_l * pref

        end if

    end subroutine calc_primitive_quantum_ekin



    subroutine calc_virial_quantum_ekin(atoms, ekin_p, ekin_l)

        type(universe), intent(in)  :: atoms
        real(dp)      , intent(out) :: ekin_p, ekin_l

        real(dp) :: cents(3, atoms%natoms)
        integer  :: i, b

        ekin_p = 0.0_dp
        ekin_l = 0.0_dp

        if (atoms%nbeads > 1) then

            cents = calc_centroid_positions(atoms)

            do b = 1, atoms%nbeads
                if (atoms%is_proj(atoms%idx(i))) then
                    ekin_p = ekin_p + sum((atoms%r(:,b,:) - cents) * -atoms%f(:,b,:))
                else
                    ekin_l = ekin_l + sum((atoms%r(:,b,:) - cents) * -atoms%f(:,b,:))
                end if
            end do

            ekin_p = 0.5_dp * ekin_p / atoms%nbeads
            ekin_l = 0.5_dp * ekin_l / atoms%nbeads

        end if

    end subroutine calc_virial_quantum_ekin




    subroutine do_ring_polymer_step(atoms)

        type(universe), intent(inout) :: atoms

        real(dp), dimension(3, atoms%nbeads, atoms%natoms)  :: p, q, newP, newQ

        real(dp) :: poly(4, atoms%nbeads), p_new(3)
        real(dp) :: twown, wk, wt, wm, cos_wt, sin_wt
        real(dp) :: betaN, piN, mass

        integer :: i, b, k

        if (.not. allocated(cjk)) call build_cjk(atoms%nbeads)

        betaN = 1.0_dp / (kB * calc_instant_temperature(atoms) * atoms%nbeads)

        ! Transform to normal mode space
        call calc_momentum_all(atoms, p)
        q = atoms%r

        newP = 0.0_dp
        newQ = 0.0_dp
        do b = 1, atoms%nbeads
            do k = 1, atoms%nBeads
                newP(:,b,:) = newP(:,b,:) + p(:,k,:)*cjk(k,b)
                newQ(:,b,:) = newQ(:,b,:) + q(:,k,:)*cjk(k,b)
            end do
        end do
        p = newP
        q = newQ

        piN = pi / atoms%nbeads
        do i = 1, atoms%natoms
            mass = atoms%m(atoms%idx(i))
            poly(1, 1) = 1.0_dp
            poly(2, 1) = 0.0_dp
            poly(3, 1) = simparams%step / mass
            poly(4, 1) = 1.0_dp

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
                q(:,b,i) = p(:,b,i) * poly(3,b) + q(:,b,i) * poly(4,b)
                p(:,b,i) = p_new
            end do
        end do

        ! Transform back to Cartesian space
        where(.not. atoms%is_fixed)
            atoms%r = 0.0_dp
            atoms%v = 0.0_dp
        end where

        do i = 1, atoms%natoms
            mass = atoms%m(atoms%idx(i))
            do b = 1, atoms%nbeads
                do k = 1, atoms%nbeads
                    where (.not. atoms%is_fixed(:,b,i))
                        atoms%r(:,b,i) = atoms%r(:,b,i) + q(:,k,i)*cjk(b,k)
                        atoms%v(:,b,i) = atoms%v(:,b,i) + p(:,k,i)*cjk(b,k)/mass
                    end where
                end do
            end do
        end do

    end subroutine do_ring_polymer_step


end module rpmd
