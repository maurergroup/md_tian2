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


    subroutine ring_polymer_step(atoms)

        type(universe), intent(inout) :: atoms

        real(dp), dimension(3, atoms%nbeads, atoms%natoms)  :: p, q, newP, newQ

        real(dp) :: poly(4, atoms%nbeads)
        real(dp), dimension(3) :: p_new
        real(dp) :: twown, wk, wt, wm, cos_wt, sin_wt
        real(dp)                                      :: betaN, piN

        integer :: i, b, k

        if (.not. allocated(cjk)) call build_cjk(atoms%nbeads)

        betaN = 1.0_dp / (kB * simparams%Tsurf * atoms%nbeads)


        !        ! Transform to normal mode space
        !        do i = 1, 3
        !            do j = 1, nParticles
        !                call rfft(p(i,j,:), nBeads)
        !                call rfft(q(i,j,:), nBeads)
        !            end do
        !        end do

        ! Transform to normal mode space
        call calc_momentum_all(atoms, p)
        q = atoms%r

        newP = 0.0_dp
        newQ = 0.0_dp
        do b = 1, atoms%nbeads
            do k = 1, atoms%nBeads
                newP(:,b,:) = newP(:,b,:) + p(:,b,:)*cjk(k,b)
                newQ(:,b,:) = newQ(:,b,:) + q(:,b,:)*cjk(k,b)
            end do
        end do
        p = newP
        q = newQ


        piN = pi / atoms%nbeads
        do i = 1, atoms%natoms
            poly(1, 1) = 1.0_dp
            poly(2, 1) = 0.0_dp
            poly(3, 1) = simparams%step / atoms%m(atoms%idx(i))
            poly(4, 1) = 1.0_dp

            if (atoms%nbeads > 1) then
                twown = 2.0_dp / betaN
                do b = 1, atoms%nbeads / 2
                    wk = twown * sin(b * piN)
                    wt = wk * simparams%step
                    wm = wk * atoms%m(atoms%idx(i))
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
            do b = 1, atoms%nbeads
                do k = 1, atoms%nbeads
                    !                newP(:,:,j) = newP(:,:,j) + p(:,:,k)*cjk(j,k)
                    !                newQ(:,:,j) = newQ(:,:,j) + q(:,:,k)*cjk(j,k)
                    where (.not. atoms%is_fixed(:,b,i))
                        atoms%r(:,b,i) = atoms%r(:,b,i) + q(:,k,i)*cjk(b,k)
                        atoms%v(:,b,i) = atoms%v(:,b,i) + p(:,k,i)*cjk(b,k)/atoms%m(atoms%idx(i))
                    end where
                end do
            end do
        end do
    !        p = newP
    !        q = newQ

        ! Transform back to Cartesian space
    !        do i = 1, 3
    !            do j = 1, nParticles
    !                call irfft(p(i,j,:), nBeads)
    !                call irfft(q(i,j,:), nBeads)
    !            end do
    !        end do



    end subroutine ring_polymer_step


end module rpmd
