module pes_emt_mod

    use constants
    use useful_things, only : split_string, lower_case
    use universe_mod

    implicit none

    type emt_pes
        private
        real(dp), allocatable :: eta2(:)
        real(dp), allocatable :: n0(:)
        real(dp), allocatable :: e0(:)
        real(dp), allocatable :: lambda(:)
        real(dp), allocatable :: v0(:)
        real(dp), allocatable :: kappa(:)
        real(dp), allocatable :: s0(:)
    end type emt_pes

    type(emt_pes) :: pes_emt

    real(dp), allocatable :: dens(:,:)

contains

    subroutine read_emt(atoms, inp_unit)

        use run_config, only : simparams

        type(universe), intent(inout) :: atoms
        integer, intent(in) :: inp_unit

        integer :: nwords, ios = 0
        character(len=max_string_length) :: buffer
        character(len=max_string_length) :: words(100)
        integer  :: idx1, idx2, ntypes, param_counter
        character(len=*), parameter :: err = "Error in read_lj: "

        ntypes = simparams%nprojectiles+simparams%nlattices

        if (.not. allocated(pes_emt%e0)) then
            allocate(pes_emt%eta2  (ntypes), &
                pes_emt%n0    (ntypes), &
                pes_emt%e0    (ntypes), &
                pes_emt%lambda(ntypes), &
                pes_emt%v0    (ntypes), &
                pes_emt%kappa (ntypes), &
                pes_emt%s0    (ntypes), &
                dens(atoms%nbeads, atoms%natoms)&
                )

            pes_emt%eta2   = default_real
            pes_emt%n0     = default_real
            pes_emt%e0     = default_real
            pes_emt%lambda = default_real
            pes_emt%v0     = default_real
            pes_emt%kappa  = default_real
            pes_emt%s0     = default_real
            dens           = default_real

        end if

        ! line should read something like "H   H   proj    proj"
        read(inp_unit, '(A)', iostat=ios) buffer
        call split_string(buffer, words, nwords)

        if (nwords /= 4) stop err // "need four entries in interaction-defining lines"

        if (words(3) == "proj" .and. words(4) == "proj" .or. &
            words(3) == "proj" .and. words(4) == "latt" .or. &
            words(3) == "latt" .and. words(4) == "proj" .or. &
            words(3) == "latt" .and. words(4) == "latt") then

            idx1 = get_idx_from_name(atoms, words(1), is_proj=(words(3)=="proj"))
            idx2 = get_idx_from_name(atoms, words(2), is_proj=(words(4)=="proj"))

            if (atoms%pes(idx1,idx2) /= default_int) then
                print *, err // "pes already defined for atoms", words(1), words(3), words(2), words(4)
                stop
            end if

        else
            stop err // "interaction must be defined via 'proj' and 'latt' keywords"
        end if

        ! set the pes type in the atoms object
        atoms%pes(idx1,idx2) = pes_id_emt
        atoms%pes(idx2,idx1) = pes_id_emt

        param_counter = 1
        do
            read(inp_unit, '(A)', iostat=ios) buffer
            call split_string(buffer, words, nwords)

            ! pes block terminated, exit
            if (nwords == 0 .or. ios /= 0) then
                exit

            ! something went wrong
            else if (nwords /= 2) then
                stop "Error in the PES file: PES parameters must consist of key value pairs. A parameter block must be terminated by a blank line."
            end if

            call lower_case(words(1))

            select case (words(1))

                case ('eta2')   ! A^-1
                    if (param_counter <= nparams_emt) then
                        read(words(2), *) pes_emt%eta2(idx1)
                    else
                        read(words(2), *) pes_emt%eta2(idx2)
                    end if

                case ('n0')     ! A^-3
                    if (param_counter <= nparams_emt) then
                        read(words(2), *) pes_emt%n0(idx1)
                    else
                        read(words(2), *) pes_emt%n0(idx2)
                    end if

                case ('e0')     ! eV
                    if (param_counter <= nparams_emt) then
                        read(words(2), *) pes_emt%e0(idx1)
                    else
                        read(words(2), *) pes_emt%e0(idx2)
                    end if

                case ('lambda') ! A^-1
                    if (param_counter <= nparams_emt) then
                        read(words(2), *) pes_emt%lambda(idx1)
                    else
                        read(words(2), *) pes_emt%lambda(idx2)
                    end if

                case ('v0')     ! eV
                    if (param_counter <= nparams_emt) then
                        read(words(2), *) pes_emt%v0(idx1)
                    else
                        read(words(2), *) pes_emt%v0(idx2)
                    end if

                case ('kappa')  ! A^-1
                    if (param_counter <= nparams_emt) then
                        read(words(2), *) pes_emt%kappa(idx1)
                    else
                        read(words(2), *) pes_emt%kappa(idx2)
                    end if

                case ('s0')     ! A
                    if (param_counter <= nparams_emt) then
                        read(words(2), *) pes_emt%s0(idx1)
                    else
                        read(words(2), *) pes_emt%s0(idx2)
                    end if

                case default
                    print *, "Error in the PES file: unknown LJ parameter", words(1)
                    stop

            end select

            param_counter = param_counter + 1
        end do

    end subroutine read_emt



    subroutine compute_emt(atoms, flag)

        !   Calculates energy and forces with EMT potential

        type(universe), intent(inout) :: atoms
        integer, intent(in)           :: flag

        integer :: i, j, k, b, n, idx_i, idx_j

        real(dp), parameter :: cutoff = 1.5_dp

        real(dp) :: beta, rcut, rr,  acut, r, theta, rtemp

        ! Energies
        real(dp), dimension(3) :: nneighs, vec
        real(dp), dimension(atoms%nbeads) ::  Ecoh
        real(dp), dimension(atoms%ntypes) :: betaeta2, kappadbeta,  betas0, igamma1, igamma2

        real(dp), dimension(atoms%ntypes, 3)            :: rnn, r3temp, x
        real(dp), dimension(atoms%ntypes, atoms%ntypes) :: chi
        real(dp), dimension(atoms%ntypes, atoms%nbeads) :: V_ref
        real(dp), dimension(atoms%nbeads, atoms%natoms) :: s

        real(dp), dimension(atoms%ntypes, atoms%ntypes, atoms%nbeads) :: V
        real(dp), dimension(atoms%nbeads, atoms%natoms, atoms%natoms) :: sigma

        ! Forces
        real(dp) :: rtemp1

        real(dp), dimension(3) :: dtheta

        real(dp), dimension(3, atoms%nbeads, atoms%ntypes, atoms%natoms, atoms%natoms) :: ds
        real(dp), dimension(3, atoms%nbeads, atoms%ntypes, atoms%ntypes, &
            atoms%natoms, atoms%natoms) :: dsigma, dV
        real(dp), dimension(3, atoms%nbeads, atoms%ntypes, atoms%ntypes, atoms%natoms) :: dsigmat
        real(dp), dimension(3, atoms%nbeads, atoms%ntypes, atoms%natoms, atoms%natoms) :: dEcoh
        real(dp), dimension(3, atoms%nbeads, atoms%ntypes, atoms%natoms) :: dVref
        real(dp), dimension(3, atoms%nbeads, atoms%natoms) :: dVf



        !----------------------VALUES OF FREQUENT USE ---------------------------------

        beta = select_beta('fcc')
        nneighs = select_nneighs('fcc')

        ! beta * s0
        betas0 = beta * pes_emt%s0

!        print *, "betas0"
!        print '(2f12.8)', betas0

        ! beta * eta2
        betaeta2 = beta * pes_emt%eta2

!        print *, "betaeta2"
!        print '(2f12.8)', betaeta2
        ! kappa / beta
        kappadbeta = pes_emt%kappa / beta

        ! 'coupling' parameters between species
        do i = 1, atoms%ntypes
            do j = 1, atoms%ntypes
                chi(i,j) = pes_emt%n0(i) / pes_emt%n0(j) * &
                    exp(0.5_dp/bohr2ang*(pes_emt%s0(j)-pes_emt%s0(i)))
            end do
        end do

        print *, "chi"
        do i = 1, atoms%ntypes
        do j = 1, atoms%ntypes
        print  '(2i4, 2f19.15)', i,j,chi(i,j)
        end do
        end do

        ! Distances to the nearest, next-nearest and next-next-nearest neighbours
        rnn(:,1) = betas0
        rnn(:,2) = rnn(:,1) * sqrt2
        rnn(:,3) = rnn(:,1) * sqrt3

        !------------------------------------------------------------------------------
        !                                  CUT-OFF
        !                                  =======
        !------------------------------------------------------------------------------
        ! We use the distance to the next-next-nearest neighbours as cut-off.
        ! We only need one cut-off and we choose one of the lattice atoms since s0
        ! is usually larger for them.

        rcut = maxval(betas0) * sqrt3

!        print *, "rcut"
!        print  '(2f19.15)', rcut

        !rcut = a_lat * sqrt3 * isqrt2
        rr = select_rr('fcc', rcut)
        acut = 9.21034037197618_dp/(rr - rcut) ! ln(10000)

!        print *, "acut"
!        print '(2f19.15)', acut

        do k = 1, 3
            x(:,k) = nneighs(k) * twelfth / (1.0_dp + exp(acut*(rnn(:,k)-rcut)))
            r3temp(:,k) = rnn(:,k) - betas0
        end do

!        print *, "x"
!        do i = 1, 3
!            print '(2f19.15)', x(:,i)
!        end do
        !-----------------------------------GAMMA--------------------------------------
        ! Gamma enforces the cut-off together with theta (see below)
        ! Gamma is defined as inverse.


        igamma1 = 0.0_dp
        igamma2 = 0.0_dp
        do k = 1, 3
            igamma1 = igamma1 + x(:,k)*exp(-pes_emt%eta2 * r3temp(:,k))
            igamma2 = igamma2 + x(:,k)*exp(-kappadbeta   * r3temp(:,k))
        end do
        igamma1 = 1.0_dp / igamma1
        igamma2 = 1.0_dp / igamma2

!        print *, "igamma1"
!        print '(2f19.15)', igamma1
!        print *, "igamma2"
!        print '(2f19.15)', igamma2




        !------------------------------------------------------------------------------
        !                          Sigma and Pair-wise Contributions
        !                          =================================
        !------------------------------------------------------------------------------

        ! initialize accumulators
        sigma = 0.0_dp
        V     = 0.0_dp
        V_ref = 0.0_dp
        s     = 0.0_dp
        Ecoh  = 0.0_dp
        dens  = 0.0_dp

        dsigma = 0.0_dp
        dsigmat  = 0.0_dp
        dV = 0.0_dp
        dVf = 0.0_dp
        ds = 0.0_dp
        dEcoh = 0.0_dp
        dVref = 0.0_dp

        do i = 1, atoms%natoms-1
            idx_i = atoms%idx(i)

            do j = i+1, atoms%natoms
                idx_j = atoms%idx(j)

                if (i == j .or. atoms%pes(idx_i,idx_j) /= pes_id_emt) cycle

                do b = 1, atoms%nbeads


                    ! Applying PBCs
                    call minimg_one(atoms, i, j, b, r, vec)

                    ! drops atoms outside (cutoff*rcut)-sphere
                    !if (r > cutoff*rcut) cycle
                    if (r > cutoff*rcut) cycle

                    vec = vec/r                       ! unit vector j -> i

                    ! cut-off function
                    rtemp = exp(acut*(r - rcut))
                    theta = 1.0_dp / (1.0_dp + rtemp)
                    rtemp1 = acut * rtemp * theta

                    ! sigma_1  i -> j
                    rtemp = theta*exp(-pes_emt%eta2(idx_i) * (r - betas0(idx_i)))
                    sigma(b,i,j) = sigma(b,i,j) + rtemp
                    dtheta = (pes_emt%eta2(idx_i) + rtemp1)*rtemp*vec
                    dsigmat(:,b,idx_j,idx_i,j) = dsigmat(:,b,idx_j,idx_i,j) - dtheta
                    dsigma(:,b,idx_j,idx_i,j,i) = dtheta


                    ! sigma_1  j <- i
                    rtemp = theta*exp(-pes_emt%eta2(idx_j) * (r - betas0(idx_j)))
                    sigma(b,j,i) = sigma(b,j,i) + rtemp
                    dtheta = (pes_emt%eta2(idx_j) + rtemp1)*rtemp*vec
                    dsigmat(:,b,idx_i,idx_j,i) = dsigmat(:,b,idx_i,idx_j,i) + dtheta
                    dsigma(:,b,idx_i,idx_j,i,j) =  -dtheta


                    ! sigma_2
                    rtemp = theta*exp(-kappadbeta(idx_i)   * (r - betas0(idx_i)))
                    V(idx_i,idx_j,b) = V(idx_i,idx_j,b) + rtemp
                    dtheta = (kappadbeta(idx_i) + rtemp1) * rtemp * vec
                    dV(:,b,idx_j,idx_i,i,j) = dV(:,b,idx_j,idx_i,i,j) + dtheta
                    dV(:,b,idx_j,idx_i,j,i) = dV(:,b,idx_j,idx_i,j,i) - dtheta

                    rtemp = theta*exp(-kappadbeta(idx_j)   * (r - betas0(idx_j)))
                    V(idx_j,idx_i,b) = V(idx_j,idx_i,b) + rtemp
                    if (idx_i /= idx_j) then
                        dtheta = (kappadbeta(idx_j) + rtemp1) * rtemp * vec
                        dV(:,b,idx_i,idx_j,j,i) = dV(:,b,idx_i,idx_j,j,i) + dtheta
                        dV(:,b,idx_i,idx_j,i,j) = dV(:,b,idx_i,idx_j,i,j) - dtheta
                    end if

                end do
            end do
        end do


!        print *, "dV"
!        do i = 1, atoms%natoms
!            do j = 1, atoms%natoms
!                do k = 1, atoms%ntypes
!                    do n = 1, atoms%ntypes
!                        print '(4i4, 3f12.5)', k, n, i, j, dV(:,:,k,n,i,j)
!                    end do
!                end do
!            end do
!        end do
        !        print *, "sigma"
        !        print '(3f15.8)', sigma

        !        print *, "dsigma"
        !        do i = 1, atoms%natoms
        !            do j = 1, atoms%natoms
        !                do k = 1, atoms%ntypes
        !                    do n = 1, atoms%ntypes
        !                        print '(4i4, 3f12.5)',k, n, i, j, dsigma(:,1,k,n,i,j)
        !                    end do
        !                end do
        !            end do
        !        end do
        !
        !        print *, "dsigma"
        !        do i = 1, atoms%natoms
        !            do j = 1, atoms%natoms
        !                print '(2i4, 3f12.5)', i, j, sum(sum(dsigma(:,:,:,:,i,j), dim=3), dim=3)
        !            end do
        !        end do
        !
        !        print *, "dsigmat"
        !        do i = 1, atoms%natoms
        !            do k = 1, atoms%ntypes
        !                do n = 1, atoms%ntypes
        !                    print '(3i4, 3f12.5)', k,n, i, dsigmat(:,:,k,n,i)
        !                end do
        !            end do
        !        end do

        !        print *, "V", sum(V)
        !        print '(2f15.8)', V(:,:,1)


        ! divide by cut-off scaling factors

        do b = 1, atoms%nbeads
            do j = 1, atoms%ntypes
                V(:,j,b) = V(:,j,b) * igamma2(j) * pes_emt%v0(j) * chi(:,j)
            end do
        end do

        do j = 1, atoms%natoms
            idx_j = atoms%idx(j)
            sigma(:,:,j) = sigma(:,:,j) * igamma1(idx_j)
        end do

        do n = 1, atoms%ntypes
            dsigma (:,:,n,:,:,:) = dsigma (:,:,n,:,:,:) * igamma1(n)
            dsigmat(:,:,n,:,:)   = dsigmat(:,:,n,:,:)   * igamma1(n)
        end do

        do k = 1, atoms%ntypes
            do n = 1, atoms%ntypes
                dV(:,:,n,k,:,:) = dV(:,:,n,k,:,:) * igamma2(n) *  pes_emt%v0(n) * chi(k,n)
            end do
        end do


        ! merge into force-suited structure (basically reverse the 1st loop)
        do i = 1, atoms%natoms-1
        idx_i = atoms%idx(i)
            do j = i+1, atoms%natoms
            idx_j = atoms%idx(j)

            if (idx_i /= idx_j) then
            dVf(:,:,i) = dVf(:,:,i) - 0.5*dV(:,:,idx_j,idx_i,i,j)
            dVf(:,:,j) = dVf(:,:,j) - 0.5*dV(:,:,idx_j,idx_i,j,i)

            dVf(:,:,i) = dVf(:,:,i) - 0.5*dV(:,:,idx_i,idx_j,j,i)
            dVf(:,:,j) = dVf(:,:,j) - 0.5*dV(:,:,idx_i,idx_j,i,j)
            else
            dVf(:,:,i) = dVf(:,:,i) - dV(:,:,idx_j,idx_i,i,j)
            dVf(:,:,j) = dVf(:,:,j) - dV(:,:,idx_j,idx_i,j,i)
            end if
            end do
        end do

        print '(3f12.5)' ,dVf
!        print *, "dV scaled"
!        do i = 1, atoms%natoms
!            do j = 1, atoms%natoms
!                do k = 1, atoms%ntypes
!                    do n = 1, atoms%ntypes
!                        print '(4i4, 3f12.5)',n, k, i, j, dV(:,:,n,k,i,j)
!                    end do
!                end do
!            end do
!        end do


        !        do j = 1, atoms%natoms
        !            idx_j = atoms%idx(j)
        !            dsigma(:,:,:,j) = dsigma(:,:,:,j) * igamma1(idx_j)
        !        end do

        !        print *, "sigma scaled"
        !        print '(3f15.8)', sigma



!        print *, "dsigma scaled"
!        do i = 1, atoms%natoms
!            do j = 1, atoms%natoms
!                do k = 1, atoms%ntypes
!                    do n = 1, atoms%ntypes
!                        print '(4i4, 3f12.5)',k, n, i, j, dsigma(:,1,k,n,i,j)
!                    end do
!                end do
!            end do
!        end do

        print *, "dsigmat scaled"
        do i = 1, atoms%natoms
            do k = 1, atoms%ntypes
                do n = 1, atoms%ntypes
                    print '(3i4, 3f12.5)', k,n, i, dsigmat(:,:,k,n,i)
                end do
            end do
        end do



        !        print *, "V scaled"
        !        print '(2f15.8)', V(:,:,1)



        !-----------------------------NEUTRAL SPHERE RADIUS----------------------------

        !        print *, "chi"
        !        print *, chi
        do i = 1, atoms%natoms
            idx_i = atoms%idx(i)
            do j = 1, atoms%natoms
                idx_j = atoms%idx(j)
                s(:,i) = s(:,i) + chi(idx_j,idx_i)*sigma(:,j,i)
            end do
        end do

        do j = 1, atoms%natoms
            idx_j = atoms%idx(j)
            do i = 1, atoms%natoms
                idx_i = atoms%idx(i)
                if (idx_i == idx_j) then
                    ds(:,:,idx_i,i,j) = -dsigma (:,:,idx_i,idx_j,i,j)
                    ds(:,:,idx_j,j,j) = -dsigmat(:,:,idx_j,idx_j,j)
                end if
            end do
        end do

        !        print *, "s"
        !        print '(3f15.8)',   s
        !
        !        print *, "betaeta2"
        !        print '(2f15.8)', betaeta2

!        print *, "ds"
!        do j = 1, atoms%natoms
!            do i = 1, atoms%natoms
!                do n = 1, atoms%ntypes
!                    print '(3i, 3f12.5)', n, i, j,  ds(:,:,n,i,j)
!                end do
!            end do
!        end do

        do i = 1, atoms%natoms

            idx_i = atoms%idx(i)
            do n = 1, atoms%ntypes
                if (idx_i /= n) then
                    ds(:,:,idx_i,i,i) = ds(:,:,idx_i,i,i) - chi(n,idx_i)*dsigmat(:,:,idx_i,n,i)
                end if

            end do

            do k = 1, 3
                ds(k,:,idx_i,:,i) = ds(k,:,idx_i,:,i)/(betaeta2(idx_i)*s)
            end do

            do j = 1, atoms%natoms
                idx_j = atoms%idx(j)
                if (idx_i /= idx_j) then
                    do k = 1, 3
                        ds(k,:,idx_i,i,j) = -chi(idx_j,idx_i)*dsigma(k,:,idx_i,idx_j,i,j)/(betaeta2(idx_i)*s(:,i))
                    end do
                end if
            end do

        end do

        print *, "ds"
        do j = 1, atoms%natoms
            do i = 1, atoms%natoms
                do n = 1, atoms%ntypes
                    print '(3i, 3f12.5)', n, i, j,  ds(:,:,n,i,j)
                end do
            end do
        end do


        do i = 1, atoms%natoms
            idx_i = atoms%idx(i)
            do b = 1, atoms%nbeads
                if (s(b,i) > 0) s(b,i) = -log(s(b,i)*twelfth)/betaeta2(idx_i)
            end do
        end do

!        print *, "s2"
!        print '(3f15.8)',  s


        !----------------------EMBEDDED ELECTRON DENSITY-------------------------------

        do i = 1, atoms%natoms
            idx_i = atoms%idx(i)
            dens(:,i) = pes_emt%n0(idx_i) * exp((0.5_dp/bohr2ang - betaeta2(idx_i))*s(:,i))
        end do

!        print *, "dens"
!        print '(3f15.8)', dens




        !---------------------------COHESIVE FUNCTION-----------------------------------

        do i = 1, atoms%natoms
            idx_i = atoms%idx(i)

            Ecoh = Ecoh + ((1.0_dp + pes_emt%lambda(idx_i)*s(:,i)) * &
                exp(-pes_emt%lambda(idx_i)*s(:,i)) - 1.0_dp) * pes_emt%e0(idx_i)
        end do

!        print *, "Ecoh"
!        print '(2f15.8)', Ecoh


        do i = 1, atoms%natoms
            do j = 1, atoms%natoms
                idx_j = atoms%idx(j)
                do k = 1, 3
                    dEcoh(k,:,idx_j,j,i) = s(:,j) * &
                        exp(-pes_emt%lambda(idx_j)*s(:,j))*ds(k,:,idx_j,j,i)

                end do

!                dEcoh(k,:,i) = sum(s * exp(-pes_emt%lambda(idx_i) * s) * ds(k,:,idx_i,:,i), dim=2)
            end do
        end do

        do n = 1, atoms%ntypes
            dEcoh(:,:,n,:,:) = dEcoh(:,:,n,:,:) * pes_emt%e0(n) * pes_emt%lambda(n) * pes_emt%lambda(n)
        end do

!        print *, "dEcoh"
!        do i = 1, atoms%natoms
!            do  j = 1, atoms%natoms
!                do n = 1, atoms%ntypes
!                    print '(3i4, 3f12.5)', n,i, j, dEcoh(:,:,n,i,j)
!                end do
!            end do
!        end do


         !----------------REFERENCE PAIR POTENTIAL CONTRIBUTIONS------------------------

        do i=1,atoms%natoms
            idx_i = atoms%idx(i)
            V_ref(idx_i,:) = V_ref(idx_i,:) + exp(-pes_emt%kappa(idx_i)*s(:,i))
        end do


        do i = 1, atoms%natoms
        do j = 1, atoms%natoms
        do k = 1, 3
            idx_j = atoms%idx(j)
            dVref(k,:,idx_j,i) = dVref(k,:,idx_j,i) + exp(-pes_emt%kappa(idx_j)*s(:,j))*ds(k,:,idx_j,j,i)
        end do
        end do
        end do

!        print *, "dVref"
!        do i = 1, atoms%natoms
!        do n = 1, atoms%ntypes
!            print '(2i4, 3f12.5)', n, i, dVref(:,:,n,i)
!        end do
!        end do
!
!stop
!        print *, "V_ref"
!        print '(2f15.8)',  V_ref


        do b = 1, atoms%nbeads
            V_ref(:,b) = V_ref(:,b) * 12.0_dp * pes_emt%v0
        end do

        do n = 1, atoms%ntypes
            dVref(:,:,n,:) = dVref(:,:,n,:) * 12.0_dp * pes_emt%v0(n) * pes_emt%kappa(n)
        end do

!        print *, "V_ref_2"
!        print '(2f15.8)',  V_ref
!        print *, "dVref"
!        do i = 1, atoms%natoms
!        do n = 1, atoms%ntypes
!            print '(2i4, 3f12.5)', n, i, dVref(:,:,n,i)
!        end do
!        end do

        print *, "dV"
        do i = 1, atoms%natoms
        do j = 1, atoms%natoms
        do n = 1, atoms%ntypes
        do k = 1, atoms%ntypes
            print '(4i4, 3f12.5)',  n,k,i, j,dV(:,:,n,k,i,j)

        end do
        end do
        end do
        end do



        !-------------------------------TOTAL ENERGY---------------------------------

        atoms%epot = Ecoh - 0.5_dp * (sum(sum(V, dim=1), dim=1) - sum(V_ref, dim=1))

        atoms%f = sum(sum(dEcoh, dim=3), dim=3) - dVf
        print *, "atoms%f"
        print '(3f12.5)', atoms%f

        print *, "Epot", atoms%epot


     !   slab%f = dEcoh_l_l + dEcoh_p_l - dV_ll_l - 0.50d0*(dV_lp_l + dV_pl_l - dvref_l_l - dvref_p_l)

        stop
    end subroutine compute_emt



    function select_nneighs(geom) result(neighs)

        character(len=3), intent(in) :: geom
        real(dp) :: neighs(3)

        select case (geom)
            case ('fcc')
                neighs = [12.0_dp, 6.0_dp, 24.0_dp]
            case default
                print *, "Error in select_nneighs(): unknown lattice structure", geom ; call abort
        end select

    end function select_nneighs



    real(dp) function select_beta(geom) result(beta)

        character(len=3), intent(in) :: geom

        select case (geom)
            case ('fcc')
                beta = (16.0_dp*pi/3)**(1.0_dp/3)*isqrt2
            case ('bcc')
                beta = (3.0_dp*pi*pi)**(1.0_dp/6.0_dp)
            case default
                print *, "Error in select_beta(): unknown lattice structure", geom ; call abort
        end select

    end function select_beta


    function select_rr(geom, rcut) result(rr)

        character(len=3), intent(in) :: geom
        real(dp), intent(in) :: rcut
        real(dp)             :: rr

        select case (geom)
            case ('fcc')
                rr = 4.0_dp * rcut / (sqrt3 + 2.0_dp)
            case ('bcc')
                rr = 4.0_dp * rcut / (sqrt(8.0_dp/3.0_dp) + sqrt(11.0_dp/3.0_dp))
            case default
                print *, "Error in select_rr(): unknown lattice structure", geom ; call abort
        end select

    end function select_rr


end module pes_emt_mod
