module pes_emt_mod

    use constants
    use useful_things, only : split_string, lower_case
    use universe_mod

    implicit none

    type emt_pes
        private
        real(dp), allocatable :: eta2(:,:)
        real(dp), allocatable :: n0(:,:)
        real(dp), allocatable :: e0(:,:)
        real(dp), allocatable :: lambda(:,:)
        real(dp), allocatable :: v0(:,:)
        real(dp), allocatable :: kappa(:,:)
        real(dp), allocatable :: s0(:,:)
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
            allocate(pes_emt%eta2  (ntypes, ntypes), &
                pes_emt%n0    (ntypes, ntypes), &
                pes_emt%e0    (ntypes, ntypes), &
                pes_emt%lambda(ntypes, ntypes), &
                pes_emt%v0    (ntypes, ntypes), &
                pes_emt%kappa (ntypes, ntypes), &
                pes_emt%s0    (ntypes, ntypes), &
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
                        read(words(2), *) pes_emt%eta2(idx1, idx2)
                    else
                        read(words(2), *) pes_emt%eta2(idx2, idx1)
                    end if

                case ('n0')     ! A^-3
                    if (param_counter <= nparams_emt) then
                        read(words(2), *) pes_emt%n0(idx1, idx2)
                    else
                        read(words(2), *) pes_emt%n0(idx2, idx1)
                    end if

                case ('e0')     ! eV
                    if (param_counter <= nparams_emt) then
                        read(words(2), *) pes_emt%e0(idx1, idx2)
                    else
                        read(words(2), *) pes_emt%e0(idx2, idx1)
                    end if

                case ('lambda') ! A^-1
                    if (param_counter <= nparams_emt) then
                        read(words(2), *) pes_emt%lambda(idx1, idx2)
                    else
                        read(words(2), *) pes_emt%lambda(idx2, idx1)
                    end if

                case ('v0')     ! eV
                    if (param_counter <= nparams_emt) then
                        read(words(2), *) pes_emt%v0(idx1, idx2)
                    else
                        read(words(2), *) pes_emt%v0(idx2, idx1)
                    end if

                case ('kappa')  ! A^-1
                    if (param_counter <= nparams_emt) then
                        read(words(2), *) pes_emt%kappa(idx1, idx2)
                    else
                        read(words(2), *) pes_emt%kappa(idx2, idx1)
                    end if

                case ('s0')     ! A
                    if (param_counter <= nparams_emt) then
                        read(words(2), *) pes_emt%s0(idx1, idx2)
                    else
                        read(words(2), *) pes_emt%s0(idx2, idx1)
                    end if

                case default
                    print *, "Error in the PES file: unknown LJ parameter", words(1)
                    stop

            end select

            param_counter = param_counter + 1
        end do

    end subroutine read_emt


    !    subroutine compute_emt(atoms, flag)
    !
    !        !   Calculates energy and forces with EMT potential
    !
    !        type(universe), intent(inout) :: atoms
    !        integer, intent(in)           :: flag
    !
    !        integer :: i, j, k, b, idx_i, idx_j
    !
    !        real(dp), parameter :: cutoff = 1.5_dp
    !
    !        real(dp) :: beta
    !        real(dp), dimension(atoms%ntypes,atoms%ntypes) :: betas0, betaeta2, kappadbeta, chi, &
    !            rcut, rr,  acut, igamma1, igamma2
    !
    !        real(dp), dimension(3,atoms%ntypes,atoms%ntypes) :: rnn, x, r3temp
    !        real(dp), dimension(3) :: vec
    !        real(dp) :: r, theta, rtemp, rtemp1
    !
    !
    !        real(dp), dimension(3) :: dtheta
    !
    !        real(dp), dimension(atoms%nbeads) :: V_ij, V_ji
    !
    !        real(dp), dimension(atoms%nbeads,atoms%natoms) :: sigma_ij, sigma_ji
    !
    !        real(dp), dimension(3,atoms%nbeads,atoms%natoms) :: dsigma_ij_i, dsigma_ji_j
    !        real(dp), dimension(3,atoms%nbeads,atoms%natoms,atoms%natoms) :: dsigma_ij_j, dsigma_ji_i
    !
    !        real(dp), dimension(3,atoms%nbeads,atoms%natoms) :: dV_ij_i,dV_ij_j, dV_ji_i, dV_ji_j
    !
    !
    !        real(dp) :: Ecoh_i, Ecoh_j, vref_i, vref_j
    !
    !        real(dp), dimension(3) :: nneighs
    !
    !        real(dp) :: s
    !
    !        real(dp) :: ds_i_i(3), ds_j_j(3), ds_i_j(3), ds_j_i(3)
    !
    !        real(dp) :: dEcoh_i_i(3), dEcoh_i_j(3), dEcoh_j_i(3), dEcoh_j_j(3)
    !
    !
    !        real(dp), dimension(3) :: dvref_i_i, dvref_i_j, dvref_j_i, dvref_j_j
    !
    !
    !
    !        !----------------------VALUES OF FREQUENT USE ---------------------------------
    !
    !        beta = select_beta('fcc')
    !        nneighs = select_nneighs('fcc')
    !
    !        ! beta * s0
    !        betas0 = beta * pes_emt%s0
    !
    !        ! beta * eta2
    !        betaeta2 = beta * pes_emt%eta2
    !        ! kappa / beta
    !        kappadbeta = pes_emt%kappa / beta
    !
    !        ! 'coupling' parameters between p and l
    !        chi = pes_emt%n0 / transpose(pes_emt%n0) &
    !            * exp(0.5_dp/bohr2ang*(pes_emt%s0-transpose(pes_emt%s0)))
    !
    !        ! Distances to the nearest, next-nearest and next-next-nearest neighbours
    !        rnn(1,:,:) = betas0
    !        rnn(2,:,:) = rnn(1,:,:) * sqrt2
    !        rnn(3,:,:) = rnn(1,:,:) * sqrt3
    !
    !        !------------------------------------------------------------------------------
    !        !                                  CUT-OFF
    !        !                                  =======
    !        !------------------------------------------------------------------------------
    !        ! We use the distance to the next-next-nearest neighbours as cut-off.
    !        ! We only need one cut-off and we choose one of the lattice atoms since s0
    !        ! is usually larger for them.
    !
    !        rcut = betas0 * sqrt3
    !        !rcut = a_lat * sqrt3 * isqrt2
    !        rr = 4 * rcut / (sqrt3 + 2.0d0)
    !        acut = 9.21034037197618_dp/(rr -rcut) ! ln(10000)
    !
    !        do k = 1, 3
    !            x(k,:,:) = nneighs(k) * twelfth / (1.0_dp + exp(acut*(rnn(k,:,:)-rcut)))
    !            r3temp(k,:,:) = rnn(k,:,:) - betas0
    !        end do
    !
    !        !-----------------------------------GAMMA--------------------------------------
    !        ! Gamma enforces the cut-off together with theta (see below)
    !        ! Gamma is defined as inverse.
    !
    !
    !        igamma1 = 0.0_dp
    !        igamma2 = 0.0_dp
    !        do k = 1, 3
    !            igamma1 = igamma1 + x(k,:,:)*exp(-pes_emt%eta2 * r3temp(k,:,:))
    !            igamma2 = igamma2 + x(k,:,:)*exp(-kappadbeta   * r3temp(k,:,:))
    !        end do
    !        igamma1 = 1.0_dp / igamma1
    !        igamma2 = 1.0_dp / igamma2
    !
    !
    !
    !
    !        !------------------------------------------------------------------------------
    !        !                          Sigma and Pair-wise Contributions
    !        !                          =================================
    !        !------------------------------------------------------------------------------
    !
    !        ! initialize accumulators
    !        sigma_ij       = 0.0_dp
    !        sigma_ji       = 0.0_dp
    !        V_ij           = 0.0_dp
    !        V_ji     = 0.0_dp
    !        dV_ij_i          = 0.0_dp
    !        dV_ij_j   = 0.0_dp
    !        dV_ji_i   = 0.0_dp
    !        dV_ji_j   = 0.0_dp
    !        dsigma_ij_i      = 0.0_dp
    !        dsigma_ij_j      = 0.0_dp
    !        dsigma_ji_i      = 0.0_dp
    !        dsigma_ji_j      = 0.0_dp
    !        vref_i      = 0.0_dp
    !        vref_j      = 0.0_dp
    !        dvref_i_i   = 0.0_dp
    !        dvref_i_j   = 0.0_dp
    !        dvref_j_i   = 0.0_dp
    !        dvref_j_j   = 0.0_dp
    !
    !
    !        do i = 1, atoms%natoms
    !            do j = 1, atoms%natoms
    !
    !                idx_i = atoms%idx(i)
    !                idx_j = atoms%idx(j)
    !
    !                if (i == j .or. atoms%pes(idx_i,idx_j) /= pes_id_emt) cycle
    !
    !                do b = 1, atoms%nbeads
    !
    !                    ! Applying PBCs
    !                    call minimg_one(atoms, i, j, b, r, vec)
    !
    !                    ! drops atoms outside (cutoff*rcut)-sphere
    !                    if (r > cutoff*rcut(idx_i,idx_j)) cycle
    !
    !                    vec = vec/r                       ! unit vector j -> i
    !
    !                    ! cut-off function
    !                    rtemp = exp(acut(idx_i,idx_j)*(r - rcut(idx_i,idx_j)))
    !                    theta = 1.0_dp / (1.0_dp + rtemp)
    !                    rtemp1 = acut(idx_i,idx_j)*rtemp*theta
    !
    !                    ! sigma_lp
    !                    rtemp = theta*exp(-pes_emt%eta2(idx_i,idx_j) * (r-betas0(idx_i,idx_j)))     ! sigma_ij*gamma1
    !                    sigma_ij(b,j) = sigma_ij(b,j) + rtemp
    !                    dtheta = (pes_emt%eta2(idx_j,idx_i) + rtemp1)*rtemp*vec
    !                    dsigma_ij_i(:,b,j) = dsigma_ij_i(:,b,j) + dtheta ! dsigma_i/dr_i
    !                    dsigma_ij_j(:,b,i,j) = -dtheta
    !
    !                    ! sigma_pl
    !                    rtemp = theta*exp(-pes_emt%eta2(idx_j,idx_i) * (r - betas0(idx_j,idx_i)) )     ! sigma_ij*gamma1
    !                    sigma_ji(b,i) = sigma_ji(b,j) + rtemp
    !                    dtheta = (pes_emt%eta2(idx_i,idx_j) + rtemp1)*rtemp*vec
    !                    dsigma_ji_j(:,b,i) = dsigma_ji_j(:,b,i) - dtheta
    !                    dsigma_ji_i(:,b,j,i) = dtheta ! dsigma_i/dr_i
    !
    !                    ! V_ij*gamma2*V_0
    !                    rtemp = theta*exp(-kappadbeta(idx_j,idx_i)*(r - betas0(idx_j,idx_i)))
    !                    V_ij(b) = V_ij(b) + rtemp
    !                    dtheta = (kappadbeta(idx_j,idx_i) + rtemp1)*rtemp*vec
    !                    dV_ij_j(:,b,i) = dV_ij_j(:,b,i) + dtheta
    !                    dV_ij_i(:,b,j) = dV_ij_i(:,b,j) - dtheta
    !
    !                    rtemp = theta*exp(-kappadbeta(idx_j,idx_i)*(r - betas0(idx_j,idx_i)))
    !                    V_ji(b) = V_ji(b) + rtemp
    !                    dtheta = (kappadbeta(idx_j,idx_i) + rtemp1)*rtemp*vec
    !                    dV_ji_j(:,b,i) = dV_ji_j(:,b,i) + dtheta
    !                    dV_ji_i(:,b,j) = dV_ji_i(:,b,j) - dtheta
    !                end do
    !            end do
    !        end do
    !
    !
    !        ! divide by cut-off scaling factors
    !        do i = 1, atoms%natoms
    !            do b = 1, atoms%nabeads
    !
    !                idx_i = atoms%idx(i)
    !
    !                sigma_ij(b,i) = sigma_ij(b,i) / igamma1(idx_i,idx_i)
    !                sigma_ij = sigma_ij*igamma1_i
    !                sigma_ji = sigma_ji*igamma1_j
    !                V_ij     =     V_ij*igamma2_i*pes_emt%v0(idx_i,idx_j)*chi_ij
    !                V_ji     =     V_ji*igamma2_j*pes_emt%v0(idx_j,idx_i)*chi_ji
    !
    !                dsigma_ii   = dsigma_ii  *igamma1_i
    !                dsigma_jj   = dsigma_jj  *igamma1_j
    !                dsigma_ij_i = dsigma_ij_i*igamma1_i
    !                dsigma_ij_j = dsigma_ij_j*igamma1_i
    !                dsigma_ji_j = dsigma_ji_j*igamma1_j
    !                dsigma_ji_i = dsigma_ji_i*igamma1_j
    !                dV_ii_i     =     dV_ii_i*igamma2_i*pes_emt%v0(idx_i,idx_j)
    !                dV_jj_j     =     dV_jj_j*igamma2_j*pes_emt%v0(idx_j,idx_i)
    !                dV_ij_i     =     dV_ij_i*igamma2_i*pes_emt%v0(idx_i,idx_j)*chi_ij
    !                dV_ij_j     =     dV_ij_j*igamma2_i*pes_emt%v0(idx_i,idx_j)*chi_ij
    !                dV_ji_i     =     dV_ji_i*igamma2_j*pes_emt%v0(idx_j,idx_i)*chi_ji
    !                dV_ji_j     =     dV_ji_j*igamma2_j*pes_emt%v0(idx_j,idx_i)*chi_ji
    !
    !                !-----------------------------NEUTRAL SPHERE RADIUS----------------------------
    !                s_i = sigma_ii + chi_ij*sigma_ij
    !                s_j = sigma_jj + chi_ji*sigma_ji
    !
    !                ! atom i
    !                ds_i_i = -dsigma_ii
    !                ds_i_i = ds_i_i - chi_ij*dsigma_ij_i
    !                ds_i_i = ds_i_i/(betaeta2_i*s_i)
    !                ds_j_i = -chi_ji*dsigma_ji_i/(betaeta2_j*s_j)
    !
    !                ! atom j
    !                ds_j_j = -dsigma_jj
    !                ds_j_j = ds_j_j - chi_ji*dsigma_ji_j
    !                ds_j_j = ds_j_j/(betaeta2_j*s_j)
    !                ds_i_j = -chi_ij*dsigma_ij_j/(betaeta2_i*s_i)
    !
    !            end do
    !        end do
    !    end do
    !    s_i = -log(s_i*twelfth)/betaeta2_i
    !    s_j = -log(s_j*twelfth)/betaeta2_j
    !
    !    !----------------------EMBEDDED ELECTRON DENSITY-------------------------------
    !
    !    rtemp = 0.5d0/bohr2ang - betaeta2_i     ! -eta_l
    !    dens(b,i) = dens(b,i) + pes_emt%n0(idx_i, idx_j)*exp(rtemp*s_i)
    !    rtemp = 0.5d0/bohr2ang - betaeta2_j     ! -eta_l
    !    dens(b,j) = dens(b,j) + pes_emt%n0(idx_j, idx_i)*exp(rtemp*s_j)
    !
    !
    !    !---------------------------COHESIVE FUNCTION-----------------------------------
    !
    !    Ecoh_i = ((1.0_dp + pes_emt%lambda(idx_i,idx_j)*s_i) * &
    !        exp(-pes_emt%lambda(idx_i,idx_j)*s_i) - 1.0_dp)*pes_emt%e0(idx_i,idx_j)
    !
    !    Ecoh_j = ((1.0_dp + pes_emt%lambda(idx_j,idx_i)*s_j) * &
    !        exp(-pes_emt%lambda(idx_j,idx_i)*s_j) - 1.0_dp)*pes_emt%e0(idx_j,idx_i)
    !
    !    ! dEcoh_l_l and dEcoh_p_l
    !    do i = 1, slab%n_atoms
    !
    !        dEcoh_l_l(1,i) = sum(s_l*exp(-pars_l(4)*s_l)*ds_l_l(1,i,:))
    !        dEcoh_l_l(2,i) = sum(s_l*exp(-pars_l(4)*s_l)*ds_l_l(2,i,:))
    !        dEcoh_l_l(3,i) = sum(s_l*exp(-pars_l(4)*s_l)*ds_l_l(3,i,:))
    !
    !        dEcoh_p_l(1,i) = sum(s_p*exp(-pars_p(4)*s_p)*ds_p_l(1,i,:))
    !        dEcoh_p_l(2,i) = sum(s_p*exp(-pars_p(4)*s_p)*ds_p_l(2,i,:))
    !        dEcoh_p_l(3,i) = sum(s_p*exp(-pars_p(4)*s_p)*ds_p_l(3,i,:))
    !
    !    end do
    !    dEcoh_l_l = pars_l(3)*pars_l(4)*pars_l(4)*dEcoh_l_l
    !    dEcoh_p_l = pars_p(3)*pars_p(4)*pars_p(4)*dEcoh_p_l
    !
    !    ! dEcoh_l_p and dEcoh_p_p
    !    do i = 1, teil%n_atoms
    !
    !        dEcoh_l_p(1,i) = sum(s_l*exp(-pars_l(4)*s_l)*ds_l_p(1,i,:))
    !        dEcoh_l_p(2,i) = sum(s_l*exp(-pars_l(4)*s_l)*ds_l_p(2,i,:))
    !        dEcoh_l_p(3,i) = sum(s_l*exp(-pars_l(4)*s_l)*ds_l_p(3,i,:))
    !
    !        dEcoh_p_p(1,i) = sum(s_p*exp(-pars_p(4)*s_p)*ds_p_p(1,i,:))
    !        dEcoh_p_p(2,i) = sum(s_p*exp(-pars_p(4)*s_p)*ds_p_p(2,i,:))
    !        dEcoh_p_p(3,i) = sum(s_p*exp(-pars_p(4)*s_p)*ds_p_p(3,i,:))
    !
    !    end do
    !    dEcoh_l_p = pars_l(3)*pars_l(4)*pars_l(4)*dEcoh_l_p
    !    dEcoh_p_p = pars_p(3)*pars_p(4)*pars_p(4)*dEcoh_p_p
    !
    !    !----------------REFERENCE PAIR POTENTIAL CONTRIBUTIONS------------------------
    !
    !    do i=1,slab%n_atoms
    !
    !        rtemp = exp(-pars_l(6)*s_l(i))
    !        vref_l = vref_l + rtemp
    !
    !        dvref_l_l = dvref_l_l + rtemp*ds_l_l(:,:,i)
    !        dvref_l_p = dvref_l_p + rtemp*ds_l_p(:,:,i)
    !
    !    end do
    !
    !    do i=1,teil%n_atoms
    !
    !        rtemp = exp(-pars_p(6)*s_p(i))
    !        vref_p = vref_p + rtemp
    !
    !        dvref_p_p = dvref_p_p + rtemp*ds_p_p(:,:,i)
    !        dvref_p_l = dvref_p_l + rtemp*ds_p_l(:,:,i)
    !
    !    end do
    !
    !    rtemp = 12.0d0 * pars_l(5)
    !    vref_l    =    vref_l*rtemp
    !    dvref_l_l = dvref_l_l*rtemp*pars_l(6)
    !    dvref_l_p = dvref_l_p*rtemp*pars_l(6)
    !
    !    rtemp = 12.0d0 * pars_p(5)
    !    vref_p    =    vref_p*rtemp
    !    dvref_p_l = dvref_p_l*rtemp*pars_p(6)
    !    dvref_p_p = dvref_p_p*rtemp*pars_p(6)
    !
    !
    !    !-------------------------------TOTAL ENERGY---------------------------------
    !
    !    Epot = Ecoh_l + Ecoh_p - V_ll - V_pp - 0.50d0*(V_lp + V_pl - vref_l - vref_p)
    !
    !    ! minus sign was taken into account in calculation of separate contributions
    !    slab%f = dEcoh_l_l + dEcoh_p_l - dV_ll_l &
    !        - 0.50d0*(dV_lp_l + dV_pl_l - dvref_l_l - dvref_p_l)
    !    teil%f = dEcoh_l_p + dEcoh_p_p - dV_pp_p &
    !        - 0.50d0*(dV_lp_p + dV_pl_p - dvref_l_p - dvref_p_p)
    !
    !
    !    deallocate(dvref_l_p, dvref_p_l, dvref_p_p, dvref_l_l)
    !    deallocate(dV_pl_p, dV_pl_l, dV_lp_p, dV_lp_l, dV_pp_p, dV_ll_l)
    !    deallocate(dEcoh_p_l, dEcoh_l_p, dEcoh_p_p, dEcoh_l_l)
    !    deallocate(ds_l_p, ds_p_l, ds_p_p, ds_l_l)
    !    deallocate(dsigma_pl_p, dsigma_pl_l, dsigma_lp_p, dsigma_lp_l)
    !    deallocate(dsigma_pp, dsigma_ll)
    !    deallocate( s_p,  s_l,  sigma_pl,  sigma_lp, sigma_pp,  sigma_ll)
    !end subroutine compute_emt



    subroutine compute_emt(atoms, flag)

        !   Calculates energy and forces with EMT potential

        type(universe), intent(inout) :: atoms
        integer, intent(in)           :: flag

        integer :: i, j, k, b, idx_i, idx_j

        real(dp), parameter :: cutoff = 1.5_dp

        real(dp) :: beta
        real(dp), dimension(atoms%ntypes, atoms%ntypes) :: igamma1, igamma2

        real(dp), dimension(atoms%ntypes,atoms%ntypes) :: betas0, betaeta2, kappadbeta, chi, &
            rcut, rr,  acut

        real(dp), dimension(3,atoms%ntypes,atoms%ntypes) ::  r3temp, x, rnn
        real(dp), dimension(3) :: vec
        real(dp) :: r, theta, rtemp!, rtemp1


        !       real(dp), dimension(3) :: dtheta




        real(dp), dimension(atoms%ntypes,atoms%ntypes,atoms%nbeads) :: V, V_ref
        real(dp), dimension(atoms%ntypes,atoms%ntypes,atoms%nbeads,atoms%natoms) ::  Ecoh
        real(dp), dimension(atoms%nbeads,atoms%natoms,atoms%natoms) :: sigma

        real(dp), dimension(atoms%nbeads,atoms%natoms,atoms%natoms) :: s

        real(dp), dimension(3) :: nneighs



        !----------------------VALUES OF FREQUENT USE ---------------------------------

        beta = select_beta('fcc')
        nneighs = select_nneighs('fcc')

        ! beta * s0
        betas0 = beta * pes_emt%s0
        print *, "betas0"
        print '(2f12.8)', betas0

        ! beta * eta2
        betaeta2 = beta * pes_emt%eta2

        print *, "betaeta2"
        print '(2f12.8)', betaeta2
        ! kappa / beta
        kappadbeta = pes_emt%kappa / beta

        ! 'coupling' parameters between p and l
        chi = pes_emt%n0 / transpose(pes_emt%n0) &
            * exp(0.5_dp/bohr2ang*(transpose(pes_emt%s0)-pes_emt%s0))

        print *, "chi"
        print  '(2f19.15)', chi

        ! Distances to the nearest, next-nearest and next-next-nearest neighbours
        rnn(1,:,:) = betas0
        rnn(2,:,:) = rnn(1,:,:) * sqrt2
        rnn(3,:,:) = rnn(1,:,:) * sqrt3

        !------------------------------------------------------------------------------
        !                                  CUT-OFF
        !                                  =======
        !------------------------------------------------------------------------------
        ! We use the distance to the next-next-nearest neighbours as cut-off.
        ! We only need one cut-off and we choose one of the lattice atoms since s0
        ! is usually larger for them.

        rcut = max(betas0, transpose(betas0)) * sqrt3
        print *, "rcut"
        print  '(2f19.15)', rcut

        !rcut = a_lat * sqrt3 * isqrt2
        rr = select_rr('fcc', rcut)
        acut = 9.21034037197618_dp/(rr - rcut) ! ln(10000)

        print *, "acut"
        print '(2f19.15)', acut

        do k = 1, 3
            x(k,:,:) = nneighs(k) * twelfth / (1.0_dp + exp(acut*(rnn(k,:,:)-rcut)))
            r3temp(k,:,:) = rnn(k,:,:) - betas0
        end do

        print *, "x"
        do i = 1, 3
            print '(2f19.15)', x(i,:,:)
        end do
        !-----------------------------------GAMMA--------------------------------------
        ! Gamma enforces the cut-off together with theta (see below)
        ! Gamma is defined as inverse.


        igamma1 = 0.0_dp
        igamma2 = 0.0_dp
        do k = 1, 3
            igamma1(:,:) = igamma1(:,:) + x(k,:,:)*exp(-pes_emt%eta2 * r3temp(k,:,:))
            igamma2(:,:) = igamma2(:,:) + x(k,:,:)*exp(-kappadbeta   * r3temp(k,:,:))
        end do
        igamma1 = 1.0_dp / transpose(igamma1)
        igamma2 = 1.0_dp / transpose(igamma2)

        print *, "igamma1"
        print '(2f19.15)', igamma1
        print *, "igamma2"
        print '(2f19.15)', igamma2




        !------------------------------------------------------------------------------
        !                          Sigma and Pair-wise Contributions
        !                          =================================
        !------------------------------------------------------------------------------

        ! initialize accumulators
        sigma       = 0.0_dp
        V           = 0.0_dp
        V_ref      = 0.0_dp
        s = 0.0_dp

        do i = 1, atoms%natoms-1
            do j = i+1, atoms%natoms


                idx_i = atoms%idx(i)
                idx_j = atoms%idx(j)

                if (i == j .or. atoms%pes(idx_i,idx_j) /= pes_id_emt) cycle

                do b = 1, atoms%nbeads


                    ! Applying PBCs
                    call minimg_one(atoms, i, j, b, r, vec)

                    ! drops atoms outside (cutoff*rcut)-sphere
                    if (r > cutoff*rcut(idx_i,idx_j)) cycle

                    vec = vec/r                       ! unit vector j -> i

                    ! cut-off function
                    rtemp = exp(acut(idx_i,idx_j)*(r - rcut(idx_i,idx_j)))
                    theta = 1.0_dp / (1.0_dp + rtemp)
                    !rtemp1 = acut(idx_i,idx_j)*rtemp*theta

                    ! sigma_1
                    !                    rtemp = theta*exp(-pes_emt%eta2(idx_i,idx_j) * (r - betas0(idx_i,idx_j)))   ! sigma_ij*gamma1
                    !                    sigma(idx_i,idx_j,b,i) = sigma(idx_i,idx_j,b,i) + rtemp
                    !
                    !                    rtemp = theta*exp(-pes_emt%eta2(idx_j,idx_i) * (r - betas0(idx_j,idx_i)))
                    !                    sigma(idx_j,idx_i,b,j) = sigma(idx_j,idx_i,b,j) + rtemp
                    rtemp = theta*exp(-pes_emt%eta2(idx_i,idx_j) * (r - betas0(idx_i,idx_j)))   ! sigma_ij*gamma1
                    sigma(b,i,j) = sigma(b,i,j) + rtemp

                    rtemp = theta*exp(-pes_emt%eta2(idx_j,idx_i) * (r - betas0(idx_j,idx_i)))
                    sigma(b,j,i) = sigma(b,j,i) + rtemp

                    ! sigma_2
                    rtemp = theta*exp(-kappadbeta(idx_i,idx_j)   * (r - betas0(idx_i,idx_j)))
                    V(idx_i,idx_j,b) = V(idx_i,idx_j,b) + rtemp

                    if (idx_i /= idx_j) then
                        rtemp = theta*exp(-kappadbeta(idx_j,idx_i)   * (r - betas0(idx_j,idx_i)))
                        V(idx_j,idx_i,b) = V(idx_j,idx_i,b) + rtemp
                    end if


                end do
            end do
        end do

        print *, "sigma"
        print '(3f15.8)', sigma

        print *, "V", sum(V)
        print '(2f15.8)', V(:,:,1)


        ! divide by cut-off scaling factors

        do b = 1, atoms%nbeads
            V(:,:,b) = V(:,:,b) * igamma2 * transpose(pes_emt%v0) * chi
            do i = 1, atoms%natoms
                idx_i = atoms%idx(i)
                do j = 1, atoms%natoms
                    idx_j = atoms%idx(j)
                    sigma(b,i,j) = sigma(b,i,j) * igamma1(idx_i,idx_j)
                end do
            end do
        end do


        print *, "sigma scaled"
        print '(3f15.8)', sigma


        print *, "V scaled"
        print '(2f15.8)', sum(V, dim=3)



         !-----------------------------NEUTRAL SPHERE RADIUS----------------------------
        do b = 1, atoms%nbeads
            do i = 1, atoms%natoms
                idx_i = atoms%idx(i)
                do j = 1, atoms%natoms
                    idx_j = atoms%idx(j)

                    s(b,i,j) = chi(idx_i,idx_j)*sigma(b,i,j)
                end do
            end do
        end do

        ! print *, "s"
        ! print '(2f15.8)',   s

        print *, "s"
        do i = 1, atoms%natoms
            print *, i
            print '(2f15.8)', s(1,i,:)
        end do


        do b = 1, atoms%nbeads
            do i = 1, atoms%natoms
                idx_i = atoms%idx(i)
                do j= 1, atoms%natoms
                    idx_j = atoms%idx(j)

                    if (s(b,i,j) > tolerance) then
                        s(b,i,j) = -log(s(b,i,j)*twelfth)/betaeta2(idx_i,idx_j)
                    end if
                end do
            end do
        end do




        print *, "s2"
        do i = 1, atoms%natoms
            print *, i
            print '(2f15.8)', s(1,i,:)
        end do
        print *, sum(sum(s, dim=1), dim=1)
         print *, sum(sum(s, dim=1), dim=2)

        stop
        !----------------------EMBEDDED ELECTRON DENSITY-------------------------------

    !        dens = 0.0_dp
    !        do i = 1, atoms%natoms
    !            do j = 1, atoms%natoms
    !                if (i == j) cycle
    !                idx_i = atoms%idx(i)
    !                idx_j = atoms%idx(j)
    !                do b = 1, atoms%nbeads
    !                    dens(b,i) = dens(b,i) + pes_emt%n0(idx_i,idx_j) * &
    !                        exp((0.5_dp/bohr2ang - betaeta2(idx_i,idx_j))*s(idx_j,idx_i,b,j))
    !                end do
    !            end do
    !        end do
    !
    !        print *, "dens"
    !        print '(2f15.8)', dens
    !
    !
    !                !---------------------------COHESIVE FUNCTION-----------------------------------
    !
    !        print '(2f15.8)',pes_emt%lambda
    !        print '(2f15.8)',pes_emt%e0
    !
    !        Ecoh = 0.0_dp
    !        do i = 1, atoms%natoms
    !            do j = 1, atoms%natoms
    !                if (i == j) cycle
    !                idx_i = atoms%idx(i)
    !                idx_j = atoms%idx(j)
    !                do b = 1, atoms%nbeads
    !                    Ecoh(idx_i,idx_j,b,i) = Ecoh(idx_i,idx_j,b,i) + &
    !                        ((1.0_dp + pes_emt%lambda(idx_i,idx_j)*s(idx_j,idx_i,b,j)) * &
    !                        exp(-pes_emt%lambda(idx_i,idx_j)*s(idx_j,idx_i,b,j)) - 1.0_dp) * &
    !                        pes_emt%e0(idx_i,idx_j)
    !                end do
    !            end do
    !        end do
    !        !sum((1.0d0 + pars_l(4)*s_l)*exp(-pars_l(4)*s_l) - 1.0d0)*pars_l(3)
    !
    !        print *, "Ecoh"
    !        print '(2f15.8)', sum(sum(sum(Ecoh, dim=1), dim=1), dim=1)
    !
    !         !----------------REFERENCE PAIR POTENTIAL CONTRIBUTIONS------------------------
    !
    !        do i=1,atoms%natoms
    !            do j = 1, atoms%natoms
    !                if (i == j) cycle
    !                idx_i = atoms%idx(i)
    !                idx_j = atoms%idx(j)
    !                do b = 1, atoms%nbeads
    !                    V_ref(idx_i,idx_j,b) = V_ref(idx_i,idx_j,b) + exp(-pes_emt%kappa(idx_i,idx_j)*s(idx_j,idx_i,b,j))
    !                end do
    !            end do
    !        end do
    !
    !        print *, "V_ref"
    !        print '(2f15.8)',  V_ref
    !
    !        do b = 1, atoms%nbeads
    !            V_ref(:,:,b) = V_ref(:,:,b) * 12.0_dp * pes_emt%v0
    !        end do
    !
    !        print *, "V_ref_2"
    !        print '(2f15.8)',  V_ref
    !
    !
    !        !-------------------------------TOTAL ENERGY---------------------------------
    !
    !        atoms%epot = sum(sum(sum(Ecoh, dim=1), dim=1), dim=2) - &
    !                0.5_dp*(sum(sum(V, dim=1), dim=1) - sum(sum(V_ref, dim=1), dim=1))
    !
    !        print *, "Epot", atoms%epot


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
        real(dp), intent(in) :: rcut(:,:)
        real(dp)             :: rr(size(rcut),size(rcut))

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
