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
                pes_emt%s0    (ntypes, ntypes)  &
                )

            pes_emt%eta2   = default_real
            pes_emt%n0     = default_real
            pes_emt%e0     = default_real
            pes_emt%lambda = default_real
            pes_emt%v0     = default_real
            pes_emt%kappa  = default_real
            pes_emt%s0     = default_real

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


    subroutine compute_emt(atoms, atom_i, atom_j, flag)

        !   Calculates energy and forces with EMT potential

        type(universe), intent(inout) :: atoms
        integer, intent(in)           :: atom_i, atom_j, flag

        integer :: i, j, idx_i, idx_j

        real(dp) :: beta
        real(dp) :: betas0_i, betaeta2_i, kappadbeta_i, chi_ij
        real(dp) :: betas0_j, betaeta2_j, kappadbeta_j, chi_ji
        real(dp) :: r, rcut, rr, acut, theta, rtemp, rtemp1
        real(dp) :: igamma1_i, igamma2_i, igamma1_j, igamma2_j
        real(dp) :: V_ji, V_ij, V_ii, V_jj, Ecoh_i, Ecoh_j, vref_i, vref_j

        real(dp), dimension(3) :: rnn_i, rnn_j         ! nnn-distances
        real(dp), dimension(3) :: nneighs
        real(dp), dimension(3) :: x_i, x_j, r3temp
        real(dp), dimension(3) :: dtheta

        real(dp) :: sigma_ii, sigma_ij, s_i
        real(dp) :: sigma_jj, sigma_ji, s_j

        real(dp) :: dsigma_ii(3), dsigma_ij_i(3), dsigma_ij_j(3)
        real(dp) :: dsigma_jj(3), dsigma_ji_i(3), dsigma_ji_j(3)

        real(dp) :: ds_i_i(3), ds_j_j(3), ds_i_j(3), ds_j_i(3)

        real(dp) :: dEcoh_i_i(3), dEcoh_i_j(3), dEcoh_j_i(3), dEcoh_j_j(3)


        real(dp), dimension(3) :: dV_ii_i, dV_jj_j, dV_ij_i, dV_ij_j, dV_ji_j, dV_ji_i
        real(dp), dimension(3) :: dvref_i_i, dvref_i_j, dvref_j_i, dvref_j_j




        !----------------------VALUES OF FREQUENT USE ---------------------------------

        idx_i = atoms%idx(i)
        idx_j = atoms%idx(j)

        beta = select_beta('fcc')
        nneighs = select_nneighs('fcc')

        ! beta * s0
        betas0_i = beta * pes_emt%s0(idx_i,idx_j)
        betas0_j = beta * pes_emt%s0(idx_j,idx_i)
        ! beta * eta2
        betaeta2_i = beta * pes_emt%eta2(idx_i,idx_j)
        betaeta2_j = beta * pes_emt%eta2(idx_i,idx_j)
        ! kappa / beta
        kappadbeta_i = pes_emt%kappa(idx_i,idx_j) / beta
        kappadbeta_j = pes_emt%kappa(idx_i,idx_j) / beta

        ! 'coupling' parameters between p and l
        chi_ij = pes_emt%n0(idx_j,idx_i) / pes_emt%n0(idx_i,idx_j) &
            *exp(0.5_dp/bohr2ang*(pes_emt%s0(idx_i,idx_j)-pes_emt%s0(idx_j,idx_i)))

        chi_ji = 1.0_dp / chi_ij

        ! Distances to the nearest, next-nearest and next-next-nearest neighbours
        rnn_i(1) = betas0_i
        rnn_i(2) = rnn_i(1) * sqrt2
        rnn_i(3) = rnn_i(1) * sqrt3
        rnn_j(1) = betas0_j
        rnn_j(2) = rnn_j(1) * sqrt2
        rnn_j(3) = rnn_j(1) * sqrt3

        !------------------------------------------------------------------------------
        !                                  CUT-OFF
        !                                  =======
        !------------------------------------------------------------------------------
        ! We use the distance to the next-next-nearest neighbours as cut-off.
        ! We only need one cut-off and we choose one of the lattice atoms since s0
        ! is usually larger for them.

        rcut = betas0_i * sqrt3
        !rcut = a_lat * sqrt3 * isqrt2
        rr = 4 * rcut / (sqrt3 + 2.0d0)
        acut = 9.21034037197618_dp/(rr -rcut) ! ln(10000)

        x_i = nneighs * twelfth / (1.0_dp + exp(acut*(rnn_i-rcut)))
        x_j = nneighs * twelfth / (1.0_dp + exp(acut*(rnn_j-rcut)))

        !-----------------------------------GAMMA--------------------------------------
        ! Gamma enforces the cut-off together with theta (see below)
        ! Gamma is defined as inverse.

        r3temp = rnn_i - betas0_i
        igamma1_i = 1.0_dp / sum(x_i*exp(-pes_emt%eta2(idx_i,idx_j) * r3temp))
        igamma2_i = 1.0_dp / sum(x_i*exp(-kappadbeta_i * r3temp))

        r3temp = rnn_j - betas0_j
        igamma1_j = 1.0_dp / sum(x_j*exp(-pes_emt%eta2(idx_j,idx_i) * r3temp))
        igamma2_j = 1.0_dp / sum(x_j*exp(-kappadbeta_j * r3temp))

        !------------------------------------------------------------------------------
        !                          Sigma and Pair-wise Contributions
        !                          =================================
        !------------------------------------------------------------------------------

        ! initialize accumulators
        sigma_ii    = 0.0_dp
        sigma_jj    = 0.0_dp
        sigma_ji    = 0.0_dp
        sigma_ij    = 0.0_dp
        V_ii        = 0.0_dp
        V_jj        = 0.0_dp
        V_ij        = 0.0_dp
        V_ji        = 0.0_dp
        vref_i      = 0.0_dp
        vref_j      = 0.0_dp
        dsigma_ii   = 0.0_dp
        dsigma_jj   = 0.0_dp
        dsigma_ij_i = 0.0_dp
        dsigma_ji_j = 0.0_dp
        dsigma_ji_i = 0.0_dp
        dsigma_ij_j = 0.0_dp
        dV_ii_i     = 0.0_dp
        dV_jj_j     = 0.0_dp
        dV_ij_i     = 0.0_dp
        dV_ij_j     = 0.0_dp
        dV_ji_i     = 0.0_dp
        dV_ji_j     = 0.0_dp
        dvref_i_i   = 0.0_dp
        dvref_i_j   = 0.0_dp
        dvref_j_i   = 0.0_dp
        dvref_j_j   = 0.0_dp


        ! Applying PBCs
        r3temp = teil%r(:,i) - slab%r(:,j)   ! distance vector
        r3temp = matmul(cell_imat, r3temp)       ! transform to direct coordinates

        r3temp(1) = r3temp(1) - Anint(r3temp(1))! imaging
        r3temp(2) = r3temp(2) - Anint(r3temp(2))
        r3temp(3) = r3temp(3) - Anint(r3temp(3))
        r3temp    = matmul(cell_mat, r3temp)    ! back to cartesian coordinates

        r =  sqrt(sum(r3temp**2))               ! distance

        ! drops atoms outside (cutoff*rcut)-sphere
        if (r > cutoff*rcut) cycle

        r3temp = r3temp/r                       ! unit vector j -> i

        ! cut-off function
        rtemp = exp(acut*(r - rcut))
        theta = 1.0d0 / (1.0d0 + rtemp)
        rtemp1 = acut*rtemp*theta

        ! sigma_lp
        rtemp = theta*exp(-pars_p(1) * (r - betas0_p) )     ! sigma_ij*gamma1
        sigma_lp(j) = sigma_lp(j) + rtemp
        dtheta = (pars_p(1) + rtemp1)*rtemp*r3temp
        dsigma_lp_l(:,j) = dsigma_lp_l(:,j) + dtheta ! dsigma_i/dr_i
        dsigma_lp_p(:,i,j) = - dtheta

        ! sigma_pl
        rtemp = theta*exp(-pars_l(1) * (r - betas0_l) )     ! sigma_ij*gamma1


        sigma_pl(i) = sigma_pl(i) + rtemp
        dtheta = (pars_l(1) + rtemp1)*rtemp*r3temp
        dsigma_pl_p(:,i) = dsigma_pl_p(:,i) - dtheta ! dsigma_i/dr_i
        dsigma_pl_l(:,j,i) = dtheta

        ! V_ij*gamma2*V_0
        rtemp = theta*exp(-kappadbeta_p*(r - betas0_p))
        V_lp = V_lp + rtemp
        dtheta = (kappadbeta_p + rtemp1)*rtemp*r3temp
        dV_lp_p(:,i) = dV_lp_p(:,i) + dtheta
        dV_lp_l(:,j) = dV_lp_l(:,j) - dtheta

        rtemp = theta*exp(-kappadbeta_l*(r - betas0_l))
        V_pl = V_pl + rtemp
        dtheta = (kappadbeta_l + rtemp1)*rtemp*r3temp
        dV_pl_p(:,i) = dV_pl_p(:,i) + dtheta
        dV_pl_l(:,j) = dV_pl_l(:,j) - dtheta


        ! divide by cut-off scaling factors
        sigma_ll = sigma_ll*igamma1l
        V_ll     =     V_ll*igamma2l*pars_l(5)
        sigma_pp = sigma_pp*igamma1p
        V_pp     =     V_pp*igamma2p*pars_p(5)
        sigma_lp = sigma_lp*igamma1l
        sigma_pl = sigma_pl*igamma1p
        V_lp     =     V_lp*igamma2l*pars_l(5)*chilp
        V_pl     =     V_pl*igamma2p*pars_p(5)*chipl

        dsigma_ll   = dsigma_ll  *igamma1l
        dsigma_pp   = dsigma_pp  *igamma1p
        dsigma_lp_l = dsigma_lp_l*igamma1l
        dsigma_lp_p = dsigma_lp_p*igamma1l
        dsigma_pl_p = dsigma_pl_p*igamma1p
        dsigma_pl_l = dsigma_pl_l*igamma1p
        dV_ll_l     =     dV_ll_l*igamma2l*pars_l(5)
        dV_pp_p     =     dV_pp_p*igamma2p*pars_p(5)
        dV_lp_l     =     dV_lp_l*igamma2l*pars_l(5)*chilp
        dV_lp_p     =     dV_lp_p*igamma2l*pars_l(5)*chilp
        dV_pl_l     =     dV_pl_l*igamma2p*pars_p(5)*chipl
        dV_pl_p     =     dV_pl_p*igamma2p*pars_p(5)*chipl

        !-----------------------------NEUTRAL SPHERE RADIUS----------------------------

        s_l = sigma_ll + chilp*sigma_lp
        s_p = sigma_pp + chipl*sigma_pl

        ds_l_l =-dsigma_ll
        do i = 1, slab%n_atoms

            ds_l_l(:,i,i) = ds_l_l(:,i,i) - chilp*dsigma_lp_l(:,i)
            ds_l_l(1,i,:) = ds_l_l(1,i,:)/(betaeta2_l*s_l)
            ds_l_l(2,i,:) = ds_l_l(2,i,:)/(betaeta2_l*s_l)
            ds_l_l(3,i,:) = ds_l_l(3,i,:)/(betaeta2_l*s_l)

            ds_p_l(1,i,:) =-chipl*dsigma_pl_l(1,i,:)/(betaeta2_p*s_p)
            ds_p_l(2,i,:) =-chipl*dsigma_pl_l(2,i,:)/(betaeta2_p*s_p)
            ds_p_l(3,i,:) =-chipl*dsigma_pl_l(3,i,:)/(betaeta2_p*s_p)

        end do

        ds_p_p =-dsigma_pp
        do i = 1, teil%n_atoms

            ds_p_p(:,i,i) = ds_p_p(:,i,i) - chipl*dsigma_pl_p(:,i)
            ds_p_p(1,i,:) = ds_p_p(1,i,:)/(betaeta2_p*s_p)
            ds_p_p(2,i,:) = ds_p_p(2,i,:)/(betaeta2_p*s_p)
            ds_p_p(3,i,:) = ds_p_p(3,i,:)/(betaeta2_p*s_p)

            ds_l_p(1,i,:) =-chilp*dsigma_lp_p(1,i,:)/(betaeta2_l*s_l)
            ds_l_p(2,i,:) =-chilp*dsigma_lp_p(2,i,:)/(betaeta2_l*s_l)
            ds_l_p(3,i,:) =-chilp*dsigma_lp_p(3,i,:)/(betaeta2_l*s_l)

        end do

        s_l = -log(s_l*twelfth)/betaeta2_l
        s_p = -log(s_p*twelfth)/betaeta2_p

        !----------------------EMBEDDED ELECTRON DENSITY-------------------------------

        rtemp = 0.5d0/bohr2ang - betaeta2_l     ! -eta_l
        slab%dens = pars_l(2)*exp(rtemp*s_l)
        rtemp = 0.5d0/bohr2ang - betaeta2_p     ! -eta_l
        teil%dens = pars_p(2)*exp(rtemp*s_p)


        !---------------------------COHESIVE FUNCTION-----------------------------------

        Ecoh_l = sum((1.0d0 + pars_l(4)*s_l)*exp(-pars_l(4)*s_l) - 1.0d0)*pars_l(3)
        Ecoh_p = sum((1.0d0 + pars_p(4)*s_p)*exp(-pars_p(4)*s_p) - 1.0d0)*pars_p(3)

        ! dEcoh_l_l and dEcoh_p_l
        do i = 1, slab%n_atoms

            dEcoh_l_l(1,i) = sum(s_l*exp(-pars_l(4)*s_l)*ds_l_l(1,i,:))
            dEcoh_l_l(2,i) = sum(s_l*exp(-pars_l(4)*s_l)*ds_l_l(2,i,:))
            dEcoh_l_l(3,i) = sum(s_l*exp(-pars_l(4)*s_l)*ds_l_l(3,i,:))

            dEcoh_p_l(1,i) = sum(s_p*exp(-pars_p(4)*s_p)*ds_p_l(1,i,:))
            dEcoh_p_l(2,i) = sum(s_p*exp(-pars_p(4)*s_p)*ds_p_l(2,i,:))
            dEcoh_p_l(3,i) = sum(s_p*exp(-pars_p(4)*s_p)*ds_p_l(3,i,:))

        end do
        dEcoh_l_l = pars_l(3)*pars_l(4)*pars_l(4)*dEcoh_l_l
        dEcoh_p_l = pars_p(3)*pars_p(4)*pars_p(4)*dEcoh_p_l

        ! dEcoh_l_p and dEcoh_p_p
        do i = 1, teil%n_atoms

            dEcoh_l_p(1,i) = sum(s_l*exp(-pars_l(4)*s_l)*ds_l_p(1,i,:))
            dEcoh_l_p(2,i) = sum(s_l*exp(-pars_l(4)*s_l)*ds_l_p(2,i,:))
            dEcoh_l_p(3,i) = sum(s_l*exp(-pars_l(4)*s_l)*ds_l_p(3,i,:))

            dEcoh_p_p(1,i) = sum(s_p*exp(-pars_p(4)*s_p)*ds_p_p(1,i,:))
            dEcoh_p_p(2,i) = sum(s_p*exp(-pars_p(4)*s_p)*ds_p_p(2,i,:))
            dEcoh_p_p(3,i) = sum(s_p*exp(-pars_p(4)*s_p)*ds_p_p(3,i,:))

        end do
        dEcoh_l_p = pars_l(3)*pars_l(4)*pars_l(4)*dEcoh_l_p
        dEcoh_p_p = pars_p(3)*pars_p(4)*pars_p(4)*dEcoh_p_p

        !----------------REFERENCE PAIR POTENTIAL CONTRIBUTIONS------------------------

        do i=1,slab%n_atoms

            rtemp = exp(-pars_l(6)*s_l(i))
            vref_l = vref_l + rtemp

            dvref_l_l = dvref_l_l + rtemp*ds_l_l(:,:,i)
            dvref_l_p = dvref_l_p + rtemp*ds_l_p(:,:,i)

        end do

        do i=1,teil%n_atoms

            rtemp = exp(-pars_p(6)*s_p(i))
            vref_p = vref_p + rtemp

            dvref_p_p = dvref_p_p + rtemp*ds_p_p(:,:,i)
            dvref_p_l = dvref_p_l + rtemp*ds_p_l(:,:,i)

        end do

        rtemp = 12.0d0 * pars_l(5)
        vref_l    =    vref_l*rtemp
        dvref_l_l = dvref_l_l*rtemp*pars_l(6)
        dvref_l_p = dvref_l_p*rtemp*pars_l(6)

        rtemp = 12.0d0 * pars_p(5)
        vref_p    =    vref_p*rtemp
        dvref_p_l = dvref_p_l*rtemp*pars_p(6)
        dvref_p_p = dvref_p_p*rtemp*pars_p(6)


        !-------------------------------TOTAL ENERGY---------------------------------

        Epot = Ecoh_l + Ecoh_p - V_ll - V_pp - 0.50d0*(V_lp + V_pl - vref_l - vref_p)

        ! minus sign was taken into account in calculation of separate contributions
        slab%f = dEcoh_l_l + dEcoh_p_l - dV_ll_l &
            - 0.50d0*(dV_lp_l + dV_pl_l - dvref_l_l - dvref_p_l)
        teil%f = dEcoh_l_p + dEcoh_p_p - dV_pp_p &
            - 0.50d0*(dV_lp_p + dV_pl_p - dvref_l_p - dvref_p_p)


        deallocate(dvref_l_p, dvref_p_l, dvref_p_p, dvref_l_l)
        deallocate(dV_pl_p, dV_pl_l, dV_lp_p, dV_lp_l, dV_pp_p, dV_ll_l)
        deallocate(dEcoh_p_l, dEcoh_l_p, dEcoh_p_p, dEcoh_l_l)
        deallocate(ds_l_p, ds_p_l, ds_p_p, ds_l_l)
        deallocate(dsigma_pl_p, dsigma_pl_l, dsigma_lp_p, dsigma_lp_l)
        deallocate(dsigma_pp, dsigma_ll)
        deallocate( s_p,  s_l,  sigma_pl,  sigma_lp, sigma_pp,  sigma_ll)
    end subroutine compute_emt



    function select_nneighs(geom) result(neighs)

        character(len=3), intent(in) :: geom
        real(dp) :: neighs(3)

        select case (geom)
            case ('fcc')
                neighs = [12, 6, 24]
            case default
                print *, "Error in select_nneighs(): unknown lattice structure", geom ; call abort
        end select

    end function select_nneighs



    real(dp) function select_beta(geom) result(beta)

        character(len=3), intent(in) :: geom

        select case (geom)
            case ('fcc')
                beta = (16*pi/3)**(1.0_dp/3)*isqrt2
            case default
                print *, "Error in select_beta(): unknown lattice structure", geom ; call abort
        end select

    end function select_beta

end module pes_emt_mod
