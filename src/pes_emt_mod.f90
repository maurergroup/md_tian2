!############################################################################
! This routine is part of
! md_tian2 (Molecular Dynamics Tian Xia 2)
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

module pes_emt_mod

    use constants
    use universe_mod
    use useful_things, only : split_string, lower_case


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
            print *, err // "interaction must be defined via 'proj' and 'latt' keywords"
            stop
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

        type(universe), intent(inout) :: atoms
        integer, intent(in)           :: flag

        integer :: i, nEMT

        ! determine number of species which are interacting via EMT potential
        nEMT = 0
        do i = 1, atoms%ntypes
            if (atoms%pes(i,i) == pes_id_emt) nEMT = nEMT + 1
        end do

        select case (nEMT)

            case (1)
                call compute_emt_nspecies(atoms, flag)  ! TODO: code a one-species version

            case (2)
                call compute_emt_2species(atoms, flag)

            case (3:)
                call compute_emt_nspecies(atoms, flag)

            case default
                print *, "Error in EMT energy subroutine: EMT compute requested, though there is no EMT PES speciefied"
                stop

        end select


    end subroutine compute_emt



    subroutine compute_emt_2species(atoms, flag)

        type(universe), intent(inout) :: atoms
        integer, intent(in)           :: flag

        real(dp), parameter :: cutoff = 1.5_dp

        integer :: i, j, b, type1, type2, idx_i, idx_j

        real(dp) :: betas0_1, betaeta2_1, kappadbeta_1, chi_21
        real(dp) :: betas0_2, betaeta2_2, kappadbeta_2, chi_12

        real(dp) :: rcut, rr, acut, theta, rtemp, rtemp1, rtemp2
        real(dp) :: igamma1_1, igamma2_1, igamma1_2, igamma2_2, beta

        real(dp), dimension(3) :: rnn_2, rnn_1         ! nnn-distances
        real(dp), dimension(3) :: x_1, x_2
        real(dp), dimension(3) :: dtheta, nneighs, r3temp

        real(dp), dimension(atoms%nbeads) :: V_12, V_21, V_22, V_11, Ecoh_2, Ecoh_1
        real(dp), dimension(atoms%nbeads) :: vref_2, vref_1, rijmag, rtemp_beads
        real(dp), dimension(3, atoms%nbeads) :: rij

        real(dp), dimension(atoms%nbeads,atoms%natoms) :: sigma_22, sigma_21, s_2
        real(dp), dimension(atoms%nbeads,atoms%natoms) :: sigma_11, sigma_12, s_1
        real(dp), dimension(atoms%nbeads,atoms%natoms) :: exp_lambda_s1, exp_lambda_s2

        real(dp), dimension(3,atoms%nbeads,atoms%natoms) :: dsigma_21_2, dsigma_12_1
        real(dp), dimension(3,atoms%nbeads,atoms%natoms) :: dEcoh_2_2, dEcoh_2_1, dEcoh_1_2, dEcoh_1_1
        real(dp), dimension(3,atoms%nbeads,atoms%natoms) :: dV_22_2, dV_11_1, dV_21_2, dV_21_1, dV_12_1, dV_12_2
        real(dp), dimension(3,atoms%nbeads,atoms%natoms) :: dvref_2_2, dvref_2_1, dvref_1_2, dvref_1_1

        real(dp), dimension(3,atoms%nbeads,atoms%natoms,atoms%natoms) :: dsigma_22, dsigma_21_1
        real(dp), dimension(3,atoms%nbeads,atoms%natoms,atoms%natoms) :: dsigma_11, dsigma_12_2
        real(dp), dimension(3,atoms%nbeads,atoms%natoms,atoms%natoms) :: ds_2_2, ds_1_1, ds_2_1, ds_1_2

        !----------------------VALUES OF FREQUENT USE ---------------------------------

        ! find the two EMT indices
        type1 = default_int
        do i = 1, atoms%ntypes
            if (atoms%pes(i,i) == pes_id_emt) then
                if (type1 == default_int) then
                    type1  = i
                else
                    type2  = i
                end if
            end if
        end do

        beta = select_beta('fcc')
        nneighs = select_nneighs('fcc')

        ! beta * s0
        betas0_1 = beta * pes_emt%s0(type1)
        betas0_2 = beta * pes_emt%s0(type2)

        ! beta * eta2
        betaeta2_1 = beta * pes_emt%eta2(type1)
        betaeta2_2 = beta * pes_emt%eta2(type2)

        ! print *, betaeta2_1, betaeta2_2

        ! kappa / beta
        kappadbeta_1 = pes_emt%kappa(type1) / beta
        kappadbeta_2 = pes_emt%kappa(type2) / beta

        ! 'coupling' parameters between 1 and 2

        chi_21 = pes_emt%n0(type1) / pes_emt%n0(type2) *exp(0.5/bohr2ang*(pes_emt%s0(type2)-pes_emt%s0(type1)))
        chi_12 = 1.0 / chi_21

        !print *,chi_12, chi_21
        !stop

        ! Distances to the nearest, next-nearest and next-next-nearest neighbours
        rnn_1(1) = betas0_1
        rnn_1(2) = rnn_1(1) * sqrt2
        rnn_1(3) = rnn_1(1) * sqrt3
        rnn_2(1) = betas0_2
        rnn_2(2) = rnn_2(1) * sqrt2
        rnn_2(3) = rnn_2(1) * sqrt3

        !------------------------------------------------------------------------------
        !                                  CUT-OFF
        !                                  =======
        !------------------------------------------------------------------------------
        ! We use the distance to the next-next-nearest neighbours as cut-off.
        ! We only need one cut-off and we choose the one of the lattice atoms since s0
        ! is usually larger for them.

        rcut = max(betas0_1, betas0_2) * sqrt3
        !rcut = a_lat * sqrt3 * isqrt2
        rr = 4.0 * rcut / (sqrt3 + 2)
        acut = 9.21034037197618_dp/(rr -rcut) ! ln(10000)

        x_1 = nneighs * twelfth / (1 + exp(acut*(rnn_1-rcut)))
        x_2 = nneighs * twelfth / (1 + exp(acut*(rnn_2-rcut)))

        !-----------------------------------GAMMA--------------------------------------
        ! Gamma enforces the cut-off together with theta (see below)
        ! Gamma is defined as inverse.

        r3temp = rnn_1 - betas0_1
        igamma1_1 = 1.0 / sum(x_1*exp(-pes_emt%eta2(type1) * r3temp))
        igamma2_1 = 1.0 / sum(x_1*exp(-kappadbeta_1 * r3temp))

        r3temp = rnn_2 - betas0_2
        igamma1_2 = 1.0 / sum(x_2*exp(-pes_emt%eta2(type2) * r3temp))
        igamma2_2 = 1.0 / sum(x_2*exp(-kappadbeta_2 * r3temp))



        !------------------------------------------------------------------------------
        !                          Sigma and Pair-wise Contributions
        !                          =================================
        !------------------------------------------------------------------------------


        ! initialize accumulators
        sigma_22    = 0.0
        sigma_11    = 0.0
        sigma_12    = 0.0
        sigma_21    = 0.0
        dsigma_22   = 0.0
        dsigma_11   = 0.0
        dsigma_21_2 = 0.0
        dsigma_12_1 = 0.0
        dsigma_12_2 = 0.0
        dsigma_21_1 = 0.0

        s_1         = 0.0
        s_2         = 0.0
        ds_1_1      = 0.0
        ds_1_2      = 0.0
        ds_2_1      = 0.0
        ds_2_2      = 0.0

        Ecoh_1      = 0.0
        Ecoh_2      = 0.0
        dEcoh_1_1   = 0.0
        dEcoh_1_2   = 0.0
        dEcoh_2_1   = 0.0
        dEcoh_2_2   = 0.0

        V_22        = 0.0
        V_11        = 0.0
        V_21        = 0.0
        V_12        = 0.0
        dV_22_2     = 0.0
        dV_11_1     = 0.0
        dV_21_2     = 0.0
        dV_21_1     = 0.0
        dV_12_2     = 0.0
        dV_12_1     = 0.0

        vref_2      = 0.0
        vref_1      = 0.0
        dvref_2_2   = 0.0
        dvref_2_1   = 0.0
        dvref_1_2   = 0.0
        dvref_1_1   = 0.0

        do i = 1, atoms%natoms-1
            idx_i = atoms%idx(i)

            do j = i+1, atoms%natoms
                idx_j = atoms%idx(j)

                if (atoms%pes(idx_i,idx_j) /= pes_id_emt) cycle

                rijmag = atoms%distances(:,i,j)
                rij    = atoms%vectors(:,:,i,j)

                if (idx_i == type2 .and. idx_j == type2) then

                    do b = 1, atoms%nbeads

                        if (rijmag(b) > cutoff*rcut) cycle

                        ! cut-off function
                        rtemp = exp(acut*(rijmag(b) - rcut))
                        theta = 1.0 / (1 + rtemp)
                        rtemp1 = acut*rtemp*theta

                        rtemp = theta*exp(-pes_emt%eta2(type2)*(rijmag(b) - betas0_2))    ! sigma_ij*gamma1
                        sigma_22(b,i) = sigma_22(b,i) + rtemp
                        sigma_22(b,j) = sigma_22(b,j) + rtemp

                        dtheta = (pes_emt%eta2(type2) + rtemp1)*rtemp * rij(:,b)/rijmag(b)
                        dsigma_22(:,b,i,i) = dsigma_22(:,b,i,i) - dtheta    ! dsigma_i/dr_i
                        dsigma_22(:,b,j,j) = dsigma_22(:,b,j,j) + dtheta
                        dsigma_22(:,b,j,i) = dtheta                       ! dsigma_i/dr_j
                        dsigma_22(:,b,i,j) =-dtheta                       ! dsigma_j/dr_i

                        rtemp = theta*exp(-kappadbeta_2*(rijmag(b) - betas0_2)) ! V_ij*gamma2*V_0
                        V_22(b) = V_22(b) + rtemp

                        dtheta = (kappadbeta_2 + rtemp1)*rtemp * rij(:,b)/rijmag(b)
                        dV_22_2(:,b,i) = dV_22_2(:,b,i) + dtheta
                        dV_22_2(:,b,j) = dV_22_2(:,b,j) - dtheta

                    end do


                else if (idx_i == type1 .and. idx_j == type1) then

                    do b = 1, atoms%nbeads

                        if (rijmag(b) > cutoff*rcut) cycle

                        ! cut-off function
                        rtemp = exp(acut*(rijmag(b) - rcut))
                        theta = 1.0 / (1 + rtemp)
                        rtemp1 = acut*rtemp*theta

                        rtemp = theta*exp(-pes_emt%eta2(type1)*(rijmag(b) - betas0_1))     ! sigma_ij*gamma1
                        sigma_11(b,i) = sigma_11(b,i) + rtemp
                        sigma_11(b,j) = sigma_11(b,j) + rtemp

                        dtheta = (pes_emt%eta2(type2) + rtemp1)*rtemp * rij(:,b)/rijmag(b)
                        dsigma_11(:,b,i,i) = dsigma_11(:,b,i,i) - dtheta    ! dsigma_i/dr_i
                        dsigma_11(:,b,j,j) = dsigma_11(:,b,j,j) + dtheta
                        dsigma_11(:,b,j,i) = dtheta                       ! dsigma_i/dr_j
                        dsigma_11(:,b,i,j) =-dtheta                       ! dsigma_j/dr_i

                        rtemp = theta*exp(-kappadbeta_1 * (rijmag(b) - betas0_1))   ! V_ij*gamma2*V_0
                        V_11(b) = V_11(b) + rtemp

                        dtheta = (kappadbeta_1 + rtemp1)*rtemp * rij(:,b)/rijmag(b)
                        dV_11_1(:,b,i) = dV_11_1(:,b,i) + dtheta
                        dV_11_1(:,b,j) = dV_11_1(:,b,j) - dtheta

                    end do


                else if (idx_i /= idx_j) then

                    do b = 1, atoms%nbeads

                        if (rijmag(b) > cutoff*rcut) cycle

                        ! cut-off function
                        rtemp = exp(acut*(rijmag(b) - rcut))
                        theta = 1.0 / (1 + rtemp)
                        rtemp1 = acut*rtemp*theta

                        ! sigma_21
                        rtemp = theta*exp(-pes_emt%eta2(type1) * (rijmag(b) - betas0_1) )     ! sigma_ij*gamma1
                        sigma_21(b,j) = sigma_21(b,j) + rtemp
                        dtheta = (pes_emt%eta2(type1) + rtemp1)*rtemp * rij(:,b)/rijmag(b)
                        dsigma_21_2(:,b,j) = dsigma_21_2(:,b,j) + dtheta ! dsigma_i/dr_i
                        dsigma_21_1(:,b,i,j) = -dtheta

                        ! sigma_12
                        rtemp = theta*exp(-pes_emt%eta2(type2) * (rijmag(b) - betas0_2) )     ! sigma_ij*gamma1
                        sigma_12(b,i) = sigma_12(b,i) + rtemp
                        dtheta = (pes_emt%eta2(type2) + rtemp1)*rtemp * rij(:,b)/rijmag(b)
                        dsigma_12_1(:,b,i) = dsigma_12_1(:,b,i) - dtheta ! dsigma_i/dr_i
                        dsigma_12_2(:,b,j,i) = dtheta

                        ! V_ij*gamma2*V_0
                        rtemp = theta*exp(-kappadbeta_1*(rijmag(b) - betas0_1))
                        V_21(b) = V_21(b) + rtemp
                        dtheta = (kappadbeta_1 + rtemp1)*rtemp * rij(:,b)/rijmag(b)
                        dV_21_1(:,b,i) = dV_21_1(:,b,i) + dtheta
                        dV_21_2(:,b,j) = dV_21_2(:,b,j) - dtheta

                        rtemp = theta*exp(-kappadbeta_2*(rijmag(b) - betas0_2))
                        V_12(b) = V_12(b) + rtemp
                        dtheta = (kappadbeta_2 + rtemp1)*rtemp * rij(:,b)/rijmag(b)
                        dV_12_1(:,b,i) = dV_12_1(:,b,i) + dtheta
                        dV_12_2(:,b,j) = dV_12_2(:,b,j) - dtheta

                    end do

                end if
            end do
        end do



        ! divide by cut-off scaling factors
        sigma_22 = sigma_22*igamma1_2
        V_22     =     V_22*igamma2_2*pes_emt%v0(type2)
        sigma_11 = sigma_11*igamma1_1
        V_11     =     V_11*igamma2_1*pes_emt%v0(type1)
        sigma_21 = sigma_21*igamma1_2
        sigma_12 = sigma_12*igamma1_1
        V_21     =     V_21*igamma2_2*pes_emt%v0(type2)*chi_21
        V_12     =     V_12*igamma2_1*pes_emt%v0(type1)*chi_12

        dsigma_22   = dsigma_22  *igamma1_2
        dsigma_11   = dsigma_11  *igamma1_1
        dsigma_21_2 = dsigma_21_2*igamma1_2
        dsigma_21_1 = dsigma_21_1*igamma1_2
        dsigma_12_1 = dsigma_12_1*igamma1_1
        dsigma_12_2 = dsigma_12_2*igamma1_1
        dV_22_2     =     dV_22_2*igamma2_2*pes_emt%v0(type2)
        dV_11_1     =     dV_11_1*igamma2_1*pes_emt%v0(type1)
        dV_21_2     =     dV_21_2*igamma2_2*pes_emt%v0(type2)*chi_21
        dV_21_1     =     dV_21_1*igamma2_2*pes_emt%v0(type2)*chi_21
        dV_12_2     =     dV_12_2*igamma2_1*pes_emt%v0(type1)*chi_12
        dV_12_1     =     dV_12_1*igamma2_1*pes_emt%v0(type1)*chi_12

        !-----------------------------NEUTRAL SPHERE RADIUS----------------------------

        ! is being used for division
        s_1 = max(1e-30, sigma_11 + chi_12*sigma_12)
        s_2 = max(1e-30, sigma_22 + chi_21*sigma_21)

!        s_1 = sigma_11 + chi_12*sigma_12
!        s_2 = sigma_22 + chi_21*sigma_21

        ds_1_1 = -dsigma_11
        ds_2_2 = -dsigma_22

        do i = 1, atoms%natoms
            if (atoms%idx(i) == type2) then

                ds_2_2(:,:,i,i) = ds_2_2(:,:,i,i) - chi_21*dsigma_21_2(:,:,i)
                do j = 1, atoms%natoms
                    if (atoms%idx(j) == type2) then

                        do b = 1, atoms%nbeads
                            ds_2_2(:,b,i,j) = ds_2_2(:,b,i,j)/(betaeta2_2*s_2(b,j))
                        end do

                    else if (atoms%idx(j) == type1) then

                        do b = 1, atoms%nbeads
                            ds_1_2(:,b,i,j) = -chi_12*dsigma_12_2(:,b,i,j)/(betaeta2_1*s_1(b,j))
                        end do

                    end if
                end do

            else if (atoms%idx(i) == type1) then

                ds_1_1(:,:,i,i) = ds_1_1(:,:,i,i) - chi_12*dsigma_12_1(:,:,i)
                do j = 1, atoms%natoms
                    if (atoms%idx(j) == type1) then

                        do b = 1, atoms%nbeads
                            ds_1_1(:,b,i,j) = ds_1_1(:,b,i,j)/(betaeta2_1*s_1(b,j))
                        end do

                    else if (atoms%idx(j) == type2) then

                        do b = 1, atoms%nbeads
                            ds_2_1(:,b,i,j) = -chi_21*dsigma_21_1(:,b,i,j)/(betaeta2_2*s_2(b,j))
                        end do

                    end if
                end do

            end if
        end do

        s_1 = -log(max(tolerance,s_1)*twelfth)/betaeta2_1
        s_2 = -log(max(tolerance,s_2)*twelfth)/betaeta2_2

        !----------COHESIVE FUNCTION AND EMBEDDED ELECTRON DENSITY--------------------------

        rtemp1 = 0.5/bohr2ang - betaeta2_1
        rtemp2 = 0.5/bohr2ang - betaeta2_2

        do i = 1, atoms%natoms
            if (atoms%idx(i) == type1) then

                dens(:,i) = pes_emt%n0(type1)*exp(rtemp1*s_1(:,i))
                exp_lambda_s1(:,i) = exp(-pes_emt%lambda(type1)*s_1(:,i))
                Ecoh_1 = Ecoh_1 + ((1 + pes_emt%lambda(type1)*s_1(:,i)) * exp_lambda_s1(:,i) - 1)

            else if (atoms%idx(i) == type2) then

                dens(:,i) = pes_emt%n0(type2)*exp(rtemp2*s_2(:,i))
                exp_lambda_s2(:,i) = exp(-pes_emt%lambda(type2)*s_2(:,i))
                Ecoh_2 = Ecoh_2 + ((1 + pes_emt%lambda(type2)*s_2(:,i)) * exp_lambda_s2(:,i) - 1)

            end if
        end do

        Ecoh_1 = Ecoh_1*pes_emt%e0(type1)
        Ecoh_2 = Ecoh_2*pes_emt%e0(type2)

        do i = 1, atoms%natoms
            if (atoms%idx(i) == type2) then

                do j = 1, atoms%natoms
                    if (atoms%idx(j) == type2) then

                        do b = 1, atoms%nbeads
                            dEcoh_2_2(:,b,i) = dEcoh_2_2(:,b,i) + s_2(b,j)*exp_lambda_s2(b,j)*ds_2_2(:,b,i,j)
                        end do

                    else if (atoms%idx(j) == type1) then

                        do b = 1, atoms%nbeads
                            dEcoh_1_2(:,b,i) = dEcoh_1_2(:,b,i) + s_1(b,j)*exp_lambda_s1(b,j)*ds_1_2(:,b,i,j)
                        end do

                    end if
                end do

            else if (atoms%idx(i) == type1) then

                do j = 1, atoms%natoms
                    if (atoms%idx(j) == type2) then

                        do b = 1, atoms%nbeads
                            dEcoh_2_1(:,b,i) = dEcoh_2_1(:,b,i) + s_2(b,j)*exp_lambda_s2(b,j)*ds_2_1(:,b,i,j)
                        end do

                    else if (atoms%idx(j) == type1) then

                        do b = 1, atoms%nbeads
                            dEcoh_1_1(:,b,i) = dEcoh_1_1(:,b,i) + s_1(b,j)*exp_lambda_s1(b,j)*ds_1_1(:,b,i,j)
                        end do

                    end if
                end do

            end if
        end do

        dEcoh_2_2 = pes_emt%e0(type2)*pes_emt%lambda(type2)*pes_emt%lambda(type2)*dEcoh_2_2
        dEcoh_1_2 = pes_emt%e0(type1)*pes_emt%lambda(type1)*pes_emt%lambda(type1)*dEcoh_1_2
        dEcoh_2_1 = pes_emt%e0(type2)*pes_emt%lambda(type2)*pes_emt%lambda(type2)*dEcoh_2_1
        dEcoh_1_1 = pes_emt%e0(type1)*pes_emt%lambda(type1)*pes_emt%lambda(type1)*dEcoh_1_1

        !----------------REFERENCE PAIR POTENTIAL CONTRIBUTIONS------------------------

        do i = 1, atoms%natoms
            if (atoms%idx(i) == type2) then

                rtemp_beads = exp(-pes_emt%kappa(type2)*s_2(:,i))
                vref_2 = vref_2 + rtemp_beads

                do j = 1, atoms%natoms
                    if (atoms%idx(j) == type2) then

                        do b = 1, atoms%nbeads
                            dvref_2_2(:,b,j) = dvref_2_2(:,b,j) + rtemp_beads(b)*ds_2_2(:,b,j,i)
                        end do

                    else if (atoms%idx(j) == type1) then

                        do b = 1, atoms%nbeads
                            dvref_2_1(:,b,j) = dvref_2_1(:,b,j) + rtemp_beads(b)*ds_2_1(:,b,j,i)
                        end do

                    end if
                end do

            else if (atoms%idx(i) == type1) then

                rtemp_beads = exp(-pes_emt%kappa(type1)*s_1(:,i))
                vref_1 = vref_1 + rtemp_beads

                do j = 1, atoms%natoms
                    if (atoms%idx(j) == type1) then

                        do b = 1, atoms%nbeads
                            dvref_1_1(:,b,j) = dvref_1_1(:,b,j) + rtemp_beads(b)*ds_1_1(:,b,j,i)
                        end do

                    else if (atoms%idx(j) == type2) then

                        do b = 1, atoms%nbeads
                            dvref_1_2(:,b,j) = dvref_1_2(:,b,j) + rtemp_beads(b)*ds_1_2(:,b,j,i)
                        end do

                    end if
                end do

            end if
        end do

        rtemp     = 12 * pes_emt%v0(type2)
        vref_2    =    vref_2*rtemp
        dvref_2_2 = dvref_2_2*rtemp*pes_emt%kappa(type2)
        dvref_2_1 = dvref_2_1*rtemp*pes_emt%kappa(type2)

        rtemp     = 12 * pes_emt%v0(type1)
        vref_1    =    vref_1*rtemp
        dvref_1_2 = dvref_1_2*rtemp*pes_emt%kappa(type1)
        dvref_1_1 = dvref_1_1*rtemp*pes_emt%kappa(type1)

        !-------------------------------TOTAL ENERGY---------------------------------

        atoms%epot = atoms%epot + (Ecoh_2 + Ecoh_1 - V_22 - V_11 - 0.5*(V_21 + V_12 - vref_2 - vref_1))

        !        print *, "Ecoh", Ecoh_2 + Ecoh_1
        !        print *, "V", - V_22 - V_11 - 0.5*(V_21 + V_12)
        !        print *, "vref",   0.5*(vref_2 + vref_1)

        !        print *, "dEcoh 11, 12, 21, 22"
        !        print *, dEcoh_1_1(:,:,1)
        !        print *,  dEcoh_1_2(:,:,1)
        !        print *,  dEcoh_2_1(:,:,1)
        !        print *,  dEcoh_2_2(:,:,1)


        !    print *, "dV_22_2, dV_21_2, dV_12_2"
        !    print *, dV_22_2(:,:,2)
        !    print *, dV_21_2(:,:,2)
        !    print *, dV_12_2(:,:,2)

        !print *, "dvref_2_2, dvref_1_2"
        !print *, dvref_2_2(:,:,2)
        !print *, dvref_1_2(:,:,2)
        !
        !
        !stop 123
        ! minus sign was taken into account in calculation of separate contributions
        do i = 1, atoms%natoms

            if (atoms%idx(i) == type2) then
                atoms%f(:,:,i) = atoms%f(:,:,i) + (dEcoh_2_2(:,:,i) + dEcoh_1_2(:,:,i) - dV_22_2(:,:,i) &
                    - 0.5*(dV_21_2(:,:,i) + dV_12_2(:,:,i) - dvref_2_2(:,:,i) - dvref_1_2(:,:,i)))
            else if (atoms%idx(i) == type1) then
                atoms%f(:,:,i) = atoms%f(:,:,i) + (dEcoh_2_1(:,:,i) + dEcoh_1_1(:,:,i) - dV_11_1(:,:,i) &
                    - 0.5*(dV_21_1(:,:,i) + dV_12_1(:,:,i) - dvref_2_1(:,:,i) - dvref_1_1(:,:,i)))
            end if

        end do

    !       print *,  "dEcoh", dEcoh_2_2(:,:,2)+ dEcoh_1_2(:,:,2)
    !       print *,  "dVf"  , -dV_22_2(:,:,2)-0.5*(dV_21_2(:,:,2) + dV_12_2(:,:,2))
    !       print *,  "dVref", 0.5*(dvref_2_2(:,:,2)+dvref_1_2(:,:,2))
    !       print *,  "atoms%f", atoms%f(:,:,2)


    !       print *,  "dEcoh", dEcoh_2_1(:,:,1)+dEcoh_1_1(:,:,1)
    !       print *,  "dVf"  , -dV_11_1(:,:,1)-0.5*(dV_21_1(:,:,1) + dV_12_1(:,:,1))
    !       print *,  "dVref", 0.5*(dvref_2_1(:,:,1)+dvref_1_1(:,:,1))
    !       print *,  "atoms%f", atoms%f(:,:,1)

    !        print *, "Epot2", atoms%epot

    end subroutine compute_emt_2species



        subroutine compute_emt_nspecies(atoms, flag)

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
            real(dp), dimension(3, atoms%nbeads, atoms%natoms) :: dVf
            real(dp), dimension(3, atoms%nbeads, atoms%ntypes, atoms%natoms) :: dVref
            real(dp), dimension(3, atoms%nbeads, atoms%ntypes, atoms%natoms, atoms%natoms) :: ds, dEcoh
            real(dp), dimension(3, atoms%nbeads, atoms%ntypes, atoms%ntypes, atoms%natoms) :: dsigmat
            real(dp), dimension(3, atoms%nbeads, atoms%ntypes, atoms%ntypes, &
                atoms%natoms, atoms%natoms) :: dsigma, dV


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


            !        print *, "chi"
            !        do i = 1, atoms%ntypes
            !        do j = 1, atoms%ntypes
            !        print  '(2i4, 2f19.15)', i,j,chi(i,j)
            !        end do
            !        end do

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

            dsigma  = 0.0_dp
            dsigmat = 0.0_dp
            dV      = 0.0_dp
            dVf     = 0.0_dp
            ds      = 0.0_dp
            dEcoh   = 0.0_dp
            dVref   = 0.0_dp

            do i = 1, atoms%natoms-1
                idx_i = atoms%idx(i)

                do j = i+1, atoms%natoms
                    idx_j = atoms%idx(j)

                    if (i == j .or. atoms%pes(idx_i,idx_j) /= pes_id_emt) cycle

                    do b = 1, atoms%nbeads

                        ! Applying PBCs
                        call minimg_one(atoms, i, j, b, b, r, vec)

                        if (r < tolerance) then
                            print *, "Error in compute_emt: distance too small between &
                            beads number"                            , b, "of atoms", i, "and", j
                            call abort
                        end if

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
                        dVf(:,:,i) = dVf(:,:,i) - 0.5_dp*dV(:,:,idx_j,idx_i,i,j)
                        dVf(:,:,j) = dVf(:,:,j) - 0.5_dp*dV(:,:,idx_j,idx_i,j,i)

                        dVf(:,:,i) = dVf(:,:,i) - 0.5_dp*dV(:,:,idx_i,idx_j,j,i)
                        dVf(:,:,j) = dVf(:,:,j) - 0.5_dp*dV(:,:,idx_i,idx_j,i,j)
                    else
                        dVf(:,:,i) = dVf(:,:,i) - dV(:,:,idx_j,idx_i,i,j)
                        dVf(:,:,j) = dVf(:,:,j) - dV(:,:,idx_j,idx_i,j,i)
                    end if
                end do
            end do

            !        print '(3f12.5)' ,dVf
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

            !        print *, "dsigmat scaled"
            !        do i = 1, atoms%natoms
            !            do k = 1, atoms%ntypes
            !                do n = 1, atoms%ntypes
            !                    print '(3i4, 3f12.5)', k,n, i, dsigmat(:,:,k,n,i)
            !                end do
            !            end do
            !        end do



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

                !            do k = 1, 3
                !                ds(k,:,idx_i,:,i) = ds(k,:,idx_i,:,i)/(betaeta2(idx_i)*s)
                !            end do

                do k = 1, 3
                    where (s > tolerance)
                        ds(k,:,idx_i,:,i) = ds(k,:,idx_i,:,i)/(betaeta2(idx_i)*s)
                    elsewhere
                        ds(k,:,idx_i,:,i) = ds(k,:,idx_i,:,i)/(betaeta2(idx_i)*tolerance)
                    end where
                end do

                do j = 1, atoms%natoms
                    idx_j = atoms%idx(j)
                    if (idx_i /= idx_j) then
                        do k = 1, 3
                            where (s(:,i) > tolerance)
                                ds(k,:,idx_i,i,j) = -chi(idx_j,idx_i) * &
                                    dsigma(k,:,idx_i,idx_j,i,j)/(betaeta2(idx_i)*s(:,i))
                            elsewhere
                                ds(k,:,idx_i,i,j) = -chi(idx_j,idx_i) * &
                                    dsigma(k,:,idx_i,idx_j,i,j)/(betaeta2(idx_i)*tolerance)
                            end where
                        end do
                    end if
                end do

            end do

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
                do b = 1, atoms%nbeads
                    s(b,i) = -log(max(s(b,i), tolerance)*twelfth)/betaeta2(idx_i)
                end do
            end do

            !              print *, "s"
            !              print '(3f15.8)',  s
            !              stop


            !----------------------EMBEDDED ELECTRON DENSITY-------------------------------

            do i = 1, atoms%natoms
                idx_i = atoms%idx(i)
                dens(:,i) = pes_emt%n0(idx_i) * exp((0.5_dp/bohr2ang - betaeta2(idx_i))*s(:,i))
            end do

            !        print '(3f15.8)', dens




            !---------------------------COHESIVE FUNCTION-----------------------------------

            do i = 1, atoms%natoms
                idx_i = atoms%idx(i)

                Ecoh = Ecoh + ((1.0_dp + pes_emt%lambda(idx_i)*s(:,i)) * &
                    exp(-pes_emt%lambda(idx_i)*s(:,i)) - 1.0_dp) * pes_emt%e0(idx_i)
            end do

            !             print *, "Ecoh"
            !             print '(2f15.8)', Ecoh
            !             stop

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

                    !        print *, "dV"
                    !        do i = 1, atoms%natoms
                    !        do j = 1, atoms%natoms
                    !        do n = 1, atoms%ntypes
                    !        do k = 1, atoms%ntypes
                    !            print '(4i4, 3f12.5)',  n,k,i, j,dV(:,:,n,k,i,j)
                    !
                    !        end do
                    !        end do
                    !        end do
                    !        end do

                    !        print *, "dVref"
                    !        do i = 1, atoms%natoms
                    !        do n = 1, atoms%ntypes
                    !        print '(2i4, 3f12.5)',  n,  i,  dVref(:,:,n,i)
                    !        end do
                    !        end do


            !-------------------------------TOTAL ENERGY---------------------------------

            !        print *, "dEcoh"
            !        atoms%f = sum(sum(dEcoh, dim=3), dim=3)
            !        print *, atoms%f(:,:,2)

            !        print *, "dVf"
            !        print *, dVf(:,:,2)

            !        print *, "dVref"
            !        atoms%f = sum(dVref, dim=3)
            !        print *, atoms%f(:,:,2)
            !        stop 234

            !    print *, "Ecoh", Ecoh
            !    print *, "V", -0.5*sum(sum(V, dim=1), dim=1)
            !    print *, "vref", 0.5*(sum(V_ref, dim=1))

            atoms%epot = atoms%epot + Ecoh - 0.5_dp * (sum(sum(V, dim=1), dim=1) - sum(V_ref, dim=1))
            !print *, atoms%epot
            atoms%f = atoms%f + sum(sum(dEcoh, dim=3), dim=3) - dVf + 0.5_dp * sum(dVref, dim=3)
        !        print *, "atoms%f"
        !        print '(3f12.5)', atoms%f

        !    print *, "dEcoh", sum(sum(dEcoh(:,:,:,:,2), dim=3), dim=3)
        !    print *, "dVf", -dVf(:,:,2)
        !    print *, "dVref", 0.5*sum(dVref(:,:,:,2), dim=3)
        !    print *, "atoms%f", atoms%f(:,:,2)

        !    print *, "dEcoh", sum(sum(dEcoh(:,:,:,:,1), dim=3), dim=3)
        !    print *, "dVf", -dVf(:,:,1)
        !    print *, "dVref", 0.5*sum(dVref(:,:,:,1), dim=3)
        !    print *, "atoms%f", atoms%f(:,:,1)

        !            print *, "Epot2", atoms%epot


         !   slab%f = dEcoh_l_l + dEcoh_p_l - dV_ll_l - 0.50d0*(dV_lp_l + dV_pl_l - dvref_l_l - dvref_p_l)

        end subroutine compute_emt_nspecies



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
                    beta = (16.0_dp*pi/3.0_dp)**(1.0_dp/3.0_dp)*isqrt2
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
