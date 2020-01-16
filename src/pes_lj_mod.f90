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

module pes_lj_mod

    use constants
    use universe_mod
    use useful_things, only : split_string, lower_case

    implicit none


    type lj_pes

        private
        real(dp), allocatable :: eps(:,:)
        real(dp), allocatable :: sigma(:,:)
        real(dp), allocatable :: shift(:,:)
        real(dp)              :: cutoff

    end type lj_pes

    type(lj_pes) :: pes_lj

contains

    subroutine read_lj(atoms, inp_unit)

        use run_config, only : simparams

        type(universe), intent(inout) :: atoms
        integer, intent(in) :: inp_unit

        integer :: nwords, ios = 0, i
        character(len=max_string_length) :: buffer
        character(len=max_string_length) :: words(100)
        integer  :: idx1, idx2, ntypes
        real(dp) :: temp3(3)
        character(len=*), parameter :: err = "Error in read_lj: "

        ntypes = simparams%nprojectiles+simparams%nlattices

        if (.not. allocated(pes_lj%eps)) then
            allocate(pes_lj%eps(ntypes, ntypes),   &
                pes_lj%sigma(ntypes, ntypes), &
                pes_lj%shift(ntypes, ntypes))
            pes_lj%eps   = default_real
            pes_lj%sigma = default_real
            pes_lj%shift = default_real

            forall (i = 1 : 3) temp3(i) = sqrt(sum(atoms%simbox(:,i)*atoms%simbox(:,i)))
            pes_lj%cutoff = minval(temp3(:))/2
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
        atoms%pes(idx1,idx2) = pes_id_lj
        atoms%pes(idx2,idx1) = pes_id_lj

        do
            read(inp_unit, '(A)', iostat=ios) buffer
            call split_string(buffer, words, nwords)

            ! pes block terminated, set shift
            if (nwords == 0 .or. ios /= 0) then

                pes_lj%shift(idx1,idx2) = 4.0_dp * pes_lj%eps(idx1,idx2) * &
                    ( (pes_lj%sigma(idx1,idx2)/pes_lj%cutoff)**12 &
                    - (pes_lj%sigma(idx1,idx2)/pes_lj%cutoff)**6 )

                pes_lj%shift(idx2,idx1)  = pes_lj%shift(idx1,idx2)

                exit

            ! something went wrong
            else if (nwords /= 2) then
                stop "Error in the PES file: PES parameters must consist of key value pairs. A parameter block must be terminated by a blank line."
            end if

            call lower_case(words(1))

            select case (words(1))

                case ('sigma')      ! given in angstrom
                    read(words(2), *) pes_lj%sigma(idx1, idx2)
                    read(words(2), *) pes_lj%sigma(idx2, idx1)

                case ('epsilon')    ! given in kelvin
                    read(words(2), *) pes_lj%eps(idx1, idx2)
                    read(words(2), *) pes_lj%eps(idx2, idx1)
                    pes_lj%eps(idx1, idx2) = pes_lj%eps(idx1, idx2) * kelvin2ev
                    pes_lj%eps(idx2, idx1) = pes_lj%eps(idx2, idx1) * kelvin2ev

                case default
                    print *, "Error in the PES file: unknown LJ parameter", words(1)
                    stop

            end select
        end do

    end subroutine read_lj



    subroutine read_simple_lj(atoms, inp_unit)

        use run_config, only : simparams

        type(universe), intent(inout) :: atoms
        integer, intent(in) :: inp_unit

        integer :: nwords, ios = 0
        character(len=max_string_length) :: buffer
        character(len=max_string_length) :: words(100)
        integer  :: idx1, idx2, ntypes
        character(len=*), parameter :: err = "Error in read_simple_lj: "

        ntypes = simparams%nprojectiles+simparams%nlattices

        if (.not. allocated(pes_lj%eps)) then
            allocate(pes_lj%eps(ntypes, ntypes), pes_lj%sigma(ntypes, ntypes))
            pes_lj%eps   = default_real
            pes_lj%sigma = default_real
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
        atoms%pes(idx1,idx2) = pes_id_simple_lj
        atoms%pes(idx2,idx1) = pes_id_simple_lj

        do
            read(inp_unit, '(A)', iostat=ios) buffer
            call split_string(buffer, words, nwords)

            ! pes block terminated, set shift
            if (nwords == 0 .or. ios /= 0) then
                exit

            ! something went wrong
            else if (nwords /= 2) then
                stop "Error in the PES file: PES parameters must consist of key value pairs. A parameter block must be terminated by a blank line."
            end if

            select case (words(1))

                case ('sigma')      ! given in angstrom
                    read(words(2), *) pes_lj%sigma(idx1, idx2)
                    read(words(2), *) pes_lj%sigma(idx2, idx1)

                case ('epsilon')    ! given in kelvin
                    read(words(2), *) pes_lj%eps(idx1, idx2)
                    read(words(2), *) pes_lj%eps(idx2, idx1)
                    pes_lj%eps(idx1, idx2) = pes_lj%eps(idx1, idx2) * kelvin2ev
                    pes_lj%eps(idx2, idx1) = pes_lj%eps(idx2, idx1) * kelvin2ev

                case default
                    print *, "Error in the PES file: unknown LJ parameter", words(1)
                    stop

            end select
        end do
    end subroutine


    subroutine compute_lj(atoms, flag)

        type(universe), intent(inout) :: atoms
        integer, intent(in)           :: flag

        integer ::  idx_i, idx_j, i, j, b
        real(dp) :: sig_rc, sig_rc_2, sig_rc_6, sig_rc_12
        real(dp), dimension(atoms%nbeads) :: sig_r, sig_r_2, sig_r_6, sig_r_12, r
        real(dp), dimension(atoms%nbeads) :: nrg, vdr
        real(dp), dimension(3, atoms%nbeads) :: f, vec

        nrg = 0.0_dp

        do i = 1, atoms%natoms-1
            do j = i+1, atoms%natoms

                idx_i = atoms%idx(i)
                idx_j = atoms%idx(j)

                if (atoms%pes(idx_i,idx_j) /= pes_id_lj) cycle

                call minimg_beads(atoms, i, j, r, vec)

                if (any(r < tolerance)) then
                    print *, "Error in compute_lj: distance too small between &
                        beads number", minloc(r), "of atoms", i, "and", j
                    call abort
                end if

                sig_rc = pes_lj%sigma(idx_i,idx_j) / pes_lj%cutoff
                sig_rc_2  = sig_rc   * sig_rc
                sig_rc_6  = sig_rc_2 * sig_rc_2 * sig_rc_2
                sig_rc_12 = sig_rc_6 * sig_rc_6

                sig_r = pes_lj%sigma(idx_i,idx_j) / r
                sig_r_2  = sig_r   * sig_r
                sig_r_6  = sig_r_2 * sig_r_2 * sig_r_2
                sig_r_12 = sig_r_6 * sig_r_6

                nrg = 4*pes_lj%eps(idx_i,idx_j) * ( (sig_r_12-sig_r_6)  &
                    + (6*sig_rc_12-3*sig_rc_6)*(r/pes_lj%cutoff)**2 &
                    - 7*sig_rc_12 + 4*sig_rc_6 )

                atoms%epot = atoms%epot + nrg


                if (flag == energy_and_force) then

                    vdr = (24/pes_lj%cutoff**2)*r*pes_lj%eps(idx_i,idx_j)*sig_rc_6*(2*sig_rc_6-1) &
                        - (24/r)*pes_lj%eps(idx_i,idx_j)*sig_r_6*(2*sig_r_6-1)

                    do b = 1, atoms%nbeads
                        f(:,b) = vdr(b) * vec(:,b)/r(b)
                    end do

                    atoms%f(:,:,i) = atoms%f(:,:,i) + f
                    atoms%f(:,:,j) = atoms%f(:,:,j) - f

                end if

            end do
        end do


    end subroutine compute_lj



    subroutine compute_simple_lj(atoms, flag)

        type(universe), intent(inout) :: atoms
        integer, intent(in)           :: flag

        integer ::  idx_i, idx_j, i, j, b

        real(dp), dimension(3, atoms%nbeads) :: f, vec
        real(dp), dimension(atoms%nbeads) :: r2, fr2, fr6, fpr, nrg, r

        nrg = 0.0_dp

        do i = 1, atoms%natoms-1
            do j = i+1, atoms%natoms

                idx_i = atoms%idx(i)
                idx_j = atoms%idx(j)

                if (atoms%pes(idx_i,idx_j) /= pes_id_simple_lj) cycle

                call minimg_beads(atoms, i, j, r, vec)

                if (any(r < tolerance)) then
                    print *, "Error in compute_simple_lj: distance too small between &
                        beads number", minloc(r), "of atoms", i, "and", j
                    call abort
                end if

                r2  = sum(vec*vec, dim=1)

                fr2 = pes_lj%sigma(idx_i, idx_j)**2 / r2
                fr6 = fr2 * fr2 * fr2

                nrg =  4 * pes_lj%eps(idx_i, idx_j) * fr6 * (fr6 - 1)
                atoms%Epot = atoms%Epot + nrg

                if (flag == energy_and_force) then
                    fpr = 24 * pes_lj%eps(idx_i, idx_j) * fr6 * (1 - 2*fr6) / r2

                    do b = 1, atoms%nbeads
                        f(:,b) = fpr(b) * vec(:,b)
                    end do

                    atoms%f(:,:,i) = atoms%f(:,:,i) + f
                    atoms%f(:,:,j) = atoms%f(:,:,j) - f
                end if
            end do
        end do
    end subroutine compute_simple_lj



end module pes_lj_mod




!pes lj
!H   H   proj    proj
!sigma   2.4
!epsilon 3.0
!
!pes lj
!H   D   proj    proj
!sigma   2.4
!epsilon 3.0
!sigma   2.4
!epsilon 3.0
!
!pes lj
!D   D   proj    proj
!sigma   2.4
!epsilon 3.0
