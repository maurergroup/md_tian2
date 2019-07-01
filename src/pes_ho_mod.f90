!############################################################################
! This routine is part of
! md_tian2 (Molecular Dynamics Xia Tian 2)
! (c) 2014-2019 Dan J. Auerbach, Sascha Kandratsenka, Svenja M. Janke, Marvin
! Kammler, Sebastian Wille
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

module pes_ho_mod

    use constants
    use useful_things, only : split_string, lower_case
    use universe_mod

    implicit none


    type ho_pes

        private
        real(dp), allocatable :: r0(:,:)
        real(dp), allocatable :: k(:,:)

    end type ho_pes

    type(ho_pes) :: pes_ho

contains

    subroutine read_ho(atoms, inp_unit)

        use run_config, only : simparams

        type(universe), intent(inout) :: atoms
        integer, intent(in) :: inp_unit

        integer :: nwords, ios = 0, i
        character(len=max_string_length) :: buffer
        character(len=max_string_length) :: words(100)
        integer  :: idx1, idx2, ntypes
        real(dp) :: dummy, reduced_mass
        character(len=*), parameter :: err = "Error in read_ho: "

        ntypes = simparams%nprojectiles+simparams%nlattices

        if (.not. allocated(pes_ho%r0)) then
            allocate(pes_ho%r0(ntypes, ntypes), pes_ho%k(ntypes, ntypes))
            pes_ho%r0 = default_real
            pes_ho%k  = default_real
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
        atoms%pes(idx1,idx2) = pes_id_ho
        atoms%pes(idx2,idx1) = pes_id_ho

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

                case ('r0')      ! given in angstrom
                    read(words(2), *) pes_ho%r0(idx1, idx2)
                    read(words(2), *) pes_ho%r0(idx2, idx1)

                ! possibilities to define the force constant
                case ('t')       ! period, given in fs
                    read(words(2), *) dummy
                    reduced_mass = atoms%m(idx1)*atoms%m(idx2) / (atoms%m(idx1)+atoms%m(idx2))
                    pes_ho%k(idx1, idx2) = (2.0_dp * pi / dummy)**2 * reduced_mass
                    pes_ho%k(idx2, idx1) = (2.0_dp * pi / dummy)**2 * reduced_mass

                case ('nu')      ! frequency, given in 1/fs
                    read(words(2), *) dummy
                    reduced_mass = atoms%m(idx1)*atoms%m(idx2) / (atoms%m(idx1)+atoms%m(idx2))
                    pes_ho%k(idx1, idx2) = (2.0_dp * pi * dummy)**2 * reduced_mass
                    pes_ho%k(idx2, idx1) = (2.0_dp * pi * dummy)**2 * reduced_mass

                case ('k')       ! force constant, given in eV/A²
                    read(words(2), *) pes_ho%k(idx1, idx2)
                    read(words(2), *) pes_ho%k(idx2, idx1)

                case default
                    print *, "Error in the PES file: unknown HO parameter", words(1)
                    stop

            end select
        end do

    end subroutine read_ho



    subroutine compute_ho(atoms, flag)

        type(universe), intent(inout) :: atoms
        integer, intent(in)           :: flag

        integer ::  idx_i, idx_j, i, j, b

        real(dp), dimension(3, atoms%nbeads) :: f, vec
        real(dp), dimension(atoms%nbeads) :: r

        do i = 1, atoms%natoms-1
            do j = i+1, atoms%natoms

                idx_i = atoms%idx(i)
                idx_j = atoms%idx(j)

                if (atoms%pes(idx_i,idx_j) /= pes_id_ho) cycle

                call minimg_beads(atoms, i, j, r, vec)

                if (any(r > 1.5*pes_ho%r0(idx_i,idx_j))) cycle

                if (any(r < tolerance)) then
                    print *, "Error in compute_ho: distance too small between &
                        beads number", minloc(r), "of atoms", i, "and", j
                    call abort
                end if

                atoms%Epot = atoms%Epot + 0.5_dp * pes_ho%k(idx_i,idx_j)*(r-pes_ho%r0(idx_i,idx_j))**2

                if (flag == energy_and_force) then
                    do b = 1, atoms%nbeads
                        f(:,b) = pes_ho%k(idx_i,idx_j) * (r(b) - pes_ho%r0(idx_i,idx_j)) * vec(:,b)/r(b)
                    end do

                    atoms%f(:,:,i) = atoms%f(:,:,i) + f
                    atoms%f(:,:,j) = atoms%f(:,:,j) - f
                end if


            end do
        end do

    end subroutine compute_ho

end module pes_ho_mod
