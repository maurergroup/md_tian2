
module pes_lj_mod

    use constants
    use useful_things, only : split_string
    use universe_mod

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

        integer :: nwords = 0, ios = 0, i
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
            pes_lj%cutoff = minval(temp3(:))
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


    subroutine compute_lj(atoms, atom_i, atom_j, flag)

        type(universe), intent(inout) :: atoms
        integer, intent(in)           :: atom_i, atom_j
        integer, intent(in)           :: flag

        integer ::  idx_i, idx_j, b
        real(dp) :: sig_rc, sig_rc_2, sig_rc_6, sig_rc_12
        real(dp), dimension(atoms%nbeads) :: sig_r, sig_r_2, sig_r_6, sig_r_12, r
        real(dp), dimension(atoms%nbeads) :: nrg, vdr
        real(dp), dimension(3, atoms%nbeads) :: f, vec
        !        real(dp) :: tmp

        nrg = 0.0_dp
        vdr = 0.0_dp

        idx_i = atoms%idx(atom_i)
        idx_j = atoms%idx(atom_j)

        call minimg_beads(atoms, atom_i, atom_j, r, vec)

        sig_rc = pes_lj%sigma(idx_i,idx_j) / pes_lj%cutoff
        sig_rc_2 = sig_rc * sig_rc
        sig_rc_6 = sig_rc_2 * sig_rc_2 * sig_rc_2
        sig_rc_12 = sig_rc_6 * sig_rc_6

        where (r < tolerance) r = tolerance


        sig_r = pes_lj%sigma(idx_i,idx_j) / r
        sig_r_2 = sig_r * sig_r
        sig_r_6 = sig_r_2 * sig_r_2 * sig_r_2
        sig_r_12 = sig_r_6 * sig_r_6

        nrg = 4*pes_lj%eps(idx_i,idx_j) * ( (sig_r_12-sig_r_6)  &
            + (6*sig_rc_12-3*sig_rc_6)*(r/pes_lj%cutoff)**2 &
            - 7*sig_rc_12 + 4*sig_rc_6 )

        atoms%epot = atoms%epot + nrg


        if (flag == energy_and_force) then



            vdr = (24/pes_lj%cutoff**2)*r*pes_lj%eps(idx_i,idx_j)*sig_rc_6*(2*sig_rc_6-1) &
                - (24/r)*pes_lj%eps(idx_i,idx_j)*sig_r_6*(2*sig_r_6-1)
            !
            !
            !                        r2 = r*r
            !                        fr2 = pes_lj%sigma(idx_i,idx_j) / r2
            !                        fr6 = fr2 * fr2 * fr2
            !                        fpr = 48._dp * pes_lj%eps(idx_i,idx_j) * fr6 * (fr6 - 0.5_dp) / r2   ! f/r

            forall(b = 1 : atoms%nbeads) f(:,b) = vdr(b) * vec(:,b)/r(b)

            atoms%f(:,:,atom_i) = atoms%f(:,:,atom_i) + f
            atoms%f(:,:,atom_j) = atoms%f(:,:,atom_j) - f

        end if


    end subroutine compute_lj



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
