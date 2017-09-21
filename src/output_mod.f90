module output_mod

    use universe_mod
    use useful_things
    use open_file
    use constants
    use run_config, only : simparams

    implicit none

    logical :: overwrite_nrg = .true.
    logical :: overwrite_xyz = .true.
    integer, parameter :: out_unit = 86
    integer :: out_id = 1


contains

    subroutine output(atoms, itraj, istep)

        type(universe), intent(in) :: atoms
        integer, intent(in)        :: itraj, istep

        integer :: i
        character(len=*), parameter :: err = "Error in output(): "

        do i = 1, size(simparams%output_type)
            if (modulo(istep, simparams%output_interval(i)) == 0) then
                select case (simparams%output_type(i))
                case (output_id_xyz)
                    call output_xyz(atoms, itraj)

                case (output_id_energy)
                    call output_nrg(atoms, itraj, istep)

                case (output_id_poscar)
                    call output_poscar(atoms)
                    out_id = out_id + 1

                case default
                    stop err // "unknown output format"
                end select
            end if
        end do

    end subroutine output



    subroutine output_xyz(atoms, itraj)

        type(universe), intent(in) :: atoms
        integer, intent(in) :: itraj

        character(len=*), parameter :: xyz_format = '(a, 3e18.8e3)'

        character(len=8) :: traj_id
        character(len=max_string_length) :: fname
        integer :: i, j
        real(dp), dimension(3, atoms%nbeads, atoms%natoms) :: dir_coords, cart_coords


        ! XXX: change system() to execute_command_line() when new compiler is available
        if (.not. dir_exists('conf')) call system('mkdir conf')

        ! to direct, fold into simbox, to cartesian

        forall (i = 1 : atoms%natoms) dir_coords(:,:,i) = matmul(atoms%isimbox, atoms%r(:,:,i))

        dir_coords = dir_coords - anint(dir_coords-0.5_dp)

        forall (i = 1 : atoms%natoms) cart_coords(:,:,i) = matmul(atoms%simbox, dir_coords(:,:,i))

        write(traj_id,'(i8.8)') itraj
        fname = 'conf/mxt_conf'//traj_id//'.xyz'

        if (overwrite_xyz) then
            call open_for_write(out_unit, fname)
            overwrite_xyz = .false.
        else
            call open_for_append(out_unit,fname)
        end if

        write(out_unit,*) atoms%natoms*atoms%nbeads
        write(out_unit,*)

        do i = 1, atoms%natoms
            do j = 1, atoms%nbeads
                write(out_unit, xyz_format) atoms%name(atoms%idx(i)), cart_coords(:,j,i)
            end do
        end do

        close(out_unit)

    end subroutine output_xyz



    subroutine output_nrg(atoms, itraj, istep)

        ! time, epot, ekin_p, ekin_l, etotal,
        use rpmd

        type(universe), intent(in) :: atoms
        integer, intent(in)        :: itraj, istep

        character(len=8)                 :: traj_id
        character(len=max_string_length) :: fname
        real(dp) :: atom_temp, bead_temp, gyr_rad_sq, a_ekin_p, a_ekin_l
        real(dp) :: b_ekin_p, b_ekin_l, atom_epot, bead_epot, etotal, b2b_dist


        ! XXX: change system() to execute_command_line() when new compiler is available
        if (.not. dir_exists('conf')) call system('mkdir conf')

        write(traj_id,'(i8.8)') itraj
        fname = 'conf/mxt_trj'//traj_id//'.dat'

        if (overwrite_nrg) then
            call open_for_write(out_unit, fname)
            write(out_unit, '(11a17)') '#time/fs', 'atom T/K', 'bead T/K', 'gyr_squared/A²', &
                'atom ekin/eV', 'bead ekin/eV', 'atom epot/eV', 'bead epot/eV','e_total/eV', 'f_total', 'b2b dist/A'
            overwrite_nrg = .false.
        else
            call open_for_append(out_unit,fname)
        end if

        atom_temp = calc_atom_temperature(atoms)
        bead_temp = calc_bead_temperature(atoms)

        gyr_rad_sq = calc_radius_of_gyration(atoms)

        call atom_ekin(atoms, a_ekin_p, a_ekin_l)
        call bead_ekin(atoms, b_ekin_p, b_ekin_l)

        atom_epot = sum(atoms%epot)/atoms%nbeads
        bead_epot = calc_bead_epot(atoms)
        b2b_dist = sum(calc_inter_bead_distances(atoms))/atoms%natoms/atoms%nbeads

        etotal = 0

        write(out_unit, '(11e17.8e2)') istep*simparams%step, atom_temp, bead_temp, &
            gyr_rad_sq, a_ekin_l, b_ekin_l, atom_epot, bead_epot, etotal, sum(atoms%f), b2b_dist

        close(out_unit)

    end subroutine output_nrg


    subroutine output_poscar(atoms)

        type(universe), intent(in) :: atoms

        character(len=max_string_length) :: fname
        character(len=8)                 :: fid
        integer :: time_vals(8), noccurrences(atoms%ntypes), i, j

        ! XXX: change system() to execute_command_line() when new compiler is available
        if (.not. dir_exists('conf')) call system('mkdir conf')

        ! prepare arrays
        noccurrences = 0
        do i = 1, atoms%natoms
            noccurrences(atoms%idx(i)) = noccurrences(atoms%idx(i)) + 1
        end do

        ! open file conf/mxt_%08d.POSCAR
        write(fid,'(i8.8)') out_id
        fname = 'conf/mxt_'//fid//'.POSCAR'
        call open_for_write(out_unit, fname)

        ! write date and time as comment line
        call date_and_time(values=time_vals)
        write(out_unit, '(i4, a, i2.2, a, i2.2, a, i2.2, a, i2.2)') &
            time_vals(1), "-", time_vals(2), "-",time_vals(3), " - ",time_vals(5), ".",time_vals(6)

        write(out_unit, '(a)')       '1.0'                    ! scaling constant
        write(out_unit, '(3f23.15)') atoms%simbox
        write(out_unit, *)           atoms%name
        write(out_unit, '(a2, i0, a1, 100i6)') ":", atoms%nbeads, ":" , (noccurrences(i), i=1,atoms%ntypes)

        if (atoms%is_cart) then
            write(out_unit, '(a)') "Cartesian"
        else
            write(out_unit, '(a)') "Direct"
        end if

        ! positions and velocities
        write(out_unit, '(3f23.15, 3l)') ((atoms%r(:,j,i), &
            .not.atoms%is_fixed(:,j,i), j=1,atoms%nbeads), i=1,atoms%natoms)
        write(out_unit, '(a)') ""
        write(out_unit, '(3f23.15)') atoms%v


        close(out_unit)

    end subroutine output_poscar

end module output_mod
