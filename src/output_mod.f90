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
    integer :: out_id_poscar = 0
    integer :: out_id_mxt    = 0


contains

    subroutine output(atoms, itraj, istep, flag)

        type(universe), intent(in) :: atoms
        integer, intent(in)        :: itraj, istep
        character(len=*), intent(in), optional :: flag

        integer :: i
        character(len=*), parameter :: err = "Error in output(): "

        ! use only for start of md trajectory
        if (present(flag)) then
            call output_scatter(atoms, itraj, istep, flag)
        end if

        do i = 1, size(simparams%output_type)
            if (modulo(istep, simparams%output_interval(i)) == 0) then

                select case (simparams%output_type(i))

                    case (output_id_xyz)
                        call output_xyz(atoms, itraj)

                    case (output_id_energy)
                        call output_nrg(atoms, itraj, istep)

                    case (output_id_poscar)
                        call output_poscar(atoms)
                        out_id_poscar = out_id_poscar + 1

                    case (output_id_mxt)
                        call output_mxt(atoms)
                        out_id_mxt = out_id_mxt + 1

                    case (output_id_scatter)
                        ! pass

                    case default
                        print *, err // "unknown output format", simparams%output_type(i)
                        stop
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
        dir_coords = dir_coords - anint(dir_coords-0.5)
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
        integer  :: i
        real(dp) :: atom_temp, bead_temp, rgyr_p, rgyr_l , a_ekin_p, a_ekin_l, &
            q_ekin_l, q_ekin_p, b_ekin_p, b_ekin_l, atom_epot, bead_epot, etotal


        ! XXX: change system() to execute_command_line() when new compiler is available
        if (.not. dir_exists('conf')) call system('mkdir conf')

        write(traj_id,'(i8.8)') itraj
        fname = 'conf/mxt_trj'//traj_id//'.dat'

        if (overwrite_nrg) then
            call open_for_write(out_unit, fname)
            write(out_unit, '(a)') '# All energies are given in eV.'
            write(out_unit, '(a1, i8, 14i17)') "#", (i, i = 1, 15)
            write(out_unit, '(15a17)') '#time/fs', 'atom T/K', 'bead T/K', 'r_gyr_sq_p/A²', &
                'r_gyr_sq_l/A²', 'atom ekin proj', 'atom ekin latt', 'bead ekin proj', &
                'bead ekin latt', 'qntm ekin proj', 'qntm ekin latt', 'atom epot', &
                'bead epot', 'e_total', 'f_total/(eV/A)'
            overwrite_nrg = .false.
        else
            call open_for_append(out_unit,fname)
        end if

        atom_temp = calc_atom_temperature(atoms)
        bead_temp = calc_bead_temperature(atoms)

        call radius_of_gyration(atoms, rgyr_p, rgyr_l)
        call atom_ekin(atoms, a_ekin_p, a_ekin_l)
        call bead_ekin(atoms, b_ekin_p, b_ekin_l)
        call virial_quantum_ekin(atoms, q_ekin_p, q_ekin_l)

        atom_epot = sum(atoms%epot)/atoms%nbeads
        bead_epot = calc_bead_epot(atoms)

        etotal = b_ekin_p + b_ekin_l + atom_epot + bead_epot

        ! primitive quantum kinetic energy estimator, do not use, just for completeness
        ! source: J. Chem. Phys. 123, 134502 (2005)
        ! q_ekin = (0.5 * atoms%dof * kb * bead_temp - bead_epot)/atoms%nbeads

        write(out_unit, '(15e17.8e2)') istep*simparams%step, atom_temp, bead_temp, &
            rgyr_p, rgyr_l, a_ekin_p, a_ekin_l, b_ekin_p, b_ekin_l, q_ekin_p, q_ekin_l, &
            atom_epot, bead_epot, etotal, sum(atoms%f)

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

        ! open file conf/poscar_%08d.dat
        write(fid,'(i8.8)') out_id_poscar+simparams%start
        fname = 'conf/poscar_'//fid//'.dat'
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




    subroutine output_mxt(atoms)

        type(universe), intent(in) :: atoms

        character(len=max_string_length) :: fname
        character(len=8)                 :: fid

        ! XXX: change system() to execute_command_line() when new compiler is available
        if (.not. dir_exists('conf')) call system('mkdir conf')

        ! open file conf/mxt_%08d.dat
        write(fid,'(i8.8)') out_id_mxt+simparams%start
        fname = 'conf/mxt_'//fid//'.dat'
        open(out_unit, file=fname, form="unformatted", status="replace")

        write(out_unit) atoms%natoms, atoms%nbeads, atoms%ntypes
        write(out_unit) atoms%r, atoms%v, atoms%is_cart, atoms%is_fixed, &
            atoms%idx, atoms%name, atoms%is_proj, atoms%simbox, atoms%isimbox

        close(out_unit)

    end subroutine output_mxt




    subroutine output_scatter(atoms, itraj, istep, flag)

        type(universe), intent(in) :: atoms
        integer, intent(in) :: itraj, istep
        character(len=*), intent(in) :: flag


        real(dp) :: ekin_p, ekin_l, atom_epot, proj_r(3), proj_v(3), proj_polar, proj_azi, time
        integer  :: turning_points
        character(len=max_string_length) :: fname
        character(len=8)                 :: fid
        character(len=*), parameter :: err = "Error in output_scatter: "

        if (count(atoms%is_proj) > 1) stop err // "output_scatter does not work for more that one projectile, please implement"

        ! XXX: change system() to execute_command_line() when new compiler is available
        if (.not. dir_exists('traj')) call system('mkdir traj')

        write(fid,'(i8.8)') simparams%start
        fname = 'traj/mxt_fin'//fid//'.dat'

        if (flag == "scatter_initial") then

            call atom_ekin(atoms, ekin_p, ekin_l)
            atom_epot = sum(atoms%epot)/atoms%nbeads
            proj_r = sum(atoms%r(:,:,1), dim=2)/atoms%nbeads
            proj_v = sum(atoms%v(:,:,1), dim=2)/atoms%nbeads

            proj_polar = 180-acos(proj_v(3)/sqrt(sum(proj_v*proj_v)))*rad2deg
            proj_azi = atan2(proj_v(2), proj_v(1))*rad2deg

            if (file_exists(fname)) then
                open(out_unit, file=fname, status="old", position="append", action="write")
            else
                open(out_unit, file=fname, status="new", action="write")
                write(out_unit, '(19a14)') "# traj_index", "ekin_p_i/eV", "ekin_l_i/eV", "epot/eV", &
                    "p_x_i/A", "p_y_i/A", "p_z_i/A", "polar_i/deg", "azi_i/deg", "ekin_p_f/eV", &
                    "ekin_l_f/eV", "epot/eV", "p_x_f/A", "p_y_f/A", "p_z_f/A", "polar_f/deg", &
                    "azi_f/deg", "time/fs", "turn_pnts"
            end if

            write(out_unit, '(i14, 8f14.7$)', advance="no") itraj, ekin_p, ekin_l, atom_epot, proj_r, proj_polar, proj_azi

            close (out_unit)

        else if (flag == "scatter_final") then

            call atom_ekin(atoms, ekin_p, ekin_l)
            atom_epot = sum(atoms%epot)/atoms%nbeads
            proj_r = sum(atoms%r(:,:,1), dim=2)/atoms%nbeads
            proj_v = sum(atoms%v(:,:,1), dim=2)/atoms%nbeads

            proj_polar = acos(proj_v(3)/sqrt(sum(proj_v*proj_v)))*rad2deg
            proj_azi = atan2(proj_v(2), proj_v(1))*rad2deg

            time = (istep-1) * simparams%step

            turning_points = calc_turning_points()

            call open_for_append(out_unit, fname)
            write(out_unit, '(9f14.7, i14)') ekin_p, ekin_l, atom_epot, proj_r, proj_polar, proj_azi, time, turning_points

            close (out_unit)

        else
            print *, err, "unknown flag", flag
            stop
        end if

    end subroutine output_scatter



    subroutine output_pes(atoms)

        type(universe), intent(in) :: atoms

        integer :: i, j
        character(len=28) :: fname
        character(len=5) :: pes
        character(len=4) :: role_i, role_j

        write(fname, '(a19, i5.5, a4)') "fit/out_params/fit_", simparams%start, ".pes"

        ! XXX: change system() to execute_command_line() when new compiler is available
        if (.not. dir_exists('fit/out_params/')) call system('mkdir fit/out_params/')

        call open_for_write(out_unit, fname)

        do i = 1, atoms%ntypes
            do j = 1, atoms%ntypes
                pes = pes_id_to_name(atoms%pes(i,j))
                role_i = universe_id_to_role(atoms, i)
                role_j = universe_id_to_role(atoms, j)

                write(out_unit, '(2a5)'), "pes", pes
                write(out_unit, '(4a8)'), atoms%name(i), atoms%name(j), role_i, role_j

            end do
        end do


    end subroutine output_pes

end module output_mod
