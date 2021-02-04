!############################################################################
! This routine is part of
! md_tian2 (Molecular Dynamics Tian Xia 2)
! (c) 2014-2020 Daniel J. Auerbach, Svenja M. Janke, Marvin Kammler,
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

module output_mod

    use universe_mod
    use useful_things
    use open_file
    use constants
    use trajectory_info
    use run_config, only : simparams
    use rpmd, only : calc_centroid_positions

    implicit none

    logical :: overwrite_nrg       = .true.
    logical :: overwrite_ads       = .true.
    logical :: overwrite_runner    = .true.
    integer, parameter :: out_unit = 86
    integer, parameter :: fin_unit = 87
    integer :: out_id_poscar       = 0
    integer :: out_id_mxt          = 0

    ! nene pes
    integer :: total_count_extrapolation_warnings = 0


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

        if (simparams%loutput) then ! this should help to reduce if conditions per step

            do i = 1, size(simparams%output_type)
                if (modulo(istep, simparams%output_interval(i)) == 0) then

                    select case (simparams%output_type(i))

                        case (output_id_xyz)
                            call output_xyz(atoms, itraj, istep)

                        case (output_id_energy)
                            call output_nrg(atoms, itraj, istep)

                        case (output_id_poscar)
                            call output_poscar(atoms)
                            out_id_poscar = out_id_poscar + 1

                        case (output_id_vasp)
                            call output_vasp(atoms, itraj, istep)

                        case (output_id_mxt)
                            call output_mxt(atoms)
                            out_id_mxt = out_id_mxt + 1

                        case (output_id_scatter)
                            ! pass

                        case (output_id_is_adsorbed)
                            call output_adsorption_status(atoms, itraj, istep)

                        case (output_id_nene)
                            print *, "MD step ", istep

                        case (output_id_aims)
                            call output_aims(atoms, itraj, istep)

                        case (output_id_runner)
                            call output_runner(atoms, itraj, istep)

                        case (output_id_beads)
                            ! pass

                        case default
                            print *, err // "unknown output format", simparams%output_type(i)
                            stop
                    end select

                end if ! interval

            end do ! requested output format

        end if ! loutput

    end subroutine output


    subroutine output_xyz(atoms, itraj, istep)

        type(universe), intent(in) :: atoms
        integer, intent(in) :: itraj, istep

        character(len=*), parameter :: xyz_format = '(a, 3e18.8e3)'

        character(len=8) :: traj_id
        character(len=max_string_length) :: fname
        integer :: i, j
        real(dp), dimension(3, atoms%nbeads, atoms%natoms) :: dir_coords, cart_coords

        if (.not. dir_exists('xyz')) call execute_command_line('mkdir xyz')

        ! to direct, fold into simbox, to cartesian
        forall (i = 1 : atoms%natoms) dir_coords(:,:,i) = matmul(atoms%isimbox, atoms%r(:,:,i))
        dir_coords = dir_coords - anint(dir_coords-0.5)
        forall (i = 1 : atoms%natoms) cart_coords(:,:,i) = matmul(atoms%simbox, dir_coords(:,:,i))
        where (cart_coords(3,:,:) > 0.5*atoms%simbox(3,3)) &
            cart_coords(3,:,:) = cart_coords(3,:,:) - atoms%simbox(3,3)

        write(traj_id,'(i8.8)') itraj
        fname = 'xyz/traj_'//traj_id//'.xyz'

        ! if this is the first call to this subroutine during a trajectory, find ID index
        ! super ugly, but the findloc intrinsic is not implemented in ifort 2013
        do i=1, size(simparams%output_type)
            if (simparams%output_type(i) == output_id_xyz) j = i
        end do

        if (simparams%output_interval(j) == istep) then
            call open_for_write(out_unit, fname)
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

        use rpmd

        type(universe), intent(in) :: atoms
        integer, intent(in)        :: itraj, istep

        character(len=8)                 :: traj_id
        character(len=max_string_length) :: fname
        integer  :: i
        real(dp) :: atom_temp, bead_temp, rgyr_p, rgyr_l , a_ekin_p, a_ekin_l, &
            q_ekin_l, q_ekin_p, b_ekin_p, b_ekin_l, atom_epot, bead_epot, etotal


        if (.not. dir_exists('energy')) call execute_command_line('mkdir energy')

        write(traj_id,'(i8.8)') itraj
        fname = 'energy/traj_'//traj_id//'.dat'

        if (overwrite_nrg) then
            call open_for_write(out_unit, fname)
            write(out_unit, '(a)') '# All energies are given in eV.'
            write(out_unit, '(a1, i8, 20i17)') "#", (i, i = 1, 21)
            write(out_unit, '(21a17)') '#time/fs', 'atom T/K', 'bead T/K', 'r_gyr_sq_p/A²', &
                'r_gyr_sq_l/A²', 'atom ekin proj', 'atom ekin latt', 'bead ekin proj', &
                'bead ekin latt', 'qntm ekin proj', 'qntm ekin latt', 'atom epot', &
                'bead epot', 'e_total', 'f_total/(eV/A)', 'r_x', 'r_y', 'r_z', 'v_x', 'v_y', 'v_z'
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



        write(out_unit, '(21e17.8e2)') istep*simparams%step, atom_temp, bead_temp, &
            rgyr_p, rgyr_l, a_ekin_p, a_ekin_l, b_ekin_p, b_ekin_l, q_ekin_p, q_ekin_l, &
            atom_epot, bead_epot, etotal, sum(atoms%f), sum(atoms%r(:,:,1), dim=2)/atoms%nbeads, &
            sum(atoms%v(:,:,1), dim=2)/atoms%nbeads

        close(out_unit)

    end subroutine output_nrg


    subroutine output_poscar(atoms)

        type(universe), intent(in) :: atoms

        character(len=max_string_length) :: fname
        character(len=8)                 :: fid
        integer :: time_vals(8), noccurrences(atoms%ntypes), i, j

        if (.not. dir_exists('poscar')) call execute_command_line('mkdir poscar')

        ! prepare arrays
        noccurrences = 0
        do i = 1, atoms%natoms
            noccurrences(atoms%idx(i)) = noccurrences(atoms%idx(i)) + 1
        end do

        ! open file poscar/poscar_%08d.dat
        write(fid,'(i8.8)') out_id_poscar+simparams%start
        fname = 'poscar/poscar_'//fid//'.dat'
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


    ! proper poscar format
    subroutine output_vasp(atoms, itraj, istep)

        type(universe), intent(in)          :: atoms
        integer, intent(in)                 :: itraj, istep

        integer                             :: zero_step
        integer                             :: bead

        character(len=max_string_length)    :: fname
        character(len=8)                    :: traj_id_vasp
        character(len=7)                    :: step_id_vasp
        character(len=2)                    :: step_id_bead

        integer                             :: time_vals(8), noccurrences(atoms%ntypes), i, j
        real(dp)                            :: cents_r(3, atoms%natoms) ! for the beads to get center of mass for positions
        real(dp)                            :: cents_v(3, atoms%natoms) ! for the beads to get center of mass for velocities
        logical                             :: cents_fix(3, atoms%natoms) ! for the beads to get fixed atoms

        zero_step = 0

        if (.not. dir_exists('vasp')) call execute_command_line('mkdir vasp')

        ! prepare arrays
        noccurrences = 0
        do i = 1, atoms%natoms
            noccurrences(atoms%idx(i)) = noccurrences(atoms%idx(i)) + 1
        end do

        ! open file
        if (istep == -1) then
             write(step_id_vasp,'(i7.7)') zero_step
        else
             write(step_id_vasp,'(i7.7)') istep
        end if
        write(traj_id_vasp,'(i8.8)') itraj


        if (simparams%lbead_output_format) then

            do bead = 1,atoms%nbeads

                write(step_id_bead,'(i2.2)') bead
                fname = 'vasp/'//traj_id_vasp//'step'//step_id_vasp//'b'//step_id_bead//'.POSCAR'
                call open_for_write(out_unit, fname)

                ! write date and time as comment line
                call date_and_time(values=time_vals)
                write(out_unit, '(i4, a, i2.2, a, i2.2, a, i2.2, a, i2.2)') &
                    time_vals(1), "-", time_vals(2), "-",time_vals(3), " - ",time_vals(5), ".",time_vals(6)

                write(out_unit, '(a)')       '1.0'                    ! scaling constant
                write(out_unit, '(3f23.15)') atoms%simbox
                write(out_unit, *)           atoms%name
                write(out_unit, '(100i6)') (noccurrences(i), i=1,atoms%ntypes)

                if (atoms%is_cart) then
                    write(out_unit, '(a)') "Cartesian"
                else
                    write(out_unit, '(a)') "Direct"
                end if

                write(out_unit, '(3f23.15, 3l)') (atoms%r(:,bead,i), &
                    .not.atoms%is_fixed(:,bead,i), i=1,atoms%natoms)
                write(out_unit, '(a)') ""
                write(out_unit, '(3f23.15)') (atoms%v(:,bead,i), i=1,atoms%natoms)

                close(out_unit)

            end do

        else

            fname = 'vasp/'//traj_id_vasp//'step'//step_id_vasp//'.POSCAR'
            call open_for_write(out_unit, fname)

            ! write date and time as comment line
            call date_and_time(values=time_vals)
            write(out_unit, '(i4, a, i2.2, a, i2.2, a, i2.2, a, i2.2)') &
                time_vals(1), "-", time_vals(2), "-",time_vals(3), " - ",time_vals(5), ".",time_vals(6)

            write(out_unit, '(a)')       '1.0'                    ! scaling constant
            write(out_unit, '(3f23.15)') atoms%simbox
            write(out_unit, *)           atoms%name
            write(out_unit, '(100i6)') (noccurrences(i), i=1,atoms%ntypes)

            if (atoms%is_cart) then
                write(out_unit, '(a)') "Cartesian"
            else
                write(out_unit, '(a)') "Direct"
            end if


            cents_r   = calc_centroid_positions(atoms)
            cents_v   = calc_centroid_velocities(atoms)
            ! in case different beads have different notations of F and T we set F if at least one entry in all beads is F
            cents_fix = calc_centroid_is_fixed(atoms)

            ! positions and velocities
            write(out_unit, '(3f23.15, 3l)') (cents_r(:,i), .not. cents_fix(:,i), i=1,atoms%natoms)
            write(out_unit, '(a)') ""
            write(out_unit, '(3f23.15)') (cents_v(:,i), i=1,atoms%natoms)

            close(out_unit)

        end if

    end subroutine output_vasp


    subroutine output_aims(atoms, itraj, istep)

        type(universe), intent(in)          :: atoms
        integer, intent(in)                 :: itraj, istep

        integer                             :: zero_step
        integer                             :: bead

        character(len=max_string_length)    :: fname
        character(len=8)                    :: traj_id_aims
        character(len=7)                    :: step_id_aims
        character(len=2)                    :: step_id_beads

        integer                             :: time_vals(8), noccurrences(atoms%ntypes), i, j, k, l
        real(dp)                            :: cents_r(3, atoms%natoms) ! for the beads to get center of mass for positions
        real(dp)                            :: cents_v(3, atoms%natoms) ! for the beads to get center of mass for velocities

        logical, dimension(3)               :: lskip

        lskip(:) = .FALSE.

        zero_step = 0 ! write the initial structure to step 0 file

        if (.not. dir_exists('aims')) call execute_command_line('mkdir aims')

        noccurrences = 0
        do i = 1, atoms%natoms
            noccurrences(atoms%idx(i)) = noccurrences(atoms%idx(i)) + 1
        end do

        if (istep == -1) then ! initial structure
             write(step_id_aims,'(i7.7)') zero_step
        else
             write(step_id_aims,'(i7.7)') istep
        end if
        write(traj_id_aims,'(i8.8)') itraj

        if (simparams%lbead_output_format) then

            do bead = 1,atoms%nbeads

                write(step_id_beads,'(i2.2)') bead
                fname = 'aims/'//traj_id_aims//'step'//step_id_aims//'b'//step_id_beads//'.in'
                call open_for_write(out_unit, fname)

                ! write date and time as comment line
                call date_and_time(values=time_vals)
                write(out_unit, '(a,i4, a, i2.2, a, i2.2, a, i2.2, a, i2.2)') &
                    "# ", time_vals(1), "-", time_vals(2), "-",time_vals(3), " - ",time_vals(5), ".",time_vals(6)

                ! lattice
                do l=1,3
                    write(out_unit, '(a,3f23.15)') "lattice_vector",atoms%simbox(:,l)
                end do

                ! add empty line, for asthetics
                write(out_unit, '(a)') ""

                ! positions and velocities
                do i=1,atoms%natoms ! atoms
                    lskip(:) = .false.
                    write(out_unit, '(a,3f23.15,x,a)') ("atom",(atoms%r(:,bead,i), atoms%name(atoms%idx(i))))
                    do j=1,3
                        if (atoms%is_fixed(j,bead,i) == .true.) then
                            lskip(j) = .true.
                        end if
                    end do
                    if (all(lskip) == .true.) then
                        write(out_unit, '(a)') "constrain_relaxation .true."
                    else
                        do j=1,3
                            if (lskip(j) == .true.) then
                                if (j == 1) then
                                    write(out_unit, '(a)') "constrain_relaxation x"
                                else if (j == 2) then
                                    write(out_unit, '(a)') "constrain_relaxation y"
                                else ! j == 3, no other option possible
                                    write(out_unit, '(a)') "constrain_relaxation z"
                                end if
                            end if
                        end do
                        write(out_unit, '(a)',advance='no') "velocity"
                        do j=1,3
                            if (lskip(j) == .false.) then
                                write(out_unit, '(f23.15)',advance='no') atoms%v(j,bead,i)
                            end if
                        end do
                        write(out_unit, '(a)') ""
                    end if

                end do ! atoms

                close(out_unit)

            end do ! beads

        else

            fname = 'aims/'//traj_id_aims//'step'//step_id_aims//'.in'
            call open_for_write(out_unit, fname)

            ! write date and time as comment line
            call date_and_time(values=time_vals)
            write(out_unit, '(a,i4, a, i2.2, a, i2.2, a, i2.2, a, i2.2)') &
                "# ", time_vals(1), "-", time_vals(2), "-",time_vals(3), " - ",time_vals(5), ".",time_vals(6)

            ! lattice
            do l=1,3
                write(out_unit, '(a,3f23.15)') "lattice_vector",atoms%simbox(:,l)
            end do

            ! add empty line, for asthetics
            write(out_unit, '(a)') ""

            cents_r = calc_centroid_positions(atoms)
            cents_v = calc_centroid_velocities(atoms)

            ! positions and velocities
            do i=1,atoms%natoms ! atoms
                lskip(:) = .false.
                write(out_unit, '(a,3f23.15,x,a)') ("atom",(cents_r(:,i), atoms%name(atoms%idx(i))))
                do k=1,atoms%nbeads
                    do j=1,3
                        if (atoms%is_fixed(j,k,i) == .true.) then
                            lskip(j) = .true.
                        end if
                    end do
                end do
                if (all(lskip) == .true.) then
                    write(out_unit, '(a)') "constrain_relaxation .true."
                else
                    do j=1,3
                        if (lskip(j) == .true.) then
                            if (j == 1) then
                                write(out_unit, '(a)') "constrain_relaxation x"
                            else if (j == 2) then
                                write(out_unit, '(a)') "constrain_relaxation y"
                            else ! j == 3, no other option possible
                                write(out_unit, '(a)') "constrain_relaxation z"
                            end if
                        end if
                    end do
                    write(out_unit, '(a)',advance='no') "velocity"
                    do j=1,3
                        if (lskip(j) == .false.) then
                            write(out_unit, '(f23.15)',advance='no') cents_v(j,i)
                        end if
                    end do
                    write(out_unit, '(a)') ""
                end if

            end do ! atoms

            close(out_unit)

        end if

    end subroutine output_aims


    subroutine output_runner(atoms, itraj, istep)

        type(universe), intent(in)          :: atoms
        integer, intent(in)                 :: itraj, istep

        integer :: zero_step

        character(len=max_string_length)    :: fname

        integer                             :: time_vals(8), noccurrences(atoms%ntypes)
        integer                             :: j, k
        real(dp)                            :: cents_r(3, atoms%natoms) ! for the beads to get center of mass for positions
        real(dp)                            :: cents_f(3, atoms%natoms) ! for the beads to get center of mass for forces

        real(dp)                            :: dummy_ce, epot

        character(len=36)   :: lattice_format = '(A7, X, F11.8, X, F11.8, X, F11.8)'
        character(len=93)   :: atom_format = '(A4, X, F16.9, X, F16.9, X, F16.9, X, A3, X, F11.8, X, F11.8, X, F11.8, X, F11.8, X, F11.8)'
        character(len=20)   :: ce_format = '(A6 ,X , F11.8)'


        zero_step = 0 ! write the initial structure to step 0 file
        dummy_ce = 0.0_dp

        if (.not. dir_exists('runner')) call execute_command_line('mkdir runner')

        fname = 'runner/input.data'

        if (overwrite_runner) then
            call open_for_write(out_unit, fname)
            ! first comment line with program name
            write(out_unit, '(a)') '# This file was generated by MDT2.'
            ! second comment line with date and time
            call date_and_time(values=time_vals)
            write(out_unit, '(a,i4, a, i2.2, a, i2.2, a, i2.2, a, i2.2)') &
            "# ", time_vals(1), "-", time_vals(2), "-",time_vals(3), " - ",time_vals(5), ".",time_vals(6)
            overwrite_runner = .false.
        else
            call open_for_append(out_unit,fname)
        end if

        epot = sum(atoms%epot)/atoms%nbeads
        cents_r = calc_centroid_positions(atoms)
        cents_f = calc_centroid_forces(atoms)

        write (out_unit,'(A5)') 'begin'
        do k = 1,3
            write(out_unit,lattice_format) 'lattice', atoms%simbox(:,k) * ang2bohr
        end do
        do j = 1,atoms%natoms
            write (out_unit,atom_format) 'atom', cents_r(:,j) * ang2bohr, atoms%name(atoms%idx(j)), dummy_ce, dummy_ce, &
            cents_f(:,j) * evang2habohr
        end do
        write (out_unit,ce_format) 'charge', dummy_ce
        write (out_unit,ce_format) 'energy', epot * ev2ha
        write (out_unit,'(A3)') 'end'

        close(out_unit)

    end subroutine output_runner


    subroutine output_mxt(atoms)

        type(universe), intent(in) :: atoms

        character(len=max_string_length) :: fname
        character(len=8)                 :: fid

        if (.not. dir_exists('mxt')) call execute_command_line('mkdir mxt')

        ! open file mxt/mxt_%08d.dat
        write(fid,'(i8.8)') out_id_mxt+simparams%start
        fname = 'mxt/mxt_'//fid//'.dat'
        open(out_unit, file=fname, form="unformatted", status="replace")

        write(out_unit) atoms%natoms, atoms%nbeads, atoms%ntypes
        write(out_unit) atoms%r, atoms%v, atoms%is_cart, atoms%is_fixed, &
            atoms%idx, atoms%name, atoms%is_proj, atoms%simbox, atoms%isimbox

        close(out_unit)

    end subroutine output_mxt




    subroutine output_scatter(atoms, itraj, istep, flag)

        use rpmd
        use pes_nene_mod,  only : max_count_extrapolation_warnings, count_extrapolation_warnings_traj_energy, count_extrapolation_warnings_traj_symfunc

        type(universe), intent(in) :: atoms
        integer, intent(in) :: itraj, istep
        character(len=*), intent(in) :: flag

        real(dp) :: ekin_p, ekin_l, b_ekin_p, b_ekin_l, proj_v(3), bead_epot
        integer  :: turning_points
        character(len=max_string_length) :: fname, fin_name
        character(len=8)                 :: fin_id
        character(len=*), parameter :: err = "Error in output_scatter: "

        if (count(atoms%is_proj) > 1) stop err // "output_scatter does not work for more that one projectile, please implement"

        if (.not. dir_exists('traj')) call execute_command_line('mkdir traj')

        write(fname, '(a12, i8.8, a4)') 'traj/mxt_fin', itraj, '.dat'

        call atom_ekin(atoms, ekin_p, ekin_l)
        proj_v = sum(atoms%v(:,:,1), dim=2)/atoms%nbeads
        bead_epot = calc_bead_epot(atoms)
        call bead_ekin(atoms, b_ekin_p, b_ekin_l)
        call finalize_bounces()

        if (flag == "scatter_initial") then

            call open_for_write(out_unit, fname)

            write(out_unit, '(a11, f14.7)') "ekin_p_i = ", ekin_p
            write(out_unit, '(a11, f14.7)') "ekin_l_i = ", ekin_l
            write(out_unit, '(a11, e21.10)') "epot_i   = ", sum(atoms%epot)/atoms%nbeads
            write(out_unit, '(a11, e21.10)') "etotal_i = ", b_ekin_p + b_ekin_l + sum(atoms%epot)/atoms%nbeads + bead_epot
            write(out_unit, '(a11, 3f14.7)')"r_i      = ", sum(atoms%r(:,:,1), dim=2)/atoms%nbeads
            write(out_unit, '(a11, f14.7)') "polar_i  = ", 180-acos(proj_v(3)/sqrt(sum(proj_v*proj_v)))*rad2deg
            write(out_unit, '(a11, f14.7)') "azi_i    = ", atan2(proj_v(2), proj_v(1))*rad2deg
            write(out_unit, '(a)') ""

            if (any(atoms%pes == pes_id_nene)) then
                count_extrapolation_warnings_traj_energy  = 0 ! re-initialize for each trajectory
                count_extrapolation_warnings_traj_symfunc = 0 ! re-initialize for each trajectory
            end if

        else if (flag == "scatter_final") then

            call open_for_append(out_unit, fname)

            write(out_unit, '(a11, f14.7)') "ekin_p_f = ", ekin_p
            write(out_unit, '(a11, f14.7)') "ekin_l_f = ", ekin_l
            write(out_unit, '(a11, e21.10)') "epot_f   = ", sum(atoms%epot)/atoms%nbeads
            write(out_unit, '(a11, e21.10)') "etotal_f = ", b_ekin_p + b_ekin_l + sum(atoms%epot)/atoms%nbeads + bead_epot
            write(out_unit, '(a11, 3f14.7)')"r_f      = ", sum(atoms%r(:,:,1), dim=2)/atoms%nbeads
            write(out_unit, '(a11, f14.7)') "polar_f  = ", acos(proj_v(3)/sqrt(sum(proj_v*proj_v)))*rad2deg
            write(out_unit, '(a11, f14.7)') "azi_f    = ", atan2(proj_v(2), proj_v(1))*rad2deg
            write(out_unit, '(a)') ""
            write(out_unit, '(a11, f14.7)') "time     = ", (istep-1) * simparams%step
            write(out_unit, '(a11, i)')     "turn_pnts = ", nturning_points
            write(out_unit, '(a11, f14.7)') "cl_appr  = ", closest_approach
            write(out_unit, '(a11, 3f14.7)')"r_min_p  = ", lowest_z
            !write(out_unit, '(a11, i)')     "bounces  = ", bounces
            !write(out_unit, '(a11, f14.7)') "time_int = ", interaction_time

            ! to make sure that a configuration file is written even if no format keyword is given:
!            if (.not. dir_exists('mxt')) call execute_command_line('mkdir mxt')

            ! open file mxt/mxt_%08d.dat
!            write(fin_id,'(i8.8)') itraj
!            fin_name = 'mxt/mxt_fin_'//fin_id//'.dat'
!            open(fin_unit, file=fin_name, form="unformatted", status="replace")

!            write(fin_unit) atoms%natoms, atoms%nbeads, atoms%ntypes
!            write(fin_unit) atoms%r, atoms%v, atoms%is_cart, atoms%is_fixed, &
!                atoms%idx, atoms%name, atoms%is_proj, atoms%simbox, atoms%isimbox

!            close(fin_unit)

            if (any(atoms%pes == pes_id_nene)) then
                total_count_extrapolation_warnings = total_count_extrapolation_warnings + sum(count_extrapolation_warnings_traj_energy) + sum(count_extrapolation_warnings_traj_symfunc)
                if (sum(count_extrapolation_warnings_traj_energy) > 0)  print *, "Extrapolation warnings (energy) this traj: ", count_extrapolation_warnings_traj_energy
                if (sum(count_extrapolation_warnings_traj_symfunc) > 0) print *, "Extrapolation warnings (symmetry functions) this traj: ", count_extrapolation_warnings_traj_symfunc
                if (total_count_extrapolation_warnings > max_count_extrapolation_warnings) stop "Error: reached max number of extrapolation warnings!"
            end if

        else
            print *, err, "unknown flag", flag
            stop
        end if

        close (out_unit)

    end subroutine output_scatter



    subroutine output_pes(atoms)

        use pes_rebo_mod, only : to_string_rebo

        type(universe), intent(in) :: atoms

        integer :: i, j
        character(len=35) :: fname
        character(len=5) :: pes
        character(len=4) :: role_i, role_j
        character(len=*), parameter :: err = "Error in output_pes(): "

        write(fname, '(a9,i5.5,a12,i5.5,a4)') "fit/fits/", simparams%start, "/out_params/", simparams%start, ".pes"
        call open_for_write(out_unit, fname)

        do i = 1, atoms%ntypes
            do j = i, atoms%ntypes
                pes = pes_id_to_name(atoms%pes(i,j))
                role_i = universe_id_to_role(atoms, i)
                role_j = universe_id_to_role(atoms, j)

                write(out_unit, '(2a5)'), "pes ", pes
                write(out_unit, '(4a8)'), atoms%name(i), atoms%name(j), role_i, role_j

                select case (atoms%pes(i,j))
                    case(pes_id_rebo)
                        call to_string_rebo(out_unit, i, j)
                    case(pes_id_no_interaction)
                        ! no interaction, no parameters

                    ! TODO: implement outputs
                    !        case (pes_id_lj)
                    !            name = pes_name_lj
                    !        case (pes_id_morse)
                    !            name = pes_name_morse
                    !        case (pes_id_emt)
                    !            name = pes_name_emt
                    !        case (pes_id_ho)
                    !            name = pes_name_ho
                    !        case (pes_id_simple_lj)
                    !            name = pes_name_simple_lj
                    !	     case (pes_id_nene)
                    !		 name = pes_name_nene
                    case default
                        print *, err, "pes output not implemented for ", pes
                        stop
                end select

                ! if this is not the last pes
                if ( i /= atoms%ntypes .or. j /= atoms%ntypes) write(out_unit, '(a)') ""

            end do
        end do

        close(out_unit)

    end subroutine output_pes




    subroutine output_adsorption_status(atoms, itraj, istep)

        type(universe), intent(in) :: atoms
        integer, intent(in)        :: itraj, istep

        character(len=max_string_length) :: fname
        character(len=8)                 :: traj_id
        integer                          :: ads_output_step

        ! XXX: change system() to execute_command_line() when new compiler is available
        if (.not. dir_exists('conf')) call system('mkdir conf')

        ! open file conf/ads_%08d.dat
        write(traj_id,'(i8.8)') itraj
        fname = 'conf/ads_'//traj_id//'.dat'

        ! overwrite file on first write
        if (overwrite_ads) then
            call open_for_write(out_unit, fname)
            overwrite_ads = .false.
        else
            call open_for_append(out_unit, fname)
        end if

        write(out_unit, '(e13.6e2, L)') (istep-1)*simparams%step, is_adsorbed

        close(out_unit)

    end subroutine output_adsorption_status

end module output_mod
