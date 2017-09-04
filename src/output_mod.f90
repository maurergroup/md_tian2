module output_mod

    use universe_mod
    use useful_things
    use open_file
    use constants
    use run_config, only : simparams

    implicit none

    logical :: overwrite = .true.
    integer, parameter :: out_unit = 86


contains

    subroutine output(atoms, itraj, istep, format)

        type(universe), intent(in) :: atoms
        integer, intent(in)        :: itraj, istep, format

        character(len=*), parameter :: err = "Error in output(): "

        select case (format)
            case (format_xyz)
                call output_xyz(atoms, itraj)

            case (format_nrg)
                call output_nrg(atoms, itraj, istep)

            case (default_int)
                print *,  err // "output format not specified"
                stop

            case default
                print *, err // "unknown output format", format
                stop
        end select

        overwrite = .false.

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

        dir_coords(:,:,:) = dir_coords - anint(dir_coords-0.5)

        forall (i = 1 : atoms%natoms) cart_coords(:,:,i) = matmul(atoms%simbox, dir_coords(:,:,i))

        write(traj_id,'(I8.8)') itraj
        fname = 'conf/mxt_conf'//traj_id//'.xyz'

        if (overwrite) then
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

        ! time, epot, ekin_p, ekin_l, etotal,

        type(universe), intent(in) :: atoms
        integer, intent(in)        :: itraj, istep

        character(len=8)                 :: traj_id
        character(len=max_string_length) :: fname
        real(dp) :: avg_epot, ekin_p, ekin_l, etotal


        ! XXX: change system() to execute_command_line() when new compiler is available
        if (.not. dir_exists('conf')) call system('mkdir conf')

        write(traj_id,'(I8.8)') itraj
        fname = 'conf/mxt_trj'//traj_id//'.dat'

        if (overwrite) then
            call open_for_write(out_unit, fname)
            write(out_unit, '(6a17)') 'time/fs', 'epot/eV', 'ekin_p/eV', 'ekin_l/eV', 'e_total/eV', 'f_total'
        else
            call open_for_append(out_unit,fname)
        end if

        avg_epot = sum(atoms%epot)/atoms%nbeads

        call simple_ekin(atoms, ekin_p, ekin_l)
        etotal = ekin_p + ekin_l + avg_epot

        write(out_unit, '(6e17.8e2)') istep*simparams%step, avg_epot, ekin_p, ekin_l, etotal, sum(atoms%f)


        close(out_unit)


    end subroutine output_nrg

end module output_mod
