module output_mod

    use universe_mod
    use useful_things
    use open_file
    use constants

    implicit none

    logical :: overwrite = .true.

contains

    subroutine output(atoms, itraj, format)

        type(universe), intent(in) :: atoms
        integer, intent(in)           :: itraj, format

        character(len=*), parameter :: err = "Error in output(): "

        select case (format)
            case (format_xyz)
                call output_xyz(atoms, itraj)

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

        integer,          parameter :: xyz_unit   = 86
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
            call open_for_write(xyz_unit, fname)
        else
            call open_for_append(xyz_unit,fname)
        end if

        write(xyz_unit,*) atoms%natoms*atoms%nbeads
        write(xyz_unit,*)

        do i = 1, atoms%natoms
            do j = 1, atoms%nbeads
                write(xyz_unit, xyz_format) atoms%name(atoms%idx(i)), cart_coords(:,j,i)
            end do
        end do

        close(xyz_unit)

    end subroutine output_xyz
end module output_mod
