module output_mod

    use universe_mod
    use useful_things
    use open_file
    use constants

    implicit none

contains

    subroutine output(atoms, itraj, format)

        type(universe), intent(inout) :: atoms
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

    end subroutine output





    subroutine output_xyz(atoms, itraj)

        type(universe), intent(inout) :: atoms
        integer, intent(in) :: itraj

        integer, parameter :: xyz_unit = 86
        character(len=8) :: traj_id
        character(len=max_string_length) :: fname
        integer :: i, j
        real(dp) :: coords(3)

        ! XXX: change system() to execute_command_line() when new compiler is available
        if (.not. dir_exists('conf')) call system('mkdir conf')

        write(traj_id,'(I8.8)') itraj
        fname = 'conf/mxt_conf'//traj_id//'.xyz'

        call to_direct(atoms)

        call open_for_append(xyz_unit,fname)

        write(xyz_unit,*) atoms%natoms*atoms%nbeads
        write(xyz_unit,*)

        do i = 1, atoms%natoms
            do j = 1, atoms%nbeads

                coords = atoms%r(:,j,i)-anint(atoms%r(:,j,i)-0.5)
                write(xyz_unit,'(a5, 3f9.6)') atoms%name(atoms%idx(i)), coords

            end do
        end do

        call to_cartesian(atoms)

        close(xyz_unit)

    end subroutine output_xyz
end module output_mod
