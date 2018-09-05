module pes_nene_mod

    use constants
    use useful_things, only : split_string, lower_case
    use universe_mod

    implicit none


    character(len=max_string_length)   :: prog_path, inp_path
    character(len=max_string_length)   :: input_string, outforces_string, outenergy_string, prog_string, grep_string


    contains

    !later: Read in here the keywords for the high-dimensional neural network potentials (HDNNPs)
    subroutine read_nene(atoms, inp_unit)

        use run_config, only : simparams

        type(universe), intent(inout) :: atoms
        integer, intent(in) :: inp_unit

        integer :: nwords, ios = 0
        character(len=max_string_length) :: buffer
        character(len=max_string_length) :: words(100)
        integer  :: idx1, idx2, ntypes
        character(len=*), parameter :: err = "Error in read_nene_ext: "

        ntypes = simparams%nprojectiles+simparams%nlattices

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
        atoms%pes(idx1,idx2) = pes_id_nene
        atoms%pes(idx2,idx1) = pes_id_nene


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

            ! insert here a readout from the pes file based on keywords in file input.nn for RuNNer when all the source code is implemented
            select case (words(1))

                case ('inp_dir')

                    read(words(2), '(A)') inp_path

                case ('prog_dir')

                    read(words(2), '(A)') prog_path


                case default
                    print *, "Error in the PES file: unknown nene_ext parameter", words(1)
                    stop

            end select

        end do


        input_string = trim(inp_path) // "input.data"
        outforces_string = trim(inp_path) // "nnforces.out"
        outenergy_string = trim(inp_path) // "energy.out"
        prog_string = "cd " // trim(inp_path) // " && " // trim(prog_path) // " > temp"
        grep_string = "grep WARNING " // trim(inp_path) // "temp"


    end subroutine read_nene


    subroutine compute_nene(atoms, flag)

        ! Calculates energy and forces with HDNNPs

        type(universe), intent(inout)   :: atoms
        integer, intent(in)             :: flag

        integer                         :: ios, eios, cios, j, k
        real(dp)                        :: dummy_ce

        ios = 0
        eios = 0
        cios = 0
        dummy_ce = 0.00000000_dp

        !write input.data for RuNNer
        open(unit=21,file=input_string,status='replace',action='write',iostat=ios)

            if (ios == 0) then

                write (21,'(A5)')'begin'

                do k = 1,3

                    write(21,lattice_format) 'lattice', atoms%simbox(:,k) * ang2bohr

                end do

                do j = 1,atoms%natoms


                    write (21,atom_format) 'atom', atoms%r(:,:,j) * ang2bohr, atoms%name(atoms%idx(j)), dummy_ce, dummy_ce, atoms%f(:,:,j)

                end do

                write (21,ce_format) 'charge', dummy_ce
                write (21,ce_format) 'energy', dummy_ce
                write (21,'(A3)')'end'

            else

                write (*,*) 'Error writing input.data for RuNNer: Check if directory is correct! iostat = ', ios
                stop

            end if

        close(unit=21)

        !execute RuNNer 
        call execute_command_line(prog_string, exitstat=eios, cmdstat=cios)
        call execute_command_line(grep_string)

            if (eios == 0 .and. cios == 0) then

                !Read out Forces from nnforces.out calculated by RuNNer
                open(unit=22,file=outforces_string,status='old',action='read',form='formatted',iostat=ios)

                    if (ios == 0) then

                        read(22,*) atoms%f(:,:,:)

                        atoms%f(:,:,:) = atoms%f(:,:,:) * forceconv

                    else

                        write (*,*) 'Error reading out forces from RuNNer: iostat = ', ios
                        stop

                    end if

                close(unit=22)

                !Read out Epot from energy.out calculated by RuNNer
                open(unit=23,file=outenergy_string,status='old',action='read',form='formatted',iostat=ios)

                    if (ios == 0) then

                        read(23,*) atoms%epot

                        atoms%epot = atoms%epot * ha2ev

                    else

                        write(*,*) 'Error reading out energy from RuNNer: iostat = ', ios
                        stop

                    end if

                close(unit=23)

            else

                write (*,*) 'Error executing RuNNer: Check if directory is correct! ', 'exitstat: ', eios, ' cmdstat: ', cios
                stop

            end if

    end subroutine compute_nene

end module pes_nene_mod
