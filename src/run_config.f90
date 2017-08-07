module run_config
        !
        ! Purpose:
        !    Initialize simulation parameters
        !
        ! Date          	Author          	History of Revison
        ! ====          	======          	==================
        ! 31.03.2017    	Marvin Kammler		    original
        !                   Sascha Kandratsenka

    use constants
    use atom_class
    use useful_things, only : lower_case, split_string
    use open_file, only : open_for_read

    implicit none

    type(simulation_parameters) :: simparams


contains

    subroutine read_input_file(input_file)

        character(len=*), intent(in) :: input_file

        integer :: i, ios = 0, line = 0, nwords
        character(len=max_string_length) :: buffer
        character(len=max_string_length) :: words(100)

        simparams = new_simulation_parameters()

        !------------------------------------------------------------------------------
        !                       READ IN INPUT FILE
        !                       ==================
        !------------------------------------------------------------------------------

        call open_for_read(38, input_file)
        ! ios < 0: end of record condition encountered or endfile condition detected
        ! ios > 0: an error is detected
        ! ios = 0  otherwise

        do while (ios == 0)

            read(38, '(A)', iostat=ios) buffer
            if (ios == 0) then

                line = line + 1
                ! Split an input string
                call split_string(buffer, words, nwords)

                select case (words(1))


                    case ('run')

                        if (simparams%run /= default_string) stop 'Error in the input file: Multiple use of the run key'
                        if (nwords == 2) then
                            read(words(2),'(A)') simparams%run
                        else
                            stop 'Error in the input file: run key needs a single argument'
                        end if


                    case('start')

                        if (simparams%start /= default_int) stop 'Error in the input file: Multiple use of the start key'
                        if (nwords == 2) then
                            read(words(2),'(i1000)', iostat=ios) simparams%start
                            if (ios /= 0)  stop 'Error in the input file: start value must be integer'
                        else
                            stop 'Error in the input file: start key needs a single argument'
                        end if


                    case('ntrajs')

                        if (simparams%ntrajs /= default_int)   stop 'Error in the input file: Multiple use of the ntrajs key'
                        if (nwords == 2) then
                            read(words(2),'(i1000)', iostat=ios) simparams%ntrajs
                            if (ios /= 0) stop 'Error in the input file: ntrajs value must be integer'
                        else
                            stop 'Error in the input file: ntrajs key needs a single argument'
                        end if


                    case('nsteps')

                        if (simparams%nsteps /= default_int)   stop 'Error in the input file: Multiple use of the nsteps key'
                        if (nwords == 2) then
                            read(words(2),'(i1000)', iostat=ios) simparams%nsteps
                            if (ios /= 0) stop 'Error in the input file: nsteps value must be integer'
                        else
                            stop 'Error in the input file: nsteps key needs a single argument'
                        end if


                    case('step')

                        if (simparams%step /= default_real)   stop 'Error in the input file: Multiple use of the step key'
                        if (nwords == 2) then
                            read(words(2),*, iostat=ios) simparams%step
                            if (ios /= 0) stop 'Error in the input file: step value must be a number'
                        else
                            stop 'Error in the input file: step key needs a single argument'
                        end if


                    case ('lattice')

                        if (simparams%nlattices /= default_int)   stop 'Error in the input file: Multiple use of the lattice key'
                        if (nwords > 1) then
                            read(words(2),'(i1000)', iostat=ios) simparams%nlattices
                            if (ios /= 0) stop 'Error in the input file: lattice key 1st argument value must be integer'
                            if (nwords == 3*simparams%nlattices + 2) then
                                allocate(simparams%name_l(simparams%nlattices), &
                                    simparams%mass_l(simparams%nlattices), &
                                    simparams%md_algo_l(simparams%nlattices))
                                simparams%name_l = words(3:nwords:3)
                                do i=1,simparams%nlattices
                                    read(words(1+3*i),*, iostat=ios) simparams%mass_l(i)
                                    if (ios /= 0) stop 'Error in the input file: lattice key - mass value is not a number'
                                    call lower_case(words(2+3*i))
                                    select case (words(2+3*i))
                                        case ('ver')                ! verlet
                                            simparams%md_algo_l(i) = prop_verlet
                                        case ('bee')                ! beeman
                                            simparams%md_algo_l(i) = prop_beeman
                                        case ('lan')                ! langevin
                                            simparams%md_algo_l(i) = prop_langevin
                                        case ('sla')                ! langevin (series)
                                            simparams%md_algo_l(i) = prop_langevin_series
                                        case default
                                            print *, 'algorithm ', trim(words(2+3*i)),&
                                                ' for lattice species ', trim(simparams%name_l(i)), ' unknown'
                                            stop
                                    end select
                                end do
                            else
                                stop 'Error in the input file: lattice key has wrong number of arguments'
                            end if

                        else
                            stop 'Error in the input file: lattice key needs arguments'
                        end if


                    case ('projectile')

                        if (simparams%nprojectiles /= default_int)   stop 'Error in the input file: Multiple use of the projectile key'
                        if (nwords > 1) then
                            read(words(2),'(i1000)', iostat=ios) simparams%nprojectiles
                            if (ios /= 0) stop 'Error in the input file: projectile key 1st argument value must be integer'
                            if (nwords == 3*simparams%nprojectiles + 2) then
                                allocate(simparams%name_p(simparams%nprojectiles),&
                                    simparams%mass_p(simparams%nprojectiles),&
                                    simparams%md_algo_p(simparams%nprojectiles))
                                simparams%name_p = words(3:nwords:3)
                                do i=1,simparams%nprojectiles
                                    read(words(1+3*i),*, iostat=ios) simparams%mass_p(i)
                                    if (ios /= 0) stop 'Error in the input file: projectile key - mass value is not a number'
                                    call lower_case(words(2+3*i))
                                    select case (words(2+3*i))
                                        case ('ver')                ! verlet
                                            simparams%md_algo_p(i) = prop_verlet
                                        case ('bee')                ! beeman
                                            simparams%md_algo_p(i) = prop_beeman
                                        case ('lan')                ! langevin
                                            simparams%md_algo_p(i) = prop_langevin
                                        case ('sla')                ! langevin (series)
                                            simparams%md_algo_p(i) = prop_langevin_series
                                        case default
                                            print *, 'algorithm ', trim(words(2+3*i)),&
                                                ' for projectile species ', trim(simparams%name_p(i)), ' unknown'
                                            stop
                                    end select
                                end do
                            else
                                stop 'Error in the input file: projectile key has wrong number of arguments'
                            end if

                        else
                            stop 'Error in the input file: projectile key needs arguments'
                        end if


                    case('incidence')

                        if (allocated(simparams%einc)) stop 'Error in the input file: Multiple use of the incidence key'
                        if (nwords > 1 .and. mod(nwords - 1, 3) == 0) then
                            allocate(   simparams%einc((nwords-1)/3),&
                                simparams%inclination((nwords-1)/3),&
                                simparams%azimuth((nwords-1)/3)      )
                            do i=0,(nwords-1)/3 -1
                                read(words(2+3*i),*, iostat=ios) simparams%einc(i+1)
                                if (ios /= 0) stop 'Error in the input file: incidence key - einc value must be a number'
                                read(words(3+3*i),*, iostat=ios) simparams%inclination(i+1)
                                if (ios /= 0) stop 'Error in the input file: incidence key - inclination value must be a number'
                                read(words(4+3*i),*, iostat=ios) simparams%azimuth(i+1)
                                if (ios /= 0) stop 'Error in the input file: incidence key - azimuth value must be a number'
                            end do
                        else
                            stop 'Error in the input file: incidence key needs (3 x n_projectiles) arguments'
                        end if


                    case('Tsurf')

                        if (simparams%Tsurf /= default_real) stop 'Error in the input file: Multiple use of the Tsurf key'
                        if (nwords == 2) then
                            read(words(2),*, iostat=ios) simparams%Tsurf
                            if (ios /= 0) stop 'Error in the input file: Tsurf value must be a number'
                        else
                            stop 'Error in the input file: Tsurf key needs a single argument'
                        end if


                    case('annealing')

                        if (simparams%sa_Tmax /= default_real) stop 'Error in the input file: Multiple use of the annealing key'
                        if (nwords == 4) then
                            read(words(2),*, iostat=ios) simparams%sa_Tmax
                            if (ios /= 0) stop 'Error in the input file: annealing key - sa_Tmax value must be a number'
                            read(words(3),'(i1000)', iostat=ios) simparams%sa_nsteps
                            if (ios /= 0) stop 'Error in the input file: annealing key - sa_nsteps value must be integer'
                            read(words(4),'(i1000)', iostat=ios) simparams%sa_interval
                            if (ios /= 0) stop 'Error in the input file: annealing key - sa_interval value must be integer'
                        else
                            stop 'Error in the input file: annealing key needs 3 arguments'
                        end if


                    case ('conf')

                        if (simparams%confname /= default_string)   stop 'Error in the input file: Multiple use of the conf key'
                        if (nwords > 2) then
                            read(words(2),'(A)') simparams%confname
                            call lower_case(simparams%confname)
                            read(words(3),'(A)') simparams%confname_file

                            select case (simparams%confname)

                                case ('poscar')

                                    if (nwords > 5) stop 'Error in the input file: conf key - poscar argument number is too large'
                                    if (words(4) /= "") then
                                        read(words(4),'(i1000)',iostat=ios) simparams%rep(1)
                                        if (ios /= 0) stop 'Error in the input file: conf key - poscar repetition arguments must be integer'
                                        if (words(5) /= "") then
                                            read(words(5),'(i1000)',iostat=ios) simparams%rep(2)
                                            if (ios /= 0) stop 'Error in the input file: conf key - poscar repetition arguments must be integer'
                                        else
                                            simparams%rep(2) = simparams%rep(1)
                                        end if
                                    end if

                                case ('mxt')
                                    ! conf mxt <mxt_file.dat>: use this mxt file as system config
                                    ! conf mxt <mxt_folder> <n>: randomly select mxt files 1<=x<=n from folder
                                    if (nwords == 4) then
                                        read(words(4),'(i1000)',iostat=ios) simparams%nconfs
                                        if (ios /= 0) stop 'Error in the input file: conf key - mxt argument must be integer'
                                    end if
                                    if (nwords < 4) stop 'Error in the input file: conf key - too few mxt arguments'
                                    if (nwords > 4) stop 'Error in the input file: conf key - too many mxt arguments'


                                case default
                                    stop 'Error in the input file: conf key - unknown conf name'
                            end select
                        else
                            stop 'Error in the input file: conf key needs at least 2 arguments'
                        end if


                    case ('output')

                        if (any(simparams%output /= default_int))   stop 'Error in the input file: Multiple use of the output key'
                        if (words(2) /= "") then
                            read(words(2),'(i1000)',iostat=ios) simparams%output(1)
                            if (ios /= 0) stop 'Error in the input file: 1st output key must be integer'
                            if (words(3) /= "") then
                                read(words(3),'(i1000)',iostat=ios) simparams%output(2)
                                if (ios /= 0) stop 'Error in the input file: 2nd output key must be integer'
                            end if
                        else
                            stop 'Error in the input file: output key has no arguments'
                        end if
                        if (nwords > 3) stop 'Error in the input file: output key has too many arguments'


                    case ('pip')

                        if (any(simparams%pip /= default_string)) stop 'Error in the input file: Multiple use of the pip key'
                        if (nwords == 4) then
                            read (words(2:4), '(A)') simparams%pip
                        else
                            stop 'Error in the input file: pip key needs 3 arguments'
                        end if


                    case ('pes')

                        if (simparams%pes_file /= default_string) stop 'Error in the input file: Multiple use of the pes key'
                        read(words(2), '(A)', iostat=ios) simparams%pes_file
                        if (ios /= 0) stop 'Error in the input file: pes file not specified'
                    !TODO: add keywords for fit


                    case default
                        if (trim(words(1)) /= '' .and. words(1)(1:1) /= '!') &
                            print *, 'Skipping invalid label ',trim(words(1)),' in line', line

                end select
            end if
        end do ! ios

    end subroutine read_input_file


    function get_atomic_indices(ident_line)

        ! Each atom has a unique index 1, 2, ..., n, where n is the number of species in
        !  the system. It starts with the projectiles and ends with the lattice atoms.
        !  The order is given in the *.inp file. It has to match the order in the geometry file.

        character(len=max_string_length), intent(in) :: ident_line(4)
        integer :: get_atomic_indices(2)

        integer :: i, j

        get_atomic_indices = -1

        do i = 1, 2
            if (ident_line(i+2) == "proj") then
                do j = 1, simparams%nprojectiles
                    if (ident_line(i) == simparams%name_p(j)) then
                        get_atomic_indices(i) = j
                    end if
                end do
            else if (ident_line(i+2) == "latt") then
                do j = 1, simparams%nlattices
                    if (ident_line(i) == simparams%name_l(j)) then
                        get_atomic_indices(i) = j + simparams%nprojectiles
                    end if
                end do
            end if
        end do

        if (any(get_atomic_indices == -1)) stop "Mismatch between atomic types in .inp file and PES file."

    end function get_atomic_indices

end module run_config
