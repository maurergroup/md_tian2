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
    use universe_mod, only : universe
    use useful_things, only : lower_case, split_string
    use open_file, only : open_for_read

    implicit none

    type simulation_parameters

        integer :: start                                            ! a trajectory to start with
        integer :: ntrajs                                           ! number of trajectories
        integer :: nsteps                                           ! number of steps
        real(dp):: step                                             ! time step in fs
        integer :: nlattices                                        ! number of lattice species
        integer :: nprojectiles                                     ! number of incident species
        character(len=3), allocatable :: name_l(:), name_p(:)       ! atomic names
        real(dp),         allocatable :: mass_l(:), mass_p(:)       ! atomic masses
        integer,          allocatable :: md_algo_l(:), md_algo_p(:) ! and respective key
        real(dp) :: einc, sigma_einc                                ! incidence energy (eV)
        real(dp) :: polar                                     ! incidence polar angle (degree)
        character(len=15) :: azimuth                                ! incidence azimuthal angle (degree) or 'r' for random
        real(dp) :: sigma_azimuth                                   ! standard dev of incidence azimuthal angle (degree)
        real(dp) :: Tsurf                                            ! surface temperature in K
        real(dp) :: Tproj                                            ! projectile temperature in K
        real(dp) :: sa_Tmax                                          ! max. Tsurf for simulated annealing in K
        integer  :: sa_nsteps                                        ! number of steps per simulated annealing cycle
        integer  :: sa_interval                                      ! number of steps per temperature interval
        character(len=7)    :: confname                             ! configuration key
        character(len=max_string_length) :: confname_file           ! name of the system configuration file or folder
        integer :: nconfs                                           ! number of configurations to read in
        character(len=max_string_length) :: merge_proj_file         ! name of folder containing projecile configuration
        integer :: merge_proj_nconfs                                ! number of projectile configurations
        integer :: rep(2)                                           ! defines in-plane repetitions
        character(len=max_string_length) :: pes_file                ! name of the file that stores the potential parameters
        character(len=3)  :: run                                    ! what to do
        integer, allocatable :: output_type(:)                      ! what to save
        integer, allocatable :: output_interval(:)                  ! when to save
        character(len=15) :: pip(3)                                 ! determine initial projectile position
        real(dp) :: andersen_time                                   ! average time between collsions per atom
        real(dp) :: pile_tau                                        ! PILE thermostat centroid mode thermostat time constant
        integer  :: force_beads                                     ! inititalizes all atoms with this many beads
        real(dp) :: proj_ul                                         ! stop trajectory if projectile's z-coordinate is higher that this value
        character(len=max_string_length) :: fit_training_folder     ! path to training data
        character(len=max_string_length) :: fit_validation_folder   ! path to validation data
        integer  :: fit_training_data                               ! number of configuration/energy pair files in fit_training_folder
        integer  :: fit_validation_data                             ! number of configuration/energy pair files in fit_validation_folder
        real(dp) :: evasp                                           ! reference energy for fit
        integer  :: maxit                                           ! maximum number of iteration during fit
        integer  :: nthreads                                        ! number of threads used for fitting

    end type

    type(simulation_parameters) :: simparams


contains

    function new_simulation_parameters()

        type(simulation_parameters) new_simulation_parameters

        new_simulation_parameters%start  = default_int
        new_simulation_parameters%ntrajs = default_int
        new_simulation_parameters%nsteps = default_int
        new_simulation_parameters%step   = default_real
        new_simulation_parameters%nlattices = default_int
        new_simulation_parameters%nprojectiles = default_int
        new_simulation_parameters%Tsurf   = default_real
        new_simulation_parameters%Tproj   = default_real
        new_simulation_parameters%einc  = default_real
        new_simulation_parameters%sigma_einc = default_real
        new_simulation_parameters%polar = default_real
        new_simulation_parameters%azimuth = default_string
        new_simulation_parameters%sigma_azimuth = default_real
        new_simulation_parameters%sa_Tmax   = default_real
        new_simulation_parameters%sa_nsteps = default_int
        new_simulation_parameters%sa_interval = default_int
        new_simulation_parameters%confname = default_string
        new_simulation_parameters%confname_file = default_string
        new_simulation_parameters%merge_proj_file = default_string
        new_simulation_parameters%merge_proj_nconfs = default_int
        new_simulation_parameters%rep = [0,0]
        new_simulation_parameters%nconfs  = default_int
        new_simulation_parameters%pes_file = default_string
        new_simulation_parameters%run = default_string
        new_simulation_parameters%pip = default_string
        new_simulation_parameters%andersen_time = 30.0_dp
        new_simulation_parameters%pile_tau = 200.0_dp
        new_simulation_parameters%force_beads = default_int
        new_simulation_parameters%proj_ul = default_real
        new_simulation_parameters%fit_training_data = default_int
        new_simulation_parameters%fit_validation_data = default_int
        new_simulation_parameters%fit_training_folder = default_string
        new_simulation_parameters%fit_validation_folder = default_string
        new_simulation_parameters%evasp = default_real
        new_simulation_parameters%maxit = 30
        new_simulation_parameters%nthreads = 1

    end function

    subroutine read_input_file(input_file)

        character(len=*), intent(in) :: input_file

        integer :: i, ios = 0, line = 0, nwords, randk
        real(dp) :: rnd
        character(len=max_string_length) :: buffer
        character(len=max_string_length) :: words(100)
        integer, parameter :: inp_unit = 38
        character(len=*), parameter :: err = "Error in the input file: "

        simparams = new_simulation_parameters()

        !------------------------------------------------------------------------------
        !                       READ IN INPUT FILE
        !                       ==================
        !------------------------------------------------------------------------------

        call open_for_read(inp_unit, input_file)
        ! ios < 0: end of record condition encountered or endfile condition detected
        ! ios > 0: an error is detected
        ! ios = 0  otherwise

        ! find the start key, to rotate random number generator
        do while (ios == 0)
            read(inp_unit, '(A)', iostat=ios) buffer
            if (ios == 0) then
                ! Split an input string
                call split_string(buffer, words, nwords)
                select case (words(1))
                    case('start')
                        if (simparams%start /= default_int) stop 'Error in the input file: Multiple use of the start key'
                        if (nwords == 2) then
                            read(words(2),'(i1000)', iostat=ios) simparams%start
                            if (ios /= 0) stop err // 'start value must be integer'
                        else
                            print *, err, 'start key needs a single argument'; stop
                        end if
                end select
            end if
        end do

        close(inp_unit)

        ! size of seed for random number generator
        randk=size(randseed)
        call random_seed(size=randk)
        call random_seed(put=randseed)

        ! rotate rng
        do i = 1, 100*simparams%start
            call random_number(rnd)
        end do

        ! read rest of the input file
        call open_for_read(inp_unit, input_file); ios = 0

        do while (ios == 0)

            read(inp_unit, '(A)', iostat=ios) buffer
            if (ios == 0) then
                line = line + 1
                ! Split an input string
                call split_string(buffer, words, nwords)

                select case (words(1))

                    case ('start')
                        ! pass

                    case ('run')

                        if (simparams%run /= default_string) stop err // "Multiple use of the run key"
                        if (nwords == 2) then
                            read(words(2),'(A)') simparams%run
                        else
                            print *, err, "run key needs a single argument"; stop
                        end if


                    case('ntrajs')

                        if (simparams%ntrajs /= default_int)   stop 'Error in the input file: Multiple use of the ntrajs key'
                        if (nwords == 2) then
                            read(words(2),'(i1000)', iostat=ios) simparams%ntrajs
                            if (ios /= 0) stop 'Error in the input file: ntrajs value must be integer'
                        else
                            print *, err, 'ntrajs key needs a single argument'; stop
                        end if


                    case('nsteps')

                        if (simparams%nsteps /= default_int)   stop 'Error in the input file: Multiple use of the nsteps key'
                        if (nwords == 2) then
                            read(words(2),'(i1000)', iostat=ios) simparams%nsteps
                            if (ios /= 0) stop 'Error in the input file: nsteps value must be integer'
                        else
                            print *, err, 'nsteps key needs a single argument'
                        end if


                    case('step')

                        if (simparams%step /= default_real) stop err // 'Multiple use of the step key'
                        if (nwords == 2) then
                            read(words(2),*, iostat=ios) simparams%step
                            if (ios /= 0) stop err // 'step value must be a number'
                        else
                            print *, err // 'step key needs a single argument'; stop
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
                                            simparams%md_algo_l(i) = prop_id_verlet
                                        !                                        case ('bee')                ! beeman
                                        !                                            simparams%md_algo_l(i) = prop_id_beeman
                                        case ('lan')                ! langevin
                                            simparams%md_algo_l(i) = prop_id_langevin
                                        !                                        case ('sla')                ! langevin (series)
                                        !                                            simparams%md_algo_l(i) = prop_id_langevin_series
                                        case ('and')
                                            simparams%md_algo_l(i) = prop_id_andersen
                                        case ('pil')
                                            simparams%md_algo_l(i) = prop_id_pile
                                        case default
                                            print *, 'algorithm ', trim(words(2+3*i)),&
                                                ' for lattice species ', trim(simparams%name_l(i)), ' unknown'
                                            stop
                                    end select
                                end do
                            else
                                print *, err, 'lattice key has wrong number of arguments'; stop
                            end if

                        else
                            print *, err, 'lattice key needs arguments'; stop
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
                                            simparams%md_algo_p(i) = prop_id_verlet
                                        !                                        case ('bee')                ! beeman
                                        !                                            simparams%md_algo_p(i) = prop_id_beeman
                                        case ('lan')                ! langevin
                                            simparams%md_algo_p(i) = prop_id_langevin
                                        !                                        case ('sla')                ! langevin (series)
                                        !                                            simparams%md_algo_p(i) = prop_id_langevin_series
                                        case ('and')
                                            simparams%md_algo_p(i) = prop_id_andersen
                                        case ('pil')
                                            simparams%md_algo_p(i) = prop_id_pile
                                        case default
                                            print *, 'algorithm ', trim(words(2+3*i)),&
                                                ' for projectile species ', trim(simparams%name_p(i)), ' unknown'
                                            stop
                                    end select
                                end do
                            else
                                print *, err, 'projectile key has wrong number of arguments'; stop
                            end if

                        else
                            print *, err, 'projectile key needs arguments'; stop
                        end if


                    case ('Einc')

                        if (simparams%einc /= default_real) stop err // 'Multiple use of the Einc key'
                        if (nwords == 2) then
                            read(words(2), *, iostat=ios) simparams%einc
                            if (ios /= 0) stop err // 'incidence energy must be a number'
                            simparams%sigma_einc = 0.0_dp
                        else if (nwords == 3) then
                            read(words(2:3), *, iostat=ios) simparams%einc, simparams%sigma_einc
                            if (ios /= 0) stop err // 'einc keyword only takes numerical arguments'
                        else
                            print *, err // 'einc key needs one or two arguments'; stop
                        end if


                    case ('polar')

                        if (simparams%polar /= default_real) stop err // 'Multiple use of the polar key'
                        if (nwords == 2) then
                            read(words(2),*, iostat=ios) simparams%polar
                            if (ios /= 0) stop err // 'polar value must be a number'
                            simparams%polar = simparams%polar * deg2rad
                        else
                            print *, err // 'polar key needs a single argument'; stop
                        end if


                    case ('azimuth')

                        if (simparams%azimuth /= default_string) stop err // 'Multiple use of the azimuth key'
                        if (nwords == 2) then
                            read(words(2),'(A)') simparams%azimuth
                            simparams%sigma_azimuth = 0.0_dp
                        else if (nwords == 3) then
                            read(words(2),'(A)') simparams%azimuth
                            read(words(3), *, iostat=ios) simparams%sigma_azimuth
                            if (ios /= 0) stop err // 'standard dev of azimuth angle must be a number'
                        else
                            print *, err // 'azimuth key needs one or two arguments'; stop
                        end if


                    case ('Tsurf')

                        if (simparams%Tsurf /= default_real) stop err // 'Multiple use of the Tsurf key'
                        if (nwords == 2) then
                            read(words(2),*, iostat=ios) simparams%Tsurf
                            if (ios /= 0) stop err // 'Tsurf value must be a number'
                        else
                            print *, err, 'Tsurf key needs a single argument'; stop
                        end if


                    case ('Tproj')

                        if (simparams%Tproj /= default_real) stop err // 'Multiple use of the Tproj key'
                        if (nwords == 2) then
                            read(words(2),*, iostat=ios) simparams%Tproj
                            if (ios /= 0) stop err // 'Tproj value must be a number'
                        else
                            print *, err, 'Tproj key needs a single argument'; stop
                        end if


                    case('annealing')

                        if (simparams%sa_Tmax /= default_real) stop err // 'Multiple use of the annealing key'
                        if (nwords == 4) then
                            read(words(2),*, iostat=ios) simparams%sa_Tmax
                            if (ios /= 0) stop err // 'annealing key - sa_Tmax value must be a number'
                            read(words(3),'(i1000)', iostat=ios) simparams%sa_nsteps
                            if (ios /= 0) stop err // 'annealing key - sa_nsteps value must be integer'
                            read(words(4),'(i1000)', iostat=ios) simparams%sa_interval
                            if (ios /= 0) stop err // 'annealing key - sa_interval value must be integer'
                        else
                            stop 'Error in the input file: annealing key needs 3 arguments'
                        end if


                    case ('conf')

                        if (simparams%confname /= default_string) stop err // 'Multiple use of the conf key'
                        if (nwords > 2) then
                            read(words(2),'(A)') simparams%confname
                            call lower_case(simparams%confname)

                            select case (simparams%confname)

                                case ('poscar')
                                    ! conf poscar <poscar_file> <x_rep> <y_rep>
                                    read(words(3),'(A)') simparams%confname_file
                                    if (nwords > 5) stop err // 'conf key - poscar argument number is too large'
                                    if (words(4) /= "") then
                                        read(words(4),'(i1000)',iostat=ios) simparams%rep(1)
                                        if (ios /= 0) stop err // 'conf key - poscar repetition arguments must be integer'
                                        if (words(5) /= "") then
                                            read(words(5),'(i1000)',iostat=ios) simparams%rep(2)
                                            if (ios /= 0) stop err // 'conf key - poscar repetition arguments must be integer'
                                        else
                                            simparams%rep(2) = simparams%rep(1)
                                        end if
                                    end if

                                case ('mxt')
                                    ! conf mxt <mxt_file.dat>: use this mxt file as system config
                                    ! conf mxt <mxt_folder> <n>: randomly select mxt files 1<=x<=n from folder
                                    read(words(3),'(A)') simparams%confname_file
                                    if (nwords == 4) then
                                        read(words(4),'(i1000)',iostat=ios) simparams%nconfs
                                        if (ios /= 0) stop err // 'conf key - mxt argument must be integer'
                                        call random_number(rnd)
                                        write(simparams%confname_file, '(2a, i8.8, a)') trim(simparams%confname_file), "/mxt_", int(rnd*simparams%nconfs)+1, ".dat"
                                    end if
                                    if (nwords < 3) stop err // 'conf key - too few mxt arguments'
                                    if (nwords > 4) stop err // 'conf key - too many mxt arguments'

                                case ('merge')
                                    ! conf mergewith <mxt_folder> <n>: randomly select mxt files 1<=x<=n from folder
                                    if (nwords /= 6) stop err // "conf merge needs projectile mxt folder, # of projectile configurations therein &
                                        lattice mxt folder and # of lattice configurations therein"
                                    ! projectile
                                    read(words(3),'(A)') simparams%merge_proj_file
                                    read(words(4),'(i1000)',iostat=ios) simparams%merge_proj_nconfs
                                    if (ios /= 0) stop err // "conf key - number of configurations must be integer"
                                    call random_number(rnd)
                                    write(simparams%merge_proj_file, '(2a, i8.8, a)') &
                                        trim(simparams%merge_proj_file), "/mxt_", int(rnd*simparams%merge_proj_nconfs)+1, ".dat"

                                    ! slab
                                    read(words(5),'(A)') simparams%confname_file
                                    read(words(6),'(i1000)',iostat=ios) simparams%nconfs
                                    if (ios /= 0) stop err // "conf key - number of configurations must be integer"
                                    call random_number(rnd)
                                    write(simparams%confname_file, '(2a, i8.8, a)') &
                                        trim(simparams%confname_file), "/mxt_", int(rnd*simparams%nconfs)+1, ".dat"


                                case default
                                    stop 'Error in the input file: conf key - unknown conf name'

                            end select
                        else
                            stop 'Error in the input file: conf key needs at least 2 arguments'
                        end if


                    case ('output')

                        if (allocated(simparams%output_type) .or. allocated(simparams%output_interval)) &
                            stop 'Error in the input file: Multiple use of the output key'
                        if (modulo(nwords, 2) /= 1) stop 'Error in the input file: number of output arguments must be even'
                        allocate(simparams%output_type((nwords-1)/2), simparams%output_interval((nwords-1)/2))
                        do i = 1, (nwords-1)/2
                            select case (words(2*i))
                                case (output_key_xyz)
                                    simparams%output_type(i) = output_id_xyz
                                case (output_key_energy)
                                    simparams%output_type(i) = output_id_energy
                                case (output_key_poscar_bead)
                                    simparams%output_type(i) = output_id_poscar_bead
                                case (output_key_poscar_true)
                                    simparams%output_type(i) = output_id_poscar_true
                                case (output_key_mxt)
                                    simparams%output_type(i) = output_id_mxt
                                case (output_key_scatter)
                                    simparams%output_type(i) = output_id_scatter
                                case default
                                    print *, 'Error in the input file: output format ', trim(words(2*i)), ' unknown'
                                    stop
                            end select

                            read (words(2*i+1), '(i1000)', iostat=ios) simparams%output_interval(i)
                            if (ios /= 0) stop 'Error in the input file: output interval must be integer'
                        end do


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


                    case ('andersen_time')

                        read(words(2), *, iostat=ios) simparams%andersen_time
                        if (ios /= 0) stop 'Error in the input file: Error reading Andersen collision time'


                    case ('pile_tau')

                        read(words(2), *, iostat=ios) simparams%pile_tau
                        if (ios /= 0) stop 'Error in the input file: Error reading PILE tau'


                    case ('force_beads')

                        if (simparams%force_beads /= default_int) stop 'Error in the input file: Multiple use of the force_beads key'
                        read(words(2), '(i1000)', iostat=ios) simparams%force_beads
                        if (ios /= 0) stop 'Error in the input file: Error reading force_beads'


                    case ('pul')

                        if (simparams%proj_ul /= default_real) stop err // "projectile upper limit set multiple times"
                        read(words(2), *, iostat=ios) simparams%proj_ul
                        if (ios /= 0) stop err // "proj_upper_limit"


                    case ('fit_training_data')

                        if (simparams%fit_training_data /= default_int) stop err // "fit training data set multiple times"
                        if (nwords /= 3) stop err // "fit_training_data key needs 2 arguments"
                        read(words(2), '(i)', iostat=ios) simparams%fit_training_data
                        if (ios /= 0) stop err // "Error reading number of training data points"
                        read(words(3), '(a)', iostat=ios) simparams%fit_training_folder
                        if (ios /= 0) stop err // "Error reading path to training data points"


                    case ('fit_validation_data')

                        if (simparams%fit_validation_data /= default_int) stop err // "fit validation data set multiple times"
                        if (nwords /= 3) stop err // "fit_validation_data key needs 2 arguments"
                        read(words(2), '(i)', iostat=ios) simparams%fit_validation_data
                        if (ios /= 0) stop err // "Error reading number of validation data points"
                        read(words(3), '(a)', iostat=ios) simparams%fit_validation_folder
                        if (ios /= 0) stop err // "Error reading path to validation data points"

                    case ('evasp')

                        if (simparams%evasp /= default_real) stop err // "reference energy for fit set multiple times"
                        read(words(2), *, iostat=ios) simparams%evasp
                        if (ios /= 0) stop err // "Error reading evasp"


                    case ('maxit')

                        read(words(2), *, iostat=ios) simparams%maxit
                        if (ios /= 0) stop err // "Error reading maxit"


                    case ('nthreads')

                        read(words(2), *, iostat=ios) simparams%nthreads
                        if (ios /= 0) stop err // "Error reading number of threads"


                    case default
                        if (trim(words(1)) /= '' .and. words(1)(1:1) /= '!') &
                            print *, 'Skipping invalid label ',trim(words(1)),' in line', line

                end select
            end if
        end do ! ios
        close(inp_unit)

    end subroutine read_input_file




end module run_config
