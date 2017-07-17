module run_config
    !
    ! Purpose:
    !    Initialize simulation parameters
    !
    ! Date          	Author          	History of Revison
    ! ====          	======          	==================
    ! 31.03.2017    	Marvin Kammler		    original
    !                   Sascha Kandratsenka

use atom_class, only : dp, amu2mass, deg2rad
use useful_things, only : lower_case

implicit none

    save

    integer  :: start_tr     = 1     ! a trajectory to start with
    integer  :: ntrajs       = 10    ! number of trajectories
    real(dp) :: step         =-1.0_dp! time step in fs
    integer  :: nsteps       = 100   ! number of steps
    integer  :: nlattices    = 0     ! number of lattice species
    integer  :: nprojectiles = 0     ! number of incident species
    character(len=3), dimension(:), allocatable :: name_l, name_p ! atomic names
    real(dp),         dimension(:), allocatable :: mass_l, mass_p ! atomic masses
    character(len=3), dimension(:), allocatable :: mdpa_name_l, mdpa_name_p ! propagation algs
    integer,          dimension(:), allocatable :: md_algo_l, md_algo_p     ! and respective key

    logical :: start_init=.false., ntrajs_init=.false., nsteps_init=.false.

    real(dp) :: einc         = 0     ! incidence energy (eV)
    real(dp) :: inclination  = 0     ! incidence polar angle (degree)
    real(dp) :: azimuth      = 0     ! incidence azimuthal angle (degree)

    character(len=1000) :: confname_file
    character(len=7)    :: confname
    integer :: n_confs = 1          ! Number of configurations to read in
    integer :: conf_nr = 1          ! number in the name of configurational file to read in.
    character(len=10)  :: fitnum    ! number of fit
    integer, dimension(2) :: rep = (/0,0/)  ! number of in-plane repetitions

    real(dp) :: Tsurf        = 300_dp ! surface temperature (Kelvin)
    logical  :: Tsurf_key = .false.
    real(dp) :: Tmax = 0_dp           ! Annealing parameters
    integer :: sasteps = 0

    character(len=1000) :: pes_file=''


contains

subroutine read_input_file(input_file)

character(*), intent(in) :: input_file

integer :: i, ios = 0, line = 0, pos1
character(len=1000) :: buffer, label!, biffer


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
            ! Find the first instance of whitespace.  Split label and data.
            buffer = adjustl(buffer)
            pos1 = scan(buffer, ' ')
            label = buffer(1:pos1)
            buffer = buffer(pos1+1:)

            select case (label)

            case('start')
                if (start_init) stop 'Error in the input file: Multiple use of the start key'
                start_init = .true.
                read(buffer,*,iostat=ios) start_tr

            case('ntrajs')
                if (ntrajs_init) stop 'Error in the input file: Multiple use of the ntrajs key'
                ntrajs_init = .true.
                read(buffer,*,iostat=ios) ntrajs

            case('step')
                if (step > 0) stop 'Error in the input file: Multiple use of the step key'
                read(buffer,*,iostat=ios) step

            case('nsteps')
                if (nsteps_init) stop 'Error in the input file: Multiple use of the nsteps key'
                nsteps_init = .true.
                read(buffer,*,iostat=ios) nsteps

            case ('lattice')
                if (nlattices > 0) stop 'Error in the input file: Multiple use of the lattice key'
                read(buffer, *, iostat=ios) nlattices
                allocate(name_l(nlattices), mass_l(nlattices), &
                         mdpa_name_l(nlattices), md_algo_l(nlattices))
                read(buffer, *, iostat=ios) nlattices, &
                                (name_l(i), mass_l(i), mdpa_name_l(i), i=1,nlattices)
                if (ios /= 0) then
                    print*, 'Error in the input file: inconsistent lattice definition'
                    print*, '                         nlattices not equal to the number of species'
                    stop
                end if
                mass_l = mass_l*amu2mass
                do i=1,nlattices
                    call lower_case(mdpa_name_l(i))
                    select case (mdpa_name_l(i))
                        case ('ver')                ! verlet
                            md_algo_l(i) = 1
                        case ('bee')                ! beeman
                            md_algo_l(i) = 2
                        case ('lan')                ! langevin
                            md_algo_l(i) = 3
                        case ('sla')                ! langevin (series)
                            md_algo_l(i) = 4
                        case ( 'pef' )              ! verlet + post Electronic friction
                            md_algo_l(i) = 5
                       case default
                            print *, 'algorithm ', mdpa_name_l(i),&
                                     ' for lattice species ', name_l(i), ' unknown'
                            stop
                    end select
                end do

            case ('projectile')
                if (nprojectiles > 0) stop 'Error in the input file: Multiple use of the projectile key'
                read(buffer, *, iostat=ios) nprojectiles
                allocate(name_p(nprojectiles), mass_p(nprojectiles), &
                         mdpa_name_p(nprojectiles), md_algo_p(nprojectiles))
                read(buffer, *, iostat=ios) nprojectiles, &
                                (name_p(i), mass_p(i), mdpa_name_p(i), i=1,nprojectiles)
                if (ios /= 0) then
                    print*, 'Error in the input file: inconsistent projectile definition'
                    print*, '                         nprojectiles not equal to the number of species'
                    stop
                end if
                mass_p = mass_p*amu2mass
                do i=1,nprojectiles
                    call lower_case(mdpa_name_p(i))
                    select case (mdpa_name_p(i))
                        case ('ver')                ! verlet
                            md_algo_p(i) = 1
                        case ('bee')                ! beeman
                            md_algo_p(i) = 2
                        case ('lan')                ! langevin
                            md_algo_p(i) = 3
                        case ('sla')                ! langevin (series)
                            md_algo_p(i) = 4
                        case ( 'pef' )              ! verlet + post Electronic friction
                            md_algo_p(i) = 5
                       case default
                            print *, 'algorithm ', mdpa_name_p(i),&
                                     ' for projectile species ', name_p(i), ' unknown'
                            stop
                    end select
                end do

            case('incidence')
                read(buffer,*,iostat=ios) einc, inclination, azimuth
                if (ios.ne.0) stop 'Error in the input file: incidence conditions'
                inclination = inclination*deg2rad
                azimuth = azimuth*deg2rad

            case ('conf')
                read(buffer, *, iostat=ios) confname, confname_file
                call lower_case(confname)
                select case (confname)
                    case ('geo')
                        read(buffer, *, iostat=ios) confname, confname_file, conf_nr
                    case ('mxt')
                        read(buffer, *, iostat=ios) confname, confname_file, n_confs
                    case ('fit')
                        read(buffer, *, iostat=ios) confname, confname_file, n_confs, fitnum
                    case ('poscar')
                    case default
                        stop 'Error in the input file: unknown conf keyword.'
                end select
                if (ios.ne.0) stop 'Error in the input file: ill-defined configuration'
                call lower_case(confname)

            case ('rep')

                read(buffer, *, iostat=ios) rep
                if (ios /= 0) then

                    read(buffer, *, iostat=ios) rep(1)
                    if (ios == 0) then
                        rep(2) = rep(1)
                        print *, 'Warning: You have specified a single number for cell &
                                           repetitions.'
                        print *, '         Using the same repetition in both directions.'
                    else
                        stop 'Error in the input file: cell repetitions ill-defined.'
                   end if

                end if

            case('Tsurf')
                read(buffer,*,iostat=ios) Tsurf
                if (ios /= 0) stop 'Error in the input file: Tsurf not properly set'
                Tsurf_key = .true.

            case('anneal')
                read(buffer, *, iostat=ios) Tmax, sasteps
                if (ios /= 0) stop 'Error in the input file: anneal not properly set'

            case ('pes')
                read(buffer, *, iostat=ios) pes_file
                if (ios /= 0) stop 'Error in the input file: pes file not specified'

            case default
                if (trim(label) /= '' .and. label(1:1) /= '!') &
                        print *, 'Skipping invalid label at line', line, trim(label)

            end select
        end if
    end do ! ios

if (pes_file == '') stop 'Error in the input file: keyword pes absent'

stop 101

end subroutine read_input_file

end module run_config
