module atom_class
    !
    ! Purpose:
    !           This module contains definitions of user types and all constants
    !
    ! Date          	Author          	    History of Revison
    ! ====          	======          	    ==================
    ! 30.03.2017    	Marvin Kammler		    new data structure
    !                   Sascha Kandratsenka
    ! 18.02.2014    	Svenja M. Janke		    Original
    !			        Sascha Kandratsenka
    !			        Dan J. Auerbach

    use constants

    implicit none

    real(dp) ::  simbox(3,3) = default_real     ! simulation cell
    real(dp) :: isimbox(3,3) = default_real     ! inverse simulation cell



    !  Type atoms
    !   structure to hold the position, velocity, force etc. for multiple atoms
    !       use rank2 array so positions, velocities, forces etc. are stored
    !       in sequentional memory locations for efficient access
    !       each array should be allocated (3, n_beads, n_atom)
    !       mass array has length of n_atoms
    type atoms

        integer                       :: natoms          ! number of atoms
        integer                       :: nbeads          ! number of beads per atom
        real(dp),         allocatable :: m(:)            ! mass
        character(len=3), allocatable :: name(:)         ! atomic name
        real(dp), allocatable         :: r(:,:,:)        ! positions
        real(dp), allocatable         :: v(:,:,:)        ! velocities
        real(dp), allocatable         :: f(:,:,:)        ! forces
        logical,  allocatable         :: is_fixed(:,:,:) ! mask array defining frozen atoms (T is frozen)
        integer,  allocatable         :: idx(:)          ! index of atom type used in program
        integer,  allocatable         :: pes(:,:)        ! defines idx-dependent pes
        logical                       :: is_cart            ! if geometry is cartesian or direct

    end type atoms

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
        real(dp),         allocatable :: einc(:)                    ! incidence energy (eV)
        real(dp),         allocatable :: inclination(:)             ! incidence polar angle (degree)
        real(dp),         allocatable :: azimuth(:)                 ! incidence azimuthal angle (degree)
        real(dp):: Tsurf                                            ! surface temperature in K
        real(dp):: sa_Tmax                                          ! max. Tsurf for simulated annealing in K
        integer :: sa_nsteps                                        ! number of steps per simulated annealing cycle
        integer :: sa_interval                                      ! number of steps per temperature interval
        character(len=7)    :: confname                             ! configuration key
        character(len=max_string_length) :: confname_file           ! name of the system configuration file or folder
        integer :: rep(2)                                           ! defines in-plane repetitions
        integer :: nconfs                                           ! number of configurations to read in
        character(len=max_string_length) :: pes_file                ! name of the file that stores the potential parameters
        character(len=3)  :: run                                    ! what to do
        integer           :: output(2)                              ! what to save
        character(len=15) :: pip(3)                                 ! determine initial projectile position

    end type

contains

    ! Constructor for type atoms
    !    input:  n_beads, n_atoms
    !    allocates arrays as (3,n_beads,n_atom)

    function new_atoms(nbeads, natoms)
        integer, intent(in) :: nbeads, natoms
        type(atoms) new_atoms

        allocate(new_atoms%name(natoms))
        allocate(new_atoms%r(3,nbeads,natoms))
        allocate(new_atoms%v(3,nbeads,natoms))
        allocate(new_atoms%f(3,nbeads,natoms))
        allocate(new_atoms%is_fixed(3,nbeads,natoms))
        allocate(new_atoms%idx(natoms))

        new_atoms%nbeads = nbeads
        new_atoms%natoms = natoms
        new_atoms%name  = default_string
        new_atoms%r     = default_real
        new_atoms%v     = default_real
        new_atoms%f     = default_real
        new_atoms%is_fixed = .false.
        new_atoms%idx   = default_int


    end function

    ! Constructor for type simulation_parameters

    function new_simulation_parameters()

        type(simulation_parameters) new_simulation_parameters

        new_simulation_parameters%start  = default_int
        new_simulation_parameters%ntrajs = default_int
        new_simulation_parameters%nsteps = default_int
        new_simulation_parameters%step   = default_real
        new_simulation_parameters%nlattices = default_int
        new_simulation_parameters%nprojectiles = default_int
        new_simulation_parameters%Tsurf   = default_real
        new_simulation_parameters%sa_Tmax   = default_real
        new_simulation_parameters%sa_nsteps = default_int
        new_simulation_parameters%sa_interval = default_int
        new_simulation_parameters%confname = default_string
        new_simulation_parameters%rep = [0,0]
        new_simulation_parameters%nconfs  = default_int
        new_simulation_parameters%pes_file = default_string
        new_simulation_parameters%run = default_string
        new_simulation_parameters%output = [default_int,default_int]
        new_simulation_parameters%pip = default_string

    end function


    subroutine to_cartesian(this)

        type(atoms), intent(inout) :: this
        integer :: i,j

        if (this%is_cart) then
            print *,"Cannot convert cartesian to cartesian coordinates" ; call abort
        end if

        do i = 1, this%natoms
            do j = 1, this%nbeads
                this%r(:,j,i) = matmul(simbox, this%r(:,j,i))
            end do
        end do

        this%is_cart = .not. this%is_cart

    end subroutine


    subroutine to_direct(this)

        type(atoms), intent(inout) :: this
        integer :: i,j

        if (.not. this%is_cart) then
            print *,"Cannot convert direct to direct coordinates" ; call abort
        end if

        do i = 1, this%natoms
            do j = 1, this%nbeads
                this%r(:,j,i) = matmul(isimbox, this%r(:,j,i))
            end do
        end do

        this%is_cart = .not. this%is_cart

    end subroutine


    integer function get_idx_from_name(this, name) result(idx)

        type(atoms), intent(in) :: this
        character(len=3), intent(in) :: name

        integer :: i

        idx = default_int
        do i = 1, this%natoms
            if (this%name(i) == name) then
                idx = this%idx(i)
                exit
            end if
        end do

        if (idx == default_int) stop "Error in get_idx_from_name(): make sure you & 
                        correctly assign element names to projectile and slab in both &
                        *.inp and *.pes files."

    end function get_idx_from_name


    subroutine create_repetitions(inslab, rep)

        use useful_things, only : invert_matrix

        type(atoms), intent(inout) :: inslab
        integer, intent(in) :: rep(2)

        type(atoms) :: out_slab
        integer :: nreps
        integer :: i, j, start, end

        if (.not. allocated(inslab%m)) stop "Error in create_repetitions(): mass attribute not allocated"

        ! new atoms
        nreps = (1+2*rep(1)) * (1+2*rep(2))
        out_slab = new_atoms(inslab%nbeads, nreps*inslab%natoms)
        allocate(out_slab%m(size(inslab%m)))
        out_slab%m = inslab%m

        do i = -rep(1), rep(1)
            do j = -rep(2), rep(2)
                ! This map -i, -i+1, ... , i-1, i and -j, -j+1, ... , j-1, j
                !  to 1, 2, 3, ... , #total_repetitions
                start = (rep(1) * (i+rep(1)) * (1+2*rep(2)) + rep(2) + j)     * inslab%natoms + 1
                end   = (rep(1) * (i+rep(1)) * (1+2*rep(2)) + rep(2) + j + 1) * inslab%natoms

                out_slab%r( 1 , :, start : end ) = inslab%r(1,:,:) + i
                out_slab%r( 2 , :, start : end ) = inslab%r(2,:,:) + j
                out_slab%r( 3 , :, start : end ) = inslab%r(3,:,:)

                out_slab%name(     start : end ) = inslab%name(:)
                out_slab%idx(      start : end ) = inslab%idx(:)
                out_slab%v( : , :, start : end ) = inslab%v(:,:,:)
                out_slab%f( : , :, start : end ) = inslab%v(:,:,:)
                out_slab%is_fixed( : , :, start : end ) = inslab%is_fixed(:,:,:)

            end do
        end do

        ! bring back to cartesian
        call to_cartesian(out_slab)

        ! enlarge the cell
        simbox(:,1) = simbox(:,1) * (1+2*rep(1))
        simbox(:,2) = simbox(:,2) * (1+2*rep(2))

        ! obtain new inverse cell matrix
        isimbox = invert_matrix(simbox)

    end subroutine create_repetitions





end module
