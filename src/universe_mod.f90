module universe_mod
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


    !  Type atoms
    !   structure to hold the position, velocity, force etc. for multiple atoms
    !       use rank2 array so positions, velocities, forces etc. are stored
    !       in sequentional memory locations for efficient access
    !       each array should be allocated (3, n_beads, n_atom)
    !       mass array has length of n_atoms
    type universe

        integer                       :: natoms          ! number of atoms
        integer                       :: nbeads          ! number of beads per atom
        integer                       :: ntypes          ! proj species + latt species
        real(dp),         allocatable :: m(:)            ! mass
        character(len=3), allocatable :: name(:)         ! atomic name
        real(dp), allocatable         :: r(:,:,:)        ! positions
        real(dp), allocatable         :: v(:,:,:)        ! velocities
        real(dp), allocatable         :: f(:,:,:)        ! forces
        real(dp), allocatable         :: a(:,:,:)        ! acceleration
        logical,  allocatable         :: is_fixed(:,:,:) ! mask array defining frozen atoms (T is frozen)

        integer, allocatable          :: idx(:)          ! index of atom type used in program
        integer, allocatable          :: pes(:,:)        ! defines idx-dependent pes
        logical, allocatable          :: is_proj(:)      ! defines projectile and lattice via idx
        real(dp), allocatable         :: epot(:)

        real(dp)                      :: simbox(3,3)     ! simulation cell
        real(dp)                      :: isimbox(3,3)    ! inverse simulation cell
        logical                       :: is_cart         ! if geometry is cartesian or direct

    end type universe



contains

    ! Constructor for type atoms
    !    input:  n_beads, n_atoms
    !    allocates arrays as (3,n_beads,n_atom)

    function new_atoms(nbeads, natoms, ntypes)
        integer, intent(in) :: nbeads, natoms, ntypes
        type(universe) new_atoms

        allocate(new_atoms%r(3,nbeads,natoms))
        allocate(new_atoms%v(3,nbeads,natoms))
        allocate(new_atoms%f(3,nbeads,natoms))
        allocate(new_atoms%a(3,nbeads,natoms))
        allocate(new_atoms%is_fixed(3,nbeads,natoms))

        allocate(new_atoms%idx(natoms))

        allocate(new_atoms%epot(nbeads))

        allocate(new_atoms%pes(ntypes,ntypes))
        allocate(new_atoms%name(ntypes))
        allocate(new_atoms%m(ntypes))
        allocate(new_atoms%is_proj(ntypes))

        new_atoms%nbeads = nbeads
        new_atoms%natoms = natoms
        new_atoms%ntypes = ntypes
        new_atoms%name  = default_string
        new_atoms%r     = default_real
        new_atoms%v     = default_real
        new_atoms%f     = default_real
        new_atoms%m     = default_real
        new_atoms%a     = default_real
        new_atoms%is_fixed = .false.
        new_atoms%is_proj  = .false.
        new_atoms%is_cart  = .false.
        new_atoms%idx   = default_int
        new_atoms%pes   = default_int
        new_atoms%simbox  = default_real
        new_atoms%isimbox = default_real
    end function


    subroutine to_cartesian(this)

        type(universe), intent(inout) :: this
        integer :: i,j

        if (this%is_cart) then
            print *,"Cannot convert cartesian to cartesian coordinates" ; call abort
        end if

        do i = 1, this%natoms
            do j = 1, this%nbeads
                this%r(:,j,i) = matmul(this%simbox, this%r(:,j,i))
            end do
        end do

        this%is_cart = .not. this%is_cart

    end subroutine


    subroutine to_direct(this)

        type(universe), intent(inout) :: this
        integer :: i,j

        if (.not. this%is_cart) then
            print *,"Cannot convert direct to direct coordinates" ; call abort
        end if

        do i = 1, this%natoms
            do j = 1, this%nbeads
                this%r(:,j,i) = matmul(this%isimbox, this%r(:,j,i))
            end do
        end do

        this%is_cart = .not. this%is_cart

    end subroutine


    integer function get_idx_from_name(this, name, is_proj) result(idx)

        type(universe), intent(in) :: this
        character(len=3), intent(in) :: name
        logical, intent(in) :: is_proj

        integer :: i

        idx = default_int
        do i = 1, this%natoms
            if (this%name(i) == name .and. this%is_proj(this%idx(i)) == is_proj) then
                idx = this%idx(i)
                exit
            end if
        end do

        if (idx == default_int) stop "Error in get_idx_from_name(): make sure you & 
                        correctly assign element names to projectile and slab in both &
                        *.inp and *.pes files."        

    end function get_idx_from_name


    subroutine create_repetitions(this, rep)

        use useful_things, only : invert_matrix

        type(universe), intent(inout) :: this
        integer, intent(in) :: rep(2)

        type(universe) :: other
        integer :: nreps
        integer :: i, j, start, end

        if (this%is_cart) call to_direct(this)

        ! new atoms
        nreps = (1+2*rep(1)) * (1+2*rep(2))
        other = new_atoms(this%nbeads, nreps*this%natoms, this%ntypes)
        other%m = this%m
        other%name = this%name
        other%is_proj = this%is_proj

        do i = -rep(1), rep(1)
            do j = -rep(2), rep(2)
                ! This map -i, -i+1, ... , i-1, i and -j, -j+1, ... , j-1, j
                !  to 1, 2, 3, ... , #total_repetitions
                start = (rep(1) * (i+rep(1)) * (1+2*rep(2)) + rep(2) + j)     * this%natoms + 1
                end   = (rep(1) * (i+rep(1)) * (1+2*rep(2)) + rep(2) + j + 1) * this%natoms

                other%r( 1 , :, start : end ) = this%r(1,:,:) + i
                other%r( 2 , :, start : end ) = this%r(2,:,:) + j
                other%r( 3 , :, start : end ) = this%r(3,:,:)


                other%idx(      start : end ) = this%idx(:)
                other%v( : , :, start : end ) = this%v(:,:,:)
                other%f( : , :, start : end ) = this%v(:,:,:)
                other%is_fixed( : , :, start : end ) = this%is_fixed(:,:,:)

            end do
        end do

        ! bring back to cartesian
        call to_cartesian(other)

        ! enlarge the cell
        other%simbox(:,1) = this%simbox(:,1) * (1+2*rep(1))
        other%simbox(:,2) = this%simbox(:,2) * (1+2*rep(2))

        ! obtain new inverse cell matrix
        other%isimbox = invert_matrix(other%simbox)

        ! replace old atoms
        this = other

    end subroutine create_repetitions


    subroutine pbc_distdir(atoms, i, j, b, r, vec)
            !
            ! Purpose: Distance between atoms a and b and vector a-->b
            !          with taking into account the periodic boundary conditions
            !

        type(universe), intent(in)         :: atoms
        integer, intent(in)                :: i, j, b
        real(8),               intent(out) :: r
        real(8), dimension(3), intent(out) :: vec

        vec = atoms%r(:,b,j)-atoms%r(:,b,i)   ! distance vector from a to b

        vec = matmul(atoms%isimbox, vec)   ! transform to direct coordinates
        vec = vec - Anint(vec)       ! imaging
        vec = matmul(atoms%simbox, vec)    ! back to cartesian coordinates

        r =  sqrt(sum(vec*vec))      ! distance

    end subroutine pbc_distdir


    subroutine set_acceleration(atoms)

        type(universe), intent(inout) :: atoms
        integer :: i

        do i = 1, atoms%natoms
            atoms%a(:,:,i) = atoms%f(:,:,i)/atoms%m(atoms%idx(i))
        end do

    end subroutine set_acceleration





end module universe_mod
