module md_init
    !
    ! Purpose:
    !    Initialize simulations.
    !    Should be able to do the following things:
    !      1.
    !
    ! Date          	Author          	History of Revison
    ! ====          	======          	==================
    ! 30.03.2017    	Marvin Kammler		    new data structure
    !                   Sascha Kandratsenka     and propagation mathods
    !
    ! 18.02.2014    	Svenja M. Janke		    Original
    !			        Sascha Kandratsenka
    !			        Dan J. Auerbach

    use universe_mod
    use useful_things
    use run_config
    use constants

    use pes_lj_mod
    use pes_ho_mod
    use pes_emt_mod
    use pes_non_mod
    implicit none

contains

    subroutine simbox_init(atoms)

        type(universe), intent(out) :: atoms

        character(len=:), allocatable :: input_file
        integer :: input_file_length, input_file_status
        integer :: randk

        ! Read in name of input file

        if (command_argument_count() == 0) stop " I need an input file"

        call get_command(length=input_file_length)
        allocate(character(input_file_length) :: input_file)
        call get_command_argument(1, input_file, input_file_length, input_file_status)
        if (input_file_status /= 0) stop " Error by reading the command line"


        ! size of seed for random number generator
        randk=size(randseed)
        call random_seed(size=randk)

        call read_input_file(input_file)        ! build the simparams object
        call ensure_input_sanity()

        call read_geometry(atoms, simparams%confname_file)
        call ensure_geometry_sanity(atoms)

        call read_pes(atoms)
        call remove_com_velocity(atoms)
        call post_process(atoms)


    !TODO: Do not forget to convert the angles from degrees to program units

    end subroutine simbox_init


    subroutine read_pes(atoms)

        use open_file, only : open_for_read

        type(universe), intent(inout) :: atoms

        integer :: ios = 0, nwords, i
        integer, parameter :: pes_unit = 38
        character(len=max_string_length) :: buffer
        character(len=max_string_length) :: words(100)
        character(len=*), parameter :: err = "Error in read_pes(): "

        if (.not. file_exists(simparams%pes_file)) stop err // "PES file does not exist"

        call open_for_read(pes_unit, simparams%pes_file)
        ! ios < 0: end of record condition encountered or endfile condition detected
        ! ios > 0: an error is detected
        ! ios = 0  otherwise

        do while (ios == 0)

            read(pes_unit, '(A)', iostat=ios) buffer
            if (ios == 0) then

                ! Split an input string
                call split_string(buffer, words, nwords)

                if (words(1) == "pes" .and. nwords == 2) then

                    select case (words(2))

                        case ('lj')
                            call read_lj(atoms, pes_unit)

                        case ('slj')
                            call read_simple_lj(atoms, pes_unit)

                        !                        case ('morse')
                        !                            call read_morse(pes_unit)

                        case ('emt')
                            call read_emt(atoms, pes_unit)
                        !
                        !                        case ('rebo')
                        !                            call read_rebo(pes_unit)

                        case ('non')
                            call read_non_interacting(atoms, pes_unit)

                        case ('ho')
                            call read_ho(atoms, pes_unit)

                        case default
                            print *, err // "unknown potential in PES file:", words(2)
                            stop

                    end select

                end if
            end if

        end do

        close(pes_unit)

        if (any(atoms%pes == default_int)) then
            print *, err // "pes matrix incomplete"
            print *, (atoms%pes(:,i), i=1,atoms%ntypes)
            stop
        end if

    end subroutine read_pes



    ! After this subroutine, the system will be initialized with all user-specified geometry
    !   and velocity information.
    subroutine read_geometry(atoms, infile)

        use useful_things, only : file_exists

        character(len=*), intent(in) :: infile
        type(universe), intent(inout) :: atoms

        if (.not. file_exists(infile)) stop "Error: geometry file does not exist"

        select case (simparams%confname)

            case ("poscar")
                call read_poscar(atoms, infile)

            case default
                stop "Error: conf keyword unknown"

        end select

    end subroutine read_geometry


    subroutine read_poscar(atoms, infile)

        use open_file, only : open_for_read

        type(universe), intent(out) :: atoms
        character(len=*), intent(in) :: infile

        character(len=*), parameter :: err_read_poscar = "Error in read_poscar: "
        real(dp) :: scaling_const
        real(dp) :: simbox(3,3)
        integer, allocatable :: natoms(:)
        integer, parameter :: geo_unit = 55
        integer :: ios, nwords, ntypes, i, j
        integer :: nbeads = 1, pos1 = 0, pos2 = 0
        character(len=3), allocatable :: elements(:)
        character(len=max_string_length) :: buffer, c_system
        character(len=max_string_length) :: words(100)
        character(len=1), allocatable :: can_move(:,:,:)


        ntypes = simparams%nprojectiles + simparams%nlattices
        allocate(natoms(ntypes), elements(ntypes))

        call open_for_read(geo_unit, infile)

        read(geo_unit, '(A)', iostat=ios) buffer    ! first line is comment
        if (ios /= 0) stop err_read_poscar // "reading comment line"

        read(geo_unit, *, iostat=ios) scaling_const             ! second line is scaling constant
        if (ios /= 0) stop err_read_poscar // "reading scaling constant"

        !        print *, buffer
        !        print *, scaling_const

        ! read the cell vectors and multiply with scaling constant
        read(geo_unit, *, iostat=ios) simbox
        if (ios /= 0) stop err_read_poscar // "reading simulation box"

        ! read atom types
        read(geo_unit, '(A)', iostat=ios) buffer
        if (ios /= 0) stop err_read_poscar // "reading atom type"
        call split_string(buffer, elements, nwords)

        ! read number of atoms and number of beads
        read(geo_unit, '(A)', iostat=ios) buffer
        if (ios /= 0) stop err_read_poscar // "reading number of atoms and beads"
        pos1 = scan(string=buffer, set=":", back=.false.)
        if (pos1 > 0) then  ! Number of beads is stated in the file
            pos2 = scan(string=buffer, set=":", back=.true.)
            read(buffer(pos1+1:pos2-1), '(i10)', iostat=ios) nbeads
            if (ios /= 0 .or. nbeads < 1) stop err_read_poscar // "reading number of beads"
        end if
        read(buffer(pos2+1:), *, iostat=ios) natoms
        if (ios /= 0) stop err_read_poscar // "reading number of atoms"

        ! now we have everything to construct atom types
        atoms = new_atoms(nbeads, sum(natoms), ntypes)

        ! Direct or Cartesian coordinates, needed later
        read(geo_unit, '(A)', iostat=ios) c_system
        if (ios /= 0) stop err_read_poscar // "reading coordinate system type (cartesian/direct)"

        ! if coordinate section has 3 entries per line, no atoms are fixed
        ! if coordinate section has 6 entries per line, some might be fixed
        read(geo_unit, '(A)', iostat=ios) buffer
        if (ios /= 0) stop err_read_poscar // "reading start of coordinate section"
        backspace(geo_unit)
        call split_string(buffer, words, nwords)

        if (nwords == 3) then

            read(geo_unit, *, iostat=ios) atoms%r
            if (ios /= 0) stop err_read_poscar // "reading coordinates"

        else if (nwords == 6) then

            if (atoms%natoms > 0) then
                allocate (can_move(3, atoms%nbeads, atoms%natoms))
                read(geo_unit, *, iostat=ios) ((atoms%r(:,j,i), can_move(:,j,i), j=1,atoms%nbeads), &
                    i=1,atoms%natoms)
                if (ios /= 0) stop err_read_poscar // "reading atom coordinates and mobility flags"
                where (can_move == "F") atoms%is_fixed = .true.
                deallocate (can_move)
            else
                stop err_read_poscar // "reading projectile coordinates and mobility flags"
            end if

        else
        stop err_read_poscar // "cannot read coordinate section"

    end if

    ! Check if there is a velocity section
    read(geo_unit, '(A)', iostat=ios) buffer
    if (ios == iostat_end) then
        ! eof of file, don't read anything more
        atoms%v = 0.0_dp

    else if (ios == 0 .and. len_trim(buffer) == 0) then
        ! blank line after coordinate section and no eof yet -> read velocities
        read(geo_unit, *, iostat=ios) atoms%v
        if (ios /= 0) stop err_read_poscar // "cannot read initial atom velocities"

    else
        stop err_read_poscar // "cannot read initial atom velocities. Did you specify &
            all species in the *.inp file?"        
    end if

    ! Set the simulation box
    simbox = simbox * scaling_const
    atoms%simbox  = simbox
    atoms%isimbox = invert_matrix(simbox)

    ! Set the projectile identification
    if (simparams%nprojectiles > 0) atoms%is_proj(1:simparams%nprojectiles)  = .true.
    if (simparams%nlattices    > 0) atoms%is_proj(simparams%nprojectiles+1:) = .false.

    ! If system was given in cartesian coordinates, convert to direct
    call lower_case(c_system)
    if (trim(adjustl(c_system)) == "d" .or. &
        trim(adjustl(c_system)) == "direct") then

        atoms%is_cart = .false.

    else if (trim(adjustl(c_system)) == "c" .or. &
        trim(adjustl(c_system)) == "cart" .or. &
        trim(adjustl(c_system)) == "cartesian") then

        atoms%is_cart = .true.
        atoms%r = atoms%r * scaling_const
        call to_direct(atoms)

    else
        stop err_read_poscar // "reading coordinate system type (cartesian/direct)"

    end if
    call set_atomic_dofs(atoms)
    call set_atomic_indices(atoms, natoms)
    call set_atomic_masses(atoms, natoms)
    call set_atomic_names(atoms, elements)
    call set_prop_algos(atoms)
    call create_repetitions(atoms, simparams%rep)

end subroutine read_poscar



subroutine ensure_input_sanity()

    character(len=*), parameter :: err  = "input sanity check error: "
    character(len=*), parameter :: warn = "input sanity check warning: "

    ! Check the *.inp file
    !!! Perform MD
    if (simparams%run == "md") then
        if (simparams%start == default_int)  stop err // "start key must be present and larger than zero."
        if (simparams%ntrajs == default_int) stop err // "ntrajs key must be present and larger than zero."
        if (simparams%nsteps == default_int) stop err // "nsteps key must be present and larger than zero."
        if (simparams%step == default_real)  stop err // "nsteps key must be present and larger than zero."
        if (simparams%nlattices == default_int .and. simparams%nprojectiles == default_int) stop err // "either lattice and/or projectile key must be present."
        if (simparams%nprojectiles == default_int) simparams%nprojectiles = 0
        if (simparams%nlattices == default_int) simparams%nlattices = 0

        if (simparams%nprojectiles > 0 .and. .not. allocated(simparams%name_p)) stop err // "projectile names not set."
        if (simparams%nprojectiles > 0 .and. .not. allocated(simparams%mass_p)) stop err // "projectile masses not set."
        if (simparams%nprojectiles > 0 .and. .not. allocated(simparams%mass_p)) stop err // "projectile masses not set."
        if (simparams%nprojectiles > 0 .and. .not. allocated(simparams%md_algo_p)) stop err // "projectile propagation algorithm not set."

        if (simparams%nlattices > 0 .and. .not. allocated(simparams%mass_l)) stop err // "lattice masses not set."
        if (simparams%nlattices > 0 .and. .not. allocated(simparams%name_l)) stop err // "lattice names not set."
        if (simparams%nlattices > 0 .and. .not. allocated(simparams%mass_l)) stop err // "lattice masses not set."
        if (simparams%nlattices > 0 .and. .not. allocated(simparams%md_algo_l)) stop err // "lattice propagation algorithm not set."

        if (.not. allocated(simparams%output_type) .or. .not. allocated(simparams%output_interval)) &
            stop err // "output options not set."
        if (any(simparams%output_interval > simparams%nsteps)) print *, warn // "one or more outputs will not show, since the output &
	    interval is larger than the total number of simulation steps."        
        

        ! If one of them is set, the others must be set as well.
        if ( (allocated(simparams%einc) .neqv. allocated(simparams%inclination)) .or. &
            (allocated(simparams%einc) .neqv. allocated(simparams%azimuth)) ) &
            stop err // "incidence conditions ambiguous."



            ! TODO: to be completed



    !!! Perform fit
    else if (simparams%run == "fit") then


    end if



end subroutine ensure_input_sanity



subroutine ensure_geometry_sanity(atoms)

    type(universe), intent(in) :: atoms
    character(len=*), parameter :: err = "geometry sanity check error: "


    !if (atoms%nbeads > 1 .and. simparams%Tsurf == default_real) stop err // "RPMD needs system temperature."
    if (simparams%force_beads < 1) stop err // "Cannot force zero or less beads."
    if (simparams%force_beads /= default_int .and. atoms%nbeads /= 1) stop err // "Cannot force any beads when system already contains multiple beads."

end subroutine ensure_geometry_sanity




subroutine post_process(this)

    use run_config, only : simparams

    type(universe), intent(inout) :: this

    type(universe) :: other
    integer :: b

    ! possibly force the system to have multiple beads
    if (simparams%force_beads /= default_int) then

        other = new_atoms(simparams%force_beads, this%natoms, this%ntypes)

        do b = 1, simparams%force_beads

            other%r(:,b,:) = this%r(:,1,:)
            other%v(:,b,:) = this%v(:,1,:)
            other%a(:,b,:) = this%a(:,1,:)
            other%f(:,b,:) = this%f(:,1,:)
            other%epot(b)  = this%epot(1)
            other%is_fixed(:,b,:) = this%is_fixed(:,1,:)

        end do

        other%idx     = this%idx
        other%m       = this%m
        other%name    = this%name
        other%pes     = this%pes
        other%algo    = this%algo
        other%is_proj = this%is_proj

        other%dof     = simparams%force_beads * this%dof
        other%simbox  = this%simbox
        other%isimbox = this%isimbox
        other%is_cart = this%is_cart

        this = other

    end if

end subroutine post_process





subroutine set_atomic_indices(atoms, natoms)

    ! Each atom has a unique index 1, 2, ..., n, where n is the number of species in
    !  the system.
    !  The order is given in the *.inp file. It has to match the order in the geometry file.

    type(universe), intent(inout) :: atoms
    integer, intent(in) :: natoms(:)

    integer :: i, j, counter

    counter = 1
    !print *, natoms
    do i = 1, size(natoms)
        do j = 1, natoms(i)
            atoms%idx(counter) = i
            counter = counter + 1
        end do
    end do
    if (any(atoms%idx == default_int))  stop "Error in set_atomic_indices: not all atom indices set."
    if (size(natoms) > size(atoms%idx)) stop "Error in set_atomic_indices: more atoms than indices."

end subroutine set_atomic_indices


subroutine set_atomic_names(atoms, elements)

    type(universe), intent(inout) :: atoms
    character(len=3), intent(in) :: elements(:)
    character(len=*), parameter :: err = "Error in set_atomic_names: "

    ! Check if element names from *.inp file and poscar file are in agreement

    if (allocated(simparams%name_p) .and. allocated(simparams%name_l)) then
        if (any(elements /= [simparams%name_p, simparams%name_l])) stop err // "atomic names in *.inp and poscar files differ"

    else if (allocated(simparams%name_p)) then
        if (any(elements /= simparams%name_p)) stop err // "atomic names in *.inp and poscar files differ"

    else if (allocated(simparams%name_l)) then
        if (any(elements /= simparams%name_l)) stop err // "atomic names in *.inp and poscar files differ"

    else
        stop err // "neither projectile nor slab exist"
    end if

    atoms%name = elements

end subroutine set_atomic_names


subroutine set_atomic_masses(atoms, natoms)

    type(universe), intent(inout) :: atoms
    integer, intent(in) :: natoms(:)

    integer :: i
    real(dp) :: masses(size(natoms))
    character(len=*), parameter :: err = "Error in set_atomic_masses: "

    if (allocated(simparams%mass_p) .and. allocated(simparams%mass_l)) then
        masses = [simparams%mass_p, simparams%mass_l]
    else if (allocated(simparams%mass_p)) then
        masses = simparams%mass_p
    else if (allocated(simparams%mass_l)) then
        masses = simparams%mass_l
    else
        stop err // "neither projectile nor slab exist"
    end if

    do i = 1, size(natoms)
        atoms%m(sum(natoms(:i-1))+1:sum(natoms(:i))) = masses(i)
    end do

    ! To program units
    atoms%m = atoms%m * amu2mass

end subroutine set_atomic_masses




subroutine set_prop_algos(atoms)

    type(universe), intent(inout) :: atoms

    character(len=*), parameter :: err = "Error in set_prop_algos: "

    ! Either (md_algo_p or md_algo_l) or both are allocated

    if (allocated(simparams%md_algo_p) .and. allocated(simparams%md_algo_l)) then
        atoms%algo = [simparams%md_algo_p, simparams%md_algo_l]
        if (size(atoms%algo) /= size([simparams%md_algo_p, simparams%md_algo_l])) stop err // "more algorithms than species"

    else if (allocated(simparams%md_algo_p)) then
        atoms%algo = simparams%md_algo_p
        if (size(atoms%algo) /= size(simparams%md_algo_p)) stop err // "more algorithms than species"

    else if (allocated(simparams%md_algo_l)) then
        atoms%algo = simparams%md_algo_l
        if (size(atoms%algo) /= size(simparams%md_algo_l)) stop err // "more algorithms than species"

    else
        stop err // "neither projectile nor slab exist"
    end if

end subroutine set_prop_algos


subroutine set_atomic_dofs(atoms)

    type(universe), intent(inout) :: atoms

    atoms%dof = 3*atoms%natoms*atoms%nbeads - count(atoms%is_fixed)

end subroutine set_atomic_dofs

!
!subroutine read_conf(nr_at_layer, nlnofix, nlno, n_p, n_l, n_p0, &
!                   slab, teil, start_l, c_matrix)
!    !
!    ! Purpose:
!    !           Read in geometries from POSCAR file and multiply them
!    !
!type(atoms), intent(out) :: teil, slab
!
!character(len=80) :: buffer
!real(8) :: cellscale
!real(8), dimension(3,3) :: c_matrix, d_matrix
!real(8), dimension(:,:), allocatable :: start_l, start_p
!real(8), dimension(:,:), allocatable :: pos_l, vel_l, pos_p
!integer :: n_l0, n_l, n_p, n_p0
!integer :: nlnofix, nlno
!integer :: ios
!character(len= 1) :: coord_sys
!integer, dimension(:), allocatable :: nr_at_layer
!
!! Read in the POSCAR-file.
!call open_for_read(38, confname_file)
!
!read(38,*) buffer
!read(38,*) cellscale
!read(38,*) c_matrix
!read(38,'(A)') buffer
!read(buffer, *, iostat=ios) n_l0, n_p
!read(38,*) coord_sys
!a_lat = c_matrix(1,1)/celldim(1)*sqrt2
!
!!------------------------------------------------------------------------------
!!           Construct simulation cell matrix and its inverse
!!------------------------------------------------------------------------------
!
!d_matrix = 0.0d0
!d_matrix(1,1) = c_matrix(2,2) / (c_matrix(1,1)*c_matrix(2,2) - c_matrix(1,2)*c_matrix(2,1))
!d_matrix(2,2) = c_matrix(1,1) / (c_matrix(1,1)*c_matrix(2,2) - c_matrix(1,2)*c_matrix(2,1))
!d_matrix(3,3) = 1.0d0/c_matrix(3,3)
!d_matrix(1,2) = c_matrix(1,2) / (c_matrix(1,2)*c_matrix(2,1) - c_matrix(1,1)*c_matrix(2,2))
!d_matrix(2,1) = c_matrix(2,1) / (c_matrix(1,2)*c_matrix(2,1) - c_matrix(1,1)*c_matrix(2,2))
!
!! Size of cell-matrix changes if slab repeated with rep
!cell_mat = 0.0d0
!cell_mat(1:2,1) = c_matrix(1:2,1)*(2*rep(1) + 1)
!cell_mat(1:2,2) = c_matrix(1:2,2)*(2*rep(2) + 1)
!cell_mat(3  ,3) = c_matrix(3  ,3)
!
!cell_imat = 0.0d0
!cell_imat(1,1:2) = d_matrix(1,1:2)/(2*rep(1) + 1)
!cell_imat(2,1:2) = d_matrix(2,1:2)/(2*rep(2) + 1)
!cell_imat(3,  3) = d_matrix(3,  3)
!
!!------------------------------------------------------------------------------
!!                           Read in coordinates
!!------------------------------------------------------------------------------
!allocate(start_l(3,n_l0))
!read(38,*) start_l
!! Read in coordinates of particle if present
!if (n_p > 0) then
!    allocate(start_p(3,n_p))
!    read(38,*) start_p
!    if (confname == 'fit') &
!       start_p(:,1:n_p) = start_p(:,1:n_p) - spread(start_p(:,1),2,n_p)
!end if
!
!! Transform the read in coordinates into direct if they are cartesians:
!if (coord_sys == 'C' .or. coord_sys == 'c') then
!    start_l=matmul(d_matrix,start_l)
!    if (n_p > 0) start_p=matmul(d_matrix,start_p)
!end if
!! Call subroutine prepare slab
!! With repetition of cell.
!call prep_slab(n_l0, start_l, c_matrix, nr_at_layer, nlnofix, nlno,&
!               n_l, pos_l, vel_l)
!
!! Create slab objects
!slab = atoms(n_l)
!! Assign slab positions and velocities
!slab%r = pos_l
!slab%v = vel_l
!
!! Assign the number of non-fixed atom
!slab%nofix = nlnofix
!!        call neighbour_list(pos_l, n_l, neigh)
!!        allocate(neigh_l(n_l,nneighs(pes_nigh)))
!!        neigh_l = neigh
!!        deallocate(neigh)
!
!! Check if particle present.
!if (md_algo_p > 0 .and. n_p == 0 .and. n_p0 == 0) then
!    md_algo_p = 0
!    print *, "Warning: Number of projectiles both in POSCAR and input file is 0."
!    print *, "         Calculations will continue without projectile."
!end if
!
!if (md_algo_p > 0) then
!    call prep_teil(n_p, n_p0, start_p, c_matrix, pos_p)
!end if
!n_p = n_p0
!! Create projectile objects
!if (n_p > 0) then
!    teil = atoms(n_p)
!    teil%r = pos_p(:,1:n_p)
!!            call neighbour_list(pos_p,n_p, neigh)
!!            allocate(neigh_p(n_p,nneighs(pes_nigh)))
!    deallocate(start_p, pos_p)
!else
!    teil = atoms(1)
!    teil%n_atoms = 0
!    teil%nofix = 0
!end if
!
!close(38)
!
!deallocate(pos_l, vel_l)
!
!end subroutine read_conf
!
!subroutine prep_slab(n_l0, start_l, c_matrix, nr_at_layer, &
!                     nlnofix, nlno, n_l, pos_l, vel_l)
!    !
!    ! Purpose:
!    !------------------------------------------------------------------------------
!    !                    Replicate input cell
!    !                    ====================
!    !------------------------------------------------------------------------------
!    ! The final array size for the lattice and the particle is determined by the
!    ! repetitions of the cell image.
!    ! 1. The new array is build by repetitions of the old one if rep is specified
!    !    in the input file. For example:
!    !       for example,    rep =(/0,0/) yields  1x1-1=  0       new images
!    !                       rep =(/0,1/)         1x3-1=  2
!    !                       rep =(/1,1/)         3x3-1=  8
!    !                       rep =(/1,2/)         3x5-1= 14
!    !                       rep =(/2,2/)         5x5-1= 24
!    !
!    ! 2. 2*(rep(1)+1)x(2*rep(2)+1) - 1 new images are formed.
!    ! 3. Every layer is repeated independently. That means layer can hold a
!    !    different number of atoms.
!    ! 4. If rep is not given, then no repetition takes place.
!    ! 5. Set velocities.
!
!integer, dimension(:), allocatable :: nr_at_layer
!integer :: i, j, k, l, itemp, itemp2
!integer :: n_l, n_l0
!integer :: nlnofix, nlno
!real(8), dimension(:,:), allocatable :: pos_l, d_l, start_l, vel_l
!real(8), dimension(3,3) :: c_matrix
!real(8) :: v_pdof
!! calculate new cellsize
!if (allocated(nr_at_layer) .eqv. .true.) then
!    if (rep(1) > 0 .or. rep(2) > 0 ) then
!        itemp = 0
!        do i = 1,celldim(3)
!            itemp = itemp + nr_at_layer(i)
!        end do
!        n_l = itemp * (2*rep(1)+1)*(2*rep(2)+1)
!    end if
!else if (allocated(nr_at_layer) .eqv. .false. ) then
!    if (rep(1) > 0 .or. rep(2) > 0 ) then
!        itemp=celldim(1)*celldim(2)
!        n_l=itemp*celldim(3)*(2*rep(1)+1)*(2*rep(2)+1)
!    end if
!end if
!! Check if Repetition necessary at all
!if (rep(1) == 0 .and. rep(2) == 0 ) then
!    n_l = n_l0
!end if
!
!allocate(d_l(3,n_l))
!d_l = 0.0d0
!! Replication. Every layer is repeated independently.
!i = 1
!itemp = 0
!if (rep(1) > 0 .or. rep(2) > 0 ) then
!if (allocated(nr_at_layer) == .true.) then
!    itemp2 = 1
!    do l = 1, celldim(3)
!        do j = -rep(1), rep(1)
!            do k = -rep(2), rep(2)
!                itemp = nr_at_layer(l)
!                d_l(1,i:i+itemp-1) = start_l(1,itemp2:itemp2+nr_at_layer(l)-1)+j
!                d_l(2,i:i+itemp-1) = start_l(2,itemp2:itemp2+nr_at_layer(l)-1)+k
!                d_l(3,i:i+itemp-1) = start_l(3,itemp2:itemp2+nr_at_layer(l)-1)
!                i = i + itemp
!            end do
!        end do
!        itemp2 = itemp2+nr_at_layer(l)
!    end do
!else if (allocated(nr_at_layer) == .false.) then
!
!    itemp=celldim(1)*celldim(2)
!    i = 1
!    do l = 1, celldim(3)
!        do j =-rep(1), rep(1)
!            do k=-rep(2), rep(2)
!                d_l(1,i:i+itemp-1) = start_l(1,(l-1)*itemp+1:l*itemp)+j
!                d_l(2,i:i+itemp-1) = start_l(2,(l-1)*itemp+1:l*itemp)+k
!                d_l(3,i:i+itemp-1) = start_l(3,(l-1)*itemp+1:l*itemp)
!                i = i+itemp
!            end do
!        end do
!    end do
!end if
!else
!    d_l = start_l
!end if
!
!! Transform back into cartesian coordinates and set velocities
!allocate(pos_l(3,n_l), vel_l(3,n_l))
!pos_l = matmul(c_matrix,d_l)
!
!! Exclude fixed atoms
!if (nlno < 0 .and. allocated(nr_at_layer) .eqv. .false.) then
!    nlnofix = n_l/celldim(3)*(celldim(3)+nlno)    ! lowest layer fixed
!else if (nlno < 0 .and. allocated(nr_at_layer) .eqv. .true.) then
!    ! calculate number of remaining atoms when lowest layer fixed
!    nlnofix = n_l + nr_at_layer(celldim(3))*(2*rep(1)+1)*(2*rep(2)+1)*nlno
!else
!    nlnofix = n_l - nlno                       ! whatever number nlno of atoms fixed
!end if
!
!
!! Sample velocities of lattice atoms from thermal distribution
!! assuming the minimum energy configuration
!v_pdof = sqrt(2.0d0*kB*Tsurf/mass_l)
!
!vel_l = 0.0d0
!do i=1,nlnofix ! set velocities for atoms that aren't fixed
!    vel_l(1,i) = normal(0.0d0,v_pdof)
!    vel_l(2,i) = normal(0.0d0,v_pdof)
!    vel_l(3,i) = normal(0.0d0,v_pdof)
!enddo
!! Set c.-of-m. velocity to zero
!vel_l(1,1:nlnofix) = vel_l(1,1:nlnofix) - sum(vel_l(1,1:nlnofix))/nlnofix
!vel_l(2,1:nlnofix) = vel_l(2,1:nlnofix) - sum(vel_l(2,1:nlnofix))/nlnofix
!vel_l(3,1:nlnofix) = vel_l(3,1:nlnofix) - sum(vel_l(3,1:nlnofix))/nlnofix
!
!deallocate(d_l)
!
!end subroutine prep_slab
!
!subroutine prep_teil(n_p, n_p0, start_p, c_matrix, pos_p)
!    !
!    ! Purpose:
!    !           Prepares and redublicates particle
!    !           Particle does not have a crystal structure, therefore no layers
!integer :: n_p, n_p0
!real(8), dimension(:,:), allocatable :: start_p, pos_p, d_p
!real(8), dimension(3,3) :: c_matrix
!integer :: itemp, i, j, l, r, s
!
!itemp = n_p*(2*rep(1)+1)*(2*rep(2)+1)
!if (n_p0 == 0) n_p0 = itemp
!
!allocate(d_p(3,n_p0))
!d_p = 0.0d0
!if (itemp < n_p0) then
!    print *, "Warning: Number of projectiles larger than can be produced"
!    print *, "         from repititions of POSCAR file."
!    print *, "         All projectile positions set to zero."
!else
!    i=1
!    j=1
!    l=1
!
!    do r =-rep(1), rep(1)
!    do s =-rep(2), rep(2)
!    do l =   1, n_p
!
!        if (j > n_p0) exit
!        d_p(1,j) = start_p(1,l)+r
!        d_p(2,j) = start_p(2,l)+s
!        d_p(3,j) = start_p(3,l)
!
!        j=j+1
!    end do
!    end do
!    end do
!end if
!
!allocate(pos_p(3,n_p0))
!pos_p = matmul(c_matrix,d_p)
!
!deallocate(d_p)
!
!end subroutine prep_teil
!
!subroutine read_mxt(nspec, teil, slab, n_p0)
!    !
!    ! Purpose:
!    !           Read in configurations from the mxt-file
!
!    type(atoms), intent(out) :: teil, slab
!    integer :: nspec, n_p0
!    integer :: n_l, nlnofix, n_p, npnofix
!    integer :: itemp
!    real(8) :: rtemp
!    character(len=8) :: str
!
!    if (confname == 'geo') then
!        write(str,'(I8.8)') conf_nr
!        open(unit=38, file=trim(confname_file)//'/mxt_conf'//str//'.bin', &
!            form='unformatted' )
!    else
!        open(unit=38, file=trim(confname_file)//'/mxt_conf00000001.bin', &
!        form='unformatted' )
!    end if
!
!    read(38) itemp
!    if (step < 0.0d0) then
!        read(38) step
!    else
!        read(38) rtemp
!    end if
!    read(38) rtemp, rtemp
!    read(38) rtemp
!    if (.not. Tsurf_key) then
!        Tsurf = rtemp
!    else if ( rtemp .ne. Tsurf) then
!        print '(A,f5.1,A)', 'Warning: The surface temperature from input (',Tsurf,' K) '
!        print '(A,f5.1,A)', '         and configuration (',rtemp, ' K) files do not match.'
!        print *,            '         Temperature set to that from input file.'
!    end if
!    read(38) nspec
!    ! potential and neighbouring
!    read(38) pes_name!, pes_nigh
!
!    read(38) name_l, n_l, nlnofix, mass_l, npars_l, key_l
!
!    if (.not. allocated(pars_l)) allocate(pars_l(npars_l))
!
!    read(38) pars_l, itemp
!    if (.not. md_algo_l_key) md_algo_l = itemp
!
!    read(38) a_lat        ! lattice constant makes us happy
!    read(38) cell_mat     ! Cell matrix
!    read(38) cell_imat    ! inverse cell matrix
!
!    slab = atoms(n_l)
!    slab%nofix = nlnofix
!
!    read(38) slab%r, slab%v, slab%a, slab%dens
!
!    if (.not. md_algo_p_key) then
!        if (nspec > 1) then
!            ! Projectile hasn't been defined in input file.
!            ! Projectile defined from configuration file
!            read(38) name_p, n_p, npnofix, mass_p, npars_p, key_p
!            if (.not. allocated(pars_p)) allocate(pars_p(npars_p))
!            read(38) pars_p, itemp
!            if (.not. md_algo_p_key) md_algo_p = itemp
!            teil = atoms(n_p)
!            read(38) teil%r, teil%v, teil%a, teil%dens
!            teil%nofix = npnofix
!        else
!            ! Projectile hasn't been defined in input file.
!            ! It does not exist in configuration file.
!            n_p = 1
!            teil = atoms(n_p)
!            teil%n_atoms = 0
!            teil%nofix = 0
!        end if
!    else
!        ! Projectile has been defined in input file.
!        ! Definition in input file always outrules other options.
!        teil = atoms(n_p0)
!    end if
!    close(38)
!
!end subroutine read_mxt
!
!subroutine read_fit(fracaimd, n_p, n_p0, n_l, n_l0, teil, slab, &
!                    start_l, c_matrix)
!
!    !
!    ! Purpose:
!    !           Read in data for fit.
!    !
!type(atoms) :: slab, teil
!
!integer :: n_p, n_l, n_p0, n_l0
!logical :: exists
!integer, dimension(2) :: fracaimd, npts, ea_fraction
!real(8), dimension(3,3):: c_matrix
!
!real(8) :: rtemp
!integer :: i, j, ios, itemp, q, l, k, r, s
!real(8), dimension(3) :: empty3
!real(8), dimension(:), allocatable :: y_eq, E_dft1, E_dft2
!real(8), dimension(:,:), allocatable :: dfix, start_l
!real(8), dimension(:,:,:),allocatable :: x_eq, aimd_l, aimd_l1, aimd_p
!real(8), dimension(:,:,:),allocatable :: d_l3, d_p3, aimd_p1
!character(len=80) :: str, buffer
!!
!real(8), dimension(3) :: r3temp
!real(8) :: minr, maxr ! distance and counter for the largest/smallest distance
!real(8), dimension(:,:), allocatable :: cpaimd
!
!! Check if Equilibrium-data exists.
!inquire(directory=trim(fit_dir),exist=exists)
!if (.not. exists) stop 'The folder with fit data does not exist.'
!npts = 0
!! read in energies for Equilibrium positions
!! and select only those that correspond to the chosen sites
!if (fracaimd(1) > 0 ) then
!    call open_for_read(39,trim(fit_dir)//fit_eq)
!    i=1
!    do
!        read(39,*,iostat=ios) itemp, rtemp, rtemp, rtemp, rtemp
!        if(ios <0) exit
!        do j = 1, nsites
!            if ((rtemp - evasp)<=dftlu(2) .and. (rtemp - evasp)>=dftlu(1)&
!            .and. itemp == ssites(j)) i=i+1
!        end do
!    end do
!    npts(1) = i - 1
!
!    if (npts(1) < fracaimd(1)) then
!        fracaimd(1) = npts(1)
!        print *, 'Your desired number of equilibrium configurations'
!        print *, 'is too large. Setting the number to', fracaimd(1)
!    end if
!    rewind(39)
!end if
!
!allocate(y_eq(npts(1)+1),x_eq(npts(1)+1,3,n_l+n_p))
!! The first element in the position arrays contains the equilibrium coordinates
!! of the slab-atoms and the particle 6.0 A above the slab.
!! This element is used to calculate the reference energy and always needs to be present.
!do j = 1, n_p
!    x_eq(1,:,j) = teil%r(:,j) + (/0.0d0,0.0d0,6.0d0/)
!end do
!x_eq(1,:,n_p+1:n_p+n_l) = slab%r
!y_eq(1) = evasp
!
!! If not only data from AIMD fitted
!if (fracaimd(1) > 0 ) then
!    i=2
!    do
!        read(39,*,iostat=ios) itemp, empty3(1), empty3(2), empty3(3), rtemp
!        if(ios <0) exit
!        do k = 1, nsites
!        if ((rtemp - evasp)<=dftlu(2) .and. (rtemp - evasp)>=dftlu(1)&
!            .and. itemp == ssites(k)) then
!            y_eq(i) = rtemp
!            do j = 1, n_p
!                x_eq(i,:,j) = teil%r(:,j) + empty3
!            end do
!            x_eq(i,:,n_p+1:n_p+n_l) = slab%r
!            i = i + 1
!        end if
!        end do
!    end do
!    close(39)
!    ea_fraction(1) = npts(1)/fracaimd(1)
!end if
!
!!------------------------------------------------------------------------------
!!                           READ IN AIMD GEOMETRIES
!!------------------------------------------------------------------------------
!if (fracaimd(2) > 0 ) then
!
!    write(str,'(A,I3.3,A)') 'traj',n_confs,'/'  ! name of directory with AIMD-data
!    ! read in energies and timesteps of the AIMD trajectory
!    ! The geometries and energies are in different files (which is why we are
!    ! opening two here.)
!    call open_for_read(17, trim(fit_dir)//trim(str)//aimd_e) ! analyse.dat for energies
!    call open_for_read(18, trim(fit_dir)//trim(str)//aimd_pos) !XDATCAR-file for positions
!
!    ! The coefficients of the transformation matrix.
!    read(18,'(A,/////)') buffer
!    n_p0 = 0
!    read(18,'(A)') buffer
!    read(buffer, *) n_l0, n_p0
!
!    ! According to the energy file, we define number of configurations
!    i = 0
!    do
!        read(17,*,iostat=ios)
!        if(ios <0) exit
!        i = i + 1
!    end do
!    rewind(17)
!    npts(2) = i - 3
!
!    allocate(E_dft2(npts(2)))
!    allocate(aimd_l1(npts(2),3,n_l0),cpaimd(3,n_l0+n_p0))
!    allocate(dfix(3,n_l0))
!    allocate(aimd_p1(npts(2),3,n_p0))
!
!! Read in positions and energies together
!    read(17,'(A)') buffer
!!    read(17,*) rtemp, rtemp, E_dft2(1)
!!    read(18,*) aimd_l1(1,:,:)
!!    if (n_p0 > 0) read(18,*) aimd_p1(1,:,:)
!
!    j=1
!    do i=1,npts(2)
!        read(17,*) rtemp, rtemp, E_dft2(j)
!        read(18,*) aimd_l1(j,:,:)
!
!        if (n_p0 > 0) then
!            read(18,*) aimd_p1(j,:,:)
!            cpaimd(:,1:n_l0) = aimd_l1(j,:,:)
!            cpaimd(:,n_l0+1:n_l0+n_p0) = aimd_p1(j,:,:)
!            ! Calculate smallest distance to Au-atoms
!            r=1000.d0
!            minr=1000.d0
!            maxr=-1000.d0
!            cpaimd = matmul(cell_mat,cpaimd(:,:))
!!            print *, cpaimd(3,1)- cpaimd(3,17)
!            do k=1,n_p0
!            do l=1,n_l0
!                r3temp= cpaimd(:,l)-cpaimd(:,n_l0+k)
!                r3temp = matmul(cell_imat, r3temp)   ! transform to direct coordinates
!                r3temp(1) = r3temp(1) - Anint(r3temp(1))! imaging
!                r3temp(2) = r3temp(2) - Anint(r3temp(2))
!                r3temp(3) = r3temp(3) - Anint(r3temp(3))
!                r3temp    = matmul(cell_mat, r3temp)    ! back to cartesian coordinates
!
!                r =  sqrt(sum(r3temp*r3temp))               ! distance
!                if (r < minr) minr = r
!                if(l==1) maxr = -r3temp(3)
!            end do
!            end do
!        end if
!
!        if((E_dft2(j) - evasp)<=aimdlu(2) .and. (E_dft2(j) - evasp)>=aimdlu(1)&
!            .and. maxr < distmima(2) .and. minr > distmima(1)) then
!            j = j + 1
!        end if
!
!    end do
!    npts(2) = j - 1
!
!! The procedure above might have excluded some of the points from the trajectories;
!! Therefore, new arrays now have to be allocated that have the correct length and
!! contain no empty spots
!allocate(E_dft1(npts(2)),aimd_l(npts(2),3,n_l0),aimd_p(npts(2),3,n_p0))
!
!do i=1,npts(2)
!        E_dft1(i)=E_dft2(i)
!        aimd_l(i,:,:) = aimd_l1(i,:,:)
!        aimd_p(i,:,:) = aimd_p1(i,:,:)
!end do
!deallocate(E_dft2,aimd_l1,aimd_p1)
!
!
!    if (npts(2) < fracaimd(2)) then
!        fracaimd(2) = npts(2)
!        print *, 'Your desired number of non-equilibrium configurations'
!        print *, 'is too large. Setting the number to', fracaimd(2)
!    end if
!    close(17)
!    close(18)
!    ! Place AIMD-atoms close to corresponding equilibrium positions
!     do i=1,npts(2)
!        dfix = aimd_l(i,:,:) - start_l
!        do j=1,n_l0
!               aimd_l(i,1,j)=aimd_l(i,1,j) - ANINT(dfix(1,j))
!               aimd_l(i,2,j)=aimd_l(i,2,j) - ANINT(dfix(2,j))
!               aimd_l(i,3,j)=aimd_l(i,3,j) - ANINT(dfix(3,j))
!        end do
!    end do
!
!    itemp=celldim(1)*celldim(2)
!    allocate(d_l3(npts(2),3,n_l))
!
!    ! Replication. We assume that the AIMD-trajectories always have the same
!    ! Number of atoms per layer.
!    do q = 1,npts(2)
!        i = 1
!        do l = 1, celldim(3)
!        do j =-rep(1), rep(1)
!        do k =-rep(2), rep(2)
!            d_l3(q,1,i:i+itemp-1) = aimd_l(q,1,(l-1)*itemp+1:l*itemp)+j
!            d_l3(q,2,i:i+itemp-1) = aimd_l(q,2,(l-1)*itemp+1:l*itemp)+k
!            d_l3(q,3,i:i+itemp-1) = aimd_l(q,3,(l-1)*itemp+1:l*itemp)
!            i = i+itemp
!        end do
!        end do
!        end do
!        d_l3(q,:,:)= matmul(c_matrix,d_l3(q,:,:))
!    end do
!
!    ! Projectile initialization and replication
!    if (n_p0 > 0) then    ! projectile existence justified
!
!        allocate(d_p3(npts(2),3,n_p))
!        j=1
!        do q = 1,npts(2)
!            j = 1
!            do r =-rep(1), rep(1)
!            do s =-rep(2), rep(2)
!            do l =   1, n_p0
!                d_p3(q,1,j) = aimd_p(q,1,l)+r
!                d_p3(q,2,j) = aimd_p(q,2,l)+s
!                d_p3(q,3,j) = aimd_p(q,3,l)
!                j=j+1
!            end do
!            end do
!            end do
!            d_p3(q,:,:)= matmul(c_matrix,d_p3(q,:,:))
!        end do
!
!    end if
!    ea_fraction(2) = npts(2)/fracaimd(2)
!
!end if
!
!
!    !------------------------------------------------------------------------------
!    !                   HOW MUCH AIMD CONTRIBUTION DO WE WANT?
!    !                   ======================================
!    !------------------------------------------------------------------------------
!
!        allocate(x_all(fracaimd(1)+fracaimd(2)+1,3,n_p+n_l))
!        allocate(y_all(fracaimd(1)+fracaimd(2)+1))
!        y_all = 0.0d0
!        x_all = 0.0d0
!
!        x_all(1,:,:) = x_eq(1,:,:)
!        x_all(2:fracaimd(1)+1,:,:) = x_eq(2:npts(1)+1:ea_fraction(1),:,:)
!        x_all(fracaimd(1)+2:,:,1:n_p) = d_p3(1:npts(2):ea_fraction(2),:,:)
!        x_all(fracaimd(1)+2:,:,n_p+1:) = d_l3(1:npts(2):ea_fraction(2),:,:)
!
!        y_all(1)  = y_eq(1)
!        y_all(2:fracaimd(1)+1)  = y_eq(2:npts(1)+1:ea_fraction(1))
!        y_all(fracaimd(1)+2:) = E_dft1(1:npts(2):ea_fraction(2))
!        y_all = y_all - evasp
!
!        if (fracaimd(2) > 0 )&
!            deallocate(d_p3, d_l3, aimd_l, aimd_p, E_dft1, dfix)
!        if (fracaimd(1) > 0 ) deallocate(y_eq, x_eq)
!
!        if (teil%n_atoms > 0) then
!            teil%r(3,:) = 6.0d0
!            np_atoms = n_p
!        end if
!        nl_atoms = n_l
!
!end subroutine read_fit
!
!subroutine traj_init(slab, teil, Eref)
!!
!! Purpose:
!!           Initialise the entire system:
!!               1. Geometry
!!               2. Interaction Potential
!!               3. Velocities
!!
!    implicit none
!
!    type(atoms), intent(inout) :: slab, teil   ! hold r, v and f for atoms in the box
!
!    integer :: traj_no, nspec, i
!    real(8) :: dum, Eref
!    integer :: ymm, nof
!    character(len=80) :: ddd
!    character(len=8) :: str
!
!
!
!!------------------------------------------------------------------------------
!!                       READ IN CONFIGURATION
!!                       =====================
!!------------------------------------------------------------------------------
!
!
!    ! Get random integer to open random trajectory for configurations
!    if (confname == 'mxt') then
!        write(str,'(I8.8)') int(ran1()*n_confs)+1
!    elseif (confname == 'geo') then
!        write(str,'(I8.8)') conf_nr
!    end if
!
!    open(unit=38, file=trim(confname_file)//'/mxt_conf'//str//'.bin', form='unformatted' )
!
!    read(38) traj_no
!    read(38) dum
!    read(38) dum, Eref
!    read(38) dum
!    read(38) nspec
!    ! potential and neighbouring
!    read(38) ddd!, ymm
!
!    read(38) ddd, ymm, nof, dum, ymm, ddd
!    read(38) (dum, i=1,ymm), ymm
!
!    read(38) dum        ! lattice constant makes us happy
!    read(38) (dum, i=1,9)
!    read(38) (dum, i=1,9)
!
!    read(38) slab%r, slab%v, slab%a, slab%dens
!
!    if (.not.md_algo_p_key .and. nspec > 1) then
!        read(38) ddd, ymm, ymm, dum, ymm, ddd
!        read(38) (dum, i=1,ymm), ymm
!        read(38) teil%r, teil%v, teil%a, teil%dens
!    end if
!
!    close(38)
!
!end subroutine traj_init
!
!subroutine particle_init(s)
!    !
!    ! Purpose:
!    !           Assign random positions to the particle atom
!    !           Furthermore set velocities
!    !
!
!    type(atoms) :: s
!    integer :: i, n_p
!    real(8) :: vinc
!    real(8), dimension(2) :: cc1, cc2
!    real(8) :: v_pdof = 0.0d0
!
!    cc1 = (/a_lat*isqrt2,0.0d0/)
!    cc2 = 0.5d0*cc1(1)*(/-1.0d0,sqrt3/)
!
!    if (confname == 'mxt' .or. confname == 'geo') then !Why do we have an if here, anyway?
!
!        select case(pip_sign)
!
!        case(0) ! Random positions
!
!            do i=1,s%n_atoms
!                s%r(1:2,i) = matmul((/ran1(), ran1()/),cell_mat(1:2,1:2))
!                s%r(3,i)   = height
!            end do
!
!        case(1) ! Put atom in designated positions
!
!            call lower_case(key_p_pos)
!            select case (key_p_pos)
!            case('top')
!                s%r(1:2,1) = (/0.,0./)
!            case('fcc')
!                s%r(1:2,1) = (cc1 + 2.0d0*cc2)/3.0d0
!            case('hcp')
!                s%r(1:2,1) = (2.0d0*cc1 + cc2)/3.0d0
!            case('bri')
!                s%r(1:2,1) = (cc1 + cc2)*0.5d0
!            end select
!
!            s%r(1,1) = s%r(1,1) - height*tan(inclination)*cos(azimuth)
!            s%r(2,1) = s%r(2,1) - height*tan(inclination)*sin(azimuth)
!            s%r(3,1) = height
!
!        case(2) ! Read in projectile positions from other file.
!
!            call open_for_read(44,key_p_pos)
!            read(44,*) n_p
!
!            if (n_p < s%n_atoms) then
!
!                print *, 'Error: Number of atoms in projectile configuration file'
!                print *, '       is smaller than that in input file.'
!                stop ' subroutine: particle_init()'
!
!            else if (n_p > s%n_atoms) then
!
!                print *, 'Warning: Number of atoms in projectile configuration file'
!                print *, '         is larger than that in input file.'
!                print '(A,I4,A)', '         Reading only first ', s%n_atoms, ' positions.'
!            end if
!                read(44,*) s%r
!
!            close(44)
!        case(3) ! Put atom in designated positions on surface
!
!            call lower_case(key_p_pos)
!            select case (key_p_pos)
!            case('top')
!                s%r(1:2,1) = (/0.,0./)
!            case('fcc')
!                s%r(1:2,1) = (cc1 + 2.0d0*cc2)/3.0d0
!            case('hcp')
!                s%r(1:2,1) = (2.0d0*cc1 + cc2)/3.0d0
!            case('bri')
!                s%r(1:2,1) = (cc1 + cc2)*0.5d0
!            end select
!            s%r(3,1) = height
!            v_pdof = sqrt(2.0d0*kB*Tsurf/mass_p)
!
!            s%v = 0.0d0
!            s%v(1,1) = normal(0.0d0,v_pdof)
!            s%v(2,1) = normal(0.0d0,v_pdof)
!            s%v(3,1) = normal(0.0d0,v_pdof)
!
!        case default
!
!        end select
!    end if
!
!    if (pip_sign .ne. -1 .or. pip_sign .ne. 3) then
!    !     Assign projectile velocities
!        vinc = sqrt(2.0d0*einc/mass_p)
!        s%v(1,:) =  vinc*sin(inclination)*cos(azimuth)
!        s%v(2,:) =  vinc*sin(inclination)*sin(azimuth)
!        s%v(3,:) = -vinc*cos(inclination)
!    end if
!
!end subroutine particle_init
!
end module md_init
