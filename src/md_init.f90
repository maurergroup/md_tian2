!############################################################################
! This routine is part of
! md_tian2 (Molecular Dynamics Tian Xia 2)
! (c) 2014-2020 Dan J. Auerbach, Svenja M. Janke, Marvin Kammler,
!               Sascha Kandratsenka, Sebastian Wille
! Dynamics at Surfaces Department
! MPI for Biophysical Chemistry Goettingen, Germany
! Georg-August-Universitaet Goettingen, Germany
!
! This program is free software: you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by the
! Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
! or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
! for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program. If not, see http://www.gnu.org/licenses.
!############################################################################

! Initialize molecular dynamics simulation
module md_init

    use universe_mod
    use useful_things
    use run_config
    use constants
    use open_file

    use pes_lj_mod,   only : read_lj, read_simple_lj
    use pes_ho_mod,   only : read_ho
    use pes_emt_mod,  only : read_emt
    use pes_non_mod,  only : read_non_interacting
    use pes_rebo_mod, only : read_rebo
    use pes_nene_mod, only : read_nene

    implicit none

    integer, parameter :: pes_unit = 38
    integer, parameter :: geo_unit = 55

contains

    subroutine simbox_init(atoms)

        type(universe), intent(out) :: atoms
        integer*8 :: seed

        character(len=:), allocatable :: input_file
        integer :: input_file_length, input_file_status
        integer :: i

        ! Read in name of input file
        if (command_argument_count() == 0) stop " I need an input file"
        if (command_argument_count() .gt. 1) stop "Only one input file is allowed"

        call get_command(length=input_file_length)
        allocate(character(input_file_length) :: input_file)
        call get_command_argument(1, input_file, input_file_length, input_file_status)
        if (input_file_status /= 0) stop " Error by reading the command line"


        ! build the simparams object
        call read_input_file(input_file)
        call ensure_input_sanity()
        
        call read_geometry(atoms, simparams%confname_file, seed)
        call ensure_geometry_sanity(atoms)

        call read_pes(atoms)
        !call remove_com_velocity(atoms)
        call post_process(atoms)



    end subroutine simbox_init




    subroutine read_pes(atoms)

        type(universe), intent(inout) :: atoms

        integer :: ios = 0, nwords, i
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

                        case (pes_name_lj)
                            call read_lj(atoms, pes_unit)

                        case (pes_name_simple_lj)
                            call read_simple_lj(atoms, pes_unit)

                        !case ('morse')
                        !    call read_morse(pes_unit)

                        case (pes_name_emt)
                            call read_emt(atoms, pes_unit)

                        case (pes_name_rebo)
                            call read_rebo(atoms, pes_unit)

                        case (pes_name_no_interaction)
                            call read_non_interacting(atoms, pes_unit)

                        case (pes_name_ho)
                            call read_ho(atoms, pes_unit)

                        case (pes_name_nene)
                            call read_nene(atoms, pes_unit)

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
    subroutine read_geometry(atoms, infile, seed)

        use useful_things, only : file_exists

        character(len=*), intent(in) :: infile
        type(universe), intent(inout) :: atoms
        integer*8 :: seed

        type(universe) :: proj

        if (.not. file_exists(infile)) stop "Error: geometry file does not exist"

        select case (simparams%confname)

            case ("poscar")
                call read_poscar(atoms, infile)

            case ("mxt", "merge")
                call read_mxt(atoms, infile)

            case default
                stop "Error: conf keyword unknown"

        end select

        ! merge projectile file if present
        if (simparams%confname == "merge") then
            call read_mxt(proj, simparams%merge_proj_file)
            call merge_universes(atoms, proj, seed)
        end if

        ! prepare initial projectile velocity vector (if provided)
        if (any(atoms%is_proj)) call set_proj_incidence(atoms, seed)

    end subroutine read_geometry




    subroutine read_poscar(atoms, infile)

        type(universe), intent(out) :: atoms
        character(len=*), intent(in) :: infile

        character(len=*), parameter :: err = "Error in read_poscar: "
        real(dp) :: scaling_const
        real(dp) :: simbox(3,3)
        integer, allocatable :: natoms(:)
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
        if (ios /= 0) stop err // "reading comment line"

        read(geo_unit, *, iostat=ios) scaling_const             ! second line is scaling constant
        if (ios /= 0) stop err // "reading scaling constant"

        ! read the cell vectors and multiply with scaling constant
        read(geo_unit, *, iostat=ios) simbox
        if (ios /= 0) stop err // "reading simulation box"

        ! read atom types
        read(geo_unit, '(A)', iostat=ios) buffer
        if (ios /= 0) stop err // "reading atom type"
        call split_string(buffer, elements, nwords)

        ! read number of atoms and number of beads
        read(geo_unit, '(A)', iostat=ios) buffer
        if (ios /= 0) stop err // "reading number of atoms and beads"
        pos1 = scan(string=buffer, set=":", back=.false.)
        if (pos1 > 0) then  ! Number of beads is stated in the file
            pos2 = scan(string=buffer, set=":", back=.true.)
            read(buffer(pos1+1:pos2-1), '(i10)', iostat=ios) nbeads
            if (ios /= 0 .or. nbeads < 1) stop err // "reading number of beads"
        end if
        read(buffer(pos2+1:), *, iostat=ios) natoms
        if (any(natoms <= 0)) stop err // "number of atoms could not be read"
        if (ios /= 0) stop err // "reading number of atoms"

        ! now we have everything to construct atom types
        atoms = new_atoms(nbeads, sum(natoms), ntypes)

        ! Direct or Cartesian coordinates, needed later
        read(geo_unit, '(A)', iostat=ios) c_system
        if (ios /= 0) stop err // "reading coordinate system type (cartesian/direct)"

        ! if coordinate section has 3 entries per line, no atoms are fixed
        ! if coordinate section has 6 entries per line, some might be fixed
        read(geo_unit, '(A)', iostat=ios) buffer
        if (ios /= 0) stop err // "reading start of coordinate section"
        backspace(geo_unit)
        call split_string(buffer, words, nwords)

        if (nwords == 3) then

            !read(geo_unit, *, iostat=ios) ((atoms%r(:,j,i), j=1,atoms%nbeads), i=1,atoms%natoms)
            read(geo_unit, *, iostat=ios) atoms%r
            if (ios /= 0) stop err // "reading coordinates"

        else if (nwords == 6) then

            if (atoms%natoms > 0) then
                allocate (can_move(3, atoms%nbeads, atoms%natoms))
                read(geo_unit, *, iostat=ios) ((atoms%r(:,j,i), can_move(:,j,i), j=1,atoms%nbeads), &
                    i=1,atoms%natoms)
                if (ios /= 0) stop err // "reading atom coordinates and mobility flags"
                where (can_move == "F") atoms%is_fixed = .true.
                deallocate (can_move)
            else
                print *, err // "reading projectile coordinates and mobility flags"
                stop
            end if

        else
            print *, err // "cannot read coordinate section"
            stop
        end if

        ! Check if there is a velocity section
        read(geo_unit, '(A)', iostat=ios) buffer
        if (ios == iostat_end) then
            ! eof of file, don't read anything more
            atoms%v = 0.0_dp

        else if (ios == 0 .and. len_trim(buffer) == 0) then
            ! blank line after coordinate section and no eof yet -> read velocities
            read(geo_unit, *, iostat=ios) atoms%v
            if (ios /= 0) stop err // "cannot read initial atom velocities"

        else
            print *,  err // "cannot read initial atom velocities. Did you specify &
            all species in the *.inp file?"
            stop
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
            print *, err // "reading coordinate system type (cartesian/direct)"
            stop
        end if
        call set_atomic_dofs(atoms)
        call set_atomic_indices(atoms, natoms)
        call set_atomic_masses(atoms, natoms)
        call set_atomic_names(atoms, elements)
        call set_prop_algos(atoms)
        call create_repetitions(atoms, simparams%rep)

    end subroutine read_poscar




    subroutine read_mxt(atoms, infile)

        type(universe), intent(out)  :: atoms
        character(len=*), intent(in) :: infile

        character(len=*), parameter  :: err = "Error in read_mxt(): "
        integer :: natoms, nbeads, ntypes, ios, i
        integer, allocatable :: atom_list(:)

        if (.not. file_exists(infile)) then
            print *, err, "file ", infile, " does not exist"
            stop
        end if

        open(geo_unit, file=infile, form="unformatted", iostat=ios)

        read(geo_unit, iostat=ios) natoms, nbeads, ntypes
        if (ios /= 0) stop err // "cannot read universe type constructor arguments"
        atoms = new_atoms(nbeads, natoms, ntypes)

        read(geo_unit, iostat=ios) atoms%r, atoms%v, atoms%is_cart, atoms%is_fixed, atoms%idx, &
            atoms%name, atoms%is_proj, atoms%simbox, atoms%isimbox
        if (ios /= 0) stop err // "cannot read file"

        close(geo_unit)

        allocate(atom_list(atoms%ntypes))
        atom_list = 0
        !print *, atoms%ntypes
        do i = 1, atoms%natoms
            atom_list(atoms%idx(i)) = atom_list(atoms%idx(i)) + 1
        end do

        call set_atomic_dofs(atoms)
        call set_atomic_masses(atoms, atom_list)
        call set_prop_algos(atoms)

    end subroutine read_mxt




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
            if (simparams%nlattices == default_int .and. simparams%nprojectiles == default_int) &
                stop err // "either lattice and/or projectile key must be present."
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

!            if (.not. allocated(simparams%output_type) .or. .not. allocated(simparams%output_interval)) &
!                stop err // "output options not set."
            if (any(simparams%output_interval > simparams%nsteps)) print *, warn // "one or more outputs will not show, &
                since the output interval is larger than the total number of simulation steps."

            ! If one of them is set, the others must be set as well.
            if (any([simparams%einc, simparams%polar] == default_real) .and.  &
                .not. all([simparams%einc, simparams%polar] == default_real)) &
                stop err // "incidence conditions ambiguous."

            if (simparams%azimuth == default_string .and. simparams%nprojectiles > 0) stop err // "incidence azimuth angle not set"

            if (any(simparams%pip == default_string) /= all(simparams%pip == default_string)) stop err // "you have set all pip arguments"

            if (simparams%force_beads /= default_int .and. simparams%confname == "mxt") stop err // "cannot force beads on mxt files"

            if (simparams%confname == "merge" .and. .not. (allocated(simparams%mass_l) .and. allocated(simparams%mass_p))) &
                stop err // "both projectile and slab have to be present in the input file to use <conf merge>"
            if (simparams%nthreads /= 1) print *, warn, "only one thread supported for md"

            if (simparams%force_beads /= default_int .and. &
                simparams%Tsurf == default_real .and. simparams%Tproj == default_real) &
                stop err // "projectile and/or slab temperature required for rpmd"


                ! TODO: to be completed

        else if (simparams%run == "min") then
            if (simparams%nlattices == default_int .and. &
                simparams%nprojectiles == default_int) stop err // &
                "either lattice and/or projectile key must be present."
            if (simparams%nprojectiles == default_int) simparams%nprojectiles = 0
            if (simparams%nlattices == default_int) simparams%nlattices = 0
            if (simparams%nthreads /= 1) print *, warn, "only one thread is supported for minimization"



        !!! Perform fit
        else if (simparams%run == "fit") then
            if (simparams%nlattices == default_int .and. &
                simparams%nprojectiles == default_int) stop err // &
                "either lattice and/or projectile key must be present."
            if (simparams%nprojectiles == default_int) simparams%nprojectiles = 0
            if (simparams%nlattices == default_int) simparams%nlattices = 0
            if (simparams%fit_training_data == default_int) stop err // "training data input missing"
            if (simparams%fit_validation_data == default_int) stop err // "validation data input missing"
            if (simparams%evasp == default_real) stop err // "reference energy evasp missing"
            if (simparams%start == default_int) stop err // "fit number not set"
            if (simparams%fit_training_folder == simparams%fit_validation_folder) stop err // "cannot use the same folder for training and validation"

        else
            print *, err // "unknown run command", simparams%run
            stop

        end if



    end subroutine ensure_input_sanity



    subroutine ensure_geometry_sanity(atoms)

        type(universe), intent(in) :: atoms
        character(len=*), parameter :: err = "geometry sanity check error: "

        if (simparams%force_beads < 1) stop err // "Cannot force zero or less beads."
        if (simparams%force_beads /= default_int .and. atoms%nbeads /= 1) stop err // "Cannot force any beads when system already contains multiple beads."
        if (atoms%nbeads > 1 .and. any(atoms%is_proj) .and. simparams%Tproj == default_real) stop err // "RPMD requires projectile temperature"
        if (atoms%nbeads > 1 .and. any(.not. atoms%is_proj) .and. simparams%Tsurf == default_real) stop err // "RPMD requires surface temperature"

    end subroutine ensure_geometry_sanity




    subroutine apply_pip(atoms, seed)

        use rpmd, only : calc_centroid_positions

        type(universe), intent(inout) :: atoms
        integer*8 :: seed
        real(dp) :: target_position, cents(3, atoms%natoms)

        integer  :: k, b, ios
        character(len=*), parameter :: err = "Error in apply_pip(): "
        character(len=*), parameter :: random_position = "r"

        call to_cartesian(atoms)
        cents = calc_centroid_positions(atoms)

        if (.not. all(atoms%is_proj)) stop err // "there are lattice atoms affected by apply_pip() subroutine"

        do k = 1, dimensionality
            ! draw random number from direct coordinates, then convert to cartesian
            if (simparams%pip(k) == random_position) then
                target_position = ranx(simparams%nran,seed,1)
                target_position = target_position * sum(atoms%simbox(:,k))

            else
                read(simparams%pip(k), *, iostat=ios) target_position
                if (ios /= 0) then
                    print *, err, "unknown pip format at position", k, "of pip array"
                    stop
                end if
            end if

            do b = 1, atoms%nbeads
                atoms%r(k,b,:) = atoms%r(k,b,:) - cents(k,:) + target_position
            end do

        end do

    end subroutine apply_pip




    subroutine merge_universes(this, other, seed)
        ! add other universe to this universe
        ! projectile must be other

        use rpmd, only : calc_centroid_positions

        type(universe), intent(inout) :: this
        type(universe), intent(inout) :: other
        integer*8 :: seed

        type(universe) :: new

        integer :: b
        real(dp), dimension(3, this%nbeads, this%natoms+other%natoms) :: dir_coords
        real(dp), dimension(3, this%natoms+other%natoms) :: centroids
        character(len=*), parameter :: err = "Error in merge_universes: "

        if (this%nbeads /= other%nbeads) stop err // "mxt and merge files differ in number of beads"
        if (any(this%is_proj))           stop err // "mxt file contains projectile atoms"
        if (any(.not. other%is_proj))    stop err // "merge file contains lattice atoms"
        if (other%natoms > 1)            stop err // "multiple projectiles not implemented for use with merge option"


        if (.not. this%is_cart)  call to_cartesian(this)
        if (.not. other%is_cart) call to_cartesian(other)

        new = new_atoms(this%nbeads, this%natoms+other%natoms, this%ntypes+other%ntypes)

        if (.not. all(simparams%pip == default_string)) then
            call apply_pip(other, seed)
        end if

        new%r(:,:,:other%natoms)        = other%r
        new%v(:,:,:other%natoms)        = other%v
        new%f(:,:,:other%natoms)        = other%f
        new%a(:,:,:other%natoms)        = other%a
        new%m(:other%natoms)            = other%m
        new%idx(:other%natoms)          = other%idx
        new%name(:other%ntypes)         = other%name
        new%is_proj(:other%ntypes)      = other%is_proj
        new%is_fixed(:,:,:other%natoms) = other%is_fixed

        new%r(:,:,other%natoms+1:)        = this%r
        new%v(:,:,other%natoms+1:)        = this%v
        new%f(:,:,other%natoms+1:)        = this%f
        new%a(:,:,other%natoms+1:)        = this%a
        new%m(other%natoms+1:)            = this%m
        new%idx(other%natoms+1:)          = this%idx+maxval(other%idx)
        new%name(other%ntypes+1:)         = this%name
        new%is_proj(other%ntypes+1:)      = this%is_proj
        new%is_fixed(:,:,other%natoms+1:) = this%is_fixed

        new%algo = [simparams%md_algo_p, simparams%md_algo_l]
        new%dof     = 3*new%natoms*new%nbeads - count(new%is_fixed)
        new%simbox  = this%simbox
        new%isimbox = this%isimbox
        new%is_cart = .true.

        this = new

        ! fold all centroids into simulation box
        centroids = calc_centroid_positions(this)
        centroids = matmul(this%isimbox, centroids)
        call to_direct(this)
        forall (b = 1 : this%nbeads) this%r(:,b,:) = this%r(:,b,:) - int(centroids)
        call to_cartesian(this)

    end subroutine merge_universes




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




    subroutine remove_atomic_com(atoms, i)

        type(universe), intent(inout) :: atoms
        integer, intent(in)           :: i

        integer  :: b
        real(dp) :: v_cm(3)

        v_cm = sum(atoms%v(:,:,i), dim=2) / atoms%nbeads

        do b = 1, atoms%nbeads
            atoms%v(:,b,i) = atoms%v(:,b,i) - v_cm
        end do

    end subroutine remove_atomic_com




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
            print *, err // "neither projectile nor slab exist"
            stop
        end if

        atoms%name = elements

    end subroutine set_atomic_names




    subroutine set_atomic_masses(atoms, natoms)

        type(universe), intent(inout) :: atoms
        integer, intent(in) :: natoms(:)

        integer :: i
        real(dp) :: masses(size(natoms))
        character(len=*), parameter :: err = "Error in set_atomic_masses: "

        if (simparams%confname == "poscar" .or. simparams%confname == "mxt") then
            if (allocated(simparams%mass_p) .and. allocated(simparams%mass_l)) then
                masses = [simparams%mass_p, simparams%mass_l]
            else if (allocated(simparams%mass_p)) then
                masses = simparams%mass_p
            else if (allocated(simparams%mass_l)) then
                masses = simparams%mass_l
            else
                print *, err // "neither projectile nor slab exist"
                stop
            end if

        else if (simparams%confname == "merge") then
            if (all(atoms%is_proj)) then
                masses = simparams%mass_p
            else if (all(.not. atoms%is_proj)) then
                masses = simparams%mass_l
            else
                print *, err, "cannot set masses, conf merge requires seperate projectile and slab mxt files"
                stop
            end if

        else
            print *, err, "unknown confname", simparams%confname
            stop
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
            if (size(atoms%algo) /= size([simparams%md_algo_p, simparams%md_algo_l]) .and. simparams%confname /= "merge") stop err // "more algorithms than species"

        else if (allocated(simparams%md_algo_p)) then
            atoms%algo = simparams%md_algo_p
            if (size(atoms%algo) /= size(simparams%md_algo_p)) stop err // "more algorithms than species"

        else if (allocated(simparams%md_algo_l)) then
            atoms%algo = simparams%md_algo_l
            if (size(atoms%algo) /= size(simparams%md_algo_l)) stop err // "more algorithms than species"

        else
            print *, err // "neither projectile nor slab exist"
            stop
        end if

    end subroutine set_prop_algos




    subroutine set_proj_incidence(atoms, seed)

        type(universe), intent(inout) :: atoms
        integer*8 :: seed

        real(dp) :: vinc, drawn_einc, proj_mass, vx, vy, vz, chosen_azimuth
        integer  :: i, ios
        character(len=*), parameter :: err = "Error in set_proj_incidence: "
        character(len=*), parameter :: random_position = "r"

        ! return if no information is provided
        if (all([simparams%einc, simparams%polar] == default_real)) return

        ! draw incidence energy from normal distribution
        call normal_deviate(simparams%einc, simparams%sigma_einc, drawn_einc)

        ! convert to velocity
        proj_mass = 0.0_dp
        do i = 1, atoms%natoms
            if (atoms%is_proj(atoms%idx(i))) proj_mass = proj_mass + atoms%m(i)
        end do

        ! set the azimuth angle
        if (simparams%azimuth == random_position) then
            chosen_azimuth = ranx(simparams%nran,seed,1)
            chosen_azimuth = 2 * pi * chosen_azimuth
        else
            read(simparams%azimuth, *, iostat=ios) chosen_azimuth
            if (ios /= 0) stop err // "first azimuth argument must be a number or r"
            call normal_deviate_0d(chosen_azimuth, simparams%sigma_azimuth, chosen_azimuth)
            chosen_azimuth = chosen_azimuth * deg2rad
        end if

        vinc =  sqrt(2*drawn_einc/proj_mass)
        vx   =  vinc * sin(simparams%polar) * cos(chosen_azimuth)
        vy   =  vinc * sin(simparams%polar) * sin(chosen_azimuth)
        vz   = -vinc * cos(simparams%polar)

        do i = 1, atoms%natoms
            if (atoms%is_proj(atoms%idx(i))) then
                call remove_atomic_com(atoms, i)
                atoms%v(1,:,i) = atoms%v(1,:,i) + vx
                atoms%v(2,:,i) = atoms%v(2,:,i) + vy
                atoms%v(3,:,i) = atoms%v(3,:,i) + vz
            else
                return    ! projectiles are always at the beginning of universe object
            end if
        end do

    end subroutine set_proj_incidence




    subroutine set_atomic_dofs(atoms)

        type(universe), intent(inout) :: atoms

        atoms%dof = 3*atoms%natoms*atoms%nbeads - count(atoms%is_fixed)

    end subroutine set_atomic_dofs




    subroutine prepare_next_traj(atoms, seed)

        type(universe), intent(inout) :: atoms
        integer*8 :: seed

        character(len=*), parameter :: err = "Error in prepare_next_traj(): "
        integer :: mxt_idx, dat_idx, new_conf
        real(dp) :: rnd
        type(universe) :: proj, slab

        if (simparams%confname == "merge") then

            ! renew simparams projectile geometry file
            mxt_idx = index(simparams%merge_proj_file, "mxt_")
            dat_idx = index(simparams%merge_proj_file, ".dat")
            if (mxt_idx /= 0 .and. dat_idx /= 0 .and. dat_idx == mxt_idx+12) then
                rnd = ranx(simparams%nran,seed,1)
                new_conf = int(rnd*simparams%merge_proj_nconfs)+1
                write(simparams%merge_proj_file, '(a, i8.8, a)') simparams%merge_proj_file(:mxt_idx+3), new_conf, ".dat"
            else
                print *, err, "projectile confname_file has wrong format. It needs to end on /mxt_%08d.dat"; stop
            end if
                
            ! renew simparams slab geometry file
            mxt_idx = index(simparams%confname_file, "mxt_")
            dat_idx = index(simparams%confname_file, ".dat")
            if (mxt_idx /= 0 .and. dat_idx /= 0 .and. dat_idx == mxt_idx+12) then
                rnd = ranx(simparams%nran,seed,1)
                new_conf = int(rnd*simparams%nconfs)+1
                write(simparams%confname_file, '(a, i8.8, a)') simparams%confname_file(:mxt_idx+3), new_conf, ".dat"
            else
                print *, err, "lattice confname_file has wrong format. It needs to end on /mxt_%08d.dat"; stop
            end if

            ! read the geometry files and combine
            call read_mxt(proj, simparams%merge_proj_file)
            call read_mxt(slab, simparams%confname_file)

            call merge_universes(slab, proj, seed)
            slab%pes = atoms%pes
            atoms = slab

        else
            print *, err, "multiple trajectories only available for <conf merge> option"; stop
        end if

        if (any(atoms%is_proj)) call set_proj_incidence(atoms, seed)

    end subroutine prepare_next_traj



end module md_init
