!############################################################################
! This routine is part of
! md_tian2 (Molecular Dynamics Tian Xia 2)
! (c) 2014-2021 Daniel J. Auerbach, Svenja M. Janke, Marvin Kammler,
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

! This module contains definitions of user types and all constants
module universe_mod

    use constants

    implicit none


    !  Type atoms
    !   structure to hold the position, velocity, force etc. for multiple atoms
    !       use rank2 array so positions, velocities, forces etc. are stored
    !       in sequential memory locations for efficient access
    !       each array should be allocated (3, n_beads, n_atom)
    !       mass array has length of n_atoms
    type universe

        integer                       :: natoms          ! number of atoms
        integer                       :: nbeads          ! number of beads per atom
        integer                       :: ntypes          ! proj species + latt species

        real(dp), allocatable         :: r(:,:,:)        ! positions
        real(dp), allocatable         :: v(:,:,:)        ! velocities
        real(dp), allocatable         :: f(:,:,:)        ! forces
        real(dp), allocatable         :: a(:,:,:)        ! acceleration
        real(dp), allocatable         :: m(:)            ! mass
        logical,  allocatable         :: is_fixed(:,:,:) ! mask array defining frozen atoms (T means atom can move, F means atom is frozen)

        integer, allocatable          :: idx(:)          ! index of atom type used in program
        character(len=3), allocatable :: name(:)         ! atomic name
        integer, allocatable          :: pes(:,:)        ! defines idx-dependent pes
        integer, allocatable          :: algo(:)         ! defines idx-dependent propagation algorithm
        logical, allocatable          :: is_proj(:)      ! defines projectile and lattice via idx

        real(dp), allocatable         :: epot(:)
        integer                       :: dof             ! degrees of freedom

        real(dp)                      :: simbox(3,3)     ! simulation cell
        real(dp)                      :: isimbox(3,3)    ! inverse simulation cell
        logical                       :: is_cart         ! if geometry is cartesian or direct

        real(dp), allocatable         :: distances(:,:,:)
        real(dp), allocatable         :: vectors(:,:,:,:)

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
        allocate(new_atoms%m(natoms))
        allocate(new_atoms%is_fixed(3,nbeads,natoms))

        allocate(new_atoms%idx(natoms))

        allocate(new_atoms%epot(nbeads))

        allocate(new_atoms%pes(ntypes,ntypes))
        allocate(new_atoms%name(ntypes))
        allocate(new_atoms%algo(ntypes))
        allocate(new_atoms%is_proj(ntypes))

        new_atoms%nbeads = nbeads
        new_atoms%natoms = natoms
        new_atoms%ntypes = ntypes
        new_atoms%name  = default_string
        new_atoms%r     = default_real
        new_atoms%v     = default_real
        new_atoms%f     = 0.0_dp
        new_atoms%m     = default_real
        new_atoms%a     = 0.0_dp
        new_atoms%is_fixed = .false.
        new_atoms%is_proj  = .false.
        new_atoms%is_cart  = .false.
        new_atoms%idx   = default_int
        new_atoms%dof   = default_int
        new_atoms%pes   = default_int
        new_atoms%algo  = default_int
        new_atoms%simbox  = default_real
        new_atoms%isimbox = default_real
    end function

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!  procedures that change the object !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




    subroutine create_repetitions(this, rep)

        use useful_things, only : invert_matrix

        type(universe), intent(inout) :: this
        integer, intent(in) :: rep(2)

        type(universe) :: other
        integer :: nreps, nprojs = 0
        integer :: i, j, start, end

        !        print *, "before"
        !        print *, this%r

        !        ! XXX: ONLY WORKS FOR rep(1) == rep(2)
        !        if (rep(1) /= rep(2)) then
        !            print *, "Repetition does not work yet for asymmetrical repetitions."
        !            call abort
        !        end if

        ! convert to direct coordinates
        if (this%is_cart) call to_direct(this)
        nreps = (1+2*rep(1)) * (1+2*rep(2))

        ! if projectile exists, put at beginning of new structure
        if (any(this%is_proj)) then
            do i = 1, this%natoms
                if (this%is_proj(this%idx(i))) then
                    nprojs = nprojs + 1
                end if
            end do
            !            print *, "nprojs", nprojs
            other = new_atoms(this%nbeads, nreps*(this%natoms-nprojs)+nprojs, this%ntypes)
            other%r(:,:,:nprojs) = this%r(:,:,:nprojs)
            other%v(:,:,:nprojs) = this%v(:,:,:nprojs)
            other%a(:,:,:nprojs) = this%a(:,:,:nprojs)
            other%f(:,:,:nprojs) = this%f(:,:,:nprojs)
            other%m(:nprojs) = this%m(:nprojs)
            other%is_fixed(:,:,:nprojs) = this%is_fixed(:,:,:nprojs)
            other%idx(:nprojs) = this%idx(:nprojs)

        else
            other = new_atoms(this%nbeads, nreps*this%natoms, this%ntypes)
        end if

        other%pes     = this%pes
        other%algo    = this%algo
        other%name    = this%name
        other%simbox  = this%simbox
        other%is_proj = this%is_proj

        start = nprojs + 1
        end = this%natoms
        do i = -rep(1), rep(1)
            do j = -rep(2), rep(2)
                ! This map -i, -i+1, ... , i-1, i and -j, -j+1, ... , j-1, j
                !  to 1, 2, 3, ... , #total_repetitions

                other%r( 1 , : , start : end ) = this%r(1,:,nprojs+1:) + i
                other%r( 2 , : , start : end ) = this%r(2,:,nprojs+1:) + j
                other%r( 3 , : , start : end ) = this%r(3,:,nprojs+1:)

                other%m(         start : end ) = this%m(nprojs+1:)
                other%idx(       start : end ) = this%idx(nprojs+1:)
                other%v( : , : , start : end ) = this%v(:,:,nprojs+1:)
                other%f( : , : , start : end ) = this%f(:,:,nprojs+1:)
                other%is_fixed( : , : , start : end ) = this%is_fixed(:,:,nprojs+1:)

                start = end + 1
                end = end + this%natoms - nprojs
            end do
        end do


        !        print *, "after"
        !        print *, other%r

        ! bring back to cartesian
        call to_cartesian(other)

        ! enlarge the cell
        other%simbox(:,1) = this%simbox(:,1) * (1+2*rep(1))
        other%simbox(:,2) = this%simbox(:,2) * (1+2*rep(2))

        ! obtain new inverse cell matrix
        other%isimbox = invert_matrix(other%simbox)

        ! set DOFs
        other%dof = 3*other%natoms*other%nbeads - count(other%is_fixed)

        ! replace old atoms
        this = other

    end subroutine create_repetitions




    subroutine remove_com_velocity(this)

        type(universe), intent(inout) :: this
        real(dp) :: v_cm(3)
        integer :: i, b, dummy

        !!! FIXME: Somehow this needs to be applied very often to work

        do dummy = 1, 1000

            v_cm = calc_com_velocity(this)

            do i = 1, this%natoms
                do b = 1, this%nbeads
                    this%v(:,b,i) = this%v(:,b,i) - v_cm
                end do
            end do
        end do

    end subroutine remove_com_velocity




    subroutine to_cartesian(this)

        type(universe), intent(inout) :: this
        integer :: i,j

        if (this%is_cart) return

        do i = 1, this%natoms
            do j = 1, this%nbeads
                this%r(:,j,i) = matmul(this%simbox, this%r(:,j,i))
            end do
        end do

        this%is_cart = .not. this%is_cart

    end subroutine to_cartesian




    subroutine to_direct(this)

        type(universe), intent(inout) :: this
        integer :: i,j

        if (.not. this%is_cart) return

        do i = 1, this%natoms
            do j = 1, this%nbeads
                this%r(:,j,i) = matmul(this%isimbox, this%r(:,j,i))
            end do
        end do

        this%is_cart = .not. this%is_cart

    end subroutine to_direct







     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!  procedures that derive quantities !!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    function calc_com_velocity(this) result(v_cm)

        ! v_cm = sum(m_i*v_i) / sum(m_i)

        type(universe), intent(in) :: this
        real(dp), dimension(3) :: v_cm, num
        integer :: i

        num = 0.0_dp
        do i = 1, this%natoms
            num = num + sum(this%v(:,:,i), dim=2) * this%m(i)
        end do

        v_cm = num / sum(this%m) / this%nbeads / this%natoms

    end function calc_com_velocity




    real(dp) function calc_atom_temperature(this) result(temperature)

        type(universe), intent(in) :: this
        real(dp)                   :: temp
        integer :: i

        temp = 0.0_dp
        do i = 1, this%natoms
            temp = temp + this%m(i) * sum(sum(this%v(:,:,i), dim=2)**2)
        end do

        temperature = temp / kB / this%dof / this%nbeads

    end function calc_atom_temperature





    subroutine calc_momentum_all(this, p)

        type(universe), intent(in)  :: this
        real(dp), intent(out) :: p(3, this%nbeads, this%natoms)

        integer :: i

        do i = 1, this%natoms
            p(:,:,i) = this%m(i) * this%v(:,:,i)
        end do

    end subroutine calc_momentum_all




    function calc_momentum_one(this, i) result(momentum)

        type(universe), intent(in)  :: this
        integer       , intent(in)  :: i
        real(dp)                    :: momentum(3, this%nbeads)

        momentum = this%m(i) * this%v(:,:,i)

    end function calc_momentum_one




    integer function get_idx_from_name(this, name, is_proj) result(idx)

        type(universe), intent(in) :: this
        character(len=3), intent(in) :: name
        logical, intent(in) :: is_proj

        integer :: i

        idx = default_int
        do i = 1, this%natoms
            !                print *, i, this%name(this%idx(i)), name, this%is_proj(this%idx(i)), is_proj, this%idx(i)
            if ((this%name(this%idx(i)) == name) .and. (this%is_proj(this%idx(i)) .eqv. is_proj)) then
                idx = this%idx(i)
                exit
            end if
        end do

        if (idx == default_int) then
            write(*,*) 'Err: ', name, is_proj
            stop "Error in get_idx_from_name(): make sure you &
                        correctly assign element names to projectile and slab in both &
                        *.inp and *.pes files."
        endif

    end function get_idx_from_name




    subroutine minimg_one(this, i, j, bi, bj, r, vec)

        ! Calculates the minimum image vector and distance
        ! between beads b of atoms i and j

        type(universe), intent(in)         :: this
        integer, intent(in)                :: i, j, bi, bj
        real(dp),               intent(out) :: r
        real(dp), dimension(3), intent(out) :: vec

        vec = this%r(:,bj,j)-this%r(:,bi,i)   ! distance vector from a to b

        vec = matmul(this%isimbox, vec)   ! transform to direct coordinates
        vec = vec - anint(vec)            ! imaging
        vec = matmul(this%simbox, vec)    ! back to cartesian coordinates

        r =  sqrt(sum(vec*vec))      ! distance

    end subroutine minimg_one




    pure subroutine minimg_beads(this, i, j, r, vec)

        ! Calculates the minimum image vector and distance
        ! between all beads of atoms i and j

        type(universe), intent(in) :: this
        integer, intent(in)        :: i, j
        real(dp), intent(out)       :: r(this%nbeads)
        real(dp), intent(out)       :: vec(3, this%nbeads)

        vec = this%r(:,:,j)-this%r(:,:,i)   ! distance vector from a to b

        vec = matmul(this%isimbox, vec)   ! transform to direct coordinates
        vec = vec - anint(vec)            ! imaging
        vec = matmul(this%simbox, vec)    ! back to cartesian coordinates
        r =  sqrt(sum(vec*vec, dim=1))    ! distance

    end subroutine minimg_beads




    subroutine atom_ekin(this, ekin_p, ekin_l)

        type(universe), intent(in) :: this
        real(dp), intent(out)      :: ekin_p, ekin_l
        real(dp)                   :: nrg
        integer :: i

        ekin_p = 0.0_dp
        ekin_l = 0.0_dp

        do i = 1, this%natoms

            nrg = this%m(i)*sum(sum(this%v(:,:,i), dim=2)**2)

            if (this%is_proj(this%idx(i))) then
                ekin_p = ekin_p + nrg
            else
                ekin_l = ekin_l + nrg
            end if

        end do

        ekin_p = 0.5_dp * ekin_p / this%nbeads / this%nbeads
        ekin_l = 0.5_dp * ekin_l / this%nbeads / this%nbeads

    end subroutine atom_ekin


    function pes_id_to_name(id) result(name)

        integer, intent(in) :: id
        character(len=5) :: name
        character(len=*), parameter :: err = "Error in pes_id_to_name(): "

        select case(id)
            case (pes_id_lj)
                name = pes_name_lj
            case (pes_id_morse)
                name = pes_name_morse
            case (pes_id_emt)
                name = pes_name_emt
            case (pes_id_rebo)
                name = pes_name_rebo
            case (pes_id_ho)
                name = pes_name_ho
            case (pes_id_simple_lj)
                name = pes_name_simple_lj
            case (pes_id_no_interaction)
                name = pes_name_no_interaction
            case (pes_id_nene)
                name = pes_name_nene
            case default
                print *, err // "unknown id", id
                stop
        end select

    end function pes_id_to_name



    function universe_id_to_role(this, id) result(role)

        type(universe), intent(in) :: this
        integer, intent(in) :: id
        character(len=4) :: role

        if (this%is_proj(this%idx(id))) then
            role = "proj"
        else
            role = "latt"
        end if

    end function universe_id_to_role


end module universe_mod
