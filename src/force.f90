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

module force

    use constants
    use universe_mod, only : universe, minimg_beads, new_atoms, to_direct

    use pes_lj_mod,   only : compute_lj, compute_simple_lj
    use pes_emt_mod,  only : compute_emt
    use pes_ho_mod,   only : compute_ho
    use pes_rebo_mod, only : compute_rebo
    use rpmd,         only : do_ring_polymer_step
    use pes_nene_mod, only : compute_nene
    
    use md_init, only: set_atomic_indices, set_atomic_dofs, read_pes
    use useful_things, only: invert_matrix
    use run_config

    use, intrinsic :: iso_c_binding

    implicit none

contains

    subroutine full_energy_force_wrapper(natoms, nbeads, ntypes, &
                                r, simbox, f, epot, pes_file, pes_file_length, natoms_list, &
                                projectile_element, surface_element, is_proj)


        integer(kind=c_int), intent(IN) :: natoms, nbeads, ntypes                     
        real(kind=c_double), intent(IN) :: r(3,nbeads,natoms), simbox(3,3)
        !real(kind=c_double), intent(IN) :: m(natoms)
        real(kind=c_double), intent(OUT) :: f(3,nbeads,natoms)
        real(kind=c_double), intent(OUT) :: epot(nbeads)
        type(universe) :: atoms
        integer(kind=c_int), intent(IN) :: natoms_list(ntypes)  
        integer(kind=c_int), intent(IN) :: pes_file_length
        character(len=pes_file_length), intent(IN) :: pes_file
        character(len=1), intent(IN) ::projectile_element
        character(len=2), intent(IN) ::surface_element
        logical*8, intent(IN) :: is_proj(ntypes)
        ! Think mass is not needed if not running dynamics
        ! f and epot are output


        ! 1. build default atoms in universe type (does allocations)
        atoms = new_atoms(nbeads, natoms, ntypes)
        atoms%r(:,:,:)=r        ! positions
        atoms%simbox = simbox
        atoms%isimbox = invert_matrix(simbox)
        atoms%is_cart=.true.
        call set_atomic_indices(atoms, natoms_list) !atoms%indx 
        atoms%name(1) = projectile_element
        atoms%name(2) = surface_element
        atoms%is_proj = is_proj       
        

        ! 2. Build required simulation parameters
        simparams = new_simulation_parameters()
        simparams%pes_file = pes_file
        simparams%nprojectiles = 1
        simparams%nlattices = 1

        ! 3. read pes
        if (atoms%pes(1,1) == default_int) then
            call read_pes(atoms)
        endif


        ! 4. calculate energy and force
        call calc_force(atoms,energy_and_force)

        ! 5. return epot and force
        f = atoms%f
        epot = atoms%epot

        ! 6. deallocations
        deallocate(atoms%r)
        deallocate(atoms%v)
        deallocate(atoms%f)
        deallocate(atoms%a)
        deallocate(atoms%m)
        deallocate(atoms%is_fixed)
        deallocate(atoms%idx)
        deallocate(atoms%epot)
        ! deallocate(atoms%pes)
        deallocate(atoms%name)
        deallocate(atoms%algo)
        deallocate(atoms%is_proj)


    end subroutine full_energy_force_wrapper

    subroutine calc_force(atoms, flag)

        type(universe), intent(inout) :: atoms
        integer, intent(in)           :: flag
        character(len=*), parameter   :: err = "Error in calc_force(): "

        real(dp) :: temp_distance(atoms%nbeads), temp_vector(3, atoms%nbeads)
        integer :: i, j

        atoms%f    = 0.0_dp
        atoms%a    = 0.0_dp
        atoms%epot = 0.0_dp

        ! XXX remove for production
        if (atoms%nbeads > 1) call do_ring_polymer_step(atoms)

        ! gather distance and vector information before calling energy subroutines
        ! calculate only one half and then add to other half (with changed sign)
        if (.not. allocated(atoms%distances)) then
            allocate(atoms%distances(atoms%nbeads, atoms%natoms, atoms%natoms))
            allocate(atoms%vectors(3, atoms%nbeads, atoms%natoms, atoms%natoms))
        end if

        atoms%distances = 0.0_dp
        atoms%vectors   = 0.0_dp

        do j = 1, atoms%natoms-1
            do i = j+1, atoms%natoms
                call minimg_beads(atoms, i, j, temp_distance, temp_vector)
                atoms%distances(:,i,j) = temp_distance
                atoms%distances(:,j,i) = temp_distance
                atoms%vectors(:,:,i,j) = -temp_vector
                atoms%vectors(:,:,j,i) = +temp_vector
            end do
        end do

        ! calculate energy (and forces if flag is set)
        if (any(atoms%pes == pes_id_lj))        call compute_lj       (atoms, flag)
        if (any(atoms%pes == pes_id_simple_lj)) call compute_simple_lj(atoms, flag)
        if (any(atoms%pes == pes_id_emt))       call compute_emt      (atoms, flag)
        if (any(atoms%pes == pes_id_ho))        call compute_ho       (atoms, flag)
        if (any(atoms%pes == pes_id_rebo))      call compute_rebo     (atoms, flag)
        if (any(atoms%pes == pes_id_nene))      call compute_nene     (atoms, flag)

        if (flag == energy_and_force) call set_acceleration(atoms)

    end subroutine calc_force


    pure subroutine set_acceleration(atoms)

        type(universe), intent(inout) :: atoms
        integer :: i

        do i = 1, atoms%natoms
            where(.not. atoms%is_fixed(:,:,i)) atoms%a(:,:,i) = atoms%f(:,:,i)/atoms%m(i)
        end do

    end subroutine set_acceleration



end module force
