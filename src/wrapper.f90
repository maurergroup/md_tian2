module wrapper

    use constants
    use universe_mod, only : universe, minimg_beads, new_atoms, to_direct

    use pes_lj_mod,   only : compute_lj, compute_simple_lj
    use pes_emt_mod,  only : compute_emt
    use pes_ho_mod,   only : compute_ho
    use pes_rebo_mod, only : compute_rebo
    use pes_nene_mod, only : compute_nene
    
    use md_init, only: set_atomic_indices, set_atomic_dofs, read_pes
    use useful_things, only: invert_matrix
    use run_config
    use force

    use, intrinsic :: iso_c_binding

    type(universe) :: atoms

    contains

    subroutine wrapper_read_pes(natoms, nbeads, ntypes, &
        simbox, pes_file, pes_file_length, natoms_list, &
        projectile_element, surface_element, is_proj)

        implicit none

        integer(kind=c_int), intent(IN) :: natoms, nbeads, ntypes
        real(kind=c_double), intent(IN) :: simbox(3,3)
        integer(kind=c_int), intent(IN) :: natoms_list(ntypes)  
        integer(kind=c_int), intent(IN) :: pes_file_length
        character(len=pes_file_length), intent(IN) :: pes_file
        character(len=1), intent(IN) ::projectile_element
        character(len=2), intent(IN) ::surface_element
        logical*8, intent(IN) :: is_proj(ntypes)


        !write(*,*) simbox
         ! 1. build default atoms in universe type (does allocations)
        atoms = new_atoms(nbeads, natoms, ntypes)
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
        call read_pes(atoms)

    end subroutine wrapper_read_pes




    subroutine wrapper_energy_force(natoms, nbeads, r, f, epot)

        integer(kind=c_int), intent(IN) :: natoms, nbeads
        real(kind=c_double), intent(IN) :: r(3,nbeads,natoms)
        real(kind=c_double), intent(OUT) :: f(3,nbeads,natoms)
        real(kind=c_double), intent(OUT) :: epot(nbeads)

        ! 1. build default atoms in universe type (does allocations
        atoms%r(:,:,:)=r        ! positions        

        ! 4. calculate energy and force
        call calc_force(atoms,energy_and_force)

        ! 5. return epot and force
        f = atoms%f
        epot = atoms%epot

    end subroutine wrapper_energy_force


    subroutine wrapper_deallocations

        implicit none

        ! 6. deallocations
        deallocate(atoms%r)
        deallocate(atoms%v)
        deallocate(atoms%f)
        deallocate(atoms%a)
        deallocate(atoms%m)
        deallocate(atoms%is_fixed)
        deallocate(atoms%idx)
        deallocate(atoms%epot)
        deallocate(atoms%pes)
        deallocate(atoms%name)
        deallocate(atoms%algo)
        deallocate(atoms%is_proj)

    end subroutine wrapper_deallocations

end module wrapper