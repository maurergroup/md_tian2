program md_tian
        ! Purpose:
        !       Do molecular dynamics, Langevin dynamics, Ring Polymer dynamics
        !
        ! Date          	Author          	    History of Revison
        ! ====          	======          	    ==================
        ! 30.03.2017    	Marvin Kammler		    new data structure
        !                   Sascha Kandratsenka     and propagation methods
        !
        ! 18.02.2014    	Svenja M. Janke		    Original
        !			        Sascha Kandratsenka
        !			        Dan J. Auerbach
        !

    use force
    use md_algo
    use md_init
    use pes_lj_mod
    use constants
    use rpmd
    use universe_mod
    use fit
    use output_mod, only : output
    use trajectory_info
    use geometry_opt

    implicit none

    integer :: itraj, istep, i
    real(dp) :: tmp
    type(universe) :: atoms

    real(dp), allocatable :: cents(:,:)
    real(dp) :: vec(3)


    call simbox_init(atoms)

    select case (simparams%run)

        case ('min')

            call optimize_geometry(atoms, geometry_opt_fire)
            call output(atoms, 1, 1)

        case ('fit')

            call perform_fit(atoms)

        case ('md')

            do itraj = simparams%start, simparams%start+simparams%ntrajs-1

                call calc_force(atoms, energy_and_force)
                if (.not. allocated(cents)) allocate(cents(3, atoms%natoms))

                if (any(simparams%output_type == output_id_scatter)) then
                    call output(atoms, itraj, -1, "scatter_initial")
                end if



                print *, "Eref", atoms%epot

                do istep = 1, simparams%nsteps

                    ! core propagation
                    call propagate_1(atoms)
                    if (atoms%nbeads > 1) call do_ring_polymer_step(atoms)
                    call calc_force(atoms, energy_and_force)
                    call propagate_2(atoms)

                    ! output and exit conditionns
                    if (any(mod(istep, simparams%output_interval) == 0)) call output(atoms, itraj, istep)

                    if ( all(atoms%r(3,:,1) > simparams%proj_ul) .and. atoms%is_proj(atoms%idx(1)) &
                        .and. any(simparams%output_type == output_id_scatter)) exit

                    ! record bounces, lowest position, etc.
                    call collect_trajectory_characteristics(atoms, istep)

                    !if (mod(istep, 10)) print *, sum(sum(atoms%r(:,:,:), dim=2), dim=2)/atoms%natoms/atoms%nbeads

                end do
                close(78)

                if (any(simparams%output_type == output_id_scatter)) call output(atoms, itraj, istep, "scatter_final")

                if (itraj < simparams%start+simparams%ntrajs-1) then
                    call random_seed(put=randseed)  ! Seed random number generator
                    do i = 1, 100*(itraj+1)
                        call random_number(tmp)   ! rotate it according to trajectory number
                    end do

                    call prepare_next_traj(atoms)
                end if

            end do

        case default
            call abort

    end select


end program md_tian
