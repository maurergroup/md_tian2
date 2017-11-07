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
    use output_mod, only : output
    use geometry_opt

    implicit none

    integer :: itraj, istep, i
    real(dp) :: tmp, proj_v(3), dotprod, old_v(3) = [0, 0, 0]
    type(universe) :: atoms

    call simbox_init(atoms)

    select case (simparams%run)

        case ('min')

            call optimize_geometry(atoms, geometry_opt_fire)
            call output(atoms, 1, 1)

        case ('md')
            !open(unit=78, file="bounces.dat")

            do itraj = simparams%start, simparams%start+simparams%ntrajs-1


                call calc_force(atoms, energy_and_force)
                if (any(simparams%output_type == output_id_scatter)) call output(atoms, itraj, -1, "scatter_initial")
                print *, "Eref", atoms%epot

                do istep = 1, simparams%nsteps

                    call propagate_1(atoms)
                    if (atoms%nbeads > 1) call do_ring_polymer_step(atoms)
                    call calc_force(atoms, energy_and_force)
                    call propagate_2(atoms)

                    if (any(mod(istep, simparams%output_interval) == 0)) call output(atoms, itraj, istep)

                    if ( sum(atoms%r(3,:,1))/atoms%nbeads > simparams%proj_ul .and. atoms%is_proj(atoms%idx(1)) &
                        .and. any(simparams%output_type == output_id_scatter)) exit


!                    !!! test bounce recognition !!!
!                    proj_v = sum(atoms%v(:,:,1), dim=2)/atoms%nbeads
!                    if (istep > 1) then
!                        dotprod = sum(old_v/sqrt(sum(old_v**2)) * proj_v/sqrt(sum(proj_v**2)))
!                        write(78, '(i, 3f12.6)') istep, &   ! time step
!                            180-acos(proj_v(3)/sqrt(sum(proj_v*proj_v)))*rad2deg, & ! polar
!                            atan2(proj_v(2), proj_v(1))*rad2deg, & ! azi
!                            dotprod
!                    end if
!                    old_v = proj_v


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
