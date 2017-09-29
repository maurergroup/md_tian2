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
    use rpmd
    use universe_mod
    use output_mod, only : output
    use geometry_opt

    implicit none

    integer :: itraj, istep
    real(dp) :: tmp
    type(universe) :: atoms

    select case (simparams%run)

        case ("min")

            call simbox_init(atoms)
            call optimize_geometry(atoms, geometry_opt_fire)
            call output(atoms, 1, 1)

        case ('md')

            do itraj = simparams%start, simparams%start+simparams%ntrajs-1

                call simbox_init(atoms)
                call calc_force(atoms, energy_and_force)
                print *, "Eref", atoms%epot

                do istep = 1, simparams%nsteps

                    call propagate_1(atoms)
                    if (atoms%nbeads > 1) call do_ring_polymer_step(atoms)
                    call calc_force(atoms, energy_and_force)
                    call propagate_2(atoms)

                    if (any(mod(istep, simparams%output_interval) == 0)) call output(atoms, itraj, istep)

                end do
            end do

        case default
            call abort

    end select


end program md_tian
