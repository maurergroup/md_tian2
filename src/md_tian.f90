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


    use md_init
    use universe_mod
    use pes_lj_mod
    use md_algo
    use rpmd

    implicit none

    integer :: itraj, istep
    type(universe) :: atoms


    call simbox_init(atoms)
    !call output_run_details()


    do itraj = simparams%start, simparams%start+simparams%ntrajs

        do istep = 1, simparams%nsteps

            call propagate_1(atoms)
            if (atoms%nbeads > 1) call ring_polymer_step(atoms)
            !call calc_forces(atoms)
            call propagate_2(atoms)

        end do
    end do


end program md_tian
