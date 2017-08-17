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

    implicit none

    integer :: itraj, istep, i, j, b
    type(universe) :: atoms


    real(dp) :: tmp, vcm(3), nom(3), denom(3), tinit, tinter


    call simbox_init(atoms)
    !call output_run_details()


    do itraj = simparams%start, simparams%start+simparams%ntrajs-1

        call calc_force(atoms)

!        print *, atoms%r
!        print *, ""
!        print *, atoms%f
!        print *, ""
!        print *, atoms%a
!                print *, ""


        do istep = 1, simparams%nsteps

            call propagate_1(atoms)
!            print *, atoms%r

            if (atoms%nbeads > 1) call ring_polymer_step(atoms)
            call calc_force(atoms)
            call propagate_2(atoms)

!            print *, atoms%r
!            print *, ""
!            print *, atoms%f
!        print *, ""
!        print *, atoms%a

            if (mod(istep, simparams%output(2)) == 0) call output(atoms, itraj, istep, simparams%output(1))



        end do
    end do


end program md_tian
