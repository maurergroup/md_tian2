
module geometry_opt

    use constants
    use universe_mod
    use run_config
    use force
    use md_algo
    use rpmd

    implicit none

contains

    subroutine optimize_geometry(atoms, method)

        type(universe), intent(inout) :: atoms
        integer, intent(in)           :: method
        character(len=*), parameter   :: err = "Error in optimize_geometry(): "

        select case (method)

        case (geometry_opt_fire)
            call fire_geometry_optimization(atoms)

        case default
            print *, err // "unknown optimization algorithm"
            stop
        end select

    end subroutine optimize_geometry



    subroutine fire_geometry_optimization(atoms)

        use output_mod

        type(universe), intent(inout) :: atoms

        ! Parameters as listed in E. Blitzek, PRL, 97, 170201 (2006)
        integer,  parameter :: n_min         = 5
        integer,  parameter :: max_opt_steps = 10000
        real(dp), parameter :: f_inc         = 1.1_dp
        real(dp), parameter :: f_dec         = 0.5_dp
        real(dp), parameter :: alpha_start   = 0.1_dp
        real(dp), parameter :: f_alpha       = 0.99_dp
        real(dp), parameter :: force_limit   = 1e-5_dp

        real(dp) :: t_max, alpha, power
        real(dp) :: f_norm, v_norm, f_unit(3)
        real(dp) :: ring_polymer_forces(3, atoms%nbeads, atoms%natoms)
        integer  :: i, b, steps_since_power_was_negative, opt_steps

        ! From a typical MD simulation time step t_MD one can obtain an initial
        ! rough estimate of t_max = 10*t_MD
        t_max = 10*simparams%step

        ! All the velocities should be on the same scale, which for heteronuclear
        ! systems can be roughly achieved by setting all the atom masses equal
        atoms%m = amu2mass

        ! set initial variables
        alpha = alpha_start
        steps_since_power_was_negative = 0
        opt_steps = 0
        ring_polymer_forces = 0.0_dp
        power = 100.0_dp


        print '(a2, a17, 2a18)', "# ", "optimization step", "power/eV*fs", "Epot/eV"

        do while (.true.)

            ! MD
            call propagate_1(atoms)
            call calc_force(atoms, energy_and_force)
            if (atoms%nbeads > 1) call do_ring_polymer_step_with_forces(atoms, ring_polymer_forces)

            call propagate_2(atoms)

            ! exit normally
            if (opt_steps >= max_opt_steps     .or. &
                maxval(atoms%f) <= force_limit .or. &
                abs(power) < 1e-7) exit

            ! exit due to error
            if (abs(power) > 1e4) then
                atoms%epot = default_real
                exit
            end if

            ! F1
            atoms%f = atoms%f + ring_polymer_forces
            power = sum(atoms%f * atoms%v)
            print '(i19, 2f18.11)', opt_steps, power, atoms%epot

            ! F2
            do i = 1, atoms%natoms
                do b = 1, atoms%nbeads
                    f_norm = sum(sqrt(atoms%f(:,b,i)*atoms%f(:,b,i)))
                    v_norm = sum(sqrt(atoms%v(:,b,i)*atoms%v(:,b,i)))
                    f_unit = atoms%f(:,b,i)/f_norm

                    where (.not. atoms%is_fixed(:,b,i))
                        atoms%v(:,b,i) = (1-alpha)*atoms%v(:,b,i) + alpha*f_unit*v_norm
                    elsewhere
                        atoms%v(:,b,i) = 0
                    end where

                end do
            end do

            ! F3
            if (power > 0 .and. steps_since_power_was_negative > n_min) then
                simparams%step = min(simparams%step*f_inc, t_max)
                alpha = alpha*f_alpha
            end if

            ! F4
            if (power <= 0) then
                simparams%step = simparams%step*f_dec
                atoms%v = 0
                alpha = alpha_start
            end if

            ! increment and/or reset counters
            if (power > 0) then
                steps_since_power_was_negative = steps_since_power_was_negative +1
            else
                steps_since_power_was_negative = 0
            end if
            opt_steps = opt_steps + 1

        end do

    end subroutine fire_geometry_optimization



end module geometry_opt

