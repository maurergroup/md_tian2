
module geometry_opt

    use constants
    use universe_mod
    use run_config
    use force

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
            stop err // "unknown optimization algorithm"

        end select

    end subroutine optimize_geometry



    subroutine fire_geometry_optimization(atoms)

        type(universe), intent(inout) :: atoms

        ! Parameters as listed in E. Blitzek, PRL, 97, 170201 (2006)
        integer,  parameter :: n_min       = 5
        real(dp), parameter :: f_inc       = 1.1_dp
        real(dp), parameter :: f_dec       = 0.5_dp
        real(dp), parameter :: alpha_start = 0.1_dp
        real(dp), parameter :: f_alpha     = 0.99_dp

        real(dp) :: t_max, alpha
        !real(dp)

        ! From a typical MD simulation time step t_MD one can obtain an initial
        ! rough estimate of t_max = 10*t_MD
        t_max = 10*simparams%step

        ! All the velocities should be on the same scale, which for heteronuclear
        ! systems can be roughly achieved by setting all the atom masses equal
        atoms%m = amu2mass






    end subroutine fire_geometry_optimization



end module geometry_opt

