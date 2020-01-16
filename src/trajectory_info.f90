!############################################################################
! This routine is part of
! md_tian2 (Molecular Dynamics Tian Xia 2)
! (c) 2014-2020 Dan J. Auerbach, Svenja M. Janke, Marvin Kammler,
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

module trajectory_info

    use run_config, only : simparams
    use universe_mod
    use constants
    use rpmd

    implicit none

    logical, private :: was_moving_down = .true.

    integer  :: nturning_points
    real(dp) :: lowest_z(3), closest_approach

contains

    subroutine collect_trajectory_characteristics(atoms, istep)

        type(universe), intent(in) :: atoms
        integer,        intent(in) :: istep

        ! reset at beginning of new trajctory
        if (istep == 1) then
            closest_approach = default_real
            lowest_z         = default_real
            nturning_points  = 0
            was_moving_down  = .true.
        end if

        call check_for_bounce(atoms)
        call check_for_lowest_position(atoms)
        call check_for_closest_approach(atoms)

    end subroutine collect_trajectory_characteristics




    subroutine check_for_closest_approach(atoms)

        type(universe), intent(in) :: atoms

        real(dp) :: vec(3), cents(3, atoms%natoms), r
        integer  :: i, j

        cents = calc_centroid_positions(atoms)

        do i = 2, atoms%natoms
            vec = cents(:,1)-cents(:,i)      ! distance vector from a to b
            vec = matmul(atoms%isimbox, vec) ! transform to direct coordinates
            vec = vec - anint(vec)           ! imaging
            vec = matmul(atoms%simbox, vec)  ! back to cartesian coordinates
            r =  sqrt(sum(vec*vec))          ! distance
            if (r < closest_approach) closest_approach = r
        end do

    end subroutine check_for_closest_approach




    subroutine check_for_lowest_position(atoms)

        type(universe), intent(in) :: atoms
        real(dp) :: proj_pos(3)

        proj_pos = sum(atoms%r(:,:,1), dim=2)/atoms%nbeads

        if (proj_pos(3) < lowest_z(3)) lowest_z = proj_pos

    end subroutine check_for_lowest_position




    subroutine check_for_bounce(atoms)

        type(universe), intent(in) :: atoms
        logical :: now_moving_down

        ! calculate projectiles centroid z position
        now_moving_down = (sum(atoms%v(3,:,1))/atoms%nbeads < 0)

        if (now_moving_down .neqv. was_moving_down) then
            nturning_points = nturning_points + 1
            was_moving_down = now_moving_down
        end if

    end subroutine check_for_bounce

end module trajectory_info
