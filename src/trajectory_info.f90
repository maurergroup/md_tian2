module trajectory_info

    use run_config, only : simparams
    use universe_mod
    use constants
    use rpmd
    use useful_things
    use open_file

    implicit none

    integer  :: nearest_surface_idx, bounces, col_partner_change
    integer  :: last_surface_idx, rFcut, vFcut, r_dot_v_ptm
    real(dp) :: lowest_z(3), closest_approach, interaction_time
    real(dp) :: r_dot_v(2), r_dot_F(2), v_dot_F(2)      ! needed for collision detection
    logical  :: is_adsorbed
    real(dp), allocatable :: centroid_positions(:,:)

contains

    subroutine collect_trajectory_characteristics(atoms, itraj, istep)

        type(universe), intent(in) :: atoms
        integer,        intent(in) :: itraj, istep

        ! reset at beginning of new trajctory
        if (istep == 1) then
            closest_approach = default_real
            lowest_z         = default_real
            interaction_time = 0
            bounces          = 0
            r_dot_v          = default_real
            r_dot_F          = default_real
            v_dot_F          = default_real
            col_partner_change = 1
            rFcut            = 0
            vFcut            = 0
            r_dot_v_ptm      = 0 ! ptm = plus to minus
        end if

        if (.not. allocated(centroid_positions)) allocate(centroid_positions(3,atoms%natoms))
        centroid_positions = calc_centroid_positions(atoms)
        call check_for_lowest_position(atoms)
        call check_for_closest_approach(atoms)
        call check_for_proj_surf_interaction(atoms)
        call check_for_bounce(atoms, itraj, istep)

    end subroutine collect_trajectory_characteristics




    subroutine check_for_closest_approach(atoms)

        type(universe), intent(in) :: atoms

        real(dp) :: vec(3), r
        real(dp) :: cl_appr_in_this_step
        integer  :: i

        cl_appr_in_this_step = default_real

        do i = 2, atoms%natoms
            ! closest approach overall in trajectory
            vec = centroid_positions(:,1)-centroid_positions(:,i)      ! distance vector from a to b
            vec = matmul(atoms%isimbox, vec) ! transform to direct coordinates
            vec = vec - anint(vec)           ! imaging
            vec = matmul(atoms%simbox, vec)  ! back to cartesian coordinates
            r =  sqrt(sum(vec*vec))          ! distance
            if (r < closest_approach) closest_approach = r
            if (r < cl_appr_in_this_step) then
                cl_appr_in_this_step = r
                nearest_surface_idx = i      ! needed for collision detection
            end if

        end do

        ! check for adsorption
        if (.not. is_adsorbed .and. cl_appr_in_this_step < simparams%adsorption_start) is_adsorbed = .true.
        if (      is_adsorbed .and. cl_appr_in_this_step > simparams%adsorption_end  ) is_adsorbed = .false.

    end subroutine check_for_closest_approach




    subroutine check_for_lowest_position(atoms)

        type(universe), intent(in) :: atoms
        real(dp) :: proj_pos(3)

        proj_pos = sum(atoms%r(:,:,1), dim=2)/atoms%nbeads

        if (proj_pos(3) < lowest_z(3)) lowest_z = proj_pos

    end subroutine check_for_lowest_position



    subroutine check_for_proj_surf_interaction(atoms)

        type(universe), intent(in) :: atoms

        if (any(atoms%f(:,:,1) /= 0)) then
            interaction_time = interaction_time + simparams%step
        end if

    end subroutine check_for_proj_surf_interaction



    subroutine check_for_bounce(atoms, itraj, istep)

        type(universe), intent(in) :: atoms
        integer       , intent(in) :: itraj, istep

        real(dp) :: r(3), v(3), F(3)
        integer  :: last_idx, this_idx
        character(len=max_string_length) :: fname
        character(len=8)                 :: traj_id


        ! distance vector from projectile to nearst surface atom (r)
        r = centroid_positions(:,nearest_surface_idx)-centroid_positions(:,1)
        r = matmul(atoms%isimbox, r)
        r = r - anint(r)
        r = matmul(atoms%simbox, r)
        call normalize(r)

        ! projectile velocity (v)
        v = calc_centroid_velocity_one(atoms, 1)
        call normalize(v)

        ! projectile force (F)
        F = calc_centroid_force_one(atoms, 1)
        call normalize(F)

        ! calculate dot products, Fortran is one-based
        last_idx          = mod(istep + 1, 2) + 1   ! idx at last time step
        this_idx          = mod(istep + 0, 2) + 1   ! idx at this time step
        r_dot_v(this_idx) = sum(r*v)
        r_dot_F(this_idx) = sum(r*F)
        v_dot_F(this_idx) = sum(v*F)

        ! collision axioms:
        ! collision is only possible with nearest surface atom
        ! an approach must turn into a departure

        ! detection:
        ! split trajectory into segments in which a bounce could have occured
        !   that is each time the nearest surface atoms changes
        ! evaluate the segments including knowledge about prior events
        ! in each collision r.v must have changed sign from + to -
        ! # bounces in this segment is the smallest number from the set
        !   r.v + to - sign changes (it really cannot be more than that)
        !   r.v intersects r.F
        !   r.v intersects v.F

!        write(traj_id,'(i8.8)') itraj
!        fname = 'conf/bounce_'//traj_id//'.dat'
!        if (.not. dir_exists('conf')) call system('mkdir conf')
        if (istep == 1) then
!            call open_for_write(34, fname)
            last_surface_idx = nearest_surface_idx
!        else
!            call open_for_append(34, fname)
        end if

        ! reset stats if nearest surface atom changed
        if (last_surface_idx /= nearest_surface_idx) then
            call finalize_bounces()
            col_partner_change = -1 * col_partner_change
            rFcut = 0
            vFcut = 0
            r_dot_v_ptm = 0     ! r.v + to -
        ! avoid crossings when entering/leaving interaction region
        else if (r_dot_F(this_idx) /= 0) then
            ! r.v cut r.F while decreasing
            if (r_dot_v(last_idx) > r_dot_F(last_idx) .and. r_dot_v(this_idx) < r_dot_F(this_idx)) then
                rFcut = rFcut + 1
!                print *, "traj", itraj, "rF cut dec at", istep, "rF = ", rFcut
            ! r.v cut r.F while increasing
            else if (r_dot_v(last_idx) < r_dot_F(last_idx) .and. r_dot_v(this_idx) > r_dot_F(this_idx)) then
                rFcut = rFcut + 1
!                print *, "traj", itraj, "rF cut inc at", istep, "rF = ", rFcut
            end if
            ! r.v cut v.F while decreasing
            if (r_dot_v(last_idx) > v_dot_F(last_idx) .and. r_dot_v(this_idx) < v_dot_F(this_idx)) then
                vFcut = vFcut + 1
!                print *, "traj", itraj, "vF cut dec at", istep, "vF = ", vFcut
            ! r.v cut v.F while increasing
            else if (r_dot_v(last_idx) < v_dot_F(last_idx) .and. r_dot_v(this_idx) > v_dot_F(this_idx)) then
                vFcut = vFcut + 1
!                print *, "traj", itraj, "vF cut inc at", istep, "vF = ", vFcut
            end if
            ! projectile turned its approach into departure
            if (r_dot_v(last_idx) > 0 .and. r_dot_v(this_idx) < 0) then
                r_dot_v_ptm = r_dot_v_ptm + 1
!                print *, "traj", itraj, "r.v ptm", istep, "r.v_ptm = ", r_dot_v_ptm
            end if
        end if

!        write(34, '(i, 3f14.7, i, f14.7)') istep, r_dot_v(this_idx), r_dot_F(this_idx), v_dot_F(this_idx), &
!                                    col_partner_change
!        close(34)

        last_surface_idx = nearest_surface_idx

    end subroutine check_for_bounce



    subroutine finalize_bounces()

        bounces = bounces + min(r_dot_v_ptm, rFcut, vFcut)

    end subroutine finalize_bounces


end module trajectory_info
