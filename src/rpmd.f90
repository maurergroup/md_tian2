module rpmd

    use universe_mod
    use run_config


    implicit none

    real(dp), allocatable :: cjk(:,:)

contains

    subroutine build_cjk(atoms)
        type(universe), intent(in) :: atoms

        integer :: j, k

        allocate(cjk(atoms%nbeads, atoms%nbeads))

        do j = 1, atoms%nbeads
            do k = 0, atoms%nbeads-1
                if (k.eq.0) then
                    cjk(j,k+1) = sqrt(1.0_dp/atoms%nbeads)
                else if (1 <= k .and. k <= atoms%nbeads/2 - 1) then
                    cjk(j,k+1) = sqrt(2.0_dp/atoms%nbeads) * cos(2.0_dp*PI*j*k/atoms%nbeads)
                else if (k .eq. atoms%nbeads/2) then
                    cjk(j,k+1) = sqrt(1.0_dp/atoms%nbeads)*(-1.0_dp)**j
                else if (atoms%nbeads/2+1 <= k .and. k <= atoms%nbeads-1) then
                    cjk(j,k+1) = sqrt(2.0_dp/atoms%nbeads) * sin(2.0_dp*PI*j*k/atoms%nbeads)
                else
                    stop "Error in build_cjk()"
                end if
            end do
        end do

    end subroutine build_cjk


    subroutine ring_polymer_step(atoms)

        type(universe), intent(inout) :: atoms

    !            real(8), intent(in)                                      :: betaN, dt
    !    integer, intent(in)                                      :: nParticles, nBeads
    !    real(8), dimension(3, nParticles, nBeads), intent(inout) :: p, q
    !    real(8), dimension(3, nParticles, nBeads)                :: newP, newQ
    !    real(8), dimension(nParticles), intent(in)               :: mass
    !    real(8), dimension(nBeads, nBeads)                       :: cjk
    !
    !    real(8) :: poly(4,nBeads)
    !    real(8), dimension(3) :: p_new
    !    real(8) :: twown, wk, wt, wm, cos_wt, sin_wt
    !    real(8), parameter :: PI = acos(-1.0d0)
    !    real(8) :: PI_N

    if (.not. allocated(cjk)) call build_cjk(atoms)


    end subroutine ring_polymer_step


end module rpmd
