module pes_rebo_mod

    use universe_mod
    use useful_things

    implicit none

    include "mkl_lapack.fi"

    ! Rebo PES object
    type rebo_pes

        private
        real(dp), dimension(:,:),   allocatable :: Dmin, Dmax, Dmaxp, Q, A, alpha, rho
        real(dp), dimension(:,:,:), allocatable :: B, beta

    end type rebo_pes

    type(rebo_pes) :: pes_rebo

    ! Fitting
    logical, dimension(:,:),   allocatable, private :: fit_Dmin, fit_Dmax, fit_Dmaxp
    logical, dimension(:,:),   allocatable, private :: fit_Q, fit_A, fit_alpha, fit_rho
    logical, dimension(:,:,:), allocatable, private :: fit_B, fit_beta


    ! Declarations
    real(dp), parameter :: TETRAHEDRON_ANGLE = -1.0_dp/3

    real(dp), parameter :: CSF_LOW  = 2.0_dp ! conjugation screening function F(x_ik), page 792, eq. 16
    real(dp), parameter :: CSF_HIGH = 3.0_dp
    real(dp), parameter :: NSF_LOW  = 3.2_dp ! neighbor screening function Q_i, page 791, eq. 12
    real(dp), parameter :: NSF_HIGH = 3.7_dp

    ! Initialization flag
    logical, public :: REBO_INIT_SUCCESS = .false.

    integer, private :: i, j, k


    ! P_CC (2D)
    integer, parameter, private                         :: MIN_I_P_CC = 0
    integer, parameter, private                         :: MAX_I_P_CC = 4
    integer, parameter, private                         :: MIN_J_P_CC = 0
    integer, parameter, private                         :: MAX_J_P_CC = 4
    real(dp), dimension(MIN_I_P_CC:MAX_I_P_CC,&
        MIN_J_P_CC:MAX_J_P_CC), private  :: P_CC, di_P_CC, dj_P_CC, cd_P_CC

    ! P_CH (2D)
    integer, parameter, private                         :: MIN_I_P_CH = 0
    integer, parameter, private                         :: MAX_I_P_CH = 4
    integer, parameter, private                         :: MIN_J_P_CH = 0
    integer, parameter, private                         :: MAX_J_P_CH = 4
    real(dp), dimension(MIN_I_P_CH:MAX_I_P_CH,&
        MIN_J_P_CH:MAX_J_P_CH), private  :: P_CH, di_P_CH, dj_P_CH, cd_P_CH

    ! F_CC (3D)
    integer, parameter, private                         :: MIN_I_F_CC = 0
    integer, parameter, private                         :: MAX_I_F_CC = 4
    integer, parameter, private                         :: MIN_J_F_CC = 0
    integer, parameter, private                         :: MAX_J_F_CC = 4
    integer, parameter, private                         :: MIN_K_F_CC = 0
    integer, parameter, private                         :: MAX_K_F_CC = 9
    real(dp), dimension(MIN_I_F_CC:MAX_I_F_CC,&
        MIN_J_F_CC:MAX_J_F_CC,&
        MIN_K_F_CC:MAX_K_F_CC), private  :: F_CC, di_F_CC, dj_F_CC, dk_F_CC,&
        dij_F_CC, dik_F_CC, djk_F_CC, cd_F_CC

    ! T_CC (3D)
    integer, parameter, private                         :: MIN_I_T_CC = 1
    integer, parameter, private                         :: MAX_I_T_CC = 3
    integer, parameter, private                         :: MIN_J_T_CC = 1
    integer, parameter, private                         :: MAX_J_T_CC = 3
    integer, parameter, private                         :: MIN_K_T_CC = 1
    integer, parameter, private                         :: MAX_K_T_CC = 9
    real(dp), dimension(MIN_I_T_CC:MAX_I_T_CC,&
        MIN_J_T_CC:MAX_J_T_CC,&
        MIN_K_T_CC:MAX_K_T_CC), private  :: T_CC, di_T_CC, dj_T_CC, dk_T_CC,&
        dij_T_CC, dik_T_CC, djk_T_CC, cd_T_CC

    ! F_CH (3D)
    integer, parameter, private                         :: MIN_I_F_CH = 0
    integer, parameter, private                         :: MAX_I_F_CH = 4
    integer, parameter, private                         :: MIN_J_F_CH = 0
    integer, parameter, private                         :: MAX_J_F_CH = 4
    integer, parameter, private                         :: MIN_K_F_CH = 0
    integer, parameter, private                         :: MAX_K_F_CH = 9
    real(dp), dimension(MIN_I_F_CH:MAX_I_F_CH,&
        MIN_J_F_CH:MAX_J_F_CH,&
        MIN_K_F_CH:MAX_K_F_CH), private  :: F_CH, di_F_CH, dj_F_CH, dk_F_CH,&
        dij_F_CH, dik_F_CH, djk_F_CH, cd_F_CH

    ! F_HH (3D)
    integer, parameter, private                         :: MIN_I_F_HH = 0
    integer, parameter, private                         :: MAX_I_F_HH = 2
    integer, parameter, private                         :: MIN_J_F_HH = 0
    integer, parameter, private                         :: MAX_J_F_HH = 2
    integer, parameter, private                         :: MIN_K_F_HH = 0
    integer, parameter, private                         :: MAX_K_F_HH = 2
    real(dp), dimension(MIN_I_F_HH:MAX_I_F_HH,&
        MIN_J_F_HH:MAX_J_F_HH,&
        MIN_K_F_HH:MAX_K_F_HH), private  :: F_HH, di_F_HH, dj_F_HH, dk_F_HH,&
        dij_F_HH, dik_F_HH, djk_F_HH, cd_F_HH


    ! Arrays to hold the spline coefficients
    real(dp), dimension(MIN_I_P_CC:MAX_I_P_CC-1,&
        MIN_J_P_CC:MAX_J_P_CC-1,16), private :: P_CC_params
    real(dp), dimension(MIN_I_P_CH:MAX_I_P_CH-1,&
        MIN_J_P_CH:MAX_J_P_CH-1,16), private :: P_CH_params
    real(dp), dimension(MIN_I_F_CC:MAX_I_F_CC-1,&
        MIN_J_F_CC:MAX_J_F_CC-1,&
        MIN_K_F_CC:MAX_K_F_CC-1,64), private :: F_CC_params
    real(dp), dimension(MIN_I_T_CC:MAX_I_T_CC-1,&
        MIN_J_T_CC:MAX_J_T_CC-1,&
        MIN_K_T_CC:MAX_K_T_CC-1,64), private :: T_CC_params
    real(dp), dimension(MIN_I_F_CH:MAX_I_F_CH-1,&
        MIN_J_F_CH:MAX_J_F_CH-1,&
        MIN_K_F_CH:MAX_K_F_CH-1,64), private :: F_CH_params
    real(dp), dimension(MIN_I_F_HH:MAX_I_F_HH-1,&
        MIN_J_F_HH:MAX_J_F_HH-1,&
        MIN_K_F_HH:MAX_K_F_HH-1,64), private :: F_HH_params

    ! gSpline
    real(dp), dimension(5), private   :: gCdom
    real(dp), dimension(4), private   :: gHdom
    real(dp), dimension(4,6), private :: gC1, gC2
    real(dp), dimension(3,6), private :: gH



    ! Subroutine and function accessabililty
    private :: init_spline2d, init_spline3d, gSpline, Sp5th, rebosi2d, rebosi3d
    public  :: rebosi_initialize, cufu

contains

    subroutine read_rebo(atoms, inp_unit)

        use run_config, only : simparams

        type(universe), intent(inout) :: atoms
        integer, intent(in) :: inp_unit

        integer :: nwords, ios = 0, i
        character(len=max_string_length) :: buffer
        character(len=max_string_length) :: words(100)
        integer  :: idx1, idx2, ntypes
        real(dp) :: temp3(3)
        character(len=*), parameter :: err = "Error in read_rebo: "

        call rebosi_initialize()

        ntypes = simparams%nprojectiles+simparams%nlattices

        if (.not. allocated(pes_rebo%Q)) then
            allocate(pes_rebo%Dmin(ntypes, ntypes),  &
                pes_rebo%Dmax(ntypes, ntypes),  &
                pes_rebo%Dmaxp(ntypes, ntypes), &
                pes_rebo%Q(ntypes, ntypes),     &
                pes_rebo%A(ntypes, ntypes),     &
                pes_rebo%alpha(ntypes, ntypes), &
                pes_rebo%rho(ntypes, ntypes),   &
                pes_rebo%B(3, ntypes, ntypes),  &
                pes_rebo%beta(3, ntypes, ntypes), &

                fit_Dmin(ntypes, ntypes),  &
                fit_Dmax(ntypes, ntypes),  &
                fit_Dmaxp(ntypes, ntypes), &
                fit_Q(ntypes, ntypes),     &
                fit_A(ntypes, ntypes),     &
                fit_alpha(ntypes, ntypes), &
                fit_rho(ntypes, ntypes),   &
                fit_B(3, ntypes, ntypes),  &
                fit_beta(3, ntypes, ntypes))

            pes_rebo%Dmin  = default_real
            pes_rebo%Dmax  = default_real
            pes_rebo%Dmaxp = default_real
            pes_rebo%Q     = default_real
            pes_rebo%A     = default_real
            pes_rebo%alpha = default_real
            pes_rebo%rho   = default_real
            pes_rebo%B     = default_real
            pes_rebo%beta  = default_real

            fit_Dmin  = .False.
            fit_Dmax  = .False.
            fit_Dmaxp = .False.
            fit_Q     = .False.
            fit_A     = .False.
            fit_alpha = .False.
            fit_rho   = .False.
            fit_B     = .False.
            fit_beta  = .False.

        end if

        ! line should read something like "H   H   proj    proj"
        read(inp_unit, '(A)', iostat=ios) buffer
        call split_string(buffer, words, nwords)

        if (nwords /= 4) stop err // "need four entries in interaction-defining lines"

        if (words(3) == "proj" .and. words(4) == "proj" .or. &
            words(3) == "proj" .and. words(4) == "latt" .or. &
            words(3) == "latt" .and. words(4) == "proj" .or. &
            words(3) == "latt" .and. words(4) == "latt") then

            idx1 = get_idx_from_name(atoms, words(1), is_proj=(words(3)=="proj"))
            idx2 = get_idx_from_name(atoms, words(2), is_proj=(words(4)=="proj"))

            if (atoms%pes(idx1,idx2) /= default_int) then
                print *, err // "pes already defined for atoms", words(1), words(3), words(2), words(4)
                stop
            end if

        else
            print *, err // "interaction must be defined via 'proj' and 'latt' keywords"
            stop
        end if

        ! set the pes type in the atoms object
        atoms%pes(idx1,idx2) = pes_id_rebo
        atoms%pes(idx2,idx1) = pes_id_rebo

        do
            read(inp_unit, '(A)', iostat=ios) buffer
            call split_string(buffer, words, nwords)

            ! pes block terminated, set shift
            if (nwords == 0 .or. ios /= 0) then
                exit

            ! something went wrong
            else if (nwords /= 2) then
                if (nwords /= 3 .or. words(3) /= "fit") then
                    print *,  err // "Error in the PES file: PES parameters must &
                        consist of key value pairs. A parameter block must be terminated by a blank line."
                    stop
                end if
            end if

            call lower_case(words(1))

            select case (words(1))

                case ('dmin')
                    read(words(2), *) pes_rebo%Dmin(idx1, idx2)
                    read(words(2), *) pes_rebo%Dmin(idx2, idx1)
                    call check_and_set_fit(words(3), idx1, idx2, fit_Dmin)

                case ('dmax')
                    read(words(2), *) pes_rebo%Dmax(idx1, idx2)
                    read(words(2), *) pes_rebo%Dmax(idx2, idx1)
                    call check_and_set_fit(words(3), idx1, idx2, fit_Dmax)

                case ('dmaxp')
                    read(words(2), *) pes_rebo%Dmaxp(idx1, idx2)
                    read(words(2), *) pes_rebo%Dmaxp(idx2, idx1)
                    call check_and_set_fit(words(3), idx1, idx2, fit_Dmaxp)

                case ('q')
                    read(words(2), *) pes_rebo%Q(idx1, idx2)
                    read(words(2), *) pes_rebo%Q(idx2, idx1)
                    call check_and_set_fit(words(3), idx1, idx2, fit_Q)

                case ('a')
                    read(words(2), *) pes_rebo%A(idx1, idx2)
                    read(words(2), *) pes_rebo%A(idx2, idx1)
                    call check_and_set_fit(words(3), idx1, idx2, fit_A)

                case ('alpha')
                    read(words(2), *) pes_rebo%alpha(idx1, idx2)
                    read(words(2), *) pes_rebo%alpha(idx2, idx1)
                    call check_and_set_fit(words(3), idx1, idx2, fit_alpha)

                case ('rho')
                    read(words(2), *) pes_rebo%rho(idx1, idx2)
                    read(words(2), *) pes_rebo%rho(idx2, idx1)
                    call check_and_set_fit(words(3), idx1, idx2, fit_rho)

                case ('b1')
                    read(words(2), *) pes_rebo%B(1, idx1, idx2)
                    read(words(2), *) pes_rebo%B(1, idx2, idx1)
                    call check_and_set_fit(words(3), idx1, idx2, fit_B(1,:,:))

                case ('b2')
                    read(words(2), *) pes_rebo%B(2, idx1, idx2)
                    read(words(2), *) pes_rebo%B(2, idx2, idx1)
                    call check_and_set_fit(words(3), idx1, idx2, fit_B(2,:,:))

                case ('b3')
                    read(words(2), *) pes_rebo%B(3, idx1, idx2)
                    read(words(2), *) pes_rebo%B(3, idx2, idx1)
                    call check_and_set_fit(words(3), idx1, idx2, fit_B(3,:,:))

                case ('beta1')
                    read(words(2), *) pes_rebo%beta(1, idx1, idx2)
                    read(words(2), *) pes_rebo%beta(1, idx2, idx1)
                    call check_and_set_fit(words(3), idx1, idx2, fit_beta(1,:,:))

                case ('beta2')
                    read(words(2), *) pes_rebo%beta(2, idx1, idx2)
                    read(words(2), *) pes_rebo%beta(2, idx2, idx1)
                    call check_and_set_fit(words(3), idx1, idx2, fit_beta(2,:,:))

                case ('beta3')
                    read(words(2), *) pes_rebo%beta(3, idx1, idx2)
                    read(words(2), *) pes_rebo%beta(3, idx2, idx1)
                    call check_and_set_fit(words(3), idx1, idx2, fit_beta(3,:,:))

                case default
                    print *, "Error in the PES file: unknown REBO parameter", words(1)
                    stop

            end select
        end do

    end subroutine read_rebo


    ! returns parameters set to fit
    subroutine get_fit_params_rebo(fit_params)

        real(dp), allocatable, intent(out) :: fit_params(:)
        integer :: i, j, k, n, nfit, ntypes

        ntypes = size(fit_Q, dim=1)

        ! count the number of fit parameters
        nfit = 0
        do i = 1, ntypes
            do j = i, ntypes
                nfit = nfit + count([fit_Dmin(j,i), fit_Dmax(j,i), fit_rho(j,i), &
                    fit_Dmaxp(j,i), fit_Q(j,i), fit_A(j,i), fit_alpha(j,i)])
                do k = 1, 3
                    nfit = nfit + count([fit_B(k,j,i), fit_beta(k,j,i)])
                end do
            end do
        end do

        allocate(fit_params(nfit))

        n = 1
        do i = 1, ntypes
            do j = i, ntypes
                if (fit_A(j,i)) then
                    fit_params(n) = pes_rebo%A(j,i)
                    n = n + 1
                end if
                if (fit_alpha(j,i)) then
                    fit_params(n) = pes_rebo%alpha(j,i)
                    n = n + 1
                end if
                do k = 1, 3
                    if (fit_B(k,j,i)) then
                        fit_params(n) = pes_rebo%B(k,j,i)
                        n = n + 1
                    end if
                    if (fit_beta(k,j,i)) then
                        fit_params(n) = pes_rebo%beta(k,j,i)
                        n = n + 1
                    end if
                end do
                if (fit_Dmax(j,i)) then
                    fit_params(n) = pes_rebo%Dmax(j,i)
                    n = n + 1
                end if
                if (fit_Dmaxp(j,i)) then
                    fit_params(n) = pes_rebo%Dmaxp(j,i)
                    n = n + 1
                end if
                if (fit_Dmin(j,i)) then
                    fit_params(n) = pes_rebo%Dmin(j,i)
                    n = n + 1
                end if
                if (fit_Q(j,i)) then
                    fit_params(n) = pes_rebo%Q(j,i)
                    n = n + 1
                end if
                if (fit_rho(j,i)) then
                    fit_params(n) = pes_rebo%rho(j,i)
                    n = n + 1
                end if
            end do
        end do

    end subroutine get_fit_params_rebo




    subroutine set_fit_params_rebo(x, pos)

        real(dp), intent(in)   :: x(:)
        integer, intent(inout) :: pos

        integer :: i, j, k, ntypes

        ntypes = size(fit_Q, dim=1)

        do i = 1, ntypes
            do j = i, ntypes
                if (fit_A(j,i)) then
                    pes_rebo%A(j,i) = x(pos)
                    pes_rebo%A(i,j) = x(pos)
                    pos = pos + 1
                end if
                if (fit_alpha(j,i)) then
                    pes_rebo%alpha(j,i) = x(pos)
                    pes_rebo%alpha(i,j) = x(pos)
                    pos = pos + 1
                end if
                do k = 1, 3
                    if (fit_B(k,j,i)) then
                        pes_rebo%B(k,j,i) = x(pos)
                        pes_rebo%B(k,i,j) = x(pos)
                        pos = pos + 1
                    end if
                    if (fit_beta(k,j,i)) then
                        pes_rebo%beta(k,j,i) = x(pos)
                        pes_rebo%beta(k,i,j) = x(pos)
                        pos = pos + 1
                    end if
                end do
                if (fit_Dmax(j,i)) then
                    pes_rebo%Dmax(j,i) = x(pos)
                    pes_rebo%Dmax(i,j) = x(pos)
                    pos = pos + 1
                end if
                if (fit_Dmaxp(j,i)) then
                    pes_rebo%Dmaxp(j,i) = x(pos)
                    pes_rebo%Dmaxp(i,j) = x(pos)
                    pos = pos + 1
                end if
                if (fit_Dmin(j,i)) then
                    pes_rebo%Dmin(j,i) = x(pos)
                    pes_rebo%Dmin(i,j) = x(pos)
                    pos = pos + 1
                end if
                if (fit_Q(j,i)) then
                    pes_rebo%Q(j,i) = x(pos)
                    pes_rebo%Q(i,j) = x(pos)
                    pos = pos + 1
                end if
                if (fit_rho(j,i)) then
                    pes_rebo%rho(j,i) = x(pos)
                    pes_rebo%rho(i,j) = x(pos)
                    pos = pos + 1
                end if
            end do
        end do

    end subroutine set_fit_params_rebo



    subroutine to_string_rebo(out_unit, idx_i, idx_j)

        integer, intent(in) :: out_unit, idx_i, idx_j

        write(out_unit, '(a6, e18.11)') "Dmin ", pes_rebo%Dmin(idx_i, idx_j)
        write(out_unit, '(a6, e18.11)') "Dmax ", pes_rebo%Dmax(idx_i, idx_j)
        write(out_unit, '(a6, e18.11)') "Dmaxp", pes_rebo%Dmaxp(idx_i, idx_j)
        write(out_unit, '(a6, e18.11)') "B1   ", pes_rebo%B(1, idx_i, idx_j)
        write(out_unit, '(a6, e18.11)') "B2   ", pes_rebo%B(2, idx_i, idx_j)
        write(out_unit, '(a6, e18.11)') "B3   ", pes_rebo%B(3, idx_i, idx_j)
        write(out_unit, '(a6, e18.11)') "beta1", pes_rebo%beta(1, idx_i, idx_j)
        write(out_unit, '(a6, e18.11)') "beta2", pes_rebo%beta(2, idx_i, idx_j)
        write(out_unit, '(a6, e18.11)') "beta3", pes_rebo%beta(3, idx_i, idx_j)
        write(out_unit, '(a6, e18.11)') "Q    ", pes_rebo%Q(idx_i, idx_j)
        write(out_unit, '(a6, e18.11)') "A    ", pes_rebo%A(idx_i, idx_j)
        write(out_unit, '(a6, e18.11)') "alpha", pes_rebo%alpha(idx_i, idx_j)
        write(out_unit, '(a6, e18.11)') "rho  ", pes_rebo%rho(idx_i, idx_j)

    end subroutine to_string_rebo



    subroutine rebosi_initialize()

        ! Do not initialize twice
        if (REBO_INIT_SUCCESS) return

        !!! P_CC initialization !!!
        ! fill P_CC(i,j) array containing values at knots
        P_CC = 0.0_dp
        P_CC(0,2) = -0.0005_dp
        P_CC(0,3) = 0.016125364564267_dp
        P_CC(1,1) = -0.01096_dp
        P_CC(1,2) = 0.006326248241119_dp
        P_CC(2,0) = 0.0_dp !-0.027603_dp
        P_CC(2,1) = 0.003179530830731_dp


        ! fill dP_CC(i,j) array containing derivatives and cross derivatives at knots
        di_P_CC = 0.0_dp
        dj_P_CC = 0.0_dp
        cd_P_CC = 0.0_dp


        !!! P_CH initialization !!!
        ! fill P_CH(i,j) array containing values at knots
        P_CH = 0.0_dp
        P_CH(0,1) =  0.2093367328250380_dp
        P_CH(0,2) = -0.064449615432525_dp
        P_CH(0,3) = -0.303927546346162_dp
        P_CH(1,0) =  0.01_dp
        P_CH(1,1) = -0.1251234006287090_dp
        P_CH(1,2) = -0.298905245783_dp
        P_CH(2,0) = -0.1220421462782555_dp
        P_CH(2,1) = -0.3005291724067579_dp
        P_CH(3,0) = -0.307584705066_dp


        ! fill dP_CH(i,j) array containing derivatives and cross derivatives at knots
        di_P_CH = 0.0_dp
        dj_P_CH = 0.0_dp
        cd_P_CH = 0.0_dp


        !!! F_CC initialization !!!
        ! fill F_CC(i,j,k) array containing values at knots
        F_CC(0,0,3:9) = 0.0049586079_dp
        F_CC(1,0,1) = 0.021693495_dp
        F_CC(0,1,1) = 0.021693495_dp
        F_CC(1,0,2:9) = 0.0049586079_dp
        F_CC(0,1,2:9) = 0.0049586079_dp
        F_CC(1,1,1) = 0.05250_dp
        F_CC(1,1,2) = -0.002088750_dp
        F_CC(1,1,3:9) = -0.00804280_dp
        F_CC(2,0,1) = 0.024698831850_dp
        F_CC(0,2,1) = 0.024698831850_dp
        F_CC(2,0,2) = -0.00597133450_dp
        F_CC(0,2,2) = -0.00597133450_dp
        F_CC(2,0,3:9) = 0.0049586079_dp
        F_CC(0,2,3:9) = 0.0049586079_dp
        F_CC(2,1,1) = 0.00482478490_dp
        F_CC(1,2,1) = 0.00482478490_dp
        F_CC(2,1,2) = 0.0150_dp
        F_CC(1,2,2) = 0.0150_dp
        F_CC(2,1,3) = -0.010_dp
        F_CC(1,2,3) = -0.010_dp
        F_CC(2,1,4) = -0.01168893870_dp
        F_CC(1,2,4) = -0.01168893870_dp
        F_CC(2,1,5) = -0.013377877400_dp
        F_CC(1,2,5) = -0.013377877400_dp
        F_CC(2,1,6) = -0.015066816000_dp
        F_CC(1,2,6) = -0.015066816000_dp
        F_CC(2,1,7:9) = -0.015066816000_dp
        F_CC(1,2,7:9) = -0.015066816000_dp
        F_CC(2,2,1) = 0.0472247850_dp
        F_CC(2,2,2) = 0.0110_dp
        F_CC(2,2,3) = 0.0198529350_dp
        F_CC(2,2,4) = 0.01654411250_dp
        F_CC(2,2,5) = 0.013235290_dp
        F_CC(2,2,6) = 0.00992646749999_dp
        F_CC(2,2,7) = 0.006617644999_dp
        F_CC(2,2,8) = 0.00330882250_dp
        F_CC(3,0,1) = -0.05989946750_dp
        F_CC(0,3,1) = -0.05989946750_dp
        F_CC(3,0,2) = -0.05989946750_dp
        F_CC(0,3,2) = -0.05989946750_dp
        F_CC(3,0,3:9) = 0.0049586079_dp
        F_CC(0,3,3:9) = 0.0049586079_dp
        F_CC(3,1,2) = -0.0624183760_dp
        F_CC(1,3,2) = -0.0624183760_dp
        F_CC(3,1,3:9) = -0.0624183760_dp
        F_CC(1,3,3:9) = -0.0624183760_dp
        F_CC(3,2,1) = -0.02235469150_dp
        F_CC(2,3,1) = -0.02235469150_dp
        F_CC(3,2,2:9) = -0.02235469150_dp
        F_CC(2,3,2:9) = -0.02235469150_dp


        ! make top end of piCC flat instead of zero
        F_CC(4,0:3,1:9) = F_CC(3,0:3,1:9)

        ! also enforces some symmetry
        do i = 0, 3
            do j = i+1, 4
                F_CC(i,j,1:9) = F_CC(j,i,1:9)
            end do
        end do
        F_CC(4,4,1:9) = F_CC(3,4,1:9)


        ! fill di_F_CC(i,j,k) array containing derivatives at knots
        di_F_CC = 0.0_dp
        di_F_CC(2,1,1)   = -0.026250_dp
        di_F_CC(2,1,5:9) = -0.0271880_dp
        di_F_CC(1,3,2)   =  0.0187723882_dp
        di_F_CC(2,3,2:9) =  0.031209_dp

        ! fill dj_F_CC(i,j,k) array containing derivatives at knots
        dj_F_CC  = 0.0_dp
        dj_F_CC(1,2,1)   = -0.026250_dp
        dj_F_CC(1,2,5:9) = -0.0271880_dp
        dj_F_CC(3,1,2)   = 0.0187723882_dp
        dj_F_CC(3,2,2:9) =  0.031209_dp

        ! fill dk_F_CC(i,j,k) array containing derivatives at knots
        dk_F_CC = 0.0_dp
        dk_F_CC(1,1,2) = -0.0302715_dp
        dk_F_CC(2,1,4:5) = -0.0100220_dp
        dk_F_CC(1,2,4:5) = -0.0100220_dp
        dk_F_CC(2,2,4:8) = -0.0033090_dp

        ! fill rest with zeros
        dij_F_CC = 0.0_dp
        dik_F_CC = 0.0_dp
        djk_F_CC = 0.0_dp
        cd_F_CC = 0.0_dp


        !!! T_CC initialization !!!
        ! fill T_CC(i,j,k) array containing values at knots
        T_CC = 0.0_dp
        T_CC(2,2,1) = -0.035140_dp
        T_CC(2,2,2:9) = -0.0040480_dp

        ! fill rest with zeros
        di_T_CC = 0.0_dp
        dj_T_CC = 0.0_dp
        dk_T_CC = 0.0_dp
        dij_T_CC = 0.0_dp
        dik_T_CC = 0.0_dp
        djk_T_CC = 0.0_dp
        cd_T_CC = 0.0_dp


        !!! F_CH initialization !!!
        ! fill F_CH(i,j,k) array containing values at knots
        F_CH = 0.0_dp

        F_CH(1,1,1:2) = -0.05_dp
        F_CH(1,1,3) = -0.3_dp
        F_CH(1,1,4:9) = -0.05_dp

        F_CH(0,2,5:9) = -0.004523893758064_dp
        F_CH(2,0,5:9) = -0.004523893758064_dp
        F_CH(1,2,2:3) = -0.250_dp
        F_CH(2,1,2:3) = -0.250_dp
        F_CH(3,1,1) = -0.1_dp
        F_CH(1,3,1) = -0.1_dp
        F_CH(3,1,2:3) = -0.125_dp
        F_CH(1,3,2:3) = -0.125_dp
        F_CH(3,1,4:9) = -0.1_dp
        F_CH(1,3,4:9) = -0.1_dp

        ! make top end of F_CH flat instead of zero
        ! also enforces some symmetry
        F_CH(4,0:3,1:9) = F_CH(3,0:3,1:9)
        do i = 0, 3
            do j = i+1, 4
                F_CH(i,j,1:9) = F_CH(j,i,1:9)
            end do
        end do
        F_CH(4,4,1:9) = F_CH(3,4,1:9)

        ! fill rest with zeros
        di_F_CH = 0.0_dp
        dj_F_CH = 0.0_dp
        dk_F_CH = 0.0_dp
        dij_F_CH = 0.0_dp
        dik_F_CH = 0.0_dp
        djk_F_CH = 0.0_dp
        cd_F_CH = 0.0_dp


        !!! F_HH initialization !!!
        ! fill F_HH(i,j,k) array containing values at knots
        F_HH = 0.0_dp
        F_HH(1,1,1) = 0.124915958_dp

        ! fill rest with zeros
        di_F_HH = 0.0_dp
        dj_F_HH = 0.0_dp
        dk_F_HH = 0.0_dp
        dij_F_HH = 0.0_dp
        dik_F_HH = 0.0_dp
        djk_F_HH = 0.0_dp
        cd_F_HH = 0.0_dp



        ! do spline interpolation
        call init_spline2d(P_CC, di_P_CC, dj_P_CC, cd_P_CC, MIN_I_P_CC, MAX_I_P_CC, MIN_J_P_CC, MAX_J_P_CC, P_CC_params)
        call init_spline2d(P_CH, di_P_CH, dj_P_CH, cd_P_CH, MIN_I_P_CH, MAX_I_P_CH, MIN_J_P_CH, MAX_J_P_CH, P_CH_params)

        call init_spline3d(F_CC, di_F_CC, dj_F_CC, dk_F_CC, dij_F_CC, dik_F_CC, djk_F_CC, cd_F_CC,&
            MIN_I_F_CC, MAX_I_F_CC, MIN_J_F_CC, MAX_J_F_CC, MIN_K_F_CC, MAX_K_F_CC, F_CC_params)
        call init_spline3d(T_CC, di_T_CC, dj_T_CC, dk_T_CC, dij_T_CC, dik_T_CC, djk_T_CC, cd_T_CC,&
            MIN_I_T_CC, MAX_I_T_CC, MIN_J_T_CC, MAX_J_T_CC, MIN_K_T_CC, MAX_K_T_CC, T_CC_params)
        call init_spline3d(F_CH, di_F_CH, dj_F_CH, dk_F_CH, dij_F_CH, dik_F_CH, djk_F_CH, cd_F_CH,&
            MIN_I_F_CH, MAX_I_F_CH, MIN_J_F_CH, MAX_J_F_CH, MIN_K_F_CH, MAX_K_F_CH, F_CH_params)
        call init_spline3d(F_HH, di_F_HH, dj_F_HH, dk_F_HH, dij_F_HH, dik_F_HH, djk_F_HH, cd_F_HH,&
            MIN_I_F_HH, MAX_I_F_HH, MIN_J_F_HH, MAX_J_F_HH, MIN_K_F_HH, MAX_K_F_HH, F_HH_params)


        ! gSpline
        gCdom = [ -1.0_dp, -2.0_dp/3, -0.5_dp, -1.0_dp/3, 1.0_dp ]
        gHdom = [ -1.0_dp, -5.0_dp/6, -0.5_dp,            1.0_dp ]

        gC1(1,1) = 0.281695_dp
        gC1(1,2) = 1.062743_dp
        gC1(1,3) = 2.1363075_dp
        gC1(1,4) = 2.5334145_dp
        gC1(1,5) = 1.5544035_dp
        gC1(1,6) = 0.3862485_dp
        gC1(2,1) = 0.282739_dp
        gC1(2,2) = 1.071877_dp
        gC1(2,3) = 2.1681365_dp
        gC1(2,4) = 2.588571_dp
        gC1(2,5) = 1.60191_dp
        gC1(2,6) = 0.402516_dp
        gC1(3,1) = 0.690025_dp
        gC1(3,2) = 5.46016_dp
        gC1(3,3) = 23.0108_dp
        gC1(3,4) = 54.90864_dp
        gC1(3,5) = 68.6124_dp
        gC1(3,6) = 34.705152_dp
        gC1(4,1) = 0.2718560918_dp
        gC1(4,2) = 0.4892740137_dp
        gC1(4,3) = -0.4328177539_dp
        gC1(4,4) = -0.5616817383_dp
        gC1(4,5) = 1.270870225_dp
        gC1(4,6) = -0.0375008379_dp

        gC2(1,1) = 0.281695_dp
        gC2(1,2) = 1.062743_dp
        gC2(1,3) = 2.1363075_dp
        gC2(1,4) = 2.5334145_dp
        gC2(1,5) = 1.5544035_dp
        gC2(1,6) = 0.3862485_dp
        gC2(2,1) = 0.282739_dp
        gC2(2,2) = 1.071877_dp
        gC2(2,3) = 2.1681365_dp
        gC2(2,4) = 2.588571_dp
        gC2(2,5) = 1.60191_dp
        gC2(2,6) = 0.402516_dp
        gC2(3,1) = 0.690025_dp
        gC2(3,2) = 5.46016_dp
        gC2(3,3) = 23.0108_dp
        gC2(3,4) = 54.90864_dp
        gC2(3,5) = 68.6124_dp
        gC2(3,6) = 34.705152_dp
        gC2(4,1) = 0.3754514434_dp
        gC2(4,2) = 1.407269131_dp
        gC2(4,3) = 2.255132012_dp
        gC2(4,4) = 2.028874746_dp
        gC2(4,5) = 1.426920732_dp
        gC2(4,6) = 0.5063519355_dp

        gH(1,1) = 270.4568000026_dp
        gH(1,2) = 1549.63580001436_dp
        gH(1,3) = 3781.77190003166_dp
        gH(1,4) = 4582.15440003486_dp
        gH(1,5) = 2721.43080001916_dp
        gH(1,6) = 630.63360000426_dp
        gH(2,1) = 16.9534406256_dp
        gH(2,2) = -21.08238756_dp
        gH(2,3) = -102.46836_dp
        gH(2,4) = -210.643236_dp
        gH(2,5) = -229.847136_dp
        gH(2,6) = -94.994646_dp
        gH(3,1) = 19.06502493216_dp
        gH(3,2) = 2.0177562846_dp
        gH(3,3) = -2.56642191986_dp
        gH(3,4) = 3.29133223466_dp
        gH(3,5) = -2.65356150626_dp
        gH(3,6) = 0.83766997536_dp


        ! inform module about initialization
        REBO_INIT_SUCCESS = .true.

    end subroutine rebosi_initialize




    ! Number of C and H neighbors are in wrong order in the paper :@
    subroutine rebosi2d(inp_i, inp_j, spec, pij, dN2)

        ! arguments
        character(len=4), intent(in)        :: spec
        real(dp), intent(in)                :: inp_i, inp_j
        real(dp), intent(out)               :: pij
        real(dp), dimension(2), intent(out) :: dN2

        ! function variables
        integer :: m, n
        real(dp) :: coeff
        real(dp) :: pos_i, pos_j
        integer :: idx_i, idx_j

        pij = 0.0_dp
        dN2 = 0.0_dp


        ! copy so that input statements are not being changed
        pos_i = inp_i
        pos_j = inp_j


        if (spec == "P_CC") then

            ! remap inside bounds
            if (pos_i <  MIN_I_P_CC) pos_i = MIN_I_P_CC
            if (pos_j <  MIN_J_P_CC) pos_j = MIN_J_P_CC
            if (pos_i >= MAX_I_P_CC) pos_i = MAX_I_P_CC - TOLERANCE
            if (pos_j >= MAX_J_P_CC) pos_j = MAX_J_P_CC - TOLERANCE

            idx_i = floor(pos_i)
            idx_j = floor(pos_j)

            ! check if close to knot
            if (abs(pos_i - idx_i) < TOLERANCE .and. &
                abs(pos_j - idx_j) < TOLERANCE) then

                pij    = P_CC(   idx_i, idx_j)
                dN2(1) = di_P_CC(idx_i, idx_j)
                dN2(2) = dj_P_CC(idx_i, idx_j)

            else
                do m=0,3
                    do n=0,3
                        coeff = P_CC_params(idx_i, idx_j, m*4+n+1)
                        pij = pij + coeff * pos_i**m * pos_j**n
                        if (m>0) dN2(1) = dN2(1) + coeff * m*pos_i**(m-1) * pos_j**n
                        if (n>0) dN2(2) = dN2(2) + coeff * pos_i**m       * n*pos_j**(n-1)
                    end do
                end do
            end if


        else if (spec == "P_CH") then

            ! remap inside bounds
            if (pos_i <  MIN_I_P_CH) pos_i = MIN_I_P_CH
            if (pos_j <  MIN_J_P_CH) pos_j = MIN_J_P_CH
            if (pos_i >= MAX_I_P_CH) pos_i = MAX_I_P_CH - TOLERANCE
            if (pos_j >= MAX_J_P_CH) pos_j = MAX_J_P_CH - TOLERANCE

            idx_i = floor(pos_i)
            idx_j = floor(pos_j)

            ! check if close to knot
            if (abs(pos_i - idx_i) < TOLERANCE .and. &
                abs(pos_j - idx_j) < TOLERANCE) then

                pij    = P_CH(   idx_i, idx_j)
                dN2(1) = di_P_CH(idx_i, idx_j)
                dN2(2) = dj_P_CH(idx_i, idx_j)

            else
                do m=0,3
                    do n=0,3
                        coeff = P_CH_params(idx_i, idx_j, m*4+n+1)
                        pij = pij + coeff * pos_i**m * pos_j**n
                        if (m>0) dN2(1) = dN2(1) + coeff * m*pos_i**(m-1) * pos_j**n
                        if (n>0) dN2(2) = dN2(2) + coeff * pos_i**m       * n*pos_j**(n-1)
                    end do
                end do
            end if

        else
            print *, "ERROR: Unknown 2D spline matrix:", spec
            stop
        end if

    end subroutine rebosi2d




    subroutine rebosi3d(inp_i, inp_j, inp_k, spec, fij, dN3)

        ! arguments
        character(len=4), intent(in)       :: spec
        real(dp), intent(in)                :: inp_i, inp_j, inp_k
        real(dp), intent(out)               :: fij
        real(dp), dimension(3), intent(out) :: dN3

        ! function variables
        real(dp)                :: pos_i, pos_j, pos_k
        real(dp), dimension(64) :: coeff
        integer                :: idx_i, idx_j, idx_k

        fij   = 0.0_dp
        dN3   = 0.0_dp

        ! copy so that input statements are not being changed
        pos_i = inp_i
        pos_j = inp_j
        pos_k = inp_k

        ! F_CH(i,j,k>9) = F_CH(i,j,9), page 797
        if (spec == "F_CH") then

            ! remap inside bounds
            if (pos_i < MIN_I_F_CH) pos_i = MIN_I_F_CH
            if (pos_j < MIN_J_F_CH) pos_j = MIN_J_F_CH
            if (pos_k < MIN_K_F_CH) pos_k = MIN_K_F_CH

            if (pos_i >= MAX_I_F_CH) pos_i = MAX_I_F_CH - TOLERANCE
            if (pos_j >= MAX_J_F_CH) pos_j = MAX_J_F_CH - TOLERANCE
            if (pos_k >= MAX_K_F_CH) pos_k = MAX_K_F_CH - TOLERANCE

            idx_i = floor(pos_i)
            idx_j = floor(pos_j)
            idx_k = floor(pos_k)

            if (abs(pos_i-idx_i) < TOLERANCE .and. &
                abs(pos_j-idx_j) < TOLERANCE .and. &
                abs(pos_k-idx_k) < TOLERANCE) then
                fij    = F_CH(   idx_i, idx_j, idx_k)
                dN3(1) = di_F_CH(idx_i, idx_j, idx_k)
                dN3(2) = dj_F_CH(idx_i, idx_j, idx_k)
                dN3(3) = dk_F_CH(idx_i, idx_j, idx_k)

            else
                coeff = F_CH_params(idx_i, idx_j, idx_k, :)
                call Sptricubic(pos_i, pos_j, pos_k, coeff, fij, dN3)

            end if


        else if (spec == "F_CC") then

            ! remap inside bounds
            if (pos_i < MIN_I_F_CC) pos_i = MIN_I_F_CC
            if (pos_j < MIN_J_F_CC) pos_j = MIN_J_F_CC
            if (pos_k < MIN_K_F_CC) pos_k = MIN_K_F_CC

            if (pos_i >= MAX_I_F_CC) pos_i = MAX_I_F_CC - TOLERANCE
            if (pos_j >= MAX_J_F_CC) pos_j = MAX_J_F_CC - TOLERANCE
            if (pos_k >= MAX_K_F_CC) pos_k = MAX_K_F_CC - TOLERANCE

            idx_i = floor(pos_i)
            idx_j = floor(pos_j)
            idx_k = floor(pos_k)

            if (abs(pos_i-idx_i) < TOLERANCE .and. &
                abs(pos_j-idx_j) < TOLERANCE .and. &
                abs(pos_k-idx_k) < TOLERANCE) then
                fij    = F_CC(   idx_i, idx_j, idx_k)
                dN3(1) = di_F_CC(idx_i, idx_j, idx_k)
                dN3(2) = dj_F_CC(idx_i, idx_j, idx_k)
                dN3(3) = dk_F_CC(idx_i, idx_j, idx_k)

            else
                coeff = F_CC_params(idx_i, idx_j, idx_k, :)
                call Sptricubic(pos_i, pos_j, pos_k, coeff, fij, dN3)

            end if


        ! T_CC(i,j,k>9) = T_CC(i,j,9); T_CC(i>3,j,k) = T_CC(3,j,k); T_CC(i,j>3,k) = T_CC(i,3,k), page 794 is
        else if (spec == "T_CC") then

            ! remap inside bounds
            if (pos_i < MIN_I_T_CC) pos_i = MIN_I_T_CC
            if (pos_j < MIN_J_T_CC) pos_j = MIN_J_T_CC
            if (pos_k < MIN_K_T_CC) pos_k = MIN_K_T_CC

            ! constraints
            if (pos_i >= MAX_I_T_CC) pos_i = MAX_I_T_CC - TOLERANCE
            if (pos_j >= MAX_J_T_CC) pos_j = MAX_J_T_CC - TOLERANCE
            if (pos_k >= MAX_K_T_CC) pos_k = MAX_K_T_CC - TOLERANCE

            idx_i = floor(pos_i)
            idx_j = floor(pos_j)
            idx_k = floor(pos_k)

            if (abs(pos_i-idx_i) < TOLERANCE .and. &
                abs(pos_j-idx_j) < TOLERANCE .and. &
                abs(pos_k-idx_k) < TOLERANCE) then
                fij    = T_CC(   idx_i, idx_j, idx_k)
                dN3(1) = di_T_CC(idx_i, idx_j, idx_k)
                dN3(2) = dj_T_CC(idx_i, idx_j, idx_k)
                dN3(3) = dk_T_CC(idx_i, idx_j, idx_k)

            else
                coeff = T_CC_params(idx_i, idx_j, idx_k, :)
                call Sptricubic(pos_i, pos_j, pos_k, coeff, fij, dN3)

            end if

        ! F_HH(i,j,k)
        else if (spec == "F_HH") then

            ! remap inside bounds
            if (pos_i < MIN_I_F_HH) pos_i = MIN_I_F_HH
            if (pos_j < MIN_J_F_HH) pos_j = MIN_J_F_HH
            if (pos_k < MIN_K_F_HH) pos_k = MIN_K_F_HH

            if (pos_i >= MAX_I_F_HH) pos_i = MAX_I_F_HH - TOLERANCE
            if (pos_j >= MAX_J_F_HH) pos_j = MAX_J_F_HH - TOLERANCE
            if (pos_k >= MAX_K_F_HH) pos_k = MAX_K_F_HH - TOLERANCE

            idx_i = floor(pos_i)
            idx_j = floor(pos_j)
            idx_k = floor(pos_k)

            if (abs(pos_i-idx_i) < TOLERANCE .and. &
                abs(pos_j-idx_j) < TOLERANCE .and. &
                abs(pos_k-idx_k) < TOLERANCE) then
                fij    = F_HH(   idx_i, idx_j, idx_k)
                dN3(1) = di_F_HH(idx_i, idx_j, idx_k)
                dN3(2) = dj_F_HH(idx_i, idx_j, idx_k)
                dN3(3) = dk_F_HH(idx_i, idx_j, idx_k)

            else
                coeff = F_HH_params(idx_i, idx_j, idx_k, :)
                call Sptricubic(pos_i, pos_j, pos_k, coeff, fij, dN3)

            end if

        else
            print *, "ERROR: Unknown 3D spline matrix:", spec
        end if

    end subroutine rebosi3d


    subroutine gSpline(atoms, costh, Nij, itype, spline, dgdc, dgdN)

        type(universe), intent(in) :: atoms
        real(dp), intent(in)    :: Nij
        integer, intent(in)     :: itype
        real(dp), intent(inout) :: costh
        real(dp), intent(out)   :: spline, dgdc, dgdN

        real(dp), dimension(6) :: coeffs
        real(dp) :: dS, g1, g2, dg1, dg2, cut
        integer :: i


        ! central atom is Carbon
        if (is_carbon(atoms, itype)) then

            if (costh < gCdom(1)) costh = gCdom(1)
            if (costh > gCdom(5)) costh = gCdom(5)

            if (Nij >= NSF_HIGH) then
                do i = 1, 4
                    if (costh >= gCdom(i) .and. costh <= gCdom(i+1)) coeffs = gC2(i,:)
                end do

                call Sp5th(costh, coeffs, spline, dgdc)
                dgdN = 0.0_dp

            else if (Nij <= NSF_LOW) then
                do i = 1, 4
                    if (costh >= gCdom(i) .and. costh <= gCdom(i+1)) coeffs = gC1(i,:)
                end do

                call Sp5th(costh, coeffs, spline, dgdc)
                dgdN = 0.0_dp

            else if (Nij > NSF_LOW .and. Nij < NSF_HIGH) then
                do i = 1, 4
                    if (costh >= gCdom(i) .and. costh <= gCdom(i+1)) coeffs = gC1(i,:)
                end do
                call Sp5th(costh, coeffs, g1, dg1)

                do i = 1, 4
                    if (costh >= gCdom(i) .and. costh <= gCdom(i+1)) coeffs = gC2(i,:)
                end do

                call Sp5th(costh, coeffs, g2, dg2)

                call cufu_one(Nij, NSF_LOW, NSF_HIGH, cut, dS)
                spline = g2 + cut * (g1-g2)
                dgdc = dg2 + cut * (dg1-dg2)
                dgdN = dS * (g1-g2)
            end if

        ! central atom is Hydrogen
        else if (is_hydrogen(atoms, itype)) then

            if (costh < gHdom(1)) costh = gHdom(1)
            if (costh > gHdom(4)) costh = gHdom(4)

            do i = 1, 3
                if (costh >= gHdom(i) .and. costh <= gHdom(i+1)) coeffs = gH(i,:)
            end do

            call Sp5th(costh, coeffs, spline, dgdc)
            dgdN = 0.0_dp
        end if

    end subroutine gSpline


    subroutine Sptricubic(i, j, k, coeffs, f, df)

        real(dp),                intent(in)  :: i, j, k
        real(dp), dimension(64), intent(in)  :: coeffs
        real(dp),                intent(out) :: f
        real(dp), dimension(3),  intent(out) :: df

        real(dp) :: in, in1, jn, jn1, kn, kn1, c
        integer :: ii, jj, kk

        f = 0.0_dp
        df = 0.0_dp

        in = 1.0_dp
        do ii = 0, 3

            jn = 1.0_dp
            do jj = 0, 3

                kn = 1.0_dp
                do kk = 0, 3

                    c = coeffs(16*ii + 4*jj + kk + 1)
                    f = f + c * in * jn * kn
                    if (ii > 0) df(1) = df(1) + c * ii * in1 * jn  * kn
                    if (jj > 0) df(2) = df(2) + c * jj * in  * jn1 * kn
                    if (kk > 0) df(3) = df(3) + c * kk * in  * jn  * kn1
                    kn1 = kn
                    kn  = kn * k
                end do
                jn1 = jn
                jn  = jn * j
            end do
            in1 = in
            in  = in * i
        end do
    end subroutine Sptricubic


    subroutine Sp5th(x, coeffs, spline, dspline)

        real(dp), intent(in) :: x
        real(dp), dimension(6), intent(in) :: coeffs
        real(dp), intent(out) :: spline, dspline

        real(dp) :: x2, x3

        x2 = x*x
        x3 = x2*x

        spline  = coeffs(1) + coeffs(2)*x + coeffs(3)*x2 + coeffs(4)*x3 &
            + coeffs(5)*x2*x2 + coeffs(6)*x2*x3

        dspline = coeffs(2) + 2.0_dp*coeffs(3)*x + 3.0_dp*coeffs(4)*x2 &
            + 4.0_dp*coeffs(5)*x3 + 5.0_dp*coeffs(6)*x2*x2

    end subroutine Sp5th





    subroutine init_spline2d(v_mat2d, di_mat2d, dj_mat2d, cr_mat2d, i_min, i_max, &
        j_min, j_max, spline_2d_params)

        ! declarations
        integer, intent(in) :: i_min, i_max, j_min, j_max    ! define array bounds
        real(dp), dimension(i_min:i_max,j_min:j_max), intent(in) :: v_mat2d, di_mat2d, dj_mat2d, cr_mat2d ! values and derivates of and at knots
        real(dp), dimension(i_min:i_max-1,j_min:j_max-1,16), intent(inout) :: spline_2d_params

        real(dp), dimension(16,16) :: coeff_mat ! coefficient matrix
        real(dp), dimension(16)    :: rhs_vec   ! right hand site vector of A*x=B
        integer                    :: i_knot_low, i_knot_high, j_knot_low, j_knot_high  ! choose corners to interpolate between

        integer, dimension(16)     :: ipiv ! needed for sys. lin. eqns. solver
        integer                    :: info ! needed for sys. lin. eqns. solver
        external                   :: dgesv


        integer                    :: m, n, i_loop, j_loop ! loop variables

        spline_2d_params = 0.0_dp

        ! loop over the grid
        do i_loop=i_min, i_max-1
            do j_loop=j_min, j_max-1

                ! find four knots to interpolate between
                i_knot_low = i_loop
                i_knot_high = i_loop+1
                j_knot_low = j_loop
                j_knot_high = j_loop+1


                ! generate coefficient matrix
                coeff_mat = 0.0_dp

                ! f(x,y)
                do m=0,3
                    do n=0,3
                        coeff_mat(1,4*m+n+1) = i_knot_low**m  * j_knot_low**n
                        coeff_mat(2,4*m+n+1) = i_knot_low**m  * j_knot_high**n
                        coeff_mat(3,4*m+n+1) = i_knot_high**m * j_knot_low**n
                        coeff_mat(4,4*m+n+1) = i_knot_high**m * j_knot_high**n
                    end do
                end do

                ! df(x,y)/dx
                do m=1,3
                    do n=0,3
                        coeff_mat(5,4*m+n+1) = m*i_knot_low**(m-1)  * j_knot_low**n
                        coeff_mat(6,4*m+n+1) = m*i_knot_low**(m-1)  * j_knot_high**n
                        coeff_mat(7,4*m+n+1) = m*i_knot_high**(m-1) * j_knot_low**n
                        coeff_mat(8,4*m+n+1) = m*i_knot_high**(m-1) * j_knot_high**n
                    end do
                end do

                ! df(x,y)/dy
                do m=0,3
                    do n=1,3
                        coeff_mat(9,4*m+n+1) = i_knot_low**m  * n*j_knot_low**(n-1)
                        coeff_mat(10,4*m+n+1) = i_knot_low**m  * n*j_knot_high**(n-1)
                        coeff_mat(11,4*m+n+1) = i_knot_high**m * n*j_knot_low**(n-1)
                        coeff_mat(12,4*m+n+1) = i_knot_high**m * n*j_knot_high**(n-1)
                    end do
                end do

                ! d²f(x,y)/dxdy
                do m=1,3
                    do n=1,3
                        coeff_mat(13,4*m+n+1) = m*i_knot_low**(m-1)  * n*j_knot_low**(n-1)
                        coeff_mat(14,4*m+n+1) = m*i_knot_low**(m-1)  * n*j_knot_high**(n-1)
                        coeff_mat(15,4*m+n+1) = m*i_knot_high**(m-1) * n*j_knot_low**(n-1)
                        coeff_mat(16,4*m+n+1) = m*i_knot_high**(m-1) * n*j_knot_high**(n-1)
                    end do
                end do

                ! generate rhs vector
                rhs_vec(1)  = v_mat2d(i_knot_low, j_knot_low)
                rhs_vec(2)  = v_mat2d(i_knot_low, j_knot_high)
                rhs_vec(3)  = v_mat2d(i_knot_high, j_knot_low)
                rhs_vec(4)  = v_mat2d(i_knot_high, j_knot_high)
                rhs_vec(5)  = di_mat2d(i_knot_low, j_knot_low)
                rhs_vec(6)  = di_mat2d(i_knot_low, j_knot_high)
                rhs_vec(7)  = di_mat2d(i_knot_high, j_knot_low)
                rhs_vec(8)  = di_mat2d(i_knot_high, j_knot_high)
                rhs_vec(9)  = dj_mat2d(i_knot_low, j_knot_low)
                rhs_vec(10) = dj_mat2d(i_knot_low, j_knot_high)
                rhs_vec(11) = dj_mat2d(i_knot_high, j_knot_low)
                rhs_vec(12) = dj_mat2d(i_knot_high, j_knot_high)
                rhs_vec(13) = cr_mat2d(i_knot_low, j_knot_low)
                rhs_vec(14) = cr_mat2d(i_knot_low, j_knot_high)
                rhs_vec(15) = cr_mat2d(i_knot_high, j_knot_low)
                rhs_vec(16) = cr_mat2d(i_knot_high, j_knot_high)

                ! call dgesv( n, nrhs, a, lda, ipiv, b, ldb, info )
                call dgesv(16, 1, coeff_mat, 16, ipiv, rhs_vec, 16, info)
                if ( info > 0 ) then
                    print *, 'The diagonal element of the triangular factor of A,'
                    print *, 'U(', info, ',', info, ') is zero, so that'
                    print *, 'A is singular; the solution could not be computed.'
                else if (info < 0) then
                    print *, -1*info, 'th parameter has illegal value.'
                end if

                ! write coefficients into holder
                spline_2d_params(i_loop,j_loop,1:16) = rhs_vec

            end do
        end do

    end subroutine init_spline2d


    subroutine init_spline3d(v_mat3d, di_mat3d, dj_mat3d, dk_mat3d,&
        dij_mat3d, dik_mat3d, djk_mat3d, cr_mat3d,&
        i_min, i_max, j_min, j_max, k_min, k_max, spline_3d_params)

        ! declarations
        integer, intent(in) :: i_min, i_max, j_min, j_max, k_min, k_max    ! define array bounds
        real(dp), intent(in), dimension(i_min:i_max,j_min:j_max,k_min:k_max) &
            :: v_mat3d, di_mat3d, dj_mat3d, dk_mat3d, dij_mat3d, dik_mat3d, &
            djk_mat3d, cr_mat3d ! values and derivates of and at knots

        real(dp), intent(out) :: spline_3d_params(i_min:i_max-1,j_min:j_max-1,k_min:k_max-1,64)

        integer                    :: i_knot_low, i_knot_high, j_knot_low, j_knot_high,&
            k_knot_low, k_knot_high ! choose corners to interpolate between
        real(dp), dimension(64,64) :: coeff_mat ! coefficient matrix
        real(dp), dimension(64)    :: rhs_vec   ! right hand site vector of A*x=B

        real(dp), dimension(64)    :: ipiv ! needed for sys. lin. eqns. solver
        integer                    :: info ! needed for sys. lin. eqns. solver
        external                   :: dgesv


        integer                    :: m, n, o, i_loop, j_loop, k_loop ! loop variables


        spline_3d_params = 0.0_dp

        do i_loop = i_min, i_max-1
            do j_loop = j_min, j_max-1
                do k_loop = k_min, k_max-1

                    ! find six knots to interpolate between
                    i_knot_low = i_loop
                    i_knot_high = i_loop+1
                    j_knot_low = j_loop
                    j_knot_high = j_loop+1
                    k_knot_low = k_loop
                    k_knot_high = k_loop+1

                    ! generate coefficient matrix
                    coeff_mat = 0.0_dp

                    ! f(x,y,z)
                    do m=0,3
                        do n=0,3
                            do o=0,3
                                coeff_mat(1,16*m+4*n+o+1) = i_knot_low**m  * j_knot_low**n  * k_knot_low**o
                                coeff_mat(2,16*m+4*n+o+1) = i_knot_low**m  * j_knot_low**n  * k_knot_high**o
                                coeff_mat(3,16*m+4*n+o+1) = i_knot_low**m  * j_knot_high**n * k_knot_low**o
                                coeff_mat(4,16*m+4*n+o+1) = i_knot_low**m  * j_knot_high**n * k_knot_high**o
                                coeff_mat(5,16*m+4*n+o+1) = i_knot_high**m * j_knot_low**n  * k_knot_low**o
                                coeff_mat(6,16*m+4*n+o+1) = i_knot_high**m * j_knot_low**n  * k_knot_high**o
                                coeff_mat(7,16*m+4*n+o+1) = i_knot_high**m * j_knot_high**n * k_knot_low**o
                                coeff_mat(8,16*m+4*n+o+1) = i_knot_high**m * j_knot_high**n * k_knot_high**o
                            end do
                        end do
                    end do

                    ! df(x,y,z)/dx
                    do m=1,3
                        do n=0,3
                            do o=0,3
                                coeff_mat(9,16*m+4*n+o+1)  = m*i_knot_low**(m-1)  * j_knot_low**n  * k_knot_low**o
                                coeff_mat(10,16*m+4*n+o+1) = m*i_knot_low**(m-1)  * j_knot_low**n  * k_knot_high**o
                                coeff_mat(11,16*m+4*n+o+1) = m*i_knot_low**(m-1)  * j_knot_high**n * k_knot_low**o
                                coeff_mat(12,16*m+4*n+o+1) = m*i_knot_low**(m-1)  * j_knot_high**n * k_knot_high**o
                                coeff_mat(13,16*m+4*n+o+1) = m*i_knot_high**(m-1) * j_knot_low**n  * k_knot_low**o
                                coeff_mat(14,16*m+4*n+o+1) = m*i_knot_high**(m-1) * j_knot_low**n  * k_knot_high**o
                                coeff_mat(15,16*m+4*n+o+1) = m*i_knot_high**(m-1) * j_knot_high**n * k_knot_low**o
                                coeff_mat(16,16*m+4*n+o+1) = m*i_knot_high**(m-1) * j_knot_high**n * k_knot_high**o
                            end do
                        end do
                    end do

                    ! df(x,y,z)/dy
                    do m=0,3
                        do n=1,3
                            do o=0,3
                                coeff_mat(17,16*m+4*n+o+1) = i_knot_low**m  * n*j_knot_low**(n-1)  * k_knot_low**o
                                coeff_mat(18,16*m+4*n+o+1) = i_knot_low**m  * n*j_knot_low**(n-1)  * k_knot_high**o
                                coeff_mat(19,16*m+4*n+o+1) = i_knot_low**m  * n*j_knot_high**(n-1) * k_knot_low**o
                                coeff_mat(20,16*m+4*n+o+1) = i_knot_low**m  * n*j_knot_high**(n-1) * k_knot_high**o
                                coeff_mat(21,16*m+4*n+o+1) = i_knot_high**m * n*j_knot_low**(n-1)  * k_knot_low**o
                                coeff_mat(22,16*m+4*n+o+1) = i_knot_high**m * n*j_knot_low**(n-1)  * k_knot_high**o
                                coeff_mat(23,16*m+4*n+o+1) = i_knot_high**m * n*j_knot_high**(n-1) * k_knot_low**o
                                coeff_mat(24,16*m+4*n+o+1) = i_knot_high**m * n*j_knot_high**(n-1) * k_knot_high**o
                            end do
                        end do
                    end do

                    ! df(x,y,z)/dz
                    do m=0,3
                        do n=0,3
                            do o=1,3
                                coeff_mat(25,16*m+4*n+o+1) = i_knot_low**m  * j_knot_low**n  * o*k_knot_low**(o-1)
                                coeff_mat(26,16*m+4*n+o+1) = i_knot_low**m  * j_knot_low**n  * o*k_knot_high**(o-1)
                                coeff_mat(27,16*m+4*n+o+1) = i_knot_low**m  * j_knot_high**n * o*k_knot_low**(o-1)
                                coeff_mat(28,16*m+4*n+o+1) = i_knot_low**m  * j_knot_high**n * o*k_knot_high**(o-1)
                                coeff_mat(29,16*m+4*n+o+1) = i_knot_high**m * j_knot_low**n  * o*k_knot_low**(o-1)
                                coeff_mat(30,16*m+4*n+o+1) = i_knot_high**m * j_knot_low**n  * o*k_knot_high**(o-1)
                                coeff_mat(31,16*m+4*n+o+1) = i_knot_high**m * j_knot_high**n * o*k_knot_low**(o-1)
                                coeff_mat(32,16*m+4*n+o+1) = i_knot_high**m * j_knot_high**n * o*k_knot_high**(o-1)
                            end do
                        end do
                    end do

                    ! d²f(x,y,z)/dxdy
                    do m=1,3
                        do n=1,3
                            do o=0,3
                                coeff_mat(33,16*m+4*n+o+1) = m*i_knot_low**(m-1)  * n*j_knot_low**(n-1)  * k_knot_low**o
                                coeff_mat(34,16*m+4*n+o+1) = m*i_knot_low**(m-1)  * n*j_knot_low**(n-1)  * k_knot_high**o
                                coeff_mat(35,16*m+4*n+o+1) = m*i_knot_low**(m-1)  * n*j_knot_high**(n-1) * k_knot_low**o
                                coeff_mat(36,16*m+4*n+o+1) = m*i_knot_low**(m-1)  * n*j_knot_high**(n-1) * k_knot_high**o
                                coeff_mat(37,16*m+4*n+o+1) = m*i_knot_high**(m-1) * n*j_knot_low**(n-1)  * k_knot_low**o
                                coeff_mat(38,16*m+4*n+o+1) = m*i_knot_high**(m-1) * n*j_knot_low**(n-1)  * k_knot_high**o
                                coeff_mat(39,16*m+4*n+o+1) = m*i_knot_high**(m-1) * n*j_knot_high**(n-1) * k_knot_low**o
                                coeff_mat(40,16*m+4*n+o+1) = m*i_knot_high**(m-1) * n*j_knot_high**(n-1) * k_knot_high**o
                            end do
                        end do
                    end do

                    ! d²f(x,y,z)/dxdz
                    do m=1,3
                        do n=0,3
                            do o=1,3
                                coeff_mat(41,16*m+4*n+o+1) = m*i_knot_low**(m-1)  * j_knot_low**n  * o*k_knot_low**(o-1)
                                coeff_mat(42,16*m+4*n+o+1) = m*i_knot_low**(m-1)  * j_knot_low**n  * o*k_knot_high**(o-1)
                                coeff_mat(43,16*m+4*n+o+1) = m*i_knot_low**(m-1)  * j_knot_high**n * o*k_knot_low**(o-1)
                                coeff_mat(44,16*m+4*n+o+1) = m*i_knot_low**(m-1)  * j_knot_high**n * o*k_knot_high**(o-1)
                                coeff_mat(45,16*m+4*n+o+1) = m*i_knot_high**(m-1) * j_knot_low**n  * o*k_knot_low**(o-1)
                                coeff_mat(46,16*m+4*n+o+1) = m*i_knot_high**(m-1) * j_knot_low**n  * o*k_knot_high**(o-1)
                                coeff_mat(47,16*m+4*n+o+1) = m*i_knot_high**(m-1) * j_knot_high**n * o*k_knot_low**(o-1)
                                coeff_mat(48,16*m+4*n+o+1) = m*i_knot_high**(m-1) * j_knot_high**n * o*k_knot_high**(o-1)
                            end do
                        end do
                    end do

                    ! d²f(x,y,z)/dydz
                    do m=0,3
                        do n=1,3
                            do o=1,3
                                coeff_mat(49,16*m+4*n+o+1) = i_knot_low**m  * n*j_knot_low**(n-1)  * o*k_knot_low**(o-1)
                                coeff_mat(50,16*m+4*n+o+1) = i_knot_low**m  * n*j_knot_low**(n-1)  * o*k_knot_high**(o-1)
                                coeff_mat(51,16*m+4*n+o+1) = i_knot_low**m  * n*j_knot_high**(n-1) * o*k_knot_low**(o-1)
                                coeff_mat(52,16*m+4*n+o+1) = i_knot_low**m  * n*j_knot_high**(n-1) * o*k_knot_high**(o-1)
                                coeff_mat(53,16*m+4*n+o+1) = i_knot_high**m * n*j_knot_low**(n-1)  * o*k_knot_low**(o-1)
                                coeff_mat(54,16*m+4*n+o+1) = i_knot_high**m * n*j_knot_low**(n-1)  * o*k_knot_high**(o-1)
                                coeff_mat(55,16*m+4*n+o+1) = i_knot_high**m * n*j_knot_high**(n-1) * o*k_knot_low**(o-1)
                                coeff_mat(56,16*m+4*n+o+1) = i_knot_high**m * n*j_knot_high**(n-1) * o*k_knot_high**(o-1)
                            end do
                        end do
                    end do

                    ! d³f(x,y,z)/dxdydz
                    do m=1,3
                        do n=1,3
                            do o=1,3
                                coeff_mat(57,16*m+4*n+o+1) = m*i_knot_low**(m-1)  * n*j_knot_low**(n-1)  * o*k_knot_low**(o-1)
                                coeff_mat(58,16*m+4*n+o+1) = m*i_knot_low**(m-1)  * n*j_knot_low**(n-1)  * o*k_knot_high**(o-1)
                                coeff_mat(59,16*m+4*n+o+1) = m*i_knot_low**(m-1)  * n*j_knot_high**(n-1) * o*k_knot_low**(o-1)
                                coeff_mat(60,16*m+4*n+o+1) = m*i_knot_low**(m-1)  * n*j_knot_high**(n-1) * o*k_knot_high**(o-1)
                                coeff_mat(61,16*m+4*n+o+1) = m*i_knot_high**(m-1) * n*j_knot_low**(n-1)  * o*k_knot_low**(o-1)
                                coeff_mat(62,16*m+4*n+o+1) = m*i_knot_high**(m-1) * n*j_knot_low**(n-1)  * o*k_knot_high**(o-1)
                                coeff_mat(63,16*m+4*n+o+1) = m*i_knot_high**(m-1) * n*j_knot_high**(n-1) * o*k_knot_low**(o-1)
                                coeff_mat(64,16*m+4*n+o+1) = m*i_knot_high**(m-1) * n*j_knot_high**(n-1) * o*k_knot_high**(o-1)
                            end do
                        end do
                    end do


                    ! generate rhs vector
                    ! f(x,y,z)
                    rhs_vec(1)  = v_mat3d(i_knot_low, j_knot_low, k_knot_low)
                    rhs_vec(2)  = v_mat3d(i_knot_low, j_knot_low, k_knot_high)
                    rhs_vec(3)  = v_mat3d(i_knot_low, j_knot_high, k_knot_low)
                    rhs_vec(4)  = v_mat3d(i_knot_low, j_knot_high, k_knot_high)
                    rhs_vec(5)  = v_mat3d(i_knot_high, j_knot_low, k_knot_low)
                    rhs_vec(6)  = v_mat3d(i_knot_high, j_knot_low, k_knot_high)
                    rhs_vec(7)  = v_mat3d(i_knot_high, j_knot_high, k_knot_low)
                    rhs_vec(8)  = v_mat3d(i_knot_high, j_knot_high, k_knot_high)

                    ! df(x,y,z)/dx
                    rhs_vec(9)  = di_mat3d(i_knot_low, j_knot_low, k_knot_low)
                    rhs_vec(10) = di_mat3d(i_knot_low, j_knot_low, k_knot_high)
                    rhs_vec(11) = di_mat3d(i_knot_low, j_knot_high, k_knot_low)
                    rhs_vec(12) = di_mat3d(i_knot_low, j_knot_high, k_knot_high)
                    rhs_vec(13) = di_mat3d(i_knot_high, j_knot_low, k_knot_low)
                    rhs_vec(14) = di_mat3d(i_knot_high, j_knot_low, k_knot_high)
                    rhs_vec(15) = di_mat3d(i_knot_high, j_knot_high, k_knot_low)
                    rhs_vec(16) = di_mat3d(i_knot_high, j_knot_high, k_knot_high)

                    ! df(x,y,z)/dy
                    rhs_vec(17) = dj_mat3d(i_knot_low, j_knot_low, k_knot_low)
                    rhs_vec(18) = dj_mat3d(i_knot_low, j_knot_low, k_knot_high)
                    rhs_vec(19) = dj_mat3d(i_knot_low, j_knot_high, k_knot_low)
                    rhs_vec(20) = dj_mat3d(i_knot_low, j_knot_high, k_knot_high)
                    rhs_vec(21) = dj_mat3d(i_knot_high, j_knot_low, k_knot_low)
                    rhs_vec(22) = dj_mat3d(i_knot_high, j_knot_low, k_knot_high)
                    rhs_vec(23) = dj_mat3d(i_knot_high, j_knot_high, k_knot_low)
                    rhs_vec(24) = dj_mat3d(i_knot_high, j_knot_high, k_knot_high)

                    ! df(x,y,z)/dz
                    rhs_vec(25) = dk_mat3d(i_knot_low, j_knot_low, k_knot_low)
                    rhs_vec(26) = dk_mat3d(i_knot_low, j_knot_low, k_knot_high)
                    rhs_vec(27) = dk_mat3d(i_knot_low, j_knot_high, k_knot_low)
                    rhs_vec(28) = dk_mat3d(i_knot_low, j_knot_high, k_knot_high)
                    rhs_vec(29) = dk_mat3d(i_knot_high, j_knot_low, k_knot_low)
                    rhs_vec(30) = dk_mat3d(i_knot_high, j_knot_low, k_knot_high)
                    rhs_vec(31) = dk_mat3d(i_knot_high, j_knot_high, k_knot_low)
                    rhs_vec(32) = dk_mat3d(i_knot_high, j_knot_high, k_knot_high)

                    ! d²f(x,y,z)/dxdy
                    rhs_vec(33) = dij_mat3d(i_knot_low, j_knot_low, k_knot_low)
                    rhs_vec(34) = dij_mat3d(i_knot_low, j_knot_low, k_knot_high)
                    rhs_vec(35) = dij_mat3d(i_knot_low, j_knot_high, k_knot_low)
                    rhs_vec(36) = dij_mat3d(i_knot_low, j_knot_high, k_knot_high)
                    rhs_vec(37) = dij_mat3d(i_knot_high, j_knot_low, k_knot_low)
                    rhs_vec(38) = dij_mat3d(i_knot_high, j_knot_low, k_knot_high)
                    rhs_vec(39) = dij_mat3d(i_knot_high, j_knot_high, k_knot_low)
                    rhs_vec(40) = dij_mat3d(i_knot_high, j_knot_high, k_knot_high)

                    ! d²f(x,y,z)/dxdz
                    rhs_vec(41) = dik_mat3d(i_knot_low, j_knot_low, k_knot_low)
                    rhs_vec(42) = dik_mat3d(i_knot_low, j_knot_low, k_knot_high)
                    rhs_vec(43) = dik_mat3d(i_knot_low, j_knot_high, k_knot_low)
                    rhs_vec(44) = dik_mat3d(i_knot_low, j_knot_high, k_knot_high)
                    rhs_vec(45) = dik_mat3d(i_knot_high, j_knot_low, k_knot_low)
                    rhs_vec(46) = dik_mat3d(i_knot_high, j_knot_low, k_knot_high)
                    rhs_vec(47) = dik_mat3d(i_knot_high, j_knot_high, k_knot_low)
                    rhs_vec(48) = dik_mat3d(i_knot_high, j_knot_high, k_knot_high)

                    ! d²f(x,y,z)/dydz
                    rhs_vec(49) = djk_mat3d(i_knot_low, j_knot_low, k_knot_low)
                    rhs_vec(50) = djk_mat3d(i_knot_low, j_knot_low, k_knot_high)
                    rhs_vec(51) = djk_mat3d(i_knot_low, j_knot_high, k_knot_low)
                    rhs_vec(52) = djk_mat3d(i_knot_low, j_knot_high, k_knot_high)
                    rhs_vec(53) = djk_mat3d(i_knot_high, j_knot_low, k_knot_low)
                    rhs_vec(54) = djk_mat3d(i_knot_high, j_knot_low, k_knot_high)
                    rhs_vec(55) = djk_mat3d(i_knot_high, j_knot_high, k_knot_low)
                    rhs_vec(56) = djk_mat3d(i_knot_high, j_knot_high, k_knot_high)

                    ! d³f(x,y,z)/dxdydz
                    rhs_vec(57) = cr_mat3d(i_knot_low, j_knot_low, k_knot_low)
                    rhs_vec(58) = cr_mat3d(i_knot_low, j_knot_low, k_knot_high)
                    rhs_vec(59) = cr_mat3d(i_knot_low, j_knot_high, k_knot_low)
                    rhs_vec(60) = cr_mat3d(i_knot_low, j_knot_high, k_knot_high)
                    rhs_vec(61) = cr_mat3d(i_knot_high, j_knot_low, k_knot_low)
                    rhs_vec(62) = cr_mat3d(i_knot_high, j_knot_low, k_knot_high)
                    rhs_vec(63) = cr_mat3d(i_knot_high, j_knot_high, k_knot_low)
                    rhs_vec(64) = cr_mat3d(i_knot_high, j_knot_high, k_knot_high)


                    !    DO m = 1, 64
                    !         WRITE(*,9998) ( coeff_mat(m,n), n = 1,64 )
                    !    END DO

                    !    print *, "NEXST"
                    !    WRITE(*,9999) (rhs_vec(n), n = 1,64 )


9998                FORMAT( 64(:,1X,F6.0) )
9999                FORMAT( 64(:,1X,F11.8) )


                    ! call dgesv( n, nrhs, a, lda, ipiv, b, ldb, info )
                    !    print *, "calling dgesv"
                    call dgesv(64, 1, coeff_mat, 64, ipiv, rhs_vec, 64, info)
                    !    print *, info
                    !    print *, "after dgesv"
                    if ( info > 0 ) then
                        print *, 'The diagonal element of the triangular factor of A,'
                        print *, 'U(', info, ',', info, ') is zero, so that'
                        print *, 'A is singular; the solution could not be computed.'
                    else if (info < 0) then
                        print *, -1*info, 'th parameter has illegal value.'
                    !    else
                    !        print *, "Solve successful."
                    end if

                    spline_3d_params(i_loop, j_loop, k_loop, 1:64) = rhs_vec

                end do
            end do
        end do
    end subroutine init_spline3d


    ! THEY FORGOT PI IN THE PAPER :@
    subroutine cufu_one(xij, xmin, xmax, cutoff, dx)

        real(dp), intent(in)  :: xij, xmin, xmax
        real(dp), intent(out) :: cutoff, dx

        real(dp) :: ratio

        ratio = (xij-xmin) / (xmax-xmin)

        if (ratio <= 0.0_dp) then
            cutoff = 1.0_dp
            dx     = 0.0_dp
        else if (ratio >= 1.0_dp) then
            cutoff = 0.0_dp
            dx     = 0.0_dp
        else
            cutoff = 0.5_dp * (1.0_dp+cos(ratio*PI))
            dx     = -0.5_dp*PI*sin(ratio*PI) / (xmax-xmin)
        end if
    end subroutine cufu_one



    subroutine cufu(xij, xmin, xmax, cutoff, dx)

        real(dp), intent(in)  :: xij(:)
        real(dp), intent(in)  :: xmin, xmax
        real(dp), intent(out) :: cutoff(:), dx(:)

        real(dp) :: ratio(size(xij))

        ratio = (xij-xmin) / (xmax-xmin)

        do i = 1, size(xij)
            if (ratio(i) <= 0.0_dp) then
                cutoff(i) = 1.0_dp
                dx(i)     = 0.0_dp
            else if (ratio(i) >= 1.0_dp) then
                cutoff(i) = 0.0_dp
                dx(i)     = 0.0_dp
            else
                cutoff(i) = 0.5_dp * (1.0_dp+cos(ratio(i)*PI))
                dx(i)     = -0.5_dp*PI*sin(ratio(i)*PI) / (xmax-xmin)
            end if
        end do


    end subroutine cufu


    integer function is_carbon(atoms, type)

        type(universe), intent(in) :: atoms
        integer, intent(in)        :: type

        is_carbon = 0
        if (atoms%name(type) == "C") is_carbon = 1

    end function is_carbon


    integer function is_hydrogen(atoms, type)

        type(universe), intent(in) :: atoms
        integer, intent(in)        :: type

        is_hydrogen = 0
        if (atoms%name(type) == "H" .or. &
            atoms%name(type) == "D" .or. &
            atoms%name(type) == "T") is_hydrogen = 1

    end function is_hydrogen




    integer function kronecker(a,b)

        character(len=3), intent(in) :: a, b

        kronecker = 0
        if (a==b) kronecker = 1

    end function kronecker



    subroutine rebo_bondorder(atoms, i, j, nC, nH, distances, vectors, VA, bij, flag)

        ! arguments
        type(universe), intent(inout) :: atoms
        integer, intent(in) :: i, j
        integer, intent(in) :: flag
        real(dp), dimension(atoms%nbeads, atoms%natoms), intent(in) :: nH, nC
        real(dp), dimension(atoms%nbeads, atoms%natoms, atoms%natoms), intent(in) :: distances
        real(dp), dimension(3, atoms%nbeads, atoms%natoms, atoms%natoms), intent(in) :: vectors
        real(dp), dimension(atoms%nbeads), intent(in) :: VA
        real(dp), dimension(atoms%nbeads), intent(out) :: bij

        ! local variables
        integer :: k, l, n, b
        integer :: itype, jtype, ktype, ltype, ntype
        character(len=*), parameter :: err = "Error in rebo_bondorder(): "

        real(dp), dimension(atoms%nbeads) :: NijC, NijH, NjiC, NjiH, Nki, Nlj
        real(dp), dimension(atoms%nbeads) :: rikmag, rjlmag, rijmag, rjimag, rknmag, rlnmag
        real(dp), dimension(atoms%nbeads) :: r23mag, r21mag, r32mag, r34mag
        real(dp), dimension(atoms%nbeads) :: wij, wik, wjl, w21, w34

        real(dp), dimension(3, atoms%nbeads) :: rij, rik, ril, rji, rjk, rjl, rkn, rln
        real(dp), dimension(3, atoms%nbeads) :: r32, r23, r34, r21
        real(dp), dimension(3, atoms%nbeads) :: cross321, cross234

        real(dp), dimension(atoms%nbeads) :: g, Tij, piRC, dummy, prefactor
        real(dp), dimension(atoms%nbeads) :: PijS, PjiS, pji, pij
        real(dp), dimension(atoms%nbeads) :: Nijconj, NconjtmpI, NconjtmpJ
        real(dp), dimension(atoms%nbeads) :: lamdajik, lamdaijl
        real(dp), dimension(atoms%nbeads) :: rr, rij2, rik2, ril2, rjk2, rjl2
        real(dp), dimension(atoms%nbeads) :: rik2i, rjl2i
        real(dp), dimension(atoms%nbeads) :: rijrik, rijrjl

        real(dp), dimension(atoms%nbeads) :: cos321, sin321, cos234, sin234, cosjik, cosijl
        real(dp), dimension(atoms%nbeads) :: cw, cwnum, cwnom, om1234, Etmp
        real(dp), dimension(atoms%nbeads) :: tmp, tmp2, tmp3
        real(dp), dimension(atoms%nbeads) :: sink2i, sinl2i
        real(dp), dimension(atoms%nbeads) :: aa, at2, costmp
        real(dp), dimension(atoms%nbeads) :: SpN

        ! derivatives
        real(dp), dimension(atoms%nbeads) :: dgdc, dgdN
        real(dp), dimension(atoms%nbeads) :: dwij, dwik, dwjl, dwkn, dwln, dw21, dw34
        real(dp), dimension(atoms%nbeads) :: dNki, dNlj
        real(dp), dimension(atoms%nbeads) :: dctik, dctij, dctil, dctji, dctjk, dctjl
        real(dp), dimension(atoms%nbeads) :: dt1dij, dt1dik, dt1dil, dt1djk, dt1djl
        real(dp), dimension(2, atoms%nbeads) :: dN2
        real(dp), dimension(3, atoms%nbeads) :: dN3

        real(dp), dimension(3, atoms%nbeads) :: dt2dij, dt2dik, dt2djl
        real(dp), dimension(3, atoms%nbeads) :: dcosjikdri, dcosjikdrj, dcosjikdrk
        real(dp), dimension(3, atoms%nbeads) :: dcosijldri, dcosijldrj, dcosijldrl

        ! forces
        real(dp), dimension(atoms%nbeads) :: fcijpc, fcikpc, fcilpc, fcjkpc, fcjlpc
        real(dp), dimension(3, atoms%nbeads) :: fi, fj, fk, fl
        real(dp), dimension(3, atoms%nbeads) :: f1, f2, f3, f4
        real(dp), dimension(3, atoms%nbeads) :: f12, f23, f24, f31, f34

        ! new
        ! real(dp), dimension(3, atoms%nbeads) :: eij, eji, eik, ejl, eijl, ejik, rijl, rjik
        ! real(dp), dimension(atoms%nbeads) :: rijlmag, rjikmag


        logical :: DEBUG = .false.

        ! real(dp) :: temp1, temp2, temp3, temp4, temp5, temp6, temp7



        !        do k = -10, 50
        !        do l = -10, 50
        !        do n = -10, 100
        !            call rebosi3d(real(k, kind=8)/10, real(l, kind=8)/10, real(n, kind=8)/10, "F_CC", piRC, dN3)
        !            print "(4f23.15)", real(k, kind=8)/10,  real(l, kind=8)/10, real(n, kind=8)/10, piRC
        !        end do
        !        end do
        !        end do
        !        stop

        !        do k = -10, 50
        !        do l = -10, 50
        !            call rebosi2d(real(k, kind=8)/10, real(l, kind=8)/10, "P_CH", piRC, dN2)
        !            print "(4f23.15)", real(k, kind=8)/10,  real(l, kind=8)/10, piRC !dN3(1), dN3(2), dN3(3)
        !        end do
        !        end do
        !        stop


        ! determine i and j types
        itype = atoms%idx(i)
        jtype = atoms%idx(j)

        ! find the distance between them
        rijmag = distances(:,i,j)
        rjimag = distances(:,j,i)

        ! determine vector between them
        rij = vectors(:,:,i,j)
        rji = vectors(:,:,j,i)

        ! evaluate interaction
        call cufu(rijmag, pes_rebo%Dmin(itype,jtype), pes_rebo%Dmax(itype,jtype), wij, dwij)

        ! determine neighbors without i-j interaction
        NijC = nC(:,i) - wij*is_carbon(atoms, jtype)
        NijH = nH(:,i) - wij*is_hydrogen(atoms, jtype)
        NjiC = nC(:,j) - wij*is_carbon(atoms, itype)
        NjiH = nH(:,j) - wij*is_hydrogen(atoms, itype)

        ! initialize accumulators to zero
        bij = 0.0_dp
        tmp = 0.0_dp
        tmp2 = 0.0_dp
        tmp3 = 0.0_dp
        dgdc = 0.0_dp
        dgdN = 0.0_dp
        piRC = 0.0_dp
        NconjtmpI = 0.0_dp
        NconjtmpJ = 0.0_dp
        Etmp = 0.0_dp

        do k = 1, atoms%natoms
            ktype = atoms%idx(k)
            if (k /= i .and. k /= j .and. atoms%pes(ktype, itype) == pes_id_rebo) then

                rikmag = distances(:,i,k)
                call cufu(rikmag, pes_rebo%Dmin(itype,ktype), pes_rebo%Dmax(itype,ktype), wik, dwik)
                if (all(wik < TOLERANCE)) cycle

                rik = vectors(:,:,i,k)

                lamdajik = 4.0_dp*is_hydrogen(atoms,itype) * &
                    ((pes_rebo%rho(itype,ktype)-rikmag) - (pes_rebo%rho(itype,jtype)-rijmag))

                Nki = nC(:,k) + nH(:,k) - wik

                cosjik = sum(rij*rik, dim=1) / (rijmag*rikmag)
                cosjik = min(cosjik, 1.0_dp)
                cosjik = max(cosjik,-1.0_dp)

                ! evaluate splines g and derivatives dg
                do b = 1, atoms%nbeads
                    call gSpline(atoms, cosjik(b), NijC(b)+NijH(b), itype, g(b), dgdc(b), dgdN(b))
                end do

                Etmp = Etmp + (wik*g*exp(lamdajik))
                tmp3 = tmp3 + (wik*dgdN*exp(lamdajik))
                call cufu(Nki, CSF_LOW, CSF_HIGH, SpN, dummy)
                NconjtmpI = NconjtmpI + is_carbon(atoms, ktype) * wik * SpN

            end if
        end do

        if (is_carbon(atoms, itype) .and. is_carbon(atoms, jtype)) then
            do b = 1, atoms%nbeads
                call rebosi2d(NijC(b), NijH(b), "P_CC", PijS(b), dN2(:,b))
            end do
        else if (is_carbon(atoms, itype) .and. is_hydrogen(atoms, jtype)) then
            do b = 1, atoms%nbeads
                call rebosi2d(NijC(b), NijH(b), "P_CH", PijS(b), dN2(:,b))
            end do
        else
            PijS = 0.0_dp
            dN2 = 0.0_dp
        end if

        pij = 1.0_dp / sqrt(1.0_dp+Etmp+PijS)
        tmp = -0.5_dp*pij*pij*pij



        ! pij forces
        if (flag == ENERGY_AND_FORCE) then
            do k = 1, atoms%natoms
                ktype = atoms%idx(k)
                if (k /= i .and. k /= j .and. atoms%pes(ktype, itype) == pes_id_rebo) then

                    rikmag = distances(:,i,k)
                    call cufu(rikmag, pes_rebo%Dmin(itype,ktype), pes_rebo%Dmax(itype,ktype), wik, dwik)
                    if (all(wik < TOLERANCE)) cycle

                    rik = vectors(:,:,i,k)

                    lamdajik = 4.0_dp*is_hydrogen(atoms,itype) * &
                        ((pes_rebo%rho(itype,ktype)-rikmag) - (pes_rebo%rho(itype,jtype)-rijmag))

                    cosjik = sum(rij*rik, dim=1) / (rijmag*rikmag)
                    cosjik = min(cosjik, 1.0_dp)
                    cosjik = max(cosjik,-1.0_dp)

                    do b = 1, atoms%nbeads
                        dcosjikdri(:,b) = (rij(:,b)+rik(:,b))/(rijmag(b)*rikmag(b)) &
                            - cosjik(b)*((rij(:,b)/(rijmag(b)*rijmag(b)))+(rik(:,b)/(rikmag(b)*rikmag(b))))
                        dcosjikdrj(:,b) = -rik(:,b)/(rijmag(b)*rikmag(b)) &
                            + cosjik(b)*( rij(:,b)/(rijmag(b)*rijmag(b)))
                        dcosjikdrk(:,b) = -rij(:,b)/(rijmag(b)*rikmag(b)) &
                            + cosjik(b)*( rik(:,b)/(rikmag(b)*rikmag(b)))
                    end do

                    ! evaluate splines g and derivatives dg
                    do b = 1, atoms%nbeads
                        call gSpline(atoms, cosjik(b), NijC(b)+NijH(b), itype, g(b), dgdc(b), dgdN(b))
                    end do

                    tmp2 = VA * 0.5_dp * tmp * wik * dgdc * exp(lamdajik)
                    do b = 1, atoms%nbeads
                        fi(:,b) = -tmp2(b) * dcosjikdri(:,b)
                        fj(:,b) = -tmp2(b) * dcosjikdrj(:,b)
                        fk(:,b) = -tmp2(b) * dcosjikdrk(:,b)
                    end do

                    tmp2 = VA * 0.5_dp * tmp * wik * g * exp(lamdajik) * 4.0_dp * is_hydrogen(atoms, itype)
                    do b = 1, atoms%nbeads
                        fi(:,b) = fi(:,b) - tmp2(b) *((-rik(:,b)/rikmag(b))+(rij(:,b)/rijmag(b)))
                        fj(:,b) = fj(:,b) - tmp2(b) * (-rij(:,b)/rijmag(b))
                        fk(:,b) = fk(:,b) - tmp2(b) * ( rik(:,b)/rikmag(b))
                    end do


                    ! coordination forces
                    ! dwik forces
                    tmp2 = VA * 0.5_dp * tmp * dwik * g * exp(lamdajik)/rikmag
                    do b = 1, atoms%nbeads
                        fi(:,b) = fi(:,b) - tmp2(b)*rik(:,b)
                        fk(:,b) = fk(:,b) + tmp2(b)*rik(:,b)
                    end do

                    ! pij forces
                    if (is_carbon(atoms, ktype)) then
                        tmp2 = VA * 0.5_dp * tmp * dN2(1,:) * dwik/rikmag
                    else if (is_hydrogen(atoms, ktype)) then
                        tmp2 = VA * 0.5_dp * tmp * dN2(2,:) * dwik/rikmag
                    else
                        print *, err // "unknown REBO type", atoms%name(ktype)
                        stop
                    end if

                    do b = 1, atoms%nbeads
                        fi(:,b) = fi(:,b) - tmp2(b)*rik(:,b)
                        fk(:,b) = fk(:,b) + tmp2(b)*rik(:,b)
                    end do

                    ! dgdN forces
                    tmp2 = VA * 0.5_dp * tmp * tmp3 * dwik/rikmag
                    do b = 1, atoms%nbeads
                        fi(:,b) = fi(:,b) - tmp2(b)*rik(:,b)
                        fk(:,b) = fk(:,b) + tmp2(b)*rik(:,b)
                    end do

                    ! add to accumulator
                    atoms%f(:,:,i) = atoms%f(:,:,i) + fi
                    atoms%f(:,:,j) = atoms%f(:,:,j) + fj
                    atoms%f(:,:,k) = atoms%f(:,:,k) + fk
                end if
            end do
        end if

        tmp  = 0.0_dp
        tmp2 = 0.0_dp
        tmp3 = 0.0_dp
        Etmp = 0.0_dp


        do l = 1, atoms%natoms
            ltype = atoms%idx(l)
            if (l /= j .and. l /= i .and. atoms%pes(ltype, jtype) == pes_id_rebo) then

                rjlmag = distances(:,j,l)

                call cufu(rjlmag, pes_rebo%Dmin(jtype,ltype), pes_rebo%Dmax(jtype,ltype), wjl, dwjl)
                if (all(wjl < TOLERANCE)) cycle

                rjl = vectors(:,:,j,l)

                lamdaijl = 4.0_dp*is_hydrogen(atoms,jtype) * &
                    ((pes_rebo%rho(jtype,ltype)-rjlmag) - (pes_rebo%rho(jtype,itype)-rjimag))

                Nlj = nC(:,l) + nH(:,l) - wjl

                cosijl = sum(rji*rjl, dim=1) / (rjimag*rjlmag)
                cosijl = min(cosijl, 1.0_dp)
                cosijl = max(cosijl,-1.0_dp)

                ! evaluate splines g and derivatives dg
                do b = 1, atoms%nbeads
                    call gSpline(atoms, cosijl(b), NjiC(b)+NjiH(b), jtype, g(b), dgdc(b), dgdN(b))
                end do

                Etmp = Etmp + (wjl*g*exp(lamdaijl))
                tmp3 = tmp3 + (wjl*dgdN*exp(lamdaijl))
                call cufu(Nlj, CSF_LOW, CSF_HIGH, SpN, dummy)
                NconjtmpJ = NconjtmpJ + is_carbon(atoms, ltype) * wjl * SpN

            end if
        end do


        if (is_carbon(atoms, jtype) .and. is_carbon(atoms, itype)) then
            do b = 1, atoms%nbeads
                call rebosi2d(NjiC(b), NjiH(b), "P_CC", PjiS(b), dN2(:,b))
            end do
        else if (is_carbon(atoms, jtype) .and. is_hydrogen(atoms, itype)) then
            do b = 1, atoms%nbeads
                call rebosi2d(NjiC(b), NjiH(b), "P_CH", PjiS(b), dN2(:,b))
            end do
        else
            PjiS = 0.0_dp
            dN2 = 0.0_dp
        end if

        pji = 1.0_dp / sqrt(1.0_dp+Etmp+PjiS)
        tmp = -0.5_dp*pji*pji*pji


        ! pji forces
        if (flag == ENERGY_AND_FORCE) then
            do l = 1, atoms%natoms
                ltype = atoms%idx(l)
                if (l /= j .and. l /= i .and. atoms%pes(ltype, jtype) == pes_id_rebo) then

                    rjlmag = distances(:,j,l)
                    call cufu(rjlmag, pes_rebo%Dmin(jtype,ltype), pes_rebo%Dmax(jtype,ltype), wjl, dwjl)
                    if (all(wjl < TOLERANCE)) cycle

                    rjl = vectors(:,:,j,l)

                    lamdaijl = 4.0_dp*is_hydrogen(atoms,jtype) * &
                        ((pes_rebo%rho(jtype,ltype)-rjlmag) - (pes_rebo%rho(jtype,itype)-rjimag))

                    cosijl = sum(rji*rjl, dim=1) / (rjimag*rjlmag)
                    cosijl = min(cosijl, 1.0_dp)
                    cosijl = max(cosijl,-1.0_dp)

                    do b = 1, atoms%nbeads
                        dcosijldri(:,b) = -rjl(:,b)/(rijmag(b)*rjlmag(b)) &
                            - cosijl(b)*rij(:,b)/(rijmag(b)*rijmag(b))
                        dcosijldrj(:,b) = (-rij(:,b)+rjl(:,b))/(rijmag(b)*rjlmag(b)) &
                            + (cosijl(b)*((rij(:,b)/(rijmag(b)*rijmag(b))) - (rjl(:,b)/(rjlmag(b)*rjlmag(b)))))
                        dcosijldrl(:,b) = rij(:,b)/(rijmag(b)*rjlmag(b)) &
                            + cosijl(b)*rjl(:,b)/(rjlmag(b)*rjlmag(b))
                    end do

                    ! evaluate splines g and derivatives dg
                    do b = 1, atoms%nbeads
                        call gSpline(atoms, cosijl(b), NjiC(b)+NjiH(b), jtype, g(b), dgdc(b), dgdN(b))
                    end do

                    tmp2 = VA * 0.5_dp * tmp * wjl * dgdc * exp(lamdaijl)
                    do b = 1, atoms%nbeads
                        fi(:,b) = -tmp2(b) * dcosijldri(:,b)
                        fj(:,b) = -tmp2(b) * dcosijldrj(:,b)
                        fl(:,b) = -tmp2(b) * dcosijldrl(:,b)
                    end do

                    tmp2 = VA * 0.5_dp * tmp * wjl * g * exp(lamdaijl) * 4.0_dp * is_hydrogen(atoms, jtype)
                    do b = 1, atoms%nbeads
                        fi(:,b) = fi(:,b) - tmp2(b) * rij(:,b)/rijmag(b)
                        fj(:,b) = fj(:,b) - tmp2(b) * (-rjl(:,b)/rjlmag(b) - rij(:,b)/rijmag(b))
                        fl(:,b) = fl(:,b) - tmp2(b) * rjl(:,b)/rjlmag(b)
                    end do


                    ! coordination forces
                    ! dwjl forces
                    tmp2 = VA * 0.5_dp * tmp * dwjl * g * exp(lamdaijl)/rjlmag
                    do b = 1, atoms%nbeads
                        fj(:,b) = fj(:,b) - tmp2(b)*rjl(:,b)
                        fl(:,b) = fl(:,b) + tmp2(b)*rjl(:,b)
                    end do

                    ! pij forces
                    if (is_carbon(atoms, ltype)) then
                        tmp2 = VA * 0.5_dp * tmp * dN2(1,:) * dwjl/rjlmag
                    else if (is_hydrogen(atoms, ltype)) then
                        tmp2 = VA * 0.5_dp * tmp * dN2(2,:) * dwjl/rjlmag
                    else
                        print *, err // "unknown REBO type", atoms%name(ltype)
                        stop
                    end if

                    do b = 1, atoms%nbeads
                        fj(:,b) = fj(:,b) - tmp2(b)*rjl(:,b)
                        fl(:,b) = fl(:,b) + tmp2(b)*rjl(:,b)
                    end do

                    ! dgdN forces
                    tmp2 = VA * 0.5_dp * tmp * tmp3 * dwjl / rjlmag
                    do b = 1, atoms%nbeads
                        fj(:,b) = fj(:,b) - tmp2(b)*rjl(:,b)
                        fl(:,b) = fl(:,b) + tmp2(b)*rjl(:,b)
                    end do

                    ! add to accumulator
                    atoms%f(:,:,i) = atoms%f(:,:,i) + fi
                    atoms%f(:,:,j) = atoms%f(:,:,j) + fj
                    atoms%f(:,:,l) = atoms%f(:,:,l) + fl
                end if
            end do
        end if


        Nijconj = 1.0_dp + NconjtmpI*NconjtmpI + NconjtmpJ*NconjtmpJ

        if (is_carbon(atoms,itype) .and. is_carbon(atoms,jtype)) then
            do b = 1, atoms%nbeads
                call rebosi3d(NijC(b)+NijH(b), NjiC(b)+NjiH(b), Nijconj(b), "F_CC", piRC(b), dN3(:,b))
            end do
        else if (is_carbon(atoms,itype) .and. is_hydrogen(atoms,jtype) .or. &
            is_hydrogen(atoms,itype) .and. is_carbon(atoms,jtype)) then
            do b = 1, atoms%nbeads
                call rebosi3d(NijC(b)+NijH(b), NjiC(b)+NjiH(b), Nijconj(b), "F_CH", piRC(b), dN3(:,b))
            end do
        else if (is_hydrogen(atoms,itype) .and. is_hydrogen(atoms,jtype)) then
            do b = 1, atoms%nbeads
                call rebosi3d(NijC(b)+NijH(b), NjiC(b)+NjiH(b), Nijconj(b), "F_HH", piRC(b), dN3(:,b))
            end do
        else
            print *, "ERROR: Unknown piRC interaction in REBO involving particles", itype, jtype
            stop
        end if

        !spline and dspline test
        !        temp1 = 1e-5
        !        do k = 0, 100
        !            do l = 0, 500
                !do n = 9000, 10000

                    ! left n right
        !                    call rebosi3d( 3.0_dp, real(l, kind=8)/1000-temp1, 9.5_dp, "F_CC", temp2, dN3)
        !                    call rebosi3d( 3.0_dp, real(l, kind=8)/1000+temp1, 9.5_dp, "F_CC", temp3, dN3)
        !
        !                    ! reference
        !                    call rebosi3d( 3.0_dp, real(l, kind=8)/1000, 9.5_dp, "F_CC", piRC, dN3)
        !
        !                    print '(3f, 2f23.15)', 3.0, real(l, kind=8)/1000, 9.5, (temp3-temp2) / (2.0_dp*temp1), dN3(2)
        !        end do
        !        end do
        !        end do
        !        stop
        !        print '(2i, 7f23.15)', i, j, NijC+NijH, NjiC+NjiH, Nijconj, piRC, dN3
        ! piRC forces
        if (flag == ENERGY_AND_FORCE) then
            do k = 1, atoms%natoms
                ktype = atoms%idx(k)
                if (k /= i .and. k /= j .and. atoms%pes(ktype, itype) == pes_id_rebo ) then
                    rikmag = distances(:,i,k)
                    call cufu(rikmag, pes_rebo%Dmin(itype,ktype), pes_rebo%Dmax(itype,ktype), wik, dwik)
                    rik = vectors(:,:,i,k)
                    Nki = nC(:,k) + nH(:,k) - wik
                    call cufu(Nki, CSF_LOW, CSF_HIGH, SpN, dNki)

                    tmp2 = VA * dN3(1,:) * dwik/rikmag
                    do b = 1, atoms%nbeads
                        atoms%f(:,b,i) = atoms%f(:,b,i) - tmp2(b)*rik(:,b)
                        atoms%f(:,b,k) = atoms%f(:,b,k) + tmp2(b)*rik(:,b)
                    end do

                    if (.not. is_carbon(atoms, ktype)) cycle

                    tmp2 = VA * dN3(3,:) * 2.0_dp * NconjtmpI * dwik * SpN/rikmag
                    do b = 1, atoms%nbeads
                        atoms%f(:,b,i) = atoms%f(:,b,i) - tmp2(b)*rik(:,b)
                        atoms%f(:,b,k) = atoms%f(:,b,k) + tmp2(b)*rik(:,b)
                    end do

                    if (any(abs(dNki) > TOLERANCE)) then
                        do n = 1, atoms%natoms
                            ntype = atoms%idx(n)
                            if (n /= i .and. n /= k .and. atoms%pes(ktype, ntype) == pes_id_rebo) then

                                rkn = vectors(:,:,k,n)
                                rknmag = distances(:,k,n)
                                call cufu(rknmag, pes_rebo%Dmin(ktype,ntype), &
                                    pes_rebo%Dmax(ktype,ntype), dummy, dwkn)

                                tmp2 = VA * dN3(3,:) * 2.0_dp * NconjtmpI * wik * dNki * dwkn/rknmag
                                do b = 1, atoms%nbeads
                                    atoms%f(:,b,k) = atoms%f(:,b,k) - tmp2(b)*rkn(:,b)
                                    atoms%f(:,b,n) = atoms%f(:,b,n) + tmp2(b)*rkn(:,b)
                                end do
                            end if
                        end do
                    end if
                end if
            end do

            ! piRC forces
            do l = 1, atoms%natoms
                ltype = atoms%idx(l)
                if (l /= i .and. l /= j .and. atoms%pes(ltype, jtype) == pes_id_rebo) then
                    rjlmag = distances(:,j,l)
                    call cufu(rjlmag, pes_rebo%Dmin(jtype,ltype), pes_rebo%Dmax(jtype,ltype), wjl, dwjl)

                    rjl = vectors(:,:,j,l)
                    Nlj = nC(:,l) + nH(:,l) - wjl
                    call cufu(Nlj, CSF_LOW, CSF_HIGH, SpN, dNlj)

                    tmp2 = VA * dN3(2,:) * dwjl/rjlmag
                    do b = 1, atoms%nbeads
                        atoms%f(:,b,j) = atoms%f(:,b,j) - tmp2(b)*rjl(:,b)
                        atoms%f(:,b,l) = atoms%f(:,b,l) + tmp2(b)*rjl(:,b)
                    end do

                    if (.not. is_carbon(atoms, ltype)) cycle

                    tmp2 = VA * dN3(3,:) * 2.0_dp * NconjtmpJ * dwjl * SpN/rjlmag
                    do b = 1, atoms%nbeads
                        atoms%f(:,b,j) = atoms%f(:,b,j) - tmp2(b)*rjl(:,b)
                        atoms%f(:,b,l) = atoms%f(:,b,l) + tmp2(b)*rjl(:,b)
                    end do

                    if (any(abs(dNlj) > TOLERANCE)) then
                        do n = 1, atoms%natoms
                            ntype = atoms%idx(n)
                            if (n /= j .and. n /= l .and. atoms%pes(ltype, ntype) == pes_id_rebo) then
                                rln = vectors(:,:,l,n)
                                rlnmag = distances(:,l,n)
                                call cufu(rlnmag, pes_rebo%Dmin(ltype,ntype), &
                                    pes_rebo%Dmax(ltype,ntype), dummy, dwln)

                                tmp2 = VA * dN3(3,:) * 2.0_dp * NconjtmpJ * wjl * dNlj * dwln/rlnmag
                                do b = 1, atoms%nbeads
                                    atoms%f(:,b,l) = atoms%f(:,b,l) - tmp2(b)*rln(:,b)
                                    atoms%f(:,b,n) = atoms%f(:,b,n) + tmp2(b)*rln(:,b)
                                end do
                            end if
                        end do
                    end if
                end if
            end do
        end if

        Tij = 0.0_dp
        Etmp = 0.0_dp
        if (is_carbon(atoms,itype) .and. is_carbon(atoms,jtype)) then
            do b = 1, atoms%nbeads
                call rebosi3d(NijC(b)+NijH(b), NjiC(b)+NjiH(b), Nijconj(b), "T_CC", Tij(b), dN3(:,b))
            end do
        else
            Tij = 0.0_dp
            dN3 = 0.0_dp
        end if


        !        do rl = 0.0_dp, 4.0_dp, 0.01_dp
        !            print *, rl, rebosi2d(0.0_dp, rl, "P_CH")
        !        end do
        !        stop



        if (any(abs(Tij) > TOLERANCE)) then
            ! atom2 = i
            ! atom3 = j
            r32 = vectors(:,:,j,i)
            r32mag = distances(:,j,i)
            r23 = -r32
            r23mag = r32mag

            do k = 1, atoms%natoms    ! atom1 = k
                if (k /= i .and. k /= j .and. atoms%pes(ktype, itype) == pes_id_rebo) then

                    ktype = atoms%idx(k)
                    r21mag = distances(:,i,k)
                    call cufu(r21mag, pes_rebo%Dmin(itype,ktype), pes_rebo%Dmaxp(itype,ktype), w21, dw21)

                    r21 = vectors(:,:,i,k)
                    cos321 = -1.0_dp * sum(r21*r32, dim=1) / (r21mag*r32mag)
                    cos321 = min(cos321,  1.0_dp)
                    cos321 = max(cos321, -1.0_dp)
                    sin321 = sqrt(1.0_dp - cos321*cos321)
                    if (all(sin321 > TOLERANCE) .and. all(r21mag > TOLERANCE)) then

                        sink2i = 1.0_dp / (sin321*sin321)
                        rik2i  = 1.0_dp / (r21mag*r21mag)
                        rr = r23mag*r23mag - r21mag*r21mag
                        !                    rjk = r21 - r23
                        rjk = vectors(:,:,j,k)
                        rjk2 = sum(rjk*rjk, dim=1)
                        rijrik = 2.0_dp * r23mag*r21mag
                        rik2 = r21mag*r21mag
                        dctik = (-rr+rjk2)/(rijrik*rik2)
                        dctij = (rr+rjk2)/(rijrik*r23mag*r23mag)
                        dctjk = -2.0_dp/rijrik

                        rijmag = r23mag
                        rikmag = r21mag
                        rij2 = r23mag*r23mag
                        rik2 = r21mag*r21mag
                        costmp = 0.5_dp*(rij2+rik2-rjk2)/rijmag/rikmag
                        !            tspjik = 0.0_dp
                        !            dtsjik = 0.0_dp

                        do l = 1, atoms%natoms    ! atom4 = l
                            if (l /= i .and. l /= j .and. l /= k &
                                .and. atoms%pes(ltype, jtype) == pes_id_rebo) then

                                ltype = atoms%idx(l)
                                r34mag = distances(:,j,l)
                                call cufu(r34mag, pes_rebo%Dmin(jtype,ltype), &
                                    pes_rebo%Dmaxp(jtype,ltype), w34, dw34)

                                r34 = vectors(:,:,j,l)
                                cos234 = sum(r32*r34, dim=1) / (r32mag*r34mag)
                                cos234 = min(cos234,  1.0_dp)
                                cos234 = max(cos234, -1.0_dp)
                                sin234 = sqrt(1.0_dp - cos234*cos234)
                                if (all(sin234 > TOLERANCE) .and. all(r34mag > TOLERANCE)) then

                                    sinl2i = 1.0_dp / (sin234*sin234)
                                    rjl2i  = 1.0_dp / (r34mag*r34mag)
                                    rr = r23mag*r23mag - r34mag*r34mag
                                    !                            ril = r23 + r34
                                    ril = vectors(:,:,i,l)
                                    ril2 = sum(ril*ril, dim=1)
                                    rijrjl = 2.0_dp * r23mag * r34mag
                                    rjl2 = r34mag * r34mag
                                    dctjl = (-rr+ril2) / (rijrjl*rjl2)
                                    dctji = (rr+ril2) / (rijrjl*r23mag*r23mag)
                                    dctil = -2.0_dp / rijrjl
                                    rjlmag = r34mag
                                    rjl2 = r34mag*r34mag
                                    costmp = 0.5_dp * (rij2+rjl2-ril2)/rijmag/rjlmag
                                    !                tspijl = 0.0_dp
                                    !                dtsijl = 0.0_dp
                                    prefactor = VA * Tij

                                    !                            cross321 = cross_product(r32, r21)
                                    !                            cross234 = cross_product(r23, r34)
                                    cross321(1,:) = r32(2,:)*r21(3,:) - r32(3,:)*r21(2,:)
                                    cross321(2,:) = r32(3,:)*r21(1,:) - r32(1,:)*r21(3,:)
                                    cross321(3,:) = r32(1,:)*r21(2,:) - r32(2,:)*r21(1,:)
                                    cross234(1,:) = r23(2,:)*r34(3,:) - r23(3,:)*r34(2,:)
                                    cross234(2,:) = r23(3,:)*r34(1,:) - r23(1,:)*r34(3,:)
                                    cross234(3,:) = r23(1,:)*r34(2,:) - r23(2,:)*r34(1,:)

                                    cwnum = sum(cross321*cross234, dim=1)
                                    cwnom = r21mag*r34mag*r23mag*r23mag*sin321*sin234
                                    om1234 = cwnum/cwnom
                                    cw = om1234
                                    Etmp = Etmp + (1.0-om1234*om1234) * w21 * w34

                                    if (flag == ENERGY_ONLY) cycle

                                    dt1dik = rik2i - dctik * sink2i * cos321
                                    dt1djk = -dctjk * sink2i * cos321
                                    dt1djl = rjl2i - dctjl * sinl2i * cos234
                                    dt1dil = -dctil * sinl2i * cos234
                                    dt1dij = 2.0_dp/(r23mag*r23mag) - dctij*sink2i*cos321 - dctji*sinl2i*cos234

                                    !                            dt2dik = cross_product( r23, cross234)
                                    !                            dt2djl = cross_product(-r23, cross321)
                                    dt2dik(1,:) = -r23(3,:)*cross234(2,:) + r23(2,:)*cross234(3,:)
                                    dt2dik(2,:) = -r23(1,:)*cross234(3,:) + r23(3,:)*cross234(1,:)
                                    dt2dik(3,:) = -r23(2,:)*cross234(1,:) + r23(1,:)*cross234(2,:)

                                    dt2djl(1,:) = -r23(2,:)*cross321(3,:) + r23(3,:)*cross321(2,:)
                                    dt2djl(2,:) = -r23(3,:)*cross321(1,:) + r23(1,:)*cross321(3,:)
                                    dt2djl(3,:) = -r23(1,:)*cross321(2,:) + r23(2,:)*cross321(1,:)

                                    dt2dij(1,:) = r21(3,:)*cross234(2,:) - r34(3,:)*cross321(2,:) - r21(2,:)*cross234(3,:) + r34(2,:)*cross321(3,:)
                                    dt2dij(2,:) = r21(1,:)*cross234(3,:) - r34(1,:)*cross321(3,:) - r21(3,:)*cross234(1,:) + r34(3,:)*cross321(1,:)
                                    dt2dij(3,:) = r21(2,:)*cross234(1,:) - r34(2,:)*cross321(1,:) - r21(1,:)*cross234(2,:) + r34(1,:)*cross321(2,:)

                                    aa = (prefactor*2.0_dp*cw/cwnom)*w21*w34
                                    at2 = aa*cwnum

                                    fcijpc = -dt1dij*at2
                                    fcikpc = -dt1dik*at2
                                    fcjlpc = -dt1djl*at2
                                    fcjkpc = -dt1djk*at2
                                    fcilpc = -dt1dil*at2

                                    do b = 1, atoms%nbeads
                                        f23(:,b) = fcijpc(b)*r23(:,b) + aa(b)*dt2dij(:,b)
                                        f12(:,b) = fcikpc(b)*r21(:,b) + aa(b)*dt2dik(:,b)
                                        f34(:,b) = fcjlpc(b)*r34(:,b) + aa(b)*dt2djl(:,b)
                                        f31(:,b) = fcjkpc(b)*rjk(:,b)
                                        f24(:,b) = fcilpc(b)*ril(:,b)
                                    end do

                                    f1 = -f12 - f31
                                    f2 =  f23 + f12 + f24
                                    f3 = -f23 + f34 + f31
                                    f4 = -f34 - f24

                                    ! coordination forces
                                    tmp2 = VA * Tij * (1.0_dp-om1234*om1234) * dw21 * w34 / r21mag
                                    do b = 1, atoms%nbeads
                                        f2(:,b) = f2(:,b) - tmp2(b)*r21(:,b)
                                        f1(:,b) = f1(:,b) + tmp2(b)*r21(:,b)
                                    end do

                                    tmp2 = VA * Tij * (1.0_dp-om1234*om1234) * w21 * dw34 / r34mag
                                    do b = 1, atoms%nbeads
                                        f3(:,b) = f3(:,b) - tmp2(b)*r34(:,b)
                                        f4(:,b) = f4(:,b) + tmp2(b)*r34(:,b)
                                    end do

                                    atoms%f(:,:,k) = atoms%f(:,:,k) + f1
                                    atoms%f(:,:,i) = atoms%f(:,:,i) + f2
                                    atoms%f(:,:,j) = atoms%f(:,:,j) + f3
                                    atoms%f(:,:,l) = atoms%f(:,:,l) + f4



!                                    ! alternative forces (TODO when time available)
!                                    rij = vectors(:,:,i,j)
!                                    rji = vectors(:,:,j,i)
!                                    rik = vectors(:,:,i,k)
!                                    rjl = vectors(:,:,j,l)
!
!                                    rijmag = distances(:,i,j)
!                                    rjimag = distances(:,j,i)
!                                    rikmag = distances(:,i,k)
!                                    rjlmag = distances(:,j,l)
!
!                                    do b = 1, atoms%nbeads
!                                        eij(:,b) = rij(:,b)/rijmag(b)
!                                        eji(:,b) = rji(:,b)/rjimag(b)
!                                        eik(:,b) = rik(:,b)/rikmag(b)
!                                        ejl(:,b) = rjl(:,b)/rjlmag(b)
!                                    end do
!
!                                    rijl = cro_pro(rji, rik)
!                                    rjik = cro_pro(rij, rjl)
!
!                                    rijlmag = sqrt(sum(rijl*rijl, dim=1))
!                                    rjikmag = sqrt(sum(rjik*rjik, dim=1))
!
!                                    do b = 1, atoms%nbeads
!                                        eijl(:,b) = rijl(:,b)/rijlmag(b)
!                                        ejik(:,b) = rjik(:,b)/rjikmag(b)
!                                    end do



                                    !print *, "kforce", -f12 -f31

                                end if  ! abs(sin234) > TOLERANCE
                            end if      ! (l /= i .and. l /= j .and. l /= k)
                        end do          ! l loop
                    end if              ! abs(sin321) > TOLERANCE
                end if                  ! (k /= i .and. k /= j)
            end do                      ! k loop

            !            print *, force


            ! Tij forces now that we have Etmp
            if (flag == ENERGY_AND_FORCE) then
                do k = 1, atoms%natoms
                    if (k /= i .and. k /= j .and. atoms%pes(ktype, itype) == pes_id_rebo) then

                        ktype = atoms%idx(k)
                        rikmag = distances(:,i,k)
                        call cufu(rikmag, pes_rebo%Dmin(itype,ktype), &
                            pes_rebo%Dmax(itype,ktype), wik, dwik)

                        rik = vectors(:,:,i,k)
                        Nki = nC(:,k)  + nH(:,k) - wik
                        call cufu(Nki, CSF_LOW, CSF_HIGH, SpN, dNki)

                        tmp2 = VA * dN3(1,:) * dwik * Etmp / rikmag
                        do b = 1, atoms%nbeads
                            atoms%f(:,b,i) = atoms%f(:,b,i) - tmp2(b)*rik(:,b)
                            atoms%f(:,b,k) = atoms%f(:,b,k) + tmp2(b)*rik(:,b)
                        end do
                        !                   if ( (i==19 .or. k==19) .and. abs(sum(tmp2*rik))>0) print '(a, 6f23.15, 3i)', "2 1", atoms%f(:,19), tmp2*rik, i, j, k

                        if (.not. is_carbon(atoms, ktype)) cycle

                        tmp2 = VA * dN3(3,:) * 2.0_dp * NconjtmpI * dwik * SpN * Etmp / rikmag
                        do b = 1, atoms%nbeads
                            atoms%f(:,b,i) = atoms%f(:,b,i) - tmp2(b)*rik(:,b)
                            atoms%f(:,b,k) = atoms%f(:,b,k) + tmp2(b)*rik(:,b)
                        end do
                        !                    if ( (i==19 .or. k==19) .and. abs(sum(tmp2*rik))>0) print '(a, 6f23.15)', "2 2", atoms%f(:,19), tmp2*rik

                        if (any(abs(dNki) > TOLERANCE)) then
                            do n = 1, atoms%natoms
                                if (n /= i .and. n /= k .and. atoms%pes(ktype, ntype) == pes_id_rebo) then
                                    ntype = atoms%idx(n)
                                    rkn = vectors(:,:,k,n)
                                    rknmag = distances(:,k,n)
                                    call cufu(rknmag, pes_rebo%Dmin(ktype,ntype), &
                                        pes_rebo%Dmax(ktype,ntype), dummy, dwkn)

                                    tmp2 = VA * dN3(3,:) * 2.0_dp * NconjtmpI * wik * dNki * dwkn * Etmp / rknmag
                                    do b = 1, atoms%nbeads
                                        atoms%f(:,b,k) = atoms%f(:,b,k) - tmp2(b)*rkn(:,b)
                                        atoms%f(:,b,n) = atoms%f(:,b,n) + tmp2(b)*rkn(:,b)
                                    end do
                                !                            if ( (n==19 .or. k==19) .and. abs(sum(tmp2*rkn))>0) print '(a, 6f23.15)', "2 3", force(:,19), tmp2*rkn
                                end if
                            end do
                        end if
                    end if
                end do


                ! Tij forces
                do l = 1, atoms%natoms
                    if (l /= i .and. l /= j .and. atoms%pes(ltype, jtype) == pes_id_rebo) then

                        ltype = atoms%idx(l)
                        rjlmag = distances(:,j,l)
                        call cufu(rjlmag, pes_rebo%Dmin(jtype,ltype), &
                            pes_rebo%Dmax(jtype,ltype), wjl, dwjl)

                        rjl = vectors(:,:,j,l)
                        Nlj = nC(:,l) + nH(:,l) - wjl
                        call cufu(Nlj, CSF_LOW, CSF_HIGH, SpN, dNlj)

                        tmp2 = VA * dN3(2,:) * dwjl * Etmp / rjlmag
                        do b = 1, atoms%nbeads
                            atoms%f(:,b,j) = atoms%f(:,b,j) - tmp2(b)*rjl(:,b)
                            atoms%f(:,b,l) = atoms%f(:,b,l) + tmp2(b)*rjl(:,b)
                        end do
                        !                    if ( (j==19 .or. l==19) .and. abs(sum(tmp2*rjl)>0) ) print '(a, 6f23.15)', "3 1", atoms%f(:,19) , tmp2*rjl

                        if (.not. is_carbon(atoms, ltype)) cycle

                        tmp2 = VA * dN3(3,:) * 2.0_dp * NconjtmpJ * dwjl * SpN * Etmp / rjlmag
                        do b = 1, atoms%nbeads
                            atoms%f(:,b,j) = atoms%f(:,b,j) - tmp2(b)*rjl(:,b)
                            atoms%f(:,b,l) = atoms%f(:,b,l) + tmp2(b)*rjl(:,b)
                        end do
                        !                    if ( (j==19 .or. l==19) .and. abs(sum(tmp2*rjl))>0) print '(a, 6f23.15)', "3 2", atoms%f(:,19), tmp2*rjl

                        if (any(abs(dNlj) > TOLERANCE)) then
                            do n = 1, atoms%natoms
                                if (n /= l .and. n /= j .and. atoms%pes(ltype, ntype) == pes_id_rebo) then

                                    ntype = atoms%idx(n)
                                    rln = vectors(:,:,l,n)
                                    rlnmag = distances(:,l,n)
                                    call cufu(rlnmag, pes_rebo%Dmin(ltype,ntype), &
                                        pes_rebo%Dmax(ltype,ntype), dummy, dwln)

                                    tmp2 = VA * dN3(3,:) * 2.0_dp * NconjtmpJ * wjl * dNlj * dwln * Etmp / rlnmag
                                    do b = 1, atoms%nbeads
                                        atoms%f(:,b,l) = atoms%f(:,b,l) - tmp2(b)*rln(:,b)
                                        atoms%f(:,b,n) = atoms%f(:,b,n) + tmp2(b)*rln(:,b)
                                    end do
                                !                            if ( (j==19 .or. l==19) .and. abs(sum(tmp2*rln))>0) print '(a, 6f23.15)', "3 3", force(:,19), tmp2*rln
                                end if
                            end do
                        end if
                    end if
                end do
            end if  ! flag ENERGY_AND_FORCE
        end if  ! abs(Tij) > TOLERANCE


        !        bij = 0.5_dp*pij
        !        bij = 0.5_dp*(pij+pji) + piRC
        !bij = 0.5_dp*(pij+pji)
        !        bij = (Tij*Etmp)
        bij = 0.5_dp*(pij+pji) + Tij*Etmp + piRC !+ (Tij*Etmp)
        !        bij = Tij


        if (DEBUG ) then ! .and. (neigh_c_i < 2 .or. neigh_c_j < 2).and. j.eq.natoms
            print *, i, j
            print *, "species i ", atoms%idx(i)
            print *, "species j ", atoms%idx(j)
            print *, "atom i has", NijC, "carbon neighbors"
            print *, "atom i has", NijH, "hydrogen neighbors"
            print *, "atom j has", NjiC, "carbon neighbors"
            print *, "atom j has", NjiH, "hydrogen neighbors"
            print *, "dist_ij", distances(:,i,j)
            print *, "num_conj", Nijconj
            print *, "first_bracket", NconjtmpI
            print *, "second_bracket", NconjtmpJ
            print *, "radical", pirc
            print *, "dihedral Tspline", Tij
            print *, "dihedral DS", Etmp
            print *, "bosp_ij", pij
            print *, "bosp_ji", pji
            !    print *, "potential_attractive", potential_attractive
            print *, "bond_order", bij
        !    print *, "potential_repulsive", potential_repulsive
        end if

    end subroutine rebo_bondorder







    subroutine compute_rebo(atoms, flag)

        !   Calculates energy and forces with REBO potential

        type(universe), intent(inout) :: atoms
        integer, intent(in) :: flag

        integer ::  i, j, k, l, b, itype, jtype, ktype, ltype
        real(dp), dimension(atoms%nbeads) :: rsq, rij, wij, fpair, wkl, dummy, dwij
        real(dp), dimension(atoms%nbeads) :: VR, pre, dVRdi, VA, term, bij, dVAdi, dVA, temp_distance

        real(dp), dimension(3, atoms%nbeads) :: temp_vector
        real(dp), dimension(atoms%nbeads, atoms%natoms) :: nH, nC
        real(dp), dimension(atoms%nbeads, atoms%natoms, atoms%natoms) :: distances
        real(dp), dimension(3, atoms%nbeads, atoms%natoms, atoms%natoms) :: vectors


        ! gather distance and vector information
        ! calculate only one half and then add to other half (with changed sign)
        distances = 0.0_dp
        vectors = 0.0_dp
        do j = 1, atoms%natoms-1
            do i = j+1, atoms%natoms
                call minimg_beads(atoms, i, j, temp_distance, temp_vector)
                distances(:,i,j) = temp_distance
                vectors(:,:,i,j) = -temp_vector
            end do
        end do

        do b = 1, atoms%nbeads
            distances(b,:,:) = distances(b,:,:) + transpose(distances(b,:,:))
            vectors(1,b,:,:) = vectors(1,b,:,:) - transpose(vectors(1,b,:,:))
            vectors(2,b,:,:) = vectors(2,b,:,:) - transpose(vectors(2,b,:,:))
            vectors(3,b,:,:) = vectors(3,b,:,:) - transpose(vectors(3,b,:,:))
        end do

        ! determine neighbors
        nC = 0.0_dp
        nH = 0.0_dp

        do k = 1, atoms%natoms - 1
            ktype = atoms%idx(k)

            do l = k+1, atoms%natoms
                ltype = atoms%idx(l)

                if (atoms%pes(ktype, ltype) == pes_id_rebo) then

                    call cufu(distances(:,k,l), pes_rebo%Dmin(ktype,ltype), &
                        pes_rebo%Dmax(ktype,ltype), wkl, dummy)

                    if (is_carbon(  atoms, ktype)) nC(:,l) = nC(:,l) + wkl
                    if (is_hydrogen(atoms, ktype)) nH(:,l) = nH(:,l) + wkl
                    if (is_carbon(  atoms, ltype)) nC(:,k) = nC(:,k) + wkl
                    if (is_hydrogen(atoms, ltype)) nH(:,k) = nH(:,k) + wkl
                end if
            end do
        end do

        ! two-body interactions
        do i = 1, atoms%natoms-1
            itype = atoms%idx(i)

            do j = i+1, atoms%natoms
                jtype = atoms%idx(j)

                if (atoms%pes(itype, jtype) /= pes_id_rebo) cycle

                rij = distances(:,i,j)
                call cufu(rij, pes_rebo%Dmin(itype,jtype), pes_rebo%Dmax(itype,jtype), wij, dwij)

                if (all(wij < tolerance)) cycle

                rsq = rij*rij

                VR = wij * (1.0_dp + (pes_rebo%Q(itype,jtype)/rij)) * &
                    pes_rebo%A(itype,jtype) * exp(-pes_rebo%alpha(itype,jtype)*rij)
                pre = wij * pes_rebo%A(itype,jtype) * exp(-pes_rebo%alpha(itype,jtype)*rij)
                dVRdi = pre * ( (-pes_rebo%alpha(itype,jtype)) - (pes_rebo%Q(itype,jtype)/rsq) &
                    - (pes_rebo%Q(itype,jtype)*pes_rebo%alpha(itype,jtype)/rij) )

                where (wij > tolerance) dVRdi = dVRdi + VR/wij * dwij

                VA = 0.0_dp
                dVA = 0.0_dp
                do k = 1, 3
                    term = -wij * pes_rebo%B(k,itype,jtype) * exp(-pes_rebo%beta(k,itype,jtype)*rij)
                    VA   = VA + term
                    dVA  = dVA - pes_rebo%beta(k,itype,jtype) * term
                end do

                where (wij > tolerance) dVA = dVA + VA/wij * dwij

                call rebo_bondorder(atoms, i, j, nC, nH, distances, vectors, VA, bij, flag)
                dVAdi = bij*dVA

                fpair = -(dVRdi+dVAdi) / rij

                do b = 1, atoms%nbeads
                    atoms%f(:,b,i) = atoms%f(:,b,i) + vectors(:,b,i,j)*fpair(b)
                    atoms%f(:,b,j) = atoms%f(:,b,j) + vectors(:,b,j,i)*fpair(b)
                end do

                !            if (i==19) print '(a, 6f23.15)', "mother", all_force(:,19), del*fpair
                !            if (j==19) print '(a, 6f23.15)', "mother", all_force(:,19), del*fpair

                atoms%epot = atoms%epot + (VR + bij*VA)
            end do
        end do

    end subroutine compute_rebo

end module pes_rebo_mod




