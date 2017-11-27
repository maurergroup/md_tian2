!** NONLINEAR LEAST SQUARE PROBLEM WITHOUT BOUNDARY CONSTRAINTS
    include 'mkl_rci.f90'
module fit

    use force
    use universe_mod
    use run_config, only : simparams
    use open_file
    use useful_things

    implicit none

    character(len=*), parameter :: train_folder = "fit/train/"
    character(len=*), parameter :: valid_folder = "fit/valid/"
    integer :: training_data_points, validation_data_points

contains

subroutine perform_fit(reference)

    type(universe), intent(inout) :: reference     ! holds the reference configuration

    type(universe), allocatable :: train_data(:), valid_data(:)
    real(dp) :: train_rmse, valid_rmse
    integer :: alloc_stat, i
    character(len=*), parameter :: err = "Error in perform_fit(): "


    training_data_points   = 0
    validation_data_points = 0

    call check_fitting_data(train_folder)
    call check_fitting_data(valid_folder)

    allocate(train_data(training_data_points), valid_data(validation_data_points), stat=alloc_stat)
    if (alloc_stat /= 0) then
        print *, err // "allocation error"; stop
    end if

    do i = 1, training_data_points
        train_data(i) = new_atoms(1, reference%natoms, reference%ntypes)
    end do
    do i = 1, validation_data_points
        valid_data(i) = new_atoms(1, reference%natoms, reference%ntypes)
    end do

    call read_fitting_data(train_folder, train_data)
    call read_fitting_data(valid_folder, valid_data)

    call calc_force(reference, flag=energy_only)


end subroutine perform_fit




subroutine read_fitting_data(folder, data)

    character(len=*), intent(in)  :: folder
    type(universe), intent(inout) :: data(:)

    character(len=*), parameter :: err = "Error in read_fitting_data(): "
    character(len=max_string_length) :: words(100), buffer
    character(len=23) :: inp_file
    integer, parameter :: inp_unit = 67
    integer :: i, ios, nwords, nfiles, nrg_pos, pos_pos, atom_pos
    character(len=23) :: pos_file, nrg_file
    logical :: config_section


    ! select folder
    if (folder == train_folder) then
        nfiles = simparams%fit_training_data
    else if (folder == valid_folder) then
        nfiles = simparams%fit_validation_data
    else
        print *, err // "internal error"; stop
    end if

    nrg_pos = 1; pos_pos = 0
    config_section = .False.

    do i = 1, nfiles

        ! select folder
        if (folder == train_folder) then
            write(inp_file, '(a10, a4, i5.5, a4)') train_folder, "nrg_", i, ".dat"     ! fit/train/nrg_00581.dat
        else if (folder == valid_folder) then
            write(inp_file, '(a10, a4, i5.5, a4)') valid_folder, "nrg_", i, ".dat"     ! fit/valid/nrg_00581.dat
        end if

        ! read energy file, and write into universe array
        call open_for_read(inp_unit, inp_file)
        ios = 0
        do while (ios == 0)

            ! read line
            read(inp_unit, '(A)', iostat=ios) buffer
            if (ios /= 0) then
                exit    ! eof reached
            end if

            ! extract energy
            call split_string(buffer, words, nwords)
            if (folder == train_folder) then
                read(words(9), *, iostat=ios) data(nrg_pos)%epot
            else if (folder == valid_folder) then
                read(words(9), *, iostat=ios) data(nrg_pos)%epot
            end if
            nrg_pos = nrg_pos + 1
            if (ios /= 0) then
                print *, err // "reading energy in file", i, "in ", folder
                stop
            end if

        end do
        close(inp_unit)


        ! select folder
        if (folder == train_folder) then
            write(inp_file, '(a10, a4, i5.5, a4)') train_folder, "pos_", i, ".dat"     ! fit/train/pos_00581.dat
        else if (folder == valid_folder) then
            write(inp_file, '(a10, a4, i5.5, a4)') valid_folder, "pos_", i, ".dat"     ! fit/valid/pos_00581.dat
        end if

        ! read energy file, and write into universe array
        call open_for_read(inp_unit, inp_file)
        ios = 0
        do while (ios == 0)

            ! read line
            read(inp_unit, '(A)', iostat=ios) buffer
            if (ios /= 0) then
                exit    ! eof reached
            end if

            ! extract energy
            call split_string(buffer, words, nwords)
            if (.not. config_section) then
                if (words(1) == "Direct" .and. words(2) == "configuration=") then
                    config_section = .True.
                else
                    continue
                end if
            end if

            ! only when in config section
            if (words(1) == "Direct" .and. words(2) == "configuration=") then
                pos_pos = pos_pos + 1
                atom_pos = 1
            else
                if (folder == train_folder) then
                    read(words(9), *, iostat=ios) data(pos_pos)%r(:,1,atom_pos)
                else if (folder == valid_folder) then
                    read(words(9), *, iostat=ios) data(pos_pos)%r(:,1,atom_pos)
                end if
                atom_pos = atom_pos + 1
            end if
            if (ios /= 0) then
                print *, err // "reading energy in file", i, "in ", folder
                stop
            end if

        end do
        close(inp_unit)




    end do





end subroutine read_fitting_data




subroutine check_fitting_data(folder)

    character(len=*), intent(in) :: folder

    character(len=23) :: inp_file
    character(len=max_string_length) :: buffer, words(100)
    character(len=*), parameter :: err = "Error in check_fitting_data(): "
    integer, parameter :: inp_unit = 67
    integer :: ndata_points ! counts all config/energy pairs
    integer :: this_data_points ! counts config/energy pairs in one file
    integer :: i, ios, nwords, nfiles

    ! select folder
    if (folder == train_folder) then
        nfiles = simparams%fit_training_data
    else if (folder == valid_folder) then
        nfiles = simparams%fit_validation_data
    else
        print *, err // "internal error"; stop
    end if

    do i = 1, nfiles

        ! select folder
        if (folder == train_folder) then
            write(inp_file, '(a10, a4, i5.5, a4)') train_folder, "nrg_", i, ".dat"     ! fit/train/nrg_00581.dat
        else if (folder == valid_folder) then
            write(inp_file, '(a10, a4, i5.5, a4)') valid_folder, "nrg_", i, ".dat"     ! fit/valid/nrg_00581.dat
        else
            print *, err // "internal error"; stop
        end if

        ! scan energy file, count lines in file and add to module variable
        call open_for_read(inp_unit, inp_file)
        ios = 0
        this_data_points = -1   ! increment also happens if ios /= 0
        do while (ios == 0)
            read(inp_unit, '(A)', iostat=ios) buffer
            this_data_points = this_data_points + 1
        end do
        close(inp_unit)

        ! select folder
        if (folder == train_folder) then
            training_data_points = training_data_points + this_data_points
            write(inp_file, '(a10, a4, i5.5, a4)') train_folder, "pos_", i, ".dat"     ! fit/train/pos_00581.dat
        else if (folder == valid_folder) then
            validation_data_points = validation_data_points + this_data_points
            write(inp_file, '(a10, a4, i5.5, a4)') valid_folder, "pos_", i, ".dat"     ! fit/valid/pos_00581.dat
        end if

        ! scan configuration file
        call open_for_read(inp_unit, inp_file)
        ios = 0
        do while (ios == 0)
            read(inp_unit, '(A)', iostat=ios) buffer
            call split_string(buffer, words, nwords)
            if (words(1) == "Direct" .and. words(2) == "configuration=") then
                ! subtract one each found time step
                this_data_points = this_data_points - 1
            end if
        end do
        close(inp_unit)
        ! perform sanity check
        if (this_data_points /= 0) then
            print *, err // "number of energy/configuration pairs does not match for file", i, "in ", folder
            stop
        end if

    end do

end subroutine check_fitting_data



subroutine nllsqwbc
    use MKL_RCI
    use MKL_RCI_TYPE
    implicit none
!** HEADER-FILE WITH DEFINITIONS (CONSTANTS, EXTERNALS)

!** USER'S OBJECTIVE FUNCTION
    external            :: EXTENDED_POWELL
!** N - NUMBER OF FUNCTION VARIABLES
    integer, parameter  :: N = 40
!** M - DIMENSION OF FUNCTION VALUE
    integer, parameter  :: M = 40
!** SOLUTION VECTOR. CONTAINS VALUES X FOR F(X)
    double precision    :: X (N)
!** PRECISIONS FOR STOP-CRITERIA (SEE MANUAL FOR MORE DETAILS)
    double precision    :: EPS (6)
!** JACOBI CALCULATION PRECISION
    double precision    :: JAC_EPS
!** REVERSE COMMUNICATION INTERFACE PARAMETER
    integer             :: RCI_REQUEST
!** FUNCTION (F(X)) VALUE VECTOR
    double precision    :: FVEC (M)
!** JACOBI MATRIX
    double precision    :: FJAC (M, N)
!** NUMBER OF ITERATIONS
    integer             :: ITER
!** NUMBER OF STOP-CRITERION
    integer             :: ST_CR
!** CONTROLS OF RCI CYCLE
    integer             :: SUCCESSFUL
!** MAXIMUM NUMBER OF ITERATIONS
    integer             :: ITER1
!** MAXIMUM NUMBER OF ITERATIONS OF CALCULATION OF TRIAL-STEP
    integer             :: ITER2
!** INITIAL STEP BOUND
    double precision    :: RS
!** INITIAL AND FINAL RESIDUALS
    double precision    :: R1, R2
!** TR SOLVER HANDLE
    TYPE(HANDLE_TR) :: HANDLE
!** CYCLES COUNTERS
    integer             :: I, J
!** RESULTS OF INPUT PARAMETER CHECKING
    integer :: INFO(6)
!** SET PRECISIONS FOR STOP-CRITERIA
    EPS  = 1.D-5
!** SET MAXIMUM NUMBER OF ITERATIONS
    ITER1 = 1000
!** SET MAXIMUM NUMBER OF ITERATIONS OF CALCULATION OF TRIAL-STEP
    ITER2 = 100
!** SET INITIAL STEP BOUND
    RS = 100.D0
!** PRECISIONS FOR JACOBI CALCULATION
    JAC_EPS = 1.D-8
!** SET THE INITIAL GUESS
    do I = 1, N/4
        X (4*I - 3) =  3.D0
        X (4*I - 2) = -1.D0
        X (4*I - 1) =  0.D0
        X (4*I)     =  1.D0
    end do
!** SET INITIAL VALUES
    do I = 1, M
        FVEC (I) = 0.D0
        do J = 1, N
            FJAC (I, J) = 0.D0
        end do
    end do
!** INITIALIZE SOLVER (ALLOCATE MEMORY, SET INITIAL VALUES)
!**     HANDLE    IN/OUT: TR SOLVER HANDLE
!**     N         IN:     NUMBER OF FUNCTION VARIABLES
!**     M         IN:     DIMENSION OF FUNCTION VALUE
!**     X         IN:     SOLUTION VECTOR. CONTAINS VALUES X FOR F(X)
!**     EPS       IN:     PRECISIONS FOR STOP-CRITERIA
!**     ITER1     IN:     MAXIMUM NUMBER OF ITERATIONS
!**     ITER2     IN:     MAXIMUM NUMBER OF ITERATIONS OF CALCULATION OF TRIAL-STEP
!**     RS        IN:     INITIAL STEP BOUND
    if (DTRNLSP_INIT (HANDLE, N, M, X, EPS, ITER1, ITER2, RS) /= TR_SUCCESS) THEN
!** IF FUNCTION DOES NOT COMPLETE SUCCESSFULLY THEN PRINT ERROR MESSAGE
        print *, '| ERROR IN DTRNLSP_INIT'
!** RELEASE INTERNAL Intel(R) MKL MEMORY THAT MIGHT BE USED FOR COMPUTATIONS.
!** NOTE: IT IS IMPORTANT TO CALL THE ROUTINE BELOW TO AVOID MEMORY LEAKS
!** UNLESS YOU DISABLE Intel(R) MKL MEMORY MANAGER
        call MKL_FREE_BUFFERS
!** AND STOP
        stop 1
    end if
!** CHECKS THE CORRECTNESS OF HANDLE AND ARRAYS CONTAINING JACOBIAN MATRIX,
!** OBJECTIVE FUNCTION, LOWER AND UPPER BOUNDS, AND STOPPING CRITERIA.
    if (DTRNLSP_CHECK (HANDLE, N, M, FJAC, FVEC, EPS, INFO) /= TR_SUCCESS) THEN
!** IF FUNCTION DOES NOT COMPLETE SUCCESSFULLY THEN PRINT ERROR MESSAGE
        print *, '| ERROR IN DTRNLSPBC_INIT'
!** RELEASE INTERNAL Intel(R) MKL MEMORY THAT MIGHT BE USED FOR COMPUTATIONS.
!** NOTE: IT IS IMPORTANT TO CALL THE ROUTINE BELOW TO AVOID MEMORY LEAKS
!** UNLESS YOU DISABLE Intel(R) MKL MEMORY MANAGER
        call MKL_FREE_BUFFERS
!** AND STOP
        stop 1
    else
!** THE HANDLE IS NOT VALID.
        if ( INFO(1) /= 0 .or. INFO(2) /= 0 .or. INFO(3) /= 0 .or. INFO(4) /= 0 ) THEN
!** THE FJAC ARRAY IS NOT VALID.
!     +     INFO(2) /= 0 .or.
!** THE FVEC ARRAY IS NOT VALID.
!     +     INFO(3) /= 0 .or.
!** THE EPS ARRAY IS NOT VALID.
!     +     INFO(4) /= 0 ) THEN
            print *, '| INPUT PARAMETERS ARE NOT VALID'
!** RELEASE INTERNAL Intel(R) MKL MEMORY THAT MIGHT BE USED FOR COMPUTATIONS.
!** NOTE: IT IS IMPORTANT TO CALL THE ROUTINE BELOW TO AVOID MEMORY LEAKS
!** UNLESS YOU DISABLE Intel(R) MKL MEMORY MANAGER
            call MKL_FREE_BUFFERS
!** AND STOP
            stop 1
        end if
    end if
!** SET INITIAL RCI CYCLE VARIABLES
    RCI_REQUEST = 0
    SUCCESSFUL = 0
!** RCI CYCLE
    do while (SUCCESSFUL == 0)
!** CALL TR SOLVER
!**   HANDLE        IN/OUT: TR SOLVER HANDLE
!**   FVEC          IN:     VECTOR
!**   FJAC          IN:     JACOBI MATRIX
!**   RCI_REQUEST   IN/OUT: RETURN NUMBER WHICH DENOTE NEXT STEP FOR PERFORMING
        if (DTRNLSP_SOLVE (HANDLE, FVEC, FJAC, RCI_REQUEST) /= TR_SUCCESS) THEN
!** IF FUNCTION DOES NOT COMPLETE SUCCESSFULLY THEN PRINT ERROR MESSAGE
            print *, '| ERROR IN DTRNLSP_SOLVE'
!** RELEASE INTERNAL Intel(R) MKL MEMORY THAT MIGHT BE USED FOR COMPUTATIONS.
!** NOTE: IT IS IMPORTANT TO CALL THE ROUTINE BELOW TO AVOID MEMORY LEAKS
!** UNLESS YOU DISABLE Intel(R) MKL MEMORY MANAGER
            call MKL_FREE_BUFFERS
!** AND STOP
            stop 1
        end if
!** ACCORDING WITH RCI_REQUEST VALUE WE DO NEXT STEP
        select case (RCI_REQUEST)
        case (-1, -2, -3, -4, -5, -6)
!**   STOP RCI CYCLE
            SUCCESSFUL = 1
        case (1)
!**   RECALCULATE FUNCTION VALUE
!**     M               IN:     DIMENSION OF FUNCTION VALUE
!**     N               IN:     NUMBER OF FUNCTION VARIABLES
!**     X               IN:     SOLUTION VECTOR
!**     FVEC            OUT:    FUNCTION VALUE F(X)
            call EXTENDED_POWELL (M, N, X, FVEC)
        case (2)
!**   COMPUTE JACOBI MATRIX
!**     EXTENDED_POWELL IN:     EXTERNAL OBJECTIVE FUNCTION
!**     N               IN:     NUMBER OF FUNCTION VARIABLES
!**     M               IN:     DIMENSION OF FUNCTION VALUE
!**     FJAC            OUT:    JACOBI MATRIX
!**     X               IN:     SOLUTION VECTOR
!**     JAC_EPS         IN:     JACOBI CALCULATION PRECISION
            if (DJACOBI (EXTENDED_POWELL, N, M, FJAC, X, JAC_EPS) /= TR_SUCCESS) THEN
!** IF FUNCTION DOES NOT COMPLETE SUCCESSFULLY THEN PRINT ERROR MESSAGE
                print *, '| ERROR IN DJACOBI'
!** RELEASE INTERNAL Intel(R) MKL MEMORY THAT MIGHT BE USED FOR COMPUTATIONS.
!** NOTE: IT IS IMPORTANT TO CALL THE ROUTINE BELOW TO AVOID MEMORY LEAKS
!** UNLESS YOU DISABLE Intel(R) MKL MEMORY MANAGER
                call MKL_FREE_BUFFERS
!** AND STOP
                stop 1
            end if
        end select
    end do
!** GET SOLUTION STATUSES
!**   HANDLE            IN: TR SOLVER HANDLE
!**   ITER              OUT: NUMBER OF ITERATIONS
!**   ST_CR             OUT: NUMBER OF STOP CRITERION
!**   R1                OUT: INITIAL RESIDUALS
!**   R2                OUT: FINAL RESIDUALS
    if (DTRNLSP_GET (HANDLE, ITER, ST_CR, R1, R2) /= TR_SUCCESS) THEN
!** IF FUNCTION DOES NOT COMPLETE SUCCESSFULLY THEN PRINT ERROR MESSAGE
        print *, '| ERROR IN DTRNLSP_GET'
!** RELEASE INTERNAL Intel(R) MKL MEMORY THAT MIGHT BE USED FOR COMPUTATIONS.
!** NOTE: IT IS IMPORTANT TO CALL THE ROUTINE BELOW TO AVOID MEMORY LEAKS
!** UNLESS YOU DISABLE Intel(R) MKL MEMORY MANAGER
        call MKL_FREE_BUFFERS
!** AND STOP
        stop 1
    end if
!** FREE HANDLE MEMORY
    if (DTRNLSP_DELETE (HANDLE) /= TR_SUCCESS) then
!** IF FUNCTION DOES NOT COMPLETE SUCCESSFULLY THEN PRINT ERROR MESSAGE
        print *, '| ERROR IN DTRNLSP_DELETE'
!** RELEASE INTERNAL Intel(R) MKL MEMORY THAT MIGHT BE USED FOR COMPUTATIONS.
!** NOTE: IT IS IMPORTANT TO CALL THE ROUTINE BELOW TO AVOID MEMORY LEAKS
!** UNLESS YOU DISABLE Intel(R) MKL MEMORY MANAGER
        call MKL_FREE_BUFFERS
!** AND STOP
        stop 1
    end if

!** RELEASE INTERNAL Intel(R) MKL MEMORY THAT MIGHT BE USED FOR COMPUTATIONS.
!** NOTE: IT IS IMPORTANT TO CALL THE ROUTINE BELOW TO AVOID MEMORY LEAKS
!** UNLESS YOU DISABLE Intel(R) MKL MEMORY MANAGER
    call MKL_FREE_BUFFERS
!** IF FINAL RESIDUAL LESS THEN REQUIRED PRECISION THEN PRINT PASS
    if (R2 < 1.D-5) then
        print *, '|         DTRNLSP POWELL............PASS'
        stop 0
!** ELSE PRINT FAILED
    else
        print *, '|         DTRNLSP POWELL............FAILED'
        stop 1
    end if

end subroutine nllsqwbc

end module fit

!** ROUTINE FOR EXTENDED POWELL FUNCTION CALCULATION
!**   M     IN:     DIMENSION OF FUNCTION VALUE
!**   N     IN:     NUMBER OF FUNCTION VARIABLES
!**   X     IN:     VECTOR FOR FUNCTION CALCULATING
!**   F     OUT:    FUNCTION VALUE F(X)
subroutine EXTENDED_POWELL (M, N, X, F)
    implicit none
    integer M, N
    double precision X (*), F (*)
    integer I

    do I = 1, N/4
        F (4*I-3) = X(4*I - 3) + 10.D0 * X(4*I - 2)
        F (4*I-2) = 2.2360679774998D0 * (X(4*I-1) - X(4*I))
        F (4*I-1) = ( X(4*I-2) - 2.D0*X(4*I-1) )**2
        F (4*I)   = 3.1622776601684D0 * (X(4*I-3) - X(4*I))**2
    end do
end subroutine EXTENDED_POWELL
