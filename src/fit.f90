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

        ! check the data for consistency, count training and validation data sets to
        ! allocate storage of array of universe objects
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

        ! # fitting params

        ! do the fit
        call nllsqwbc(train_data, valid_data)

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

        nrg_pos = 0; pos_pos = 0

        do i = 1, nfiles

            !print *, "reading", folder, i

            ! select folder
            write(inp_file, '(a10, a4, i5.5, a4)') folder, "nrg_", i, ".dat"     ! ex: fit/train/nrg_00581.dat

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
                nrg_pos = nrg_pos + 1
                read(words(9), *, iostat=ios) data(nrg_pos)%epot
                if (ios /= 0) then
                    print *, err // "reading energy in file", i, "in ", folder
                    stop
                end if

            end do
            close(inp_unit)


            ! select folder
            write(inp_file, '(a10, a4, i5.5, a4)') folder, "pos_", i, ".dat"     ! ex: fit/train/pos_00581.dat

            ! read position file, and write into universe array
            call open_for_read(inp_unit, inp_file)
            ios = 0
            config_section = .False.
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
                        cycle
                    end if
                end if

                ! only when in config section
                if (words(1) == "Direct" .and. words(2) == "configuration=") then
                    pos_pos = pos_pos + 1
                    atom_pos = 0
                else
                    atom_pos = atom_pos + 1
                    read(words(1:3), *, iostat=ios) data(pos_pos)%r(:,1,atom_pos)
                end if
                if (ios /= 0) then
                    print *, err // "reading energy in file", i, "in ", folder
                    stop
                end if

            end do
            close(inp_unit)

            if (pos_pos /= nrg_pos) then
                print *, err // "discrepancy between energies and positions in file ", i, " in ", folder
                stop
            else if (atom_pos /= data(i)%natoms) then
                print *, err, atom_pos, " atoms found in file ", i, " in ", folder, " but data structure requires", data(i)%natoms
                stop
            end if

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

            !print *, "checking", folder, i

            ! select folder
            write(inp_file, '(a10, a4, i5.5, a4)') folder, "nrg_", i, ".dat"     ! ex: fit/train/nrg_00581.dat

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
            else if (folder == valid_folder) then
                validation_data_points = validation_data_points + this_data_points
            end if

            write(inp_file, '(a10, a4, i5.5, a4)') folder, "pos_", i, ".dat"     ! ex: fit/valid/pos_00581.dat


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



    function from_pes_params_to_x(ref) result(x)

        type(universe), intent(in) :: ref
        real(dp), allocatable :: x(:)

        ! TODO: add other PESs

        if (any(ref%pes == pes_id_rebo)) call inquire_nfit_params_rebo()




    end function from_pes_params_to_x





    subroutine nllsqwbc(td, vd)
        use mkl_rci
        use mkl_rci_type
        implicit none

        type(universe), intent(inout) :: td(:), vd(:)

        real(dp) :: train_nrg(training_data_points)
        real(dp) :: valid_nrg(validation_data_points)

        !** user's objective function
        external            :: obj_func
        !** n - number of function variables
        integer  :: n
        !** m - dimension of function value
        integer  :: m
        !** solution vector. contains flattened potential parameters
        real(dp), allocatable :: x(:)
        !** precisions for stop-criteria (see manual for more details)
        real(dp)            :: eps(6)
        !** jacobi calculation precision
        real(dp)            :: jac_eps
        !** reverse communication interface parameter
        integer             :: rci_request
        !** function (f(x)) value vector
        real(dp), allocatable :: fvec(:)
        !** jacobi matrix
        real(dp), allocatable :: fjac (:,:)
        !** number of iterations
        integer             :: iter
        !** number of stop-criterion
        integer             :: st_cr
        !** controls of rci cycle
        integer             :: successful
        !** maximum number of iterations
        integer             :: iter1
        !** maximum number of iterations of calculation of trial-step
        integer             :: iter2
        !** initial step bound
        real(dp)            :: rs
        !** initial and final residuals
        real(dp)            :: r1, r2
        !** tr solver handle
        type(handle_tr) :: handle
        !** cycles counters
        integer             :: i, j
        !** results of input parameter checking
        integer :: info(6)
        !** set precisions for stop-criteria
        eps  = 1.0e-5_dp
        !** set maximum number of iterations
        iter1 = 1000
        !** set maximum number of iterations of calculation of trial-step
        iter2 = 100
        !** set initial step bound
        rs = 100.0_dp
        !** precisions for jacobi calculation
        jac_eps = 1.0e-8_dp

        m = training_data_points



        !** set the initial guess
        do i = 1, n/4
            x (4*i - 3) =  3.d0
            x (4*i - 2) = -1.d0
            x (4*i - 1) =  0.d0
            x (4*i)     =  1.d0
        end do
        !** set initial values
        fvec = 0.0_dp
        fjac = 0.0_dp

        !** initialize solver (allocate memory, set initial values)
        !**     handle    in/out: tr solver handle
        !**     n         in:     number of function variables
        !**     m         in:     dimension of function value
        !**     x         in:     solution vector. contains values x for f(x)
        !**     eps       in:     precisions for stop-criteria
        !**     iter1     in:     maximum number of iterations
        !**     iter2     in:     maximum number of iterations of calculation of trial-step
        !**     rs        in:     initial step bound
        if (dtrnlsp_init (handle, n, m, x, eps, iter1, iter2, rs) /= tr_success) then
            !** if function does not complete successfully then print error message
            print *, '| error in dtrnlsp_init'
            !** release internal intel(r) mkl memory that might be used for computations.
            !** note: it is important to call the routine below to avoid memory leaks
            !** unless you disable intel(r) mkl memory manager
            call mkl_free_buffers
            !** and stop
            stop 1
        end if
        !** checks the correctness of handle and arrays containing jacobian matrix,
        !** objective function, lower and upper bounds, and stopping criteria.
        if (dtrnlsp_check (handle, n, m, fjac, fvec, eps, info) /= tr_success) then
            !** if function does not complete successfully then print error message
            print *, '| error in dtrnlspbc_init'
            !** release internal intel(r) mkl memory that might be used for computations.
            !** note: it is important to call the routine below to avoid memory leaks
            !** unless you disable intel(r) mkl memory manager
            call mkl_free_buffers
            !** and stop
            stop 1
        else
            !** the handle is not valid.
            if ( info(1) /= 0 .or. info(2) /= 0 .or. info(3) /= 0 .or. info(4) /= 0 ) then
                !** the fjac array is not valid.
                !     +     info(2) /= 0 .or.
                !** the fvec array is not valid.
                !     +     info(3) /= 0 .or.
                !** the eps array is not valid.
                !     +     info(4) /= 0 ) then
                print *, '| input parameters are not valid'
                !** release internal intel(r) mkl memory that might be used for computations.
                !** note: it is important to call the routine below to avoid memory leaks
                !** unless you disable intel(r) mkl memory manager
                call mkl_free_buffers
                !** and stop
                stop 1
            end if
        end if
        !** set initial rci cycle variables
        rci_request = 0
        successful = 0
        !** rci cycle
        do while (successful == 0)
            !** call tr solver
            !**   handle        in/out: tr solver handle
            !**   fvec          in:     vector
            !**   fjac          in:     jacobi matrix
            !**   rci_request   in/out: return number which denote next step for performing
            if (dtrnlsp_solve (handle, fvec, fjac, rci_request) /= tr_success) then
                !** if function does not complete successfully then print error message
                print *, '| error in dtrnlsp_solve'
                !** release internal intel(r) mkl memory that might be used for computations.
                !** note: it is important to call the routine below to avoid memory leaks
                !** unless you disable intel(r) mkl memory manager
                call mkl_free_buffers
                !** and stop
                stop 1
            end if
            !** according with rci_request value we do next step
            select case (rci_request)
                case (-1, -2, -3, -4, -5, -6)
                    !**   stop rci cycle
                    successful = 1
                case (1)
                    !**   recalculate function value
                    !**     m               in:     dimension of function value
                    !**     n               in:     number of function variables
                    !**     x               in:     solution vector
                    !**     fvec            out:    function value f(x)
                    call obj_func (m, n, x, fvec)
                case (2)
                    !**   compute jacobi matrix
                    !**     extended_powell in:     external objective function
                    !**     n               in:     number of function variables
                    !**     m               in:     dimension of function value
                    !**     fjac            out:    jacobi matrix
                    !**     x               in:     solution vector
                    !**     jac_eps         in:     jacobi calculation precision
                    if (djacobi (obj_func, n, m, fjac, x, jac_eps) /= tr_success) then
                        !** if function does not complete successfully then print error message
                        print *, '| error in djacobi'
                        !** release internal intel(r) mkl memory that might be used for computations.
                        !** note: it is important to call the routine below to avoid memory leaks
                        !** unless you disable intel(r) mkl memory manager
                        call mkl_free_buffers
                        !** and stop
                        stop 1
                    end if
            end select
        end do
        !** get solution statuses
        !**   handle            in: tr solver handle
        !**   iter              out: number of iterations
        !**   st_cr             out: number of stop criterion
        !**   r1                out: initial residuals
        !**   r2                out: final residuals
        if (dtrnlsp_get (handle, iter, st_cr, r1, r2) /= tr_success) then
            !** if function does not complete successfully then print error message
            print *, '| error in dtrnlsp_get'
            !** release internal intel(r) mkl memory that might be used for computations.
            !** note: it is important to call the routine below to avoid memory leaks
            !** unless you disable intel(r) mkl memory manager
            call mkl_free_buffers
            !** and stop
            stop 1
        end if
        !** free handle memory
        if (dtrnlsp_delete (handle) /= tr_success) then
            !** if function does not complete successfully then print error message
            print *, '| error in dtrnlsp_delete'
            !** release internal intel(r) mkl memory that might be used for computations.
            !** note: it is important to call the routine below to avoid memory leaks
            !** unless you disable intel(r) mkl memory manager
            call mkl_free_buffers
            !** and stop
            stop 1
        end if

        !** release internal intel(r) mkl memory that might be used for computations.
        !** note: it is important to call the routine below to avoid memory leaks
        !** unless you disable intel(r) mkl memory manager
        call mkl_free_buffers
        !** if final residual less then required precision then print pass
        if (r2 < 1.d-5) then
            print *, '|         dtrnlsp powell............pass'
            stop 0
        !** else print failed
        else
            print *, '|         dtrnlsp powell............failed'
            stop 1
        end if

    end subroutine nllsqwbc

end module fit

!** ROUTINE FOR EXTENDED POWELL FUNCTION CALCULATION
!**   M     IN:     DIMENSION OF FUNCTION VALUE
!**   N     IN:     NUMBER OF FUNCTION VARIABLES
!**   X     IN:     VECTOR FOR FUNCTION CALCULATING
!**   F     OUT:    FUNCTION VALUE F(X)
subroutine obj_func (M, N, X, F)
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
end subroutine obj_func
