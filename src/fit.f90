!** NONLINEAR LEAST SQUARE PROBLEM WITHOUT BOUNDARY CONSTRAINTS
    include 'mkl_rci.f90'
module fit

    use force
    use universe_mod
    use run_config, only : simparams
    use open_file
    use useful_things
    use pes_rebo_mod

    implicit none

    character(len=*), parameter :: train_folder = "fit/train/"
    character(len=*), parameter :: valid_folder = "fit/valid/"
    integer :: training_data_points, validation_data_points
    real(dp), allocatable :: train_nrg(:), valid_nrg(:)
    type(universe), allocatable :: train_data(:), valid_data(:)
    integer :: nfit_params

contains

    subroutine perform_fit(reference)

        type(universe), intent(inout) :: reference     ! holds the reference configuration


        real(dp) :: train_rmse, valid_rmse
        integer :: alloc_stat, i, values(8)
        real(dp), allocatable :: x(:)
        real(dp), allocatable :: train_dev(:), valid_dev(:)
        character(len=*), parameter :: err = "Error in perform_fit(): "


        training_data_points   = 0
        validation_data_points = 0

        ! check the data for consistency, count training and validation data sets to
        ! allocate storage of array of universe objects
        call check_fitting_data(train_folder)
        call check_fitting_data(valid_folder)

        ! 1st entry is the reference minimum energy structure
        allocate(train_data(training_data_points+1), &
                 valid_data(validation_data_points), &
                 train_nrg(training_data_points+1), &
                 valid_nrg(validation_data_points), &
                 train_dev(training_data_points+1), &
                 valid_dev(validation_data_points), stat=alloc_stat)
        if (alloc_stat /= 0) then
            print *, err // "allocation error"; stop
        end if

        ! build training and validation data sets
        train_data(1) = reference
        train_data(1)%epot(1) = simparams%evasp
        do i = 1, training_data_points
            train_data(i+1)     = new_atoms(1, reference%natoms, reference%ntypes)
            train_data(i+1)%pes = reference%pes
            train_data(i+1)%idx = reference%idx
            train_data(i+1)%name = reference%name
            train_data(i+1)%simbox = reference%simbox
            train_data(i+1)%isimbox = reference%isimbox
        end do
        do i = 1, validation_data_points
            valid_data(i) = new_atoms(1, reference%natoms, reference%ntypes)
            valid_data(i)%pes = reference%pes
            valid_data(i)%idx = reference%idx
            valid_data(i)%name = reference%name
            valid_data(i)%simbox = reference%simbox
            valid_data(i)%isimbox = reference%isimbox
        end do

        call read_fitting_data(train_folder, train_data)
        call read_fitting_data(valid_folder, valid_data)


        ! subtract DFT reference energy
        do i = 1, training_data_points + 1
            train_nrg(i) = train_data(i)%epot(1) - simparams%evasp
        end do
        do i = 1, validation_data_points
            valid_nrg(i) = valid_data(i)%epot(1) - simparams%evasp
        end do

        x = from_pes_params_to_x(reference)



        ! do the fit
        call nllsqwbc(train_data, valid_data, x, train_dev, train_rmse)

        valid_dev = default_real
        do i = 1, validation_data_points
            print *, i
            call calc_force(valid_data(i), energy_only)
            print *, valid_nrg(i), valid_data(i)%epot(1)-train_data(1)%epot(1)
            valid_dev(i) = valid_nrg(i) - (valid_data(i)%epot(1)-train_data(1)%epot(1))
        end do

        ! evaluate validation data set
        valid_rmse = sqrt(sum(valid_dev*valid_dev)/validation_data_points)
        call date_and_time(VALUES=values)
        print '(i4, a, i2.2, a, i2.2, a, i2.2, a, i2.2, a5, f10.5, a12, f10.5)', &
            values(1), "-", values(2), "-",values(3), " - ",values(5), &
            ".",values(6), "max:", maxval(abs(valid_dev)), "valid rmse:", valid_rmse

    end subroutine perform_fit




    subroutine read_fitting_data(folder, data)

        character(len=*), intent(in)  :: folder
        type(universe), intent(inout) :: data(:)

        character(len=*), parameter :: err = "Error in read_fitting_data(): "
        character(len=max_string_length) :: words(100), buffer
        character(len=23) :: inp_file
        integer, parameter :: inp_unit = 67
        integer :: i, ios, nwords, nfiles, nrg_pos, pos_pos, atom_pos
        logical :: config_section


        ! select folder
        if (folder == train_folder) then
            nfiles = simparams%fit_training_data
            nrg_pos = 1; pos_pos = 1
        else if (folder == valid_folder) then
            nfiles = simparams%fit_validation_data
            nrg_pos = 0; pos_pos = 0
        else
            print *, err // "internal error"; stop
        end if

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
                read(words(9), *, iostat=ios) data(nrg_pos)%epot(1)
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

        real(dp), allocatable :: params_rebo(:)


        ! TODO: add other PESs
        if (any(ref%pes == pes_id_rebo)) params_rebo = get_fit_params_rebo()

        ! add other sizes when implementing other PESs
        nfit_params = size([params_rebo])

        allocate(x(nfit_params))

        ! add other PES param array when implementing other PESs
        x = [params_rebo]

    end function from_pes_params_to_x



    subroutine from_x_to_pes_params(x)

        real(dp), intent(in)       :: x(:)

        integer :: pos

        pos = 1 ! counts the position in the x array

        if (REBO_INIT_SUCCESS) call set_fit_params_rebo(x, pos)

    end subroutine from_x_to_pes_params


    subroutine nllsqwbc(td, vd, x, fvec, train_rmse)
        use mkl_rci
        use mkl_rci_type
        implicit none

        !** traning and validation data
        type(universe), intent(inout) :: td(:), vd(:)

        !** solution vector. contains flattened potential parameters
        real(dp), intent(inout) :: x(:)

        !** vector containing deviations from target energies
        real(dp), intent(out) :: fvec(:), train_rmse

        !** user's objective function
        !external            :: obj_func
        !** n - number of function variables
        integer  :: n
        !** m - dimension of function value
        integer  :: m
        !** precisions for stop-criteria (see manual for more details)
        real(dp)            :: eps(6)
        !** jacobi calculation precision
        real(dp)            :: jac_eps
        !** reverse communication interface parameter
        integer             :: rci_request
        !** jacobi matrix
        real(dp), allocatable :: fjac (:,:)
        !** number of iterations
        integer             :: iter = 0
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
        integer :: values(8)
        real(dp) :: max_val
        !** set precisions for stop-criteria
        eps  = 1.0e-5_dp
        !** set maximum number of iterations
        iter1 = simparams%maxit
        !** set maximum number of iterations of calculation of trial-step
        iter2 = 100
        !** set initial step bound
        rs = 100.0_dp
        !** precisions for jacobi calculation
        jac_eps = 1.0e-8_dp

        m = training_data_points+1
        n = size(x)

        allocate(fjac(m,n))

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
        !print *, "dtrnlsp_init success"
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
        !print *, "dtrnlsp_check success"
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
            !print *, "starting dtrnlsp_solve"
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
                    !print *, "RCI request = 1: calculate obj_func"
                    call obj_func (m, n, x, fvec)
                case (2)
                    !**   compute jacobi matrix
                    !**     extended_powell in:     external objective function
                    !**     n               in:     number of function variables
                    !**     m               in:     dimension of function value
                    !**     fjac            out:    jacobi matrix
                    !**     x               in:     solution vector
                    !**     jac_eps         in:     jacobi calculation precision
                    !print *, "RCI request = 2: calculate djacobi"
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

                    do i=2,M
                        if (mod(i, M/10) == 0) then
                            write(*,"(2i6, 2f12.4)") iter, i-1, train_nrg(i), train_nrg(i)-fvec(i)
                        end if
                    end do

                    train_rmse = sqrt(sum(fvec(2:)*fvec(2:))/training_data_points)
                    max_val = maxval(abs(fvec)) ! max deviation

                    call date_and_time(VALUES=values)

                    print '(i4, a, i2.2, a, i2.2, a, i2.2, a, i2.2, a6, i6, a5, f10.5, a12, f10.5)', &
                        values(1), "-", values(2), "-",values(3), " - ",values(5), &
                        ".",values(6), "iter:", iter, "max:", max_val, "train rmse:", train_rmse

                    iter = iter + 1
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

        call from_x_to_pes_params(x)

        call mkl_free_buffers

        !print *, 'initial residual', r1
        !print *, 'final residual', r2

    end subroutine nllsqwbc


    !**   M     IN:     DIMENSION OF FUNCTION VALUE
    !**   N     IN:     NUMBER OF FUNCTION VARIABLES
    !**   X     IN:     VECTOR FOR FUNCTION CALCULATING
    !**   F     OUT:    FUNCTION VALUE F(X)
    subroutine obj_func (m, n, x, f)

        implicit none
        integer,  intent(in)  :: m, n
        real(dp), intent(in)  :: x(n)
        real(dp), intent(out) :: f(m)

        integer :: i

        call from_x_to_pes_params(x)

        ! Eref is train_data(1)%epot(1)
        do i = 1, m
            call calc_force(train_data(i), energy_only)
            f(i) = train_nrg(i) - (train_data(i)%epot(1)-train_data(1)%epot(1))
        end do

    end subroutine obj_func

end module fit


