!############################################################################
! This routine is part of
! md_tian2 (Molecular Dynamics Tian Xia 2)
! (c) 2014-2021 Daniel J. Auerbach, Svenja M. Janke, Marvin Kammler,
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

!** NONLINEAR LEAST SQUARE PROBLEM WITHOUT BOUNDARY CONSTRAINTS
    include 'mkl_rci.f90'
    include 'mkl_service.f90'
    include 'omp_lib.f90'

module fit

    use force
    use universe_mod
    use run_config, only : simparams
    use open_file
    use useful_things
    use pes_rebo_mod
    use output_mod

    implicit none

    integer, parameter :: fit_out = 89
    integer :: training_data_points, validation_data_points
    real(dp), allocatable :: train_nrg(:), valid_nrg(:)
    type(universe), allocatable :: train_data(:), valid_data(:)
    type(universe)              :: reference
    integer :: nfit_params

contains

    subroutine perform_fit(dummy)

        type(universe), intent(in) :: dummy     ! holds the reference configuration

        real(dp) :: train_rmse, valid_rmse
        integer :: alloc_stat, i
        real(dp), allocatable :: x(:)
        real(dp), allocatable :: train_dev(:), valid_dev(:)
        character(len=*), parameter :: err = "Error in perform_fit(): "

        ! make dummmy available to module
        reference = dummy

        ! set up folder structure and open file handles
        call fit_setup()

        training_data_points   = 0
        validation_data_points = 0

        ! check the data for consistency, count training and validation data sets to
        ! allocate storage of array of universe objects
        call check_fitting_data(simparams%fit_training_folder)
        call timestamp(fit_out)
        write(fit_out, '(a, i)') "training data points found:  ", training_data_points

        call check_fitting_data(simparams%fit_validation_folder)
        call timestamp(fit_out)
        write(fit_out, '(a, i)') "validation data points found:", validation_data_points

        ! 1st entry is the reference minimum energy structure
        allocate(train_data(training_data_points), &
                 valid_data(validation_data_points), &
                 train_nrg(training_data_points), &
                 valid_nrg(validation_data_points), &
                 train_dev(training_data_points), &
                 valid_dev(validation_data_points), stat=alloc_stat)
        if (alloc_stat /= 0) then
            print *, err // "allocation error"; stop
        end if

        ! build training and validation data sets
        do i = 1, training_data_points
            train_data(i)     = new_atoms(1, reference%natoms, reference%ntypes)
            train_data(i)%pes = reference%pes
            train_data(i)%idx = reference%idx
            train_data(i)%name = reference%name
            train_data(i)%simbox = reference%simbox
            train_data(i)%isimbox = reference%isimbox
        end do
        do i = 1, validation_data_points
            valid_data(i) = new_atoms(1, reference%natoms, reference%ntypes)
            valid_data(i)%pes = reference%pes
            valid_data(i)%idx = reference%idx
            valid_data(i)%name = reference%name
            valid_data(i)%simbox = reference%simbox
            valid_data(i)%isimbox = reference%isimbox
        end do

        call timestamp(fit_out); write(fit_out, '(a)') "allocation complete"

        call read_fitting_data(simparams%fit_training_folder, train_data)
        call timestamp(fit_out)
        write(fit_out, '(a)') "training data points read"

        call read_fitting_data(simparams%fit_validation_folder, valid_data)
        call timestamp(fit_out)
        write(fit_out, '(a)') "validation data points read"

        ! subtract DFT reference energy
        do i = 1, training_data_points
            train_nrg(i) = train_data(i)%epot(1) - simparams%evasp
        end do
        do i = 1, validation_data_points
            valid_nrg(i) = valid_data(i)%epot(1) - simparams%evasp
        end do

        call from_pes_params_to_x(reference, x)

        if (size(x) > size(train_nrg)) stop err // "number of parameters to fit is larger than training data set"

        ! do the fit
        call timestamp(fit_out)

!        do i = 1, training_data_points
!            print *, train_data(i)%nbeads
!        end do

        if (simparams%maxit > 0) then
            write(fit_out,'(a, i, a)') "starting fit with max ", simparams%maxit, " iterations"
            call nllsqwbc(x, train_dev, train_rmse)
        else
            write(fit_out,'(a)') "comparison to training data set"
            call calc_force(reference, energy_only)
            call obj_func(training_data_points, nfit_params, x, train_dev)
            if (training_data_points > 0) then
                train_rmse = sqrt(sum(train_dev*train_dev)/training_data_points)
            else
                train_rmse = 0.0_dp
            end if
        end if

        ! compare to validation data set
        call timestamp(fit_out)
        write(fit_out, '(a)') "comparison to validation data set"
        !$omp parallel do private(i)
        do i = 1, validation_data_points
            call calc_force(valid_data(i), energy_only)
            valid_dev(i) = valid_nrg(i) - (valid_data(i)%epot(1)-reference%epot(1))
        end do
        !$omp end parallel do

        ! evaluate validation data set
        if (validation_data_points > 0) then
            valid_rmse = sqrt(sum(valid_dev*valid_dev)/validation_data_points)
        else
            valid_rmse = 0.0_dp
        end if

        write(fit_out, '(a, f13.1, a)') "training rmse:  ", 1000*train_rmse, " meV"
        write(fit_out, '(a, f13.1, a)') "validation rmse:", 1000*valid_rmse, " meV"

        close(fit_out)

        ! output results of the fit
        call write_deviations_of_fit()
        call output_pes(reference)

    end subroutine perform_fit




    subroutine read_fitting_data(folder, data)

        character(len=max_string_length), intent(in)  :: folder
        type(universe), intent(inout) :: data(:)

        character(len=*), parameter :: err = "Error in read_fitting_data(): "
        character(len=max_string_length) :: words(100), buffer
        character(len=max_string_length) :: inp_file
        integer, parameter :: inp_unit = 67
        integer :: i, ios, nwords, nfiles, nrg_pos, pos_pos, atom_pos
        logical :: config_section


        ! select folder
        if (folder == simparams%fit_training_folder) then
            nfiles = simparams%fit_training_data
        else if (folder == simparams%fit_validation_folder) then
            nfiles = simparams%fit_validation_data
        else
            print *, err // "internal error"; stop
        end if

        nrg_pos = 0; pos_pos = 0

        do i = 1, nfiles

            ! build input file string
            write(inp_file, '(a4, i5.5, a4)') "nrg_", i, ".dat"     ! ex: fit/train/nrg_00581.dat
            inp_file =  trim(folder) // "/" // inp_file

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

            ! build input file string
            write(inp_file, '(a4, i5.5, a4)') "pos_", i, ".dat"     ! ex: fit/train/nrg_00581.dat
            inp_file =  trim(folder) // "/" // inp_file

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

        character(len=max_string_length), intent(in) :: folder

        character(len=max_string_length) :: inp_file
        character(len=max_string_length) :: buffer, words(100)
        character(len=*), parameter :: err = "Error in check_fitting_data(): "
        integer, parameter :: inp_unit = 67
        integer :: this_data_points ! counts config/energy pairs in one file
        integer :: i, ios, nwords, nfiles

        ! select folder
        if (folder == simparams%fit_training_folder) then
            nfiles = simparams%fit_training_data
        else if (folder == simparams%fit_validation_folder) then
            nfiles = simparams%fit_validation_data
        else
            print *, err // "internal error in check_fitting_data"; stop
        end if

        do i = 1, nfiles

            ! build input file string
            write(inp_file, '(a4, i5.5, a4)') "nrg_", i, ".dat"     ! ex: fit/train/nrg_00581.dat
            inp_file =  trim(folder) // "/" // inp_file

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
            if (folder == simparams%fit_training_folder) then
                training_data_points = training_data_points + this_data_points
            else if (folder == simparams%fit_validation_folder) then
                validation_data_points = validation_data_points + this_data_points
            end if

            ! build input file string
            write(inp_file, '(a4, i5.5, a4)') "pos_", i, ".dat"     ! ex: fit/valid/pos_00581.dat
            inp_file =  trim(folder) // "/" // inp_file

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



    subroutine fit_setup()

        character(len=14) :: folder
        character(len=25) :: param_folder
        character(len=26) :: dev_folder
        character(len=30) :: fname

        ! XXX: change system() to execute_command_line() when new compiler is available
        ! create the fits directory
        if (.not. dir_exists("fit/"))      call system("mkdir fit/")
        if (.not. dir_exists("fit/fits/")) call system("mkdir fit/fits/")

        ! create folder for output of this fit
        write(folder, '(a9, i5.5)') "fit/fits/", simparams%start
        if (.not. dir_exists(folder)) call system("mkdir " // folder)

        ! create folder for fitted parameters
        write(param_folder, '(a14, a11)') folder, "/out_params"
        if (.not. dir_exists(param_folder)) call system("mkdir " // param_folder)

        ! create folder for deviations
        write(dev_folder, '(a14, a12)') folder, "/deviations/"
        if (.not. dir_exists(dev_folder)) call system("mkdir " // dev_folder)

        ! file that contains details to the fit
        write(fname, '(a14, a16)') folder, "/fit_details.txt"
        call open_for_write(fit_out, fname)

        call timestamp(fit_out); write(fit_out, '(a)') "fit started"

        call omp_set_num_threads(simparams%nthreads)
        call mkl_set_num_threads(1)

    end subroutine fit_setup


    subroutine write_deviations_of_fit()

        character(len=38) :: t_dev_file
        character(len=40) :: v_dev_file
        integer :: i

        write(t_dev_file, '(a9, i5.5, a24)') "fit/fits/", simparams%start, "/deviations/training.dat"
        call open_for_write(fit_out, t_dev_file)
        write(fit_out, '(a26, a18)') "# target epot/eV", "this epot/eV"
        do i = 1, training_data_points
            write(fit_out, '(i6, 2f18.5)') i, train_nrg(i), train_data(i)%epot(1)-reference%epot(1)
        end do
        close(fit_out)

        write(v_dev_file, '(a9, i5.5, a26)') "fit/fits/", simparams%start, "/deviations/validation.dat"
        call open_for_write(fit_out, v_dev_file)
        write(fit_out, '(a26, a18)') "# target epot/eV", "this epot/eV"
        do i = 1, validation_data_points
            write(fit_out, '(i6, 2f18.5)') i, valid_nrg(i), valid_data(i)%epot(1)-reference%epot(1)
        end do
        close(fit_out)

    end subroutine write_deviations_of_fit



    subroutine from_pes_params_to_x(ref, x)

        type(universe), intent(in) :: ref
        real(dp), allocatable, intent(out) :: x(:)

        real(dp), allocatable :: params_rebo(:)


        ! TODO: add other PESs
        if (any(ref%pes == pes_id_rebo)) call get_fit_params_rebo(params_rebo)

        ! add other sizes when implementing other PESs
        nfit_params = size([params_rebo])

        allocate(x(nfit_params))

        ! add other PES param array when implementing other PESs
        x = [params_rebo]

    end subroutine from_pes_params_to_x



    subroutine from_x_to_pes_params(x)

        real(dp), intent(in)       :: x(:)

        integer :: pos

        pos = 1 ! counts the position in the x array

        if (REBO_INIT_SUCCESS) call set_fit_params_rebo(x, pos)

    end subroutine from_x_to_pes_params


    subroutine nllsqwbc(x, fvec, train_rmse)
        use mkl_rci
        use mkl_rci_type
        implicit none


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

        !** results of input parameter checking
        integer :: info(6)
        real(dp) :: max_val, train_rmse_old, mean_dev
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

        m = training_data_points
        n = size(x)

        allocate(fjac(m,n))

        !** set initial values
        fvec = 0.0_dp
        fjac = 0.0_dp
        train_rmse_old = 0.0_dp

        !** write to fit details file
        write(fit_out, '(a)') ""
        write(fit_out, '(a30, 3a13, a14)') "iter", "max err/meV", "avg err/meV", "rmse/meV", "Δrmse/meV"

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

            ! calculate reference energy
            call from_x_to_pes_params(x)
            call calc_force(reference, energy_only)

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

                    call output_pes(reference)

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

                    max_val = maxval(abs(fvec)) ! max deviation
                    mean_dev = sum(fvec)/training_data_points
                    train_rmse = sqrt(sum(fvec*fvec)/training_data_points)

                    call timestamp(fit_out)
                    write(fit_out, '(i6, 4f13.1)') iter, 1000*max_val, 1000*mean_dev, &
                        1000*train_rmse, 1000*(train_rmse-train_rmse_old)

                    train_rmse_old = train_rmse
                    iter = iter + 1

!                    do i=2,M
!                        if (mod(i, M/10) == 0) then
!                            write(*,"(2i6, 2f12.4)") iter, i-1, train_nrg(i), train_nrg(i)-fvec(i)
!                        end if
!                    end do


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

        max_val = maxval(abs(fvec)) ! max deviation
        mean_dev = sum(fvec)/training_data_points
        train_rmse = sqrt(sum(fvec*fvec)/training_data_points)
        call timestamp(fit_out)
        write(fit_out, '(i6, 4f13.1)') iter, 1000*max_val, 1000*mean_dev, &
                        1000*train_rmse, 1000*(train_rmse-train_rmse_old)

        call mkl_free_buffers

        write(fit_out, "(a)") ""
        call timestamp(fit_out); write(fit_out, "(a)") "fit complete"

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

        call calc_force(reference, energy_only)

        !$omp parallel do private(i)
        do i = 1, m
            !print *, "point", i, "on thread", omp_get_thread_num(), "of", omp_get_num_threads()
            call calc_force(train_data(i), energy_only)
            f(i) = train_nrg(i) - (train_data(i)%epot(1)-reference%epot(1))
        end do
        !$omp end parallel do

    end subroutine obj_func

end module fit
