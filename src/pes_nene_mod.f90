module pes_nene_mod

    use constants
    use useful_things, only : split_string, lower_case, file_exists
    use universe_mod
    use open_file, only : open_for_read

    implicit none


    contains

    type runner_input_parameters

        ! input.nn
        integer :: nn_type_short
        integer :: mode
        logical :: lshort
        integer :: maxnum_layers_short_atomic
        integer :: nelem
        integer :: npairs
        integer :: global_hidden_layers_short
        character(len=3), dimension(atoms%ntypes) :: element

        

        ! shared
        integer :: ndim ! in
        
        ! weights.XXX.data
        ! readweights.f90
        integer :: maxnum_weights_local ! in
        integer :: num_weights_local(ndim) ! in
        real(dp)  :: weights_local, dimension(maxnum_weights_local, ndim) ! out

        ! scaling.data
        ! readscale.f90
        integer :: maxnum_funcvalues_local                        ! in
        integer :: num_funcvalues_local(ndim)                     ! in
        real(dp)  :: avvalue_local(ndim,maxnum_funcvalues_local)  ! out
        real(dp)  :: maxvalue_local(ndim,maxnum_funcvalues_local) ! out
        real(dp)  :: minvalue_local(ndim,maxnum_funcvalues_local) ! out
        real(dp)  :: eshortmin                                    ! out
        real(dp)  :: eshortmax                                    ! out


    end type

    type(runner_input_parameters) :: inpnnparam

    function new_runner_input_parameters

        type(runner_input_parameters) new_runner_input_parameters

        new_runner_input_parameters%nn_type_short               = 1 ! default, no pair available
        new_runner_input_parameters%mode                        = 3 ! default, only RuNNer mode 3 implemented
        new_runner_input_parameters%lshort                      = .true. ! default, only short
        new_runner_input_parameters%maxnum_layers_short_atomic  = 0
        new_runner_input_parameters%nelem                       = default_int
        new_runner_input_parameters%npairs                      = 0
        new_runner_input_parameters%global_hidden_layers_short  = default_int
        new_runner_input_parameters%element                     = default_string
        new_runner_input_parameters%weights_local               = default_real



    end function


    ! Here all necessary files and keywords are read in for the high-dimensional neural network potentials (HDNNPs)
    subroutine read_nene(atoms, inp_unit)

        use run_config, only : simparams

        type(universe), intent(inout) :: atoms
        integer, intent(in) :: inp_unit

        integer :: nwords, ios = 0, line = 0
        character(len=max_string_length) :: buffer
        character(len=max_string_length) :: words(100)
        character(len=max_string_length) :: inp_path

        character(len=max_string_length)   :: filename_input_nn, filename_scaling
        character(len=max_string_length), dimension(atoms%ntypes)   :: weight_names, weight_names_list

        logical, dimension(2) :: lshort

        integer  :: idx1, idx2, ntypes, weight_counter, nelem_counter_1, nelem_counter_2, element_counter
        character(len=*), parameter :: err = "Error in read_nene: "

        ntypes = simparams%nprojectiles+simparams%nlattices

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
        atoms%pes(idx1,idx2) = pes_id_nene
        atoms%pes(idx2,idx1) = pes_id_nene


        do
            read(inp_unit, '(A)', iostat=ios) buffer
            call split_string(buffer, words, nwords)

            ! pes block terminated, exit
            if (nwords == 0 .or. ios /= 0) then
                exit

            ! something went wrong
            else if (nwords /= 2 .and. words(2) /= "weights") then
                print *,  err // "Error in the PES file: PES parameters must &
                        consist of key value pairs. A parameter block must be terminated by a blank line."
                stop

            else if (words(2) == "weights" .and. nwords .neqv. atoms%ntypes+1) then ! check for valid number of weight file names given
                print *,  err // "Number of weight files given and number of elements in the structure does not match!"
                stop
            end if

            call lower_case(words(1))

            ! readout of path related to the RuNNer framework and additional the weight file names
            select case (words(1))

                case ('inp_dir')

                    read(words(2), '(A)') inp_path

                case ('weights')

                    do weight_counter = 1,atoms%ntypes

                        read(words(weight_counter+1), '(A)') weight_names(weight_counter)

                    end do

                case default

                    print *, "Error in the PES file: unknown nene parameter", words(1)
                    stop

            end select

        end do

        ! set name strings for RuNNer related files
        filename_input_nn   = trim(inp_path) // "input.nn"
        filename_scaling = trim(inp_path) // "scaling.data"

        ! loop over weight file names
        do weight_counter = 1,atoms%ntypes
            weight_names_list(weight_counter)  = trim(inp_path) // trim(weight_names(weight_counter))
        end do


        ! check existance of input.nn
        if (.not. file_exists(filename_input_nn)) stop err // "input.nn file does not exist"

        ! read all input keywords from input.nn
        call open_for_read(input_nn_unit, filename_input_nn); ios = 0

        do while (ios == 0) ! analog to read_input_file subroutine in run_config.f90
            read(input_nn_unit, '(A)', iostat=ios) buffer
            if (ios == 0) then
                line = line + 1
                call split_string(buffer, words, nwords)

                ! read in all necessary keywords
                ! only read in necessary one, only look for multiple setting for relevant keywords to avoid multiple variable settings
                select case (words(1))

                case ('nn_type_short') ! no pair type available, skip keyword
                        if (inpnnparam%nn_type_short /= default_int) stop err // 'Multiple use of the nn_type_short key'
                        if (nwords == 2) then
                            read(words(2),'(A)') inpnnparam%nn_type_short
                        else
                            print *, err, "nn_type_short key needs a single argument"; stop
                        end if

                    case ('runner_mode') ! mode 3 default, skip keyword
                        if (inpnnparam%mode /= default_int) stop err // 'Multiple use of the runner_mode key'
                        if (nwords == 2) then
                            read(words(2),'(A)') inpnnparam%mode
                            if (words(2) .neqv. 3) then
                                print *, err, "mode ", words(2), " not supported"; stop
                        else
                            print *, err, "runner_mode key needs a single argument"; stop
                        end if

                    case ('parallel_mode') ! other than default -> skip with message?
                        if (inpnnparam%paramode /= default_) stop err // 'Multiple use of the parallel_mode key'
                        if (nwords == 2) then
                            read(words(2),'(A)') inpnnparam%paramode
                            if (words(2) .neqv. 1) then
                                print *, err, "parallel_mode ", words(2), " not supported"; stop
                        else
                            print *, err, "parallel_mode key needs a single argument"; stop
                        end if

                    case ('number_of_elements') ! check with md_tian.inp
                        if (inpnnparam%nelem /= default_int) stop err // 'Multiple use of the number_of_elements key'
                        if (nwords == 2) then
                            read(words(2),'(A)') inpnnparam%nelem
                            do nelem_counter_1=1,inpnnparam%nelem
                                do nelem_counter_2=1,inpnnparam%nelem
                                    inpnnparam%npairs = inpnnparam%npairs + 1
                                end do
                            end do
                        else
                            print *, err, "number_of_elements key needs a single argument"; stop
                        end if

                    case ('elements') ! check with md_tian.inp
                        if (inpnnparam%element /= default_string) stop err // 'Multiple use of the elements key'
                        if (nwords == atoms%ntypes+1) then
                            do element_counter = 1,atoms%ntypes
                                read(words(element_counter+1),'(A)') inpnnparam%element(element_counter)
                            end do
                        else
                            print *, err, "elements key does not match with number of element types"; stop
                        end if

                    case ('random_seed')
                        if (inpnnparam% /= default_) stop err // 'Multiple use of the random_seed key'

                    case ('random_number_type')
                        if (inpnnparam% /= default_) stop err // 'Multiple use of the random_number_type key'

                    case ('remove_atom_energies')
                        if (inpnnparam% /= default_) stop err // 'Multiple use of the remove_atom_energies key'

                    case ('atom_energy')
                        if (inpnnparam% /= default_) stop err // 'Multiple use of the atom_energy key'

                    case ('energy_threshold')
                        if (inpnnparam% /= default_) stop err // 'Multiple use of the energy_threshold key'

                    case ('bond_threshold')
                        if (inpnnparam% /= default_) stop err // 'Multiple use of the bond_threshold key'

                    case ('use_short_nn')
                        if (inpnnparam% /= default_bool) stop err // 'Multiple use of the use_short_nn key'
                        if (nwords == 1) then
                            inpnnparam% = .true.
                        else
                            print *, err, "use_short_nn key has additional arguments"; stop
                        end if


                    case ('global_hidden_layers_short')
                        if (inpnnparam%global_hidden_layers_short /= default_int) stop err // 'Multiple use of the global_hidden_layers_short key'
                        if (nwords == 2) then
                            inpnnparam%global_hidden_layers_short = inpnnparam%global_hidden_layers_short + 1
                        else
                            print *, err, "global_hidden_layers_short key needs a single argument"; stop
                        end if

                    case ('global_nodes_short')
                        if (inpnnparam% /= default_) stop err // 'Multiple use of the global_nodes_short key'

                    case ('global_activation_short')
                        if (inpnnparam% /= default_) stop err // 'Multiple use of the global_activation_short key'

                    case ('element_hidden_layers_short')
                        if (inpnnparam% /= default_) stop err // 'Multiple use of the element_hidden_layers_short key'

                    case ('element_nodes_short')
                        if (inpnnparam% /= default_) stop err // 'Multiple use of the element_nodes_short key'

                    case ('element_activation_short')
                        if (inpnnparam% /= default_) stop err // 'Multiple use of the element_activation_short key'

                    case ('cutoff_type')
                        if (inpnnparam% /= default_) stop err // 'Multiple use of the cutoff_type key'

                    case ('symfunction_short')
                        if (inpnnparam% /= default_) stop err // 'Multiple use of the symfunction_short key'

                    case ('points_in_memory')
                        if (inpnnparam% /= default_) stop err // 'Multiple use of the points_in_memory key'

                    case ('scale_symmetry_functions')
                        if (inpnnparam% /= default_) stop err // 'Multiple use of the scale_symmetry_functions key'

                    case ('center_symmetry_functions')
                        if (inpnnparam% /= default_) stop err // 'Multiple use of the center_symmetry_functions key'

                    case ('fitting_unit')
                        if (inpnnparam% /= default_) stop err // 'Multiple use of the fitting_unit key'

                    case ('use_old_scaling')
                        if (inpnnparam% /= default_) stop err // 'Multiple use of the use_old_scaling key'

                    case ('scale_min_short')
                        if (inpnnparam% /= default_) stop err // 'Multiple use of the scale_min_short key'

                    case ('scale_max_short')
                        if (inpnnparam% /= default_) stop err // 'Multiple use of the scale_max_short key'

                    case ('calculate_forces')
                        if (inpnnparam% /= default_) stop err // 'Multiple use of the calculate_forces key'

                    case ('check_forces')
                        if (inpnnparam% /= default_) stop err // 'Multiple use of the check_forces key'

                    case ('calculate_stress')
                        if (inpnnparam% /= default_) stop err // 'Multiple use of the calculate_stress key'

                    case ('use_atom_energies')
                        if (inpnnparam% /= default_) stop err // 'Multiple use of the  key'

                    case ('')
                        if (inpnnparam% /= default_) stop err // 'Multiple use of the  key'

                    case ('')
                        if (inpnnparam% /= default_) stop err // 'Multiple use of the  key'

                    case ('')
                        if (inpnnparam% /= default_) stop err // 'Multiple use of the  key'

                    case ('')
                        if (inpnnparam% /= default_) stop err // 'Multiple use of the  key'

                    case ('')
                        if (inpnnparam% /= default_) stop err // 'Multiple use of the  key'




                    case default ! for all unnecessary keywords make a dummy case?
                        if (trim(words(1)) /= '' .and. words(1)(1:1) /= '#') & !check for empty and comment lines
                            print *, 'Skipping invalid label ',trim(words(1)),' in line', line

                end select

            else

                write(*,*) err // 'Error reading out from input.nn: iostat = ', ios
                stop

            end if
        close(input_nn_unit)


        ! check existance of scaling.data
        if (.not. file_exists(filename_scaling)) stop err // 'scaling.data file does not exist'

        ! read in all data from scaling.data
        call open_for_read(scaling_unit, filename_scaling); ios = 0

        do while (ios == 0) 
            read(scaling_unit, '(A)', iostat=ios) buffer
            if (ios == 0) then

                ! readscale.f90
                ! do i1=1,ndim
                !     do i2=1,num_funcvalues_local(i1)
                !         read(scaleunit,*)i3,i3,minvalue_local(i1,i2),maxvalue_local(i1,i2),avvalue_local(i1,i2)
                !     enddo ! i2
                ! enddo ! i1

                ! read(scaleunit,*)eshortmin,eshortmax
           
            else

                    write(*,*) err // 'Error reading out from scaling.data: iostat = ', ios
                    stop

            end if
        end do ! ios
        close(scaling_unit)


        ! loop over weight.XXX.data files and read in all data
        do weight_counter = 1,atoms%ntypes ! weights

            ! check existance of each weight file before reading
            if (.not. file_exists(weight_names_list(weight_counter))) stop err // weight_names_list(weight_counter), ' file does not exist'

            call open_for_read(weight_unit, weight_names_list(weight_counter)); ios = 0

            do while (ios == 0) ! ios
                read(scaling_unit, '(A)', iostat=ios) buffer
                if (ios == 0) then

                        ! readweights.f90
                        ! do i1=1,ndim
                        !     do i2=1,num_weights_local(i1)
                        !         read(wunit,*)weights_local(i2,i1)
                        !     end do
                        ! end do

                else

                        write(*,*) err // 'Error reading out from',weight_names_list(weight_counter),': iostat = ', ios
                        stop

                end if

            end do ! ios
            close(weight_unit)

        end do ! weights

!       ! from RuNNer main.f90
!       ! start mpi routines
!        call mpi_init(mpierror)
!            if(mpierror.ne.0)then
!                print *, 'Error in mpi_init ',mpierror
!                stop
!            endif
!       ! get number of processes mpisize
!        call mpi_comm_size(mpi_comm_world,mpisize,mpierror)
        ! get process id mpirank
!        call mpi_comm_rank(mpi_comm_world,mpirank,mpierror)


    end subroutine read_nene


    subroutine compute_nene(atoms, flag)

        ! Calculates energy and forces with HDNNPs

        type(universe), intent(inout)   :: atoms
        integer, intent(in)             :: flag



    end subroutine compute_nene

!    subroutine cleanup_nene

        ! from RuNNer cleanup.f90
        ! deallocate everything
!        if(lshort.and.(nn_type_short.eq.1))then
!            deallocate(weights_short_atomic)
!            deallocate(symfunction_short_atomic_list)
!            deallocate(num_funcvalues_short_atomic)
!            deallocate(windex_short_atomic)
!            deallocate(num_layers_short_atomic)
!            deallocate(actfunc_short_atomic)
!            deallocate(nodes_short_atomic)
!            deallocate(num_weights_short_atomic)
!            deallocate(function_type_short_atomic)
!            deallocate(symelement_short_atomic)
!            deallocate(funccutoff_short_atomic)
!            deallocate(eta_short_atomic)
!            deallocate(zeta_short_atomic)
!            deallocate(lambda_short_atomic)
!            deallocate(rshift_short_atomic)
!        endif

!        call mpi_barrier(mpi_comm_world,mpierror)

        ! from RuNNer main.f90
        ! shutdown mpi routines
!        call mpi_finalize(mpierror)
!        if(mpierror.ne.0)then
!            print *, 'mpierror finalize ', mpierror
!        endif

!    end subroutine cleanup_nene

end module pes_nene_mod
