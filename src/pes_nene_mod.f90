!############################################################################
! This routine is part of
! md_tian2 (Molecular Dynamics Tian Xia 2)
! (c) 2014-2019 Dan J. Auerbach, Svenja M. Janke, Marvin Kammler,
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

module pes_nene_mod

    !use constants, only : max_string_length, pes_id_nene, default_string, default_int, default_real, default_bool, inpnn_unit, scaling_unit, weight_unit,
    !use constants, disabled1 => pi, disabled2 => rad2deg ! with as the opposite of "use module, only :"
    use universe_mod
    use get_defaults

    implicit none

    contains

    !type runner_input_parameters

        !private
        !character(len=max_string_length)                :: filename_inpnn, filename_scaling
        !character(len=max_string_length), allocatable   :: filename_weights(:)

    !end type

    !type(runner_input_parameters) :: rinpparam



!   2do:
!   use RuNNer modules if necessary, otherwise make own ones
!   change how the seeed for the random number generator will be (add new variable?)
!   variable declarations concerning RuNNer in the corresponding modules, but set to (our) default values has to be done before reading out keywords (own subroutine called in compute_nene)
!   move RuNNer related files to folder and change the makefile
!   check how many mpi routines has to stay in the code, at least set the few default values so that no error will occur due to wrong default mpi settings, therefore the mpi_dummy_routines.f90 file makes sense
!   don't explicitly give weight file names, use RuNNer routine instead


    ! Here all necessary files and keywords are read in for the high-dimensional neural network potentials (HDNNPs)
    subroutine read_nene(atoms, inp_unit)

        use constants, only : max_string_length, pes_id_nene, default_string, default_int, default_real, default_bool
        use open_file, only : lower_case, open_for_read, split_string
        !use run_config, only : simparams
        use useful_things, only : file_exists
        use get_defaults


        type(universe), intent(inout) :: atoms
        integer, intent(in) :: inp_unit

        integer :: nwords, ios = 0
        character(len=max_string_length) :: buffer
        character(len=max_string_length) :: words(100)
        character(len=max_string_length) :: inp_path

        character(len=max_string_length)                :: filename_inpnn, filename_scaling
        character(len=max_string_length), allocatable   :: weights_path(:), filename_weights(:)

        integer, parameter  :: inpnn_unit       = 61
        integer, parameter  :: scaling_unit     = 62
        integer, parameter  :: weight_unit      = 63

        integer  :: idx1, idx2, weight_counter
        integer  :: npairs_counter_1, npairs_counter_2, element_counter, nodes_counter
        integer  :: general_counter_1, general_counter_2, general_counter_3

        character(len=*), parameter :: err = "Error in read_nene: "
        character(len=*), parameter :: err_pes = "Error in the PES file: "

        character(len=*), parameter :: err_inpnn = "Error when reading input.nn: "
        character(len=*), parameter :: err_scaling = "Error when reading scaling.data: "
        character(len=*), parameter :: err_weight = "Error when reading the following weight file: "
        character(len=*), parameter :: warn_inpnn = "Warning when reading input.nn: "

        ! first read the pes file:
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
        atoms%pes(idx1,idx2) = pes_id_nene
        atoms%pes(idx2,idx1) = pes_id_nene

        ! initialize readout variables
        if (.not. allocated(weights_path)) then

            allocate(weights_path(atoms%ntypes))
            allocate(filename_weights(atoms%ntypes))

            weights_path        = default_string
            filename_weights    = default_string
            inp_path            = default_string

        end if


        do
            read(inp_unit, '(A)', iostat=ios) buffer
            call split_string(buffer, words, nwords)

            ! pes block terminated, exit
            if (nwords == 0 .or. ios /= 0) then
                exit

            ! something went wrong
            else if (nwords /= 2 .and. words(1) /= "weights") then
                print *,  err // err_pes // "PES parameters must &
                        consist of key value pairs. A parameter block must be terminated by a blank line."
                stop

            ! maybe add comparison of element type and number in weight filename?
            else if (words(1) == "weights" .and. nwords /= atoms%ntypes+1) then ! check for valid number of weight file names given
                print *,  err // "Number of weight files given and number of elements in the structure file do not match!"
                stop
            end if

            call lower_case(words(1))

            ! readout of path related to the RuNNer framework and additional the weight file names
            select case (words(1))

                case ('inp_dir')

                    if (inp_path /= default_string) stop err // err_pes // 'Multiple use of the inp_dir key'
                    read(words(2), '(A)') inp_path

                case ('weights')

                    if (weights_path /= default_string) stop err // err_pes // 'Multiple use of the weights key'
                    do weight_counter = 1,atoms%ntypes ! check if ntypes = number of elements in proj AND latt (elements per structure NOT latt/proj types)!!

                        read(words(weight_counter+1), '(A)') weights_path(weight_counter)

                    end do

                case default

                    print *, err // err_pes // "unknown nene parameter ", words(1)
                    stop

            end select

        end do

        ! set name strings for RuNNer related files
        filename_inpnn   = trim(inp_path) // "input.nn"
        filename_scaling = trim(inp_path) // "scaling.data"

        ! loop over weight file names
        do weight_counter = 1,atoms%ntypes
            filename_weights(weight_counter)  = trim(inp_path) // trim(weights_path(weight_counter))
        end do

        ! in case of the HDNNPs several additional input files have to be read

        ! read all input keywords from input.nn several times to respect dependencies

!       read in keywords related to input.nn according to the following files from RuNNer (chronologically)
!       1) getdimensions.f90
!       2) paircount.f90
!       3) readkeywords.f90
!       4) readinput.f90

        ! start readout according to main.f90
        call mpi_init(mpierror)
        call mpi_comm_size(mpi_comm_world,mpisize,mpierror)
        call mpi_comm_rank(mpi_comm_world,mpirank,mpierror)

        ! set all variables to default values -> rethink this subroutine!!
        call set_defaults()

        !call initnn(iseed)

        ! start readout of input.nn according to initnn.f90
        call get_nnconstants()
        !call writeheader()

        !call initialization(ielem,lelement) -> here only getdimensions, paircount and checkstructures are needed

        ! start readout according to initialization.f90

        ! start readout according to getdimensions.f90

        ! check existance of input.nn
        if (.not. file_exists(filename_inpnn)) stop err // err_inpnn // "file does not exist"

        call open_for_read(inpnn_unit, filename_inpnn); ios = 0

        do while (ios == 0)
            read(inpnn_unit, '(A)', iostat=ios) buffer
            if (ios == 0) then
                !line = line + 1
                call split_string(buffer, words, nwords)

                select case (words(1))

                    case ('nn_type_short')
                        if (rinpparam%nn_type_short /= default_int) stop err // err_inpnn // 'Multiple use of the nn_type_short key'
                        if (nwords == 2) then
                            read(words(2),'(i1000)', iostat=ios) rinpparam%nn_type_short
                            if (ios /= 0) stop err // err_inpnn // "nn_type_short value must be integer"

                            select case (words(2))

                                case (1)
                                    ! Gandalf says: you shall pass

                                case (2)
                                    print *, err, err_inpnn, "nn_type_short 2 not supported, Pair NN not implemented!"
                                    stop

                                case default
                                    print *, err, err_inpnn, "Error in nn_type_short key value, ", words(2), " not implemented"
                                    stop

                            end select

                        else
                            print *, err, err_inpnn, "nn_type_short key needs a single argument"; stop
                        end if

                    case ('runner_mode')
                        if (rinpparam%mode /= default_int) stop err // err_inpnn // 'Multiple use of the runner_mode key'
                        if (nwords == 2) then
                            read(words(2),'(i1000)', iostat=ios) rinpparam%mode
                            if (ios /= 0) stop err // err_inpnn // "runner_mode value must be integer"
                            if (rinpparam%mode /= 3) then
                                print *, err, err_inpnn, "Only mode 3 (prediction mode) available"; stop
                        else
                            print *, err, err_inpnn, "runner_mode key needs a single argument"; stop
                        end if

                    case ('use_short_nn')
                        if (rinpparam%lshort /= default_bool) stop err // err_inpnn // 'Multiple use of the use_short_nn key'
                        if (nwords == 1) then
                            rinpparam%lshort = .true.
                        else
                            print *, err, err_inpnn, "use_short_nn key needs no argument(s)"; stop
                        end if

                    case ('use_electrostatics')
                        if (rinpparam%lelec /= default_bool) stop err // err_inpnn // 'Multiple use of the use_electrostatics key'
                        if (nwords == 1) then
                            rinpparam%lelec = .true.
                        else
                            print *, err, err_inpnn, "use_electrostatics key needs no argument(s)"; stop
                        end if

                    case ('electrostatic_type', 'nn_type_elec')
                        if (rinpparam%nn_type_elec /= default_int) stop err // err_inpnn // 'Multiple use of the electrostatic_type/nn_type_elec key'
                        if (nwords == 2) then
                            read(words(2),'(i1000)', iostat=ios) rinpparam%nn_type_elec
                            if (ios /= 0) stop err // err_inpnn // "electrostatic_type/nn_type_elec value must be integer"
                        else
                            print *, err, err_inpnn, "electrostatic_type/nn_type_elec key needs a single argument"; stop
                        end if

                    case default
                        ! for every other keyword pass here, check for unrecognized keywords later

                end select

            else
                write(*,*) err // err_inpnn // 'iostat = ', ios
                stop
            end if
        end do

        close(inpnn_unit)


        call open_for_read(inpnn_unit, filename_inpnn); ios = 0

        do while (ios == 0)
            read(inpnn_unit, '(A)', iostat=ios) buffer
            if (ios == 0) then
                line = line + 1
                call split_string(buffer, words, nwords)

                select case (words(1))

                    case ('debug_mode')
                        if (rinpparam%lfounddebug /= default_bool) stop err // err_inpnn // 'Multiple use of the debug_mode key'
                        if (nwords == 1) then
                            rinpparam%lfounddebug = .true.
                            rinpparam%ldebug = .true.
                        else
                            print *, err, err_inpnn, "debug_mode key needs no argument(s)"; stop
                        end if

                    case ('global_hidden_layers_short')
                        if (rinpparam%lfound_num_layersshort /= default_bool) stop err // err_inpnn // 'Multiple use of the global_hidden_layers_short key'
                        if (nwords == 2) then
                            rinpparam%lfound_num_layersshort = .true.
                            read(words(2),'(i1000)', iostat=ios) rinpparam%maxnum_layers_short_atomic
                            if (ios /= 0) stop err // err_inpnn // "global_hidden_layers_short value must be integer"
                            rinpparam%maxnum_layers_short_atomic = rinpparam%maxnum_layers_short_atomic + 1
                        else
                            print *, err, err_inpnn, "global_hidden_layers_short key needs a single argument"; stop
                        end if

                    case ('global_hidden_layers_electrostatic')
                        if (rinpparam%lfound_num_layersewald /= default_bool) stop err // err_inpnn // 'Multiple use of the global_hidden_layers_electrostatic key'
                        if (nwords == 2) then
                            rinpparam%lfound_num_layersewald = .true.
                            read(words(2),'(i1000)', iostat=ios) rinpparam%maxnum_layers_elec
                            if (ios /= 0) stop err // err_inpnn // "global_hidden_layers_electrostatic value must be integer"
                            rinpparam%maxnum_layers_elec = rinpparam%maxnum_layers_elec + 1
                        else
                            print *, err, err_inpnn, "global_hidden_layers_electrostatic key needs a single argument"; stop
                        end if

                    case ('global_hidden_layers_pair')
                        print *, err, err_inpnn, "global_hidden_layers_pair key not supported, Pair NN not implemented"; stop

                    case ('use_atom_energies')
                        if (rinpparam%lfound_luseatomenergies /= default_bool) stop err // err_inpnn // 'Multiple use of the use_atom_energies key'
                        if (nwords == 1) then
                            rinpparam%lfound_luseatomenergies = .true.
                            rinpparam%luseatomenergies = .true.
                        else
                            print *, err, err_inpnn, "use_atom_energies key needs no argument(s)"; stop
                        end if

                    case ('use_atom_charges')
                        if (rinpparam%lfound_luseatomcharges /= default_bool) stop err // err_inpnn // 'Multiple use of the use_atom_charges key'
                        if (nwords == 1) then
                            rinpparam%lfound_luseatomcharges = .true.
                            rinpparam%luseatomcharges = .true.
                        else
                            print *, err, err_inpnn, "use_atom_charges key needs no argument(s)"; stop
                        end if

                    case ('number_of_elements')
                        if (rinpparam%lfound_nelem /= default_bool) stop err // err_inpnn // 'Multiple use of the number_of_elements key'
                        if (nwords == 2) then
                            rinpparam%lfound_nelem = .true.
                            read(words(2),'(i1000)', iostat=ios) rinpparam%nelem
                            if (ios /= 0) stop err // err_inpnn // "number_of_elements value must be integer"
                            if (rinpparam%nelem /= atoms%ntypes) stop err // err_inpnn // "number of elements in input.nn and in structure file differ"
                            rinpparam%npairs = 0
                            do npairs_counter_1 = 1,rinpparam%nelem
                                do npairs_counter_2 = npairs_counter_1,rinpparam%nelem
                                    rinpparam%npairs = rinpparam%npairs + 1
                                end do
                            end do
                        else
                            print *, err, err_inpnn, "number_of_elements key needs a single argument"; stop
                        end if

                    case default
                        ! for every other keyword pass here, check for unrecognized keywords later

                end select

            else
                write(*,*) err // err_inpnn // 'iostat = ', ios
                stop
            end if
        end do

        close(inpnn_unit)

        if (rinpparam%lshort .and. (rinpparam%nn_type_short == 1) .and. (rinpparam%maxnum_layers_short_atomic == default_int)) stop err // err_inpnn // 'global_hidden_layers_short key not set'
        !if (rinpparam%lshort .and. (rinpparam%nn_type_short == 2) .and. (rinpparam%maxnum_layers_short_pair == default_int)) stop err // err_inpnn // 'global_hidden_layers_pairs key not set'
        if (rinpparam%lelec .and. (rinpparam%nn_type_elec == 1) .and. (rinpparam%maxnum_layers_elec == default_int)) stop err // err_inpnn // 'global_hidden_layers_electrostatic key not set'

        allocate(rinpparam%nucelem(rinpparam%nelem))
        allocate(rinpparam%element(rinpparam%nelem))
        allocate(rinpparam%dmin_element(rinpparam%nelem*(rinpparam%nelem+1)/2))


        call open_for_read(inpnn_unit, filename_inpnn); ios = 0 ! maybe move to readout before since nelem is given by atoms%ntypes

        do while (ios == 0)
            read(inpnn_unit, '(A)', iostat=ios) buffer
            if (ios == 0) then
                !line = line + 1
                call split_string(buffer, words, nwords)

                select case (words(1))

                    case ('elements')
                        if (any(rinpparam%element /= default_string)) stop err // err_inpnn // 'Multiple use of the elements key'
                        !if (nwords == atoms%nnelem+1) then
                        if (nwords == atoms%ntypes+1) then ! check if ntypes is equal to number of elements!!
                        !if (nwords == rinpparam%nelem+1) then
                            !do element_counter = 1,atoms%nelem
                            do element_counter = 1,atoms%ntypes
                            !do element_counter = 1,rinpparam%nelem
                                read(words(element_counter+1),'(A)') rinpparam%element(element_counter) !! check for valid name will be done later!
                            end do
                            if (any(rinpparam%element /= atoms%name)) stop err // err_inpnn // "element names in input.nn and in structure file differ"
                        else
                            print *, err, err_inpnn, "elements key does not match with number of element types"; stop
                        end if

                    case default
                        ! for every other keyword pass here, check for unrecognized keywords later

                end select

            else
                write(*,*) err // err_inpnn // 'iostat = ', ios
                stop
            end if
        end do

        close(inpnn_unit)

        do nuc_counter=1,atoms%ntypes
            call nuccharge(rinpparam%element(i),rinpparam%nucelem(i))
        end do

        call sortelements()

        if (rinpparam%maxnum_layers_short_atomic == default_int) then
            rinpparam%maxnum_layers_short_atomic = 0
        end if
        if (rinpparam%maxnum_layers_elec == default_int) then
            rinpparam%maxnum_layers_elec = 0
        end if
        !if (rinpparam%maxnum_layers_short_pair == default_int) then ! not needed
        !    rinpparam%maxnum_layers_short_pair = 0
        !end if
        if (rinpparam%lfound_nelem == default_bool) stop err // err_inpnn // "number_of_elements key not found"

        allocate(rinpparam%num_funcvalues_local(102))
        allocate(rinpparam%num_funcvaluese_local(102))
        !allocate(rinpparam%num_funcvaluesp_local(102,102)) ! not needed
        rinpparam%num_funcvalues_local(:) = 0
        rinpparam%num_funcvaluese_local(:) = 0
        !rinpparam%num_funcvaluesp_local(:) = 0 ! not needed

        if (rinpparam%maxnum_layers_short_atomic .gt. 0) then
            allocate(rinpparam%nodes_short_local(0:rinpparam%maxnum_layers_short_atomic))
            rinpparam%nodes_short_local(:) = default_int ! = 0 in getdimensions.f90
        end if
        if (rinpparam%maxnum_layers_selec .gt. 0) then
            allocate(rinpparam%nodes_ewald_local(0:rinpparam%maxnum_layers_elec))
            rinpparam%nodes_ewald_local(:) = default_int ! = 0 in getdimensions.f90
        end if
        if (rinpparam%maxnum_layers_short_pair .gt. 0) then ! not needed
            allocate(rinpparam%nodes_pair_local(0:rinpparam%maxnum_layers_short_pair))
            rinpparam%nodes_pair_local(:) = default_int ! = 0 in getdimensions.f90
        end if


        call open_for_read(inpnn_unit, filename_inpnn); ios = 0

        do while (ios == 0)
            read(inpnn_unit, '(A)', iostat=ios) buffer
            if (ios == 0) then
                !line = line + 1
                call split_string(buffer, words, nwords)

                select case (words(1))

                    case ('global_nodes_short')
                        if (rinpparam%nodes_short_local /= default_int) stop err // err_inpnn // 'Multiple use of the global_nodes_short key'
                        if (nwords == rinpparam%maxnum_layers_short_atomic+1) then
                            do nodes_counter = 1,rinpparam%maxnum_layers_short_atomic-1
                                read(words(nodes_counter+1),'(i1000)', iostat=ios) rinpparam%nodes_short_local(nodes_counter)
                                if (ios /= 0) stop err // err_inpnn // "global_nodes_short value", nodes_counter, " must be integer"
                            end do
                        else
                            print *, err, err_inpnn, "global_nodes_short argument number does not match with global_hidden_layers_short value"; stop
                        end if

                    case ('global_nodes_electrostatic')
                        if (rinpparam%nodes_ewald_local /= default_int) stop err // err_inpnn // 'Multiple use of the global_nodes_electrostatic key'
                        if (nwords == rinpparam%maxnum_layers_elec+1) then
                            do nodes_counter = 1,rinpparam%maxnum_layers_elec-1
                                read(words(nodes_counter+1),'(i1000)', iostat=ios) rinpparam%nodes_ewald_local(nodes_counter)
                                if (ios /= 0) stop err // err_inpnn // "global_nodes_electrostatic value", nodes_counter, " must be integer"
                            end do
                        else
                            print *, err, err_inpnn, "global_nodes_electrostatic argument number ", nwords, " does not match with global_hidden_layers_electrostatic value ", maxnum_layers_elec-1; stop
                        end if

!                   case ('global_nodes_pair') ! not needed - make dummy
!                       if (rinpparam%nodes_pair_local /= default_int) stop err // err_inpnn // 'Multiple use of the global_nodes_pair key'
!                       if (nwords == rinpparam%maxnum_layers_short_pair+1) then
!                           do nodes_counter = 1,rinpparam%maxnum_layers_short_pair-1
!                               read(words(nodes_counter+1),'(i1000)', iostat=ios) rinpparam%nodes_pair_local(nodes_counter)
!                               if (ios /= 0) stop err // err_inpnn // "global_nodes_pair value", nodes_counter, " must be integer"
!                           end do
!                       else
!                           print *, err, err_inpnn, "global_nodes_pair argument number does not match with global_hidden_layers_pair value"; stop
!                       end if

                    case ('global_nodes_pair')
                        print *, err, err_inpnn, "global_nodes_pair key not supported, Pair NN not implemented"; stop

                    case ('element_symfunction_short')
                        !if (rinpparam% /= default_int) stop err // err_inpnn // 'Multiple use of the element_symfunction_short key'
                        !if (nwords == 2) then
                        read(words(2),'(A)', iostat=ios) rinpparam%elementtemp
                        read(words(3),'(i1000)', iostat=ios) rinpparam%function_type_local
                        if (ios /= 0) stop err // err_inpnn // "element_symfunction_short second argument value must be integer"
                        call nuccharge(rinpparam%elementtemp, rinpparam%ztemp)
                        if (rinpparam%num_funcvalues_local(rinpparam%ztemp) /= default_int) then
                            print *, err // err_inpnn // 'Error in element_symfunction_short: Element with atomic number ', rinpparam%num_funcvalues_local(rinpparam%ztemp), 'already set, check for multiple use of key'
                            stop
                        end if
                        rinpparam%num_funcvalues_local(rinpparam%ztemp) = 0
                        call lower_case(words(3))

                        select case (words(3))

                            case (1,2,4)
                                rinpparam%num_funcvalues_local(rinpparam%ztemp) = rinpparam%num_funcvalues_local(rinpparam%ztemp) + rinpparam%nelem

                            case (3,8,9)
                                rinpparam%num_funcvalues_local(rinpparam%ztemp) = rinpparam%num_funcvalues_local(rinpparam%ztemp) + rinpparam%nelem
                                if (rinpparam%nelem .gt. 1) then
                                    do general_counter_1 = 1,rinpparam$nelem-1
                                        rinpparam%num_funcvalues_local(rinpparam%ztemp) = rinpparam%num_funcvalues_local(rinpparam%ztemp) + general_counter_1
                                    end do
                                end if

                            case (5,6) ! only for Pair NN!
                                rinpparam%num_funcvalues_local(rinpparam%ztemp) = rinpparam%num_funcvalues_local(rinpparam%ztemp) + 1

                            case default
                                print *, err, err_inpnn, "Error in element_symfunction_short key, symfunction type ", words(3), " not implemented"
                                stop

                        end select
                        !else
                            !print *, err, err_inpnn, "element_symfunction_short key needs a single argument"; stop
                        !end if

                    case ('element_symfunction_electrostatic')
                        read(words(2),'(A)', iostat=ios) rinpparam%elementtemp
                        read(words(3),'(i1000)', iostat=ios) rinpparam%function_type_local
                        if (ios /= 0) stop err // err_inpnn // "element_symfunction_electrostatic second argument value must be integer"
                        call nuccharge(rinpparam%elementtemp, rinpparam%ztemp)
                        call lower_case(words(3))

                        select case (words(3))

                            case (1,2,4)
                                rinpparam%num_funcvaluese_local(rinpparam%ztemp) = rinpparam%num_funcvaluese_local(rinpparam%ztemp) + rinpparam%nelem

                            case (3,8,9)
                                rinpparam%num_funcvaluese_local(rinpparam%ztemp) = rinpparam%num_funcvaluese_local(rinpparam%ztemp) + rinpparam%nelem
                                if (rinpparam%nelem .gt. 1) then
                                    do general_counter_1 = 1,rinpparam%nelem-1
                                        rinpparam%num_funcvaluese_local(rinpparam%ztemp) = rinpparam%num_funcvaluese_local(rinpparam%ztemp) + general_counter_1
                                    end do
                                end if

                            case (5,6) ! only for Pair NN!
                                rinpparam%num_funcvalues_local(rinpparam%ztemp) = rinpparam%num_funcvalues_local(rinpparam%ztemp) + 1

                            case default
                                print *, err, err_inpnn, "Error in element_symfunction_electrostatic key, symfunction type ", words(3), " not implemented"
                                stop

                        end select

                    case ('global_symfunction_short')
                        read(words(2),'(i1000)', iostat=ios) rinpparam%function_type_local !!check if it will read only the function type and not the symbol!!
                        if (ios /= 0) stop err // err_inpnn // "global_symfunction_short (second?) argument value must be integer"
                        !call nuccharge(rinpparam%elementtemp, rinpparam%ztemp)
                        call lower_case(words(2))

                        select case (words(2))

                            case (1,2,4)
                                do general_counter_1 = 1,rinpparam%nelem
                                    rinpparam%num_funcvalues_local(rinpparam%nucelem(general_counter_1)) = rinpparam%num_funcvalues_local(rinpparam%nucelem(general_counter_1)) + rinpparam%nelem
                                end do

                            case (3,8,9)
                                do general_counter_1 = 1,rinpparam%nelem
                                    rinpparam%num_funcvalues_local(rinpparam%nucelem(general_counter_1)) = rinpparam%num_funcvalues_local(rinpparam%nucelem(general_counter_1)) + rinpparam%nelem
                                end do
                                do general_counter_1 = 1,rinpparam%nelem
                                    if (rinpparam%nelem .gt. 1) then
                                        do general_counter_2 = 1,rinpparam%nelem-1
                                            rinpparam%num_funcvalues_local(rinpparam%nucelem(general_counter_1)) = rinpparam%num_funcvalues_local(rinpparam%nucelem(general_counter_1)) + general_counter_2
                                        end do
                                    end if
                                end do

                            case (5,6) ! only for Pair NN!
                                do general_counter_1 = 1,rinpparam%nelem
                                    rinpparam%num_funcvalues_local(rinpparam%nucelem(general_counter_1)) = rinpparam%num_funcvalues_local(rinpparam%nucelem(general_counter_1)) + 1
                                end do

                            case default
                                print *, err, err_inpnn, "Error in global_symfunction_short key, symfunction type ", words(2), " not implemented"
                                stop

                        end select

                    case ('global_symfunction_electrostatic')
                        read(words(2),'(i1000)', iostat=ios) rinpparam%function_type_local !!check if it will read only the function type and not the symbol!!
                        if (ios /= 0) stop err // err_inpnn // "global_symfunction_electrostatic (second?) argument value must be integer"
                        !call nuccharge(rinpparam%elementtemp, rinpparam%ztemp)
                        call lower_case(words(2))

                        select case (words(2))

                            case (1,2,4)
                                do general_counter_1 = 1,rinpparam%nelem
                                    rinpparam%num_funcvaluese_local(rinpparam%nucelem(general_counter_1)) = rinpparam%num_funcvaluese_local(rinpparam%nucelem(general_counter_1)) + rinpparam%nelem
                                end do

                            case (3,8,9)
                                do general_counter_1 = 1,rinpparam%nelem
                                    rinpparam%num_funcvaluese_local(rinpparam%nucelem(general_counter_1)) = rinpparam%num_funcvaluese_local(rinpparam%nucelem(general_counter_1)) + rinpparam%nelem
                                end do
                                do general_counter_1 = 1,rinpparam%nelem
                                    if (rinpparam%nelem .gt. 1) then
                                        do general_counter_2 = 1,rinpparam%nelem-1
                                            rinpparam%num_funcvaluese_local(rinpparam%nucelem(general_counter_1)) = rinpparam%num_funcvaluese_local(rinpparam%nucelem(general_counter_1)) + general_counter_2
                                        end do
                                    end if
                                end do

                            case (5,6) ! only for Pair NN
                                do general_counter_1 = 1,rinpparam%nelem
                                    rinpparam%num_funcvaluese_local(rinpparam%nucelem(general_counter_1)) = rinpparam%num_funcvaluese_local(rinpparam%nucelem(general_counter_1)) + 1
                                end do

                            case default
                                print *, err, err_inpnn, "Error in global_symfunction_electrostatic key, symfunction type ", words(2), " not implemented"
                                stop

                        end select

                    case ('symfunction_short')
                        read(words(2),'(A)', iostat=ios) rinpparam%elementtemp
                        call nuccharge(rinpparam%elementtemp, rinpparam%ztemp)
                        rinpparam%num_funcvalues_local(rinpparam%ztemp) = rinpparam%num_funcvalues_local(rinpparam%ztemp) + 1

                    case ('symfunction_electrostatic')
                        read(words(2),'(A)', iostat=ios) rinpparam%elementtemp
                        call nuccharge(rinpparam%elementtemp, rinpparam%ztemp)
                        rinpparam%num_funcvaluese_local(rinpparam%ztemp) = rinpparam%num_funcvaluese_local(rinpparam%ztemp) + 1

                    case ('pairsymfunction_short') ! not needed
                        print *, err, err_inpnn, "pairsymfunction_short key is not supported, Pair NN not implemented"

                    case ('element_pairsymfunction_short') ! not needed
                        print *, err, err_inpnn, "element_pairsymfunction_short key is not supported, Pair NN not implemented"

                    case ('global_pairsymfunction_short') ! not needed
                        print *, err, err_inpnn, "global_pairsymfunction_short key is not supported, Pair NN not implemented"

                    case ('global_symfunction_short_pair') ! not needed
                        print *, err, err_inpnn, "global_symfunction_short_pair key is not supported, Pair NN not implemented"

                    case default
                        ! for every other keyword pass here, check for unrecognized keywords later

                end select

            else
                write(*,*) err // err_inpnn // 'iostat = ', ios
                stop
            end if
        end do

        close(inpnn_unit)

        if (rinpparam%maxnum_layers_short_atomic .gt. 0) then
            do general_counter_1 = 1,rinpparam%maxnum_layers_short_atomic
                rinpparam%maxnodes_short_atomic = max(rinpparam%maxnodes_short_atomic, rinpparam%nodes_short_local(general_counter_1))
            end do
        end if

        if (rinpparam%maxnum_layers_elec .gt. 0) then
            do general_counter_1 = 1,rinpparam%maxnum_layers_elec
                rinpparam%maxnodes_elec = max(rinpparam%maxnodes_elec, rinpparam%nodes_ewald_local(general_counter_1))
            end do
        end if

        if (allocated(rinpparam%nodes_ewald_local)) deallocate(rinpparam%nodes_ewald_local)
        if (allocated(rinpparam%nodes_short_local)) deallocate(rinpparam%nodes_short_local)
        !if (allocated(rinpparam%nodes_pair_local)) deallocate(rinpparam%nodes_pair_local)

        do general_counter_1 = 1,102
            rinpparam%maxnum_funcvalues_short_atomic = max(rinpparam%maxnum_funcvalues_short_atomic, rinpparam%num_funcvalues_local(general_counter_1))
            rinpparam%maxnum_funcvalues_elec = max(rinpparam%maxnum_funcvalues_elec, rinpparam%num_funcvaluese_local(i))
        end do

        deallocate(rinpparam%num_funcvalues_local)
        deallocate(rinpparam%num_funcvaluese_local)
        !deallocate(rinpparam%num_funcvaluesp_local) ! not needed

        deallocate(rinpparam%nucelem)
        deallocate(rinpparam%element)
        !end readout according to getdimensions.f90

        !start readout according to paircount.f90
        if (rinpparam%nn_type_short == 1) then

            call open_for_read(inpnn_unit, filename_inpnn); ios = 0

            do while (ios == 0)
                read(inpnn_unit, '(A)', iostat=ios) buffer
                if (ios == 0) then
                    !line = line + 1
                    call split_string(buffer, words, nwords)

                    select case (words(1))

                        case ('global_symfunction_short')
                            read(words(2),'(i1000)', iostat=ios) rinpparam%function_type_temp
                            if (ios /= 0) stop err // err_inpnn // "global_symfunction_short first argument value must be integer"
                            call lower_case(words(2))
                            select case (words(2))

                                case (1)
                                    if (nwords == 3)
                                        read(words(3),*, iostat=ios) rinpparam%funccutoff_local
                                        if (ios /= 0) stop err // err_inpnn // "global_symfunction_short type ", words(2), " argument ", nwords-1, " must be a number"
                                    else
                                        print *, err, err_inpnn, "global_symfunction_short type ", words(2), " needs 2 arguments"; stop
                                    end if

                                case (2)
                                    if (nwords == 5)
                                        read(words(5),*, iostat=ios) rinpparam%funccutoff_local
                                        if (ios /= 0) stop err // err_inpnn // "global_symfunction_short type ", words(2), " argument ", nwords-1, " must be a number"
                                    else
                                        print *, err, err_inpnn, "global_symfunction_short type ", words(2), " needs 4 arguments"; stop
                                    end if

                                case (3)
                                    if (nwords == 6)
                                        read(words(6),*, iostat=ios) rinpparam%funccutoff_local
                                        if (ios /= 0) stop err // err_inpnn // "global_symfunction_short type ", words(2), " argument ", nwords-1, " must be a number"
                                    else
                                        print *, err, err_inpnn, "global_symfunction_short type ", words(2), " needs 5 arguments"; stop
                                    end if

                                case (4)
                                    if (nwords == 4)
                                        read(words(4),*, iostat=ios) rinpparam%funccutoff_local
                                        if (ios /= 0) stop err // err_inpnn // "global_symfunction_short type ", words(2), " argument ", nwords-1, " must be a number"
                                    else
                                        print *, err, err_inpnn, "global_symfunction_short type ", words(2), " needs 3 arguments"; stop
                                    end if

                                case (8)
                                    if (nwords == 5)
                                        read(words(5),*, iostat=ios) rinpparam%funccutoff_local
                                        if (ios /= 0) stop err // err_inpnn // "global_symfunction_short type ", words(2), " argument ", nwords-1, " must be a number"
                                    else
                                        print *, err, err_inpnn, "global_symfunction_short type ", words(2), " needs 4 arguments"; stop
                                    end if

                                case (9)
                                    if (nwords == 6)
                                        read(words(6),*, iostat=ios) rinpparam%funccutoff_local
                                        if (ios /= 0) stop err // err_inpnn // "global_symfunction_short type ", words(2), " argument ", nwords-1, " must be a number"
                                    else
                                        print *, err, err_inpnn, "global_symfunction_short type ", words(2), " needs 5 arguments"; stop
                                    end if

                                case default
                                    print *, err, err_inpnn, "Error in global_symfunction_short key, symfunction type ", words(2), " not implemented"
                                    stop

                            end select

                        case ('element_symfunction_short')
                            read(words(2),'(A)') rinpparam%elementtemp1
                            read(words(3),'(i1000)', iostat=ios) rinpparam%function_type_temp
                            if (ios /= 0) stop err // err_inpnn // "element_symfunction_short second argument value must be integer"
                            call lower_case(words(2))

                            select case (words(3))

                                case (1)
                                    if (nwords == 4)
                                        read(words(4),*, iostat=ios) rinpparam%funccutoff_local
                                        if (ios /= 0) stop err // err_inpnn // "element_symfunction_short type ", words(3), " argument ", nwords-1, " must be a number"
                                    else
                                        print *, err, err_inpnn, "element_symfunction_short type ", words(3), " needs 3 arguments"; stop
                                    end if

                                case (2)
                                    if (nwords == 6)
                                        read(words(6),*, iostat=ios) rinpparam%funccutoff_local
                                        if (ios /= 0) stop err // err_inpnn // "element_symfunction_short type ", words(3), " argument ", nwords-1, " must be a number"
                                    else
                                        print *, err, err_inpnn, "element_symfunction_short type ", words(3), " needs 5 arguments"; stop
                                    end if

                                case (3)
                                    if (nwords == 7)
                                        read(words(7),*, iostat=ios) rinpparam%funccutoff_local
                                        if (ios /= 0) stop err // err_inpnn // "element_symfunction_short type ", words(3), " argument ", nwords-1, " must be a number"
                                    else
                                        print *, err, err_inpnn, "element_symfunction_short type ", words(3), " needs 6 arguments"; stop
                                    end if

                                case (4)
                                    if (nwords == 5)
                                        read(words(5),*, iostat=ios) rinpparam%funccutoff_local
                                        if (ios /= 0) stop err // err_inpnn // "element_symfunction_short type ", words(3), " argument ", nwords-1, " must be a number"
                                    else
                                        print *, err, err_inpnn, "element_symfunction_short type ", words(3), " needs 4 arguments"; stop
                                    end if

                                case (8)
                                    if (nwords == 6)
                                        read(words(6),*, iostat=ios) rinpparam%funccutoff_local
                                        if (ios /= 0) stop err // err_inpnn // "element_symfunction_short type ", words(3), " argument ", nwords-1, " must be a number"
                                    else
                                        print *, err, err_inpnn, "element_symfunction_short type ", words(3), " needs 5 arguments"; stop
                                    end if

                                case (9)
                                    if (nwords == 7)
                                        read(words(7),*, iostat=ios) rinpparam%funccutoff_local
                                        if (ios /= 0) stop err // err_inpnn // "element_symfunction_short type ", words(3), " argument ", nwords-1, " must be a number"
                                    else
                                        print *, err, err_inpnn, "element_symfunction_short type ", words(3), " needs 6 arguments"; stop
                                    end if

                                case default
                                    print *, err, err_inpnn, "Error in element_symfunction_short key, symfunction type ", words(3), " not implemented"
                                    stop

                            end select

                        case ('symfunction_short')
                            read((words(2),'(A)') rinpparam%elementtemp1
                            read(words(3),'(i1000)', iostat=ios) rinpparam%function_type_temp
                            if (ios /= 0) stop err // err_inpnn // "symfunction_short second argument value must be integer"
                            call lower_case(words(2))

                            select case (words(3))

                                case (1)
                                    if (nwords == 5)
                                        read(words(4),'(A)') rinpparam%elementtemp2
                                        read(words(5),*, iostat=ios) rinpparam%funccutoff_local
                                        if (ios /= 0) stop err // err_inpnn // "symfunction_short type ", words(3), " argument ", nwords-1, " must be a number"
                                    else
                                        print *, err, err_inpnn, "symfunction_short type ", words(3), " needs 4 arguments"; stop
                                    end if

                                case (2)
                                    if (nwords == 7)
                                        read(words(4),'(A)') rinpparam%elementtemp2
                                        read(words(7),*, iostat=ios) rinpparam%funccutoff_local
                                        if (ios /= 0) stop err // err_inpnn // "symfunction_short type ", words(3), " argument ", nwords-1, " must be a number"
                                    else
                                        print *, err, err_inpnn, "symfunction_short type ", words(3), " needs 6 arguments"; stop
                                    end if

                                case (3)
                                    if (nwords == 9)
                                        read(words(4),'(A)') rinpparam%elementtemp2
                                        read(words(5),'(A)') rinpparam%elementtemp3
                                        read(words(9),*, iostat=ios) rinpparam%funccutoff_local
                                        if (ios /= 0) stop err // err_inpnn // "symfunction_short type ", words(3), " argument ", nwords-1, " must be a number"
                                    else
                                        print *, err, err_inpnn, "symfunction_short type ", words(3), " needs 8 arguments"; stop
                                    end if

                                case (4)
                                    if (nwords == 6)
                                        read(words(4),'(A)') rinpparam%elementtemp2
                                        read(words(6),*, iostat=ios) rinpparam%funccutoff_local
                                        if (ios /= 0) stop err // err_inpnn // "symfunction_short type ", words(3), " argument ", nwords-1, " must be a number"
                                    else
                                        print *, err, err_inpnn, "symfunction_short type ", words(3), " needs 5 arguments"; stop
                                    end if

                                case (8)
                                    if (nwords == 8)
                                        read(words(4),'(A)') rinpparam%elementtemp2
                                        read(words(5),'(A)') rinpparam%elementtemp3
                                        read(words(8),*, iostat=ios) rinpparam%funccutoff_local
                                        if (ios /= 0) stop err // err_inpnn // "symfunction_short type ", words(3), " argument ", nwords-1, " must be a number"
                                    else
                                        print *, err, err_inpnn, "symfunction_short type ", words(3), " needs 7 arguments"; stop
                                    end if

                                case (9)
                                    if (nwords == 9)
                                        read(words(4),'(A)') rinpparam%elementtemp2
                                        read(words(5),'(A)') rinpparam%elementtemp3
                                        read(words(9),*, iostat=ios) rinpparam%funccutoff_local
                                        if (ios /= 0) stop err // err_inpnn // "symfunction_short type ", words(3), " argument ", nwords-1, " must be a number"
                                    else
                                        print *, err, err_inpnn, "symfunction_short type ", words(3), " needs 8 arguments"; stop
                                    end if

                                case default
                                    print *, err, err_inpnn, "Error in symfunction_short key, symfunction type ", words(3), " not implemented"
                                    stop

                            end select

                        ! for every other keyword pass here, check for unrecognized keywords later
                        case default
                            !if (trim(words(1)) /= '' .and. words(1)(1:1) /= '#') & ! check for empty and comment lines
                                !print *, warn_inpnn, 'Skipping invalid label ', trim(words(1)),' in line ', line

                    end select

                    rinpparam%maxcutoff_local = max(rinpparam%maxcutoff_local, rinpparam%funccutoff_local)

                else
                    write(*,*) err // err_inpnn // 'iostat = ', ios
                    stop
                end if
            end do

            close(inpnn_unit)

        end if

        if (rinpparam%nn_type_elec == 1) .or. (rinpparam%nn_type_elec == 3) .or. (rinpparam%nn_type_elec == 4) then

            call open_for_read(inpnn_unit, filename_inpnn); ios = 0

            do while (ios == 0)
                read(inpnn_unit, '(A)', iostat=ios) buffer
                if (ios == 0) then
                    !line = line + 1
                    call split_string(buffer, words, nwords)

                    select case (words(1))

                        case ('global_symfunction_electrostatic')
                            read(words(2),'(i1000)', iostat=ios) rinpparam%function_type_temp
                            if (ios /= 0) stop err // err_inpnn // "global_symfunction_electrostatic first argument value must be integer"
                            call lower_case(words(2))

                            select case (words(2))

                                case (1)
                                    if (nwords == 3)
                                        read(words(3),*, iostat=ios) rinpparam%funccutoff_local
                                        if (ios /= 0) stop err // err_inpnn // "global_symfunction_electrostatic type ", words(2), " argument ", nwords-1, " must be a number"
                                    else
                                        print *, err, err_inpnn, "global_symfunction_electrostatic type ", words(2), " needs 2 arguments"; stop
                                    end if

                                case (2)
                                    if (nwords == 5)
                                        read(words(5),*, iostat=ios) rinpparam%funccutoff_local
                                        if (ios /= 0) stop err // err_inpnn // "global_symfunction_electrostatic type ", words(2), " argument ", nwords-1, " must be a number"
                                    else
                                        print *, err, err_inpnn, "global_symfunction_electrostatic type ", words(2), " needs 4 arguments"; stop
                                    end if

                                case (3)
                                    if (nwords == 6)
                                        read(words(6),*, iostat=ios) rinpparam%funccutoff_local
                                        if (ios /= 0) stop err // err_inpnn // "global_symfunction_electrostatic type ", words(2), " argument ", nwords-1, " must be a number"
                                    else
                                        print *, err, err_inpnn, "global_symfunction_electrostatic type ", words(2), " needs 5 arguments"; stop
                                    end if

                                case (4)
                                    if (nwords == 4)
                                        read(words(4),*, iostat=ios) rinpparam%funccutoff_local
                                        if (ios /= 0) stop err // err_inpnn // "global_symfunction_electrostatic type ", words(2), " argument ", nwords-1, " must be a number"
                                    else
                                        print *, err, err_inpnn, "global_symfunction_electrostatic type ", words(2), " needs 3 arguments"; stop
                                    end if

                                case (8)
                                    if (nwords == 5)
                                        read(words(5),*, iostat=ios) rinpparam%funccutoff_local
                                        if (ios /= 0) stop err // err_inpnn // "global_symfunction_electrostatic type ", words(2), " argument ", nwords-1, " must be a number"
                                    else
                                        print *, err, err_inpnn, "global_symfunction_electrostatic type ", words(2), " needs 4 arguments"; stop
                                    end if

                                case (9)
                                    if (nwords == 6)
                                        read(words(6),*, iostat=ios) rinpparam%funccutoff_local
                                        if (ios /= 0) stop err // err_inpnn // "global_symfunction_electrostatic type ", words(2), " argument ", nwords-1, " must be a number"
                                    else
                                        print *, err, err_inpnn, "global_symfunction_electrostatic type ", words(2), " needs 6 arguments"; stop
                                    end if

                                case default
                                    print *, err, err_inpnn, "Error in global_symfunction_electrostatic key, symfunction type ", words(2), " not implemented"
                                    stop

                            end select

                        case ('element_symfunction_electrostatic')
                            read(words(2),'(A)') rinpparam%elementtemp1
                            read(words(3),'(i1000)', iostat=ios) rinpparam%function_type_temp
                            if (ios /= 0) stop err // err_inpnn // "element_symfunction_electrostatic second argument value must be integer"
                            call lower_case(words(2))

                            select case (words(3))

                                case (1)
                                    if (nwords == 4)
                                        read(words(4),*, iostat=ios) rinpparam%funccutoff_local
                                        if (ios /= 0) stop err // err_inpnn // "element_symfunction_electrostatic type ", words(3), " argument ", nwords-1, " must be a number"
                                    else
                                        print *, err, err_inpnn, "element_symfunction_electrostatic type ", words(3), " needs 3 arguments"; stop
                                    end if

                                case (2)
                                    if (nwords == 6)
                                        read(words(6),*, iostat=ios) rinpparam%funccutoff_local
                                        if (ios /= 0) stop err // err_inpnn // "element_symfunction_electrostatic type ", words(3), " argument ", nwords-1, " must be a number"
                                    else
                                        print *, err, err_inpnn, "element_symfunction_electrostatic type ", words(3), " needs 5 arguments"; stop
                                    end if

                                case (3)
                                    if (nwords == 7)
                                        read(words(7),*, iostat=ios) rinpparam%funccutoff_local
                                        if (ios /= 0) stop err // err_inpnn // "element_symfunction_electrostatic type ", words(3), " argument ", nwords-1, " must be a number"
                                    else
                                        print *, err, err_inpnn, "element_symfunction_electrostatic type ", words(3), " needs 6 arguments"; stop
                                    end if

                                case (4)
                                    if (nwords == 5)
                                        read(words(5),*, iostat=ios) rinpparam%funccutoff_local
                                        if (ios /= 0) stop err // err_inpnn // "element_symfunction_electrostatic type ", words(3), " argument ", nwords-1, " must be a number"
                                    else
                                        print *, err, err_inpnn, "element_symfunction_electrostatic type ", words(3), " needs 4 arguments"; stop
                                    end if

                                case (8)
                                    if (nwords == 6)
                                        read(words(6),*, iostat=ios) rinpparam%funccutoff_local
                                        if (ios /= 0) stop err // err_inpnn // "element_symfunction_electrostatic type ", words(3), " argument ", nwords-1, " must be a number"
                                    else
                                        print *, err, err_inpnn, "element_symfunction_electrostatic type ", words(3), " needs 5 arguments"; stop
                                    end if

                                case (9)
                                    if (nwords == 7)
                                        read(words(7),*, iostat=ios) rinpparam%funccutoff_local
                                        if (ios /= 0) stop err // err_inpnn // "element_symfunction_electrostatic type ", words(3), " argument ", nwords-1, " must be a number"
                                    else
                                        print *, err, err_inpnn, "element_symfunction_electrostatic type ", words(3), " needs 6 arguments"; stop
                                    end if

                                case default
                                    print *, err, err_inpnn, "Error in element_symfunction_electrostatic key, symfunction type ", words(3), " not implemented"
                                    stop

                            end select

                        case ('symfunction_electrostatic')
                            read((words(2),'(A)') rinpparam%elementtemp1
                            read(words(3),'(i1000)', iostat=ios) rinpparam%function_type_temp
                            if (ios /= 0) stop err // err_inpnn // "symfunction_electrostatic second argument value must be integer"
                            call lower_case(words(2))

                            select case (words(3))

                                case (1)
                                    if (nwords == 5)
                                        read(words(4),'(A)') rinpparam%elementtemp2
                                        read(words(5),*, iostat=ios) rinpparam%funccutoff_local
                                        if (ios /= 0) stop err // err_inpnn // "symfunction_electrostatic type ", words(3), " argument ", nwords-1, " must be a number"
                                    else
                                        print *, err, err_inpnn, "symfunction_electrostatic type ", words(3), " needs 4 arguments"; stop
                                    end if

                                case (2)
                                    if (nwords == 7)
                                        read(words(4),'(A)') rinpparam%elementtemp2
                                        read(words(7),*, iostat=ios) rinpparam%funccutoff_local
                                        if (ios /= 0) stop err // err_inpnn // "symfunction_electrostatic type ", words(3), " argument ", nwords-1, " must be a number"
                                    else
                                        print *, err, err_inpnn, "symfunction_electrostatic type ", words(3), " needs 6 arguments"; stop
                                    end if

                                case (3)
                                    if (nwords == 9)
                                        read(words(4),'(A)') rinpparam%elementtemp2
                                        read(words(5),'(A)') rinpparam%elementtemp3
                                        read(words(9),*, iostat=ios) rinpparam%funccutoff_local
                                        if (ios /= 0) stop err // err_inpnn // "symfunction_electrostatic type ", words(3), " argument ", nwords-1, " must be a number"
                                    else
                                        print *, err, err_inpnn, "symfunction_electrostatic type ", words(3), " needs 8 arguments"; stop
                                    end if

                                case (4)
                                    if (nwords == 7)
                                        read(words(4),'(A)') rinpparam%elementtemp2
                                        read(words(7),*, iostat=ios) rinpparam%funccutoff_local
                                        if (ios /= 0) stop err // err_inpnn // "symfunction_electrostatic type ", words(3), " argument ", nwords-1, " must be a number"
                                    else
                                        print *, err, err_inpnn, "symfunction_electrostatic type ", words(3), " needs 6 arguments"; stop
                                    end if

                                case (8)
                                    if (nwords == 8)
                                        read(words(4),'(A)') rinpparam%elementtemp2
                                        read(words(5),'(A)') rinpparam%elementtemp3
                                        read(words(8),*, iostat=ios) rinpparam%funccutoff_local
                                        if (ios /= 0) stop err // err_inpnn // "symfunction_electrostatic type ", words(3), " argument ", nwords-1, " must be a number"
                                    else
                                        print *, err, err_inpnn, "symfunction_electrostatic type ", words(3), " needs 7 arguments"; stop
                                    end if

                                case (9)
                                    if (nwords == 9)
                                        read(words(4),'(A)') rinpparam%elementtemp2
                                        read(words(5),'(A)') rinpparam%elementtemp3
                                        read(words(9),*, iostat=ios) rinpparam%funccutoff_local
                                        if (ios /= 0) stop err // err_inpnn // "symfunction_electrostatic type ", words(3), " argument ", nwords-1, " must be a number"
                                    else
                                        print *, err, err_inpnn, "symfunction_electrostatic type ", words(3), " needs 8 arguments"; stop
                                    end if

                                case default
                                    print *, err, err_inpnn, "Error in symfunction_electrostatic key, symfunction type ", words(3), " not implemented"; stop

                            end select

                        ! for every other keyword pass here, check for unrecognized keywords later
                        case default
                            !if (trim(words(1)) /= '' .and. words(1)(1:1) /= '#') & ! check for empty and comment lines
                                !print *, warn_inpnn, 'Skipping invalid label ', trim(words(1)),' in line ', line

                    end select

                    rinpparam%maxcutoff_local = max(rinpparam%maxcutoff_local, rinpparam%funccutoff_local)

                else
                    write(*,*) err // err_inpnn // 'iostat = ', ios
                    stop
                end if
            end do

            close(inpnn_unit)

        end if

        if (rinpparam%maxcutoff_local == default_real) then
            print *, err, err_inpnn, "maxcutoff_local is not set, specify symmetry functions"
            stop
        end if

        rinpparam%max_num_pairs = 0
        ! end readout according to paircount.f90

        ! start checkstructures.f90

        !call checkonestructure(i1,lelement)

        ! according to checkonestructure.f90
        atoms%simbox(3,3) -> read(dataunit,*,err=90)keyword,(lattice(nlattice,i),i=1,3) (nlattice = 1) ! don't forget to convert
        atoms%r(:,:,:) -> ! don't forget to convert

        if(keyword.eq.'lattice') then
            nlattice=nlattice+1
            backspace(dataunit)
            read(dataunit,*,err=90)keyword,(lattice(nlattice,i),i=1,3)
        endif

        if(keyword.eq.'atom') then
            num_atoms=num_atoms+1
            backspace(dataunit)
            read(dataunit,*,err=91)keyword,(xyzstruct(i,num_atoms),i=1,3),&
                elementsymbol(num_atoms),atomcharge(num_atoms),&
                atomenergy(num_atoms),(totalforce(i,num_atoms),i=1,3)
            call nuccharge(elementsymbol(num_atoms),zelem(num_atoms))
            lelement(zelem(num_atoms))=.true. ! element found
        endif

        !! check if lattice vectors make sense
        if(lperiodic)then
            call getvolume(lattice,volume)
            if(volume.lt.0.0000001d0)then
                write(ounit,*)'ERROR: volume of a periodic structure is very small ',volume
                stop
            endif
        endif

        ielem=0
        do i1=1,102
            if(lelement(i1)) ielem=ielem+1
        enddo
        ! end checkstructures.f90

        ! further readout according to initnn.f90

        call distribute_nnflags() ! check if this is really needed

        if(rinpparam%lshort.and.(rinpparam%nn_type_short.eq.1))then
        allocate (rinpparam%num_funcvalues_short_atomic(rinpparam%nelem))
        rinpparam%num_funcvalues_short_atomic(:)=0
        allocate (rinpparam%windex_short_atomic(2*rinpparam%maxnum_layers_short_atomic,rinpparam%nelem))
        allocate (rinpparam%num_layers_short_atomic(rinpparam%nelem))
        rinpparam%num_layers_short_atomic(:)=rinpparam%maxnum_layers_short_atomic
        allocate (rinpparam%actfunc_short_atomic(rinpparam%maxnodes_short_atomic,rinpparam%maxnum_layers_short_atomic,rinpparam%nelem))
        allocate (rinpparam%nodes_short_atomic(0:rinpparam%maxnum_layers_short_atomic,rinpparam%nelem))
        rinpparam%nodes_short_atomic(:,:)=0
        allocate (rinpparam%num_weights_short_atomic(rinpparam%nelem))
        rinpparam%num_weights_short_atomic(:)=0
        end if

        if(rinpparam%lelec.and.(rinpparam%nn_type_elec.eq.1))then
        allocate (rinpparam%num_funcvalues_elec(rinpparam%nelem))
        rinpparam%num_funcvalues_elec(:)=0
        allocate (rinpparam%windex_elec(2*rinpparam%maxnum_layers_elec,rinpparam%nelem))
        allocate (rinpparam%num_layers_elec(rinpparam%nelem))
        rinpparam%num_layers_elec(:)=rinpparam%maxnum_layers_elec
        allocate (rinpparam%actfunc_elec(rinpparam%maxnodes_elec,rinpparam%maxnum_layers_elec,rinpparam%nelem))
        allocate (rinpparam%nodes_elec(0:rinpparam%maxnum_layers_elec,rinpparam%nelem))
        rinpparam%nodes_elec(:,:)=0
        allocate (rinpparam%num_weights_elec(rinpparam%nelem))
        rinpparam%num_weights_elec(:)=0
        endif

        allocate (rinpparam%fixedcharge(rinpparam%nelem))
        rinpparam%fixedcharge(:)=0.0d0
        allocate (rinpparam%nucelem(rinpparam%nelem))
        allocate (rinpparam%element(rinpparam%nelem))
        allocate (rinpparam%atomrefenergies(rinpparam%nelem))
        allocate (rinpparam%elempair(rinpparam%npairs,2))
        rinpparam%elempair(:,:)=0

        call allocatesymfunctions()

        !call readinput(ielem,iseed,lelement) !ielem iseed defined in main.f90/initnn.f90 -> I defined it in get_defaults.f90

        ! start readout of input.nn according to readinput.f90


            call initializecounters() ! even if we use default values, this should stay!!

            if(rinpparam%lshort.and.(rinpparam%nn_type_short.eq.1))then
                rinpparam%nodes_short_atomic_temp(:)   =0
                rinpparam%actfunc_short_atomic_dummy(:)=' '
            endif
            endif
            if(rinpparam%lelec.and.(rinpparam%nn_type_elec.eq.1))then
                rinpparam%nodes_elec_temp(:)           =0
                rinpparam%actfunc_elec_dummy(:)        =' '
            endif
            rinpparam%kalmanlambda_local =0.98000d0
            rinpparam%kalmanlambdae_local=0.98000d0
            rinpparam%iseed=200

            call inputnndefaults() ! own subroutine in pes_nene_mod_supply.f90

            if(rinpparam%lshort.and.(rinpparam%nn_type_short.eq.1))then
                rinpparam%windex_short_atomic(:,:)    =0
                rinpparam%num_weights_short_atomic(:) =0
                rinpparam%maxnum_weights_short_atomic =0
            endif
            if(rinpparam%lelec.and.(rinpparam%nn_type_elec.eq.1))then
                rinpparam%windex_elec(:,:)            =0
                rinpparam%num_weights_elec(:)         =0
                rinpparam%maxnum_weights_elec         =0
            endif

            !call readkeywords(rinpparam%iseed, rinpparam%nodes_short_atomic_temp,rinpparam%nodes_elec_temp,rinpparam%nodes_short_pair_temp, rinpparam%kalmanlambda_local,rinpparam%kalmanlambdae_local)


            ! start readout according to readkeywords.f90
            call open_for_read(inpnn_unit, filename_inpnn); ios = 0

            do while (ios == 0)
                read(inpnn_unit, '(A)', iostat=ios) buffer
                if (ios == 0) then
                    line = line + 1
                    call split_string(buffer, words, nwords)

                    select case (words(1))

                        ! add already read keywords here so that a keyword which is NOT related to RuNNer will be recognized!!
                        case ('nn_type_short') !in readkeywords.f90
                            ! do nothing here, just let it pass

                        case ('runner_mode') !in readkeywords.f90
                            ! do nothing here, just let it pass

                        case ('use_short_nn') !in readkeywords.f90
                            ! do nothing here, just let it pass

                        case ('use_electrostatics') !in readkeywords.f90
                            ! do nothing here, just let it pass

                        case ('electrostatic_type', 'nn_type_elec') !in readkeywords.f90
                            ! do nothing here, just let it pass

                        case ('debug_mode') !in readkeywords.f90
                            ! do nothing here, just let it pass

                        case ('global_hidden_layers_short') !in readkeywords.f90
                            ! do nothing here, just let it pass

                        case ('global_hidden_layers_electrostatic') !in readkeywords.f90
                            ! do nothing here, just let it pass

                        !case ('global_hidden_layers_pair') !in readkeywords.f90
                            ! do nothing here, just let it pass

                        !case ('use_atom_energies') !in readkeywords.f90
                            ! do nothing here, just let it pass

                        !case ('use_atom_charges') !in readkeywords.f90
                            ! do nothing here, just let it pass

                        case ('number_of_elements') !in readkeywords.f90
                            ! do nothing here, just let it pass

                        case ('elements') !in readkeywords.f90
                            ! do nothing here, just let it pass

                        !case ('global_nodes_short')
                            ! do nothing here, just let it pass

                        !case ('global_nodes_electrostatic')
                            ! do nothing here, just let it pass

                        !case ('global_nodes_pair')
                            ! do nothing here, just let it pass

                        case ('element_symfunction_short') !in readkeywords.f90
                            ! do nothing here, just let it pass

                        case ('element_symfunction_electrostatic') !in readkeywords.f90
                            ! do nothing here, just let it pass

                        case ('global_symfunction_short', 'global_symfunction_short_atomic') !in readkeywords.f90
                            ! do nothing here, just let it pass

                        case ('global_symfunction_electrostatic', 'global_symfunction_elec') !in readkeywords.f90
                            ! do nothing here, just let it pass

                        case ('symfunction_short') !in readkeywords.f90
                            ! do nothing here, just let it pass

                        case ('symfunction_electrostatic') !in readkeywords.f90
                            ! do nothing here, just let it pass

                        !case ('pairsymfunction_short') !in readkeywords.f90
                            ! done before, not needed here anymore

                        !case ('element_pairsymfunction_short') !in readkeywords.f90
                            ! done before, not needed here anymore

                        !case ('global_pairsymfunction_short') !in readkeywords.f90
                            ! done before, not needed here anymore


                        case ('check_input_forces')
                            if (rinpparam%inputforcethreshold /= default_real) stop err // err_inpnn // 'Multiple use of the check_input_forces key'
                            if (nwords == 2) then
                                rinpparam%lcheckinputforces = .true.
                                read(words(2),*, iostat=ios) rinpparam%inputforcethreshold
                                if (ios /= 0) stop err // err_inpnn // "check_input_forces value must be a number"
                            else
                                print *, err, err_inpnn, "check_input_forces key needs a single argument"; stop
                            end if

                        case ('print_force_components')
                            if (rinpparam%lprintforcecomponents /= default_bool) stop err // err_inpnn // 'Multiple use of the print_force_components key'
                            if (nwords == 1) then
                                rinpparam%lprintforcecomponents = .true.
                            else
                                print *, err, err_inpnn, "print_force_components key needs no argument(s)"; stop
                            end if

                        case ('ion_forces_only')
                            if (rinpparam%lionforcesonly /= default_int) stop err // err_inpnn // 'Multiple use of the ion_forces_only key'
                            if (nwords == 1) then
                                rinpparam%lionforcesonly = .true.
                            else
                                print *, err, err_inpnn, "ion_forces_only key needs no argument(s)"; stop
                            end if

                        case ('use_electrostatic_nn')
                            print *, err, err_inpnn, "use_electrostatic_nn key is obsolete, please use electrostatic_type and use_electrostatics instead"; stop

!                       case ('debug_mode')
!                           if (rinpparam%ldebug /= default_bool) stop err // err_inpnn // 'Multiple use of the debug_mode key'
!                           if (nwords == 1) then
!                               rinpparam%ldebug = .true.
!                           else
!                               print *, err, err_inpnn, "debug_mode key needs no argument(s)"; stop
!                           end if

                        case ('cutoff_type')
                            if (rinpparam%cutoff_type /= default_int) stop err // err_inpnn // 'Multiple use of the cutoff_type key'
                            if (nwords == 2) then
                                read(words(2),'(i1000)', iostat=ios) rinpparam%cutoff_type
                                if (ios /= 0) stop err // err_inpnn // "cutoff_type value must be integer"
                            else
                                print *, err, err_inpnn, "cutoff_type key needs a single argument"; stop
                            end if

                        case ('cutoff_alpha')
                            if (rinpparam%cutoff_alpha /= default_real) stop err // err_inpnn // 'Multiple use of the cutoff_alpha key'
                            if (nwords == 2) then
                                read(words(2),*, iostat=ios) rinpparam%cutoff_alpha
                                if (ios /= 0) stop err // err_inpnn // "cutoff_alpha value must be a number"
                            else
                                print *, err, err_inpnn, "cutoff_alpha key needs a single argument"; stop
                            end if

                        case ('dynamic_force_grouping')
                            if (rinpparam%ldynforcegroup /= default_bool) stop err // err_inpnn // 'Multiple use of the dynamic_force_grouping key'
                            if (nwords == 3) then
                                rinpparam%ldynforcegroup = .true.
                                read(words(2),'(i1000)', iostat=ios) rinpparam%dynforcegroup_start
                                if (ios /= 0) stop err // err_inpnn // "dynamic_force_grouping first argument value must be integer"
                                read(words(3),'(i1000)', iostat=ios) rinpparam%dynforcegroup_step
                                if (ios /= 0) stop err // err_inpnn // "dynamic_force_grouping second argument value must be integer"
                            else
                                print *, err, err_inpnn, "dynamic_force_grouping key needs a single argument"; stop
                            end if

                        case ('detect_saturation')
                            if (rinpparam%ldetect_saturation /= default_bool) stop err // err_inpnn // 'Multiple use of the detect_saturation key'
                            if (nwords == 2) then
                                rinpparam%ldetect_saturation = .true.
                                read(words(2),*, iostat=ios) rinpparam%saturation_threshold
                                if (ios /= 0) stop err // err_inpnn // "detect_saturation value must be a number"
                            else
                                print *, err, err_inpnn, "detect_saturation key needs a single argument"; stop
                            end if

                        case ('data_clustering')
                            if (rinpparam%ldataclustering /= default_bool) stop err // err_inpnn // 'Multiple use of the data_clustering key'
                            if (nwords == 3) then
                                rinpparam%ldataclustering=.true.
                                read(words(2),*, iostat=ios) rinpparam%dataclusteringthreshold1
                                if (ios /= 0) stop err // err_inpnn // "data_clustering value must be a number"
                                read(words(3),*, iostat=ios) rinpparam%dataclusteringthreshold2
                                if (ios /= 0) stop err // err_inpnn // "data_clustering value must be a number"
                            else
                                print *, err, err_inpnn, "data_clustering key needs two arguments"; stop
                            end if

                        case ('analyze_error_energy_step')
                            if (rinpparam%analyze_error_energy_step /= default_real) stop err // err_inpnn // 'Multiple use of the analyze_error_energy_step key'
                            if (nwords == 2) then
                                read(words(2),*, iostat=ios) rinpparam%analyze_error_energy_step
                                if (ios /= 0) stop err // err_inpnn // "analyze_error_energy_step value must be a number"
                            else
                                print *, err, err_inpnn, "analyze_error_energy_step key needs a single argument"; stop
                            end if

                        case ('analyze_error_force_step')
                            if (rinpparam%analyze_error_force_step /= default_real) stop err // err_inpnn // 'Multiple use of the analyze_error_force_step key'
                            if (nwords == 2) then
                                read(words(2),*, iostat=ios) rinpparam%analyze_error_force_step
                                if (ios /= 0) stop err // err_inpnn // "analyze_error_force_step value must be a number"
                            else
                                print *, err, err_inpnn, "analyze_error_force_step key needs a single argument"; stop
                            end if

                        case ('analyze_error_charge_step')
                            if (rinpparam%analyze_error_charge_step /= default_real) stop err // err_inpnn // 'Multiple use of the analyze_error_charge_step key'
                            if (nwords == 2) then
                                read(words(2),*, iostat=ios) rinpparam%analyze_error_charge_step
                                if (ios /= 0) stop err // err_inpnn // "analyze_error_charge_step value must be a number"
                            else
                                print *, err, err_inpnn, "analyze_error_charge_step key needs a single argument"; stop
                            end if

                        case ('parallel_mode')
                            if (rinpparam%paramode /= default_int) stop err // err_inpnn // 'Multiple use of the parallel_mode key'
                            if (nwords == 2) then
                                read(words(2),'(i1000)', iostat=ios) rinpparam%paramode
                                if (ios /= 0) stop err // err_inpnn // "parallel_mode value must be integer"
                                if (rinpparam%paramode /= 1) stop err // err_inpnn // "parallel_mode ", rinpparam%paramode, " not implemented, only parallel_mode 1 (serial version) available"
                            else
                                print *, err, err_inpnn, "parallel_mode key needs a single argument"; stop
                            end if

                        case ('symfunction_correlation')
                            if (rinpparam%lpearson_correlation /= default_bool) stop err // err_inpnn // 'Multiple use of the symfunction_correlation key'
                            if (nwords == 1) then
                                rinpparam%lpearson_correlation = .true.
                            else
                                print *, err, err_inpnn, "symfunction_correlation key needs no argument(s)"; stop
                            end if

                        case ('weight_analysis')
                            if (rinpparam%lweightanalysis /= default_bool) stop err // err_inpnn // 'Multiple use of the weight_analysis key'
                            if (nwords == 1) then
                                rinpparam%lweightanalysis = .true.
                            else
                                print *, err, err_inpnn, "weight_analysis key needs no argument(s)"; stop
                            end if

                        case ('environment_analysis')
                            if (rinpparam%lenvironmentanalysis /= default_bool) stop err // err_inpnn // 'Multiple use of the environment_analysis key'
                            if (nwords == 1) then
                                rinpparam%lenvironmentanalysis = .true.
                            else
                                print *, err, err_inpnn, "environment_analysis key needs no argument(s)"; stop
                            end if

                        case ('find_contradictions')
                            if (rinpparam%lfindcontradictions /= default_bool) stop err // err_inpnn // 'Multiple use of the find_contradictions key'
                            if (nwords == 3) then
                                rinpparam%lfindcontradictions = .true.
                                read(words(2),*, iostat=ios) rinpparam%deltagthres
                                if (ios /= 0) stop err // err_inpnn // "find_contradictions first argument value must be a number"
                                read(words(3),*, iostat=ios) rinpparam%deltafthres
                                if (ios /= 0) stop err // err_inpnn // "find_contradictions second argument value must be a number"
                            else
                                print *, err, err_inpnn, "find_contradictions key needs two arguments"; stop
                            end if

                        case ('use_old_scaling')
                            if (rinpparam%luseoldscaling /= default_bool) stop err // err_inpnn // 'Multiple use of the use_old_scaling key'
                            if (nwords == 1) then
                                rinpparam%luseoldscaling = .true.
                            else
                                print *, err, err_inpnn, "use_old_scaling key needs no argument(s)"; stop
                            end if

                        case ('md_mode')
                            if (rinpparam%lmd /= default_bool) stop err // err_inpnn // 'Multiple use of the md_mode key'
                            if (nwords == 1) then
                                rinpparam%lmd = .true.
                            else
                                print *, err, err_inpnn, "md_mode key needs no argument(s)"; stop
                            end if

                        !case ('global_nodes_short', 'global_nodes_short_atomic')
                        !case ('global_nodes_short_atomic')
                            !print *, err, err_inpnn, "global_nodes_short_atomic key is obsolete, please use global_nodes_short instead"; stop

                        case ('global_nodes_short')
                            if(rinpparam%lshort .and. (rainpparam%nn_type_short == 1)) then
                                if (nwords == rinpparam%maxnum_layers_short_atomic) then
                                    do general_counter_1 = 1,rinpparam%maxnum_layers_short_atomic-1
                                        read(words(general_counter_1+1),'(i1000)', iostat=ios) rinpparam%nodes_short_atomic_temp(general_counter_1)
                                        if (ios /= 0) stop err // err_inpnn // "global_nodes_short argument ", general_counter_1, " value must be integer"
                                    end do
                                    do general_counter_1 = 1,rinpparam%nelem
                                        do general_counter_2 = 1,rinpparam%maxnum_layers_short_atomic
                                            rinpparam%nodes_short_atomic(general_counter_2,general_counter_1) = rinpparam%nodes_short_atomic_temp(general_counter_2)
                                        end do
                                    end do
                                else
                                    print *, err, err_inpnn, "global_nodes_short key needs ", rinpparam%maxnum_layers_short_atomic-1, " arguments"; stop
                                end if
                            end if

                        case ('global_nodes_electrostatic')
                            if(rinpparam%lelec .and. (rainpparam%nn_type_elec == 1)) then
                                if (nwords == rinpparam%maxnum_layers_elec) then
                                    do general_counter_1 = 1,rinpparam%maxnum_layers_elec-1
                                        read(words(general_counter_1+1),'(i1000)', iostat=ios) rinpparam%nodes_elec_temp(general_counter_1)
                                        if (ios /= 0) stop err // err_inpnn // "global_nodes_electrostatic argument ", general_counter_1, " value must be integer"
                                    end do
                                    do general_counter_1 = 1,rinpparam%nelem
                                        do general_counter_2 = 1,rinpparam%maxnum_layers_elec
                                            rinpparam%nodes_elec(general_counter_2,general_counter_1) = rinpparam%nodes_elec_temp(general_counter_2)
                                        end do
                                    end do
                                else
                                    print *, err, err_inpnn, "global_nodes_electrostatic key needs ", rinpparam%maxnum_layers_elec-1, " arguments"; stop
                                end if
                            end if

                        !case ('global_nodes_short_pair')
                            !print *, err, err_inpnn, "global_nodes_short_pair key is obsolete, please use global_nodes_pair instead"; stop

                        !case ('global_nodes_pair')
                            !print *, err, err_inpnn, "global_nodes_pair key not supported, Pair NN not implemented"; stop

                        case ('global_output_nodes_short')
                            print *, err, err_inpnn, "global_output_nodes_short key is obsolete, please remove it"; stop

                        case ('global_output_nodes_electrostatic')
                            print *, err, err_inpnn, "global_output_nodes_electrostatic key is obsolete, please remove it"; stop

                        case ('global_output_nodes_pair')
                            print *, err, err_inpnn, "global_output_nodes_pair key is obsolete, please remove it"; stop

                        case ('ewald_alpha')
                            if (rinpparam%ewaldalpha /= default_real) stop err // err_inpnn // 'Multiple use of the ewald_alpha key'
                            if (nwords == 2) then
                                read(words(2),*, iostat=ios) rinpparam%ewaldalpha
                                if (ios /= 0) stop err // err_inpnn // "ewald_alpha value must be a number"
                            else
                                print *, err, err_inpnn, "ewald_alpha key needs a single argument"; stop
                            end if

                        case ('ewald_cutoff')
                            if (rinpparam%ewaldcutoff /= default_real) stop err // err_inpnn // 'Multiple use of the ewald_cutoff key'
                            if (nwords == 2) then
                                read(words(2),*, iostat=ios) rinpparam%ewaldcutoff
                                if (ios /= 0) stop err // err_inpnn // "ewald_cutoff value must be a number"
                            else
                                print *, err, err_inpnn, "ewald_cutoff key needs a single argument"; stop
                            end if

                        case ('ewald_kmax')
                            if (rinpparam%ewaldkmax /= default_int) stop err // err_inpnn // 'Multiple use of the ewald_kmax key'
                            if (nwords == 2) then
                                read(words(2),'(i1000)', iostat=ios) rinpparam%ewaldkmax
                                if (ios /= 0) stop err // err_inpnn // "ewald_kmax value must be integer"
                            else
                                print *, err, err_inpnn, "ewald_kmax key needs a single argument"; stop
                            end if

                        case ('precondition_weights')
                            if (rinpparam%lprecond /= default_bool) stop err // err_inpnn // 'Multiple use of the precondition_weights key'
                            if (nwords == 1) then
                                rinpparam%lprecond = .true.
                            else
                                print *, err, err_inpnn, "precondition_weights key needs no argument(s)"; stop
                            end if

                        case ('initialization_only')
                            if (rinpparam%linionly /= default_bool) stop err // err_inpnn // 'Multiple use of the initialization_only key'
                            if (nwords == 1) then
                                rinpparam%linionly = .true.
                            else
                                print *, err, err_inpnn, "initialization_only key needs no argument(s)"; stop
                            end if

                        case ('force_grouping_by_structure')
                            if (rinpparam%lfgroupbystruct /= default_bool) stop err // err_inpnn // 'Multiple use of the force_grouping_by_structure key'
                            if (nwords == 1) then
                                rinpparam%lfgroupbystruct = .true.
                            else
                                print *, err, err_inpnn, "force_grouping_by_structure key needs no argument(s)"; stop
                            end if

                        case ('charge_grouping_by_structure')
                            if (rinpparam%lqgroupbystruct /= default_bool) stop err // err_inpnn // 'Multiple use of the charge_grouping_by_structure key'
                            if (nwords == 1) then
                                rinpparam%lqgroupbystruct = .true.
                            else
                                print *, err, err_inpnn, "charge_grouping_by_structure key needs no argument(s)"; stop
                            end if

                        case ('mix_all_points')
                            if (rinpparam%lmixpoints /= default_bool) stop err // err_inpnn // 'Multiple use of the mix_all_points key'
                            if (nwords == 1) then
                                rinpparam%lmixpoints = .true.
                            else
                                print *, err, err_inpnn, "mix_all_points key needs no argument(s)"; stop
                            end if

                        case ('print_convergence_vector')
                            if (rinpparam%lprintconv /= default_bool) stop err // err_inpnn // 'Multiple use of the print_convergence_vector key'
                            if (nwords == 1) then
                                rinpparam%lprintconv = .true.
                            else
                                print *, err, err_inpnn, "print_convergence_vector key needs no argument(s)"; stop
                            end if

                        case ('print_mad')
                            if (rinpparam%lprintmad /= default_bool) stop err // err_inpnn // 'Multiple use of the print_mad key'
                            if (nwords == 1) then
                                rinpparam%lprintmad = .true.
                            else
                                print *, err, err_inpnn, "print_mad key needs no argument(s)"; stop
                            end if

                        case ('noise_energy')
                            if (rinpparam%noisee /= default_real) stop err // err_inpnn // 'Multiple use of the noise_energy key'
                            if (nwords == 2) then
                                read(words(2),*, iostat=ios) rinpparam%noisee
                                if (ios /= 0) stop err // err_inpnn // "noise_energy value must be a number"
                            else
                                print *, err, err_inpnn, "noise_energy key needs a single argument"; stop
                            end if

                        case ('noise_force')
                            if (rinpparam%noisef /= default_real) stop err // err_inpnn // 'Multiple use of the noise_force key'
                            if (nwords == 2) then
                                read(words(2),*, iostat=ios) rinpparam%noisef
                                if (ios /= 0) stop err // err_inpnn // "noise_force value must be a number"
                            else
                                print *, err, err_inpnn, "noise_force key needs a single argument"; stop
                            end if

                        case ('noise_charge')
                            if (rinpparam%noiseq /= default_real) stop err // err_inpnn // 'Multiple use of the noise_charge key'
                            if (nwords == 2) then
                                read(words(2),*, iostat=ios) rinpparam%noiseq
                                if (ios /= 0) stop err // err_inpnn // "noise_charge value must be a number"
                            else
                                print *, err, err_inpnn, "noise_charge key needs a single argument"; stop
                            end if

                        case ('short_energy_group')
                            if (rinpparam%nenergygroup /= default_int) stop err // err_inpnn // 'Multiple use of the short_energy_group key'
                            if (nwords == 2) then
                                read(words(2),'(i1000)', iostat=ios) rinpparam%nenergygroup
                                if (ios /= 0) stop err // err_inpnn // "short_energy_group value must be integer"
                            else
                                print *, err, err_inpnn, "short_energy_group key needs a single argument"; stop
                            end if

                        case ('short_force_group')
                            if (rinpparam%nforcegroup /= default_int) stop err // err_inpnn // 'Multiple use of the short_force_group key'
                            if (nwords == 2) then
                                read(words(2),'(i1000)', iostat=ios) rinpparam%nforcegroup
                                if (ios /= 0) stop err // err_inpnn // "short_force_group value must be integer"
                            else
                                print *, err, err_inpnn, "short_force_group key needs a single argument"; stop
                            end if

                        case ('charge_group')
                            if (rinpparam%nchargegroup /= default_int) stop err // err_inpnn // 'Multiple use of the charge_group key'
                            if (nwords == 2) then
                                read(words(2),'(i1000)', iostat=ios) rinpparam%nchargegroup
                                if (ios /= 0) stop err // err_inpnn // "charge_group value must be integer"
                            else
                                print *, err, err_inpnn, "charge_group key needs a single argument"; stop
                            end if

                       case ('use_short_forces')
                            if (rinpparam%luseforces /= default_bool) stop err // err_inpnn // 'Multiple use of the use_short_forces key'
                            if (nwords == 1) then
                                rinpparam%luseforces = .true.
                            else
                                print *, err, err_inpnn, "use_short_forces key needs no argument(s)"; stop
                            end if

                        case ('short_energy_fraction')
                            if (rinpparam%energyrnd /= default_real) stop err // err_inpnn // 'Multiple use of the short_energy_fraction key'
                            if (nwords == 2) then
                                read(words(2),*, iostat=ios) rinpparam%energyrnd
                                if (ios /= 0) stop err // err_inpnn // "short_energy_fraction value must be a number"
                            else
                                print *, err, err_inpnn, "short_energy_fraction key needs a single argument"; stop
                            end if

                        case ('short_force_fraction')
                            if (rinpparam%forcernd /= default_real) stop err // err_inpnn // 'Multiple use of the short_force_fraction key'
                            if (nwords == 2) then
                                read(words(2),*, iostat=ios) rinpparam%forcernd
                                if (ios /= 0) stop err // err_inpnn // "short_force_fraction value must be a number"
                            else
                                print *, err, err_inpnn, "short_force_fraction key needs a single argument"; stop
                            end if

                        case ('charge_fraction')
                            if (rinpparam%chargernd /= default_real) stop err // err_inpnn // 'Multiple use of the charge_fraction key'
                            if (nwords == 2) then
                                read(words(2),*, iostat=ios) rinpparam%chargernd
                                if (ios /= 0) stop err // err_inpnn // "charge_fraction value must be a number"
                            else
                                print *, err, err_inpnn, "charge_fraction key needs a single argument"; stop
                            end if

                        case ('use_atom_charges')
                            ! let it pass, since done before
                            !if (rinpparam%luseatomcharges /= default_bool) stop err // err_inpnn // 'Multiple use of the use_atom_charges key'
                            !if (nwords == 1) then
                            !    rinpparam%luseatomcharges = .true.
                            !else
                            !    print *, err, err_inpnn, "use_atom_charges key needs no argument(s)"; stop
                            !end if

                        case ('use_atom_energies')
                            ! let it pass, since done before
                            !if (rinpparam%luseatomenergies /= default_bool) stop err // err_inpnn // 'Multiple use of the use_atom_energies key'
                            !if (nwords == 1) then
                            !    rinpparam%luseatomenergies = .true.
                            !else
                            !    print *, err, err_inpnn, "use_atom_energies key needs no argument(s)"; stop
                            !end if

                        case ('remove_atom_energies')
                            if (rinpparam%lremoveatomenergies /= default_bool) stop err // err_inpnn // 'Multiple use of the remove_atom_energies key'
                            if (nwords == 1) then
                                rinpparam%lremoveatomenergies = .true.
                            else
                                print *, err, err_inpnn, "remove_atom_energies key needs no argument(s)"; stop
                            end if

                        case ('analyze_error')
                            if (rinpparam%lanalyzeerror /= default_bool) stop err // err_inpnn // 'Multiple use of the analyze_error key'
                            if (nwords == 1) then
                                rinpparam%lanalyzeerror = .true.
                            else
                                print *, err, err_inpnn, "analyze_error key needs no argument(s)"; stop
                            end if

                        case ('use_charge_constraint')
                            if (rinpparam%lchargeconstraint /= default_bool) stop err // err_inpnn // 'Multiple use of the use_charge_constraint key'
                            if (nwords == 1) then
                                rinpparam%lchargeconstraint = .true.
                            else
                                print *, err, err_inpnn, "use_charge_constraint key needs no argument(s)"; stop
                            end if

                        case ('fitmode')
                            if (rinpparam%fitmode /= default_int) stop err // err_inpnn // 'Multiple use of the fitmode key'
                            if (nwords == 2) then
                                read(words(2),'(i1000)', iostat=ios) rinpparam%fitmode
                                if (ios /= 0) stop err // err_inpnn // "fitmode value must be integer"

                                select case (words(2))

                                    case (1, 2)
                                        ! Just let it pass

                                    case default
                                        print *, err, err_inpnn, "Error in fitmode key value, only 1 and 2 possible"
                                        stop

                                end select

                            else
                                print *, err, err_inpnn, "fitmode key needs a single argument"; stop
                            end if

                        case ('energy_threshold')
                            if (rinpparam%fitethres /= default_real) stop err // err_inpnn // 'Multiple use of the energy_threshold key'
                            if (nwords == 2) then
                                rinpparam%lfitethres = .true.
                                read(words(2),*, iostat=ios) rinpparam%fitethres
                                if (ios /= 0) stop err // err_inpnn // "energy_threshold value must be a number"
                            else
                                print *, err, err_inpnn, "energy_threshold key needs a single argument"; stop
                            end if

                        case ('force_threshold')
                            if (rinpparam%fitfthres /= default_real) stop err // err_inpnn // 'Multiple use of the force_threshold key'
                            if (nwords == 2) then
                                rinpparam%lfitfthres = .true.
                                read(words(2),*, iostat=ios) rinpparam%fitfthres
                                if (ios /= 0) stop err // err_inpnn // "force_threshold value must be a number"
                            else
                                print *, err, err_inpnn, "force_threshold key needs a single argument"; stop
                            end if

                       case ('bond_threshold')
                            if (rinpparam%rmin /= default_real) stop err // err_inpnn // 'Multiple use of the bond_threshold key'
                            if (nwords == 2) then
                                read(words(2),*, iostat=ios) rinpparam%rmin
                                if (ios /= 0) stop err // err_inpnn // "bond_threshold value must be a number"
                            else
                                print *, err, err_inpnn, "bond_threshold key needs a single argument"; stop
                            end if

                        case ('optmode_short_energy')
                            if (rinpparam%optmodee /= default_int) stop err // err_inpnn // 'Multiple use of the optmode_short_energy key'
                            if (nwords == 2) then
                                read(words(2),'(i1000)', iostat=ios) rinpparam%optmodee
                                if (ios /= 0) stop err // err_inpnn // "optmode_short_energy value must be integer"
                            else
                                print *, err, err_inpnn, "optmode_short_energy key needs a single argument"; stop
                            end if

                        case ('optmode_short_force')
                            if (rinpparam%optmodef /= default_int) stop err // err_inpnn // 'Multiple use of the optmode_short_force key'
                            if (nwords == 2) then
                                read(words(2),'(i1000)', iostat=ios) rinpparam%optmodef
                                if (ios /= 0) stop err // err_inpnn // "optmode_short_force value must be integer"
                            else
                                print *, err, err_inpnn, "optmode_short_force key needs a single argument"; stop
                            end if

                        case ('optmode_charge')
                            if (rinpparam%optmodeq /= default_int) stop err // err_inpnn // 'Multiple use of the optmode_charge key'
                            if (nwords == 2) then
                                read(words(2),'(i1000)', iostat=ios) rinpparam%optmodeq
                                if (ios /= 0) stop err // err_inpnn // "optmode_charge value must be integer"
                            else
                                print *, err, err_inpnn, "optmode_charge key needs a single argument"; stop
                            end if

                        case ('random_seed')
                            if (rinpparam%iseed /= default_int) stop err // err_inpnn // 'Multiple use of the random_seed key'
                            if (nwords == 2) then
                                read(words(2),'(i1000)', iostat=ios) rinpparam%iseed
                                if (ios /= 0) stop err // err_inpnn // "random_seed value must be integer"
                            else
                                print *, err, err_inpnn, "random_seed key needs a single argument"; stop
                            end if

                        case ('points_in_memory', 'nblock') ! think about to set it according to number of atoms from structure file
                            if (rinpparam%nblock /= default_int) stop err // err_inpnn // 'Multiple use of the points_in_memory/nblock key'
                            if (nwords == 2) then
                                read(words(2),'(i1000)', iostat=ios) rinpparam%nblock
                                if (ios /= 0) stop err // err_inpnn // "points_in_memory/nblock value must be integer"
                            else
                                print *, err, err_inpnn, "points_in_memory/nblock key needs a single argument"; stop
                            end if

                        case ('epochs')
                            if (rinpparam%nepochs /= default_int) stop err // err_inpnn // 'Multiple use of the epochs key'
                            if (nwords == 2) then
                                read(words(2),'(i1000)', iostat=ios) rinpparam%nepochs
                                if (ios /= 0) stop err // err_inpnn // "epochs value must be integer"
                            else
                                print *, err, err_inpnn, "epochs key needs a single argument"; stop
                            end if

                        case ('write_weights_epoch')
                            if (rinpparam%iwriteweight /= default_int) stop err // err_inpnn // 'Multiple use of the write_weights_epoch key'
                            if (nwords == 2) then
                                read(words(2),'(i1000)', iostat=ios) rinpparam%iwriteweight
                                if (ios /= 0) stop err // err_inpnn // "write_weights_epoch value must be integer"
                            else
                                print *, err, err_inpnn, "write_weights_epoch key needs a single argument"; stop
                            end if

                        case ('write_temporary_weights')
                            if (rinpparam%lwritetmpweights /= default_bool) stop err // err_inpnn // 'Multiple use of the write_temporary_weights key'
                            if (nwords == 1) then
                                rinpparam%lwritetmpweights = .true.
                            else
                                print *, err, err_inpnn, "write_temporary_weights key needs no argument(s)"; stop
                            end if

                        case ('write_symfunctions')
                            if (rinpparam%lwritesymfunctions /= default_bool) stop err // err_inpnn // 'Multiple use of the write_symfunctions key'
                            if (nwords == 1) then
                                rinpparam%lwritesymfunctions = .true.
                            else
                                print *, err, err_inpnn, "write_symfunctions key needs no argument(s)"; stop
                            end if

                        case ('test_fraction')
                            if (rinpparam%splitthres /= default_real) stop err // err_inpnn // 'Multiple use of the test_fraction key'
                            if (nwords == 2) then
                                read(words(2),*, iostat=ios) rinpparam%splitthres
                                if (ios /= 0) stop err // err_inpnn // "test_fraction value must be a number"
                            else
                                print *, err, err_inpnn, "test_fraction key needs a single argument"; stop
                            end if

                        case ('scale_min_short_atomic')
                            if (rinpparam%scmin_short_atomic /= default_real) stop err // err_inpnn // 'Multiple use of the scale_min_short_atomic key'
                            if (nwords == 2) then
                                read(words(2),*, iostat=ios) rinpparam%scmin_short_atomic
                                if (ios /= 0) stop err // err_inpnn // "scale_min_short_atomic value must be a number"
                            else
                                print *, err, err_inpnn, "scale_min_short_atomic key needs a single argument"; stop
                            end if

                        case ('scale_max_short_atomic')
                            if (rinpparam%scmax_short_atomic /= default_real) stop err // err_inpnn // 'Multiple use of the scale_max_short_atomic key'
                            if (nwords == 2) then
                                read(words(2),*, iostat=ios) rinpparam%scmax_short_atomic
                                if (ios /= 0) stop err // err_inpnn // "scale_max_short_atomic value must be a number"
                            else
                                print *, err, err_inpnn, "scale_max_short_atomic key needs a single argument"; stop
                            end if

                        case ('scale_min_short_pair')
                            print *, err, err_inpnn, "scale_min_short_pair key not supported, Pair NN not implemented"; stop

                        case ('scale_max_short_pair')
                            print *, err, err_inpnn, "scale_max_short_pair key not supported, Pair NN not implemented"; stop

                        case ('scale_min_elec')
                            if (rinpparam%scmin_elec /= default_real) stop err // err_inpnn // 'Multiple use of the scale_min_elec key'
                            if (nwords == 2) then
                                read(words(2),*, iostat=ios) rinpparam%scmin_elec
                                if (ios /= 0) stop err // err_inpnn // "scale_min_elec value must be a number"
                            else
                                print *, err, err_inpnn, "scale_min_elec key needs a single argument"; stop
                            end if

                        case ('scale_max_elec')
                            if (rinpparam%scmax_elec /= default_real) stop err // err_inpnn // 'Multiple use of the scale_max_elec key'
                            if (nwords == 2) then
                                read(words(2),*, iostat=ios) rinpparam%scmax_elec
                                if (ios /= 0) stop err // err_inpnn // "scale_max_elec value must be a number"
                            else
                                print *, err, err_inpnn, "scale_max_elec key needs a single argument"; stop
                            end if

                        case ('short_energy_error_threshold')
                            if (rinpparam%kalmanthreshold /= default_real) stop err // err_inpnn // 'Multiple use of the short_energy_error_threshold key'
                            if (nwords == 2) then
                                read(words(2),*, iostat=ios) rinpparam%kalmanthreshold
                                if (ios /= 0) stop err // err_inpnn // "short_energy_error_threshold value must be a number"
                            else
                                print *, err, err_inpnn, "short_energy_error_threshold key needs a single argument"; stop
                            end if

                        case ('short_force_error_threshold')
                            if (rinpparam%kalmanthresholdf /= default_real) stop err // err_inpnn // 'Multiple use of the short_force_error_threshold key'
                            if (nwords == 2) then
                                read(words(2),*, iostat=ios) rinpparam%kalmanthresholdf
                                if (ios /= 0) stop err // err_inpnn // "short_force_error_threshold value must be a number"
                            else
                                print *, err, err_inpnn, "short_force_error_threshold key needs a single argument"; stop
                            end if

                        case ('charge_error_threshold')
                            if (rinpparam%kalmanthresholde /= default_real) stop err // err_inpnn // 'Multiple use of the charge_error_threshold key'
                            if (nwords == 2) then
                                read(words(2),*, iostat=ios) rinpparam%kalmanthresholde
                                if (ios /= 0) stop err // err_inpnn // "charge_error_threshold value must be a number"
                            else
                                print *, err, err_inpnn, "charge_error_threshold key needs a single argument"; stop
                            end if

                        case ('total_charge_error_threshold')
                            if (rinpparam%kalmanthresholdc /= default_real) stop err // err_inpnn // 'Multiple use of the total_charge_error_threshold key'
                            if (nwords == 2) then
                                read(words(2),*, iostat=ios) rinpparam%kalmanthresholdc
                                if (ios /= 0) stop err // err_inpnn // "total_charge_error_threshold value must be a number"
                            else
                                print *, err, err_inpnn, "total_charge_error_threshold key needs a single argument"; stop
                            end if

                        case ('kalman_damp_short')
                            if (rinpparam%kalman_dampe /= default_real) stop err // err_inpnn // 'Multiple use of the kalman_damp_short key'
                            if (nwords == 2) then
                                read(words(2),*, iostat=ios) rinpparam%kalman_dampe
                                if (ios /= 0) stop err // err_inpnn // "kalman_damp_short value must be a number"
                            else
                                print *, err, err_inpnn, "kalman_damp_short key needs a single argument"; stop
                            end if

                        case ('kalman_damp_force')
                            if (rinpparam%kalman_dampf /= default_real) stop err // err_inpnn // 'Multiple use of thekalman_damp_force  key'
                            if (nwords == 2) then
                                read(words(2),*, iostat=ios) rinpparam%kalman_dampf
                                if (ios /= 0) stop err // err_inpnn // "kalman_damp_force value must be a number"
                            else
                                print *, err, err_inpnn, "kalman_damp_force key needs a single argument"; stop
                            end if

                        case ('kalman_damp_charge')
                            if (rinpparam%kalman_dampq /= default_real) stop err // err_inpnn // 'Multiple use of the kalman_damp_charge key'
                            if (nwords == 2) then
                                read(words(2),*, iostat=ios) rinpparam%kalman_dampq
                                if (ios /= 0) stop err // err_inpnn // "kalman_damp_charge value must be a number"
                            else
                                print *, err, err_inpnn, "kalman_damp_charge key needs a single argument"; stop
                            end if

                        case ('kalman_lambda_short')
                            if (rinpparam%kalmanlambda_local /= default_real) stop err // err_inpnn // 'Multiple use of the kalman_lambda_short key'
                            if (nwords == 2) then
                                read(words(2),*, iostat=ios) rinpparam%kalmanlambda_local
                                if (ios /= 0) stop err // err_inpnn // "kalman_lambda_short value must be a number"
                            else
                                print *, err, err_inpnn, "kalman_lambda_short key needs a single argument"; stop
                            end if

                        case ('kalman_lambda_charge')
                            if (rinpparam%kalmanlambdae_local /= default_real) stop err // err_inpnn // 'Multiple use of the kalman_lambda_charge key'
                            if (nwords == 2) then
                                read(words(2),*, iostat=ios) rinpparam%kalmanlambdae_local
                                if (ios /= 0) stop err // err_inpnn // "kalman_lambda_charge value must be a number"
                            else
                                print *, err, err_inpnn, "kalman_lambda_charge key needs a single argument"; stop
                            end if

                        case ('kalman_lambda_charge_constraint')
                            if (rinpparam%kalmanlambdac /= default_real) stop err // err_inpnn // 'Multiple use of the kalman_lambda_charge_constraint key'
                            if (nwords == 2) then
                                read(words(2),*, iostat=ios) rinpparam%kalmanlambdac
                                if (ios /= 0) stop err // err_inpnn // "kalman_lambda_charge_constraint value must be a number"
                            else
                                print *, err, err_inpnn, "kalman_lambda_charge_constraint key needs a single argument"; stop
                            end if

                        case ('kalman_nue_short')
                            if (rinpparam%kalmannue /= default_real) stop err // err_inpnn // 'Multiple use of the kalman_nue_short key'
                            if (nwords == 2) then
                                read(words(2),*, iostat=ios) rinpparam%kalmannue
                                if (ios /= 0) stop err // err_inpnn // "kalman_nue_short value must be a number"
                            else
                                print *, err, err_inpnn, "kalman_nue_short key needs a single argument"; stop
                            end if

                        case ('kalman_nue_charge')
                            if (rinpparam%kalmannuee /= default_real) stop err // err_inpnn // 'Multiple use of the kalman_nue_charge key'
                            if (nwords == 2) then
                                read(words(2),*, iostat=ios) rinpparam%kalmannuee
                                if (ios /= 0) stop err // err_inpnn // "kalman_nue_charge value must be a number"
                            else
                                print *, err, err_inpnn, "kalman_nue_charge key needs a single argument"; stop
                            end if

                        case ('kalman_nue_charge_constraint')
                            if (rinpparam%kalmannuec /= default_real) stop err // err_inpnn // 'Multiple use of the kalman_nue_charge_constraint key'
                            if (nwords == 2) then
                                read(words(2),*, iostat=ios) rinpparam%kalmannuec
                                if (ios /= 0) stop err // err_inpnn // "kalman_nue_charge_constraint value must be a number"
                            else
                                print *, err, err_inpnn, "kalman_nue_charge_constraint key needs a single argument"; stop
                            end if

                        case ('use_noisematrix')
                            if (rinpparam%lusenoisematrix /= default_bool) stop err // err_inpnn // 'Multiple use of the use_noisematrix key'
                            if (nwords == 1) then
                                rinpparam%lusenoisematrix = .true.
                            else
                                print *, err, err_inpnn, "use_noisematrix key needs no argument(s)"; stop
                            end if

                        case ('kalman_q0')
                            if (rinpparam%kalman_q0 /= default_real) stop err // err_inpnn // 'Multiple use of the kalman_q0 key'
                            if (nwords == 2) then
                                read(words(2),*, iostat=ios) rinpparam%kalman_q0
                                if (ios /= 0) stop err // err_inpnn // "kalman_q0 value must be a number"
                            else
                                print *, err, err_inpnn, "kalman_q0 key needs a single argument"; stop
                            end if

                        case ('kalman_qtau')
                            if (rinpparam%kalman_qtau /= default_real) stop err // err_inpnn // 'Multiple use of the kalman_qtau key'
                            if (nwords == 2) then
                                read(words(2),*, iostat=ios) rinpparam%kalman_qtau
                                if (ios /= 0) stop err // err_inpnn // "kalman_qtau value must be a number"
                            else
                                print *, err, err_inpnn, "kalman_qtau key needs a single argument"; stop
                            end if

                       case ('kalman_qmin')
                            if (rinpparam%kalman_qmin /= default_real) stop err // err_inpnn // 'Multiple use of the kalman_qmin key'
                            if (nwords == 2) then
                                read(words(2),*, iostat=ios) rinpparam%kalman_qmin
                                if (ios /= 0) stop err // err_inpnn // "kalman_qmin value must be a number"
                            else
                                print *, err, err_inpnn, "kalman_qmin key needs a single argument"; stop
                            end if

                        case ('kalman_epsilon')
                            if (rinpparam%kalman_epsilon /= default_real) stop err // err_inpnn // 'Multiple use of the kalman_epsilon key'
                            if (nwords == 2) then
                                read(words(2),*, iostat=ios) rinpparam%kalman_epsilon
                                if (ios /= 0) stop err // err_inpnn // "kalman_epsilon value must be a number"
                            else
                                print *, err, err_inpnn, "kalman_epsilon key needs a single argument"; stop
                            end if

                        case ('steepest_descent_step_energy_short')
                            if (rinpparam%steepeststepe /= default_real) stop err // err_inpnn // 'Multiple use of the steepest_descent_step_energy_short key'
                            if (nwords == 2) then
                                read(words(2),*, iostat=ios) rinpparam%steepeststepe
                                if (ios /= 0) stop err // err_inpnn // "steepest_descent_step_energy_short value must be a number"
                            else
                                print *, err, err_inpnn, "steepest_descent_step_energy_short key needs a single argument"; stop
                            end if

                       case ('steepest_descent_step_force_short')
                            if (rinpparam%steepeststepf /= default_real) stop err // err_inpnn // 'Multiple use of the steepest_descent_step_force_short key'
                            if (nwords == 2) then
                                read(words(2),*, iostat=ios) rinpparam%steepeststepf
                                if (ios /= 0) stop err // err_inpnn // "steepest_descent_step_force_short value must be a number"
                            else
                                print *, err, err_inpnn, "steepest_descent_step_force_short key needs a single argument"; stop
                            end if

                        case ('steepest_descent_step_charge')
                            if (rinpparam%steepeststepq /= default_real) stop err // err_inpnn // 'Multiple use of the steepest_descent_step_charge key'
                            if (nwords == 2) then
                                read(words(2),*, iostat=ios) rinpparam%steepeststepq
                                if (ios /= 0) stop err // err_inpnn // "steepest_descent_step_charge value must be a number"
                            else
                                print *, err, err_inpnn, "steepest_descent_step_charge key needs a single argument"; stop
                            end if

                        case ('force_update_scaling')
                            if (rinpparam%scalefactorf /= default_real) stop err // err_inpnn // 'Multiple use of the force_update_scaling key'
                            if (nwords == 2) then
                                read(words(2),*, iostat=ios) rinpparam%scalefactorf
                                if (ios /= 0) stop err // err_inpnn // "force_update_scaling value must be a number"
                            else
                                print *, err, err_inpnn, "force_update_scaling key needs a single argument"; stop
                            end if

                        case ('charge_update_scaling')
                            if (rinpparam%scalefactorq /= default_real) stop err // err_inpnn // 'Multiple use of the charge_update_scaling key'
                            if (nwords == 2) then
                                read(words(2),*, iostat=ios) rinpparam%scalefactorq
                                if (ios /= 0) stop err // err_inpnn // "charge_update_scaling value must be a number"
                            else
                                print *, err, err_inpnn, "charge_update_scaling key needs a single argument"; stop
                            end if

                        case ('random_order_training')
                            print *, err, err_inpnn, "random_order_training key is obsolete, please use mix_all_points instead"; stop

                        case ('scale_symmetry_functions')
                            if (rinpparam%lscalesym /= default_bool) stop err // err_inpnn // 'Multiple use of the scale_symmetry_functions key'
                            if (nwords == 1) then
                                rinpparam%lscalesym = .true.
                            else
                                print *, err, err_inpnn, "scale_symmetry_functions key needs no argument(s)"; stop
                            end if

                        case ('center_symmetry_functions')
                            if (rinpparam%lcentersym /= default_bool) stop err // err_inpnn // 'Multiple use of the center_symmetry_functions key'
                            if (nwords == 1) then
                                rinpparam%lcentersym = .true.
                            else
                                print *, err, err_inpnn, "center_symmetry_functions key needs no argument(s)"; stop
                            end if

                      case ('use_old_weights_short')
                            if (rinpparam%luseoldweightsshort /= default_bool) stop err // err_inpnn // 'Multiple use of the use_old_weights_short key'
                            if (nwords == 1) then
                                rinpparam%luseoldweightsshort = .true.
                            else
                                print *, err, err_inpnn, "use_old_weights_short key needs no argument(s)"; stop
                            end if

                        case ('use_old_weights_charge')
                            if (rinpparam%luseoldweightscharge /= default_bool) stop err // err_inpnn // 'Multiple use of the use_old_weights_charge key'
                            if (nwords == 1) then
                                rinpparam%luseoldweightscharge = .true.
                            else
                                print *, err, err_inpnn, "use_old_weights_charge key needs no argument(s)"; stop
                            end if

                        case ('save_kalman_matrices')
                            if (rinpparam%lsavekalman /= default_bool) stop err // err_inpnn // 'Multiple use of the save_kalman_matrices key'
                            if (nwords == 1) then
                                rinpparam%lsavekalman = .true.
                            else
                                print *, err, err_inpnn, "save_kalman_matrices key needs no argument(s)"; stop
                            end if

                        case ('read_kalman_matrices')
                            if (rinpparam%lrestkalman /= default_bool) stop err // err_inpnn // 'Multiple use of the read_kalman_matrices key'
                            if (nwords == 1) then
                                rinpparam%lrestkalman = .true.
                            else
                                print *, err, err_inpnn, "read_kalman_matrices key needs no argument(s)"; stop
                            end if

                        case ('update_single_element')
                            if (rinpparam%lupdatebyelement /= default_bool) stop err // err_inpnn // 'Multiple use of the update_single_element key'
                            if (nwords == 2) then
                                rinpparam%lupdatebyelement = .true.
                                read(words(2),'(i1000)', iostat=ios) rinpparam%elemupdate
                                if (ios /= 0) stop err // err_inpnn // "update_single_element value must be integer"
                            else
                                print *, err, err_inpnn, "update_single_element key needs a single argument"; stop
                            end if

                        case ('update_worst_short_energies')
                            if (rinpparam%luseworste /= default_bool) stop err // err_inpnn // 'Multiple use of the update_worst_short_energies key'
                            if (nwords == 2) then
                                rinpparam%luseworste = .true.
                                read(words(2),*, iostat=ios) rinpparam%worste
                                if (ios /= 0) stop err // err_inpnn // "update_worst_short_energies value must be a number"
                            else
                                print *, err, err_inpnn, "update_worst_short_energies key needs a single argument"; stop
                            end if

                        case ('update_worst_short_forces')
                            if (rinpparam%luseworstf /= default_bool) stop err // err_inpnn // 'Multiple use of the update_worst_short_forces key'
                            if (nwords == 2) then
                                rinpparam%luseworstf = .true.
                                read(words(2),*, iostat=ios) rinpparam%worstf
                                if (ios /= 0) stop err // err_inpnn // "update_worst_short_forces value must be a number"
                            else
                                print *, err, err_inpnn, "update_worst_short_forces key needs a single argument"; stop
                            end if

                        case ('update_worst_charges')
                            if (rinpparam%luseworstq /= default_bool) stop err // err_inpnn // 'Multiple use of the update_worst_charges key'
                            if (nwords == 2) then
                                rinpparam%luseworstq = .true.
                                read(words(2),*, iostat=ios) rinpparam%worstq
                                if (ios /= 0) stop err // err_inpnn // "update_worst_charges value must be a number"
                            else
                                print *, err, err_inpnn, "update_worst_charges key needs a single argument"; stop
                            end if

                        case ('growth_mode')
                            if (rinpparam%lgrowth /= default_bool) stop err // err_inpnn // 'Multiple use of the growth_mode key'
                            if (nwords == 3) then
                                rinpparam%lgrowth = .true.
                                read(words(2),'(i1000)', iostat=ios) rinpparam%ngrowth
                                if (ios /= 0) stop err // err_inpnn // "growth_mode first argument value must be integer"
                                read(words(3),'(i1000)', iostat=ios) rinpparam%growthstep
                                if (ios /= 0) stop err // err_inpnn // "growth_mode second argument must be integer"
                            else
                                print *, err, err_inpnn, "growth_mode key needs a single argument"; stop
                            end if

                        case ('use_damping')
                            if (rinpparam%ldampw /= default_bool) stop err // err_inpnn // 'Multiple use of the use_damping key'
                            if (nwords == 2) then
                                rinpparam%ldampw = .true.
                                read(words(2),*, iostat=ios) rinpparam%dampw
                                if (ios /= 0) stop err // err_inpnn // "use_damping value must be a number"
                            else
                                print *, err, err_inpnn, "use_damping key needs a single argument"; stop
                            end if

                        case ('fix_weights')
                            if (rinpparam%lfixweights /= default_bool) stop err // err_inpnn // 'Multiple use of the fix_weights key'
                            if (nwords == 1) then
                                rinpparam%lfixweights = .true.
                            else
                                print *, err, err_inpnn, "fix_weights key needs no argument(s)"; stop
                            end if

                        case ('calculate_forces')
                            if (rinpparam%ldoforces /= default_bool) stop err // err_inpnn // 'Multiple use of the calculate_forces key'
                            if (nwords == 1) then
                                rinpparam%ldoforces = .true.
                            else
                                print *, err, err_inpnn, "calculate_forces key needs no argument(s)"; stop
                            end if

                        case ('calculate_hessian')
                            if (rinpparam%ldohessian /= default_bool) stop err // err_inpnn // 'Multiple use of the calculate_hessian key'
                            if (nwords == 1) then
                                rinpparam%ldohessian = .true.
                            else
                                print *, err, err_inpnn, "calculate_hessian key needs no argument(s)"; stop
                            end if

                        case ('calculate_stress')
                            if (rinpparam%ldostress /= default_bool) stop err // err_inpnn // 'Multiple use of the calculate_stress key'
                            if (nwords == 1) then
                                rinpparam%ldostress = .true.
                            else
                                print *, err, err_inpnn, "calculate_stress key needs no argument(s)"; stop
                            end if

                        case ('enforce_max_num_neighbors_atomic')
                            if (rinpparam%max_num_neighbors_atomic_input /= default_int) stop err // err_inpnn // 'Multiple use of the enforce_max_num_neighbors_atomic key'
                            if (nwords == 2) then
                                rinpparam%lenforcemaxnumneighborsatomic = .true.
                                read(words(2),'(i1000)', iostat=ios) rinpparam%max_num_neighbors_atomic_input
                                if (ios /= 0) stop err // err_inpnn // "enforce_max_num_neighbors_atomic value must be integer"
                            else
                                print *, err, err_inpnn, "enforce_max_num_neighbors_atomic key needs a single argument"; stop
                            end if

                        case ('detailed_timing')
                            if (rinpparam%lfinetime /= default_bool) stop err // err_inpnn // 'Multiple use of the detailed_timing key'
                            if (nwords == 1) then
                                rinpparam%lfinetime = .true.
                            else
                                print *, err, err_inpnn, "detailed_timing key needs no argument(s)"; stop
                            end if

                        case ('detailed_timing_epoch')
                            if (rinpparam%lfinetimeepoch /= default_bool) stop err // err_inpnn // 'Multiple use of the detailed_timing_epoch key'
                            if (nwords == 1) then
                                rinpparam%lfinetimeepoch = .true.
                            else
                                print *, err, err_inpnn, "detailed_timing_epoch key needs no argument(s)"; stop
                            end if

                        case ('write_pdb')
                            print *, err, err_inpnn, "write_pdb key is obsolete, please remove it"; stop

                        case ('write_xyz')
                            print *, err, err_inpnn, "write_xyz key is obsolete, please remove it"; stop

                        case ('write_pov')
                            print *, err, err_inpnn, "write_pov key is obsolete, please remove it"; stop

                        case ('write_pwscf')
                            print *, err, err_inpnn, "write_pwscf key is obsolete, please remove it"; stop

                        case ('write_trainpoints')
                            if (rinpparam%lwritetrainpoints /= default_bool) stop err // err_inpnn // 'Multiple use of the write_trainpoints key'
                            if (nwords == 1) then
                                rinpparam%lwritetrainpoints = .true.
                            else
                                print *, err, err_inpnn, "write_trainpoints key needs no argument(s)"; stop
                            end if

                        case ('write_trainforces')
                            if (rinpparam%lwritetrainforces /= default_bool) stop err // err_inpnn // 'Multiple use of the write_trainforces key'
                            if (nwords == 1) then
                                rinpparam%lwritetrainforces = .true.
                            else
                                print *, err, err_inpnn, "write_trainforces key needs no argument(s)"; stop
                            end if

                        case ('write_traincharges')
                            if (rinpparam%lwritetraincharges /= default_bool) stop err // err_inpnn // 'Multiple use of the write_traincharges key'
                            if (nwords == 1) then
                                rinpparam%lwritetraincharges = .true.
                            else
                                print *, err, err_inpnn, "write_traincharges key needs no argument(s)"; stop
                            end if

                        case ('max_force')
                            if (rinpparam%maxforce /= default_real) stop err // err_inpnn // 'Multiple use of the max_force key'
                            if (nwords == 2) then
                                read(words(2),*, iostat=ios) rinpparam%maxforce
                                if (ios /= 0) stop err // err_inpnn // "max_force value must be a number"
                            else
                                print *, err, err_inpnn, "max_force key needs a single argument"; stop
                            end if

                        case ('max_energy')
                            if (rinpparam%maxenergy /= default_real) stop err // err_inpnn // 'Multiple use of the max_energy key'
                            if (nwords == 2) then
                                read(words(2),*, iostat=ios) rinpparam%maxenergy
                                if (ios /= 0) stop err // err_inpnn // "max_energy value must be a number"
                            else
                                print *, err, err_inpnn, "max_energy key needs a single argument"; stop
                            end if

                        case ('nn_type')
                            print *, err, err_inpnn, "nn_type key is obsolete, please use nn_type_short instead"; stop

                        case ('random_number_type')
                            if (rinpparam%nran /= default_int) stop err // err_inpnn // 'Multiple use of the random_number_type key'
                            if (nwords == 2) then
                                read(words(2),'(i1000)', iostat=ios) rinpparam%nran
                                if (ios /= 0) stop err // err_inpnn // "random_number_type value must be integer"
                            else
                                print *, err, err_inpnn, "random_number_type key needs a single argument"; stop
                            end if

                        case ('calculate_final_force')
                            if (rinpparam%lfinalforce /= default_bool) stop err // err_inpnn // 'Multiple use of the calculate_final_force key'
                            if (nwords == 1) then
                                rinpparam%lfinalforce = .true.
                            else
                                print *, err, err_inpnn, "calculate_final_force key needs no argument(s)"; stop
                            end if

                        case ('normalize_nodes')
                            if (rinpparam%lnormnodes /= default_bool) stop err // err_inpnn // 'Multiple use of the normalize_nodes key'
                            if (nwords == 1) then
                                rinpparam%lnormnodes = .true.
                            else
                                print *, err, err_inpnn, "normalize_nodes key needs no argument(s)"; stop
                            end if

                        case ('atom_energy')
                            ! just let it pass

                        case ('weight_constraint')
                            count_wconstraint=count_wconstraint+1
                            ! just let it pass

                        case ('weighte_constraint')
                            count_wconstraint=count_wconstraint+1 ! same as previous?
                            ! just let it pass

                        case ('weights_min')
                            if (rinpparam%weights_min /= default_real) stop err // err_inpnn // 'Multiple use of the weights_min key'
                            if (nwords == 2) then
                                read(words(2),*, iostat=ios) rinpparam%weights_min
                                if (ios /= 0) stop err // err_inpnn // "weights_min value must be a number"
                            else
                                print *, err, err_inpnn, "weights_min key needs a single argument"; stop
                            end if

                        case ('weights_max')
                            if (rinpparam%weights_max /= default_real) stop err // err_inpnn // 'Multiple use of the weights_max key'
                            if (nwords == 2) then
                                read(words(2),*, iostat=ios) rinpparam%weights_max
                                if (ios /= 0) stop err // err_inpnn // "weights_max value must be a number"
                            else
                                print *, err, err_inpnn, "weights_max key needs a single argument"; stop
                            end if

                        case ('separate_bias_ini_short')
                            if (rinpparam%lseparatebiasini /= default_bool) stop err // err_inpnn // 'Multiple use of the separate_bias_ini_short key'
                            if (nwords == 1) then
                                rinpparam%lseparatebiasini = .true.
                            else
                                print *, err, err_inpnn, "separate_bias_ini_short key needs no argument(s)"; stop
                            end if

                        case ('biasweights_min')
                            if (rinpparam%biasweights_min /= default_real) stop err // err_inpnn // 'Multiple use of the biasweights_min key'
                            if (nwords == 2) then
                                read(words(2),*, iostat=ios) rinpparam%biasweights_min
                                if (ios /= 0) stop err // err_inpnn // "biasweights_min value must be a number"
                            else
                                print *, err, err_inpnn, "biasweights_min key needs a single argument"; stop
                            end if

                        case ('biasweights_max')
                            if (rinpparam%biasweights_max /= default_real) stop err // err_inpnn // 'Multiple use of the biasweights_max key'
                            if (nwords == 2) then
                                read(words(2),*, iostat=ios) rinpparam%biasweights_max
                                if (ios /= 0) stop err // err_inpnn // "biasweights_max value must be a number"
                            else
                                print *, err, err_inpnn, "biasweights_max key needs a single argument"; stop
                            end if

                        case ('weightse_min')
                            if (rinpparam%weightse_min /= default_real) stop err // err_inpnn // 'Multiple use of the weightse_min key'
                            if (nwords == 2) then
                                read(words(2),*, iostat=ios) rinpparam%weightse_min
                                if (ios /= 0) stop err // err_inpnn // "weightse_min value must be a number"
                            else
                                print *, err, err_inpnn, "weightse_min key needs a single argument"; stop
                            end if

                        case ('weightse_max')
                            if (rinpparam%weightse_max /= default_real) stop err // err_inpnn // 'Multiple use of the weightse_max key'
                            if (nwords == 2) then
                                read(words(2),*, iostat=ios) rinpparam%weightse_max
                                if (ios /= 0) stop err // err_inpnn // "weightse_max value must be a number"
                            else
                                print *, err, err_inpnn, "weightse_max key needs a single argument"; stop
                            end if

                        case ('use_systematic_weights_short')
                            if (rinpparam%lsysweights /= default_bool) stop err // err_inpnn // 'Multiple use of the use_systematic_weights_short key'
                            if (nwords == 1) then
                                rinpparam%lsysweights = .true.
                            else
                                print *, err, err_inpnn, "use_systematic_weights_short key needs no argument(s)"; stop
                            end if

                        case ('use_systematic_weights_electrostatic')
                            if (rinpparam%lsysweightse /= default_bool) stop err // err_inpnn // 'Multiple use of the use_systematic_weights_electrostatic key'
                            if (nwords == 1) then
                                rinpparam%lsysweightse = .true.
                            else
                                print *, err, err_inpnn, "use_systematic_weights_electrostatic key needs no argument(s)"; stop
                            end if

                        case ('print_sensitivity')
                            if (rinpparam%lsens /= default_bool) stop err // err_inpnn // 'Multiple use of the print_sensitivity key'
                            if (nwords == 1) then
                                rinpparam%lsens = .true.
                            else
                                print *, err, err_inpnn, "print_sensitivity key needs no argument(s)"; stop
                            end if

                        case ('read_unformatted')
                            if (rinpparam%lreadunformatted /= default_bool) stop err // err_inpnn // 'Multiple use of the read_unformatted key'
                            if (nwords == 1) then
                                rinpparam%lreadunformatted = .true.
                            else
                                print *, err, err_inpnn, "read_unformatted key needs no argument(s)"; stop
                            end if

                        case ('write_unformatted')
                            if (rinpparam%lwriteunformatted /= default_bool) stop err // err_inpnn // 'Multiple use of the write_unformatted key'
                            if (nwords == 1) then
                                rinpparam%lwriteunformatted = .true.
                            else
                                print *, err, err_inpnn, "write_unformatted key needs no argument(s)"; stop
                            end if

                        case ('reset_kalman')
                            if (rinpparam%lresetkalman /= default_bool) stop err // err_inpnn // 'Multiple use of the reset_kalman key'
                            if (nwords == 1) then
                                rinpparam%lresetkalman = .true.
                            else
                                print *, err, err_inpnn, "reset_kalman key needs no argument(s)"; stop
                            end if

                        case ('separate_kalman_short')
                            if (rinpparam%lsepkalman /= default_bool) stop err // err_inpnn // 'Multiple use of the separate_kalman_short key'
                            if (nwords == 1) then
                                rinpparam%lsepkalman = .true.
                            else
                                print *, err, err_inpnn, "separate_kalman_short key needs no argument(s)"; stop
                            end if

                        case ('repeated_energy_update')
                            if (rinpparam%lrepeate /= default_bool) stop err // err_inpnn // 'Multiple use of the repeated_energy_update key'
                            if (nwords == 1) then
                                rinpparam%lrepeate = .true.
                            else
                                print *, err, err_inpnn, "repeated_energy_update key needs no argument(s)"; stop
                            end if

                        case ('enforce_totcharge')
                            if (rinpparam%enforcetotcharge /= default_int) stop err // err_inpnn // 'Multiple use of the enforce_totcharge key'
                            if (nwords == 2) then
                                read(words(2),'(i1000)', iostat=ios) rinpparam%enforcetotcharge
                                if (ios /= 0) stop err // err_inpnn // "enforce_totcharge value must be integer"
                            else
                                print *, err, err_inpnn, "enforce_totcharge key needs a single argument"; stop
                            end if

                        case ('shuffle_weights_short_atomic')
                            if (rinpparam%lshuffle_weights_short_atomic /= default_bool) stop err // err_inpnn // 'Multiple use of the shuffle_weights_short_atomic key'
                            if (nwords == 3) then
                                rin pparam%lshuffle_weights_short_atomic = .true.
                                read(words(2),'(i1000)', iostat=ios) rinpparam%nshuffle_weights_short_atomic
                                if (ios /= 0) stop err // err_inpnn // "shuffle_weights_short_atomic first argument value must be integer"
                                read(words(3),*, iostat=ios) rinpparam%shuffle_weights_short_atomic
                                if (ios /= 0) stop err // err_inpnn // "shuffle_weights_short_atomic second argument value must be a number"
                            else
                                print *, err, err_inpnn, "shuffle_weights_short_atomic key needs two arguments"; stop
                            end if

                        case ('check_forces')
                            if (rinpparam%lcheckf /= default_bool) stop err // err_inpnn // 'Multiple use of the check_forces key'
                            if (nwords == 1) then
                                rinpparam%lcheckf = .true.
                            else
                                print *, err, err_inpnn, "check_forces key needs no argument(s)"; stop
                            end if

                        case ('write_fit_statistics')
                            if (rinpparam%lfitstats /= default_bool) stop err // err_inpnn // 'Multiple use of the write_fit_statistics key'
                            if (nwords == 1) then
                                rinpparam%lfitstats = .true.
                            else
                                print *, err, err_inpnn, "write_fit_statistics key needs no argument(s)"; stop
                            end if

                        case ('fixed_short_energy_error_threshold')
                            if (rinpparam%lfixederrore /= default_bool) stop err // err_inpnn // 'Multiple use of the fixed_short_energy_error_threshold key'
                            if (nwords == 2) then
                                rinpparam%lfixederrore = .true.
                                read(words(2),*, iostat=ios) rinpparam%fixederrore
                                if (ios /= 0) stop err // err_inpnn // "fixed_short_energy_error_threshold value must be a number"
                            else
                                print *, err, err_inpnn, "fixed_short_energy_error_threshold key needs a single argument"; stop
                            end if

                        case ('fixed_short_force_error_threshold')
                            if (rinpparam%lfixederrorf /= default_bool) stop err // err_inpnn // 'Multiple use of the fixed_short_force_error_threshold key'
                            if (nwords == 2) then
                                rinpparam%lfixederrorf = .true.
                                read(words(2),*, iostat=ios) rinpparam%fixederrorf
                                if (ios /= 0) stop err // err_inpnn // "fixed_short_force_error_threshold value must be a number"
                            else
                                print *, err, err_inpnn, "fixed_short_force_error_threshold key needs a single argument"; stop
                            end if

                        case ('restrict_weights')
                            if (rinpparam%restrictw /= default_real) stop err // err_inpnn // 'Multiple use of the restrict_weights key'
                            if (nwords == 2) then
                                read(words(2),*, iostat=ios) rinpparam%restrictw
                                if (ios /= 0) stop err // err_inpnn // "restrict_weights value must be a number"
                            else
                                print *, err, err_inpnn, "restrict_weights key needs a single argument"; stop
                            end if

                        case ('screen_electrostatics')
                            if (rinpparam%lscreen /= default_bool) stop err // err_inpnn // 'Multiple use of the screen_electrostatics key'
                            if (nwords == 3) then
                                rinpparam%lscreen = .true.
                                read(words(2),*, iostat=ios) rinpparam%rscreen_onset
                                if (ios /= 0) stop err // err_inpnn // "screen_electrostatics first argument value must be a number"
                                read(words(3),*, iostat=ios) rinpparam%rscreen_cut
                                if (ios /= 0) stop err // err_inpnn // "screen_electrostatics second argument value must be a number"
                            else
                                print *, err, err_inpnn, " key needs two arguments"; stop
                            end if

                        case ('silent_mode')
                            if (rinpparam%lsilent /= default_bool) stop err // err_inpnn // 'Multiple use of the silent_mode key'
                            if (nwords == 1) then
                                rinpparam%lsilent = .true.
                            else
                                print *, err, err_inpnn, "silent_mode key needs no argument(s)"; stop
                            end if

                        case ('prepare_md')
                            if (rinpparam%lpreparemd /= default_bool) stop err // err_inpnn // 'Multiple use of the prepare_md key'
                            if (nwords == 1) then
                                rinpparam%lpreparemd = .true.
                            else
                                print *, err, err_inpnn, "prepare_md key needs no argument(s)"; stop
                            end if

                        case ('fitting_unit')
                            if (rinpparam%fitting_unit /= default_int) stop err // err_inpnn // 'Multiple use of the fitting_unit'
                            if (nwords == 2) then
                                call lower_case(words(2))

                                select case (words(2))

                                    case ('eV')
                                        fitting_unit = 1

                                    case ('Ha')
                                        fitting_unit = 2

                                    case default
                                        print *, err, err_inpnn, "Error when reading fitting_unit: unknown energy unit specified ", words(2); stop

                                end select

                            else
                                print *, err, err_inpnn, "fitting_unit key needs a single"; stop
                            end if

                        case ('global_activation_short')
                            ! just let it pass

                        case ('global_activation_electrostatic')
                            ! just let it pass

                        case ('global_activation_pair')
                            print *, err, err_inpnn, "global_activation_pair key not supported, Pair NN not implemented"; stop

                        case ('element_hidden_layers_short')
                            ! just let it pass

                        case ('element_hidden_layers_electrostatic')
                            ! just let it pass

                        case ('element_hidden_layers_pair')
                            print *, err, err_inpnn, "element_hidden_layers_pair key not supported, Pair NN not implemented"; stop

                        case ('element_nodes_short')
                            ! just let it pass

                        case ('element_nodes_electrostatic')
                            ! just let it pass

                        case ('element_nodes_pair')
                            print *, err, err_inpnn, "element_nodes_pair key not supported, Pair NN not implemented"; stop

                        case ('element_activation_short')
                            ! just let it pass

                        case ('element_activation_electrostatic')
                            ! just let it pass

                        case ('element_activation_pair')
                            print *, err, err_inpnn, "element_activation_pair key not supported, Pair NN not implemented"; stop

                        case ('joint_energy_force_update')
                            if (rinpparam%ljointefupdate /= default_bool) stop err // err_inpnn // 'Multiple use of the joint_energy_force_update key'
                            if (nwords == 1) then
                                rinpparam%ljointefupdate = .true.
                            else
                                print *, err, err_inpnn, "joint_energy_force_update key needs no argument(s)"; stop
                            end if

                        case ('use_fixed_charges')
                            print *, err, err_inpnn, "use_fixed_charges key is obsolete, please use electrostatic_type 3 instead"; stop

                        case ('use_omp_mkl')
                            if (rinpparam%lompmkl /= default_bool) stop err // err_inpnn // 'Multiple use of the use_omp_mkl key'
                            if (nwords == 1) then
                                rinpparam%lompmkl = .true.
                            else
                                print *, err, err_inpnn, "use_omp_mkl key needs no argument(s)"; stop
                            end if

                        case ('nguyen_widrow_weights_short')
                            if (rinpparam%lnwweights /= default_bool) stop err // err_inpnn // 'Multiple use of the nguyen_widrow_weights_short key'
                            if (nwords == 1) then
                                rinpparam%lnwweights = .true.
                            else
                                print *, err, err_inpnn, "nguyen_widrow_weights_short key needs no argument(s)"; stop
                            end if

                        case ('nguyen_widrow_weights_ewald')
                            if (rinpparam%lnwweightse /= default_bool) stop err // err_inpnn // 'Multiple use of the nguyen_widrow_weights_ewald key'
                            if (nwords == 1) then
                                rinpparam%lnwweightse = .true.
                            else
                                print *, err, err_inpnn, "nguyen_widrow_weights_ewald key needs no argument(s)"; stop
                            end if

                        case ('print_date_and_time')
                            if (rinpparam%lprintdateandtime /= default_bool) stop err // err_inpnn // 'Multiple use of the print_date_and_time key'
                            if (nwords == 1) then
                                rinpparam%lprintdateandtime = .true.
                            else
                                print *, err, err_inpnn, "print_date_and_time key needs no argument(s)"; stop
                            end if

                        case ('enable_on_the_fly_input')
                            if (rinpparam%lenableontheflyinput /= default_bool) stop err // err_inpnn // 'Multiple use of the enable_on_the_fly_input key'
                            if (nwords == 1) then
                                rinpparam%lenableontheflyinput = .true.
                            else
                                print *, err, err_inpnn, "enable_on_the_fly_input key needs no argument(s)"; stop
                            end if

                        case ('element_decoupled_kalman')
                            if (rinpparam%luseedkalman /= default_bool) stop err // err_inpnn // 'Multiple use of the element_decoupled_kalman key'
                            if (nwords == 1) then
                                rinpparam%luseedkalman = .true.
                            else
                                print *, err, err_inpnn, "element_decoupled_kalman key needs no argument(s)"; stop
                            end if

                        case ('element_decoupled_forces_v2')
                            if (rinpparam%ledforcesv2 /= default_bool) stop err // err_inpnn // 'Multiple use of the element_decoupled_forces_v2 key'
                            if (nwords == 1) then
                                rinpparam%ledforcesv2 = .true.
                            else
                                print *, err, err_inpnn, "element_decoupled_forces_v2 key needs no argument(s)"; stop
                            end if

                        case ('analyze_composition')
                            if (rinpparam%lanalyzecomposition /= default_bool) stop err // err_inpnn // 'Multiple use of the analyze_composition key'
                            if (nwords == 1) then
                                rinpparam%lanalyzecomposition = .true.
                            else
                                print *, err, err_inpnn, "analyze_composition key needs no argument(s)"; stop
                            end if

                        case ('fixed_charge')
                            ! just let it pass

                        case ('print_all_short_weights')
                            ! write(pstring(1:1),'(a1)')'1'
                            ! just let it pass

                        case ('print_all_electrostatic_weights')
                            ! write(pstring(2:2),'(a1)')'1'
                            ! just let it pass

                        case ('print_all_deshortdw')
                            ! write(pstring(3:3),'(a1)')'1'
                            ! just let it pass

                        case ('print_all_dfshortdw')
                            ! write(pstring(4:4),'(a1)')'1'
                            ! just let it pass my precious

                        ! check for keyword which is not related to RuNNer keywords
                        case default
                            if (trim(words(1)) /= '' .and. words(1)(1:1) /= '#') & ! check for empty and comment lines
                                print *, err, err_inpnn, 'The keyword ', trim(words(1)),' in line ', line, ' was not recognized, check the spelling or look at the manual'; stop

                    end select

                else
                    write(*,*) err // err_inpnn // 'iostat = ', ios
                    stop
                end if
            end do

            close(inpnn_unit)
            ! end of readout according to readkeywords.f90

            ! further readout according to readinput.f90
            if (lshort .and. (nn_type_short == 1)) then
                do general_counter=1,rinpparam%nelem
                    rinpparam%nodes_short_atomic(rinpparam%maxnum_layers_short_atomic,general_counter)=1
                    if(rinpparam%lelec.and.(rinpparam%nn_type_elec.eq.2))then
                        rinpparam%nodes_short_atomic(rinpparam%maxnum_layers_short_atomic,general_counter) = rinpparam%nodes_short_atomic(rinpparam%maxnum_layers_short_atomic,general_counter)+1
                    endif
                enddo
            endif
            if(rinpparam%lshort.and.(rinpparam%nn_type_short.eq.2))then
                do general_counter=1,rinpparam%npairs
                    rinpparam%nodes_short_pair(rinpparam%maxnum_layers_short_pair,general_counter)=1
                enddo
            endif
            if(rinpparam%lelec.and.(rinpparam%nn_type_elec.eq.1))then
                do general_counter=1,rinpparam%nelem
                    rinpparam%nodes_elec(rinpparam%maxnum_layers_elec,general_counter)=1
                enddo
            endif

            call open_for_read(inpnn_unit, filename_inpnn); ios = 0

            do while (ios == 0)
                read(inpnn_unit, '(A)', iostat=ios) buffer
                if (ios == 0) then
                    call split_string(buffer, words, nwords)

                    select case (words(1))

                        case ('global_activation_short')
                            if (any(actfunc_short_atomic /= default_string)) stop err // err_inpnn // 'Multiple use of the global_activation_short key'
                            if (nwords == maxnum_layers_short_atomic+1) then
                                do general_counter_1 = 1,maxnum_layers_short_atomic
                                    do general_counter_3 = 1,nelem
                                        do general_counter_2 = 1,nodes_short_atomic
                                            read(words(general_counter_1+1),'(A)', iostat=ios) actfunc_short_atomic_dummy(general_counter_1)
                                            actfunc_short_atomic(general_counter_2, general_counter_1, general_counter_3) = actfunc_short_atomic_dummy(general_counter_1)
                                            !actfunc_short_atomic(general_counter_2, general_counter_1, general_counter_3) = words(general_counter_1+1)
                                        end do
                                        if(nodes_short_atomic(general_counter_1, general_counter_3) .lt. maxnodes_short_atomic)then
                                            do general_counter_2 = nodes_short_atomic(general_counter_1, general_counter_3)+1,maxnodes_short_atomic
                                                actfunc_short_atomic(general_counter_2, general_counter_1, general_counter_3) = ' '
                                            enddo
                                    end do
                                end do
                            else
                                print *, err, err_inpnn, "global_activation_short key needs ", maxnum_layers_short_atomic, " arguments according to maxnum_layers_short_atomic value"; stop
                            end if

                        case ('global_activation_electrostatic')
                            if (any(actfunc_elec /= default_string)) stop err // err_inpnn // 'Multiple use of the global_activation_electrostatic key'
                            if (nwords == maxnum_layers_elec+1) then
                                do general_counter_1 = 1,maxnum_layers_elec
                                    do general_counter_3 = 1,nelem
                                        do general_counter_2 = 1,nodes_elec
                                            read(words(general_counter_1+1),'(A)', iostat=ios) actfunc_elec_dummy(general_counter_1)
                                            actfunc_elec(general_counter_2, general_counter_1, general_counter_3) = actfunc_elec_dummy(general_counter_1)
                                            !actfunc_elec(general_counter_2, general_counter_1, general_counter_3) = words(general_counter_1+1)
                                        end do
                                        if (nodes_elec(general_counter_1, general_counter_3) .lt. maxnodes_elec) then
                                            do general_counter_2 = nodes_elec(general_counter_1, general_counter_3)+1,maxnodes_elec
                                                actfunc_elec(general_counter_2, general_counter_1, general_counter_3) = ' '
                                            enddo
                                    end do
                                end do
                            else
                                print *, err, err_inpnn, "global_activation_electrostatic key needs ", maxnum_layers_elec, " arguments according to maxnum_layers_elec value"; stop
                            end if

                        case ('global_activation_pair')
                            print *, err, err_inpnn, "global_activation_pair key not supported, Pair NN not implemented"; stop

                        case default
                            ! Smeagol shows them secret ways that nobody else could find

                    end select

                else
                    write(*,*) err // err_inpnn // 'iostat = ', ios
                    stop
                end if
            end do

            close(inpnn_unit)

            do i=1,nelem
                call nuccharge(element(i),nucelem(i))
            enddo
            call sortelements()

            if(lelec.and.(nn_type_elec.eq.3))then
                call open_for_read(inpnn_unit, filename_inpnn); ios = 0

                do while (ios == 0)
                    read(inpnn_unit, '(A)', iostat=ios) buffer
                    if (ios == 0) then

                        call split_string(buffer, words, nwords)

                        select case (words(1))

                            case ('fixed_charge')
                                !if ((rinpparam%elementtemp /= default_string) .and. (rinpparam%chargetemp /= default_real)) stop err // err_inpnn // 'Multiple use of the fixed_charge key'
                                if (nwords == 3) then
                                    read(words(2),'(A)', iostat=ios) elementtemp
                                    read(words(3),*, iostat=ios) chargetemp
                                    if (ios /= 0) stop err // err_inpnn // "fixed_charge second argument value must be a number"
                                    call nuccharge(elementtemp,ztemp)
                                    fixedcharge(elementindex(ztemp))=chargetemp
                                else
                                    print *, err, err_inpnn, "fixed_charge key needs a single argument"; stop
                                end if

                            case default
                                ! just let it pass

                        end select

                    else
                        write(*,*) err // err_inpnn // 'iostat = ', ios
                        stop
                    end if
                end do

                close(inpnn_unit)

                do general_counter_1 = 1,nelem
                    if (fixedcharge(general_counter_1) .gt. 10.0d0) then
                        print *, err, err_inpnn, "Error when reading fixed_charge: No fixed charge specified for element ",element(general_counter_1); stop
                    endif
                enddo
            endif

            call open_for_read(inpnn_unit, filename_inpnn); ios = 0

            do while (ios == 0)
                read(inpnn_unit, '(A)', iostat=ios) buffer
                if (ios == 0) then
                    call split_string(buffer, words, nwords)

                    select case (words(1))

                        case ('element_hidden_layers_short')
                            !if (rinpparam%elementtemp /= default_string) stop err // err_inpnn // 'Multiple use of the element_hidden_layers_short key'
                            if (nwords == 3) then
                                read(words(2),'(A)', iostat=ios) rinpparam%elementtemp
                                call checkelement(elementtemp)
                                call nuccharge(elementtemp,ztemp)
                                read(words(3),'(i1000)', iostat=ios) rinpparam%num_layers_short_atomic(elementindex(ztemp))
                                if (ios /= 0) stop err // err_inpnn // "element_hidden_layers_short second argument value for element ", ztemp," must be integer"
                                num_layers_short_atomic(elementindex(ztemp)) = num_layers_short_atomic(elementindex(ztemp)) + 1
                                if (num_layers_short_atomic(elementindex(ztemp)) .gt. maxnum_layers_short_atomic) then
                                    print *, err, err_inpnn, "Error when reading element_hidden_layers_short: element ", ztemp, " has too many hidden layers"; stop
                                endif
                                nodes_short_atomic(num_layers_short_atomic(elementindex(ztemp)),elementindex(ztemp)) = 1
                                do general_counter_1=2,maxnodes_short_atomic
                                    actfunc_short_atomic(general_counter_1,num_layers_short_atomic(elementindex(ztemp)),elementindex(ztemp)) = ' '
                                enddo
                            else
                                print *, err, err_inpnn, "element_hidden_layers_short key for element ", ztemp, " needs two arguments"; stop
                            end if

                        case ('element_hidden_layers_electrostatic')
                            !if (rinpparam%elementtemp /= default_string) stop err // err_inpnn // 'Multiple use of the element_hidden_layers_short key'
                            if (nwords == 3) then
                                read(words(2),'(A)', iostat=ios) rinpparam%elementtemp
                                call checkelement(elementtemp)
                                call nuccharge(elementtemp,ztemp)
                                read(words(3),'(i1000)', iostat=ios) rinpparam%num_layers_elec(elementindex(ztemp))
                                if (ios /= 0) stop err // err_inpnn // "element_hidden_layers_electrostatic second argument value for element ", ztemp," must be integer"
                                num_layers_elec(elementindex(ztemp)) = num_layers_elec(elementindex(ztemp)) + 1
                                if (num_layers_elec(elementindex(ztemp)) .gt. maxnum_layers_elec) then
                                    print *, err, err_inpnn, "Error when reading element_hidden_layers_electrostatic: element ", ztemp, " has too many hidden layers"; stop
                                endif
                                nodes_elec(num_layers_elec(elementindex(ztemp)),elementindex(ztemp)) = 1
                                do general_counter_1=2,maxnodes_elec
                                    actfunc_elec(general_counter_1,num_layers_elec(elementindex(ztemp)),elementindex(ztemp)) = ' '
                                enddo
                            else
                                print *, err, err_inpnn, "element_hidden_layers_electrostatic key for element ", ztemp, " needs two arguments"; stop
                            end if

                        case ('element_hidden_layers_pair')
                            print *, err, err_inpnn, "element_hidden_layers_pair key not supported, Pair NN not implemented"; stop

                        case default
                            ! just let it pass

                    end select

                else
                    write(*,*) err // err_inpnn // 'iostat = ', ios
                    stop
                end if
            end do

            close(inpnn_unit)

            call open_for_read(inpnn_unit, filename_inpnn); ios = 0

            do while (ios == 0)
                read(inpnn_unit, '(A)', iostat=ios) buffer
                if (ios == 0) then
                    call split_string(buffer, words, nwords)

                    select case (words(1))

                        case ('element_nodes_short')
                            if (nwords == 4) then
                                read(words(2),'(A)', iostat=ios) elementtemp
                                call checkelement(elementtemp)
                                call nuccharge(elementtemp,ztemp)
                                read(words(3),'(i1000)', iostat=ios) layer
                                if (ios /= 0) stop err // err_inpnn // "element_nodes_short second argument value for element ", element(elementindex(ztemp)), " must be integer"
                                read(words(4),'(i1000)', iostat=ios) node
                                if (ios /= 0) stop err // err_inpnn // "element_nodes_short third argument value for element ", element(elementindex(ztemp)), " must be integer"
                                if (layer .eq. num_layers_short_atomic(elementindex(ztemp))) then
                                    print *, err, err_inpnn, "Error when reading element_nodes_short: do not modifiy the number of output nodes for element ", element(elementindex(ztemp)); stop
                                endif
                                if (node .gt. maxnodes_short_atomic) then
                                   print *, err, err_inpnn, "Error when reading element_nodes_short: too many nodes requested for element ", element(elementindex(ztemp)); stop
                                end if
                                nodes_short_atomic(layer,elementindex(ztemp)) = node
                                do general_counter_1 = nodes_short_atomic(layer,elementindex(ztemp))+1,maxnodes_short_atomic
                                    actfunc_short_atomic(general_counter_1,layer,elementindex(ztemp)) = ' '
                                end do
                            else
                                print *, err, err_inpnn, "element_nodes_short key for element ", element(elementindex(ztemp)), " needs three arguments"; stop
                            end if

                        case ('element_nodes_electrostatic')
                            if (nwords == 4) then
                                read(words(2),'(A)', iostat=ios) elementtemp
                                call checkelement(elementtemp)
                                call nuccharge(elementtemp,ztemp)
                                read(words(3),'(i1000)', iostat=ios) layer
                                if (ios /= 0) stop err // err_inpnn // "element_nodes_electrostatic second argument value for element ", element(elementindex(ztemp)), " must be integer"
                                read(words(4),'(i1000)', iostat=ios) node
                                if (ios /= 0) stop err // err_inpnn // "element_nodes_electrostatic third argument value for element ", element(elementindex(ztemp)), " must be integer"
                                if (layer .eq. num_layers_elec(elementindex(ztemp))) then
                                    print *, err, err_inpnn, "Error when reading element_nodes_electrostatic: do not modifiy the number of output nodes for element ", element(elementindex(ztemp)); stop
                                endif
                                if (node .gt. maxnodes_elec) then
                                   print *, err, err_inpnn, "Error when reading element_nodes_electrostatic: too many nodes requested for element ", element(elementindex(ztemp)); stop
                                end if
                                nodes_elec(layer,elementindex(ztemp)) = node
                                do general_counter_1 = nodes_elec(layer,elementindex(ztemp))+1,maxnodes_elec
                                    actfunc_elec(general_counter_1,layer,elementindex(ztemp)) = ' '
                                end do
                            else
                                print *, err, err_inpnn, "element_nodes_electrostatic key for element ", element(elementindex(ztemp)), " needs three arguments"; stop
                            end if

                        case ('element_nodes_pair')
                            print *, err, err_inpnn, "element_nodes_pair key not supported, Pair NN not implemented"; stop

                        case default
                            ! just let it pass

                    end select

                else
                    write(*,*) err // err_inpnn // 'iostat = ', ios
                    stop
                end if
            end do

            close(inpnn_unit)

            call open_for_read(inpnn_unit, filename_inpnn); ios = 0

            do while (ios == 0)
                read(inpnn_unit, '(A)', iostat=ios) buffer
                if (ios == 0) then
                    call split_string(buffer, words, nwords)

                    select case (words(1))

                        case ('element_activation_short')
                            if (nwords == 5) then
                                read(words(2),'(A)', iostat=ios) rinpparam%elementtemp
                                call checkelement(elementtemp)
                                call nuccharge(elementtemp,ztemp)
                                read(words(3),'(i1000)', iostat=ios) rinpparam%layer
                                if (ios /= 0) stop err // err_inpnn // "element_activation_short second argument value for element ", element(elementindex(ztemp)), " must be integer"
                                read(words(4),'(i1000)', iostat=ios) rinpparam%node
                                if (ios /= 0) stop err // err_inpnn // "element_activation_short third argument value for element ", element(elementindex(ztemp)), " must be integer"
                                read(words(5),'(A)', iostat=ios) rinpparam%actfunc
                                if (layer .gt. num_layers_short_atomic(elementindex(ztemp))) then
                                    print *, err, err_inpnn, "Error when reading element_activation_short: layer is too large for element ", element(elementindex(ztemp)); stop
                                endif
                                if (node .gt. nodes_short_atomic(layer,elementindex(ztemp))) then
                                    print *, err, err_inpnn, "Error when reading element_activation_short: node is too large for element ", element(elementindex(ztemp)); stop
                                endif
                                actfunc_short_atomic(node,layer,elementindex(ztemp))=actfunc
                            else
                                print *, err, err_inpnn, "element_activation_short key for element ", element(elementindex(ztemp)), " needs four arguments"; stop
                            end if

                        case ('element_activation_electrostatic')
                            if (nwords == 5) then
                                read(words(2),'(A)', iostat=ios) rinpparam%elementtemp
                                call checkelement(elementtemp)
                                call nuccharge(elementtemp,ztemp)
                                read(words(3),'(i1000)', iostat=ios) rinpparam%layer
                                if (ios /= 0) stop err // err_inpnn // "element_activation_electrostatic second argument value for element ", element(elementindex(ztemp)), " must be integer"
                                read(words(4),'(i1000)', iostat=ios) rinpparam%node
                                if (ios /= 0) stop err // err_inpnn // "element_activation_electrostatic third argument value for element ", element(elementindex(ztemp)), " must be integer"
                                read(words(5),'(A)', iostat=ios) rinpparam%actfunc
                                if (layer .gt. num_layers_elec(elementindex(ztemp))) then
                                    print *, err, err_inpnn, "Error when reading element_activation_electrostatic: layer is too large for element ", element(elementindex(ztemp)); stop
                                endif
                                if (node .gt. nodes_elec(layer,elementindex(ztemp))) then
                                    print *, err, err_inpnn, "Error when reading element_activation_electrostatic: node is too large for element ", element(elementindex(ztemp)); stop
                                endif
                                actfunc_elec(node,layer,elementindex(ztemp))=actfunc
                            else
                                print *, err, err_inpnn, "element_activation_electrostatic key for element ", element(elementindex(ztemp)), " needs four arguments"; stop
                            end if

                        case ('element_activation_pair')
                            print *, err, err_inpnn, "element_activation_pair key not supported, Pair NN not implemented"; stop

                        case default
                            ! just let it pass

                    end select

                else
                    write(*,*) err // err_inpnn // 'iostat = ', ios
                    stop
                end if
            end do

            close(inpnn_unit)

            if (lshort .and. (nn_type_short == 1)) then
                sym_short_atomic_count(:)=0
                num_funcvalues_short_atomic(:)=0
            endif
            if (lelec .and. (nn_type_elec == 1)) then
                sym_elec_count(:)=0
                num_funcvalues_elec(:)=0
            endif

            call open_for_read(inpnn_unit, filename_inpnn); ios = 0

            do while (ios == 0)
                read(inpnn_unit, '(A)', iostat=ios) buffer
                if (ios == 0) then
                    call split_string(buffer, words, nwords)

                    select case (words(1))

                        case ('symfunction_short') ! allocation of arrays is done in allocatesymfunctions() from module symfunctions.f90 when called in initnn.f90! -> no additional allocation subroutine needed!!
                            if (lshort .and. (nn_type_short == 1)) then
                                !call allocate_readsymfunctionatomic(maxnum_funcvalues_short_atomic, sym_short_atomic_count, function_type_short_atomic, symelement_short_atomic, &
                                     !funccutoff_short_atomic, eta_short_atomic, zeta_short_atomic, rshift_short_atomic, lambda_short_atomic) ! maybe there is a better way to set the dimensions of variables needed for the readout?
                                !if (nwords == 5) then
                                read(words(2),'(A)', iostat=ios) rinpparam%elementtemp1
                                call checkelement(elementtemp1)
                                call nuccharge(elementtemp1,ztemp1)
                                sym_short_atomic_count(elementindex(ztemp1)) = sym_short_atomic_count(elementindex(ztemp1)) + 1
                                read(words(3),'(i1000)', iostat=ios) rinpparam%function_type_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),elementindex(ztemp1))
                                if (ios /= 0) stop err // err_inpnn // "symfunction_short argument 2 value must be integer"

                                select case(words(3))

                                    case ('1')
                                        if (nwords == 5) then
                                            read(words(4),'(A)') rinpparam%elementtemp2
                                            call checkelement(elementtemp2)
                                            call nuccharge(elementtemp2,ztemp2)
                                            read(words(5),*, iostat=ios) rinpparam%funccutoff_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),elementindex(ztemp1))
                                            if (ios /= 0) stop err // err_inpnn // "symfunction_short type ", words(3), " argument 4 must be a number"
                                            symelement_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),1,elementindex(ztemp1))=ztemp2
                                        else
                                            print *, err, err_inpnn, "symfunction_short type ", words(3), " needs 4 arguments"; stop
                                        end if

                                    case ('2')
                                        if (nwords == 7) then
                                            read(words(4),'(A)') rinpparam%elementtemp2
                                            call checkelement(elementtemp2)
                                            call nuccharge(elementtemp2,ztemp2)
                                            read(words(5),*, iostat=ios) rinpparam%eta_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),elementindex(ztemp1))
                                            if (ios /= 0) stop err // err_inpnn // "symfunction_short type ", words(3), " argument 4 must be a number"
                                            read(words(6),*, iostat=ios) rinpparam%rshift_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),elementindex(ztemp1))
                                            if (ios /= 0) stop err // err_inpnn // "symfunction_short type ", words(3), " argument 5 must be a number"
                                            read(words(7),*, iostat=ios) rinpparam%funccutoff_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),elementindex(ztemp1))
                                            if (ios /= 0) stop err // err_inpnn // "symfunction_short type ", words(3), " argument 6 must be a number"
                                            symelement_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),1,elementindex(ztemp1))=ztemp2
                                        else
                                            print *, err, err_inpnn, "symfunction_short type ", words(3), " needs 6 arguments"; stop
                                        end if

                                    case ('3')
                                        if (nwords == 9) then
                                            read(words(4),'(A)') rinpparam%elementtemp2
                                            read(words(5),'(A)') rinpparam%elementtemp3
                                            call checkelement(elementtemp2)
                                            call checkelement(elementtemp3)
                                            call nuccharge(elementtemp2,ztemp2)
                                            call nuccharge(elementtemp2,ztemp3)
                                            if (ztemp2 .gt. ztemp3) then
                                                itemp = ztemp2
                                                ztemp2 = ztemp3
                                                ztemp3 = itemp
                                            endif
                                            read(words(6),*, iostat=ios) rinpparam%eta_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),elementindex(ztemp1))
                                            if (ios /= 0) stop err // err_inpnn // "symfunction_short type ", words(3), " argument 5 must be a number"
                                            read(words(7),*, iostat=ios) rinpparam%lambda_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),elementindex(ztemp1))
                                            if (ios /= 0) stop err // err_inpnn // "symfunction_short type ", words(3), " argument 6 must be a number"
                                            read(words(8),*, iostat=ios) rinpparam%zeta_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),elementindex(ztemp1))
                                            if (ios /= 0) stop err // err_inpnn // "symfunction_short type ", words(3), " argument 7 must be a number"
                                            read(words(9),*, iostat=ios) rinpparam%funccutoff_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),elementindex(ztemp1))
                                            if (ios /= 0) stop err // err_inpnn // "symfunction_short type ", words(3), " argument 8 must be a number"
                                            symelement_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),1,elementindex(ztemp1))=ztemp2
                                            symelement_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),2,elementindex(ztemp1))=ztemp3
                                        else
                                            print *, err, err_inpnn, "symfunction_short type ", words(3), " needs 8 arguments"; stop
                                        end if

                                    case ('4')
                                        if (nwords == 6) then
                                            read(words(4),'(A)') rinpparam%elementtemp2
                                            call checkelement(elementtemp2)
                                            call nuccharge(elementtemp2,ztemp2)
                                            read(words(5),*, iostat=ios) rinpparam%eta_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),elementindex(ztemp1))
                                            if (ios /= 0) stop err // err_inpnn // "symfunction_short type ", words(3), " argument 4 must be a number"
                                            read(words(6),*, iostat=ios) rinpparam%funccutoff_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),elementindex(ztemp1))
                                            if (ios /= 0) stop err // err_inpnn // "symfunction_short type ", words(3), " argument 5 must be a number"
                                            symelement_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),1,elementindex(ztemp1))=ztemp2
                                        else
                                            print *, err, err_inpnn, "symfunction_short type ", words(3), " needs 5 arguments"; stop
                                        end if

                                    case ('5')
                                        if (nwords == 4) then
                                            read(words(4),*, iostat=ios) rinpparam%eta_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),elementindex(ztemp1))
                                            if (ios /= 0) stop err // err_inpnn // "symfunction_short type ", words(3), " argument 3 must be a number"
                                        else
                                            print *, err, err_inpnn, "symfunction_short type ", words(3), " needs 3 arguments"; stop
                                        end if

                                    case ('6')
                                        if (nwords == 5) then
                                            read(words(4),'(A)') rinpparam%elementtemp2
                                            call checkelement(elementtemp2)
                                            call nuccharge(elementtemp2,ztemp2)
                                            read(words(5),*, iostat=ios) rinpparam%funccutoff_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),elementindex(ztemp1))
                                            if (ios /= 0) stop err // err_inpnn // "symfunction_short type ", words(3), " argument 4 must be a number"
                                            symelement_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),1,elementindex(ztemp1))=ztemp2
                                        else
                                            print *, err, err_inpnn, "symfunction_short type ", words(3), " needs 4 arguments"; stop
                                        end if

                                    case ('8')
                                        if (nwords == 8) then
                                            read(words(4),'(A)') rinpparam%elementtemp2
                                            read(words(5),'(A)') rinpparam%elementtemp3
                                            call checkelement(elementtemp2)
                                            call checkelement(elementtemp3)
                                            call nuccharge(elementtemp2,ztemp2)
                                            call nuccharge(elementtemp2,ztemp3)
                                            if (ztemp2 .gt. ztemp3) then
                                                itemp = ztemp2
                                                ztemp2 = ztemp3
                                                ztemp3 = itemp
                                            endif
                                            read(words(6),*, iostat=ios) rinpparam%eta_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),elementindex(ztemp1))
                                            if (ios /= 0) stop err // err_inpnn // "symfunction_short type ", words(3), " argument 5 must be a number"
                                            read(words(7),*, iostat=ios) rinpparam%rshift_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),elementindex(ztemp1))
                                            if (ios /= 0) stop err // err_inpnn // "symfunction_short type ", words(3), " argument 6 must be a number"
                                            read(words(8),*, iostat=ios) rinpparam%funccutoff_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),elementindex(ztemp1))
                                            if (ios /= 0) stop err // err_inpnn // "symfunction_short type ", words(3), " argument 7 must be a number"
                                            symelement_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),1,elementindex(ztemp1))=ztemp2
                                            symelement_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),2,elementindex(ztemp1))=ztemp3
                                        else
                                            print *, err, err_inpnn, "symfunction_short type ", words(3), " needs 7 arguments"; stop
                                        end if

                                    case ('9')
                                        if (nwords == 9) then
                                            read(words(4),'(A)') rinpparam%elementtemp2
                                            read(words(5),'(A)') rinpparam%elementtemp3
                                            call checkelement(elementtemp2)
                                            call checkelement(elementtemp3)
                                            call nuccharge(elementtemp2,ztemp2)
                                            call nuccharge(elementtemp2,ztemp3)
                                            if (ztemp2 .gt. ztemp3) then
                                                itemp = ztemp2
                                                ztemp2 = ztemp3
                                                ztemp3 = itemp
                                            endif
                                            read(words(6),*, iostat=ios) rinpparam%eta_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),elementindex(ztemp1))
                                            if (ios /= 0) stop err // err_inpnn // "symfunction_short type ", words(3), " argument 5 must be a number"
                                            read(words(7),*, iostat=ios) rinpparam%lambda_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),elementindex(ztemp1))
                                            if (ios /= 0) stop err // err_inpnn // "symfunction_short type ", words(3), " argument 6 must be a number"
                                            read(words(8),*, iostat=ios) rinpparam%zeta_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),elementindex(ztemp1))
                                            if (ios /= 0) stop err // err_inpnn // "symfunction_short type ", words(3), " argument 7 must be a number"
                                            read(words(9),*, iostat=ios) rinpparam%funccutoff_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),elementindex(ztemp1))
                                            if (ios /= 0) stop err // err_inpnn // "symfunction_short type ", words(3), " argument 8 must be a number"
                                            symelement_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),1,elementindex(ztemp1))=ztemp2
                                            symelement_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),2,elementindex(ztemp1))=ztemp3
                                        else
                                            print *, err, err_inpnn, "symfunction_short type ", words(3), " needs 8 arguments"; stop
                                        end if

                                    case default
                                        print *, err, err_inpnn, "Error in symfunction_short key, symfunction type ", words(3), " not implemented"
                                        stop

                                end select

                            end if

                        case ('element_symfunction_short')
                            if (lshort .and. (nn_type_short == 1)) then
                                read(words(2),'(A)', iostat=ios) elementtemp1
                                call checkelement(elementtemp1)
                                call nuccharge(elementtemp1,ztemp1)
                                read(words(3),'(i1000)', iostat=ios) function_type_temp
                                if (ios /= 0) stop err // err_inpnn // "element_symfunction_short argument 2 value must be integer"

                                select case(words(3))

                                    case ('1')
                                        if (nwords == 4) then
                                            read(words(4),*, iostat=ios) funccutoff_temp
                                            if (ios /= 0) stop err // err_inpnn // "element_symfunction_short type ", words(3), " argument 3 must be a number"
                                            do general_counter_1 = 1,nelem
                                                sym_short_atomic_count(elementindex(ztemp1) = sym_short_atomic_count(elementindex(ztemp1) + 1
                                                function_type_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),elementindex(ztemp1)) = function_type_temp
                                                funccutoff_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),elementindex(ztemp1)) = funccutoff_temp
                                                symelement_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),1,elementindex(ztemp1)) = nucelem(general_counter_1)
                                            end do
                                        else
                                            print *, err, err_inpnn, "element_symfunction_short type ", words(3), " needs 3 arguments"; stop
                                        end if

                                    case ('2')
                                        if (nwords == 6) then
                                            read(words(4),*, iostat=ios) eta_temp
                                            if (ios /= 0) stop err // err_inpnn // "element_symfunction_short type ", words(3), " argument 3 must be a number"
                                            read(words(5),*, iostat=ios) rshift_temp
                                            if (ios /= 0) stop err // err_inpnn // "element_symfunction_short type ", words(3), " argument 4 must be a number"
                                            read(words(6),*, iostat=ios) funccutoff_temp
                                            if (ios /= 0) stop err // err_inpnn // "element_symfunction_short type ", words(3), " argument 5 must be a number"
                                            do general_counter_1 = 1,nelem
                                                sym_short_atomic_count(elementindex(ztemp1) = sym_short_atomic_count(elementindex(ztemp1) + 1
                                                function_type_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),elementindex(ztemp1)) = function_type_temp
                                                eta_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),elementindex(ztemp1)) = eta_temp
                                                rshift_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),elementindex(ztemp1)) = rshift_temp
                                                funccutoff_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),elementindex(ztemp1)) = funccutoff_temp
                                                symelement_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),1,elementindex(ztemp1)) = nucelem(general_counter_1)
                                            end do
                                        else
                                            print *, err, err_inpnn, "element_symfunction_short type ", words(3), " needs 5 arguments"; stop
                                        end if

                                    case ('3')
                                        if (nwords == 7) then
                                            read(words(4),*, iostat=ios) eta_temp
                                            if (ios /= 0) stop err // err_inpnn // "element_symfunction_short type ", words(3), " argument 3 must be a number"
                                            read(words(5),*, iostat=ios) lambda_temp
                                            if (ios /= 0) stop err // err_inpnn // "element_symfunction_short type ", words(3), " argument 4 must be a number"
                                            read(words(6),*, iostat=ios) zeta_temp
                                            if (ios /= 0) stop err // err_inpnn // "element_symfunction_short type ", words(3), " argument 5 must be a number"
                                            read(words(7),*, iostat=ios) funccutoff_temp
                                            if (ios /= 0) stop err // err_inpnn // "element_symfunction_short type ", words(3), " argument 6 must be a number"
                                            do general_counter_1 = 1,nelem
                                                sym_short_atomic_count(elementindex(ztemp1) = sym_short_atomic_count(elementindex(ztemp1) + 1
                                                function_type_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),elementindex(ztemp1)) = function_type_temp
                                                eta_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),elementindex(ztemp1)) = eta_temp
                                                lambda_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),elementindex(ztemp1)) = lambda_temp
                                                zeta_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),elementindex(ztemp1)) = zeta_temp
                                                funccutoff_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),elementindex(ztemp1)) = funccutoff_temp
                                                symelement_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),1,elementindex(ztemp1)) = nucelem(general_counter_1)
                                                symelement_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),2,elementindex(ztemp1)) = nucelem(general_counter_1)
                                            end do
                                            do general_counter_1 = 1,nelem
                                                if (nelem .gt. 1) then
                                                    do general_counter_2 = 1,general_counter_1-1
                                                        sym_short_atomic_count(elementindex(ztemp1) = sym_short_atomic_count(elementindex(ztemp1) + 1
                                                        function_type_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),elementindex(ztemp1)) = function_type_temp
                                                        eta_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),elementindex(ztemp1)) = eta_temp
                                                        lambda_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),elementindex(ztemp1)) = lambda_temp
                                                        zeta_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),elementindex(ztemp1)) = zeta_temp
                                                        funccutoff_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),elementindex(ztemp1)) = funccutoff_temp
                                                        symelement_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),1,elementindex(ztemp1)) = nucelem(general_counter_2)
                                                        symelement_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),2,elementindex(ztemp1)) = nucelem(general_counter_1)
                                                    end do
                                                end if
                                            end do
                                        else
                                            print *, err, err_inpnn, "element_symfunction_short type ", words(3), " needs 6 arguments"; stop
                                        end if

                                    case ('4')
                                        if (nwords == 5) then
                                            read(words(4),*, iostat=ios) eta_temp
                                            if (ios /= 0) stop err // err_inpnn // "element_symfunction_short type ", words(3), " argument 3 must be a number"
                                            read(words(5),*, iostat=ios) funccutoff_temp
                                            if (ios /= 0) stop err // err_inpnn // "element_symfunction_short type ", words(3), " argument 4 must be a number"
                                            do general_counter_1 = 1,nelem
                                                sym_short_atomic_count(elementindex(ztemp1) = sym_short_atomic_count(elementindex(ztemp1) + 1
                                                function_type_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),elementindex(ztemp1)) = function_type_temp
                                                eta_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),elementindex(ztemp1)) = eta_temp
                                                funccutoff_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),elementindex(ztemp1)) = funccutoff_temp
                                                symelement_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),1,elementindex(ztemp1)) = nucelem(general_counter_1)
                                            end do
                                        else
                                            print *, err, err_inpnn, "element_symfunction_short type ", words(3), " needs 4 arguments"; stop
                                        end if

                                    case ('5')
                                        if (nwords == 4) then
                                            read(words(4),*, iostat=ios) eta_temp
                                            if (ios /= 0) stop err // err_inpnn // "element_symfunction_short type ", words(3), " argument 3 must be a number"
                                            sym_short_atomic_count(elementindex(ztemp1) = sym_short_atomic_count(elementindex(ztemp1) + 1
                                            function_type_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),elementindex(ztemp1)) = function_type_temp
                                            eta_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),elementindex(ztemp1)) = eta_temp
                                        else
                                            print *, err, err_inpnn, "element_symfunction_short type ", words(3), " needs 3 arguments"; stop
                                        end if

                                    case ('6')
                                        if (nwords == 4) then
                                            read(words(4),*, iostat=ios) funccutoff_temp
                                            if (ios /= 0) stop err // err_inpnn // "element_symfunction_short type ", words(3), " argument 3 must be a number"
                                            sym_short_atomic_count(elementindex(ztemp1) = sym_short_atomic_count(elementindex(ztemp1) + 1
                                            function_type_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),elementindex(ztemp1)) = function_type_temp
                                        else
                                            print *, err, err_inpnn, "element_symfunction_short type ", words(3), " needs 3 arguments"; stop
                                        end if

                                    case ('8')
                                        if (nwords == 6) then
                                            read(words(4),*, iostat=ios) eta_temp
                                            if (ios /= 0) stop err // err_inpnn // "element_symfunction_short type ", words(3), " argument 3 must be a number"
                                            read(words(5),*, iostat=ios) rshift_temp
                                            if (ios /= 0) stop err // err_inpnn // "element_symfunction_short type ", words(3), " argument 4 must be a number"
                                            read(words(6),*, iostat=ios) funccutoff_temp
                                            if (ios /= 0) stop err // err_inpnn // "element_symfunction_short type ", words(3), " argument 5 must be a number"
                                            do general_counter_1 = 1,nelem
                                                sym_short_atomic_count(elementindex(ztemp1) = sym_short_atomic_count(elementindex(ztemp1) + 1
                                                function_type_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),elementindex(ztemp1)) = function_type_temp
                                                eta_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),elementindex(ztemp1)) = eta_temp
                                                rshift_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),elementindex(ztemp1)) = rshift_temp
                                                funccutoff_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),elementindex(ztemp1)) = funccutoff_temp
                                                symelement_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),1,elementindex(ztemp1)) = nucelem(general_counter_1)
                                                symelement_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),2,elementindex(ztemp1)) = nucelem(general_counter_1)
                                            end do
                                            do general_counter_1 = 1,nelem
                                                if (nelem .gt. 1) then
                                                    do general_counter_2 = 1,general_counter_1-1
                                                        sym_short_atomic_count(elementindex(ztemp1) = sym_short_atomic_count(elementindex(ztemp1) + 1
                                                        function_type_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),elementindex(ztemp1)) = function_type_temp
                                                        eta_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),elementindex(ztemp1)) = eta_temp
                                                        rshift_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),elementindex(ztemp1)) = rshift_temp
                                                        funccutoff_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),elementindex(ztemp1)) = funccutoff_temp
                                                        symelement_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),1,elementindex(ztemp1)) = nucelem(general_counter_2)
                                                        symelement_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),2,elementindex(ztemp1)) = nucelem(general_counter_1)
                                                    end do
                                                end if
                                            end do
                                        else
                                            print *, err, err_inpnn, "element_symfunction_short type ", words(3), " needs 5 arguments"; stop
                                        end if

                                    case ('9')
                                        if (nwords == 7) then
                                            read(words(4),*, iostat=ios) eta_temp
                                            if (ios /= 0) stop err // err_inpnn // "element_symfunction_short type ", words(3), " argument 3 must be a number"
                                            read(words(5),*, iostat=ios) lambda_temp
                                            if (ios /= 0) stop err // err_inpnn // "element_symfunction_short type ", words(3), " argument 4 must be a number"
                                            read(words(6),*, iostat=ios) zeta_temp
                                            if (ios /= 0) stop err // err_inpnn // "element_symfunction_short type ", words(3), " argument 5 must be a number"
                                            read(words(7),*, iostat=ios) funccutoff_temp
                                            if (ios /= 0) stop err // err_inpnn // "element_symfunction_short type ", words(3), " argument 6 must be a number"
                                            do general_counter_1 = 1,nelem
                                                sym_short_atomic_count(elementindex(ztemp1) = sym_short_atomic_count(elementindex(ztemp1) + 1
                                                function_type_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),elementindex(ztemp1)) = function_type_temp
                                                eta_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),elementindex(ztemp1)) = eta_temp
                                                lambda_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),elementindex(ztemp1)) = lambda_temp
                                                zeta_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),elementindex(ztemp1)) = zeta_temp
                                                funccutoff_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),elementindex(ztemp1)) = funccutoff_temp
                                                symelement_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),1,elementindex(ztemp1)) = nucelem(general_counter_1)
                                            end do
                                            do general_counter_1 = 1,nelem
                                                if (nelem .gt. 1) then
                                                    do general_counter_2 = 1,general_counter_1-1
                                                        sym_short_atomic_count(elementindex(ztemp1) = sym_short_atomic_count(elementindex(ztemp1) + 1
                                                        function_type_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),elementindex(ztemp1)) = function_type_temp
                                                        eta_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),elementindex(ztemp1)) = eta_temp
                                                        lambda_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),elementindex(ztemp1)) = lambda_temp
                                                        zeta_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),elementindex(ztemp1)) = zeta_temp
                                                        funccutoff_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),elementindex(ztemp1)) = funccutoff_temp
                                                        symelement_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),1,elementindex(ztemp1)) = nucelem(general_counter_2)
                                                        symelement_short_atomic(sym_short_atomic_count(elementindex(ztemp1)),2,elementindex(ztemp1)) = nucelem(general_counter_1)
                                                    end do
                                                end if
                                            end do
                                        else
                                            print *, err, err_inpnn, "element_symfunction_short type ", words(3), " needs 6 arguments"; stop
                                        end if

                                    case default
                                        print *, err, err_inpnn, "Error in element_symfunction_short key, symfunction type ", words(3), " not implemented"
                                        stop

                                end select

                            end if

                        case ('global_symfunction_short')
                            if (lshort .and. (nn_type_short == 1)) then
                                read(words(2),'(i1000)', iostat=ios) function_type_temp
                                if (ios /= 0) stop err // err_inpnn // "global_symfunction_short argument 1 value must be integer"

                                select case(words(2))

                                    case ('1')
                                        if (nwords == 3) then
                                            read(words(3),*, iostat=ios) funccutoff_temp
                                            if (ios /= 0) stop err // err_inpnn // "global_symfunction_short type ", words(2), " argument 2 must be a number"
                                            do general_counter_1 = 1,nelem
                                                do general_counter_2 = 1,nelem
                                                    sym_short_atomic_count(general_counter_1) = sym_short_atomic_count(general_counter_1) + 1
                                                    function_type_short_atomic(sym_short_atomic_count(general_counter_1),general_counter_1) = function_type_temp
                                                    funccutoff_short_atomic(sym_short_atomic_count(general_counter_1),general_counter_1) = funccutoff_temp
                                                    symelement_short_atomic(sym_short_atomic_count(general_counter_1),1,general_counter_1) = nucelem(general_counter_2)
                                                end do
                                            end do
                                        else
                                            print *, err, err_inpnn, "global_symfunction_short type ", words(2), " needs 2 arguments"; stop
                                        end if

                                    case ('2')
                                        if (nwords == 5) then
                                            read(words(3),*, iostat=ios) eta_temp
                                            if (ios /= 0) stop err // err_inpnn // "global_symfunction_short type ", words(2), " argument 2 must be a number"
                                            read(words(4),*, iostat=ios) rshift_temp
                                            if (ios /= 0) stop err // err_inpnn // "global_symfunction_short type ", words(2), " argument 3 must be a number"
                                            read(words(5),*, iostat=ios) funccutoff_temp
                                            if (ios /= 0) stop err // err_inpnn // "global_symfunction_short type ", words(2), " argument 4 must be a number"
                                            do general_counter_1 = 1,nelem
                                                do general_counter_2 = 1,nelem
                                                    sym_short_atomic_count(general_counter_1) = sym_short_atomic_count(general_counter_1) + 1
                                                    function_type_short_atomic(sym_short_atomic_count(general_counter_1),general_counter_1) = function_type_temp
                                                    eta_short_atomic(sym_short_atomic_count(general_counter_1),general_counter_1) = eta_temp
                                                    rshift_short_atomic(sym_short_atomic_count(general_counter_1),general_counter_1) = rshift_temp
                                                    funccutoff_short_atomic(sym_short_atomic_count(general_counter_1),general_counter_1) = funccutoff_temp
                                                    symelement_short_atomic(sym_short_atomic_count(general_counter_1),1,general_counter_1) = nucelem(general_counter_2)
                                                end do
                                            end do
                                        else
                                            print *, err, err_inpnn, "global_symfunction_short type ", words(2), " needs 4 arguments"; stop
                                        end if

                                    case ('3')
                                        if (nwords == 6) then
                                            read(words(3),*, iostat=ios) eta_temp
                                            if (ios /= 0) stop err // err_inpnn // "global_symfunction_short type ", words(2), " argument 2 must be a number"
                                            read(words(4),*, iostat=ios) lambda_temp
                                            if (ios /= 0) stop err // err_inpnn // "global_symfunction_short type ", words(2), " argument 3 must be a number"
                                            read(words(5),*, iostat=ios) zeta_temp
                                            if (ios /= 0) stop err // err_inpnn // "global_symfunction_short type ", words(2), " argument 4 must be a number"
                                            read(words(6),*, iostat=ios) funccutoff_temp
                                            if (ios /= 0) stop err // err_inpnn // "global_symfunction_short type ", words(2), " argument 5 must be a number"
                                            do general_counter_3 = 1,nelem
                                                do general_counter_1 = 1,nelem
                                                    sym_short_atomic_count(general_counter_3) = sym_short_atomic_count(general_counter_3) + 1
                                                    function_type_short_atomic(sym_short_atomic_count(general_counter_3),general_counter_3) = function_type_temp
                                                    eta_short_atomic(sym_short_atomic_count(general_counter_3),general_counter_3) = eta_temp
                                                    lambda_short_atomic(sym_short_atomic_count(general_counter_3),general_counter_3) = lambda_temp
                                                    zeta_short_atomic(sym_short_atomic_count(general_counter_3),general_counter_3) = zeta_temp
                                                    funccutoff_short_atomic(sym_short_atomic_count(general_counter_3),general_counter_3) = funccutoff_temp
                                                    symelement_short_atomic(sym_short_atomic_count(general_counter_3),1,general_counter_3) = nucelem(general_counter_1)
                                                    symelement_short_atomic(sym_short_atomic_count(general_counter_3),2,general_counter_3) = nucelem(general_counter_1)
                                                end do
                                                do general_counter_1 = 1,nelem
                                                    if (nelem .gt. 1) then
                                                        do general_counter_2 = 1,general_counter_1-1
                                                            sym_short_atomic_count(general_counter_3) = sym_short_atomic_count(general_counter_3) + 1
                                                            function_type_short_atomic(sym_short_atomic_count(general_counter_3),general_counter_3) = function_type_temp
                                                            eta_short_atomic(sym_short_atomic_count(general_counter_3),general_counter_3) = eta_temp
                                                            lambda_short_atomic(sym_short_atomic_count(general_counter_3),general_counter_3) = lambda_temp
                                                            zeta_short_atomic(sym_short_atomic_count(general_counter_3),general_counter_3) = zeta_temp
                                                            funccutoff_short_atomic(sym_short_atomic_count(general_counter_3),general_counter_3) = funccutoff_temp
                                                            symelement_short_atomic(sym_short_atomic_count(general_counter_3),1,general_counter_3) = nucelem(general_counter_2)
                                                            symelement_short_atomic(sym_short_atomic_count(general_counter_3),2,general_counter_3) = nucelem(general_counter_1)
                                                        end do
                                                    end if
                                                end do
                                            end do
                                        else
                                            print *, err, err_inpnn, "global_symfunction_short type ", words(2), " needs 5 arguments"; stop
                                        end if

                                    case ('4')
                                        if (nwords == 4) then
                                            read(words(3),*, iostat=ios) eta_temp
                                            if (ios /= 0) stop err // err_inpnn // "global_symfunction_short type ", words(2), " argument 2 must be a number"
                                            read(words(4),*, iostat=ios) funccutoff_temp
                                            if (ios /= 0) stop err // err_inpnn // "global_symfunction_short type ", words(2), " argument 3 must be a number"
                                            do general_counter_3 = 1,nelem
                                                do general_counter_1 = 1,nelem
                                                    sym_short_atomic_count(general_counter_3) = sym_short_atomic_count(general_counter_3) + 1
                                                    function_type_short_atomic(sym_short_atomic_count(general_counter_3),general_counter_3) = function_type_temp
                                                    eta_short_atomic(sym_short_atomic_count(general_counter_3),general_counter_3) = eta_temp
                                                    funccutoff_short_atomic(sym_short_atomic_count(general_counter_3),general_counter_3) = funccutoff_temp
                                                    symelement_short_atomic(sym_short_atomic_count(general_counter_3),1,general_counter_3) = nucelem(general_counter_1)
                                                end do
                                            end do
                                        else
                                            print *, err, err_inpnn, "global_symfunction_short type ", words(2), " needs 3 arguments"; stop
                                        end if

                                    case ('5')
                                        if (nwords == 3) then
                                            read(words(3),*, iostat=ios) eta_temp
                                            if (ios /= 0) stop err // err_inpnn // "global_symfunction_short type ", words(2), " argument 2 must be a number"
                                            do general_counter_3 = 1,nelem
                                                sym_short_atomic_count(general_counter_3) = sym_short_atomic_count(general_counter_3) + 1
                                                function_type_short_atomic(sym_short_atomic_count(general_counter_3),general_counter_3) = function_type_temp
                                                eta_short_atomic(sym_short_atomic_count(general_counter_3),general_counter_3) = eta_temp
                                                symelement_short_atomic(sym_short_atomic_count(general_counter_3),1,general_counter_3) = nucelem(general_counter_3)
                                            end do
                                        else
                                            print *, err, err_inpnn, "global_symfunction_short type ", words(2), " needs 2 arguments"; stop
                                        end if

                                    case ('6')
                                        if (nwords == 3) then
                                            read(words(3),*, iostat=ios) funccutoff_temp
                                            if (ios /= 0) stop err // err_inpnn // "global_symfunction_short type ", words(2), " argument 2 must be a number"
                                            do general_counter_3 = 1,nelem
                                                sym_short_atomic_count(general_counter_3) = sym_short_atomic_count(general_counter_3) + 1
                                                function_type_short_atomic(sym_short_atomic_count(general_counter_3),general_counter_3) = function_type_temp
                                                funccutoff_short_atomic(sym_short_atomic_count(general_counter_3),general_counter_3) = funccutoff_temp
                                                symelement_short_atomic(sym_short_atomic_count(general_counter_3),1,general_counter_3) = nucelem(general_counter_3)
                                            end do
                                        else
                                            print *, err, err_inpnn, "global_symfunction_short type ", words(2), " needs 2 arguments"; stop
                                        end if

                                    case ('8')
                                        if (nwords == 5) then
                                            read(words(3),*, iostat=ios) eta_temp
                                            if (ios /= 0) stop err // err_inpnn // "global_symfunction_short type ", words(2), " argument 2 must be a number"
                                            read(words(4),*, iostat=ios) rshift_temp
                                            if (ios /= 0) stop err // err_inpnn // "global_symfunction_short type ", words(2), " argument 3 must be a number"
                                            read(words(5),*, iostat=ios) funccutoff_temp
                                            if (ios /= 0) stop err // err_inpnn // "global_symfunction_short type ", words(2), " argument 4 must be a number"
                                            do general_counter_3 = 1,nelem
                                                do general_counter_1 = 1,nelem
                                                    sym_short_atomic_count(general_counter_3) = sym_short_atomic_count(general_counter_3) + 1
                                                    function_type_short_atomic(sym_short_atomic_count(general_counter_3),general_counter_3) = function_type_temp
                                                    eta_short_atomic(sym_short_atomic_count(general_counter_3),general_counter_3) = eta_temp
                                                    rshift_short_atomic(sym_short_atomic_count(general_counter_3),general_counter_3) = rshift_temp
                                                    funccutoff_short_atomic(sym_short_atomic_count(general_counter_3),general_counter_3) = funccutoff_temp
                                                    symelement_short_atomic(sym_short_atomic_count(general_counter_3),1,general_counter_3) = nucelem(general_counter_1)
                                                    symelement_short_atomic(sym_short_atomic_count(general_counter_3),2,general_counter_3) = nucelem(general_counter_1)
                                                end do
                                                do general_counter_1 = 1,nelem
                                                    if (nelem .gt. 1) then
                                                        do general_counter_2 = 1,general_counter_1-1
                                                            sym_short_atomic_count(general_counter_3) = sym_short_atomic_count(general_counter_3) + 1
                                                            function_type_short_atomic(sym_short_atomic_count(general_counter_3),general_counter_3) = function_type_temp
                                                            eta_short_atomic(sym_short_atomic_count(general_counter_3),general_counter_3) = eta_temp
                                                            rshift_short_atomic(sym_short_atomic_count(general_counter_3),general_counter_3) = rshift_temp
                                                            funccutoff_short_atomic(sym_short_atomic_count(general_counter_3),general_counter_3) = funccutoff_temp
                                                            symelement_short_atomic(sym_short_atomic_count(general_counter_3),1,general_counter_3) = nucelem(general_counter_2)
                                                            symelement_short_atomic(sym_short_atomic_count(general_counter_3),2,general_counter_3) = nucelem(general_counter_1)
                                                        end do
                                                    end if
                                                end do
                                            end do
                                        else
                                            print *, err, err_inpnn, "global_symfunction_short type ", words(2), " needs 4 arguments"; stop
                                        end if

                                    case ('9')
                                        if (nwords == 6) then
                                            read(words(3),*, iostat=ios) eta_temp
                                            if (ios /= 0) stop err // err_inpnn // "global_symfunction_short type ", words(2), " argument 2 must be a number"
                                            read(words(4),*, iostat=ios) lambda_temp
                                            if (ios /= 0) stop err // err_inpnn // "global_symfunction_short type ", words(2), " argument 3 must be a number"
                                            read(words(5),*, iostat=ios) zeta_temp
                                            if (ios /= 0) stop err // err_inpnn // "global_symfunction_short type ", words(2), " argument 4 must be a number"
                                            read(words(6),*, iostat=ios) funccutoff_temp
                                            if (ios /= 0) stop err // err_inpnn // "global_symfunction_short type ", words(2), " argument 5 must be a number"
                                            do general_counter_3 = 1,nelem
                                                do general_counter_1 = 1,nelem
                                                    sym_short_atomic_count(general_counter_3) = sym_short_atomic_count(general_counter_3) + 1
                                                    function_type_short_atomic(sym_short_atomic_count(general_counter_3),general_counter_3) = function_type_temp
                                                    eta_short_atomic(sym_short_atomic_count(general_counter_3),general_counter_3) = eta_temp
                                                    lambda_short_atomic(sym_short_atomic_count(general_counter_3),general_counter_3) = lambda_temp
                                                    zeta_short_atomic(sym_short_atomic_count(general_counter_3),general_counter_3) = zeta_temp
                                                    funccutoff_short_atomic(sym_short_atomic_count(general_counter_3),general_counter_3) = funccutoff_temp
                                                    symelement_short_atomic(sym_short_atomic_count(general_counter_3),1,general_counter_3) = nucelem(general_counter_1)
                                                    symelement_short_atomic(sym_short_atomic_count(general_counter_3),2,general_counter_3) = nucelem(general_counter_1)
                                                end do
                                                do general_counter_1 = 1,nelem
                                                    if (nelem .gt. 1) then
                                                        do general_counter_2 = 1,general_counter_1-1
                                                            sym_short_atomic_count(general_counter_3) = sym_short_atomic_count(general_counter_3) + 1
                                                            function_type_short_atomic(sym_short_atomic_count(general_counter_3),general_counter_3) = function_type_temp
                                                            eta_short_atomic(sym_short_atomic_count(general_counter_3),general_counter_3) = eta_temp
                                                            lambda_short_atomic(sym_short_atomic_count(general_counter_3),general_counter_3) = lambda_temp
                                                            zeta_short_atomic(sym_short_atomic_count(general_counter_3),general_counter_3) = zeta_temp
                                                            funccutoff_short_atomic(sym_short_atomic_count(general_counter_3),general_counter_3) = funccutoff_temp
                                                            symelement_short_atomic(sym_short_atomic_count(general_counter_3),1,general_counter_3) = nucelem(general_counter_2)
                                                            symelement_short_atomic(sym_short_atomic_count(general_counter_3),2,general_counter_3) = nucelem(general_counter_1)
                                                        end do
                                                    end if
                                                end do
                                            end do
                                        else
                                            print *, err, err_inpnn, "global_symfunction_short type ", words(2), " needs 5 arguments"; stop
                                        end if

                                    case default
                                        print *, err, err_inpnn, "Error in global_symfunction_short key, symfunction type ", words(2), " not implemented"
                                        stop

                                end select

                            end if


                        case ('symfunction_electrostatic')
                            if (lelec .and. (nn_type_elec == 1)) then
                                read(words(2),'(A)', iostat=ios) rinpparam%elementtemp1
                                call checkelement(elementtemp1)
                                call nuccharge(elementtemp1,ztemp1)
                                sym_elec_count(elementindex(ztemp1)) = sym_elec_count(elementindex(ztemp1)) + 1
                                read(words(3),'(i1000)', iostat=ios) function_type_elec(sym_elec_count(elementindex(ztemp1)),elementindex(ztemp1))
                                if (ios /= 0) stop err // err_inpnn // "symfunction_electrostatic argument 2 value must be integer"

                                select case(words(3))

                                    case ('1')
                                        if (nwords == 5) then
                                            read(words(4),'(A)') rinpparam%elementtemp2
                                            call checkelement(elementtemp2)
                                            call nuccharge(elementtemp2,ztemp2)
                                            read(words(5),*, iostat=ios) rinpparam%funccutoff_elec(sym_elec_count(elementindex(ztemp1)),elementindex(ztemp1))
                                            if (ios /= 0) stop err // err_inpnn // "symfunction_electrostatic type ", words(3), " argument 4 must be a number"
                                            symelement_elec(sym_elec_count(elementindex(ztemp1)),1,elementindex(ztemp1))=ztemp2
                                        else
                                            print *, err, err_inpnn, "symfunction_electrostatic type ", words(3), " needs 4 arguments"; stop
                                        end if

                                    case ('2')
                                        if (nwords == 7) then
                                            read(words(4),'(A)') elementtemp2
                                            call checkelement(elementtemp2)
                                            call nuccharge(elementtemp2,ztemp2)
                                            read(words(5),*, iostat=ios) rinpparam%eta_elec(sym_elec_count(elementindex(ztemp1)),elementindex(ztemp1))
                                            if (ios /= 0) stop err // err_inpnn // "symfunction_electrostatic type ", words(3), " argument 4 must be a number"
                                            read(words(6),*, iostat=ios) rinpparam%rshift_elec(sym_elec_count(elementindex(ztemp1)),elementindex(ztemp1))
                                            if (ios /= 0) stop err // err_inpnn // "symfunction_electrostatic type ", words(3), " argument 5 must be a number"
                                            read(words(7),*, iostat=ios) rinpparam%funccutoff_elec(sym_elec_count(elementindex(ztemp1)),elementindex(ztemp1))
                                            if (ios /= 0) stop err // err_inpnn // "symfunction_electrostatic type ", words(3), " argument 6 must be a number"
                                            symelement_elec(sym_elec_count(elementindex(ztemp1)),1,elementindex(ztemp1))=ztemp2
                                        else
                                            print *, err, err_inpnn, "symfunction_electrostatic type ", words(3), " needs 6 arguments"; stop
                                        end if

                                    case ('3')
                                        if (nwords == 9) then
                                            read(words(4),'(A)') rinpparam%elementtemp2
                                            read(words(5),'(A)') rinpparam%elementtemp3
                                            call checkelement(elementtemp2)
                                            call checkelement(elementtemp3)
                                            call nuccharge(elementtemp2,ztemp2)
                                            call nuccharge(elementtemp2,ztemp3)
                                            if (ztemp2 .gt. ztemp3) then
                                                itemp = ztemp2
                                                ztemp2 = ztemp3
                                                ztemp3 = itemp
                                            endif
                                            read(words(6),*, iostat=ios) rinpparam%eta_elec(sym_elec_count(elementindex(ztemp1)),elementindex(ztemp1))
                                            if (ios /= 0) stop err // err_inpnn // "symfunction_electrostatic type ", words(3), " argument 5 must be a number"
                                            read(words(7),*, iostat=ios) rinpparam%lambda_elec(sym_elec_count(elementindex(ztemp1)),elementindex(ztemp1))
                                            if (ios /= 0) stop err // err_inpnn // "symfunction_electrostatic type ", words(3), " argument 6 must be a number"
                                            read(words(8),*, iostat=ios) rinpparam%zeta_elec(sym_elec_count(elementindex(ztemp1)),elementindex(ztemp1))
                                            if (ios /= 0) stop err // err_inpnn // "symfunction_electrostatic type ", words(3), " argument 7 must be a number"
                                            read(words(9),*, iostat=ios) rinpparam%funccutoff_elec(sym_elec_count(elementindex(ztemp1)),elementindex(ztemp1))
                                            if (ios /= 0) stop err // err_inpnn // "symfunction_electrostatic type ", words(3), " argument 8 must be a number"
                                            symelement_elec(sym_elec_count(elementindex(ztemp1)),1,elementindex(ztemp1))=ztemp2
                                            symelement_elec(sym_elec_count(elementindex(ztemp1)),2,elementindex(ztemp1))=ztemp3
                                        else
                                            print *, err, err_inpnn, "symfunction_electrostatic type ", words(3), " needs 8 arguments"; stop
                                        end if

                                    case ('4')
                                        if (nwords == 6) then
                                            read(words(4),'(A)') rinpparam%elementtemp2
                                            call checkelement(elementtemp2)
                                            call nuccharge(elementtemp2,ztemp2)
                                            read(words(5),*, iostat=ios) eta_elec(sym_elec_count(elementindex(ztemp1)),elementindex(ztemp1))
                                            if (ios /= 0) stop err // err_inpnn // "symfunction_electrostatic type ", words(3), " argument 4 must be a number"
                                            read(words(6),*, iostat=ios) rinpparam%funccutoff_elecc(sym_elec_count(elementindex(ztemp1)),elementindex(ztemp1))
                                            if (ios /= 0) stop err // err_inpnn // "symfunction_electrostatic type ", words(3), " argument 5 must be a number"
                                            symelement_elec(sym_elec_count(elementindex(ztemp1)),1,elementindex(ztemp1))=ztemp2
                                        else
                                            print *, err, err_inpnn, "symfunction_electrostatic type ", words(3), " needs 5 arguments"; stop
                                        end if

                                    case ('5')
                                        if (nwords == 4) then
                                            read(words(4),*, iostat=ios) rinpparam%eta_elec(sym_elec_count(elementindex(ztemp1)),elementindex(ztemp1))
                                            if (ios /= 0) stop err // err_inpnn // "symfunction_electrostatic type ", words(3), " argument 3 must be a number"
                                        else
                                            print *, err, err_inpnn, "symfunction_electrostatic type ", words(3), " needs 3 arguments"; stop
                                        end if

                                    case ('6')
                                        if (nwords == 5) then
                                            read(words(4),'(A)') rinpparam%elementtemp2
                                            call checkelement(elementtemp2)
                                            call nuccharge(elementtemp2,ztemp2)
                                            read(words(5),*, iostat=ios) rinpparam%funccutoff_elec(sym_elec_count(elementindex(ztemp1)),elementindex(ztemp1))
                                            if (ios /= 0) stop err // err_inpnn // "symfunction_electrostatic type ", words(3), " argument 4 must be a number"
                                            symelement_elec(sym_elec_count(elementindex(ztemp1)),1,elementindex(ztemp1))=ztemp2
                                        else
                                            print *, err, err_inpnn, "symfunction_electrostatic type ", words(3), " needs 4 arguments"; stop
                                        end if

                                    case ('8')
                                        if (nwords == 8) then
                                            read(words(4),'(A)') rinpparam%elementtemp2
                                            read(words(5),'(A)') rinpparam%elementtemp3
                                            call checkelement(elementtemp2)
                                            call checkelement(elementtemp3)
                                            call nuccharge(elementtemp2,ztemp2)
                                            call nuccharge(elementtemp2,ztemp3)
                                            if (ztemp2 .gt. ztemp3) then
                                                itemp = ztemp2
                                                ztemp2 = ztemp3
                                                ztemp3 = itemp
                                            endif
                                            read(words(6),*, iostat=ios) rinpparam%eta_elec(sym_elec_count(elementindex(ztemp1)),elementindex(ztemp1))
                                            if (ios /= 0) stop err // err_inpnn // "symfunction_electrostatic type ", words(3), " argument 5 must be a number"
                                            read(words(7),*, iostat=ios) rinpparam%rshift_elec(sym_elec_count(elementindex(ztemp1)),elementindex(ztemp1))
                                            if (ios /= 0) stop err // err_inpnn // "symfunction_electrostatic type ", words(3), " argument 6 must be a number"
                                            read(words(8),*, iostat=ios) rinpparam%funccutoff_elec(sym_elec_count(elementindex(ztemp1)),elementindex(ztemp1))
                                            if (ios /= 0) stop err // err_inpnn // "symfunction_electrostatic type ", words(3), " argument 7 must be a number"
                                            symelement_elec(sym_elec_count(elementindex(ztemp1)),1,elementindex(ztemp1))=ztemp2
                                            symelement_elec(sym_elec_count(elementindex(ztemp1)),2,elementindex(ztemp1))=ztemp3
                                        else
                                            print *, err, err_inpnn, "symfunction_electrostatic type ", words(3), " needs 7 arguments"; stop
                                        end if

                                    case ('9')
                                        if (nwords == 9) then
                                            read(words(4),'(A)') rinpparam%elementtemp2
                                            read(words(5),'(A)') rinpparam%elementtemp3
                                            call checkelement(elementtemp2)
                                            call checkelement(elementtemp3)
                                            call nuccharge(elementtemp2,ztemp2)
                                            call nuccharge(elementtemp2,ztemp3)
                                            if (ztemp2 .gt. ztemp3) then
                                                itemp = ztemp2
                                                ztemp2 = ztemp3
                                                ztemp3 = itemp
                                            endif
                                            read(words(6),*, iostat=ios) rinpparam%eta_elec(sym_elec_count(elementindex(ztemp1)),elementindex(ztemp1))
                                            if (ios /= 0) stop err // err_inpnn // "symfunction_electrostatic type ", words(3), " argument 5 must be a number"
                                            read(words(7),*, iostat=ios) rinpparam%lambda_elec(sym_elec_count(elementindex(ztemp1)),elementindex(ztemp1))
                                            if (ios /= 0) stop err // err_inpnn // "symfunction_electrostatic type ", words(3), " argument 6 must be a number"
                                            read(words(8),*, iostat=ios) rinpparam%zeta_elec(sym_elec_count(elementindex(ztemp1)),elementindex(ztemp1))
                                            if (ios /= 0) stop err // err_inpnn // "symfunction_electrostatic type ", words(3), " argument 7 must be a number"
                                            read(words(9),*, iostat=ios) rinpparam%funccutoff_elec(sym_elec_count(elementindex(ztemp1)),elementindex(ztemp1))
                                            if (ios /= 0) stop err // err_inpnn // "symfunction_electrostatic type ", words(3), " argument 8 must be a number"
                                            symelement_elec(sym_elec_count(elementindex(ztemp1)),1,elementindex(ztemp1))=ztemp2
                                            symelement_elec(sym_elec_count(elementindex(ztemp1)),2,elementindex(ztemp1))=ztemp3
                                        else
                                            print *, err, err_inpnn, "symfunction_electrostatic type ", words(3), " needs 8 arguments"; stop
                                        end if

                                    case default
                                        print *, err, err_inpnn, "Error in symfunction_electrostatict key, symfunction type ", words(3), " not implemented"
                                        stop

                                end select

                            end if

                        case ('element_symfunction_electrostatic')
                            if (lelec .and. (nn_type_elec == 1)) then
                                read(words(2),'(A)', iostat=ios) elementtemp1
                                call checkelement(elementtemp1)
                                call nuccharge(elementtemp1,ztemp1)
                                read(words(3),'(i1000)', iostat=ios) function_type_temp
                                if (ios /= 0) stop err // err_inpnn // "element_symfunction_electrostatic argument 2 value must be integer"

                                select case(words(3))

                                    case ('1')
                                        if (nwords == 4) then
                                            read(words(4),*, iostat=ios) funccutoff_temp
                                            if (ios /= 0) stop err // err_inpnn // "element_symfunction_electrostatic type ", words(3), " argument 3 must be a number"
                                            do general_counter_1 = 1,nelem
                                                sym_elec_count(elementindex(ztemp1) = sym_elec_count(elementindex(ztemp1) + 1
                                                function_type_elec(sym_elec_count(elementindex(ztemp1)),elementindex(ztemp1)) = function_type_temp
                                                funccutoff_elec(sym_elec_count(elementindex(ztemp1)),elementindex(ztemp1)) = funccutoff_temp
                                                symelement_elec(sym_elec_count(elementindex(ztemp1)),1,elementindex(ztemp1)) = nucelem(general_counter_1)
                                            end do
                                        else
                                            print *, err, err_inpnn, "element_symfunction_electrostatic type ", words(3), " needs 3 arguments"; stop
                                        end if

                                    case ('2')
                                        if (nwords == 6) then
                                            read(words(4),*, iostat=ios) eta_temp
                                            if (ios /= 0) stop err // err_inpnn // "element_symfunction_electrostatic type ", words(3), " argument 3 must be a number"
                                            read(words(5),*, iostat=ios) rshift_temp
                                            if (ios /= 0) stop err // err_inpnn // "element_symfunction_electrostatic type ", words(3), " argument 4 must be a number"
                                            read(words(6),*, iostat=ios) funccutoff_temp
                                            if (ios /= 0) stop err // err_inpnn // "element_symfunction_electrostatic type ", words(3), " argument 5 must be a number"
                                            do general_counter_1 = 1,nelem
                                                sym_elec_count(elementindex(ztemp1) = sym_elec_count(elementindex(ztemp1) + 1
                                                function_type_elec(sym_elec_count(elementindex(ztemp1)),elementindex(ztemp1)) = function_type_temp
                                                eta_elec(sym_elec_count(elementindex(ztemp1)),elementindex(ztemp1)) = eta_temp
                                                rshift_elec(sym_elec_count(elementindex(ztemp1)),elementindex(ztemp1)) = rshift_temp
                                                funccutoff_elec(sym_elec_count(elementindex(ztemp1)),elementindex(ztemp1)) = funccutoff_temp
                                                symelement_elec(sym_elec_count(elementindex(ztemp1)),1,elementindex(ztemp1)) = nucelem(general_counter_1)
                                            end do
                                        else
                                            print *, err, err_inpnn, "element_symfunction_electrostatic type ", words(3), " needs 5 arguments"; stop
                                        end if

                                    case ('3')
                                        if (nwords == 7) then
                                            read(words(4),*, iostat=ios) eta_temp
                                            if (ios /= 0) stop err // err_inpnn // "element_symfunction_electrostatic type ", words(3), " argument 3 must be a number"
                                            read(words(5),*, iostat=ios) lambda_temp
                                            if (ios /= 0) stop err // err_inpnn // "element_symfunction_electrostatic type ", words(3), " argument 4 must be a number"
                                            read(words(6),*, iostat=ios) zeta_temp
                                            if (ios /= 0) stop err // err_inpnn // "element_symfunction_electrostatic type ", words(3), " argument 5 must be a number"
                                            read(words(7),*, iostat=ios) funccutoff_temp
                                            if (ios /= 0) stop err // err_inpnn // "element_symfunction_electrostatic type ", words(3), " argument 6 must be a number"
                                            do general_counter_1 = 1,nelem
                                                sym_elec_count(elementindex(ztemp1) = sym_elec_count(elementindex(ztemp1) + 1
                                                function_type_elec(sym_elec_count(elementindex(ztemp1)),elementindex(ztemp1)) = function_type_temp
                                                eta_elec(sym_elec_count(elementindex(ztemp1)),elementindex(ztemp1)) = eta_temp
                                                lambda_elec(sym_elec_count(elementindex(ztemp1)),elementindex(ztemp1)) = lambda_temp
                                                zeta_elec(sym_elec_count(elementindex(ztemp1)),elementindex(ztemp1)) = zeta_temp
                                                funccutoff_elec(sym_elec_count(elementindex(ztemp1)),elementindex(ztemp1)) = funccutoff_temp
                                                symelement_elec(sym_elec_count(elementindex(ztemp1)),1,elementindex(ztemp1)) = nucelem(general_counter_1)
                                                symelement_elec(sym_elec_count(elementindex(ztemp1)),2,elementindex(ztemp1)) = nucelem(general_counter_1)
                                            end do
                                            do general_counter_1 = 1,nelem
                                                if (nelem .gt. 1) then
                                                    do general_counter_2 = 1,general_counter_1-1
                                                        sym_elec_count(elementindex(ztemp1) = sym_elec_count(elementindex(ztemp1) + 1
                                                        function_type_elec(sym_elec_count(elementindex(ztemp1)),elementindex(ztemp1)) = function_type_temp
                                                        eta_elec(sym_elec_count(elementindex(ztemp1)),elementindex(ztemp1)) = eta_temp
                                                        lambda_elec(sym_elec_count(elementindex(ztemp1)),elementindex(ztemp1)) = lambda_temp
                                                        zeta_elec(sym_elec_count(elementindex(ztemp1)),elementindex(ztemp1)) = zeta_temp
                                                        funccutoff_elec(sym_elec_count(elementindex(ztemp1)),elementindex(ztemp1)) = funccutoff_temp
                                                        symelement_elec(sym_elec_count(elementindex(ztemp1)),1,elementindex(ztemp1)) = nucelem(general_counter_2)
                                                        symelement_elec(sym_elec_count(elementindex(ztemp1)),2,elementindex(ztemp1)) = nucelem(general_counter_1)
                                                    end do
                                                end if
                                            end do
                                        else
                                            print *, err, err_inpnn, "element_symfunction_electrostatic type ", words(3), " needs 6 arguments"; stop
                                        end if

                                    case ('4')
                                        if (nwords == 5) then
                                            read(words(4),*, iostat=ios) eta_temp
                                            if (ios /= 0) stop err // err_inpnn // "element_symfunction_electrostatic type ", words(3), " argument 3 must be a number"
                                            read(words(5),*, iostat=ios) funccutoff_temp
                                            if (ios /= 0) stop err // err_inpnn // "element_symfunction_electrostatic type ", words(3), " argument 4 must be a number"
                                            do general_counter_1 = 1,nelem
                                                sym_elec_count(elementindex(ztemp1) = sym_elec_count(elementindex(ztemp1) + 1
                                                function_type_elec(sym_elec_count(elementindex(ztemp1)),elementindex(ztemp1)) = function_type_temp
                                                eta_elec(sym_elec_count(elementindex(ztemp1)),elementindex(ztemp1)) = eta_temp
                                                funccutoff_elec(sym_elec_count(elementindex(ztemp1)),elementindex(ztemp1)) = funccutoff_temp
                                                symelement_elec(sym_elec_count(elementindex(ztemp1)),1,elementindex(ztemp1)) = nucelem(general_counter_1)
                                            end do
                                        else
                                            print *, err, err_inpnn, "element_symfunction_electrostatic type ", words(3), " needs 4 arguments"; stop
                                        end if

                                    case ('5')
                                        if (nwords == 4) then
                                            read(words(4),*, iostat=ios) eta_temp
                                            if (ios /= 0) stop err // err_inpnn // "element_symfunction_electrostatic type ", words(3), " argument 3 must be a number"
                                            sym_elec_count(elementindex(ztemp1) = sym_elec_count(elementindex(ztemp1) + 1
                                            function_type_elec(sym_elec_count(elementindex(ztemp1)),elementindex(ztemp1)) = function_type_temp
                                            eta_elec(sym_elec_count(elementindex(ztemp1)),elementindex(ztemp1)) = eta_temp
                                        else
                                            print *, err, err_inpnn, "element_symfunction_electrostatic type ", words(3), " needs 3 arguments"; stop
                                        end if

                                    case ('6')
                                        if (nwords == 4) then
                                            read(words(4),*, iostat=ios) funccutoff_temp
                                            if (ios /= 0) stop err // err_inpnn // "element_symfunction_electrostatic type ", words(3), " argument 3 must be a number"
                                            sym_elec_count(elementindex(ztemp1) = sym_elec_count(elementindex(ztemp1) + 1
                                            function_type_elec(sym_elec_count(elementindex(ztemp1)),elementindex(ztemp1)) = function_type_temp
                                        else
                                            print *, err, err_inpnn, "element_symfunction_electrostatic type ", words(3), " needs 3 arguments"; stop
                                        end if

                                    case ('8')
                                        if (nwords == 6) then
                                            read(words(4),*, iostat=ios) eta_temp
                                            if (ios /= 0) stop err // err_inpnn // "element_symfunction_electrostatic type ", words(3), " argument 3 must be a number"
                                            read(words(5),*, iostat=ios) rshift_temp
                                            if (ios /= 0) stop err // err_inpnn // "element_symfunction_electrostatic type ", words(3), " argument 4 must be a number"
                                            read(words(6),*, iostat=ios) funccutoff_temp
                                            if (ios /= 0) stop err // err_inpnn // "element_symfunction_electrostatic type ", words(3), " argument 5 must be a number"
                                            do general_counter_1 = 1,nelem
                                                sym_elec_count(elementindex(ztemp1) = sym_elec_count(elementindex(ztemp1) + 1
                                                function_type_elec(sym_elec_count(elementindex(ztemp1)),elementindex(ztemp1)) = function_type_temp
                                                eta_elec(sym_elec_count(elementindex(ztemp1)),elementindex(ztemp1)) = eta_temp
                                                rshift_elec(sym_elec_count(elementindex(ztemp1)),elementindex(ztemp1)) = rshift_temp
                                                funccutoff_elec(sym_elec_count(elementindex(ztemp1)),elementindex(ztemp1)) = funccutoff_temp
                                                symelement_elec(sym_elec_count(elementindex(ztemp1)),1,elementindex(ztemp1)) = nucelem(general_counter_1)
                                                symelement_elec(sym_elec_count(elementindex(ztemp1)),2,elementindex(ztemp1)) = nucelem(general_counter_1)
                                            end do
                                            do general_counter_1 = 1,nelem
                                                if (nelem .gt. 1) then
                                                    do general_counter_2 = 1,general_counter_1-1
                                                        sym_elec_count(elementindex(ztemp1) = sym_elec_count(elementindex(ztemp1) + 1
                                                        function_type_elec(sym_elec_count(elementindex(ztemp1)),elementindex(ztemp1)) = function_type_temp
                                                        eta_elec(sym_elec_count(elementindex(ztemp1)),elementindex(ztemp1)) = eta_temp
                                                        rshift_elec(sym_elec_count(elementindex(ztemp1)),elementindex(ztemp1)) = rshift_temp
                                                        funccutoff_elec(sym_elec_count(elementindex(ztemp1)),elementindex(ztemp1)) = funccutoff_temp
                                                        symelement_elec(sym_elec_count(elementindex(ztemp1)),1,elementindex(ztemp1)) = nucelem(general_counter_2)
                                                        symelement_elec(sym_elec_count(elementindex(ztemp1)),2,elementindex(ztemp1)) = nucelem(general_counter_1)
                                                    end do
                                                end if
                                            end do
                                        else
                                            print *, err, err_inpnn, "element_symfunction_electrostatic type ", words(3), " needs 5 arguments"; stop
                                        end if

                                    case ('9')
                                        if (nwords == 7) then
                                            read(words(4),*, iostat=ios) eta_temp
                                            if (ios /= 0) stop err // err_inpnn // "element_symfunction_electrostatic type ", words(3), " argument 3 must be a number"
                                            read(words(5),*, iostat=ios) lambda_temp
                                            if (ios /= 0) stop err // err_inpnn // "element_symfunction_electrostatic type ", words(3), " argument 4 must be a number"
                                            read(words(6),*, iostat=ios) zeta_temp
                                            if (ios /= 0) stop err // err_inpnn // "element_symfunction_electrostatic type ", words(3), " argument 5 must be a number"
                                            read(words(7),*, iostat=ios) funccutoff_temp
                                            if (ios /= 0) stop err // err_inpnn // "element_symfunction_electrostatic type ", words(3), " argument 6 must be a number"
                                            do general_counter_1 = 1,nelem
                                                sym_elec_count(elementindex(ztemp1) = sym_elec_count(elementindex(ztemp1) + 1
                                                function_type_elec(sym_elec_count(elementindex(ztemp1)),elementindex(ztemp1)) = function_type_temp
                                                eta_elec(sym_elec_count(elementindex(ztemp1)),elementindex(ztemp1)) = eta_temp
                                                lambda_elec(sym_elec_count(elementindex(ztemp1)),elementindex(ztemp1)) = lambda_temp
                                                zeta_elec(sym_elec_count(elementindex(ztemp1)),elementindex(ztemp1)) = zeta_temp
                                                funccutoff_elec(sym_elec_count(elementindex(ztemp1)),elementindex(ztemp1)) = funccutoff_temp
                                                symelement_elec(sym_elec_count(elementindex(ztemp1)),1,elementindex(ztemp1)) = nucelem(general_counter_1)
                                            end do
                                            do general_counter_1 = 1,nelem
                                                if (nelem .gt. 1) then
                                                    do general_counter_2 = 1,general_counter_1-1
                                                        sym_elec_count(elementindex(ztemp1) = sym_elec_count(elementindex(ztemp1) + 1
                                                        function_type_elec(sym_elec_count(elementindex(ztemp1)),elementindex(ztemp1)) = function_type_temp
                                                        eta_elec(sym_elec_count(elementindex(ztemp1)),elementindex(ztemp1)) = eta_temp
                                                        lambda_elec(sym_elec_count(elementindex(ztemp1)),elementindex(ztemp1)) = lambda_temp
                                                        zeta_elec(sym_elec_count(elementindex(ztemp1)),elementindex(ztemp1)) = zeta_temp
                                                        funccutoff_elec(sym_elec_count(elementindex(ztemp1)),elementindex(ztemp1)) = funccutoff_temp
                                                        symelement_elec(sym_elec_count(elementindex(ztemp1)),1,elementindex(ztemp1)) = nucelem(general_counter_2)
                                                        symelement_elec(sym_elec_count(elementindex(ztemp1)),2,elementindex(ztemp1)) = nucelem(general_counter_1)
                                                    end do
                                                end if
                                            end do
                                        else
                                            print *, err, err_inpnn, "element_symfunction_electrostatic type ", words(3), " needs 6 arguments"; stop
                                        end if

                                    case default
                                        print *, err, err_inpnn, "Error in element_symfunction_electrostatic key, symfunction type ", words(3), " not implemented"
                                        stop

                                end select

                            end if

                        case ('global_symfunction_electrostatic')
                            if (lelec .and. (nn_type_elec == 1)) then
                                read(words(2),'(i1000)', iostat=ios) function_type_temp
                                if (ios /= 0) stop err // err_inpnn // "global_symfunction_electrostatic argument 1 value must be integer"

                                select case(words(2))

                                    case ('1')
                                        if (nwords == 3) then
                                            read(words(3),*, iostat=ios) funccutoff_temp
                                            if (ios /= 0) stop err // err_inpnn // "global_symfunction_electrostatic type ", words(2), " argument 2 must be a number"
                                            do general_counter_1 = 1,nelem
                                                do general_counter_2 = 1,nelem
                                                    sym_elec_count(general_counter_1) = sym_elec_count(general_counter_1) + 1
                                                    function_type_elec(sym_elec_count(general_counter_1),general_counter_1) = function_type_temp
                                                    funccutoff_elec(sym_elec_count(general_counter_1),general_counter_1) = funccutoff_temp
                                                    symelement_elec(sym_elec_count(general_counter_1),1,general_counter_1) = nucelem(general_counter_2)
                                                end do
                                            end do
                                        else
                                            print *, err, err_inpnn, "global_symfunction_electrostatic type ", words(2), " needs 2 arguments"; stop
                                        end if

                                    case ('2')
                                        if (nwords == 5) then
                                            read(words(3),*, iostat=ios) eta_temp
                                            if (ios /= 0) stop err // err_inpnn // "global_symfunction_electrostatic type ", words(2), " argument 2 must be a number"
                                            read(words(4),*, iostat=ios) rshift_temp
                                            if (ios /= 0) stop err // err_inpnn // "global_symfunction_electrostatic type ", words(2), " argument 3 must be a number"
                                            read(words(5),*, iostat=ios) funccutoff_temp
                                            if (ios /= 0) stop err // err_inpnn // "global_symfunction_electrostatic type ", words(2), " argument 4 must be a number"
                                            do general_counter_1 = 1,nelem
                                                do general_counter_2 = 1,nelem
                                                    sym_elec_count(general_counter_1) = sym_elec_count(general_counter_1) + 1
                                                    function_type_elec(sym_elec_count(general_counter_1),general_counter_1) = function_type_temp
                                                    eta_elec(sym_elec_count(general_counter_1),general_counter_1) = eta_temp
                                                    rshift_elec(sym_elec_count(general_counter_1),general_counter_1) = rshift_temp
                                                    funccutoff_elec(sym_elec_count(general_counter_1),general_counter_1) = funccutoff_temp
                                                    symelement_elec(sym_elec_count(general_counter_1),1,general_counter_1) = nucelem(general_counter_2)
                                                end do
                                            end do
                                        else
                                            print *, err, err_inpnn, "global_symfunction_electrostatic type ", words(2), " needs 4 arguments"; stop
                                        end if

                                    case ('3')
                                        if (nwords == 6) then
                                            read(words(3),*, iostat=ios) eta_temp
                                            if (ios /= 0) stop err // err_inpnn // "global_symfunction_electrostatic type ", words(2), " argument 2 must be a number"
                                            read(words(4),*, iostat=ios) lambda_temp
                                            if (ios /= 0) stop err // err_inpnn // "global_symfunction_electrostatic type ", words(2), " argument 3 must be a number"
                                            read(words(5),*, iostat=ios) zeta_temp
                                            if (ios /= 0) stop err // err_inpnn // "global_symfunction_electrostatic type ", words(2), " argument 4 must be a number"
                                            read(words(6),*, iostat=ios) funccutoff_temp
                                            if (ios /= 0) stop err // err_inpnn // "global_symfunction_electrostatic type ", words(2), " argument 5 must be a number"
                                            do general_counter_3 = 1,nelem
                                                do general_counter_1 = 1,nelem
                                                    sym_elec_count(general_counter_3) = sym_elec_count(general_counter_3) + 1
                                                    function_type_elec(sym_elec_count(general_counter_3),general_counter_3) = function_type_temp
                                                    eta_elec(sym_elec_count(general_counter_3),general_counter_3) = eta_temp
                                                    lambda_elec(sym_elec_count(general_counter_3),general_counter_3) = lambda_temp
                                                    zeta_elec(sym_elec_count(general_counter_3),general_counter_3) = zeta_temp
                                                    funccutoff_elec(sym_elec_count(general_counter_3),general_counter_3) = funccutoff_temp
                                                    symelement_elec(sym_elec_count(general_counter_3),1,general_counter_3) = nucelem(general_counter_1)
                                                    symelement_elec(sym_elec_count(general_counter_3),2,general_counter_3) = nucelem(general_counter_1)
                                                end do
                                                do general_counter_1 = 1,nelem
                                                    if (nelem .gt. 1) then
                                                        do general_counter_2 = 1,general_counter_1-1
                                                            sym_elec_count(general_counter_3) = sym_elec_count(general_counter_3) + 1
                                                            function_type_elec(sym_elec_count(general_counter_3),general_counter_3) = function_type_temp
                                                            eta_elec(sym_elec_count(general_counter_3),general_counter_3) = eta_temp
                                                            lambda_elec(sym_elec_count(general_counter_3),general_counter_3) = lambda_temp
                                                            zeta_elec(sym_elec_count(general_counter_3),general_counter_3) = zeta_temp
                                                            funccutoff_elec(sym_elec_count(general_counter_3),general_counter_3) = funccutoff_temp
                                                            symelement_elec(sym_elec_count(general_counter_3),1,general_counter_3) = nucelem(general_counter_2)
                                                            symelement_elec(sym_elec_count(general_counter_3),2,general_counter_3) = nucelem(general_counter_1)
                                                        end do
                                                    end if
                                                end do
                                            end do
                                        else
                                            print *, err, err_inpnn, "global_symfunction_electrostatic type ", words(2), " needs 5 arguments"; stop
                                        end if

                                    case ('4')
                                        if (nwords == 4) then
                                            read(words(3),*, iostat=ios) eta_temp
                                            if (ios /= 0) stop err // err_inpnn // "global_symfunction_electrostatic type ", words(2), " argument 2 must be a number"
                                            read(words(4),*, iostat=ios) funccutoff_temp
                                            if (ios /= 0) stop err // err_inpnn // "global_symfunction_electrostatic type ", words(2), " argument 3 must be a number"
                                            do general_counter_3 = 1,nelem
                                                do general_counter_1 = 1,nelem
                                                    sym_elec_count(general_counter_3) = sym_elec_count(general_counter_3) + 1
                                                    function_type_elec(sym_elec_count(general_counter_3),general_counter_3) = function_type_temp
                                                    eta_elec(sym_elec_count(general_counter_3),general_counter_3) = eta_temp
                                                    funccutoff_elec(sym_elec_count(general_counter_3),general_counter_3) = funccutoff_temp
                                                    symelement_elec(sym_elec_count(general_counter_3),1,general_counter_3) = nucelem(general_counter_1)
                                                end do
                                            end do
                                        else
                                            print *, err, err_inpnn, "global_symfunction_electrostatic type ", words(2), " needs 3 arguments"; stop
                                        end if

                                    case ('5')
                                        if (nwords == 3) then
                                            read(words(3),*, iostat=ios) eta_temp
                                            if (ios /= 0) stop err // err_inpnn // "global_symfunction_electrostatic type ", words(2), " argument 2 must be a number"
                                            do general_counter_3 = 1,nelem
                                                sym_elec_count(general_counter_3) = sym_elec_count(general_counter_3) + 1
                                                function_type_elec(sym_elec_count(general_counter_3),general_counter_3) = function_type_temp
                                                eta_elec(sym_elec_count(general_counter_3),general_counter_3) = eta_temp
                                                symelement_elec(sym_elec_count(general_counter_3),1,general_counter_3) = nucelem(general_counter_3)
                                            end do
                                        else
                                            print *, err, err_inpnn, "global_symfunction_electrostatic type ", words(2), " needs 2 arguments"; stop
                                        end if

                                    case ('6')
                                        if (nwords == 3) then
                                            read(words(3),*, iostat=ios) funccutoff_temp
                                            if (ios /= 0) stop err // err_inpnn // "global_symfunction_electrostatic type ", words(2), " argument 2 must be a number"
                                            do general_counter_3 = 1,nelem
                                                sym_elec_count(general_counter_3) = sym_elec_count(general_counter_3) + 1
                                                function_type_elec(sym_elec_count(general_counter_3),general_counter_3) = function_type_temp
                                                funccutoff_elec(sym_elec_count(general_counter_3),general_counter_3) = funccutoff_temp
                                                symelement_elec(sym_elec_count(general_counter_3),1,general_counter_3) = nucelem(general_counter_3)
                                            end do
                                        else
                                            print *, err, err_inpnn, "global_symfunction_electrostatic type ", words(2), " needs 2 arguments"; stop
                                        end if

                                    case ('8')
                                        if (nwords == 5) then
                                            read(words(3),*, iostat=ios) eta_temp
                                            if (ios /= 0) stop err // err_inpnn // "global_symfunction_electrostatic type ", words(2), " argument 2 must be a number"
                                            read(words(4),*, iostat=ios) rshift_temp
                                            if (ios /= 0) stop err // err_inpnn // "global_symfunction_electrostatic type ", words(2), " argument 3 must be a number"
                                            read(words(5),*, iostat=ios) funccutoff_temp
                                            if (ios /= 0) stop err // err_inpnn // "global_symfunction_electrostatic type ", words(2), " argument 4 must be a number"
                                            do general_counter_3 = 1,nelem
                                                do general_counter_1 = 1,nelem
                                                    sym_elec_count(general_counter_3) = sym_elec_count(general_counter_3) + 1
                                                    function_type_elec(sym_elec_count(general_counter_3),general_counter_3) = function_type_temp
                                                    eta_elec(sym_elec_count(general_counter_3),general_counter_3) = eta_temp
                                                    rshift_elec(sym_elec_count(general_counter_3),general_counter_3) = rshift_temp
                                                    funccutoff_elec(sym_elec_count(general_counter_3),general_counter_3) = funccutoff_temp
                                                    symelement_elec(sym_elec_count(general_counter_3),1,general_counter_3) = nucelem(general_counter_1)
                                                    symelement_elec(sym_elec_count(general_counter_3),2,general_counter_3) = nucelem(general_counter_1)
                                                end do
                                                do general_counter_1 = 1,nelem
                                                    if (nelem .gt. 1) then
                                                        do general_counter_2 = 1,general_counter_1-1
                                                            sym_elec_count(general_counter_3) = sym_elec_count(general_counter_3) + 1
                                                            function_type_elec(sym_elec_count(general_counter_3),general_counter_3) = function_type_temp
                                                            eta_elec(sym_elec_count(general_counter_3),general_counter_3) = eta_temp
                                                            rshift_elec(sym_elec_count(general_counter_3),general_counter_3) = rshift_temp
                                                            funccutoff_elec(sym_elec_count(general_counter_3),general_counter_3) = funccutoff_temp
                                                            symelement_elec(sym_elec_count(general_counter_3),1,general_counter_3) = nucelem(general_counter_2)
                                                            symelement_elec(sym_elec_count(general_counter_3),2,general_counter_3) = nucelem(general_counter_1)
                                                        end do
                                                    end if
                                                end do
                                            end do
                                        else
                                            print *, err, err_inpnn, "global_symfunction_electrostatic type ", words(2), " needs 4 arguments"; stop
                                        end if

                                    case ('9')
                                        if (nwords == 6) then
                                            read(words(3),*, iostat=ios) eta_temp
                                            if (ios /= 0) stop err // err_inpnn // "global_symfunction_electrostatic type ", words(2), " argument 2 must be a number"
                                            read(words(4),*, iostat=ios) lambda_temp
                                            if (ios /= 0) stop err // err_inpnn // "global_symfunction_electrostatic type ", words(2), " argument 3 must be a number"
                                            read(words(5),*, iostat=ios) zeta_temp
                                            if (ios /= 0) stop err // err_inpnn // "global_symfunction_electrostatic type ", words(2), " argument 4 must be a number"
                                            read(words(6),*, iostat=ios) funccutoff_temp
                                            if (ios /= 0) stop err // err_inpnn // "global_symfunction_electrostatic type ", words(2), " argument 5 must be a number"
                                            do general_counter_3 = 1,nelem
                                                do general_counter_1 = 1,nelem
                                                    sym_elec_count(general_counter_3) = sym_elec_count(general_counter_3) + 1
                                                    function_type_elec(sym_elec_count(general_counter_3),general_counter_3) = function_type_temp
                                                    eta_elec(sym_elec_count(general_counter_3),general_counter_3) = eta_temp
                                                    lambda_elec(sym_elec_count(general_counter_3),general_counter_3) = lambda_temp
                                                    zeta_elec(sym_elec_count(general_counter_3),general_counter_3) = zeta_temp
                                                    funccutoff_elec(sym_elec_count(general_counter_3),general_counter_3) = funccutoff_temp
                                                    symelement_elec(sym_elec_count(general_counter_3),1,general_counter_3) = nucelem(general_counter_1)
                                                    symelement_elec(sym_elec_count(general_counter_3),2,general_counter_3) = nucelem(general_counter_1)
                                                end do
                                                do general_counter_1 = 1,nelem
                                                    if (nelem .gt. 1) then
                                                        do general_counter_2 = 1,general_counter_1-1
                                                            sym_elec_count(general_counter_3) = sym_elec_count(general_counter_3) + 1
                                                            function_type_elec(sym_elec_count(general_counter_3),general_counter_3) = function_type_temp
                                                            eta_elec(sym_elec_count(general_counter_3),general_counter_3) = eta_temp
                                                            lambda_elec(sym_elec_count(general_counter_3),general_counter_3) = lambda_temp
                                                            zeta_elec(sym_elec_count(general_counter_3),general_counter_3) = zeta_temp
                                                            funccutoff_elec(sym_elec_count(general_counter_3),general_counter_3) = funccutoff_temp
                                                            symelement_elec(sym_elec_count(general_counter_3),1,general_counter_3) = nucelem(general_counter_2)
                                                            symelement_elec(sym_elec_count(general_counter_3),2,general_counter_3) = nucelem(general_counter_1)
                                                        end do
                                                    end if
                                                end do
                                            end do
                                        else
                                            print *, err, err_inpnn, "global_symfunction_electrostatic type ", words(2), " needs 5 arguments"; stop
                                        end if

                                    case default
                                        print *, err, err_inpnn, "Error in global_symfunction_electrostatic key, symfunction type ", words(2), " not implemented"
                                        stop

                                end select

                            end if

                        case ('global_pairsymfunction_short')
                            print *, err, err_inpnn, "global_pairsymfunction_short key not supported, Pair NN not implemented"; stop

                        case ('element_pairsymfunction_short')
                            print *, err, err_inpnn, "element_pairsymfunction_short key not supported, Pair NN not implemented"; stop

                        case ('pairsymfunction_short')
                            print *, err, err_inpnn, "pairsymfunction_short key not supported, Pair NN not implemented"; stop

                        case default
                            ! just let it pass

                    end select

                else
                    write(*,*) err // err_inpnn // 'iostat = ', ios
                    stop
                end if
            end do

            close(inpnn_unit)

            do i1=1,nelem
                if(lshort.and.(nn_type_short.eq.1))then
                    num_funcvalues_short_atomic(i1)=sym_short_atomic_count(i1)
                    nodes_short_atomic(0,i1)=num_funcvalues_short_atomic(i1)
                endif
                if(lelec.and.(nn_type_elec.eq.1))then
                    num_funcvalues_elec(i1)=sym_elec_count(i1)
                    nodes_elec(0,i1)=num_funcvalues_elec(i1)
                endif
            enddo

            if(lshort.and.(nn_type_short.eq.1))then
                do i1=1,nelem
                    if(num_funcvalues_short_atomic(i1).eq.0)then
                        print *, 'ERROR: No short range symfunctions specified for ',element(i1)
                        stop
                    endif
                enddo
            endif
            if(lelec.and.(nn_type_elec.eq.1))then
                do i1=1,nelem
                    if(num_funcvalues_elec(i1).eq.0)then
                        print *, 'ERROR: No electrostatic symfunctions specified for ',element(i1)
                        stop
                    endif
                enddo
            endif

            ! call set_runner_counters() ! to avoid unwanted error messages in checkinputnn() subroutine

            call checkinputnn(err, err_inpnn) ! own subroutine in pes_nene_mod_supply.f90

            !call printinputnn(iseed,ielem,& ! should be skipped completely, but ask Jorg/Sascha if anything might be useful
            !    nodes_short_atomic_temp,nodes_elec_temp,nodes_short_pair_temp,&
            !    kalmanlambda_local,kalmanlambdae_local,&
            !    actfunc_short_atomic_dummy,actfunc_elec_dummy,actfunc_short_pair_dummy)

            write(ounit,'(a15,i4,a30)')' Element pairs: ',npairs,' , shortest distance (Bohr)'
            icount=0
            do i=1,nelem
                do j=i,nelem
                    icount=icount+1
                    if(dmin_element(icount).lt.9999.d0)then
                        write(ounit,'(a6,i4,2a3,1x,f10.3)')' pair ',&
                        icount,element(i),element(j),dmin_element(icount)
                    endif
                enddo
            enddo
            print *, '============================================================='

            if(lshort.and.(nn_type_short.eq.1))then
                do i1=1,nelem
                    wcount=0
                    do i=1,num_layers_short_atomic(i1)
                        wcount=wcount+1
                        windex_short_atomic(wcount,i1)=num_weights_short_atomic(i1)+1
                        num_weights_short_atomic(i1)=num_weights_short_atomic(i1)&
                          +nodes_short_atomic(i-1,i1)*nodes_short_atomic(i,i1)
                        wcount=wcount+1
                        windex_short_atomic(wcount,i1)=num_weights_short_atomic(i1)+1
                        num_weights_short_atomic(i1)=num_weights_short_atomic(i1)&
                          +nodes_short_atomic(i,i1) ! bias weights
                    enddo
                    if((mode.eq.2).or.(mode.eq.3))then
                        write(ounit,'(a,a3,i10)')' => short range NN weights type 1                ',&
                        element(i1),num_weights_short_atomic(i1)
                    endif
                    maxnum_weights_short_atomic=max(maxnum_weights_short_atomic,num_weights_short_atomic(i1))
                enddo
            endif

            if(lelec.and.(nn_type_elec.eq.1))then
                do i1=1,nelem
                    wcount=0
                    do i=1,num_layers_elec(i1)
                        wcount=wcount+1
                        windex_elec(wcount,i1)=num_weights_elec(i1)+1
                        num_weights_elec(i1)=num_weights_elec(i1)+nodes_elec(i-1,i1)*nodes_elec(i,i1)
                        wcount=wcount+1
                        windex_elec(wcount,i1)=num_weights_elec(i1)+1
                        num_weights_elec(i1)=num_weights_elec(i1)+nodes_elec(i,i1)
                    enddo
                    write(ounit,'(a,a3,i10)')' => electrostatic NN weights                     ',element(i1),num_weights_elec(i1)
                    maxnum_weights_elec=max(maxnum_weights_elec,num_weights_elec(i1))
                enddo
            endif
            write(ounit,*)'-------------------------------------------------------------'

            if(nn_type_short.eq.1)then
                maxnum_weights_short_pair=1
            endif
            if((.not.lelec).or.(lelec.and.(nn_type_elec.ne.1)))then
                maxnum_weights_elec=1
            endif

            if(lremoveatomenergies)then
                !call readatomenergies()

                ! start readout of input.nn according to readatomenergies.f90
                call open_for_read(inpnn_unit, filename_inpnn); ios = 0

                do while (ios == 0)
                    read(inpnn_unit, '(A)', iostat=ios) buffer
                    if (ios == 0) then
                        line = line + 1
                        call split_string(buffer, words, nwords)
                        atom_energy_counter = 0

                        select case (words(1))

                            case ('atom_energy')
                                !if ( /= default_int) stop err // err_inpnn // 'Multiple use of the  key'
                                if (nwords == 3) then
                                    read(words(2),'(A)', iostat=ios) elementtemp
                                    call nuccharge(elementtemp,ztemp)
                                    do general_counter_1 = 1,nelem
                                        if (ztemp.eq.nucelem(general_counter_1)) then
                                            atom_energy_counter = atom_energy_counter + 1
                                            read(words(2),'(A)', iostat=ios) elementtemp(atom_energy_counter)
                                            read(words(3),*, iostat=ios) atomrefenergies(atom_energy_counter)
                                            if (ios /= 0) stop err // err_inpnn // "atom_energy key in line ", line, " second argument value must be a number"
                                        else
                                            print *, warn_inpnn, 'atom_energy for element ',elementtemp,' is ignored'
                                        end if
                                    end do
                                else
                                    print *, err, err_inpnn, "atom_energy key in line ", line, " needs 2 arguments"; stop
                                end if

                            case default
                                ! just let it pass

                        end select

                    else
                        write(*,*) err // err_inpnn // 'iostat = ', ios
                        stop
                    end if
                end do

                close(inpnn_unit)

                do general_counter_1 = 1,atom_energy_counter
                    call nuccharge(elementsymbol(atom_energy_counter),zelem(atom_energy_counter))
                end do

                if(nelem.gt.1)then
                    do atom_energy_counter = 1,nelem-1
                        if (zelem(atom_energy_counter) .gt. zelem(atom_energy_counter+1)) then
                            ztemp=zelem(atom_energy_counter)
                            elementtemp=elementsymbol(atom_energy_counter)
                            etemp=atomrefenergies(atom_energy_counter)
                            zelem(atom_energy_counter)=zelem(atom_energy_counter+1)
                            elementsymbol(atom_energy_counter)=elementsymbol(atom_energy_counter+1)
                            atomrefenergies(atom_energy_counter)=atomrefenergies(atom_energy_counter+1)
                            zelem(atom_energy_counter+1)=ztemp
                            elementsymbol(atom_energy_counter+1)=elementtemp
                            atomrefenergies(atom_energy_counter+1)=etemp
                        endif
                    enddo
                endif

                lfound(:)=.false.
                do atom_energy_counter=1,nelem
                    lfound(elementindex(zelem(atom_energy_counter)))=.true.
                enddo

                do atom_energy_counter=1,nelem
                    if (lfound(atom_energy_counter) .eqv. .false.) then
                        print *, err, err_inpnn, 'Error: atom_energy not found for element ', nucelem(atom_energy_counter)
                        stop
                    endif
                enddo

                print *, 'atomic reference energies read from input.nn:'

                do atom_energy_counter=1,nelem
                    write(ounit,'(a1,a2,x,f18.8)')' ',elementsymbol(atom_energy_counter),atomrefenergies(atom_energy_counter)
                enddo
                ! end readout of input.nn according to readatomenergies.f90

            endif

            call open_for_read(inpnn_unit, filename_inpnn); ios = 0

            do while (ios == 0)
                read(inpnn_unit, '(A)', iostat=ios) buffer
                if (ios == 0) then
                    call split_string(buffer, words, nwords)

                    select case (words(1))

                        case ('node_activation_short')
                            !if ( lnode_activation_short /= default_bool) stop err // err_inpnn // 'Multiple use of the node_activation_short key'
                            !lnode_activation_short = .true.
                            print *, err, err_inpnn, "node_activation_short key was found, read activation functions of individual nodes is not implemented"; stop

                        case default
                            ! just let it pass

                    end select

                else
                    write(*,*) err // err_inpnn // 'iostat = ', ios
                    stop
                end if
            end do

            close(inpnn_unit)

            if(lshort.and.(nn_type_short.eq.1).and.(mode.ne.1))then
                do i3=1,nelem
                    write(ounit,*)'-------------------------------------------------'
                    write(ounit,*)'Atomic short range NN for element: ',element(i3)
                    write(ounit,'(a,10i5)')' architecture    ',(nodes_short_atomic(i1,i3),i1=0,num_layers_short_atomic(i3))
                    write(ounit,*)'-------------------------------------------------'
                    itemp=0
                    do i1=0,num_layers_short_atomic(i3)
                        itemp=max(itemp,nodes_short_atomic(i1,i3))
                    enddo ! i1
                    do i1=1,itemp ! loop over all lines with hidden nodes
                        if(i1.le.nodes_short_atomic(0,i3))then ! still input node to be printed
                            if(i1.le.maxnodes_short_atomic)then ! still hidden nodes present
                                write(ounit,'(i4,x,9a3)')i1,'  G',(actfunc_short_atomic(i1,i2,i3),i2=1,num_layers_short_atomic(i3))
                            else
                                write(ounit,'(i4,x,a3)')i1,'  G'
                            endif
                        else ! no input node in front of hidden nodes
                            write(ounit,'(i4,4x,8a3)')i1,(actfunc_short_atomic(i1,i2,i3),i2=1,num_layers_short_atomic(i3))
                        endif
                    enddo
                enddo ! i3
            endif

            if(lelec.and.(nn_type_elec.eq.1).and.(mode.ne.1))then
                do i3=1,nelem
                    write(ounit,*)'---------------------------------------------------'
                    write(ounit,*)'Electrostatic NN for element: ',element(i3)
                    write(ounit,'(a,10i5)')' architecture    ',(nodes_elec(i1,i3),i1=0,num_layers_elec(i3))
                    write(ounit,*)'---------------------------------------------------'
                    itemp=0
                    do i1=0,num_layers_elec(i3)
                        itemp=max(itemp,nodes_elec(i1,i3))
                    enddo ! i1
                    do i1=1,itemp ! loop over all lines with hidden nodes
                        if(i1.le.nodes_elec(0,i3))then ! still input node to be printed
                            if(i1.le.maxnodes_elec)then ! still hidden nodes present
                                write(ounit,'(i4,x,9a3)')i1,'  G',(actfunc_elec(i1,i2,i3),i2=1,num_layers_elec(i3))
                            else
                                write(ounit,'(i4,x,a3)')i1,'  G'
                            endif
                        else ! no input node in front of hidden nodes
                            write(ounit,'(i4,4x,8a3)')i1,(actfunc_elec(i1,i2,i3),i2=1,num_layers_elec(i3))
                        endif
                    enddo
                enddo ! i3
            endif
            write(ounit,*)'-------------------------------------------------------------'

            if((nn_type_short.eq.1).and.lshort)then
                call sortsymfunctions(&
                  maxnum_funcvalues_short_atomic,num_funcvalues_short_atomic,function_type_short_atomic,symelement_short_atomic,eta_short_atomic,zeta_short_atomic,rshift_short_atomic,lambda_short_atomic,funccutoff_short_atomic)
            endif

            if(lelec.and.(nn_type_elec.eq.1))then
                call sortsymfunctions(&
                  maxnum_funcvalues_elec,num_funcvalues_elec,function_type_elec,symelement_elec,eta_elec,zeta_elec,rshift_elec,lambda_elec,funccutoff_elec)
            endif

            if(lshort.and.(nn_type_short.eq.1))then
        do i1=1,nelem
          write(ounit,*)'-------------------------------------------------------------'
          write(ounit,*)' short range atomic symmetry &
                          &functions element ',element(i1),' :'
          write(ounit,*)'-------------------------------------------------------------'
          do i2=1,num_funcvalues_short_atomic(i1)
            if(function_type_short_atomic(i2,i1).eq.1)then
              write(ounit,'(i5,a3,i3,x,a3,3x,24x,f8.3)')&
                i2,element(i1),function_type_short_atomic(i2,i1),&
                element(elementindex(symelement_short_atomic(i2,1,i1))),&
                funccutoff_short_atomic(i2,i1)
            elseif(function_type_short_atomic(i2,i1).eq.2)then
              write(ounit,'(i5,a3,i3,x,a3,3x,8x,3f8.3)')&
                i2,element(i1),function_type_short_atomic(i2,i1),&
                element(elementindex(symelement_short_atomic(i2,1,i1))),&
                eta_short_atomic(i2,i1),rshift_short_atomic(i2,i1),funccutoff_short_atomic(i2,i1)
            elseif(function_type_short_atomic(i2,i1).eq.3)then
              write(ounit,'(i5,a3,i3,x,2a3,4f8.3)')&
                i2,element(i1),function_type_short_atomic(i2,i1),&
                element(elementindex(symelement_short_atomic(i2,1,i1))),&
                element(elementindex(symelement_short_atomic(i2,2,i1))),&
                eta_short_atomic(i2,i1),lambda_short_atomic(i2,i1),&
                zeta_short_atomic(i2,i1),funccutoff_short_atomic(i2,i1)
            elseif(function_type_short_atomic(i2,i1).eq.4)then
              write(ounit,'(i5,a3,i3,x,a3,3x,16x,2f8.3)')&
                i2,element(i1),function_type_short_atomic(i2,i1),&
                element(elementindex(symelement_short_atomic(i2,1,i1))),&
                eta_short_atomic(i2,i1),funccutoff_short_atomic(i2,i1)
            elseif(function_type_short_atomic(i2,i1).eq.5)then
              write(ounit,'(i5,a3,i3,4x,27x,f8.3)')&
                i2,element(i1),function_type_short_atomic(i2,i1),eta_short_atomic(i2,i1)
            elseif(function_type_short_atomic(i2,i1).eq.6)then
              write(ounit,'(i5,a3,i3,x,a3,3x,24x,f8.3)')&
                i2,element(i1),function_type_short_atomic(i2,i1),&
                element(elementindex(symelement_short_atomic(i2,1,i1))),&
                funccutoff_short_atomic(i2,i1)
            elseif(function_type_short_atomic(i2,i1).eq.8)then
              write(ounit,'(i5,a3,i3,x,2a3,4f8.3)')&
                i2,element(i1),function_type_short_atomic(i2,i1),&
                element(elementindex(symelement_short_atomic(i2,1,i1))),&
                element(elementindex(symelement_short_atomic(i2,2,i1))),&
                eta_short_atomic(i2,i1),rshift_short_atomic(i2,i1),&
                funccutoff_short_atomic(i2,i1)
            elseif(function_type_short_atomic(i2,i1).eq.9)then
              write(ounit,'(i5,a3,i3,x,2a3,4f8.3)')&
                i2,element(i1),function_type_short_atomic(i2,i1),&
                element(elementindex(symelement_short_atomic(i2,1,i1))),&
                element(elementindex(symelement_short_atomic(i2,2,i1))),&
                eta_short_atomic(i2,i1),lambda_short_atomic(i2,i1),&
                zeta_short_atomic(i2,i1),funccutoff_short_atomic(i2,i1)
            else
              write(ounit,*)'Error: printing unknown symfunction in readinput '
              stop
            endif
          enddo ! i2
        enddo ! i1=1,nelem
      endif ! lshort

          if(lelec.and.(nn_type_elec.eq.1))then
        do i1=1,nelem
          write(ounit,*)'-------------------------------------------------------------'
          write(ounit,*)' electrostatic symmetry functions element ',element(i1),' :'
          write(ounit,*)'-------------------------------------------------------------'
          do i2=1,num_funcvalues_elec(i1)
            if(function_type_elec(i2,i1).eq.1)then
              write(ounit,'(i5,a3,i3,x,a3,3x,24x,f8.3)')&
                i2,element(i1),function_type_elec(i2,i1),&
                element(elementindex(symelement_elec(i2,1,i1))),&
                funccutoff_elec(i2,i1)
            elseif(function_type_elec(i2,i1).eq.2)then
              write(ounit,'(i5,a3,i3,x,a3,3x,8x,3f8.3)')&
                i2,element(i1),function_type_elec(i2,i1),&
                element(elementindex(symelement_elec(i2,1,i1))),&
                eta_elec(i2,i1),rshift_elec(i2,i1),funccutoff_elec(i2,i1)
            elseif(function_type_elec(i2,i1).eq.3)then
              write(ounit,'(i5,a3,i3,x,2a3,4f8.3)')&
                i2,element(i1),function_type_elec(i2,i1),&
                element(elementindex(symelement_elec(i2,1,i1))),&
                element(elementindex(symelement_elec(i2,2,i1))),&
                eta_elec(i2,i1),lambda_elec(i2,i1),&
                zeta_elec(i2,i1),funccutoff_elec(i2,i1)
            elseif(function_type_elec(i2,i1).eq.4)then
              write(ounit,'(i5,a3,i3,x,a3,3x,16x,2f8.3)')&
                i2,element(i1),function_type_elec(i2,i1),&
                element(elementindex(symelement_elec(i2,1,i1))),&
                eta_elec(i2,i1),funccutoff_elec(i2,i1)
            elseif(function_type_elec(i2,i1).eq.5)then
              write(ounit,'(i5,a3,i3,4x,27x,f8.3)')&
                i2,element(i1),function_type_elec(i2,i1),eta_elec(i2,i1)
            elseif(function_type_elec(i2,i1).eq.6)then
              write(ounit,'(i5,a3,i3,x,a3,3x,24x,f8.3)')&
                i2,element(i1),function_type_elec(i2,i1),&
                element(elementindex(symelement_elec(i2,1,i1))),&
                funccutoff_elec(i2,i1)
            elseif(function_type_elec(i2,i1).eq.8)then
              write(ounit,'(i5,a3,i3,x,2a3,4f8.3)')&
                i2,element(i1),function_type_elec(i2,i1),&
                element(elementindex(symelement_elec(i2,1,i1))),&
                element(elementindex(symelement_elec(i2,2,i1))),&
                eta_elec(i2,i1),rshift_elec(i2,i1),&
                funccutoff_elec(i2,i1)
            elseif(function_type_elec(i2,i1).eq.9)then
              write(ounit,'(i5,a3,i3,x,2a3,4f8.3)')&
                i2,element(i1),function_type_elec(i2,i1),&
                element(elementindex(symelement_elec(i2,1,i1))),&
                element(elementindex(symelement_elec(i2,2,i1))),&
                eta_elec(i2,i1),lambda_elec(i2,i1),&
                zeta_elec(i2,i1),funccutoff_elec(i2,i1)
            else
              write(ounit,*)'Error: printing unknown symfunctione in readinput '
              stop
            endif
          enddo ! i2
        enddo ! i1=1,nelem
      endif ! lelec
      write(ounit,*)'-------------------------------------------------------------'
      ! end of readout according to readinput.f90








        ! further readout according to initnn.f90
        call getlistdim()

        !call distribute_predictionoptions() only mpi

        !call distribute_symfunctions() in symfunctions.f90, only mpi

        !call distribute_globaloptions() only mpi

        if(rinpparam%lshort.and.(rinpparam%nn_type_short.eq.1))then
            allocate (rinpparam%weights_short_atomic(rinpparam%maxnum_weights_short_atomic,rinpparam%nelem))
            rinpparam%weights_short_atomic(:,:)=0.0d0
            allocate (rinpparam%symfunction_short_atomic_list(rinpparam%maxnum_funcvalues_short_atomic,rinpparam%max_num_atoms,rinpparam%nblock))
            rinpparam%symfunction_short_atomic_list(:,:,:)=0.0d0
        end if

        if(rinpparam%lelec.and.(rinpparam%nn_type_elec.eq.1))then
            allocate (rinpparam%weights_elec(rinpparam%maxnum_weights_elec,rinpparam%nelem))
            rinpparam%weights_elec(:,:)=0.0d0
            allocate (rinpparam%symfunction_elec_list(rinpparam%maxnum_funcvalues_elec,rinpparam%max_num_atoms,rinpparam%nblock))
            rinpparam%symfunction_elec_list(:,:,:)=0.0d0
        end if
        ! end of readout according to initnn.f90, all things have been read and set up, ready for compute_nene()!!






        ! the dummy readout
        call open_for_read(inpnn_unit, filename_inpnn); ios = 0

        do while (ios == 0)
            read(inpnn_unit, '(A)', iostat=ios) buffer
            if (ios == 0) then
                line = line + 1
                call split_string(buffer, words, nwords)

                select case (words(1))

                    case ('')
                        if ( /= default_int) stop err // err_inpnn // 'Multiple use of the  key'
                        if (nwords == 2) then
                            read(words(2),'(i1000)', iostat=ios)
                            if (ios /= 0) stop err // err_inpnn // " value must be integer"
                        else
                            print *, err, err_inpnn, " key needs a single argument"; stop
                        end if

                    case ('')
                        if ( /= default_real) stop err // err_inpnn // 'Multiple use of the  key'
                        if (nwords == 2) then
                            read(words(2),*, iostat=ios)
                            if (ios /= 0) stop err // err_inpnn // " value must be a number"
                        else
                            print *, err, err_inpnn, " key needs a single argument"; stop
                        end if

                    case ('')
                        if ( /= default_bool) stop err // err_inpnn // 'Multiple use of the  key'
                        if (nwords == 1) then
                             = .true.
                        else
                            print *, err, err_inpnn, " key needs no argument(s)"; stop
                        end if

                    case ('')
                        if ( /= default_string) stop err // err_inpnn // 'Multiple use of the  key'
                        if (nwords == 2) then
                            read(words(2),'(A)', iostat=ios)
                        else
                            print *, err, err_inpnn, " key needs a single argument"; stop
                        end if

                    case ('')
                            print *, err, err_inpnn, " key is obsolete, please remove it"; stop

                    case default
                        if (trim(words(1)) /= '' .and. words(1)(1:1) /= '#') & ! check for empty and comment lines
                            print *, warn_inpnn, 'Skipping invalid label ', trim(words(1)),' in line ', line

                end select

            else
                write(*,*) err // err_inpnn // 'iostat = ', ios
                stop
            end if
        end do

        close(inpnn_unit)





        ! here the full list of keywords, remove after implementing according to readkeywords.f90!!

                    case ('fix_weights')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the fix_weights key'

                    case ('fixed_charge')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the fixed_charge key'

                    case ('fixed_short_energy_error_threshold')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the fixed_short_energy_error_threshold key'

                    case ('fixed_short_force_error_threshold')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the fixed_short_force_error_threshold key'

                    case ('force_grouping_by_structure')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the force_grouping_by_structure key'

                    case ('force_threshold')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the force_threshold key'

                    case ('force_update_scaling')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the force_update_scaling key'

                    case ('global_activation_electrostatic')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the global_activation_electrostatic key'

                    case ('global_activation_pair')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the global_activation_pair key'

                    case ('global_activation_short')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the global_activation_short key'

                    case ('global_hidden_layers_electrostatic')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the global_hidden_layers_electrostatic key'

                    case ('global_hidden_layers_pair')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the global_hidden_layers_pair key'

                    case ('global_hidden_layers_short')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the global_hidden_layers_short key'

                    case ('global_nodes_electrostatic')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the global_nodes_electrostatic key'

                    case ('global_nodes_pair' .or. 'global_nodes_short_pair')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the global_nodes_pair/global_nodes_short_pair key'

                    case ('global_nodes_short' .or. 'global_nodes_short_atomic')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the global_nodes_short/global_nodes_short_atomic key'

                    case ('global_pairsymfunction_short' .or. 'global_symfunction_short_pair')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the global_pairsymfunction_short/global_symfunction_short_pair key'

                    case ('global_symfunction_electrostatic' .or. 'global_symfunction_elec')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the global_symfunction_electrostatic/global_symfunction_elec key'

                    case ('global_symfunction_short' .or. 'global_symfunction_short_atomic')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the global_symfunction_short/global_symfunction_short_atomic key'

                    case ('growth_mode')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the growth_mode key'

                    case ('initialization_only')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the initialization_only key'

                    case ('ion_forces_only')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the ion_forces_only key'

                    case ('joint_energy_force_update')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the joint_energy_force_update key'

                    case ('kalman_damp_charge')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the kalman_damp_charge key'

                    case ('kalman_damp_force')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the kalman_damp_force key'

                    case ('kalman_damp_short')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the kalman_damp_short key'

                    case ('kalman_epsilon')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the kalman_epsilon key'

                    case ('kalman_lambda_charge')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the kalman_lambda_charge key'

                    case ('kalman_lambda_charge_constraint')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the kalman_lambda_charge_constraint key'

                    case ('kalman_lambda_short')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the kalman_lambda_short key'

                    case ('kalman_nue_charge')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the kalman_nue_charge key'

                    case ('kalman_nue_charge_constraint')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the kalman_nue_charge_constraint key'

                    case ('kalman_nue_short')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the kalman_nue_short key'

                    case ('kalman_q0')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the kalman_q0 key'

                    case ('kalman_qtau')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the kalman_qtau key'

                    case ('kalman_qmin')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the kalman_qmin key'

                    case ('max_energy')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the max_energy key'

                    case ('max_force')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the max_force key'

                    case ('md_mode')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the md_mode key'

                    case ('mix_all_points')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the mix_all_points key'

                    case ('nguyen_widrow_weights_ewald')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the nguyen_widrow_weights_ewald key'

                    case ('nguyen_widrow_weights_short')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the nguyen_widrow_weights_short key'

                    case ('noise_charge')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the noise_charge key'

                    case ('noise_energy')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the noise_energy key'

                    case ('noise_force')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the noise_force key'

                    case ('normalize_nodes')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the normalize_nodes key'

                    case ('optmode_charge')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the optmode_charge key'

                    case ('optmode_short_energy')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the optmode_short_energy key'

                    case ('optmode_short_force')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the optmode_short_force key'

                    case ('pairsymfunction_short')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the pairsymfunction_short key'

                    case ('points_in_memory' .or. 'nblock')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the points_in_memory/nblock key'

                    case ('precondition_weights')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the precondition_weights key'

                    case ('prepare_md')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the prepare_md key'

                    case ('print_date_and_time')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the print_date_and_time key'

                    case ('print_force_components')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the print_force_components key'

                    case ('read_kalman_matrices')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the read_kalman_matrices key'

                    case ('read_unformatted')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the read_unformatted key'

                    case ('remove_atom_energies')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the remove_atom_energies key'

                    case ('repeated_energy_update')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the repeated_energy_update key'

                    case ('reset_kalman')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the reset_kalman key'

                    case ('restrict_weights')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the restrict_weights key'

                    case ('save_kalman_matrices')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the save_kalman_matrices key'

                    case ('scale_max_elec')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the scale_max_elec key'

                    case ('scale_max_short')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the scale_max_short key'

                    case ('scale_max_short_atomic')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the scale_max_short_atomic key'

                    case ('scale_max_short_pair')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the scale_max_short_pair key'

                    case ('scale_min_elec')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the scale_min_elec key'

                    case ('scale_min_short')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the scale_min_short key'

                    case ('scale_min_short_atomic')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the scale_min_short_atomic key'

                    case ('scale_min_short_pair')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the scale_min_short_pair key'

                    case ('scale_symmetry_functions')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the scale_symmetry_functions key'

                    case ('screen_electrostatics')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the screen_electrostatics key'

                    case ('separate_bias_ini_short')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the separate_bias_ini_short key'

                    case ('separate_kalman_short')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the separate_kalman_short key'

                    case ('short_energy_error_threshold')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the short_energy_error_threshold key'

                    case ('short_energy_fraction')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the short_energy_fraction key'

                    case ('short_energy_group')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the short_energy_group key'

                    case ('short_force_error_threshold')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the short_force_error_threshold key'

                    case ('short_force_fraction')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the short_force_fraction key'

                    case ('short_force_group')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the short_force_group key'

                    case ('shuffle_weights_short_atomic')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the shuffle_weights_short_atomic key'

                    case ('silent_mode')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the silent_mode key'

                    case ('steepest_descent_step_charge')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the steepest_descent_step_charge key'

                    case ('steepest_descent_step_energy_short')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the steepest_descent_step_energy_short key'

                    case ('steepest_descent_step_force_short')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the steepest_descent_step_force_short key'

                    case ('symfunction_correlation')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the symfunction_correlation key'

                    case ('symfunction_electrostatic')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the symfunction_electrostatic key'

                    case ('symfunction_short')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the symfunction_short key'

                    case ('update_single_element')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the update_single_element key'

                    case ('update_worst_charges')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the update_worst_charges key'

                    case ('update_worst_short_energies')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the update_worst_short_energies key'

                    case ('update_worst_short_forces')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the update_worst_short_forces key'

                    case ('use_atom_charges')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the use_atom_charges key'

                    case ('use_atom_energies')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the use_atom_energies key'

                    case ('use_charge_constraint')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the use_charge_constraint key'

                    case ('use_damping')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the use_damping key'

                    case ('use_electrostatic_nn')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the use_electrostatic_nn key'

                    case ('use_noisematrix')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the use_noisematrix key'

                    case ('use_old_scaling')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the use_old_scaling key'

                    case ('use_old_weights_charge')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the use_old_weights_charge key'

                    case ('use_old_weights_short')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the use_old_weights_short key'

                    case ('use_short_forces')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the use_short_forces key'

                    case ('use_systematic_weights_electrostatic')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the use_systematic_weights_electrostatic key'

                    case ('use_systematic_weights_short')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the use_systematic_weights_short key'



        ! check existance of scaling.data
        if (.not. file_exists(filename_scaling)) stop err // err_scaling // 'file does not exist'

        ! read in all data from scaling.data
        call open_for_read(scaling_unit, filename_scaling); ios = 0

        do while (ios == 0) ! ios loop
            read(scaling_unit, '(A)', iostat=ios) buffer
            if (ios == 0) then

                ! readscale.f90
                ! do i1=1,ndim
                !     do i2=1,num_funcvalues_local(i1)
                !         read(scaleunit,*)i3,i3,minvalue_local(i1,i2),maxvalue_local(i1,i2),avvalue_local(i1,i2)
                !     enddo ! i2
                ! enddo ! i1

                ! read(scaleunit,*)eshortmin,eshortmax

                        ! initmode3.f90
                        ! call readscale(nelem,1,&  --> ndim = nelem!
!                           maxnum_funcvalues_short_atomic,num_funcvalues_short_atomic,&
!                           minvalue_short_atomic,maxvalue_short_atomic,avvalue_short_atomic,&
!                           eshortmin,eshortmax,rdummy,rdummy)





            else
                    write(*,*) err // err_scaling // 'iostat = ', ios
                    stop
            end if
        end do ! ios loop

        close(scaling_unit)


        ! loop over weight.XXX.data files and read in all data
        do weight_counter = 1,atoms%ntypes ! weights loop

            ! check existance of each weight file before reading
            if (.not. file_exists(weight_names_list(weight_counter))) stop err // err_weights // weight_names_list(weight_counter), ' file does not exist!'

            call open_for_read(weight_unit, weight_names_list(weight_counter)); ios = 0

            do while (ios == 0) ! ios loop
                read(weight_unit, '(A)', iostat=ios) buffer
                if (ios == 0) then

                        ! readweights.f90
                        ! do i1=1,ndim
                        !     do i2=1,num_weights_local(i1)
                        !         read(wunit,*)weights_local(i2,i1)
                        !     end do
                        ! end do

      !initmode3.f90
      !if(mpirank.eq.0)then
      !   if(lshort.and.(nn_type_short.eq.1))then
      !     call readweights(0,nelem,&
      !       maxnum_weights_short_atomic,num_weights_short_atomic,&
      !       weights_short_atomic)
      !   endif ! lshort
      !endif ! mpirank.eq.0





                else

                        write(*,*) err // err_weight // weight_names_list(weight_counter),', iostat = ', ios
                        stop

                end if
            end do ! ios loop

            close(weight_unit)

        end do ! weights loop

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

!       from predict.f90
!       allocate(sens(nelem,maxnum_funcvalues_short_atomic))


!       from initmode3.f90
!      maxcutoff_short_atomic       =0.0d0
!! get maxcutoff_short_atomic
!     if(lshort.and.(nn_type_short.eq.1))then
!       do i2=1,nelem
!        do i1=1,num_funcvalues_short_atomic(i2)
!          maxcutoff_short_atomic=max(maxcutoff_short_atomic,funccutoff_short_atomic(i1,i2))
!        enddo ! i1
!      enddo
!     endif ! lshort

!     minvalue_short_atomic(:,:)   =0.0d0
!     maxvalue_short_atomic(:,:)   =0.0d0
!     avvalue_short_atomic(:,:)    =0.0d0
!     minvalue_short_pair(:,:)     =0.0d0
!     maxvalue_short_pair(:,:)     =0.0d0
!     avvalue_short_pair(:,:)      =0.0d0
!     minvalue_elec(:,:) =0.0d0
!     maxvalue_elec(:,:) =0.0d0
!     avvalue_elec(:,:)  =0.0d0
!     chargemin(:)       =0.0d0
!     chargemax(:)       =0.0d0
!     if(mpirank.eq.0)then
!       if(lshort.and.(nn_type_short.eq.1))then
!         call readscale(nelem,1,&
!           maxnum_funcvalues_short_atomic,num_funcvalues_short_atomic,&
!           minvalue_short_atomic,maxvalue_short_atomic,avvalue_short_atomic,&
!           eshortmin,eshortmax,rdummy,rdummy)
!       endif ! lshort
!     endif ! mpirank.eq.0

!     if(lshort.and.(nn_type_short.eq.1))then
!       call
!     mpi_bcast(minvalue_short_atomic,nelem*maxnum_funcvalues_short_atomic,mpi_real8,0,mpi_comm_world,mpierror)
!       call
!     mpi_bcast(maxvalue_short_atomic,nelem*maxnum_funcvalues_short_atomic,mpi_real8,0,mpi_comm_world,mpierror)
!       call
!     mpi_bcast(avvalue_short_atomic,nelem*maxnum_funcvalues_short_atomic,mpi_real8,0,mpi_comm_world,mpierror)
!       call mpi_bcast(eshortmin,1,mpi_real8,0,mpi_comm_world,mpierror)
!       call mpi_bcast(eshortmax,1,mpi_real8,0,mpi_comm_world,mpierror)
!     endif

    end subroutine read_nene


    subroutine compute_nene(atoms, flag)

        ! Calculates energy and forces with HDNNPs

        ! md_tian2 related modules
        use constants, only : habohr2evang, timestep_ha2ev => ha2ev ! think about better way to realize that here we need our more precise ha2ev variable!!; conflict due to ha2ev from nnconstants.f90 and constants.f90

        ! RuNNer related modules (needed for and in predictionshortatomic.f90)
        use fileunits
        use globaloptions
        use mpi_mod
        use nnflags
        use nnshort_atomic
        use predictionoptions
        use saturation
        use structures
        use symfunctions
        use timings

        ! RuNNer related modules (needed for and in initnn.f90)
        !use fittingoptions
        !use mode1options
        !use nnewald
        !use nnconstants


        type(universe), intent(inout)   :: atoms
        integer, intent(in)             :: flag

        character(len=*), parameter :: err = "Error in compute_nene: "


        ! the elements have to be sorted according to RuNNer before calling the prediction -> better way than calling sortelements in every MD step -> do that when reading the structure file in md_tian2!!

        ! according to predict.f90

        ! move allocation/deallocation to read_nene() and before closing program like with the cleanup subroutine??
        if(lshort.and.(nn_type_short.eq.1))then
          allocate(sens(nelem,maxnum_funcvalues_short_atomic))
        !elseif(lshort.and.(nn_type_short.eq.2))then
        !  allocate(sens(npairs,maxnum_funcvalues_short_pair))
        endif

        if(lelec.and.(nn_type_elec.eq.1).or.(nn_type_elec.eq.3).or.(nn_type_elec.eq.4))then
          allocate(sense(nelem,maxnum_funcvalues_elec))
        endif

        !call getstructure_mode3(i4,num_atoms,num_pairs,zelem,&
        !  num_atoms_element,lattice,xyzstruct,&
        !  totalenergy,totalcharge,totalforce,atomenergy,atomcharge,&
        !  elementsymbol,lperiodic)

        ! from getstructure_mode3.f90
        if(mpirank.eq.0)then
            if(npoints.eq.1)then
                open(dataunit,file='input.data',form='formatted',status='old')
                rewind(dataunit)
            endif
            call readonestructure(num_atoms,&
                zelem,num_atoms_element,lattice,&
                totalcharge,totalenergy,atomcharge,atomenergy,xyzstruct,&
                totalforce,elementsymbol,lperiodic)
            if(npoints.eq.totnum_structures)then
                close(dataunit)
            endif
        endif

        !call initmode3(i4,&
        !  minvalue_short_atomic,maxvalue_short_atomic,avvalue_short_atomic,&
        !  minvalue_short_pair,maxvalue_short_pair,avvalue_short_pair,&
        !  minvalue_elec,maxvalue_elec,avvalue_elec,&
        !  eshortmin,eshortmax,chargemin,chargemax)

        if(lshort .and. nn_type_short == 1) then
          !if(nn_type_short.eq.1)then
            call predictionshortatomic(&
              num_atoms,num_atoms_element,zelem,&
              lattice,xyzstruct,&
              minvalue_short_atomic,maxvalue_short_atomic,avvalue_short_atomic,&
              eshortmin,eshortmax,&
              nntotalenergy,nnshortforce,&
              nnatomenergy,nnshortenergy,nnstress_short,&
              atomenergysum,sens,lperiodic)
          !elseif(nn_type_short.eq.2)then
            !call predictionshortpair(&
              !num_atoms,num_atoms_element,zelem,&
              !lattice,xyzstruct,&
              !minvalue_short_pair,maxvalue_short_pair,avvalue_short_pair,&
              !eshortmin,eshortmax,&
              !nntotalenergy,nnshortforce,&
              !nnatomenergy,nnpairenergy,nnshortenergy,&
              !nnstress_short,pairs_charge,&
              !atomenergysum,sens,lperiodic)
          !endif
        endif

        if(lelec .and. ((nn_type_elec == 1) .or. (nn_type_elec == 3) .or. (nn_type_elec == 4))) then ! will probably be implemented later, leave in comments or make dummy subroutine!
          call predictionelectrostatic(&
            num_atoms,zelem,&
            minvalue_elec,maxvalue_elec,avvalue_elec,&
            lattice,xyzstruct,&
            nntotalcharge,nnatomcharge,&
            chargemin,chargemax,nnelecenergy,&
            nnelecforce,nnstress_elec,sense,lperiodic)
        else
          nnatomcharge(:)=0.0d0
        endif

        ! combine short range and electrostatic energies
        nntotalenergy=nnshortenergy+nnelecenergy

        ! add energies of free atoms
        if(lremoveatomenergies.and.lshort)then
          call addatoms(num_atoms,&
            zelem,num_atoms_element,&
            atomenergysum,nnatomenergy)
          nntotalenergy=nntotalenergy+atomenergysum
        endif


        ! combination of short-range and electrostatic forces
        if(ldoforces)then
          nntotalforce(:,:)=nnshortforce(:,:)+nnelecforce(:,:)
        endif

        ! calculate the volume, needed also for stress
        if(lperiodic)then
          volume=0.0d0
          call getvolume(lattice,volume)
          if((mpirank.eq.0).and.(.not.lmd))then
            write(ounit,*)'-------------------------------------------------------------'
            write(ounit,*)'volume ',volume,' Bohr^3 for configuration ', i4
          endif
        endif

        ! combination of short-range and electrostatic stress
        if(ldostress.and.lperiodic)then
          nnstress(:,:)=nnstress_short(:,:)+nnstress_elec(:,:)
          nnstress(:,:)=nnstress(:,:)/volume
        endif

        ! check sum of forces if requested -> do we need that? seems fine, but clarify
        if(lcheckf)then
            forcesum(:)=0.0d0
            do i3=1,num_atoms
              do i2=1,3
                forcesum(i2)=forcesum(i2)+nntotalforce(i2,i3)
              enddo ! i2
            enddo ! i3
            write(ounit,'(A10,3A25)')'Conf.','Sum of Fx(Ha/Bohr)', 'Sum of Fy(Ha/Bohr)','Sum of Fz(Ha/Bohr)'
            write(ounit,'(I10,3f25.8)')i1,forcesum(1),forcesum(2),forcesum(3)
            do i2=1,3
              if(abs(forcesum(i2)).gt.0.000001d0)then
                write(ounit,'(I10,A31,I10,f25.8)')i4,'Error in forces of component: ',&
                  i2,forcesum(i2)
                stop
              endif
            enddo ! i2
          endif

        if(lshort.and.(nn_type_short.eq.1))then
            deallocate(sens)
        elseif(lshort.and.(nn_type_short.eq.2))then
            deallocate(sens)
        endif
        if(lelec.and.(nn_type_elec.eq.1).or.(nn_type_elec.eq.3).or.(nn_type_elec.eq.4))then
            deallocate(sense)
        endif

        ! in every md step we had to convert the positions from Angstrom to Bohr, look for the lattice (should only be converted once) and


        ! return the two following variables only
        atoms%epot = nntotalenergy * Ha2eV
        atoms%f(:,:,:) = nntotalforce * habohr2evang ! check if dimensions match

    end subroutine compute_nene

    ! according to main.f90
!   subroutine pes_nene_cleanup() ! in the end add in the source code if pes_id_nene or pes_name_nene => move to pes_nene_mod_supply.f90 and call it with use module, only: pes_nene_cleanup
!
!       use mpi_mod
!       use fileunits
!       use timings
!       use nnshort_atomic
!       use nnewald
!       use nnshort_pair
!       use symfunctions
!       use fittingoptions
!       use nnflags
!       use globaloptions
!
!       implicit none
!
!
!       ! according to main.f90
!       !call cleanup()
!       ! according to cleanup.f90
!       if(lshort.and.(nn_type_short.eq.1))then
!           deallocate(weights_short_atomic)
!           deallocate(symfunction_short_atomic_list)
!           deallocate(num_funcvalues_short_atomic)
!           deallocate(windex_short_atomic)
!           deallocate(num_layers_short_atomic)
!           deallocate(actfunc_short_atomic)
!           deallocate(nodes_short_atomic)
!           deallocate(num_weights_short_atomic)
!           deallocate(function_type_short_atomic)
!           deallocate(symelement_short_atomic)
!           deallocate(funccutoff_short_atomic)
!           deallocate(eta_short_atomic)
!           deallocate(zeta_short_atomic)
!           deallocate(lambda_short_atomic)
!           deallocate(rshift_short_atomic)
!       endif
!
!       if(lelec.and.(nn_type_elec.eq.1))then
!           deallocate(weights_elec)
!           deallocate(symfunction_elec_list)
!           deallocate(num_funcvalues_elec)
!           deallocate(windex_elec)
!           deallocate(num_layers_elec)
!           deallocate(actfunc_elec)
!           deallocate(nodes_elec)
!           deallocate(num_weights_elec)
!           deallocate(function_type_elec)
!           deallocate(symelement_elec)
!           deallocate(funccutoff_elec)
!           deallocate(eta_elec)
!           deallocate(zeta_elec)
!           deallocate(lambda_elec)
!           deallocate(rshift_elec)
!       endif
!
!       deallocate(nucelem)
!       deallocate(element)
!       deallocate(dmin_element)
!       if(allocated(atomrefenergies))deallocate(atomrefenergies)
!       if(allocated(fixedcharge))deallocate(fixedcharge)
!       if(allocated(elempair))deallocate(elempair)
!
!       call mpi_barrier(mpi_comm_world,mpierror)
!
!       ! according to main.f90
!       call mpi_finalize(mpierror)
!
!   end subroutine pes_nene_cleanup

end module pes_nene_mod
