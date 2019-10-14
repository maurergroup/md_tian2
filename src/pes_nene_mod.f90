!############################################################################
! This routine is part of
! md_tian2 (Molecular Dynamics Xia Tian 2)
! (c) 2014-2019 Dan J. Auerbach, Sascha Kandratsenka, Svenja M. Janke, Marvin
! Kammler, Sebastian Wille
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

    use constants
    use useful_things, only : split_string, lower_case, file_exists
    use universe_mod
    use open_file, only : open_for_read
!2do:
!   use RuNNer modules if necessary, otherwise make own ones
!   include RuNNer subroutine files -> are all subroutines completely independent or are they using global variables, if not ask Jorg to change that!!

    implicit none



    contains

    type runner_input_parameters

        ! input.nn
        ! 1) getdimensions.f90, 2) readinput.f90, 3) readkeywords.f90, 4) checkinputnn.f90
        ! getdimensions.f90
        integer :: nn_type_short
        integer :: mode
        logical :: lshort
        logical :: lelec
        integer :: nn_type_elec
        logical :: lfounddebug
        logical :: ldebug
        logical :: lfound_num_layersshort
        integer :: maxnum_layers_short_atomic
        logical :: lfound_num_layersewald
        integer :: maxnum_layers_elec
        logical :: lfound_num_layerspair
        integer :: maxnum_layers_short_pair
        logical :: lfound_luseatomenergies
        logical :: luseatomenergies
        logical :: lfound_luseatomcharges
        logical :: luseatomcharges
        logical :: lfound_nelem
        integer :: nelem
        integer :: npairs

        character(len=3), dimension(atoms%ntypes)                  :: element
        integer, dimension(atoms%ntypes)                           :: nucelem
        real(dp), dimension(atoms%ntypes * (atoms%ntypes + 1) / 2) :: dmin_element

        integer, dimension(:),   allocatable :: nodes_short_local
        integer, dimension(:),   allocatable :: nodes_ewald_local
        integer, dimension(:),   allocatable :: nodes_pair_local
        integer, dimension(:),   allocatable :: num_funcvalues_local
        integer, dimension(:),   allocatable :: num_funcvaluese_local
        integer, dimension(:,:), allocatable :: num_funcvaluesp_local

        character(len=3) :: elementtemp
        integer :: ztemp
        integer :: maxnum_funcvalues_short_atomic
        integer :: maxnum_funcvalues_elec
        integer :: maxnum_funcvalues_short_pair ! not needed?

        integer :: function_type_local
        integer :: function_type_temp
	real(dp) :: funccutoff_local
	real(dp) :: maxcutoff_local
	character(len=3) :: elementtemp1, elementtemp2, elementtemp3





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

    type(runner_input_parameters) :: rinpparam

    function new_runner_input_parameters

        type(runner_input_parameters) new_runner_input_parameters

        new_runner_input_parameters%nn_type_short                  = default_int ! 1 default, no pair available
        new_runner_input_parameters%mode                           = default_int ! 3 default, only RuNNer mode 3 implemented
        new_runner_input_parameters%lshort                         = default_bool ! .true. default, only short implemented
        new_runner_input_parameters%lelec                          = default_bool
        new_runner_input_parameters%nn_type_elec                   = default_int
        new_runner_input_parameters%lfounddebug                    = default_bool
        new_runner_input_parameters%ldebug                         = default_bool
        new_runner_input_parameters%lfound_num_layersshort         = default_bool
        new_runner_input_parameters%maxnum_layers_short_atomic     = default_int
        new_runner_input_parameters%lfound_num_layersewald         = default_bool
        new_runner_input_parameters%maxnum_layers_elec             = default_int
        new_runner_input_parameters%lfound_num_layerspair          = default_bool
        new_runner_input_parameters%maxnum_layers_short_pair       = default_int
        new_runner_input_parameters%lfound_luseatomenergies        = default_bool
        new_runner_input_parameters%luseatomenergies               = default_bool
        new_runner_input_parameters%lfound_luseatomcharges         = default_bool
        new_runner_input_parameters%luseatomcharges                = default_bool
        new_runner_input_parameters%lfound_nelem                   = default_bool
        new_runner_input_parameters%nelem                          = default_int
        new_runner_input_parameters%npairs                         = default_int
        new_runner_input_parameters%element                        = default_string
        new_runner_input_parameters%nucelem                        = default_int
        new_runner_input_parameters%dmin_element                   = default_real
        new_runner_input_parameters%nodes_short_local              = default_int
        new_runner_input_parameters%nodes_ewald_local              = default_int
        new_runner_input_parameters%nodes_pair_local               = default_int
        new_runner_input_parameters%num_funcvalues_local           = default_int
        new_runner_input_parameters%num_funcvaluese_local          = default_int
        new_runner_input_parameters%num_funcvaluesp_local          = default_int
        new_runner_input_parameters%elementtemp                    = default_string
        new_runner_input_parameters%ztemp                          = default_int
        new_runner_input_parameters%maxnum_funcvalues_short_atomic = default_int
        new_runner_input_parameters%maxnum_funcvalues_elec         = default_int
        new_runner_input_parameters%maxnum_funcvalues_short_pair   = default_int ! not needed?
        new_runner_input_parameters%function_type_local            = default_int
        new_runner_input_parameters%function_type_temp             = default_int
        new_runner_input_parameters%funccutoff_local               = default_real
	new_runner_input_parameters%maxcutoff_local                = default_real
	new_runner_input_parameters%elementtemp1                   = default_string
	new_runner_input_parameters%elementtemp2                   = default_string
	new_runner_input_parameters%elementtemp3                   = default_string

        new_runner_input_parameters%           = default_


        new_runner_input_parameters%weights_local                  = default_real



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

        character(len=max_string_length) :: filename_inpnn, filename_scaling
        character(len=max_string_length), dimension(atoms%ntypes) :: weights_path, filename_weights

        !logical, dimension(2) :: lshort

        integer  :: idx1, idx2, ntypes, weight_counter
        integer  :: npairs_counter_1, npairs_counter_2, element_counter, nodes_counter, general_counter_1, general_counter_2
        character(len=*), parameter :: err = "Error in read_nene: "
        character(len=*), parameter :: err_inpnn = "Error when reading input.nn: "
        character(len=*), parameter :: err_scaling = "Error when reading scaling.data: "
        character(len=*), parameter :: err_weight = "Error when reading the following weight file: "

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
            else if (nwords /= 2 .and. words(1) /= "weights") then
                print *,  err // "Error in the PES file: PES parameters must &
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

                    read(words(2), '(A)') inp_path

                case ('weights')

                    do weight_counter = 1,atoms%ntypes

                        read(words(weight_counter+1), '(A)') weights_path(weight_counter)

                    end do

                case default

                    print *, "Error in the PES file: unknown nene parameter ", words(1)
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


        ! check existance of input.nn
        if (.not. file_exists(filename_inpnn)) stop err // err_inpnn // "file does not exist"

        ! read all input keywords from input.nn several times to respect dependencies

!       read in keywords related to input.nn according to the following files from RuNNer (chronologically)
!       1) getdimensions.f90
!       2) paircount.f90; not all, but have a further look
!       3) readkeywords.f90
!       4) readinput.f90

!       following the dummy readout
!
!       call open_for_read(inpnn_unit, filename_inpnn); ios = 0
!
!       do while (ios == 0)
!           read(inpnn_unit, '(A)', iostat=ios) buffer
!           if (ios == 0) then
!               line = line + 1
!               call split_string(buffer, words, nwords)
!
!               select case (words(1))
!
!                   case ('')
!                       if (rinpparam% /= default_int) stop err // err_inpnn // 'Multiple use of the  key'
!                       if (nwords == 2) then
!                           read(words(2),'(i1000)', iostat=ios) rinpparam%
!                           if (ios /= 0) stop err // err_inpnn // " value must be integer"
!                       else
!                           print *, err, err_inpnn, " key needs a single argument"; stop
!                       end if
!
!                   case default
!                       if (trim(words(1)) /= '' .and. words(1)(1:1) /= '#') & ! check for empty and comment lines
!                           print *, 'Skipping invalid label', trim(words(1)),'in line', line
!
!               end select
!           else
!               write(*,*) err // err_inpnn // 'iostat = ', ios
!               stop
!           end if
!
!       close(inpnn_unit)

        call open_for_read(inpnn_unit, filename_inpnn); ios = 0

        do while (ios == 0)
            read(inpnn_unit, '(A)', iostat=ios) buffer
            if (ios == 0) then
                line = line + 1
                call split_string(buffer, words, nwords)

                select case (words(1))

                    case ('nn_type_short')
                        if (rinpparam%nn_type_short /= default_int) stop err // err_inpnn // 'Multiple use of the nn_type_short key'
                        if (nwords == 2) then
                            read(words(2),'(i1000)', iostat=ios) rinpparam%nn_type_short
                            if (ios /= 0) stop err // err_inpnn // "nn_type_short value must be integer"

                            select case (words(2))

                                case (1)
                                    !pass

                                case (2)
                                    print *, err // err_inpnn // "nn_type_short 2 not supported, Pair NN not implemented!"
                                    stop

                                case default
                                    print *, err // err_inpnn // "Error in nn_type_short key value, ", words(2), " not implemented"
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
                        if (trim(words(1)) /= '' .and. words(1)(1:1) /= '#') & ! check for empty and comment lines
                            print *, 'Skipping invalid label ', trim(words(1)),' in line ', line

                end select

            else
                write(*,*) err // err_inpnn // 'iostat = ', ios
                stop
            end if

        close(inpnn_unit)


        call open_for_read(inpnn_unit, filename_inpnn); ios = 0

        do while (ios == 0)
            read(inpnn_unit, '(A)', iostat=ios) buffer
            if (ios == 0) then
                line = line + 1
                call split_string(buffer, words, nwords)

                select case (words(1))

                    case ('debug_mode') ! not needed
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

                    case ('global_hidden_layers_pair') ! not needed - make dummy and include error message that pair type is not implemented
                        if (rinpparam%lfound_num_layerspair  /= default_bool) stop err // err_inpnn // 'Multiple use of the global_hidden_layers_pair key'
                        if (nwords == 2) then
                            rinpparam%lfound_num_layerspair = .true.
                            read(words(2),'(i1000)', iostat=ios) rinpparam%maxnum_layers_short_pair
                            if (ios /= 0) stop err // err_inpnn // "global_hidden_layers_pair value must be integer"
                            rinpparam%maxnum_layers_short_pair = rinpparam%maxnum_layers_short_pair + 1
                        else
                            print *, err, err_inpnn, "global_hidden_layers_pair key needs a single argument"; stop
                        end if

                    case ('use_atom_energies')
                        if (rinpparam%lfound_luseatomenergies /= default_bool) stop err // err_inpnn // 'Multiple use of the use_atom_energies key'
                        if (nwords == 1) then
                            rinpparam%lfound_luseatomenergies = .true.
                            rinpparam%luseatomenergies = .true.
                        else
                            print *, err, err_inpnn, "use_atom_energies key needs no argument(s)"; stop
                        end if

                    case ('use_atom_charges')
                        if (rinpparam%lfound_luseatomcharges /= default_int) stop err // err_inpnn // 'Multiple use of the use_atom_charges key'
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
                            if (rinpparam%nelem /= atoms%ntypes) stop err // err_inpnn // "number of elements in input.nn and structure file differ"
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
                        if (trim(words(1)) /= '' .and. words(1)(1:1) /= '#') & ! check for empty and comment lines
                            print *, 'Skipping invalid label ', trim(words(1)),' in line ', line

                end select

            else
                write(*,*) err // err_inpnn // 'iostat = ', ios
                stop
            end if

        close(inpnn_unit)

        if (rinpparam%lshort .and. (rinpparam%nn_type_short == 1) .and. (rinpparam%maxnum_layers_short_atomic == default_int)) stop err // err_inpnn // 'global_hidden_layers_short key not set'
        if (rinpparam%lshort .and. (rinpparam%nn_type_short == 2) .and. (rinpparam%maxnum_layers_short_pair == default_int)) stop err // err_inpnn // 'global_hidden_layers_pairs key not set'
        if (rinpparam%lelec .and. (rinpparam%nn_type_elec == 1) .and. (rinpparam%maxnum_layers_elec == default_int)) stop err // err_inpnn // 'global_hidden_layers_electrostatic key not set'

        !allocate(rinpparam%nucelem(rinpparam%nelem))
        !allocate(rinpparam%element(rinpparam%nelem))
        !allocate(rinpparam%dmin_element(rinpparam%nelem*(rinpparam%nelem+1)/2))


        call open_for_read(inpnn_unit, filename_inpnn); ios = 0 ! maybe move to readout before since nelem is given by atoms%ntypes

        do while (ios == 0)
            read(inpnn_unit, '(A)', iostat=ios) buffer
            if (ios == 0) then
                line = line + 1
                call split_string(buffer, words, nwords)

                select case (words(1))

                    case ('elements')
                        if (rinpparam%element /= default_string) stop err // err_inpnn // 'Multiple use of the elements key'
                        if (nwords == atoms%ntypes+1) then
                            do element_counter = 1,atoms%ntypes
                                read(words(element_counter+1),'(A)', iostat=ios) rinpparam%element(element_counter)
                                !if (ios /= 0) stop err // err_inpnn // "element symbol ", element_counter, " must be string not a number" -> maybe check for proper string without digits -> use white list (periodic table) and forget about black list
                            end do
                        else
                            print *, err, err_inpnn, "elements key does not match with number of element types"; stop
                        end if

                    case default
                        if (trim(words(1)) /= '' .and. words(1)(1:1) /= '#') & ! check for empty and comment lines
                            print *, 'Skipping invalid label ', trim(words(1)),' in line ', line

                end select

            else
                write(*,*) err // err_inpnn // 'iostat = ', ios
                stop
            end if

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
        if (rinpparam%maxnum_layers_short_pair == default_int) then ! not needed
            rinpparam%maxnum_layers_short_pair = 0
        end if
        if (rinpparam%lfound_nelem == default_bool) stop err // err_inpnn // "number of elements not found"

        allocate(rinpparam%num_funcvalues_local(102))
        allocate(rinpparam%num_funcvaluese_local(102))
        allocate(rinpparam%num_funcvaluesp_local(102,102)) ! not needed
        rinpparam%num_funcvalues_local(:) = 0
        rinpparam%num_funcvaluese_local(:) = 0
        rinpparam%num_funcvaluesp_local(:) = 0 ! not needed

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
                line = line + 1
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
                            print *, err, err_inpnn, "global_nodes_electrostatic argument number does not match with global_hidden_layers_electrostatic value"; stop
                        end if

                    case ('global_nodes_pair') ! not needed - make dummy
                        if (rinpparam%nodes_pair_local /= default_int) stop err // err_inpnn // 'Multiple use of the global_nodes_pair key'
                        if (nwords == rinpparam%maxnum_layers_short_pair+1) then
                            do nodes_counter = 1,rinpparam%maxnum_layers_short_pair-1
                                read(words(nodes_counter+1),'(i1000)', iostat=ios) rinpparam%nodes_pair_local(nodes_counter)
                                if (ios /= 0) stop err // err_inpnn // "global_nodes_pair value", nodes_counter, " must be integer"
                            end do
                        else
                            print *, err, err_inpnn, "global_nodes_pair argument number does not match with global_hidden_layers_pair value"; stop
                        end if

                    case ('element_symfunction_short')
                        !if (rinpparam% /= default_int) stop err // err_inpnn // 'Multiple use of the element_symfunction_short key'
                        !if (nwords == 2) then
                        read(words(2),'(A)', iostat=ios) rinpparam%elementtemp
                        read(words(3),'(i1000)', iostat=ios) rinpparam%function_type_local
                        if (ios /= 0) stop err // err_inpnn // "element_symfunction_short second argument value must be integer"
                        call nuccharge(rinpparam%elementtemp, rinpparam%ztemp)
                        if (rinpparam%num_funcvalues_local(rinpparam%ztemp) /= default_int) then
                            print *, err // err_inpnn // 'Error in element_symfunction_short: Element with atomic number', rinpparam%num_funcvalues_local(rinpparam%ztemp), 'already set, check for multiple use of key'
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

                            case (5,6) ! only for Pair NN
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

                            case (5,6) ! only for Pair NN
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

                            case (5,6) ! only for Pair NN
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

                    case default
                        if (trim(words(1)) /= '' .and. words(1)(1:1) /= '#') & ! check for empty and comment lines
                            print *, 'Skipping invalid label ', trim(words(1)),' in line ', line

                end select

            else
                write(*,*) err // err_inpnn // 'iostat = ', ios
                stop
            end if

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
        if (allocated(rinpparam%nodes_pair_local)) deallocate(rinpparam%nodes_pair_local) ! not needed

        do general_counter_1 = 1,102
            rinpparam%maxnum_funcvalues_short_atomic = max(rinpparam%maxnum_funcvalues_short_atomic, )
        end do

        deallocate(rinpparam%num_funcvalues_local)
        deallocate(rinpparam%num_funcvaluese_local)
        deallocate(rinpparam%num_funcvaluesp_local) ! not needed

        deallocate(rinpparam%nucelem)
        deallocate(rinpparam%element)
        !end readout according to getdimensions.f90


        !start readout according to paircount.f90
        call open_for_read(inpnn_unit, filename_inpnn); ios = 0

        do while (ios == 0)
            read(inpnn_unit, '(A)', iostat=ios) buffer
            if (ios == 0) then
                line = line + 1
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
                                    print *, err, err_inpnn, "global_symfunction_short type ", words(2), " needs ", nwords-1, " arguments"; stop
                                end if

                            case (2)
                                if (nwords == 5)
                                    read(words(5),*, iostat=ios) rinpparam%funccutoff_local
                                    if (ios /= 0) stop err // err_inpnn // "global_symfunction_short type ", words(2), " argument ", nwords-1, " must be a number"
                                else
                                    print *, err, err_inpnn, "global_symfunction_short type ", words(2), " needs ", nwords-1, " arguments"; stop
                                end if

                            case (3)
                                if (nwords == 6)
                                    read(words(6),*, iostat=ios) rinpparam%funccutoff_local
                                    if (ios /= 0) stop err // err_inpnn // "global_symfunction_short type ", words(2), " argument ", nwords-1, " must be a number"
                                else
                                    print *, err, err_inpnn, "global_symfunction_short type ", words(2), " needs ", nwords-1, " arguments"; stop
                                end if

                            case (4)
                                if (nwords == 4)
                                    read(words(4),*, iostat=ios) rinpparam%funccutoff_local
                                    if (ios /= 0) stop err // err_inpnn // "global_symfunction_short type ", words(2), " argument ", nwords-1, " must be a number"
                                else
                                    print *, err, err_inpnn, "global_symfunction_short type ", words(2), " needs ", nwords-1, " arguments"; stop
                                end if

                            case (8)
                                if (nwords == 5)
                                    read(words(5),*, iostat=ios) rinpparam%funccutoff_local
                                    if (ios /= 0) stop err // err_inpnn // "global_symfunction_short type ", words(2), " argument ", nwords-1, " must be a number"
                                else
                                    print *, err, err_inpnn, "global_symfunction_short type ", words(2), " needs ", nwords-1, " arguments"; stop
                                end if

                            case (9)
                                if (nwords == 6)
                                    read(words(6),*, iostat=ios) rinpparam%funccutoff_local
                                    if (ios /= 0) stop err // err_inpnn // "global_symfunction_short type ", words(2), " argument ", nwords-1, " must be a number"
                                else
                                    print *, err, err_inpnn, "global_symfunction_short type ", words(2), " needs ", nwords-1, " arguments"; stop
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
                                    print *, err, err_inpnn, "element_symfunction_short type ", words(3), " needs ", nwords-1, " arguments"; stop
                                end if

                            case (2)
                                if (nwords == 6)
                                    read(words(6),*, iostat=ios) rinpparam%funccutoff_local
                                    if (ios /= 0) stop err // err_inpnn // "element_symfunction_short type ", words(3), " argument ", nwords-1, " must be a number"
                                else
                                    print *, err, err_inpnn, "element_symfunction_short type ", words(3), " needs ", nwords-1, " arguments"; stop
                                end if

                            case (3)
                                if (nwords == 7)
                                    read(words(7),*, iostat=ios) rinpparam%funccutoff_local
                                    if (ios /= 0) stop err // err_inpnn // "element_symfunction_short type ", words(3), " argument ", nwords-1, " must be a number"
                                else
                                    print *, err, err_inpnn, "element_symfunction_short type ", words(3), " needs ", nwords-1, " arguments"; stop
                                end if

                            case (4)
                                if (nwords == 5)
                                    read(words(5),*, iostat=ios) rinpparam%funccutoff_local
                                    if (ios /= 0) stop err // err_inpnn // "element_symfunction_short type ", words(3), " argument ", nwords-1, " must be a number"
                                else
                                    print *, err, err_inpnn, "element_symfunction_short type ", words(3), " needs ", nwords-1, " arguments"; stop
                                end if

                            case (8)
                                if (nwords == 6)
                                    read(words(6),*, iostat=ios) rinpparam%funccutoff_local
                                    if (ios /= 0) stop err // err_inpnn // "element_symfunction_short type ", words(3), " argument ", nwords-1, " must be a number"
                                else
                                    print *, err, err_inpnn, "element_symfunction_short type ", words(3), " needs ", nwords-1, " arguments"; stop
                                end if

                            case (9)
                                if (nwords == 7)
                                    read(words(7),*, iostat=ios) rinpparam%funccutoff_local
                                    if (ios /= 0) stop err // err_inpnn // "element_symfunction_short type ", words(3), " argument ", nwords-1, " must be a number"
                                else
                                    print *, err, err_inpnn, "element_symfunction_short type ", words(3), " needs ", nwords-1, " arguments"; stop
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
                        select case (words(2))

                            case (1)
                                if (nwords == 5)
                                    read(words(4),'(A)') rinpparam%elementtemp2
                                    read(words(5),*, iostat=ios) rinpparam%funccutoff_local
                                    if (ios /= 0) stop err // err_inpnn // "symfunction_short type ", words(3), " argument ", nwords-1, " must be a number"
                                else
                                    print *, err, err_inpnn, "symfunction_short type ", words(3), " needs ", nwords-1, " arguments"; stop
                                end if

                            case (2)
                                if (nwords == 7)
                                    read(words(4),'(A)') rinpparam%elementtemp2
                                    read(words(7),*, iostat=ios) rinpparam%funccutoff_local
                                    if (ios /= 0) stop err // err_inpnn // "symfunction_short type ", words(3), " argument ", nwords-1, " must be a number"
                                else
                                    print *, err, err_inpnn, "symfunction_short type ", words(3), " needs ", nwords-1, " arguments"; stop
                                end if

                            case (3)
                                if (nwords == 9)
                                    read(words(4),'(A)') rinpparam%elementtemp2
                                    read(words(5),'(A)') rinpparam%elementtemp3
                                    read(words(9),*, iostat=ios) rinpparam%funccutoff_local
                                    if (ios /= 0) stop err // err_inpnn // "symfunction_short type ", words(3), " argument ", nwords-1, " must be a number"
                                else
                                    print *, err, err_inpnn, "symfunction_short type ", words(3), " needs ", nwords-1, " arguments"; stop
                                end if

                            case (4)
                                if (nwords == 7)
                                    read(words(4),'(A)') rinpparam%elementtemp2
                                    read(words(7),*, iostat=ios) rinpparam%funccutoff_local
                                    if (ios /= 0) stop err // err_inpnn // "symfunction_short type ", words(3), " argument ", nwords-1, " must be a number"
                                else
                                    print *, err, err_inpnn, "symfunction_short type ", words(3), " needs ", nwords-1, " arguments"; stop
                                end if

                            case (8)
                                if (nwords == 8)
                                    read(words(4),'(A)') rinpparam%elementtemp2
                                    read(words(5),'(A)') rinpparam%elementtemp3
                                    read(words(8),*, iostat=ios) rinpparam%funccutoff_local
                                    if (ios /= 0) stop err // err_inpnn // "symfunction_short type ", words(3), " argument ", nwords-1, " must be a number"
                                else
                                    print *, err, err_inpnn, "symfunction_short type ", words(3), " needs ", nwords-1, " arguments"; stop
                                end if

                            case (9)
                                if (nwords == 9)
                                    read(words(4),'(A)') rinpparam%elementtemp2
                                    read(words(5),'(A)') rinpparam%elementtemp3
                                    read(words(9),*, iostat=ios) rinpparam%funccutoff_local
                                    if (ios /= 0) stop err // err_inpnn // "symfunction_short type ", words(3), " argument ", nwords-1, " must be a number"
                                else
                                    print *, err, err_inpnn, "symfunction_short type ", words(3), " needs ", nwords-1, " arguments"; stop
                                end if

                            case default
                                print *, err, err_inpnn, "Error in symfunction_short key, symfunction type ", words(3), " not implemented"
                                stop

                        end select

                    case default
                        if (trim(words(1)) /= '' .and. words(1)(1:1) /= '#') & ! check for empty and comment lines
                            print *, 'Skipping invalid label ', trim(words(1)),' in line ', line

                end select

                rinpparam%maxcutoff_local = max(rinpparam%maxcutoff_local, rinpparam%funccutoff_local)

            else
                write(*,*) err // err_inpnn // 'iostat = ', ios
                stop
            end if

        close(inpnn_unit)

        if (rinpparam%maxcutoff_local == default_real) then
            print *, err, err_inpnn, "maxcutoff_local is zero, specify symmetry functions"
            stop
        end if



        call open_for_read(inpnn_unit, filename_inpnn); ios = 0

        do while (ios == 0)
            read(inpnn_unit, '(A)', iostat=ios) buffer
            if (ios == 0) then
                line = line + 1
                call split_string(buffer, words, nwords)

                select case (words(1))

                    case ('')
                        if (rinpparam% /= default_int) stop err // err_inpnn // 'Multiple use of the  key'
                        if (nwords == 2) then
                            read(words(2),'(i1000)', iostat=ios) rinpparam%
                            if (ios /= 0) stop err // err_inpnn // " value must be integer"
                        else
                            print *, err, err_inpnn, " key needs a single argument"; stop
                        end if

                    case default
                        if (trim(words(1)) /= '' .and. words(1)(1:1) /= '#') & ! check for empty and comment lines
                            print *, 'Skipping invalid label ', trim(words(1)),' in line ', line

                end select
            else
                write(*,*) err // err_inpnn // 'iostat = ', ios
                stop
            end if

        close(inpnn_unit)






        call open_for_read(inpnn_unit, filename_inpnn); ios = 0

        do while (ios == 0) ! analog to read_input_file subroutine in run_config.f90
            read(inpnn_unit, '(A)', iostat=ios) buffer
            if (ios == 0) then
                line = line + 1
                call split_string(buffer, words, nwords)

                select case (words(1))

                    case ('analyze_composition')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the analyze_composition key'

                    case ('analyze_error')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the analyze_error key'

                    case ('analyze_error_charge_step')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the analyze_error_charge_step key'

                    case ('analyze_error_energy_step')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the analyze_error_energy_step key'

                    case ('analyze_error_force_step')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the analyze_error_force_step key'

                    case ('atom_energy')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the atom_energy key'

                    case ('biasweights_max')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the biasweights_max key'

                    case ('biasweights_min')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the biasweights_min key'

                    case ('bond_threshold')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the bond_threshold key'

                    case ('calculate_final_force')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the calculate_final_force key'

                    case ('calculate_forces')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the calculate_forces key'

                    case ('calculate_stress')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the calculate_stress key'

                    case ('center_symmetry_functions')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the center_symmetry_functions key'

                    case ('charge_error_threshold')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the charge_error_threshold key'

                    case ('charge_fraction')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the charge_fraction key'

                    case ('charge_group')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the charge_group key'

                    case ('charge_grouping_by_structure')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the charge_grouping_by_structure key'

                    case ('charge_update_scaling')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the charge_update_scaling key'

                    case ('check_forces')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the check_forces key'

                    case ('check_input_forces')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the check_input_forces key'

                    case ('cutoff_type')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the cutoff_type key'

                    case ('data_clustering')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the data_clustering key'

                    case ('debug_mode')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the debug_mode key'

                    case ('detailed_timing')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the detailed_timing key'

                    case ('detailed_timing_epoch')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the detailed_timing_epoch key'

                    case ('detect_saturation')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the detect_saturation key'

                    case ('dynamic_force_grouping')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the dynamic_force_grouping key'

                    case ('electrostatic_type' .or. 'nn_type_elec')
                        if (rinpparam%nn_type_elec_local /= default_) stop err // err_inpnn // 'Multiple use of the electrostatic_type/nn_type_elec key'

                    case ('element_activation_electrostatic')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the element_activation_electrostatic key'

                    case ('element_activation_pair')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the element_activation_pair key'

                    case ('element_activation_short')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the element_activation_short key'

                    case ('element_decoupled_forces_v2')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the element_decoupled_forces_v2 key'

                    case ('element_decoupled_kalman')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the element_decoupled_kalman key'

                    case ('element_hidden_layers_electrostatic')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the element_hidden_layers_electrostatic key'

                    case ('element_hidden_layers_pair')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the element_hidden_layers_pair key'

                    case ('element_hidden_layers_short')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the element_hidden_layers_short key'

                    case ('element_nodes_electrostatic')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the element_nodes_electrostatic key'

                    case ('element_nodes_pair')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the element_nodes_pair key'

                    case ('element_nodes_short')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the element_nodes_short key'

                    case ('element_pairsymfunction_short')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the element_pairsymfunction_short key'
                        ! you shall pass

                    case ('element_symfunction_electrostatic')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the element_symfunction_electrostatic key'
                        ! you shall pass

                    case ('element_symfunction_short')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the element_symfunction_short key'
                        ! you shall pass

                    case ('elements')
                        if (any(rinpparam%element /= default_string)) stop err // err_inpnn // 'Multiple use of the elements key'
                        if (nwords == atoms%nelem+1)
                            do i=1,atoms%nelem
                                read(words(i+1),'(A)', iostat=ios) element(i)
                                if (ios /= 0) stop err // err_inpnn // "elements keyword values must be string"
                                if (any(rinpparam%element /= atoms%name)) stop err // err_inpnn // "element names in input.nn and *.inp/poscar files differ"
                        else
                            print *, err, err_inpnn, "Error: number of element symbols given does not match with number of elements from structure file"; stop
                        end if


                    case ('enable_on_the_fly_input')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the enable_on_the_fly_input key'

                    case ('energy_threshold')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the energy_threshold key'

                    case ('enforce_max_num_neighbors_atomic')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the enforce_max_num_neighbors_atomic key'

                    case ('enforce_totcharge')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the enforce_totcharge key'

                    case ('environment_analysis')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the environment_analysis key'

                    case ('epochs')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the epochs key'

                    case ('ewald_alpha')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the ewald_alpha key'

                    case ('ewald_cutoff')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the ewald_cutoff key'

                    case ('ewald_kmax')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the ewald_kmax key'

                    case ('find_contradictions')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the find_contradictions key'

                    case ('fitmode')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the fitmode key'

                    case ('fitting_unit')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the fitting_unit key'

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

                    case ('global_output_nodes_electrostatic')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the global_output_nodes_electrostatic key'
                        print *, err, err_inpnn, "Error: global_output_nodes_electrostatic keyword is obsolete, please remove it"; stop

                    case ('global_output_nodes_pair')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the global_output_nodes_pair key'
                        print *, err, err_inpnn, "Error: global_output_nodes_pair keyword is obsolete, please remove it"; stop

                    case ('global_output_nodes_short')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the global_output_nodes_short key'
                        print *, err, err_inpnn, "Error: global_output_nodes_short keyword is obsolete, please remove it"; stop

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

                    case ('nn_type')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the nn_type key'
                        print *, err, err_inpnn, "Error: nn_type keyword is obsolete, use nn_type_short instead"; stop

                    case ('nn_type_short') ! only short available
                        if (rinpparam%nn_type_short /= default_int) stop err // err_inpnn // 'Multiple use of the nn_type_short key'
                        if (nwords == 2) then
                            read(words(2),'(i1000)', iostat=ios) rinpparam%nn_type_short
                            if (ios /= 0) stop err // err_inpnn // "nn_type_short value must be integer"
                            if (words(2) /= 1) then
                                print *, err, err_inpnn, "Only nn_type_short 1 (Behler-Parrinello) available!"; stop
                            end if
                        else
                            print *, err, err_inpnn, "nn_type_short key needs a single argument"; stop
                        end if

                    case ('noise_charge')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the noise_charge key'

                    case ('noise_energy')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the noise_energy key'

                    case ('noise_force')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the noise_force key'

                    case ('normalize_nodes')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the normalize_nodes key'

                    case ('number_of_elements')
                        if (rinpparam%nelem /= default_int) stop err // err_inpnn // 'Multiple use of the number_of_elements key'
                        if (nwords == 2) then
                            read(words(2),'(i1000)',iostat=ios) rinpparam%nelem
                            if (ios /= 0) stop err // err_inpnn // "number_of_elements value must be integer"
                            if (rinpparam%nelem /= atoms%ntypes) stop err // err_inpnn // "element number in input.nn and *.inp/poscar files differ"
                        else
                            print *, err, err_inpnn, "number_of_elements key needs a single argument"
                        end if




                    case ('number_of_elements') ! check with md_tian input
                        if (rinpparam%nelem /= default_int) stop err // err_inpnn // 'Multiple use of the number_of_elements key'
                        if (nwords == 2) then
                            read(words(2),'(A)',iostat=ios) rinpparam%nelem

                            do nelem_counter_1=1,rinpparam%nelem
                                do nelem_counter_2=1,rinpparam%nelem
                                    rinpparam%npairs = rinpparam%npairs + 1
                                end do
                            end do
                        else
                            print *, err, err_inpnn, "number_of_elements key needs a single argument"; stop
                        end if



                    case ('optmode_charge')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the optmode_charge key'

                    case ('optmode_short_energy')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the optmode_short_energy key'

                    case ('optmode_short_force')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the optmode_short_force key'

                    case ('pairsymfunction_short')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the pairsymfunction_short key'

                    case ('parallel_mode')
                        if (rinpparam%paramode /= default_) stop err // err_inpnn // 'Multiple use of the parallel_mode key'
                        if (nwords == 2) then
                            read(words(2),'(A)') rinpparam%paramode
                            if (words(2) .neqv. 1) then
                                print *, err, err_inpnn, "Only parallel_mode 1 (serial) available"; stop
                        else
                            print *, err, err_inpnn, "parallel_mode key needs a single argument"; stop
                        end if

                    case ('points_in_memory' .or. 'nblock')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the points_in_memory/nblock key'

                    case ('precondition_weights')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the precondition_weights key'

                    case ('prepare_md')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the prepare_md key'

                    case ('print_all_deshortdw')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the print_all_deshortdw key'

                    case ('print_all_dfshortdw')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the print_all_dfshortdw key'

                    case ('print_all_electrostatic_weights')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the print_all_electrostatic_weights key'

                    case ('print_all_short_weights')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the print_all_short_weights key'

                    case ('print_convergence_vector')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the print_convergence_vector key'

                    case ('print_date_and_time')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the print_date_and_time key'

                    case ('print_force_components')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the print_force_components key'

                    case ('print_mad')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the print_mad key'

                    case ('print_sensitivity')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the print_sensitivity key'

                    case ('random_number_type')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the random_number_type key'

                    case ('random_order_training')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the random_order_training key'
                        print *, err, err_inpnn, "Error: random_order_training keyword is obsolete, please use mix_all_points
                        instead"; stop

                    case ('random_seed')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the random_seed key'

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

                    case ('runner_mode') ! only mode 3 available
                        if (rinpparam%mode /= default_int) stop err // err_inpnn // 'Multiple use of the runner_mode key'

                        if (nwords == 2) then
                            read(words(2),'') rinpparam%mode
                            if (words(2) /= 3)
                                print *, err, err_inpnn, "Only mode 3 available"; stop
                        else
                            print *, err, err_inpnn, "runner_mode key needs a single argument"; stop
                        end if

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

                    case ('test_fraction')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the test_fraction key'

                    case ('total_charge_error_threshold')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the total_charge_error_threshold key'

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

                    case ('use_electrostatics')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the use_electrostatics key'

                    case ('use_fixed_charges')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the use_fixed_charges key'

                    case ('use_ipi')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the use_ipi key'

                    case ('use_noisematrix')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the use_noisematrix key'

                    case ('use_old_scaling')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the use_old_scaling key'

                    case ('use_old_weights_charge')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the use_old_weights_charge key'

                    case ('use_old_weights_short')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the use_old_weights_short key'

                    case ('use_omp_mkl')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the use_omp_mkl key'

                    case ('use_short_forces')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the use_short_forces key'

                    case ('use_short_nn')
                        if (rinpparam%lshort_local /= default_bool) stop err // err_inpnn // 'Multiple use of the use_short_nn key'
                        if (nwords == 1) then
                            rinpparam%lshort_local = .true.
                        else
                            print *, err, err_inpnn, "use_short_nn key has additional arguments"; stop
                        end if

                    case ('use_systematic_weights_electrostatic')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the use_systematic_weights_electrostatic key'

                    case ('use_systematic_weights_short')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the use_systematic_weights_short key'

                    case ('weight_analysis')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the weight_analysis key'

                    case ('weight_constraint')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the weight_constraint key'

                    case ('weighte_constraint')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the weighte_constraint key'

                    case ('weights_max')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the weights_max key'

                    case ('weights_min')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the weights_min key'

                    case ('weightse_max')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the weightse_max key'

                    case ('weightse_min')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the weightse_min key'

                    case ('write_fit_statistics')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the write_fit_statistics key'

                    case ('write_pdb')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the write_pdb key'
                        print *, err, err_inpnn, "Error: write_pdb keyword is no longer supported"; stop

                    case ('write_pov')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the write_pov key'
                        print *, err, err_inpnn, "Error: write_pov keyword is no longer supported"; stop

                    case ('write_pwscf')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the write_pwscf key'
                        print *, err, err_inpnn, "Error: write_pwscf keyword is no longer supported"; stop

                    case ('write_temporary_weights')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the write_temporary_weights key'

                    case ('write_traincharges')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the write_traincharges key'

                    case ('write_trainforces')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the write_trainforces key'

                    case ('write_trainpoints')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the write_trainpoints key'

                    case ('write_unformatted')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the write_unformatted key'

                    case ('write_weights_epoch')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the write_weights_epoch key'

                    case ('write_xyz')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the write_xyz key'
                        print *, err, err_inpnn, "Error: write_xyz keyword is no longer supported"; stop




















                    case ('elements') ! check with md_tian.inp
                        if (rinpparam%element /= default_string) stop err // err_inpnn // 'Multiple use of the elements key'
                        if (nwords == atoms%ntypes+1) then
                            do element_counter = 1,atoms%ntypes
                                read(words(element_counter+1),'(A)') rinpparam%element(element_counter)
                            end do
                        else
                            print *, err, "elements key does not match with number of element types"; stop
                        end if





                    case ('global_hidden_layers_short')
                        if (rinpparam%global_hidden_layers_short /= default_int) stop err // err_inpnn // 'Multiple use of the global_hidden_layers_short key'
                        if (nwords == 2) then
                            rinpparam%global_hidden_layers_short = rinpparam%global_hidden_layers_short + 1
                        else
                            print *, err, err_inpnn, "global_hidden_layers_short key needs a single argument"; stop
                        end if



                    case ('global_output_nodes_short')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the  key'
                        print *, "Error: global_output_nodes_short keyword is obsolete, please remove it"
                        stop

                    case ('global_output_nodes_electrostatic')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the  key'
                        print *, "Error: global_output_nodes_electrostatic keyword is obsolete, please remove it"
                        stop

                    case ('global_output_nodes_pair')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the  key'
                        print *, "Error: global_output_nodes_pair keyword is obsolete, please remove it"
                        stop



                    case ('random_order_training')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the  key'
                        print *, "Error: random_order_training keyword is obsolete, please use mix_all_points instead"
                        stop







                    case ('nn_type')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the  key'
                        print *, "Error: nn_type keyword is obsolete, use nn_type_short instead"
                        stop






                    case ('use_fixed_charges')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the  key'
                        print *, "use_fixed_charges keyword is obsolete, use electrostatic_type 3 instead"
                        stop



                    case ('use_electrostatic_nn')
                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the  key'
                        print *, "Error: use_electrostatic_nn keyword is obsolete, please use electrostatic_type and use_electrostatics instead"
                        stop

!                    case ('')
!                        if (rinpparam% /= default_) stop err // err_inpnn // 'Multiple use of the  key'

                    case default
                        if (trim(words(1)) /= '' .and. words(1)(1:1) /= '#') & !check for empty and comment lines
                            print *, 'Skipping invalid label ',trim(words(1)),' in line ', line

                end select

            else

                write(*,*) err // err_inpnn // 'iostat = ', ios
                stop

            end if

        close(inpnn_unit)


!       from paircount.f90
!! get the cutoff for nn_type_short_local = 1
!     if(nn_type_short_local.eq.1)then
!       open(nnunit,file='input.nn')
!       rewind(nnunit)

!110     read(nnunit,*,END=100) keyword

!         if(keyword.eq.'global_symfunction_short')then
!           backspace(nnunit)
!           read(nnunit,*)dummy,function_type_temp
!           if(function_type_temp.eq.1)then
!             backspace(nnunit)
!             read(nnunit,*)dummy,function_type_temp,funccutoff_local
!           elseif(function_type_temp.eq.2)then
!             backspace(nnunit)
!             read(nnunit,*)dummy,function_type_temp,dummy,dummy,funccutoff_local
!           elseif(function_type_temp.eq.3)then
!             backspace(nnunit)
!             read(nnunit,*)dummy,function_type_temp,dummy,dummy,dummy,funccutoff_local
!           elseif(function_type_temp.eq.4)then
!             backspace(nnunit)
!             read(nnunit,*)dummy,function_type_temp,dummy,funccutoff_local
!           elseif(function_type_temp.eq.8)then
!             backspace(nnunit)
!             read(nnunit,*)dummy,function_type_temp,dummy,dummy,funccutoff_local
!           elseif(function_type_temp.eq.9)then
!             backspace(nnunit)
!             read(nnunit,*)dummy,function_type_temp,dummy,dummy,dummy,funccutoff_local
!           else
!             write(ounit,*)'Error: unknown symfunction in paircount'
!             stop
!           endif

!         elseif(keyword.eq.'element_symfunction_short')then
!           backspace(nnunit)
!           read(nnunit,*)dummy,elementtemp1,function_type_temp
!           if(function_type_temp.eq.1)then
!             backspace(nnunit)
!             read(nnunit,*)dummy,elementtemp1,function_type_temp,&
!               funccutoff_local
!           elseif(function_type_temp.eq.2)then
!             backspace(nnunit)
!             read(nnunit,*)dummy,elementtemp1,function_type_temp,&
!               dummy,dummy,funccutoff_local
!           elseif(function_type_temp.eq.3)then
!             backspace(nnunit)
!             read(nnunit,*)dummy,elementtemp1,function_type_temp,&
!               dummy,dummy,dummy,funccutoff_local
!           elseif(function_type_temp.eq.4)then
!             backspace(nnunit)
!             read(nnunit,*)dummy,elementtemp1,function_type_temp,&
!               dummy,funccutoff_local
!           elseif(function_type_temp.eq.8)then
!             backspace(nnunit)
!             read(nnunit,*)dummy,elementtemp1,function_type_temp,&
!               dummy,dummy,funccutoff_local
!           elseif(function_type_temp.eq.9)then
!             backspace(nnunit)
!             read(nnunit,*)dummy,elementtemp1,function_type_temp,&
!               dummy,dummy,dummy,funccutoff_local
!           else
!             write(ounit,*)'Error: unknown symfunction in paircount'
!             stop
!          endif
!!
!         elseif(keyword.eq.'symfunction_short')then
!           backspace(nnunit)
!           read(nnunit,*)dummy,elementtemp1,function_type_temp
!           if(function_type_temp.eq.1)then
!             backspace(nnunit)
!             read(nnunit,*)dummy,elementtemp1,function_type_temp,&
!               elementtemp2,funccutoff_local
!           elseif(function_type_temp.eq.2)then
!             backspace(nnunit)
!             read(nnunit,*)dummy,elementtemp1,function_type_temp,&
!               elementtemp2,dummy,dummy,funccutoff_local
!           elseif(function_type_temp.eq.3)then
!             backspace(nnunit)
!             read(nnunit,*)dummy,elementtemp1,function_type_temp,&
!               elementtemp2,elementtemp3,dummy,dummy,dummy,funccutoff_local
!           elseif(function_type_temp.eq.4)then
!             backspace(nnunit)
!             read(nnunit,*)dummy,elementtemp1,function_type_temp,&
!               elementtemp2,dummy,dummy,funccutoff_local
!           elseif(function_type_temp.eq.8)then
!             backspace(nnunit)
!             read(nnunit,*)dummy,elementtemp1,function_type_temp,&
!               elementtemp2,elementtemp3,dummy,dummy,funccutoff_local
!           elseif(function_type_temp.eq.9)then
!             backspace(nnunit)
!             read(nnunit,*)dummy,elementtemp1,function_type_temp,&
!               elementtemp2,elementtemp3,dummy,dummy,dummy,funccutoff_local
!           else
!             write(ounit,*)'Error: unknown symfunction in paircount'
!             stop
!          endif
!         endif
!!
!         maxcutoff_local=max(maxcutoff_local,funccutoff_local)
!         goto 110
!100       continue
!       close(nnunit)

!! check if maxcutoff_local has been found properly
!       if(maxcutoff_local.eq.0.0d0)then
!         write(ounit,*)'Error: maxcutoff_local is zero in paircount'
!         write(ounit,*)'Did you forget to specify symmetry functions???'
!         stop !'
!       endif
!     endif ! --> nn_type_short_local = 1

!??
!       count_struct=0
!       open(dataunit,file='input.data',form='formatted')
!         backspace(dataunit)
!10       continue
!         read(dataunit,*,END=30) keyword
!         if(keyword.eq.'begin') then
!           count_struct= count_struct+1
!           nlattice  = 0
!           num_atoms = 0
!           lperiodic =.false.
!         endif
!         if(keyword.eq.'lattice') then
!           nlattice=nlattice+1
!           backspace(dataunit)
!           read(dataunit,*)keyword,(lattice(nlattice,i1),i1=1,3)
!         endif
!         if(keyword.eq.'atom') then
!           backspace(dataunit)
!           num_atoms=num_atoms+1
!           read(dataunit,*)keyword,(xyzstruct(i1,num_atoms),i1=1,3),elementsymbol(num_atoms)
!           call nuccharge(elementsymbol(num_atoms),zelem(num_atoms))
!         endif
!         if(keyword.eq.'end') then
!           if(nlattice.eq.3)then
!             lperiodic=.true.
!           endif
!           if(lperiodic)then
!              call translate(num_atoms,lattice,xyzstruct)
!           endif
!! determine num_pairs
!           call getnumpairs(num_atoms,num_pairs,zelem,&
!             maxcutoff_local,lattice,xyzstruct,dmin_temp,lperiodic)
!           max_num_pairs=max(max_num_pairs,num_pairs)
!           do i1=1,nelem*(nelem+1)/2
!!              write(ounit,*)i1,dmin_element(i1),dmin_temp(i1)
!             dmin_element(i1)=min(dmin_element(i1),dmin_temp(i1))
!           enddo
!!
!         endif
!         goto 10
!30       continue
!       close(dataunit)
!??




!     if(nn_type_short_local.ne.2)then
!       max_num_pairs=0
!     endif


        ! check existance of scaling.data
        if (.not. file_exists(filename_scaling)) stop err // err_scaling // 'file does not exist'

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

                        ! initmode3.f90
                        ! call readscale(nelem,1,&  --> ndim = nelem!
!                           maxnum_funcvalues_short_atomic,num_funcvalues_short_atomic,&
!                           minvalue_short_atomic,maxvalue_short_atomic,avvalue_short_atomic,&
!                           eshortmin,eshortmax,rdummy,rdummy)





            else

                    write(*,*) err // err_scaling // 'iostat = ', ios
                    stop

            end if
        end do ! ios
        close(scaling_unit)


        ! loop over weight.XXX.data files and read in all data
        do weight_counter = 1,atoms%ntypes ! weights

            ! check existance of each weight file before reading
            if (.not. file_exists(weight_names_list(weight_counter))) stop err // err_weights // weight_names_list(weight_counter),
            ' file does not exist!'

            call open_for_read(weight_unit, weight_names_list(weight_counter)); ios = 0

            do while (ios == 0) ! ios
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

        type(universe), intent(inout)   :: atoms
        integer, intent(in)             :: flag

!       from predict.f90
!       call predictionshortatomic(&
!             num_atoms,num_atoms_element,zelem,&
!             lattice,xyzstruct,&
!             minvalue_short_atomic,maxvalue_short_atomic,avvalue_short_atomic,&
!             eshortmin,eshortmax,&
!             nntotalenergy,nnshortforce,&
!             nnatomenergy,nnshortenergy,nnstress_short,&
!             atomenergysum,sens,lperiodic)




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


    subroutine sortelements()

        type(runner_input_parameters), intent(inout)   :: rinpparam

        integer counter
        integer ztemp
        integer nuc_counter_1,nuc_counter_2

        character*2 elementtemp

        rinpparam%elementindex(:) = 0

        if(nelem.gt.1)then
            do nuc_counter_1 = 1,rinpparam%nelem - 1
                if (rinpparam%nucelem(nuc_counter_1) .gt. rinpparam%nucelem(nuc_counter_2)) then
                    ztemp = rinpparam%nucelem(nuc_counter_1)
                    elementtemp = rinpparam%element(nuc_counter_1)
                    rinpparam%nucelem(nuc_counter_1) = rinpparam%nucelem(nuc_counter_1 + 1)
                    rinpparam%element(nuc_counter_1) = rinpparam%element(nuc_counter_1 + 1)
                    rinpparam%nucelem(nuc_counter_1 + 1) = ztemp
                    rinpparam%element(nuc_counter_1 + 1) = elementtemp
                end if
            end do
        end if

        do nuc_counter_1 = 1,102
            do nuc_counter_2 = 1,rinpparam%nelem
                if (rinpparam%nucelem(nuc_counter_2) .eq. nuc_counter_1) then
                    rinpparam%elementindex(nuc_counter_1) = nuc_counter_2
                end if
            end do
        end do

        rinpparam%pairindex(:,:) = 0
        counter = 0
        do nuc_counter_1 = 1,rinpparam%nelem
            do nuc_counter_2 = 1,rinpparam%nelem
                if (rinpparam%nucelem(nuc_counter_2) .ge. rinpparam%nucelem(nuc_counter_1)) then
                    counter = cpunter + 1
                    rinpparam%pairindex(rinpparam%nucelem(nuc_counter_1),rinpparam%nucelem(nuc_counter_2)) = counter
                    rinpparam%pairindex(rinpparam%nucelem(nuc_counter_2),rinpparam%nucelem(nuc_counter_1)) = counter
                end if
            end do
        end do

    end subroutine sortelements

end module pes_nene_mod
