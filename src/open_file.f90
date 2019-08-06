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

! This module contains routines to open files smoothly
module open_file

contains

    subroutine open_for_read(lun,file_name)
        implicit none
        integer, intent(in)           :: lun
        character(len=*), intent(in)  :: file_name

        integer                       :: ios
        character(len=120)            :: error_message

        open(unit=lun, file=file_name, status='old', action='read', iostat=ios, iomsg=error_message)
        if (ios==0) return
        print '(/ "Error on open file ", (a), " for read i/o status=", i4 )', TRIM(file_name), ios
        print '( "error message=", (a) )', error_message

        STOP 101

    end subroutine open_for_read

    subroutine open_for_write(lun,file_name)
        implicit none
        integer, intent(in)           :: lun
        character(len=*), intent(in)  :: file_name

        integer                       :: ios
        character(120)                :: error_message

        open(unit=lun, file=file_name, status='new', action='write', iostat=ios, iomsg=error_message)
        if (ios==0) return

        open(unit=lun, file=file_name, status='replace', action='write', iostat=ios, iomsg=error_message)

        if (ios==0) return
        print *, 'failed to open file ', file_name, ' for write with status=replace. i/o status =',ios
        print *, 'error message: ', error_message
        STOP 103

    end subroutine open_for_write

    subroutine open_for_append(lun,file_name)
        implicit none
        integer, intent(in)           :: lun
        character(len=*), intent(in)  :: file_name

        integer                       :: ios
        character(120)                :: error_message

        open(unit=lun, file=file_name, status='new', action='write', iostat=ios, iomsg=error_message)
        if (ios==0) return

        open(unit=lun, file=file_name, status='old', access='append', action='write', &
            iostat=ios, iomsg=error_message)

        if (ios==0) return
        print *, 'failed to open file ', file_name, ' for write with status=old. i/o status =',ios
        print *, 'error message: ', error_message
        STOP 103

    end subroutine open_for_append

end module
