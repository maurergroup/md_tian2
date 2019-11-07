!############################################################################
! This routine is part of
! md_tian2 (Molecular Dynamics Tian Xia 2)
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

subroutine write_info()

      implicit none

      print *, '--------------------------------------------------------------------------'
      print *, '                md_tian2 (Molecular Dynamics Tian Xia 2)                  '
      print *, '---------  (c) 2014-2019 Dan J. Auerbach, Sascha Kandratsenka,   ---------'
      print *, '---------    Svenja M. Janke, Marvin Kammler, Sebastian Wille    ---------'
      print *, '---------            Dynamics at Surfaces Department             ---------'
      print *, '---------    MPI for Biophysical Chemistry Goettingen, Germany   ---------'
      print *, '---------      Georg-August-Universitaet Goettingen, Germany     ---------'
      print *, '--------------------------------------------------------------------------'
      print *, '--------------------------------------------------------------------------'
      print *, 'This program is free software: you can redistribute it and/or modify it   '
      print *, 'under the terms of the GNU General Public License as published by the     '
      print *, 'Free Software Foundation, either version 3 of the License, or             '
      print *, '(at your option) any later version.                                       '
      print *, '                                                                          '
      print *, 'This program is distributed in the hope that it will be useful, but       '
      print *, 'WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY'
      print *, 'or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License   '
      print *, 'for more details.                                                         '
      print *, '                                                                          '
      print *, 'You should have received a copy of the GNU General Public License along   '
      print *, 'with this program. If not, see http://www.gnu.org/licenses.               '
      print *, '--------------------------------------------------------------------------'
      print *, '--------------------------------------------------------------------------'
      print *, 'When using md_tian2, please cite the following papers:                    '

      !2do: maybe add further machine specific details, papers to cite and further things to mention
      !2do: add RuNNer/REBO etc. related papers also in here?

end subroutine write_info()
