!######################################################################
! This routine is part of
! RuNNer - RuNNer Neural Network Energy Representation
! (c) 2008-2020 Prof. Dr. Joerg Behler 
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
!######################################################################
!! called by:
!! - checkonestructure.f90
!! - getrefatoms.f90
!! - readinput.90
!! - readonestructure.f90
!!
      subroutine nuccharge(elementtemp,zelemtemp)
!!
      implicit none
!!
      integer zelemtemp           ! out
!!
      character*2 elementtemp     ! in
!!
        if(trim(elementtemp).eq.'H') then
          zelemtemp=1
        elseif(trim(elementtemp).eq.'He') then
          zelemtemp=2
        elseif(trim(elementtemp).eq.'Li') then
          zelemtemp=3
        elseif(trim(elementtemp).eq.'Be') then
          zelemtemp=4
        elseif(trim(elementtemp).eq.'B') then
          zelemtemp=5
        elseif(trim(elementtemp).eq.'C') then
          zelemtemp=6
        elseif(trim(elementtemp).eq.'N') then
          zelemtemp=7
        elseif(trim(elementtemp).eq.'O') then
          zelemtemp=8
        elseif(trim(elementtemp).eq.'F') then
          zelemtemp=9
        elseif(trim(elementtemp).eq.'Ne') then
          zelemtemp=10
        elseif(trim(elementtemp).eq.'Na') then
          zelemtemp=11
        elseif(trim(elementtemp).eq.'Mg') then
          zelemtemp=12
        elseif(trim(elementtemp).eq.'Al') then
          zelemtemp=13
        elseif(trim(elementtemp).eq.'Si') then
          zelemtemp=14
        elseif(trim(elementtemp).eq.'P') then
          zelemtemp=15
        elseif(trim(elementtemp).eq.'S') then
          zelemtemp=16
        elseif(trim(elementtemp).eq.'Cl') then
          zelemtemp=17
        elseif(trim(elementtemp).eq.'Ar') then
          zelemtemp=18
        elseif(trim(elementtemp).eq.'K') then
          zelemtemp=19
        elseif(trim(elementtemp).eq.'Ca') then
          zelemtemp=20
        elseif(trim(elementtemp).eq.'Sc') then
          zelemtemp=21
        elseif(trim(elementtemp).eq.'Ti') then
          zelemtemp=22
        elseif(trim(elementtemp).eq.'V') then
          zelemtemp=23
        elseif(trim(elementtemp).eq.'Cr') then
          zelemtemp=24
        elseif(trim(elementtemp).eq.'Mn') then
          zelemtemp=25
        elseif(trim(elementtemp).eq.'Fe') then
          zelemtemp=26
        elseif(trim(elementtemp).eq.'Co') then
          zelemtemp=27
        elseif(trim(elementtemp).eq.'Ni') then
          zelemtemp=28
        elseif(trim(elementtemp).eq.'Cu') then
          zelemtemp=29
        elseif(trim(elementtemp).eq.'Zn') then
          zelemtemp=30
        elseif(trim(elementtemp).eq.'Ga') then
          zelemtemp=31
        elseif(trim(elementtemp).eq.'Ge') then
          zelemtemp=32
        elseif(trim(elementtemp).eq.'As') then
          zelemtemp=33
        elseif(trim(elementtemp).eq.'Se') then
          zelemtemp=34
        elseif(trim(elementtemp).eq.'Br') then
          zelemtemp=35
        elseif(trim(elementtemp).eq.'Kr') then
          zelemtemp=36
        elseif(trim(elementtemp).eq.'Rb') then
          zelemtemp=37
        elseif(trim(elementtemp).eq.'Sr') then
          zelemtemp=38
        elseif(trim(elementtemp).eq.'Y') then
          zelemtemp=39
        elseif(trim(elementtemp).eq.'Zr') then
          zelemtemp=40
        elseif(trim(elementtemp).eq.'Nb') then
          zelemtemp=41
        elseif(trim(elementtemp).eq.'Mo') then
          zelemtemp=42
        elseif(trim(elementtemp).eq.'Tc') then
          zelemtemp=43
        elseif(trim(elementtemp).eq.'Ru') then
          zelemtemp=44
        elseif(trim(elementtemp).eq.'Rh') then
          zelemtemp=45
        elseif(trim(elementtemp).eq.'Pd') then
          zelemtemp=46
        elseif(trim(elementtemp).eq.'Ag') then
          zelemtemp=47
        elseif(trim(elementtemp).eq.'Cd') then
          zelemtemp=48
        elseif(trim(elementtemp).eq.'In') then
          zelemtemp=49
        elseif(trim(elementtemp).eq.'Sn') then
          zelemtemp=50
        elseif(trim(elementtemp).eq.'Sb') then
          zelemtemp=51
        elseif(trim(elementtemp).eq.'Te') then
          zelemtemp=52
        elseif(trim(elementtemp).eq.'I') then
          zelemtemp=53
        elseif(trim(elementtemp).eq.'Xe') then
          zelemtemp=54
        elseif(trim(elementtemp).eq.'Cs') then
          zelemtemp=55
        elseif(trim(elementtemp).eq.'Ba') then
          zelemtemp=56
        elseif(trim(elementtemp).eq.'La') then
          zelemtemp=57
        elseif(trim(elementtemp).eq.'Ce') then
          zelemtemp=58
        elseif(trim(elementtemp).eq.'Pr') then
          zelemtemp=59
        elseif(trim(elementtemp).eq.'Nd') then
          zelemtemp=60
        elseif(trim(elementtemp).eq.'Pm') then
          zelemtemp=61
        elseif(trim(elementtemp).eq.'Sm') then
          zelemtemp=62
        elseif(trim(elementtemp).eq.'Eu') then
          zelemtemp=63
        elseif(trim(elementtemp).eq.'Gd') then
          zelemtemp=64
        elseif(trim(elementtemp).eq.'Tb') then
          zelemtemp=65
        elseif(trim(elementtemp).eq.'Dy') then
          zelemtemp=66
        elseif(trim(elementtemp).eq.'Ho') then
          zelemtemp=67
        elseif(trim(elementtemp).eq.'Er') then
          zelemtemp=68
        elseif(trim(elementtemp).eq.'Tm') then
          zelemtemp=69
        elseif(trim(elementtemp).eq.'Yb') then
          zelemtemp=70
        elseif(trim(elementtemp).eq.'Lu') then
          zelemtemp=71
        elseif(trim(elementtemp).eq.'Hf') then
          zelemtemp=72
        elseif(trim(elementtemp).eq.'Ta') then
          zelemtemp=73
        elseif(trim(elementtemp).eq.'W') then
          zelemtemp=74
        elseif(trim(elementtemp).eq.'Re') then
          zelemtemp=75
        elseif(trim(elementtemp).eq.'Os') then
          zelemtemp=76
        elseif(trim(elementtemp).eq.'Ir') then
          zelemtemp=77
        elseif(trim(elementtemp).eq.'Pt') then
          zelemtemp=78
        elseif(trim(elementtemp).eq.'Au') then
          zelemtemp=79
        elseif(trim(elementtemp).eq.'Hg') then
          zelemtemp=80
        elseif(trim(elementtemp).eq.'Tl') then
          zelemtemp=81
        elseif(trim(elementtemp).eq.'Pb') then
          zelemtemp=82
        elseif(trim(elementtemp).eq.'Bi') then
          zelemtemp=83
        elseif(trim(elementtemp).eq.'Po') then
          zelemtemp=84
        elseif(trim(elementtemp).eq.'At') then
          zelemtemp=85
        elseif(trim(elementtemp).eq.'Rn') then
          zelemtemp=86
        elseif(trim(elementtemp).eq.'Fr') then
          zelemtemp=87
        elseif(trim(elementtemp).eq.'Ra') then
          zelemtemp=88
        elseif(trim(elementtemp).eq.'Ac') then
          zelemtemp=89
        elseif(trim(elementtemp).eq.'Th') then
          zelemtemp=90
        elseif(trim(elementtemp).eq.'Pa') then
          zelemtemp=91
        elseif(trim(elementtemp).eq.'U') then
          zelemtemp=92
        elseif(trim(elementtemp).eq.'Np') then
          zelemtemp=93
        elseif(trim(elementtemp).eq.'Pu') then
          zelemtemp=94
        elseif(trim(elementtemp).eq.'Am') then
          zelemtemp=95
        elseif(trim(elementtemp).eq.'Cm') then
          zelemtemp=96
        elseif(trim(elementtemp).eq.'Bk') then
          zelemtemp=97
        elseif(trim(elementtemp).eq.'Cf') then
          zelemtemp=98
        elseif(trim(elementtemp).eq.'Es') then
          zelemtemp=99
        elseif(trim(elementtemp).eq.'Fm') then
          zelemtemp=100
        elseif(trim(elementtemp).eq.'Md') then
          zelemtemp=101
        elseif(trim(elementtemp).eq.'No') then
          zelemtemp=102
        else
          write(*,*)'Error: Unknown element ',elementtemp,' error'
          stop
        endif
!!
      return
      end
