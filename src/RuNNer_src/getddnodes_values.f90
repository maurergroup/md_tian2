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
!! - getdfpairdw.f90
!! - getdfpairdw_para.f90
!! - getdfshortdw.f90
!! - getdfshortdw_para.f90
!!
      subroutine getddnodes_values(maxnum_layers_local,num_layers_local,maxnodes_local,&
        nodes_local,itemp,ndim,nodes_sum_local,nodes_values_local,dnodes_values_local,&
        ddnodes_values_local,actfunc)
!!
      use fileunits
      use globaloptions
!!
      implicit none
!!
      integer maxnum_layers_local                                          ! in
      integer num_layers_local(ndim)                                       ! in
      integer maxnodes_local                                               ! in
      integer nodes_local(0:maxnum_layers_local,ndim)                      ! in
      integer i2,i3                                                        ! internal
      integer ndim                                                         ! in
      integer itemp                                                        ! in
!!
      real*8 nodes_values_local(maxnum_layers_local,maxnodes_local)        ! in
      real*8 nodes_sum_local(maxnum_layers_local,maxnodes_local)           ! in
      real*8 dnodes_values_local(maxnum_layers_local,maxnodes_local)       ! in 
      real*8 ddnodes_values_local(maxnum_layers_local,maxnodes_local)      ! out 
      real*8 alphagaussian                                                 ! internal
      real*8 norm                                                          ! internal
!!
      logical lnan
      logical isnan
!!
      character*1 actfunc(maxnodes_local,maxnum_layers_local,ndim)         ! in
!!
!! initializations
      alphagaussian     = 0.5d0
      ddnodes_values_local(:,:)= 0.0d0 ! initialization
      norm=1.0d0
!!
!!
        do i2=1,num_layers_local(itemp)
          if(lnormnodes)then
            norm=1.d0/dble(nodes_local(i2-1,itemp))
          endif
          do i3=1,nodes_local(i2,itemp)
            if(actfunc(i3,i2,itemp).eq.'l')then
              ddnodes_values_local(i2,i3)=0.d0              ! CHECK!
            elseif(actfunc(i3,i2,itemp).eq.'t')then
              ddnodes_values_local(i2,i3)&
              =-2.d0*nodes_values_local(i2,i3)*dnodes_values_local(i2,i3) ! tanh(x)''=-2*tanh(x)*(1-tanh(x)**2)
            elseif(actfunc(i3,i2,itemp).eq.'g')then
              ddnodes_values_local(i2,i3)&
              =2.d0*alphagaussian*(-1.d0+2.d0*alphagaussian*nodes_sum_local(i2,i3)**2)*nodes_values_local(i2,i3)
            elseif(actfunc(i3,i2,itemp).eq.'c')then
              ddnodes_values_local(i2,i3)&
              =-1.0d0*nodes_values_local(i2,i3) ! cos(x)''=-cos(x)
            elseif(actfunc(i3,i2,itemp).eq.'s')then
              ddnodes_values_local(i2,i3)&
                =-1.d0*dnodes_values_local(i2,i3) &
                + (2.d0*dexp(-2.d0*nodes_sum_local(i2,i3)))/((1.d0+dexp(-1.d0*nodes_sum_local(i2,i3)))**3.d0)
!! Caution: for large values of nodes_sum_local it can happen that dexp(-nodes_sum_local) becomes infinity
!! => ddnodes_values_local becomes NaN. We check this here and stop in that case
              lnan=isnan(ddnodes_values_local(i2,i3))
              if(lnan)then
                write(*,*)'ddnodes_values_local(i2,i3) ',ddnodes_values_local(i2,i3),i2,i3
!!                write(*,*)'NaN in ddnodes_values_local(i2,i3) ',i2,i3
!! dirty fix for NaN here: ddnodes_values_local is 0 for very small or very large nodes_sum_local
!! => if this causes numerical problems, we set it manually
!!                write(*,*)'ddnodes_values_local set to 0'
!!                ddnodes_values_local(i2,i3)=0.0d0
                stop
              endif
            elseif(actfunc(i3,i2,itemp).eq.'S')then
              ddnodes_values_local(i2,i3)&
                =dnodes_values_local(i2,i3) &
                - (2.d0*dexp(-2.d0*nodes_sum_local(i2,i3)))&
                /((1.d0+dexp(-1.d0*nodes_sum_local(i2,i3)))**3.d0)
            elseif(actfunc(i3,i2,itemp).eq.'e')then
              ddnodes_values_local(i2,i3)=nodes_values_local(i2,i3)
!! 2nd derivative of softplus is 1st derivative of sigmoid and can be written as
!! ( e^y - 1 ) / e^(2y), where y is softplus. Since ( e^y - 1 ) / e^y is the first derivativei f',
!! the whole expression can also be written as f'/e^y.
            elseif(actfunc(i3,i2,itemp).eq.'p')then
              ddnodes_values_local(i2,i3)= dnodes_values_local(i2,i3) / dexp( nodes_values_local(i2,i3) )
!! 2nd derivative of harmonic potential function
            elseif(actfunc(i3,i2,itemp).eq.'h')then
              ddnodes_values_local(i2,i3)= 2d0 
            else
              write(ounit,*)'Error, activation function not implemented ',actfunc(i3,i2,itemp)
              stop !'
            endif
            if(lnormnodes)then
              ddnodes_values_local(i2,i3)=ddnodes_values_local(i2,i3)*norm
            endif
!!
          enddo ! i3
        enddo ! i2
!!
!!
      return
      end
