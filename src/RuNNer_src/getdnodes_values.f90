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
!! multipurpose subroutine

!! called by:
!! - getcoulombdchargedxyz.f90
!! - getcoulombdchargedxyz_para.f90 
!! - getcoulombforces.f90
!! - getdchargedxyz.f90
!! - getdfpairdw.f90
!! - getdfpairdw_para.f90
!! - getdfshortdw.f90
!! - getdfshortdw_para.f90
!! - getonedeshortdw.f90
!! - getoneshortforce.f90
!! - getoneshortforce_para.f90
!! - getoneshortforcepair_para.f90
!! - getshortforces.f90
!! - getshortforces_para.f90
!! - getshortforces_parapair.f90
!! - getshortforcespair.f90       
!!
      subroutine getdnodes_values(maxnum_layers_local,num_layers_local,maxnodes_local,&
        nodes_local,itemp,ndim,nodes_sum_local,nodes_values_local,dnodes_values_local,actfunc)
!!
      use fileunits
      use globaloptions
!!
      implicit none
!!
      integer maxnum_layers_local                                           ! in
      integer num_layers_local(ndim)                                        ! in
      integer maxnodes_local                                                ! in
      integer nodes_local(0:maxnum_layers_local,ndim)                       ! in
      integer i2,i3                                                         ! internal
      integer ndim                                                          ! in
      integer itemp                                                         ! in
!!
      real*8 nodes_values_local(maxnum_layers_local,maxnodes_local)         ! in
      real*8 nodes_sum_local(maxnum_layers_local,maxnodes_local)            ! in
      real*8 dnodes_values_local(maxnum_layers_local,maxnodes_local)        ! out
      real*8 alphagaussian                                                  ! internal
      real*8 norm                                                           ! internal
!!
      character*1 actfunc(maxnodes_local,maxnum_layers_local,ndim)          ! in
!!
!!
!! initializations
      alphagaussian     = 0.5d0
      dnodes_values_local(:,:)= 0.0d0 ! initialization
      norm=1.0d0
!!
!! calculate part of the derivative at each node (\frac{\partial f_a(G_i)}{\partial G_i} )
        do i2=1,num_layers_local(itemp)
          if(lnormnodes)then
            norm=1.d0/dble(nodes_local(i2-1,itemp))
          endif
          do i3=1,nodes_local(i2,itemp)
            if(actfunc(i3,i2,itemp).eq.'l')then
              dnodes_values_local(i2,i3)=1.d0
            elseif(actfunc(i3,i2,itemp).eq.'t')then
              dnodes_values_local(i2,i3)&
              =1.d0-nodes_values_local(i2,i3)**2 ! tanh(x)'=1-tanh(x)**2
            elseif(actfunc(i3,i2,itemp).eq.'g')then
              dnodes_values_local(i2,i3)&
              =-2.d0*alphagaussian*nodes_sum_local(i2,i3)*nodes_values_local(i2,i3)
            elseif(actfunc(i3,i2,itemp).eq.'c')then
              dnodes_values_local(i2,i3)&
              =-1.d0*dsin(nodes_sum_local(i2,i3)) ! cos(x)'=-sin(x)
            elseif(actfunc(i3,i2,itemp).eq.'s')then
!! FIXME: for sigmoid we have: f'(x)=f(x)[1-f(x)], this should be faster to calculate, check!
              dnodes_values_local(i2,i3)&
              =dexp(-1.d0*nodes_sum_local(i2,i3))&
                /(1.d0+dexp(-1.d0*nodes_sum_local(i2,i3)))**2.0d0 ! (1/(1+exp(-x)))'=exp(-x)/(1+exp(-x))**2
            elseif(actfunc(i3,i2,itemp).eq.'S')then
              dnodes_values_local(i2,i3)&
              =-1.d0*dexp(-1.d0*nodes_sum_local(i2,i3))/(1.d0+dexp(-1.d0*nodes_sum_local(i2,i3)))**2.0d0
            elseif(actfunc(i3,i2,itemp).eq.'e')then
              dnodes_values_local(i2,i3)&
              =-1.d0*nodes_values_local(i2,i3) ! exp(-x)'=-exp(-x)
!! MG: Derivative of new softplus activation function is the sigmoid function,
!! MG: which can be written as ( e^y - 1.0 )/e^y where y is the softplus function.
            elseif(actfunc(i3,i2,itemp).eq.'p')then
              dnodes_values_local(i2,i3)&
              =( dexp(nodes_values_local(i2,i3)) - 1d0 )/dexp(nodes_values_local(i2,i3)) 
!! MG: Harmonic potential derivative
            elseif(actfunc(i3,i2,itemp).eq.'h')then
              dnodes_values_local(i2,i3) = 2d0*nodes_sum_local(i2,i3) 
            else
              write(ounit,*)'Error, activation function not implemented ',actfunc(i3,i2,itemp)
              stop !'
            endif
            if(lnormnodes)then
              dnodes_values_local(i2,i3)=norm*dnodes_values_local(i2,i3)
            endif
          enddo ! i3
        enddo ! i2
!!
      return
      end
