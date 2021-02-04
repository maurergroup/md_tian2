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
!! - calconecharge.f90
!! - calconecharge_para.f90
!! - calconeshort.f90
!! - calconeshort_para.f90
!! - calconeshortpair.f90
!! - calconeshort_parapair.f90
!! - getalphalm.f90
!! - getcoulombdchargedxyz_para.f90
!! - getcoulombforces.f90
!! - getdchargedxyz.f90
!! - getdfshortdw.f90
!! - getdfshortdw_para.f90
!! - getdfpairdw.f90
!! - getdfpairdw_para.f90
!! - getonedeshortdw.f90
!! - getoneshortforce.f90
!! - getoneshortforce_para.f90
!! - getshortforcespair.f90
!! - getoneshortforcepair_para.f90
!! - getshortforces.f90
!! - getshortforces_para.f90
!! - optimize_ewald.f90
!!
      subroutine calconenn(nnoutdim,maxnum_funcvalues_local,maxnodes_local,&
          maxnum_layers_local,num_layers_local,&
          maxnum_weights_local,nodes_local,&
          symfunction_atom,weights_local,&
          nodes_values_local,nodes_sum_local,&
          nnoutput,actfunc_local)
!!
      use globaloptions
      use saturation
      use nnflags
      use fileunits
!!
      implicit none
!!
      integer maxnum_funcvalues_local                                       ! in
      integer maxnum_layers_local                                           ! in
      integer num_layers_local                                              ! in
      integer maxnum_weights_local                                          ! in
      integer nodes_local(0:maxnum_layers_local)                            ! in
      integer maxnodes_local                                                ! in
      integer i1,i2,i3                                                      ! internal
      integer icount                                                        ! internal
      integer nnoutdim                                                      ! in
!!
      real*8 symfunction_atom(maxnum_funcvalues_local)                      ! in 
      real*8 weights_local(maxnum_weights_local)                            ! in
      real*8 nodes_sum_local(maxnum_layers_local,maxnodes_local)            ! out
      real*8 nodes_values_local(maxnum_layers_local,maxnodes_local)         ! out
!! CAUTION: just one output node is assumed here
      real*8 nnoutput(nnoutdim)                                             ! out
      real*8 alphagaussian                    
!!
      character*1 actfunc_local(maxnodes_local,maxnum_layers_local)         ! in
!!
!!
!!
!!
      icount                 = 0
      nodes_values_local(:,:)= 0.0d0
      nodes_sum_local(:,:)   = 0.0d0
      alphagaussian          = 0.5d0
!!
!!
!! loop over all hidden and output layers
      do i1=1,num_layers_local
        do i3=1,nodes_local(i1-1)
          do i2=1,nodes_local(i1)
            if(i1.eq.1)then ! i1-1 is input layer
              icount=icount+1
              nodes_values_local(i1,i2)&
              =nodes_values_local(i1,i2)+weights_local(icount)*symfunction_atom(i3) 
            else
              icount=icount+1
              nodes_values_local(i1,i2)&
              =nodes_values_local(i1,i2)+weights_local(icount)*nodes_values_local(i1-1,i3) 
            endif
          enddo ! i2
        enddo ! i3
!! add bias weight
        do i2=1,nodes_local(i1)
          icount=icount+1
          nodes_values_local(i1,i2)=nodes_values_local(i1,i2)+weights_local(icount)
!! if requested normalize by the number of nodes in previous layer
          if(lnormnodes)then
            nodes_values_local(i1,i2)=nodes_values_local(i1,i2)/dble(nodes_local(i1-1))
          endif
!! store node values before application of activation function (needed for derivatives) 
          nodes_sum_local(i1,i2)=nodes_values_local(i1,i2)
!! check for saturation of node if requested
          if(lshort.and.(nn_type_short.eq.1))then
            if(ldetect_saturation.and.(i1.lt.num_layers_local).and.(saturation_element.ne.0))then
              if(abs(nodes_values_local(i1,i2)).gt.saturation_threshold)then
                nodes_saturation(i1,i2,saturation_element)=nodes_saturation(i1,i2,saturation_element)+1
              endif
              nodes_total(i1,i2,saturation_element)=nodes_total(i1,i2,saturation_element)+1
            endif
          endif
!! apply activation function
          if(actfunc_local(i2,i1).eq."t")then
            nodes_values_local(i1,i2)=tanh(nodes_values_local(i1,i2))
          elseif(actfunc_local(i2,i1).eq."g")then
            nodes_values_local(i1,i2)=dexp(-alphagaussian*(nodes_values_local(i1,i2))**2)
          elseif(actfunc_local(i2,i1).eq."l")then
          elseif(actfunc_local(i2,i1).eq."c")then
            nodes_values_local(i1,i2)=dcos(nodes_values_local(i1,i2))
          elseif(actfunc_local(i2,i1).eq."s")then
            nodes_values_local(i1,i2)=1.d0/(1.d0+dexp(-1.d0*nodes_values_local(i1,i2)))
          elseif(actfunc_local(i2,i1).eq."S")then
            nodes_values_local(i1,i2)=1.d0-1.d0/(1.d0+dexp(-1.d0*nodes_values_local(i1,i2)))
          elseif(actfunc_local(i2,i1).eq."e")then
            nodes_values_local(i1,i2)=dexp(-1.d0*nodes_values_local(i1,i2))
!! MG: New softplus activation function ln( e^x + 1 )
          elseif(actfunc_local(i2,i1).eq."p")then
            nodes_values_local(i1,i2)=dlog( dexp(nodes_values_local(i1,i2)) + 1d0 )
!! MG: Harmonic potential as activation function
          elseif(actfunc_local(i2,i1).eq."h")then
            nodes_values_local(i1,i2)=nodes_values_local(i1,i2)**2
          else
            write(ounit,*)"Error: Unknown activation function ",actfunc_local(i2,i1)
            stop
          endif
        enddo ! i2 loop over all nodes
      enddo ! i1 loop over all hidden layers and output layers
!!
!! prepare output of the NN
      do i1=1,nodes_local(num_layers_local)
        nnoutput(i1)=nodes_values_local(num_layers_local,i1)
      enddo
!!
!!      write(ounit,*)'1nodes_sum ',nodes_sum_local
!!      write(ounit,*)'1nodes_values ',nodes_values_local
      return
      end
