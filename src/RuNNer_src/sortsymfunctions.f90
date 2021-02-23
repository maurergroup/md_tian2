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
!! Multipurpose subroutine

!! called by:
!! - readinput.f90
!!
      subroutine sortsymfunctions(&
        maxnum_funcvalues_local,num_funcvalues_local,&
        function_type_local,symelement_local,&
        eta_loca,zeta_loca,rshift_local,lambda_local,funccutoff_local)
!!
      use fileunits
      use globaloptions
!!
      implicit none
!!
      integer i1,i2
      integer maxnum_funcvalues_local                         ! in
      integer num_funcvalues_local(nelem)                     ! in
      integer function_type_local_temp                        ! internal
      integer function_type_local(maxnum_funcvalues_local,nelem)    ! in/out
      integer symelement_local(maxnum_funcvalues_local,2,nelem)     ! in/out
      integer symelement_local_temp(2)                        ! internal
      integer counter(9,nelem)                                ! internal
      integer num_type5                                       ! internal
      integer num_type6                                       ! internal
      integer num_others                                      ! internal
!!
      real*8 funccutoff_local(maxnum_funcvalues_local,nelem)        ! in/out
      real*8 eta_loca(maxnum_funcvalues_local,nelem)               ! in/out
      real*8 zeta_loca(maxnum_funcvalues_local,nelem)              ! in/out
      real*8 lambda_local(maxnum_funcvalues_local,nelem)            ! in/out
      real*8 rshift_local(maxnum_funcvalues_local,nelem)            ! in/out
      real*8 funccutoff_local_temp                            ! internal
      real*8 eta_loca_temp                                   ! internal
      real*8 zeta_loca_temp                                  ! internal
      real*8 lambda_local_temp                                ! internal
      real*8 rshift_local_temp                                ! internal
      real*8 thres                                      ! internal
!!
!!
      thres=0.0001d0
!!
!! check if there are unkown symmetry function types
!!
      do i1=1,nelem
        do i2=1,num_funcvalues_local(i1)
          if((function_type_local(i2,i1).gt.9).or.(function_type_local(i2,i1).eq.7))then
            write(ounit,*)'ERROR in sortsymfunctions'
            write(ounit,*)'unknown function_type_local ',function_type_local(i2,i1)
            write(ounit,*)'redimension also counter in sortsymfunctions.f90!'
            stop
          endif
        enddo
      enddo
!!
!! first sort according to function_type_local
!!
      do i1=1,nelem
        if(num_funcvalues_local(i1).gt.1)then
 10       continue
          do i2=1,num_funcvalues_local(i1)-1
            if(function_type_local(i2,i1).gt.function_type_local(i2+1,i1))then
              function_type_local_temp    =function_type_local(i2,i1)
              function_type_local(i2,i1)  =function_type_local(i2+1,i1)
              function_type_local(i2+1,i1)=function_type_local_temp
              symelement_local_temp(1)    =symelement_local(i2,1,i1)
              symelement_local(i2,1,i1)   =symelement_local(i2+1,1,i1)
              symelement_local(i2+1,1,i1) =symelement_local_temp(1)
              symelement_local_temp(2)    =symelement_local(i2,2,i1)
              symelement_local(i2,2,i1)   =symelement_local(i2+1,2,i1)
              symelement_local(i2+1,2,i1) =symelement_local_temp(2)
              eta_loca_temp    =eta_loca(i2,i1)
              eta_loca(i2,i1)  =eta_loca(i2+1,i1)
              eta_loca(i2+1,i1)=eta_loca_temp
              zeta_loca_temp    =zeta_loca(i2,i1)
              zeta_loca(i2,i1)  =zeta_loca(i2+1,i1)
              zeta_loca(i2+1,i1)=zeta_loca_temp
              lambda_local_temp    =lambda_local(i2,i1)
              lambda_local(i2,i1)  =lambda_local(i2+1,i1)
              lambda_local(i2+1,i1)=lambda_local_temp
              rshift_local_temp    =rshift_local(i2,i1)
              rshift_local(i2,i1)  =rshift_local(i2+1,i1)
              rshift_local(i2+1,i1)=rshift_local_temp
              funccutoff_local_temp    =funccutoff_local(i2,i1)
              funccutoff_local(i2,i1)  =funccutoff_local(i2+1,i1)
              funccutoff_local(i2+1,i1)=funccutoff_local_temp
              goto 10
            endif
          enddo ! i2
        endif
      enddo ! i1
!!
!! count symmetry functions of each type
!!
      counter(:,:)=0
      do i1=1,nelem
        do i2=1,num_funcvalues_local(i1)
          counter(function_type_local(i2,i1),i1)=counter(function_type_local(i2,i1),i1)+1
        enddo ! i2
      enddo ! i1
!!
!! sort according to cutoff
!!
      do i1=1,nelem
        if(num_funcvalues_local(i1).gt.1)then
 11       continue
          do i2=1,num_funcvalues_local(i1)-1
            if((function_type_local(i2,i1).eq.function_type_local(i2+1,i1))&
              .and.(funccutoff_local(i2,i1).gt.funccutoff_local(i2+1,i1)))then
              function_type_local_temp    =function_type_local(i2,i1)
              function_type_local(i2,i1)  =function_type_local(i2+1,i1)
              function_type_local(i2+1,i1)=function_type_local_temp
              symelement_local_temp(1)    =symelement_local(i2,1,i1)
              symelement_local(i2,1,i1)   =symelement_local(i2+1,1,i1)
              symelement_local(i2+1,1,i1) =symelement_local_temp(1)
              symelement_local_temp(2)    =symelement_local(i2,2,i1)
              symelement_local(i2,2,i1)   =symelement_local(i2+1,2,i1)
              symelement_local(i2+1,2,i1) =symelement_local_temp(2)
              eta_loca_temp    =eta_loca(i2,i1)
              eta_loca(i2,i1)  =eta_loca(i2+1,i1)
              eta_loca(i2+1,i1)=eta_loca_temp
              zeta_loca_temp    =zeta_loca(i2,i1)
              zeta_loca(i2,i1)  =zeta_loca(i2+1,i1)
              zeta_loca(i2+1,i1)=zeta_loca_temp
              lambda_local_temp    =lambda_local(i2,i1)
              lambda_local(i2,i1)  =lambda_local(i2+1,i1)
              lambda_local(i2+1,i1)=lambda_local_temp
              rshift_local_temp    =rshift_local(i2,i1)
              rshift_local(i2,i1)  =rshift_local(i2+1,i1)
              rshift_local(i2+1,i1)=rshift_local_temp
              funccutoff_local_temp    =funccutoff_local(i2,i1)
              funccutoff_local(i2,i1)  =funccutoff_local(i2+1,i1)
              funccutoff_local(i2+1,i1)=funccutoff_local_temp
              goto 11
            endif
          enddo ! i2
        endif
      enddo ! i1
!!
!! sort according to eta_loca
!!
      do i1=1,nelem
        if(num_funcvalues_local(i1).gt.1)then
 12       continue
          do i2=1,num_funcvalues_local(i1)-1
            if((function_type_local(i2,i1).eq.function_type_local(i2+1,i1))&
              .and.(funccutoff_local(i2,i1).eq.funccutoff_local(i2+1,i1))&
              .and.(eta_loca(i2,i1).gt.eta_loca(i2+1,i1)))then
              function_type_local_temp    =function_type_local(i2,i1)
              function_type_local(i2,i1)  =function_type_local(i2+1,i1)
              function_type_local(i2+1,i1)=function_type_local_temp
              symelement_local_temp(1)    =symelement_local(i2,1,i1)
              symelement_local(i2,1,i1)   =symelement_local(i2+1,1,i1)
              symelement_local(i2+1,1,i1) =symelement_local_temp(1)
              symelement_local_temp(2)    =symelement_local(i2,2,i1)
              symelement_local(i2,2,i1)   =symelement_local(i2+1,2,i1)
              symelement_local(i2+1,2,i1) =symelement_local_temp(2)
              eta_loca_temp    =eta_loca(i2,i1)
              eta_loca(i2,i1)  =eta_loca(i2+1,i1)
              eta_loca(i2+1,i1)=eta_loca_temp
              zeta_loca_temp    =zeta_loca(i2,i1)
              zeta_loca(i2,i1)  =zeta_loca(i2+1,i1)
              zeta_loca(i2+1,i1)=zeta_loca_temp
              lambda_local_temp    =lambda_local(i2,i1)
              lambda_local(i2,i1)  =lambda_local(i2+1,i1)
              lambda_local(i2+1,i1)=lambda_local_temp
              rshift_local_temp    =rshift_local(i2,i1)
              rshift_local(i2,i1)  =rshift_local(i2+1,i1)
              rshift_local(i2+1,i1)=rshift_local_temp
              funccutoff_local_temp    =funccutoff_local(i2,i1)
              funccutoff_local(i2,i1)  =funccutoff_local(i2+1,i1)
              funccutoff_local(i2+1,i1)=funccutoff_local_temp
              goto 12
            endif
          enddo ! i2
        endif
      enddo ! i1
!!
!! sort according to zeta_loca
!!
      do i1=1,nelem
        if(num_funcvalues_local(i1).gt.1)then
 13       continue
          do i2=1,num_funcvalues_local(i1)-1
            if((function_type_local(i2,i1).eq.function_type_local(i2+1,i1))&
              .and.(funccutoff_local(i2,i1).eq.funccutoff_local(i2+1,i1))&
              .and.(eta_loca(i2,i1).eq.eta_loca(i2+1,i1))&
              .and.(zeta_loca(i2,i1).gt.zeta_loca(i2+1,i1)))then
              function_type_local_temp    =function_type_local(i2,i1)
              function_type_local(i2,i1)  =function_type_local(i2+1,i1)
              function_type_local(i2+1,i1)=function_type_local_temp
              symelement_local_temp(1)    =symelement_local(i2,1,i1)
              symelement_local(i2,1,i1)   =symelement_local(i2+1,1,i1)
              symelement_local(i2+1,1,i1) =symelement_local_temp(1)
              symelement_local_temp(2)    =symelement_local(i2,2,i1)
              symelement_local(i2,2,i1)   =symelement_local(i2+1,2,i1)
              symelement_local(i2+1,2,i1) =symelement_local_temp(2)
              eta_loca_temp    =eta_loca(i2,i1)
              eta_loca(i2,i1)  =eta_loca(i2+1,i1)
              eta_loca(i2+1,i1)=eta_loca_temp
              zeta_loca_temp    =zeta_loca(i2,i1)
              zeta_loca(i2,i1)  =zeta_loca(i2+1,i1)
              zeta_loca(i2+1,i1)=zeta_loca_temp
              lambda_local_temp    =lambda_local(i2,i1)
              lambda_local(i2,i1)  =lambda_local(i2+1,i1)
              lambda_local(i2+1,i1)=lambda_local_temp
              rshift_local_temp    =rshift_local(i2,i1)
              rshift_local(i2,i1)  =rshift_local(i2+1,i1)
              rshift_local(i2+1,i1)=rshift_local_temp
              funccutoff_local_temp    =funccutoff_local(i2,i1)
              funccutoff_local(i2,i1)  =funccutoff_local(i2+1,i1)
              funccutoff_local(i2+1,i1)=funccutoff_local_temp
              goto 13
            endif
          enddo ! i2
        endif
      enddo ! i1
!!
!! sort according to lambda_local 
!!
      do i1=1,nelem
        if(num_funcvalues_local(i1).gt.1)then
 14       continue
          do i2=1,num_funcvalues_local(i1)-1
            if((function_type_local(i2,i1).eq.function_type_local(i2+1,i1))&
              .and.(funccutoff_local(i2,i1).eq.funccutoff_local(i2+1,i1))&
              .and.(eta_loca(i2,i1).eq.eta_loca(i2+1,i1))&
              .and.(zeta_loca(i2,i1).eq.zeta_loca(i2+1,i1))&
              .and.(lambda_local(i2,i1).gt.lambda_local(i2+1,i1)))then
              function_type_local_temp    =function_type_local(i2,i1)
              function_type_local(i2,i1)  =function_type_local(i2+1,i1)
              function_type_local(i2+1,i1)=function_type_local_temp
              symelement_local_temp(1)    =symelement_local(i2,1,i1)
              symelement_local(i2,1,i1)   =symelement_local(i2+1,1,i1)
              symelement_local(i2+1,1,i1) =symelement_local_temp(1)
              symelement_local_temp(2)    =symelement_local(i2,2,i1)
              symelement_local(i2,2,i1)   =symelement_local(i2+1,2,i1)
              symelement_local(i2+1,2,i1) =symelement_local_temp(2)
              eta_loca_temp    =eta_loca(i2,i1)
              eta_loca(i2,i1)  =eta_loca(i2+1,i1)
              eta_loca(i2+1,i1)=eta_loca_temp
              zeta_loca_temp    =zeta_loca(i2,i1)
              zeta_loca(i2,i1)  =zeta_loca(i2+1,i1)
              zeta_loca(i2+1,i1)=zeta_loca_temp
              lambda_local_temp    =lambda_local(i2,i1)
              lambda_local(i2,i1)  =lambda_local(i2+1,i1)
              lambda_local(i2+1,i1)=lambda_local_temp
              rshift_local_temp    =rshift_local(i2,i1)
              rshift_local(i2,i1)  =rshift_local(i2+1,i1)
              rshift_local(i2+1,i1)=rshift_local_temp
              funccutoff_local_temp    =funccutoff_local(i2,i1)
              funccutoff_local(i2,i1)  =funccutoff_local(i2+1,i1)
              funccutoff_local(i2+1,i1)=funccutoff_local_temp
              goto 14
            endif
          enddo ! i2
        endif
      enddo ! i1
!!
!! sort according to rshift_local 
!!
      do i1=1,nelem
        if(num_funcvalues_local(i1).gt.1)then
 15       continue
          do i2=1,num_funcvalues_local(i1)-1
            if((function_type_local(i2,i1).eq.function_type_local(i2+1,i1))&
              .and.(funccutoff_local(i2,i1).eq.funccutoff_local(i2+1,i1))&
              .and.(eta_loca(i2,i1).eq.eta_loca(i2+1,i1))&
              .and.(zeta_loca(i2,i1).eq.zeta_loca(i2+1,i1))&
              .and.(lambda_local(i2,i1).eq.lambda_local(i2+1,i1))&
              .and.(rshift_local(i2,i1).gt.rshift_local(i2+1,i1)))then
              function_type_local_temp    =function_type_local(i2,i1)
              function_type_local(i2,i1)  =function_type_local(i2+1,i1)
              function_type_local(i2+1,i1)=function_type_local_temp
              symelement_local_temp(1)    =symelement_local(i2,1,i1)
              symelement_local(i2,1,i1)   =symelement_local(i2+1,1,i1)
              symelement_local(i2+1,1,i1) =symelement_local_temp(1)
              symelement_local_temp(2)    =symelement_local(i2,2,i1)
              symelement_local(i2,2,i1)   =symelement_local(i2+1,2,i1)
              symelement_local(i2+1,2,i1) =symelement_local_temp(2)
              eta_loca_temp    =eta_loca(i2,i1)
              eta_loca(i2,i1)  =eta_loca(i2+1,i1)
              eta_loca(i2+1,i1)=eta_loca_temp
              zeta_loca_temp    =zeta_loca(i2,i1)
              zeta_loca(i2,i1)  =zeta_loca(i2+1,i1)
              zeta_loca(i2+1,i1)=zeta_loca_temp
              lambda_local_temp    =lambda_local(i2,i1)
              lambda_local(i2,i1)  =lambda_local(i2+1,i1)
              lambda_local(i2+1,i1)=lambda_local_temp
              rshift_local_temp    =rshift_local(i2,i1)
              rshift_local(i2,i1)  =rshift_local(i2+1,i1)
              rshift_local(i2+1,i1)=rshift_local_temp
              funccutoff_local_temp    =funccutoff_local(i2,i1)
              funccutoff_local(i2,i1)  =funccutoff_local(i2+1,i1)
              funccutoff_local(i2+1,i1)=funccutoff_local_temp
              goto 15
            endif
          enddo ! i2
        endif
      enddo ! i1
!!
!! sort according to symelement_local(:,1,:) 
!!
      do i1=1,nelem
        if(num_funcvalues_local(i1).gt.1)then
 16       continue
          do i2=1,num_funcvalues_local(i1)-1
            if((function_type_local(i2,i1).eq.function_type_local(i2+1,i1))&
              .and.(funccutoff_local(i2,i1).eq.funccutoff_local(i2+1,i1))&
              .and.(eta_loca(i2,i1).eq.eta_loca(i2+1,i1))&
              .and.(zeta_loca(i2,i1).eq.zeta_loca(i2+1,i1))&
              .and.(lambda_local(i2,i1).eq.lambda_local(i2+1,i1))&
              .and.(rshift_local(i2,i1).eq.rshift_local(i2+1,i1))&
              .and.(symelement_local(i2,1,i1).gt.symelement_local(i2+1,1,i1)))then
              function_type_local_temp    =function_type_local(i2,i1)
              function_type_local(i2,i1)  =function_type_local(i2+1,i1)
              function_type_local(i2+1,i1)=function_type_local_temp
              symelement_local_temp(1)    =symelement_local(i2,1,i1)
              symelement_local(i2,1,i1)   =symelement_local(i2+1,1,i1)
              symelement_local(i2+1,1,i1) =symelement_local_temp(1)
              symelement_local_temp(2)    =symelement_local(i2,2,i1)
              symelement_local(i2,2,i1)   =symelement_local(i2+1,2,i1)
              symelement_local(i2+1,2,i1) =symelement_local_temp(2)
              eta_loca_temp    =eta_loca(i2,i1)
              eta_loca(i2,i1)  =eta_loca(i2+1,i1)
              eta_loca(i2+1,i1)=eta_loca_temp
              zeta_loca_temp    =zeta_loca(i2,i1)
              zeta_loca(i2,i1)  =zeta_loca(i2+1,i1)
              zeta_loca(i2+1,i1)=zeta_loca_temp
              lambda_local_temp    =lambda_local(i2,i1)
              lambda_local(i2,i1)  =lambda_local(i2+1,i1)
              lambda_local(i2+1,i1)=lambda_local_temp
              rshift_local_temp    =rshift_local(i2,i1)
              rshift_local(i2,i1)  =rshift_local(i2+1,i1)
              rshift_local(i2+1,i1)=rshift_local_temp
              funccutoff_local_temp    =funccutoff_local(i2,i1)
              funccutoff_local(i2,i1)  =funccutoff_local(i2+1,i1)
              funccutoff_local(i2+1,i1)=funccutoff_local_temp
              goto 16
            endif
          enddo ! i2
        endif
      enddo ! i1
!!
!! sort according to symelement_local(:,2,:)
!!
      do i1=1,nelem
        if(num_funcvalues_local(i1).gt.1)then
 17       continue
          do i2=1,num_funcvalues_local(i1)-1
            if((function_type_local(i2,i1).eq.function_type_local(i2+1,i1))&
              .and.(funccutoff_local(i2,i1).eq.funccutoff_local(i2+1,i1))&
              .and.(eta_loca(i2,i1).eq.eta_loca(i2+1,i1))&
              .and.(zeta_loca(i2,i1).eq.zeta_loca(i2+1,i1))&
              .and.(lambda_local(i2,i1).eq.lambda_local(i2+1,i1))&
              .and.(rshift_local(i2,i1).eq.rshift_local(i2+1,i1))&
              .and.(symelement_local(i2,1,i1).eq.symelement_local(i2+1,1,i1))&
              .and.(symelement_local(i2,2,i1).gt.symelement_local(i2+1,2,i1)))then
              function_type_local_temp    =function_type_local(i2,i1)
              function_type_local(i2,i1)  =function_type_local(i2+1,i1)
              function_type_local(i2+1,i1)=function_type_local_temp
              symelement_local_temp(1)    =symelement_local(i2,1,i1)
              symelement_local(i2,1,i1)   =symelement_local(i2+1,1,i1)
              symelement_local(i2+1,1,i1) =symelement_local_temp(1)
              symelement_local_temp(2)    =symelement_local(i2,2,i1)
              symelement_local(i2,2,i1)   =symelement_local(i2+1,2,i1)
              symelement_local(i2+1,2,i1) =symelement_local_temp(2)
              eta_loca_temp    =eta_loca(i2,i1)
              eta_loca(i2,i1)  =eta_loca(i2+1,i1)
              eta_loca(i2+1,i1)=eta_loca_temp
              zeta_loca_temp    =zeta_loca(i2,i1)
              zeta_loca(i2,i1)  =zeta_loca(i2+1,i1)
              zeta_loca(i2+1,i1)=zeta_loca_temp
              lambda_local_temp    =lambda_local(i2,i1)
              lambda_local(i2,i1)  =lambda_local(i2+1,i1)
              lambda_local(i2+1,i1)=lambda_local_temp
              rshift_local_temp    =rshift_local(i2,i1)
              rshift_local(i2,i1)  =rshift_local(i2+1,i1)
              rshift_local(i2+1,i1)=rshift_local_temp
              funccutoff_local_temp    =funccutoff_local(i2,i1)
              funccutoff_local(i2,i1)  =funccutoff_local(i2+1,i1)
              funccutoff_local(i2+1,i1)=funccutoff_local_temp
              goto 17
            endif
          enddo ! i2
        endif
      enddo ! i1
!!
!! check for double symmetry functions
!!
      do i1=1,nelem
        do i2=1,num_funcvalues_local(i1)-1
          if((function_type_local(i2,i1).eq.function_type_local(i2+1,i1))&
            .and.(abs(funccutoff_local(i2,i1)-funccutoff_local(i2+1,i1)).lt.thres)&
            .and.(abs(eta_loca(i2,i1)-eta_loca(i2+1,i1)).lt.thres)&
            .and.(abs(zeta_loca(i2,i1)-zeta_loca(i2+1,i1)).lt.thres)&
            .and.(abs(lambda_local(i2,i1)-lambda_local(i2+1,i1)).lt.thres)&
            .and.(abs(rshift_local(i2,i1)-rshift_local(i2+1,i1)).lt.thres)&
            .and.(symelement_local(i2,1,i1).eq.symelement_local(i2+1,1,i1))&
            .and.(symelement_local(i2,2,i1).eq.symelement_local(i2+1,2,i1)))then
            write(ounit,*)'ERROR: a symmetry function is specified twice!'
            write(ounit,*)'element is ',element(i1)
            write(ounit,*)'function_type_local ',function_type_local(i2,i1)
            write(ounit,'(a14,f10.3)')'funccutoff_local    ',funccutoff_local(i2,i1)
            write(ounit,'(a14,f10.3)')'eta_loca           ',eta_loca(i2,i1)
            write(ounit,'(a14,f10.3)')'zeta_loca          ',zeta_loca(i2,i1)
            write(ounit,'(a14,f10.3)')'lambda_local        ',lambda_local(i2,i1)
            write(ounit,'(a14,f10.3)')'rshift_local        ',rshift_local(i2,i1)
            write(ounit,*)'symelement_local(1) ',symelement_local(i2,1,i1)
            write(ounit,*)'symelement_local(2) ',symelement_local(i2,2,i1)
            stop
          endif
        enddo ! i2
      enddo ! i1
!!
!! JB 2019: check if function types 5 or 6 are wrongly used for PES construction
!!          for this we count symmetry functions of each type, types 5 and 6 must only be used alone!
      do i1=1,nelem
        num_type5=0
        num_type6=0
        num_others=0
        do i2=1,num_funcvalues_local(i1)
           if(function_type_local(i2,i1).eq.5)then
             num_type5=num_type5+1
           elseif(function_type_local(i2,i1).eq.6)then
             num_type6=num_type6+1
           else
             num_others=num_others+1
           endif
        enddo ! i2
        if(num_type5.gt.0)then
          if((num_type6.gt.0).or.(num_others.gt.0))then
            write(ounit,*)'ERROR: symmetry function type 5 must not be combined with other types!'
            stop
          endif
        endif
        if(num_type6.gt.0)then
          if((num_type5.gt.0).or.(num_others.gt.0))then
            write(ounit,*)'ERROR: symmetry function type 6 must not be combined with other types!'
            stop
          endif
        endif
      enddo ! i1
!!
      return
      end
