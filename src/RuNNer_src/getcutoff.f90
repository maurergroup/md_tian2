!######################################################################
! This routine is part of
! RuNNer - RuNNer Neural Network Energy Representation
! (c) 2008-2019 Prof. Dr. Joerg Behler 
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
!######################################################################!!                                                                    
!!
      subroutine getcutoff(&
        cutoff_type,cutoff_alpha,maxnum_funcvalues_local,nelem,&
        i2,iindex,&
        funccutoff_local,r_local,fcut_local,temp1)
!!
      use nnconstants
      use fileunits
!!
      implicit none
!!

      integer nelem                                          ! in
      integer i2                                             ! in
      integer iindex                                         ! in
      integer cutoff_type                                    ! in
      integer maxnum_funcvalues_local                        ! in

      real*8 temp1                                           ! out 
      real*8 r_local                                         ! in
      real*8 r_temp                                          ! internal
      real*8 fcut_local                                      ! out
      real*8 funccutoff_local(maxnum_funcvalues_local,nelem) ! in
      real*8 cutoff_alpha                                    ! in and modified by kenko
      real*8 rc                                              ! internal and modified by kenko
      real*8 rcinv                                           ! internal and modified by kenko
      real*8 iw                                              ! internal and modified by kenko
      real*8 x                                               ! internal and modified by kenko
      real*8 rci                                            ! internal and modified by kenko 


      if (cutoff_alpha .gt.1.0)then
         write(ounit,*)'ERROR: cutoff_alpha < 1.0 is required '
         stop
      end if 

      rc = funccutoff_local(i2,iindex)
      r_temp=r_local/rc    ! R_ij/R_c
      rci= cutoff_alpha*rc ! inner cutoff
      rcinv = 1.0d0/rc     ! inverse cutoff
      iw = 1.0d0/(rc-rci)     

      if(cutoff_type.eq.0)then
!!      fc=1 (be aware that this usually does not make sense)
        fcut_local=1.0d0
        temp1 =0.0d0      
      elseif(cutoff_type.eq.1)then
!!      fc=0.5*(cos(pi*r/Rc)+1)
!!      fc'=0.5*(-sin(pi*r/Rc))*pi/Rc
        if (r_local .le. rci) then
            fcut_local = 1.0d0
            temp1 = 0.0d0
        else
            x = (r_local-rci)*iw
            fcut_local=0.5d0*(dcos(pi*x)+1.d0)
            temp1 =(-dsin(pi*x))*(pi/2.0d0)*iw
        end if
      elseif(cutoff_type.eq.2)then
!!      fc=(tanh(1-r/Rc))^3
!!      fc'=3*(tanh(1-r/Rc))^2*(1-(tanh(1-r/Rc))^2)*(-1/Rc)
!!         =-3/Rc * ((tanh(1-r/Rc))^2-(tanh(1-r/Rc))^4)
        if (r_local .le. rci) then
            fcut_local = (tanh(1.0d0))**3
            temp1 = 0.0d0
        else
            x = (r_local-rci)*iw
            fcut_local=(tanh(1.d0-r_temp))**3
            temp1 =(-3.d0/rc)&
                  *((tanh(1.d0-r_temp))**2 -(tanh(1.d0-r_temp))**4)
        end if 
      elseif(cutoff_type.eq.3)then
!!      fc=((e^1+1)/(e^1-1))^3*(tanh(1-r/Rc))^3
!!      fc'=((e^1+1)/(e^1-1))^3*(-3/Rc)*((tanh(1-r/Rc))^2-(tanh(1-r/Rc))^4)
        if (r_local .le. rci) then
            fcut_local =((exp(1.0d0)+1.0d0)/(exp(1.0d0)-1.0d0))**3*(tanh(1.0d0))**3
            temp1      = 0.0d0
        else  
            fcut_local=((exp(1.0d0)+1.0d0)/(exp(1.0d0)-1.0d0))**3 *(tanh(1.d0-r_temp))**3
            temp1 =((exp(1.0d0)+1.0d0)/(exp(1.0d0)-1.0d0))**3 &
                   *((-3.d0/rc)*((tanh(1.d0-r_temp))**2 -(tanh(1.d0-r_temp))**4))
        end if
      elseif(cutoff_type.eq.4)then
!!      fc=e^(1-1/(1-(r/Rc)^2))
!!      fc'=-(2*r/((1-(r/Rc)^2)^2))*(1/rc^2)*e^(-1/(1-(r/Rc)^2))
        if (r_local .le. rci) then
           fcut_local = 1.0d0
           temp1      = 0.0d0
        else
           x = (r_local-rci)*iw            
           fcut_local=exp(1.0d0-1.0d0/(1.0d0-x**2))
           temp1=(-2.0d0*iw*x/((1.0d0-x**2)**2))&
           *exp(-1.0d0/(1.0d0-x**2))*exp(1.0d0)
        end if
      elseif(cutoff_type.eq.5)then
!!      fc=(2r/Rc-3)*(r/Rc)^2+1
!!      fc'=(6r^2-6rRc)/Rc^3=(6(r/Rc)^2-6r/Rc)/Rc
        if (r_local .le. rci) then
           fcut_local = 1.0d0
           temp1      = 0.0d0
        else
           x = (r_local-rci)*iw
           fcut_local=(2.0d0*x-3.0d0)*x**2 + 1.0d0
           temp1 = iw*x*(6.0d0*x-6.0d0)
        end if
      elseif(cutoff_type.eq.6)then
!!      fc=((15-6r/Rc)r/Rc-10)*(r/Rc)^3+1
!!      fc=(r/Rc(r/Rc(20r/Rc-70)+84)-35)*(r/Rc)^4+1
        if (r_local .le. rci) then
           fcut_local = 1.0d0
           temp1      = 0.0d0 
        else
           x = (r_local-rci)*iw  
           fcut_local=((15.0d0-6.0d0*x)*x-10.0d0)*(x)**3+1.0d0
           temp1 = iw *x*x*((60.0d0-30.0d0*x)*x-30.0)
        end if
      elseif(cutoff_type.eq.7)then
!!      fc = (r/Rc(r/Rc(20r/Rc-70)+84)-35)*(r/Rc)^4+1
!!      fc'=[140(r/Rc)^3 *((r/Rc)^3-3(r/Rc)^2+3r/Rc-1)]/Rc
        if (r_local .le. rci) then
           fcut_local = 1.0d0
           temp1 = 0.0d0
        else
           x = (r_local-rci)*iw
           fcut_local=(x*(x*(20.0d0*x-70.0d0)+84.0d0)-35.0d0)&
                      *x**4+1.0d0
           temp1 =(x**3)*iw*(x*(x*(140.0d0*x-420.0d0)+420.0d0)-140.0d0)
        end if
      elseif(cutoff_type.eq.8)then
!!      fc=(r/Rc(r/Rc(r/Rc(315-70r/Rc)-540)+420)-126)*(r/Rc)^5+1
!!      fc'=[630(r/Rc)^4*(4(r/Rc)^3 - (r/Rc)^4 -6(r/Rc)^2 +4r/Rc -1)]/Rc
        if (r_local .le. rci) then
           fcut_local = 1.0d0
           temp1      = 0.0d0
        else
           x = (r_local-rci)*iw
           fcut_local=(x*(x*(x*(315.0d0-70.0d0*x)-540.0d0)+420.0d0)-126.0d0)&
                      *x**5+1.0d0
           temp1 =iw*(x**4)*(x*(x*((2520.0d0-630.d0*x)*x-3780.0d0)+2520.0d0)-630.0d0)
        end if
      else !'
        write(ounit,*)'ERROR: Unknown cutoff_type in calconefunction_para, please use value between 0 and 8 in input.nn'
        stop
      endif !'
!!
      return
      end 

