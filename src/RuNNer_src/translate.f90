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
!! - readonestructure.f90
!! - paircount.f90
!!
      subroutine translate(num_atoms,lattice,xyzstruct)
!!
      use globaloptions
!!
      implicit none
!!
      integer i
      integer num_atoms
      real*8 lattice(3,3)
      real*8 xyzstruct(3,max_num_atoms)
      real*8 axb(3),bxc(3),cxa(3)
      real*8 n1(3),n2(3),n3(3),n4(3),n5(3),n6(3)
      real*8 d1,d2,d3,d4,d5,d6
      real*8 distance1,distance2
      real*8 vorzeichen
      logical ltest
!!
!!      ltest=.true.
      ltest=.false.
!!
      if(ltest)then
        write(6,*)'lattice vectors (A)'
        do i=1,3
          write(6,'(3f14.6)')lattice(i,1),lattice(i,2),lattice(i,3)
        enddo
      endif
!!
!! calculation of the vector products
      axb(1)=lattice(1,2)*lattice(2,3)-lattice(1,3)*lattice(2,2)
      axb(2)=lattice(1,3)*lattice(2,1)-lattice(1,1)*lattice(2,3)
      axb(3)=lattice(1,1)*lattice(2,2)-lattice(1,2)*lattice(2,1)
!!
      bxc(1)=lattice(2,2)*lattice(3,3)-lattice(2,3)*lattice(3,2)
      bxc(2)=lattice(2,3)*lattice(3,1)-lattice(2,1)*lattice(3,3)
      bxc(3)=lattice(2,1)*lattice(3,2)-lattice(2,2)*lattice(3,1)
!!
      cxa(1)=lattice(3,2)*lattice(1,3)-lattice(3,3)*lattice(1,2)
      cxa(2)=lattice(3,3)*lattice(1,1)-lattice(3,1)*lattice(1,3)
      cxa(3)=lattice(3,1)*lattice(1,2)-lattice(3,2)*lattice(1,1)
!!
!! plane (0,0,0),(a,0,0),(0,b,0)
!! a x b through origin
      n1(1)=axb(1)/dsqrt(axb(1)**2+axb(2)**2+axb(3)**2)
      n1(2)=axb(2)/dsqrt(axb(1)**2+axb(2)**2+axb(3)**2)
      n1(3)=axb(3)/dsqrt(axb(1)**2+axb(2)**2+axb(3)**2)
      d1=0.d0 ! Plane through origin
!!
!! plane (0,0,c),(a,0,c),(0,b,c)
!! parallel to a x b through c
      n2(1)=axb(1)/dsqrt(axb(1)**2+axb(2)**2+axb(3)**2)
      n2(2)=axb(2)/dsqrt(axb(1)**2+axb(2)**2+axb(3)**2)
      n2(3)=axb(3)/dsqrt(axb(1)**2+axb(2)**2+axb(3)**2)
      d2=-n2(1)*lattice(3,1)-n2(2)*lattice(3,2)-n2(3)*lattice(3,3)
!!
!! plane (0,0,0),(0,b,0),(0,0,c)
!! b x c through origin
      n3(1)=bxc(1)/dsqrt(bxc(1)**2+bxc(2)**2+bxc(3)**2)
      n3(2)=bxc(2)/dsqrt(bxc(1)**2+bxc(2)**2+bxc(3)**2)
      n3(3)=bxc(3)/dsqrt(bxc(1)**2+bxc(2)**2+bxc(3)**2)
      d3=0.d0
!!
!! plane (a,0,0),(a,b,0),(a,0,c)
!! parallel to b x c trough a
      n4(1)=bxc(1)/dsqrt(bxc(1)**2+bxc(2)**2+bxc(3)**2)
      n4(2)=bxc(2)/dsqrt(bxc(1)**2+bxc(2)**2+bxc(3)**2)
      n4(3)=bxc(3)/dsqrt(bxc(1)**2+bxc(2)**2+bxc(3)**2)
      d4=-n4(1)*lattice(1,1)-n4(2)*lattice(1,2)-n4(3)*lattice(1,3)
!!
!! plane (0,0,0),(0,0,c),(a,0,0)
!! c x a through origin
      n5(1)=cxa(1)/dsqrt(cxa(1)**2+cxa(2)**2+cxa(3)**2)
      n5(2)=cxa(2)/dsqrt(cxa(1)**2+cxa(2)**2+cxa(3)**2)
      n5(3)=cxa(3)/dsqrt(cxa(1)**2+cxa(2)**2+cxa(3)**2)
      d5=0.0d0
!!
!! plane (0,b,0),(0,b,c),(a,b,0)
!! parallel to c x a through b
      n6(1)=cxa(1)/dsqrt(cxa(1)**2+cxa(2)**2+cxa(3)**2)
      n6(2)=cxa(2)/dsqrt(cxa(1)**2+cxa(2)**2+cxa(3)**2)
      n6(3)=cxa(3)/dsqrt(cxa(1)**2+cxa(2)**2+cxa(3)**2)
      d6=-n6(1)*lattice(2,1)-n6(2)*lattice(2,2)-n6(3)*lattice(2,3)
!!
      do i=1,num_atoms
 10     continue
        if(ltest) then
          write(6,'(a,3f14.6)')'atom',xyzstruct(1,i),xyzstruct(2,i),xyzstruct(3,i)
        endif
!!
!! calculate distance to plane 1
        distance1=n1(1)*xyzstruct(1,i)+n1(2)*xyzstruct(2,i)+n1(3)*xyzstruct(3,i)+d1
      if(ltest) write(6,'(a,f14.6)')'distance to plane 1 ',distance1
!!
!! calculate distance to plane 2
        distance2=n2(1)*xyzstruct(1,i)+n2(2)*xyzstruct(2,i)+n2(3)*xyzstruct(3,i)+d2
      if(ltest) write(6,'(a,f14.6)')'distance to plane 2 ',distance2
!!
        vorzeichen=n1(1)*lattice(3,1)+n1(2)*lattice(3,2)&
                 +n1(3)*lattice(3,3)
        vorzeichen=vorzeichen/abs(vorzeichen)
!!        write(6,*)'Vorzeichen',vorzeichen
!! for numerical stability 0.0d0 is not a good criterion
        if((distance1.lt.-0.00001d0).and.(distance2.lt.-0.00001d0))then
!!        if((distance1.lt.0.0d0).and.(distance2.lt.0.d0))then
          xyzstruct(1,i)=xyzstruct(1,i)+lattice(3,1)*vorzeichen
          xyzstruct(2,i)=xyzstruct(2,i)+lattice(3,2)*vorzeichen
          xyzstruct(3,i)=xyzstruct(3,i)+lattice(3,3)*vorzeichen
      if(ltest) write(6,*)'Point is below plane 1'
!!           write(6,*)'Atom is shifted into cell'
          goto 10
      elseif((distance1.gt.0.00001d0).and.(distance2.gt.0.00001d0))then
!!        elseif((distance1.gt.0.d0).and.(distance2.gt.0.d0))then
          xyzstruct(1,i)=xyzstruct(1,i)-lattice(3,1)*vorzeichen
          xyzstruct(2,i)=xyzstruct(2,i)-lattice(3,2)*vorzeichen
          xyzstruct(3,i)=xyzstruct(3,i)-lattice(3,3)*vorzeichen
      if(ltest) write(6,*)'Point is above plane 2'
!!           write(6,*)'Atom is shifted into cell'
          goto 10
        else
      if(ltest) write(6,*)'Point is between planes 1 and 2'
        endif
!!
!! calculate distance to plane 3
        distance1=n3(1)*xyzstruct(1,i)+n3(2)*xyzstruct(2,i)+n3(3)*xyzstruct(3,i)+d3
      if(ltest) write(6,'(a,f14.6)')'distance to plane 3 ',distance1
!!
!! calculate distance to plane 4
        distance2=n4(1)*xyzstruct(1,i)+n4(2)*xyzstruct(2,i)+n4(3)*xyzstruct(3,i)+d4
      if(ltest) write(6,'(a,f14.6)')'distance to plane 4 ',distance2
!!
        vorzeichen=n1(1)*lattice(3,1)+n1(2)*lattice(3,2)&
                 +n1(3)*lattice(3,3)
        vorzeichen=vorzeichen/abs(vorzeichen)
!!        write(6,*)'Vorzeichen',vorzeichen
        if((distance1.lt.-0.00001d0).and.(distance2.lt.-0.00001d0))then
!!        if((distance1.lt.0.0d0).and.(distance2.lt.0.d0))then
          xyzstruct(1,i)=xyzstruct(1,i)+lattice(1,1)*vorzeichen
          xyzstruct(2,i)=xyzstruct(2,i)+lattice(1,2)*vorzeichen
          xyzstruct(3,i)=xyzstruct(3,i)+lattice(1,3)*vorzeichen
      if(ltest) write(6,*)'Point is below plane 3'
!!           write(6,*)'Atom is shifted into cell'
          goto 10
      elseif((distance1.gt.0.00001d0).and.(distance2.gt.0.00001d0))then
!!        elseif((distance1.gt.0.d0).and.(distance2.gt.0.d0))then
          xyzstruct(1,i)=xyzstruct(1,i)-lattice(1,1)*vorzeichen
          xyzstruct(2,i)=xyzstruct(2,i)-lattice(1,2)*vorzeichen
          xyzstruct(3,i)=xyzstruct(3,i)-lattice(1,3)*vorzeichen
      if(ltest) write(6,*)'Point is above plane 4'
!!           write(6,*)'Atom is shifted into cell'
          goto 10
        else
      if(ltest) write(6,*)'Point is between planes 3 and 4'
        endif
!!
!! calculate distance to plane 5
        distance1=n5(1)*xyzstruct(1,i)+n5(2)*xyzstruct(2,i)+n5(3)*xyzstruct(3,i)+d5
      if(ltest) write(6,'(a,f14.6)')'distance to plane 5 ',distance1
!!
!! calculate distance to plane 6
        distance2=n6(1)*xyzstruct(1,i)+n6(2)*xyzstruct(2,i)+n6(3)*xyzstruct(3,i)+d6
      if(ltest) write(6,'(a,f14.6)')'distance to plane 6 ',distance2
!!
        vorzeichen=n1(1)*lattice(3,1)+n1(2)*lattice(3,2)&
                 +n1(3)*lattice(3,3)
        vorzeichen=vorzeichen/abs(vorzeichen)
!!        write(6,*)'Vorzeichen',vorzeichen
        if((distance1.lt.-0.00001d0).and.(distance2.lt.-0.00001d0))then
!!        if((distance1.lt.0.0d0).and.(distance2.lt.0.d0))then
          xyzstruct(1,i)=xyzstruct(1,i)+lattice(2,1)*vorzeichen
          xyzstruct(2,i)=xyzstruct(2,i)+lattice(2,2)*vorzeichen
          xyzstruct(3,i)=xyzstruct(3,i)+lattice(2,3)*vorzeichen
      if(ltest) write(6,*)'Point is below plane 5'
!!           write(6,*)'Atom is shifted into cell'
          goto 10
      elseif((distance1.gt.0.00001d0).and.(distance2.gt.0.00001d0))then
!!        elseif((distance1.gt.0.d0).and.(distance2.gt.0.d0))then
          xyzstruct(1,i)=xyzstruct(1,i)-lattice(2,1)*vorzeichen
          xyzstruct(2,i)=xyzstruct(2,i)-lattice(2,2)*vorzeichen
          xyzstruct(3,i)=xyzstruct(3,i)-lattice(2,3)*vorzeichen
      if(ltest) write(6,*)'Point is above plane 6'
!!           write(6,*)'Atom is shifted into cell'
          goto 10
        else
      if(ltest) write(6,*)'Point is between planes 5 and 6'
        endif
!!
        if(ltest) then
          write(6,'(a,3f14.6)')'atom at end ',xyzstruct(1,i),xyzstruct(2,i),xyzstruct(3,i)
        endif
!!
      enddo
!!
      return
      end
