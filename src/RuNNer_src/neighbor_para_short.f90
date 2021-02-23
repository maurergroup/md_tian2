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
!######################################################################!! TODO: check if this also might work for 
!!    * ewaldreal_mode3
!!    * ewaldreal_para
!! If this is the case, replace neighbor_para with this one
!! called by: 
!!
      subroutine neighbor_para_short(n_start,n_end,&
                 num_atoms,zelem,&
                 lsta,lstb,lstc,lste,&
                 cutoff,lattice,xyz,lperiodic)
!!
      use globaloptions
      use fileunits 
!!
      implicit none
!!
      integer num_atoms 
      integer i,j
      integer indc
      integer lsta(2,max_num_atoms)
      integer lstc(listdim)
      integer lste(listdim)
      integer na,nb,nc
      integer n1,n2,n3
      integer zelem(max_num_atoms)
      integer icount,jcount                  ! internal
      integer n_start                        ! in
      integer n_end                          ! in
!!
      real*8 lstb(listdim,4)
      real*8 xtemp,ytemp,ztemp
      real*8 xyz(3,max_num_atoms)
      real*8 lattice(3,3)
      real*8 cutoff,rr
      real*8 xrel1,xrel2,xrel3
      real*8 axb(3),axc(3),bxc(3) ! vector products
      real*8 absaxb,absaxc,absbxc
      real*8 proja,projb,projc
      logical lperiodic
!!
      if(lperiodic) then
!! determine how often we have to multiply the cell in which direction
      na = 0
      nb = 0
      nc = 0
!! calculate the norm vectors of the 3 planes
      axb(1)=lattice(1,2)*lattice(2,3)-lattice(1,3)*lattice(2,2)
      axb(2)=lattice(1,3)*lattice(2,1)-lattice(1,1)*lattice(2,3)
      axb(3)=lattice(1,1)*lattice(2,2)-lattice(1,2)*lattice(2,1)
      absaxb=sqrt(axb(1)**2+axb(2)**2+axb(3)**2)
      axb(1)=axb(1)/absaxb
      axb(2)=axb(2)/absaxb
      axb(3)=axb(3)/absaxb
      axc(1)=lattice(1,2)*lattice(3,3)-lattice(1,3)*lattice(3,2)
      axc(2)=lattice(1,3)*lattice(3,1)-lattice(1,1)*lattice(3,3)
      axc(3)=lattice(1,1)*lattice(3,2)-lattice(1,2)*lattice(3,1)
      absaxc=sqrt(axc(1)**2+axc(2)**2+axc(3)**2)
      axc(1)=axc(1)/absaxc
      axc(2)=axc(2)/absaxc
      axc(3)=axc(3)/absaxc
      bxc(1)=lattice(2,2)*lattice(3,3)-lattice(2,3)*lattice(3,2)
      bxc(2)=lattice(2,3)*lattice(3,1)-lattice(2,1)*lattice(3,3)
      bxc(3)=lattice(2,1)*lattice(3,2)-lattice(2,2)*lattice(3,1)
      absbxc=sqrt(bxc(1)**2+bxc(2)**2+bxc(3)**2)
      bxc(1)=bxc(1)/absbxc
      bxc(2)=bxc(2)/absbxc
      bxc(3)=bxc(3)/absbxc
!! calculate the projections
      proja=lattice(1,1)*bxc(1)+lattice(1,2)*bxc(2)+lattice(1,3)*bxc(3)
      projb=lattice(2,1)*axc(1)+lattice(2,2)*axc(2)+lattice(2,3)*axc(3)
      projc=lattice(3,1)*axb(1)+lattice(3,2)*axb(2)+lattice(3,3)*axb(3)
      proja=abs(proja)
      projb=abs(projb)
      projc=abs(projc)
!! determine how often we have to multiply the cell in which direction
      na = 0
      nb = 0
      nc = 0
      do while(dble(na)*proja.le.cutoff)
        na=na+1
      enddo
      do while(dble(nb)*projb.le.cutoff)
        nb=nb+1
      enddo
      do while(dble(nc)*projc.le.cutoff)
        nc=nc+1
      enddo
!!      write(*,*)'proja,projb,projc',proja,projb,projc
!!      write(*,*)'na,nb,nc',na,nb,nc
      endif
!!
!! determine lsta,lstb,lstc,lste
      icount=0
      if(lperiodic)then
      indc=0
      jcount=n_start-1
!! loop over atoms for which we want the neighbors
!!      do i=1,num_atoms
      do i=n_start,n_end
        jcount=jcount+1
        lsta(1,jcount)=indc+1
!! loop over all potentially neighboring atoms
        do j=1,num_atoms
          do n1=-na,na
            do n2=-nb,nb
              do n3=-nc,nc
!! avoid interaction of atom with itself
         if((n1.eq.0).and.(n2.eq.0).and.(n3.eq.0).and.(i.eq.j)) then
         else
           xtemp=xyz(1,j)&
      +dble(n1)*lattice(1,1)+dble(n2)*lattice(2,1)+dble(n3)*lattice(3,1)
           ytemp=xyz(2,j)&
      +dble(n1)*lattice(1,2)+dble(n2)*lattice(2,2)+dble(n3)*lattice(3,2)
           ztemp=xyz(3,j)&
      +dble(n1)*lattice(1,3)+dble(n2)*lattice(2,3)+dble(n3)*lattice(3,3)
!!           write(*,'(i4,3f14.6)')j,xtemp,ytemp,ztemp
           xrel1=xtemp-xyz(1,i)
           xrel2=ytemp-xyz(2,i)
           xrel3=ztemp-xyz(3,i)
           rr=xrel1**2+xrel2**2+xrel3**2
           if(rr.le.cutoff**2) then
             indc=indc+1
             if(indc.gt.listdim)then
                write(*,*)'Error: neighbor_para_short: redimension lstb'
                write(*,*)'This error is often a consequence of wrong units in the file'
                write(*,*)'input.data (must be Bohr for postions and lattice vectors'
                stop
             endif
             lstb(indc,1)=xtemp
             lstb(indc,2)=ytemp
             lstb(indc,3)=ztemp
             lstb(indc,4)=sqrt(rr)
!!             write(*,'(3i4,4f14.6)')indc,i,j,xtemp,ytemp,ztemp,sqrt(rr)
             lstc(indc)=j ! identification of atom
             lste(indc)=zelem(j) ! nuclear charge of atom
           endif
         endif
              enddo ! n3
            enddo ! n2
          enddo ! n1
        enddo ! j
        lsta(2,jcount)=indc
!!        write(*,*)'atom',i,'from',lsta(1,i),'to',lsta(2,i)
!!        write(*,*)'atom',i,'from',lsta(1,jcount),'to',lsta(2,jcount)
      enddo ! i
!!
      else ! .not. lperiodic
      indc=0
      jcount=n_start-1
!! loop over atoms for which we want the neighbors
!!    do i=1,num_atoms

      do i=n_start,n_end
        jcount=jcount+1

        lsta(1,jcount)=indc+1
!! loop over all potentially neighboring atoms
        do j=1,num_atoms
         if(i.eq.j) then
         else
           xtemp=xyz(1,j)
           ytemp=xyz(2,j)
           ztemp=xyz(3,j)
           xrel1=xtemp-xyz(1,i)
           xrel2=ytemp-xyz(2,i)
           xrel3=ztemp-xyz(3,i)
           rr=xrel1**2+xrel2**2+xrel3**2
           if(rr.le.cutoff**2) then
             indc=indc+1
             if(indc.gt.listdim)then
                write(*,*)'Error: neighbor_para_short: redimension lstb'
                write(*,*)'This error is often a consequence of wrong units in the file'
                write(*,*)'input.data (must be Bohr for postions and lattice vectors'
                stop
             endif
             lstb(indc,1)=xtemp
             lstb(indc,2)=ytemp
             lstb(indc,3)=ztemp
             lstb(indc,4)=sqrt(rr)
             lstc(indc)=j ! identification of atom
             lste(indc)=zelem(j) ! nuclear charge of atom
           endif
         endif
        enddo ! j
        lsta(2,jcount)=indc
!!        write(*,*)'atom',i,'from',lsta(1,i),'to',lsta(2,i)
!!        write(*,*)'atom',i,'from',lsta(1,jcount),'to',lsta(2,jcount)
      enddo ! i
      endif
!!
      return
      end
