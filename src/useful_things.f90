!############################################################################
! This routine is part of
! md_tian2 (Molecular Dynamics Tian Xia 2)
! (c) 2014-2020 Daniel J. Auerbach, Svenja M. Janke, Marvin Kammler,
!               Sascha Kandratsenka, Sebastian Wille
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

! Contains useful math routine
module useful_things

    use constants

    implicit none


    interface normal_deviate
        module procedure normal_deviate_0d
        module procedure normal_deviate_1d
        module procedure normal_deviate_2d
        module procedure normal_deviate_3d
    end interface

    public :: normal_deviate

contains



    function ran1()  !returns random number between 0 - 1
        real(8) ran1,x
        call random_number(x) ! built in fortran 90 random number function
        ran1 = x
    end function ran1

!    function normal(mean,sigma) !returns a normal distribution
!        real(8) normal,tmp
!        real(8) mean,sigma   ! Sigma is the velocity we want to achieve
!        integer flag
!        real(8) fac,gsave,rsq,r1,r2
!        save flag,gsave
!        data flag /0/
!        if (flag.eq.0) then
!            rsq=2.0d0
!            do while(rsq.ge.1.0d0.or.rsq.eq.0.0d0) ! new from for do
!                r1=2.0d0*ran1()-1.0d0
!                r2=2.0d0*ran1()-1.0d0
!                rsq=r1*r1+r2*r2
!            enddo
!            fac=sqrt(-2.0d0*log(rsq)/rsq)
!            gsave=r1*fac        ! shouldn't those two values be below zero?
!            tmp=r2*fac          !
!            flag=1
!        else
!            tmp=gsave
!            flag=0
!        endif
!        normal=tmp*sigma+mean
!        return
!    end function normal

    function ranx(nran, seed, switch) result(rnd)

        implicit none

        integer  :: nran ! set which rng to use
        real(dp) :: rnd
        integer*8, intent(inout) :: seed ! the seed will be the trajectory id
        integer :: switch, i
        integer, dimension(1) :: putseed

        integer :: randk


        select case (nran)

            case (1)

                select case (switch)

                    case (0) ! fixed seed, rotation according to traj id

                        ! size of seed for random number generator
                        randk=size(randseed)
                        call random_seed(size=randk)
                        call random_seed(put=randseed)

                        ! rotate rng
                        do i = 1, 100*seed
                            call random_number(rnd)
                        end do

                    case (1)

                        call random_number(rnd)

                    case default

                        print *, 'Unknown value for switch in ranx function'
                        stop
                
                end select

            case (2) ! seed is traj id, without rotation

                select case (switch)

                    case (0)

                        putseed = (/seed/)
                        call random_seed(put=putseed)

                        !do i = 1, 100*seed
                            !call random_number(rnd)
                        !end do

                    case (1)

                        call random_number(rnd)

                    case default

                        print *, 'Unknown value for switch in ranx function'
                        stop

                end select
                
            case (3) ! rng with traj id as seed

                select case (switch)

                    case (0)

                        print *, 'Warning: Reproducability of random numbers not guarenteed when recalculate single trajectory!'
                        ! just let it pass

                    case (1)

                        rnd = ran3(seed)

                    case default
                         print *, 'Unknown value for switch in ranx function'
                         stop

                end select

            case default

                print *, 'Error: Unknown random number generator type in ranx function'
                stop

        end select

    end function ranx

    function ran3(seed) result(rnd)

        implicit none

        real(dp) :: rnd
        integer*8, intent(inout) :: seed

        rnd = xorshift64star(seed)

    end function ran3

    ! xorshift64*
    ! took implementation of the algorithm in the Spire package as reference (MIT license)
    ! https://github.com/typelevel/spire/blob/master/extras/src/main/scala/spire/random/rng/XorShift64Star.scala
    ! has a period of 2**64-1, translated to fortran by JAF
    function xorshift64star(state) result(rnd)

        implicit none

        real*8 :: rnd
        integer*8, intent(inout) :: state

        integer*8 :: b(4)
        integer*8, parameter :: f(4) = [56605, 20332, 62609, 9541] ! 2685821657736338717 in base 2**16
        integer*8 :: m(4)

        state = ieor(state, ishft(state, -12))
        state = ieor(state, ishft(state,  25))
        state = ieor(state, ishft(state, -27))

        ! fortran has no unsigned ints.
        ! below would be how we convert signed to unsigned (twos complement)
        ! if (x < 0) then
        !   x = 2_16**64 + state
        ! else
        !   x = state
        ! end if
        ! doing the modulo is equivalent however
        ! BUT: no 128 bit int support from intel compiler
        ! we therefore use some trickery by transforming the number into base 2**16
        ! the version below would work with gnu compiler
        ! x = modulo(state * 2685821657736338717_16, 2_16**64)
        ! rnd = real(x, 8) / 2_16**64
        ! write(*,*) rnd
        ! this is the version for intel

        call toBase16(state, b)
        call multBase16(b, f, m)
        ! call fromBase16(m, xx)
        ! to return a double we divide by 2^64
        rnd = m(4) / 2.d0**16 + m(3) / 2.d0**32 + m(2) / 2.d0**48 + m(1) / 2.d0**64

    end function xorshift64star

    ! becaus we cannot use 128 bit ints in ifort we represent our number using 4 ints in base 2**16
    ! x = sum_i b(i) * 2**(16*(i-1))
    ! treats x as unsigned int
    subroutine toBase16(x, b)

        implicit none

        integer*8, intent(in) :: x
        integer*8, intent(out) :: b(4)
        integer :: i
        integer*8 :: t

        b(:) = 0
        if (x<0) then
            t = (9223372036854775807_8 + x) + 1
            b(4) = 2_8**15
        else
            t = x
        end if

        do i=1,4
            b(i) = b(i) + modulo(t, 2_8**16)
            t = t / 2_8**16
        end do

    end subroutine

    ! also treats x as unsigned int
    subroutine fromBase16(b, x) ! not needed, but left in for support information

        implicit none

        integer*8, intent(in) :: b(4)
        integer*8, intent(out) :: x
        integer*8 :: t
        integer :: i

        x = 0
        do i=1,3
            x = x + b(i) * 2_8**(16*(i-1))
        end do
        x = x + modulo(b(4), 2_8**16) * 2_8**(16*3)

    end subroutine

    ! does multiplication mod 2**64 of two numbers in the base of 2**16
    subroutine multBase16(a, b, c)

        implicit none

        integer*8, intent(in) :: a(4), b(4)
        integer*8, intent(out) :: c(4)

        integer :: i,j

        c(:) = 0
        do i=1,4
            do j=1,5-i
                c(i+j-1) = c(i+j-1) + a(i) * b(j)
            end do
        end do

        do i=1,3
            c(i+1) = c(i+1) + c(i) / 2_8**16
            c(i) = modulo(c(i), 2_8**16)
        end do
        c(4) = modulo(c(4), 2_8**16)

    end subroutine multBase16

    subroutine normal_deviate_0d(mu, sigma, nrml_dvt)

        real(dp), intent(in)  :: mu, sigma
        real(dp), intent(out) :: nrml_dvt
        real(dp)              :: rnd1, rnd2

        call random_number(rnd1)
        call random_number(rnd2)
        nrml_dvt = sqrt(-2*log(rnd1)) * cos(2*pi*rnd2) * sigma + mu

    end subroutine normal_deviate_0d


    subroutine normal_deviate_1d(mu, sigma, nrml_dvt)

        real(dp), intent(in)                :: mu, sigma
        real(dp), intent(out)               :: nrml_dvt(:)
        real(dp), dimension(size(nrml_dvt)) :: rnd1, rnd2

        call random_number(rnd1)
        call random_number(rnd2)
        nrml_dvt = sqrt(-2*log(rnd1)) * cos(2*pi*rnd2) * sigma + mu

    end subroutine normal_deviate_1d


    subroutine normal_deviate_2d(mu, sigma, nrml_dvt)

        real(dp), intent(in)   :: mu, sigma
        real(dp), intent(out)  :: nrml_dvt(:,:)
        real(dp), dimension(size(nrml_dvt, dim=1), &
                            size(nrml_dvt, dim=2)) :: rnd1, rnd2

        call random_number(rnd1)
        call random_number(rnd2)
        nrml_dvt = sqrt(-2*log(rnd1)) * cos(2*pi*rnd2) * sigma + mu

    end subroutine normal_deviate_2d


    subroutine normal_deviate_3d(mu, sigma, nrml_dvt)

        real(dp), intent(in)   :: mu, sigma
        real(dp), intent(out)  :: nrml_dvt(:,:,:)
        real(dp), dimension(size(nrml_dvt, dim=1), &
                            size(nrml_dvt, dim=2), &
                            size(nrml_dvt, dim=3)) :: rnd1, rnd2

        call random_number(rnd1)
        call random_number(rnd2)
        nrml_dvt = sqrt(-2*log(rnd1)) * cos(2*pi*rnd2) * sigma + mu

    end subroutine normal_deviate_3d




    subroutine lower_case(str)
        character(*), intent(in out) :: str
        integer :: i

        do i = 1, len(str)
            select case(str(i:i))
                case("A":"Z")
                    str(i:i) = achar(iachar(str(i:i))+32)
            end select
        end do
    end subroutine lower_case




    subroutine split_string ( line, words, nw )
        character(*), intent(in)  :: line
        character(*), intent(out) :: words(:)
        integer,      intent(out) :: nw
        character(len(words)) :: buf( size(words) )
        integer :: i, ios

        nw = 0 ; words(:) = ""

        do i = 1, size(words)
            read( line, *, iostat=ios ) buf( 1 : i )
            if ( ios /= 0 ) exit
            nw = i
            words( 1 : nw ) = buf( 1 : nw )
        enddo

    endsubroutine


    subroutine norm_dist(vec1, vec2, length, norm)
        !
        ! Purpose: normalised distance between 2 vectors
        !

        integer :: length
        real(8), dimension(length) :: vec1, vec2
        real(8) :: norm, n1, n2

        norm = dot_product(vec1 - vec2, vec1 - vec2)
        n1   = dot_product(vec1,vec1)
        n2   = dot_product(vec2,vec2)

        n1 = max(n1,n2)
        if (n1 .eq. 0.0d0) then
            norm = 0.0d0
        else
            norm=Sqrt(norm/n1)
        end if

    end subroutine norm_dist



    subroutine pbc_dist(a, b, cmat, cimat, r)
            !
            ! Purpose: Distance between atoms a and b
            !          with taking into account the periodic boundary conditions
            !

        real(8), dimension(3),   intent(in)  :: a, b
        real(8), dimension(3,3), intent(in)  :: cmat, cimat
        real(8),                 intent(out) :: r

        real(8), dimension(3) :: r3temp


        ! Applying PBCs
        r3temp = b - a   ! distance vector from a to b
        r3temp = matmul(cimat, r3temp)   ! transform to direct coordinates

        r3temp(1) = r3temp(1) - Anint(r3temp(1))! imaging
        r3temp(2) = r3temp(2) - Anint(r3temp(2))
        r3temp(3) = r3temp(3) - Anint(r3temp(3))
        r3temp    = matmul(cmat, r3temp)    ! back to cartesian coordinates

        r =  sqrt(sum(r3temp*r3temp))               ! distance

    end subroutine pbc_dist



    function lines_in_file(lunit, file_name)
        !
        ! Purpose: Count the number of lines in file 'file_name'.
        !          This allows for run-time determination of number sample points.
        !
        integer, intent(in)         :: lunit
        character(*), intent(in)    :: file_name
        integer                     :: ios, lines_in_file

        lines_in_file = 0
        open(lunit, file=file_name)
        do
            read(lunit,*, IOSTAT=ios)
            if (ios /= 0) exit
            lines_in_file = lines_in_file + 1
        end do
        close(lunit)
    end function


    logical function file_exists(fname) result(exists)

        character(len=*), intent(in) :: fname

        inquire(file=fname, exist=exists)

    end function file_exists



    logical function dir_exists(fname) result(exists)

        character(len=*), intent(in) :: fname

        inquire(directory=fname, exist=exists)

    end function dir_exists



    function invert_matrix(A) result(B)
        !! Performs a direct calculation of the inverse of a 3×3 matrix.
        real(dp), intent(in) :: A(3,3)   !! Matrix
        real(dp)             :: B(3,3)   !! Inverse matrix
        real(dp)             :: detinv, det

        det = (A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)&
            - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)&
            + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))

        if (det < tolerance) stop "Error in invert_matrix: matrix is singular"

        ! Calculate the inverse determinant of the matrix
        detinv = 1/det

        ! Calculate the inverse of the matrix
        B(1,1) = +detinv * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
        B(2,1) = -detinv * (A(2,1)*A(3,3) - A(2,3)*A(3,1))
        B(3,1) = +detinv * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
        B(1,2) = -detinv * (A(1,2)*A(3,3) - A(1,3)*A(3,2))
        B(2,2) = +detinv * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
        B(3,2) = -detinv * (A(1,1)*A(3,2) - A(1,2)*A(3,1))
        B(1,3) = +detinv * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
        B(2,3) = -detinv * (A(1,1)*A(2,3) - A(1,3)*A(2,1))
        B(3,3) = +detinv * (A(1,1)*A(2,2) - A(1,2)*A(2,1))
    end function



    function cro_pro(r1, r2)

        real(dp), intent(in) :: r1(:,:), r2(:,:)
        real(dp) :: cro_pro(3, size(r1, dim=2))

        if (size(r1, dim=1) /= 3 .or. size(r2, dim=1) /= 3) then
            print *, "First dimension in cross product function must be 3"
            stop

        else if (size(r1, dim=2) /= size(r2, dim=2)) then
            print *, "Arrays in cross product function differ in their 2nd dimension"
            stop

        end if

        cro_pro(1,:) = r1(2,:)*r2(3,:) - r1(3,:)*r2(2,:)
        cro_pro(2,:) = r1(3,:)*r2(1,:) - r1(1,:)*r2(3,:)
        cro_pro(3,:) = r1(1,:)*r2(2,:) - r1(2,:)*r2(1,:)

    end function cro_pro



    subroutine check_and_set_fit(vary, idx1, idx2, arr)

        character(len=3), intent(in) :: vary
        integer, intent(in) :: idx1, idx2
        logical, intent(out) :: arr(:,:)

        if (vary == "fit") then
            arr(idx1, idx2) = .True.
            arr(idx2, idx1) = .True.
        end if

    end subroutine check_and_set_fit



    subroutine timestamp(out_unit)

        integer, intent(in) :: out_unit
        integer :: val(8)
        character(len=*), parameter :: fstring = "(a, i4, a, i2.2, a, i2.2, a, i2.2, a, i2.2, a, i2.2, a2$)"

        call date_and_time(values=val)

        write(out_unit, fstring) "[", val(1), "-", val(2), "-",val(3), " - ",val(5), ".",val(6), ":", val(7), "] "

    end subroutine timestamp

end module useful_things
