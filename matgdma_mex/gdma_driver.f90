MODULE gdma_driver

TYPE gdma_input
    real(kind(1d0)) :: bigexp ! threshold for old/new algorithm
    integer, allocatable :: limit(:), & ! maximum order of multipole output, numSites by 1
                            shellNprims(:), &
                            shell2atom(:), &
                            shellNfuncs(:)
    real(kind(1d0)), allocatable :: nucleiCharges(:), &
                                    xyzSites(:,:), & ! 3 by numSites, row-major
                                    primExps(:), &
                                    primCoefs(:), &
                                    density(:,:)
END TYPE

CONTAINS

SUBROUTINE gdma_driver_routine(q_out, input_args)

!  Distributed Multipole Analysis for Gaussian Wavefunctions
!
!  Copyright (C) 2005-14  Anthony J. Stone
!
!  This program is free software; you can redistribute it and/or
!  modify it under the terms of the GNU General Public License
!  as published by the Free Software Foundation; either version 2
!  of the License, or (at your option) any later version.
!  
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!  
!  You should have received a copy of the GNU General Public License
!  along with this program; if not, write to
!  the Free Software Foundation, Inc., 51 Franklin Street,
!  Fifth Floor, Boston, MA 02110-1301, USA.

USE input

USE dma
USE atom_grids, ONLY: debug_grid => debug
USE timing, ONLY: start_timer, timer, time_and_date
IMPLICIT NONE

INTEGER, PARAMETER :: dp=kind(1d0)

CHARACTER(LEN=100) :: file
CHARACTER(LEN=80) :: buffer
CHARACTER(LEN=20) :: key
CHARACTER(LEN=8) :: whichd="SCF"
CHARACTER(LEN=24) :: datestring

!  Maximum number of sites is number of atoms + nextra
INTEGER :: nextra=16
INTEGER :: ncoorb, maxl, cmax, nprim, nx, num, n, ich, mul, nbf
INTEGER, ALLOCATABLE :: shell_type(:), shell_nfuncs(:)
INTEGER :: i, j, k, kp=0
LOGICAL :: eof, fchk, first, ok=.false.

REAL(dp), ALLOCATABLE :: densty(:,:), dtri(:), input_density(:,:)
INTEGER :: ir=5 ! Input stream

LOGICAL :: verbose=.false., debug(0:2)=.false.

INTEGER :: aok
REAL(dp) :: e, rt3v2, td(5,6), tf(7,10), tg(9,15)
REAL(dp), ALLOCATABLE :: temp(:,:), primCoefs(:)
REAL(dp), PARAMETER :: PI=3.14159265358979d0
REAL(dp), PARAMETER :: rt2=1.4142135623731d0,                          &
    rt3=1.73205080756888d0, rt5=2.23606797749979d0, rt7=2.64575131106459d0
REAL(dp), PARAMETER :: rt10=rt2*rt5, rt14=rt2*rt7, rt21=rt3*rt7,       &
    rt35=rt5*rt7
INTEGER, PARAMETER :: v400=1, v040=2, v004=3, v310=4, v301=5,          &
    v130=6, v031=7, v103=8, v013=9, v220=10, v202=11, v022=12,         &
    v211=13, v121=14, v112=15

! input arguments
type(gdma_input), intent(in) :: input_args

! output arguments
real(dp), allocatable, intent(out) :: q_out(:,:)

fchk=.false.
first=.true.


!! Start input

nat = size(input_args%limit)
! Deep copy zan = input_args%nucleiCharges
if(allocated(zan)) deallocate(zan)
allocate(zan(nat))
do i = 1, size(input_args%nucleiCharges)
    zan(i) = input_args%nucleiCharges(i)
end do

! Deep copy c = input_args%xyzSites ; cartesian coordinates of sites
if(allocated(c)) deallocate(c)
allocate(c(3, nat))
do i = 1, size(input_args%xyzSites, 1)
    do j = 1, size(input_args%xyzSites, 2)
        c(i,j) = input_args%xyzSites(i,j)
    end do
end do

! Deep copy user_given_limit = input_args%limit ; maximum output multipole order
if(allocated(user_given_limit)) deallocate(user_given_limit)
allocate(user_given_limit(nat))
do i = 1, size(input_args%limit)
    user_given_limit(i) = input_args%limit(i)
end do

!~ ich = XXX ! total charge
!~ mul = XXX ! total multiplicity; I guess they just affect output print

nshell = size(input_args%shellNfuncs)
! Deep copy kng = input_args%shellNprims
if(allocated(kng)) deallocate(kng)
allocate(kng(nshell))
do i = 1, nshell
    kng(i) = input_args%shellNprims(i)
end do

! Deep copy katom = input_args%shell2atom
if(allocated(katom)) deallocate(katom)
allocate(katom(nshell+1))
katom = 0
do i = 1, nshell
    katom(i) = input_args%shell2atom(i)
end do

! Deep copy shell_nfuncs = input_args%shellNfuncs
if(allocated(shell_nfuncs)) deallocate(shell_nfuncs)
allocate(shell_nfuncs(nshell))
do i = 1, nshell
    shell_nfuncs(i) = input_args%shellNfuncs(i)
end do

! Determine shell_type according to shell_nfuncs
if(allocated(shell_type)) deallocate(shell_type)
allocate(shell_type(nshell))
do i = 1, nshell
    select case (shell_nfuncs(i))
        case (1)
            shell_type(i) = 0
        case (3)
            shell_type(i) = 1
        case (5, 6)
            shell_type(i) = 2
        case (7, 10)
            shell_type(i) = 3
        case (9, 15)
            shell_type(i) = 4
        case default
            shell_type(i) = -1
    end select
end do

nprim = sum(input_args%shellNprims)
! Deep copy ex = input_args%primExps
if(allocated(ex)) deallocate(ex)
allocate(ex(nprim))
do i = 1, nprim
    ex(i) = input_args%primExps(i)
end do

! Deep copy primCoefs = input_args%primCoefs
if(allocated(primCoefs)) deallocate(primCoefs)
allocate(primCoefs(nprim))
do i = 1, nprim
    primCoefs(i) = input_args%primCoefs(i)
end do

nbf = sum(shell_nfuncs)
! Deep copy input_density = input_args%density
if(allocated(input_density)) deallocate(input_density)
allocate(input_density(nbf, nbf))
do i = 1, nbf
    do j = 1, nbf
        input_density(i,j) = input_args%density(i,j)
    end do
end do

bigexp = input_args%bigexp ! threshold for old/new algorithm

!! Done input


lmax = maxval(user_given_limit)
num = size(input_density, 1)

maxcen=nat
maxs=maxcen+nextra

allocate (kstart(nshell), ktype(nshell),          &
    kloc(nshell), kmin(nshell), kmax(nshell),stat=aok)
kstart = 0
ktype = 0
kloc = 0
kmin = 0
kmax = 0

allocate(cs(nprim), cp(nprim), stat=aok)
cs=0d0; cp=0d0

n=0
do i=1,nshell
  select case(shell_type(i))
  case(0)
    n=n+1
  case(1)
    n=n+3
  case(2)
    n=n+6
  case(3)
    n=n+10
  case(4)
    n=n+15
  end select
end do
maxbfn=n
allocate(iax(n+1), stat=aok)

k=1
do i=1,nshell
  kstart(i)=k
  k=k+kng(i)
end do

k = 0
cp = 0.0
cs = 0.0
do i=1,nshell
  do j=kstart(i),kstart(i)+kng(i)-1
    k = k + 1
    e = primCoefs(k)
    if (shell_type(i) .eq. 1) then
      cp(j)=e
    else
      cs(j)=e
    end if
  end do
end do
nx = n*(n+1)/2

allocate(densty(n,n),temp(n,n))
do i = 1, size(input_density, 1)
    do j = 1, size(input_density, 2)
        densty(i, j) = input_density(i, j)
    end do
end do
!~ deallocate(input_density)
!~ densty(1:size(input_density,1), 1:size(input_density,2)) = input_density

!  We use unnormalized primitive functions, so we transfer the
!  normalising factor to the contraction coefficients. This is
!  the factor for z^n exp(-e*r^2). General formula is
!  (4e)^(n/2).(2e/pi)^{3/4}/sqrt{(2n-1)!!}
do i=1,nshell
  do j=kstart(i),kstart(i)+kng(i)-1
    e=ex(j)
    select case(abs(shell_type(i)))
    case(0,1)
      cs(j)=cs(j)*sqrt(sqrt((2d0*e/pi)**3))
      cp(j)=cp(j)*sqrt(4d0*e*sqrt((2d0*e/pi)**3))
    case(5) ! h shell
      cs(j)=cs(j)*(4d0*e)**2*sqrt(4d0*e*sqrt((2d0*e/pi)**3)/945d0)
    case(4) ! g shell
      cs(j)=cs(j)*(4d0*e)**2*sqrt(sqrt((2d0*e/pi)**3)/105d0)
    case(3) ! f shell
      cs(j)=cs(j)*4d0*e*sqrt(4d0*e*sqrt((2d0*e/pi)**3)/15d0)
    case(2) ! d shell
      cs(j)=cs(j)*4d0*e*sqrt(sqrt((2d0*e/pi)**3)/3d0)
    end select
  end do
end do


!  Conversion from normalised spherical form to normalised Cartesian
!  Schlegel & Frisch, IJQC (1995) 54, 83-87.
rt3v2=rt3/2d0
td(1,:)=(/-0.5d0, -0.5d0, 1d0, 0d0, 0d0, 0d0/)
td(2,:)=(/0d0,    0d0,    0d0, 0d0, 1d0, 0d0/)
td(3,:)=(/0d0,    0d0,    0d0, 0d0, 0d0, 1d0/)
td(4,:)=(/rt3v2,  -rt3v2, 0d0, 0d0, 0d0, 0d0/)
td(5,:)=(/0d0,    0d0,    0d0, 1d0, 0d0, 0d0/)

!  f functions
!   1   2   3   4   5   6   7   8   9   10
!  xxx yyy zzz xxy xxz xyy yyz xzz yzz xyz
tf(:,:)=0d0
! 30
tf(1,3)=1d0; tf(1,5)=-1.5d0/sqrt(5d0); tf(1,7)=-1.5d0/sqrt(5d0)
! 31c ( F+1 in Gaussian notation )
tf(2,1)=-sqrt(3d0/8d0); tf(2,6)=-sqrt(1.2d0)/4d0; tf(2,8)=sqrt(1.2d0)
! 31s ( F-1 )
tf(3,2)=-sqrt(3d0/8d0); tf(3,4)=-sqrt(1.2d0)/4d0; tf(3,9)=sqrt(1.2d0)
! 32c ( F+2 )
tf(4,5)=sqrt(0.75d0); tf(4,7)=-sqrt(0.75d0)
! 32s
tf(5,10)=1d0
! 33c
tf(6,1)=sqrt(10d0)/4d0; tf(6,6)=-0.75d0*sqrt(2d0)
! 33s
tf(7,2)=-sqrt(10d0)/4d0; tf(7,4)=0.75d0*sqrt(2d0)

!  g functions
!   1    2    3    4    5    6    7    8    9   10   11   12   13   14   15
! xxxx yyyy zzzz xxxy xxxz xyyy yyyz xzzz yzzz xxyy xxzz yyzz xxyz xyyz xyzz
tg=0d0
!  40
tg(1,v400)=0.375d0; tg(1,v040)=0.375d0; tg(1,v004)=1d0
tg(1,v220)=0.75d0*rt3/rt35; tg(1,v202)=-3d0*rt3/rt35; tg(1,v022)=-3d0*rt3/rt35
!  41c
tg(2,v103)=rt10/rt7; tg(2,v301)=-0.75d0*rt10/rt7; tg(2,v121)=-0.75d0*rt2/rt7
!  41s
tg(3,v013)=rt10/rt7; tg(3,v031)=-0.75d0*rt10/rt7; tg(3,v211)=-0.75d0*rt2/rt7
!  42c
tg(4,v202)=1.5d0*rt3/rt7; tg(4,v022)=-1.5d0*rt3/rt7; tg(4,v400)=-rt5/4d0; tg(4,v040)=rt5/4d0
!  42s
tg(5,v112)=3d0/rt7; tg(5,v310)=-rt5/(2d0*rt7); tg(5,v130)=-rt5/(2d0*rt7)
!  43c
tg(6,v301)=rt10/4d0; tg(6,v121)=-0.75d0*rt2
!  43s
tg(7,v031)=-rt10/4d0; tg(7,v211)=0.75d0*rt2
!  44c
tg(8,v400)=rt35/8d0; tg(8,v040)=rt35/8d0; tg(8,v220)=-0.75*rt3
!  44s
tg(9,v310)=rt5/2d0; tg(9,v130)=-rt5/2d0

!  Deal with shell types, transforming from spherical to cartesian
!  basis if necessary
k=0
do i=1,nshell
  kloc(i)=k+1 ! First basis function for shell i
  select case(shell_type(i))
  case(0) ! s shell
    kmin(i)=1
    kmax(i)=1
    ktype(i)=1
  case(1) ! p shell
    kmin(i)=2
    kmax(i)=4
    ktype(i)=2
  case(2) ! d shell
    kmin(i)=5
    kmax(i)=10
    ktype(i)=3
    if (shell_nfuncs(i) .eq. 5) then ! Spherical d shell
      temp(1:num,1:k)=densty(1:num,1:k)
      temp(1:num,k+1:k+6)=matmul(densty(1:num,k+1:k+5),td)
      temp(1:num,k+7:num+1)=densty(1:num,k+6:num)
      num=num+1
      densty(1:k,1:num)=temp(1:k,1:num)
      densty(k+1:k+6,1:num)=matmul(transpose(td),temp(k+1:k+5,1:num))
      densty(k+7:num,1:num)=temp(k+6:num-1,1:num)
    endif
  case(3) ! f shell
    kmin(i)=11
    kmax(i)=20
    ktype(i)=4
    if (shell_nfuncs(i) .eq. 7) then ! Spherical f shell
      temp(1:num,1:k)=densty(1:num,1:k)
      temp(1:num,k+1:k+10)=matmul(densty(1:num,k+1:k+7),tf)
      if (i<nshell) temp(1:num,k+11:num+3)=densty(1:num,k+8:num)
      num=num+3
      densty(1:k,1:num)=temp(1:k,1:num)
      densty(k+1:k+10,1:num)=matmul(transpose(tf),temp(k+1:k+7,1:num))
      if (i<nshell) densty(k+11:num,1:num)=temp(k+8:num-3,1:num)
    endif
  case(4) ! g shell
    kmin(i)=21
    kmax(i)=35
    ktype(i)=5
    ! print "(a,i0,a,i0)", "num = ", num, "  k = ", k
    if (shell_type(i) .lt. 0) then ! Spherical g shell
      temp(1:num,1:k)=densty(1:num,1:k)
      temp(1:num,k+1:k+15)=matmul(densty(1:num,k+1:k+9),tg)
      if (i<nshell) temp(1:num,k+16:num+6)=densty(1:num,k+10:num)
      num=num+6
      densty(1:k,1:num)=temp(1:k,1:num)
      densty(k+1:k+15,1:num)=matmul(transpose(tg),temp(k+1:k+9,1:num))
      if (i<nshell) densty(k+16:num,1:num)=temp(k+10:num-6,1:num)
    endif
  case default
    write (buffer,"(a,i0)") "Unrecognized or unimplemented shell type ", i
    call die(trim(buffer),.false.)
  end select
  k=k+kmax(i)-kmin(i)+1
end do

allocate(dtri(nx))
k=0
do i=1,num
do j=1,i
  k=k+1
  dtri(k)=densty(i,j)
end do
end do
!~ deallocate(densty)
first=.false.
call dma_main(dtri,kp)

q_out = q

deallocate(zan)
deallocate(c)
deallocate(user_given_limit)
deallocate(kng)
deallocate(katom)
deallocate(shell_nfuncs)
deallocate(shell_type)
deallocate(ex)
deallocate(primCoefs)
deallocate(input_density)
deallocate(kstart, ktype, kloc, kmin, kmax)
deallocate(cs, cp)
deallocate(iax)
deallocate(densty, temp)
deallocate(dtri)

END SUBROUTINE gdma_driver_routine

END MODULE gdma_driver
