!--------author:jing ma,jingma@bu.edu,physics department,BU----------!
!--------------------time:16:45,09/22/2016---------------------------!
implicit none

integer(8)             :: nw,n,r0,b,i,j,x,dest
integer(8),allocatable :: matrc(:,:)
integer(8),allocatable :: num(:)
real(8)                :: p
real(8),allocatable    :: matrd(:,:),arr(:)
!----------------input and initialization below----------------------!
open(1,file='read.in',status='old')
read(1,*)nw,n,r0
close(1)

allocate(num(n))

b=64
allocate(matrc(-n:n,0:b-1))
allocate(matrd(-n:n,0:b-1))
allocate(arr(0:b-1))
matrc=0
matrd=0
arr=0
!------------------------calculation below---------------------------!
do i=1,nw
    call crea(num,r0,n)
    do j=0,b-1
        x=dest(num,j,n)
        matrc(x,j)=matrc(x,j)+1
    enddo
enddo

do i=-n,n
    do j=0,b-1
        matrd(i,j)=real(matrc(i,j),8)/real(nw,8)-p(i,n)
    enddo
enddo

do i=-n,n
    do j=0,b-1
        arr(j)=arr(j)+matrd(i,j)**2
    enddo
enddo
!--------------------------output below------------------------------!
open(2,file='d.dat',status='replace')
do j=0,b-1
    write(2,'(i6,f30.18)') j,sqrt(real(nw,8)*arr(j))
enddo
close(2)

do j=0,b-1
    call writ(matrc,j,n,b,nw)
enddo

end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                  above are the main programme                      !
!               below are subroutines and functions                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!--subroutine below creates an array of random numbers of length n,--!
!-using the seed r0 and updates r0 as the last element of array num--!
subroutine crea(num,r0,n)
implicit none

integer(8)             :: num(n),r0,n,a,c,i

a=2862933555777941757_8
c=1013904243_8

num(1)=a*r0+c
do i=2,n
    num(i)=a*num(i-1)+c
enddo

r0=num(n)

end subroutine crea

!----function below calculates the destination of a random walk------!
!---of length n, determined by a random array and the bit postion b--!
function dest(num,b,n)
implicit none

integer(8)             :: dest,num(n),b,n,i

dest=0

do i=1,n
    if (btest(num(i),b)) then
        dest=dest+1
    else
        dest=dest-1
    endif
enddo

end function dest

!----function below calculates the probability distribution of x-----!
!------------------given a random walk of length n-------------------!
function p(x,n)
implicit none

real(8)                :: p,lgfac
integer(8)             :: x,n

if (mod(x,2)==mod(n,2)) then
    p=exp(lgfac(n)-n*log(2._8)-lgfac((n+x)/2)-lgfac((n-x)/2))
else
    p=0._8
endif

end function p

!---------function below calculates log of factorial(n)--------------!
function lgfac(n)
implicit none

real(8)                   :: lgfac
integer(8)                :: n,i

lgfac=0

do i=1,n
    lgfac=lgfac+log(real(i,8))
enddo

end function lgfac

!---------subroutine below writes the full distributions for b=j-----!
!------------------------to files named pnn.dat----------------------!
subroutine writ(matrc,j,n,b,nw)
implicit none

integer(8)                :: j,n,b,nw,i
integer(8)                :: matrc(-n:n,0:b-1)
character(len=7)          :: nam
character                 :: nam1,nam2
real(8)                   :: p

nam1=achar(48+j/10)
nam2=achar(48+j-j/10*10)
nam='p'//nam1//nam2//'.dat'

open(j+3,file=nam,status='replace')
do i=-n,n
    if (matrc(i,j)/=0) then
        write(j+3,'(i6,2f12.6)') i,real(matrc(i,j),8)/real(nw,8),p(i,n)
    endif
enddo
close(j+3)
end subroutine writ
