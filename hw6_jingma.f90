!--------author:jing ma,jingma@bu.edu,physics department,BU----------!
!--------------------time:15:55,11/20/2016---------------------------!
 module systempara
 implicit none
 integer :: l,bins,reps,steps,n
 real(8) :: t
 real(8) :: pflip(-1:1,-4:4)
 integer, allocatable :: spin(:)
 character(len=3) :: naml
 character(len=15) :: namt
 end module systempara
!--------------------------------------------------------------------!
 module para1
 implicit none
 real(8), allocatable :: mag(:)
 end module para1
!--------------------------------------------------------------------!
 module para2
 implicit none
 real(8) :: tim
 end module para2
!--------------------------------------------------------------------!
 program myising2d
 use systempara; implicit none
 integer :: i

 open(1,file='read.in',status='old')
 read(1,*) l,t
 n=l**2
 write(naml,'(i3)') l
 write(namt,'(f15.13)') t
 read(1,*) bins,reps,steps
 close(1)
 do i=-4,4
    pflip(-1,i)=exp(+2.*i/t)
    pflip(+1,i)=exp(-2.*i/t)
 enddo
 allocate (spin(0:n-1))

 if (steps/=0) then
    call part1
 else
    call part2
 endif

 deallocate(spin)

 end program myising2d
!--------------------------------------------------------------------!
 subroutine part1
 use systempara; use para1; implicit none
 integer :: i,j,k
 print*, 'this is part1!'

 allocate (mag(1:steps))
 call initran(1)

 open(2,file='bindata_'//trim(adjustl(naml))//'_2.dat')
 do i=1,bins
    mag=0.d0
    do j=i,reps
       spin=1
       do k=1,steps
          call mcstep
          call measure(k)
       enddo
    enddo
    call writebindata
 enddo
 close(2)

 deallocate(mag)
 call average1

 end subroutine part1
!--------------------------------------------------------------------!
 subroutine initran(w)
 implicit none
 integer(8) :: irmax
 integer(4) :: w,nb,b
 real(8)    :: dmu64
 integer(8) :: ran64,mul64,add64
 common/bran64/dmu64,ran64,mul64,add64
      
 irmax=2_8**31
 irmax=2*(irmax**2-1)+1
 mul64=2862933555777941757_8
 add64=1013904243
 dmu64=0.5d0/dble(irmax)

 open(10,file='seed.in',status='old')
 read(10,*)ran64
 close(10)
 if (w.ne.0) then
    open(10,file='seed.in',status='unknown')
    write(10,*)abs((ran64*mul64)/5+5265361)
    close(10)
 endif

 end subroutine initran
!--------------------------------------------------------------------!
 real(8) function myran()
 implicit none
 real(8)    :: dmu64
 integer(8) :: ran64,mul64,add64
 common/bran64/dmu64,ran64,mul64,add64

 ran64=ran64*mul64+add64
 myran=0.5d0+dmu64*dble(ran64)

 end function myran
!--------------------------------------------------------------------!
 subroutine mcstep
 use systempara; use para1; implicit none
 integer :: i,s,x,y,s1,s2,s3,s4
 real(8), external :: myran

 do i=1,n
    s=int(myran()*n)
    x=mod(s,l); y=s/l
    s1=spin(mod(x+1,l)+y*l)
    s2=spin(x+mod(y+1,l)*l)
    s3=spin(mod(x-1+l,l)+y*l)
    s4=spin(x+mod(y-1+l,l)*l)
    if (myran()<pflip(spin(s),s1+s2+s3+s4)) spin(s)=-spin(s)
 enddo

 end subroutine mcstep
!--------------------------------------------------------------------!
 subroutine measure(k)
 use systempara; use para1; implicit none
 integer :: k,m

 m=sum(spin)
 mag(k)=mag(k)+real(m,8)

 end subroutine measure
!--------------------------------------------------------------------!
 subroutine writebindata
 use systempara; use para1; implicit none
 integer :: i

 mag=mag/(real(reps,8)*real(n,8))
 do i=1,steps
    write(2,1) i,mag(i)
 enddo
 1 format(i8,f18.12)

 end subroutine writebindata
!--------------------------------------------------------------------!
 subroutine average1
 use systempara; implicit none
 integer :: nr,i,j
 real(8), allocatable :: av(:,:)
 real(8) :: mg

 allocate(av(2,steps))

 nr=0
 av=0.d0
 open(1,file='bindata_'//trim(adjustl(naml))//'_2.dat',status='old')
 do
    do i=1,steps
       read(1,*,end=2)j,mg
       av(1,i)=av(1,i)+mg
       av(2,i)=av(2,i)+mg**2
    enddo
    nr=nr+1
 enddo
 2 close(1)
 write(*,*)'Number of bins read: ',nr
 av=av/real(nr,8)
 av(2,:)=sqrt((av(2,:)-av(1,:)**2)/real(nr,8))

 open(3,file='bindatacal_'//trim(adjustl(naml))//'_2.dat')
 do i=1,steps
    write(3,'(i8,2f14.8)') i,av(:,i)
 enddo
 close(3)

 deallocate(av)

 end subroutine average1
!--------------------------------------------------------------------!
 subroutine part2
 use systempara; use para2; implicit none
 integer :: i,j,flag,pos1,pos2,timetime
 real :: time0,time1
 print*, 'this is part2!'

 call cpu_time(time0)
 call initran(1)

 open(3,file='res_'//trim(adjustl(naml))//'_'//trim(adjustl(namt))//'.dat')
 do i=1,bins
    tim=0.d0
    do j=i,reps
       spin=1
       flag=0
       pos1=0
       pos2=0
       do
          pos1=pos1+1
          call mcstep2(flag,pos2)
          if (flag==1) exit
       enddo
       call measure2(pos1,pos2)
    enddo
    call writedata
 enddo
 close(3)

 call cpu_time(time1)
 timetime=int(time1-time0)
 call average2(timetime)

 end subroutine part2
!--------------------------------------------------------------------!
 subroutine mcstep2(flag,pos2)
 use systempara; use para2; implicit none
 integer :: flag,pos2,i,s,x,y,s1,s2,s3,s4,m
 real(8), external :: myran

 do i=1,n
    s=int(myran()*n)
    x=mod(s,l); y=s/l
    s1=spin(mod(x+1,l)+y*l)
    s2=spin(x+mod(y+1,l)*l)
    s3=spin(mod(x-1+l,l)+y*l)
    s4=spin(x+mod(y-1+l,l)*l)
    if (myran()<pflip(spin(s),s1+s2+s3+s4)) then
       spin(s)=-spin(s)
       m=sum(spin)
       if (m==0) then
          flag=1
          pos2=i
       endif
    endif
 enddo

 end subroutine mcstep2
!--------------------------------------------------------------------!
 subroutine measure2(pos1,pos2)
 use systempara; use para2; implicit none
 integer :: pos1,pos2

 tim=tim+real(pos1-1,8)+real(pos2,8)/n

 end subroutine measure2
!--------------------------------------------------------------------!
 subroutine writedata
 use systempara; use para2; implicit none

 tim=tim/(real(reps,8))
 write(3,2) tim
 2 format(f18.12)

 end subroutine writedata
!--------------------------------------------------------------------!
 subroutine average2(timetime)
 use systempara; implicit none
 integer :: timetime,nr
 real(8) :: mg,av(2)

 nr=0
 av=0.d0
 open(4,file='res_'//trim(adjustl(naml))//'_'//trim(adjustl(namt))//'.dat',status='old')
 do 
    read(4,*,end=5) mg
    av(1)=av(1)+mg
    av(2)=av(2)+mg**2
    nr=nr+1
 enddo
 5 close(4)
 write(*,*)'Number of bins read: ',nr
 av=av/real(nr,8)
 av(2)=sqrt((av(2)-av(1)**2)/real(nr,8))

 open(6,file='rescal_'//trim(adjustl(naml))//'_'//trim(adjustl(namt))//'.dat')
 write(6,'(2f16.4,i10)') av(:),timetime
 close(6)

 end subroutine average2
!--------------------------------------------------------------------!
