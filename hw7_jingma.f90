!--------author:jing ma,jingma@bu.edu,physics department,BU----------!
!--------------------time:23:23,12/06/2016---------------------------!
 module systempara
 implicit none
 integer :: l=20,nn,nr=10,m=200,n=1e4
 real(8) :: betal=0.d0,betar=8.d0,delta,t
 integer, allocatable :: interaction(:,:)
 integer, allocatable :: spin(:)
 real(8), allocatable :: ene(:,:),mag(:,:),ene0(:,:)
 character(len=2) :: namel
 character(len=2) :: namenr
 end module systempara
!--------------------------------------------------------------------!
 program annealing
 use systempara; implicit none
 integer :: i,k,j

 nn=l**3
 write(namel,'(i2)') l
 write(namenr,'(i2)') nr
 delta=(betar-betal)/real(n*m,8)
 call readin
 allocate(spin(nn))
 allocate(ene(m,nr))
 ene=0.0
 allocate(mag(m,nr))
 mag=0.0
 allocate(ene0(m,nr))
 ene0=0.0
 call initran(1)
 do i=1,nr
    spin=1
    do k=1,m
       do j=1,n
          t=1/(betal+((k-1)*n+j)*delta)
          call mcstep
          call measure(k,i)
       enddo
       if (ene0(k,i)>ene(k,i)) then
          ene0(k,i)=ene(k,i)
       endif
       if (k<m) then
          ene0(k+1,i)=ene0(k,i)
       endif
    enddo
    print*,'run',i,'is done!'
 enddo

 do i=1,nr
    do k=1,m
       ene(k,i)=ene(k,i)/real(n,8)
       mag(k,i)=mag(k,i)/real(n,8)
       ene0(k,i)=ene0(k,i)*real(nn,8)/real(n,8)
    enddo
 enddo

 call writedata

 deallocate(interaction)
 deallocate(spin)
 deallocate(ene)
 deallocate(mag)

 end program annealing
!--------------------------------------------------------------------!
 subroutine readin
 use systempara; implicit none
 integer :: i

 allocate(interaction(nn,7))
 open(1,file='l'//trim(adjustl(namel))//'.in',status='old')
 do i=1,nn
    read(1,*) interaction(i,:)
 enddo
 close(1)
 end subroutine readin
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
 use systempara; implicit none
 integer :: i,s,x,y,z,s1,s2,s3,s4,s5,s6
 real(8) :: myran,pflip

 do i=1,nn
    s=int(myran()*nn)+1
    x=mod(s-1,l); y=mod(s-1,l**2)/l; z=(s-1)/(l**2)
    s1=spin(1+mod(x+1,l)+y*l+z*l**2)
    s2=spin(1+mod(x-1+l,l)+y*l+z*l**2)
    s3=spin(1+x+mod(y+1,l)*l+z*l**2)
    s4=spin(1+x+mod(y-1+l,l)*l+z*l**2)
    s5=spin(1+x+y*l+mod(z+1,l)*l**2)
    s6=spin(1+x+y*l+mod(z-1+l,l)*l**2)
    if (myran()<pflip(s,s1,s2,s3,s4,s5,s6)) spin(s)=-spin(s)
 enddo

 end subroutine mcstep
!--------------------------------------------------------------------!
 function pflip(s,s1,s2,s3,s4,s5,s6)
 use systempara; implicit none
 integer :: s,s1,s2,s3,s4,s5,s6
 real(8) :: pflip,temp

 temp=interaction(s,2)*s1+interaction(s,3)*s2+interaction(s,4)*s3
 temp=temp+interaction(s,5)*s4+interaction(s,6)*s5+interaction(s,7)*s6
 pflip=exp(2.d0*spin(s)*temp/t)

 end function pflip
!--------------------------------------------------------------------!
 subroutine measure(k,i)
 use systempara; implicit none
 integer :: k,i,s,x,y,z,s1,s3,s5
 real(8) :: tempe,tempm

 tempe=0.d0
 do s=1,nn
    x=mod(s-1,l); y=mod(s-1,l**2)/l; z=(s-1)/(l**2)
    s1=spin(1+mod(x+1,l)+y*l+z*(l**2))
    s3=spin(1+x+mod(y+1,l)*l+z*(l**2))
    s5=spin(1+x+y*l+mod(z+1,l)*(l**2))
    tempe=tempe+spin(s)*(interaction(s,2)*s1+interaction(s,4)*s3+interaction(s,6)*s5)
 enddo
 ene(k,i)=ene(k,i)+tempe/real(nn,8)
 tempm=abs(sum(spin)) 
 mag(k,i)=mag(k,i)+tempm/real(nn,8)

 end subroutine measure
!--------------------------------------------------------------------!
 subroutine writedata
 use systempara; implicit none
 integer :: k

 open(2,file='e_'//trim(adjustl(namel))//'.dat')
 open(3,file='m_'//trim(adjustl(namel))//'.dat')
 open(4,file='e0_'//trim(adjustl(namel))//'.dat')
 do k=1,m
    call writeline(k)
 enddo
 close(2)
 close(3)
 close(4)

 end subroutine writedata
!--------------------------------------------------------------------!
 subroutine writeline(k)
 use systempara; implicit none
 integer :: k,i
 real(8) :: beta,ave1,ave2,avm1,avm2,ave01,ave02
 
 beta=betal+real(k,8)*n*delta

 ave1=0.d0; ave2=0.d0
 do i=1,nr
    ave1=ave1+ene(k,i)
    ave2=ave2+ene(k,i)**2
 enddo
 ave1=ave1/real(nr,8)
 ave2=sqrt((ave2/real(nr,8)-ave1**2)/real(nr,8))
 write(2,'(3f11.6,'//trim(adjustl(namenr))//'f11.6)')beta,ave1,ave2,ene(k,:)

 avm1=0.d0; avm2=0.d0
 do i=1,nr
    avm1=avm1+mag(k,i)
    avm2=avm2+mag(k,i)**2
 enddo
 avm1=avm1/real(nr,8)
 avm2=sqrt((avm2/real(nr,8)-avm1**2)/real(nr,8))
 write(3,'(3f11.6,'//trim(adjustl(namenr))//'f11.6)')beta,avm1,avm2,mag(k,:)

 ave01=0.d0; ave02=0.d0
 do i=1,nr
    ave01=ave01+ene0(k,i)
    ave02=ave02+ene0(k,i)**2
 enddo
 ave01=ave01/real(nr,8)
 ave02=sqrt((ave02/real(nr,8)-ave01**2)/real(nr,8))
 write(4,'(3f11.2,'//trim(adjustl(namenr))//'f11.2)')beta,ave01,ave02,ene0(k,:)

 end subroutine writeline
!--------------------------------------------------------------------!
















