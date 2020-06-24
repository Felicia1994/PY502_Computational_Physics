!--------author:jing ma,jingma@bu.edu,physics department,BU----------!
!--------------------time:14:43,11/05/2016---------------------------!
 module potentialparameters
 implicit none
 real(8),parameter :: v0=0.1d0,lx=5.d0,ly=10.d0,x0=2.d0,y0=2.d0,y1=5.d0
 real(8), parameter :: pi=3.14159265358979328d0,beta=13.11d0
 integer :: nbasisx,nbasisy  

 end module potentialparameters
!-----------------------------------------------!
 program variational
 use potentialparameters; implicit none

 integer, parameter :: npsix=100,npsiy=200
 integer :: pn,i,j
 real(8) :: energ0,velement,psi(0:npsix,0:npsiy)
 real(8), allocatable :: eig(:),ham(:,:)
 character(len=3) :: nam

 print*,'Number of basis states'
 read(*,*) nbasisx
 read(*,*) nbasisy
 write(nam,'(i3)') nbasisx
 allocate(eig(nbasisx*nbasisy))
 allocate(ham(nbasisx*nbasisy,nbasisx*nbasisy))

 ham=0.d0
 do j=1,nbasisx*nbasisy
   ham(j,j)=energ0(j)+velement(j,j)
   do i=1,j-1
      ham(i,j)=velement(i,j)
      ham(j,i)=ham(i,j)
   end do
 end do
 call diasym(ham,eig,nbasisx*nbasisy)
!--------------------------------------!
 open(1,file=trim(adjustl(nam))//'_eig.dat',status='replace')
 do i=1,nbasisx*nbasisy
    write(1,1)i,eig(i)
 end do
 1 format(i4,' ',f12.7)
 close (1)
!--------------------------------------!
 call fullwavefunction(ham(:,1),psi,npsix,npsiy)
 open(2,file=trim(adjustl(nam))//'_psi1.dat',status='replace')
 do i=0,npsix
   do j=0,npsiy
     write(2,2) lx*real(i,8)/real(npsix,8),ly*real(j,8)/real(npsiy,8),psi(i,j)
   enddo
 enddo
 2 format(f8.5,' ',f8.5,' 'f12.7)
 close (2)

 call fullwavefunction(ham(:,2),psi,npsix,npsiy)
 open(3,file=trim(adjustl(nam))//'_psi2.dat',status='replace')
 do i=0,npsix
   do j=0,npsiy
     write(3,3) lx*real(i,8)/real(npsix,8),ly*real(j,8)/real(npsiy,8),psi(i,j)
   enddo
 enddo
 3 format(f8.5,' ',f8.5,' 'f12.7)
 close (3)

 call fullwavefunction(ham(:,3),psi,npsix,npsiy)
 open(4,file=trim(adjustl(nam))//'_psi3.dat',status='replace')
 do i=0,npsix
   do j=0,npsiy
     write(4,4) lx*real(i,8)/real(npsix,8),ly*real(j,8)/real(npsiy,8),psi(i,j)
   enddo
 enddo
 4 format(f8.5,' ',f8.5,' 'f12.7)
 close (4)

 call fullwavefunction(ham(:,4),psi,npsix,npsiy)
 open(5,file=trim(adjustl(nam))//'_psi4.dat',status='replace')
 do i=0,npsix
   do j=0,npsiy
     write(5,5) lx*real(i,8)/real(npsix,8),ly*real(j,8)/real(npsiy,8),psi(i,j)
   enddo
 enddo
 5 format(f8.5,' ',f8.5,' 'f12.7)
 close (5)

 end program variational
!---------------------------------------------------------------------!
 real(8) function energ0(k)
 use potentialparameters; implicit none
 integer :: k,kx,ky

 kx=mod(k-1,nbasisx)+1
 ky=(k-1)/nbasisx+1
 energ0=pi**2*((real(kx,8)/lx)**2+(real(ky,8)/ly)**2)/(2*beta)

 end function energ0
!---------------------------------------!
 real(8) function velement(k1,k2)
 use potentialparameters; implicit none
 integer :: k1,k2,k1x,k1y,k2x,k2y
 real(8) :: vx,vy

 k1x=mod(k1-1,nbasisx)+1
 k1y=(k1-1)/nbasisx+1
 k2x=mod(k2-1,nbasisx)+1
 k2y=(k2-1)/nbasisx+1
 velement=-v0*vx(k1x,k2x)*vy(k1y,k2y)
! velement=-v0*vy(k1y,k2y)

 end function velement
!---------------------------------------!
 real(8) function vx(k1x,k2x)
 use potentialparameters; implicit none
 integer :: k1x,k2x
 real(8) :: a,b

 a=real(k1x,8)*pi/lx
 b=real(k2x,8)*pi/lx
 if (k1x==k2x) then
   vx=(lx-2*x0)/2.d0-(sin(2.d0*a*(lx-x0))-sin(2.d0*a*x0))/4.d0/a
 else
   vx=(sin((a-b)*(lx-x0))-sin((a-b)*x0))/2.d0/(a-b)-(sin((a+b)*(lx-x0))-sin((a+b)*x0))/2.d0/(a+b)
 endif
   vx=vx*2.d0/lx
 
 end function vx
!---------------------------------------!
 real(8) function vy(k1y,k2y)
 use potentialparameters; implicit none
 integer :: k1y,k2y
 real(8) :: a,b

 a=real(k1y,8)*pi/ly
 b=real(k2y,8)*pi/ly
 if (k1y==k2y) then
   vy=(ly-2*y0-y1)/2.d0-(sin(2.d0*a*(ly-y0))-sin(a*(ly+y1)))/4.d0/a-(sin(a*(ly-y1))-sin(2.d0*a*y0))/4.d0/a
 else
   vy=(sin((a-b)*(ly-y0))-sin((a-b)*(ly+y1)/2.d0))/2.d0/(a-b)-(sin((a+b)*(ly-y0))-sin((a+b)*(ly+y1)/2.d0))/2.d0/(a+b)
   vy=vy+(sin((a-b)*(ly-y1)/2.d0)-sin((a-b)*y0))/2.d0/(a-b)-(sin((a+b)*(ly-y1)/2.d0)-sin((a+b)*y0))/2.d0/(a+b)
 endif
   vy=vy*2.d0/ly
 
 end function vy

!---------------------------------------------------------!
!Calls the LAPACK diagonalization subroutine DSYEV        !
!input:  a(n,n) = real symmetric matrix to be diagonalized!
!            n  = size of a                               !
!output: a(n,n) = orthonormal eigenvectors of a           !
!        eig(n) = eigenvalues of a in ascending order     !
!---------------------------------------------------------!
 subroutine diasym(a,eig,n)
 implicit none
 integer :: n,l,inf
 real(8) :: a(n,n),eig(n),work(n*(3+n/2))

 l=n*(3+n/2)
 call dsyev('V','U',n,a,n,eig,work,l,inf)

 end subroutine diasym
!-----------------------------------------------------------!
 subroutine fullwavefunction(vec,psi,npsix,npsiy)
 use potentialparameters; implicit none
 integer :: npsix,npsiy,i,j,k,kx,ky
 real(8) :: vec(nbasisx*nbasisy),psi(0:npsix,0:npsiy),x,y,wavefunc

 psi=0.d0
 do i=0,npsix
   do j=0,npsiy
     x=lx*real(i,8)/real(npsix,8)
     y=ly*real(j,8)/real(npsiy,8)
     do k=1,nbasisx*nbasisy
       kx=mod(k-1,nbasisx)+1
       ky=(k-1)/nbasisx+1
       psi(i,j)=psi(i,j)+vec(k)*wavefunc(kx,ky,x,y)
     enddo
   enddo
 enddo

 end subroutine fullwavefunction
!---------------------------------------------------------------! 
 real(8) function wavefunc(kx,ky,x,y)
 use potentialparameters;implicit none
 integer :: kx,ky
 real(8) :: x,y

 wavefunc=sin(real(kx,8)*pi*x/lx)*sin(real(ky,8)*pi*y/ly)*2.d0/sqrt(lx*ly)

 end function wavefunc
!---------------------------------------------------------------!
