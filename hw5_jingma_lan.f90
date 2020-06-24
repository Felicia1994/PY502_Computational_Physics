!--------author:jing ma,jingma@bu.edu,physics department,BU----------!
!--------------------time:14:43,11/05/2016---------------------------!
 module systemparameters
 implicit none
 real(8),parameter :: v0=0.1d0,lx=5.d0,ly=10.d0,x0=2.d0,y0=2.d0,y1=5.d0
 real(8), parameter :: beta=13.11d0
 integer :: d=2
 integer :: nx,ny
 real(8), allocatable :: vpot(:)
 real(8) :: delta

 end module systemparameters
!-----------------------------------------!
 module lanczosvectors
 implicit none
 real(8), save, allocatable :: f0(:)
 real(8), save, allocatable :: f1(:)
 real(8), save, allocatable :: f2(:)
 real(8), save, allocatable :: aa(:)
 real(8), save, allocatable :: nn(:)

 end module lanczosvectors
!-------------------------------------------!
 program lanczos
 use systemparameters; use lanczosvectors; implicit none

 integer :: niter,i,j
 real(8) :: escale
 real(8), allocatable :: psi(:)
 real(8), allocatable :: eig(:)
 real(8), allocatable :: vec(:,:)
 character(len=8) :: str1,str2

 write(*,*)'Delta'
 read(*,*) delta
 nx=lx/delta
 ny=ly/delta
 escale=2.d0*delta**2
 write(str1,'(f8.6)') delta
 write(*,*)'Lanczos iterations'
 read(*,*) niter
 write(str2,'(i8)') niter

 allocate(f0(nx*ny))
 allocate(f1(nx*ny))
 allocate(f2(nx*ny))
 allocate(aa(0:niter))
 allocate(nn(0:niter))

 allocate(psi(nx*ny))
 allocate(eig(0:niter-1))
 allocate(vec(0:niter-1,0:niter-1))

 call potentialenergy(escale)
 call initstate(psi,nx*ny) 
 call lanczos1(nx*ny,niter,psi,eig,vec)
!------------------------------------------!
 open(1,file=trim(adjustl(str1))//'_'//trim(adjustl(str2))//'_'//'e.dat',status='replace')
 do i=0,niter-1
    write(1,1) i+1,eig(i)/escale/beta
 enddo
 1 format(i4,' ',f12.7)
 close(1)
!------------------------------------------!
 call lanczos2(nx*ny,niter,psi,vec(:,0))
 open(2,file=trim(adjustl(str1))//'_'//trim(adjustl(str2))//'_'//'p1.dat',status='replace')
 do i=1,nx
   do j=1,ny
     write(2,2) delta*real(i,8),delta*real(j,8),psi(i+(j-1)*nx)
   enddo
 enddo
 2 format(f8.5,' ',f8.5,' 'f12.7)
 close (2)

 call lanczos2(nx*ny,niter,psi,vec(:,1))
 open(3,file=trim(adjustl(str1))//'_'//trim(adjustl(str2))//'_'//'p2.dat',status='replace')
 do i=1,nx
   do j=1,ny
     write(3,3) delta*real(i,8),delta*real(j,8),psi(i+(j-1)*nx)
   enddo
 enddo
 3 format(f8.5,' ',f8.5,' 'f12.7)
 close (3)

 call lanczos2(nx*ny,niter,psi,vec(:,2))
 open(4,file=trim(adjustl(str1))//'_'//trim(adjustl(str2))//'_'//'p3.dat',status='replace')
 do i=1,nx
   do j=1,ny
     write(4,4) delta*real(i,8),delta*real(j,8),psi(i+(j-1)*nx)
   enddo
 enddo
 4 format(f8.5,' ',f8.5,' 'f12.7)
 close (4)

 call lanczos2(nx*ny,niter,psi,vec(:,3))
 open(5,file=trim(adjustl(str1))//'_'//trim(adjustl(str2))//'_'//'p4.dat',status='replace')
 do i=1,nx
   do j=1,ny
     write(5,5) delta*real(i,8),delta*real(j,8),psi(i+(j-1)*nx)
   enddo
 enddo
 5 format(f8.5,' ',f8.5,' 'f12.7)
 close (5)
!------------------------------------------!
 deallocate(f0)
 deallocate(f1)
 deallocate(f2)
 deallocate(aa)
 deallocate(nn)
 deallocate(psi)
 deallocate(eig)
 deallocate(vec)
 deallocate(vpot)
 end program lanczos
!-----------------------------------------------------------------!
 subroutine potentialenergy(escale)
 use systemparameters; implicit none
 real(8) :: escale
 integer :: i,j
 logical :: bool1,bool2,bool3

 allocate(vpot(nx*ny))
 vpot(:)=0.d0
 do i=1,nx
   do j=1,ny
     bool1=(delta*i>=x0).and.(delta*i<lx-x0)
     bool2=(delta*j>=y0).and.(delta*j<(ly-y1)/2.d0)
     bool3=(delta*j>=(ly+y1)/2.d0).and.(delta*j<(ly-y0))
     if (bool1.and.(bool2.or.bool3)) then
       vpot(i+(j-1)*nx)=-beta*v0
     endif
   enddo
 enddo
 vpot(:)=vpot(:)*escale
 vpot(:)=vpot(:)+2.d0*d

 end subroutine potentialenergy
!------------------------------!
 subroutine initstate(psi,n)
 implicit none
 integer :: n,i
 real(8) :: psi(n),ran,norm

 call initran(1)
 do i=1,n
    psi(i)=ran()-0.5d0
 end do
 norm=1.d0/sqrt(dot_product(psi,psi))
 do i=1,n
    psi(i)=psi(i)*norm
 end do
 
 end subroutine initstate
!-----------------------------------------------------!
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
!--------------------------------------------!
 real(8) function ran()
 implicit none
 real(8)    :: dmu64
 integer(8) :: ran64,mul64,add64
 common/bran64/dmu64,ran64,mul64,add64

 ran64=ran64*mul64+add64
 ran=0.5d0+dmu64*dble(ran64)

 end function ran
!----------------------------------------------!
 subroutine lanczos1(n,niter,p0,eig,vec)
!--------------------------------------------------------------!
! This subroutine is constructing the normalized basis states  !
! directly, with the coefficients aa(m) and nn(m) being the    !
! matrix elements of the tri-diagonal hamiltonian matrix.      !
!--------------------------------------------------------------!
 use systemparameters; use lanczosvectors; implicit none
 integer :: n,niter,m
 real(8) :: p0(n),eig(0:niter-1),vec(0:niter-1,0:niter-1),t

 f0=p0
 nn(0)=1.d0
 call hamoperation(nx*ny,f0,f1)
 aa(0)=dot_product(f0,f1)
 f1=f1-aa(0)*f0
 nn(1)=sqrt(dot_product(f1,f1))
 f1=f1/nn(1)
 do m=2,niter
    write(*,*)'Initial Lanczos iteration: ',m
    call hamoperation(nx*ny,f1,f2)
    aa(m-1)=dot_product(f1,f2)
    f2=f2-aa(m-1)*f1-nn(m-1)*f0
    nn(m)=sqrt(dot_product(f2,f2))
    f2=f2/nn(m)
    f0=f1
    f1=f2
 enddo
 call diatri(niter,eig,vec)

 end subroutine lanczos1
!----------------------------!
 subroutine diatri(n,eig,vec)
 use lanczosvectors; implicit none
 integer :: n,inf
 real(8) :: eig(n),vec(n,n),d(n),e(n),work(max(1,2*n-2))

 d=aa(0:n-1)
 e=nn(1:n)
 call dstev('V',n,d,e,vec,n,work,inf)
 eig=d

 end subroutine diatri
!------------------------------------!
 subroutine lanczos2(n,niter,psi,vec)
!-------------------------------------------------------------------!
! This subroutine re-constructs the normalized Lanczos basis states
! using the coefficients aa(m) and nn(m) previously computed in the 
! subroutine lanczos1(). The vector vec() is an eigenstate in the 
! Lanczos basis and it's tranformed to the real-space basis state 
! stored as the vector psi().
!-------------------------------------------------------------------!
 use systemparameters; use lanczosvectors; implicit none
 integer :: n,niter,m
 real(8) :: psi(n),vec(0:niter-1)
 
 f0=psi
 psi=psi*vec(0)
 call hamoperation(nx*ny,f0,f1)
 f1=(f1-aa(0)*f0)/nn(1)
 psi=psi+vec(1)*f1
 do m=2,niter-1
    write(*,*)'Second Lanczos iteration: ',m
    call hamoperation(nx*ny,f1,f2)
    f2=(f2-aa(m-1)*f1-nn(m-1)*f0)/nn(m)
    psi=psi+vec(m)*f2
    f0=f1
    f1=f2
 enddo
 
 end subroutine lanczos2
!-----------------------!

!--------------------------------!
 subroutine hamoperation(n,f1,f2)
!-------------------------------------------------------!
! Acting with H on a state vector f1(), leading to f2().!
!-------------------------------------------------------!
 use systemparameters; implicit none
 integer :: n,k,kx,ky,kz
 real(8) :: f1(n),f2(n)

 f2=vpot*f1
! if (d==1) then
!    do j=1,n
!       if (j.ne.1) f2(j)=f2(j)-f1(j-1)
!       if (j.ne.n) f2(j)=f2(j)-f1(j+1)
!    end do
! else if (d==2) then
    do k=1,nx*ny
       kx=mod(k-1,nx)+1
       ky=(k-1)/nx+1
       if (kx.ne.1) f2(k)=f2(k)-f1(k-1)
       if (kx.ne.nx) f2(k)=f2(k)-f1(k+1)
       if (ky.ne.1) f2(k)=f2(k)-f1(k-nx)
       if (ky.ne.ny) f2(k)=f2(k)-f1(k+nx)
    end do
! else if (d==3) then
!    ll=l**2
!    do j=1,n
!       x=1+mod(j-1,l)
!       y=1+mod(j-1,ll)/l
!       z=1+(j-1)/ll    
!       if (x.ne.1) f2(j)=f2(j)-f1(j-1)
!       if (x.ne.l) f2(j)=f2(j)-f1(j+1)
!       if (y.ne.1) f2(j)=f2(j)-f1(j-l)
!       if (y.ne.l) f2(j)=f2(j)-f1(j+l)
!       if (z.ne.1) f2(j)=f2(j)-f1(j-ll)
!       if (z.ne.l) f2(j)=f2(j)-f1(j+ll)
!    end do
! end if

 end subroutine hamoperation
!---------------------------!
