!--------author:jing ma,jingma@bu.edu,physics department,BU----------!
!--------------------time:22:01,10/08/2016---------------------------!
 module popularvariables
 implicit none
 integer(8)        :: nr
 real(8)           :: rmax,h,h12,a,pi=3.141592653589793d0
 real(8),parameter :: r0=0.1d0,e=-2.226d0

 end module popularvariables
!--------------------------------------------------------------------!
 program schrodinger
 use popularvariables
 implicit none
 integer(8)          :: i,j
 real(8)             :: dv=1.d0,v0,v1,v2,b0,b1,b2,eps=1.d-6,r,norm
 real(8),allocatable :: u(:)

 write(*,*)'integrating from rmax to r0; give rmax:'
 read(*,*)rmax
 write(*,*)'number of points:'
 read(*,*)nr
 allocate(u(0:nr))
 h=(rmax-r0)/real(nr,8)
 h12=h**2/12.d0

 open(1,file='vr.dat')
 do j=0,50
   a=0.5d0+(3.d0-0.5d0)*real(j,8)/50.d0
   v1=0.d0
   call integrate(v1,u)
   b1=u(0)
   do
     v2=v1+dv
     call integrate(v2,u)
     b2=u(0)
     if (b1*b2 < 0.d0) exit
     v1=v2
     b1=b2
   enddo

   do
     v0=(v1+v2)/2.d0
     call integrate(v0,u)
     b0=u(0)
     if (abs(b0) <= eps) exit
     if (b0*b1 <= 0.d0) then
       v2=v0      
       b2=b0 
     else
       v1=v0
       b1=b0
     endif
   enddo

   norm=0.d0
   do i=1,nr-3,2
     r=r0+real(i,8)*h
     norm=norm+4.d0*(u(i)**2)*(r**2)+2.d0*(u(i+1)**2)*((r+h)**2)
   enddo
   norm=norm+(u(0)**2)*(r0**2)+4.d0*(u(nr-1)**2)*((rmax-h)**2)+(u(nr)**2)*(rmax**2)
   r=sqrt(norm*pi*h/3.d0)
   write(1,1)a,v0,r
   1 format(f12.6,f12.6,f12.6)
 enddo
 close(1)

 end program schrodinger
!--------------------------------------------------------!
 subroutine integrate(v0,u)
 use popularvariables
 implicit none
 real(8)             :: v0,p1,q2,q1,q0,f,r
 real(8)             :: u(0:nr)
 integer(8)          :: i

 u(nr-1)=1.d-4
 u(nr)=u(nr-1)*exp(-sqrt(0.053667)*h)
 call setinitcond(v0,u(nr),u(nr-1),q2,q1,f,p1)
 do i=2,nr
   r=r0+real(nr-i,8)*h
   call numerovstep(r,v0,q2,q1,f,p1)
   u(nr-i)=p1
 end do
 call normalize(u)

 end subroutine integrate
!------------------------------------------------------------!
 subroutine setinitcond(v0,u2,u1,q2,q1,f,p1)
 use popularvariables
 implicit none
 real(8) :: v0,u2,u1,q2,q1,f,p1,potential

 f=-0.053667d0*(potential(rmax,v0)-1)        
 q2=u2*(1.d0-h12*f)
 f=-0.053667d0*(potential(rmax-h,v0)-1)        
 q1=u1*(1.d0-h12*f)
 p1=u1
	
 end subroutine setinitcond
!----------------------------------------------------------------!
 subroutine numerovstep(r,v0,q2,q1,f,p1)
 use popularvariables
 implicit none
 real(8)             :: r,v0,q2,q1,f,p1,q0,potential

 q0=2.d0*q1-q2+h**2*f*p1
 q2=q1; q1=q0
 f=-0.053667d0*(potential(r,v0)-1)        
 p1=q1/(1.d0-h12*f)
         
 end subroutine numerovstep
!------------------------------------------------------------------!
 real(8) function potential(r,v0)
 use popularvariables
 implicit none
 real(8) :: potential,r,v0

 if (r >= r0) then
   potential=-v0*exp(-r/a)/(e*r/a)
 else
   print*,r
   print*,'wrong radius!'
 endif 

 end function potential
!----------------------------------------------------------------!
 subroutine normalize(u)
 use popularvariables
 implicit none

 integer(8)   :: i
 real(8)      :: u(0:nr),norm

 norm=0.d0
 do i=1,nr-3,2
   norm=norm+4.d0*u(i)**2+2.d0*u(i+1)**2
 end do
 norm=norm+u(0)**2+4.d0*u(nr-1)**2+u(nr)**2
 norm=1.d0/sqrt(norm*4*pi*h/3.d0)
 do i=0,nr
   u(i)=u(i)*norm
 end do

 end subroutine normalize
!-------------------------------------------------------!
