!--------author:jing ma,jingma@bu.edu,physics department,BU----------!
!--------------------time:22:01,10/08/2016---------------------------!
 module systemparam
 implicit none

 real(8), parameter :: pi=3.141592653589793d0
 real(8), parameter :: ts = 86164d0
 real(8), parameter :: rgeo = 42168d3
 real(8), parameter :: g = 6.6743d-11
 real(8), parameter :: me = 5.9736d24
 real(8), parameter :: rm = 384400d3

 end module systemparam
!--------------------------------------------------------------------!
 program satellite

 use systemparam
 implicit none
 
 real(8)     :: alpha,delta,a,t,ax,ay,az,vx,vy,vz,x,y,z
 integer(8)  :: nt,tmax,nw,i,nu=0,nd=0,nuu=0,ndd=0
 real(8)     :: phi,r,deltaphi,deltar,theta

 print*,'inclination angle alpha in degrees'; read*,alpha
 alpha = alpha/180*pi
 print*,'time step delta as 1 nt-th of the siderial day'; read*,nt
 delta = ts/nt
 print*,'integration time tmax times of the siderial day';read*,tmax
 tmax = nt*tmax
 print*,'writing data every nw steps';read*,nw
 print*,'initial conditions set up by a';read*,a

 x = rgeo*a**(2./3); y = 0d0; z = 0d0
 vx = 0d0; vy = sqrt(g*me/x); vz = 0d0

 open(1,file='sat.dat',status='replace')
 do i=1,tmax
     t = delta*(i-1)
     if (mod(i,nw)==0) then
         phi = atan(y/x)+pi*(nu+nd)
         r = sqrt(x**2+y**2)/1d3
         deltaphi = phi-2*pi*t/ts
         deltar = r-rgeo/1d3
         theta = atan(z/r)+pi*(nuu+ndd)
         write(1,1) t,phi,r,deltaphi,deltar,theta
         1 format(f12.2,' ',f12.6,' ',f12.6,' ',f12.8,' ',f12.6,' ',f20.8)
     endif
     call acceleration(ax,ay,az,x,y,z,t,alpha)
     vx = vx+delta*ax; vy = vy+delta*ay; vz = vz+delta*az
     if (y>0 .and. x>0 .and. x+delta*vx<0) then
         nu = nu+1
     endif
     if (y<0 .and. x<0 .and. x+delta*vx>0) then
         nd = nd+1
     endif
     if (z>0 .and. x>0 .and. x+delta*vx<0) then
         nuu = nuu+1
     endif
     if (z>0 .and. x<0 .and. x+delta*vx>0) then
         nuu = nuu-1
     endif
     if (z<0 .and. x<0 .and. x+delta*vx>0) then
         ndd = ndd+1
     endif
     if (z<0 .and. x>0 .and. x+delta*vx<0) then
         ndd = ndd-1
     endif
     x=x+delta*vx; y=y+delta*vy; z=z+delta*vz
 enddo
 close(1)

!print*,deltaphi

 end program satellite
!--------------------------------------------------------------------!
 subroutine acceleration(ax,ay,az,x,y,z,t,alpha)
 
 use systemparam
 implicit none

 real(8)     :: ax,ay,az,x,y,z,t,alpha,xm,ym,zm

 ax = -g*me*x/(sqrt(x**2+y**2+z**2)**3)
 ay = -g*me*y/(sqrt(x**2+y**2+z**2)**3)
 az = -g*me*z/(sqrt(x**2+y**2+z**2)**3)
 
 xm = rm*cos(alpha)*cos(2*pi*t/(27.25d0*ts))
 ym = rm*sin(2*pi*t/(27.25d0*ts))
 zm = rm*sin(alpha)*cos(2*pi*t/(27.25d0*ts))

 ax = ax-1.23d-2*g*me*(x-xm)/(sqrt((x-xm)**2+(y-ym)**2+(z-zm)**2)**3)
 ay = ay-1.23d-2*g*me*(y-ym)/(sqrt((x-xm)**2+(y-ym)**2+(z-zm)**2)**3)
 az = az-1.23d-2*g*me*(z-zm)/(sqrt((x-xm)**2+(y-ym)**2+(z-zm)**2)**3)

 end subroutine acceleration








