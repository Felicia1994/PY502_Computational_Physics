!--------author:jing ma,jingma@bu.edu,physics department,BU----------!
!--------------------time:17:22,09/27/2016---------------------------!
implicit none

integer(8)             :: seed
integer                :: npt,nbi,i,j
real(8)                :: xyz(3),totx,totz,ax,az,sigx,sigz
real                   :: r1,r2,rho1,rho2,rho,ix,iz

!----------------input and initialization below----------------------!
open(1,file='seed.in',status='old')
read(1,*) seed
close(1)
open(2,file='read.in',status='old')
read(2,*) r1,r2,rho1,rho2,npt,nbi
close(2)
ax=0;az=0
sigx=0;sigz=0
!-----------------------calculations below---------------------------!
open(3,file='bin.dat')
do i=1,nbi
    totz=0
    totx=0
    do j=1,npt
        call posi(xyz,seed,r1)
        if ((xyz(1)**2+xyz(2)**2+xyz(3)**2)>r1**2) then
            rho=0
        elseif ((xyz(1)**2+xyz(2)**2)<r2**2) then
            rho=rho2
        else
            rho=rho1
        endif
        totz=totz+rho*(xyz(1)**2+xyz(2)**2)
        totx=totx+rho*(xyz(2)**2+xyz(3)**2)
    enddo
    ix=(totx/npt)*8*(r1**3)
    iz=(totz/npt)*8*(r1**3)
    ax=ax+ix;az=az+iz
    sigx=sigx+ix**2;sigz=sigz+iz**2
    write(3,*)i,iz,ix
enddo
close(3)
ax=ax/nbi;az=az/nbi
sigx=sqrt((sigx/nbi-ax**2)/(nbi-1));sigz=sqrt((sigz/nbi-az**2)/(nbi-1))
open(4,file='res.dat')
write(4,*)ax,sigx
write(4,*)az,sigz
close(4)
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                  above are the main programme                      !
!                     below are subroutines                          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine posi(xyz,seed,r1)
implicit none
integer(8)             :: seed,temp(3),a,c
real(8)                :: dmul,xyz(3)
real                   :: r1
a=2862933555777941757_8
c=1013904243_8
temp(1)=a*seed+c
temp(2)=a*temp(1)+c
temp(3)=a*temp(2)+c
seed=temp(3)
dmul=1.d0/dble(2*(2_8**62-1)+1)
xyz(1)=dble(temp(1))*dmul*r1
xyz(2)=dble(temp(2))*dmul*r1
xyz(3)=dble(temp(3))*dmul*r1
end subroutine posi

