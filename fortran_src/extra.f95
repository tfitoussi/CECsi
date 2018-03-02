#include "../temp/preprocessing.f95"
!=====================================================================
!     Here are summarize all the functions and routines needed 
!     for photons and leptons 
!=====================================================================

!..............................................................................
!          random number generator from Numerical Recipes (Fortran90)         .
!..............................................................................
function ran0()
  use Constantes, only : iseed, prec
  implicit none
  real(prec) :: ran0
  integer,parameter :: IA=16807,IM=2147483647,IQ=127773,IR=2836
  real, save :: am
  integer, save :: ix=-1, iy=-1,k
  !$OMP THREADPRIVATE(am,ix,iy)

  if (iseed <= 0 .or. iy < 0) then
     am = nearest(1.e0,-1.e0)/IM
     iy = ior(ieor(888889999,abs(iseed)),1)
     ix = ieor(777755555,abs(iseed))
     iseed = abs(iseed)+1
  end if
  ix = ieor(ix,ishft(ix,13))
  ix = ieor(ix,ishft(ix,-17))
  ix = ieor(ix,ishft(ix,5))
  k = iy/IQ
  iy = IA*(iy-k*IQ)-IR*k
  if (iy<0) iy = iy+IM
  ran0 = am*ior(iand(IM,ieor(ix,iy)),1)
end function ran0

! Compute the random distribution of a cone
!===========================================
subroutine isotrop(u,v,w)
   use Constantes, only : prec, pi
   implicit none
   real(prec), intent(out):: u,v,w
   real(prec):: a1,a2,ran0
   a1=ran0()
   a2=ran0()
   w=a1*2.-1.
   u=sqrt(1-w*w)*cos(2*pi*a2)
   v=sqrt(1.-w*w)*sin(2.*pi*a2)
end subroutine isotrop

!------------------------------------------------------------------------------ 
! Draw a photon from the EBL
!------------------------------------------------------------------------------
subroutine selecPhotEBL(ee,zint)
   use Constantes, only : prec, enerEBL, maxenerEBL, minenerEBL
   implicit none
   external NeInterpol 
   real(prec), intent(in):: zint
   real(prec), intent(out):: ee
   real(prec) :: a1,a2,ran0,maxEBL,Ne,NN
   real(prec), dimension(:), allocatable :: nEBLinter
   logical :: reject
   integer :: i

   ! Compute the "curve" of the EBL density vs energy for our specific redshift
   ! which is most probably not tabulated

   allocate(nEBLinter(size(enerEBL)))
   do i=1, size(enerEBL)
      Call NeInterpol(enerEBL(i),zint,Ne)
      nEBLinter(i) = Ne
   end do
   ! extraction of the min, max values
   maxEBL=maxval(enerEBL*nEBLinter)
   deallocate(nEBLinter)

   ! Rejection on the EBL values
   reject=.TRUE.
   do while (reject)
      a1=ran0()
      a2=ran0()
      ee=log(minEnerEBL)+a1*(log(maxEnerEBL/minEnerEBL))
      NN=a2*(maxEBL)
      Call NeInterpol(exp(ee),zint,Ne)
      if (NN<Ne*exp(ee)) then
         reject=.FALSE.
      end if
      ee=exp(ee)
   end do
end subroutine selecPhotEBL
 
!------------------------------------------------------------------------------ 
! Draw a photon from the CMB
!------------------------------------------------------------------------------
subroutine selecPhotCMB(ee,zint)
   use Constantes, only : prec, Tcmbp
   implicit none
   real(prec), intent(in):: zint
   real(prec), intent(out):: ee
   real(prec):: cc,d,a1,a2,a3,a4,ran0,zeta
   integer:: mm

   a2=ran0()
   a3=ran0()
   a4=ran0()
   d=a2*a3*a4

   a1=ran0()
   a1=1.202d0*a1
   mm=1
   zeta=1
   do while (a1.ge.zeta)
      mm=mm+1
      cc=mm
      zeta=zeta+1/cc**3
   end do
   cc=mm
   ee=-Tcmbp*(1+zint)/cc*log(d)
end subroutine selecPhotCMB 

subroutine genx(ee,emin,zint)
!    e,emin,kt in mec^2 units   
   use Constantes, only : prec, Tcmbp
   implicit none
   real(prec), intent(in):: zint, emin
   real(prec), intent(out):: ee
   real(prec):: x=0,x0,ran0,a1,g
   real(prec):: dx1,dx00,b,dx0
   real(prec):: ktr,xt,yt,deps

   ktr=Tcmbp*(1+zint)

   x0=emin/ktr
   xt=1
   yt=0
   do while(xt.gt.yt)
      a1=ran0()
      dx00=-log(a1)
      dx0=0
      dx1=dx00
      b=1+X0+x0*x0/2
      deps=1d-6
      do while (abs((dx1-dx0))/dx1.gt.deps)
         dx0=dx1
         x=x0+dx0
         g=(1.+x+x*x/2)/b
         dx1=dx00+log(g)
      enddo
      yt=(1-exp(-x0))/(1-exp(-x))
      xt=ran0()
   enddo
   ee=x*ktr
end subroutine genx

subroutine renorme(u,v,w)
   use Constantes, only : prec
   implicit none
   real(prec), intent(inout) :: u,v,w
   real(prec) :: norm
   norm=sqrt(u**2+v**2+w**2)
   u=u/norm
   v=v/norm
   w=w/norm
end subroutine renorme

! Rotation of theta around (u,v,w) of the vector (u1,v1,w1)
!===========================================================
subroutine rotate(u1,theta,u)
   use Constantes, only : prec
   implicit none
   real(prec), intent(in) :: theta
   real(prec), intent(inout) :: u(3)
   real(prec), intent(out) :: u1(3)
   real(prec) :: u2(3)
   real(prec) :: cx, cy, cz, dx, dy, dz, ex, ey, ez, cos_theta, sin_theta

   call renorme(u(1),u(2),u(3))
   ! Definition of the rotation matrix (given in https://en.wikipedia.org/wiki/Rotation_matrix) 
   !      / cx cy cz \
   !  R = | dx dy dz |
   !      \ ex ey ez /

   cos_theta=cos(theta)
   sin_theta=sin(theta)

   cx = cos_theta+(1-cos_theta)*u(1)**2
   cy = u(1)*u(2)*(1-cos_theta)-u(3)*sin_theta
   cz = u(1)*u(3)*(1-cos_theta)+u(2)*sin_theta
   dx = u(1)*u(2)*(1-cos_theta)+u(3)*sin_theta
   dy = cos_theta+(1-cos_theta)*u(2)**2
   dz = u(2)*u(3)*(1-cos_theta)-u(1)*sin_theta
   ex = u(1)*u(3)*(1-cos_theta)-u(2)*sin_theta
   ey = u(2)*u(3)*(1-cos_theta)+u(1)*sin_theta
   ez = cos_theta+(1-cos_theta)*u(3)**2
      
   u2(1)=u1(1)
   u2(2)=u1(2)
   u2(3)=u1(3)
   ! Apply the rotation to (u1,v1,w1)
   u1(1) = u2(1)*cx+u2(2)*cy+u2(3)*cz
   u1(2) = u2(1)*dx+u2(2)*dy+u2(3)*dz
   u1(3) = u2(1)*ex+u2(2)*ey+u2(3)*ez

   call renorme(u1(1),u1(2),u1(3))

end subroutine rotate

subroutine pdtVect (u1,u2,u3)
   use Constantes, only: prec
   implicit none
   real(prec), intent(in)  :: u1(3),u2(3)
   real(prec), intent(out) :: u3(3)
   ! u3   = u1   X u2
   ! |ux3   |ux1   |ux2
   ! |uy3 = |uy1 X |uy2
   ! |uz3   |uz1   |uz2
   u3(1) = u1(2)*u2(3)-u1(3)*u2(2)
   u3(2) = u1(3)*u2(1)-u1(1)*u2(3)
   u3(3) = u1(1)*u2(2)-u1(2)*u2(1)
end subroutine pdtVect
