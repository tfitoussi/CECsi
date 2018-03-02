#include "../temp/preprocessing.f95"
!===============================================================================
!     Here are put all the functions and routines necessary to compute
!     integrals in the program
!===============================================================================

!********************************************************************************
!* Calculation of GAUSS-LEGendRE abscissas and weights for Gaussian Quadrature
!* integration of polynomial functions.
!*      For normalized lower and upper limits of integration -1.0 & 1.0, and
!* given n, this routine calculates, arrays xabsc(1:n) and  weig(1:n) of length n,
!* containing the abscissas and weights of the Gauss-Legendre n-point quadrature
!* formula.  For detailed explanations finding weights & abscissas, see
!* "Numerical Recipes in Fortran */
!********************************************************************************
subroutine  gauleg(ngp, xabsc, weig)
   use Constantes, ONLY : prec,pi
   implicit none
   integer  i, j, m
   real(prec)  p1, p2, p3, pp, z, z1
   integer, intent(in) :: ngp            ! # of Gauss Points
   real(prec), intent(out) :: xabsc(ngp), weig(ngp)
   real(prec)  :: eps
   parameter (eps=3.0d-15)        !EPS is the relative precision

   m = (ngp + 1) / 2
   !* Roots are symmetric in the interval - so only need to find half of them  */
   do i = 1, m    ! Loop over the desired roots */
      z = cos( pi * (i-0.25) / (ngp+0.5) )
      !*   Starting with the above approximation to the ith root,
      !*   we enter the main loop of refinement by NEWTON'S method   */
      100  p1 = 1.0
      p2 = 0.0
      !*  Loop up the recurrence relation to get the Legendre
      !*  polynomial evaluated at z                 */

      do j = 1, ngp
         p3 = p2
         p2 = p1
         p1 = ((2.0*j-1.0) * z * p2 - (j-1.0)*p3) / j
      enddo

      !* p1 is now the desired Legendre polynomial. We next compute pp,
      !* its derivative, by a standard relation involving also p2, the
      !* polynomial of one lower order.      */
      pp = ngp*(z*p1-p2)/(z*z-1.0)
      z1 = z
      z = z1 - p1/pp             ! Newton's Method  */

      if (abs(z-z1) .gt. EPS) GOTO  100

      xabsc(i) =  - z                     ! Roots will be bewteen -1.0 & 1.0 */
      xabsc(ngp+1-i) =  + z                 ! and symmetric about the origin  */
      weig(i) = 2.0/((1.0-z*z)*pp*pp) ! Compute the weight and its       */
      weig(ngp+1-i) = weig(i)               ! symmetric counterpart         */

   end do     ! i loop

End subroutine gauleg

!********************************************************************************
!*    GAUSS LEGendRE inTEGRATION
!*     Returns the SinGLE integral of the function (of ONE VARIABLE) "func"
!* between x1 and x2 by N-point Gauss-Legendre integration. The function
!* is evaluated exactly N times at interior points in the range of
!* integration.       */
!********************************************************************************
recursive function qgauss(func, x1, x2) result(intgrl)
   use Constantes, only : prec, ngp, xabsc, weig
   implicit none
   real(prec)  intgrl, x1, x2, func
   real(prec)  xm, xl
   integer j

   if (x1==x2) then
      intgrl=0
   else if (abs(x1-x2) < 1.d-15) then
      intgrl= (x2-x1)*func(x1)
   else
      intgrl = 0
      xm = 0.5d0 * (x2 + x1)
      xl = 0.5d0 * (x2 - x1)
      do j = 1, ngp
         intgrl = intgrl + weig(j) * func( xm + xl*xabsc(j))
      end do
      intgrl = intgrl * xl;    !Scale the answer to the range of integration  */
   end if

end function qgauss

recursive function qgauss1(func, x1, x2, redzz) result(intgrl)
   use Constantes, only : prec, ngp, xabsc, weig
   implicit none
   real(prec) :: intgrl, x1, x2, func
   real(prec) :: xm, xl, redzz
   integer j

   if (x1==x2) then
      intgrl=0
   else if (abs(x1-x2) < 1.d-15) then
      intgrl= (x2-x1)*func(x1,redzz)
   else
      intgrl = 0
      xm = 0.5d0 * (x2 + x1)
      xl = 0.5d0 * (x2 - x1)
      do j = 1, ngp
         intgrl = intgrl + weig(j) * func( xm + xl*xabsc(j), redzz)
      end do
      intgrl = intgrl * xl;    !Scale the answer to the range of integration  */
   end if

end function qgauss1

recursive function qgauss2(func, x1, x2, z1,gamma1) result(intgrl)
   use Constantes, only : prec, ngp, xabsc, weig
   implicit none
   real(prec) :: intgrl, x1, x2, func
   real(prec) :: xm, xl,z1,gamma1
   integer j

   if (x1==x2) then
      intgrl=0
   else if (abs(x1-x2) < 1.d-15) then
      intgrl= (x2-x1)*func(x1,z1,gamma1)
   else
      intgrl = 0
      xm = 0.5d0 * (x2 + x1)
      xl = 0.5d0 * (x2 - x1)
      do j = 1, ngp
         intgrl = intgrl + weig(j) * func( xm + xl*xabsc(j), z1,gamma1)
      end do
      intgrl = intgrl * xl;    !Scale the answer to the range of integration  */
   end if

end function qgauss2


!********************************************************************************
!*    inTEGRAND useD
!********************************************************************************
!  In photon computation
!==========================
real(prec) function integranddtaudxPP (y,redzz)
   use Constantes, only : prec, k, m, c, Tcmbp, Tcmb, pi,&
                           currentEnergy, currentRedshift
   implicit none
   external NeInterpol 
   real(prec) :: eps, phiBar, n, serie, serieAvant, pres, beta0, w0, kk, s0
   real(prec) :: y, Egammaa, TcmbpC, redzz, nCMB, nEBL,theta

   eps=exp(y)
   Egammaa=currentEnergy*(1+redzz)/(1+currentRedshift) ! redz redshift of the initial particle
   s0=eps*Egammaa    ! without dimension
   beta0=sqrt(1-1/s0)
   w0=(1+beta0)/(1-beta0)
   pres=1.D-5  ! precision of the computation
   TcmbpC=Tcmbp*(1+redzz)

   ! Calculation of phiBar (cf article Gould Schreder 1967)
   serie=0
   kk=0                 
   serieAvant=10 ! start the while statement
   do while(abs(serie-serieAvant)>pres*serieAvant)
      serieAvant=serie
      kk=kk+1
      serie = serie + ((-1)**(kk-1))*(kk**(-2))*(w0**(-kk))
   end do
   phiBar=(1.+beta0**2)*s0*log(w0)-(beta0**2)*log(w0) -log(w0)**2 -4*beta0*s0 &
      +2*beta0+4*log(w0)*log(w0+1)+4*(-(0.5d0)*((log(w0))**2)-(pi**2)/12+serie)

   ! Photons distribution
   call NeInterpol(eps, redzz, nEBL)  !EBL
   if (eps/(TcmbpC)<1D-2) then
      theta=eps/(Tcmb*k*(1+redzz)/(m*c*c))
      nCMB=((eps/pi)**2)*(theta+(theta**2)/2+(theta**3)/6+(theta**4)/24)**(-1)
   else
      nCMB=((eps/pi)**2)*((exp(eps/TcmbpC)-1)**(-1)) ! CMB (blackbody)
   end if
   n=nCMB+nEBL

   integranddtaudxPP=n*phiBar*exp(-y)
   return
end function integranddtaudxPP

real(prec) function integranddtaudxPPCMB (y,redzz)
   use Constantes, only : prec, k, m, c, Tcmbp, Tcmb, pi,&
                           currentEnergy, currentRedshift
   implicit none
   real(prec) :: eps, phiBar, serie, serieAvant, pres, beta0, w0, kk, s0
   real(prec) :: y, Egammaa, TcmbpC, redzz, nCMB, theta

   eps=exp(y)
   Egammaa=currentEnergy*(1+redzz)/(1+currentRedshift)
   s0=eps*Egammaa    ! without dimension
   beta0=sqrt((1-1/s0))
   w0=(1+beta0)*(1+beta0)*s0
   pres=1.D-5  ! precision of the computation
   TcmbpC=Tcmbp*(1+redzz)

   ! Calculation of phiBar (cf article Gould Schreder 1967)
   serie=0
   kk=0                 
   serieAvant=10 ! start the while statement
   do  while(abs(serie-serieAvant)>pres*serieAvant)
      serieAvant=serie
      kk=kk+1
      serie = serie + ((-1)**(kk-1))*(kk**(-2))*(w0**(-kk))
   end do
   phiBar=(1+beta0**2)*s0*log(w0)-(beta0**2)*log(w0) -log(w0)**2 -4*beta0*s0 &
      +2*beta0+4*log(w0)*log(w0+1)+4*(-(0.5d0)*((log(w0))**2)-(pi**2)/12+serie)

   ! Photons distribution
   nCMB=((eps/pi)**2)*((exp(eps/TcmbpC)-1)**(-1)) ! CMB (blackbody)
   theta=eps/(Tcmb*k*(1+redzz)/(m*c*c))
   if (eps/(TcmbpC)<1D-2) then
      nCMB=((eps/pi)**2)*(theta+(theta**2)/2.+(theta**3)/6.+(theta**4)/24.)**(-1)
   end if
   
   integranddtaudxPPCMB=((exp(-y)))*nCMB*phiBar
   return
end function integranddtaudxPPCMB

real(prec) function integranddtaudxPPEBL (y,redzz)
   use Constantes, only : prec, Tcmbp, pi, currentEnergy, currentRedshift
   implicit none
   external NeInterpol 
   real(prec) :: eps, phiBar, serie, serieAvant, pres, beta0, w0, kk, s0
   real(prec) :: y, Egammaa, TcmbpC, redzz, nEBL

   eps=exp(y)
   Egammaa=currentEnergy*(1+redzz)/(1+currentRedshift)
   s0=eps*Egammaa    ! without dimension
   beta0=sqrt((1-1/s0))
   w0=(1+beta0)*(1+beta0)*s0
   pres=1.D-5  ! precision of the computation
   TcmbpC=Tcmbp*(1+redzz)

   ! Calculation of phiBar (cf article Gould Schreder 1967)
   serie=0
   kk=0                 
   serieAvant=10 ! start the while statement
   do  while(abs(serie-serieAvant)>pres*serieAvant)
      serieAvant=serie
      kk=kk+1
      serie = serie + ((-1)**(kk-1))*(kk**(-2))*(w0**(-kk))
   end do
   phiBar=(1+beta0**2)*s0*log(w0)-(beta0**2)*log(w0) -log(w0)**2 -4*beta0*s0 &
      +2*beta0+4*log(w0)*log(w0+1)+4*(-(0.5d0)*((log(w0))**2)-(pi**2)/12+serie)

   ! Photons distribution
   call NeInterpol(eps, redzz, nEBL)  !EBL
   
   integranddtaudxPPEBL=((exp(-y)))*nEBL*phiBar
   return
end function integranddtaudxPPEBL

real(prec) function integrandTauPP (reddz)
   use Constantes, only : prec, omegaM, omegaK, omegaL, currentEnergy, currentRedshift
   implicit none
   real(prec) qgauss1
   external qgauss1 
   real(prec) integranddtaudxPP
   external integranddtaudxPP 
   real(prec) :: aa, bb
   real(prec) :: reddz, dtaudxPP, Egammaa

   !redz = borne sup, reddz = variable 
   Egammaa=currentEnergy*(1+reddz)/(1+currentRedshift) 
   aa=log(1/Egammaa) 
   bb=log(max(2.402285564062D-5,1/Egammaa))

   if (aa==bb) then
      dtaudxPP=0
   else
      dtaudxPP=qgauss1(integranddtaudxPP, aa, bb, reddz)
   end if
   integrandTauPP=(1/(Egammaa**2))*dtaudxPP/((1+reddz)*sqrt((omegaM*((1+reddz)**3)+(omegaK)*((1+reddz)**2)+omegaL)))
   return

end function integrandTauPP

!==========================================================================================
! Integrals for the Compton interaction
!==========================================================================================
real(prec) function integranddtaudxCom (y,reddz)
   use Constantes, only : prec, hb, Tcmb, k, m, c, pi, sigmaT, &
                           xiComp, currentEnergy, currentRedshift
   implicit none
   external integrandx, NeInterpol, qgauss
   real(prec) integrandx, qgauss
   real(prec) :: a, b, n, beta, eps, nEBL, nCMB,y, theta
   real(prec) :: x, reddz,xx,w,gamma
  
   eps=exp(y)
   gamma=sqrt(1+(currentEnergy**2-1)*((1+reddz)/(1+currentRedshift))**2)
   theta=eps/(Tcmb*k*(1+reddz))
   if (theta>1D-2) then
      nCMB=((hb*c)**(-3))*((eps/pi)**2)*((exp(theta)-1)**(-1)) !pour CMB corps noir
   else 
      nCMB=((hb*c)**(-3))*((eps/pi)**2)*(theta+(theta**2)/2.+(theta**3)/6.+(theta**4)/24.)**(-1)
   end if

   eps=eps/(m*c*c)
   call NeInterpol(eps, reddz, nEBL) ! Pour photons du CIB, à partir du fichier de Dominguez et al.
   eps=eps*m*c*c
   nEBL=nEBL*(((hb*c)**(-3))*(((m*c*c))**2))

   n=(nCMB+nEBL)/xiComp
  
   w=eps/(m*c**2)
   x=2.*gamma*w
#ifdef approx_Thomson_regime 
   integranddtaudxCom=n*sigmaT*exp(y)
#else
   if ( x < 1.d-3 ) then
      xx = 1 - x + (13./10)*x**2 - (133./80)*x**3 + (143./70)*x**4
   else
      beta=sqrt(1-gamma**(-2))
      b=2*gamma*(1+beta)*w
      a=2*gamma*(1-beta)*w

      xx = qgauss(integrandx,a,b)
      xx = 3 / (32*gamma**2*beta*w**2) * xx
   end if
   integranddtaudxcom = n * sigmat * xx * exp(y)
#endif
end function integranddtaudxCom

real(prec) function integranddtaudxComCMB (y,reddz)
   use Constantes, only : prec, hb, Tcmb, k, m, c, pi, sigmaT, &
                           xiComp, currentEnergy, currentRedshift
   implicit none
   external integrandx, qgauss
   real(prec) integrandx, qgauss
   real(prec) :: a, b, n, beta, eps, nCMB,y, theta
   real(prec) :: x,reddz,xx,w,gamma
  
   eps=exp(y)
   gamma=sqrt(1+(currentEnergy**2-1)*((1+reddz)/(1+currentRedshift))**2)

   theta=eps/(Tcmb*k*(1+reddz))
   if (theta>1D-2) then
      nCMB=((hb*c)**(-3))*((eps/pi)**2)*((exp(theta)-1)**(-1)) !pour CMB corps noir
   else 
      nCMB=((hb*c)**(-3))*((eps/pi)**2)*(theta+(theta**2)/2+(theta**3)/6+((theta)**4)/24)**(-1)
   end if
   n=nCMB/xiComp
   w=eps/(m*c**2)

#ifdef approx_Thomson_regime 
   integranddtaudxComCMB=n*sigmaT*exp(y)
#else
   x=2*gamma*w
   if (x<1.D-3) then
      xx=sigmaT*(1-x+(13/10)*x**2-(133/80)*x**3+(143/70)*x**4)
      integranddtaudxComCMB=n*xx*exp(y)
   else
      beta=sqrt(1-gamma**(-2))
      b=2*gamma*(1+beta)*w
      a=2*gamma*(1-beta)*w

      xx=qgauss(integrandx,a,b)
      integranddtaudxComCMB=(3*sigmaT*n/(32*gamma**2*beta*w**2))*xx*exp(y)
   end if
#endif
   return
end function integranddtaudxComCMB

real(prec) function integranddtaudxComEBL (y,reddz)
   use Constantes, only : prec, hb, Tcmb, k, m, c, sigmaT, &
                           currentEnergy, currentRedshift,xiComp
   implicit none
   external integrandx, NeInterpol, qgauss
   real(prec) integrandx, qgauss
   real(prec) :: a, b, n, beta, eps, nEBL,y, theta
   real(prec) :: TcmbpC, x, reddz,xx,w,gamma

   eps=exp(y)
   TcmbpC=Tcmb*k*(1+reddz)
   theta=eps/(Tcmb*k*(1+reddz))

   eps=eps/(m*c*c)
   call NeInterpol(eps, reddz, nEBL) 
   eps=eps*m*c*c

   nEBL=nEBL*(((hb*c)**(-3))*(((m*c*c))**2))
   n=nEBL/xiComp

   gamma=sqrt(1+(currentEnergy**2-1)*((1+reddz)/(1+currentRedshift))**2)

#ifdef approx_Thomson_regime 
   integranddtaudxComEBL=n*sigmaT*exp(y)
#else
   w=eps/(m*c**2)
   x=2*gamma*w
   if (x<1.D-3) then
      xx=sigmaT*(1-x+(13/10)*x**2-(133/80)*x**3+(143/70)*x**4)
      integranddtaudxComEBL=n*xx*exp(y)
   else
      beta=sqrt(1-gamma**(-2))
      b=2*gamma*(1+beta)*w
      a=2*gamma*(1-beta)*w
      xx=qgauss(integrandx,a,b)
      integranddtaudxComEBL=(3*sigmaT*n/(32*gamma**2*beta*w**2))*xx*exp(y)
   end if
#endif
   return
end function integranddtaudxComEBL

real(prec) function integrandTauCom (reddz)
   use Constantes, only : prec, omegaM, omegaK, omegaL
   implicit none
   external qgauss1, integranddtaudxCom
   real(prec) qgauss1, integranddtaudxCom
   real(prec) :: reddz, dtaudxCom,aa,bb

   aa=log(1.0D-17)
   bb=log(3.0D-11)
  
   dtaudxCom=qgauss1(integranddtaudxCom,aa,bb,reddz)
   integrandTauCom=dtaudxCom/((1+reddz)*sqrt((omegaM*((1+reddz)**3)+omegaK*((1+reddz)**2)+omegaL)))
   return
end function integrandTauCom

!==========================================================================================
! Integrands distance and time
!==========================================================================================
real(prec) function integrandx(x)
   use Constantes, only : prec
   implicit none
   real(prec) :: x, xp1
   xp1 = x+1
   integrandx=(1-4/x*(1+2/x))*log(xp1)+8/x-0.5d0/xp1**2+0.5d0
   return
end function integrandx

real(prec) function integrandt(z)
   use Constantes, only : omegaM,omegaK,omegaL,prec
   implicit none
   real(prec) :: z
   integrandt=-1/( (1 + z) * sqrt(omegaM*((1+z)**3)+omegaK*((1+z)**2)+omegaL))
   return
end function integrandt

real(prec) function integrandtpro (z, z1,gamma1)
   use Constantes, only : omegaM,omegaL,prec
   implicit none
   real(prec) :: z,z1,zp1,gamma1, gamma
   zp1 = z+1
   gamma = sqrt(1+(gamma1**2-1)* (zp1/(1+z1))**2)
   integrandtpro=-1/( zp1 * gamma * sqrt(omegaM*zp1**3+omegaL))
   return
end function integrandtpro

real(prec) function integrandeta (z,z1,gamma1)
   use Constantes, only : omegaM,omegaK,omegaL,prec
   implicit none
   real(prec) :: z,z1,zp1,gamma1, gamma
   zp1 = z+1
   gamma = sqrt(1+(gamma1**2-1)* (zp1/(1+z1))**2)
   integrandeta=-zp1/(gamma*sqrt(omegaM*(zp1**3)+(omegaK)*(zp1**2)+omegaL))
   return
end function integrandeta

real(prec) function integrandR (reddz)
   use Constantes, only : omegaM,omegaK,omegaL,prec
   implicit none
   real(prec) :: reddz
   integrandR=1/(sqrt(omegaM*((1+reddz)**3)+(omegaK)*((1+reddz)**2)+omegaL))
   return
end function integrandR

real(prec) function integrandDr (x, z1,gamma1)
   use Constantes, only : omegaM,omegaK,omegaL,prec
   implicit none
   real(prec) :: x, z1,gamma1,gammasquare,beta
   gammasquare = 1+(gamma1**2-1)* ((1+x)/(1+z1))**2 
   beta = sqrt(1 - 1/gammasquare)
   integrandDr= -beta/(sqrt(omegaM*((1+x)**3)+omegaK*((1+x)**2)+omegaL))
   return
end function integrandDr

subroutine integrandDistXYZ (z, z1,gamma1, integrand) 
   use Constantes, only : omegaM,omegaK,omegaL,prec,H0
   use EGMF, only : wlarmor
   use particles, only: charge
   implicit none
   real(prec) :: z,z1,gamma1, gamma, v, theta
   real(prec) :: qgauss2, integrandtpro, dvdz
   real(prec), intent(out), dimension(3) :: integrand
   external :: qgauss2, integrandtpro

   gamma = sqrt(1+(gamma1**2-1)*((1+z)/(1+z1))**2)
   v = (gamma1/gamma)*((gamma**2-1)/(gamma1**2-1))
   theta = qgauss2(integrandtpro,z1,z,z1,gamma1)/H0
   theta = -charge*wLarmor*theta
   dvdz= -v/((1+z)*sqrt(omegaM*((1+z)**3)+omegaK*((1+z)**2)+omegaL))

   integrand(1)=dvdz
   integrand(2)=dvdz*cos(theta)
   integrand(3)=dvdz*sin(theta)

   return
end subroutine integrandDistXYZ
