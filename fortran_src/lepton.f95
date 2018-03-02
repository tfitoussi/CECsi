#include "../temp/preprocessing.f95"
!===============================================================================         
                           MODULE Lepton_Routines
!===============================================================================         
use constantes, only: prec
implicit none
!===============================================================================         
                                   CONTAINS
!===============================================================================         
subroutine Check_if_lepton_interact(lepton, Zfinal)
   !----------------------------------------------------------------------------
   ! Check if the lepton will interact (logical willInteract) 
   ! or if it reaches the detector ?
   !----------------------------------------------------------------------------
   use Constantes, only : currentEnergy, ComptonEnergyThreshold, m, c, &
                     distSourcetoEarth, H0, omegaL, omegaM, XiComp, Mpc
   use Particles, only : particle, duetoCMB, willInteract
   implicit none
   type (particle), intent(in) :: lepton
   type (particle) :: lepton_modify
   real(prec), intent(out) ::  Zfinal
   real(prec) :: distToSource,zint,redzmax,redzmin, tauint, Zlim
   real(prec) :: tau=0, qgauss, ProbaInt, ProbaLim, ran0
   real(prec) :: aa, bb, dtaudx, deltaZ, f, fprime, taulim
   real(prec) :: dtaudxPPCMB,dtaudxPPEBL,integrandTauCom,PCMB!,PEBL
   real(prec) :: integranddtaudxComCMB,integranddtaudxComEBL, qgauss1
   external qgauss,qgauss1,integrandTauCom, ran0
   external integranddtaudxComCMB,integranddtaudxComEBL,integranddtaudxCom
   integer :: n

   ! Compute Compton Accumulation if lepton energy beyond the threshold
   !----------------------------------------------------------------------------
   XiComp=1
   if (lepton%energy < ComptonEnergyThreshold+1) then
      XiComp=(ComptonEnergyThreshold+2)/(lepton%energy+1)
   end if


#ifdef z_0_limit   
   Zlim = 0
#else
   ! Compute redshift of the lepton when arrived in the detection area without interaction
   !----------------------------------------------------------------------------
   redzmin = -1
   redzmax = lepton%redshift
 
   distToSource = sqrt(sum(lepton%pos**2))
   zint = -2 ! precaution
   ! source X------------------- d -----------------------------X earth
   ! source X------------- d' ---------------X Lepton 
   ! when d' = d (comobile distance) => zlim
   ! Principle: dichotomy
   do while( abs(distToSource/distSourceToEarth-1) > (1.d-8))
      zint = (redzmin+redzmax)/2
      if (redzmin == zint .or. redzmax == zint) then
         distToSource = distSourceToEarth
      else
         ! We reinitialize a copy of the lepton in which we can modify its properties
         ! without affected the "real" lepton
         lepton_modify = lepton
         call MoveLepton0(lepton_modify,zint) 
         distToSource = sqrt(sum(lepton_modify%pos**2))
         if (redzmin == zint .or. redzmax == zint) then
            distToSource = distSourceToEarth
         end if
         if (distToSource < distSourceToEarth) then
            redzmax = zint
         else
            redzmin = zint
         end if
      end if
   end do
   Zlim = zint
#endif

   ! compute the redshift of interaction 
   !----------------------------------------------------------------------------
   ! Compute the redshift of interaction is a lot consumming of time. If the particle 
   ! interact beyond the detector (Earth) it is only waste time.  
   ! Basically, we compute the probability of interaction at the redshift limit (P0)
   ! and we draw a probability of interaction for the lepton (P). if P < P0,
   ! we assume that the lepton will not interact, otherwise we compute the
   ! redshift of interaction.
   currentEnergy = currentEnergy-1
   taulim = qgauss(integrandTauCom, Zlim, lepton%redshift)
   taulim = taulim * c/H0
#ifdef approx_lambda_ic
   if (taulim < 1) then ! No interaction before reaching the limit
      willInteract = .FALSE.
   else
      willInteract = .TRUE.
      tauint = 1
   end if 
#else
   tauint = -log(1-ran0())
   if (tauLim < tauint ) then ! No interaction before reaching the limit
      willInteract = .FALSE.
   else
      willInteract = .TRUE.
   end if
#endif

   ! We create a copy of the lepton in which we can modify its properties
   ! without affected the "real" lepton
   if (.not. willInteract) then ! No interaction before reaching the limit
      Zfinal = Zlim
   else
      aa=log(1.0D-17) !erg
      bb=log(3.0D-11) !erg

      dtaudx=qgauss1(integranddtaudxCom, aa, bb, lepton%redshift)
      deltaZ=-tauint*(H0/c)*((1+lepton%redshift)*sqrt((omegaM*((1+lepton%redshift)**3)+omegaL)))/dtaudx
      if (abs(deltaZ)<1d-6) then ! lepton really close to Earth
         Zfinal=lepton%redshift+deltaZ
      else ! we use Newton method
         zint = lepton%redshift
         n=0
         tau=taulim
         do while(abs(tau-tauint) > (1.d-4)*abs(tauint))
            n=n+1
            !zint = (redzmin+redzmax)/2.
            tau = qgauss(integrandTauCom, zint, lepton%redshift)
            tau = tau*c/H0
            f = tau - tauint
            fprime = - integrandTauCom(zint)*c/H0
            zint = zint - f/fprime
            zint=max(zint,zlim)
            !if(isnan(zint)) then
            if(zint==zlim) then
               print*, n, zint, taulim, tau, tauint, fprime
               print*, lepton%Energy, lepton%dir
               print*, lepton%time, lepton%pos
               print*, lepton%weight
               stop
            end if
            if (n > 25) then 
               ! sometimes Newton oscillate around zint without enough precision
               ! then take the last value as best estimation rather than reduce
               ! the precision for everyone
               !print*, zint, tau, tauint, f, f/tauint
               tau = tauint
            end if
            !if (redzmax == zint .or. redzmin == zint) then
            !   tau = tauint
            !end if
            !if ((tau-tauint) < 0. ) then
            !   redzmax=zint
            !else if ((tau-tauint) > 0.) then
            !   redzmin=zint
            !end if
         end do
         Zfinal = zint
         !stop
         !if(zint <1.d-4) then
         !   willInteract = .FALSE.
         !end if
      end if 
      !print*, zfinal, deltaZ,tauint, dtaudx

#ifdef approx_no_EBL_ic
      duetoCMB= .TRUE.
#else
      ! Is the interaction come from a EBL lepton or a CMB lepton?
      ! Important to compute the pair production (globale variable: logical duetoCMB)
      dtaudxPPCMB=qgauss1(integranddtaudxComCMB, aa, bb, Zfinal)
      dtaudxPPEBL=qgauss1(integranddtaudxComEBL, aa, bb, Zfinal)
      PCMB=dtaudxPPCMB/(dtaudxPPCMB+dtaudxPPEBL)
      !PEBL=dtaudxPPEBL/(dtaudxPPCMB+dtaudxPPEBL)

      ProbaInt=ran0() ! random probability to select if it comes from CMB or EBL
      if (ProbaInt < PCMB) then
         duetoCMB= .TRUE.
      else
         duetoCMB= .FALSE.
      end if
#endif
   end if
   currentEnergy = currentEnergy+1

end subroutine Check_if_lepton_interact

subroutine MoveLepton(lepton, zfinal)
   use Particles, only : particle
   use Constantes , only : prec, a0, c, H0, Mpc
   use EGMF, only : lambdaEGMF
   implicit none
   type (particle), intent(inout) :: lepton
   real(prec), intent(in) :: zfinal
   real(prec) :: Dz_max, Dz
   integer :: i, Nstep, N=1

   Dz_max = -lambdaEGMF*Mpc*H0/c !/N ! < 0
   Dz = zfinal-lepton%redshift      ! < 0
   Nstep = ceiling(Dz/Dz_max)
   do i=1, Nstep
      call MoveLepton0(lepton, lepton%redshift+Dz/Nstep)
   end do

end subroutine MoveLepton

subroutine MoveLepton0(lepton, zfinal)
   ! Compute the distance travel by the particle and modify its position and energy
   !----------------------------------------------------------------------------
   use Particles, only : particle
   use Constantes , only : prec, a0, c, H0
   use EGMF, only : selectEGMF, wLarmor
   implicit none
   type (particle), intent(inout) :: lepton
   real(prec), intent(in) :: zfinal
   real(prec) :: gammaI, gammaF, eta, zini, B(3)
   real(prec) :: VV, theta, sin2_theta2, sin_theta, sinc_theta, dl_sinc_sin_theta
   real(prec) :: qgauss2, integrandeta, e1(3), e2(3), dl_sinc_theta_minus_1
   external :: qgauss2, integrandeta

   ! 0. Choose the magnetic field normed direction (Bx, By, Bz) ------------
   ! first approx: the lepton move in the cube where he is at the beginning
   ! => constant magnetic field
   call selectEGMF(lepton%pos,B)

   zini = lepton%redshift
   gammaI = lepton%Energy
   VV = c * sqrt(gammaI**2-1) /(zini+1)

   eta=qgauss2(integrandeta,zini,zfinal,zini,gammaI)/H0     
   theta = real(lepton%charge)*wLarmor*eta
   sin_theta = sin(theta)
   sin2_theta2 = sin(theta/2)**2
   if (theta <1d-50) then
      sinc_theta = 1
   else
      sinc_theta = sin(theta)/theta
   end if

   call pdtVect(B,lepton%dir,e1)
   call pdtVect(B,e1,e2)

   if (abs(theta) > 1d-100) then
      lepton%pos=lepton%pos+eta*VV*(lepton%dir+2*sin2_theta2/theta*e1+(1-sinc_theta)*e2)
   else
      dl_sinc_sin_theta = theta/2 - (1./3)*(theta/2)**3 + (2./45)*(theta/2)**5
      dl_sinc_theta_minus_1 = -(1./6)*(theta/2)**2 + (1./120)*(theta/2)**4 
      lepton%pos=lepton%pos+eta*VV*(lepton%dir+dl_sinc_sin_theta*e1-dl_sinc_theta_minus_1*e2)
   end if
   lepton%dir = lepton%dir + sin_theta * e1 + 2*sin2_theta2 * e2
   call renorme(lepton%dir(1),lepton%dir(2),lepton%dir(3))
   
   gammaF = sqrt(1+(gammaI**2-1)*((1+zfinal)/(1+zini))**2)
   lepton%energy = gammaF
   
   lepton%redshift = zfinal

end subroutine MoveLepton0

subroutine Inverse_Compton_scattering(lepton,photon)
   use particles, only : particle, emission_redshift
   use Constantes, only : prec, errorId, nbCompton, nbErrors, pi, xiComp
   implicit none
   type (particle), intent(inout) :: lepton, photon
   real(prec) :: EleptBefore, sinalpha
   real(prec) :: alphaa, phi, ran0, deltaE
   real(prec) :: arrivingDir(3),vectPdt(3)
   external :: ran0

   !$OMP ATOMIC
   nbCompton =nbCompton+1

   EleptBefore = lepton%energy-1
   photon%charge=0
   photon%time=lepton%time
   photon%pos=lepton%pos
   photon%energy=0
   photon%dir=lepton%dir
   photon%redshift=lepton%redshift

   lepton%energy=lepton%energy-1 !kinematic energy, Elept adim
   arrivingDir=lepton%dir

   call isophotcompton(lepton%energy,lepton%dir(1),lepton%dir(2),lepton%dir(3),&
      photon%energy,photon%dir(1),photon%dir(2),photon%dir(3),lepton%redshift)

   !Compton accumulation
   vectPdt(1)=arrivingDir(2)*lepton%dir(3)-arrivingDir(3)*lepton%dir(2)
   vectPdt(2)=arrivingDir(3)*lepton%dir(1)-arrivingDir(1)*lepton%dir(3)
   vectPdt(3)=arrivingDir(1)*lepton%dir(2)-arrivingDir(2)*lepton%dir(1)
   sinalpha=sqrt(vectPdt(1)**2+vectPdt(2)**2+vectPdt(3)**2)
   alphaa=asin(sinalpha)
   lepton%dir=arrivingDir
   if ( alphaa .ne. 0 ) then
      call rotate(sqrt(XiComp)*alphaa,vectPdt(1),vectPdt(2),vectPdt(3),&
                     lepton%dir(1),lepton%dir(2),lepton%dir(3))
      phi = ran0()*2*pi
      call rotate(phi,arrivingDir(1),arrivingDir(2),arrivingDir(3),&
                     lepton%dir(1),lepton%dir(2),lepton%dir(3))
   end if
   !Compton accumulation on the energies
#ifdef approx_Eic
   photon%energy = 3.3392426539537134 *1d9/510998.918  *(EleptBefore *510998.918d-12)**2 &
                     *(1+emission_redshift)
   deltaE = photon%energy 
#else
   deltaE = EleptBefore-lepton%energy
#endif
   lepton%energy=EleptBefore-deltaE*xiComp
   lepton%energy=lepton%energy+1

   photon%weight=lepton%weight*xiComp

   photon%generation=lepton%generation+1

   !NbPhotonsProd=NbPhotonsProd+photon%weight
   !write(*,*) arrivingDir(1),arrivingDir(2),arrivingDir(3)
   !write(*,*) lepton%dir(1),lepton%dir(2),lepton%dir(3)
   !write(*,*) photon%dir(1),photon%dir(2),photon%dir(3)

   if (lepton%energy<1) then
      ! $OMP CRITICAL Errors
      nbErrors = nbErrors + 1
      write(errorID,*) 'error ', nbErrors, 'Compton accumulation threshold too high'
      write(errorID,*) lepton
      write(errorID,*) photon
      lepton%energy=1+1d-10
      ! $OMP end CRITICAL Errors
   end if
end subroutine Inverse_Compton_scattering

!--------------------------------------------------------------------
!      Compton process
!--------------------------------------------------------------------
subroutine isophotcompton(e1,u1,v1,w1,ep,up,vp,wp,zint)
! assuming isotropic photon field
! selects incident angle and energy of photon
! for compton interaction
! interaction is then simulated by a call to subroutine compton

! note: all ref change in this routine and its kids have not been rechecked
   use Constantes, only :  pi, prec
   use particles, only : particle, duetoCMB
   implicit none
   real(prec) :: e1,u1,v1,w1
   real(prec) :: e2,u2,v2,w2,zint
   real(prec) :: ep,up,vp,wp,ran0
   real(prec) :: csi,mu,smu
   logical rejec
   real(prec) :: phi,cp,sp,cx,sx
   real(prec) :: yt,xt,u3,v3,w3
   real(prec) :: gama,beta,ksi
   real(prec) :: coef,e3
   external ran0

   ! cis angle between the projetion of r in the ex ez plane and the ez vector
   gama=e1+1
   rejec=.true.
   xt=1
   yt=0 

   ! TARGET ENERGY
   do while (xt.gt.yt)
      if (duetoCMB ) then ! si compton est du au CMB
         call selecphotCMB(e2,zint)
         !! BB truncated at low energy (to check)
      else  ! si le Compton est du à l'EBL
        call selecphotEBL(e2,zint)
      end if
          
      yt=yu(e1,e2)
      ! compute the angle averaged Compton cross section,
      xt=ran0()
      ! used for rejection test
   enddo

   !write(*,*) e1, e2, e1/(e1+e2)*100, duetoCMB

   ! TARGET DIRECTION
   call genmu(e1,e2,mu)
   ! select angle between electron and photon (isotropic photons)
   phi=2*pi*ran0()

   cp=cos(phi)
   sp=sin(phi)
   smu=sqrt(1-mu**2)
   if (w1==0) then 
      csi=pi/2
   else
      csi=atan(u1/w1)                  
   end if
   cx=cos(csi)
   sx=sin(csi)      

   ! determine incident photon direction in global frame:  
   u2=smu*(cp*cx+sp*sx*v1)+mu*u1
   v2=-smu*sp*(w1*cx+u1*sx)+mu*v1
   w2=smu*(-cp*sx+sp*cx*v1)+mu*w1

   beta=sqrt(e1*(e1+2))/gama
   coef=1-mu*beta
   ksi=2*e2*gama*coef
   e3=e1
   ep=e2
   up=u2
   vp=v2
   wp=w2
   u3=u1
   v3=v1
   w3=w1

   call compton(up,vp,wp,ep,e3,ksi,u3,v3,w3)

   call renorme(up,vp,wp) 
   call renorme(u3,v3,w3) 
   e1=e3
   u1=u3
   v1=v3
   w1=w3

end subroutine isophotcompton

subroutine genmu(ec,e,mu)
! For isotropic monoenergetic photon gas,
!  generates mu=cos(theta) theta= angle between incident electron and photon. 
   use Constantes, only : prec
   implicit none
   real(prec) :: gama,e,mu,xt,y
   real(prec) :: xmin,xmax,f0,lnf0
   real(prec) :: beta,fmin,fmax,a1,ran0
   real(prec) :: x0,gfdex,a,b,lnfx,ec

   gama=ec+1
   beta=sqrt(ec*(ec+2))/gama

   if (1-beta.gt.1d-5) then
      xmin=2*gama*e*(1-beta)
   else 
      xmin=e/gama
   endif
   xmax=2*gama*e*(1+beta)

   if (xmax.lt.0.1d0) then
   !  KN effects are weak: use simple rejection technique
      y=0
      xt=1
      do while (xt.gt.y) 
         a1=ran0()
         mu=(1-sqrt((1+beta)**2-4*beta*a1))/beta
         if (abs(mu).gt.1) then
            write(6,*) 'error mu',a1,beta,mu
            mu=mu/abs(mu)
         endif
         x0=2*E*gama*(1-mu*beta)
         y=0.75d0*xsx(x0)/x0
         xt=ran0()
      enddo
   else
      ! KN effects are strong
      ! numerical inversion of an aproximation to the full distribution
      ! then rejection.
      fmin=gfx(xmin)
      fmax=gfx(xmax)

      y=0
      xt=1
      do while (xt.gt.y) 
         a1=ran0()
         F0=a1*(fmax-fmin)+fmin
         lnf0=log(f0)        
         x0=2*e*gama         
         gfdex=gfx(x0)        

!      solving gfx(x)=f0 using using iterative method
         do while (abs(gfdex/f0-1).gt.1d-8) 
            lnfx=log(gfdex)
            a=x0*fx(x0)/gfdex
            b=lnfx-lnf0-a*log(x0)
         
            x0=exp(-b/a)
            gfdex=gfx(x0)
         enddo
         y=xsx(x0)/fx(x0)
         xt=ran0()
      enddo
      mu=(1-x0/2/gama/e)/beta
   endif
end subroutine genmu

!*******************************************************************
subroutine compton(u,v,w,e,ec,xi,dex,dey,dez)
!*******************************************************************
   use Constantes, only : m,c,pi, prec
   implicit none
   real(prec) :: u,v,w,e,ran0
   real(prec) :: beta,gama,dex,dey,dez
   real(prec) :: ud,vd,wd,a1,a2,mud,q,phi
   real(prec) ::  ro,y,xi,de,eps,xdsxi,vv
   real(prec) :: gama2,e2
   real(prec) :: mu0,vv2,fact, dex2,dey2,dez2
   real(prec) :: csphi,snphi,smud,t,csphismud,snphismud
   real(prec) :: ec
   
!     common /mc2/mc2 
   gama=ec+1
   beta=sqrt(ec*(ec+2))/gama

!     direction du photon diffuse (par test de rejection)
   a1=1.
! le bug etait la :
   y=0.
! la ligne precedente avait ete supprimmee
   do while (2*a1.ge.y)
      
      a1=ran0()
      a2=ran0()
      phi=2*pi*a2
      
      mu0=2*a1-1
      !mud=min((beta+mu0)/(1+beta*mu0),1.)
      !t=gama*(1.-mud*beta)
      mud=1-(1-mu0)/(gama**2 *(1+beta)*(1+beta*mu0))
      t=1/((1+beta)*gama)*(1+(1-mu0)*beta/(1+beta*mu0))
!       min( added on 17/4/09 to avoid apparition of  NaN in some rare cases when comuting smud
      smud=sqrt(1-mud**2)
      
      csphi=cos(phi)
      snphi=sin(phi)
      csphismud=csphi*smud
      snphismud=snphi*smud
      ro=(dex**2+dey**2)**5.d-1
      
      ud=mud*dex+ &
          (dey*csphismud+dex*dez*snphismud)/ro
      vd=mud*dey+ &
         (-dex*csphismud+dey*dez*snphismud)/ro

      wd=mud*dez-snphismud*ro
      
      q=u*ud+v*vd+w*wd
      
      eps=e*(1-q)/(t)
      xdsxi=1/(1+eps)
      
      
      a1=ran0()
      
      y=4*((eps/xi)**2-eps/xi)*xdsxi**2+xdsxi+xdsxi**3
      
   enddo
!      write(6,*) 'sc:',smud,mud,mu0,beta,ro,y
!     energie du photon diffuse
   
   e2=xi*xdsxi/(2*t)
   
   de=e-e2
   gama2=de+gama

   vv=c*sqrt(1-gama**(-2))   
   vv2=c*sqrt(1-gama2**(-2))
   fact=1/(gama2*m*vv2)
   dex2=fact*(m*c*c*e/c*u+gama*m*vv*dex-m*c*c*e2/c*ud)
   dey2=fact*(m*c*c*e/c*v+gama*m*vv*dey-m*c*c*e2/c*vd)
   dez2=fact*(m*c*c*e/c*w+gama*m*vv*dez-m*c*c*e2/c*wd)

   gama=gama2
   ec=ec+de
   u=ud
   v=vd
   w=wd
   e=e2
   dex=dex2
   dey=dey2
   dez=dez2
end subroutine compton

function yu(e1,e)
   use Constantes, only : prec
   implicit none
   !   angle averaged  Compton cross-section for lepton of energy e1
   !   and photon of energy e2 
   !   (in units of sigma_T)
   real(prec) :: yu,g,e,e1,b
   real(prec) :: xx,xxm,xxp,umb,xst

   g=e1+1
   b=sqrt(e1*(e1+2))/g
   if (g.gt.1d2) then
      umb=1/2/g**2
   else
      umb=1-b
   endif
   xx=2*e*g
   xxp=xx*(1+b)
   xxm=xx*umb
   xst=(phix(xxp)-phix(xxm))/xx**2
   xst=xst*3/(8*b)
   yu=xst
end function yu

function phix(x)
   use Constantes, only : prec
   implicit none
! phix=\int_0^x (x*sigma(x))dx   (PSS83 p 319)
   real(prec) :: x,phix
   if (x.lt.0.5d0) then 
      phix=x*x/6+0.047d0*x*x*x-0.03d0*x**4+0.5d0*x*x/(1+x)
   else 
      if (x.lt.3.5d0) then 
         phix=(1+x)*log(1+x)-0.94d0*x-9.25d-3
      else 
         phix=(1+x)*log(1+x)-0.5d0*x*log(2+0.076d0*x)+9.214d0
      endif
   endif
end function phix

function xsx(x)
   use Constantes, only : prec
   implicit none
! x*sigma(x)  see PSS83, p319
   real(prec) :: x,xsx
   if (x.gt.0.5d0) then
    xsx=(1-4/x-8/x/x)*log(1+x)+0.5d0
    xsx=xsx+8/x-1/2/(1+x)**2
   else
     xsx=1/3+0.141d0*x-0.12d0*x*x+(1+0.5d0*x)/(1+x)**2
     xsx=xsx*x
   endif
end function xsx
     
function gfx(x)
   use Constantes, only : prec
   implicit none
!   primitive of fx see below
   real(prec) :: gfx,x,u
   u=1.+x
   gfx=u*log(u)+0.5d0*(exp(-x)-x-1)
end function gfx

function fx(x)
   use Constantes, only : prec
   implicit none
! 0.7<xsx/fx<1.     fx-> xsx when x-> \infty  
! fx is used to model the Compton incident angle by rejection tecnique
   real(prec) :: fx,x
   fx=log(1+x)+0.5d0*(1-exp(-x))
end function fx

!===============================================================================         
                        END MODULE Lepton_Routines
!===============================================================================         
