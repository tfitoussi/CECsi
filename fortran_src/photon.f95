#include "../temp/preprocessing.f95"
!===============================================================================         
                           MODULE Photon_Routines
!===============================================================================         
use constantes, only: prec
implicit none
!===============================================================================         
                                   CONTAINS
!===============================================================================         
subroutine Check_if_photon_interact(photon, Zfinal)
   use Constantes, ONLY : c, H0, prec, pi, r0, m, hb, omegaM, omegaL, distSourceToEarth
   use Particles, only: particle, dueToCMB, willInteract
   implicit none
   type (particle), intent(in) :: photon
   type (particle) :: photon_modify
   real(prec), intent(out) :: Zfinal
   real(prec) :: distToSource,zint,redzmax,redzmin, comp, Zlim
   real(prec) :: Norma, tau, qgauss, ProbaInt, ProbaLim, ran0
   real(prec) :: aa, bb, dtaudx, deltaZ,integranddtaudxPP
   real(prec) :: dtaudxPPCMB,dtaudxPPEBL,integrandTauPP,PCMB!,PEBL
   real(prec) :: integranddtaudxPPCMB,integranddtaudxPPEBL, qgauss1
   external qgauss,qgauss1,integrandTauPP, ran0
   external integranddtaudxPPCMB,integranddtaudxPPEBL,integranddtaudxPP

   tau=0

#ifdef z_0_limit   
   Zlim = 0
#else
   ! Compute redshift of the photon when arrived in the detection area without interaction
   !=======================================================================================
   redzmin = -1
   redzmax = photon%redshift
 
   distToSource = sqrt(photon%pos(1)**2+photon%pos(2)**2+photon%pos(3)**2)
   ! source X------------------- d -----------------------------X earth
   ! source X------------- d' ---------------X Photon 
   ! when d' = d (comobile distance) => zlim
   ! Principle: dichotomy
   do while( abs(distToSource-distSourceToEarth)/distSourceToEarth > (1.d-15))
      zint = (redzmin+redzmax)/2
      if (redzmin == zint .or. redzmax == zint) then
         distToSource = distSourceToEarth
      else
         ! We reinitialize a copy of the photon in which we can modify its properties
         ! without affected the "real" photon
         photon_modify = photon
         call MovePhoton(photon_modify,zint) 
         distToSource = sqrt(photon_modify%pos(1)**2+photon_modify%pos(2)**2+photon_modify%pos(3)**2)
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
   !=====================================
   ! Compute the redshift of interaction is a lot consumming of time. If the particle 
   ! interact beyond the detector (Earth) it is only waste time.  
   ! Basically, we compute the probability of interaction at the redshift limit (P0)
   ! and we draw a probability of interaction for the photon (P). if P < P0,
   ! we assume that the photon will not interact, otherwise we compute the
   ! redshift of interaction.
   Norma = (c/H0)*pi*(r0**2)*(m*c/hb)**3
   tau = qgauss(integrandTauPP, Zlim, photon%redshift)
   tau = tau * Norma
#ifdef approx_lambda_gg
   if (tau < 1) then ! No interaction before reaching the limit
      willInteract = .FALSE.
   else
      willInteract = .TRUE.
      comp = -1
   end if 
#else
   comp = log(1.-ran0())
   if ( tau < -comp ) then ! No interaction before reaching the limit
      willInteract = .FALSE.
   else
      willInteract = .TRUE.
   end if 
#endif

   if (.not. willinteract) then ! No interaction before reaching the limit
      Zfinal = Zlim
   else
      redzmax = photon%redshift
      redzmin = Zlim

      aa = log(1/photon%energy)
      bb = log(max(2.402285564062D-5,1/photon%energy))
      
      dtaudx=qgauss1(integranddtaudxPP, aa, bb, photon%redshift)
      deltaZ=comp*Norma*((1+photon%redshift)*sqrt((omegaM*((1.+photon%redshift)**3)+omegaL)))/(dtaudx/(photon%energy**2))

      if (abs(deltaZ)<1d-6) then ! photon really close to Earth
         Zfinal=photon%redshift+deltaZ
      else ! we use dichotomy
         do while(abs(tau+comp)  > (1.d-8)*abs(comp))
            zint = (redzmin+redzmax)/2
            tau = qgauss(integrandTauPP, zint, photon%redshift)
            tau = tau*Norma 
            if (redzmax == zint .or. redzmin == zint) then
               tau = -comp
            end if
            if ((tau+comp) < 0 ) then
               redzmax=zint
            else if ((tau+comp) > 0) then
               redzmin=zint
            end if
         end do
         Zfinal = zint
      end if 

      ! Is the interaction come from a EBL photon or a CMB photon?
      ! Important to compute the pair production (globale variable: logical duetoCMB)
      photon_modify%energy=photon%energy*(1+Zfinal)/(1+photon%redshift)
      aa = log(1/photon_modify%energy)
      bb = log(max(2.402285564062D-5,1/photon_modify%energy))
      if (aa==bb) then
         duetoCMB= .FALSE.
      else
         dtaudxPPCMB=qgauss1(integranddtaudxPPCMB, aa, bb, Zfinal)
         dtaudxPPEBL=qgauss1(integranddtaudxPPEBL, aa, bb, Zfinal)
         PCMB=dtaudxPPCMB/(dtaudxPPCMB+dtaudxPPEBL)
         !PEBL=dtaudxPPEBL/(dtaudxPPCMB+dtaudxPPEBL)

         ProbaInt=ran0() ! random probability to select if it comes from CMB or EBL
         if (ProbaInt < PCMB) then
            duetoCMB= .TRUE.
         else
            duetoCMB= .FALSE.
         end if
      end if
   end if

end subroutine Check_if_photon_interact

! Compute the distance travel by the particle and modify its position and energy
!================================================================================
subroutine MovePhoton(photon, zfinal)
   use Particles, only: particle
   use Constantes, ONLY : c, H0, a0,prec
   implicit none
   type (particle), intent(inout) :: photon
   real(prec), intent(in) :: zfinal
   real(prec) :: deltaR
   real(prec) :: integrandR,integrandt,qgauss
   external integrandR, integrandt, qgauss

   deltaR=(c/(H0*a0))*qgauss(integrandR,zfinal,photon%redshift)! Compute the distance travel by the photon
   photon%pos(1)=(photon%dir(1)*deltaR+photon%pos(1))
   photon%pos(2)=(photon%dir(2)*deltaR+photon%pos(2))
   photon%pos(3)=(photon%dir(3)*deltaR+photon%pos(3))
   photon%energy=photon%energy*(1+zfinal)/(1+photon%redshift)
   photon%redshift = zfinal

end subroutine MovePhoton

! compute the pair production 
!=============================
subroutine make_pairs(photon_gamma, electron, positron)
   use Particles, only: particle, duetoCMB
   use Constantes, only: prec, NbPairsProd, NbLeptonsProd, currentEnergy, Tcmbp
   implicit none
   type (particle), intent(inout) :: photon_gamma, electron, positron
   real (prec) :: Emin, hnu1, hnu2, cosTheta, ktr
   real (prec) :: gama, beta, betaSquare
   real (prec) :: Xsection, randomXsection, ran0
   real (prec) :: dirPhotonLowEnergy(3)
   real (prec) :: uu(2), vv(2), ww(2)
   external :: ran0
   logical :: rejection
   integer :: it
   cosTheta=0 
   rejection = .TRUE.
   Xsection=0
   randomXsection=1

   electron%charge = -1
   electron%pos = photon_gamma%pos
   electron%energy = 0
   electron%redshift = photon_gamma%redshift
   positron%charge = 1
   positron%pos = photon_gamma%pos
   positron%energy = 0
   positron%redshift = photon_gamma%redshift

   !$OMP ATOMIC
   nbPairsprod= nbpairsprod +1
   ! Compute the threshold of reaction if angle between photons is pi
   hnu1=photon_gamma%energy 

   do while (rejection)
      Emin = 1 / hnu1
      hnu2=0
      ! Select target photon: must have an energy higher than the threshold
      ! previously computed
      do while (hnu2 < Emin)
         ! Check if the interaction come from the CMB or the EBL
         if (duetoCMB) then
            emin = 1/currentEnergy
            ktr=Tcmbp*(1+photon_gamma%redshift)
            !    e,emin,kt in mec^2 units   
            if (emin/ktr.lt.1) then
               call selecPhotCMB(hnu2,photon_gamma%redshift) !redshift interaction
            else
               call genx(hnu2,emin,photon_gamma%redshift)
            end if 
         else
            call selecPhotEBL(hnu2,photon_gamma%redshift)
         end if
         ! select an angle accordig to the isotropic distribution
         costheta = 1
         do while (costheta==1)
            call isotrop(dirPhotonLowEnergy(1), dirPhotonLowEnergy(2),dirPhotonLowEnergy(3))
            ! cosTheta = cos (angle between photons)
            cosTheta = photon_gamma%dir(1)*dirPhotonLowEnergy(1) + & 
                       photon_gamma%dir(2)*dirPhotonLowEnergy(2) + &
                       photon_gamma%dir(3)*dirPhotonLowEnergy(3) 
         end do 
         Emin = 2 / photon_gamma%energy / (1-cosTheta)
         ! threshold which take into account the angle between photons
         it = it+1
      end do

      ! Lorentz factor of the created leptons in their CM frame
      gama = hnu1 * hnu2 * (1-cosTheta)/2 
      if (gama < 1) then
         write(6,*) "ERREUR GAMMA^2 < 1"
         write(6,*) photon_gamma%dir
         write(6,*) photon_gamma%energy
         STOP
      end if
      ! here we have a candidate target photon drawn from isotropic distribution
      ! we have to decide wether we keep it or not
      ! we do that with a rejection test on the amplitude of the cross section
      gama = sqrt(gama)
      betaSquare = 1-gama**(-2)
      beta = sqrt(betaSquare)
      Xsection = (3-betaSquare**2)*2*log((1+beta)*gama)
      Xsection = Xsection-2*beta*(2-betaSquare)
      Xsection = Xsection/((1+beta)*gama**2)*(1-cosTheta)/(1.36341d0*2) 
      ! 1,36341 = max of Xsession => 0 < Xsession < 1

      randomXsection = ran0()
      if (randomXsection < Xsection) then
         uu(1) = photon_gamma%dir(1)
         vv(1) = photon_gamma%dir(2)
         ww(1) = photon_gamma%dir(3)
         uu(2) = dirPhotonLowEnergy(1)
         vv(2) = dirPhotonLowEnergy(2)
         ww(2) = dirPhotonLowEnergy(3)
         !write(*,*), hnu1, hnu2, uu, vv, ww, cosTheta, gama
         call pair_production(hnu1, hnu2, uu, vv, ww, cosTheta, gama, betaSquare, beta)
         !write(*,*), hnu1,hnu2,hnu1+hnu2

#ifdef approx_Ee
         hnu1 = photon_gamma%energy/2
         hnu2 = photon_gamma%energy/2
#endif

         electron%dir(1) = uu(1)
         electron%dir(2) = vv(1)
         electron%dir(3) = ww(1)
         electron%energy = hnu1-1 ! total energy -> -1 kinematic energy
         electron%weight = photon_gamma%weight
         electron%time   = photon_gamma%time
         electron%generation = photon_gamma%generation+1
         positron%dir(1) = uu(2)
         positron%dir(2) = vv(2)
         positron%dir(3) = ww(2)
         positron%energy = hnu2-1
         positron%weight = photon_gamma%weight
         positron%time   = photon_gamma%time
         positron%generation = photon_gamma%generation+1

         rejection = .FALSE.
         
      end if
   end do

   nbLeptonsProd = nbLeptonsProd + positron%weight + electron%weight

end subroutine make_pairs

subroutine pair_production(hnu1,hnu2,uu,vv,ww,up1up2,gama,betaca,beta)
   use Constantes, only : prec, pi
   implicit none
   
   real(prec) :: hnu1,hnu2,up1,vp1,vp2,wp1,wp2,gama1
   real(prec) :: gama,beta,gambeta,betacm,gamacm
   real(prec) :: uu(2),vv(2),ww(2)
   real(prec) :: gambetacm,a,ucmup1,up1up2,up2,a1,a2,a3,b1,b2,b3,c1,c2,c3
   real(prec) :: al,ymax,y1,u,yt,x,y,zz,ran0,ap,am,phi,a11
   real(prec) :: mu,ucmu1cm,ucm,vcm,wcm,u1cm,v1cm,w1cm,shnu
   real(prec) :: y2,betaca,gama2,Etmp
   real(prec) :: u1,v1,w1,u2,v2,w2

   if (hnu1 > hnu2) then
      up1=uu(2)
      up2=uu(1)
      vp1=vv(2)
      vp2=vv(1)
      wp1=ww(2)
      wp2=ww(1)
      Etmp=hnu2
      hnu2=hnu1
      hnu1=Etmp
   else
      up1=uu(1)
      up2=uu(2)
      vp1=vv(1)
      vp2=vv(2)
      wp1=ww(1)
      wp2=ww(2)
   end if

   a=0
   ucm=0
   vcm=0
   wcm=0
   
   shnu=hnu1+hnu2
   
   gamacm=.5d0*shnu/gama
   betacm=sqrt(1-gamacm**(-2))
   gambetacm=gamacm*betacm
   gambeta=gama*beta
   
   y1=1+2*betaca*(1-betaca)
   y2=2*(betaca-1)**2-1
   yt=1
   y=0
   ymax=2
   
   do while (yt.gt.y)
      x=2*ran0()-1
      ap=(1+beta)**x
      am=(1-beta)**x

      x=(ap-am)/(ap+am)
      u=x**2
      y=y1-(y2*u+u**2)/(1-u)
      yt=ymax*ran0()
   enddo
   
   mu=x/beta
   phi=2*pi*ran0()
   y1=(1-mu**2)**.5d0
   y=cos(phi)*y1
   zz=sin(phi)*y1

   if (betacm.gt.1d-10) then

      a=shnu*betacm
      a1=hnu1/a
      a2=hnu2/a
      
      ucm=a1*up1+a2*up2
      vcm=a1*vp1+a2*vp2
      wcm=a1*wp1+a2*wp2
      
      ucmup1=a1+a2*up1up2

      a11=hnu1*(-gambetacm+(gamacm-1)*ucmup1)/gama
      a2=hnu1/gama
      
   else 
      
      a11=0
      a2=1

   endif
   
   a1=a11*ucm+a2*up1
   b1=a11*vcm+a2*vp1
   c1=a11*wcm+a2*wp1
     
   if (abs(a1).lt.7.d-1) then
      
      a2=(1-a1**2)**.5d0
      al=-a1/a2 
      b2=al*b1
      c2=al*c1
      
   else
      
      b2=sqrt(1-b1**2)
      al=-b1/b2
      a2=al*a1
      c2=al*c1
      
   endif
   
   a3=b1*c2-c1*b2
   b3=c1*a2-c2*a1
   c3=a1*b2-b1*a2
   
   u1cm=mu*a1+y*a2+zz*a3
   v1cm=mu*b1+y*b2+zz*b3
   w1cm=mu*c1+y*c2+zz*c3
   
   ucmu1cm=ucm*u1cm+vcm*v1cm+wcm*w1cm
   
   y=shnu/2
   zz=a*beta*ucmu1cm/2
   
   gama1=y+zz
   gama2=y-zz
            
   a1=gambetacm*gama
   b1=gambeta*(gamacm-1)*ucmu1cm
   c1=gambeta
   a3=sqrt(gama1**2-1)
   b3=sqrt(gama2**2-1)
   if (a3.gt.0.and.b3.gt.0) then
      u1=((a1+b1)*ucm+c1*u1cm)/a3       
      v1=((a1+b1)*vcm+c1*v1cm)/a3
      w1=((a1+b1)*wcm+c1*w1cm)/a3
      
      if (abs(sqrt(up2**2+vp2**2+wp2**2)-1).gt.1d-5) then
         
         write(6,*) 'erreur pair prod  2'
         
         write(6,*) u1
         write(6,*) v1
         write(6,*) w1
         write(6,*) sqrt(u1**2+v1**2+w1**2)
      endif
      
      u2=((a1-b1)*ucm-c1*u1cm)/b3       
      v2=((a1-b1)*vcm-c1*v1cm)/b3
      w2=((a1-b1)*wcm-c1*w1cm)/b3
      
      
      if (abs(sqrt(up2**2+vp2**2+wp2**2)-1).gt.1d-5) then
         
         write(6,*) 'erreur anil  2'
         write(6,*) u2
         write(6,*) v2
         write(6,*) w2
         write(6,*) sqrt(u2**2+v2**2+w2**2)
      endif
   else
      write(6,*) 'Erreur pair_production'
      ! no momentum conservation !!!!!!!
      call isotrop(u1,v1,w1)
      call isotrop(u2,v2,w2)
   endif
   
   !if (gama1.ge.gama2) then
      hnu1=gama1          
      hnu2=gama2
      uu(1)=u1
      vv(1)=v1
      ww(1)=w1
      uu(2)=u2
      vv(2)=v2
      ww(2)=w2
   !else
   !   hnu1=gama2
   !   hnu2=gama1
   !   uu(1)=u2
   !   vv(1)=v2
   !   ww(1)=w2
   !   uu(2)=u1
   !   vv(2)=v1
   !   ww(2)=w1
   !endif

end subroutine pair_production

!===============================================================================         
                        END MODULE Photon_Routines
!===============================================================================         
