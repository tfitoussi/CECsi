#include "../preprocessing.f95"

module test
use constantes, only: prec
integer, parameter :: nh=100
real(prec) :: e,eh(nh+1),hist(nh),dens(nh),x(nh),dx(nh)
real(prec), parameter :: emin=1.e-4, emax=2.e1
end module test

!===============================================================================
                               PROGRAM CASCADE
!===============================================================================
use Constantes, only : prec, OMP_num_threads, iseed, currentEnergy, currentRedshift , Mpc
use EGMF, only : intEGMF, lambdaEGMF, Nm, selectEGMF, initiateEGMF
#if defined B_turbulent
use EGMF, only : EGMFtab, init_B
#elif defined B_pseudo_random
use EGMF, only : EGMFtab, init_B
#elif defined B_random
use EGMF, only : EGMFtab
#endif
use test, only: nh,e,eh,hist,dens,x,dx,emin,emax
implicit none
real(kind=4) :: tbegin
integer, dimension(3) :: current_time
integer :: i, nbphot=0, nmax=1
real(prec) :: b(3), posini(3), ran0, zi, zmax, zmin
external ran0

call initiate()
! JE DEFINIS LA GRILLE POUR MON HISTOGRAMME PERSO
print*,'  lecture EBL'
eh = emin*(emax/emin)**(/((i-1.)/nh,i=1,nh+1)/)
dens = 0
do i=1,nh
   x(i) = 0.5*(eh(i)+eh(i+1))
   dx(i) = eh(i+1)-eh(i)
   call NeInterpol(x(i)/511.e3,0.01,dens(i))
enddo
dens = dens/maxval(dens)
open(unit=323,file='output/EBL_CMB_test.dat',action='write')

print*,'  Main loop...'
!$OMP PARALLEL do SCHEDULE(DYNAMIC) NUM_THreadS(OMP_num_threads)
do i=1, nmax
   iseed = -i*3591857 ! initiate the seed generator (reproductibility)
   call test_CMB_EBL_in_IC()
   if (modulo(nbphot*100, nmax) == 0) then
      call print_time()
   end if
end do
!$OMP end PARALLEL do

! JE NORMALISE MON HISTOGRAM PERSO ET J'ECRIS TOUT DANS UN FICHIER
hist = hist/dx
hist = hist/maxval(hist)
print*,'  Ecriture fichier'
open(unit=47,file='output/EBL_CMB_test2.dat',action='write')
do i=1,nh
   write(47,*) x(i),hist(i),dens(i)
enddo
close(47)
close(323)

call finish()

!===============================================================================
                                   CONTAINS
!===============================================================================
subroutine test_CMB_EBL_in_IC()
   use constantes, only : currentEnergy, currentRedshift, H0,c,Mpc
   use particles, only: charge, particle, duetoCMB, initiate_source_particle
   use Lepton_Routines, only : check_if_lepton_interact
   use test, only: nh,e,eh,hist,dens,x,dx,emin,emax
   implicit none
   type (particle) :: part
   real (prec) :: Zfinal, posini(3), dirini(3), cos_theta, e2
   integer :: n, stackSize, k

   call initiate_source_particle(part)
   currentRedshift = part%redshift
   currentEnergy= part%energy
   charge = float(part%charge)

   if (part%charge==0) then ! it's a photon
   else ! it's a lepton
      call Check_if_lepton_interact(part, Zfinal)

      if (duetoCMB) then
         call selecPhotCMB(e2,currentRedshift)
      else
         call selecPhotEBL(e2,currentRedshift)
      end if
      !$OMP CRITICAL (test)
      if ( .not. duetoCMB ) then
         k = 1+floor(nh*log(e2*510998.918/emin)/log(emax/emin))
         hist(k) = hist(k) + 1
      end if
      ! JE STOCKE AUSSI DIRECTEMENT la VALEUR DE L'ENERGIE TIREE
      write(323,*) currentRedshift, e2*510998.918
      !$OMP END CRITICAL (test)

   end if
   !$OMP CRITICAL (nbphotons)
   nbphot=nbphot+1
   !$OMP end CRITICAL (nbphotons)
end subroutine test_CMB_EBL_in_IC
!-------------------------------------------------------------------------------

subroutine print_time()
   implicit none
   real(kind=4) :: tarray(2), tend, t_in_sec, nb_sec
   integer :: nb_hours, nb_mins
   character(len=*), parameter :: erase = repeat(achar(8),80)
   call itime(current_time)
   call ETIME(tarray, tend)
   t_in_sec = tend-tbegin
   nb_hours = int(t_in_sec/3600)
   t_in_sec = t_in_sec - nb_hours * 3600
   nb_mins  = int(t_in_sec/60)
   nb_sec   = t_in_sec - nb_mins * 60
   write(*,'(a,a,i2.2,a,i2.2,a,i2.2,a,i3,a,i4,a,i2,a,f5.2,a,f10.1,a)',advance='no') &
      erase,"  [",current_time(1),":",current_time(2),":",current_time(3),"] ", &
      int(float(nbphot)/float(nmax)*100)," % done using ", &
      nb_hours, "h ", nb_mins, "min ", nb_sec, "s CPU time"
   !write(*,'(a,a,i2.2,a,i2.2,a,i2.2,a,i3,a,i4,a,i2,a,f5.2,a,f10.1,a)',advance='no') &
   !   erase,"  [",current_time(1),":",current_time(2),":",current_time(3),"] ", &
   !   int(float(i)/float(nmax)*100)," % done using ", &
   !   nb_hours, "h ", nb_mins, "min ", nb_sec, "s CPU time"
end subroutine print_time
!-------------------------------------------------------------------------------

subroutine initiate()
   use constantes, only: Tcmbp, errorId, ResultID, deltaID, alphaS, ComptonEnergyThreshold, &
      ComptonThreshold, Ethreshold, ngp, xabsc, weig
   use particles, only: Emin_source, Emax_source, emission_redshift, Nature_source, startsource
   implicit none
   real(kind=4) :: tarray(2)

   ! Read the namelist parameters
   !=============================
   namelist /source/ emission_redshift,Emin_source,Emax_source, Nature_source, nmax
   namelist /EGMF_param/ IntEGMF, lambdaEGMF, Nm
   namelist /simulation/ OMP_num_threads, alphaS, Ethreshold, ComptonThreshold
   open(unit=10, file="input_parameters.f95")
   read(10, nml=source)
   read(10, nml=EGMF_param)
   read(10, nml=simulation)
   close(10)
   !emission_redshift=4 ! sinon problème de calcul des distances d'interaction (a
   !                    ! atteint la terre!)

   call ETIME(tarray, tbegin)
   write(*,*) "=================================================================="
   write(*,*) "|                       TESTS EBL CMB                            |"
   write(*,*) "=================================================================="
   call gauleg(ngp, xabsc, weig) ! Gauss-Legendre, weight and position
   call startsource()
   Call initiateEGMF()
   call ReadEBLfiles() ! load EBL model from files

   write(*,*) "Approximations ==================================================="
#ifdef approx_motionless_leptons
   write(*,*) " > Motionless leptons"
#endif
#ifdef approx_Eic
   write(*,*) " > Eic = 3.24 (Ee/1TeV)^2 GeV"
#endif
#ifdef approx_Thomson_regime
   write(*,*) " > Inverse Compton scattering always in Thomson regime"
#endif
#ifdef approx_lambda_ic
   write(*,*) " > lambda_ic such as tau_ic = 1"
#endif
#ifdef approx_Ee
   write(*,*) " > Ee = Egamma/2"
#endif
#ifdef approx_lambda_gg
   write(*,*) " > lambda_gg such as tau_gg = 1"
#endif

   write(*,*) "Acceleration ====================================================="
   write(*,'(a,f6.0,a)') "  > Energy threshold  ", Ethreshold," MeV"
   Ethreshold=Ethreshold*1.d6/(510998.918) ! adimensionnal
   ! Compute the Compton accumulation energy threshold (only on the CMB)
   ComptonEnergyThreshold = ComptonThreshold/((4./3.)*2.7*Tcmbp*(1+emission_redshift))
   write(*,'(a,f6.2,a,f4.2,a)') "  > Compton threshold ", Comptonthreshold*1d2,' % (',ComptonEnergyThreshold*510998.918d-12,' TeV)'
   write(*,'(a,f6.2,a)') "  > alphaS            ", alphaS
   !$ write(*,'(a,i4,a)') "  > On ", OMP_num_threads, " threads"

   ! File which contains all the particles "detected"
   open (unit=ResultID,file="output/results.dat")
   open (unit=errorID,file="output/Compton_Errors.dat")
   open (unit=deltaID,file="output/lepton_deflection.dat")

   call itime(current_time)
   write(*,*) "CPU time ========================================================="
   if (Nature_source == 0) then
      write(*,'(a,i2.2,a,i2.2,a,i2.2,a,i8,a)') "  [",current_time(1),":",current_time(2), &
         ":",current_time(3),"] Start running over ", nmax, " initial photons"
   else
      write(*,'(a,i2.2,a,i2.2,a,i2.2,a,i8,a)') "  [",current_time(1),":",current_time(2), &
         ":",current_time(3),"] Start running over ", nmax, " initial leptons"
   end if
end subroutine initiate
!-------------------------------------------------------------------------------

subroutine finish()
   use constantes, only: distSourceToEarth, Mpc, nbCompton, nbPairsProd, nbLeptonsProd, &
      errorId, ResultID, deltaID, densityEBL, redshiftEBL, enerEBL
   use particles, only: Emax_source, emission_redshift
   implicit none
   ! allocate in Routines/ReadEBL.f95: ReadEBLfiles
   deallocate(densityEBL)
   deallocate(redshiftEBL)
   deallocate(enerEBL)
#ifndef B_uniform
   deallocate(EGMFtab)
#endif

   close(ResultID)
   close(errorID)
   close(deltaID)

   open (unit=111,file="output/profile.dat")
   write(unit=111,fmt=*) IntEGMF, Emax_source*(510998.918d-12), distSourceToEarth/Mpc, &
      emission_redshift, nmax, nbLeptonsProd
   close(111)

   write(*,*) ""
   write(*,*) "=================================================================="
   write(*,'(a,i8)') "  > Inverse Compton scattering computed:",nbCompton
   write(*,'(a,i8)') "  > Pairs production computed:          ",nbPairsprod
   write(*,'(a,i8)') "  > Number of leptons generated:        ",nbLeptonsProd
   write(*,*) "=================================================================="
   write(*,*) "|                       End of Simulation                        |"
   write(*,*) "=================================================================="
end subroutine finish
!===============================================================================
                              END PROGRAM CASCADE
!===============================================================================         
