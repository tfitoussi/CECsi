#include "../temp/preprocessing.f95"
!===============================================================================         
                               PROGRAM CASCADE
!===============================================================================         
use Constantes, only : prec, OMP_num_threads, iseed, seed
#ifdef file_positrons
   use constantes, only: positronsID, positronsLoosedID
#endif
use EGMF, only : intEGMF, lambdaEGMF, Nm, selectEGMF, initiateEGMF
#if defined B_turbulent
use EGMF, only : EGMFtab, init_B
#elif defined B_pseudo_random
use EGMF, only : EGMFtab, init_B
#elif defined B_random
use EGMF, only : EGMFtab
#endif
implicit none
real(kind=4) :: tbegin
integer, dimension(3) :: current_time
integer :: i, EBLevts=0, nbphot=0, nmax=1

call initiate()
! Loop over the "nmax" photons from the source
!=============================================
!$OMP PARALLEL do SCHEDULE(DYNAMIC) NUM_THreadS(OMP_num_threads) 
do i=1, nmax
   iseed = -i*seed ! initiate the seed generator (reproductibility)
   call make_cascade()
   if (modulo(nbphot*100, nmax) == 0) then
      call print_time()
   end if
end do
!$OMP end PARALLEL do
call finish()

!===============================================================================         
                                   CONTAINS
!===============================================================================         
subroutine make_cascade()
   use constantes, only : currentEnergy, currentRedshift, deltaID,H0,c,Mpc
   use particles, only: charge, willInteract,  stacksizeini, computetime, store_particle, &
      register_detected_particle, initiate_source_particle, particle, duetocmb, Esource
   use Photon_Routines, only : movePhoton, check_if_photon_interact, make_pairs
   use Lepton_Routines, only : moveLepton, check_if_lepton_interact, Inverse_Compton_scattering
   implicit none
   type (particle), dimension(:), allocatable :: tempPart, partstack
   real (prec) :: Zfinal, posini(3), posfinal(3), EGMF(3)
   integer :: n, stackSize
   real(prec) :: integrandt, qgauss
   external      integrandt, qgauss
   ! Initialize the stack of particles
   !===================================
   !  - every particle produce during the cascade will store in a stack
   !  - stack sort by increasing energy (last particle -> lower energy)
   !  - loop over the stack until the stack is empty
   stackSize = 1
   allocate(partStack(stacksizeini))    
   allocate(tempPart(2)) ! particles producted by IC and Pair prod
   ! Initialse first gamma-ray photon characteristics
   call initiate_source_particle(partStack(stackSize))

   do while (stackSize.ne.0)
      ! Take the last element of the stack (lower energy particle) and start study
      ! store the current redshift and energy (needed in some integrals)
      currentRedshift = partStack(stackSize)%redshift
      currentEnergy= partStack(stackSize)%energy
      charge = float(partStack(stackSize)%charge)
#ifdef file_cascade_traj
      posini = partStack(stackSize)%pos
#endif

      if (partStack(stackSize)%charge==0) then ! it's a photon
         call Check_if_photon_interact(partStack(stackSize), Zfinal)
         call MovePhoton(partStack(stackSize), Zfinal)
         call ComputeTime(partStack(stackSize), Zfinal)
         !if (partStack(stackSize)%generation==0) then
         !   print*, c/(H0*Mpc)* qgauss(integrandt,currentRedshift,zfinal)
         !   stop
         !end if
#ifdef file_cascade_traj
         posfinal = partStack(stackSize)%pos
#endif
      
         if ( willInteract ) then ! the photon can interact before reaching the detector
            ! compute the pair production and store the two leptons produce in the
            ! stack
            call make_pairs(partStack(stackSize),tempPart(1),tempPart(2))
#ifdef file_positrons
            !$OMP CRITICAL (positrons)
            if (tempPart(1)%charge == 1) then
               write(unit=positronsID,fmt=*) tempPart(1)%generation, tempPart(1)%weight, & 
                  tempPart(1)%redshift, tempPart(1)%energy*510.998918*1d-6, &
                  Esource*510.998918*1d-6
            else
               write(unit=positronsID,fmt=*) tempPart(2)%generation, tempPart(2)%weight, & 
                  tempPart(2)%redshift, tempPart(2)%energy*510.998918*1d-6, &
                  Esource*510.998918*1d-6
            end if
            !$OMP end CRITICAL (positrons)
#endif
            call store_particle(tempPart,partStack,stackSize)

         else ! the photon interact further than the detector so it is detected before
            ! write the final file "Result" which contains the properties of the
            ! photon (charge, Zlim, energy, position and direction)
            call register_detected_particle(partStack(stackSize),stacksize)
         end if
      
      else ! it's a lepton
         call Check_if_lepton_interact(partStack(stackSize), Zfinal)
#ifdef file_lepton_deflection
         !$OMP CRITICAL 
         call selectEGMF(partStack(stackSize)%pos,EGMF)
         write(unit=deltaID,fmt=*) partStack(stackSize)%energy*(510.998918d-6),&
            partStack(stackSize)%pos, partStack(stackSize)%dir, EGMF
         !$OMP end CRITICAL 
#endif

#ifdef approx_motionless_leptons
         posini = partStack(stackSize)%pos
#endif
         call MoveLepton(partStack(stackSize), Zfinal)
#ifdef approx_motionless_leptons
         ! change direction not the position and redshift ~ local and instantaneous cooling
         partStack(stackSize)%redshift = currentRedshift
         partStack(stackSize)%pos = posini 
#else
         call ComputeTime(partStack(stackSize), Zfinal)
#endif
#ifdef file_cascade_traj
         posfinal = partStack(stackSize)%pos
#endif

         if ( willInteract ) then ! the lepton can interact before reaching the detector
            ! compute the inverse Compton scattering and store the two leptons produce in the
            ! stack
            tempPart(1) = partStack(stackSize) ! store the initial lepton -> will be modify
            call Inverse_Compton_scattering(tempPart(1),tempPart(2))
            call store_particle(tempPart,partStack,stackSize)
            if (.not. duetocmb) then
               EBLevts = EBLevts +1
            end if
         
         else ! the lepton interact further than the detector so it is detected before
            ! write the final file "Result" which contains the properties of the
            ! lepton (charge, Zlim, energy, position and direction)
            call register_detected_particle(partStack(stackSize),stacksize)
         end if            
      end if 
#ifdef file_cascade_traj
      !$OMP CRITICAL 
      write(unit=deltaID,fmt=*) int(charge), posini/Mpc, posfinal/Mpc
      !$OMP end CRITICAL 
#endif
   end do
   deallocate(tempPart)
   deallocate(partStack)    
   !$OMP CRITICAL (nbphotons)
   nbphot=nbphot+1
   !$OMP end CRITICAL (nbphotons)
end subroutine make_cascade
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
end subroutine print_time
!-------------------------------------------------------------------------------

subroutine initiate()
   use constantes, only: Tcmbp, distSourceToEarth, errorId, ResultID, deltaID, &
      alphaS, ComptonEnergyThreshold, ComptonThreshold, Ethreshold, ngp, xabsc, &
      weig, Eth_lept 
   use particles, only: Emin_source, Emax_source, emission_redshift, Nature_source, startsource
   implicit none
   real(kind=4) :: tarray(2)
#ifdef ProbeEBL
   integer :: i
   real(prec) :: Eebl, z
#endif 

   ! Read the namelist parameters
   !=============================
   namelist /source/ emission_redshift,Emin_source,Emax_source, Nature_source, nmax
   namelist /EGMF_param/ IntEGMF, lambdaEGMF , Nm
   namelist /simulation/ OMP_num_threads, alphaS, Ethreshold, ComptonThreshold, seed
   open(unit=10, file="temp/input_parameters.f95")
   read(10, nml=source)
   read(10, nml=EGMF_param)
   read(10, nml=simulation)
   close(10)

   if (Emin_source < Emax_source ) then
      Ethreshold = Emin_source*1d3
   end if
   ! limit for a lepton to produce a gamma with an energy upper than
   ! Ethreshold analytic expression: Egamma = 3.24 (Ee/1TeV)^2 GeV 
   ! take 30% (security)
   !Eth_lept = 555.6*sqrt(Ethreshold*1d-3)
   Eth_lept = Ethreshold
      
   call ETIME(tarray, tbegin) 
   write(*,*) "=================================================================="
   write(*,*) "|        Cosmological Electromagnetic Cascade Simulation         |"
   write(*,*) "=================================================================="
   call gauleg(ngp, xabsc, weig) ! Gauss-Legendre, weight and position
   call startsource()
   call ReadEBLfiles() ! load EBL model from files
   Call initiateEGMF()

   write(*,*) "Approximations ==================================================="
#ifdef approx_motionless_leptons
   write(*,*) " > Motionless leptons"
#endif
#ifdef approx_Eic
   write(*,*) " > Eic = 3.24 (Ee/1TeV)^2 GeV"
#endif
#ifdef approx_no_EBL_ic
   write(*,*) " > No EBL photons target in inverse Compton scattering"
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
   write(*,'(a,es6.0,a,es8.2,a)') "  > Energy threshold  ",Ethreshold," GeV (lept: ",Eth_lept," GeV)"
   Ethreshold = Ethreshold*1.d9/(510998.918) ! adimensionnal
   Eth_lept = Eth_lept*1.d9/(510998.918) ! adimensionnal
   ! Compute the Compton accumulation energy threshold (only on the CMB)
   ComptonEnergyThreshold = ComptonThreshold/((4./3.)*2.7*Tcmbp*(1+emission_redshift))
   write(*,'(a,f6.2,a,f4.2,a)') "  > Compton threshold ", Comptonthreshold*1d2,' % (',ComptonEnergyThreshold*510998.918d-12,' TeV)'
   write(*,'(a,f6.2,a)') "  > alphaS            ", alphaS
   !$ write(*,'(a,i4,a)') "  > On ", OMP_num_threads, " threads"

   ! File which contains all the particles "detected"
#ifdef z_0_limit   
   open (unit=ResultID,file="temp/results_z_0.dat")     
#else
   open (unit=ResultID,file="temp/results.dat")     
#endif
   open (unit=errorID,file="temp/Compton_Errors.dat")     
#ifdef file_lepton_deflection
   open (unit=deltaID,file="temp/lepton_deflection.dat")     
#endif
#ifdef file_cascade_traj
   open (unit=deltaID,file="temp/cascade_traj.dat")     
#endif
#ifdef file_positrons
   open (unit=positronsID,file="temp/positrons.dat")     
   open (unit=positronsLoosedID,file="temp/positrons_loosed.dat")     
#endif

   call itime(current_time)
   write(*,*) "CPU time ========================================================="
   if (Nature_source == 0) then 
      write(*,'(a,i2.2,a,i2.2,a,i2.2,a,i9,a)')  "  [",current_time(1),":",current_time(2), &
         ":",current_time(3),"] Start running over ", nmax, " initial photons"
   else 
      write(*,'(a,i2.2,a,i2.2,a,i2.2,a,i9,a)')  "  [",current_time(1),":",current_time(2), &
         ":",current_time(3),"] Start running over ", nmax, " initial leptons"
   end if
end subroutine initiate
!-------------------------------------------------------------------------------

subroutine finish()
   use constantes, only: distSourceToEarth, Mpc, nbCompton, nbPairsProd, nbLeptonsProd, &
      errorId, ResultID, deltaID, densityEBL, redshiftEBL, enerEBL 
   use particles, only: Emax_source, emission_redshift, nbLeptonsSaved, nbPhotonsSaved
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
#ifdef file_positrons
   close(positronsID)
   close(positronsLoosedID)
#endif

   open (unit=111,file="temp/profile.dat")     
   write(unit=111,fmt=*) IntEGMF, Emax_source*(510998.918d-12), distSourceToEarth/Mpc, emission_redshift, &
      nmax, nbLeptonsProd
   close(111)

   call itime(current_time)

   write(*,*) ""
   write(*,*) "=================================================================="
   write(*,'(a,i8)') "  > Inverse Compton scattering computed:",nbCompton
   write(*,'(a,i8)') "  > Pairs production computed:          ",nbPairsprod
   write(*,'(a,i8)') "  > Number of leptons generated:        ",nbLeptonsProd
   write(*,'(a,i8)') "  > Number of lepton events saved:      ",nbLeptonsSaved
   write(*,'(a,i8)') "  > Number of photon events saved:      ",nbPhotonsSaved
   write(*,*) "=================================================================="
   write(*,*) "|                       End of Simulation                        |"
   write(*,*) "=================================================================="
end subroutine finish
!===============================================================================         
                              END PROGRAM CASCADE
!===============================================================================         
