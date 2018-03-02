#include "../temp/preprocessing.f95"
!===============================================================================         
                              MODULE Particles
!===============================================================================         
use constantes, only: prec
implicit none

type particle
   integer :: charge    ! -1 = electron, 0 = part, 1 = positron
   real(prec) :: time
   real(prec) :: pos(3) ! position -> cartesian (x, y, z)
   real(prec) :: energy
   real(prec) :: dir(3) ! direction -> cartesian (u, v, w)
   real(prec) :: redshift, weight
   integer :: generation
end type particle

! source
integer :: Nature_source
real(prec) :: emission_redshift, Emin_source, Emax_source, Esource, dir_source
real(prec) :: zEarth=0.

real(prec) :: charge
logical :: duetoCMB, willInteract
integer, parameter :: stacksizeini=10

integer :: nbLeptonsSaved=0, nbPhotonsSaved=0

!$OMP THREADPRIVATE(duetoCMB, willInteract,charge,Esource,dir_source)

!===============================================================================         
                                   CONTAINS
!===============================================================================         
subroutine startsource()
   use constantes, only: c, h0, mpc, distSourceToEarth
   implicit none
   real (prec) :: qgauss, integrandR
   external qgauss, integrandR
   ! Compute the distance from the source to the Earth (Global variable)
   distSourceToEarth = qgauss(integrandR, zEarth, emission_redshift) 
   distSourceToEarth = (c/H0)*distSourceToEarth
   write(*,*) "Source ==========================================================="
   write(*,'(a,es10.2)')   "  > redshift          ", emission_redshift
   write(*,'(a,f10.2,a)')  "  > distance          ", distSourceToEarth/Mpc, " Mpc"
   write(*,'(a,es10.1,a)') "  > Energy min        ", Emin_source, " TeV"
   write(*,'(a,es10.1,a)') "  > Energy max        ", Emax_source, " TeV"
   Emin_source=Emin_source*1.d12 /(510998.918) ! cgs unit
   Emax_source=Emax_source*1.d12 /(510998.918) ! cgs unit
end subroutine startsource

subroutine initiate_source_particle(part)
   use constantes, only: c, h0, mpc
   implicit none
   type (particle), intent(out) :: part
   real(kind=prec) :: logE, ran0, integrandt, qgauss
   external ran0, integrandt, qgauss
   
   ! uniform distribution (in log) between Emin and Emax
   logE = log10(Emin_source) + ran0() * log10(Emax_source/Emin_source)
   Esource = 10**logE
   part%energy = Esource

   part%charge = Nature_source
   part%redshift = emission_redshift ! Redshift emission
   part%pos=(/0,0,0/)
   part%dir=(/0,0,1/)
   part%weight=1.
   part%generation=0

   part%time = -qgauss(integrandt,emission_redshift,zearth)

end subroutine initiate_source_particle

subroutine ComputeTime(part, zfinal)
   use Constantes, only : prec, currentRedshift
   implicit none
   type (particle), intent(inout) :: part
   real(prec), intent(in) :: zfinal
   real(prec) :: tauf, integrandt, qgauss
   external integrandt, qgauss

   tauf = qgauss(integrandt,currentRedshift,zfinal)
   part%time = part%time + tauf

end subroutine ComputeTime

! add a particle to stack                                              
!=========================
subroutine store_particle(part,partstack,StackSize)
   use Constantes, only : prec, currentEnergy, Ethreshold, Eth_lept, alphaS, positronsLoosedID
   implicit none
   type (particle), dimension(:), allocatable, intent(in) :: part
   type (particle), allocatable, dimension(:) :: partstack
   type (particle), allocatable, dimension(:) :: tempstack
   real (prec) :: frac, Emin, ran0
   logical :: found
   integer :: i, j, n, StackSize
   external ran0

   stackSize = stackSize-1
   do n=1, size(part)
      ! store only particle with a sufficiently high energy -> speed it up
      frac=(part(n)%energy/currentEnergy)**alphaS
      if ((part(n)%charge == 0))  then
         Emin = Ethreshold
      else
         Emin = Eth_lept
      end if

      if ((ran0()<frac) .and. (part(n)%energy > Emin)) then

         if (stacksize == size(partstack)) then
            allocate(tempstack(stacksize))
            tempstack(:) = partstack(:)
            deallocate(partstack)
            allocate(partstack(stacksize+stacksizeini))
            partstack(1:stacksize)=tempstack(1:stacksize)
            deallocate(tempstack)
            print*, " Increasing stack size from ",stacksize,"to",size(partstack)
         end if

         ! particle sort by decreasing energy, sampling weight for produced particle
         found = .false.
         StackSize=StackSize+1  ! enlarge stack
         if (Stacksize==1) then ! 1st particle in stack
            partStack(1)= part(n)
            partStack(1)%weight=partStack(1)%weight/frac
         else ! StackSize > 1: add particle to E-ordered stack
            i=1
            do while ( (i<StackSize) .and. (.not. found) ) 
               j=StackSize-i
               if(partStack(j)%energy > part(n)%energy) then 
                  partStack(j+1)= part(n)                    ! add particle to stack
                  partStack(j+1)%weight=partStack(j+1)%weight/frac
                  found = .true.
               else
                  partStack(j+1)=partStack(j)                ! re-arrange stack
               endif
               i=i+1
            end do
         end if
#ifdef file_positrons
      else if (part(n)%charge == 1) then
      !$OMP CRITICAL (positrons_loosed)
         write(unit=positronsLoosedID,fmt=*) part(n)%generation, part(n)%weight, & 
            part(n)%redshift, part(n)%energy*510.998918*1d-6, Esource*510.998918*1d-6
      !$OMP end CRITICAL (positrons_loosed)
#endif
      end if
   end do
end subroutine store_particle

! Register particles detected
!=============================
subroutine register_detected_particle(part,stacksize)
   use Constantes, only : prec, ResultID, H0, pi 
   type (particle) :: part
   real(prec) :: theta_p, cos_theta_d, sin_phi_d, cos_phi_d
   real(prec) :: ud(3), vd(3)
   integer :: stacksize

   stackSize = stackSize-1

   ! normalize arrival position vector
   part%pos =part%pos/sqrt(sum(part%pos**2)) 
   ! compute arrival position (spherical coord. theta)
   theta_p = acos(part%pos(3))

   ! compute arrival direction (spherical coord. theta: cos)
   cos_theta_d = dot_product(part%pos,part%dir)
   if (cos_theta_d > 1) then
      cos_theta_d = 1
   else if (cos_theta_d < -1) then
      cos_theta_d = -1
   end if

   !$OMP CRITICAL (register)
   ! compute arrival position (spherical coord. phi: sin and cos)
   if (theta_p == 0) then
      write(unit=ResultID,fmt=*) part%generation, part%weight, &
         part%energy*510.998918*1d-6, part%time/H0, theta_p, acos(cos_theta_d), 0, &
         Esource*510.998918*1d-6, part%charge
   else
      vd = (/-part%pos(2),part%pos(1),0./)/sin(theta_p)
      cos_phi_d = dot_product(part%dir,vd)
      call pdtVect(vd,part%pos,ud)
      sin_phi_d = dot_product(part%dir,ud)
      cos_phi_d = cos_phi_d / sqrt(sin_phi_d**2+cos_phi_d**2)
      sin_phi_d = sin_phi_d / sqrt(sin_phi_d**2+cos_phi_d**2)

      ! Energy [GeV], time [s], position, direction and arrival angle [rad]
      write(unit=ResultID,fmt=*) part%generation, part%weight, &
         part%energy*510.998918*1d-6, part%time/H0, theta_p, acos(cos_theta_d), &
         atan2(sin_phi_d,cos_phi_d), Esource*510.998918*1d-6, part%charge
   end if
  
   if (charge == 0) then
      nbPhotonsSaved=nbPhotonsSaved+1
   else
      nbLeptonsSaved=nbLeptonsSaved+1 
   end if
   !$OMP end CRITICAL (register)

end subroutine register_detected_particle

subroutine spherical_coord(pos)
   use Constantes, only : pi
   real(prec), intent(inout) :: pos(3)
   real(prec) :: hyp, theta, phi

   hyp = sqrt(sum(pos**2))
   theta = acos(pos(3)/hyp)
   if (pos(1) == 0) then 
      phi = 0 ! could be whatever
   else
      phi = atan(pos(2)/pos(1))
   end if

   pos(1) = hyp
   pos(2) = theta
   pos(3) = phi

end subroutine spherical_coord

!===============================================================================         
                           END MODULE Particles
!===============================================================================         
