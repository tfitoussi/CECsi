#include "../temp/preprocessing.f95"
! ============================================================================         
!                      Extragalactic magnetic field
                              MODULE EGMF 
! ============================================================================         
use constantes, only: DistSourceToEarth, Mpc, prec, pi
implicit none
!### Namelist parameters ###
real(prec) ::IntEGMF, lambdaEGMF ! magnetic field strength and correlation
real(prec) :: wLarmor
integer :: Nm ! turbulent magnetic field

integer :: nCube, sizeini  = 50000

#if defined(B_pseudo_random)
real(prec) :: sizeEGMF_r,ncube_r  ! local varialbes (must be real and not integers
                                ! because they can become very large)
real(prec) :: scalef
integer, parameter :: sizeEGMF = 1000000
real(prec), dimension(:), allocatable :: EGMFtab
#elif defined(B_random)
type EGMFdir
   real(prec) :: B(3)
   integer :: i, j, k
end type EGMFdir
integer :: sizeEGMF = 0
type (EGMFdir), dimension(:), allocatable :: EGMFtab
#elif defined(B_turbulent)
type EGMFturb
   real(prec) :: Xi_n(3)
   real(prec) :: A_kn, k_n, phi_n, theta_n, phase_n
end type EGMFturb
type (EGMFturb), dimension(:), allocatable :: EGMFtab
#endif

!==============================================================================         
                                   CONTAINS
!==============================================================================         
subroutine initiateEGMF()
   use constantes, only: e, m, c
   implicit none
   wLarmor = (IntEGMF*e)/(m*c)
   nCube = int(2*(DistSourceToEarth/Mpc)/lambdaEGMF)+1
   write(*,*) "EGMF ============================================================="
#if defined B_uniform
   write(*,*) " > Magnetic field uniform in the whole Universe"
#elif defined(B_turbulent)
   write(*,*) " > Turbulent magnetic field"
   write(*,*) " > Number of modes:", Nm
   allocate(EGMFtab(Nm))
   call init_b()
#elif defined(B_pseudo_random)
   write(*,*) " > Pseudo-random distribution of the magnetic field"
   allocate(EGMFtab(sizeEGMF))
   call init_b()
#elif defined(B_random)
   write(*,*) " > Fully random distribution of the magnetic field"
   allocate(EGMFtab(sizeini))
#endif
   write(*,'(a,es10.1,a)') "  > Strenght          ", intEGMF,    " Gauss"
   write(*,'(a,es10.1,a)') "  > Coherence length  ", lambdaEGMF, " Mpc"
   write(*,'(a,es10.1,a)') "  > Larmor pulsation  ", wLarmor,    " Hz"
end subroutine initiateEGMF

!===================== EGMF UNIFORM ===========================
#if defined(B_uniform)
subroutine selectEGMF(leptonpos,B)
   implicit none
   real(prec), intent(out) :: B(3)
   real(prec), dimension(3) , intent(in) :: leptonpos
   b=(/0,1,0/)
end subroutine selectEGMF

!===================== EGMF PSEUdo-RANDOM ===========================
#elif defined(B_pseudo_random)
subroutine init_b()
   ! ------------------------------------------------------------------------------
   ! This routine initialize the reference array once for the entire simulation
   ! with random values.
   ! ------------------------------------------------------------------------------
   !$ use constantes, only: OMP_num_threads
   implicit none
   integer :: i
   real(prec) :: ran0
   external ran0
   scalef = real(huge(1),kind=prec)/(real(nCube,kind=prec)*(real(nCube,kind=prec)+10))

   !$OMP PARALLEL do SCHEDULE(DYNAMIC) NUM_THREADS(OMP_num_threads)
   do i=1,sizeEGMF
       EGMFtab(i) = ran0()
   enddo
   !$OMP end PARALLEL do
end subroutine init_b

subroutine selectEGMF(leptonpos,B)
   implicit none
   real(prec), intent(out) :: B(3)
   real(prec), dimension(3) , intent(in) :: leptonpos

   real(prec) :: cos_t,sin_t,phi   ! local angle variables
   real(prec) :: ratio
   real(prec) :: irand
   integer :: n
   integer, dimension(3) :: IndEGMF

   ! find position
   do n=1, 3 ! x, y, z
      ratio = leptonpos(n)/DistSourceToEarth
      if (abs(ratio) <= 1) then
         IndEGMF(n)=int((ratio+1.)*(nCube/2))+1
      else if (ratio > 0) then
         IndEGMF(n) = nCube/2 +1
      else
         IndEGMF(n) = 1
      end if
   end do

   ! ------------------------------------------------------------------------------
   ! This routine computes the magnetic 3-vector in cell (i,j,k) by drawing values 
   ! from the reference array in a pseudo-random manner.
   ! The size sizeEGMF of the reference array must be much larger then the length nCube of 
   ! the cube in one direction: sizeEGMF>>nCube (better: sizeEGMF>>nCube**2) to minimize 
   ! the formation of non-random patterns.
   ! The result is perfectly random for sizeEGMF>=nCube**3.
   ! ------------------------------------------------------------------------------

   ! Compute angles
   phi = 0
   cos_t = 2 * ran1(IndEGMF(1),IndEGMF(2),IndEGMF(3)) - 1
   sin_t = sqrt(1-cos_t*cos_t)
   phi = 2*pi * ran2(IndEGMF(1),IndEGMF(2),IndEGMF(3))

   ! Compute components of the B-field
   b(1) = sin_t * cos(phi)
   b(2) = sin_t * sin(phi)
   b(3) = cos_t
end subroutine selectEGMF

function ran1(i,j,k) result(xrand)
   use constantes, only: pi,prec
   implicit none
   real(prec), parameter :: sqrt2 = sqrt(2.)
   integer, intent(in) :: i,j,k
   integer ::  irand,seed
   real(prec) :: xrand

   ! Get array index --------------------------------------
   seed = int((1.*i+j*real(k,kind=prec))*scalef)
   xrand = seededran(seed)

   seed = int((2.*j+i*real(k,kind=prec))*scalef)
   xrand = xrand + pi * seededran(seed)

   seed = int((3.*k+i*real(j,kind=prec))*scalef)
   xrand = xrand + sqrt2 * seededran(seed)

   irand = int(xrand/10.*huge(1))
   irand = modulo(irand,sizeEGMF)

   ! Retrieve random value --------------------------------
   xrand = EGMFtab(irand)

end function ran1

function ran2(i,j,k) result(xrand)
   use constantes, only: pi,prec
   implicit none
   real(prec), parameter :: sqrt2=sqrt(2.),sqrt5=sqrt(5.)
   integer, intent(in) :: i,j,k
   integer ::  irand,seed
   real(prec) :: xrand

   ! Get array index --------------------------------------
   seed = int((7.*i+j*real(k,kind=prec))*scalef)
   xrand = sqrt2 * seededran(seed)

   seed = int((5.*j+i*real(k,kind=prec))*scalef)
   xrand = xrand + sqrt5 * seededran(seed)

   seed = int((3.*k+i*real(j,kind=prec))*scalef)
   xrand = xrand + pi * seededran(seed)

   irand = int(xrand/50.*huge(1))
   irand = modulo(irand,sizeEGMF)

   ! Retrieve random value --------------------------------
   xrand = EGMFtab(irand)
end function ran2

function seededran(iseed) result(x)
! Custom random generator.
! It generates a random number with a uniform distribution between 0 and 1
! 
! INPUT:
!       none
! OUTPUT:
!       ran0: the generated value for the random varaible
   implicit none
   real(prec) :: x
   integer,parameter :: IA=16807,IM=2147483647,IQ=127773,IR=2836
   real, parameter :: am=nearest(1.e0,-1.e0)/IM
   integer :: ix, iy, k
   integer, intent(in) :: iseed

   ix = ieor(777755555,iseed)
   ix = ieor(ix,ishft(ix,13))
   ix = ieor(ix,ishft(ix,-17))
   ix = ieor(ix,ishft(ix,5))

   iy = ior(ieor(888889999,iseed),1)
   k = iy/IQ
   iy = IA*(iy-k*IQ)-IR*k
   if (iy<0) iy = iy+IM

   x = am*ior(iand(IM,ieor(ix,iy)),1)

end function seededran

!===================== EGMF FULLY RANdoM ===========================
#elif defined(B_random)
subroutine selectEGMF(leptonpos,B)
   implicit none
   real(prec), intent(out) :: B(3)
   real(prec), dimension(3), intent(in) :: leptonpos
   logical :: find_or_insert

   real(prec) :: ratio
   integer :: n, i
   integer, dimension(3) :: IndEGMF

   ! find position
   do n=1, 3 ! x, y, z
      ratio = leptonpos(n)/DistSourceToEarth
      if (abs(ratio) <= 1) then
         IndEGMF(n)=int((ratio+1.)*(nCube/2))+1
      else if (ratio > 0) then
         IndEGMF(n) = nCube/2 +1
      else
         IndEGMF(n) = 1
      end if
   end do

   !$OMP CRITICAL (egmf)
   if (sizeEGMF==0) then 
      call insert_value(1,indEGMF,b)
   else
      find_or_insert = .false.
      do n=1, sizeEGMF !while ((indEGMF(1)/=EGMFtab(n)%i) .and. (n<=sizeEGMF))
         if (indEGMF(1) == EGMFtab(n)%i) then 
            if (indEGMF(2) == EGMFtab(n)%j) then
               if(indEGMF(3) == EGMFtab(n)%k) then
                  B = EGMFtab(n)%b
                  find_or_insert = .true.
                  exit
               else if (indEGMF(3) < EGMFtab(n)%k) then
                  call insert_value(n,indEGMF,b)
                  find_or_insert = .true.
                  exit
               end if 
            else if (indEGMF(2) < EGMFtab(n)%j) then
               call insert_value(n,indEGMF,b)
               find_or_insert = .true.
               exit
            end if
         else if (indEGMF(1) < EGMFtab(n)%i) then
            call insert_value(n,indEGMF,b)
            find_or_insert = .true.
            exit
         end if
      end do
      if (.not. find_or_insert) then
         call insert_value(sizeEGMF+1,indEGMF,b)
      end if
   end if
   !$OMP end CRITICAL (egmf)
end subroutine selectEGMF

subroutine insert_value(n,ind,b)
   implicit none
   integer, intent(in) :: n,ind(3)
   real(prec), intent(out) :: b(3)
   type (EGMFdir), dimension(:), allocatable :: EGMFsave
   integer :: p, rp

   if (sizeEGMF == size(EGMFtab)) then
      allocate(EGMFsave(sizeEGMF))
      EGMFsave(1:sizeEGMF)=EGMFtab(1:sizeEGMF)
      deallocate(EGMFtab)
      allocate(EGMFtab(sizeEGMF+sizeini))
      EGMFtab(1:sizeEGMF)=EGMFsave(1:sizeEGMF)
      deallocate(EGMFsave)
      print*, " Increasing EGMF table to",size(EGMFtab)
   endif

   do p=n, sizeEGMF
      rp = sizeEGMF+n-p
      EGMFtab(rp+1) = EGMFtab(rp)
   end do
   sizeEGMF=sizeEGMF+1
   EGMFtab(n)%i = ind(1)
   EGMFtab(n)%j = ind(2)
   EGMFtab(n)%k = ind(3)
   call isotrop(B(1),B(2),B(3))
   call renorme(B(1),B(2),B(3))
   EGMFtab(n)%B = b
end subroutine insert_value

!===================== EGMF TURBULENT ===========================
#elif defined(B_turbulent)
subroutine init_b()
   ! ------------------------------------------------------------------------------
   ! This routine initialize the reference array once for the entire simulation
   ! with random values.
   ! ------------------------------------------------------------------------------
   !$ use constantes, only: OMP_num_threads
   use Constantes, only : Ethreshold
   implicit none
   integer :: i
   real(prec) :: k_min, k_max, alpha, alpha_n, G_kn, s, q, kn_norm, delta_kn
   real(prec) :: ran0
   external ran0

   s = 5./3.
   q = 0

   open(unit=111, file="Ak.txt")
   k_min = 1./lambdaEGMF ! Mpc-1
   ! Caution: k_max = 1/lambda_min > 1/R_L (using approx of R_L)
   k_max = 1d-2 *(intEGMF/1d-17)/(Ethreshold*510998.918*1d-12) ! Mpc-1
   alpha = (k_max/k_min)**(1/float(Nm)) ! log distrib norm
   !!$OMP PARALLEL do SCHEDULE(DYNAMIC) NUM_THREADS(OMP_num_threads)
   do i=1,Nm
      kn_norm = alpha**(i-1)
      delta_kn = (1-1/alpha)*kn_norm
      G_kn = kn_norm**q / (1+kn_norm**2)**((s+q)/2.)
      EGMFtab(i)%A_kn = G_kn*delta_kn
      EGMFtab(i)%k_n = kn_norm * k_min ! Mpc-1
      alpha_n = ran0()*2*pi
      EGMFtab(i)%phi_n = ran0()*2*pi
      EGMFtab(i)%theta_n = ran0()*pi
      EGMFtab(i)%Xi_n(1) =  cos(EGMFtab(i)%theta_n)*cos(EGMFtab(i)%phi_n)*cos(alpha_n) &
                            + sin(EGMFtab(i)%phi_n)*sin(alpha_n)
      EGMFtab(i)%Xi_n(2) =  cos(EGMFtab(i)%theta_n)*sin(EGMFtab(i)%phi_n)*cos(alpha_n) &
                            - cos(EGMFtab(i)%phi_n)*sin(alpha_n)
      EGMFtab(i)%Xi_n(3) =  -sin(EGMFtab(i)%theta_n)*cos(alpha_n)
      EGMFtab(i)%phase_n = ran0()*2*pi
      write(unit=111,fmt=*) EGMFtab(i)%k_n, delta_kn, G_kn, EGMFtab(i)%A_kn
   end do
   !!$OMP end PARALLEL do
   EGMFtab%A_kn = sqrt(EGMFtab%A_kn / sum(EGMFtab%A_kn))
   close(111)
End subroutine init_b

subroutine selectEGMF(leptonpos,B)
   implicit none
   real(prec), intent(out) :: B(3)
   real(prec), dimension(3) , intent(in) :: leptonpos
   integer :: i
   real(prec) :: z_prime, cos_kn_z

   b(1)=0
   b(2)=0
   b(3)=0
   do i=1, Nm
      z_prime = sin(EGMFtab(i)%theta_n)*cos(EGMFtab(i)%phi_n)*leptonpos(1) &
                + sin(EGMFtab(i)%theta_n)*sin(EGMFtab(i)%phi_n)*leptonpos(2) &
                + cos(EGMFtab(i)%theta_n)*leptonpos(3)
      z_prime = z_prime/Mpc
      cos_kn_z = cos(EGMFtab(i)%k_n*z_prime +EGMFtab(i)%phase_n)
      b=b+EGMFtab(i)%Xi_n*EGMFtab(i)%A_kn*cos_kn_z
   end do
   b= b*sqrt(2.)
   !call renorme(B(1),B(2),B(3))
end subroutine selectEGMF
#endif

!==============================================================================         
                              END MODULE EGMF
!==============================================================================         
