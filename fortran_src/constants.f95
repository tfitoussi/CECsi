!===============================================================================         
                           MODULE Constantes
!===============================================================================         
implicit none
! Precision
integer, parameter :: prec=8

! Physical constantes
real(prec), parameter :: e=4.8032068D-10 ! esu (cgs units)
real(prec), parameter :: m=9.1093897D-28 ! g
real(prec), parameter :: c=2.99792458D+10 ! cm.s-1
real(prec), parameter :: k=1.380658E-16 ! erg.K-1
real(prec), parameter :: h=6.62606957D-27 ! erg.s
real(prec), parameter :: hb=1.05457266E-27 ! erg.s
real(prec), parameter :: Tcmb=2.725 ! K
real(prec), parameter :: Tcmbp=(k*Tcmb)/(m*c*c) ! Adimentionnal thermal energy
real(prec), parameter :: pi=4*atan(1.) 
real(prec), parameter :: alpha=1./137.035999679
real(prec), parameter :: lambda=3.86E-11 ! cm 
real(prec), parameter :: lambdaC=hb/(m*c)
real(prec), parameter :: Mpc=(3.0856776D+16)*1d8 ! Mpc to cm
real(prec), parameter :: eV_to_erg=1.602176565e-12
real(prec), parameter :: s_to_yr=3600*24*365.25

! Cosmology
real(prec), parameter :: a0=1.
real(prec), parameter :: H0=67.8*1d5/(Mpc) ! s-1
real(prec), parameter :: omegaM=0.3
real(prec), parameter :: omegaK=0.
real(prec), parameter :: omegaL=0.7

! Particles physics
real(prec), parameter :: r0=e*e/(m*c*c)
real(prec), parameter :: sigmaT=8.*pi*(r0**2)/3. ! Thomson cross section

! Gauss-Legendre integration
integer, parameter :: ngp=30 !number of points 
real(prec) :: xabsc(ngp), weig(ngp) ! position and weight 
integer, parameter :: n_bin = 81

integer :: iseed, seed ! for the random number generation
integer :: OMP_num_threads
real(prec) :: alphaS, Ethreshold, Eth_lept

! data from ELB files
real(prec),dimension(:), allocatable :: redshiftEBL,enerEBL 
real(prec),dimension(:,:), allocatable :: densityEBL 
real(prec) :: maxEnerEBL, minEnerEBL, maxRedshiftEBL, minRedshiftEBL

! Some usefull globales variables
real(prec) :: DistSourceToEarth
real(prec) :: currentRedshift, currentEnergy

! Need for the Compton accumulation
real(prec) :: xiComp=1., ComptonThreshold, ComptonEnergyThreshold
integer :: NbErrors=0 ! number of errors due to Compton Acc. threshold too high
integer :: NbCompton=0, NbPairsprod=0 
integer :: NbPhotonsProd=0, NbLeptonsProd=0 

! Files IDs
integer, parameter :: ResultID=11 
integer, parameter :: errorID=12
integer, parameter :: deltaID=13
integer, parameter :: positronsID=14, positronsLoosedID=15 

!$OMP THREADPRIVATE(currentRedshift,currentEnergy,xiComp,iseed)
!===============================================================================         
                           END MODULE Constantes
!===============================================================================         
