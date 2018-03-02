#include "../temp/preprocessing.f95"
!===============================================================
subroutine  ReadEBLfiles()
   use Constantes, only : redshiftEBL, densityEBL, enerEBL, prec, c, h, hb, pi, m, &
                 maxEnerEBL, minEnerEBL, maxRedshiftEBL, minRedshiftEBL, eV_to_erg 
   implicit none
   real(prec), dimension(:), allocatable :: hv
   real(prec), dimension(:,:), allocatable :: mat, nezz
   integer :: i, j, k, n_EBL1, n_EBL2

   write(*,*) "EBL model ========================================================"
#if defined(EBL_model_Dominguez)
   write(*,*) " > Model from Dominguez et Al; "
   write(*,*) "   data taken from http://side.iaa.es/EBL/"
   n_EBL1 = 50
   n_EBL2 = 18
#elif defined(EBL_model_Finke)
   write(*,*) " > Model C from Finke et al. (2010), arXiv:0905.1115; "
   write(*,*) "   data taken from http://www.phy.ohiou.edu/~finke/EBL/"
   n_EBL1 = 397
   n_EBL2 = 500
#elif defined(EBL_model_Franceschini)
   write(*,*) " > Franceschini model, arXiv:0805.1841; "
   write(*,*) "   data taken directly from paper"
   n_EBL1 = 31
   n_EBL2 = 11
#elif defined(EBL_model_Gilmore)
   write(*,*) " > Model from Gilmore et al. (2012), arXiv:1104.0671; "
   write(*,*) "   data taken from physics.ucsc.edu/~joel/EBLdata-Gilmore2012/"
   n_EBL1 = 101
   n_EBL2 = 18
#elif defined(EBL_model_Lower_limit)
   write(*,*) " > 'lower limit' model of Kneiske and dole"
   n_EBL1 = 160
   n_EBL2 = 200
#elif defined(EBL_model_Best_fit)
   write(*,*) " > 'best-fit' model of Kneiske and dole"
   n_EBL1 = 201
   n_EBL2 = n_EBL1
#endif
   allocate (redshiftEBL(n_EBL2))
   allocate (enerEBL(n_EBL1))
   allocate (densityEBL(n_EBL1,n_EBL2))
   allocate(nezz(n_EBL1,n_EBL2))
   allocate(hv(n_EBL1))

#if defined(EBL_model_Dominguez)
   ! file z_EBL contains redshift only
   open(unit=60,file="EBL_files/z_Dominguez.dat")
   do i=1,size(redshiftEBL)
      read(60,*) redshiftEBL(i)
   end do
   close(unit=60)

   ! file EBL contains:
   !  - 1st column: lambda
   !  - other columns: lambda * I_lambda
   open(unit=61,file="EBL_files/lambdaI_Dominguez.dat")
   allocate(mat(n_EBL1,n_EBL2+1))
   do i=1,size(mat,1)
      read(61,*)  mat(i,:)
      ! Convertion: lambda -> erg
      hv(i) = h*c/(mat(i,1)*1.D-4)
      do j=2,size(mat,2)
         ! Convertion lambda*I_lambda -> photon / cm^3 / erg
         nezz(i,j-1) = (1.D-6)*(4*pi/c)*mat(i,j)/(hv(i)**2)
      end do
   end do
   close(unit=61)
   deallocate(mat)
#elif defined(EBL_model_Finke)
   do i=1,n_EBL2
      redshiftEBL(i)=0.01d0*(i-1)
   enddo    

   open(1,file='EBL_files/n_Finke.dat',status='old')
   do i=1,n_EBL1
      read(1,*) hv(i),nezz(i,:)  
   end do
   close(unit=1)
#elif defined(EBL_model_Franceschini)
   do i=1,n_EBL2
      redshiftEBL(i)=0.2d0*(i-1)
   enddo    

   open(1,file='EBL_files/n_Fra_2017.dat',status='old')
   do i=1,n_EBL1
      read(1,*) hv(i),nezz(i,:)  
   end do
   close(unit=1)
#elif defined(EBL_model_Gilmore)
   open(1,file='EBL_files/z-IR_Gil.dat',status='old')
   do i=1,n_EBL2
      read(1,*) k,redshiftEBL(i)
   enddo
   close(1)     

   open(1,file='EBL_files/n_Gil.dat',status='old')
   do i=1,n_EBL1
      read(1,*) hv(i),nezz(i,:)  
   enddo
   close(1)   
#elif defined(EBL_model_Lower_limit)
   open(1,file='EBL_files/z-IR.dat',status='old')
   do i=1,n_EBL2
      read(1,*) k,redshiftEBL(i)
   enddo
   close(1)     

   open(1,file='EBL_files/n_lowerlimit10.dat',status='old')
   do i=1,n_EBL1
      read(1,*) hv(i),nezz(i,:)  
   enddo
   close(1)   
#elif defined(EBL_model_Best_fit)
   do i=1,n_EBL2
      redshiftEBL(i)=0.025d0*(i-1)
   enddo    

   open(1,file='EBL_files/n_bestfit10.dat',status='old')
   do i=1,n_EBL1
      read(1,*) enerEBL(i),densityEBL(i,:)  
   enddo
   close(1)     
#endif

#ifndef EBL_model_Best_fit
   ! inversion to have hnu increasing: i -> size(densityEBL,1)-i+1
   do i=1,size(densityEBL,1)
      enerEBL(i)=hv(size(nezz,1)-i+1)
      densityEBL(i,:)=nezz(size(nezz,1)-i+1,:)
   end do
#endif

   deallocate(nezz)
   deallocate(hv)

   ! adimenssionnal 
#ifdef EBL_model_Dominguez
   enerEBL = enerEBL /(m*c*c)
   densityEBL = densityEBL *((hb*c)**3)/((m*c*c)**2)
#else
   enerEBL = enerEBL *eV_to_erg/(m*c*c)
   densityEBL = densityEBL *((hb*c)**3)/((m*c*c)**2)/eV_to_erg 
#endif

   maxEnerEBL = maxval(enerEBL)
   minEnerEBL = minval(enerEBL)
   maxRedshiftEBL = maxval(redshiftEBL)
   minRedshiftEBL = minval(redshiftEBL)
  
end subroutine ReadEBLfiles

! interpolation of EBL data to find EBL density
subroutine  NeInterpol(epss,redddz, EBLinterpol)
   use Constantes
   implicit none
   real(prec) :: fx, fy
   real(prec), intent(in) :: epss, redddz 
   real(prec), intent(out) :: EBLinterpol 
   integer ::  leftE, rightE, leftZ, rightZ

   !Si le point est hors grille des données, on met n=0
   if (epss>maxenerEBL .OR. epss<minenerEBL) then
      EBLinterpol=0. 

   else if ( redddz<minredshiftEBL ) then
      CAll r8vec_bracket(size(enerEBL), enerEBL, epss, leftE, rightE )
      fx=(epss-enerEBL(leftE))/(enerEBL(rightE)-enerEBL(leftE))
      EBLinterpol=((1.-fx)*densityEBL(leftE,1)+ fx*densityEBL(rightE,1))

   else if ( redddz>maxredshiftEBL ) then
      CAll r8vec_bracket(size(enerEBL), enerEBL, epss, leftE, rightE )
      fx=(epss-enerEBL(leftE))/(enerEBL(rightE)-enerEBL(leftE))
      leftZ = size(redshiftEBL)
      EBLinterpol=((1.-fx)*densityEBL(leftE,leftZ)+ fx*densityEBL(rightE,leftZ))

   else !recherche des (leftE,rightE) et (leftZ,rightZ) qui encadrent les valeurs 
      !de epss et reddz demandées
      CAll r8vec_bracket(size(enerEBL), enerEBL, epss, leftE, rightE )
      CAll r8vec_bracket(size(redshiftEBL), redshiftEBL, redddz, leftZ, rightZ )

      fx=(epss-enerEBL(leftE))/(enerEBL(rightE)-enerEBL(leftE))
      fy=(redddz-redshiftEBL(leftZ))/(redshiftEBL(rightZ)-redshiftEBL(leftZ))
      EBLinterpol=(1.-fx)*(1.-fy)*densityEBL(leftE,leftZ)+&
         (1.-fx)*fy*densityEBL(leftE,rightZ)+&
         fx*(1.-fy)*densityEBL(rightE,leftZ)+&
         fx*fy*densityEBL(rightE,rightZ)
   end if
#ifndef EBL_model_Gilmore
   EBLinterpol=((1.+redddz)**3)*EBLinterpol 
#endif
end  subroutine NeInterpol 

subroutine r8vec_bracket ( n, x, xval, left, right )
!*****************************************************************************80
!
!! R8VEC_BRACKET searches a sorted R8VEC for successive brackets of a value.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    If the values in the vector are thought of as defining intervals
!    on the real line, then this routine searches for the interval
!    nearest to or containing the given value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Constantes:
!
!    Input, integer ( kind = 4 ) N, length of input array.
!
!    Input, real ( kind = 8 ) X(N), an array sorted into ascending order.
!
!    Input, real ( kind = 8 ) XVAL, a value to be bracketed.
!
!    Output, integer ( kind = 4 ) LEFT, RIGHT, the results of the search.
!    Either:
!      XVAL < X(1), when LEFT = 1, RIGHT = 2;
!      X(N) < XVAL, when LEFT = N-1, RIGHT = N;
!    or
!      X(LEFT) <= XVAL <= X(RIGHT).

   use Constantes, only : prec
   implicit none
   integer ( kind = 4 ) n
   integer ( kind = 4 ) i
   integer ( kind = 4 ) left
   integer ( kind = 4 ) right
   real ( kind = prec ) x(n)
   real ( kind = prec) xval

   do i = 2, n - 1
      if ( xval < x(i) ) then
         left = i - 1
         right = i
         return
      end if
   end do

   left = n - 1
   right = n

   return
end
