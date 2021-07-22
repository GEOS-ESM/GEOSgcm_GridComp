module mcica_subcol_gen_sw

   use cloud_condensate_inhomogeneity, only: &
      condensate_inhomogeneous, zcw_lookup
   use parrrsw, only : nbndsw, ngptsw
   use rrsw_con, only: grav
   use rrsw_wvn, only: ngb
   use rrsw_vsn
    
   implicit none

   public :: mcica_sw      
      
contains

   !-------------------------------------------------------------------------------------------------
   subroutine mcica_sw(ncol, nlay, nsubcol, play, cld, clwp, ciwp, &
                       cld_stoch, clwp_stoch, ciwp_stoch, &
                       CDF, CDF2, CDF3, alpha, zm, alat, DOY, rdl, adl) 
   !-------------------------------------------------------------------------------------------------
   ! Contact: Cecile Hannay (hannay@ucar.edu)
   ! 
   ! Original code: Based on Raisanen et al., QJRMS, 2004.
   !
   ! Modifications: Generalized for use with RRTMG and added Mersenne Twister as the default
   !   random number generator, which can be changed to the optional kissvec random number generator
   !   with flag 'irng'. Some extra functionality has been commented or removed.  
   !   Michael J. Iacono, AER, Inc., February 2007
   !
   ! Given a profile of cloud fraction, cloud water and cloud ice, we produce a set of subcolumns.
   ! Each layer within each subcolumn is homogeneous, with cloud fraction equal to zero or one 
   ! and uniform cloud liquid and cloud ice concentration.
   ! The ensemble as a whole reproduces the probability function of cloud liquid and ice within each layer 
   ! and obeys an overlap assumption in the vertical.   
   ! 
   ! Overlap assumption:
   !  The cloud are consistent with 4 overlap assumptions: random, maximum, maximum-random and exponential. 
   !  The default option is maximum-random (option 3)
   !  The options are: 1=random overlap, 2=max/random, 3=maximum overlap, 4=exponential overlap
   !  This is set with the variable "overlap" 
   !mji - Exponential overlap option (overlap=4) has been deactivated in this version
   !  The exponential overlap uses also a length scale, Zo. (real,    parameter  :: Zo = 2500. ) 
   ! 
   ! PDF assumption:
   !  We can use arbitrary complicated PDFS. 
   !  In the present version, we produce homogeneuous clouds (the simplest case).  
   !  Future developments include using the PDF scheme of Ben Johnson. 
   !
   ! History file:
   !  Option to add diagnostics variables in the history file. (using FINCL in the namelist)
   !  nsubcol = number of subcolumns
   !  overlap = overlap type (1-3)
   !  Zo = length scale 
   !  CLOUD_S = mean of the subcolumn cloud fraction ('_S" means Stochastic)
   !  CLDLIQ_S = mean of the subcolumn cloud water
   !  CLDICE_S = mean of the subcolumn cloud ice 
   !
   ! Note:
   !   Here: we force that the cloud condensate to be consistent with the cloud fraction 
   !   i.e we only have cloud condensate when the cell is cloudy. 
   !   In CAM: The cloud condensate and the cloud fraction are obtained from 2 different equations 
   !   and the 2 quantities can be inconsistent (i.e. CAM can produce cloud fraction 
   !   without cloud condensate or the opposite).
   !--------------------------------------------------------------------------------------------------------

      use mcica_random_numbers

      integer, intent(in) :: ncol            ! number of layers
      integer, intent(in) :: nlay            ! number of layers
      integer, intent(in) :: nsubcol         ! number of sub-columns (g-point intervals)

      ! Column state (cloud fraction, cloud water, cloud ice) + variables needed to read physics state 
      real , intent(in) :: play(:,:)          ! layer pressure (Pa)
                                                      !    Dimensions: (ncol,nlay)
      real , intent(in) :: cld(:,:)           ! cloud fraction 
                                                      !    Dimensions: (ncol,nlay)
      real , intent(in) :: clwp(:,:)          ! in-cloud liquid water path (g/m2)
                                                      !    Dimensions: (ncol,nlay)
      real , intent(in) :: ciwp(:,:)          ! in-cloud ice water path (g/m2)
                                                      !    Dimensions: (ncol,nlay)
      real , intent(out), dimension(:,:,:) :: CDF       
      real , intent(out), dimension(:,:,:) :: CDF2
      real , intent(out), dimension(:,:,:) :: CDF3
      real , intent(out), dimension(:,:) :: alpha
      
      real , intent(in), dimension(:) :: ALAT
      real , intent(in), dimension(:,:) :: zm
      integer , intent(in) :: DOY
      
      real , intent(out) :: cld_stoch(:,:,:)  ! subcolumn cloud fraction 
                                                      !    Dimensions: (ngptsw,ncol,nlay)
      real , intent(out) :: clwp_stoch(:,:,:) ! subcolumn in-cloud liquid water path
                                                      !    Dimensions: (ngptsw,ncol,nlay)
      real , intent(out) :: ciwp_stoch(:,:,:) ! subcolumn in-cloud ice water path
                                                      !    Dimensions: (ngptsw,ncol,nlay)
      
      ! -- Local variables

      real, intent(out) :: adl(:)                
      real, intent(out) :: rdl(:)
    
      ! Constants (min value for cloud fraction and cloud water and ice)
      real, parameter :: cldmin = 1.0e-20    ! min cloud fraction

      ! Variables related to random number and seed 
      integer :: seed1, seed2, seed3, seed4  ! seed (kissvec)
      real :: rand_num                       ! random number (kissvec)

      ! Indices
      integer :: ilev, isubcol, i, n, icol   ! indices
      real    :: am1,am2,am3,am4,amr
      real    :: sigma_qcw, zcw
      logical :: cond_inhomo
      !------------------------------------------------------------------------------------------ 

      ! save for speed
      cond_inhomo = condensate_inhomogeneous()

      !$acc kernels
      CDF  = 0.
      CDF2 = 0.
      CDF3 = 0.
      alpha = 0.
      !$acc end kernels

      ! Exponential overlap: weighting between maximum and random overlap increases with the distance. 
      ! The random numbers for exponential overlap verify:
      ! j=1   RAN(j)=RND1
      ! j>1   if RND1 < alpha(j,j-1) => RAN(j) = RAN(j-1)
      !                                 RAN(j) = RND2
      ! alpha is obtained from the equation
      ! alpha = exp(- (Zi-Zj-1)/Zo) where Zo is a characteristic length scale    

      ! -----------------------------------
      ! compute decorrelation length scales
      ! -----------------------------------

      ! cloud presence decorrelation length scale
      !$acc kernels
      am1 = 1.4315
      am2 = 2.1219
      am4 = -25.584
      amr = 7.
      if (doy .gt. 181) then
         am3 = -4.*amr/365*(doy-272)
      else
         am3 = 4.*amr/365.*(doy-91)
      endif
      do i = 1, ncol
         adl(i) = (am1+am2*exp(-(alat(i)*180./3.141592-am3)**2/am4**2))*1.e3  ! [m]
      end do

      ! condensate decorrelation length scale
      am1 = 0.72192
      am2 = 0.78996
      am4 = 40.404
      amr = 8.5
      if (doy .gt. 181) then
         am3 = -4.*amr/365*(doy-272)
      else
         am3 = 4.*amr/365.*(doy-91)
      endif
      do i = 1, ncol
         rdl(i) = (am1+am2*exp(-(alaT(i)*180./3.141592-am3)**2/am4**2))*1.e3  ! [m]
      end do
      !$acc end kernels      
   
      !$acc kernels
      do i = 1, ncol
         alpha(i,1) = 0.
         do ilev = 2, nlay
            alpha(i,ilev) = exp( -(zm(i,ilev) - zm(i,ilev-1)) / adl(i))
         end do
      end do
      !$acc end kernels
     
      !$acc kernels
      do ilev = 1,nlay
         do i = 1, ncol
            seed1 = (play(i,1) - int(play(i,1)))  * 100000000 - ilev
            seed2 = (play(i,2) - int(play(i,2)))  * 100000000 + ilev
            seed3 = (play(i,3) - int(play(i,3)))  * 100000000 + ilev * 6.2
            seed4 = (play(i,4) - int(play(i,4)))  * 100000000           
            do isubcol = 1,nsubcol
               call kissvec(seed1,seed2,seed3,seed4,rand_num)
               CDF(i,ilev,isubcol) = rand_num 
               call kissvec(seed1,seed2,seed3,seed4,rand_num)
               CDF2(i,ilev,isubcol) = rand_num 
            end do
         end do
      end do
      !$acc end kernels        
    
      !$acc kernels
      do icol = 1, ncol
         do isubcol = 1, nsubcol
            do ilev = 2, nlay
               if (CDF2(icol, ilev, isubcol) < alpha(icol,ilev)) then
                  CDF(icol, ilev, isubcol) = CDF(icol, ilev-1, isubcol) 
               end if
            end do
         end do
      end do
      !$acc end kernels
    
      if (cond_inhomo) then
        
         !$acc kernels  
         do i = 1, ncol
            alpha(i,1)=0.
            do ilev = 2, nlay
               alpha(i, ilev) = exp( -( zm(i, ilev) - zm(i, ilev-1))/rdl(i))
            end do
         end do
         !$acc end kernels
        
         !$acc kernels
         do ilev = 1,nlay
            do i = 1, ncol
               seed1 = (play(i,1) - int(play(i,1)))  * 100000000 - ilev
               seed2 = (play(i,2) - int(play(i,2)))  * 100000000 + ilev
               seed3 = (play(i,3) - int(play(i,3)))  * 100000000 + ilev * 6.2
               seed4 = (play(i,4) - int(play(i,4)))  * 100000000           
               do isubcol = 1,nsubcol
                  call kissvec(seed1,seed2,seed3,seed4,rand_num)
                  CDF2(i,ilev,isubcol) = rand_num 
                  call kissvec(seed1,seed2,seed3,seed4,rand_num)
                  CDF3(i,ilev,isubcol) = rand_num 
               end do
            end do
         end do
         !$acc end kernels

         !$acc kernels
         do icol = 1, ncol
            do isubcol = 1, nsubcol
               do ilev = 2, nlay
                  if (CDF2(icol, ilev, isubcol) < alpha(icol,ilev)) then
                     CDF3(icol, ilev, isubcol) = CDF3(icol, ilev-1, isubcol)
                  end if
               end do
            end do
         end do
         !$acc end kernels        
      end if

      ! where the subcolumn is cloudy, the subcolumn cloud fraction is 1;
      ! where the subcolumn is not cloudy, the subcolumn cloud fraction is 0;
      ! where there is a cloud, define the subcolumn cloud properties,
      ! otherwise set these to zero

      !$acc kernels 
      do ilev = 1,nlay
         do i = 1, ncol
            do isubcol = 1, nsubcol
               if (cond_inhomo .and. cdf(i,ilev,isubcol) >= (1. - cld(i,ilev))) then

                  cld_stoch(i,ilev,isubcol) = 1. 
                  
                  ! Cloud fraction sets level of inhomogeneity
                  if (cld(i,ilev) .gt. 0.99) then
                     sigma_qcw = 0.5
                  elseif (cld(i,ilev) .gt. 0.9) then
                     sigma_qcw = 0.71
                  else  
                     sigma_qcw = 1.0
                  endif
                  
                  ! horizontally variable clouds
                  zcw = zcw_lookup(cdf3(i,ilev,isubcol),sigma_qcw)
                  clwp_stoch(i,ilev,isubcol) = clwp(i,ilev) * zcw
                  ciwp_stoch(i,ilev,isubcol) = ciwp(i,ilev) * zcw
                
               elseif (cdf(i,ilev,isubcol) >= (1. - cld(i,ilev))) then

                   cld_stoch(i,ilev,isubcol) = 1. 
                  clwp_stoch(i,ilev,isubcol) = clwp(i,ilev)
                  ciwp_stoch(i,ilev,isubcol) = ciwp(i,ilev)

               else

                  ! a clear subcolumn
                   cld_stoch(i,ilev,isubcol) = 0. 
                  clwp_stoch(i,ilev,isubcol) = 0. 
                  ciwp_stoch(i,ilev,isubcol) = 0. 

               endif
            enddo
         enddo
      enddo
      !$acc end kernels

   end subroutine mcica_sw

   !------------------------------------------------------
   subroutine kissvec(seed1, seed2, seed3, seed4, ran_num)
   !------------------------------------------------------

      real, intent(inout) :: ran_num
      integer, intent(inout) :: seed1, seed2, seed3, seed4
!     integer :: m, k, n, kiss
      integer :: kiss

!     ! inline function
!     m(k, n) = ieor (k, ishft (k, n) )

     !seed1 = 69069 * seed1 + 1327217885 ! in LW
      seed1 = 69069 * seed1 + 132721785  ! in SW
     !seed2 = m (m (m (seed2, 13), - 17), 5)                    ! in LW
      seed2 = 11002 * iand (seed2, 65535) + ishft (seed2, -16)  ! in SW
      seed3 = 18000 * iand (seed3, 65535) + ishft (seed3, -16)
      seed4 = 30903 * iand (seed4, 65535) + ishft (seed4, -16)
      kiss = seed1 + seed2 + ishft (seed3, 16) + seed4
      ran_num = kiss * 2.328306e-10 + 0.5

   end subroutine kissvec

end module mcica_subcol_gen_sw
