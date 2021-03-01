!     path:      $Source$
!     author:    $Author$
!     revision:  $Revision$
!     created:   $Date$
!

module mcica_subcol_gen_sw
      use parrrsw, only : nbndsw, ngptsw
      use rrsw_con, only: grav
      use rrsw_wvn, only: ngb
      use rrsw_vsn
      use tab_xcw

    
      implicit none

      public :: mcica_sw      
      
      contains
!-------------------------------------------------------------------------------------------------
      subroutine mcica_sw(ncol, nlay, nsubcol, icld, irng, play, cld, clwp, ciwp, &
                               tauc, ssac, asmc, fsfc, cld_stoch, clwp_stoch, ciwp_stoch, &
                               tauc_stoch, ssac_stoch, asmc_stoch, fsfc_stoch, changeSeed, CDF,&
                               CDF2, CDF3, alpha, zm, alat, DOY, rdl, adl) 
!-------------------------------------------------------------------------------------------------

  !----------------------------------------------------------------------------------------------------------------
  ! ---------------------
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
  ! Seed:
  !  If the stochastic cloud generator is called several times during the same timestep, 
  !  one should change the seed between the call to insure that the subcolumns are different.
  !  This is done by changing the argument 'changeSeed'
  !  For example, if one wants to create a set of columns for the shortwave and another set for the longwave ,
  !  use 'changeSeed = 1' for the first call and'changeSeed = 2' for the second call 
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
  !---------------------------------------------------------------------------------------------------------------

      use mcica_random_numbers
! The Mersenne Twister random number engine
      !use MersenneTwister, only: randomNumberSequence, &   
      !                           new_RandomNumberSequence, getRandomReal

      !type(randomNumberSequence) :: randomNumbers

! -- Arguments

      integer , intent(in) :: ncol            ! number of layers
      integer , intent(in) :: nlay            ! number of layers
      integer , intent(in) :: icld            ! clear/cloud, cloud overlap flag
      integer , intent(in) :: irng         ! flag for random number generator
                                                      !  0 = kissvec
                                                      !  1 = Mersenne Twister
      integer , intent(in) :: nsubcol         ! number of sub-columns (g-point intervals)
      integer , optional, intent(in) :: changeSeed     ! allows permuting seed

! Column state (cloud fraction, cloud water, cloud ice) + variables needed to read physics state 
      real , intent(in) :: play(:,:)          ! layer pressure (Pa)
                                                      !    Dimensions: (ncol,nlay)
      real , intent(in) :: cld(:,:)           ! cloud fraction 
                                                      !    Dimensions: (ncol,nlay)
      real , intent(in) :: clwp(:,:)          ! in-cloud liquid water path (g/m2)
                                                      !    Dimensions: (ncol,nlay)
      real , intent(in) :: ciwp(:,:)          ! in-cloud ice water path (g/m2)
                                                      !    Dimensions: (ncol,nlay)
      real , intent(in) :: tauc(:,:,:)        ! in-cloud optical depth (non-delta scaled)
                                                      !    Dimensions: (nbndsw,ncol,nlay)
      real , intent(in) :: ssac(:,:,:)        ! in-cloud single scattering albedo (non-delta scaled)
                                                      !    Dimensions: (nbndsw,ncol,nlay)
      real , intent(in) :: asmc(:,:,:)        ! in-cloud asymmetry parameter (non-delta scaled)
                                                      !    Dimensions: (nbndsw,ncol,nlay)
      real , intent(in) :: fsfc(:,:,:)        ! in-cloud forward scattering fraction (non-delta scaled)
                                                      !    Dimensions: (nbndsw,ncol,nlay)
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
      real , intent(out) :: tauc_stoch(:,:,:) ! subcolumn in-cloud optical depth
                                                      !    Dimensions: (ngptsw,ncol,nlay)
      real , intent(out) :: ssac_stoch(:,:,:) ! subcolumn in-cloud single scattering albedo
                                                      !    Dimensions: (ngptsw,ncol,nlay)
      real , intent(out) :: asmc_stoch(:,:,:) ! subcolumn in-cloud asymmetry parameter
                                                      !    Dimensions: (ngptsw,ncol,nlay)
      real , intent(out) :: fsfc_stoch(:,:,:) ! subcolumn in-cloud forward scattering fraction
                                                      !    Dimensions: (ngptsw,ncol,nlay)
      
! -- Local variables
      real, intent(out)  :: adl(:)                
      real, intent(out)  :: rdl(:)
    

! Set overlap
      integer  :: overlap                     ! 1 = random overlap, 2 = maximum/random,


! Constants (min value for cloud fraction and cloud water and ice)
      real , parameter :: cldmin = 1.0e-20  ! min cloud fraction


! Variables related to random number and seed 
     
      integer :: seed1, seed2, seed3, seed4  ! seed to create random number
    
      integer  :: iseed                        ! seed to create random number (Mersenne Twister)
      real  :: rand_num_mt                     ! random number (Mersenne Twister)
      real  :: kiss

! Indices
      integer  :: ilev, isubcol, i, n, ngbm, iplon   ! indices
      integer  :: inhm
      real     :: am1,am2,am3,am4,amr
      real     :: SIGMA_QCW
      real     :: rind1, rind2
      integer  :: ind1, ind2
      real     :: zcw
!------------------------------------------------------------------------------------------ 

! Check that irng is in bounds; if not, set to default
   

! Pass input cloud overlap setting to local variable
      overlap = icld
      inhm = 0

!$acc kernels
   CDF = 0.0
   CDF2 = 0.0
   CDF3 = 0.0
   alpha = 0.0
!$acc end kernels

! ------ Apply overlap assumption --------

! generate the random numbers  

 
if (icld==1) then
!$acc kernels 
   do ilev = 1,nlay
      do i = 1, ncol
         seed1 = (play(i,1) - int(play(i,1)))  * 100000000 - ilev
         seed2 = (play(i,2) - int(play(i,2)))  * 100000000 + ilev
         seed3 = (play(i,3) - int(play(i,3)))  * 100000000 + ilev * 6.2
         seed4 = (play(i,4) - int(play(i,4)))  * 100000000           
         do isubcol = 1,nsubcol
            seed1 = 69069  * seed1 + 132721785 
            seed2 = 11002  * iand (seed2, 65535 ) + ishft (seed2, - 16 )
            seed3 = 18000  * iand (seed3, 65535 ) + ishft (seed3, - 16 )
            seed4 = 30903  * iand (seed4, 65535 ) + ishft (seed4, - 16 )
            kiss = seed1 + seed2 + ishft (seed3, 16 ) + seed4
            CDF(i,ilev,isubcol) = kiss*2.328306e-10  + 0.5 
         end do
      end do
   end do
!$acc end kernels      
elseif (icld==4) then
      ! Exponential overlap: weighting between maximum and random overlap increases with the distance. 
      ! The random numbers for exponential overlap verify:
      ! j=1   RAN(j)=RND1
      ! j>1   if RND1 < alpha(j,j-1) => RAN(j) = RAN(j-1)
      !                                 RAN(j) = RND2
      ! alpha is obtained from the equation
      ! alpha = exp(- (Zi-Zj-1)/Zo) where Zo is a characteristic length scale    

!$acc kernels
      AM1 = 1.4315
      AM2 = 2.1219
      AM4 = -25.584
      AMR = 7.
      IF ( DOY .GT. 181 ) THEN
         AM3 = -4.*AMR/365*(DOY-272)
      ELSE
         AM3 = 4.*AMR/365.*(DOY-91)
      ENDIF
      do i = 1, ncol
         adl(i) = AM1+AM2*EXP(-(ALAT(i)*180./3.141592-AM3)**2/AM4**2)
         adl(i) = adl(i)*1.e3
      end do
      AM1 = 0.72192
      AM2 = 0.78996
      AM4 = 40.404
      AMR = 8.5
      IF ( DOY .GT. 181 ) THEN
         AM3 = -4.*AMR/365*(DOY-272)
      ELSE
         AM3 = 4.*AMR/365.*(DOY-91)
      ENDIF
      do i = 1, ncol
         rdl(i) = AM1+AM2*EXP(-(ALAT(i)*180./3.141592-AM3)**2/AM4**2)
         rdl(i) = rdl(i)*1.e3
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
               seed1 = 69069  * seed1 + 132721785 
               seed2 = 11002  * iand (seed2, 65535 ) + ishft (seed2, - 16 )
               seed3 = 18000  * iand (seed3, 65535 ) + ishft (seed3, - 16 )
               seed4 = 30903  * iand (seed4, 65535 ) + ishft (seed4, - 16 )
               kiss = seed1 + seed2 + ishft (seed3, 16 ) + seed4
               CDF(i,ilev,isubcol) = kiss*2.328306e-10  + 0.5 
               ! for testing purposes to generate determanistic sequence
               !CDF(i,ilev,isubcol) = real(MOD((isubcol * 37 + ilev * 71) , 10000))/10000.0
               
               seed1 = 69069  * seed1 + 132721785 
               seed2 = 11002  * iand (seed2, 65535 ) + ishft (seed2, - 16 )
               seed3 = 18000  * iand (seed3, 65535 ) + ishft (seed3, - 16 )
               seed4 = 30903  * iand (seed4, 65535 ) + ishft (seed4, - 16 )
               kiss = seed1 + seed2 + ishft (seed3, 16 ) + seed4
               CDF2(i,ilev,isubcol) = kiss*2.328306e-10  + 0.5 
               ! for testing purposes to generate determanistic sequence
               !CDF2(i,ilev,isubcol) =  real(MOD((isubcol * 37 + ilev * 73) , 10000))/10000.0
            end do
         end do
      end do
!$acc end kernels        
       
!$acc kernels
   do iplon = 1, ncol
      do isubcol = 1, nsubcol
         do ilev = 2, nlay
            if (CDF2(iplon, ilev, isubcol) < alpha(iplon,ilev)) then
               CDF(iplon, ilev, isubcol) = CDF(iplon, ilev-1, isubcol) 
            end if
         end do
      end do
   end do
!$acc end kernels
   inhm = 1
    
end if

    if (inhm>=1) then
        
      
        
   
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
               seed1 = 69069  * seed1 + 132721785 
               seed2 = 11002  * iand (seed2, 65535 ) + ishft (seed2, - 16 )
               seed3 = 18000  * iand (seed3, 65535 ) + ishft (seed3, - 16 )
               seed4 = 30903  * iand (seed4, 65535 ) + ishft (seed4, - 16 )
               kiss = seed1 + seed2 + ishft (seed3, 16 ) + seed4
               CDF2(i,ilev,isubcol) = kiss*2.328306e-10  + 0.5 
               ! for testing purposes to generate determanistic sequence
               !CDF2(i,ilev,isubcol) = real(MOD((isubcol * 7 + ilev * 101) , 10000))/10000.0
               
               seed1 = 69069  * seed1 + 132721785 
               seed2 = 11002  * iand (seed2, 65535 ) + ishft (seed2, - 16 )
               seed3 = 18000  * iand (seed3, 65535 ) + ishft (seed3, - 16 )
               seed4 = 30903  * iand (seed4, 65535 ) + ishft (seed4, - 16 )
               kiss = seed1 + seed2 + ishft (seed3, 16 ) + seed4
               CDF3(i,ilev,isubcol) = kiss*2.328306e-10  + 0.5 
               ! for testing purposes to generate determanistic sequence
               !CDF3(i,ilev,isubcol) =  real(MOD((isubcol * 41 + ilev * 89) , 10000))/10000.0
            end do
         end do
      end do
!$acc end kernels

!$acc kernels
      do iplon = 1, ncol
         do isubcol = 1, nsubcol
            do ilev = 2, nlay
               if (CDF2(iplon, ilev, isubcol) < alpha(iplon,ilev)) then
                  CDF3(iplon, ilev, isubcol) = CDF3(iplon, ilev-1, isubcol)
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

      ngbm = ngb(1) - 1
!$acc kernels 
      do ilev = 1,nlay
         do i = 1, ncol
            do isubcol = 1, nsubcol
               if ( inhm >= 1 .and. CDF(i,ilev,isubcol)>=(1.0 - cld(i,ilev))) then
                  cld_stoch(i,ilev,isubcol) = 1.0 
                  
                  if (cld(i,ilev) .gt. 0.99 ) then
                    SIGMA_QCW = 0.5
                  elseif (cld(i,ilev) .gt. 0.9 ) then
                    SIGMA_QCW = 0.71
                    else  
                    SIGMA_QCW = 1
                  endif
                  
                  RIND1 = CDF3(i,ilev,isubcol) * (N1 - 1) + 1.0
                  IND1  = MAX(1, MIN(INT(RIND1), N1-1))
                  RIND1 = RIND1 - IND1
                  RIND2 = 40.0 * SIGMA_QCW    - 3.0
                  IND2  = MAX(1, MIN(INT(RIND2), N2-1))
                  RIND2 = RIND2 - IND2

                   ZCW = (1.0-RIND1) * (1.0-RIND2) * XCW(IND1,IND2) &
                       + (1.0-RIND1) * RIND2       * XCW(IND1,IND2+1) &
                       + RIND1 * (1.0-RIND2)       * XCW(IND1+1,IND2) &
                       + RIND1 * RIND2             * XCW(IND1+1,IND2+1)
                   
                 

                  
                  clwp_stoch(i,ilev,isubcol) = clwp(i,ilev) * ZCW
                  ciwp_stoch(i,ilev,isubcol) = ciwp(i,ilev) * ZCW
                  n = ngb(isubcol) - ngbm
                  tauc_stoch(i,ilev,isubcol) = tauc(i,ilev,n)
                  ssac_stoch(i,ilev,isubcol) = ssac(i,ilev,n)
                  asmc_stoch(i,ilev,isubcol) = asmc(i,ilev,n)
                  fsfc_stoch(i,ilev,isubcol) = fsfc(i,ilev,n)
                
                
               elseif ( CDF(i,ilev,isubcol)>=(1.0 - cld(i,ilev)) ) then
                  cld_stoch(i,ilev,isubcol) = 1.0 
                  clwp_stoch(i,ilev,isubcol) = clwp(i,ilev)
                  ciwp_stoch(i,ilev,isubcol) = ciwp(i,ilev)
                  n = ngb(isubcol) - ngbm
                  tauc_stoch(i,ilev,isubcol) = tauc(i,ilev,n)
                  ssac_stoch(i,ilev,isubcol) = ssac(i,ilev,n)
                  asmc_stoch(i,ilev,isubcol) = asmc(i,ilev,n)
                  fsfc_stoch(i,ilev,isubcol) = fsfc(i,ilev,n)
               else
                  cld_stoch(i,ilev,isubcol) = 0. 
                  clwp_stoch(i,ilev,isubcol) = 0. 
                  ciwp_stoch(i,ilev,isubcol) = 0. 
                  tauc_stoch(i,ilev,isubcol) = 0. 
                  ssac_stoch(i,ilev,isubcol) = 1. 
                  asmc_stoch(i,ilev,isubcol) = 0. 
                  fsfc_stoch(i,ilev,isubcol) = 0. 
               endif
            enddo
         enddo
      enddo
!$acc end kernels



      end subroutine mcica_sw


      end module mcica_subcol_gen_sw


