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

!-------------------------------------------------------------------------------------------------
! Original code: Based on Raisanen et al., QJRMS, 2004.
! Original Contact: Cecile Hannay (hannay@ucar.edu)
!
! Modifications:
! ...
! Peter Norris, GMAO, Apr 2021:
!   Tidy up code and comments
!   Keep only exponential (generalized) overlap.
!
! Given a profile of cloud fraction, cloud water and cloud ice, produce a set of subcolumns.
! For ih == 0, each layer within each subcolumn is homogeneous, with cloud frac zero or one
!   and uniform cloud liquid and cloud ice concentration.
! For ih > 0, each layer has horizontal condensate variability.
! The ensemble of subcolumns statistically reproduces the cloud fraction within each layer
!   and its prescribed vertical overlap, and, for ih > 0, the PDF of cloud liquid and ice
!   within each layer and its prescribed vertical correlation.
!
! Overlap assumption:
! Exponential (generalized) overlap (Raisannen, Pincus) using pre-calculated correlations
! alphad and rcorrd based on decorrelation length scales from Oreopoulos et al.
!
!------------------------------------------------------------------------------------------------

! ----- Arguments -----

      integer , intent(in) :: ncol            ! number of layers
      integer , intent(in) :: nlay            ! number of layers
      integer , intent(in) :: icld            ! clear/cloud, cloud overlap flag
      integer , intent(in) :: irng            ! flag for random number generator
                                                 !  0 = kissvec (ONLY ONE WORKING)
                                                 !  1 = Mersenne Twister
      integer , intent(in) :: nsubcol         ! number of sub-columns (g-point intervals)

      integer , optional, intent(in) :: changeSeed     ! allows permuting seed

      ! grid-column state
      real , intent(in) :: play(:,:)          ! layer pressure (Pa)             (ncol,nlay)
      real , intent(in) :: cld(:,:)           ! cld fraction                    (ncol,nlay)
      real , intent(in) :: clwp(:,:)          ! in-cld liquid water path (g/m2) (ncol,nlay)
      real , intent(in) :: ciwp(:,:)          ! in-cld ice    water path (g/m2) (ncol,nlay)
      ! non-delta scaled cloud optical parameters 
      real , intent(in) :: tauc(:,:,:)        ! in-cld optical depth            (nbndsw,ncol,nlay)
      real , intent(in) :: ssac(:,:,:)        ! in-cld single scattering albedo (nbndsw,ncol,nlay)
      real , intent(in) :: asmc(:,:,:)        ! in-cld asymmetry parameter      (nbndsw,ncol,nlay)
      real , intent(in) :: fsfc(:,:,:)        ! in-cld forward scattering frac  (nbndsw,ncol,nlay)

      ! random numbers used for overlap
      real , intent(out), dimension(:,:,:) :: CDF       
      real , intent(out), dimension(:,:,:) :: CDF2
      real , intent(out), dimension(:,:,:) :: CDF3

      ! vertical decorrelation length calcs
      real,    intent(in),  dimension(:)   :: ALAT
      integer, intent(in)                  :: DOY
      real,    intent(in),  dimension(:,:) :: zm
      real,    intent(out), dimension(:)   :: adl
      real,    intent(out), dimension(:)   :: rdl
      real,    intent(out), dimension(:,:) :: alpha
      
      ! output stochastic subcolumns
      real , intent(out) :: cld_stoch(:,:,:)  ! subcol cld fraction                    (ngptsw,ncol,nlay)
      real , intent(out) :: clwp_stoch(:,:,:) ! subcol in-cld liquid water path        (ngptsw,ncol,nlay)
      real , intent(out) :: ciwp_stoch(:,:,:) ! subcol in-cld ice water path           (ngptsw,ncol,nlay)
      real , intent(out) :: tauc_stoch(:,:,:) ! subcol in-cld optical depth            (ngptsw,ncol,nlay)
      real , intent(out) :: ssac_stoch(:,:,:) ! subcol in-cld single scattering albedo (ngptsw,ncol,nlay)
      real , intent(out) :: asmc_stoch(:,:,:) ! subcol in-cld asymmetry parameter      (ngptsw,ncol,nlay)
      real , intent(out) :: fsfc_stoch(:,:,:) ! subcol in-cld forward scattering frac  (ngptsw,ncol,nlay)
      
! ----- Local variables -----

      ! Variables related to random number and seed
      integer :: seed1, seed2, seed3, seed4  ! seed (kissvec)
      real  :: kiss
      real :: rand_num                       ! random number (kissvec)

      ! Locals
      real :: am1,am2,am3,am4,amr
      real :: RIND1, RIND2, ZCW, SIGMA_QCW
      integer :: IND1, IND2
      integer :: inhm

      ! Indices
      integer :: iplon, isubcol, ilev, ibnd, ngbm

!------------------------------------------------------------------------------------------ 

!$acc kernels
      CDF = 0.0
      CDF2 = 0.0
      CDF3 = 0.0
      alpha = 0.0
!$acc end kernels

! ------ Apply overlap assumption --------

! generate the random numbers  

      ! Exponential overlap: weighting between maximum and random overlap
      !   increases with the distance. 
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
      do iplon = 1, ncol
         adl(iplon) = AM1+AM2*EXP(-(ALAT(iplon)*180./3.141592-AM3)**2/AM4**2)
         adl(iplon) = adl(iplon)*1.e3
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
      do iplon = 1, ncol
         rdl(iplon) = AM1+AM2*EXP(-(ALAT(iplon)*180./3.141592-AM3)**2/AM4**2)
         rdl(iplon) = rdl(iplon)*1.e3
      end do
!$acc end kernels      
      
!$acc kernels
      do iplon = 1, ncol
         alpha(iplon,1) = 0.
         do ilev = 2, nlay
            alpha(iplon,ilev) = exp( -(zm(iplon,ilev) - zm(iplon,ilev-1)) / adl(iplon))
         end do
      end do
!$acc end kernels
        
!$acc kernels
      do ilev = 1,nlay
         do iplon = 1, ncol
            seed1 = (play(iplon,1) - int(play(iplon,1)))  * 100000000 - ilev
            seed2 = (play(iplon,2) - int(play(iplon,2)))  * 100000000 + ilev
            seed3 = (play(iplon,3) - int(play(iplon,3)))  * 100000000 + ilev * 6.2
            seed4 = (play(iplon,4) - int(play(iplon,4)))  * 100000000           
            do isubcol = 1,nsubcol
               seed1 = 69069  * seed1 + 132721785 
               seed2 = 11002  * iand (seed2, 65535 ) + ishft (seed2, - 16 )
               seed3 = 18000  * iand (seed3, 65535 ) + ishft (seed3, - 16 )
               seed4 = 30903  * iand (seed4, 65535 ) + ishft (seed4, - 16 )
               kiss = seed1 + seed2 + ishft (seed3, 16 ) + seed4
               CDF(iplon,ilev,isubcol) = kiss*2.328306e-10  + 0.5 
               ! for testing purposes to generate determanistic sequence
               !CDF(iplon,ilev,isubcol) = real(MOD((isubcol * 37 + ilev * 71) , 10000))/10000.0
               
               seed1 = 69069  * seed1 + 132721785 
               seed2 = 11002  * iand (seed2, 65535 ) + ishft (seed2, - 16 )
               seed3 = 18000  * iand (seed3, 65535 ) + ishft (seed3, - 16 )
               seed4 = 30903  * iand (seed4, 65535 ) + ishft (seed4, - 16 )
               kiss = seed1 + seed2 + ishft (seed3, 16 ) + seed4
               CDF2(iplon,ilev,isubcol) = kiss*2.328306e-10  + 0.5 
               ! for testing purposes to generate determanistic sequence
               !CDF2(iplon,ilev,isubcol) =  real(MOD((isubcol * 37 + ilev * 73) , 10000))/10000.0
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
    
   if (inhm >= 1) then
        
!$acc kernels  
      do iplon = 1, ncol
         alpha(iplon,1)=0.
         do ilev = 2, nlay
             alpha(iplon, ilev) = exp( -( zm(iplon, ilev) - zm(iplon, ilev-1)) / rdl(iplon) )
         end do
      end do
!$acc end kernels
        
!$acc kernels
      do ilev = 1,nlay
         do iplon = 1, ncol
            seed1 = (play(iplon,1) - int(play(iplon,1)))  * 100000000 - ilev
            seed2 = (play(iplon,2) - int(play(iplon,2)))  * 100000000 + ilev
            seed3 = (play(iplon,3) - int(play(iplon,3)))  * 100000000 + ilev * 6.2
            seed4 = (play(iplon,4) - int(play(iplon,4)))  * 100000000           
            do isubcol = 1,nsubcol
               seed1 = 69069  * seed1 + 132721785 
               seed2 = 11002  * iand (seed2, 65535 ) + ishft (seed2, - 16 )
               seed3 = 18000  * iand (seed3, 65535 ) + ishft (seed3, - 16 )
               seed4 = 30903  * iand (seed4, 65535 ) + ishft (seed4, - 16 )
               kiss = seed1 + seed2 + ishft (seed3, 16 ) + seed4
               CDF2(iplon,ilev,isubcol) = kiss*2.328306e-10  + 0.5 
               ! for testing purposes to generate determanistic sequence
               !CDF2(iplon,ilev,isubcol) = real(MOD((isubcol * 7 + ilev * 101) , 10000))/10000.0
               
               seed1 = 69069  * seed1 + 132721785 
               seed2 = 11002  * iand (seed2, 65535 ) + ishft (seed2, - 16 )
               seed3 = 18000  * iand (seed3, 65535 ) + ishft (seed3, - 16 )
               seed4 = 30903  * iand (seed4, 65535 ) + ishft (seed4, - 16 )
               kiss = seed1 + seed2 + ishft (seed3, 16 ) + seed4
               CDF3(iplon,ilev,isubcol) = kiss*2.328306e-10  + 0.5 
               ! for testing purposes to generate determanistic sequence
               !CDF3(iplon,ilev,isubcol) =  real(MOD((isubcol * 41 + ilev * 89) , 10000))/10000.0
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

   end if  ! horiz variability



! where the subcolumn is cloudy, the subcolumn cloud fraction is 1;
! where the subcolumn is not cloudy, the subcolumn cloud fraction is 0;
! where there is a cloud, define the subcolumn cloud properties,
! otherwise set these to zero

      ngbm = ngb(1) - 1

!$acc kernels 
      do ilev = 1,nlay
         do iplon = 1, ncol
            do isubcol = 1, nsubcol
               if ( inhm >= 1 .and. CDF(iplon,ilev,isubcol) >= (1.0 - cld(iplon,ilev))) then
                  cld_stoch(iplon,ilev,isubcol) = 1.0 
                  
                  if (cld(iplon,ilev) .gt. 0.99 ) then
                    SIGMA_QCW = 0.5
                  elseif (cld(iplon,ilev) .gt. 0.9 ) then
                    SIGMA_QCW = 0.71
                  else  
                    SIGMA_QCW = 1
                  endif
                  
                  RIND1 = CDF3(iplon,ilev,isubcol) * (N1 - 1) + 1.0
                  IND1  = MAX(1, MIN(INT(RIND1), N1-1))
                  RIND1 = RIND1 - IND1
                  RIND2 = 40.0 * SIGMA_QCW - 3.0
                  IND2  = MAX(1, MIN(INT(RIND2), N2-1))
                  RIND2 = RIND2 - IND2

                  ZCW =  (1.0-RIND1) * (1.0-RIND2) * XCW(IND1,IND2) &
                       + (1.0-RIND1) * RIND2       * XCW(IND1,IND2+1) &
                       + RIND1 * (1.0-RIND2)       * XCW(IND1+1,IND2) &
                       + RIND1 * RIND2             * XCW(IND1+1,IND2+1)
                   
                  clwp_stoch(iplon,ilev,isubcol) = clwp(iplon,ilev) * ZCW
                  ciwp_stoch(iplon,ilev,isubcol) = ciwp(iplon,ilev) * ZCW

                  ibnd = ngb(isubcol) - ngbm
                  tauc_stoch(iplon,ilev,isubcol) = tauc(iplon,ilev,ibnd)
                  ssac_stoch(iplon,ilev,isubcol) = ssac(iplon,ilev,ibnd)
                  asmc_stoch(iplon,ilev,isubcol) = asmc(iplon,ilev,ibnd)
                  fsfc_stoch(iplon,ilev,isubcol) = fsfc(iplon,ilev,ibnd)
                
                
               elseif ( CDF(iplon,ilev,isubcol) >= (1.0 - cld(iplon,ilev)) ) then

                  cld_stoch(iplon,ilev,isubcol) = 1.0 
                  clwp_stoch(iplon,ilev,isubcol) = clwp(iplon,ilev)
                  ciwp_stoch(iplon,ilev,isubcol) = ciwp(iplon,ilev)

                  ibnd = ngb(isubcol) - ngbm
                  tauc_stoch(iplon,ilev,isubcol) = tauc(iplon,ilev,ibnd)
                  ssac_stoch(iplon,ilev,isubcol) = ssac(iplon,ilev,ibnd)
                  asmc_stoch(iplon,ilev,isubcol) = asmc(iplon,ilev,ibnd)
                  fsfc_stoch(iplon,ilev,isubcol) = fsfc(iplon,ilev,ibnd)

               else

                  cld_stoch(iplon,ilev,isubcol) = 0. 
                  clwp_stoch(iplon,ilev,isubcol) = 0. 
                  ciwp_stoch(iplon,ilev,isubcol) = 0. 

                  tauc_stoch(iplon,ilev,isubcol) = 0. 
                  ssac_stoch(iplon,ilev,isubcol) = 1. 
                  asmc_stoch(iplon,ilev,isubcol) = 0. 
                  fsfc_stoch(iplon,ilev,isubcol) = 0. 

               endif
            enddo
         enddo
      enddo
!$acc end kernels

      end subroutine mcica_sw


   end module mcica_subcol_gen_sw


