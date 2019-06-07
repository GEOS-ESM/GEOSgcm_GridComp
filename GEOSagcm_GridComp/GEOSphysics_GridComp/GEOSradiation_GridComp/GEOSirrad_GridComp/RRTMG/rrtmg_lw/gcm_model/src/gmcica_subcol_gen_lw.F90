!     path:      $Source$
!     author:    $Author$
!     revision:  $Revision$
!     created:   $Date$
!
#include "_gpudef.inc"
    
      module gpu_mcica_subcol_gen_lw

!  --------------------------------------------------------------------------
! |                                                                          |
! |  Copyright 2006-2009, Atmospheric & Environmental Research, Inc. (AER).  |
! |  This software may be used, copied, or redistributed as long as it is    |
! |  not sold and this copyright notice is reproduced on each copy made.     |
! |  This model is provided as is without any express or implied warranties. |
! |                       (http://www.rtweb.aer.com/)                        |
! |                                                                          |
!  --------------------------------------------------------------------------

! Purpose: Create McICA stochastic arrays for cloud physical or optical properties.
! Two options are possible:
! 1) Input cloud physical properties: cloud fraction, ice and liquid water
!    paths, ice fraction, and particle sizes.  Output will be stochastic
!    arrays of these variables.  (inflag = 1)
! 2) Input cloud optical properties directly: cloud optical depth, single
!    scattering albedo and asymmetry parameter.  Output will be stochastic
!    arrays of these variables.  (inflag = 0; longwave scattering is not
!    yet available, ssac and asmc are for future expansion)

! --------- Modules ----------

       !use parkind, only : im => kind , rb => kind 
      use parrrtm, only : nbndlw, ngptlw
      use rrlw_con, only: grav
      use rrlw_wvn, only: ngb
      use rrlw_vsn
      use WaterDistributionMod, only: n1, n2, tabulate_xcw_beta, tabulate_xcw_gamma

#ifdef _CUDA
      use cudafor
      use cudadevice
#endif

      implicit none

      real , dimension(n1,n2) :: XCW

      real  _gpudev, allocatable :: pmidd(:, :),xcwd(:,:)
      real  _gpudev, allocatable :: cldfracd(:,:), clwpd(:,:), ciwpd(:,:), taucd(:,:,:)
      real  _gpudev, allocatable :: alphad(:,:), rcorrd(:,:)
      integer  _gpudev, allocatable :: cloudFlagd(:,:)


! public interfaces/functions/subroutines
      !public :: mcica_subcol_lwg, generate_stochastic_cloudsg 

      contains

!------------------------------------------------------------------
! Public subroutines
!------------------------------------------------------------------

      subroutine mcica_subcol_lwg(iplon, ncol, nlay, icld, permuteseed, irng, play, &
                       cldfrac, ciwp, clwp, tauc, ngbd, cldfmcl, &
                       ciwpmcl, clwpmcl,taucmcl, cloudFlag,cloudMH, cloudHH, zmd, alatd, doy, inhm)

! ----- Input -----
! Control
      integer , intent(in) :: iplon           ! column/longitude index
      integer , intent(in) :: ncol            ! number of columns
      integer , intent(in) :: nlay            ! number of model layers
      integer , intent(in) :: icld            ! clear/cloud, cloud overlap flag
      integer , intent(in) :: permuteseed     ! if the cloud generator is called multiple times, 
                                                      ! permute the seed between each call.
                                                      ! between calls for LW and SW, recommended
                                                      ! permuteseed differes by 'ngpt'
      integer , intent(in) :: irng         ! flag for random number generator
                                                      !  0 = kissvec
                                                      !  1 = Mersenne Twister
      integer , intent(in) :: cloudMH, cloudHH

! Atmosphere
      real , intent(in) :: play(:,:)          ! layer pressures (mb) 
                                                      !    Dimensions: (ncol,nlay)

! Atmosphere/clouds - cldprop
      real , intent(in) :: cldfrac(:,:)       ! layer cloud fraction
                                                      !    Dimensions: (ncol,nlay)
      real , intent(in) :: tauc(:,:,:)        ! in-cloud optical depth
                                                      !    Dimensions: (nbndlw,ncol,nlay)
      real , intent(in) :: ciwp(:,:)          ! in-cloud ice water path
                                                      !    Dimensions: (ncol,nlay)
      real , intent(in) :: clwp(:,:)          ! in-cloud liquid water path
                                                      !    Dimensions: (ncol,nlay)
      integer  _gpudev, intent(in) :: ngbd(:)
! ----- Output -----
! Atmosphere/clouds - cldprmc [mcica]
      real  _gpudev, intent(out) :: cldfmcl(:,:,:)    ! cloud fraction [mcica]
                                                      !    Dimensions: (ngptlw,ncol,nlay)
      real  _gpudev, intent(out) :: ciwpmcl(:,:,:)    ! in-cloud ice water path [mcica]
                                                      !    Dimensions: (ngptlw,ncol,nlay)
      real  _gpudev, intent(out) :: clwpmcl(:,:,:)    ! in-cloud liquid water path [mcica]
                                                      !    Dimensions: (ngptlw,ncol,nlay)
      real  _gpudev, intent(out) :: taucmcl(:,:,:)    ! in-cloud optical depth [mcica]
                                                      !    Dimensions: (ngptlw,ncol,nlay)
      integer , intent(out) :: cloudFlag(:,:)
      

      real  _gpudev, intent(in) :: zmd(:,:)  ! Height of midpoints (above surface)
                                                      !    Dimensions: (ncol,nlay)
      real  _gpudev, intent(in) :: alatd(:)  ! Latitude of column
                                                      !    Dimensions: (ncol)
      integer , intent(in) :: doy           ! Day of year
      integer , intent(in)  :: inhm

! ----- Local -----

! Stochastic cloud generator variables [mcica]
      integer , parameter :: nsubclw = ngptlw ! number of sub-columns (g-point intervals)
      integer  :: ilev                        ! loop index

      real  :: pmid(ncol, nlay)               ! layer pressures (Pa) 
      real  :: alpha(ncol, nlay), rcorr(ncol, nlay)

      real  :: rdl(ncol),adl(ncol)
      real  :: am1,am2,am3,am4,amr
#ifdef _CUDA
      type(dim3) :: dimGrid, dimBlock
#endif
      integer, save :: counter = 0
      !real , parameter  :: Zo = 2000.         ! length scale (m) 
      !real , parameter  :: Zc = 1000.         ! length scale (m) 
      integer :: i,j,k,tk
      real :: t1, t2
  
! Return if clear sky; or stop if icld out of range
      if (icld.eq.0) then 
        cldfmcl = 0.0
        ciwpmcl = 0.0
        clwpmcl = 0.0
        taucmcl = 0.0
        cloudFlag = 0.0

        return
      end if 
      if (icld.lt.0.or.icld.gt.4) then 
         stop 'MCICA_SUBCOL: INVALID ICLD'
      endif 
   
! NOTE: For GCM mode, permuteseed must be offset between LW and SW by at least the number of subcolumns


! Pass particle sizes to new arrays, no subcolumns for these properties yet
! Convert pressures from mb to Pa

       pmid(:ncol,:nlay) = play(:ncol,:nlay)*1.e2 

           if     (inhm == 1) then
            call tabulate_xcw_beta(xcw)
         elseif (inhm == 2) then
            call tabulate_xcw_gamma(xcw)
         endif

      allocate( pmidd(ncol, nlay), cldfracd(ncol, nlay+1))
      allocate( clwpd(ncol, nlay+1), ciwpd(ncol, nlay+1), taucd(ncol, nbndlw, nlay))
      allocate( xcwd( n1, n2))
      allocate( alphad( ncol, nlay), rcorrd(ncol, nlay))
      allocate( cloudFlagd( ncol, 4))

      ! compute alpha and rcorr

      !$acc kernels
      am1 = 1.4315
      am2 = 2.1219
      am4 = -25.584
      amr = 7.
      if ( doy .gt. 181 ) then
         am3 = -4.*amr/365*(doy-272)
      else
         am3 = 4.*amr/365.*(doy-91)
      endif
      adl(:) = am1+am2*exp(-(alatd(:)*180./3.141592-am3)**2/am4**2)
      adl(:) = adl(:)*1.e3
      am1 = 0.72192
      am2 = 0.78996
      am4 = 40.404
      amr = 8.5
      if ( doy .gt. 181 ) then
         am3 = -4.*amr/365*(doy-272)
      else
         am3 = 4.*amr/365.*(doy-91)
      endif
      rdl(:) = am1+am2*exp(-(alatd(:)*180./3.141592-am3)**2/am4**2)
      rdl(:) = rdl(:)*1.e3
      !$acc end kernels      

        
      !$acc kernels 
      do j = 1, ncol
        alphad(j, 1) = 0.0
        rcorrd(j, 1) = 0.0
        do i = 2, nlay
             alphad(j,i) = exp( -( zmd (j, i) -  zmd (j, i-1)) / adl(j))
             rcorrd(j,i) = exp( -( zmd (j, i) -  zmd (j, i-1)) / rdl(j))
         end do
        end do
      !$acc end kernels


      pmidd = pmid
    
      cldfracd = cldfrac
      clwpd = clwp
      ciwpd = ciwp
      taucd = tauc
      xcwd = xcw


 end subroutine mcica_subcol_lwg

 subroutine mcica_cldc(ncol, icld, nlay, cloudmh, cloudhh, cloudflag, cldfmcl)


     integer :: i,j,k,tk
     integer , parameter :: nsubclw = ngptlw ! number of sub-columns (g-point intervals)
     integer , intent(in) :: ncol            ! number of columns
     integer , intent(in) :: nlay            ! number of model layers
     integer , intent(in) :: cloudMH, cloudHH
     integer , intent(out) :: cloudFlag(:,:)
     integer , intent(in) :: icld
     real  _gpudev, intent(in) :: cldfmcl(:,:,:)    ! cloud fraction [mcica]
     if (icld.eq.0) return                                             


     !$acc kernels
     do i = 1, ncol
        cloudFlagd(i,1)=0        
        cloudFlagd(i,2)=0
        cloudFlagd(i,3)=0
        cloudFlagd(i,4)=0
       
        do j = 1, nsubclw
            tk=0
            do k = 1, nlay
                 if (cldfmcl(i,j,k) > 0.5) then 
                    tk=1
                 end if
            end do
            if (tk .eq. 0) then 
                cloudFlagd(i,1) = cloudFlagd(i,1)+1
            end if
            tk=0
            do k = cloudHH, nlay
                 if (cldfmcl(i,j,k) > 0.5) then 
                    tk=1
                 end if
            end do
            if (tk .eq. 0) then 
                cloudFlagd(i,2) = cloudFlagd(i,2)+1
            end if
              tk=0
            do k = cloudMH, cloudHH+1
                 if (cldfmcl(i,j,k) > 0.5) then 
                    tk=1
                 end if
            end do
            if (tk .eq. 0) then 
                cloudFlagd(i,3) = cloudFlagd(i,3)+1
            end if
              tk=0
            do k = 1, cloudMH+1
                 if (cldfmcl(i,j,k) > 0.5) then 
                    tk=1
                 end if
            end do
            if (tk .eq. 0) then 
                cloudFlagd(i,4) = cloudFlagd(i,4)+1
            end if 
        end do
     end do
     !$acc end kernels
   
     cloudFlag = cloudFlagd
     
     deallocate( pmidd, cldfracd)
     deallocate( clwpd, ciwpd, taucd, xcwd, cloudFlagd)
     deallocate( alphad, rcorrd)
  
    

     end subroutine mcica_cldc


!-------------------------------------------------------------------------------------------------
       _gpuker subroutine generate_stochastic_cloudsg(ncol, nlay,icld, &
                                   ngbd, cld_stoch, clwp_stoch, ciwp_stoch, tauc_stoch, ih, zm) 
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


 
! -- Arguments

      integer , intent(in)  :: ncol            ! number of columns
      integer , intent(in):: nlay            ! number of layers
      integer , intent(in) :: icld            ! clear/cloud, cloud overlap flag
  
       integer  _gpudev, intent(in) :: ngbd(:)
!      real , intent(in) :: ssac(:,:,:)       ! in-cloud single scattering albedo
                                                      !    Dimensions: (nbndlw,ncol,nlay)
                                                      !   inactive - for future expansion
!      real , intent(in) :: asmc(:,:,:)       ! in-cloud asymmetry parameter
                                                      !    Dimensions: (nbndlw,ncol,nlay)
                                                      !   inactive - for future expansion

      real  _gpudev, intent(out) :: cld_stoch(:,:,:)  ! subcolumn cloud fraction 
                                                      !    Dimensions: (ngptlw,ncol,nlay)
      real  _gpudev, intent(out) :: clwp_stoch(:,:,:) ! subcolumn in-cloud liquid water path
                                                      !    Dimensions: (ngptlw,ncol,nlay)
      real  _gpudev, intent(out) :: ciwp_stoch(:,:,:) ! subcolumn in-cloud ice water path
                                                      !    Dimensions: (ngptlw,ncol,nlay)
      real  _gpudev, intent(out) :: tauc_stoch(:,:,:) ! subcolumn in-cloud optical depth
                                                      !    Dimensions: (ngptlw,ncol,nlay)
!      real , intent(out) :: ssac_stoch(:,:,:)! subcolumn in-cloud single scattering albedo
                                                      !    Dimensions: (ngptlw,ncol,nlay)
                                                      !   inactive - for future expansion
!      real , intent(out) :: asmc_stoch(:,:,:)! subcolumn in-cloud asymmetry parameter
                                                      !    Dimensions: (ngptlw,ncol,nlay)
                                                      !   inactive - for future expansion
      real  _gpudev, intent(in) :: zm(:,:)  ! Height of midpoints (above surface)
                                                      !    Dimensions: (ncol,nlay)
     
      integer , value, intent(in) :: ih      
      !integer, value, intent(in) :: counter
   


!GPU These are needed to allocate space for temporary automatic local
!GPU arrays on the GPU. On the CPU this ifndef converts them to nlay

#ifndef GPU_MAXLEVS
#define GPU_MAXLEVS nlay
#endif

      
     
       
! Cloud condensate
      
       real  :: RIND1, RIND2, ZCW, SIGMA_QCW
       integer  :: IND1, IND2
     
       real  :: CDF3(GPU_MAXLEVS)      ! random numbers

       integer  :: inhm_ 
       real  :: cfs
       integer, parameter :: nsubcol = 140
       
! Constants (min value for cloud fraction and cloud water and ice)
     ! real , parameter :: cldmin = 1.0e-20  ! min cloud fraction
!      real , parameter :: qmin   = 1.0e-10    ! min cloud water and cloud ice (not used)

! Variables related to random number and seed 
      real  :: CDF(GPU_MAXLEVS), CDF2(GPU_MAXLEVS)      ! random numbers
      integer  :: seed1, seed2, seed3, seed4 ! seed to create random number (kissvec)
      real  :: rand_num      ! random number (kissvec)
      integer  :: iseed                       ! seed to create random number (Mersenne Teister)
      real  :: rand_num_mt                    ! random number (Mersenne Twister)

! Flag to identify cloud fraction in subcolumns
   !   logical :: iscloudy(GPU_MAXLEVS)   ! flag that says whether a gridbox is cloudy

! Indices
      integer  :: ilev, isubcol, i, n         ! indices
      
      integer :: iplon, gp
#ifdef _CUDA
      iplon = (blockidx%x-1) * blockdim%x + threadidx%x
      gp = (blockidx%y-1) * blockdim%y + threadidx%y

!------------------------------------------------------------------------------------------ 
     !   print *, "ppp ", iplon, gp
      if (iplon <= ncol .and. gp <= nsubcol) then
#else
      do iplon = 1, ncol
      do gp = 1, nsubcol
#endif



! ----- Create seed  --------
   
! Advance randum number generator by changeseed values
   
! For kissvec, create a seed that depends on the state of the columns. Maybe not the best way, but it works.  
! Must use pmid from bottom four layers. 
      
           seed1 = (pmidd(iplon,1) - int(pmidd(iplon,1)))  * 1000000000 + (gp) * 11 
           seed3 = (pmidd(iplon,3) - int(pmidd(iplon,3)))  * 1000000000 + (gp) * 13 
           seed2 = seed1 + gp
           seed4 = seed3 - gp 
  

! ------ Apply overlap assumption --------

! generate the random numbers  



       select case (icld)

       case(1) 

           do ilev = 1,nlay
             call kissvec(seed1, seed2, seed3, seed4, rand_num)
             CDF(ilev) = rand_num
           end do
        

       case(2)

           do ilev = 1,nlay
             call kissvec(seed1, seed2, seed3, seed4, rand_num)
             CDF(ilev) = rand_num
           end do
          

           do ilev = 2,nlay
             if (CDF(ilev-1) > 1.  - cldfracd(iplon, ilev-1)) then 
                CDF(ilev) = CDF(ilev-1)
             else
                 CDF(ilev) = CDF(ilev) * (1. - cldfracd(iplon, ilev-1))
             end if
           end do
            
       case(3)

           call kissvec(seed1, seed2, seed3, seed4, rand_num)
           do ilev = 1,nlay
            CDF(ilev) = rand_num
           end do


       case(4)

           do ilev = 1,nlay
              call kissvec(seed1, seed2, seed3, seed4, rand_num)
              CDF(ilev) = rand_num
              call kissvec(seed1, seed2, seed3, seed4, rand_num)
              CDF2(ilev) = rand_num
              call kissvec(seed1, seed2, seed3, seed4, rand_num)
              CDF3(ilev) = rand_num 
           end do

        do ilev = 2,nlay

    

           if (CDF2(ilev) < alphad(iplon,ilev)) then
            CDF(ilev) = CDF(ilev-1)
           end if
           
        end do

        end select 


        if (ih > 0) then

        do ilev = 1,nlay
              call kissvec(seed1, seed2, seed3, seed4, rand_num)
              CDF2(ilev) = rand_num
              call kissvec(seed1, seed2, seed3, seed4, rand_num)
              CDF3(ilev) = rand_num 
        end do


        do ilev = 2,nlay
           if (CDF2(ilev) < rcorrd (iplon,ilev)) then
               CDF3(ilev) = CDF3(ilev-1) 
            end if
        end do

        end if

     
   
      
         
! -- generate subcolumns for homogeneous clouds -----
     ! do ilev = 1,nlay
     !    iscloudy(ilev) = (CDF(ilev) >= 1.  - cldfracd(iplon,ilev))
     ! enddo

! where the subcolumn is cloudy, the subcolumn cloud fraction is 1;
! where the subcolumn is not cloudy, the subcolumn cloud fraction is 0;
! where there is a cloud, define the subcolumn cloud properties, 
! otherwise set these to zero
          n = ngbd(gp)
      do ilev = 1,nlay
        cfs = cldfracd(iplon, ilev)
         !  do gp = 1, nsubcol
               if (CDF(ilev) >=1.  - cfs) then
                  cld_stoch(iplon,gp,ilev) = 1. 
                  clwp_stoch(iplon,gp,ilev) = clwpd(iplon,ilev)
                  ciwp_stoch(iplon,gp,ilev) = ciwpd(iplon,ilev)
                
                  tauc_stoch(iplon,gp,ilev) = taucd(iplon,n,ilev)
                  
                   if (ih > 0) then  
                     ! Cloud fraction ~ inhomogeneity
                     if (cfs .gt. 0.99 ) then
                        SIGMA_QCW = 0.5
                     elseif (cfs .gt. 0.9 ) then
                        SIGMA_QCW = 0.71
                     else  
                        SIGMA_QCW = 1
                     endif

                     cld_stoch(iplon,gp,ilev) = 1. 

                     ! Horizontally veriable clouds:
                     ! Determine ZCW = ratio of cloud condensate miximg ratio QC for this cell to
                     ! its mean value for all cloudy cells in this layer.
                     ! Use bilinear interpolation of ZCW tabulated in array XCW as a function of
                     !    * cumulative probability Y (= CDF3)
                     !    * relative standard deviation SIGMA
                     ! Take care that the definition of RIND2 is consistent with subroutine
                     ! TABULATE_XCW
                     RIND1 = CDF3(ilev) * (N1 - 1) + 1.0
                     IND1  = MAX(1, MIN(INT(RIND1), N1-1))
                     RIND1 = RIND1 - IND1
                     RIND2 = 40.0 * SIGMA_QCW    - 3.0
                     IND2  = MAX(1, MIN(INT(RIND2), N2-1))
                     RIND2 = RIND2 - IND2

                     ZCW = (1.0-RIND1) * (1.0-RIND2) * XCWD(IND1,IND2) &
                          + (1.0-RIND1) * RIND2       * XCWD(IND1,IND2+1) &
                          + RIND1 * (1.0-RIND2)       * XCWD(IND1+1,IND2) &
                          + RIND1 * RIND2             * XCWD(IND1+1,IND2+1)

                     clwp_stoch(iplon,gp,ilev) = clwpd(iplon,ilev) * ZCW
                     ciwp_stoch(iplon,gp,ilev) = ciwpd(iplon,ilev) * ZCW
                  
                     tauc_stoch(iplon,gp,ilev) = taucd(iplon,n,ilev)
!                  ssac_stoch(isubcol,i,ilev) = ssac(n,i,ilev)
!                  asmc_stoch(isubcol,i,ilev) = asmc(n,i,ilev)
                   end if
               else
                  cld_stoch(iplon,gp,ilev) = 0. 
                  clwp_stoch(iplon,gp,ilev) = 0. 
                  ciwp_stoch(iplon,gp,ilev) = 0. 
                  tauc_stoch(iplon,gp,ilev) = 0. 
!                  ssac_stoch(isubcol,i,ilev) = 1. 
!                  asmc_stoch(isubcol,i,ilev) = 1. 
               endif
           
       !end do
      enddo

#ifdef _CUDA
      endif
#else
      end do
      end do
#endif


      end subroutine generate_stochastic_cloudsg

      _gpuked  subroutine kissvec(seed1,seed2,seed3,seed4,ran_arr)
!-------------------------------------------------------------------------------------------------- 

! public domain code
! made available from http://www.fortran.com/
! downloaded by pjr on 03/16/04 for NCAR CAM
! converted to vector form, functions inlined by pjr,mvr on 05/10/2004

! The  KISS (Keep It Simple Stupid) random number generator. Combines:
! (1) The congruential generator x(n)=69069*x(n-1)+1327217885, period 2^32.
! (2) A 3-shift shift-register generator, period 2^32-1,
! (3) Two 16-bit multiply-with-carry generators, period 597273182964842497>2^59
!  Overall period>2^123; 
!
      real , intent(inout)  :: ran_arr
      integer , intent(inout) :: seed1,seed2,seed3,seed4
      integer  :: i,sz,kiss
      integer  :: m, k, n

! inline function 
      m(k, n) = ieor (k, ishft (k, n) )

     
     
         seed1 = 69069 * seed1 + 1327217885
         seed2 = m (m (m (seed2, 13), - 17), 5)
         seed3 = 18000 * iand (seed3, 65535) + ishft (seed3, - 16)
         seed4 = 30903 * iand (seed4, 65535) + ishft (seed4, - 16)
         kiss = seed1 + seed2 + ishft (seed3, 16) + seed4
         ran_arr = kiss*2.328306e-10  + 0.5 
     
    
      end subroutine kissvec


      end module gpu_mcica_subcol_gen_lw
