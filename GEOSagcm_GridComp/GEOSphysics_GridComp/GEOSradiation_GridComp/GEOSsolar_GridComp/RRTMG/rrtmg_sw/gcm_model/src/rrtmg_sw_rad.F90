!  --------------------------------------------------------------------------
! |                                                                          |
! |  Copyright 2002-2009, Atmospheric & Environmental Research, Inc. (AER).  |
! |  This software may be used, copied, or redistributed as long as it is    |
! |  not sold and this copyright notice is reproduced on each copy made.     |
! |  This model is provided as is without any express or implied warranties. |
! |                       (http://www.rtweb.aer.com/)                        |
! |                                                                          |
!  --------------------------------------------------------------------------
!
! ****************************************************************************
! *                                                                          *
! *                             RRTMG_SW                                     *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                 a rapid radiative transfer model                         *
! *                  for the solar spectral region                           *
! *           for application to general circulation models                  *
! *                                                                          *
! *                                                                          *
! *           Atmospheric and Environmental Research, Inc.                   *
! *                       131 Hartwell Avenue                                *
! *                       Lexington, MA 02421                                *
! *                                                                          *
! *                                                                          *
! *                          Eli J. Mlawer                                   *
! *                       Jennifer S. Delamere                               *
! *                        Michael J. Iacono                                 *
! *                        Shepard A. Clough                                 *
! *                       David M. Berthiaume                                *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                      email:  miacono@aer.com                             *
! *                      email:  emlawer@aer.com                             *
! *                      email:  jdelamer@aer.com                            *
! *                                                                          *
! *       The authors wish to acknowledge the contributions of the           *
! *       following people:  Steven J. Taubman, Patrick D. Brown,            *
! *       Ronald E. Farren, Luke Chen, Robert Bergstrom.                     *
! *                                                                          *
! ****************************************************************************

    
    
#ifdef _CUDA
#define gpu_device ,device
#else
#define gpu_device 
#endif
    
      module rrtmg_sw_rad

! --------- Modules ---------

      use rrsw_vsn
      use mcica_subcol_gen_sw, only: mcica_sw
      use rrtmg_sw_cldprmc, only: cldprmc_sw
      use rrtmg_sw_setcoef, only: setcoef_sw
      use rrtmg_sw_spcvmc, only: spcvmc_sw

      implicit none


      public :: rrtmg_sw,  earth_sun


    contains

      subroutine rrtmg_sw &
            (rpart, ncol    ,nlay    ,icld    , iaer, &
             play    ,plev    ,tlay    ,tlev    ,tsfc   , &
             h2ovmr , o3vmr   ,co2vmr  ,ch4vmr  ,n2ovmr ,o2vmr , &
             asdir   ,asdif   ,aldir   ,aldif   , &
             coszen  ,adjes   ,dyofyr  ,scon    ,isolvar, &
             inflgsw ,iceflgsw,liqflgsw,cld, &
             tauc ,ssac ,asmc ,fsfc , &
             ciwp ,clwp ,rei ,rel , &
             tauaer  ,ssaaer  ,asmaer  ,ecaer   , &
             swuflx  ,swdflx  ,swhr    ,swuflxc ,swdflxc ,swhrc, &
             nirr    ,nirf    ,parr    ,parf    ,uvrr    ,uvrf , normFlx, &
             zm, alat, numCPUs ,&
! optional I/O
             bndsolvar,indsolvar,solcycfrac)




      use parrrsw, only : nbndsw, ngptsw, naerec, nstr, nmol, mxmol, &
                          jpband, jpb1, jpb2, rrsw_scon
      use rrsw_aer, only : rsrtaua, rsrpiza, rsrasya
      use rrsw_con, only : heatfac, oneminus, pi,  grav, avogad
      use rrsw_wvn, only : wavenum1, wavenum2
      use rrsw_cld, only : extliq1, ssaliq1, asyliq1, &
                           extice2, ssaice2, asyice2, &
                           extice3, ssaice3, asyice3, fdlice3, &
                           extice4, ssaice4, asyice4, &
                           abari, bbari, cbari, dbari, ebari, fbari
      use rrsw_wvn, only : wavenum2, ngb
      use rrsw_ref, only : preflog, tref
      use iso_fortran_env, only : error_unit
#ifdef _CUDA
      use cudafor
#endif 

      

! ------- Declarations

      integer , intent(in) :: rpart
      integer, intent(in) :: ncol            ! Number of horizontal columns     
      integer, intent(in) :: nlay            ! Number of model layers
      integer, intent(inout) :: icld         ! Cloud overlap method
                                             !    0: Clear only
                                             !    1: Random
                                             !    2: Maximum/random
                                             !    3: Maximum
      integer , intent(in) :: iaer
      real, intent(in) :: play(:,:)          ! Layer pressures (hPa, mb)
                                             !    Dimensions: (ncol,nlay)
      real, intent(in) :: plev(:,:)          ! Interface pressures (hPa, mb)
                                             !    Dimensions: (ncol,nlay+1)
      real, intent(in) :: tlay(:,:)          ! Layer temperatures (K)
                                             !    Dimensions: (ncol,nlay)
      real, intent(in) :: tlev(:,:)          ! Interface temperatures (K)
                                             !    Dimensions: (ncol,nlay+1)
      real, intent(in) :: tsfc(:)            ! Surface temperature (K)
                                             !    Dimensions: (ncol)
      real, intent(in) :: h2ovmr(:,:)        ! H2O volume mixing ratio
                                             !    Dimensions: (ncol,nlay)
      real, intent(in) :: o3vmr(:,:)         ! O3 volume mixing ratio
                                             !    Dimensions: (ncol,nlay)
      real, intent(in) :: co2vmr(:,:)        ! CO2 volume mixing ratio
                                             !    Dimensions: (ncol,nlay)
      real, intent(in) :: ch4vmr(:,:)        ! Methane volume mixing ratio
                                             !    Dimensions: (ncol,nlay)
      real, intent(in) :: n2ovmr(:,:)        ! Nitrous oxide volume mixing ratio
                                             !    Dimensions: (ncol,nlay)
      real, intent(in) :: o2vmr(:,:)         ! Oxygen volume mixing ratio
                                             !    Dimensions: (ncol,nlay)
      real, intent(in) :: asdir(:)           ! UV/vis surface albedo direct rad
                                             !    Dimensions: (ncol)
      real, intent(in) :: aldir(:)           ! Near-IR surface albedo direct rad
                                             !    Dimensions: (ncol)
      real, intent(in) :: asdif(:)           ! UV/vis surface albedo: diffuse rad
                                             !    Dimensions: (ncol)
      real, intent(in) :: aldif(:)           ! Near-IR surface albedo: diffuse rad
                                             !    Dimensions: (ncol)

      integer, intent(in) :: dyofyr          ! Day of the year (used to get Earth/Sun
                                             !  distance if adjflx not provided)
      real, intent(in) :: adjes              ! Flux adjustment for Earth/Sun distance
      real, intent(in) :: coszen(:)          ! Cosine of solar zenith angle
                                             !    Dimensions: (ncol)
      real, intent(in) :: scon               ! Solar constant (W/m2)
                                             !    Total solar irradiance averaged 
                                             !    over the solar cycle.
                                             !    If scon = 0.0, the internal solar 
                                             !    constant, which depends on the  
                                             !    value of isolvar, will be used. 
                                             !    For isolvar=-1, scon=1368.22 Wm-2,
                                             !    For isolvar=0,1,3, scon=1360.85 Wm-2,
                                             !    If scon > 0.0, the internal solar
                                             !    constant will be scaled to the 
                                             !    provided value of scon.
      integer, intent(in) :: isolvar         ! Flag for solar variability method
                                             !   -1 = (when scon .eq. 0.0): No solar variability
                                             !        and no solar cycle (Kurucz solar irradiance
                                             !        of 1368.22 Wm-2 only);
                                             !        (when scon .ne. 0.0): Kurucz solar irradiance
                                             !        scaled to scon and solar variability defined
                                             !        (optional) by setting non-zero scale factors
                                             !        for each band in bndsolvar
                                             !    0 = (when SCON .eq. 0.0): No solar variability 
                                             !        and no solar cycle (NRLSSI2 solar constant of 
                                             !        1360.85 Wm-2 for the 100-50000 cm-1 spectral 
                                             !        range only), with facular and sunspot effects 
                                             !        fixed to the mean of Solar Cycles 13-24;
                                             !        (when SCON .ne. 0.0): No solar variability 
                                             !        and no solar cycle (NRLSSI2 solar constant of 
                                             !        1360.85 Wm-2 for the 100-50000 cm-1 spectral 
                                             !        range only), is scaled to SCON
                                             !    1 = Solar variability (using NRLSSI2  solar
                                             !        model) with solar cycle contribution
                                             !        determined by fraction of solar cycle
                                             !        with facular and sunspot variations
                                             !        fixed to their mean variations over the
                                             !        average of Solar Cycles 13-24;
                                             !        two amplitude scale factors allow
                                             !        facular and sunspot adjustments from
                                             !        mean solar cycle as defined by indsolvar 
                                             !    2 = Solar variability (using NRLSSI2 solar
                                             !        model) over solar cycle determined by 
                                             !        direct specification of Mg (facular)
                                             !        and SB (sunspot) indices provided
                                             !        in indsolvar (scon = 0.0 only)
                                             !    3 = (when scon .eq. 0.0): No solar variability
                                             !        and no solar cycle (NRLSSI2 solar irradiance
                                             !        of 1360.85 Wm-2 only);
                                             !        (when scon .ne. 0.0): NRLSSI2 solar irradiance
                                             !        scaled to scon and solar variability defined
                                             !        (optional) by setting non-zero scale factors
                                             !        for each band in bndsolvar
      real, intent(in), optional :: indsolvar(:) ! Facular and sunspot amplitude 
                                                 ! scale factors (isolvar=1), or
                                                 ! Mg and SB indices (isolvar=2)
                                                 !    Dimensions: (2)
      real, intent(in), optional :: bndsolvar(:) ! Solar variability scale factors 
                                                 ! for each shortwave band
                                                 !    Dimensions: (nbndsw=14)
      real, intent(in), optional :: solcycfrac   ! Fraction of averaged 11-year solar cycle (0-1)
                                                 !    at current time (isolvar=1)
                                                 !    0.0 represents the first day of year 1
                                                 !    1.0 represents the last day of year 11

      integer, intent(in) :: inflgsw         ! Flag for cloud optical properties
      integer, intent(in) :: iceflgsw        ! Flag for ice particle specification
      integer, intent(in) :: liqflgsw        ! Flag for liquid droplet specification

      real , intent(in) :: cld(:,:)           ! Cloud fraction
                                              !    Dimensions: (ncol,nlay)
      real , intent(in) :: tauc(:,:,:)        ! In-cloud optical depth
                                              !    Dimensions: (ncol,nlay,nbndsw)
      real , intent(in) :: ssac(:,:,:)        ! In-cloud single scattering albedo
                                              !    Dimensions: (ncol,nlay,nbndsw)
      real , intent(in) :: asmc(:,:,:)        ! In-cloud asymmetry parameter
                                              !    Dimensions: (ncol,nlay,nbndsw)
      real , intent(in) :: fsfc(:,:,:)        ! In-cloud forward scattering fraction
                                              !    Dimensions: (ncol,nlay,nbndsw)
      real , intent(in) :: ciwp(:,:)          ! In-cloud ice water path (g/m2)
                                              !    Dimensions: (ncol, nlay)
      real , intent(in) :: clwp(:,:)          ! In-cloud liquid water path (g/m2)
                                              !    Dimensions: (ncol, nlay)
      real , intent(in) :: rei(:,:)           ! Cloud ice effective radius (microns)
                                              !    Dimensions: (ncol, nlay)

      real , intent(in) :: rel(:,:)           ! Cloud water drop effective radius (microns)
                                              !    Dimensions: (ncol,nlay)
      real, intent(in) :: tauaer(:,:,:)       ! Aerosol optical depth (iaer=10 only)
                                              !    Dimensions: (ncol,nlay,nbndsw)
                                              ! (non-delta scaled)      
      real, intent(in) :: ssaaer(:,:,:)       ! Aerosol single scattering albedo (iaer=10 only)
                                              !    Dimensions: (ncol,nlay,nbndsw)
                                              ! (non-delta scaled)      
      real, intent(in) :: asmaer(:,:,:)       ! Aerosol asymmetry parameter (iaer=10 only)
                                              !    Dimensions: (ncol,nlay,nbndsw)
                                              ! (non-delta scaled)      
      real, intent(in) :: ecaer(:,:,:)        ! Aerosol optical depth at 0.55 micron (iaer=6 only)
                                              !    Dimensions: (ncol,nlay,naerec)
                                              ! (non-delta scaled)      
      integer , intent(in) :: normFlx         ! Normalize fluxes flag
                                              !  0 = no normalization
                                              !  1 = normalize fluxes ( / (scon * coszen) )
      real , intent(in) :: zm(:,:)            ! Heights of level midpoints
                                              !    Dimensions: (ncol,nlay)
      real , intent(in) :: alat(:)            ! Latitude of column
                                              !    Dimensions: (ncol)
      integer , intent(in) :: numCPUs         ! Number of cores per node
                                              

! ----- Output -----

      real, intent(out) :: swuflx(:,:)       ! Total sky shortwave upward flux (W/m2)
                                             !    Dimensions: (ncol,nlay+1)
      real, intent(out) :: swdflx(:,:)       ! Total sky shortwave downward flux (W/m2)
                                             !    Dimensions: (ncol,nlay+1)
      real, intent(out) :: swhr(:,:)         ! Total sky shortwave radiative heating rate (K/d)
                                             !    Dimensions: (ncol,nlay)
      real, intent(out) :: swuflxc(:,:)      ! Clear sky shortwave upward flux (W/m2)
                                             !    Dimensions: (ncol,nlay+1)
      real, intent(out) :: swdflxc(:,:)      ! Clear sky shortwave downward flux (W/m2)
                                             !    Dimensions: (ncol,nlay+1)
      real, intent(out) :: swhrc(:,:)        ! Clear sky shortwave radiative heating rate (K/d)
                                             !    Dimensions: (ncol,nlay)
   ! Output added for Land/Surface process
      real , intent(out) :: nirr(:)           ! Near-IR direct downward shortwave flux (w/m2)
                                              !    Dimensions: (ncol)
      real , intent(out) :: nirf(:)           ! Near-IR diffuse downward shortwave flux (w/m2)
                                              !    Dimensions: (ncol)
      real , intent(out) :: parr(:)           ! Visible direct downward shortwave flux (w/m2)
                                              !    Dimensions: (ncol)
      real , intent(out) :: parf(:)           ! Visible diffuse downward shortwave flux (w/m2)
                                              !    Dimensions: (ncol)
      real , intent(out) :: uvrr(:)           ! UV direct downward shortwave flux (w/m2)
                                              !    Dimensions: (ncol)
      real , intent(out) :: uvrf(:)           ! UV diffuse downward shortwave flux (w/m2)
                                              !    Dimensions: (ncol)


      integer :: npart, pncol
      
      ! ASSERTs to catch unphysical inputs
      if (any(play   < 0.)) then
        write(error_unit,*) 'file:', __FILE__, ', line:', __LINE__
        error stop 'negative values in input: play'
      end if
      if (any(plev   < 0.)) then
        write(error_unit,*) 'file:', __FILE__, ', line:', __LINE__
        error stop 'negative values in input: plev'
      end if
      if (any(tlay   < 0.)) then
        write(error_unit,*) 'file:', __FILE__, ', line:', __LINE__
        error stop 'negative values in input: tlay'
      end if
      if (any(tlev   < 0.)) then
        write(error_unit,*) 'file:', __FILE__, ', line:', __LINE__
        error stop 'negative values in input: tlev'
      end if
      if (any(tsfc   < 0.)) then
        write(error_unit,*) 'file:', __FILE__, ', line:', __LINE__
        error stop 'negative values in input: tsfc'
      end if
      if (any(h2ovmr < 0.)) then
        write(error_unit,*) 'file:', __FILE__, ', line:', __LINE__
        error stop 'negative values in input: h2ovmr'
      end if
      if (any(o3vmr  < 0.)) then
        write(error_unit,*) 'file:', __FILE__, ', line:', __LINE__
        error stop 'negative values in input: o3vmr'
      end if
      if (any(co2vmr < 0.)) then
        write(error_unit,*) 'file:', __FILE__, ', line:', __LINE__
        error stop 'negative values in input: co2vmr'
      end if
      if (any(ch4vmr < 0.)) then
        write(error_unit,*) 'file:', __FILE__, ', line:', __LINE__
        error stop 'negative values in input: ch4vmr'
      end if
      if (any(n2ovmr < 0.)) then
        write(error_unit,*) 'file:', __FILE__, ', line:', __LINE__
        error stop 'negative values in input: n2ovmr'
      end if
      if (any(o2vmr  < 0.)) then
        write(error_unit,*) 'file:', __FILE__, ', line:', __LINE__
        error stop 'negative values in input: o2vmr'
      end if
      if (any(asdir  < 0.)) then
        write(error_unit,*) 'file:', __FILE__, ', line:', __LINE__
        error stop 'negative values in input: asdir'
      end if
      if (any(aldir  < 0.)) then
        write(error_unit,*) 'file:', __FILE__, ', line:', __LINE__
        error stop 'negative values in input: aldir'
      end if
      if (any(asdif  < 0.)) then
        write(error_unit,*) 'file:', __FILE__, ', line:', __LINE__
        error stop 'negative values in input: asdif'
      end if
      if (any(aldif  < 0.)) then
        write(error_unit,*) 'file:', __FILE__, ', line:', __LINE__
        error stop 'negative values in input: aldif'
      end if
      if (any(cld    < 0.)) then
        write(error_unit,*) 'file:', __FILE__, ', line:', __LINE__
        error stop 'negative values in input: cld'
      end if
      if (any(ciwp   < 0.)) then
        write(error_unit,*) 'file:', __FILE__, ', line:', __LINE__
        error stop 'negative values in input: ciwp'
      end if
      if (any(clwp   < 0.)) then
        write(error_unit,*) 'file:', __FILE__, ', line:', __LINE__
        error stop 'negative values in input: clwp'
      end if
      if (any(rei    < 0.)) then
        write(error_unit,*) 'file:', __FILE__, ', line:', __LINE__
        error stop 'negative values in input: rei'
      end if
      if (any(rel    < 0.)) then
        write(error_unit,*) 'file:', __FILE__, ', line:', __LINE__
        error stop 'negative values in input: rel'
      end if
      if (any(tauaer < 0.)) then
        write(error_unit,*) 'file:', __FILE__, ', line:', __LINE__
        error stop 'negative values in input: tauaer'
      end if
      if (any(ssaaer < 0.)) then
        write(error_unit,*) 'file:', __FILE__, ', line:', __LINE__
        error stop 'negative values in input: ssaaer'
      end if
!     if (any(asmaer < 0.)) then
!       write(error_unit,*) 'file:', __FILE__, ', line:', __LINE__
!       error stop 'negative values in input: asmaer'
!     end if

#ifdef _CUDA
      type(cudadeviceprop) :: prop
      real :: gmem
      integer :: err
      integer :: munits
      integer :: numDevices, numCPUsPerGPU
      real :: maxmem
#endif
      
      if (rpart > 0) then
         pncol = rpart
      else

#ifdef _CUDA
 
      err = cudaGetDeviceProperties( prop, 0)
      gmem = prop%totalGlobalMem / (1024.0 * 1024.0)
      !print *, "total GPU global memory is ", gmem , "MB"

      err = cudaGetDeviceCount(numDevices)
      !print *, "total number of GPUs is ", numDevices

      numCPUsPerGPU = ceiling( real(numCPUs) / real(numDevices) )
      !print *, "number of CPUs per GPU is ", numCPUsPerGPU

      maxmem = gmem/real(numCPUsPerGPU)
      !print *, "available GPU global memory per CPU is ", maxmem , "MB"
      
      ! dmb 2013
      ! Here 
      ! The optimal partition size is determined by the following conditions
      ! 1. Powers of 2 are the most efficient.
      ! 2. The second to largest power of 2 that can fit on 
      !    the GPU is most efficient.
      ! 3. Having a small remainder for the final partiion is inefficient.
      
      if (gmem > 5000) then
         pncol = 4096
      else if (gmem > 3000) then
         pncol = 2048
      else if (gmem > 1000) then
         pncol = 1024
      else 
         pncol = 512
      end if

      !print *,"pncol based on gmem is: ", pncol 

      pncol = pncol / numCPUsPerGPU

      !print *,"pncol based on gmem per numCPUsPerGPU is: ", pncol 

      ! the smallest allowed partition size is 32
      do err = 1, 6
          if (pncol > ncol .and. pncol>32) then 
              pncol = pncol/2
          end if
      end do
      
      ! if we have a very large number of columns, account for the 
      ! static ncol memory requirement 
      if (ncol>29000 .and. pncol>4000) then
          pncol = pncol/2
      end if

#else
      pncol = 2
      
#endif 


      !print *, "Final partition size is ", pncol
      end if
      
      
                                                      
      call rrtmg_sw_sub &
            (pncol, ncol, nlay    ,icld    , iaer, &
             play    ,plev    ,tlay    ,tlev    ,tsfc   , &
             h2ovmr , o3vmr   ,co2vmr  ,ch4vmr  ,n2ovmr ,o2vmr , &
             asdir   ,asdif   ,aldir   ,aldif   , &
             coszen  ,adjes   ,dyofyr  ,scon    , isolvar, &
             inflgsw ,iceflgsw,liqflgsw,cld , &
             tauc ,ssac ,asmc ,fsfc , &
             ciwp ,clwp ,rei ,rel , &
             tauaer  ,ssaaer  ,asmaer  ,ecaer   , &
             swuflx  ,swdflx  ,swhr    ,swuflxc ,swdflxc ,swhrc, &
             nirr    ,nirf    ,parr    ,parf    ,uvrr    ,uvrf , normFlx, zm, alat ,&
! optional I/O
             bndsolvar,indsolvar,solcycfrac)
                               


                                                      
      end subroutine rrtmg_sw                                                     


      subroutine rrtmg_sw_sub &
            (ncol ,gncol,  nlay    ,icld    , iaer, &
             gplay    ,gplev    ,gtlay    ,gtlev    ,gtsfc   , &
             gh2ovmr , go3vmr   ,gco2vmr  ,gch4vmr  ,gn2ovmr ,go2vmr , &
             gasdir   ,gasdif   ,galdir   ,galdif   , &
             gcoszen  ,adjes   ,dyofyr  ,scon    , isolvar, &
             inflgsw ,iceflgsw,liqflgsw,gcld , &
             gtauc ,gssac ,gasmc ,gfsfc , &
             gciwp ,gclwp ,grei ,grel , &
             gtauaer  ,gssaaer  ,gasmaer  ,gecaer   , &
             swuflx  ,swdflx  ,swhr    ,swuflxc ,swdflxc ,swhrc, &
             nirr    ,nirf    ,parr    ,parf    ,uvrr    ,uvrf  , normFlx, gzm, galat, &
! optional I/O
             bndsolvar,indsolvar,solcycfrac)



      use parrrsw, only : nbndsw, ngptsw, naerec, nstr, nmol, mxmol, &
                          jpband, jpb1, jpb2, rrsw_scon
      use rrsw_aer, only : rsrtaua, rsrpiza, rsrasya
      use rrsw_con, only : heatfac, oneminus, pi,  grav, avogad
      use rrsw_wvn, only : wavenum1, wavenum2
      use rrsw_cld, only : extliq1, ssaliq1, asyliq1, &
                           extice2, ssaice2, asyice2, &
                           extice3, ssaice3, asyice3, fdlice3, &
                           extice4, ssaice4, asyice4, &
                           abari, bbari, cbari, dbari, ebari, fbari
      use rrsw_wvn, only : wavenum2, ngb, icxa, nspa, nspb
      use rrsw_ref, only : preflog, tref
      use tab_xcw

      use rrsw_kg16, kao16 => kao, kbo16 => kbo, selfrefo16 => selfrefo, forrefo16 => forrefo, sfluxrefo16 => sfluxrefo
      use rrsw_kg16, ka16 => ka, kb16 => kb, selfref16 => selfref, forref16 => forref, sfluxref16 => sfluxref
      use rrsw_kg16, irradnceo16 => irradnceo, facbrghto16 => facbrghto, snsptdrko16 => snsptdrko

      use rrsw_kg17, kao17 => kao, kbo17 => kbo, selfrefo17 => selfrefo, forrefo17 => forrefo, sfluxrefo17 => sfluxrefo
      use rrsw_kg17, ka17 => ka, kb17 => kb, selfref17 => selfref, forref17 => forref, sfluxref17 => sfluxref
      use rrsw_kg17, irradnceo17 => irradnceo, facbrghto17 => facbrghto, snsptdrko17 => snsptdrko

      use rrsw_kg18, kao18 => kao, kbo18 => kbo, selfrefo18 => selfrefo, forrefo18 => forrefo, sfluxrefo18 => sfluxrefo
      use rrsw_kg18, ka18 => ka, kb18 => kb, selfref18 => selfref, forref18 => forref, sfluxref18 => sfluxref
      use rrsw_kg18, irradnceo18 => irradnceo, facbrghto18 => facbrghto, snsptdrko18 => snsptdrko

      use rrsw_kg19, kao19 => kao, kbo19 => kbo, selfrefo19 => selfrefo, forrefo19 => forrefo, sfluxrefo19 => sfluxrefo
      use rrsw_kg19, ka19 => ka, kb19 => kb, selfref19 => selfref, forref19 => forref, sfluxref19 => sfluxref
      use rrsw_kg19, irradnceo19 => irradnceo, facbrghto19 => facbrghto, snsptdrko19 => snsptdrko

      use rrsw_kg20, kao20 => kao, kbo20 => kbo, selfrefo20 => selfrefo, forrefo20 => forrefo, &
         sfluxrefo20 => sfluxrefo, absch4o20 => absch4o
      use rrsw_kg20, ka20 => ka, kb20 => kb, selfref20 => selfref, forref20 => forref, &
        sfluxref20 => sfluxref, absch420 => absch4
      use rrsw_kg20, irradnceo20 => irradnceo, facbrghto20 => facbrghto, snsptdrko20 => snsptdrko

      use rrsw_kg21, kao21 => kao, kbo21 => kbo, selfrefo21 => selfrefo, forrefo21 => forrefo, sfluxrefo21 => sfluxrefo
      use rrsw_kg21, ka21 => ka, kb21 => kb, selfref21 => selfref, forref21 => forref, sfluxref21 => sfluxref
      use rrsw_kg21, irradnceo21 => irradnceo, facbrghto21 => facbrghto, snsptdrko21 => snsptdrko

      use rrsw_kg22, kao22 => kao, kbo22 => kbo, selfrefo22 => selfrefo, forrefo22 => forrefo, sfluxrefo22 => sfluxrefo
      use rrsw_kg22, ka22 => ka, kb22 => kb, selfref22 => selfref, forref22 => forref, sfluxref22 => sfluxref
      use rrsw_kg22, irradnceo22 => irradnceo, facbrghto22 => facbrghto, snsptdrko22 => snsptdrko

      use rrsw_kg23, kao23 => kao, selfrefo23 => selfrefo, forrefo23 => forrefo, sfluxrefo23 => sfluxrefo, raylo23 => raylo
      use rrsw_kg23, ka23 => ka, selfref23 => selfref, forref23 => forref, sfluxref23 => sfluxref, rayl23 => rayl
      use rrsw_kg23, irradnceo23 => irradnceo, facbrghto23 => facbrghto, snsptdrko23 => snsptdrko

      use rrsw_kg24, kao24 => kao, kbo24 => kbo, selfrefo24 => selfrefo, forrefo24 => forrefo, sfluxrefo24 => sfluxrefo
      use rrsw_kg24, abso3ao24 => abso3ao, abso3bo24 => abso3bo, raylao24 => raylao, raylbo24 => raylbo
      use rrsw_kg24, ka24 => ka, kb24 => kb, selfref24 => selfref, forref24 => forref, sfluxref24 => sfluxref
      use rrsw_kg24, abso3a24 => abso3a, abso3b24 => abso3b, rayla24 => rayla, raylb24 => raylb
      use rrsw_kg24, irradnceo24 => irradnceo, facbrghto24 => facbrghto, snsptdrko24 => snsptdrko

      use rrsw_kg25, kao25 => kao, sfluxrefo25=>sfluxrefo
      use rrsw_kg25, abso3ao25 => abso3ao, abso3bo25 => abso3bo, raylo25 => raylo
      use rrsw_kg25, ka25 => ka, sfluxref25=>sfluxref
      use rrsw_kg25, abso3a25 => abso3a, abso3b25 => abso3b, rayl25 => rayl
      use rrsw_kg25, irradnceo25 => irradnceo, facbrghto25 => facbrghto, snsptdrko25 => snsptdrko
     
      use rrsw_kg26, sfluxrefo26 => sfluxrefo
      use rrsw_kg26, sfluxref26 => sfluxref
      use rrsw_kg26, irradnceo26 => irradnceo, facbrghto26 => facbrghto, snsptdrko26 => snsptdrko

      use rrsw_kg27, kao27 => kao, kbo27 => kbo, sfluxrefo27 => sfluxrefo, rayl27=>rayl
      use rrsw_kg27, ka27 => ka, kb27 => kb, sfluxref27 => sfluxref, raylo27=>raylo
      use rrsw_kg27, irradnceo27 => irradnceo, facbrghto27 => facbrghto, snsptdrko27 => snsptdrko

      use rrsw_kg28, kao28 => kao, kbo28 => kbo, sfluxrefo28 => sfluxrefo
      use rrsw_kg28, ka28 => ka, kb28 => kb, sfluxref28 => sfluxref
      use rrsw_kg28, irradnceo28 => irradnceo, facbrghto28 => facbrghto, snsptdrko28 => snsptdrko

      use rrsw_kg29, kao29 => kao, kbo29 => kbo, selfrefo29 => selfrefo, forrefo29 => forrefo, sfluxrefo29 => sfluxrefo
      use rrsw_kg29, absh2oo29 => absh2oo, absco2o29 => absco2o
      use rrsw_kg29, ka29 => ka, kb29 => kb, selfref29 => selfref, forref29 => forref, sfluxref29 => sfluxref
      use rrsw_kg29, absh2o29 => absh2o, absco229 => absco2
      use rrsw_kg29, irradnceo29 => irradnceo, facbrghto29 => facbrghto, snsptdrko29 => snsptdrko

! ------- Declarations



      integer , intent(in) :: ncol
      integer , intent(in) :: gncol          ! Number of horizontal columns     
      integer , intent(in) :: nlay           ! Number of model layers
      integer , intent(inout) :: icld        ! Cloud overlap method
                                             !    0: Clear only
                                             !    1: Random
                                             !    2: Maximum/random
                                             !    3: Maximum
      integer , intent(in) :: iaer
      integer , intent(in) :: dyofyr         ! Day of the year (used to get Earth/Sun
                                             !  distance if adjflx not provided)                                                      
      real , intent(in) :: adjes             ! Flux adjustment for Earth/Sun distance

      real, intent(in) :: scon               ! Solar constant (W/m2)
                                             !    Total solar irradiance averaged 
                                             !    over the solar cycle.
                                             !    If scon = 0.0, the internal solar 
                                             !    constant, which depends on the  
                                             !    value of isolvar, will be used. 
                                             !    For isolvar=-1, scon=1368.22 Wm-2,
                                             !    For isolvar=0,1,3, scon=1360.85 Wm-2,
                                             !    If scon > 0.0, the internal solar
                                             !    constant will be scaled to the 
                                             !    provided value of scon.
      integer, intent(in) :: isolvar         ! Flag for solar variability method
                                             !   -1 = (when scon .eq. 0.0): No solar variability
                                             !        and no solar cycle (Kurucz solar irradiance
                                             !        of 1368.22 Wm-2 only);
                                             !        (when scon .ne. 0.0): Kurucz solar irradiance
                                             !        scaled to scon and solar variability defined
                                             !        (optional) by setting non-zero scale factors
                                             !        for each band in bndsolvar
                                             !    0 = (when SCON .eq. 0.0): No solar variability 
                                             !        and no solar cycle (NRLSSI2 solar constant of 
                                             !        1360.85 Wm-2 for the 100-50000 cm-1 spectral 
                                             !        range only), with facular and sunspot effects 
                                             !        fixed to the mean of Solar Cycles 13-24;
                                             !        (when SCON .ne. 0.0): No solar variability 
                                             !        and no solar cycle (NRLSSI2 solar constant of 
                                             !        1360.85 Wm-2 for the 100-50000 cm-1 spectral 
                                             !        range only), is scaled to SCON
                                             !    1 = Solar variability (using NRLSSI2  solar
                                             !        model) with solar cycle contribution
                                             !        determined by fraction of solar cycle
                                             !        with facular and sunspot variations
                                             !        fixed to their mean variations over the
                                             !        average of Solar Cycles 13-24;
                                             !        two amplitude scale factors allow
                                             !        facular and sunspot adjustments from
                                             !        mean solar cycle as defined by indsolvar 
                                             !    2 = Solar variability (using NRLSSI2 solar
                                             !        model) over solar cycle determined by 
                                             !        direct specification of Mg (facular)
                                             !        and SB (sunspot) indices provided
                                             !        in indsolvar (scon = 0.0 only)
                                             !    3 = (when scon .eq. 0.0): No solar variability
                                             !        and no solar cycle (NRLSSI2 solar irradiance
                                             !        of 1360.85 Wm-2 only);
                                             !        (when scon .ne. 0.0): NRLSSI2 solar irradiance
                                             !        scaled to scon and solar variability defined
                                             !        (optional) by setting non-zero scale factors
                                             !        for each band in bndsolvar
      real, intent(in), optional :: indsolvar(2)    ! Facular and sunspot amplitude 
                                                    ! scale factors (isolvar=1), or
                                                    ! Mg and SB indices (isolvar=2)
                                                    !    Dimensions: (2)
      real, intent(in), optional :: bndsolvar(nbndsw) ! Solar variability scale factors 
                                                      ! for each shortwave band
                                                      !    Dimensions: (nbndsw=14)
      real, intent(in), optional :: solcycfrac   ! Fraction of averaged 11-year solar cycle (0-1)
                                                 !    at current time (isolvar=1)
                                                 !    0.0 represents the first day of year 1
                                                 !    1.0 represents the last day of year 11
      integer , intent(in) :: inflgsw            ! Flag for cloud optical properties
      integer , intent(in) :: iceflgsw           ! Flag for ice particle specification
      integer , intent(in) :: liqflgsw           ! Flag for liquid droplet specification
      
      real , intent(in) :: gcld(gncol, nlay)          ! Cloud fraction
                                                      !    Dimensions: (ncol,nlay)
      real , intent(in) :: gtauc(gncol,nlay,nbndsw)   ! In-cloud optical depth
                                                      !    Dimensions: (ncol,nlay,nbndsw)
      real , intent(in) :: gssac(gncol,nlay,nbndsw)   ! In-cloud single scattering albedo
                                                      !    Dimensions: (ncol,nlay,nbndsw)
      real , intent(in) :: gasmc(gncol,nlay,nbndsw)   ! In-cloud asymmetry parameter
                                                      !    Dimensions: (ncol,nlay,nbndsw)
      real , intent(in) :: gfsfc(gncol,nlay,nbndsw)   ! In-cloud forward scattering fraction
                                                      !    Dimensions: (ncol,nlay,nbndsw)
      real , intent(in) :: gciwp(gncol, nlay)         ! In-cloud ice water path (g/m2)
                                                      !    Dimensions: (ncol,nlay)
      real , intent(in) :: gclwp(gncol, nlay)         ! In-cloud liquid water path (g/m2)
                                                      !    Dimensions: (ncol,nlay)
                                                      
      real , intent(in) :: grei(gncol, nlay)          ! Cloud ice effective radius (microns)
                                                      !    Dimensions: (ncol,nlay)

      real , intent(in) :: grel(gncol, nlay)          ! Cloud water drop effective radius (microns)
                                                      !    Dimensions: (ncol,nlay)
                                                      
      
      real , intent(in) :: gplay(gncol,nlay)          ! Layer pressures (hPa, mb)
                                                      !    Dimensions: (ncol,nlay)
      real , intent(in) :: gplev(gncol,nlay+1)        ! Interface pressures (hPa, mb)
                                                      !    Dimensions: (ncol,nlay+1)
      real , intent(in) :: gtlay(gncol,nlay)          ! Layer temperatures (K)
                                                      !    Dimensions: (ncol,nlay)
      real , intent(in) :: gtlev(gncol,nlay+1)        ! Interface temperatures (K)
                                                      !    Dimensions: (ncol,nlay+1)
      real , intent(in) :: gtsfc(gncol)               ! Surface temperature (K)
                                                      !    Dimensions: (ncol)
      real , intent(in) :: gh2ovmr(gncol,nlay)        ! H2O volume mixing ratio
                                                      !    Dimensions: (ncol,nlay)
      real , intent(in) :: go3vmr(gncol,nlay)         ! O3 volume mixing ratio
                                                      !    Dimensions: (ncol,nlay)
      real , intent(in) :: gco2vmr(gncol,nlay)        ! CO2 volume mixing ratio
                                                      !    Dimensions: (ncol,nlay)
      real , intent(in) :: gch4vmr(gncol,nlay)        ! Methane volume mixing ratio
                                                      !    Dimensions: (ncol,nlay)
      real , intent(in) :: gn2ovmr(gncol,nlay)        ! Nitrous oxide volume mixing ratio
                                                      !    Dimensions: (ncol,nlay)
      real , intent(in) :: go2vmr(gncol,nlay)         ! Oxygen volume mixing ratio
                                                      !    Dimensions: (ncol,nlay)
      real , intent(in) :: gasdir(gncol)              ! UV/vis surface albedo direct rad
                                                      !    Dimensions: (ncol)
      real , intent(in) :: galdir(gncol)              ! Near-IR surface albedo direct rad
                                                      !    Dimensions: (ncol)
      real , intent(in) :: gasdif(gncol)              ! UV/vis surface albedo: diffuse rad
                                                      !    Dimensions: (ncol)
      real , intent(in) :: galdif(gncol)              ! Near-IR surface albedo: diffuse rad
                                                      !    Dimensions: (ncol)

      
      real , intent(in) :: gcoszen(gncol)             ! Cosine of solar zenith angle
                                                      !    Dimensions: (ncol)
    
      real , intent(in) :: gtauaer(gncol,nlay,nbndsw) ! Aerosol optical depth (iaer=10 only)
                                                      !    Dimensions: (ncol,nlay,nbndsw)
                                                      ! (non-delta scaled)      
      real , intent(in) :: gssaaer(gncol,nlay,nbndsw) ! Aerosol single scattering albedo (iaer=10 only)
                                                      !    Dimensions: (ncol,nlay,nbndsw)
                                                      ! (non-delta scaled)      
      real , intent(in) :: gasmaer(gncol,nlay,nbndsw) ! Aerosol asymmetry parameter (iaer=10 only)
                                                      !    Dimensions: (ncol,nlay,nbndsw)
                                                      ! (non-delta scaled)      
      real , intent(in) :: gecaer(gncol,nlay,naerec)  ! Aerosol optical depth at 0.55 micron (iaer=6 only)
                                                      !    Dimensions: (ncol,nlay,naerec)
                                                      ! (non-delta scaled)      
      integer , intent(in) :: normFlx                 ! Normalize fluxes flag
                                                      !  0 = no normalization
                                                      !  1 = normalize fluxes ( / (scon * coszen) )
      real, intent(in) :: gzm(gncol, nlay)            ! Heights of level midpoints
                                                      !    Dimensions: (ncol,nlay)
      real, intent(in) :: galat(gncol)                ! Latitudes of columns
                                                      !    Dimensions: (ncol)
                                              
! ----- Output -----

      real , intent(out) :: swuflx(:,:)               ! Total sky shortwave upward flux (W/m2)
                                                      !    Dimensions: (ncol,nlay+1)
      real , intent(out) :: swdflx(:,:)               ! Total sky shortwave downward flux (W/m2)
                                                      !    Dimensions: (ncol,nlay+1)
      real , intent(out) :: swhr(:,:)                 ! Total sky shortwave radiative heating rate (K/d)
                                                      !    Dimensions: (ncol,nlay)
      real , intent(out) :: swuflxc(:,:)              ! Clear sky shortwave upward flux (W/m2)
                                                      !    Dimensions: (ncol,nlay+1)
      real , intent(out) :: swdflxc(:,:)              ! Clear sky shortwave downward flux (W/m2)
                                                      !    Dimensions: (ncol,nlay+1)
      real , intent(out) :: swhrc(:,:)                ! Clear sky shortwave radiative heating rate (K/d)
                                                      !    Dimensions: (ncol,nlay)


   ! Output added for Land/Surface process
      real , intent(out) :: nirr(:)                   ! Near-IR direct downward shortwave flux (w/m2)
                                                      !    Dimensions: (ncol)
      real , intent(out) :: nirf(:)                   ! Near-IR diffuse downward shortwave flux (w/m2)
                                                      !    Dimensions: (ncol)
      real , intent(out) :: parr(:)                   ! Visible direct downward shortwave flux (w/m2)
                                                      !    Dimensions: (ncol)
      real , intent(out) :: parf(:)                   ! Visible diffuse downward shortwave flux (w/m2)
                                                      !    Dimensions: (ncol)
      real , intent(out) :: uvrr(:)                   ! UV direct downward shortwave flux (w/m2)
                                                      !    Dimensions: (ncol)
      real , intent(out) :: uvrf(:)                   ! UV diffuse downward shortwave flux (w/m2)
                                                      !    Dimensions: (ncol)

! ----- Local -----

! Control
     
      integer :: istart              ! beginning band of calculation
      integer :: iend                ! ending band of calculation
      integer :: icpr                ! cldprop/cldprmc use flag
      integer :: iout                ! output option flag
  
      integer :: idelm               ! delta-m scaling flag
                                     ! [0 = direct and diffuse fluxes are unscaled]
                                     ! [1 = direct and diffuse fluxes are scaled]
                                     ! (total downward fluxes are always delta scaled)
      integer :: isccos              ! instrumental cosine response flag (inactive)
      integer :: iplon               ! column loop index
      integer :: i                   ! layer loop index                       ! jk
      integer :: ib                  ! band loop index                        ! jsw
      integer :: ia, ig              ! indices
      integer :: k                   ! layer loop index
      integer :: ims                 ! value for changing mcica permute seed
      integer :: imca                ! flag for mcica [0=off, 1=on]

      real :: zepsec, zepzen         ! epsilon
      real :: zdpgcp                 ! flux to heating conversion ratio

! Atmosphere


      real  :: coldry(ncol,nlay+1)          ! dry air column amount
      real  :: wkl(ncol,mxmol,nlay)         ! molecular amounts (mol/cm-2)

      real  :: cossza(ncol)                 ! Cosine of solar zenith angle
      real  :: adjflux(jpband)              ! adjustment for current Earth/Sun distance

                                            !  default value of 1368.22 Wm-2 at 1 AU
      real  :: albdir(ncol,nbndsw)          ! surface albedo, direct          ! zalbp
      real  :: albdif(ncol,nbndsw)          ! surface albedo, diffuse         ! zalbd
      
      real  :: rdl(ncol), adl(ncol)

      real  :: CDF(ncol,nlay,ngptsw)
      real  :: CDF2(ncol,nlay,ngptsw)
      real  :: CDF3(ncol,nlay,ngptsw)
      real  :: alpha(ncol,nlay)


! Atmosphere - setcoef
      integer  :: laytrop(ncol)             ! tropopause layer index
      integer  :: layswtch(ncol)            ! tropopause layer index
      integer  :: laylow(ncol)              ! tropopause layer index
      integer  :: jp(ncol,nlay+1)           ! 
      integer  :: jt(ncol,nlay+1)           !
      integer  :: jt1(ncol,nlay+1)          !

      real  :: colh2o(ncol,nlay+1)          ! column amount (h2o)
      real  :: colco2(ncol,nlay+1)          ! column amount (co2)
      real  :: colo3(ncol,nlay+1)           ! column amount (o3)
      real  :: coln2o(ncol,nlay+1)          ! column amount (n2o)
      real  :: colch4(ncol,nlay+1)          ! column amount (ch4)
      real  :: colo2(ncol,nlay+1)           ! column amount (o2)
      real  :: colmol(ncol,nlay+1)          ! column amount
      real  :: co2mult(ncol,nlay+1)         ! column amount 

      integer  :: indself(ncol,nlay+1) 
      integer  :: indfor(ncol,nlay+1) 
      real  :: selffac(ncol,nlay+1) 
      real  :: selffrac(ncol,nlay+1) 
      real  :: forfac(ncol,nlay+1) 
      real  :: forfrac(ncol,nlay+1) 

      real  :: fac00(ncol,nlay+1) , fac01(ncol,nlay+1) , &
               fac10(ncol,nlay+1) , fac11(ncol,nlay+1)  
      
      real :: play(ncol,nlay)               ! Layer pressures (hPa, mb)
                                            !    Dimensions: (ncol,nlay)
      real :: plev(ncol,nlay+1)             ! Interface pressures (hPa, mb)
                                            !    Dimensions: (ncol,nlay+1)
      real :: tlay(ncol,nlay)               ! Layer temperatures (K)
                                            !    Dimensions: (ncol,nlay)
      real :: tlev(ncol,nlay+1)             ! Interface temperatures (K)
                                            !    Dimensions: (ncol,nlay+1)
      real :: tsfc(ncol)                    ! Surface temperature (K)
                                            !    Dimensions: (ncol)
                                                      
      real :: coszen(ncol)   
      real :: swdflx_at_top(gncol)           ! swdflx at TOA (ncol)

! Atmosphere/clouds - cldprop
      integer :: ncbands             ! number of cloud spectral bands
 

      real   :: cld(ncol,nlay)              ! Cloud fraction
      real   :: tauc(ncol,nlay,nbndsw)      ! In-cloud optical depth
      real   :: ssac(ncol,nlay,nbndsw)      ! In-cloud single scattering 
      real   :: asmc(ncol,nlay,nbndsw)      ! In-cloud asymmetry parameter
      real   :: fsfc(ncol,nlay,nbndsw)      ! In-cloud forward scattering fraction
      real   :: ciwp(ncol,nlay)             ! In-cloud ice water path (g/m2)
      real   :: clwp(ncol,nlay)             ! In-cloud liquid water path (g/m2)
      real   :: rei(ncol,nlay)              ! Cloud ice effective radius (microns)
      real   :: rel(ncol,nlay)              ! Cloud water drop effective radius (microns)
      
      real   :: alat(ncol)
      real   :: zm(ncol, nlay)
                                                      
      real, dimension(ncol) :: znirr,znirf,zparr,zparf,zuvrr,zuvrf
      
      real  :: taucmc(ncol,nlay+1,ngptsw)    ! in-cloud optical depth [mcica]
      real  :: taormc(ncol,nlay+1,ngptsw)    ! unscaled in-cloud optical depth [mcica]
      real  :: ssacmc(ncol,nlay+1,ngptsw)    ! in-cloud single scattering albedo [mcica]
      real  :: asmcmc(ncol,nlay+1,ngptsw)    ! in-cloud asymmetry parameter [mcica]
      real  :: fsfcmc(ncol,nlay+1,ngptsw)    ! in-cloud forward scattering fraction [mcica]
      
      
      real :: cldfmcl(ncol,nlay+1,ngptsw)    ! cloud fraction [mcica]
      real :: ciwpmcl(ncol,nlay+1,ngptsw)    ! in-cloud ice water path [mcica]
      real :: clwpmcl(ncol,nlay+1,ngptsw)    ! in-cloud liquid water path [mcica]
                                                     


! Atmosphere/clouds/aerosol - spcvrt,spcvmc
      real  :: ztauc(ncol,nlay+1,nbndsw)     ! cloud optical depth
      real  :: ztaucorig(ncol,nlay+1,nbndsw) ! unscaled cloud optical depth
      real  :: zasyc(ncol,nlay+1,nbndsw)     ! cloud asymmetry parameter 
                                             !  (first moment of phase function)
      real  :: zomgc(ncol,nlay+1,nbndsw)     ! cloud single scattering albedo
   
      real  :: taua(ncol, nlay+1, nbndsw)
      real  :: asya(ncol, nlay+1, nbndsw)
      real  :: omga(ncol, nlay+1, nbndsw)
   

      real  :: zbbfu(ncol,nlay+2)            ! temporary upward shortwave flux (w/m2)
      real  :: zbbfd(ncol,nlay+2)            ! temporary downward shortwave flux (w/m2)
      real  :: zbbcu(ncol,nlay+2)            ! temporary clear sky upward shortwave flux (w/m2)
      real  :: zbbcd(ncol,nlay+2)            ! temporary clear sky downward shortwave flux (w/m2)
      real  :: zbbfddir(ncol,nlay+2)         ! temporary downward direct shortwave flux (w/m2)
      real  :: zbbcddir(ncol,nlay+2)         ! temporary clear sky downward direct shortwave flux (w/m2)
      real  :: zuvfd(ncol,nlay+2)            ! temporary UV downward shortwave flux (w/m2)
      real  :: zuvcd(ncol,nlay+2)            ! temporary clear sky UV downward shortwave flux (w/m2)
      real  :: zuvfddir(ncol,nlay+2)         ! temporary UV downward direct shortwave flux (w/m2)
      real  :: zuvcddir(ncol,nlay+2)         ! temporary clear sky UV downward direct shortwave flux (w/m2)
      real  :: znifd(ncol,nlay+2)            ! temporary near-IR downward shortwave flux (w/m2)
      real  :: znicd(ncol,nlay+2)            ! temporary clear sky near-IR downward shortwave flux (w/m2)
      real  :: znifddir(ncol,nlay+2)         ! temporary near-IR downward direct shortwave flux (w/m2)
      real  :: znicddir(ncol,nlay+2)         ! temporary clear sky near-IR downward direct shortwave flux (w/m2)

! Optional output fields 
      real  :: swnflx(ncol,nlay+2)           ! Total sky shortwave net flux (W/m2)
      real  :: swnflxc(ncol,nlay+2)          ! Clear sky shortwave net flux (W/m2)
      real  :: dirdflux(ncol,nlay+2)         ! Direct downward shortwave surface flux
      real  :: difdflux(ncol,nlay+2)         ! Diffuse downward shortwave surface flux
      real  :: uvdflx(ncol,nlay+2)           ! Total sky downward shortwave flux, UV/vis  
      real  :: nidflx(ncol,nlay+2)           ! Total sky downward shortwave flux, near-IR 
      real  :: dirdnuv(ncol,nlay+2)          ! Direct downward shortwave flux, UV/vis
      real  :: difdnuv(ncol,nlay+2)          ! Diffuse downward shortwave flux, UV/vis
      real  :: dirdnir(ncol,nlay+2)          ! Direct downward shortwave flux, near-IR
      real  :: difdnir(ncol,nlay+2)          ! Diffuse downward shortwave flux, near-IR

! Solar variability
      real :: svar_f                 ! Solar variability facular multiplier
      real :: svar_s                 ! Solar variability sunspot multiplier
      real :: svar_i                 ! Solar variability baseline irradiance multiplier
      real :: svar_f_bnd(jpband)     ! Solar variability facular multiplier (by band)
      real :: svar_s_bnd(jpband)     ! Solar variability sunspot multiplier (by band)
      real :: svar_i_bnd(jpband)     ! Solar variability baseline irradiance multiplier (by band)


      
      real gpu_device :: zgco(ncol,ngptsw,nlay+1), zomco(ncol,ngptsw,nlay+1)  
      real gpu_device :: zrdnd(ncol,ngptsw,nlay+1) 
      real gpu_device :: zref(ncol,ngptsw,nlay+1)  , zrefo(ncol,ngptsw,nlay+1)  
      real gpu_device :: zrefd(ncol,ngptsw,nlay+1)  , zrefdo(ncol,ngptsw,nlay+1)  
      real gpu_device :: ztauo(ncol,ngptsw,nlay)  
      real gpu_device :: zdbt(ncol,ngptsw,nlay+1)  ,ztdbt(ncol,ngptsw,nlay+1)   
      real gpu_device :: ztra(ncol,ngptsw,nlay+1)  , ztrao(ncol,ngptsw,nlay+1)  
      real gpu_device :: ztrad(ncol,ngptsw,nlay+1)  , ztrado(ncol,ngptsw,nlay+1)  
      real gpu_device :: zfd(ncol,ngptsw,nlay+1)  , zfu(ncol,ngptsw,nlay+1)  
      real gpu_device :: zsflxzen(ncol,ngptsw)
      real gpu_device :: ssi(ncol,ngptsw)
      real gpu_device :: ztaur(ncol,nlay,ngptsw), ztaug(ncol,nlay,ngptsw) 

      integer :: npartc, npart, npartb, cldflag(gncol), profic(gncol), profi(gncol)

      real , parameter :: amd = 28.9660     ! Effective molecular weight of dry air (g/mol)
      real , parameter :: amw = 18.0160     ! Molecular weight of water vapor (g/mol)


! Set molecular weight ratios (for converting mmr to vmr)
!  e.g. h2ovmr = h2ommr * amdw)
      real , parameter :: amdw = 1.607793   ! Molecular weight of dry air / water vapor
      real , parameter :: amdc = 0.658114   ! Molecular weight of dry air / carbon dioxide
      real , parameter :: amdo = 0.603428   ! Molecular weight of dry air / ozone
      real , parameter :: amdm = 1.805423   ! Molecular weight of dry air / methane
      real , parameter :: amdn = 0.658090   ! Molecular weight of dry air / nitrous oxide
      real , parameter :: amdo2 = 0.905140  ! Molecular weight of dry air / oxygen

      real , parameter :: sbc = 5.67e-08    ! Stefan-Boltzmann constant (W/m2K4)

      integer  :: isp, l, ix, n, imol       ! Loop indices
      real  :: amm, summol                  ! 
      real  :: adjflx                       ! flux adjustment for Earth/Sun distance
      integer :: prt
      integer :: piplon
      
      integer :: ipart, cols, cole, colr, ncolc, ncolb
      integer :: irng, cc, ncolst
      real :: tt1, tt2

      real :: solvar(jpband)                   ! solar constant scaling factor by band
                                               !  Dimension(jpband=29)
      real :: indsolvar_scl(2)                 ! Adjusted facular and sunspot amplitude 
                                               ! scale factors (isolvar=1)
      real :: indsolvar_ndx(2)                 ! Facular and sunspot indices (isolvar=2)
      real :: solcycfr                         ! Local solar cycle fraction (default = 0.0
                                               ! unless solcycfrac is present)

      real, parameter ::  solcycfrac_min = 0.0189    ! Solar cycle fraction at solar minimum
      real, parameter ::  solcycfrac_max = 0.3750    ! Solar cycle fraction at solar maximum
      real, parameter ::  fracdiff_min2max = 0.3561  ! 0.3750 - 0.0189
      real, parameter ::  fracdiff_max2min = 0.6439  ! 1.0189 - 0.3750
      real :: wgt                              ! Weighting factor for amplitude scale factor adjustment
      real :: svar_f_0, svar_s_0               ! Solar variability indices for current fractional
                                               !  position in typical solar cycle, interpolated
                                               !  from lookup table of values over solar cycle
      real :: svar_cprim                       ! Solar variability intermediate value
      real :: svar_r                           ! Solar variability intermediate value
      integer :: sfid                          ! Solar variability solar cycle fraction index
      real :: tmp_f_0, tmp_s_0                 ! Solar variability temporary quantities
      real :: fraclo, frachi, intfrac          ! Solar variability interpolation factors

! Mean quiet sun, facular brightening, and sunspot dimming coefficient terms (NRLSSI2, 100-50000 cm-1), 
! spectrally integrated (from hi-res values after mapping to g-point space)
      real, parameter :: Iint = 1360.37     ! Solar quiet sun irradiance term, integrated
      real, parameter :: Fint = 0.996047    ! Solar facular brightening term (index-offset), integrated
      real, parameter :: Sint = -0.511590   ! Solar sunspot dimming term (index-offset), integrated
      real, parameter :: Foffset = 0.14959542    ! Solar variability facular offset
      real, parameter :: Soffset = 0.00066696    ! Solar variability sunspot offset

! Mg and SB indices for average solar cycle integrated over solar cycle
      real, parameter :: svar_f_avg = 0.1567652  ! Solar variability NRLSSI2 Mg "Bremen" index 
                                                 !  time-averaged over Solar Cycles 13-24
                                                 !  and averaged over solar cycle (132 values
                                                 !  excluding end points in Mg and SB arrays)
      real, parameter :: svar_s_avg = 909.71260  ! Solar variability NRLSSI2 SB "SPOT67" index 
                                                 !  time-averaged over Solar Cycles 13-24
                                                 !  and averaged over solar cycle (132 values
                                                 !  excluding end points in Mg and SB arrays)
      integer, parameter :: nsolfrac = 134       ! Number of elements in solar arrays 
                                                 !  132 values (excluding end points) represent
                                                 !  the center dates of the 12 months per year 
                                                 !  over the mean 11-year solar cycle;
                                                 !  2 end points represent the first day of the 
                                                 !  first month of year 1 and the last day of
                                                 !  the last month of year 11
      real :: intrvl_len                         !  Fractional interval length of mgavgcyc
                                                 !  and sbavgcyc
      real :: intrvl_len_hf                      !  Fractional half interval length of mgavgcyc
                                                 !  and sbavgcyc

! Mg and SB index look-up tables for average solar cycle as a function of solar cycle
      real :: mgavgcyc(nsolfrac)               ! Facular index from NRLSSI2 Mg "Bremen" index 
                                               !  time-averaged over Solar Cycles 13-24
      real :: sbavgcyc(nsolfrac)               ! Sunspot index from NRLSSI2 SB "SPOT67" index 
                                               !  time-averaged over Solar Cycles 13-24
      mgavgcyc(:) = (/ &
        &   0.150737,  0.150746,  0.150733,  0.150718,  0.150725,  0.150762, &
        &   0.150828,  0.150918,  0.151017,  0.151113,  0.151201,  0.151292, &
        &   0.151403,  0.151557,  0.151766,  0.152023,  0.152322,  0.152646, &
        &   0.152969,  0.153277,  0.153579,  0.153899,  0.154252,  0.154651, &
        &   0.155104,  0.155608,  0.156144,  0.156681,  0.157178,  0.157605, &
        &   0.157971,  0.158320,  0.158702,  0.159133,  0.159583,  0.160018, &
        &   0.160408,  0.160725,  0.160960,  0.161131,  0.161280,  0.161454, &
        &   0.161701,  0.162034,  0.162411,  0.162801,  0.163186,  0.163545, &
        &   0.163844,  0.164029,  0.164054,  0.163910,  0.163621,  0.163239, &
        &   0.162842,  0.162525,  0.162344,  0.162275,  0.162288,  0.162369, &
        &   0.162500,  0.162671,  0.162878,  0.163091,  0.163251,  0.163320, &
        &   0.163287,  0.163153,  0.162927,  0.162630,  0.162328,  0.162083, &
        &   0.161906,  0.161766,  0.161622,  0.161458,  0.161266,  0.161014, &
        &   0.160666,  0.160213,  0.159690,  0.159190,  0.158831,  0.158664, &
        &   0.158634,  0.158605,  0.158460,  0.158152,  0.157691,  0.157152, &
        &   0.156631,  0.156180,  0.155827,  0.155575,  0.155406,  0.155280, &
        &   0.155145,  0.154972,  0.154762,  0.154554,  0.154388,  0.154267, &
        &   0.154152,  0.154002,  0.153800,  0.153567,  0.153348,  0.153175, &
        &   0.153044,  0.152923,  0.152793,  0.152652,  0.152510,  0.152384, &
        &   0.152282,  0.152194,  0.152099,  0.151980,  0.151844,  0.151706, &
        &   0.151585,  0.151496,  0.151437,  0.151390,  0.151347,  0.151295, &
        &   0.151220,  0.151115,  0.150993,  0.150883,  0.150802,  0.150752, &
        &   0.150729,  0.150737/)
      sbavgcyc(:) = (/ &
        &    50.3550,   44.1322,   52.0179,   59.2231,   66.3702,   71.7545, &
        &    76.8671,   83.4723,   91.1574,   98.4915,  105.3173,  115.1791, &
        &   130.9432,  155.0483,  186.5379,  221.5456,  256.9212,  291.5276, &
        &   325.2953,  356.4789,  387.2470,  422.8557,  466.1698,  521.5139, &
        &   593.2833,  676.6234,  763.6930,  849.1200,  928.4259,  994.9705, &
        &  1044.2605, 1087.5703, 1145.0623, 1224.3491, 1320.6497, 1413.0979, &
        &  1472.1591, 1485.7531, 1464.1610, 1439.1617, 1446.2449, 1496.4323, &
        &  1577.8394, 1669.5933, 1753.0408, 1821.9296, 1873.2789, 1906.5240, &
        &  1920.4482, 1904.6881, 1861.8397, 1802.7661, 1734.0215, 1665.0562, &
        &  1608.8999, 1584.8208, 1594.0162, 1616.1486, 1646.6031, 1687.1962, &
        &  1736.4778, 1787.2419, 1824.9084, 1835.5236, 1810.2161, 1768.6124, &
        &  1745.1085, 1748.7762, 1756.1239, 1738.9929, 1700.0656, 1658.2209, &
        &  1629.2925, 1620.9709, 1622.5157, 1623.4703, 1612.3083, 1577.3031, &
        &  1516.7953, 1430.0403, 1331.5112, 1255.5171, 1226.7653, 1241.4419, &
        &  1264.6549, 1255.5559, 1203.0286, 1120.2747, 1025.5101,  935.4602, &
        &   855.0434,  781.0189,  718.0328,  678.5850,  670.4219,  684.1906, &
        &   697.0376,  694.8083,  674.1456,  638.8199,  602.3454,  577.6292, &
        &   565.6213,  553.7846,  531.7452,  503.9732,  476.9708,  452.4296, &
        &   426.2826,  394.6636,  360.1086,  324.9731,  297.2957,  286.1536, &
        &   287.4195,  288.9029,  282.7594,  267.7211,  246.6594,  224.7318, &
        &   209.2318,  204.5217,  204.1653,  200.0440,  191.0689,  175.7699, &
        &   153.9869,  128.4389,  103.8445,   85.6083,   73.6264,   64.4393, &
        &    56.5779,   50.3550/)


! Initializations

      zepsec = 1.e-06
      zepzen = 1.e-10
      oneminus = 1.0 - zepsec
      pi = 2. * asin(1.)
      irng = 0

      istart = jpb1
      iend = jpb2
      iout = 0
      icpr = 1
      ims = 2


      solvar(:) = 1.0
      adjflux(:) = 1.0
      svar_f = 1.0 
      svar_s = 1.0 
      svar_i = 1.0 
      svar_f_bnd(:) = 1.0 
      svar_s_bnd(:) = 1.0 
      svar_i_bnd(:) = 1.0 

! Adjust amplitude scaling of mean solar cycle to be 1.0 at solar minimum (solcycfrac_min=0.0189),
! to be the requested indsolvar at solar maximum (solcycfrac_max=0.3750), and to vary between 
! those values at intervening values of solcycfrac. 
      if (isolvar .eq. 1) then 
! Check for presence of indsolvar and solcycfrac when isolvar = 1. 
! Use a solar cycle fraction of 0.0 and no scaling by default unless both indsolvar and solcycfrac are present. 
         solcycfr = 0.0
         indsolvar_scl(1:2) = 1.0
         if (present(indsolvar) .and. present(solcycfrac)) then 
            solcycfr = solcycfrac
            if (indsolvar(1).ne.1.0.or.indsolvar(2).ne.1.0) then 
               if (solcycfrac .ge. 0.0 .and. solcycfrac .lt. solcycfrac_min) then
                  wgt = (solcycfrac+1.0-solcycfrac_max)/fracdiff_max2min
                  indsolvar_scl(1) = indsolvar(1) + wgt * (1.0-indsolvar(1))
                  indsolvar_scl(2) = indsolvar(2) + wgt * (1.0-indsolvar(2))
               endif
               if (solcycfrac .ge. solcycfrac_min .and. solcycfrac .le. solcycfrac_max) then
                  wgt = (solcycfrac-solcycfrac_min)/fracdiff_min2max
                  indsolvar_scl(1) = 1.0 + wgt * (indsolvar(1)-1.0)
                  indsolvar_scl(2) = 1.0 + wgt * (indsolvar(2)-1.0)
               endif
               if (solcycfrac .gt. solcycfrac_max .and. solcycfrac .le. 1.0) then
                  wgt = (solcycfrac-solcycfrac_max)/fracdiff_max2min
                  indsolvar_scl(1) = indsolvar(1) + wgt * (1.0-indsolvar(1))
                  indsolvar_scl(2) = indsolvar(2) + wgt * (1.0-indsolvar(2))
               endif
            endif
         endif
      endif

! Check for presence of indsolvar when isolvar = 2. 
      if (isolvar .eq. 2) then 
! Use mean solar cycle facular and sunspot indices by default unless indsolvar is present
         indsolvar_ndx(1) = svar_f_avg
         indsolvar_ndx(2) = svar_s_avg
         if (present(indsolvar)) then 
            indsolvar_ndx(1) = indsolvar(1)
            indsolvar_ndx(2) = indsolvar(2)
         endif
      endif

! Set flux adjustment for current Earth/Sun distance (two options).
! 1) Use Earth/Sun distance flux adjustment provided by GCM (input as adjes);
      adjflx = adjes
!
! 2) Calculate Earth/Sun distance from DYOFYR, the cumulative day of the year.
!    (Set adjflx to 1. to use constant Earth/Sun distance of 1 AU). 

! MATMAT We supply dyofyr for MCICA exponential cloud
!        overlap purposes. We are passing in the MAPL
!        DIST as ADJES in the Solar Grid Comp. 

!     if (dyofyr .gt. 0) then
!        adjflx = earth_sun(dyofyr)
!     endif

! Apply selected solar variability option based on ISOLVAR and input 
! solar constant.
! For scon = 0, use internally defined solar constant, which is
! 1368.22 Wm-2 (for ISOLVAR=-1) and 1360.85 Wm-2 (for ISOLVAR=0,3;
! options ISOLVAR=1,2 model solar cycle variations from 1360.85 Wm-2)
!
! SCON = 0 
! Use internal TSI value
      SCON_IS_0: if (scon .eq. 0.0) then 

!   No solar cycle and no solar variability (Kurucz solar source function)
!   Apply constant scaling by band if first element of bndsolvar specified
         if (isolvar .eq. -1) then
            solvar(jpb1:jpb2) = 1.0
            if (present(bndsolvar)) solvar(jpb1:jpb2) = bndsolvar(:)
         endif 

!   Mean solar cycle with no solar variability (NRLSSI2 model solar irradiance)
!   Quiet sun, facular, and sunspot terms averaged over the mean solar cycle 
!   (defined as average of Solar Cycles 13-24).
         if (isolvar .eq. 0) then
            svar_f = 1.0
            svar_s = 1.0
            svar_i = 1.0
         endif 

!   Mean solar cycle with solar variability (NRLSSI2 model)
!   Facular and sunspot terms interpolated from LUTs to input solar cycle 
!   fraction for mean solar cycle. Scalings defined below to convert from 
!   averaged Mg and SB terms to Mg and SB terms interpolated here.
!   (Includes optional facular and sunspot amplitude scale factors)
         if (isolvar .eq. 1) then
!   Interpolate svar_f_0 and svar_s_0 from lookup tables using provided solar cycle fraction
            if (solcycfr .le. 0.0) then
               tmp_f_0 = mgavgcyc(1)
               tmp_s_0 = sbavgcyc(1)
            elseif (solcycfr .ge. 1.0) then
               tmp_f_0 = mgavgcyc(nsolfrac)
               tmp_s_0 = sbavgcyc(nsolfrac)
            else
               intrvl_len = 1.0 / (nsolfrac-2)
               intrvl_len_hf = 0.5 * intrvl_len
!   Initial half interval (1)
               if (solcycfr .le. intrvl_len_hf) then 
                  sfid = 1
                  fraclo = 0.0
                  frachi = intrvl_len_hf
               endif
!   Main whole intervals (131)
               if (solcycfr .gt. intrvl_len_hf .and. solcycfr .lt. 1.0-intrvl_len_hf) then 
                  sfid = floor((solcycfr-intrvl_len_hf) * (nsolfrac-2)) + 2
                  fraclo = (sfid-2) * intrvl_len + intrvl_len_hf
                  frachi = fraclo + intrvl_len
               endif
!   Final half interval (1)
               if (solcycfr .ge. 1.0-intrvl_len_hf) then 
                  sfid = (nsolfrac-2) + 1
                  fraclo = 1.0 - intrvl_len_hf
                  frachi = 1.0
               endif
               intfrac = (solcycfr - fraclo) / (frachi - fraclo)
               tmp_f_0 = mgavgcyc(sfid) + intfrac * (mgavgcyc(sfid+1) - mgavgcyc(sfid))
               tmp_s_0 = sbavgcyc(sfid) + intfrac * (sbavgcyc(sfid+1) - sbavgcyc(sfid))
            endif
            svar_f_0 = tmp_f_0
            svar_s_0 = tmp_s_0
            svar_f = indsolvar_scl(1) * (svar_f_0 - Foffset) / (svar_f_avg - Foffset)
            svar_s = indsolvar_scl(2) * (svar_s_0 - Soffset) / (svar_s_avg - Soffset)
            svar_i = 1.0
         endif 

!   Specific solar cycle with solar variability (NRLSSI2 model)
!   Facular and sunspot index terms input directly to model specific 
!   solar cycle.  Scalings defined below to convert from averaged
!   Mg and SB terms to specified Mg and SB terms. 
         if (isolvar .eq. 2) then
            svar_f = (indsolvar_ndx(1) - Foffset) / (svar_f_avg - Foffset)
            svar_s = (indsolvar_ndx(2) - Soffset) / (svar_s_avg - Soffset)
            svar_i = 1.0
         endif 

!   Mean solar cycle with no solar variability (NRLSSI2 model)
!   Averaged facular, sunspot and quiet sun terms from mean solar cycle 
!   (derived as average of Solar Cycles 13-24). This information is built
!   into coefficient terms specified by g-point elsewhere. Separate
!   scaling by spectral band is applied as defined by bndsolvar. 
         if (isolvar .eq. 3) then
            solvar(jpb1:jpb2) = 1.0
            if (present(bndsolvar)) solvar(jpb1:jpb2) = bndsolvar(:)
            do ib = jpb1,jpb2
               svar_f_bnd(ib) = solvar(ib)
               svar_s_bnd(ib) = solvar(ib)
               svar_i_bnd(ib) = solvar(ib)
            enddo
         endif 

      endif SCON_IS_0

! SCON > 0 
! Scale from internal TSI to externally specified TSI value (scon)
      SCON_GT_0: if (scon .gt. 0.0) then 

!   No solar cycle and no solar variability (Kurucz solar source function)
!   Scale from internal solar constant to requested solar constant.
!   Apply optional constant scaling by band if first element of bndsolvar > 0.0
         if (isolvar .eq. -1) then
            if (.not. present(bndsolvar)) solvar(jpb1:jpb2) = scon / rrsw_scon 
            if (present(bndsolvar)) solvar(jpb1:jpb2) = bndsolvar(:) * scon / rrsw_scon 
         endif 

!   Mean solar cycle with no solar variability (NRLSSI2 model solar irradiance)
!   Quiet sun, facular, and sunspot terms averaged over the mean solar cycle 
!   (defined as average of Solar Cycles 13-24).
!   Scale internal solar constant to requested solar constant. 
!!   Fint is provided as the product of (svar_f_avg-Foffset) and Fint, 
!!   Sint is provided as the product of (svar_s_avg-Soffset) and Sint
         if (isolvar .eq. 0) then
            svar_cprim = Fint + Sint + Iint
            svar_r = scon / svar_cprim
            svar_f = svar_r
            svar_s = svar_r
            svar_i = svar_r
         endif 

!   Mean solar cycle with solar variability (NRLSSI2 model)
!   Facular and sunspot terms interpolated from LUTs to input solar cycle 
!   fraction for mean solar cycle. Scalings defined below to convert from 
!   averaged Mg and SB terms to Mg and SB terms interpolated here.
!   Scale internal solar constant to requested solar constant. 
!   (Includes optional facular and sunspot amplitude scale factors)
         if (isolvar .eq. 1) then
!   Interpolate svar_f_0 and svar_s_0 from lookup tables using provided solar cycle fraction
            if (solcycfr .le. 0.0) then
               tmp_f_0 = mgavgcyc(1)
               tmp_s_0 = sbavgcyc(1)
            elseif (solcycfr .ge. 1.0) then
               tmp_f_0 = mgavgcyc(nsolfrac)
               tmp_s_0 = sbavgcyc(nsolfrac)
            else
               intrvl_len = 1.0 / (nsolfrac-2)
               intrvl_len_hf = 0.5 * intrvl_len
!   Initial half interval (1)
               if (solcycfr .le. intrvl_len_hf) then 
                  sfid = 1
                  fraclo = 0.0
                  frachi = intrvl_len_hf
               endif
!   Main whole intervals (131)
               if (solcycfr .gt. intrvl_len_hf .and. solcycfr .lt. 1.0-intrvl_len_hf) then 
                  sfid = floor((solcycfr-intrvl_len_hf) * (nsolfrac-2)) + 2
                  fraclo = (sfid-2) * intrvl_len + intrvl_len_hf
                  frachi = fraclo + intrvl_len
               endif
!   Final half interval (1)
               if (solcycfr .ge. 1.0-intrvl_len_hf) then 
                  sfid = (nsolfrac-2) + 1
                  fraclo = 1.0 - intrvl_len_hf
                  frachi = 1.0
               endif
               intfrac = (solcycfr - fraclo) / (frachi - fraclo)
               tmp_f_0 = mgavgcyc(sfid) + intfrac * (mgavgcyc(sfid+1) - mgavgcyc(sfid))
               tmp_s_0 = sbavgcyc(sfid) + intfrac * (sbavgcyc(sfid+1) - sbavgcyc(sfid))
            endif
            svar_f_0 = tmp_f_0
            svar_s_0 = tmp_s_0
!   Define Cprime 
!            svar_cprim = indsolvar(1) * svar_f_avg * Fint + indsolvar(2) * svar_s_avg * Sint + Iint
!   Fint is provided as the product of (svar_f_avg-Foffset) and Fint, 
!   Sint is provided as the product of (svar_s_avg-Soffset) and Sint
            svar_i = (scon - (indsolvar_scl(1) * Fint + indsolvar_scl(2) * Sint)) / Iint
            svar_f = indsolvar_scl(1) * (svar_f_0 - Foffset) / (svar_f_avg - Foffset)
            svar_s = indsolvar_scl(2) * (svar_s_0 - Soffset) / (svar_s_avg - Soffset)
         endif 

!   Specific solar cycle with solar variability (NRLSSI2 model)
!   (Not available for SCON > 0)
!         if (isolvar .eq. 2) then
!            scon = 0.0
!            svar_f = (indsolvar_ndx(1) - Foffset) / (svar_f_avg - Foffset)
!            svar_s = (indsolvar_ndx(2) - Soffset) / (svar_s_avg - Soffset)
!            svar_i = 1.0
!         endif 

! MAT The code below is provided by Peter Norris
         if (isolvar .eq. 2) then
            svar_f = (indsolvar_ndx(1) - Foffset) / (svar_f_avg - Foffset)
            svar_s = (indsolvar_ndx(2) - Soffset) / (svar_s_avg - Soffset)
            svar_i = (scon - (svar_f * Fint + svar_s * Sint)) / Iint 
         endif

!   Mean solar cycle with no solar variability (NRLSSI2 model)
!   Averaged facular, sunspot and quiet sun terms from mean solar cycle 
!   (derived as average of Solar Cycles 13-24). This information is built
!   into coefficient terms specified by g-point elsewhere. Separate
!   scaling by spectral band is applied as defined by bndsolvar. 
!   Scale internal solar constant (svar_cprim) to requested solar constant (scon)
!   Fint is provided as the product of (svar_f_avg-Foffset) and Fint, 
!   Sint is provided as the product of (svar_s_avg-Soffset) and Sint
         if (isolvar .eq. 3) then
            svar_cprim = Fint + Sint + Iint
            if (.not. present(bndsolvar)) solvar(jpb1:jpb2) = scon / svar_cprim
            if (present(bndsolvar)) solvar(jpb1:jpb2) = bndsolvar(:) * scon / svar_cprim
            do ib = jpb1,jpb2
               svar_f_bnd(ib) = solvar(ib)
               svar_s_bnd(ib) = solvar(ib)
               svar_i_bnd(ib) = solvar(ib)
            enddo
         endif 

      endif SCON_GT_0

! Combine Earth-Sun adjustment and solar constant scaling
! when no solar variability and no solar cycle requested
      if (isolvar .lt. 0) then
         do ib = jpb1,jpb2
            adjflux(ib) = adjflx * solvar(ib)
         enddo
      endif
! Define Earth-Sun adjustment when solar variability requested
      if (isolvar .ge. 0) then
         do ib = jpb1,jpb2
            adjflux(ib) = adjflx
         enddo
      endif
      
    

      if (icld.lt.0.or.icld.gt.4) icld = 1
      
      
    ! determine cloud profile
    cldflag=0
    do iplon = 1, gncol
        if (any(gcld(iplon,:) > 0)) cldflag(iplon)=1
    end do



    ! build profile separation
    cols = 0
    cole = 0


    do iplon = 1, gncol
        if (cldflag(iplon)==1) then
            cole=cole+1
            profi(cole) = iplon
        else
            cols=cols+1
            profic(cols) = iplon
        end if
    end do
    

if (icld==4) then
    call TABULATE_XCW_BETA
end if

        
    
!$acc data copyout(swuflxc, swdflxc, swuflx, swdflx, swnflxc, swnflx, swhrc, swhr) &
!$acc create(laytrop, layswtch, laylow, jp, jt, jt1, &
!$acc co2mult, colch4, colco2, colh2o, colmol, coln2o, &
!$acc colo2, colo3, fac00, fac01, fac10, fac11, &
!$acc selffac, selffrac, indself, forfac, forfrac, indfor, &
!$acc zbbfu, zbbfd, zbbcu, zbbcd,zbbfddir, zbbcddir, zuvfd, zuvcd, zuvfddir, &
!$acc zuvcddir, znifd, znicd, znifddir,znicddir, &
!$acc cldfmcl, ciwpmcl, clwpmcl,  &
!$acc taormc, taucmc, ssacmc, asmcmc, fsfcmc) &
!$acc deviceptr(zref,zrefo,zrefd,zrefdo,&
!$acc ztauo,ztdbt,&
!$acc ztra,ztrao,ztrad,ztrado,&
!$acc zfd,zfu,zdbt,zgco,&
!$acc zomco,zrdnd,ztaug, ztaur,zsflxzen,ssi)&
!$acc create(ciwp, clwp, cld, tauc, ssac, asmc, fsfc, rei, rel, rdl, adl) &
!$acc create(play, tlay, plev, tlev, tsfc, cldflag, coszen, swdflx_at_top) &
!$acc create(coldry, wkl, znirr,znirf,zparr,zparf,zuvrr,zuvrf) &
!$acc create(extliq1, ssaliq1, asyliq1, extice2, ssaice2, asyice2) &
!$acc create(extice3, ssaice3, asyice3, fdlice3, abari, bbari, cbari, dbari, ebari, fbari) &
!$acc create(taua, asya, omga,gtauaer,gssaaer,gasmaer, zm, alat) &
!$acc create(CDF, CDF2, CDF3, alpha) &
!$acc copyin(wavenum2, ngb) &
!$acc copyin(tref, preflog, albdif, albdir, cossza)&
!$acc copyin(icxa, adjflux, nspa, nspb)&
!$acc copyin(svar_f, svar_s, svar_i)&
!$acc copyin(svar_f_bnd, svar_s_bnd, svar_i_bnd)&
!$acc copyin(kao16,kbo16,selfrefo16,forrefo16,sfluxrefo16)&
!$acc copyin(ka16,kb16,selfref16,forref16,sfluxref16)&
!$acc copyin(irradnceo16,facbrghto16,snsptdrko16)&
!$acc copyin(kao17,kbo17,selfrefo17,forrefo17,sfluxrefo17)&
!$acc copyin(ka17,kb17,selfref17,forref17,sfluxref17)&
!$acc copyin(irradnceo17,facbrghto17,snsptdrko17)&
!$acc copyin(kao18,kbo18,selfrefo18,forrefo18,sfluxrefo18)&
!$acc copyin(ka18,kb18,selfref18,forref18,sfluxref18)&
!$acc copyin(irradnceo18,facbrghto18,snsptdrko18)&
!$acc copyin(kao19,kbo19,selfrefo19,forrefo19,sfluxrefo19)&
!$acc copyin(ka19,kb19,selfref19,forref19,sfluxref19)&
!$acc copyin(irradnceo19,facbrghto19,snsptdrko19)&
!$acc copyin(kao20,kbo20,selfrefo20,forrefo20,sfluxrefo20,absch4o20)&
!$acc copyin(ka20,kb20,selfref20,forref20,sfluxref20,absch420)&
!$acc copyin(irradnceo20,facbrghto20,snsptdrko20)&
!$acc copyin(kao21,kbo21,selfrefo21,forrefo21,sfluxrefo21)&
!$acc copyin(ka21,kb21,selfref21,forref21,sfluxref21)&
!$acc copyin(irradnceo21,facbrghto21,snsptdrko21)&
!$acc copyin(kao22,kbo22,selfrefo22,forrefo22,sfluxrefo22)&
!$acc copyin(ka22,kb22,selfref22,forref22,sfluxref22)&
!$acc copyin(irradnceo22,facbrghto22,snsptdrko22)&
!$acc copyin(kao23,selfrefo23,forrefo23,sfluxrefo23,raylo23)&
!$acc copyin(ka23,selfref23,forref23,sfluxref23,rayl23)&
!$acc copyin(irradnceo23,facbrghto23,snsptdrko23)&
!$acc copyin(kao24,kbo24,selfrefo24,forrefo24,sfluxrefo24,abso3ao24,abso3bo24,raylao24,raylbo24)&
!$acc copyin(ka24,kb24,selfref24,forref24,sfluxref24,abso3a24,abso3b24,rayla24,raylb24)&
!$acc copyin(irradnceo24,facbrghto24,snsptdrko24)&
!$acc copyin(kao25,sfluxrefo25,abso3ao25,abso3bo25,raylo25)&
!$acc copyin(ka25,sfluxref25,abso3a25,abso3b25,rayl25)&
!$acc copyin(irradnceo25,facbrghto25,snsptdrko25)&
!$acc copyin(sfluxrefo26)&
!$acc copyin(sfluxref26,gzm,galat)&
!$acc copyin(irradnceo26,facbrghto26,snsptdrko26)&
!$acc copyin(kao27,kbo27,sfluxrefo27, raylo27)&
!$acc copyin(ka27,kb27,sfluxref27, rayl27)&
!$acc copyin(irradnceo27,facbrghto27,snsptdrko27)&
!$acc copyin(kao28,kbo28,sfluxrefo28)&
!$acc copyin(ka28,kb28,sfluxref28,gtauc, gssac, gasmc, gfsfc)&
!$acc copyin(irradnceo28,facbrghto28,snsptdrko28)&
!$acc copyin(kao29,kbo29,selfrefo29,forrefo29,sfluxrefo29,absh2oo29,absco2o29)&
!$acc copyin(ka29,kb29,selfref29,forref29,sfluxref29,absh2o29,absco229)&
!$acc copyin(irradnceo29,facbrghto29,snsptdrko29)&
!$acc copyin(gh2ovmr, gco2vmr, go3vmr, gn2ovmr, gch4vmr, go2vmr)&
!$acc copyin(gcld, gciwp, gclwp, grei, grel, gplay, gplev, gtlay, gtlev, gtsfc)&
!$acc copyin(gasdir, galdir, gasdif, galdif,profi,profic,gcoszen)&
!$acc copyout(nirr,nirf,parr,parf,uvrr,uvrf)


!$acc data copyin(XCW) if(icld==4)

!$acc update device(extliq1, ssaliq1, asyliq1, extice2, ssaice2, asyice2) &
!$acc device(extice3, ssaice3, asyice3, fdlice3, abari, bbari, cbari, dbari, ebari, fbari) &
!$acc device(preflog)


ncolc = cols
ncolb = cole

npartc = ceiling( real(ncolc) / real(ncol) )
npartb = ceiling( real(ncolb) / real(ncol) )

!$acc kernels    
    cldfmcl = 0.0
    ciwpmcl = 0.0
    clwpmcl = 0.0     
!$acc end kernels
  
      idelm = 1


!$acc kernels
taua = 0.0
asya = 0.0
omga = 1.0
!$acc end kernels

if (iaer==10) then

!$acc update device(gtauaer,gssaaer,gasmaer)

end if





      


! PARTITION LOOP ----------------------------------------------------------------------------
do cc = 1, 2

     if (cc==1) then 
         
         npart = npartc
         ncolst = ncolc
     else
        
         npart = npartb
         ncolst = ncolb
         
     end if



      do ipart = 0,npart-1
        cols = ipart * ncol + 1
        cole = (ipart + 1) * ncol
        if (cole>ncolst) cole=ncolst
        colr = cole - cols + 1

!$acc kernels            
            taormc = 0.0 
            taucmc = 0.0
            ssacmc = 1.0
            asmcmc = 0.0
            fsfcmc = 0.0
!$acc end kernels            

! Clear cases
      if (cc==1) then    
 !$acc kernels loop private(piplon)
 do iplon = 1, colr
      piplon = profic(iplon + cols - 1)

      do ib=1,8
         albdir(iplon,ib)  = galdir(piplon)
         albdif(iplon,ib)  = galdif(piplon)
         enddo
         albdir(iplon,nbndsw)  = galdir(piplon)
         albdif(iplon,nbndsw)  = galdif(piplon)
!  UV/visible bands 25-28 (10-13), 16000-50000 cm-1, 0.200-0.625 micron
     
         do ib=10,13
         albdir(iplon,ib)  = gasdir(piplon)
         albdif(iplon,ib)  = gasdif(piplon)
      enddo

!  Transition band 9, 12850-16000 cm-1, 0.625-0.778 micron, Take average, dmlee
       albdir(iplon, 9) = (gasdir(piplon)+galdir(piplon))/2.
       albdif(iplon, 9) = (gasdif(piplon)+galdif(piplon))/2.

         enddo
!$acc end kernels      

!$acc kernels 
do iplon = 1, colr

     piplon = profic(iplon + cols - 1)
    
     play(iplon,:) = gplay(piplon, 1:nlay)
     plev(iplon,:) = gplev(piplon, 1:nlay+1)
     tlay(iplon,:) = gtlay(piplon, 1:nlay)
     tlev(iplon,:) = gtlev(piplon, 1:nlay+1)
     tsfc(iplon)   = gtsfc(piplon)

               enddo
!$acc end kernels

if (iaer==10) then
    
!$acc kernels
    do iplon = 1, colr
     piplon = profic(iplon + cols - 1)
     taua(iplon, 1:nlay, :) = gtauaer(piplon, 1:nlay, :)
     asya(iplon, 1:nlay, :) = gasmaer(piplon, 1:nlay, :)
     omga(iplon, 1:nlay, :) = gssaaer(piplon, 1:nlay, :)
   
            enddo
!$acc end kernels

         endif   


     
!$acc kernels
do iplon = 1, colr
     piplon = profic(iplon + cols - 1)
     wkl(iplon,1,:) = gh2ovmr(piplon,1:nlay)
     wkl(iplon,2,:) = gco2vmr(piplon,1:nlay)
     wkl(iplon,3,:) = go3vmr(piplon,1:nlay)
     wkl(iplon,4,:) = gn2ovmr(piplon,1:nlay)
     wkl(iplon,5,:) = 0.0
     wkl(iplon,6,:) = gch4vmr(piplon,1:nlay)
     wkl(iplon,7,:) = go2vmr(piplon,1:nlay)   
     coszen(iplon)  = gcoszen(piplon)
     
     
  
   
end do
!$acc end kernels
!************** cloudy cases ***************
                  else
          
 !$acc kernels loop private(piplon)
 do iplon = 1, colr
      piplon = profi(iplon + cols - 1)
     
      do ib=1,8
         albdir(iplon,ib)  = galdir(piplon)
         albdif(iplon,ib)  = galdif(piplon)
               enddo
         albdir(iplon,nbndsw)  = galdir(piplon)
         albdif(iplon,nbndsw)  = galdif(piplon)
        !  UV/visible bands 25-28 (10-13), 16000-50000 cm-1, 0.200-0.625 micron
     
      do ib=10,13
         albdir(iplon,ib)  = gasdir(piplon)
         albdif(iplon,ib)  = gasdif(piplon)
            enddo

!  Transition band 9, 12850-16000 cm-1, 0.625-0.778 micron, Take average, dmlee
       albdir(iplon, 9) = (gasdir(piplon)+galdir(piplon))/2.
       albdif(iplon, 9) = (gasdif(piplon)+galdif(piplon))/2.

               enddo
!$acc end kernels               
          
!$acc kernels 
do iplon = 1, colr
   
     piplon = profi(iplon + cols - 1)
     
     play(iplon,:) = gplay(piplon, 1:nlay)
     plev(iplon,:) = gplev(piplon, 1:nlay+1)
     tlay(iplon,:) = gtlay(piplon, 1:nlay)
     tlev(iplon,:) = gtlev(piplon, 1:nlay+1)
     tsfc(iplon) = gtsfc(piplon)
     cld(iplon,:) = gcld(piplon, 1:nlay)
     ciwp(iplon,:) = gciwp(piplon, 1:nlay)
     clwp(iplon,:) = gclwp(piplon, 1:nlay)
     rei(iplon,:) = grei(piplon, 1:nlay) 
     rel(iplon,:) = grel(piplon, 1:nlay)
     zm(iplon,:) = gzm(piplon, 1:nlay)
     alat(iplon) = galat(piplon)
            enddo
!$acc end kernels

if (iaer==10) then

!$acc kernels    
    do iplon = 1, colr
     piplon = profi(iplon + cols - 1)
     taua(iplon, 1:nlay, :) = gtauaer(piplon, 1:nlay, :)
     asya(iplon, 1:nlay, :) = gasmaer(piplon, 1:nlay, :)
     omga(iplon, 1:nlay, :) = gssaaer(piplon, 1:nlay, :)
   
    end do
!$acc end kernels

         endif


! Copy the direct cloud optical properties over to the temp arrays
! and then onto the GPU
! We are on the CPU here


!$acc kernels 
do iplon = 1, colr
    piplon = profi(iplon + cols - 1)
     tauc(iplon, 1:nlay, :) = gtauc(piplon, 1:nlay, :)
     ssac(iplon, 1:nlay, :) = gssac(piplon, 1:nlay, :)
     asmc(iplon, 1:nlay, :) = gasmc(piplon, 1:nlay, :)
     fsfc(iplon, 1:nlay, :) = gfsfc(piplon, 1:nlay, :)
         enddo
!$acc end kernels



!$acc kernels
do iplon = 1, colr
     piplon = profi(iplon + cols - 1)
     wkl(iplon,1,:) = gh2ovmr(piplon,1:nlay)
     wkl(iplon,2,:) = gco2vmr(piplon,1:nlay)
     wkl(iplon,3,:) = go3vmr(piplon,1:nlay)
     wkl(iplon,4,:) = gn2ovmr(piplon,1:nlay)
     wkl(iplon,5,:) = 0.0
     wkl(iplon,6,:) = gch4vmr(piplon,1:nlay)
     wkl(iplon,7,:) = go2vmr(piplon,1:nlay)  
     coszen(iplon)  = gcoszen(piplon)

         enddo
!$acc end kernels
end if


!$acc kernels
do iplon = 1, colr
     cossza(iplon) = max(zepzen,coszen(iplon))
         enddo
!$acc end kernels  


!$acc kernels
  do iplon = 1,colr
      
      do l = 1,nlay

         coldry(iplon, l) = (plev(iplon, l)-plev(iplon, l+1)) * 1.e3  * avogad / &
                     (1.e2  * grav * ((1.  - wkl(iplon, 1,l)) * amd + wkl(iplon, 1,l) * amw) * &
                     (1.  + wkl(iplon, 1,l)))
    
      end do
      enddo
!$acc end kernels

!$acc kernels
  do iplon = 1,colr

      do l = 1,nlay
        do imol = 1, nmol
           wkl(iplon,imol,l) = coldry(iplon,l) * wkl(iplon,imol,l)
        end do
       end do
    end do
!$acc end kernels


   if (cc==2) then
   call mcica_sw(colr, nlay, ngptsw, icld,irng, play, &
                       cld, clwp, ciwp, tauc, ssac, asmc, fsfc, &
                       cldfmcl, clwpmcl, ciwpmcl, &
                       taucmc, ssacmc, asmcmc, fsfcmc,1,CDF, CDF2, CDF3, alpha, zm, &
                       alat, dyofyr, rdl, adl)
   end if   



   if (cc==2) then
   call cldprmc_sw(colr, nlay, inflgsw, iceflgsw, liqflgsw,  &
                         cldfmcl , ciwpmcl , clwpmcl , rei, rel, &
                         taormc, taucmc, ssacmc, asmcmc, fsfcmc)
   end if



   call setcoef_sw(colr, nlay, play , tlay , plev , tlev , tsfc , &
                        coldry , wkl , &
                         laytrop, layswtch, laylow, jp , jt , jt1 , &
                         co2mult , colch4 , colco2 , colh2o , colmol , coln2o , &
                         colo2 , colo3 , fac00 , fac01 , fac10 , fac11 , &
                         selffac , selffrac , indself , forfac , forfrac , indfor )




   call spcvmc_sw(cc,ncol, colr, nlay, istart, iend, icpr, idelm, iout, &
              play, tlay, plev, tlev, &
              tsfc, albdif, albdir, &
              cldfmcl, taucmc, asmcmc, ssacmc, taormc, &
              taua, asya, omga,cossza, coldry, adjflux, &
              isolvar, svar_f, svar_s, svar_i, &
              svar_f_bnd, svar_s_bnd, svar_i_bnd, &
              laytrop, layswtch, laylow, jp, jt, jt1, &
              co2mult, colch4, colco2, colh2o, colmol, &
              coln2o, colo2, colo3, &
              fac00, fac01, fac10, fac11, &
              selffac, selffrac, indself, forfac, forfrac, indfor, &
              zbbfd, zbbfu, zbbcd, zbbcu, zuvfd, &
              zuvcd, znifd, znicd, &
              zbbfddir, zbbcddir, zuvfddir, zuvcddir, znifddir, znicddir,&
              zgco,zomco,zrdnd,zref,zrefo,zrefd,zrefdo,ztauo,zdbt,ztdbt,&
              ztra,ztrao,ztrad,ztrado,zfd,zfu,ztaug, ztaur, zsflxzen, ssi,&
              znirr,znirf,zparr,zparf,zuvrr,zuvrf)




! Transfer up and down, clear and total sky fluxes to output arrays.
! Vertical indexing goes from bottom to top; reverse here for GCM if necessary.



if (cc==1) then
!$acc kernels loop independent
    do iplon = 1, colr
         piplon = profic(iplon + cols - 1)
        
         do i = 1, nlay+1


            swuflxc(piplon,i) = zbbcu(iplon,i) 
            swdflxc(piplon,i) = zbbcd(iplon,i) 
            swuflx(piplon,i) = zbbfu(iplon,i) 
            swdflx(piplon,i) = zbbfd(iplon,i) 

         enddo

!  Total and clear sky net fluxes

         do i = 1, nlay+1
            swnflxc(iplon,i)  = swdflxc(piplon,i) - swuflxc(piplon,i)
            swnflx(iplon,i)  = swdflx(piplon,i) - swuflx(piplon,i)
         enddo

!  Total and clear sky heating rates

         do i = 1, nlay
            zdpgcp = heatfac / (plev(iplon, i) - plev(iplon, i+1))
            swhrc(piplon,i) = (swnflxc(iplon,i+1)  - swnflxc(iplon,i) ) * zdpgcp
            swhr(piplon,i) = (swnflx(iplon,i+1)  - swnflx(iplon,i) ) * zdpgcp
         enddo
         swhrc(piplon,nlay) = 0. 
         swhr(piplon,nlay) = 0. 

! End longitude loop
      enddo
!$acc end kernels 

!$acc kernels loop independent
do iplon = 1, colr
         piplon = profic(iplon + cols - 1)
         nirr(piplon) = znirr(iplon)
         nirf(piplon) = znirf(iplon) - znirr(iplon)
         parr(piplon) = zparr(iplon)
         parf(piplon) = zparf(iplon) - zparr(iplon)
         uvrr(piplon) = zuvrr(iplon)
         uvrf(piplon) = zuvrf(iplon) - zuvrr(iplon)

end do
!$acc end kernels 
            else
!$acc kernels loop independent
    do iplon = 1, colr
         piplon = profi(iplon + cols - 1)

         do i = 1, nlay+1
             

            swuflxc(piplon,i) = zbbcu(iplon,i) 
            swdflxc(piplon,i) = zbbcd(iplon,i) 
            swuflx(piplon,i) = zbbfu(iplon,i) 
            swdflx(piplon,i) = zbbfd(iplon,i) 

            enddo

!  Total and clear sky net fluxes

         do i = 1, nlay+1
            swnflxc(iplon,i)  = swdflxc(piplon,i) - swuflxc(piplon,i)
            swnflx(iplon,i)  = swdflx(piplon,i) - swuflx(piplon,i)
         enddo

!  Total and clear sky heating rates

         do i = 1, nlay
            zdpgcp = heatfac / (plev(iplon, i) - plev(iplon, i+1))
            swhrc(piplon,i) = (swnflxc(iplon,i+1)  - swnflxc(iplon,i) ) * zdpgcp
            swhr(piplon,i) = (swnflx(iplon,i+1)  - swnflx(iplon,i) ) * zdpgcp
         enddo
         swhrc(piplon,nlay) = 0. 
         swhr(piplon,nlay) = 0. 

! End longitude loop
            enddo
!$acc end kernels 

!$acc kernels loop independent
do iplon = 1, colr
         piplon = profi(iplon + cols - 1)
         nirr(piplon) = znirr(iplon)
         nirf(piplon) = znirf(iplon) - znirr(iplon)
         parr(piplon) = zparr(iplon)
         parf(piplon) = zparf(iplon) - zparr(iplon)
         uvrr(piplon) = zuvrr(iplon)
         uvrf(piplon) = zuvrf(iplon) - zuvrr(iplon)

         enddo
!$acc end kernels 
      endif




      enddo


     
            enddo

! If the user requested 'normalized' fluxes, then here we
! divide the fluxes by the solar constant divided by coszen

! MAT This requires only lit points passed in

if (normFlx==1) then


!$acc kernels
   swdflx_at_top(:) = max(swdflx(:,nlay+1),1e-7)

   do k = 1, nlay+1
      swuflxc(:,k)=swuflxc(:,k)/swdflx_at_top(:)
      swdflxc(:,k)=swdflxc(:,k)/swdflx_at_top(:)
      swuflx (:,k)=swuflx (:,k)/swdflx_at_top(:)
      swdflx (:,k)=swdflx (:,k)/swdflx_at_top(:)
   enddo

   nirr(:)=nirr(:)/swdflx_at_top(:)
   nirf(:)=nirf(:)/swdflx_at_top(:)
   parr(:)=parr(:)/swdflx_at_top(:)
   parf(:)=parf(:)/swdflx_at_top(:)
   uvrr(:)=uvrr(:)/swdflx_at_top(:)
   uvrf(:)=uvrf(:)/swdflx_at_top(:)
!$acc end kernels
endif

! end of data statement for xcw with icld==4
!$acc end data

!$acc end data
      
end subroutine rrtmg_sw_sub

!*************************************************************************
      real  function earth_sun(idn)
!*************************************************************************
!
!  Purpose: Function to calculate the correction factor of Earth's orbit
!  for current day of the year

!  idn        : Day of the year
!  earth_sun  : square of the ratio of mean to actual Earth-Sun distance

! ------- Modules -------

      use rrsw_con, only : pi

      integer , intent(in) :: idn

      real  :: gamma

      gamma = 2. *pi*(idn-1)/365. 

! Use Iqbal's equation 1.2.1

      earth_sun = 1.000110  + .034221  * cos(gamma) + .001289  * sin(gamma) + &
                   .000719  * cos(2. *gamma) + .000077  * sin(2. *gamma)

      end function earth_sun

      end module rrtmg_sw_rad


