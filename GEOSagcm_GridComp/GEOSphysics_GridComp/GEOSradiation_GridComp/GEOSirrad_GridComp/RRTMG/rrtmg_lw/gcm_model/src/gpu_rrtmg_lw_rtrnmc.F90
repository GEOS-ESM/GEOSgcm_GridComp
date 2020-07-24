! (dmb 2012) This is the GPU version of the rtrnmc subroutine.  This has been greatly
! modified to be efficiently run on the GPU.  Originally, there was a g-point loop within
! this subroutine to perform the summation of the fluxes over the g-points.  This has been
! modified so that this subroutine can be run in parallel across the g-points.  This was
! absolutely critical because of two reasons.
! 1. For a relatively low number of profiles, there wouldn't be enough threads to keep
!    the GPU busy enough to run at full potential.  As a result of this, this subroutine
!    would end up being a bottleneck.
! 2. The memory access for the GPU arrays would be innefient because there would be very
!    little coalescing which is critical for obtaining optimal performance.
#include "_gpudef.inc"


      module gpu_rrtmg_lw_rtrnmc

!  --------------------------------------------------------------------------
! |                                                                          |
! |  Copyright 2002-2009, Atmospheric & Environmental Research, Inc. (AER).  |
! |  This software may be used, copied, or redistributed as long as it is    |
! |  not sold and this copyright notice is reproduced on each copy made.     |
! |  This model is provided as is without any express or implied warranties. |
! |                       (http://www.rtweb.aer.com/)                        |
! |                                                                          |
!  --------------------------------------------------------------------------

! --------- Modules ----------

       !use parkind, only : im => kind , rb => kind 
      use parrrtm, only : mg, nbndlw, ngptlw
      use rrlw_con, only: fluxfac, heatfac
      use rrlw_tbl, only: bpade, tblint, tau_tbl, exp_tbl, tfn_tbl

#ifdef _CUDA
      use cudafor
#endif      
   
    
      implicit none 
      
      integer(kind=4), parameter :: ntbl = 10000
      integer  _gpucon :: ngsd(nbndlw)      

      ! (dmb 2012) I moved most GPU variables so that they are module level variables.
      ! PGI Fortran seems to sometimes have trouble passing arrays into kernels correctly.
      ! Using module level variables bypasses this issue and allows for cleaner code.
      ! Atmosphere
      
      real , allocatable _gpudev :: taucmcd(:,:,:)
   
      real , allocatable _gpudev, dimension(:,:) :: pzd            ! level (interface) pressures (hPa, mb)

                                                      !    Dimensions: (0:nlayers)
      real , allocatable _gpudev, dimension(:) :: pwvcmd              ! precipitable water vapor (cm)
      real , allocatable _gpudev, dimension(:,:) :: semissd          ! lw surface emissivity
                                                      !    Dimensions: (nbndlw)
      real , allocatable _gpudev, dimension(:,:,:) :: planklayd      ! 
                                                      !    Dimensions: (nlayers,nbndlw)
      real , allocatable _gpudev, dimension(:,:,:) :: planklevd     ! 
                                                      !    Dimensions: (0:nlayers,nbndlw)
      real, allocatable _gpudev, dimension(:,:) :: plankbndd        ! 
                                                      !    Dimensions: (nbndlw)
   
      real , allocatable _gpudev :: gurad(:,:,:)     ! upward longwave flux (w/m2)
      real , allocatable _gpudev :: gdrad(:,:,:)     ! downward longwave flux (w/m2)
      real , allocatable _gpudev :: gclrurad(:,:,:)     ! clear sky upward longwave flux (w/m2)
      real , allocatable _gpudev :: gclrdrad(:,:,:)     ! clear sky downward longwave flux (w/m2)

      real, allocatable  _gpudev :: gdtotuflux_dtd(:,:,:) ! change in upward longwave flux (w/m2/k)
                                              ! with respect to surface temperature

      real, allocatable  _gpudev :: gdtotuclfl_dtd(:,:,:) ! change in clear sky upward longwave flux (w/m2/k)
                                              ! with respect to surface temperature
  
! Clouds
      integer  _gpudev :: idrvd            ! flag for calculation of dF/dt from 
                                                      ! Planck derivative [0=off, 1=on]
      real  _gpucon :: bpaded
      real  _gpucon :: heatfacd
      real  _gpucon :: fluxfacd
      real  _gpucon :: a0d(nbndlw), a1d(nbndlw), a2d(nbndlw)
      integer  _gpucon :: delwaved(nbndlw)
      real , allocatable _gpudev :: totufluxd(:,:)     ! upward longwave flux (w/m2)
      real , allocatable _gpudev :: totdfluxd(:,:)     ! downward longwave flux (w/m2)
      real , allocatable _gpudev :: fnetd(:,:)         ! net longwave flux (w/m2)
      real , allocatable _gpudev :: htrd(:,:)          ! longwave heating rate (k/day)
      real , allocatable _gpudev :: totuclfld(:,:)     ! clear sky upward longwave flux (w/m2)
      real , allocatable _gpudev :: totdclfld(:,:)     ! clear sky downward longwave flux (w/m2)
      real , allocatable _gpudev :: fnetcd(:,:)        ! clear sky net longwave flux (w/m2)
      real , allocatable _gpudev :: htrcd(:,:)         ! clear sky longwave heating rate (k/day)
      real , allocatable _gpudev :: dtotuflux_dtd(:,:) ! change in upward longwave flux (w/m2/k)
                                              ! with respect to surface temperature
      real , allocatable _gpudev :: dtotuclfl_dtd(:,:) ! change in clear sky upward longwave flux (w/m2/k)
                                              ! with respect to surface temperature
      real , allocatable _gpudev :: dplankbnd_dtd(:,:) 
     
      ! TOA OLR in bands 9-11 and their derivatives with surface temp
      real , allocatable _gpudev :: olrb09d(:)     ! (W/m2)
      real , allocatable _gpudev :: dolrb09_dtd(:) ! (W/m2/K)
      real , allocatable _gpudev :: olrb10d(:)     ! (W/m2)
      real , allocatable _gpudev :: dolrb10_dtd(:) ! (W/m2/K)
      real , allocatable _gpudev :: olrb11d(:)     ! (W/m2)
      real , allocatable _gpudev :: dolrb11_dtd(:) ! (W/m2/K)

      contains

!-----------------------------------------------------------------------------
      ! TODO add _gpuker 
      _gpuker subroutine  rtrnzero(ncol, nlayers)
!-----------------------------------------------------------------------------

! MAT This is a new routine that moves some of the zeroing out of rtrnmcg. The
!     odd loop structure that is good for GPUs is bad for CPUs. A move to
!     array syntax on CPUs allows the compiler to use fancy techniques to save
!     a rather significant amount of time.

! ----- Input -----

      integer(kind=4), value, intent(in) :: nlayers        ! total number of layers
      integer(kind=4), value, intent(in) :: ncol           ! total number of columns
     
      integer :: iplon, igc, lay

      ! (dmb 2012) Here we compute the index for the column and band dimensions

#ifdef _CUDA
      iplon = (blockidx%x-1) * blockdim%x + threadidx%x
      igc   = (blockidx%y-1) * blockdim%y + threadidx%y

      ! (dmb 2012) Make sure that the column and bands are within the proper ranges
      if (iplon <= ncol .and. igc<=ngptlw) then
    
            gurad(iplon, igc, 0) = 0.0 
            gdrad(iplon, igc, 0) = 0.0 
            !totuflux(iplon,igc,0) = 0.0 
            !totdflux(iplon,igc,0) = 0.0 
            gclrurad(iplon, igc, 0) = 0.0 
            gclrdrad(iplon, igc, 0) = 0.0 
            !totuclfl(iplon,igc,0) = 0.0 
            !totdclfl(iplon,igc,0) = 0.0 
            if (idrvd .eq. 1) then
               gdtotuflux_dtd(iplon,igc,0) = 0.0 
               gdtotuclfl_dtd(iplon,igc,0) = 0.0 
            endif

            do lay = 1, nlayers
               gurad(iplon, igc, lay) = 0.0 
               gdrad(iplon, igc, lay) = 0.0 
               gclrurad(iplon, igc, lay) = 0.0 
               gclrdrad(iplon, igc, lay) = 0.0 
               
               ! (dmb 2012) I removed the band loop here because it was terribly inefficient
               ! I now set the required variables outside of the kernel

               if (idrvd .eq. 1) then
                  gdtotuflux_dtd(iplon,igc,lay) = 0.0 
                  gdtotuclfl_dtd(iplon,igc,lay) = 0.0 
               endif
            enddo
      endif

#else
      gurad = 0.0 
      gdrad = 0.0 
      gclrurad = 0.0 
      gclrdrad = 0.0 
      if (idrvd .eq. 1) then
         gdtotuflux_dtd = 0.0 
         gdtotuclfl_dtd = 0.0 
      endif
#endif

       end subroutine rtrnzero

    !-----------------------------------------------------------------------------
      ! TODO add _gpuker 
      _gpuker subroutine  rtrnmcg(ncol, nlayers, istart, iend, iout, &
                        ngb,icldlyr, taug, fracsd, cldfmcd)
!-----------------------------------------------------------------------------
!
!  Original version:   E. J. Mlawer, et al. RRTM_V3.0
!  Revision for GCMs:  Michael J. Iacono; October, 2002
!  Revision for F90:  Michael J. Iacono; June, 2006
!  Revision for dFdT option: M. J. Iacono and E. J. Mlawer, November 2009
!
!  This program calculates the upward fluxes, downward fluxes, and
!  heating rates for an arbitrary clear or cloudy atmosphere.  The input
!  to this program is the atmospheric profile, all Planck function
!  information, and the cloud fraction by layer.  A variable diffusivity 
!  angle (SECDIFF) is used for the angle integration.  Bands 2-3 and 5-9 
!  use a value for SECDIFF that varies from 1.50 to 1.80 as a function of 
!  the column water vapor, and other bands use a value of 1.66.  The Gaussian 
!  weight appropriate to this angle (WTDIFF=0.5) is applied here.  Note that 
!  use of the emissivity angle for the flux integration can cause errors of 
!  1 to 4 W/m2 within cloudy layers.  
!  Clouds are treated with the McICA stochastic approach and maximum-random
!  cloud overlap. 
!  This subroutine also provides the optional capability to calculate
!  the derivative of upward flux respect to surface temperature using
!  the pre-tabulated derivative of the Planck function with respect to 
!  temperature integrated over each spectral band.
!***************************************************************************

! ------- Declarations -------

! ----- Input -----

       integer(kind=4), value, intent(in) :: nlayers         ! total number of layers
      integer(kind=4), value, intent(in) :: ncol           ! total number of columns
      integer(kind=4), value, intent(in) :: istart          ! beginning band of calculation
      integer(kind=4), value, intent(in) :: iend            ! ending band of calculation
      integer(kind=4), value, intent(in) :: iout            ! output option flag
      integer , intent(in) :: ngb(:)                ! band index
     
      integer , intent(in) :: icldlyr(:,:)
      real  _gpudev :: taug(:,:,:)
      real  _gpudev :: fracsd(:,:,:)
      real  _gpudev :: cldfmcd(:,:,:)
     
!GPU These are needed to allocate space for temporary automatic local
!GPU arrays on the GPU. On the CPU this ifndef converts them to nlayers

#ifndef GPU_MAXLEVS
#define GPU_MAXLEVS nlayers
#endif

   
         
   ! ----- Local -----
! Declarations for radiative transfer
   
      real  :: atot(GPU_MAXLEVS)
      real  :: atrans(GPU_MAXLEVS)
      real  :: bbugas(GPU_MAXLEVS)
      real  :: bbutot(GPU_MAXLEVS)
     
      real  :: uflux(0:GPU_MAXLEVS)
      real  :: dflux(0:GPU_MAXLEVS)
      real  :: uclfl(0:GPU_MAXLEVS)
      real  :: dclfl(0:GPU_MAXLEVS)
    
      real  :: odclds
      real  :: efclfracs
      real  :: absclds

      real  :: secdiff                 ! secant of diffusivity angle
      real  :: transcld, radld, radclrd, plfrac, blay, dplankup, dplankdn
      real  :: odepth, odtot, odepth_rec, odtot_rec, gassrc
      real  :: tblind, tfactot, bbd, bbdtot, tfacgas, transc, tausfac
      real  :: rad0, reflect, radlu, radclru
      real  :: d_rad0_dt, d_radlu_dt, d_radclru_dt

   
      integer  :: ibnd, ib, lay, lev, l, ig  ! loop indices
      integer  :: igc                               ! g-point interval counter
      integer  :: iclddn                            ! flag for cloud in down path
      integer  :: ittot, itgas, itr                 ! lookup table indices
   

  
! ------- Definitions -------
! input
!    nlayers                      ! number of model layers
!    ngptlw                       ! total number of g-point subintervals
!    nbndlw                       ! number of longwave spectral bands
!    ncbands                      ! number of spectral bands for clouds
!    secdiff                      ! diffusivity angle
!    wtdiff                       ! weight for radiance to flux conversion
!    pavel                        ! layer pressures (mb)
!    pz                           ! level (interface) pressures (mb)
!    tavel                        ! layer temperatures (k)
!    tz                           ! level (interface) temperatures(mb)
!    tbound                       ! surface temperature (k)
!    cldfrac                      ! layer cloud fraction
!    taucloud                     ! layer cloud optical depth
!    itr                          ! integer look-up table index
!    icldlyr                      ! flag for cloudy layers
!    iclddn                       ! flag for cloud in column at any layer
!    semiss                       ! surface emissivities for each band
!    reflect                      ! surface reflectance
!    bpade                        ! 1/(pade constant)
!    tau_tbl                      ! clear sky optical depth look-up table
!    exp_tbl                      ! exponential look-up table for transmittance
!    tfn_tbl                      ! tau transition function look-up table

! local
!    atrans                       ! gaseous absorptivity
!    abscld                       ! cloud absorptivity
!    atot                         ! combined gaseous and cloud absorptivity
!    odclr                        ! clear sky (gaseous) optical depth
!    odcld                        ! cloud optical depth
!    odtot                        ! optical depth of gas and cloud
!    tfacgas                      ! gas-only pade factor, used for planck fn
!    tfactot                      ! gas and cloud pade factor, used for planck fn
!    bbdgas                       ! gas-only planck function for downward rt
!    bbugas                       ! gas-only planck function for upward rt
!    bbdtot                       ! gas and cloud planck function for downward rt
!    bbutot                       ! gas and cloud planck function for upward calc.
!    gassrc                       ! source radiance due to gas only
!    efclfrac                     ! effective cloud fraction
!    radlu                        ! spectrally summed upward radiance 
!    radclru                      ! spectrally summed clear sky upward radiance 
!    urad                         ! upward radiance by layer
!    clrurad                      ! clear sky upward radiance by layer
!    radld                        ! spectrally summed downward radiance 
!    radclrd                      ! spectrally summed clear sky downward radiance 
!    drad                         ! downward radiance by layer
!    clrdrad                      ! clear sky downward radiance by layer
!    d_radlu_dt                   ! spectrally summed upward radiance 
!    d_radclru_dt                 ! spectrally summed clear sky upward radiance 
!    d_urad_dt                    ! upward radiance by layer
!    d_clrurad_dt                 ! clear sky upward radiance by layer

! output
!    totuflux                     ! upward longwave flux (w/m2)
!    totdflux                     ! downward longwave flux (w/m2)
!    fnet                         ! net longwave flux (w/m2)
!    htr                          ! longwave heating rate (k/day)
!    totuclfl                     ! clear sky upward longwave flux (w/m2)
!    totdclfl                     ! clear sky downward longwave flux (w/m2)
!    fnetc                        ! clear sky net longwave flux (w/m2)
!    htrc                         ! clear sky longwave heating rate (k/day)
!    dtotuflux_dt                 ! change in upward longwave flux (w/m2/k)
!                                 ! with respect to surface temperature
!    dtotuclfl_dt                 ! change in clear sky upward longwave flux (w/m2/k)
!    olrb10                       ! TOA OLR in band 10
!    dolrb10_dt                   !   and its derivative wrty surface temp
!    
   



  
! This secant and weight corresponds to the standard diffusivity 
! angle.  This initial value is redefined below for some bands.
     real , parameter :: wtdiff = 0.5       ! epsilon
     real , parameter :: rec_6 = 0.166667   ! minimum value for cloud quantities

! Reset diffusivity angle for Bands 2-3 and 5-9 to vary (between 1.50
! and 1.80) as a function of total column water vapor.  The function
! has been defined to minimize flux and cooling rate errors in these bands
! over a wide range of precipitable water values.


      
      integer :: iplon
      real :: bbb
      ! (dmb 2012) Here we compute the index for the column and band dimensions
#ifdef _CUDA
      iplon = (blockidx%x-1) * blockdim%x + threadidx%x
      igc = (blockidx%y-1) * blockdim%y + threadidx%y
       ! (dmb 2012) Make sure that the column and bands are within the proper ranges
      if (iplon <= ncol .and. igc<=ngptlw) then
    
#else
     do iplon = 1, ncol
        do igc = 1, ngptlw 

#endif


       
      ibnd = ngb(igc)
      
         if (ibnd.eq.1 .or. ibnd.eq.4 .or. ibnd.ge.10) then
           secdiff = 1.66 
         else
           secdiff = a0d(ibnd) + a1d(ibnd)*exp(a2d(ibnd)*pwvcmd(iplon))
           if (secdiff .gt. 1.80 ) secdiff = 1.80 
           if (secdiff .lt. 1.50 ) secdiff = 1.50 
         endif
    
! Radiative transfer starts here.
         radld = 0. 
         radclrd = 0. 
         iclddn = 0

! Downward radiative transfer loop.  

         do lev = nlayers, 1, -1
               plfrac = fracsd(iplon,lev,igc)
               blay = planklayd(iplon,lev,ibnd)
               dplankup = planklevd(iplon,lev,ibnd) - blay
               dplankdn = planklevd(iplon,lev-1,ibnd) - blay
               odepth = secdiff * taug(iplon,lev,igc)
               if (odepth .lt. 0.0 ) odepth = 0.0 
!  Cloudy layer
               if (icldlyr(iplon, lev).eq.1) then
                  iclddn = 1
                  ! (dmb 2012) Here instead of using the lookup tables to compute 
                  ! the optical depth and related quantities, I compute them on the 
                  ! fly because this is actually much more efficient on the GPU.
                  odclds = secdiff * taucmcd(iplon,igc,lev)
                  absclds = 1.  - exp(-odclds)
                  efclfracs = absclds * cldfmcd(iplon, igc,lev)
                  odtot = odepth + odclds
                
#ifdef _CUDA
                     tblind = odepth/(bpaded+odepth)
                     itgas = tblint*tblind+0.5 
                     bbb = itgas / float(tblint)
                     odepth = bpaded * bbb / (1.  - bbb)

                     atrans(lev) = exp( -odepth)
                     atrans(lev) = 1  -atrans(lev)
                     ! (dmb 2012) Compute tfacgas on the fly.  Even though this is an expensive operation,
                     ! it is more efficient to do the calculation within the kernel on the GPU. 
                     if (odepth < 0.06) then
                     tfacgas = odepth/6. 
                     else
                     tfacgas = 1. -2. *((1. /odepth)-((1.  - atrans(lev))/(atrans(lev))))
                     endif
                     gassrc = atrans(lev) * plfrac * (blay + tfacgas*dplankdn)

                     odtot = odepth + odclds
                     tblind = odtot/(bpaded+odtot)
                     ittot = tblint*tblind + 0.5 
                     bbb = ittot / float(tblint)
                     bbb = bpaded * bbb / (1.  - bbb)
                     atot(lev) = 1.  - exp(-bbb)
                     if (bbb < 0.06) then
                     tfactot = bbb/6. 
                     else
                     tfactot = 1. -2. *((1. /bbb)-((1-atot(lev))/(atot(lev))))
                     endif
                     bbdtot = plfrac * (blay + tfactot*dplankdn)
                     bbd = plfrac*(blay+tfacgas*dplankdn)
#else
                     tblind = odepth/(bpade+odepth)
                     itgas = tblint*tblind+0.5 
                     odepth = tau_tbl(itgas)
                     atrans(lev) = 1.  - exp_tbl(itgas)
                     tfacgas = tfn_tbl(itgas)
                     gassrc = atrans(lev) * plfrac * (blay + tfacgas*dplankdn)

                     odtot = odepth + odclds
                     tblind = odtot/(bpade+odtot)
                     ittot = tblint*tblind + 0.5 
                     tfactot = tfn_tbl(ittot)
                     bbdtot = plfrac * (blay + tfactot*dplankdn)
                     bbd = plfrac*(blay+tfacgas*dplankdn)
                     atot(lev) = 1.  - exp_tbl(ittot)
#endif

                    radld = radld - radld * (atrans(lev) + &
                    efclfracs * (1.  - atrans(lev))) + &
                    gassrc + cldfmcd(iplon, igc,lev) * &
                    (bbdtot * atot(lev) - gassrc)
                  gdrad(iplon, igc, lev-1) = gdrad(iplon, igc, lev-1) + radld
                  bbugas(lev) = plfrac * (blay + tfacgas * dplankup)
                  bbutot(lev) = plfrac * (blay + tfactot * dplankup)
              

                 
!  Clear layer
               else

#ifdef _CUDA
                     tblind = odepth/(bpaded+odepth)
                     itr = tblint*tblind+0.5 
                     ! (dmb 2012) Compute the atrans and related values on the fly instead
                     ! of using the lookup tables.
                     bbb = itr/float(tblint)
                     bbb = bpaded * bbb / (1.  - bbb)
                     transc = exp( -bbb )
                     if (transc < 1.e-20 ) transc = 1.e-20 
                     atrans(lev) = 1. -transc

                     if (bbb < 0.06 ) then
                     tausfac = bbb/6. 
                     else
                     tausfac = 1. -2. *((1. /bbb)-(transc/(1.-transc)))
                     endif 


                     bbd = plfrac*(blay+tausfac*dplankdn)
                     bbugas(lev) = plfrac * (blay + tausfac * dplankup)
                      radld = radld + (bbd-radld)*atrans(lev)
#else
                
                      tblind = odepth/(bpade+odepth)
                     itr = tblint*tblind+0.5 
                     transc = exp_tbl(itr)
                     atrans(lev) = 1. -transc
                     tausfac = tfn_tbl(itr)
                     bbd = plfrac*(blay+tausfac*dplankdn)
                     bbugas(lev) = plfrac * (blay + tausfac * dplankup)
                      radld = radld + (bbd-radld)*atrans(lev)
                  
               
#endif
                  gdrad(iplon, igc, lev-1) = gdrad(iplon, igc, lev-1) + radld
                endif
!  Set clear sky stream to total sky stream as long as layers
!  remain clear.  Streams diverge when a cloud is reached (iclddn=1),
!  and clear sky stream must be computed separately from that point.
                  if (iclddn.eq.1) then
                     radclrd = radclrd + (bbd-radclrd) * atrans(lev) 
                     ! (dmb 2012) Rather than summing up the results and then computing the 
                     ! total fluxes, I store the g-point specific values in GPU arrays to be 
                     ! summed up later in a new kernel.  This ensures that we can parallelize 
                     ! across enough dimensions so that the GPU remains busy.
                     gclrdrad(iplon, igc, lev-1) = gclrdrad(iplon, igc, lev-1) + radclrd
                  else
                     radclrd = radld
                     gclrdrad(iplon, igc, lev-1) = gdrad(iplon, igc, lev-1)
                  endif
            enddo

! Spectral emissivity & reflectance
!  Include the contribution of spectrally varying longwave emissivity
!  and reflection from the surface to the upward radiative transfer.
!  Note: Spectral and Lambertian reflection are identical for the
!  diffusivity angle flux integration used here.
!  Note: The emissivity is applied to plankbnd and dplankbnd_dt when 
!  they are defined in subroutine setcoef. 
    
         rad0 = fracsd(iplon,1,igc) * plankbndd(iplon,ibnd)
         if (idrvd .eq. 1) then
            d_rad0_dt = fracsd(iplon,1,igc) * dplankbnd_dtd(iplon,ibnd)
         endif

!  Add in specular reflection of surface downward radiance.
         reflect = 1.  - semissd(iplon,ibnd)
         radlu = rad0 + reflect * radld
         radclru = rad0 + reflect * radclrd


! Upward radiative transfer loop.
         gurad(iplon, igc, 0) = gurad(iplon, igc, 0) + radlu
         gclrurad(iplon, igc, 0) = gclrurad(iplon, igc, 0) + radclru
         if (idrvd .eq. 1) then
            ! (dmb 2012) For the optional derivatives, we likewise store the intermediate
            ! g-point values for later summation.
            d_radlu_dt = d_rad0_dt
            gdtotuflux_dtd(iplon, igc, 0) = d_radlu_dt
            d_radclru_dt = d_rad0_dt
            gdtotuclfl_dtd(iplon, igc, 0) = d_radclru_dt
         endif

         do lev = 1, nlayers
!  Cloudy layer
            if (icldlyr(iplon, lev) .eq. 1) then
               gassrc = bbugas(lev) * atrans(lev)
               odclds = secdiff * taucmcd(iplon,igc,lev)
               absclds = 1.  - exp(-odclds)
               efclfracs = absclds * cldfmcd(iplon, igc,lev)
               radlu = radlu - radlu * (atrans(lev) + &
                   efclfracs * (1.  - atrans(lev))) + &
                   gassrc + cldfmcd(iplon, igc,lev) * &
                   (bbutot(lev) * atot(lev) - gassrc)
               gurad(iplon, igc, lev) = gurad(iplon, igc, lev) + radlu
               if (idrvd .eq. 1) then
                  d_radlu_dt = d_radlu_dt * cldfmcd(iplon,igc,lev) * (1.0  - atot(lev)) + &
                         d_radlu_dt * (1.0  - cldfmcd(iplon,igc,lev)) * (1.0  - atrans(lev))
                  
                  gdtotuflux_dtd( iplon, igc, lev)  = gdtotuflux_dtd(iplon, igc, lev) + d_radlu_dt
               endif
!  Clear layer
            else
               radlu = radlu + (bbugas(lev)-radlu)*atrans(lev)
               gurad(iplon, igc, lev) = gurad(iplon, igc, lev) + radlu
               if (idrvd .eq. 1) then
                  d_radlu_dt = d_radlu_dt * (1.0  - atrans(lev))
                  !d_urad_dt(lev) = d_urad_dt(lev) + d_radlu_dt
                  gdtotuflux_dtd( iplon, igc, lev) = gdtotuflux_dtd( iplon, igc, lev) + d_radlu_dt
               endif
            endif

!  Set clear sky stream to total sky stream as long as all layers
!  are clear (iclddn=0).  Streams must be calculated separately at 
!  all layers when a cloud is present (ICLDDN=1), because surface 
!  reflectance is different for each stream.
            if (iclddn.eq.1) then
               radclru = radclru + (bbugas(lev)-radclru)*atrans(lev) 
               gclrurad(iplon, igc, lev) = gclrurad(iplon, igc, lev) + radclru
            else
               radclru = radlu
               gclrurad(iplon, igc, lev) = gurad(iplon, igc, lev)
            endif
            if (idrvd .eq. 1) then
               if (iclddn.eq.1) then
                  d_radclru_dt = d_radclru_dt * (1.0  - atrans(lev))
                  gdtotuclfl_dtd(iplon, igc, lev) = gdtotuclfl_dtd(iplon, igc, lev) + d_radclru_dt
               else
                  d_radclru_dt = d_radlu_dt
                  gdtotuclfl_dtd(iplon, igc, lev) = gdtotuflux_dtd(iplon, igc, lev)
               endif
            endif

         enddo
          
          
         tblind = wtdiff * delwaved(ibnd) * fluxfacd
         ! (dmb 2012) Now that the g-points values were created, we modify them 
         ! so that later summation (integration) will be simpler.  
         do lev = 0, nlayers  
           gurad(iplon, igc, lev) = gurad(iplon, igc, lev) * tblind
           gdrad(iplon, igc, lev) = gdrad(iplon, igc, lev) * tblind
           gclrurad(iplon, igc, lev) = gclrurad(iplon, igc, lev) * tblind
           gclrdrad(iplon, igc, lev) = gclrdrad(iplon, igc, lev) * tblind
         end do

         if (idrvd .eq. 1) then
           ! (dmb 2012) we do the same for the optional derivatives.
           do lev = 0, nlayers
             gdtotuflux_dtd(iplon, igc, lev) = gdtotuflux_dtd(iplon, igc, lev) * tblind
             gdtotuclfl_dtd(iplon, igc, lev) = gdtotuclfl_dtd(iplon, igc, lev) * tblind
           end do

         endif

#ifdef _CUDA
      endif
#else
        end do
      end do
#endif

       end subroutine rtrnmcg

     ! (dmb 2012) This subroutine adds up the indivial g-point fluxes to arrive at a 
     ! final upward and downward flux value for each column and layer.  This subroutine 
     ! is parallelized across the column and layer dimensions.  As long as we parallelize 
     ! across two of the three dimesnions, we should usually have enough GPU saturation.
     _gpuker subroutine rtrnadd(ncol, nlay, ngpt, drvf, ngb)

        integer, intent(in), value :: ncol
        integer, intent(in), value :: nlay
        integer, intent(in), value :: ngpt
        integer, intent(in), value :: drvf
        integer, intent(in)        :: ngb(:) ! band index

        
        
        integer :: iplon, ilay, igp
        !real :: d(ngptlw)
        ! (dmb 2012) compute the column and layer indices from the grid and block 
        ! configurations. 

#ifdef _CUDA
        iplon = (blockidx%x-1) * blockdim%x + threadidx%x
        ilay = (blockidx%y-1) * blockdim%y + threadidx%y - 1
       
        ! (dmb 2012) make sure that the column and layer are within range
        if (ilay <= nlay .and. iplon <= ncol) then
#else

        do iplon = 1, ncol
        do ilay = 0, nlay

#endif

        do igp = 1, ngpt
              
               totufluxd(iplon, ilay)=totufluxd(iplon, ilay)+gurad(iplon, igp, ilay)
               totdfluxd(iplon, ilay)=totdfluxd(iplon, ilay)+gdrad(iplon, igp, ilay)
               totuclfld(iplon, ilay)=totuclfld(iplon, ilay)+gclrurad(iplon, igp, ilay)
               totdclfld(iplon, ilay)=totdclfld(iplon, ilay)+gclrdrad(iplon, igp, ilay)
               if (ilay .eq. nlay) then
                 if (ngb(igp) .eq.  9) olrb09d(iplon) = olrb09d(iplon) + gurad(iplon, igp, nlay)
                 if (ngb(igp) .eq. 10) olrb10d(iplon) = olrb10d(iplon) + gurad(iplon, igp, nlay)
                 if (ngb(igp) .eq. 11) olrb11d(iplon) = olrb11d(iplon) + gurad(iplon, igp, nlay)
               end if

        end do
        if (drvf .eq. 1) then

            do igp = 1, ngpt
                
                dtotuflux_dtd(iplon, ilay) = dtotuflux_dtd(iplon, ilay) + gdtotuflux_dtd( iplon, igp, ilay)
                dtotuclfl_dtd(iplon, ilay) = dtotuclfl_dtd(iplon, ilay) + gdtotuclfl_dtd( iplon, igp, ilay)
                if (ilay .eq. nlay) then
                  if (ngb(igp) .eq.  9) dolrb09_dtd(iplon) = dolrb09_dtd(iplon) + gdtotuflux_dtd(iplon, igp, nlay)
                  if (ngb(igp) .eq. 10) dolrb10_dtd(iplon) = dolrb10_dtd(iplon) + gdtotuflux_dtd(iplon, igp, nlay)
                  if (ngb(igp) .eq. 11) dolrb11_dtd(iplon) = dolrb11_dtd(iplon) + gdtotuflux_dtd(iplon, igp, nlay)
                end if

            end do

        end if

#ifdef _CUDA
        end if
#else
        end do
        end do
#endif

        end subroutine


        ! (dmb 2012) This kernel computes the heating rates separately.  It is parallelized across the 
        ! columnn and layer dimensions.
        _gpuker subroutine rtrnheatrates(ncol, nlay)

         integer, intent(in), value :: ncol
        integer, intent(in), value :: nlay
      
        
        real :: t2
        integer :: iplon, ilay
#ifdef _CUDA        
        iplon = (blockidx%x-1) * blockdim%x + threadidx%x
        ilay = (blockidx%y-1) * blockdim%y + threadidx%y - 1
    
             if (ilay<nlay .and. iplon<=ncol) then
#else
        do iplon = 1, ncol
        do ilay = 0, nlay - 1
#endif
       
               t2 = pzd(iplon, ilay ) - pzd(iplon, ilay + 1)
               htrd(iplon, ilay) = heatfacd * ((totufluxd(iplon, ilay) - totdfluxd(iplon, ilay)) &
                - (totufluxd(iplon, ilay+1) - totdfluxd(iplon, ilay+1)))/t2
               htrcd(iplon, ilay) = heatfacd * ((totuclfld(iplon, ilay) - totdclfld(iplon, ilay)) &
               - (totuclfld(iplon, ilay+1) - totdclfld(iplon, ilay+1)))/t2
#ifdef _CUDA
            end if
#else
        end do
        end do
#endif
       
        end subroutine

        ! (dmb 2012) Copy needed variables over to the GPU.  These arrays are pretty small so simple 
        ! stream 0 assignment operators suffice.
        subroutine copyGPUrtrnmcg(pz, pwvcm, idrv)
            
            real , intent(in) :: pz(:,0:)             ! level (interface) pressures (hPa, mb)
            integer , intent(in) :: idrv             ! flag for calculation of dF/dt from 
            real , intent(in) :: pwvcm(:)

   
            
            pzd(:,:) = pz(:, 0:ubound(pzd,2))
            pwvcmd(:) = pwvcm
            idrvd = idrv
            bpaded = bpade
            heatfacd = heatfac
            fluxfacd = fluxfac
         
       end subroutine

       ! (dmb 2012) Allocate the arrays for the rtrnmc routine on the GPU.  Some of these arrays are 
       ! quite large as they contain all 3 dimensions.  Luckily, for the gurad arrays, no copying of data
       ! from the CPU is needed because they are only stored on the GPU.
       subroutine allocateGPUrtrnmcg(ncol, nlay, ngptlw, drvf)

          integer , intent(in) :: ncol, nlay, ngptlw, drvf

          allocate( taucmcd(ncol, ngptlw, nlay))
          allocate( pzd(ncol, 0:nlay))
          allocate( pwvcmd(ncol))
          allocate( semissd(ncol, nbndlw))
          allocate( planklayd(ncol,nlay,nbndlw))
          allocate( planklevd(ncol, 0:nlay, nbndlw))
          allocate( plankbndd(ncol,nbndlw))
          allocate ( gurad(ncol,ngptlw,0:nlay))     ! upward longwave flux (w/m2)
          allocate ( gdrad(ncol,ngptlw,0:nlay))     ! downward longwave flux (w/m2)
          allocate ( gclrurad(ncol,ngptlw,0:nlay))     ! clear sky upward longwave flux (w/m2)
          allocate ( gclrdrad(ncol,ngptlw,0:nlay))     ! clear sky downward longwave flux (w/m2)
          ! (dmb 2012) Only allocate the optional derivative arrays if the flag is set
          if (drvf .eq. 1) then
            
            allocate( gdtotuflux_dtd( ncol, ngptlw, 0:nlay))
            allocate( gdtotuclfl_dtd( ncol, ngptlw, 0:nlay))

          endif
          
          
          allocate (totufluxd(ncol, 0:nlay))     ! upward longwave flux (w/m2)
          allocate (totdfluxd(ncol, 0:nlay))     ! downward longwave flux (w/m2)
          allocate (fnetd(ncol, 0:nlay))         ! net longwave flux (w/m2)
          allocate (htrd(ncol, 0:nlay))          ! longwave heating rate (k/day)
          allocate (totuclfld(ncol, 0:nlay))     ! clear sky upward longwave flux (w/m2)
          allocate (totdclfld(ncol, 0:nlay))     ! clear sky downward longwave flux (w/m2)
          allocate (fnetcd(ncol, 0:nlay))        ! clear sky net longwave flux (w/m2)
          allocate (htrcd(ncol, 0:nlay))         ! clear sky longwave heating rate (k/day)
          allocate (dtotuflux_dtd(ncol, 0:nlay)) ! change in upward longwave flux (w/m2/k)
          allocate (dtotuclfl_dtd(ncol, 0:nlay))
          allocate (dplankbnd_dtd(ncol,nbndlw)) 

          allocate (olrb09d(ncol))
          allocate (dolrb09_dtd(ncol))
          allocate (olrb10d(ncol))
          allocate (dolrb10_dtd(ncol))
          allocate (olrb11d(ncol))
          allocate (dolrb11_dtd(ncol))

       end subroutine 

       ! (dmb 2012) This subroutine deallocates rtrnmc related GPU arrays.
       subroutine deallocateGPUrtrnmcg( drvf )
        integer , intent(in) :: drvf
          

          deallocate( taucmcd)
          deallocate( pzd)
          deallocate( pwvcmd)
          deallocate( semissd)
          deallocate( planklayd)
          deallocate( planklevd)
          deallocate( plankbndd)
          deallocate ( gurad)     ! upward longwave flux (w/m2)
          deallocate ( gdrad)     ! downward longwave flux (w/m2)
          deallocate ( gclrurad)     ! clear sky upward longwave flux (w/m2)
          deallocate ( gclrdrad)     ! clear sky downward longwave flux (w/m2)
          deallocate (totufluxd)     ! upward longwave flux (w/m2)
          deallocate (totdfluxd)     ! downward longwave flux (w/m2)
          deallocate (fnetd)         ! net longwave flux (w/m2)
          deallocate (htrd)          ! longwave heating rate (k/day)
          deallocate (totuclfld)     ! clear sky upward longwave flux (w/m2)
          deallocate (totdclfld)     ! clear sky downward longwave flux (w/m2)
          deallocate (fnetcd)        ! clear sky net longwave flux (w/m2)
          deallocate (htrcd)         ! clear sky longwave heating rate (k/day)
          deallocate (dtotuflux_dtd) ! change in upward longwave flux (w/m2/k)
          deallocate (dtotuclfl_dtd)
          deallocate (dplankbnd_dtd) 
          if ( drvf .eq. 1) then
            deallocate( gdtotuflux_dtd, gdtotuclfl_dtd )
          end if

          deallocate (olrb09d, dolrb09_dtd)
          deallocate (olrb10d, dolrb10_dtd)
          deallocate (olrb11d, dolrb11_dtd)

       end subroutine 
      end module gpu_rrtmg_lw_rtrnmc

