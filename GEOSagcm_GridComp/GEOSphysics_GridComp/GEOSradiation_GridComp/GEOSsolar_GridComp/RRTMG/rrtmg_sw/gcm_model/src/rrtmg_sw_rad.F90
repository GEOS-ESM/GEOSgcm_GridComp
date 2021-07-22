!?pmn: usage of pi from modules only? (as per earth_sun but later: want to keep zero diff for now)
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
      use iso_fortran_env, only : error_unit

      implicit none

      public :: rrtmg_sw, earth_sun

   contains

      subroutine rrtmg_sw ( &
         rpart, ncol, nlay, &
         scon, adjes, coszen, isolvar, &
         play, plev, tlay, &
         h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr, &
         iceflgsw, liqflgsw, &
         cld, ciwp, clwp, rei, rel, &
         dyofyr, zm, alat, &
         iaer, tauaer, ssaaer, asmaer, &
         asdir, asdif, aldir, aldif, &
         normFlx, numCPUs, &
         swuflx, swdflx, swhr, swuflxc, swdflxc, swhrc, &
         nirr, nirf, parr, parf, uvrr, uvrf, &

         ! optional inputs
         bndscl, indsolvar, solcycfrac)

!pmn: put heating rates outside??

      use parrrsw,  only : nbndsw
#ifdef _CUDA
      use cudafor
#endif 

      
      ! ----- Inputs -----

      ! dimensions
      ! ----------
      integer, intent(in) :: rpart                   ! Number of columns in a partition
      integer, intent(in) :: ncol                    ! Number of horizontal columns     
      integer, intent(in) :: nlay                    ! Number of model layers

      ! orbit
      ! -----
      real, intent(in) :: scon                       ! Solar constant (W/m2)

      ! Notes: SCON is total solar irradiance averaged over the solar cycle. If scon == 0,
      ! an internal solar constant will be used. This depends on the value of isolvar. For
      ! isolvar=-1, scon=1368.22 Wm-2 (Kurucz), while for isolvar>=0, scon=1360.85 Wm-2
      ! (NRLSSI2). If scon > 0.0, this internal solar constant will be scaled to the value
      ! of SCON provided.

      real, intent(in) :: adjes                      ! Flux adjustment for Earth/Sun distance
                                                     !    i.e., ~ 1/dist(Earth-Sun in AU)^2
      real, intent(in) :: coszen (ncol)              ! Cosine of solar zenith angle

      ! solar variability
      ! -----------------
      integer, intent(in) :: isolvar                 ! Flag for solar variability method

      ! Notes: isolvar = -1 uses the Kurucz source function, while isolvar >= 0 use the NRLSSI2 
      ! solar model. First, the behavior for SCON = 0: isolvar in {-1,0,3} all have a fixed solar
      ! input at 1AU. For isolvar = -1 it is the Kurucz solar constant of 1368.22 Wm-2, while for
      ! isolvar = {0,3} it is the NRLSSI2 solar constant of 1360.85 Wm-2 (for the 100-50000 cm-1
      ! spectral range only, and based on the mean of solar cycle 13-24). The spectral breakdown
      ! is hardcoded into the model, including the "quiet sun, faculae brightening and sunspot
      ! darkening terms" for NRLSSI2, which are based on the means over solar cycles 13-24. But
      ! an optional scaling by band (via multiplier bndscl(:)) can be applied for isolvar {-1,3}.
      ! For SCON > 0, the behavior of these three cases (isolvar {-1,0,3}) is simple: the solar
      ! input is just scaled to the provided SCON, uniformly across the spectrum (and for each of
      ! the quiet sun, faculae and sunpot terms for NRLSSI2). This scaling is applied before the
      ! optional bndscl scaling.
      !    The isolvar = {1,2} cases require more detailed explanation. Both use the NRLSSI2 model,
      ! and but allow the mix of the quiet sun, faculae, and sunspot terms to vary with time. Since
      ! each of these terms has a different spectral response, this gives a solar input spectra that
      ! varies with time. This is achieved by setting the factors svar_{i,f,s}, which are multipliers
      ! to the average spectral response of the quiet sun, faculae, and sunspot terms (with integral
      ! values {I,F,S}int W/m2.)
      !    The isolvar = 1 case is designed to allow simulations over a generic solar cycle. This
      ! is an 11-year long average solar cycle derived from actual solar cycles 13-24. We denote
      ! this cycle "AvgCyc11". This cycle is implemented as follows: an 11-year cycle is provided
      ! for the faculae and sunpot indices Mg and SB via the {mg,sb}avgcyc arrays in the NRLSSI2
      ! module, which tabulate the index values per month for the 11 years. In conjunction with
      ! this, a linear relationship is provided between Mg and svar_f and between SB and svar_s.
      ! The main input for isolvar = 1 is solcycfrac in [0,1], which is the normalized input
      ! position in the 11-year average cycle. This solcycfrac is used to interpolate into the
      ! {mg,sb}avgcyc arrays to values Mg and SB, which are then converted to svar_f and svar_s
      ! via the linear relationship provided. These linear relationships are desined such that
      ! svar_f = 1 at Mg = <Mg>, the time average Mg index over AvgCyc11, and similarly for svar_s
      ! and SB, such that <svar_{s,f}> are both unity. So, if solcyclfr is uniformly cycled in
      ! [0,1] by the caller of rrtmg_sw_rad(), then each of the faculae and sunspot terms will
      ! also cycle, providing a time and spectrally varying solar cycle, but still with average
      ! flux contributions over the cycle of Fint and Sint (for SCON.eq.0). An extra facility is
      ! provided via the optional indsolar(2) array, which allows for a further multiplication
      ! of the svar_{f,s} terms after they are formed by the linearization above. This multiplier
      ! is time varying, being one at solar minimum, and the value of indsolvar(1,2) at solar
      ! maximum (for 1=Mg, 2=SB respecively), and to vary linearly with solcycfrac between those
      ! extrema (see NRLSSI2 module for further details).
      !    Still discussing isolvar = 1, for SCON.eq.0 we take the hint from not explicitly
      ! setting SCON to let an indsolvar.ne.1 choice cause a deviation from the internal solar 
      ! constant, because, while the svar_{f,s} average to unity over a cycle *without* a time-
      ! varying indsolvar multiplier, they do not do so with it. If, on the other hand, a value
      ! SCON > 0 is provided, we ASSUME that it is a value that we should honor as a MEAN over
      ! the AvgCyc11 cycle. We do this by adjusting a time-invariant quiet sun svar_i factor
      ! to honor the provided mean SCON, even in the presence of indsolvar.ne.1. Further
      ! details are found in the code.
      !    Finally, the insolvar = 2 case. This is a data-driven case provided for when values
      ! of TSI (SCON) and Mg and SB indices are available from data. In this case, the AvgCyc11
      ! cycle is not used explicitly, but the index-to-svar linearizations ARE used. The input
      ! indsolvar in this case provides the (1=Mg,2=SB) indices and produces concomitant mult-
      ! ipliers svar_{f,s} by these linear relationships. Again, for SCON.eq.0 we take the hint
      ! from not explicitly setting SCON to let the solar constant vary naturally as a result
      ! of the above data driven svar_{f,s} values (and by keeping svar_i = 1). BUT, for input
      ! SCON > 0, we take this as the TIME-VARYING, DATA-SUPPLIED TSI and honor it at every time
      ! (not just in a mean sense as for isolvar = 1). We do this by adjusting the quiet sun
      ! svar_i term (see the code).

      real, intent(in), optional :: indsolvar (2)    ! Facular and sunspot amplitude scale facs (isolvar=1),
                                                     !    or Mg and SB indices (isolvar=2)
      real, intent(in), optional :: bndscl (nbndsw)  ! Scale factors for each band
      real, intent(in), optional :: solcycfrac       ! Fraction of averaged 11-year solar cycle (0-1)
                                                     !    at current time (isolvar=1)
                                                     !    0. represents the first day of year 1
                                                     !    1. represents the last day of year 11

      ! profile
      ! -------
      real, intent(in) :: play   (ncol,nlay)         ! Layer pressures (hPa)
      real, intent(in) :: plev   (ncol,nlay+1)       ! Interface pressures (hPa)
      real, intent(in) :: tlay   (ncol,nlay)         ! Layer temperatures (K)

      ! gases
      ! -----
      real, intent(in) :: h2ovmr (ncol,nlay)         ! H2O volume mixing ratio
      real, intent(in) :: o3vmr  (ncol,nlay)         ! O3 volume mixing ratio
      real, intent(in) :: co2vmr (ncol,nlay)         ! CO2 volume mixing ratio
      real, intent(in) :: ch4vmr (ncol,nlay)         ! Methane volume mixing ratio
      real, intent(in) :: n2ovmr (ncol,nlay)         ! Nitrous oxide volume mixing ratio
      real, intent(in) :: o2vmr  (ncol,nlay)         ! Oxygen volume mixing ratio

      ! cloud optics flags
      ! ------------------
      integer, intent(in) :: iceflgsw                ! Flag for ice particle specification
      integer, intent(in) :: liqflgsw                ! Flag for liquid droplet specification

      ! clouds
      ! ------
      real, intent(in) :: cld    (ncol,nlay)         ! Cloud fraction
      real, intent(in) :: ciwp   (ncol,nlay)         ! In-cloud ice water path (g/m2)
      real, intent(in) :: clwp   (ncol,nlay)         ! In-cloud liquid water path (g/m2)
      real, intent(in) :: rei    (ncol,nlay)         ! Cloud ice effective radius (microns)
      real, intent(in) :: rel    (ncol,nlay)         ! Cloud water drop effective radius (microns)

      ! cloud overlap (exponential)
      ! ---------------------------
      integer, intent(in) :: dyofyr                  ! Day of the year
      real, intent(in) :: zm     (ncol,nlay)         ! Heights of level midpoints
      real, intent(in) :: alat   (ncol)              ! Latitude of column

      ! aerosols (optical props, non-delta-scaled)
      ! ------------------------------------------
      integer, intent(in) :: iaer                    ! aerosol flag (0=off, 10=on)
      real, intent(in) :: tauaer (ncol,nlay,nbndsw)  ! aer optical depth    (iaer=10 only)
      real, intent(in) :: ssaaer (ncol,nlay,nbndsw)  ! aer single scat albedo (iaer=10 only)
      real, intent(in) :: asmaer (ncol,nlay,nbndsw)  ! aer asymmetry param    (iaer=10 only)

      ! surface albedos
      ! ---------------
      real, intent(in) :: asdir  (ncol)              ! UV/vis  surface albedo: direct rad
      real, intent(in) :: asdif  (ncol)              ! UV/vis  surface albedo: diffuse rad
      real, intent(in) :: aldir  (ncol)              ! Near-IR surface albedo: direct rad
      real, intent(in) :: aldif  (ncol)              ! Near-IR surface albedo: diffuse rad

      ! etc
      ! ---
      integer, intent(in) :: normFlx                 ! Normalize fluxes?
                                                     !   0 = no normalization
                                                     !   1 = normalize (by scon*coszen)

      integer, intent(in) :: numCPUs                 ! Number of cores per node
                                              
      ! ----- Outputs -----

      real, intent(out) :: swuflx  (ncol,nlay+1)     ! Total sky SW up   flux (W/m2)
      real, intent(out) :: swdflx  (ncol,nlay+1)     ! Total sky SW down flux (W/m2)
      real, intent(out) :: swuflxc (ncol,nlay+1)     ! Clear sky SW up   flux (W/m2)
      real, intent(out) :: swdflxc (ncol,nlay+1)     ! Clear sky SW down flux (W/m2)
      real, intent(out) :: swhr    (ncol,nlay)       ! Total sky SW heating rate (K/d)
      real, intent(out) :: swhrc   (ncol,nlay)       ! Clear sky SW heating rate (K/d)

      ! Output added for Land/Surface process
      real, intent(out) :: nirr    (ncol)            ! Near-IR direct  down SW flux (W/m2)
      real, intent(out) :: nirf    (ncol)            ! Near-IR diffuse down SW flux (W/m2)
      real, intent(out) :: parr    (ncol)            ! Visible direct  down SW flux (W/m2)
      real, intent(out) :: parf    (ncol)            ! Visible diffuse down SW flux (W/m2)
      real, intent(out) :: uvrr    (ncol)            ! UV      direct  down SW flux (W/m2)
      real, intent(out) :: uvrf    (ncol)            ! UV      diffuse down SW flux (W/m2)

      ! ----- Locals -----

      integer :: pncol
      
#ifdef _CUDA
      type(cudadeviceprop) :: prop
      real :: gmem
      integer :: err
      integer :: munits
      integer :: numDevices, numCPUsPerGPU
      real :: maxmem
#endif
      
      ! ASSERTs to catch unphysical or invalid inputs
      ! ----------------------------------------------

      if (any(play   < 0.)) then
        write(error_unit,*) 'file:', __FILE__, ', line:', __LINE__
        write(error_unit,*) 'minval(play):', minval(play)
        error stop 'negative values in input: play'
      end if
      if (any(plev   < 0.)) then
        write(error_unit,*) 'file:', __FILE__, ', line:', __LINE__
        write(error_unit,*) 'minval(plev):', minval(plev)
        error stop 'negative values in input: plev'
      end if
      if (any(tlay   < 0.)) then
        write(error_unit,*) 'file:', __FILE__, ', line:', __LINE__
        write(error_unit,*) 'minval(tlay):', minval(tlay)
        error stop 'negative values in input: tlay'
      end if
      if (any(h2ovmr < 0.)) then
        write(error_unit,*) 'file:', __FILE__, ', line:', __LINE__
        write(error_unit,*) 'minval(h2ovmr):', minval(h2ovmr)
        error stop 'negative values in input: h2ovmr'
      end if
      if (any(o3vmr  < 0.)) then
        write(error_unit,*) 'file:', __FILE__, ', line:', __LINE__
        write(error_unit,*) 'minval(o3vmr):', minval(o3vmr)
        error stop 'negative values in input: o3vmr'
      end if
      if (any(co2vmr < 0.)) then
        write(error_unit,*) 'file:', __FILE__, ', line:', __LINE__
        write(error_unit,*) 'minval(co2vmr):', minval(co2vmr)
        error stop 'negative values in input: co2vmr'
      end if
      if (any(ch4vmr < 0.)) then
        write(error_unit,*) 'file:', __FILE__, ', line:', __LINE__
        write(error_unit,*) 'minval(ch4vmr):', minval(ch4vmr)
        error stop 'negative values in input: ch4vmr'
      end if
      if (any(n2ovmr < 0.)) then
        write(error_unit,*) 'file:', __FILE__, ', line:', __LINE__
        write(error_unit,*) 'minval(n2ovmr):', minval(n2ovmr)
        error stop 'negative values in input: n2ovmr'
      end if
      if (any(o2vmr  < 0.)) then
        write(error_unit,*) 'file:', __FILE__, ', line:', __LINE__
        write(error_unit,*) 'minval(o2vmr):', minval(o2vmr)
        error stop 'negative values in input: o2vmr'
      end if
      if (any(asdir  < 0.)) then
        write(error_unit,*) 'file:', __FILE__, ', line:', __LINE__
        write(error_unit,*) 'minval(asdir):', minval(asdir)
        error stop 'negative values in input: asdir'
      end if
      if (any(aldir  < 0.)) then
        write(error_unit,*) 'file:', __FILE__, ', line:', __LINE__
        write(error_unit,*) 'minval(aldir):', minval(aldir)
        error stop 'negative values in input: aldir'
      end if
      if (any(asdif  < 0.)) then
        write(error_unit,*) 'file:', __FILE__, ', line:', __LINE__
        write(error_unit,*) 'minval(asdif):', minval(asdif)
        error stop 'negative values in input: asdif'
      end if
      if (any(aldif  < 0.)) then
        write(error_unit,*) 'file:', __FILE__, ', line:', __LINE__
        write(error_unit,*) 'minval(aldif):', minval(aldif)
        error stop 'negative values in input: aldif'
      end if
      if (any(cld    < 0.)) then
        write(error_unit,*) 'file:', __FILE__, ', line:', __LINE__
        write(error_unit,*) 'minval(cld):', minval(cld)
        error stop 'negative values in input: cld'
      end if
      if (any(ciwp   < 0.)) then
        write(error_unit,*) 'file:', __FILE__, ', line:', __LINE__
        write(error_unit,*) 'minval(ciwp):', minval(ciwp)
        error stop 'negative values in input: ciwp'
      end if
      if (any(clwp   < 0.)) then
        write(error_unit,*) 'file:', __FILE__, ', line:', __LINE__
        write(error_unit,*) 'minval(clwp):', minval(clwp)
        error stop 'negative values in input: clwp'
      end if
      if (any(rei    < 0.)) then
        write(error_unit,*) 'file:', __FILE__, ', line:', __LINE__
        write(error_unit,*) 'minval(rei):', minval(rei)
        error stop 'negative values in input: rei'
      end if
      if (any(rel    < 0.)) then
        write(error_unit,*) 'file:', __FILE__, ', line:', __LINE__
        write(error_unit,*) 'minval(rel):', minval(rel)
        error stop 'negative values in input: rel'
      end if
      if (any(tauaer < 0.)) then
        write(error_unit,*) 'file:', __FILE__, ', line:', __LINE__
        write(error_unit,*) 'minval(tauaer):', minval(tauaer)
        error stop 'negative values in input: tauaer'
      end if
      if (any(ssaaer < 0.)) then
        write(error_unit,*) 'file:', __FILE__, ', line:', __LINE__
        write(error_unit,*) 'minval(ssaaer):', minval(ssaaer)
        error stop 'negative values in input: ssaaer'
      end if

      ! set column partition size pncol
      ! -------------------------------

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
      
      ! do a partition
      ! --------------

      call rrtmg_sw_sub ( &
         pncol, ncol, nlay, &
         scon, adjes, coszen, isolvar, &
         play, plev, tlay, &
         h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr, &
         iceflgsw, liqflgsw, &
         cld, ciwp, clwp, rei, rel, &
         dyofyr, zm, alat, &
         iaer, tauaer, ssaaer, asmaer, &
         asdir, asdif, aldir, aldif, &
         normFlx, swuflx, swdflx, swhr, swuflxc, swdflxc, swhrc, &
         nirr, nirf, parr, parf, uvrr, uvrf, &

         ! optional inputs
         bndscl, indsolvar, solcycfrac)
                                                      
      end subroutine rrtmg_sw                                                     


      subroutine rrtmg_sw_sub ( &
         pncol, gncol, nlay, &
         scon, adjes, gcoszen, isolvar, &
         gplay, gplev, gtlay, &
         gh2ovmr, go3vmr, gco2vmr, gch4vmr, gn2ovmr, go2vmr, &
         iceflgsw, liqflgsw, &
         gcld, gciwp, gclwp, grei, grel, &
         dyofyr, gzm, galat, &
         iaer, gtauaer, gssaaer, gasmaer, &
         gasdir, gasdif, galdir, galdif, &
         normFlx, swuflx, swdflx, swhr, swuflxc, swdflxc, swhrc, &
         nirr, nirf, parr, parf, uvrr, uvrf, &

         ! optional inputs
         bndscl, indsolvar, solcycfrac)


!?pmn: all needed?
      use parrrsw, only : nbndsw, ngptsw, nmol, mxmol, &
                          jpband, jpb1, jpb2, rrsw_scon
      use rrsw_con, only : heatfac, oneminus, pi, grav, avogad
      use rrsw_wvn, only : wavenum1, wavenum2
      use rrsw_cld, only : extliq1, ssaliq1, asyliq1, &
                           extice2, ssaice2, asyice2, &
                           extice3, ssaice3, asyice3, fdlice3, &
                           extice4, ssaice4, asyice4, &
                           abari, bbari, cbari, dbari, ebari, fbari
      use rrsw_wvn, only : ngb, icxa, nspa, nspb
      use rrsw_ref, only : preflog, tref
      use tab_xcw
      use NRLSSI2, only : initialize_NRLSSI2, &
                          adjust_solcyc_amplitudes, &
                          interpolate_indices, &
                          Iint, Fint, Sint, &
                          Mg_avg, Mg_0, &
                          SB_avg, SB_0, &
                          isolvar_1_mean_svar_f, &
                          isolvar_1_mean_svar_s

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

      ! ----- Inputs -----
      ! (see rrtmg_sw() for more detailed comments)

      ! dimensions
      integer, intent(in) :: pncol                     ! Nominal horiz cols in a partition
      integer, intent(in) :: gncol                     ! Global number of horizontal columns
      integer, intent(in) :: nlay                      ! Number of model layers

      ! orbit
      real, intent(in) :: scon                         ! Solar constant (W/m2)
      real, intent(in) :: adjes                        ! Flux adjustment for Earth/Sun distance
      real, intent(in) :: gcoszen (gncol)              ! Cosine of solar zenith angle

      ! solar variability
      integer, intent(in) :: isolvar
      real, intent(in), optional :: indsolvar (2)
      real, intent(in), optional :: bndscl (nbndsw)
      real, intent(in), optional :: solcycfrac

      ! profile
      real, intent(in) :: gplay   (gncol,nlay)         ! Layer pressures (hPa)
      real, intent(in) :: gplev   (gncol,nlay+1)       ! Interface pressures (hPa)
      real, intent(in) :: gtlay   (gncol,nlay)         ! Layer temperatures (K)

      ! gases
      real, intent(in) :: gh2ovmr (gncol,nlay)         ! H2O volume mixing ratio
      real, intent(in) :: go3vmr  (gncol,nlay)         ! O3 volume mixing ratio
      real, intent(in) :: gco2vmr (gncol,nlay)         ! CO2 volume mixing ratio
      real, intent(in) :: gch4vmr (gncol,nlay)         ! Methane volume mixing ratio
      real, intent(in) :: gn2ovmr (gncol,nlay)         ! Nitrous oxide volume mixing ratio
      real, intent(in) :: go2vmr  (gncol,nlay)         ! Oxygen volume mixing ratio

      ! cloud optics flags
      integer, intent(in) :: iceflgsw                  ! Flag for ice particle specifn
      integer, intent(in) :: liqflgsw                  ! Flag for liquid droplet specifn
      
      ! clouds
      real, intent(in) :: gcld    (gncol,nlay)         ! Cloud fraction
      real, intent(in) :: gciwp   (gncol,nlay)         ! In-cloud ice water path (g/m2)
      real, intent(in) :: gclwp   (gncol,nlay)         ! In-cloud liquid water path (g/m2)
      real, intent(in) :: grei    (gncol,nlay)         ! Cloud ice effective radius (um)
      real, intent(in) :: grel    (gncol,nlay)         ! Cloud drop effective radius (um)
                                                      
      ! cloud overlap
      integer, intent(in) :: dyofyr                    ! Day of the year
      real, intent(in) :: gzm     (gncol,nlay)         ! Heights of level midpoints
      real, intent(in) :: galat   (gncol)              ! Latitudes of columns
                                              
      ! aerosols (optical props, non-delta-scaled)
      integer, intent(in) :: iaer                      ! aerosol flag (0=off, 10=on)
      real, intent(in) :: gtauaer (gncol,nlay,nbndsw)  ! aer optical depth   (iaer=10 only)    
      real, intent(in) :: gssaaer (gncol,nlay,nbndsw)  ! aer single scat alb (iaer=10 only)    
      real, intent(in) :: gasmaer (gncol,nlay,nbndsw)  ! aer asymmetry param (iaer=10 only)    

      ! surface albedos
      real, intent(in) :: gasdir  (gncol)              ! UV/vis  surface albedo: direct rad
      real, intent(in) :: gasdif  (gncol)              ! UV/vis  surface albedo: diffuse rad
      real, intent(in) :: galdir  (gncol)              ! Near-IR surface albedo: direct rad
      real, intent(in) :: galdif  (gncol)              ! Near-IR surface albedo: diffuse rad

      integer, intent(in) :: normFlx                   ! Normalize fluxes flag

      ! ----- Output -----

      real, intent(out) :: swuflx  (gncol,nlay+1)      ! Total sky SW up   flux (W/m2)
      real, intent(out) :: swdflx  (gncol,nlay+1)      ! Total sky SW down flux (W/m2)
      real, intent(out) :: swuflxc (gncol,nlay+1)      ! Clear sky SW up   flux (W/m2)
      real, intent(out) :: swdflxc (gncol,nlay+1)      ! Clear sky SW down flux (W/m2)
      real, intent(out) :: swhr    (gncol,nlay)        ! Total sky SW heating rate (K/d)
      real, intent(out) :: swhrc   (gncol,nlay)        ! Clear sky SW heating rate (K/d)

      ! Output added for Land/Surface process
      real , intent(out) :: nirr   (gncol)             ! Near-IR direct  down SW flux (w/m2)
      real , intent(out) :: nirf   (gncol)             ! Near-IR diffuse down SW flux (w/m2)
      real , intent(out) :: parr   (gncol)             ! Visible direct  down SW flux (w/m2)
      real , intent(out) :: parf   (gncol)             ! Visible diffuse down SW flux (w/m2)
      real , intent(out) :: uvrr   (gncol)             ! UV      direct  down SW flux (w/m2)
      real , intent(out) :: uvrf   (gncol)             ! UV      diffuse down SW flux (w/m2)

      ! ----- Locals -----

!? pmn review all used
     
      ! Control
      integer :: istart                  ! beginning band of calculation
      integer :: iend                    ! ending band of calculation
      real :: zepsec, zepzen             ! epsilon
      real :: zdpgcp                     ! flux to heating conversion ratio

      integer :: ibnd, icol, ilay, ilev  ! various indices

      ! Atmosphere
      real :: coldry (pncol,nlay+1)      ! dry air column amount
      real :: wkl (pncol,mxmol,nlay)     ! molecular amounts (mol/cm-2)

      ! solar input
      real :: coszen (pncol)             ! Cosine of solar zenith angle
      real :: cossza (pncol)             ! Cosine of solar zenith angle
      real :: adjflux (jpband)           ! adjustment for curr Earth/Sun distance
      real :: swdflx_at_top (gncol)      ! swdflx at TOA

      ! surface albedos
      real :: albdir (pncol,nbndsw)      ! surface albedo, direct
      real :: albdif (pncol,nbndsw)      ! surface albedo, diffuse
      
      ! cloud overlap
      real :: rdl (pncol), adl (pncol)
      real :: CDF   (pncol,nlay,ngptsw)
      real :: CDF2  (pncol,nlay,ngptsw)
      real :: CDF3  (pncol,nlay,ngptsw)
      real :: alpha (pncol,nlay)

      ! Atmosphere - setcoef
      ! --------------------

      integer :: laytrop  (pncol)            ! tropopause layer index
      integer :: laylow   (pncol)            ! tropopause layer index
      integer :: jp  (pncol,nlay+1)          ! 
      integer :: jt  (pncol,nlay+1)          !
      integer :: jt1 (pncol,nlay+1)          !

      ! gasesous absorbers
      real :: colh2o  (pncol,nlay+1)         ! column amount (h2o)
      real :: colco2  (pncol,nlay+1)         ! column amount (co2)
      real :: colo3   (pncol,nlay+1)         ! column amount (o3)
      real :: coln2o  (pncol,nlay+1)         ! column amount (n2o)
      real :: colch4  (pncol,nlay+1)         ! column amount (ch4)
      real :: colo2   (pncol,nlay+1)         ! column amount (o2)
      real :: colmol  (pncol,nlay+1)         ! column amount
      real :: co2mult (pncol,nlay+1)         ! column amount 

      integer :: indself (pncol,nlay+1) 
      integer :: indfor  (pncol,nlay+1) 
      real :: selffac    (pncol,nlay+1) 
      real :: selffrac   (pncol,nlay+1) 
      real :: forfac     (pncol,nlay+1) 
      real :: forfrac    (pncol,nlay+1) 

      real, dimension (pncol,nlay+1) :: &
         fac00, fac01, fac10, fac11  
      
      ! general
      real :: play (pncol,nlay)               ! Layer pressures (hPa)
      real :: plev (pncol,nlay+1)             ! Interface pressures (hPa)
      real :: tlay (pncol,nlay)               ! Layer temperatures (K)

      ! Atmosphere/clouds - cldprop
      ! ---------------------------

      integer :: ncbands                      ! num of cloud spectral bands

      real :: cld  (pncol,nlay)               ! Cloud fraction
      real :: ciwp (pncol,nlay)               ! In-cloud ice water path (g/m2)
      real :: clwp (pncol,nlay)               ! In-cloud liq water path (g/m2)
      real :: rei  (pncol,nlay)               ! Cloud ice effective radius (um)
      real :: rel  (pncol,nlay)               ! Cloud drop effective radius (um)
      
      real :: alat (pncol)
      real :: zm (pncol,nlay)
                                                      
      real, dimension (pncol) :: &
         znirr, znirf, zparr, zparf, zuvrr, zuvrf
      
      real :: taucmc (pncol,nlay+1,ngptsw)    ! in-cloud optical depth [mcica]
      real :: taormc (pncol,nlay+1,ngptsw)    ! unscaled in-cloud optl depth [mcica]
      real :: ssacmc (pncol,nlay+1,ngptsw)    ! in-cloud single scat albedo [mcica]
      real :: asmcmc (pncol,nlay+1,ngptsw)    ! in-cloud asymmetry param [mcica]
      real :: fsfcmc (pncol,nlay+1,ngptsw)    ! in-cloud forward scat frac [mcica]
      
      real :: cldfmcl (pncol,nlay+1,ngptsw)   ! cloud fraction [mcica]
      real :: ciwpmcl (pncol,nlay+1,ngptsw)   ! in-cloud ice water path [mcica]
      real :: clwpmcl (pncol,nlay+1,ngptsw)   ! in-cloud liquid water path [mcica]

      ! Atmosphere/clouds/aerosol - spcvrt,spcvmc
      ! -----------------------------------------

!? pmn why nlay+1
      real :: ztauc (pncol,nlay+1,nbndsw)     ! cloud optical depth
      real :: ztaucorig (pncol,nlay+1,nbndsw) ! unscaled cloud optical depth
      real :: zasyc (pncol,nlay+1,nbndsw)     ! cloud asymmetry parameter 
      real :: zomgc (pncol,nlay+1,nbndsw)     ! cloud single scattering albedo
   
      real :: taua (pncol,nlay+1,nbndsw)
      real :: asya (pncol,nlay+1,nbndsw)
      real :: omga (pncol,nlay+1,nbndsw)

!? pmn why nlay+2
      real :: zbbfu    (pncol,nlay+2)         ! temporary up SW flux (w/m2)
      real :: zbbfd    (pncol,nlay+2)         ! temporary down SW flux (w/m2)
      real :: zbbcu    (pncol,nlay+2)         ! temporary clear sky up SW flux (w/m2)
      real :: zbbcd    (pncol,nlay+2)         ! temporary clear sky down SW flux (w/m2)
      real :: zbbfddir (pncol,nlay+2)         ! temporary down direct SW flux (w/m2)
      real :: zbbcddir (pncol,nlay+2)         ! temporary clear sky down direct SW flux (w/m2)
      real :: zuvfd    (pncol,nlay+2)         ! temporary UV down SW flux (w/m2)
      real :: zuvcd    (pncol,nlay+2)         ! temporary clear sky UV down SW flux (w/m2)
      real :: zuvfddir (pncol,nlay+2)         ! temporary UV down direct SW flux (w/m2)
      real :: zuvcddir (pncol,nlay+2)         ! temporary clear sky UV down direct SW flux (w/m2)
      real :: znifd    (pncol,nlay+2)         ! temporary near-IR down SW flux (w/m2)
      real :: znicd    (pncol,nlay+2)         ! temporary clear sky near-IR down SW flux (w/m2)
      real :: znifddir (pncol,nlay+2)         ! temporary near-IR down direct SW flux (w/m2)
      real :: znicddir (pncol,nlay+2)         ! temporary clear sky near-IR down direct SW flux (w/m2)

      ! Output fields 
      ! -------------

      real :: swnflx   (pncol,nlay+2)         ! Total sky SW net flux (W/m2)
      real :: swnflxc  (pncol,nlay+2)         ! Clear sky SW net flux (W/m2)
      real :: dirdflux (pncol,nlay+2)         ! Direct down SW surface flux
      real :: difdflux (pncol,nlay+2)         ! Diffuse down SW surface flux
      real :: uvdflx   (pncol,nlay+2)         ! Total sky down SW flux, UV/vis  
      real :: nidflx   (pncol,nlay+2)         ! Total sky down SW flux, near-IR 
      real :: dirdnuv  (pncol,nlay+2)         ! Direct down SW flux, UV/vis
      real :: difdnuv  (pncol,nlay+2)         ! Diffuse down SW flux, UV/vis
      real :: dirdnir  (pncol,nlay+2)         ! Direct down SW flux, near-IR
      real :: difdnir  (pncol,nlay+2)         ! Diffuse down SW flux, near-IR

      ! Solar variability multipliers
      ! -----------------------------

      real :: svar_f               ! facular multiplier
      real :: svar_s               ! sunspot multiplier
      real :: svar_i               ! baseline irradiance multiplier
      real :: svar_f_bnd (jpband)  ! facular multiplier (by band)
      real :: svar_s_bnd (jpband)  ! sunspot multiplier (by band)
      real :: svar_i_bnd (jpband)  ! baseline irradiance multiplier (by band)

!? pmn
      real gpu_device :: zgco  (pncol,ngptsw,nlay+1), zomco  (pncol,ngptsw,nlay+1)  
      real gpu_device :: zrdnd (pncol,ngptsw,nlay+1) 
      real gpu_device :: zref  (pncol,ngptsw,nlay+1), zrefo  (pncol,ngptsw,nlay+1)  
      real gpu_device :: zrefd (pncol,ngptsw,nlay+1), zrefdo (pncol,ngptsw,nlay+1)  
      real gpu_device :: ztauo (pncol,ngptsw,nlay)  
      real gpu_device :: zdbt  (pncol,ngptsw,nlay+1), ztdbt  (pncol,ngptsw,nlay+1)   
      real gpu_device :: ztra  (pncol,ngptsw,nlay+1), ztrao  (pncol,ngptsw,nlay+1)  
      real gpu_device :: ztrad (pncol,ngptsw,nlay+1), ztrado (pncol,ngptsw,nlay+1)  
      real gpu_device :: zfd   (pncol,ngptsw,nlay+1), zfu    (pncol,ngptsw,nlay+1)  
      real gpu_device :: zsflxzen(pncol,ngptsw)
      real gpu_device :: ssi   (pncol,ngptsw)
      real gpu_device :: ztaur (pncol,nlay,ngptsw), ztaug (pncol,nlay,ngptsw) 

      integer :: npart_clr, npart_cld, npart
      integer, dimension (gncol) :: &
         cldflag, gicol_clr, gicol_cld

      real, parameter :: amd = 28.9660     ! Effective molecular weight of dry air (g/mol)
      real, parameter :: amw = 18.0160     ! Molecular weight of water vapor (g/mol)

! Set molecular weight ratios (for converting mmr to vmr), e.g. h2ovmr = h2ommr * amdw

      real, parameter :: amdw  = 1.607793  ! Molecular weight of dry air / water vapor
      real, parameter :: amdc  = 0.658114  ! Molecular weight of dry air / carbon dioxide
      real, parameter :: amdo  = 0.603428  ! Molecular weight of dry air / ozone
      real, parameter :: amdm  = 1.805423  ! Molecular weight of dry air / methane
      real, parameter :: amdn  = 0.658090  ! Molecular weight of dry air / nitrous oxide
      real, parameter :: amdo2 = 0.905140  ! Molecular weight of dry air / oxygen

      real, parameter :: sbc = 5.67e-08    ! Stefan-Boltzmann constant (W/m2K4)

      integer :: n, imol, gicol            ! Loop indices
      real :: adjflx                       ! flux adjustment for Earth/Sun distance
      
      integer :: ipart, ncol_clr, ncol_cld, col_last, cols, cole, ncol, cc

      ! other solar variability locals
      ! ------------------------------
      real :: solvar (jpband)              ! solar constant scaling factor by band
      real :: indsolvar_scl (2)            ! Adjusted facular and sunspot amplitude 
                                           !   scale factors (isolvar=1)
      real :: indsolvar_ndx (2)            ! Facular and sunspot indices (isolvar=2)

      real :: solcycfr, Mg_now, SB_now
      real :: scon_int, svar_r

      ! Initializations
      ! ---------------

!? pmn
      zepsec = 1.e-06
      zepzen = 1.e-10
      oneminus = 1.0 - zepsec
      pi = 2. * asin(1.)

!? pmn
      istart = jpb1
      iend = jpb2

      ! solar variability: default values
      solvar(:) = 1.
      adjflux(:) = 1.
      svar_f = 1.
      svar_s = 1. 
      svar_i = 1. 
      svar_f_bnd(:) = 1. 
      svar_s_bnd(:) = 1. 
      svar_i_bnd(:) = 1. 

      ! isolvar == 1 specifies the position in AvgCyc11 through solcycfrac
      ! and allows scaling of solar cycle amplitudes as described in notes.
      ! ------------------------------------------------------------------

      if (isolvar .eq. 1) then 

         ! require solcycfrac present, else what's the point of using isolvar=1 ?
         if (.not.present(solcycfrac)) then
            write(error_unit,*) 'file:', __FILE__, ', line:', __LINE__
            error stop 'RRTMG_SW: isolvar == 1 requires solcycfrac present!'
         end if
         solcycfr = solcycfrac

         ! No amplitude scaling unless indsolvar is present. 
         indsolvar_scl(1:2) = 1.

         if (present(indsolvar)) then 

            ! Adjust amplitude scaling of mean solar cycle to be unity at
            ! solar minimum (solcycfrac_min), to be the requested indsolvar
            ! at solar maximum (solcycfrac_max), and to vary linearly with
            ! solcycfr between those values.

            if (indsolvar(1).ne.1. .or. indsolvar(2).ne.1.) &
               call adjust_solcyc_amplitudes(solcycfr, indsolvar, indsolvar_scl)

         endif

      endif

      ! isolvar == 2 allows direct specification of Mg and SB via indsolvar
      ! -------------------------------------------------------------------
      
      if (isolvar .eq. 2) then 

         ! default to mean indices
         indsolvar_ndx(1) = Mg_avg
         indsolvar_ndx(2) = SB_avg

         ! update to specified indices if provided
         if (present(indsolvar)) then 
            indsolvar_ndx(1) = indsolvar(1)
            indsolvar_ndx(2) = indsolvar(2)
         endif

      endif

      ! pre-calculated constants (will only do calcs once internally)
      ! -------------------------------------------------------------
      call initialize_NRLSSI2 (isolvar, indsolvar)

      ! Set flux adjustment for current Earth/Sun distance (two options)
      ! ----------------------------------------------------------------
      ! (Set adjflx to 1. to use constant Earth/Sun distance of 1 AU). 

      ! 1) Provided by GCM via ADJES (from MAPL sun factor DIST ~ 1/r^2)
      adjflx = adjes

      ! 2) Calc Earth/Sun dist adj from DYOFYR, the cumulative day of year
      ! (Turned off but DYOFYR used by MCICA exponential cloud overlap).

      ! if (dyofyr .gt. 0) then
      !    adjflx = earth_sun(dyofyr)
      ! endif

      ! --------------------------------------------------------
      ! Apply selected solar variability option based on ISOLVAR
      ! and input solar constant SCON.
      ! --------------------------------------------------------

      if (scon == 0.) then 

         ! For scon = 0, use internally defined solar constant, which is
         ! 1368.22 Wm-2 (for ISOLVAR=-1) and 1360.85 Wm-2 (For ISOLVAR=0,3;
         ! Options ISOLVAR=1,2 model sol cyc varations from 1360.85 Wm-2).

         if (isolvar .eq. -1) then

            ! Constant sun (Kurucz)
            ! Apply optional scaling by band if bndscl present.

            solvar(jpb1:jpb2) = 1.
            if (present(bndscl)) solvar(jpb1:jpb2) = bndscl(:)

         elseif (isolvar .eq. 0) then

            ! Constant sun (NRLSSI2 model)
            ! Quiet sun, facular, and sunspot terms averaged over AvgCyc11.

            svar_f = 1.
            svar_s = 1.
            svar_i = 1.

         elseif (isolvar .eq. 1) then

            ! Apply NRLSSI2 solar irradiance model at a specified solcycfr
            ! within AvgCyc11, with the additional amplitude scalings in 
            ! indsolvar_scl.

            ! interpolate mean solar cycle to solcycfr
            call interpolate_indices (solcycfr, Mg_now, SB_now)

            ! Apply linear index-to-flux-multiplier-svar relationship
            ! with the additional indsolvar_scl scaling.
            svar_f = indsolvar_scl(1) * (Mg_now - Mg_0) / (Mg_avg - Mg_0)
            svar_s = indsolvar_scl(2) * (SB_now - SB_0) / (SB_avg - SB_0)
            svar_i = 1.

         elseif (isolvar .eq. 2) then

            ! Specified solar cycle with solar variability based on NRLSSI2 model.
            ! Facular and sunspot index terms input directly.

            svar_f = (indsolvar_ndx(1) - Mg_0) / (Mg_avg - Mg_0)
            svar_s = (indsolvar_ndx(2) - SB_0) / (SB_avg - SB_0)
            svar_i = 1.

         elseif (isolvar .eq. 3) then

            ! Constant sun (NRLSSI2 model).
            ! Averaged facular, sunspot and quiet sun terms from AvgCyc11.
            ! Apply optional scaling by band if bndscl present.

            solvar(jpb1:jpb2) = 1.
            if (present(bndscl)) solvar(jpb1:jpb2) = bndscl(:)
            do ibnd = jpb1,jpb2
               svar_f_bnd(ibnd) = solvar(ibnd)
               svar_s_bnd(ibnd) = solvar(ibnd)
               svar_i_bnd(ibnd) = solvar(ibnd)
            enddo

         else
            write(error_unit,*) 'file:', __FILE__, ', line:', __LINE__
            write(error_unit,*) 'bad isolvar value:', isolvar
            error stop 'RRTMG_SW: invalid isolvar'
         endif 

      elseif (scon > 0.) then 

         ! Scale from internal to externally specified SCON.

         if (isolvar .eq. -1) then

            ! Constant sun (Kurucz)
            ! Scale from internal to requested solar constant.
            ! Apply optional scaling by band if bndscl present.

            solvar(jpb1:jpb2) = scon / rrsw_scon 
            if (present(bndscl)) &
               solvar(jpb1:jpb2) = solvar(jpb1:jpb2) * bndscl(:)

         elseif (isolvar .eq. 0) then

            ! Constant sun (NRLSSI2 model)
            ! Quiet sun, facular, and sunspot terms averaged over AvgCyc11.
            ! Scale from internal to requested solar constant. 

            scon_int = Fint + Sint + Iint
            svar_r = scon / scon_int
            svar_f = svar_r
            svar_s = svar_r
            svar_i = svar_r

         elseif (isolvar .eq. 1) then

            ! Apply NRLSSI2 solar irradiance model at a specified solcycfr
            ! within AvgCyc11, with the additional amplitude scalings in
            ! indsolvar_scl. Scale from the internal to the requested solar
            ! constant, which is treated as a required *cycle average*.

            ! interpolate mean solar cycle to solcycfr
            call interpolate_indices (solcycfr, Mg_now, SB_now)

            ! Apply linear index-to-flux-multiplier-svar relationship
            ! with the additional indsolvar_scl scaling. Select a constant
            ! svar_i such that chosen scon is the <cycle average>.
            ! scon = svar_i * Iint + <svar_f> * Fint + <svar_s> * Sint >
            ! => svar_i = [scon - (<svar_f> * Fint + <svar_s> * Sint)] / Iint

            svar_f = indsolvar_scl(1) * (Mg_now - Mg_0) / (Mg_avg - Mg_0)
            svar_s = indsolvar_scl(2) * (SB_now - SB_0) / (SB_avg - SB_0)
            svar_i = (scon - (isolvar_1_mean_svar_f * Fint + &
                              isolvar_1_mean_svar_s * Sint)) / Iint

         elseif (isolvar .eq. 2) then

            ! Specified solar cycle with solar variability based on NRLSSI2 model.
            ! Facular and sunspot index terms input directly. Scale from internal
            ! to requested solar constant by setting svar_i so that
            !   svar_i * Iint + svar_f * Fint + svar_s * Sint = scon.
            ! So, scon is honored at EACH time, because it too is assumed to be
            ! specified from time-varying data.

            svar_f = (indsolvar_ndx(1) - Mg_0) / (Mg_avg - Mg_0)
            svar_s = (indsolvar_ndx(2) - SB_0) / (SB_avg - SB_0)
            svar_i = (scon - (svar_f * Fint + svar_s * Sint)) / Iint 

         elseif (isolvar .eq. 3) then

            ! Constant sun (NRLSSI2 model).
            ! Averaged facular, sunspot and quiet sun terms from AvgCyc11.
            ! Scale from internal to requested solar constant.
            ! Apply optional scaling by band if bndscl present.

            scon_int = Fint + Sint + Iint
            solvar(jpb1:jpb2) = scon / scon_int
            if (present(bndscl)) solvar(jpb1:jpb2) = solvar(jpb1:jpb2) * bndscl(:)
            do ibnd = jpb1,jpb2
               svar_f_bnd(ibnd) = solvar(ibnd)
               svar_s_bnd(ibnd) = solvar(ibnd)
               svar_i_bnd(ibnd) = solvar(ibnd)
            enddo

         else
            write(error_unit,*) 'file:', __FILE__, ', line:', __LINE__
            write(error_unit,*) 'bad isolvar value:', isolvar
            error stop 'RRTMG_SW: invalid isolvar'
         endif 

      else
         write(error_unit,*) 'file:', __FILE__, ', line:', __LINE__
         error stop 'RRTMG_SW: scon cannot be negative!'
      endif

      ! Earth-Sun distance adjustment
      adjflux(jpb1:jpb2) = adjflx

      ! Combine with solar constant scaling for Kurucz
      ! (done separately via svar_ for NRLSSI2)
      if (isolvar < 0) then
         adjflux(jpb1:jpb2) = adjflux(jpb1:jpb2) * solvar(jpb1:jpb2)
      endif
      
      ! determine cloud profile
      cldflag = 0
      do gicol = 1,gncol
         if (any(gcld(gicol,:) > 0)) cldflag(gicol) = 1
      end do

      ! build profile separation (clear/cloudy)
      ncol_clr = 0
      ncol_cld = 0

      do gicol = 1,gncol
         if (cldflag(gicol)==1) then
            ncol_cld = ncol_cld + 1
            gicol_cld(ncol_cld) = gicol
         else
            ncol_clr = ncol_clr + 1
            gicol_clr(ncol_clr) = gicol
         end if
      end do

!? pmn do as for lw eventually
      call TABULATE_XCW_BETA

!$acc data copyout(swuflxc, swdflxc, swuflx, swdflx, swnflxc, swnflx, swhrc, swhr) &
!$acc create(laytrop, laylow, jp, jt, jt1, &
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
!$acc create(ciwp, clwp, cld, rei, rel, rdl, adl) &
!$acc create(play, tlay, plev, cldflag, coszen, swdflx_at_top) &
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
!$acc copyin(ka28,kb28,sfluxref28)&
!$acc copyin(irradnceo28,facbrghto28,snsptdrko28)&
!$acc copyin(kao29,kbo29,selfrefo29,forrefo29,sfluxrefo29,absh2oo29,absco2o29)&
!$acc copyin(ka29,kb29,selfref29,forref29,sfluxref29,absh2o29,absco229)&
!$acc copyin(irradnceo29,facbrghto29,snsptdrko29)&
!$acc copyin(gh2ovmr, gco2vmr, go3vmr, gn2ovmr, gch4vmr, go2vmr)&
!$acc copyin(gcld, gciwp, gclwp, grei, grel, gplay, gplev, gtlay)&
!$acc copyin(gasdir, galdir, gasdif, galdif,gicol_cld,gicol_clr,gcoszen)&
!$acc copyout(nirr,nirf,parr,parf,uvrr,uvrf)

!$acc data copyin(XCW)

!$acc update device(extliq1, ssaliq1, asyliq1, extice2, ssaice2, asyice2) &
!$acc device(extice3, ssaice3, asyice3, fdlice3, abari, bbari, cbari, dbari, ebari, fbari) &
!$acc device(preflog)

      ! number of pncol partitions needed for each of clear and cloudy profiles
      npart_clr = ceiling( real(ncol_clr) / real(pncol) )
      npart_cld = ceiling( real(ncol_cld) / real(pncol) )

      ! zero McICA cloud physical props
      !$acc kernels    
      cldfmcl = 0.
      ciwpmcl = 0.
      clwpmcl = 0.     
      !$acc end kernels
  
      ! zero aerosols
      !$acc kernels
      taua = 0.
      asya = 0.
      omga = 1.
      !$acc end kernels

      ! aerosols requested
      if (iaer==10) then
         !$acc update device(gtauaer,gssaaer,gasmaer)
      end if

      ! partitioning over clear (cc=1) and cloudy (cc=2) columns
      ! --------------------------------------------------------

      do cc = 1,2  ! outer loop over clear then cloudy columns

         if (cc==1) then 
         
            ! clear
            npart = npart_clr
            col_last = ncol_clr

         else
        
            ! cloudy
            npart = npart_cld
            col_last = ncol_cld
         
         end if

         ! loop over partitions
         do ipart = 0,npart-1

            ! partition dimensions
            cols = ipart * pncol + 1
            cole = (ipart + 1) * pncol
            if (cole > col_last) cole = col_last
            ncol = cole - cols + 1

!?pmn defaults for clear cc==1 I thinjk ... add comment
            ! zero McICA cloud optical props
            !$acc kernels            
            taormc = 0.
            taucmc = 0.
            ssacmc = 1.
            asmcmc = 0.
!?pmn needed --- dont think so
            fsfcmc = 0.
            !$acc end kernels            

            ! copy inputs into partition
            ! --------------------------

            if (cc==1) then    

               ! -------------
               ! Clear columns
               ! -------------

               !$acc kernels loop private(gicol)
               do icol = 1,ncol
                  gicol = gicol_clr(icol + cols - 1)

                  ! assign surface albedos to bands

                  ! near IR bands 14=nbndsw and 1-8
                  ! 820-12850 cm-1, 0.778-12.2 um
                  do ibnd=1,8
                     albdir(icol,ibnd) = galdir(gicol)
                     albdif(icol,ibnd) = galdif(gicol)
                  enddo
                  albdir(icol,nbndsw) = galdir(gicol)
                  albdif(icol,nbndsw) = galdif(gicol)

                  ! UV/Vis bands 10-13
                  ! 16000-50000 cm-1, 0.200-0.625 um
                  do ibnd=10,13
                     albdir(icol,ibnd) = gasdir(gicol)
                     albdif(icol,ibnd) = gasdif(gicol)
                  enddo

                  ! Transition band 9
                  ! 12850-16000 cm-1, 0.625-0.778 um
                  ! Take average, dmlee
                  albdir(icol,9) = (gasdir(gicol)+galdir(gicol))/2.
                  albdif(icol,9) = (gasdif(gicol)+galdif(gicol))/2.

               enddo
               !$acc end kernels      

               ! copy in partition (general)
               !$acc kernels 
               do icol = 1,ncol
                  gicol = gicol_clr(icol + cols - 1)
    
                  play(icol,:) = gplay(gicol,1:nlay)
                  plev(icol,:) = gplev(gicol,1:nlay+1)
                  tlay(icol,:) = gtlay(gicol,1:nlay)

               enddo
               !$acc end kernels

               ! copy in partition (aerosols)
               if (iaer==10) then
                  !$acc kernels
                  do icol = 1,ncol
                     gicol = gicol_clr(icol + cols - 1)
                     taua(icol,1:nlay,:) = gtauaer(gicol,1:nlay,:)
                     asya(icol,1:nlay,:) = gasmaer(gicol,1:nlay,:)
                     omga(icol,1:nlay,:) = gssaaer(gicol,1:nlay,:)
!?pmn this ordering is very inefficient
                  enddo
                  !$acc end kernels
               endif   

               ! copy in partition (gases)
               !$acc kernels
               do icol = 1,ncol
                  gicol = gicol_clr(icol + cols - 1)
                  wkl(icol,1,:) = gh2ovmr(gicol,1:nlay)
                  wkl(icol,2,:) = gco2vmr(gicol,1:nlay)
                  wkl(icol,3,:) = go3vmr (gicol,1:nlay)
                  wkl(icol,4,:) = gn2ovmr(gicol,1:nlay)
                  wkl(icol,5,:) = 0.
                  wkl(icol,6,:) = gch4vmr(gicol,1:nlay)
                  wkl(icol,7,:) = go2vmr (gicol,1:nlay)   
                  coszen(icol)  = gcoszen(gicol)
                end do
                !$acc end kernels

            else

               ! --------------
               ! Cloudy columns
               ! --------------
          
               !$acc kernels loop private(gicol)
               do icol = 1,ncol
                  gicol = gicol_cld(icol + cols - 1)
     
                  ! assign surface albedos to bands

                  ! near IR bands 14=nbndsw and 1-8
                  ! 820-12850 cm-1, 0.778-12.2 um
                  do ibnd=1,8
                     albdir(icol,ibnd) = galdir(gicol)
                     albdif(icol,ibnd) = galdif(gicol)
                  enddo
                  albdir(icol,nbndsw) = galdir(gicol)
                  albdif(icol,nbndsw) = galdif(gicol)

                  ! UV/Vis bands 10-13
                  ! 16000-50000 cm-1, 0.200-0.625 um
                  do ibnd=10,13
                     albdir(icol,ibnd) = gasdir(gicol)
                     albdif(icol,ibnd) = gasdif(gicol)
                  enddo

                  ! Transition band 9
                  ! 12850-16000 cm-1, 0.625-0.778 um
                  ! Take average, dmlee
                  albdir(icol,9) = (gasdir(gicol)+galdir(gicol))/2.
                  albdif(icol,9) = (gasdif(gicol)+galdif(gicol))/2.

               enddo
               !$acc end kernels               
          
               ! copy in partition (general and cloud physical props)
               !$acc kernels 
               do icol = 1,ncol
                  gicol = gicol_cld(icol + cols - 1)
     
                  play(icol,:) = gplay(gicol,1:nlay)
                  plev(icol,:) = gplev(gicol,1:nlay+1)
                  tlay(icol,:) = gtlay(gicol,1:nlay)
                  cld (icol,:) = gcld (gicol,1:nlay)
                  ciwp(icol,:) = gciwp(gicol,1:nlay)
                  clwp(icol,:) = gclwp(gicol,1:nlay)
                  rei (icol,:) = grei (gicol,1:nlay) 
                  rel (icol,:) = grel (gicol,1:nlay)
                  zm  (icol,:) = gzm  (gicol,1:nlay)
                  alat(icol)   = galat(gicol)
               enddo
               !$acc end kernels

               ! copy in partition (aerosols)
               if (iaer==10) then
                  !$acc kernels    
                  do icol = 1,ncol
                     gicol = gicol_cld(icol + cols - 1)
                     taua(icol,1:nlay,:) = gtauaer(gicol,1:nlay,:)
                     asya(icol,1:nlay,:) = gasmaer(gicol,1:nlay,:)
                     omga(icol,1:nlay,:) = gssaaer(gicol,1:nlay,:)
                  end do
                  !$acc end kernels
               endif

               ! copy in partition (gases)
               !$acc kernels
               do icol = 1,ncol
                  gicol = gicol_cld(icol + cols - 1)
                  wkl(icol,1,:) = gh2ovmr(gicol,1:nlay)
                  wkl(icol,2,:) = gco2vmr(gicol,1:nlay)
                  wkl(icol,3,:) = go3vmr(gicol,1:nlay)
                  wkl(icol,4,:) = gn2ovmr(gicol,1:nlay)
                  wkl(icol,5,:) = 0.
                  wkl(icol,6,:) = gch4vmr(gicol,1:nlay)
                  wkl(icol,7,:) = go2vmr(gicol,1:nlay)  
                  coszen(icol)  = gcoszen(gicol)
               enddo
               !$acc end kernels

            end if  ! clear or cloudy columns

            ! limit tiny cosine zenith angles
            !$acc kernels
            do icol = 1,ncol
               cossza(icol) = max(zepzen,coszen(icol))
            enddo
            !$acc end kernels  

            ! evaluate dry air molecules/cm^2
            ! (see details in rrtmg_lw_rad())
            !$acc kernels
            do icol = 1,ncol
               do ilay = 1,nlay
                  coldry(icol,ilay) = (plev(icol,ilay)-plev(icol,ilay+1)) * 1.e3 * avogad / &
                     (1.e2 * grav * ((1.-wkl(icol,1,ilay)) * amd + wkl(icol,1,ilay) * amw) * &
                     (1. + wkl(icol,1,ilay)))
               enddo
            enddo
            !$acc end kernels

            ! gases also to molecules/cm^2
            !$acc kernels
            do icol = 1,ncol
               do ilay = 1,nlay
                  do imol = 1,nmol
                     wkl(icol,imol,ilay) = coldry(icol,ilay) * wkl(icol,imol,ilay)
                  end do
               end do
            end do
            !$acc end kernels

!issue that pncol is full partition(=dimension) while ncol is actual used
! so either input pncol as well to dimension or else only send in actual 
! needed like in LW ... study
!pmn needed working on and the replacing by abstract as per LW
            ! McICA subcolumn generation
            if (cc==2) then
               call mcica_sw( &
                  ncol, nlay, ngptsw, play, &
                  cld, clwp, ciwp, &
                  cldfmcl, clwpmcl, ciwpmcl, &
                  CDF, CDF2, CDF3, alpha, zm, &
                  alat, dyofyr, rdl, adl)
            end if   

            ! cloud optical property generation
            if (cc==2) then
               call cldprmc_sw( &
                  ncol, nlay, iceflgsw, liqflgsw,  &
                  cldfmcl, ciwpmcl, clwpmcl, rei, rel, &
                  taormc, taucmc, ssacmc, asmcmc)
            end if

            ! Calculate information needed by the radiative transfer routine
            ! that is specific to this atmosphere, especially some of the
            ! coefficients and indices needed to compute the optical depths
            ! by interpolating data from stored reference atmospheres.

            call setcoef_sw( &
               ncol, nlay, play, tlay, coldry, wkl, &
               laytrop, laylow, jp, jt, jt1, &
               co2mult, colch4, colco2, colh2o, colmol, coln2o, &
               colo2, colo3, fac00, fac01, fac10, fac11, &
               selffac, selffrac, indself, forfac, forfrac, indfor)

            ! compute sw radiative fluxes
            call spcvmc_sw( &
               cc, pncol, ncol, nlay, istart, iend, &
               play, tlay, albdif, albdir, &
               cldfmcl, taucmc, asmcmc, ssacmc, taormc, &
               taua, asya, omga,cossza, coldry, adjflux, &
               isolvar, svar_f, svar_s, svar_i, &
               svar_f_bnd, svar_s_bnd, svar_i_bnd, &
               laytrop, laylow, jp, jt, jt1, &
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

            ! Copy out up and down, clear and total sky fluxes to output arrays.
            ! Vertical indexing goes from bottom to top; reverse here for GCM if necessary.

            if (cc==1) then  ! clear columns

               !$acc kernels loop independent
               do icol = 1,ncol
                  gicol = gicol_clr(icol + cols - 1)
        
                  ! up and down fluxes
                  do ilev = 1,nlay+1
                     swuflxc(gicol,ilev) = zbbcu(icol,ilev) 
                     swdflxc(gicol,ilev) = zbbcd(icol,ilev) 
                     swuflx (gicol,ilev) = zbbfu(icol,ilev) 
                     swdflx (gicol,ilev) = zbbfd(icol,ilev) 
                  enddo

                  ! net fluxes
                  do ilev = 1,nlay+1
                     swnflxc(icol,ilev)  = swdflxc(gicol,ilev) - swuflxc(gicol,ilev)
                     swnflx (icol,ilev)  = swdflx (gicol,ilev) - swuflx (gicol,ilev)
                  enddo

                  ! heating rates
                  do ilay = 1,nlay
                     zdpgcp = heatfac / (plev(icol,ilay) - plev(icol,ilay+1))
                     swhrc(gicol,ilay) = (swnflxc(icol,ilay+1) - swnflxc(icol,ilay) ) * zdpgcp
                     swhr (gicol,ilay) = (swnflx (icol,ilay+1) - swnflx (icol,ilay) ) * zdpgcp
                  enddo
                  swhrc(gicol,nlay) = 0. 
                  swhr (gicol,nlay) = 0. 

               enddo
               !$acc end kernels 

               ! surface broadband fluxes
               !$acc kernels loop independent
               do icol = 1,ncol
                  gicol = gicol_clr(icol + cols - 1)
                  nirr(gicol) = znirr(icol)
                  nirf(gicol) = znirf(icol) - znirr(icol)
                  parr(gicol) = zparr(icol)
                  parf(gicol) = zparf(icol) - zparr(icol)
                  uvrr(gicol) = zuvrr(icol)
                  uvrf(gicol) = zuvrf(icol) - zuvrr(icol)
               end do
               !$acc end kernels 

            else ! cloudy columns

               !$acc kernels loop independent
               do icol = 1,ncol
                  gicol = gicol_cld(icol + cols - 1)
                  do ilev = 1,nlay+1
                     swuflxc(gicol,ilev) = zbbcu(icol,ilev) 
                     swdflxc(gicol,ilev) = zbbcd(icol,ilev) 
                     swuflx (gicol,ilev) = zbbfu(icol,ilev) 
                     swdflx (gicol,ilev) = zbbfd(icol,ilev) 
                  enddo

                  do ilev = 1,nlay+1
                     swnflxc(icol,ilev)  = swdflxc(gicol,ilev) - swuflxc(gicol,ilev)
                     swnflx (icol,ilev)  = swdflx (gicol,ilev) - swuflx (gicol,ilev)
                  enddo

                  do ilay = 1,nlay
                     zdpgcp = heatfac / (plev(icol,ilay) - plev(icol,ilay+1))
                     swhrc(gicol,ilay) = (swnflxc(icol,ilay+1) - swnflxc(icol,ilay)) * zdpgcp
                     swhr (gicol,ilay) = (swnflx (icol,ilay+1) - swnflx (icol,ilay)) * zdpgcp
                  enddo
                  swhrc(gicol,nlay) = 0. 
                  swhr (gicol,nlay) = 0. 

               enddo
               !$acc end kernels 

               !$acc kernels loop independent
               do icol = 1,ncol
                  gicol = gicol_cld(icol + cols - 1)
                  nirr(gicol) = znirr(icol)
                  nirf(gicol) = znirf(icol) - znirr(icol)
                  parr(gicol) = zparr(icol)
                  parf(gicol) = zparf(icol) - zparr(icol)
                  uvrr(gicol) = zuvrr(icol)
                  uvrf(gicol) = zuvrf(icol) - zuvrr(icol)
               enddo
               !$acc end kernels 

            endif  ! clear/cloudy

         enddo  ! over partitions

      enddo  ! outer loop (cc) over clear then cloudy columns

      ! If the user requests 'normalized' fluxes, divide
      ! the fluxes by the solar constant times coszen
      ! MAT This requires only lit points passed in

      if (normFlx == 1) then

         !$acc kernels
         swdflx_at_top(:) = max(swdflx(:,nlay+1),1e-7)

         do ilev = 1,nlay+1
            swuflxc(:,ilev) = swuflxc(:,ilev) / swdflx_at_top(:)
            swdflxc(:,ilev) = swdflxc(:,ilev) / swdflx_at_top(:)
            swuflx (:,ilev) = swuflx (:,ilev) / swdflx_at_top(:)
            swdflx (:,ilev) = swdflx (:,ilev) / swdflx_at_top(:)
         enddo

         nirr(:) = nirr(:) / swdflx_at_top(:)
         nirf(:) = nirf(:) / swdflx_at_top(:)
         parr(:) = parr(:) / swdflx_at_top(:)
         parf(:) = parf(:) / swdflx_at_top(:)
         uvrr(:) = uvrr(:) / swdflx_at_top(:)
         uvrf(:) = uvrf(:) / swdflx_at_top(:)
         !$acc end kernels

      endif

      !$acc end data ! for XCW

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

      integer, intent(in) :: idn

      real :: gamma

      gamma = 2. * pi * (idn-1)/365. 

      ! Use Iqbal's equation 1.2.1

      earth_sun = 1.000110  + .034221  * cos(gamma) + .001289  * sin(gamma) + &
                   .000719  * cos(2. *gamma) + .000077  * sin(2. *gamma)

      end function earth_sun

   end module rrtmg_sw_rad


