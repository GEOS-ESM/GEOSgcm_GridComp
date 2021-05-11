module rrtmg_lw_rad

   use cloud_condensate_inhomogeneity, only: &
      initialize_inhomogeneity, release_inhomogeneity
   use rrtmg_lw_cldprmc, only : cldprmc
   use rrtmg_lw_setcoef, only : setcoef, setcoef_free
   use rrtmg_lw_taumol, only : taumol
   use rrtmg_lw_rtrnmc, only : rtrnmc

   use iso_fortran_env, only : error_unit  ! for debugging

   implicit none

contains

   ! -----------------------------------------------------------------------------
   subroutine rrtmg_lw( &
      ncol, nlay, psize, idrv, &
      play, plev, tlay, tlev, tsfc, emis, & 
      h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr, &
      cfc11vmr, cfc12vmr, cfc22vmr, ccl4vmr, &
      cldf, ciwp, clwp, rei, rel, iceflglw, liqflglw, &
      tauaer, zm, alat, dyofyr, cloudLM, cloudMH, clearCounts, &
      uflx, dflx, uflxc, dflxc, duflx_dt, duflxc_dt, &
      olrb06, dolrb06_dt, olrb09, dolrb09_dt, &
      olrb10, dolrb10_dt, olrb11, dolrb11_dt)
   ! -----------------------------------------------------------------------------
   !
   ! This program is the driver subroutine for RRTMG_LW, the AER LW radiation 
   ! model for application to GCMs, that has been adapted from RRTM_LW for
   ! improved efficiency.
   !
   ! NOTE: all layering in RRTMG is ordered from surface to toa. 
   !
   ! NOTE: The call to RRTMG_LW_INI should be moved to the GCM initialization
   !  area, since this has to be called only once. 
   !
   ! This routine:
   !    a) calls GENERATE_STOCHASTIC_CLOUDS to generate subcolumn binary
   !       clouds for McICA (Monte Carlo Independent Column Approximation,
   !       Pincus et al., JC, 2003) treatment of cloudy radiative transfer;
   !    b) calls CLDPRMC to set cloud optical depth for McICA based 
   !       on input cloud physical properties;
   !    c) calls SETCOEF to calculate various quantities needed for 
   !       the radiative transfer algorithm;
   !    d) calls TAUMOL to calculate gaseous optical depths for each 
   !       of the 16 spectral bands;
   !    e) calls RTRNMC (for both clear and cloudy profiles) to perform the
   !       radiative transfer calculation;
   !    f) passes the necessary fluxes back to GCM.
   !
   ! Random number generation is by KISSVEC.
   !
   ! Cloud fraction and cloud physical properties are input and cloud optical
   ! properties are calculated by cldprmc based on iceflglw and liqflglw.
   ! Ice particle size provided must be appropriately defined for the ice
   ! parameterization selected. 
   !
   ! Aerosol optical depth is input directly by layer and spectral band
   ! as band average optical depth at mid-point of each spectral band.
   ! RRTMG_LW currently treats only aerosol absorption; scattering
   ! capability is not presently available.
   !
   ! The optional calculation of the change in upward flux as a function of
   ! surface temperature is available (controlled by input flag idrv). This
   ! can be utilized to approximate adjustments to the upward flux profile
   ! caused only by a change in surface temperature between full radiation
   ! calls. This feature uses the pre-calculated derivative of the Planck
   ! function with respect to surface temperature. 
   ! (1) Normal forward calculation for the input profile (idrv=0)
   ! (2) Normal forward calculation with optional calculation of the change
   !     in upward flux as a function of surface temperature for clear sky
   !     and total sky flux. Flux partial derivatives are provided in arrays
   !     duflx_dt and duflxc_dt for total and clear sky. (idrv=1)
   !
   ! ------- Modifications -------
   !
   ! This version of RRTMG_LW has been modified from RRTM_LW to use a reduced 
   ! set of g-points for application to GCMs.  
   !
   !-- Original version (derived from RRTM_LW), reduction of g-points, other
   !   revisions for use with GCMs.  
   !     1999: M. J. Iacono, AER, Inc.
   !-- Adapted for use with NCAR/CAM.
   !     May 2004: M. J. Iacono, AER, Inc.
   !-- Revised to add McICA capability. 
   !     Nov 2005: M. J. Iacono, AER, Inc.
   !-- Conversion to F90 formatting for consistency with rrtmg_sw.
   !     Feb 2007: M. J. Iacono, AER, Inc.
   !-- Modifications to formatting to use assumed-shape arrays.
   !     Aug 2007: M. J. Iacono, AER, Inc.
   !-- Modified to add longwave aerosol absorption.
   !     Apr 2008: M. J. Iacono, AER, Inc.
   !-- Added capability to calculate derivative of upward flux wrt surface temperature. 
   !     Nov 2009: M. J. Iacono, E. J. Mlawer, AER, Inc.
   !-- Added capability to run on GPU
   !     Aug 2012: David Berthiaume, AER, Inc.
   !-- Added alat and dyofyr to echo cloud-overlap scheme in RRTMGPU_SW
   !     Dec 2013: Matt Thompson, NASA/GMAO
   !-- Added numCPUs to allow for multple CPUs per GPU.
   !     Jan 2014: Matt Thompson, NASA/GMAO
   !-- Added new routine rtrnzero to help speed up zeroing on CPUs.
   !     Mar 2014: Matt Thompson, NASA/GMAO
   !-- Added band 10 OLR and d/dT for water vapor channel simulation
   !     Mar 2019: Peter Norris, NASA/GMAO
   !-- Major cleanup, re-org, and de-GPU
   !     Apr 2020: Peter Norris, NASA/GMAO
   ! -----------------------------------------------------------------------------

      use parrrtm, only: nbndlw

      ! ----- Input -----

      ! dimensions
      integer, intent(in) :: ncol   ! Number of horizontal columns
      integer, intent(in) :: nlay   ! Number of model layers
      integer, intent(in) :: psize  ! column partitioning size

      ! for requesting upwards flux derivatives wrt surface tempeerature
      integer, intent(in) :: idrv   ! 0: normal, 1: adds duflx_dt and duflxc_dt

      ! column characterization
      real,    intent(in) :: play (ncol,nlay)    ! Layer pressures [hPa] 
      real,    intent(in) :: plev (ncol,0:nlay)  ! Interface pressures [hPa]
      real,    intent(in) :: tlay (ncol,nlay)    ! Layer temperatures [K]
      real,    intent(in) :: tlev (ncol,0:nlay)  ! Interface temperatures [K]

      ! surface characterization
      real,    intent(in) :: tsfc (ncol)         ! Surface temperature [K]
      real,    intent(in) :: emis (ncol,nbndlw)  ! Surface emissivity

      ! column gaseous properties
      ! Note: All volume mixing ratios are in dimensionless units of mole fraction
      ! obtained by scaling mass mixing ratio (g/g) with the appropriate molecular
      ! weights (g/mol). The mole fraction is with respect to *dry air*.
      real,    intent(in) :: h2ovmr   (ncol,nlay)  ! H2O volume mixing ratio
      real,    intent(in) :: o3vmr    (ncol,nlay)  ! O3 volume mixing ratio
      real,    intent(in) :: co2vmr   (ncol,nlay)  ! CO2 volume mixing ratio
      real,    intent(in) :: ch4vmr   (ncol,nlay)  ! Methane volume mixing ratio
      real,    intent(in) :: n2ovmr   (ncol,nlay)  ! Nitrous oxide vol mixing ratio
      real,    intent(in) :: o2vmr    (ncol,nlay)  ! Oxygen volume mixing ratio
      real,    intent(in) :: cfc11vmr (ncol,nlay)  ! CFC11 volume mixing ratio
      real,    intent(in) :: cfc12vmr (ncol,nlay)  ! CFC12 volume mixing ratio
      real,    intent(in) :: cfc22vmr (ncol,nlay)  ! CFC22 volume mixing ratio
      real,    intent(in) :: ccl4vmr  (ncol,nlay)  ! CCL4 volume mixing ratio

      ! cloud physical properties
      real,    intent(in) :: cldf     (ncol,nlay)  ! Cloud fraction
      real,    intent(in) :: ciwp     (ncol,nlay)  ! In-cloud ice water path [g/m2]
      real,    intent(in) :: clwp     (ncol,nlay)  ! In-cloud liquid water path [g/m2]
      real,    intent(in) :: rei      (ncol,nlay)  ! Cloud ice crystal effective size [um]
      real,    intent(in) :: rel      (ncol,nlay)  ! Cloud droplet effective radius [um]

      ! -- NOTE: specific definition of rei depends on setting of iceflglw -------
      ! iceflglw = 0: ice effective radius, r_ec, (Ebert and Curry, 1992),
      !               r_ec must be >= 10.0 microns
      ! iceflglw = 1: ice effective radius, r_ec, (Ebert and Curry, 1992),
      !               r_ec range is limited to 13.0 to 130.0 microns
      ! iceflglw = 2: ice effective radius, r_k, (Key, Streamer Ref. Manual, 1996)
      !               r_k range is limited to 5.0 to 131.0 microns
      ! iceflglw = 3: generalized effective size, dge, (Fu, 1996),
      !               dge range is limited to 5.0 to 140.0 microns
      !               [dge = 1.0315 * r_ec]
      ! --------------------------------------------------------------------------

      ! flags for conversion of physical to optical cloud properties
      integer, intent(in) :: iceflglw  ! for ice crystals
      integer, intent(in) :: liqflglw  ! for liquid droplets

      ! Aerosol optical depth at mid-point of LW spectral bands
      real,    intent(in) :: tauaer (ncol,nlay,nbndlw)

      ! for cloud overlap calculations
      real,    intent(in) :: zm   (ncol,nlay)  ! Heights of level midpoints
      real,    intent(in) :: alat (ncol)       ! Latitude of column
      integer, intent(in) :: dyofyr            ! Day of the year

      ! pressure super-layer interface levels for cloud fractions
      integer, intent(in) :: cloudLM  ! Low-mid
      integer, intent(in) :: cloudMH  ! Mid-high

      ! ----- Output -----

      ! subcolumn clear counts for Tot|High|Mid|Low bands
      integer, intent(out) :: clearCounts(ncol,4)

      real, intent(out) :: uflx  (ncol,nlay+1)  ! Total sky LW upward flux [W/m2]
      real, intent(out) :: dflx  (ncol,nlay+1)  ! Total sky LW downward flux [W/m2]
      real, intent(out) :: uflxc (ncol,nlay+1)  ! Clear sky LW upward flux [W/m2]
      real, intent(out) :: dflxc (ncol,nlay+1)  ! Clear sky LW downward flux [W/m2]

      ! change in upward longwave flux wrt surface temperature [W/m2/K]
      real, intent(out) :: duflx_dt  (ncol,nlay+1)  ! total sky
      real, intent(out) :: duflxc_dt (ncol,nlay+1)  ! clear sky

      ! OLR for bands 9-11 and temperature derivatives [W/m2, W/m2/K]
      real, intent(out), dimension(ncol) :: olrb06, dolrb06_dt
      real, intent(out), dimension(ncol) :: olrb09, dolrb09_dt
      real, intent(out), dimension(ncol) :: olrb10, dolrb10_dt
      real, intent(out), dimension(ncol) :: olrb11, dolrb11_dt

      ! ----- Locals -----
      integer :: n, nparts

      ! ----------------------------------
      ! ASSERTs to catch unphysical inputs
      ! ----------------------------------
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
      if (any(tlev   < 0.)) then
        write(error_unit,*) 'file:', __FILE__, ', line:', __LINE__
        write(error_unit,*) 'minval(tlev):', minval(tlev)
        error stop 'negative values in input: tlev'
      end if
      if (any(tsfc   < 0.)) then
        write(error_unit,*) 'file:', __FILE__, ', line:', __LINE__
        write(error_unit,*) 'minval(tsfc):', minval(tsfc)
        error stop 'negative values in input: tsfc'
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
      if (any(cfc11vmr < 0.)) then
        write(error_unit,*) 'file:', __FILE__, ', line:', __LINE__
        write(error_unit,*) 'minval(cfc11vmr):', minval(cfc11vmr)
        error stop 'negative values in input: cfc11vmr'
      end if
      if (any(cfc12vmr < 0.)) then
        write(error_unit,*) 'file:', __FILE__, ', line:', __LINE__
        write(error_unit,*) 'minval(cfc12vmr):', minval(cfc12vmr)
        error stop 'negative values in input: cfc12vmr'
      end if
      if (any(cfc22vmr < 0.)) then
        write(error_unit,*) 'file:', __FILE__, ', line:', __LINE__
        write(error_unit,*) 'minval(cfc22vmr):', minval(cfc22vmr)
        error stop 'negative values in input: cfc22vmr'
      end if
      if (any(ccl4vmr < 0.)) then
        write(error_unit,*) 'file:', __FILE__, ', line:', __LINE__
        write(error_unit,*) 'minval(ccl4vmr):', minval(ccl4vmr)
        error stop 'negative values in input: ccl4vmr'
      end if
      if (any(emis    < 0.)) then
        write(error_unit,*) 'file:', __FILE__, ', line:', __LINE__
        write(error_unit,*) 'minval(emis):', minval(emis)
        error stop 'negative values in input: emis'
      end if
      if (any(cldf < 0.)) then
        write(error_unit,*) 'file:', __FILE__, ', line:', __LINE__
        write(error_unit,*) 'minval(cldf):', minval(cldf)
        error stop 'negative values in input: cldf'
      end if
      if (any(ciwp    < 0.)) then
        write(error_unit,*) 'file:', __FILE__, ', line:', __LINE__
        write(error_unit,*) 'minval(ciwp):', minval(ciwp)
        error stop 'negative values in input: ciwp'
      end if
      if (any(clwp    < 0.)) then
        write(error_unit,*) 'file:', __FILE__, ', line:', __LINE__
        write(error_unit,*) 'minval(clwp):', minval(clwp)
        error stop 'negative values in input: clwp'
      end if
      if (any(rei     < 0.)) then
        write(error_unit,*) 'file:', __FILE__, ', line:', __LINE__
        write(error_unit,*) 'minval(rei):', minval(rei)
        error stop 'negative values in input: rei'
      end if
      if (any(rel     < 0.)) then
        write(error_unit,*) 'file:', __FILE__, ', line:', __LINE__
        write(error_unit,*) 'minval(rel):', minval(rel)
        error stop 'negative values in input: rel'
      end if
      if (any(tauaer < 0.)) then
        write(error_unit,*) 'file:', __FILE__, ', line:', __LINE__
        write(error_unit,*) 'minval(tauaer):', minval(tauaer)
        error stop 'negative values in input: tauaer'
      end if

      ! set up condensate inhomogeneity tables
      call initialize_inhomogeneity(1)
! pmn: put in GCM init eventually

      ! ---------------------------------
      ! partition columns for performance
      ! ---------------------------------
      nparts = ceiling(real(ncol)/real(psize))

      do n = 0,nparts-1

         call rrtmg_lw_part (nparts, ncol, &
            n * psize + 1, min(psize, ncol - n * psize), &
            nlay, idrv, play, plev, tlay, tlev, tsfc, emis, & 
            h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr, &
            cfc11vmr, cfc12vmr, cfc22vmr, ccl4vmr, &
            cldf, ciwp, clwp, rei, rel, iceflglw, liqflglw, &
            tauaer, zm, alat, dyofyr, cloudLM, cloudMH, clearCounts, &
            uflx, dflx, uflxc, dflxc, duflx_dt, duflxc_dt, &
            olrb06, dolrb06_dt, olrb09, dolrb09_dt, &
            olrb10, dolrb10_dt, olrb11, dolrb11_dt)

      end do

      ! release condensate inhomogeneity resources
      call release_inhomogeneity

   end subroutine rrtmg_lw


   ! ------------------------------------------------------------------
   subroutine rrtmg_lw_part( &
      nparts, ncol, colstart, pncol, nlay, idrv, &
      play, plev, tlay, tlev, tsfc, emis, & 
      h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr, &
      cfc11vmr, cfc12vmr, cfc22vmr, ccl4vmr, &
      cldf, ciwp, clwp, rei, rel, iceflglw, liqflglw, &
      tauaer, zm, alat, dyofyr, cloudLM, cloudMH, clearCounts, &
      uflx, dflx, uflxc, dflxc, duflx_dt, duflxc_dt, &
      olrb06, dolrb06_dt, olrb09, dolrb09_dt, &
      olrb10, dolrb10_dt, olrb11, dolrb11_dt)
   ! ------------------------------------------------------------------

      use parrrtm, only : nbndlw, ngptlw
      use rrlw_con, only: fluxfac, oneminus, pi
      use rrlw_wvn, only: ngb
      use cloud_subcol_gen, only: generate_stochastic_clouds, clearCounts_threeBand

      ! ----- Input -----

      ! dimensions etc.
      integer, intent(in) :: nparts    ! Number of partitions
      integer, intent(in) :: ncol      ! Number of horizontal columns
      integer, intent(in) :: pncol     ! Number of columns in this partition
      integer, intent(in) :: colstart  ! Starting column of this partition
      integer, intent(in) :: nlay      ! Number of model layers
      integer, intent(in) :: idrv      ! 0: normal, 1: adds duflx_dt and duflxc_dt

      ! column characterization
      real,    intent(in) :: play (ncol,nlay)    ! Layer pressures [hPa] 
      real,    intent(in) :: plev (ncol,0:nlay)  ! Interface pressures [hPa]
      real,    intent(in) :: tlay (ncol,nlay)    ! Layer temperatures [K]
      real,    intent(in) :: tlev (ncol,0:nlay)  ! Interface temperatures [K]

      ! surface characterization
      real,    intent(in) :: tsfc (ncol)         ! Surface temperature [K]
      real,    intent(in) :: emis (ncol,nbndlw)  ! Surface emissivity

      ! column gaseous properties
      real,    intent(in) :: h2ovmr   (ncol,nlay)  ! H2O volume mixing ratio
      real,    intent(in) :: o3vmr    (ncol,nlay)  ! O3 volume mixing ratio
      real,    intent(in) :: co2vmr   (ncol,nlay)  ! CO2 volume mixing ratio
      real,    intent(in) :: ch4vmr   (ncol,nlay)  ! Methane volume mixing ratio
      real,    intent(in) :: n2ovmr   (ncol,nlay)  ! Nitrous oxide vol mixing ratio
      real,    intent(in) :: o2vmr    (ncol,nlay)  ! Oxygen volume mixing ratio
      real,    intent(in) :: cfc11vmr (ncol,nlay)  ! CFC11 volume mixing ratio
      real,    intent(in) :: cfc12vmr (ncol,nlay)  ! CFC12 volume mixing ratio
      real,    intent(in) :: cfc22vmr (ncol,nlay)  ! CFC22 volume mixing ratio
      real,    intent(in) :: ccl4vmr  (ncol,nlay)  ! CCL4 volume mixing ratio

      ! cloud physical properties
      real,    intent(in) :: cldf     (ncol,nlay)  ! Cloud fraction
      real,    intent(in) :: ciwp     (ncol,nlay)  ! In-cloud ice water path [g/m2]
      real,    intent(in) :: clwp     (ncol,nlay)  ! In-cloud liquid water path [g/m2]
      real,    intent(in) :: rei      (ncol,nlay)  ! Cloud ice crystal effective size [um]
      real,    intent(in) :: rel      (ncol,nlay)  ! Cloud droplet effective radius [um]

      ! flags for conversion of physical to optical cloud properties
      integer, intent(in) :: iceflglw  ! For ice crystals
      integer, intent(in) :: liqflglw  ! For liquid droplets

      ! Aerosol optical depth at mid-point of LW spectral bands
      real,    intent(in) :: tauaer (ncol,nlay,nbndlw)

      ! for cloud overlap calculations
      real,    intent(in) :: zm   (ncol,nlay)  ! Heights of level midpoints
      real,    intent(in) :: alat (ncol)       ! Latitude of column
      integer, intent(in) :: dyofyr            ! Day of the year

      ! pressure super-layer interface levels for cloud fractions
      integer, intent(in) :: cloudLM  ! Low-mid
      integer, intent(in) :: cloudMH  ! Mid-high

      ! ----- Output -----

      ! subcolumn clear counts for Tot|High|Mid|Low bands
      integer, intent(out) :: clearCounts(ncol,4)

      real, intent(out) :: uflx  (ncol,nlay+1)  ! Total sky LW upward flux [W/m2]
      real, intent(out) :: dflx  (ncol,nlay+1)  ! Total sky LW downward flux [W/m2]
      real, intent(out) :: uflxc (ncol,nlay+1)  ! Clear sky LW upward flux [W/m2]
      real, intent(out) :: dflxc (ncol,nlay+1)  ! Clear sky LW downward flux [W/m2]

      ! change in upward longwave flux wrt surface temperature [W/m2/K]
      real, intent(out) :: duflx_dt  (ncol,nlay+1)  ! total sky
      real, intent(out) :: duflxc_dt (ncol,nlay+1)  ! clear sky

      ! OLR for bands 9-11 and temperature derivatives [W/m2, W/m2/K]
      real, intent(out), dimension(ncol) :: olrb06, dolrb06_dt
      real, intent(out), dimension(ncol) :: olrb09, dolrb09_dt
      real, intent(out), dimension(ncol) :: olrb10, dolrb10_dt
      real, intent(out), dimension(ncol) :: olrb11, dolrb11_dt

      ! local partitioned variables ...
      ! ===============================

      ! general
      real :: p_zm   (pncol,nlay)    ! mid-layer heights
      real :: p_alat (pncol)         ! latitudes
      real :: p_play (pncol,nlay)    ! layer pressures [hPa]
      real :: p_tlay (pncol,nlay)    ! layer temperatures [K]
      real :: p_plev (pncol,0:nlay)  ! level (interface) pressures [hPa]
      real :: p_tlev (pncol,0:nlay)  ! level (interface) temperatures [K]
      real :: p_tsfc (pncol)         ! surface temperature [K]
      real :: p_emis (pncol,nbndlw)  ! lw surface emissivity

      ! molecular volume mixing ratios
      real, dimension (pncol,nlay) :: &
         p_h2ovmr, p_o3vmr, p_co2vmr, p_ch4vmr, p_n2ovmr, p_o2vmr, p_covmr, &
         p_cfc11vmr, p_cfc12vmr, p_cfc22vmr, p_ccl4vmr

      ! cloud input profiles
      real :: p_cldf (pncol,nlay)  ! Cloud fraction
      real :: p_ciwp (pncol,nlay)  ! In-cloud ice water path [g/m2]
      real :: p_clwp (pncol,nlay)  ! In-cloud liquid water path [g/m2]
      real :: p_rei  (pncol,nlay)  ! Cloud ice particle effective size [um]
      real :: p_rel  (pncol,nlay)  ! Cloud water drop effective radius [um]

      ! Aerosol optical depth at mid-point of LW spectral bands
      real :: p_tauaer (pncol,nlay,nbndlw)

      ! gas optical depths and Planck fractions
      real :: taug   (pncol,nlay,ngptlw)  ! gas + aerosol optical depth
      real :: pfracs (pncol,nlay,ngptlw)  ! Planck fractions

      ! mcica generated clouds
      real :: cldfmc (pncol,ngptlw,nlay)  ! cloud fraction
      real :: ciwpmc (pncol,ngptlw,nlay)  ! cloud ice water path [g/m2]
      real :: clwpmc (pncol,ngptlw,nlay)  ! cloud liq water path [g/m2]
      real :: taucmc (pncol,ngptlw,nlay)  ! cloud optical depth
      integer :: p_clearCounts(pncol,4)   ! for super-band cld fractions

      ! cloudy for ANY subcol of column? [0=no,1=yes]
      integer :: icldlyr (pncol,nlay)

      ! spectrally summed fluxes and upward flux derivatives wrt Tsurf
      real :: totuflux     (pncol,0:nlay)  ! upward longwave flux (W/m2)
      real :: totdflux     (pncol,0:nlay)  ! downward longwave flux (W/m2)
      real :: totuclfl     (pncol,0:nlay)  ! clrsky upward lw flux (W/m2)
      real :: totdclfl     (pncol,0:nlay)  ! clrsky downward lw flux (W/m2)
      real :: dtotuflux_dt (pncol,0:nlay)  ! d/d(Tsurf) (W/m2/K)
      real :: dtotuclfl_dt (pncol,0:nlay)  ! d/d(Tsurf) (W/m2/K)

      ! TOA OLR in bands 6 & 9-11 and their derivatives wrt Tsurf
      real :: p_olrb06     (pncol)  ! (W/m2)
      real :: p_olrb09     (pncol)  ! (W/m2)
      real :: p_olrb10     (pncol)  ! (W/m2)
      real :: p_olrb11     (pncol)  ! (W/m2)
      real :: p_dolrb06_dt (pncol)  ! (W/m2/K)
      real :: p_dolrb09_dt (pncol)  ! (W/m2/K)
      real :: p_dolrb10_dt (pncol)  ! (W/m2/K)
      real :: p_dolrb11_dt (pncol)  ! (W/m2/K)

      ! keep at one (handle to unused special treatment of band 16
      ! when it is set to 16 (see rrtmg_lw_setcoef() for details))
      integer, parameter :: istart = 1

      ! set some globals
      oneminus = 1. - 1.e-6 
      pi = 2. * asin(1. )
      fluxfac = pi * 2.e4                   ! orig:   fluxfac = pi * 2.d4  
! pmn: put this in init of make parameters in _con

      ! copy partition
      p_zm       = zm       (colstart:(colstart+pncol-1),:)
      p_alat     = alat     (colstart:(colstart+pncol-1))
      p_play     = play     (colstart:(colstart+pncol-1),:)
      p_tlay     = tlay     (colstart:(colstart+pncol-1),:)
      p_plev     = plev     (colstart:(colstart+pncol-1),:)
      p_tlev     = tlev     (colstart:(colstart+pncol-1),:)
      p_tsfc     = tsfc     (colstart:(colstart+pncol-1))
      p_emis     = emis     (colstart:(colstart+pncol-1),:)
      p_h2ovmr   = h2ovmr   (colstart:(colstart+pncol-1),:)
      p_o3vmr    = o3vmr    (colstart:(colstart+pncol-1),:)
      p_co2vmr   = co2vmr   (colstart:(colstart+pncol-1),:)
      p_ch4vmr   = ch4vmr   (colstart:(colstart+pncol-1),:)
      p_n2ovmr   = n2ovmr   (colstart:(colstart+pncol-1),:)
      p_o2vmr    = o2vmr    (colstart:(colstart+pncol-1),:)
      p_covmr    = 0.
      p_cfc11vmr = cfc11vmr (colstart:(colstart+pncol-1),:)
      p_cfc12vmr = cfc12vmr (colstart:(colstart+pncol-1),:)
      p_cfc22vmr = cfc22vmr (colstart:(colstart+pncol-1),:)
      p_ccl4vmr  = ccl4vmr  (colstart:(colstart+pncol-1),:)
! pmn: consider adding CO absorption since it is calculable !!!!!!!!!!!!!
! pmn: i.e., pass through from RRTMG_LW
      p_cldf     = cldf     (colstart:(colstart+pncol-1),:)
      p_ciwp     = ciwp     (colstart:(colstart+pncol-1),:)
      p_clwp     = clwp     (colstart:(colstart+pncol-1),:)
      p_rei      = rei      (colstart:(colstart+pncol-1),:)
      p_rel      = rel      (colstart:(colstart+pncol-1),:)
      p_tauaer   = tauaer   (colstart:(colstart+pncol-1),:,:)

      ! Call model and data initialization, compute lookup tables, perform
      ! In a GCM this call should be placed in the model initialization
      ! area, since this has to be called only once.  
      ! call rrtmg_lw_ini()

      ! Generate stochastic subcolumns of cloud physical properties

      call generate_stochastic_clouds( &
         pncol, ngptlw, nlay, &
         p_zm, p_alat, dyofyr, &
         p_play, p_cldf, p_ciwp, p_clwp, &
         cldfmc, clwpmc, ciwpmc)

      ! for super-band cloud fractions

      call clearCounts_threeBand( &
         pncol, ngptlw, nlay, cloudLM, cloudMH, cldfmc, &
         p_clearCounts)

      ! cloud physical to physical properties

      call cldprmc (pncol, nlay, &
        cldfmc, ciwpmc, clwpmc, p_rei, p_rel, &
        iceflglw, liqflglw, ngb, &
        taucmc, icldlyr)

      ! Calculate information needed by the radiative transfer routine
      ! that is specific to this atmosphere, especially some of the 
      ! coefficients and indices needed to compute the optical depths
      ! by interpolating data from stored reference atmospheres. 

      call setcoef (pncol, nlay, istart, idrv, &
         p_play, p_tlay, p_plev, p_tlev, p_tsfc, p_emis, &
         p_h2ovmr, p_o3vmr, p_co2vmr, p_ch4vmr, p_n2ovmr, p_o2vmr, p_covmr, &
         p_cfc11vmr, p_cfc12vmr, p_cfc22vmr, p_ccl4vmr)

      ! Calculate the gaseous optical depths and Planck fractions for 
      ! each longwave spectral band. Also adds in aerosol optical depths.

      call taumol (pncol, nlay, ngb, p_play, p_tauaer, taug, pfracs)

      ! Call the radiative transfer routine

      call rtrnmc (pncol, nlay, idrv, ngb, &
         p_emis, taug, pfracs, icldlyr, cldfmc, taucmc, &
         totuflux, totdflux, totuclfl, totdclfl, &
         dtotuflux_dt, dtotuclfl_dt, &
         p_olrb06, p_olrb09, p_olrb10, p_olrb11, &
         p_dolrb06_dt, p_dolrb09_dt, p_dolrb10_dt, p_dolrb11_dt)

      ! copy the partitioned results back
      clearCounts (colstart:(colstart+pncol-1),:) = p_clearCounts
      uflx  (colstart:(colstart+pncol-1),1:(nlay+1)) = totuflux(:,0:nlay)
      dflx  (colstart:(colstart+pncol-1),1:(nlay+1)) = totdflux(:,0:nlay)
      uflxc (colstart:(colstart+pncol-1),1:(nlay+1)) = totuclfl(:,0:nlay)
      dflxc (colstart:(colstart+pncol-1),1:(nlay+1)) = totdclfl(:,0:nlay)
      olrb06(colstart:(colstart+pncol-1)) = p_olrb06(:)
      olrb09(colstart:(colstart+pncol-1)) = p_olrb09(:)
      olrb10(colstart:(colstart+pncol-1)) = p_olrb10(:)
      olrb11(colstart:(colstart+pncol-1)) = p_olrb11(:)
      if (idrv == 1) then
         duflx_dt (colstart:(colstart+pncol-1),1:(nlay+1)) = dtotuflux_dt(:,0:nlay)
         duflxc_dt(colstart:(colstart+pncol-1),1:(nlay+1)) = dtotuclfl_dt(:,0:nlay)
         dolrb06_dt(colstart:(colstart+pncol-1)) = p_dolrb06_dt(:)
         dolrb09_dt(colstart:(colstart+pncol-1)) = p_dolrb09_dt(:)
         dolrb10_dt(colstart:(colstart+pncol-1)) = p_dolrb10_dt(:)
         dolrb11_dt(colstart:(colstart+pncol-1)) = p_dolrb11_dt(:)
      end if
! pmn: these copies can be cleaned up RHS without ()

      call setcoef_free

   end subroutine rrtmg_lw_part

end module rrtmg_lw_rad
