module rrtmg_lw_rad

   use cloud_condensate_inhomogeneity_lw, only: &
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
      ncol, nlay, psize, dudTs, &
      play, plev, tlay, tlev, tsfc, emis, & 
      h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr, &
      cfc11vmr, cfc12vmr, cfc22vmr, ccl4vmr, &
      cldf, ciwp, clwp, rei, rel, iceflglw, liqflglw, &
      tauaer, zm, alat, dyofyr, cloudLM, cloudMH, clearCounts, &
      uflx, dflx, uflxc, dflxc, duflx_dTs, duflxc_dTs, &
      band_output, olrb, dolrb_dTs)
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
   ! The optional calculation of the change in upward flux as a function
   ! of surface temperature is controlled by input flag dudTs. This can be
   ! utilized to approximate adjustments to the upward flux profile caused
   ! only by a change in surface temperature between full radiation calls.
   ! This feature uses the pre-calculated derivative of the Planck function
   ! with respect to surface temperature: 
   ! (dudTs false) Normal forward calculation;
   ! (dudTs true ) Normal forward calculation with optional calculation
   !    of the change in upward flux as a function of surface temperature
   !    for total-sky and clear-sky flux. Flux partial derivatives are
   !    provided in arrays duflx_dTs and duflxc_dTs respectively.
   ! If dudTs false, DO NOT use any d_dTs, as values will be undefined.
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

      ! for requesting upwards flux derivatives wrt Tsurf
      logical, intent(in) :: dudTs

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
      real, intent(out) :: duflx_dTs  (ncol,nlay+1)  ! total sky
      real, intent(out) :: duflxc_dTs (ncol,nlay+1)  ! clear sky

      ! ----- band OLRs -----

      ! which band OLRs to calculate?
      logical, intent(in) :: band_output (nbndlw)

      ! band OLRs and d/dTs
      real, intent(out) :: olrb      (nbndlw,ncol)  ! [W/m2]
      real, intent(out) :: dolrb_dTs (nbndlw,ncol)  ! [W/m2/K]

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

      ! Call model and data initialization, compute lookup tables, perform
      ! In a GCM this call should be placed in the model initialization
      ! area, since this has to be called only once.  
      ! call rrtmg_lw_ini()

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
            nlay, dudTs, play, plev, tlay, tlev, tsfc, emis, & 
            h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr, &
            cfc11vmr, cfc12vmr, cfc22vmr, ccl4vmr, &
            cldf, ciwp, clwp, rei, rel, iceflglw, liqflglw, &
            tauaer, zm, alat, dyofyr, cloudLM, cloudMH, clearCounts, &
            uflx, dflx, uflxc, dflxc, duflx_dTs, duflxc_dTs, &
            band_output, olrb, dolrb_dTs)

      end do

      ! release condensate inhomogeneity resources
      call release_inhomogeneity

   end subroutine rrtmg_lw


   ! ------------------------------------------------------------------
   subroutine rrtmg_lw_part( &
      nparts, ncol, colstart, pncol, nlay, dudTs, &
      play, plev, tlay, tlev, tsfc, emis, & 
      h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr, &
      cfc11vmr, cfc12vmr, cfc22vmr, ccl4vmr, &
      cldf, ciwp, clwp, rei, rel, iceflglw, liqflglw, &
      tauaer, zm, alat, dyofyr, cloudLM, cloudMH, clearCounts, &
      uflx, dflx, uflxc, dflxc, duflx_dTs, duflxc_dTs, &
      band_output, olrb, dolrb_dTs)
   ! ------------------------------------------------------------------

      use parrrtm, only : nbndlw, ngptlw
      use cloud_subcol_gen_lw, only: &
         generate_stochastic_clouds, clearCounts_threeBand

      ! ----- Input -----

      ! dimensions etc.
      integer, intent(in) :: nparts    ! Number of partitions
      integer, intent(in) :: ncol      ! Number of horizontal columns
      integer, intent(in) :: pncol     ! Number of columns in this partition
      integer, intent(in) :: colstart  ! Starting column of this partition
      integer, intent(in) :: nlay      ! Number of model layers
      logical, intent(in) :: dudTs     ! true adds d/dTs derivs

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
      real, intent(out) :: duflx_dTs  (ncol,nlay+1)  ! total sky
      real, intent(out) :: duflxc_dTs (ncol,nlay+1)  ! clear sky

      ! ----- band OLRs -----

      ! which band OLRs to calculate?
      logical, intent(in) :: band_output (nbndlw)

      ! band OLRs and d/dTs
      real, intent(out) :: olrb      (nbndlw,ncol)  ! [W/m2]
      real, intent(out) :: dolrb_dTs (nbndlw,ncol)  ! [W/m2/K]

      ! local partitioned variables ...
      ! ===============================

      ! general
      real :: p_zm     (nlay,pncol)  ! mid-layer heights [m]
      real :: p_alat        (pncol)  ! latitudes
      real :: p_play   (nlay,pncol)  ! layer pressures [hPa]
      real :: p_tlay   (nlay,pncol)  ! layer temperatures [K]
      real :: p_plev (0:nlay,pncol)  ! level (interface) pressures [hPa]
      real :: p_tlev (0:nlay,pncol)  ! level (interface) temperatures [K]
      real :: p_tsfc        (pncol)  ! surface temperature [K]
      real :: p_emis (nbndlw,pncol)  ! lw surface emissivity

      ! molecular volume mixing ratios
      real, dimension (nlay,pncol) :: &
         p_h2ovmr, p_o3vmr, p_co2vmr, p_ch4vmr, p_n2ovmr, p_o2vmr, p_covmr, &
         p_cfc11vmr, p_cfc12vmr, p_cfc22vmr, p_ccl4vmr

      ! cloud input profiles
      real :: p_cldf (nlay,pncol)  ! Cloud fraction
      real :: p_ciwp (nlay,pncol)  ! In-cloud ice water path [g/m2]
      real :: p_clwp (nlay,pncol)  ! In-cloud liquid water path [g/m2]
      real :: p_rei  (nlay,pncol)  ! Cloud ice particle effective size [um]
      real :: p_rel  (nlay,pncol)  ! Cloud water drop effective radius [um]

      ! Aerosol optical depth at mid-point of LW spectral bands
      real :: p_tauaer (nlay,nbndlw,pncol)

      ! gas optical depths and Planck fractions
      real :: taug   (nlay,ngptlw,pncol)  ! gas + aerosol optical depth
      real :: pfracs (nlay,ngptlw,pncol)  ! Planck fractions

      ! mcica generated clouds
      logical :: cldymc (nlay,ngptlw,pncol)  ! cloudy or not?
      real    :: ciwpmc (nlay,ngptlw,pncol)  ! cloud ice water path [g/m2]
      real    :: clwpmc (nlay,ngptlw,pncol)  ! cloud liq water path [g/m2]
      real    :: taucmc (nlay,ngptlw,pncol)  ! cloud optical depth
      integer :: p_clearCounts (4,pncol)     ! for super-band cld fractions

      ! cloudy for ANY subcol/gpoint of column?
      logical :: cloudy (nlay,pncol)

      ! spectrally summed fluxes and upward flux derivatives wrt Tsurf
      real :: totuflux      (0:nlay,pncol)  ! upward longwave flux [W/m2]
      real :: totdflux      (0:nlay,pncol)  ! downward longwave flux [W/m2]
      real :: totuclfl      (0:nlay,pncol)  ! clrsky upward lw flux [W/m2]
      real :: totdclfl      (0:nlay,pncol)  ! clrsky downward lw flux [W/m2]
      real :: dtotuflux_dTs (0:nlay,pncol)  ! d/d(Tsurf) [W/m2/K]
      real :: dtotuclfl_dTs (0:nlay,pncol)  ! d/d(Tsurf) [W/m2/K]

      ! band OLRs and d/dTs
      real :: p_olrb      (nbndlw,pncol)  ! [W/m2]
      real :: p_dolrb_dTs (nbndlw,pncol)  ! [W/m2/K]

      ! keep at one (handle to unused special treatment of band 16
      ! when it is set to 16 (see rrtmg_lw_setcoef() for details))
      integer, parameter :: istart = 1

      ! indices
      integer :: ilev, ilay, ibnd, n

      ! copy partition and reorder for speed
      p_alat = alat (colstart:(colstart+pncol-1))
      p_tsfc = tsfc (colstart:(colstart+pncol-1))
      do ilay = 1,nlay
         p_zm       (ilay,:) = zm       (colstart:(colstart+pncol-1),ilay)
         p_play     (ilay,:) = play     (colstart:(colstart+pncol-1),ilay)
         p_tlay     (ilay,:) = tlay     (colstart:(colstart+pncol-1),ilay)
         p_cldf     (ilay,:) = cldf     (colstart:(colstart+pncol-1),ilay)
         p_ciwp     (ilay,:) = ciwp     (colstart:(colstart+pncol-1),ilay)
         p_clwp     (ilay,:) = clwp     (colstart:(colstart+pncol-1),ilay)
         p_rei      (ilay,:) = rei      (colstart:(colstart+pncol-1),ilay)
         p_rel      (ilay,:) = rel      (colstart:(colstart+pncol-1),ilay)
         p_h2ovmr   (ilay,:) = h2ovmr   (colstart:(colstart+pncol-1),ilay)
         p_o3vmr    (ilay,:) = o3vmr    (colstart:(colstart+pncol-1),ilay)
         p_co2vmr   (ilay,:) = co2vmr   (colstart:(colstart+pncol-1),ilay)
         p_ch4vmr   (ilay,:) = ch4vmr   (colstart:(colstart+pncol-1),ilay)
         p_n2ovmr   (ilay,:) = n2ovmr   (colstart:(colstart+pncol-1),ilay)
         p_o2vmr    (ilay,:) = o2vmr    (colstart:(colstart+pncol-1),ilay)
         p_covmr    (ilay,:) = 0.
! pmn: consider adding CO absorption since it is calculable !!!!!!!!!!!!!
! pmn: i.e., pass through from RRTMG_LW
         p_cfc11vmr (ilay,:) = cfc11vmr (colstart:(colstart+pncol-1),ilay)
         p_cfc12vmr (ilay,:) = cfc12vmr (colstart:(colstart+pncol-1),ilay)
         p_cfc22vmr (ilay,:) = cfc22vmr (colstart:(colstart+pncol-1),ilay)
         p_ccl4vmr  (ilay,:) = ccl4vmr  (colstart:(colstart+pncol-1),ilay)
      end do
      do ilev = 0,nlay
         p_plev (ilev,:) = plev (colstart:(colstart+pncol-1),ilev)
         p_tlev (ilev,:) = tlev (colstart:(colstart+pncol-1),ilev)
      end do
      do ibnd = 1,nbndlw
         p_emis (ibnd,:) = emis (colstart:(colstart+pncol-1),ibnd)
         do ilay = 1,nlay
            p_tauaer (ilay,ibnd,:) = tauaer (colstart:(colstart+pncol-1),ilay,ibnd)
        end do
      end do

      ! Generate stochastic subcolumns of cloud physical properties

      call generate_stochastic_clouds( &
         pncol, pncol, ngptlw, nlay, &
         p_zm, p_alat, dyofyr, &
         p_play, p_cldf, p_ciwp, p_clwp, 1.e-20, &
         cldymc, ciwpmc, clwpmc)

      ! for super-band cloud fractions

      call clearCounts_threeBand( &
         pncol, pncol, ngptlw, nlay, cloudLM, cloudMH, cldymc, &
         p_clearCounts)
      do n = 1,4
         clearCounts (colstart:(colstart+pncol-1),n) = p_clearCounts(n,:)
      end do

      ! cloud physical to physical properties

      call cldprmc (pncol, nlay, &
        cldymc, ciwpmc, clwpmc, p_rei, p_rel, &
        iceflglw, liqflglw, taucmc, cloudy)

      ! Calculate information needed by the radiative transfer routine
      ! that is specific to this atmosphere, especially some of the 
      ! coefficients and indices needed to compute the optical depths
      ! by interpolating data from stored reference atmospheres. 

      call setcoef (pncol, nlay, istart, dudTs, &
         p_play, p_tlay, p_plev, p_tlev, p_tsfc, p_emis, &
         p_h2ovmr, p_o3vmr, p_co2vmr, p_ch4vmr, p_n2ovmr, p_o2vmr, p_covmr, &
         p_cfc11vmr, p_cfc12vmr, p_cfc22vmr, p_ccl4vmr)

      ! Calculate the gaseous optical depths and Planck fractions for 
      ! each longwave spectral band. Also adds in aerosol optical depths.

      call taumol (pncol, nlay, p_play, p_tauaer, taug, pfracs)

      ! Call the radiative transfer routine

      call rtrnmc (pncol, nlay, dudTs, &
         p_emis, taug, pfracs, cloudy, taucmc, &
         totuflux, totdflux, totuclfl, totdclfl, &
         dtotuflux_dTs, dtotuclfl_dTs, &
         band_output, p_olrb, p_dolrb_dTs)

      ! copy the partitioned fluxes back
      do ilev = 0,nlay
         uflx  (colstart:(colstart+pncol-1),ilev+1) = totuflux(ilev,:)
         dflx  (colstart:(colstart+pncol-1),ilev+1) = totdflux(ilev,:)
         uflxc (colstart:(colstart+pncol-1),ilev+1) = totuclfl(ilev,:)
         dflxc (colstart:(colstart+pncol-1),ilev+1) = totdclfl(ilev,:)
         if (dudTs) then
            duflx_dTs (colstart:(colstart+pncol-1),ilev+1) = dtotuflux_dTs(ilev,:)
            duflxc_dTs(colstart:(colstart+pncol-1),ilev+1) = dtotuclfl_dTs(ilev,:)
         end if
      end do

      ! band OLRs
      do ibnd = 1,nbndlw
         if (band_output(ibnd)) then
            olrb(ibnd,colstart:(colstart+pncol-1)) = p_olrb(ibnd,:)
            if (dudTs) &
               dolrb_dTs(ibnd,colstart:(colstart+pncol-1)) = p_dolrb_dTs(ibnd,:)
         end if
      end do

      ! free internal stae of setcoef
      call setcoef_free

   end subroutine rrtmg_lw_part

end module rrtmg_lw_rad
