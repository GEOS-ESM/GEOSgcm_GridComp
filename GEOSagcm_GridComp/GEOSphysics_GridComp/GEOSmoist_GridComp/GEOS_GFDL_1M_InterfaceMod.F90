! $Id$

#include "MAPL.h"

!=============================================================================
!BOP

! !MODULE: GEOS_GFDL_1M_InterfaceMod -- A Module to interface with the
!   GFDL_1M cloud microphysics

module GEOS_GFDL_1M_InterfaceMod

   use ESMF
   use MAPL
   use gFTL2_StringVector, only: StringVector
   use GEOS_UtilsMod
   use GEOS_RadarMod
   use GEOSmoist_Process_Library
   use Aer_Actv_Single_Moment
   use aer_cloud, only : aer_cloud_init
   use gfdl2_cloud_microphys_mod, only : gfdl_cloud_microphys_init, &
                                         gfdl_cloud_microphys_driver, &
                                         ICE_LSC_VFALL_PARAM, &
                                         ICE_CNV_VFALL_PARAM
   use gfdl_mp_mod, only : gfdl_mp_init, &
                           gfdl_mp_driver, &
                           do_ref, &
                           do_hail, &
                           do_sedi_heat, &
                           do_sedi_melt_qi, &
                           do_sedi_melt_qs, &
                           do_sedi_melt_qg, &
                           ifflag

   implicit none

   private

   integer :: STATUS

   ! Local resource variables
   real :: TURNRHCRIT_TOP, TURNRHCRIT_SFC
   real :: MIN_RH_CRIT, MAX_RH_CRIT, MIN_RH_UNSTABLE, MIN_RH_STABLE
   real :: TAU_EVAP, CCW_EVAP_EFF
   real :: TAU_SUBL, CCI_EVAP_EFF
   real :: ANV_ICEFALL
   real :: LS_ICEFALL
   real :: FAC_RL
   real :: MIN_RL
   real :: MAX_RL
   real :: FAC_RI
   real :: MIN_RI
   real :: MAX_RI
   logical :: LHYDROSTATIC
   logical :: LPHYS_HYDROSTATIC
   logical :: LMELTFRZ_CLDMACRO
   logical :: LMELTFRZ_CLDMICRO
   real :: GFDL_MP_KLID

   logical :: REPORT_GFDL_1M_NEGATIVES

   logical :: GFDL_MP3
   logical :: USE_PYMOIST_GFDL1M

   public :: GFDL_1M_Setup, GFDL_1M_Initialize, GFDL_1M_Run, GFDL_1M_Finalize

contains

   subroutine GFDL_1M_Setup(GC, ADVECTED, RC)
      type(ESMF_GridComp), intent(inout) :: GC ! gridded component
      type(StringVector), intent(inout) :: ADVECTED ! internals to subscribe to advection (TRADV)
      integer, optional :: RC ! return code

      ! !INTERNAL STATE:

      ! Advected by dynamics (MAPL2 FRIENDLYTO=DYNAMICS). MAPL2 also had
      ! Q, QLLS, QLCN, QILS, QICN, QRAIN, QSNOW, QGRAUPEL friendly to TURBULENCE
      ! and Q to CHEMISTRY:ANALYSIS; no MAPL3 service for those yet.
      call ADVECTED%push_back('Q')
      call ADVECTED%push_back('CLLS')
      call ADVECTED%push_back('CLCN')
      call ADVECTED%push_back('QLLS')
      call ADVECTED%push_back('QLCN')
      call ADVECTED%push_back('QILS')
      call ADVECTED%push_back('QICN')
      call ADVECTED%push_back('QRAIN')
      call ADVECTED%push_back('QSNOW')
      call ADVECTED%push_back('QGRAUPEL')

      !BOS

      call MAPL_GridCompAddSpec(GC, &
           STATE_INTENT=ESMF_STATEINTENT_INTERNAL, &
           SHORT_NAME='Q', &
           STANDARD_NAME='specific_humidity', &
           UNITS='kg kg-1', &
           DIMS='xyz', &
           VERTICAL_STAGGER=MAPL_VERTICAL_STAGGER_CENTER, &
           RESTART_MODE=MAPL_RESTART_REQUIRED, &
           FILL_VALUE=1.0e-6, &
           _RC)

      call MAPL_GridCompAddSpec(GC, &
           STATE_INTENT=ESMF_STATEINTENT_INTERNAL, &
           SHORT_NAME='QLLS', &
           STANDARD_NAME='mass_fraction_of_large_scale_cloud_liquid_water', &
           UNITS='kg kg-1', &
           DIMS='xyz', &
           VERTICAL_STAGGER=MAPL_VERTICAL_STAGGER_CENTER, &
           FILL_VALUE=0.0, &
           _RC)

      call MAPL_GridCompAddSpec(GC, &
           STATE_INTENT=ESMF_STATEINTENT_INTERNAL, &
           SHORT_NAME='QLCN', &
           STANDARD_NAME='mass_fraction_of_convective_cloud_liquid_water', &
           UNITS='kg kg-1', &
           DIMS='xyz', &
           VERTICAL_STAGGER=MAPL_VERTICAL_STAGGER_CENTER, &
           FILL_VALUE=0.0, &
           _RC)

      call MAPL_GridCompAddSpec(GC, &
           STATE_INTENT=ESMF_STATEINTENT_INTERNAL, &
           SHORT_NAME='CLLS', &
           STANDARD_NAME='large_scale_cloud_area_fraction', &
           UNITS='1', &
           DIMS='xyz', &
           VERTICAL_STAGGER=MAPL_VERTICAL_STAGGER_CENTER, &
           FILL_VALUE=0.0, &
           _RC)

      call MAPL_GridCompAddSpec(GC, &
           STATE_INTENT=ESMF_STATEINTENT_INTERNAL, &
           SHORT_NAME='CLCN', &
           STANDARD_NAME='convective_cloud_area_fraction', &
           UNITS='1', &
           DIMS='xyz', &
           VERTICAL_STAGGER=MAPL_VERTICAL_STAGGER_CENTER, &
           FILL_VALUE=0.0, &
           _RC)

      call MAPL_GridCompAddSpec(GC, &
           STATE_INTENT=ESMF_STATEINTENT_INTERNAL, &
           SHORT_NAME='QILS', &
           STANDARD_NAME='mass_fraction_of_large_scale_cloud_ice_water', &
           UNITS='kg kg-1', &
           DIMS='xyz', &
           VERTICAL_STAGGER=MAPL_VERTICAL_STAGGER_CENTER, &
           FILL_VALUE=0.0, &
           _RC)

      call MAPL_GridCompAddSpec(GC, &
           STATE_INTENT=ESMF_STATEINTENT_INTERNAL, &
           SHORT_NAME='QICN', &
           STANDARD_NAME='mass_fraction_of_convective_cloud_ice_water', &
           UNITS='kg kg-1', &
           DIMS='xyz', &
           VERTICAL_STAGGER=MAPL_VERTICAL_STAGGER_CENTER, &
           FILL_VALUE=0.0, &
           _RC)

      call MAPL_GridCompAddSpec(GC, &
           STATE_INTENT=ESMF_STATEINTENT_INTERNAL, &
           SHORT_NAME='QRAIN', &
           STANDARD_NAME='mass_fraction_of_rain', &
           UNITS='kg kg-1', &
           DIMS='xyz', &
           VERTICAL_STAGGER=MAPL_VERTICAL_STAGGER_CENTER, &
           FILL_VALUE=0.0, &
           _RC)

      call MAPL_GridCompAddSpec(GC, &
           STATE_INTENT=ESMF_STATEINTENT_INTERNAL, &
           SHORT_NAME='QSNOW', &
           STANDARD_NAME='mass_fraction_of_snow', &
           UNITS='kg kg-1', &
           DIMS='xyz', &
           VERTICAL_STAGGER=MAPL_VERTICAL_STAGGER_CENTER, &
           FILL_VALUE=0.0, &
           _RC)

      call MAPL_GridCompAddSpec(GC, &
           STATE_INTENT=ESMF_STATEINTENT_INTERNAL, &
           SHORT_NAME='QGRAUPEL', &
           STANDARD_NAME='mass_fraction_of_graupel', &
           UNITS='kg kg-1', &
           DIMS='xyz', &
           VERTICAL_STAGGER=MAPL_VERTICAL_STAGGER_CENTER, &
           FILL_VALUE=0.0, &
           _RC)

      call MAPL_GridCompAddSpec(GC, &
           STATE_INTENT=ESMF_STATEINTENT_INTERNAL, &
           SHORT_NAME='NACTL', &
           STANDARD_NAME='activ aero # conc liq phase for 1-mom', &
           UNITS='m-3', &
           DIMS='xyz', &
           VERTICAL_STAGGER=MAPL_VERTICAL_STAGGER_CENTER, &
           RESTART_MODE=MAPL_RESTART_SKIP, &
           FILL_VALUE=0.0, &
           _RC)

      call MAPL_GridCompAddSpec(GC, &
           STATE_INTENT=ESMF_STATEINTENT_INTERNAL, &
           SHORT_NAME='NACTI', &
           STANDARD_NAME='activ aero # conc ice phase for 1-mom', &
           UNITS='m-3', &
           DIMS='xyz', &
           VERTICAL_STAGGER=MAPL_VERTICAL_STAGGER_CENTER, &
           RESTART_MODE=MAPL_RESTART_SKIP, &
           FILL_VALUE=0.0, &
           _RC)

      call MAPL_GridCompAddSpec(GC, &
           STATE_INTENT=ESMF_STATEINTENT_EXPORT, &
           SHORT_NAME='REF_DBZ', &
           STANDARD_NAME='Simulated_gfdl_radar_reflectivity', &
           UNITS='dBZ', &
           DIMS='xyz', &
           VERTICAL_STAGGER=MAPL_VERTICAL_STAGGER_CENTER, &
           _RC)

      call MAPL_GridCompAddSpec(GC, &
           STATE_INTENT=ESMF_STATEINTENT_EXPORT, &
           SHORT_NAME='REF_DBZ_MAX', &
           STANDARD_NAME='Maximum_composite_gfdl_radar_reflectivity', &
           UNITS='dBZ', &
           DIMS='xy', &
           VERTICAL_STAGGER=MAPL_VERTICAL_STAGGER_NONE, &
           _RC)

      call MAPL_GridCompAddSpec(GC, &
           STATE_INTENT=ESMF_STATEINTENT_EXPORT, &
           SHORT_NAME='REF_DBZ_1KM', &
           STANDARD_NAME='Base_1KM_AGL_gfdl_radar_reflectivity', &
           UNITS='dBZ', &
           DIMS='xy', &
           VERTICAL_STAGGER=MAPL_VERTICAL_STAGGER_NONE, &
           _RC)

      call MAPL_GridCompAddSpec(GC, &
           STATE_INTENT=ESMF_STATEINTENT_EXPORT, &
           SHORT_NAME='REF_DBZ_TOP', &
           STANDARD_NAME='Echo_top_gfdl_radar_reflectivity', &
           UNITS='m', &
           DIMS='xy', &
           VERTICAL_STAGGER=MAPL_VERTICAL_STAGGER_NONE, &
           _RC)

      call MAPL_GridCompAddSpec(GC, &
           STATE_INTENT=ESMF_STATEINTENT_EXPORT, &
           SHORT_NAME='REF_DBZ_M10C', &
           STANDARD_NAME='Minus_10C_gfdl_radar_reflectivity', &
           UNITS='dBZ', &
           DIMS='xy', &
           VERTICAL_STAGGER=MAPL_VERTICAL_STAGGER_NONE, &
           _RC)

   end subroutine GFDL_1M_Setup

   subroutine GFDL_1M_Initialize(GC, CLOCK, IMPORT, EXPORT, RC)
      type(ESMF_GridComp), intent(inout) :: GC ! gridded component
      type(ESMF_Clock), intent(inout) :: CLOCK ! The clock
      integer, optional :: RC ! return code

      type(ESMF_State) :: INTERNAL
      type(ESMF_State), intent(inout) :: IMPORT
      type(ESMF_State), intent(inout) :: EXPORT

      real(kind=ESMF_KIND_R8) :: DT_R8
      real :: DT_MOIST

      real, pointer, dimension(:, :, :) :: Q, QLLS, QLCN, QILS, QICN, QRAIN, QSNOW, QGRAUPEL

      character(len=ESMF_MAXSTR) :: errmsg

      ! Moist's timestep (MAPL2: the RUNALARM ring interval)
      call MAPL_ClockGet(CLOCK, dt=DT_MOIST, _RC)
      DT_R8 = DT_MOIST

      call MAPL_GridCompGetResource(GC, 'PHYS_HYDROSTATIC', LPHYS_HYDROSTATIC, default=.true., _RC)
      LHYDROSTATIC = LPHYS_HYDROSTATIC
      call MAPL_GridCompGetResource(GC, 'MELTFRZ_CLDMACRO', LMELTFRZ_CLDMACRO, default=.true., _RC)
      call MAPL_GridCompGetResource(GC, 'MELTFRZ_CLDMICRO', LMELTFRZ_CLDMICRO, default=.false., _RC)

      call MAPL_GridCompGetInternalState(GC, INTERNAL, _RC)

      call MAPL_StateGetPointer(INTERNAL, Q, 'Q', _RC)
      call MAPL_StateGetPointer(INTERNAL, QRAIN, 'QRAIN', _RC)
      call MAPL_StateGetPointer(INTERNAL, QSNOW, 'QSNOW', _RC)
      call MAPL_StateGetPointer(INTERNAL, QGRAUPEL, 'QGRAUPEL', _RC)
      call MAPL_StateGetPointer(INTERNAL, QLLS, 'QLLS', _RC)
      call MAPL_StateGetPointer(INTERNAL, QLCN, 'QLCN', _RC)
      call MAPL_StateGetPointer(INTERNAL, QILS, 'QILS', _RC)
      call MAPL_StateGetPointer(INTERNAL, QICN, 'QICN', _RC)

      call MAPL_GridCompGetResource(GC, 'REPORT_GFDL_1M_NEGATIVES', REPORT_GFDL_1M_NEGATIVES, default=.false., _RC)

      call MAPL_GridCompGetResource(GC, 'GFDL_MP3', GFDL_MP3, default=.true., _RC)
      if (DT_R8 <= 150.0) do_hail = .true.
      if (DT_R8 <= 150.0) do_sedi_heat = .true.
      if (DT_R8 <= 150.0) do_sedi_melt_qi = .true.
      if (DT_R8 <= 150.0) do_sedi_melt_qs = .true.
      if (DT_R8 <= 150.0) do_sedi_melt_qg = .true.

      if (GFDL_MP3) then
         call gfdl_mp_init(LHYDROSTATIC, DT_MOIST)
         call WRITE_PARALLEL("INITIALIZED GFDL_1M gfdl_mp v3 in non-generic GC INIT")
         call MAPL_GridCompGetResource(GC, 'DO_GFDL_REFLECTIVITY', do_ref, default=.true., _RC)
      else
         call gfdl_cloud_microphys_init()
         call WRITE_PARALLEL("INITIALIZED GFDL_1M gfdl_cloud_microphys in non-generic GC INIT")
         do_ref = .false. ! Force to false so MAPL DBZ Calc triggers, as older driver has no DBZ3D
      end if

      ! DSL entry point
      call MAPL_GridCompGetResource(GC, 'USE_PYMOIST_GFDL1M', USE_PYMOIST_GFDL1M, default=.false., _RC)

      if (USE_PYMOIST_GFDL1M) then
         call MAPL_ConfigSetAttribute(CF, DT_MOIST, 'DSL__GFLD1M_DT:', _RC)
         call MAPL_pybridge_gcinit("pyMoist.fortran.param_interfaces.microphysics.GFDL1M_interface", MAPL, IMPORT, &
              EXPORT)
      end if

      call MAPL_GridCompGetResource(GC, 'SH_MD_DP', SH_MD_DP, default=.true., _RC)

      ! -----------------------------------------------------------------------------------------
      ! MACROPHYSICS: CRITICAL RELATIVE HUMIDITY (RHCRIT) PROFILES
      ! -----------------------------------------------------------------------------------------
      ! Controls the grid-mean relative humidity threshold required for sub-grid cloud formation
      ! and maintenance. A higher RHCRIT assumes a highly homogeneous grid box (requiring near
      ! 100% RH to form clouds), while a lower RHCRIT assumes high sub-grid moisture variance.
      ! The vertical profile transitions from a high, well-mixed boundary layer plateau to a
      ! lower free-troposphere floor, modulated by the Estimated Inversion Strength (EIS).
      !
      !   - TURNRHCRIT_SFC  : Height [m] of the boundary layer top. If <= 0.0, dynamically tracks
      !                       the prognostic PBL height (KPBLSC). If > 0.0, acts as a static height.
      !   - MAX_RH_CRIT     : Absolute ceiling [fraction] applied near the surface (well-mixed BL).
      !   - MIN_RH_UNSTABLE : Target RHCRIT at PBL top for unstable regimes (e.g., Trade Cumulus).
      !   - MIN_RH_STABLE   : Target RHCRIT at PBL top for highly stable regimes (e.g., Stratocumulus).
      !   - MIN_RH_CRIT     : Absolute floor [fraction] reached in the heterogeneous free troposphere.
      !   - TURNRHCRIT_TOP  : Height [m] at which the profile fully relaxes to the MIN_RH_CRIT floor.
      call MAPL_GridCompGetResource(GC, 'TURNRHCRIT_SFC', TURNRHCRIT_SFC, default=-1., _RC)
      call MAPL_GridCompGetResource(GC, 'MAX_RH_CRIT', MAX_RH_CRIT, default=0.9900, _RC)
      call MAPL_GridCompGetResource(GC, 'MIN_RH_UNSTABLE', MIN_RH_UNSTABLE, default=0.9750, _RC)
      call MAPL_GridCompGetResource(GC, 'MIN_RH_STABLE', MIN_RH_STABLE, default=0.8750, _RC)
      call MAPL_GridCompGetResource(GC, 'MIN_RH_CRIT', MIN_RH_CRIT, default=0.7250, _RC)
      call MAPL_GridCompGetResource(GC, 'TURNRHCRIT_TOP', TURNRHCRIT_TOP, default=3250., _RC)
      call MAPL_GridCompGetResource(GC, 'PDFSHAPE', PDFSHAPE, default=1, _RC)

      ! -----------------------------------------------------------------------------------------
      ! ICE SETTLING & FALL SPEED MULTIPLIERS
      ! -----------------------------------------------------------------------------------------
      ! Controls the prognostic fall speeds of cloud ice, which directly dictates the residence
      ! time of cirrus anvils and large-scale stratiform ice in the upper troposphere.
      !   - ICE_*_VFALL_PARAM : Selects the empirical fall-speed relationship (1 = Standard).
      !   - ANV_ICEFALL       : Linear multiplier on convective anvil ice fall speed.
      !   - LS_ICEFALL        : Linear multiplier on large-scale (grid-scale) ice fall speed.
      if (.not. GFDL_MP3) then ! GFDL_MP3 has internal parameters for this
         call MAPL_GridCompGetResource(GC, 'ICE_LSC_VFALL_PARAM', ICE_LSC_VFALL_PARAM, default=1, _RC)
         call MAPL_GridCompGetResource(GC, 'ICE_CNV_VFALL_PARAM', ICE_CNV_VFALL_PARAM, default=1, _RC)
         call MAPL_GridCompGetResource(GC, 'ANV_ICEFALL', ANV_ICEFALL, default=1.0, _RC)
         call MAPL_GridCompGetResource(GC, 'LS_ICEFALL', LS_ICEFALL, default=1.0, _RC)
      else
         ! use gfdl_mp namelist parameters
      end if

      ! -----------------------------------------------------------------------------------------
      ! CLOUD OPTICS & EFFECTIVE RADII [m]
      ! -----------------------------------------------------------------------------------------
      ! LIQUID RADIUS OPTIONS (LIQ_RADII_PARAM):
      !   1 = Empirical power-law (legacy)
      !   2 = Liu & Daum empirical with dispersion (legacy)
      !   3 = Physical: Uses NNL from aerosol activation (RECOMMENDED)
      !   4+ = Morrison-Gettelman gamma closure (two-moment schemes)
      !
      ! ICE RADIUS OPTIONS (ICE_RADII_PARAM):
      !   1 = Wyser temperature-dependent (legacy)
      !   2 = Sun temperature-dependent with hexagonal geometry (legacy)
      !   3 = Physical: Uses NNI from aerosol activation with tight bounds (RECOMMENDED)
      !   4+ = Morrison-Gettelman power-law (two-moment schemes)
      !
      ! DISPERSION FACTORS: Scale volume mean → effective radius (account for distribution width)
      ! FAC_R*: Tuning multipliers | MIN/MAX_R*: Hard bounds [m]
      ! -----------------------------------------------------------------------------------------
      call MAPL_GridCompGetResource(GC, 'LIQ_RADII_PARAM', LIQ_RADII_PARAM, default=3, _RC)
      call MAPL_GridCompGetResource(GC, 'ICE_RADII_PARAM', ICE_RADII_PARAM, default=3, _RC)
      call MAPL_GridCompGetResource(GC, 'FAC_RI', FAC_RI, default=1.0, _RC)
      call MAPL_GridCompGetResource(GC, 'MIN_RI', MIN_RI, default=5.e-6, _RC)
      call MAPL_GridCompGetResource(GC, 'MAX_RI', MAX_RI, default=150.e-6, _RC)
      call MAPL_GridCompGetResource(GC, 'FAC_RL', FAC_RL, default=1.0, _RC)
      call MAPL_GridCompGetResource(GC, 'MIN_RL', MIN_RL, default=5.0e-6, _RC)
      call MAPL_GridCompGetResource(GC, 'MAX_RL', MAX_RL, default=60.0e-6, _RC)

      ! -----------------------------------------------------------------------------------------
      ! GRID-SCALE EVAPORATION / SUBLIMATION EFFICIENCIES [s^-1]
      ! -----------------------------------------------------------------------------------------
      ! Controls the baseline efficiency (A_EFF) of clear-air entrainment and sub-grid evaporation
      ! when grid-mean relative humidity drops below RHCRIT.
      !   - CCW_EVAP_EFF : Efficiency for Cloud Liquid Water.
      !   - CCI_EVAP_EFF : Efficiency for Cloud Ice.
      CCW_EVAP_EFF = 4.e-3
      call MAPL_GridCompGetResource(GC, 'CCW_EVAP_EFF', CCW_EVAP_EFF, default=(CCW_EVAP_EFF), _RC)
      CCI_EVAP_EFF = 4.e-3
      call MAPL_GridCompGetResource(GC, 'CCI_EVAP_EFF', CCI_EVAP_EFF, default=(CCI_EVAP_EFF), _RC)

      ! -----------------------------------------------------------------------------------------
      ! CONVECTIVE CLOUD FRACTION (CAPE THRESHOLDS) [J/kg]
      ! -----------------------------------------------------------------------------------------
      ! Diagnoses macrophysical convective area fraction derived from Convective Available Potential
      ! Energy (CAPE). Scales linearly from 0.0 at MIN to 1.0 at MAX.
      call MAPL_GridCompGetResource(GC, 'CNV_FRACTION_MIN', CNV_FRACTION_MIN, default=500.0, _RC)
      call MAPL_GridCompGetResource(GC, 'CNV_FRACTION_MAX', CNV_FRACTION_MAX, default=1500.0, _RC)

      ! -----------------------------------------------------------------------------------------
      ! MACROPHYSICS & PHASE SWITCHES
      ! -----------------------------------------------------------------------------------------
      ! ICE_FRACTION_POLYNOMIAL : Selects the temperature-dependent liquid/ice partitioning curve.
      ! USE_AEROSOL_NN          : Activates aerosol coupling for nucleation.
      ! USE_BERGERON            : Activates macrophysical Wegener-Bergeron-Findeisen logic
      !                           (e.g., modifying mixed-phase cloud fraction and phase overlap).
      !                           Note: Microphysical WBF mass transfer is controlled separately.
      call MAPL_GridCompGetResource(GC, 'ICE_FRACTION_POLYNOMIAL', ICE_FRACTION_POLYNOMIAL, default=V12_ICE_POLYNOMIAL,&
           & _RC)
      call MAPL_GridCompGetResource(GC, 'USE_AEROSOL_NN', USE_AEROSOL_NN, default=(USE_AEROSOL_NN), _RC)
      call MAPL_GridCompGetResource(GC, 'USE_BERGERON', USE_BERGERON, default=.false., _RC)

      ! -----------------------------------------------------------------------------------------
      ! AEROSOL-ACTIVATED NUMBER CONCENTRATION BOUNDS
      ! -----------------------------------------------------------------------------------------
      ! Clamp limits (min/max) applied to the activated hydrometeor number concentrations
      ! returned from the aerosol activation scheme (Aer_Activation). All values are in
      ! SI number density [# m^-3]. These bounds act as physical guardrails and as the
      ! effective default when the activation scheme returns values outside a plausible range;
      ! they are exposed as resources so they can be tuned without recompilation.
      !
      ! Radiative relevance: the liquid number NN feeds the droplet effective radius via
      ! Reff ~ (rho*qc / N)^(1/3) (see LDRADIUS4), so an overly high NN_MIN_LIQ forces
      ! droplets too small over clean marine air, brightening low cloud in shortwave.
      ! The ice bounds constrain the DeMott (2010) ice-nucleus number and thus high/anvil
      ! cloud amount and OLR.
      !
      ! NN_MIN_LIQ : Minimum cloud-droplet number concentration [# m^-3].
      ! NN_MAX_LIQ : Maximum cloud-droplet number concentration [# m^-3].
      ! NN_FAC_LIQ : Global scale factor for the CCN
      !
      ! NN_MIN_ICE : Minimum ice-crystal / ice-nucleus number concentration [# m^-3].
      ! NN_MAX_ICE : Maximum ice-crystal / ice-nucleus number concentration [# m^-3].
      ! NN_FAC_ICE : Global scale factor for the CIN
      ! -----------------------------------------------------------------------------------------
      call MAPL_GridCompGetResource(GC, 'NN_MIN_LIQ', NN_MIN_LIQ, default=50.0e6, _RC)
      call MAPL_GridCompGetResource(GC, 'NN_MAX_LIQ', NN_MAX_LIQ, default=1000.0e6, _RC)
      call MAPL_GridCompGetResource(GC, 'NN_FAC_LIQ', NN_FAC_LIQ, default=1.0, _RC)
      call MAPL_GridCompGetResource(GC, 'NN_MIN_ICE', NN_MIN_ICE, default=1.0e0, _RC)
      call MAPL_GridCompGetResource(GC, 'NN_MAX_ICE', NN_MAX_ICE, default=50.0e6, _RC)
      call MAPL_GridCompGetResource(GC, 'NN_FAC_ICE', NN_FAC_ICE, default=1.0, _RC)

      if (USE_AEROSOL_NN) then
         ! NOTE: For now we hard code in .false. for use_wnet as that is only an option with MG and will be handled
         ! there
         call aer_cloud_init(use_wnet=.false.)
         call WRITE_PARALLEL("INITIALIZED aer_cloud_init")
      end if

      call MAPL_GridCompGetResource(GC, 'GFDL_MP_KLID', GFDL_MP_KLID, default=-999.0, _RC)

   end subroutine GFDL_1M_Initialize

   subroutine GFDL_1M_Run(GC, IMPORT, EXPORT, CLOCK, RC)
      type(ESMF_GridComp), intent(inout) :: GC ! Gridded component
      type(ESMF_State), intent(inout) :: IMPORT ! Import state
      type(ESMF_State), intent(inout) :: EXPORT ! Export state
      type(ESMF_Clock), intent(inout) :: CLOCK ! The clock
      integer, optional, intent(out) :: RC ! Error code:

      ! Local derived type aliases

      type(ESMF_Grid) :: GRID
      type(ESMF_State) :: INTERNAL
      real(kind=ESMF_KIND_R8) :: DT_R8
      real :: DT_MOIST
      logical :: alarm_is_ringing

      ! Internals
      real, pointer, dimension(:, :, :) :: Q, QLLS, QLCN, CLLS, CLCN, QILS, QICN, QRAIN, QSNOW, QGRAUPEL
      real, pointer, dimension(:, :, :) :: NACTL, NACTI
      ! Imports
      real, pointer, dimension(:, :, :) :: ZLE, PLE, T, U, V, W, KH
      real, pointer, contiguous, dimension(:, :, :) :: p3d ! scratch for the edge-field bounds remap
      real, pointer, dimension(:, :) :: AREA, PHIS, FRLAND, TS, DTSX, SH, EVAP, KPBLSC
      real, pointer, dimension(:, :, :) :: SL2, SL3, QT2, QT3, W2, W3, SLQT, WQT, WSL, PDF_A
      real, pointer, dimension(:, :, :) :: OMEGA
      ! Local
      real, allocatable, dimension(:, :, :) :: U0, V0
      real, allocatable, dimension(:, :, :) :: PLEmb, ZLE0
      real, allocatable, dimension(:, :, :) :: PLmb, ZL0
      real, allocatable, dimension(:, :, :) :: PL
      real, allocatable, dimension(:, :, :) :: DZ, DZET, DP, MASS, iMASS
      real, allocatable, dimension(:, :, :) :: QSliq, QSice
      real, allocatable, dimension(:, :, :) :: DBZ3D
      real, allocatable, dimension(:, :, :) :: DQVDTmic, DQLDTmic, DQRDTmic, DQIDTmic, &
           DQSDTmic, DQGDTmic, DQADTmic, &
           DUDTmic, DVDTmic, DTDTmic, DWDTmic
      real, allocatable, dimension(:, :, :) :: QHAIL_, QGRAUPEL_
      integer, allocatable, dimension(:, :) :: KLCL
      integer :: KLID
      real, allocatable, dimension(:, :, :) :: TMP3D
      real, allocatable, dimension(:, :) :: TMP2D
      real, allocatable, dimension(:) :: TMP1D
      ! Exports
      real, allocatable, dimension(:, :) :: LONS, LATS
      real, pointer, dimension(:, :, :) :: NACTR
      ! Re-exported internals (were ALLOC'd exports shared with the driver and the other schemes)
      real, pointer, dimension(:, :) :: LS_PRCP, LS_SNR, ICE, FRZR, CNV_FRC, SRF_TYPE
      real, pointer, dimension(:, :, :) :: DQVDT_macro, DQIDT_macro, DQLDT_macro, DQADT_macro, DQRDT_macro, DQSDT_macro&
           , DQGDT_macro
      real, pointer, dimension(:, :, :) :: DUDT_macro, DVDT_macro, DTDT_macro
      real, pointer, dimension(:, :, :) :: DQVDT_micro, DQIDT_micro, DQLDT_micro, DQADT_micro, DQRDT_micro, DQSDT_micro&
           , DQGDT_micro
      real, pointer, dimension(:, :, :) :: DUDT_micro, DVDT_micro, DTDT_micro
      real, pointer, dimension(:, :, :) :: CLDREFFL, CLDREFFI
      real, pointer, dimension(:, :, :) :: PFL_LS, PFL_AN
      real, pointer, dimension(:, :, :) :: PFI_LS, PFI_AN
      real, pointer, dimension(:, :) :: EIS, LTS
      ! Exports computed in local arrays (were ALLOC'd), copied to the export if it is connected
      real, allocatable, dimension(:, :) :: PRCP_WATER, PRCP_RAIN, PRCP_SNOW, PRCP_ICE, PRCP_GRAUPEL
      real, allocatable, dimension(:, :, :) :: RAD_CF, RAD_QV, RAD_QL, RAD_QI, RAD_QR, RAD_QS, RAD_QG
      real, allocatable, dimension(:, :, :) :: EVAPC, SUBLC
      real, allocatable, dimension(:, :, :) :: RHX, REV_LS, RSU_LS
      real, allocatable, dimension(:, :, :) :: VFALL_ICE, VFALL_SNOW, VFALL_GRAUPEL, VFALL_RAIN
      real, allocatable, dimension(:, :, :) :: PFR_LS, PFS_LS, PFG_LS
      real, allocatable, dimension(:, :, :) :: WTHV2, WQL, PDFITERS
      real, allocatable, dimension(:, :, :) :: RHCRIT3D
      ! Their export pointers
      real, pointer, dimension(:, :) :: PRCP_WATER_ptr, PRCP_RAIN_ptr, PRCP_SNOW_ptr, PRCP_ICE_ptr, PRCP_GRAUPEL_ptr
      real, pointer, dimension(:, :, :) :: RAD_CF_ptr, RAD_QV_ptr, RAD_QL_ptr, RAD_QI_ptr, RAD_QR_ptr, RAD_QS_ptr, &
           RAD_QG_ptr
      real, pointer, dimension(:, :, :) :: EVAPC_ptr, SUBLC_ptr
      real, pointer, dimension(:, :, :) :: RHX_ptr, REV_LS_ptr, RSU_LS_ptr
      real, pointer, dimension(:, :, :) :: VFALL_ICE_ptr, VFALL_SNOW_ptr, VFALL_GRAUPEL_ptr, VFALL_RAIN_ptr
      real, pointer, dimension(:, :, :) :: PFR_LS_ptr, PFS_LS_ptr, PFG_LS_ptr
      real, pointer, dimension(:, :, :) :: WTHV2_ptr, WQL_ptr, PDFITERS_ptr
      real, pointer, dimension(:, :, :) :: RHCRIT3D_ptr
      ! Exports
      real, pointer, dimension(:, :, :) :: CNV_PRC3
      real, pointer, dimension(:, :) :: DBZ_MAX, DBZ_1KM, DBZ_TOP, DBZ_M10C
      real, pointer, dimension(:, :, :) :: DBZ
      real, pointer, dimension(:, :, :) :: DQVDT_FILL
      real, pointer, dimension(:, :, :) :: DQLLSDT_FILL
      real, pointer, dimension(:, :, :) :: DQLCNDT_FILL
      real, pointer, dimension(:, :, :) :: DQILSDT_FILL
      real, pointer, dimension(:, :, :) :: DQICNDT_FILL
      real, pointer, dimension(:, :, :) :: DQRDT_FILL
      real, pointer, dimension(:, :, :) :: DQSDT_FILL
      real, pointer, dimension(:, :, :) :: DQGDT_FILL
      real, pointer, dimension(:, :, :) :: PTR3D
      real, pointer, dimension(:, :) :: PTR2D

      ! Local variables
      real :: tmp_val, rand1
      real :: x_norm, transition_thickness, start_height, upper_scale
      real :: ALPHA, RHCRIT
      real :: fraction_hail

      real, allocatable :: facEIS_2d(:, :), minrhcrit_2d(:, :), turnrhcrit_2d(:, :)
      real, allocatable :: min_rh_free_2d(:, :)
      real, allocatable :: qg_col(:), qh_col(:), prs_col(:), dbz_col(:)

      integer :: IM, JM, LM
      integer :: I, J, L
      type(ESMF_VM) :: VMG

      call ESMF_GridCompGet(GC, VM=VMG, _RC)

      call MAPL_GridCompTimerStart(GC, "--GFDL_1M", _RC)

      ! Get parameters from the component
      !-----------------------------------

      call MAPL_GridCompGet(GC, num_levels=LM, GRID=GRID, _RC)
      call MAPL_GridGet(GRID, IM=IM, JM=JM, _RC)
      call MAPL_GridGetCoordinates(GRID, longitudes=LONS, latitudes=LATS, _RC) ! radians
      call MAPL_GridCompGetInternalState(GC, INTERNAL, _RC)

      ! Moist's timestep (MAPL2: the RUNALARM ring interval)
      call MAPL_ClockGet(CLOCK, dt=DT_MOIST, _RC)
      DT_R8 = DT_MOIST

      if (USE_PYMOIST_GFDL1M) then
         call MAPL_pybridge_gcrun_with_internal("pyMoist.fortran.param_interfaces.microphysics.GFDL1M_interface", MAPL,&
              & IMPORT, EXPORT, INTERNAL)
      else

         call MAPL_StateGetPointer(INTERNAL, Q, 'Q', _RC)
         call MAPL_StateGetPointer(INTERNAL, QRAIN, 'QRAIN', _RC)
         call MAPL_StateGetPointer(INTERNAL, QSNOW, 'QSNOW', _RC)
         call MAPL_StateGetPointer(INTERNAL, QGRAUPEL, 'QGRAUPEL', _RC)
         call MAPL_StateGetPointer(INTERNAL, QLLS, 'QLLS', _RC)
         call MAPL_StateGetPointer(INTERNAL, QLCN, 'QLCN', _RC)
         call MAPL_StateGetPointer(INTERNAL, CLCN, 'CLCN', _RC)
         call MAPL_StateGetPointer(INTERNAL, CLLS, 'CLLS', _RC)
         call MAPL_StateGetPointer(INTERNAL, QILS, 'QILS', _RC)
         call MAPL_StateGetPointer(INTERNAL, QICN, 'QICN', _RC)
         call MAPL_StateGetPointer(INTERNAL, NACTL, 'NACTL', _RC)
         call MAPL_StateGetPointer(INTERNAL, NACTI, 'NACTI', _RC)

         ! Import State
         call MAPL_StateGetPointer(IMPORT, AREA, 'AREA', _RC)
         call MAPL_StateGetPointer(IMPORT, ZLE, 'ZLE', _RC)
         call MAPL_StateGetPointer(IMPORT, PLE, 'PLE', _RC)
         call MAPL_StateGetPointer(IMPORT, T, 'T', _RC)
         call MAPL_StateGetPointer(IMPORT, U, 'U', _RC)
         call MAPL_StateGetPointer(IMPORT, V, 'V', _RC)
         call MAPL_StateGetPointer(IMPORT, W, 'W', _RC)
         call MAPL_StateGetPointer(IMPORT, FRLAND, 'FRLAND', _RC)
         call MAPL_StateGetPointer(IMPORT, KH, 'KH', _RC)
         ! Edge fields: MAPL3 pointers are 1:LM+1; MAPL2's were 0:LM, which this code assumes
         p3d => ZLE
         ZLE(1:IM, 1:JM, 0:LM) => p3d
         p3d => PLE
         PLE(1:IM, 1:JM, 0:LM) => p3d
         p3d => KH
         KH(1:IM, 1:JM, 0:LM) => p3d
         call MAPL_StateGetPointer(IMPORT, PDF_A, 'PDF_A', _RC)
         call MAPL_StateGetPointer(IMPORT, W2, 'W2', _RC)
         call MAPL_StateGetPointer(IMPORT, W3, 'W3', _RC)
         call MAPL_StateGetPointer(IMPORT, WQT, 'WQT', _RC)
         call MAPL_StateGetPointer(IMPORT, WSL, 'WSL', _RC)
         call MAPL_StateGetPointer(IMPORT, SL2, 'SL2', _RC)
         call MAPL_StateGetPointer(IMPORT, SL3, 'SL3', _RC)
         call MAPL_StateGetPointer(IMPORT, QT2, 'QT2', _RC)
         call MAPL_StateGetPointer(IMPORT, QT3, 'QT3', _RC)
         call MAPL_StateGetPointer(IMPORT, SLQT, 'SLQT', _RC)
         call MAPL_StateGetPointer(IMPORT, TS, 'TS', _RC)
         call MAPL_StateGetPointer(IMPORT, KPBLSC, 'KPBL_SC', _RC)
         call MAPL_StateGetPointer(IMPORT, SH, 'SH', _RC)
         call MAPL_StateGetPointer(IMPORT, EVAP, 'EVAP', _RC)
         call MAPL_StateGetPointer(IMPORT, OMEGA, 'OMEGA', _RC)
         call MAPL_StateGetPointer(IMPORT, PHIS, 'PHIS', _RC)

         ! Allocatables
         ! Edge variables
         allocate(ZLE0(IM, JM, 0:LM))
         allocate(PLEmb(IM, JM, 0:LM))
         ! Layer variables
         allocate(U0(IM, JM, LM))
         allocate(V0(IM, JM, LM))
         allocate(ZL0(IM, JM, LM))
         allocate(PLmb(IM, JM, LM))
         allocate(PL(IM, JM, LM))
         allocate(DZET(IM, JM, LM))
         allocate(DZ(IM, JM, LM))
         allocate(DP(IM, JM, LM))
         allocate(MASS(IM, JM, LM))
         allocate(iMASS(IM, JM, LM))
         allocate(QSliq(IM, JM, LM))
         allocate(QSice(IM, JM, LM))
         allocate(DBZ3D(IM, JM, LM))
         allocate(TMP3D(IM, JM, LM))
         ! Local tendencies
         allocate(DQVDTmic(IM, JM, LM))
         allocate(DQLDTmic(IM, JM, LM))
         allocate(DQIDTmic(IM, JM, LM))
         allocate(DQRDTmic(IM, JM, LM))
         allocate(DQSDTmic(IM, JM, LM))
         allocate(DQGDTmic(IM, JM, LM))
         allocate(DQADTmic(IM, JM, LM))
         allocate(DUDTmic(IM, JM, LM))
         allocate(DVDTmic(IM, JM, LM))
         allocate(DTDTmic(IM, JM, LM))
         allocate(DWDTmic(IM, JM, LM))
         ! 2D Variables
         allocate(TMP2D(IM, JM))
         ! 1D Variables
         allocate(TMP1D(LM))

         ! Initialize to clear DBZ
         DBZ3D = -30.0

         ! Derived States
         PLEmb = PLE * .01
         PLmb = 0.5 * (PLEmb(:, :, 0:LM - 1) + PLEmb(:, :, 1:LM))
         PL = 100.0 * PLmb
         do L = 0, LM
            ZLE0(:, :, L) = ZLE(:, :, L) - ZLE(:, :, LM) ! Edge Height (m) above the surface
         end do
         ZL0 = 0.5 * (ZLE0(:, :, 0:LM - 1) + ZLE0(:, :, 1:LM)) ! Layer Height (m) above the surface
         DZET = (ZLE0(:, :, 0:LM - 1) - ZLE0(:, :, 1:LM)) ! Layer thickness (m)
         QSliq = GEOS_QsatLQU(T, PL)
         QSice = GEOS_QsatICE(T, PL)
         DP = (PLE(:, :, 1:LM) - PLE(:, :, 0:LM - 1))
         MASS = DP / MAPL_GRAV
         iMASS = 1.0 / MASS
         U0 = U
         V0 = V
         KLCL = FIND_KLCL(T, Q, PLmb, IM, JM, LM)
         if (GFDL_MP_KLID > 0.0) then
            KLID = GFDL_MP_KLID
         else
            KLID = 1
         end if

         ! Export and/or scratch Variable: computed locally, zero-filled as MAPL2's ALLOC did
         allocate(RAD_CF(IM, JM, LM), RAD_QV(IM, JM, LM), RAD_QL(IM, JM, LM), RAD_QI(IM, JM, LM), source=0.0)
         allocate(RAD_QR(IM, JM, LM), RAD_QS(IM, JM, LM), RAD_QG(IM, JM, LM), source=0.0)
         ! Re-exported internals
         call MAPL_StateGetPointer(INTERNAL, CLDREFFL, 'RL', _RC)
         call MAPL_StateGetPointer(INTERNAL, CLDREFFI, 'RI', _RC)
         ! This export MUST have been filled in the GridComp
         call MAPL_StateGetPointer(INTERNAL, CNV_FRC, 'CNV_FRC', _RC)
         call MAPL_StateGetPointer(INTERNAL, SRF_TYPE, 'SRF_TYPE', _RC)
         ! Exports  required below: computed locally
         allocate(EVAPC(IM, JM, LM), SUBLC(IM, JM, LM), source=0.0)
         allocate(PRCP_WATER(IM, JM), PRCP_RAIN(IM, JM), PRCP_SNOW(IM, JM), PRCP_ICE(IM, JM), PRCP_GRAUPEL(IM, JM), &
              source=0.0)
         ! Exports to be filled: re-exported internals
         call MAPL_StateGetPointer(INTERNAL, LS_PRCP, 'LS_PRCP', _RC)
         call MAPL_StateGetPointer(INTERNAL, LS_SNR, 'LS_SNR', _RC)
         call MAPL_StateGetPointer(INTERNAL, ICE, 'ICE', _RC)
         call MAPL_StateGetPointer(INTERNAL, FRZR, 'FRZR', _RC)
         call MAPL_StateGetPointer(INTERNAL, PFL_AN, 'PFL_AN', _RC)
         call MAPL_StateGetPointer(INTERNAL, PFL_LS, 'PFL_LS', _RC)
         call MAPL_StateGetPointer(INTERNAL, PFI_AN, 'PFI_AN', _RC)
         call MAPL_StateGetPointer(INTERNAL, PFI_LS, 'PFI_LS', _RC)
         ! Edge fields: MAPL3 pointers are 1:LM+1; MAPL2's were 0:LM, which this code assumes
         p3d => PFL_AN
         PFL_AN(1:IM, 1:JM, 0:LM) => p3d
         p3d => PFL_LS
         PFL_LS(1:IM, 1:JM, 0:LM) => p3d
         p3d => PFI_AN
         PFI_AN(1:IM, 1:JM, 0:LM) => p3d
         p3d => PFI_LS
         PFI_LS(1:IM, 1:JM, 0:LM) => p3d
         ! Exports to be filled: computed locally
         allocate(RHX(IM, JM, LM), REV_LS(IM, JM, LM), RSU_LS(IM, JM, LM), source=0.0)
         allocate(PFR_LS(IM, JM, 0:LM), PFS_LS(IM, JM, 0:LM), PFG_LS(IM, JM, 0:LM), source=0.0) ! edge
         allocate(WTHV2(IM, JM, LM), WQL(IM, JM, LM), PDFITERS(IM, JM, LM), source=0.0)
         allocate(VFALL_ICE(IM, JM, LM), VFALL_SNOW(IM, JM, LM), VFALL_GRAUPEL(IM, JM, LM), VFALL_RAIN(IM, JM, LM), &
              source=0.0)
         ! Unused Exports (forced to 0.0)
         call MAPL_StateGetPointer(INTERNAL, PTR2D, 'CN_PRCP', _RC)
         PTR2D = 0.0
         call MAPL_StateGetPointer(INTERNAL, PTR2D, 'AN_PRCP', _RC)
         PTR2D = 0.0
         call MAPL_StateGetPointer(INTERNAL, PTR2D, 'SC_PRCP', _RC)
         PTR2D = 0.0
         call MAPL_StateGetPointer(INTERNAL, PTR2D, 'CN_SNR', _RC)
         PTR2D = 0.0
         call MAPL_StateGetPointer(INTERNAL, PTR2D, 'AN_SNR', _RC)
         PTR2D = 0.0
         call MAPL_StateGetPointer(INTERNAL, PTR2D, 'SC_SNR', _RC)
         PTR2D = 0.0
         ! Lower tropospheric stability and estimated inversion strength from MoistGC
         call MAPL_StateGetPointer(INTERNAL, LTS, 'LTS', _RC)
         call MAPL_StateGetPointer(INTERNAL, EIS, 'EIS', _RC)

         if (DEBUG_TQ_ERRORS) then
            do L = 1, LM
               do J = 1, JM
                  do I = 1, IM
                     if ((T(I, J, L) > 333.0) .or. (T(I, J, L) /= T(I, J, L)) .or. &
                          (Q(I, J, L) < 0.0) .or. (Q(I, J, L) /= Q(I, J, L)) .or. &
                          (QLLS(I, J, L) < 0.0) .or. (QLLS(I, J, L) /= QLLS(I, J, L)) .or. &
                          (QLCN(I, J, L) < 0.0) .or. (QLCN(I, J, L) /= QLCN(I, J, L)) .or. &
                          (QILS(I, J, L) < 0.0) .or. (QILS(I, J, L) /= QILS(I, J, L)) .or. &
                          (QICN(I, J, L) < 0.0) .or. (QICN(I, J, L) /= QICN(I, J, L)) .or. &
                          (QRAIN(I, J, L) < 0.0) .or. (QRAIN(I, J, L) /= QRAIN(I, J, L)) .or. &
                          (QSNOW(I, J, L) < 0.0) .or. (QSNOW(I, J, L) /= QSNOW(I, J, L)) .or. &
                          (QGRAUPEL(I, J, L) < 0.0) .or. (QGRAUPEL(I, J, L) /= QGRAUPEL(I, J, L))) then
                        print *, "T or Q  spike detected : ", T(I, J, L)
                        print *, "    On Entry to GFDL   : "
                        print *, "  Latitude       =", LATS(I, J) * 180.0 / MAPL_PI
                        print *, "  Longitude      =", LONS(I, J) * 180.0 / MAPL_PI
                        print *, "  Pressure (mb)  =", PLmb(I, J, L)
                        print *, "                        CLLS=", CLLS(I, J, L), "   CLCN=", CLCN(I, J, L)
                        print *, "    QV=", Q(I, J, L), " QLLS=", QLLS(I, J, L), "   QLCN=", QLCN(I, J, L)
                        print *, "                        QILS=", QILS(I, J, L), "   QICN=", QICN(I, J, L)
                        print *, "    QR=", QRAIN(I, J, L), "   QS=", QSNOW(I, J, L), "     QG=", QGRAUPEL(I, J, L)
                     end if
                  end do
               end do
            end do
         end if

         call MAPL_GridCompTimerStart(GC, "---CLDMACRO", _RC)
         call MAPL_StateGetPointer(INTERNAL, DQVDT_macro, 'DQVDT_macro', _RC)
         call MAPL_StateGetPointer(INTERNAL, DQIDT_macro, 'DQIDT_macro', _RC)
         call MAPL_StateGetPointer(INTERNAL, DQLDT_macro, 'DQLDT_macro', _RC)
         call MAPL_StateGetPointer(INTERNAL, DQADT_macro, 'DQADT_macro', _RC)
         call MAPL_StateGetPointer(INTERNAL, DQRDT_macro, 'DQRDT_macro', _RC)
         call MAPL_StateGetPointer(INTERNAL, DQSDT_macro, 'DQSDT_macro', _RC)
         call MAPL_StateGetPointer(INTERNAL, DQGDT_macro, 'DQGDT_macro', _RC)
         call MAPL_StateGetPointer(INTERNAL, DUDT_macro, 'DUDT_macro', _RC)
         call MAPL_StateGetPointer(INTERNAL, DVDT_macro, 'DVDT_macro', _RC)
         call MAPL_StateGetPointer(INTERNAL, DTDT_macro, 'DTDT_macro', _RC)

         !$OMP parallel do default(none) &
         !$OMP shared(LM, JM, IM, DUDT_macro, DVDT_macro, DTDT_macro, DQVDT_macro, &
         !$OMP        DQLDT_macro, DQIDT_macro, DQADT_macro, DQRDT_macro, DQSDT_macro, &
         !$OMP        DQGDT_macro, U, V, T, Q, QLCN, QLLS, QICN, QILS, CLCN, CLLS, &
         !$OMP        QRAIN, QSNOW, QGRAUPEL, EVAPC, SUBLC, PDFITERS, RHX) &
         !$OMP private(I, J, L)
         do L = 1, LM
            do J = 1, JM
               do I = 1, IM
                  DUDT_macro(I, J, L) = U(I, J, L)
                  DVDT_macro(I, J, L) = V(I, J, L)
                  DTDT_macro(I, J, L) = T(I, J, L)
                  DQVDT_macro(I, J, L) = Q(I, J, L)
                  DQLDT_macro(I, J, L) = QLCN(I, J, L) + QLLS(I, J, L)
                  DQIDT_macro(I, J, L) = QICN(I, J, L) + QILS(I, J, L)
                  DQADT_macro(I, J, L) = CLCN(I, J, L) + CLLS(I, J, L)
                  DQRDT_macro(I, J, L) = QRAIN(I, J, L)
                  DQSDT_macro(I, J, L) = QSNOW(I, J, L)
                  DQGDT_macro(I, J, L) = QGRAUPEL(I, J, L)
                  ! Clear exports
                  EVAPC(I, J, L) = 0.0
                  SUBLC(I, J, L) = 0.0
                  PDFITERS(I, J, L) = 0.0
                  RHX(I, J, L) = 0.0
               end do
            end do
         end do

         ! Include shallow precip condensates if present
         call MAPL_StateGetPointer(INTERNAL, PTR3D, 'SHLW_PRC3', _RC)
         if (associated(PTR3D)) then
            !$OMP parallel do default(none) &
            !$OMP shared(LM, JM, IM, QRAIN, PTR3D, DT_MOIST) &
            !$OMP private(I, J, L)
            do L = 1, LM
               do J = 1, JM
                  do I = 1, IM
                     QRAIN(I, J, L) = QRAIN(I, J, L) + PTR3D(I, J, L) * DT_MOIST
                  end do
               end do
            end do
         end if

         call MAPL_StateGetPointer(INTERNAL, PTR3D, 'SHLW_SNO3', _RC)
         if (associated(PTR3D)) then
            !$OMP parallel do default(none) &
            !$OMP shared(LM, JM, IM, QSNOW, PTR3D, DT_MOIST) &
            !$OMP private(I, J, L)
            do L = 1, LM
               do J = 1, JM
                  do I = 1, IM
                     QSNOW(I, J, L) = QSNOW(I, J, L) + PTR3D(I, J, L) * DT_MOIST
                  end do
               end do
            end do
         end if

         ! -----------------------------------------------------------------
         ! Precalculate 2D arrays (Hoisted out of the L loop)
         ! -----------------------------------------------------------------
         allocate(facEIS_2d(IM, JM), minrhcrit_2d(IM, JM), turnrhcrit_2d(IM, JM))
         allocate(min_rh_free_2d(IM, JM))

         !$OMP parallel do default(none) &
         !$OMP shared(IM, JM, EIS, SRF_TYPE, MIN_RH_UNSTABLE, MIN_RH_STABLE, &
         !$OMP        TURNRHCRIT_SFC, ZL0, KPBLSC, facEIS_2d, minrhcrit_2d, &
         !$OMP        turnrhcrit_2d, min_rh_free_2d, MIN_RH_CRIT) &
         !$OMP private(I, J)
         do J = 1, JM
            do I = 1, IM

               facEIS_2d(I, J) = get_fac_eis(EIS(I, J), SRF_TYPE(I, J))
               minrhcrit_2d(I, J) = MIN_RH_UNSTABLE * (1.0 - facEIS_2d(I, J)) + MIN_RH_STABLE * facEIS_2d(I, J)

               min_rh_free_2d(I, J) = MIN(minrhcrit_2d(I, J), MIN_RH_CRIT + 0.05)

               minrhcrit_2d(I, J) = max(0.7, minrhcrit_2d(I, J))

               if (TURNRHCRIT_SFC <= 0.0) then
                  turnrhcrit_2d(I, J) = ZL0(I, J, NINT(KPBLSC(I, J)))
               else
                  turnrhcrit_2d(I, J) = TURNRHCRIT_SFC
               end if
            end do
         end do

         ! evap/subl/pdf
         allocate(RHCRIT3D(IM, JM, LM), source=0.0) ! local, copied to the RHCRIT export below

         !$OMP parallel do default(none) &
         !$OMP shared(LM, JM, IM, Q, T, QLLS, QILS, CLLS, QLCN, QICN, CLCN, KLID, &
         !$OMP        facEIS_2d, minrhcrit_2d, turnrhcrit_2d, TURNRHCRIT_TOP, min_rh_free_2d, MAX_RH_CRIT, PLmb, &
         !$OMP        AREA, RHCRIT3D, DT_MOIST, PDFSHAPE, CNV_FRC, SRF_TYPE, ZL0, NACTL, &
         !$OMP        NACTI, WSL, WQT, SL2, QT2, SLQT, W3, W2, QT3, SL3, PDF_A, PDFITERS, &
         !$OMP        WTHV2, WQL, USE_BERGERON, RHX, LMELTFRZ_CLDMACRO, CCW_EVAP_EFF, &
         !$OMP        EVAPC, CCI_EVAP_EFF, SUBLC, QSliq, QSice) &
         !$OMP private(L, J, I, transition_thickness, start_height, upper_scale, &
         !$OMP         x_norm, RHCRIT, ALPHA)
         do L = 1, LM
            do J = 1, JM
               do I = 1, IM
                  ! cleanup clouds before cldmacro
                  call FIX_UP_CLOUDS(Q(I, J, L), T(I, J, L), QLLS(I, J, L), QILS(I, J, L), CLLS(I, J, L), &
                       QLCN(I, J, L), QICN(I, J, L), CLCN(I, J, L), &
                       REMOVE_CLOUDS=(L < KLID))

                  ! -----------------------------------------------------------------
                  ! Turn RHCRIT in the vertical column
                  ! -----------------------------------------------------------------
                  transition_thickness = 0.25 * turnrhcrit_2d(I, J)
                  start_height = turnrhcrit_2d(I, J) - transition_thickness

                  if (ZL0(I, J, L) >= turnrhcrit_2d(I, J)) then
                     ! 1. Free Troposphere: Sharp relaxation targeting the RH floor
                     if (ZL0(I, J, L) < TURNRHCRIT_TOP) then
                        upper_scale = (ZL0(I, J, L) - turnrhcrit_2d(I, J)) / (TURNRHCRIT_TOP - turnrhcrit_2d(I, J))
                        RHCRIT = minrhcrit_2d(I, J) - (minrhcrit_2d(I, J) - min_rh_free_2d(I, J)) * upper_scale
                     else
                        ! Lock straight onto your min_rh_free_2d(I,J) default floor at and above TURNRHCRIT_TOP
                        RHCRIT = min_rh_free_2d(I, J)
                     end if

                  else if (ZL0(I, J, L) <= start_height) then
                     ! 2. Hard plateau: Lock at MAX_RH_CRIT near the ground
                     RHCRIT = MAX_RH_CRIT

                  else
                     ! 3. Boundary Layer S-curve transition
                     x_norm = (turnrhcrit_2d(I, J) - ZL0(I, J, L)) / transition_thickness
                     x_norm = max(0.0, MIN(1.0, x_norm))

                     RHCRIT = minrhcrit_2d(I, J) + (MAX_RH_CRIT - minrhcrit_2d(I, J)) * &
                          (x_norm * x_norm * (3.0 - 2.0 * x_norm))
                  end if

                  ! -----------------------------------------------------------------
                  ! Scale-Aware Blending for RHCRIT
                  ! -----------------------------------------------------------------
                  RHCRIT = MAX_RH_CRIT + (RHCRIT - MAX_RH_CRIT) * SQRT(SQRT(AREA(I, J) / 1.e10))
                  ! ===============================================================
                  ALPHA = max(0.0, MIN(0.3, (1.0 - RHCRIT)))
                  ! fill RHCRIT export
                  RHCRIT3D(I, J, L) = 1.0 - ALPHA
                  ! Do CLOUD MACRO below the pressure lid
                  if (L >= KLID) then
                     ! Put condensates in touch with the PDF
                     call hystpdf( &
                          DT_MOIST, &
                          ALPHA, &
                          PDFSHAPE, &
                          CNV_FRC(I, J), &
                          SRF_TYPE(I, J), &
                          PLmb(I, J, L), &
                          ZL0(I, J, L), &
                          Q(I, J, L), &
                          QLLS(I, J, L), &
                          QLCN(I, J, L), &
                          QILS(I, J, L), &
                          QICN(I, J, L), &
                          T(I, J, L), &
                          CLLS(I, J, L), &
                          CLCN(I, J, L), &
                          NACTL(I, J, L), &
                          NACTI(I, J, L), &
                          WSL(I, J, L), &
                          WQT(I, J, L), &
                          SL2(I, J, L), &
                          QT2(I, J, L), &
                          SLQT(I, J, L), &
                          W3(I, J, L), &
                          W2(I, J, L), &
                          QT3(I, J, L), &
                          SL3(I, J, L), &
                          PDF_A(I, J, L), &
                          PDFITERS(I, J, L), &
                          WTHV2(I, J, L), &
                          WQL(I, J, L), &
                          .false., &
                          USE_BERGERON)
                     RHX(I, J, L) = Q(I, J, L) / GEOS_QSAT(T(I, J, L), PLmb(I, J, L))
                     if (LMELTFRZ_CLDMACRO) then
                        ! meltfrz new condensates
                        call MELTFRZ(DT_MOIST, &
                             CNV_FRC(I, J), &
                             SRF_TYPE(I, J), &
                             T(I, J, L), &
                             QLCN(I, J, L), &
                             QICN(I, J, L))
                     end if
                     ! evaporation for CN
                     if (CCW_EVAP_EFF > 0.0) then ! else evap done inside GFDL
                        EVAPC(I, J, L) = Q(I, J, L)
                        call EVAP3( &
                             DT_MOIST, &
                             CCW_EVAP_EFF, &
                             RHCRIT3D(I, J, L), &
                             PLmb(I, J, L), &
                             T(I, J, L), &
                             Q(I, J, L), &
                             QLCN(I, J, L), &
                             QICN(I, J, L), &
                             CLCN(I, J, L), &
                             NACTL(I, J, L), &
                             NACTI(I, J, L), &
                             QSliq(I, J, L), &
                             CNV_FRC(I, J))
                        EVAPC(I, J, L) = (Q(I, J, L) - EVAPC(I, J, L)) / DT_MOIST
                     end if
                     ! sublimation for CN
                     if (CCI_EVAP_EFF > 0.0) then ! else subl done inside GFDL
                        SUBLC(I, J, L) = Q(I, J, L)
                        call SUBL3( &
                             DT_MOIST, &
                             CCI_EVAP_EFF, &
                             RHCRIT3D(I, J, L), &
                             PLmb(I, J, L), &
                             T(I, J, L), &
                             Q(I, J, L), &
                             QLCN(I, J, L), &
                             QICN(I, J, L), &
                             CLCN(I, J, L), &
                             NACTL(I, J, L), &
                             NACTI(I, J, L), &
                             QSice(I, J, L), &
                             CNV_FRC(I, J))
                        SUBLC(I, J, L) = (Q(I, J, L) - SUBLC(I, J, L)) / DT_MOIST
                     end if
                  end if
                  ! cleanup clouds after cldmacro
                  call FIX_UP_CLOUDS(Q(I, J, L), T(I, J, L), QLLS(I, J, L), QILS(I, J, L), CLLS(I, J, L), &
                       QLCN(I, J, L), QICN(I, J, L), CLCN(I, J, L), &
                       REMOVE_CLOUDS=(L < KLID))
               end do ! IM loop
            end do ! JM loop
         end do ! LM loop

         ! Clean up memory allocations
         deallocate(facEIS_2d, minrhcrit_2d, turnrhcrit_2d)

         ! Get fill negative export pointers if requested
         ! ----------------------------------------------
         call MAPL_StateGetPointer(EXPORT, DQVDT_FILL, 'DQVDT_FILL_CLDMACRO', _RC)
         call MAPL_StateGetPointer(EXPORT, DQLLSDT_FILL, 'DQLLSDT_FILL_CLDMACRO', _RC)
         call MAPL_StateGetPointer(EXPORT, DQLCNDT_FILL, 'DQLCNDT_FILL_CLDMACRO', _RC)
         call MAPL_StateGetPointer(EXPORT, DQILSDT_FILL, 'DQILSDT_FILL_CLDMACRO', _RC)
         call MAPL_StateGetPointer(EXPORT, DQICNDT_FILL, 'DQICNDT_FILL_CLDMACRO', _RC)
         call MAPL_StateGetPointer(EXPORT, DQRDT_FILL, 'DQRDT_FILL_CLDMACRO', _RC)
         call MAPL_StateGetPointer(EXPORT, DQSDT_FILL, 'DQSDT_FILL_CLDMACRO', _RC)
         call MAPL_StateGetPointer(EXPORT, DQGDT_FILL, 'DQGDT_FILL_CLDMACRO', _RC)
         ! Cleanup negative water species
         ! ------------------------------
         call FILLQ2ZERO(Q, MASS, dt=DT_MOIST, DQDT=DQVDT_FILL, VM=VMG, _RC)
         call FILLQ2ZERO(QLLS, MASS, dt=DT_MOIST, DQDT=DQLLSDT_FILL, VM=VMG, _RC)
         call FILLQ2ZERO(QLCN, MASS, dt=DT_MOIST, DQDT=DQLCNDT_FILL, VM=VMG, _RC)
         call FILLQ2ZERO(QILS, MASS, dt=DT_MOIST, DQDT=DQILSDT_FILL, VM=VMG, _RC)
         call FILLQ2ZERO(QICN, MASS, dt=DT_MOIST, DQDT=DQICNDT_FILL, VM=VMG, _RC)
         call FILLQ2ZERO(QRAIN, MASS, dt=DT_MOIST, DQDT=DQRDT_FILL, VM=VMG, _RC)
         call FILLQ2ZERO(QSNOW, MASS, dt=DT_MOIST, DQDT=DQSDT_FILL, VM=VMG, _RC)
         call FILLQ2ZERO(QGRAUPEL, MASS, dt=DT_MOIST, DQDT=DQGDT_FILL, VM=VMG, _RC)

         ! Update macrophysics tendencies
         !$OMP parallel do default(none) &
         !$OMP shared(LM, JM, IM, DUDT_macro, DVDT_macro, DTDT_macro, DQVDT_macro, &
         !$OMP        DQLDT_macro, DQIDT_macro, DQADT_macro, DQRDT_macro, DQSDT_macro, &
         !$OMP        DQGDT_macro, U, V, T, Q, QLCN, QLLS, QICN, QILS, CLCN, CLLS, &
         !$OMP        QRAIN, QSNOW, QGRAUPEL, DT_MOIST) &
         !$OMP private(I, J, L)
         do L = 1, LM
            do J = 1, JM
               do I = 1, IM
                  DUDT_macro(I, J, L) = (U(I, J, L) - DUDT_macro(I, J, L)) / DT_MOIST
                  DVDT_macro(I, J, L) = (V(I, J, L) - DVDT_macro(I, J, L)) / DT_MOIST
                  DTDT_macro(I, J, L) = (T(I, J, L) - DTDT_macro(I, J, L)) / DT_MOIST
                  DQVDT_macro(I, J, L) = (Q(I, J, L) - DQVDT_macro(I, J, L)) / DT_MOIST
                  DQLDT_macro(I, J, L) = ((QLCN(I, J, L) + QLLS(I, J, L)) - DQLDT_macro(I, J, L)) / DT_MOIST
                  DQIDT_macro(I, J, L) = ((QICN(I, J, L) + QILS(I, J, L)) - DQIDT_macro(I, J, L)) / DT_MOIST
                  DQADT_macro(I, J, L) = ((CLCN(I, J, L) + CLLS(I, J, L)) - DQADT_macro(I, J, L)) / DT_MOIST
                  DQRDT_macro(I, J, L) = (QRAIN(I, J, L) - DQRDT_macro(I, J, L)) / DT_MOIST
                  DQSDT_macro(I, J, L) = (QSNOW(I, J, L) - DQSDT_macro(I, J, L)) / DT_MOIST
                  DQGDT_macro(I, J, L) = (QGRAUPEL(I, J, L) - DQGDT_macro(I, J, L)) / DT_MOIST
               end do
            end do
         end do
         call MAPL_GridCompTimerStop(GC, "---CLDMACRO", _RC)

         if (DEBUG_TQ_ERRORS) then
            do L = 1, LM
               do J = 1, JM
                  do I = 1, IM
                     if ((T(I, J, L) > 333.0) .or. (T(I, J, L) /= T(I, J, L)) .or. &
                          (Q(I, J, L) < 0.0) .or. (Q(I, J, L) /= Q(I, J, L)) .or. &
                          (QLLS(I, J, L) < 0.0) .or. (QLLS(I, J, L) /= QLLS(I, J, L)) .or. &
                          (QLCN(I, J, L) < 0.0) .or. (QLCN(I, J, L) /= QLCN(I, J, L)) .or. &
                          (QILS(I, J, L) < 0.0) .or. (QILS(I, J, L) /= QILS(I, J, L)) .or. &
                          (QICN(I, J, L) < 0.0) .or. (QICN(I, J, L) /= QICN(I, J, L)) .or. &
                          (QRAIN(I, J, L) < 0.0) .or. (QRAIN(I, J, L) /= QRAIN(I, J, L)) .or. &
                          (QSNOW(I, J, L) < 0.0) .or. (QSNOW(I, J, L) /= QSNOW(I, J, L)) .or. &
                          (QGRAUPEL(I, J, L) < 0.0) .or. (QGRAUPEL(I, J, L) /= QGRAUPEL(I, J, L))) then
                        print *, "T or Q  spike detected : ", T(I, J, L)
                        print *, "    GFDL DTDT_macro Temp Increment : ", DTDT_macro(I, J, L) * DT_MOIST
                        print *, "  Latitude       =", LATS(I, J) * 180.0 / MAPL_PI
                        print *, "  Longitude      =", LONS(I, J) * 180.0 / MAPL_PI
                        print *, "  Pressure (mb)  =", PLmb(I, J, L)
                        print *, "                        CLLS=", CLLS(I, J, L), "   CLCN=", CLCN(I, J, L)
                        print *, "    QV=", Q(I, J, L), " QLLS=", QLLS(I, J, L), "   QLCN=", QLCN(I, J, L)
                        print *, "                        QILS=", QILS(I, J, L), "   QICN=", QICN(I, J, L)
                        print *, "    QR=", QRAIN(I, J, L), "   QS=", QSNOW(I, J, L), "     QG=", QGRAUPEL(I, J, L)
                     end if
                  end do
               end do
            end do
         end if

         call MAPL_GridCompTimerStart(GC, "---CLDMICRO", _RC)
         ! Zero-out microphysics tendencies
         call MAPL_StateGetPointer(INTERNAL, DQVDT_micro, 'DQVDT_micro', _RC)
         call MAPL_StateGetPointer(INTERNAL, DQIDT_micro, 'DQIDT_micro', _RC)
         call MAPL_StateGetPointer(INTERNAL, DQLDT_micro, 'DQLDT_micro', _RC)
         call MAPL_StateGetPointer(INTERNAL, DQADT_micro, 'DQADT_micro', _RC)
         call MAPL_StateGetPointer(INTERNAL, DQRDT_micro, 'DQRDT_micro', _RC)
         call MAPL_StateGetPointer(INTERNAL, DQSDT_micro, 'DQSDT_micro', _RC)
         call MAPL_StateGetPointer(INTERNAL, DQGDT_micro, 'DQGDT_micro', _RC)
         call MAPL_StateGetPointer(INTERNAL, DUDT_micro, 'DUDT_micro', _RC)
         call MAPL_StateGetPointer(INTERNAL, DVDT_micro, 'DVDT_micro', _RC)
         call MAPL_StateGetPointer(INTERNAL, DTDT_micro, 'DTDT_micro', _RC)

         !$OMP parallel do default(none) &
         !$OMP shared(LM, JM, IM, Q, QLLS, QLCN, QILS, QICN, QRAIN, QSNOW, QGRAUPEL, &
         !$OMP        CLLS, CLCN, U, V, T, DZET, &
         !$OMP        DQVDT_micro, DQLDT_micro, DQIDT_micro, DQRDT_micro, DQSDT_micro, &
         !$OMP        DQGDT_micro, DQADT_micro, DUDT_micro, DVDT_micro, DTDT_micro, &
         !$OMP        DZ, DQVDTmic, DQLDTmic, DQRDTmic, DQIDTmic, DQSDTmic, DQGDTmic, &
         !$OMP        DQADTmic, DUDTmic, DVDTmic, DTDTmic, PFI_LS, PFL_LS, &
         !$OMP        RAD_CF, RAD_QL, RAD_QI, RAD_QV, RAD_QR, RAD_QS, RAD_QG) &
         !$OMP private(I, J, L)
         do L = 1, LM
            do J = 1, JM
               do I = 1, IM
                  DQVDT_micro(I, J, L) = Q(I, J, L)
                  DQLDT_micro(I, J, L) = QLLS(I, J, L) + QLCN(I, J, L)
                  DQIDT_micro(I, J, L) = QILS(I, J, L) + QICN(I, J, L)
                  DQRDT_micro(I, J, L) = QRAIN(I, J, L)
                  DQSDT_micro(I, J, L) = QSNOW(I, J, L)
                  DQGDT_micro(I, J, L) = QGRAUPEL(I, J, L)
                  DQADT_micro(I, J, L) = CLLS(I, J, L) + CLCN(I, J, L)
                  DUDT_micro(I, J, L) = U(I, J, L)
                  DVDT_micro(I, J, L) = V(I, J, L)
                  DTDT_micro(I, J, L) = T(I, J, L)

                  ! Delta-Z layer thickness (gfdl expects this to be negative)
                  DZ(I, J, L) = -1.0 * DZET(I, J, L)

                  ! Zero-out local microphysics tendencies
                  DQVDTmic(I, J, L) = 0.0
                  DQLDTmic(I, J, L) = 0.0
                  DQRDTmic(I, J, L) = 0.0
                  DQIDTmic(I, J, L) = 0.0
                  DQSDTmic(I, J, L) = 0.0
                  DQGDTmic(I, J, L) = 0.0
                  DQADTmic(I, J, L) = 0.0
                  DUDTmic(I, J, L) = 0.0
                  DVDTmic(I, J, L) = 0.0
                  DTDTmic(I, J, L) = 0.0

                  ! Zero-out 3D Precipitation Fluxes
                  PFI_LS(I, J, L) = 0.0
                  PFL_LS(I, J, L) = 0.0

                  ! Cloud fractions and condensates for radiation
                  ! Using Maximum-Random Overlap for horizontal combination
                  RAD_CF(I, J, L) = CLCN(I, J, L) + CLLS(I, J, L) - (CLCN(I, J, L) * CLLS(I, J, L))
                  RAD_CF(I, J, L) = MIN(max(RAD_CF(I, J, L), 0.0), 1.0) ! Safegaurd to [0,1]
                  RAD_QL(I, J, L) = QLCN(I, J, L) + QLLS(I, J, L)
                  RAD_QI(I, J, L) = QICN(I, J, L) + QILS(I, J, L)
                  RAD_QV(I, J, L) = Q(I, J, L)
                  RAD_QR(I, J, L) = QRAIN(I, J, L)
                  RAD_QS(I, J, L) = QSNOW(I, J, L)
                  RAD_QG(I, J, L) = QGRAUPEL(I, J, L)
               end do
            end do
         end do

         ! Run the driver
         if (GFDL_MP3) then
            call gfdl_mp_driver( &
            ! Input water/cloud species and liquid+ice NACTL & NACTI (#/m^3)
                 RAD_QV, RAD_QL, RAD_QR, RAD_QI, RAD_QS, RAD_QG, RAD_CF, DBZ3D, NACTL, NACTI, &
            ! Input fields
                 T, W, U, V, DZ, DP, &
            ! Other inputs
                 DT_MOIST, RHCRIT3D, PHIS, CNV_FRC, EIS, AREA, SRF_TYPE, &
            ! Output precipitates
                 PRCP_WATER, PRCP_RAIN, PRCP_ICE, PRCP_SNOW, PRCP_GRAUPEL, &
            ! constant grid/time information
                 LHYDROSTATIC, 1, IM * JM, 1, LM, KLID, &
            ! Output tendencies
                 DQADTmic, &
            ! Output rain re-evaporation and sublimation
                 REV_LS, RSU_LS, &
            ! Output mass flux during sedimentation (Pa kg/kg)
                 PFL_LS(:, :, 1:LM), PFR_LS(:, :, 1:LM), PFI_LS(:, :, 1:LM), PFS_LS(:, :, 1:LM), PFG_LS(:, :, 1:LM))
            ! Convert evap/subl/cloud/precipitation flux exports from (mm/day) to (kg m-2 s-1)
            !$OMP parallel do default(none) &
            !$OMP shared(LM, JM, IM, REV_LS, RSU_LS, PFL_LS, PFI_LS, PFR_LS, PFS_LS, PFG_LS) &
            !$OMP private(I, J, L)
            do L = 1, LM
               do J = 1, JM
                  do I = 1, IM
                     REV_LS(I, J, L) = REV_LS(I, J, L) * (1.0 / 86400.0)
                     RSU_LS(I, J, L) = RSU_LS(I, J, L) * (1.0 / 86400.0)
                     PFL_LS(I, J, L) = PFL_LS(I, J, L) * (1.0 / 86400.0)
                     PFI_LS(I, J, L) = PFI_LS(I, J, L) * (1.0 / 86400.0)
                     PFR_LS(I, J, L) = PFR_LS(I, J, L) * (1.0 / 86400.0)
                     PFS_LS(I, J, L) = PFS_LS(I, J, L) * (1.0 / 86400.0)
                     PFG_LS(I, J, L) = PFG_LS(I, J, L) * (1.0 / 86400.0)
                  end do
               end do
            end do
            if (do_ref) then
               call MAPL_GridCompTimerStart(GC, "---CLD_REF_DBZ", _RC)
               call MAPL_StateGetPointer(EXPORT, DBZ, 'REF_DBZ', _RC)
               call MAPL_StateGetPointer(EXPORT, DBZ_MAX, 'REF_DBZ_MAX', _RC)
               call MAPL_StateGetPointer(EXPORT, DBZ_1KM, 'REF_DBZ_1KM', _RC)
               call MAPL_StateGetPointer(EXPORT, DBZ_TOP, 'REF_DBZ_TOP', _RC)
               call MAPL_StateGetPointer(EXPORT, DBZ_M10C, 'REF_DBZ_M10C', _RC)
               if (associated(DBZ)) DBZ = DBZ3D
               if (associated(DBZ_MAX)) then
                  DBZ_MAX = -9999.0
                  do L = 1, LM
                     do J = 1, JM
                        do I = 1, IM
                           DBZ_MAX(I, J) = max(DBZ_MAX(I, J), DBZ3D(I, J, L))
                        end do
                     end do
                  end do
               end if
               if (associated(DBZ_1KM)) then
                  call cs_interpolator(1, IM, 1, JM, LM, DBZ3D, 1000., ZLE0, DBZ_1KM, -20.)
               end if
               if (associated(DBZ_TOP)) then
                  DBZ_TOP = MAPL_UNDEF
                  do J = 1, JM
                     do I = 1, IM
                        do L = LM, 1, -1
                           if (ZLE0(I, J, L) >= 25000.) continue
                           if (DBZ3D(I, J, L) >= 18.5) then
                              DBZ_TOP(I, J) = ZLE0(I, J, L)
                              exit
                           end if
                        end do
                     end do
                  end do
               end if
               if (associated(DBZ_M10C)) then
                  DBZ_M10C = MAPL_UNDEF
                  do J = 1, JM
                     do I = 1, IM
                        do L = LM, 1, -1
                           if (ZLE0(I, J, L) >= 25000.) continue
                           if (T(I, J, L) <= MAPL_TICE - 10.0) then
                              DBZ_M10C(I, J) = DBZ3D(I, J, L)
                              exit
                           end if
                        end do
                     end do
                  end do
               end if
               call MAPL_GridCompTimerStop(GC, "---CLD_REF_DBZ", _RC)
            end if
         else
            call gfdl_cloud_microphys_driver( &
            ! Input water/cloud species and liquid+ice NACTL & NACTI (#/m^3)
                 RAD_QV, RAD_QL, RAD_QR, RAD_QI, RAD_QS, RAD_QG, RAD_CF, NACTL, NACTI, &
            ! Output tendencies
                 DQVDTmic, DQLDTmic, DQRDTmic, DQIDTmic, &
                 DQSDTmic, DQGDTmic, DQADTmic, DTDTmic, &
            ! Input fields
                 T, W, U, V, DUDTmic, DVDTmic, DZ, DP, &
            ! constant inputs
                 AREA, DT_MOIST, FRLAND, CNV_FRC, SRF_TYPE, EIS, &
                 RHCRIT3D, ANV_ICEFALL, LS_ICEFALL, &
            ! Output rain re-evaporation and sublimation
                 REV_LS, RSU_LS, &
            ! Output fall speeds
                 VFALL_ICE, VFALL_SNOW, VFALL_GRAUPEL, VFALL_RAIN, &
            ! Output precipitates
                 PRCP_RAIN, PRCP_SNOW, PRCP_ICE, PRCP_GRAUPEL, &
            ! Output mass flux during sedimentation (Pa kg/kg)
                 PFL_LS(:, :, 1:LM), PFI_LS(:, :, 1:LM), &
            ! constant grid/time information
                 LHYDROSTATIC, LPHYS_HYDROSTATIC, &
                 1, IM, 1, JM, 1, LM, KLID, LM)
            ! Convert precipitation fluxes from (Pa kg/kg) to (kg m-2 s-1)
            PFL_LS = PFL_LS / (MAPL_GRAV * DT_MOIST)
            PFI_LS = PFI_LS / (MAPL_GRAV * DT_MOIST)
            PFR_LS = 0.0
            PFS_LS = 0.0
            PFG_LS = 0.0
            T = T + DTDTmic * DT_MOIST
            U = U + DUDTmic * DT_MOIST
            V = V + DVDTmic * DT_MOIST
            ! Apply moist/cloud species tendencies
            RAD_QV = RAD_QV + DQVDTmic * DT_MOIST
            RAD_QL = RAD_QL + DQLDTmic * DT_MOIST
            RAD_QR = RAD_QR + DQRDTmic * DT_MOIST
            RAD_QI = RAD_QI + DQIDTmic * DT_MOIST
            RAD_QS = RAD_QS + DQSDTmic * DT_MOIST
            RAD_QG = RAD_QG + DQGDTmic * DT_MOIST
         end if
         ! CleanUp Negative Water Species
         if (REPORT_GFDL_1M_NEGATIVES) then
            call neg_adj_external(IM, JM, LM, RAD_QV, RAD_QL, RAD_QR, RAD_QI, RAD_QS, RAD_QG, &
                 WARNING_LABEL="After GFDL_1M Driver", VM=VMG, _RC)
         else
            call neg_adj_external(IM, JM, LM, RAD_QV, RAD_QL, RAD_QR, RAD_QI, RAD_QS, RAD_QG, &
                 VM=VMG, _RC)
         end if

         ! Update cloud fraction, redistribute clouds, and fill vapor/precip states
         !$OMP parallel do default(none) &
         !$OMP shared(LM, JM, IM, RAD_QL, RAD_QI, RAD_CF, DQADTmic, DT_MOIST, &
         !$OMP        CLCN, CLLS, QLCN, QLLS, QICN, QILS, RAD_QV, T, &
         !$OMP        Q, QRAIN, RAD_QR, QSNOW, RAD_QS, QGRAUPEL, RAD_QG) &
         !$OMP private(I, J, L)
         do L = 1, LM
            do J = 1, JM
               do I = 1, IM
                  ! 1. Update cloud fraction (Replaces where/elsewhere)
                  if (RAD_QL(I, J, L) + RAD_QI(I, J, L) > 0.0) then
                     RAD_CF(I, J, L) = MIN(1.0, max(0.0, RAD_CF(I, J, L) + DQADTmic(I, J, L) * DT_MOIST))
                  else
                     RAD_CF(I, J, L) = 0.0
                  end if

                  ! 2. Redistribute clouds (Calls the new scalar subroutine)
                  call REDISTRIBUTE_CLOUDS_SCALAR(RAD_CF(I, J, L), RAD_QL(I, J, L), RAD_QI(I, J, L), &
                       CLCN(I, J, L), CLLS(I, J, L), QLCN(I, J, L), QLLS(I, J, L), &
                       QICN(I, J, L), QILS(I, J, L), RAD_QV(I, J, L), T(I, J, L))

                  ! 3. Fill vapor/rain/snow/graupel state (Replaces whole-array copies)
                  Q(I, J, L) = RAD_QV(I, J, L)
                  QRAIN(I, J, L) = RAD_QR(I, J, L)
                  QSNOW(I, J, L) = RAD_QS(I, J, L)
                  QGRAUPEL(I, J, L) = RAD_QG(I, J, L)
               end do
            end do
         end do

         ! Get fill negative export pointers if requested
         ! ----------------------------------------------
         call MAPL_StateGetPointer(EXPORT, DQVDT_FILL, 'DQVDT_FILL_CLDMICRO', _RC)
         call MAPL_StateGetPointer(EXPORT, DQLLSDT_FILL, 'DQLLSDT_FILL_CLDMICRO', _RC)
         call MAPL_StateGetPointer(EXPORT, DQLCNDT_FILL, 'DQLCNDT_FILL_CLDMICRO', _RC)
         call MAPL_StateGetPointer(EXPORT, DQILSDT_FILL, 'DQILSDT_FILL_CLDMICRO', _RC)
         call MAPL_StateGetPointer(EXPORT, DQICNDT_FILL, 'DQICNDT_FILL_CLDMICRO', _RC)
         call MAPL_StateGetPointer(EXPORT, DQRDT_FILL, 'DQRDT_FILL_CLDMICRO', _RC)
         call MAPL_StateGetPointer(EXPORT, DQSDT_FILL, 'DQSDT_FILL_CLDMICRO', _RC)
         call MAPL_StateGetPointer(EXPORT, DQGDT_FILL, 'DQGDT_FILL_CLDMICRO', _RC)
         ! Cleanup negative water species
         ! ------------------------------
         if (REPORT_GFDL_1M_NEGATIVES) then
            call FILLQ2ZERO(Q, MASS, dt=DT_MOIST, DQDT=DQVDT_FILL, WARNING_LABEL="QV   After GFDL_1M" // &
                 " REDISTRIBUTE_CLOUDS", VM=VMG, _RC)
            call FILLQ2ZERO(QLLS, MASS, dt=DT_MOIST, DQDT=DQLLSDT_FILL, WARNING_LABEL="QLLS After GFDL_1M" // &
                 " REDISTRIBUTE_CLOUDS", VM=VMG, _RC)
            call FILLQ2ZERO(QLCN, MASS, dt=DT_MOIST, DQDT=DQLCNDT_FILL, WARNING_LABEL="QLCN After GFDL_1M" // &
                 " REDISTRIBUTE_CLOUDS", VM=VMG, _RC)
            call FILLQ2ZERO(QILS, MASS, dt=DT_MOIST, DQDT=DQILSDT_FILL, WARNING_LABEL="QILS After GFDL_1M" // &
                 " REDISTRIBUTE_CLOUDS", VM=VMG, _RC)
            call FILLQ2ZERO(QICN, MASS, dt=DT_MOIST, DQDT=DQICNDT_FILL, WARNING_LABEL="QICN After GFDL_1M" // &
                 " REDISTRIBUTE_CLOUDS", VM=VMG, _RC)
            call FILLQ2ZERO(QRAIN, MASS, dt=DT_MOIST, DQDT=DQRDT_FILL, WARNING_LABEL="QR   After GFDL_1M" // &
                 " REDISTRIBUTE_CLOUDS", VM=VMG, _RC)
            call FILLQ2ZERO(QSNOW, MASS, dt=DT_MOIST, DQDT=DQSDT_FILL, WARNING_LABEL="QS   After GFDL_1M" // &
                 " REDISTRIBUTE_CLOUDS", VM=VMG, _RC)
            call FILLQ2ZERO(QGRAUPEL, MASS, dt=DT_MOIST, DQDT=DQGDT_FILL, WARNING_LABEL="QG   After GFDL_1M" // &
                 " REDISTRIBUTE_CLOUDS", VM=VMG, _RC)
         else
            call FILLQ2ZERO(Q, MASS, dt=DT_MOIST, DQDT=DQVDT_FILL, VM=VMG, _RC)
            call FILLQ2ZERO(QLLS, MASS, dt=DT_MOIST, DQDT=DQLLSDT_FILL, VM=VMG, _RC)
            call FILLQ2ZERO(QLCN, MASS, dt=DT_MOIST, DQDT=DQLCNDT_FILL, VM=VMG, _RC)
            call FILLQ2ZERO(QILS, MASS, dt=DT_MOIST, DQDT=DQILSDT_FILL, VM=VMG, _RC)
            call FILLQ2ZERO(QICN, MASS, dt=DT_MOIST, DQDT=DQICNDT_FILL, VM=VMG, _RC)
            call FILLQ2ZERO(QRAIN, MASS, dt=DT_MOIST, DQDT=DQRDT_FILL, VM=VMG, _RC)
            call FILLQ2ZERO(QSNOW, MASS, dt=DT_MOIST, DQDT=DQSDT_FILL, VM=VMG, _RC)
            call FILLQ2ZERO(QGRAUPEL, MASS, dt=DT_MOIST, DQDT=DQGDT_FILL, VM=VMG, _RC)
         end if

         ! -----------------------------------------------------------------
         ! 1. 2D Precip Diagnostics (Threaded + Div to Mult)
         ! -----------------------------------------------------------------
         !$OMP parallel do default(none) &
         !$OMP shared(IM, JM, PRCP_WATER, PRCP_RAIN, PRCP_SNOW, PRCP_ICE, PRCP_GRAUPEL, &
         !$OMP        LS_PRCP, LS_SNR, ICE, FRZR) &
         !$OMP private(I, J)
         do J = 1, JM
            do I = 1, IM
               PRCP_WATER(I, J) = max(PRCP_WATER(I, J) * (1.0 / 86400.0), 0.0)
               PRCP_RAIN(I, J) = max(PRCP_RAIN(I, J) * (1.0 / 86400.0), 0.0)
               PRCP_SNOW(I, J) = max(PRCP_SNOW(I, J) * (1.0 / 86400.0), 0.0)
               PRCP_ICE(I, J) = max(PRCP_ICE(I, J) * (1.0 / 86400.0), 0.0)
               PRCP_GRAUPEL(I, J) = max(PRCP_GRAUPEL(I, J) * (1.0 / 86400.0), 0.0)

               LS_PRCP(I, J) = PRCP_RAIN(I, J)
               LS_SNR(I, J) = PRCP_SNOW(I, J)
               ICE(I, J) = PRCP_ICE(I, J) + PRCP_GRAUPEL(I, J)
               FRZR(I, J) = 0.0
            end do
         end do

         ! -----------------------------------------------------------------
         ! 2. Fused Flux Redistribution, MeltFreeze, and RadCouple
         ! -----------------------------------------------------------------
         !$OMP parallel do default(none) &
         !$OMP shared(LM, JM, IM, QLCN, RAD_QL, PFL_AN, PFL_LS, PFR_LS, QICN, RAD_QI, &
         !$OMP        PFI_AN, PFI_LS, PFS_LS, PFG_LS, LMELTFRZ_CLDMICRO, DT_MOIST, &
         !$OMP        CNV_FRC, SRF_TYPE, T, QLLS, QILS, Q, CLLS, CLCN, KLID, PLmb, &
         !$OMP        QRAIN, QSNOW, QGRAUPEL, NACTL, NACTI, RAD_QV, RAD_QR, RAD_QS, &
         !$OMP        RAD_QG, RAD_CF, CLDREFFL, CLDREFFI, FAC_RL, MIN_RL, MAX_RL, &
         !$OMP        FAC_RI, MIN_RI, MAX_RI, AREA) &
         !$OMP private(I, J, L, tmp_val)
         do L = 1, LM
            do J = 1, JM
               do I = 1, IM
                  ! Redistribute precipitation fluxes for chemistry (TMP3D eliminated)
                  tmp_val = MIN(1.0, max(QLCN(I, J, L) / max(RAD_QL(I, J, L), 1.E-8), 0.0))
                  PFL_AN(I, J, L) = (PFL_LS(I, J, L) + PFR_LS(I, J, L)) * tmp_val
                  PFL_LS(I, J, L) = (PFL_LS(I, J, L) + PFR_LS(I, J, L)) - PFL_AN(I, J, L)

                  tmp_val = MIN(1.0, max(QICN(I, J, L) / max(RAD_QI(I, J, L), 1.E-8), 0.0))
                  PFI_AN(I, J, L) = (PFI_LS(I, J, L) + PFS_LS(I, J, L) + PFG_LS(I, J, L)) * tmp_val
                  PFI_LS(I, J, L) = (PFI_LS(I, J, L) + PFS_LS(I, J, L) + PFG_LS(I, J, L)) - PFI_AN(I, J, L)

                  ! MeltFreeze and FixUp
                  if (LMELTFRZ_CLDMICRO) then
                     call MELTFRZ(DT_MOIST, CNV_FRC(I, J), SRF_TYPE(I, J), T(I, J, L), QLCN(I, J, L), QICN(I, J, L))
                     call MELTFRZ(DT_MOIST, CNV_FRC(I, J), SRF_TYPE(I, J), T(I, J, L), QLLS(I, J, L), QILS(I, J, L))
                     call FIX_UP_CLOUDS(Q(I, J, L), T(I, J, L), QLLS(I, J, L), QILS(I, J, L), CLLS(I, J, L), &
                          QLCN(I, J, L), QICN(I, J, L), CLCN(I, J, L), REMOVE_CLOUDS=(L < KLID))
                  end if

                  call RADCOUPLE_SCALE_AWARE(T(I, J, L), PLmb(I, J, L), CLLS(I, J, L), CLCN(I, J, L), &
                       Q(I, J, L), QLLS(I, J, L), QILS(I, J, L), QLCN(I, J, L), QICN(I, J, L), &
                       QRAIN(I, J, L), QSNOW(I, J, L), QGRAUPEL(I, J, L), NACTL(I, J, L), NACTI(I, J, L), CNV_FRC(I, J)&
                       , &
                       RAD_QV(I, J, L), RAD_QL(I, J, L), RAD_QI(I, J, L), RAD_QR(I, J, L), RAD_QS(I, J, L), &
                       RAD_QG(I, J, L), RAD_CF(I, J, L), CLDREFFL(I, J, L), CLDREFFI(I, J, L), &
                       FAC_RL, MIN_RL, MAX_RL, FAC_RI, MIN_RI, MAX_RI)
               end do
            end do
         end do

         ! -----------------------------------------------------------------
         ! 3. Cleanup negative water species (Sequential, DO NOT THREAD)
         ! -----------------------------------------------------------------
         if (REPORT_GFDL_1M_NEGATIVES) then
            call FILLQ2ZERO(Q, MASS, WARNING_LABEL="QV   After GFDL_1M RADCOUPLE", VM=VMG, _RC)
            call FILLQ2ZERO(QLLS, MASS, WARNING_LABEL="QLLS After GFDL_1M RADCOUPLE", VM=VMG, _RC)
            call FILLQ2ZERO(QLCN, MASS, WARNING_LABEL="QLCN After GFDL_1M RADCOUPLE", VM=VMG, _RC)
            call FILLQ2ZERO(QILS, MASS, WARNING_LABEL="QILS After GFDL_1M RADCOUPLE", VM=VMG, _RC)
            call FILLQ2ZERO(QICN, MASS, WARNING_LABEL="QICN After GFDL_1M RADCOUPLE", VM=VMG, _RC)
            call FILLQ2ZERO(QRAIN, MASS, WARNING_LABEL="QR   After GFDL_1M RADCOUPLE", VM=VMG, _RC)
            call FILLQ2ZERO(QSNOW, MASS, WARNING_LABEL="QS   After GFDL_1M RADCOUPLE", VM=VMG, _RC)
            call FILLQ2ZERO(QGRAUPEL, MASS, WARNING_LABEL="QG   After GFDL_1M RADCOUPLE", VM=VMG, _RC)
         else
            call FILLQ2ZERO(Q, MASS, VM=VMG, _RC)
            call FILLQ2ZERO(QLLS, MASS, VM=VMG, _RC)
            call FILLQ2ZERO(QLCN, MASS, VM=VMG, _RC)
            call FILLQ2ZERO(QILS, MASS, VM=VMG, _RC)
            call FILLQ2ZERO(QICN, MASS, VM=VMG, _RC)
            call FILLQ2ZERO(QRAIN, MASS, VM=VMG, _RC)
            call FILLQ2ZERO(QSNOW, MASS, VM=VMG, _RC)
            call FILLQ2ZERO(QGRAUPEL, MASS, VM=VMG, _RC)
         end if

         ! -----------------------------------------------------------------
         ! 4. Update microphysics tendencies (Threaded array assignments)
         ! -----------------------------------------------------------------
         !$OMP parallel do default(none) &
         !$OMP shared(LM, JM, IM, DQVDT_micro, DQLDT_micro, DQIDT_micro, DQADT_micro, &
         !$OMP        DQRDT_micro, DQSDT_micro, DQGDT_micro, DUDT_micro, DVDT_micro, &
         !$OMP        DTDT_micro, Q, QLLS, QLCN, QILS, QICN, CLLS, CLCN, QRAIN, &
         !$OMP        QSNOW, QGRAUPEL, U, V, T, DT_MOIST) &
         !$OMP private(I, J, L)
         do L = 1, LM
            do J = 1, JM
               do I = 1, IM
                  DQVDT_micro(I, J, L) = (Q(I, J, L) - DQVDT_micro(I, J, L)) / DT_MOIST
                  DQLDT_micro(I, J, L) = ((QLLS(I, J, L) + QLCN(I, J, L)) - DQLDT_micro(I, J, L)) / DT_MOIST
                  DQIDT_micro(I, J, L) = ((QILS(I, J, L) + QICN(I, J, L)) - DQIDT_micro(I, J, L)) / DT_MOIST
                  DQADT_micro(I, J, L) = ((CLLS(I, J, L) + CLCN(I, J, L)) - DQADT_micro(I, J, L)) / DT_MOIST
                  DQRDT_micro(I, J, L) = (QRAIN(I, J, L) - DQRDT_micro(I, J, L)) / DT_MOIST
                  DQSDT_micro(I, J, L) = (QSNOW(I, J, L) - DQSDT_micro(I, J, L)) / DT_MOIST
                  DQGDT_micro(I, J, L) = (QGRAUPEL(I, J, L) - DQGDT_micro(I, J, L)) / DT_MOIST
                  DUDT_micro(I, J, L) = (U(I, J, L) - DUDT_micro(I, J, L)) / DT_MOIST
                  DVDT_micro(I, J, L) = (V(I, J, L) - DVDT_micro(I, J, L)) / DT_MOIST
                  DTDT_micro(I, J, L) = (T(I, J, L) - DTDT_micro(I, J, L)) / DT_MOIST
               end do
            end do
         end do
         call MAPL_GridCompTimerStop(GC, "---CLDMICRO", _RC)

         if (DEBUG_TQ_ERRORS) then
            do L = 1, LM
               do J = 1, JM
                  do I = 1, IM
                     if ((T(I, J, L) > 333.0) .or. (T(I, J, L) /= T(I, J, L)) .or. &
                          (Q(I, J, L) < 0.0) .or. (Q(I, J, L) /= Q(I, J, L)) .or. &
                          (QLLS(I, J, L) < 0.0) .or. (QLLS(I, J, L) /= QLLS(I, J, L)) .or. &
                          (QLCN(I, J, L) < 0.0) .or. (QLCN(I, J, L) /= QLCN(I, J, L)) .or. &
                          (QILS(I, J, L) < 0.0) .or. (QILS(I, J, L) /= QILS(I, J, L)) .or. &
                          (QICN(I, J, L) < 0.0) .or. (QICN(I, J, L) /= QICN(I, J, L)) .or. &
                          (QRAIN(I, J, L) < 0.0) .or. (QRAIN(I, J, L) /= QRAIN(I, J, L)) .or. &
                          (QSNOW(I, J, L) < 0.0) .or. (QSNOW(I, J, L) /= QSNOW(I, J, L)) .or. &
                          (QGRAUPEL(I, J, L) < 0.0) .or. (QGRAUPEL(I, J, L) /= QGRAUPEL(I, J, L))) then
                        print *, "T or Q  spike detected : ", T(I, J, L)
                        print *, "    GFDL DTDT_micro Temp Increment : ", DTDT_micro(I, J, L) * DT_MOIST
                        print *, "  Latitude       =", LATS(I, J) * 180.0 / MAPL_PI
                        print *, "  Longitude      =", LONS(I, J) * 180.0 / MAPL_PI
                        print *, "  Pressure (mb)  =", PLmb(I, J, L)
                        print *, "                        CLLS=", CLLS(I, J, L), "   CLCN=", CLCN(I, J, L)
                        print *, "                      RAD_QL=", RAD_QL(I, J, L), " RAD_QI=", RAD_QI(I, J, L)
                        print *, "    QV=", Q(I, J, L), " QLLS=", QLLS(I, J, L), "   QLCN=", QLCN(I, J, L)
                        print *, "                        QILS=", QILS(I, J, L), "   QICN=", QICN(I, J, L)
                        print *, "    QR=", QRAIN(I, J, L), "   QS=", QSNOW(I, J, L), "     QG=", QGRAUPEL(I, J, L)
                     end if
                  end do
               end do
            end do
         end if

         call MAPL_GridCompTimerStart(GC, "---CLDDIAGS", _RC)

         call MAPL_StateGetPointer(EXPORT, PTR3D, 'DQRL', _RC)
         if (associated(PTR3D)) PTR3D = DQRDT_macro + DQRDT_micro

         ! dissipative heating tendency from KE across the macro/micro physics
         call MAPL_StateGetPointer(EXPORT, PTR3D, 'DTDTFRIC', _RC)
         if (associated(PTR3D)) then
            call dissipative_ke_heating(IM, JM, LM, MASS, U0, V0, &
                 DUDT_macro + DUDT_micro,&
                 DVDT_macro + DVDT_micro, PTR3D)
         end if

         ! Call the shared radar diagnostics routine
         call MAPL_GridCompTimerStart(GC, "---RADAR_DIAGS", _RC)
         call compute_radar_diagnostics(EXPORT, CLOCK, IM, JM, LM, &
              Q, QRAIN, QSNOW, QGRAUPEL, T, PLmb, W, ZLE0, &
              STATUS)
         _VERIFY(STATUS)
         call MAPL_GridCompTimerStop(GC, "---RADAR_DIAGS", _RC)

         call MAPL_StateGetPointer(EXPORT, PTR3D, 'QRTOT', _RC)
         if (associated(PTR3D)) PTR3D = QRAIN

         call MAPL_StateGetPointer(EXPORT, PTR3D, 'QSTOT', _RC)
         if (associated(PTR3D)) PTR3D = QSNOW

         call MAPL_StateGetPointer(EXPORT, PTR3D, 'QGTOT', _RC)
         if (associated(PTR3D)) PTR3D = QGRAUPEL

         call MAPL_StateGetPointer(EXPORT, PTR2D, 'LWP', _RC)
         if (associated(PTR2D)) PTR2D = SUM((QLCN + QLLS + QRAIN) * MASS, 3)

         call MAPL_StateGetPointer(EXPORT, PTR2D, 'IWP', _RC)
         if (associated(PTR2D)) PTR2D = SUM((QICN + QILS + QSNOW + QGRAUPEL) * MASS, 3)

         call MAPL_StateGetPointer(EXPORT, PTR3D, 'NCPL_VOL', _RC)
         if (associated(PTR3D)) PTR3D = NACTL

         call MAPL_StateGetPointer(EXPORT, PTR3D, 'NCPI_VOL', _RC)
         if (associated(PTR3D)) PTR3D = NACTI

         ! Copy the local arrays to their exports, if connected
         call MAPL_StateGetPointer(EXPORT, RAD_CF_ptr, 'FCLD', _RC)
         call MAPL_StateGetPointer(EXPORT, RAD_QV_ptr, 'QV', _RC)
         call MAPL_StateGetPointer(EXPORT, RAD_QL_ptr, 'QL', _RC)
         call MAPL_StateGetPointer(EXPORT, RAD_QI_ptr, 'QI', _RC)
         call MAPL_StateGetPointer(EXPORT, RAD_QR_ptr, 'QR', _RC)
         call MAPL_StateGetPointer(EXPORT, RAD_QS_ptr, 'QS', _RC)
         call MAPL_StateGetPointer(EXPORT, RAD_QG_ptr, 'QG', _RC)
         call MAPL_StateGetPointer(EXPORT, EVAPC_ptr, 'EVAPC', _RC)
         call MAPL_StateGetPointer(EXPORT, SUBLC_ptr, 'SUBLC', _RC)
         call MAPL_StateGetPointer(EXPORT, PRCP_WATER_ptr, 'PRCP_WATER', _RC)
         call MAPL_StateGetPointer(EXPORT, PRCP_RAIN_ptr, 'PRCP_RAIN', _RC)
         call MAPL_StateGetPointer(EXPORT, PRCP_SNOW_ptr, 'PRCP_SNOW', _RC)
         call MAPL_StateGetPointer(EXPORT, PRCP_ICE_ptr, 'PRCP_ICE', _RC)
         call MAPL_StateGetPointer(EXPORT, PRCP_GRAUPEL_ptr, 'PRCP_GRAUPEL', _RC)
         call MAPL_StateGetPointer(EXPORT, RHX_ptr, 'RHX', _RC)
         call MAPL_StateGetPointer(EXPORT, REV_LS_ptr, 'REV_LS', _RC)
         call MAPL_StateGetPointer(EXPORT, RSU_LS_ptr, 'RSU_LS', _RC)
         call MAPL_StateGetPointer(EXPORT, PFR_LS_ptr, 'PFR_LS', _RC)
         call MAPL_StateGetPointer(EXPORT, PFS_LS_ptr, 'PFS_LS', _RC)
         call MAPL_StateGetPointer(EXPORT, PFG_LS_ptr, 'PFG_LS', _RC)
         call MAPL_StateGetPointer(EXPORT, WTHV2_ptr, 'WTHV2', _RC)
         call MAPL_StateGetPointer(EXPORT, WQL_ptr, 'WQL', _RC)
         call MAPL_StateGetPointer(EXPORT, PDFITERS_ptr, 'PDFITERS', _RC)
         call MAPL_StateGetPointer(EXPORT, VFALL_ICE_ptr, 'VFALL_ICE', _RC)
         call MAPL_StateGetPointer(EXPORT, VFALL_SNOW_ptr, 'VFALL_SNOW', _RC)
         call MAPL_StateGetPointer(EXPORT, VFALL_GRAUPEL_ptr, 'VFALL_GRAUPEL', _RC)
         call MAPL_StateGetPointer(EXPORT, VFALL_RAIN_ptr, 'VFALL_RAIN', _RC)
         call MAPL_StateGetPointer(EXPORT, RHCRIT3D_ptr, 'RHCRIT', _RC)
         if (associated(RAD_CF_ptr)) RAD_CF_ptr = RAD_CF
         if (associated(RAD_QV_ptr)) RAD_QV_ptr = RAD_QV
         if (associated(RAD_QL_ptr)) RAD_QL_ptr = RAD_QL
         if (associated(RAD_QI_ptr)) RAD_QI_ptr = RAD_QI
         if (associated(RAD_QR_ptr)) RAD_QR_ptr = RAD_QR
         if (associated(RAD_QS_ptr)) RAD_QS_ptr = RAD_QS
         if (associated(RAD_QG_ptr)) RAD_QG_ptr = RAD_QG
         if (associated(EVAPC_ptr)) EVAPC_ptr = EVAPC
         if (associated(SUBLC_ptr)) SUBLC_ptr = SUBLC
         if (associated(PRCP_WATER_ptr)) PRCP_WATER_ptr = PRCP_WATER
         if (associated(PRCP_RAIN_ptr)) PRCP_RAIN_ptr = PRCP_RAIN
         if (associated(PRCP_SNOW_ptr)) PRCP_SNOW_ptr = PRCP_SNOW
         if (associated(PRCP_ICE_ptr)) PRCP_ICE_ptr = PRCP_ICE
         if (associated(PRCP_GRAUPEL_ptr)) PRCP_GRAUPEL_ptr = PRCP_GRAUPEL
         if (associated(RHX_ptr)) RHX_ptr = RHX
         if (associated(REV_LS_ptr)) REV_LS_ptr = REV_LS
         if (associated(RSU_LS_ptr)) RSU_LS_ptr = RSU_LS
         if (associated(PFR_LS_ptr)) PFR_LS_ptr = PFR_LS
         if (associated(PFS_LS_ptr)) PFS_LS_ptr = PFS_LS
         if (associated(PFG_LS_ptr)) PFG_LS_ptr = PFG_LS
         if (associated(WTHV2_ptr)) WTHV2_ptr = WTHV2
         if (associated(WQL_ptr)) WQL_ptr = WQL
         if (associated(PDFITERS_ptr)) PDFITERS_ptr = PDFITERS
         if (associated(VFALL_ICE_ptr)) VFALL_ICE_ptr = VFALL_ICE
         if (associated(VFALL_SNOW_ptr)) VFALL_SNOW_ptr = VFALL_SNOW
         if (associated(VFALL_GRAUPEL_ptr)) VFALL_GRAUPEL_ptr = VFALL_GRAUPEL
         if (associated(VFALL_RAIN_ptr)) VFALL_RAIN_ptr = VFALL_RAIN
         if (associated(RHCRIT3D_ptr)) RHCRIT3D_ptr = RHCRIT3D

         call MAPL_GridCompTimerStop(GC, "---CLDDIAGS", _RC)
      end if ! USE_PYMOIST_GFDL1M

      call MAPL_GridCompTimerStop(GC, "--GFDL_1M", _RC)

   end subroutine GFDL_1M_Run

   subroutine GFDL_1M_Finalize(GC, IMPORT, EXPORT, RC)

      type(ESMF_GridComp), intent(inout) :: GC ! Gridded component
      type(ESMF_State), intent(inout) :: IMPORT ! Import state
      type(ESMF_State), intent(inout) :: EXPORT ! Export state
      integer, optional, intent(out) :: RC ! Error code

      type(MAPL_MetaComp), pointer :: MAPL

      ! Get my internal MAPL_Generic state
      !-----------------------------------
      call MAPL_GetObjectFromGC(GC, MAPL, _RC)

      if (USE_PYMOIST_GFDL1M) then
         call MAPL_pybridge_gcfinalize("pyMoist.fortran.param_interfaces.microphysics.GFDL1M_interface", MAPL, IMPORT, &
              EXPORT)
      end if

   end subroutine GFDL_1M_Finalize

end module GEOS_GFDL_1M_InterfaceMod
