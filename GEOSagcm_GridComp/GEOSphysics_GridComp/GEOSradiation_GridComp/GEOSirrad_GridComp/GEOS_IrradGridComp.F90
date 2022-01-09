
!  $Id$

#include "MAPL_Generic.h"

module GEOS_IrradGridCompMod

!BOP
! !MODULE: GEOS_Irrad -- A Module to compute longwaves radiative transfer through a cloudy atmosphere

! !DESCRIPTION:
! 
!   {\tt Irrad} is a light-weight gridded component to compute longwave 
! radiative fluxes. It operates on the ESMF grid that appears in the
! gridded component passed to its {\tt Initialize} method. Unlike
! heavier gridded components, it does not enforce its own grid.
! The only restrictions are that it be a 3-dimensional grid
! in which one dimension is aligned with the vertical coordinate and
! only the horizontal dimensions are decomposed.
!
!   The radiative transfer calculation is based on M-D Chou's IRRAD routine.
! A full documentation of the code may be found in
! "A Thermal Infrared Radiation Parameterization for Atmospheric Studies"
! M.-D. Chou et al., NASA/TM-2001-104606, Vol. 19, 55 pp, 2003.
! Based on the 1996-version of the Air Force Geophysical Laboratory HITRAN data
! base (Rothman et al., 1998), the parameterization includes the absorption due
! to major gaseous absorption (water vapor, CO2 , O3 ) and most of the minor 
! trace gases (N2O, CH4 , CFC's), as well as clouds and aerosols. The thermal
! infrared spectrum is divided into nine bands and a subband. To achieve a high
! degree of accuracy and speed, various approaches of computing the transmission
! function are applied to different spectral bands and gases. The gaseous 
! transmission function is computed either using the k-distribution method or 
! the table look-up method. To include the effect of scattering due to clouds 
! and aerosols, the optical thickness is scaled by the single-scattering albedo
! and asymmetry factor. The optical thickness, the single-scattering albedo, 
! and the asymmetry factor of clouds are parameterized as functions of the ice
! and water content and the particle size.

!   All outputs are optional and are filled only if they have been
! initialized by a coupler. 
!
!   The net (+ve downward) fluxes are returned at the layer
! interfaces, which are indexed from the top of the atmosphere (L=0)
! to the surface. It also computes the sensitivity of net downward flux to 
! surface temperature and emission by the surface.
! The full transfer calculation, including the linearization w.r.t. the surface temperature,
! is done intermitently, on the component's main time step and its results are 
! kept in the internal state. Exports are refreshed each heartbeat based on the
! latest surface temperature.
!
!   Radiation should be called either before or after thos components
!    (usually SURFACE and DYNAMICS) that use its fluxes and modify
!    its inputs. If it is called before, the intemittent refresh should
!    occur during the first step of the radiation cycle, while if it
!    is called after, it should occur during the last step. The behavior
!    of the component needs to be somewhat different in these two cases
!    and so a means is provided, through the logical attribute \texttt{CALL\_LAST} in
!    configuration, of telling the component how it is being used. The 
!    default is \texttt{CALL\_LAST = "TRUE"}. 
!
!
! !USES:

  use ESMF
  use MAPL
  use GEOS_UtilsMod

  use rrtmg_lw_rad, only: rrtmg_lw
  use rrtmg_lw_init, only: rrtmg_lw_ini
  use parrrtm, only: ngptlw, nbndlw
  use rrlw_wvn, only: wavenum1, wavenum2

  ! for RRTMGP
  use mo_gas_optics_rrtmgp, only: ty_gas_optics_rrtmgp

#ifdef _CUDA
  use cudafor
  ! NOTE: USE renames are used below to prevent name clashes with
  !       CUDA copies to the GPU.
  use rad_constants, only: &
        AIB_IR_CONST=>AIB_IR, AWB_IR_CONST=>AWB_IR, &
        AIW_IR_CONST=>AIW_IR, AWW_IR_CONST=>AWW_IR, &
        AIG_IR_CONST=>AIG_IR, AWG_IR_CONST=>AWG_IR
  use irrad_constants, only: &
        XKW_CONST=>XKW, XKE_CONST=>XKE,  MW_CONST=>MW,  &
         AW_CONST=>AW,   BW_CONST=>BW,                  &
         PM_CONST=>PM,  FKW_CONST=>FKW, GKW_CONST=>GKW, &
         CB_CONST=>CB,  DCB_CONST=>DCB,                 &
        W11_CONST=>W11, W12_CONST=>W12, W13_CONST=>W13, &
        P11_CONST=>P11, P12_CONST=>P12, P13_CONST=>P13, &
        DWE_CONST=>DWE, DPE_CONST=>DPE,                 &
         C1_CONST=>C1,   C2_CONST=>C2,   C3_CONST=>C3,  &
        OO1_CONST=>OO1, OO2_CONST=>OO2, OO3_CONST=>OO3, &
        H11_CONST=>H11, H12_CONST=>H12, H13_CONST=>H13, &
        H21_CONST=>H21, H22_CONST=>H22, H23_CONST=>H23, &
        H81_CONST=>H81, H82_CONST=>H82, H83_CONST=>H83
  use irradmod, only: &
        ! Subroutines
        IRRAD, &
        ! Parameters
        NX, NO, NC, NH, &
        ! Inputs
        PLE_DEV, TA_DEV, WA_DEV, OA_DEV, TB_DEV, &
        N2O_DEV, CH4_DEV, CFC11_DEV, CFC12_DEV, CFC22_DEV, &
        FS_DEV, TG_DEV, EG_DEV, TV_DEV, EV_DEV, &
        RV_DEV, CWC_DEV, FCLD_DEV, REFF_DEV, &
        ! Aerosol inputs
        TAUA_DEV, SSAA_DEV, ASYA_DEV, &
        ! Constant arrays in global memory
         C1,  C2,  C3, &
        OO1, OO2, OO3, &
        H11, H12, H13, &
        H21, H22, H23, &
        H81, H82, H83, &
        ! Outputs
        FLXU_DEV, FLXAU_DEV, FLCU_DEV, FLAU_DEV, &
        FLXD_DEV, FLXAD_DEV, FLCD_DEV, FLAD_DEV, &
        DFDTS_DEV, SFCEM_DEV, TAUDIAG_DEV, &
        ! Constants
        XKW, XKE, MW, AW, BW, PM, FKW, &
        GKW, AIB_IR, AWB_IR, AIW_IR, AWW_IR, AIG_IR, AWG_IR, &
        CB, DCB, W11, W12, W13, P11, P12, &
        P13, DWE, DPE
#else
  use irradmod, only: IRRAD
#endif
  
  implicit none
  private

! !PUBLIC MEMBER FUNCTIONS:

  public SetServices

!EOP

   ! -------------------------------------------
   ! Select which RRTMG bands support OLR output
   ! -------------------------------------------
   !    via OLRBbbRG and TBRBbbRG exports ...
   ! (These exports require support space in the
   ! internal state so we choose only the ones we want
   ! to offer here at compile time. In the future, if
   ! most of the bands are required, then we will change
   ! strategy and reserve space for ALL of them in the
   ! internal state. In that case we would only require
   ! runtime band selection via the EXPORTS chosen.)

   ! NOTE: band 16 should be requested with caution ...
   ! Band 16 is technically 2600-3250 cm-1. But when RT
   ! is performed across all 16 bands, as it is in GEOS-5
   ! usage, then band 16 includes the integrated Planck
   ! values from 2600 cm-1 to infinity. So, the brightness
   ! temperature (Tbr) calculations in Update_Flx(), which
   ! use specified wavenumber endpoints, may require mod-
   ! ification for band 16 (pmn: TODO). For the moment,
   ! the limits [2600,3250] are used. 

   ! Which bands are supported?
   !    (Currently RRTMG only)
   !    (actual calculation only if export is requested)
   ! Supported?    Band  Requested by (and use)
   logical, parameter :: band_output_supported (nbndlw) = [ &
      .false. , &!  01
      .false. , &!  02
      .false. , &!  03
      .false. , &!  04
      .false. , &!  05
      .true.  , &!  06   A. Collow (Window)
      .false. , &!  07
      .false. , &!  08
      .true.  , &!  09   W. Putman (Water Vapor)
      .true.  , &!  10   W. Putman (Water Vapor)
      .true.  , &!  11   W. Putman (Water Vapor)
      .false. , &!  12
      .false. , &!  13
      .false. , &!  14
      .false. , &!  15
      .false. ]  !  16

   ! PS: We may later have an RRTMG internal state like
   ! RRTMGP below with various rrtmg_lw_init data, etc.
   ! TODO

   ! -----------------------------------------------
   ! RRTMGP internal state
   ! This will be attached to the Gridded Component
   ! used to provide efficient initialization
   type ty_RRTMGP_state
     private
     logical :: initialized = .false.
     type (ty_gas_optics_rrtmgp) :: k_dist
   end type ty_RRTMGP_state

   ! wrapper to access RRTMGP internal state
   type ty_RRTMGP_wrap
     type (ty_RRTMGP_state), pointer :: ptr => null()
   end type ty_RRTMGP_wrap

contains

!BOP
! !IROUTINE: SetServices -- Sets ESMF services for this component

! !INTERFACE:
  subroutine SetServices ( GC, RC )

! !ARGUMENTS:
    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    integer, optional                  :: RC  ! return code

! !DESCRIPTION: This version uses the MAPL\_GenericSetServices. This function sets
!                the Initialize and Finalize services, as well as allocating
!   our instance of a generic state and putting it in the 
!   gridded component (GC). Here we only need to set the run method and
!   add the state variable specifications (also generic) to our instance
!   of the generic state. This is the way our true state variables get into
!   the ESMF\_State INTERNAL, which is in the MAPL\_MetaComp.

!EOP

!=============================================================================

    character(len=ESMF_MAXSTR) :: IAm
    integer                    :: STATUS
    character(len=ESMF_MAXSTR) :: COMP_NAME

    type (ESMF_Config) :: CF

    integer :: MY_STEP
    integer :: ACCUMINT
    real    :: DT

    type (ty_RRTMGP_state), pointer :: rrtmgp_state => null()
    type (ty_RRTMGP_wrap)           :: wrap

    ! for OLRBbbRG, TBRBbbRG
    real :: RFLAG
    logical :: USE_RRTMG
    integer :: ibnd
    character*2 :: bb
    character*9 :: wvn_rng  ! xxxx-yyyy

!=============================================================================

! Get my name and set-up traceback handle
! ---------------------------------------

    Iam = 'SetServices'
    call ESMF_GridCompGet(GC, NAME=COMP_NAME, __RC__)
    Iam = trim(COMP_NAME) // Iam

    ! save pointer to the wrapped RRTMGP internal state in the GC
    allocate(rrtmgp_state, __STAT__)
    wrap%ptr => rrtmgp_state
    call ESMF_UserCompSetInternalState(GC, 'RRTMGP_state', wrap, status)
    VERIFY_(status)

! Set the Run entry point
! -----------------------

    call MAPL_GridCompSetEntryPoint(GC, ESMF_METHOD_RUN, Run, __RC__)

! Get the configuration
! ---------------------

    call ESMF_GridCompGet(GC, CONFIG=CF, __RC__)

! Get the intervals; "heartbeat" must exist
! -----------------------------------------

    call ESMF_ConfigGetAttribute(CF, DT, Label="RUN_DT:", __RC__)

! Refresh interval defaults to heartbeat. This will also be read by
! MAPL_Generic and set as the component's main time step.
! -----------------------------------------------------------------

    call ESMF_ConfigGetAttribute(CF, DT, Label=trim(COMP_NAME)//"_DT:", default=DT, __RC__)
    MY_STEP = nint(DT)

! Averaging interval defaults to the refresh interval.
!-----------------------------------------------------

    call ESMF_ConfigGetAttribute(CF, DT, Label=trim(COMP_NAME)//'Avrg:', default=DT, __RC__)
    ACCUMINT = nint(DT)

! Is RRTMG LW being run?
! ----------------------

    call ESMF_ConfigGetAttribute(CF, RFLAG, LABEL='USE_RRTMG_IRRAD:', DEFAULT=0., __RC__)
    USE_RRTMG = RFLAG /= 0.

! Set the state variable specs.
! -----------------------------

!BOS

!  !IMPORT STATE:

     call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME         = 'PLE',                               &
        LONG_NAME          = 'air_pressure',                      &
        UNITS              = 'Pa',                                &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationEdge,                  &
        AVERAGING_INTERVAL = ACCUMINT,                            &
        REFRESH_INTERVAL   = MY_STEP,                      __RC__ )

     call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME         = 'T',                                 &
        LONG_NAME          = 'air_temperature',                   &
        UNITS              = 'K',                                 &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
        AVERAGING_INTERVAL = ACCUMINT,                            &
        REFRESH_INTERVAL   = MY_STEP,                      __RC__ )

     call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME         = 'QV',                                &
        LONG_NAME          = 'specific_humidity',                 &
        UNITS              = 'kg kg-1',                           &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
        AVERAGING_INTERVAL = ACCUMINT,                            &
        REFRESH_INTERVAL   = MY_STEP,                      __RC__ )

     call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME         = 'QL',                                &
        LONG_NAME          = 'mass_fraction_of_cloud_liquid_water_in_air', &
        UNITS              = 'kg kg-1',                           &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
        AVERAGING_INTERVAL = ACCUMINT,                            &
        REFRESH_INTERVAL   = MY_STEP,                      __RC__ )

     call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME         = 'QI',                                &
        LONG_NAME          = 'mass_fraction_of_cloud_ice_in_air', &
        UNITS              = 'kg kg-1',                           &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
        AVERAGING_INTERVAL = ACCUMINT,                            &
        REFRESH_INTERVAL   = MY_STEP,                      __RC__ )

     call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME         = 'QR',                                &
        LONG_NAME          = 'mass_fraction_of_rain_water_in_air',&
        UNITS              = 'kg kg-1',                           &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
        AVERAGING_INTERVAL = ACCUMINT,                            &
        REFRESH_INTERVAL   = MY_STEP,                      __RC__ )

     call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME         = 'QS',                                &
        LONG_NAME          = 'mass_fraction_of_snow_in_air',      &
        UNITS              = 'kg kg-1',                           &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
        AVERAGING_INTERVAL = ACCUMINT,                            &
        REFRESH_INTERVAL   = MY_STEP,                      __RC__ )

     call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME         = 'RL',                                &
        LONG_NAME          = 'effective_radius_of_cloud_liquid_water_particles',      &
        UNITS              = 'm',                                 &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
        AVERAGING_INTERVAL = ACCUMINT,                            &
        REFRESH_INTERVAL   = MY_STEP,                      __RC__ )

     call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME         = 'RI',                                &
        LONG_NAME          = 'effective_radius_of_cloud_ice_particles',   &
        UNITS              = 'm',                                 &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
        AVERAGING_INTERVAL = ACCUMINT,                            &
        REFRESH_INTERVAL   = MY_STEP,                      __RC__ )

     call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME         = 'RR',                                &
        LONG_NAME          = 'effective_radius_of_rain_particles',&
        UNITS              = 'm',                                 &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
        AVERAGING_INTERVAL = ACCUMINT,                            &
        REFRESH_INTERVAL   = MY_STEP,                      __RC__ )

     call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME         = 'RS',                                &
        LONG_NAME          = 'effective_radius_of_snow_particles',&
        UNITS              = 'm',                                 &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
        AVERAGING_INTERVAL = ACCUMINT,                            &
        REFRESH_INTERVAL   = MY_STEP,                      __RC__ )

     call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME         = 'O3',                                &
        LONG_NAME          = 'ozone_mass_mixing_ratio',           &
        UNITS              = 'kg kg-1',                           &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
        AVERAGING_INTERVAL = ACCUMINT,                            &
        REFRESH_INTERVAL   = MY_STEP,                      __RC__ )

     call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME         = 'CH4',                               &
        LONG_NAME          = 'methane_concentration',             &
        UNITS              = 'pppv',                              &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
        AVERAGING_INTERVAL = ACCUMINT,                            &
        REFRESH_INTERVAL   = MY_STEP,                      __RC__ )

     call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME         = 'N2O',                               &
        LONG_NAME          = 'nitrous_oxide_concentration',       &
        UNITS              = 'pppv',                              &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
        AVERAGING_INTERVAL = ACCUMINT,                            &
        REFRESH_INTERVAL   = MY_STEP,                      __RC__ )

     call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME         = 'CFC11',                             &
        LONG_NAME          = 'CFC11_concentration',               &
        UNITS              = 'pppv',                              &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
        AVERAGING_INTERVAL = ACCUMINT,                            &
        REFRESH_INTERVAL   = MY_STEP,                      __RC__ )

     call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME         = 'CFC12',                             &
        LONG_NAME          = 'CFC12_concentration',               &
        UNITS              = 'pppv',                              &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
        AVERAGING_INTERVAL = ACCUMINT,                            &
        REFRESH_INTERVAL   = MY_STEP,                      __RC__ )

     call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME         = 'HCFC22',                            &
        LONG_NAME          = 'HCFC22_concentration',              &
        UNITS              = 'pppv',                              &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
        AVERAGING_INTERVAL = ACCUMINT,                            &
        REFRESH_INTERVAL   = MY_STEP,                      __RC__ )

     call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME         = 'FCLD',                              &
        LONG_NAME          = 'cloud_area_fraction_in_atmosphere_layer', &
        UNITS              = '1',                                 &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
        AVERAGING_INTERVAL = ACCUMINT,                            &
        REFRESH_INTERVAL   = MY_STEP,                      __RC__ )

     call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME         = 'TS',                                &
        LONG_NAME          = 'surface_skin_temperature',          &
        UNITS              = 'K',                                 &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        AVERAGING_INTERVAL = ACCUMINT,                            &
        REFRESH_INTERVAL   = MY_STEP,                      __RC__ )

     call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME         = 'EMIS',                              &
        LONG_NAME          = 'surface_emissivity',                &
        UNITS              = '1',                                 &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        AVERAGING_INTERVAL = ACCUMINT,                            &
        REFRESH_INTERVAL   = MY_STEP,                      __RC__ )

     call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME         = 'PREF',                              &
        LONG_NAME          = 'reference_air_pressure',            &
        UNITS              = 'Pa',                                &
        DIMS               = MAPL_DimsVertOnly,                   &
        VLOCATION          = MAPL_VLocationEdge,           __RC__ )

! Instantaneous TS is used only for updating the IR fluxes due to TS change

     call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME         = 'TSINST',                            &
        LONG_NAME          = 'surface_skin_temperature',          &
        UNITS              = 'K',                                 &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,           __RC__ )


    call MAPL_AddImportSpec(GC,                                   &
       LONG_NAME  = 'aerosols',                                   &
       UNITS      = 'kg kg-1',                                    &
       SHORT_NAME = 'AERO',                                       &
       DIMS       = MAPL_DimsHorzVert,                            &
       VLOCATION  = MAPL_VLocationCenter,                         &
       DATATYPE   = MAPL_StateItem,                               &
       RESTART    = MAPL_RestartSkip,                      __RC__ )

!  !EXPORT STATE:

    call MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME = 'FLX',                                       &
        LONG_NAME  = 'net_downward_longwave_flux_in_air',         &
        UNITS      = 'W m-2',                                     &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationEdge,                   __RC__ )

    call MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME = 'FLXA',                                      &
        LONG_NAME  = 'net_downward_longwave_flux_in_air_and_no_aerosol', &
        UNITS      = 'W m-2',                                     &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationEdge,                   __RC__ )

    call MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME = 'FLXD',                                      &
        LONG_NAME  = 'downward_longwave_flux_in_air',             &
        UNITS      = 'W m-2',                                     &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationEdge,                   __RC__ )

    call MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME = 'FLXAD',                                     &
        LONG_NAME  = 'downward_longwave_flux_in_air_and_no_aerosol', &
        UNITS      = 'W m-2',                                     &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationEdge,                   __RC__ )

    call MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME = 'FLXU',                                      &
        LONG_NAME  = 'upward_longwave_flux_in_air',               &
        UNITS      = 'W m-2',                                     &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationEdge,                   __RC__ )

    call MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME = 'FLXAU',                                     &
        LONG_NAME  = 'upward_longwave_flux_in_air_and_no_aerosol',&
        UNITS      = 'W m-2',                                     &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationEdge,                   __RC__ )

    call MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME = 'FLC',                                       &
        LONG_NAME  = 'net_downward_longwave_flux_in_air_assuming_clear_sky', &
        UNITS      = 'W m-2',                                     &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationEdge,                   __RC__ )

    call MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME = 'FLCD',                                      &
        LONG_NAME  = 'downward_longwave_flux_in_air_assuming_clear_sky', &
        UNITS      = 'W m-2',                                     &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationEdge,                   __RC__ )
    
    call MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME = 'FLCU',                                      &
        LONG_NAME  = 'upward_longwave_flux_in_air_assuming_clear_sky', &
        UNITS      = 'W m-2',                                     &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationEdge,                   __RC__ )

    call MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME = 'FLA',                                       &
        LONG_NAME  = 'net_downward_longwave_flux_in_air_assuming_clear_sky_and_no_aerosol',&
        UNITS      = 'W m-2',                                     &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationEdge,                   __RC__ )

    call MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME = 'FLAD',                                      &
        LONG_NAME  = 'downward_longwave_flux_in_air_assuming_clear_sky_and_no_aerosol', &
        UNITS      = 'W m-2',                                     &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationEdge,                   __RC__ )

    call MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME = 'FLAU',                                      &
        LONG_NAME  = 'upward_longwave_flux_in_air_assuming_clear_sky_and_no_aerosol', &
        UNITS      = 'W m-2',                                     &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationEdge,                   __RC__ )

    call MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME = 'SFCEM',                                     &
        LONG_NAME  = 'longwave_flux_emitted_from_surface',        &
        UNITS      = 'W m-2',                                     &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,                   __RC__ )

    call MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME = 'SFCEM0',                                    &
        LONG_NAME  = 'longwave_flux_emitted_from_surface_at_reference_time',&
        UNITS      = 'W m-2',                                     &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,                   __RC__ )

    call MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME = 'LWS0',                                      &
        LONG_NAME  = 'surface_absorbed_longwave_radiation_at_reference_time',&
        UNITS      = 'W m-2',                                     &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,                   __RC__ )

    call MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME = 'DSFDTS',                                    &
        LONG_NAME  = 'sensitivity_of_longwave_flux_emitted_from_surface_to_surface_temperature', &
        UNITS      = 'W m-2 K-1',                                 &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,                   __RC__ )

    call MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME = 'DSFDTS0',                                   &
        LONG_NAME  = 'sensitivity_of_longwave_flux_emitted_from_surface_to_surface_temperature_at_reference_time', &
        UNITS      = 'W m-2 K-1',                                 &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,                   __RC__ )

    call MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME = 'TSREFF',                                    &
        LONG_NAME  = 'surface_temperature',                       &
        UNITS      = 'K',                                         &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,                   __RC__ )

    call MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME = 'OLR',                                       &
        LONG_NAME  = 'upwelling_longwave_flux_at_toa',            &
        UNITS      = 'W m-2',                                     &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,                   __RC__ )

    call MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME = 'OLRA',                                      &
        LONG_NAME  = 'upwelling_longwave_flux_at_toa_and_no_aerosol', &
        UNITS      = 'W m-2',                                     &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,                   __RC__ )

    call MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME = 'OLC',                                       &
        LONG_NAME  = 'upwelling_longwave_flux_at_toa_assuming_clear_sky',&
        UNITS      = 'W m-2',                                     &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,                   __RC__ )

    call MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME = 'OLCC5',                                     &
        LONG_NAME  = 'upwelling_longwave_flux_at_toa_assuming_clear_sky_masked_using_cldtt_LE_5',&
        UNITS      = 'W m-2',                                     &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,                   __RC__ )

    call MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME = 'OLA',                                       &
        LONG_NAME  = 'upwelling_longwave_flux_at_toa_assuming_clear_sky_and_no_aerosol',&
        UNITS      = 'W m-2',                                     &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,                   __RC__ )

    if (USE_RRTMG) then
       do ibnd = 1,nbndlw
          if (band_output_supported(ibnd)) then
             write(bb,'(I0.2)') ibnd
             write(wvn_rng,'(I0,"-",I0)') nint(wavenum1(ibnd)), nint(wavenum2(ibnd))
   
             call MAPL_AddExportSpec(GC,                                    &
                SHORT_NAME = 'OLRB'//bb//'RG',                              &
                LONG_NAME  = 'upwelling_longwave_flux_at_TOA_in_RRTMG_band' &
                                //bb//' ('//trim(wvn_rng)//' cm-1)',        &
                UNITS      = 'W m-2',                                       &
                DIMS       = MAPL_DimsHorzOnly,                             &
                VLOCATION  = MAPL_VLocationNone,                     __RC__ )

             call MAPL_AddExportSpec(GC,                                    &
                SHORT_NAME = 'TBRB'//bb//'RG',                              &
                LONG_NAME  = 'brightness_temperature_in_RRTMG_band'         &
                                //bb//' ('//trim(wvn_rng)//' cm-1)',        &
                UNITS      = 'K',                                           &
                DIMS       = MAPL_DimsHorzOnly,                             &
                VLOCATION  = MAPL_VLocationNone,                     __RC__ )

          end if
       end do
    end if

    call MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME = 'FLNS',                                      &
        LONG_NAME  = 'surface_net_downward_longwave_flux',        &
        UNITS      = 'W m-2',                                     &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,                   __RC__ )

    call MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME = 'FLNSNA',                                    &
        LONG_NAME  = 'surface_net_downward_longwave_flux_and_no_aerosol', &
        UNITS      = 'W m-2',                                     &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,                   __RC__ )

    call MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME = 'FLNSC',                                     &
        LONG_NAME  = 'surface_net_downward_longwave_flux_assuming_clear_sky',&
        UNITS      = 'W m-2',                                     &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,                   __RC__ )

    call MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME = 'FLNSA',                                     &
        LONG_NAME  = 'surface_net_downward_longwave_flux_assuming_clear_sky_and_no_aerosol',&
        UNITS      = 'W m-2',                                     &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,                   __RC__ )

    call MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME = 'LWS',                                       &
        LONG_NAME  = 'surface_absorbed_longwave_radiation',       &
        UNITS      = 'W m-2',                                     &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,                   __RC__ )

    call MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME = 'LWSA',                                      &
        LONG_NAME  = 'surface_absorbed_longwave_radiation_and_no_aerosol', &
        UNITS      = 'W m-2',                                     &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,                   __RC__ )

    call MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME = 'LCS',                                       &
        LONG_NAME  = 'surface_absorbed_longwave_radiation_assuming_clear_sky',&
        UNITS      = 'W m-2',                                     &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,                   __RC__ )

    call MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME = 'LCSC5',                                     &
        LONG_NAME  = 'surface_absorbed_longwave_radiation_assuming_clear_sky_masked_using_cldtt_LE_5',&
        UNITS      = 'W m-2',                                     &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,                   __RC__ )

    call MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME = 'LAS',                                       &
        LONG_NAME  = 'surface_absorbed_longwave_radiation_assuming_clear_sky_and_no_aerosol',&
        UNITS      = 'W m-2',                                     &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,                   __RC__ )

    call MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME = 'CLDTMP',                                    &
        LONG_NAME  = 'cloud_top_temperature',                     &
        UNITS      = 'K',                                         &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,                   __RC__ )

    call MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME = 'CLDPRS',                                    &
        LONG_NAME  = 'cloud_top_pressure',                        &
        UNITS      = 'Pa',                                        &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,                   __RC__ )

    call MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME = 'TAUIR',                                     &
        LONG_NAME  = 'longwave_cloud_optical_thickness_at_800_cm-1',&
        UNITS      = 'W m-2',                                     &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,                 __RC__ )

    call MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME = 'CLDTT'  ,                                   &
        LONG_NAME  = 'total_2D_cloud_area_fraction',              &
        UNITS      = '1',                                         &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,                   __RC__ )

    call MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME = 'CLDTTLW',                                   &
        LONG_NAME  = 'total_cloud_area_fraction_rrtmg_lw',        &
        UNITS      = '1',                                         &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,                   __RC__ )

    call MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME = 'CLDHILW',                                   &
        LONG_NAME  = 'high-level_cloud_area_fraction_rrtmg_lw',   &
        UNITS      = '1',                                         &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,                   __RC__ )

    call MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME = 'CLDMDLW',                                   &
        LONG_NAME  = 'mid-level_cloud_area_fraction_rrtmg_lw',    &
        UNITS      = '1',                                         &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,                   __RC__ )

    call MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME = 'CLDLOLW',                                   &
        LONG_NAME  = 'low_level_cloud_area_fraction_rrtmg_lw',    &
        UNITS      = '1',                                         &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,                   __RC__ )

!  Irrad does not have a "real" internal state. To update the net_longwave_flux
!  due to the change of surface temperature every time step, we keep 
!  several variables in the internal state.

!  !INTERNAL STATE:

    call MAPL_AddInternalSpec(GC,                                 &
        SHORT_NAME = 'FLX',                                       &
        LONG_NAME  = 'net_downward_longwave_flux_in_air',         &
        UNITS      = 'W m-2',                                     &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationEdge,                   __RC__ )

    call MAPL_AddInternalSpec(GC,                                 &
        SHORT_NAME = 'FLC',                                       &
        LONG_NAME  = 'net_downward_longwave_flux_in_air_for_clear_sky',&
        UNITS      = 'W m-2',                                     &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationEdge,                   __RC__ )

    call MAPL_AddInternalSpec(GC,                                 &
        SHORT_NAME = 'FLA',                                       &
        LONG_NAME  = 'net_downward_longwave_flux_in_air_for_clear_sky_and_no_aerosol',  &
        UNITS      = 'W m-2',                                     &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationEdge,                   __RC__ )

    call MAPL_AddInternalSpec(GC,                                 &
        SHORT_NAME = 'FLXD',                                      &
        LONG_NAME  = 'downward_longwave_flux_in_air',             &
        UNITS      = 'W m-2',                                     &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationEdge,                   __RC__ )

    call MAPL_AddInternalSpec(GC,                                 &
        SHORT_NAME = 'FLXU',                                      &
        LONG_NAME  = 'upward_longwave_flux_in_air',               &
        UNITS      = 'W m-2',                                     &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationEdge,                   __RC__ )

    call MAPL_AddInternalSpec(GC,                                 &
        SHORT_NAME = 'FLCD',                                      &
        LONG_NAME  = 'downward_longwave_flux_in_air_for_clear_sky',&
        UNITS      = 'W m-2',                                     &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationEdge,                   __RC__ )

    call MAPL_AddInternalSpec(GC,                                 &
        SHORT_NAME = 'FLCU',                                      &
        LONG_NAME  = 'upward_longwave_flux_in_air_for_clear_sky', &
        UNITS      = 'W m-2',                                     &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationEdge,                   __RC__ )

    call MAPL_AddInternalSpec(GC,                                 &
        SHORT_NAME = 'FLAD',                                      &
        LONG_NAME  = 'downward_longwave_flux_in_air_for_clear_sky_and_no_aerosol',&
        UNITS      = 'W m-2',                                     &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationEdge,                   __RC__ )

    call MAPL_AddInternalSpec(GC,                                 &
        SHORT_NAME = 'FLAU',                                      &
        LONG_NAME  = 'upward_longwave_flux_in_air_for_clear_sky_and_no_aerosol',&
        UNITS      = 'W m-2',                                     &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationEdge,                   __RC__ )

    call MAPL_AddInternalSpec(GC,                                 &
        SHORT_NAME = 'DFDTS',                                     &
        LONG_NAME  = 'sensitivity_of_net_downward_longwave_flux_in_air_to_surface_temperature',&
        UNITS      = 'W m-2 K-1',                                 &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationEdge,                   __RC__ )

    call MAPL_AddInternalSpec(GC,                                 &
        SHORT_NAME = 'DFDTSC',                                    &
        LONG_NAME  = 'sensitivity_of_net_downward_longwave_flux_in_air_to_surface_temperature_for_clear_sky',&
        UNITS      = 'W m-2 K-1',                                 &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationEdge,                   __RC__ )

    call MAPL_AddInternalSpec(GC,                                 &
        SHORT_NAME = 'DFDTSNA',                                   &
        LONG_NAME  = 'sensitivity_of_net_downward_longwave_flux_in_air_to_surface_temperature_no_aerosol',&
        UNITS      = 'W m-2 K-1',                                 &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationEdge,                   __RC__ )

    call MAPL_AddInternalSpec(GC,                                 &
        SHORT_NAME = 'DFDTSCNA',                                  &
        LONG_NAME  = 'sensitivity_of_net_downward_longwave_flux_in_air_to_surface_temperature_for_clear_sky_no_aerosol',&
        UNITS      = 'W m-2 K-1',                                 &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationEdge,                   __RC__ )

    call MAPL_AddInternalSpec(GC,                                 &
        SHORT_NAME = 'SFCEM',                                     &
        LONG_NAME  = 'longwave_flux_emitted_from_surface',        &
        UNITS      = 'W m-2',                                     &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,                   __RC__ )

    call MAPL_AddInternalSpec(GC,                                 &
        SHORT_NAME = 'TS',                                        &
        LONG_NAME  = 'surface_temperature',                       &
        UNITS      = 'K',                                         &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,                   __RC__ )

    call MAPL_AddInternalSpec(GC,                                 &
        SHORT_NAME = 'FLXA',                                      &
        LONG_NAME  = 'net_downward_longwave_flux_in_air_and_no_aerosol',&
        UNITS      = 'W m-2',                                     &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationEdge,                   __RC__ )

    call MAPL_AddInternalSpec(GC,                                 &
        SHORT_NAME = 'FLXAD',                                     &
        LONG_NAME  = 'downward_longwave_flux_in_air_and_no_aerosol',&
        UNITS      = 'W m-2',                                     &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationEdge,                   __RC__ )

    call MAPL_AddInternalSpec(GC,                                 &
        SHORT_NAME = 'FLXAU',                                     &
        LONG_NAME  = 'upward_longwave_flux_in_air_and_no_aerosol',&
        UNITS      = 'W m-2',                                     &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationEdge,                   __RC__ )

    if (USE_RRTMG) then
       do ibnd = 1,nbndlw
          if (band_output_supported(ibnd)) then
             write(bb,'(I0.2)') ibnd

             call MAPL_AddInternalSpec(GC,                                       &
                SHORT_NAME = 'OLRB'//bb//'RG',                                   &
                LONG_NAME  = 'upwelling_longwave_flux_at_TOA_in_RRTMG_band'//bb, &
                UNITS      = 'W m-2',                                            &
                DIMS       = MAPL_DimsHorzOnly,                                  &
                VLOCATION  = MAPL_VLocationNone,                          __RC__ )

             call MAPL_AddInternalSpec(GC,                                       &
                SHORT_NAME = 'DOLRB'//bb//'RGDT',                                &
                LONG_NAME  = 'derivative_of_upwelling_longwave_flux_at_TOA'//    &
                                '_in_RRTMG_band'//bb//'_wrt_surface_temp',       &
                UNITS      = 'W m-2 K-1',                                        &
                DIMS       = MAPL_DimsHorzOnly,                                  &
                VLOCATION  = MAPL_VLocationNone,                          __RC__ )

          end if
       end do
    end if

!EOS

! Set the Profiling timers
! ------------------------

    call MAPL_TimerAdd(GC, name="-LW_DRIVER"               , __RC__)
    call MAPL_TimerAdd(GC, name="--IRRAD"                  , __RC__)
    call MAPL_TimerAdd(GC, name="---IRRAD_RUN"             , __RC__)
    call MAPL_TimerAdd(GC, name="---IRRAD_DATA"            , __RC__)
    call MAPL_TimerAdd(GC, name="----IRRAD_DATA_DEVICE"    , __RC__)
    call MAPL_TimerAdd(GC, name="----IRRAD_DATA_CONST"     , __RC__)
    call MAPL_TimerAdd(GC, name="---IRRAD_ALLOC"           , __RC__)
    call MAPL_TimerAdd(GC, name="---IRRAD_DEALLOC"         , __RC__)
    call MAPL_TimerAdd(GC, name="--RRTMG"                  , __RC__)
    call MAPL_TimerAdd(GC, name="---RRTMG_RUN"             , __RC__)
    call MAPL_TimerAdd(GC, name="---RRTMG_INIT"            , __RC__)
    call MAPL_TimerAdd(GC, name="---RRTMG_FLIP"            , __RC__)
    call MAPL_TimerAdd(GC, name="--RRTMGP"                 , __RC__)
    call MAPL_TimerAdd(GC, name="---RRTMGP_SETUP_1"        , __RC__)
    call MAPL_TimerAdd(GC, name="---RRTMGP_SETUP_2"        , __RC__)
    call MAPL_TimerAdd(GC, name="---RRTMGP_SETUP_3"        , __RC__)
    call MAPL_TimerAdd(GC, name="---RRTMGP_SETUP_4"        , __RC__)
    call MAPL_TimerAdd(GC, name="---RRTMGP_IO_1"           , __RC__)
    call MAPL_TimerAdd(GC, name="---RRTMGP_IO_2"           , __RC__)
    call MAPL_TimerAdd(GC, name="---RRTMGP_CLOUD_OPTICS"   , __RC__)
    call MAPL_TimerAdd(GC, name="---RRTMGP_AEROSOL_SETUP"  , __RC__)
    call MAPL_TimerAdd(GC, name="---RRTMGP_SUBSET"         , __RC__)
    call MAPL_TimerAdd(GC, name="---RRTMGP_MCICA"          , __RC__)
    call MAPL_TimerAdd(GC, name="---RRTMGP_GAS_OPTICS"     , __RC__)
    call MAPL_TimerAdd(GC, name="---RRTMGP_RT"             , __RC__)
    call MAPL_TimerAdd(GC, name="---RRTMGP_POST"           , __RC__)
    call MAPL_TimerAdd(GC, name="--MISC"                   , __RC__)
    call MAPL_TimerAdd(GC, name="---AEROSOLS"              , __RC__)
    call MAPL_TimerAdd(GC, name="-UPDATE_FLX"              , __RC__)

! Set generic init and final methods
! ----------------------------------

    call MAPL_GenericSetServices(GC, __RC__)

    RETURN_(ESMF_SUCCESS)
  
  end subroutine SetServices

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!BOP
! !IROUTINE: RUN -- Run method for the LW component

! !INTERFACE:
subroutine RUN ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:
  type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
  type(ESMF_State),    intent(inout) :: IMPORT ! Import state
  type(ESMF_State),    intent(inout) :: EXPORT ! Export state
  type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
  integer, optional,   intent(  out) :: RC     ! Error code:

! !DESCRIPTION: Periodically refreshes the fluxes and their derivatives
!                w.r.t surface skin temperature. On every step it produces
!                a linear estimate of the fluxes based on the instantaneous
!                surface temperature.

!EOP


! ErrLog Variables

  character(len=ESMF_MAXSTR)         :: IAm
  integer                            :: STATUS
  character(len=ESMF_MAXSTR)         :: COMP_NAME

! Local derived type aliases

  type (MAPL_MetaComp),     pointer  :: MAPL
  type (ESMF_Grid)                   :: ESMFGRID
  type (ESMF_State)                  :: INTERNAL
  type (ESMF_Alarm)                  :: ALARM

  integer                            :: IM, JM, LM
  integer                            :: CalledLast

  type (ty_RRTMGP_state), pointer    :: rrtmgp_state => null()
  type (ty_RRTMGP_wrap)              :: wrap

! Pointers to internal

   real, pointer, dimension(:,:  )   :: SFCEM_INT
   real, pointer, dimension(:,:  )   :: TS_INT
   real, pointer, dimension(:,:,:)   :: FLX_INT
   real, pointer, dimension(:,:,:)   :: FLXA_INT
   real, pointer, dimension(:,:,:)   :: FLC_INT
   real, pointer, dimension(:,:,:)   :: FLA_INT
   real, pointer, dimension(:,:,:)   :: FLXU_INT
   real, pointer, dimension(:,:,:)   :: FLXAU_INT
   real, pointer, dimension(:,:,:)   :: FLCU_INT
   real, pointer, dimension(:,:,:)   :: FLAU_INT
   real, pointer, dimension(:,:,:)   :: FLXD_INT
   real, pointer, dimension(:,:,:)   :: FLXAD_INT
   real, pointer, dimension(:,:,:)   :: FLCD_INT
   real, pointer, dimension(:,:,:)   :: FLAD_INT

   real, pointer, dimension(:,:,:)   :: DFDTS
   real, pointer, dimension(:,:,:)   :: DFDTSNA
   real, pointer, dimension(:,:,:)   :: DFDTSC
   real, pointer, dimension(:,:,:)   :: DFDTSCNA

   real, external :: getco2

! Concerning what radiation to use

   logical :: USE_RRTMGP, USE_RRTMGP_SORAD
   logical :: USE_RRTMG , USE_RRTMG_SORAD
   logical :: USE_CHOU  , USE_CHOU_SORAD
   real    :: RFLAG

! Additional pointers for RRTMG

   real, pointer, dimension(:,:  )   :: LONS
   real, pointer, dimension(:,:  )   :: LATS

   type (ESMF_VM)                    :: VM
   integer                           :: CoresPerNode
   integer                           :: COMM

   ! which bands require OLR output?
   ! (only RRTMG currently; OLRBbbRG, TBRBbbRG)
   real, pointer, dimension(:,:) :: ptr2d
   logical :: band_output (nbndlw)
   integer :: ibnd
   character*2 :: bb

!=============================================================================

! Begin... 

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

   Iam = "Run"
   call ESMF_GridCompGet( GC, name=COMP_NAME, GRID=ESMFGRID, RC=STATUS )
   VERIFY_(STATUS)
   Iam = trim(COMP_NAME) // Iam

! Get my internal MAPL_Generic state
!-----------------------------------

   call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
   VERIFY_(STATUS)

   call MAPL_TimerOn(MAPL,"TOTAL")

! Get parameters from generic state. The RUNALARM is used to control
!  the calling of the full transfer calculation
!-------------------------------------------------------------------

   call MAPL_Get(MAPL,              &
         IM=IM, JM=JM, LM=LM,                    &
         LONS=LONS, LATS=LATS,                   &
         RUNALARM            = ALARM,            &
         INTERNAL_ESMF_STATE = INTERNAL,         &
                                       RC=STATUS )
   VERIFY_(STATUS)

! Get number of cores per node for RRTMG GPU
! ------------------------------------------

   call ESMF_VMGetCurrent(VM, RC=STATUS)
   VERIFY_(STATUS)

   call ESMF_VmGet(VM, mpiCommunicator=COMM, RC=STATUS)
   VERIFY_(STATUS)

   CoresPerNode = MAPL_CoresPerNodeGet(COMM,RC=STATUS)
   VERIFY_(STATUS)

! Decide which radiation to use
! RRTMGP dominates RRTMG dominates Chou-Suarez
! Chou-Suarez is the default if nothing else asked for in Resource file
! These USE_ flags are shared globally by contained LW_Driver() and Update_Flx()
!-------------------------------------------------------------------------------

   ! first for IRRAD
   USE_RRTMGP = .false.
   USE_RRTMG  = .false.
   USE_CHOU   = .false.
   call MAPL_GetResource( MAPL, RFLAG ,'USE_RRTMGP_IRRAD:', DEFAULT=0., __RC__)
   USE_RRTMGP = RFLAG /= 0.
   if (.not. USE_RRTMGP) then
     call MAPL_GetResource( MAPL, RFLAG ,'USE_RRTMG_IRRAD:', DEFAULT=0., __RC__)
     USE_RRTMG = RFLAG /= 0.
     USE_CHOU  = .not.USE_RRTMG
   end if

   ! then SOLAR
   USE_RRTMGP_SORAD = .false.
   USE_RRTMG_SORAD  = .false.
   USE_CHOU_SORAD   = .false.
   call MAPL_GetResource( MAPL, RFLAG ,'USE_RRTMGP_SORAD:', DEFAULT=0., __RC__)
   USE_RRTMGP_SORAD = RFLAG /= 0.
   if (.not. USE_RRTMGP_SORAD) then
     call MAPL_GetResource( MAPL, RFLAG ,'USE_RRTMG_SORAD:', DEFAULT=0., __RC__)
     USE_RRTMG_SORAD = RFLAG /= 0.
     USE_CHOU_SORAD  = .not.USE_RRTMG_SORAD
   end if

   ! select which bands require OLRB output ...
   ! ------------------------------------------
   ! Currently only available for RRTMG
   ! must be supported AND requested by export 'OLRBbbRG' OR 'TBRBbbRG'
   if (USE_RRTMG) then
      do ibnd = 1,nbndlw
         band_output(ibnd) = .false.
         if (.not. band_output_supported(ibnd)) cycle
         write(bb,'(I0.2)') ibnd
         call MAPL_GetPointer(EXPORT, ptr2d, 'OLRB'//bb//'RG', __RC__)
         if (associated(ptr2d)) then
            band_output(ibnd) = .true.
            cycle
         end if
         call MAPL_GetPointer(EXPORT, ptr2d, 'TBRB'//bb//'RG', __RC__)
         if (associated(ptr2d)) then
            band_output(ibnd) = .true.
            cycle
         end if
      end do
   end if

! Pointers to Internals; these are needed by both Update and Refresh
!-------------------------------------------------------------------

   call MAPL_GetPointer(INTERNAL, SFCEM_INT, 'SFCEM', RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(INTERNAL, FLX_INT,   'FLX',   RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(INTERNAL, FLXA_INT,  'FLXA',  RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(INTERNAL, FLC_INT,   'FLC',   RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(INTERNAL, FLA_INT,   'FLA',   RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(INTERNAL, FLXU_INT,  'FLXU',  RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(INTERNAL, FLXAU_INT, 'FLXAU', RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(INTERNAL, FLCU_INT,  'FLCU',  RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(INTERNAL, FLAU_INT,  'FLAU',  RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(INTERNAL, FLXD_INT,  'FLXD',  RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(INTERNAL, FLXAD_INT, 'FLXAD', RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(INTERNAL, FLCD_INT,  'FLCD',  RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(INTERNAL, FLAD_INT,  'FLAD',  RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(INTERNAL, TS_INT,    'TS',    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(INTERNAL, DFDTS,     'DFDTS', RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(INTERNAL, DFDTSC,    'DFDTSC',RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(INTERNAL, DFDTSNA,   'DFDTSNA', RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(INTERNAL, DFDTSCNA,  'DFDTSCNA',RC=STATUS); VERIFY_(STATUS)

! Determine calling sequence
!---------------------------

   call MAPL_GetResource(MAPL,CalledLast,'CALLED_LAST:', default=1, RC=STATUS)
   VERIFY_(STATUS)

! Fill exported fluxed based on latest Ts
!----------------------------------------

   if(CalledLast/=0) then
      call MAPL_TimerOn(MAPL,"-UPDATE_FLX")
       call Update_Flx( IM,JM,LM, RC=STATUS )
       VERIFY_(STATUS)
      call MAPL_TimerOff(MAPL,"-UPDATE_FLX")
   endif

! If it is time, refresh internal state.
!---------------------------------------
   
   if ( ESMF_AlarmIsRinging   (ALARM, RC=STATUS) ) then
      call ESMF_AlarmRingerOff(ALARM, RC=STATUS)
      VERIFY_(STATUS)

      call MAPL_TimerOn(MAPL,"-LW_DRIVER")
       call LW_Driver( IM,JM,LM,LATS,LONS,CoresPerNode, RC=STATUS )
       VERIFY_(STATUS)
      call MAPL_TimerOff(MAPL,"-LW_DRIVER")
 
   endif

! Fill exported fluxes based on latest Ts
!----------------------------------------

   if(CalledLast==0) then
      call MAPL_TimerOn(MAPL,"-UPDATE_FLX")
       call Update_Flx( IM,JM,LM, RC=STATUS )
       VERIFY_(STATUS)
      call MAPL_TimerOff(MAPL,"-UPDATE_FLX")
   endif

   call MAPL_TimerOff(MAPL,"TOTAL")

!  All done
!-----------

   RETURN_(ESMF_SUCCESS)

contains



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine Lw_Driver(IM,JM,LM,LATS,LONS,CoresPerNode,RC)

   ! RRTMGP module uses
   use mo_rte_kind,                only: wp
   use mo_gas_concentrations,      only: ty_gas_concs
   use mo_cloud_optics,            only: ty_cloud_optics
   use mo_cloud_sampling,          only: draw_samples, &
                                         sampled_mask_max_ran, sampled_mask_exp_ran
   use mo_optical_props,           only: ty_optical_props, &
                                         ty_optical_props_arry, ty_optical_props_1scl, &
                                         ty_optical_props_2str, ty_optical_props_nstr
   use mo_source_functions,        only: ty_source_func_lw
   use mo_fluxes,                  only: ty_fluxes_broadband
   use mo_rte_lw,                  only: rte_lw
   use mo_load_coefficients,       only: load_and_init
   use mo_load_cloud_coefficients, only: load_cld_lutcoeff, load_cld_padecoeff

   ! Type of MKL VSL Basic RNGs
   ! (1) Mersenne Twister types
   ! brng = VSL_BRNG_MT19937
   ! Alternatives are VSL_BRNG_SFMT19937, maybe VSL_BRNG_MT2203?
   ! (2) Counter based PRNGs (CBPRNGs)
   ! brng = VSL_BRNG_PHILOX4X32X10  ! 10-round Philox 4x32 counter, 2x32 key
   ! Alternatives are VSL_BRNG_ARS5 ! faster if AES-NI instructions hardware supported
   !
#ifdef HAVE_MKL
   use MKL_VSL_TYPE
   use mo_rng_mklvsl_plus, only: ty_rng_mklvsl_plus
#endif

   integer,                   intent(IN )    :: IM, JM, LM, CoresPerNode
   real,    dimension(IM,JM), intent(IN )    :: LATS, LONS
   integer, optional,         intent(OUT)    :: RC
     
!  Locals

   character(len=ESMF_MAXSTR)        :: IAm
   integer                           :: STATUS

! local variables

   logical, parameter :: TRACE    = .true.

   integer, parameter :: NS       = 1       ! number of sub-grid surface types 

   integer, parameter :: KICE     = 1
   integer, parameter :: KLIQUID  = 2
   integer, parameter :: KRAIN    = 3
   integer, parameter :: KSNOW    = 4

   real    :: CO2

   real    :: TAUCRIT                       ! pressure separating low and middle clouds
   real    :: PRS_LOW_MID                   ! pressure separating low and middle clouds
   real    :: PRS_MID_HIGH                  ! pressure separating low and high clouds
   integer :: LCLDMH                        ! model level separating high and middle clouds
   integer :: LCLDLM                        ! model level separating low  and middle clouds

   integer :: NUM_BANDS                     ! Holds value from AGCM.rc
   integer :: TOTAL_RAD_BANDS

   character(len=ESMF_MAXSTR), pointer :: AEROSOLS(:)

   integer :: i, j, K, L, YY, DOY

   real, dimension (IM,JM)         :: T2M   !  fractional cover of sub-grid regions
   real, dimension (IM,JM,NS)      :: FS    !  fractional cover of sub-grid regions
   real, dimension (IM,JM,NS)      :: TG    !  land or ocean surface temperature
   real, dimension (IM,JM,NS,10)   :: EG    !  land or ocean surface emissivity
   real, dimension (IM,JM,NS)      :: TV    !  vegetation temperature
   real, dimension (IM,JM,NS,10)   :: EV    !  vegetation emissivity
   real, dimension (IM,JM,NS,10)   :: RV    !  vegetation reflectivity
   real, dimension (IM,JM,LM,4)    :: CWC   !  in-cloud cloud water mixing ratio
   real, dimension (IM,JM,LM,4)    :: REFF  !  effective radius of cloud particles
   real, dimension (IM,JM,LM,10)   :: TAUDIAG
   real, dimension (IM,JM,LM)      :: RH, PL

! Local Aerosol Variables
! -----------------------

   REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: TAUA
   REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: SSAA
   REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: ASYA

   REAL :: X
   INTEGER :: IB, NA

   integer, parameter :: NB_CHOU   = 10        ! Number of bands in IRRAD calcs for Chou
   integer, parameter :: NB_RRTMG  = 16        ! Number of bands in IRRAD calcs for RRTMG
   integer, parameter :: NB_RRTMGP = 16        ! Number of bands in IRRAD calcs for RRTMGP

   integer, parameter :: NB_CHOU_SORAD   = 8   ! Number of bands in SORAD calcs for Chou
   integer, parameter :: NB_RRTMG_SORAD  = 14  ! Number of bands in SORAD calcs for RRTMG
   integer, parameter :: NB_RRTMGP_SORAD = 14  ! Number of bands in SORAD calcs for RRTMGP

   integer            :: NB_IRRAD              ! Number of bands in IRRAD calcs

   INTEGER :: OFFSET

! AERO state variables
! --------------------
   type (ESMF_State)                    :: AERO
   type (ESMF_Field)                    :: AS_FIELD
   character(len=ESMF_MAXSTR)           :: AS_FIELD_NAME   
   type (ESMF_Field)                    :: AS_FIELD_Q
   integer                              :: AS_STATUS
   real, pointer,     dimension(:,:,:)  :: AS_PTR_3D
   real, pointer,     dimension(:,:,:)  :: AS_PTR_PLE
   real, pointer,     dimension(:,:,:)  :: AS_PTR_T
   real, pointer,     dimension(:,:,:)  :: AS_PTR_Q
   real, allocatable, dimension(:,:,:)  :: AS_ARR_RH
   real, allocatable, dimension(:,:,:)  :: AS_ARR_PL

   real, allocatable, dimension(:,:,:,:):: AEROSOL_EXT
   real, allocatable, dimension(:,:,:,:):: AEROSOL_SSA
   real, allocatable, dimension(:,:,:,:):: AEROSOL_ASY

   real, pointer,     dimension(:,:,:)  :: VAR_PTR_3D

   logical                              :: implements_aerosol_optics 

   integer                              :: band

! Variables for RRTMG Code
! ------------------------

   integer :: iceflglw        ! Flag for ice particle specification
   integer :: liqflglw        ! Flag for liquid droplet specification
   logical :: Ts_derivs       ! calculate Tsurf derivatives of upward fluxes
   integer :: NN, IJ, LV

   real,    allocatable, dimension(:,:)   :: FCLD_R
   real,    allocatable, dimension(:,:)   :: TLEV_R       ! Edge Level temperature
   real,    allocatable, dimension(:,:)   :: PLE_R        ! Reverse of level pressure
   real,    allocatable, dimension(:,:)   :: ZM_R         ! Reverse of layer height
   real,    allocatable, dimension(:,:)   :: EMISS        ! Surface emissivity at 16 RRTMG bands
   real,    allocatable, dimension(:,:)   :: CLIQWP       ! Cloud liquid water path 
   real,    allocatable, dimension(:,:)   :: CICEWP       ! Cloud ice water path 
   real,    allocatable, dimension(:,:)   :: RELIQ        ! Cloud liquid effective radius 
   real,    allocatable, dimension(:,:)   :: REICE        ! Cloud ice effective radius 
   real,    allocatable, dimension(:,:,:) :: TAUAER
   real,    allocatable, dimension(:,:)   :: PL_R, T_R,  Q_R, O2_R,  O3_R
   real,    allocatable, dimension(:,:)   :: CO2_R, CH4_R, N2O_R, CFC11_R, CFC12_R, CFC22_R, CCL4_R
   real,    allocatable, dimension(:)     :: TSFC
   real,    allocatable, dimension(:,:)   :: UFLX, DFLX, UFLXC, DFLXC, DUFLX_DTS, DUFLXC_DTS
   integer, allocatable, dimension(:,:)   :: CLEARCOUNTS
   real,    allocatable, dimension(:)     :: ALAT
   real,    allocatable, dimension(:,:)   :: OLRBRG, DOLRBRG_DTS

   ! pmn: should we update these?
   real, parameter :: O2   = 0.2090029E+00 ! preexisting
   real, parameter :: N2   = 0.7906400E+00 ! approx from rrtmgp input file
   real, parameter :: CCL4 = 0.1105000E-09 ! preexisting
   real, parameter :: CO   = 0.            ! currently zero

   ! variables for RRTMGP code
   ! -------------------------

   ! conversion factor (see below)
   real(wp), parameter :: cwp_fac = real(1000./MAPL_GRAV,kind=wp)

   ! input arrays: dimensions (col, lay)
   real(wp), dimension(:,:), allocatable         :: p_lay, t_lay, p_lev, dp_wp, cf_wp
   real(wp), dimension(:,:), allocatable, target :: t_lev

   ! surface input arrays
   real(wp), dimension(:),   allocatable         :: t_sfc
   real(wp), dimension(:,:), allocatable         :: emis_sfc ! first dim is band

   ! fluxes
   real(wp), dimension(:,:), allocatable, target :: flux_up_clrsky, flux_dn_clrsky, dfupdts_clrsky, &
                                                    flux_up_clrnoa, flux_dn_clrnoa, dfupdts_clrnoa, &
                                                    flux_up_allsky, flux_dn_allsky, dfupdts_allsky, &
                                                    flux_up_allnoa, flux_dn_allnoa, dfupdts_allnoa

   ! derived types for interacting with RRTMGP
   type(ty_gas_optics_rrtmgp), pointer           :: k_dist
   type(ty_gas_concs)                            :: gas_concs, gas_concs_subset
   type(ty_cloud_optics)                         :: cloud_optics
   type(ty_source_func_lw)                       :: sources
   type(ty_fluxes_broadband)                     :: fluxes_clrsky, fluxes_clrnoa, &
                                                    fluxes_allsky, fluxes_allnoa

   ! The band-space (ncol,nlay,nbnd) aerosol and in-cloud optical properties
   ! Polymorphic with dynamic type (#streams) defined later
   class(ty_optical_props_arry), allocatable :: &
     cloud_props, cloud_props_subset, &
       aer_props,   aer_props_subset
 
   ! The g-point cloud optical properties used for mcICA
   class(ty_optical_props_arry), allocatable :: cloud_props_gpt

   ! The g-point optical properties used in RT calculations for clean|dirty exports
   ! Polymorphic with dynamic type (#streams) defined later
   class(ty_optical_props_arry), allocatable :: clean_optical_props, dirty_optical_props

   ! RRTMGP locals
   logical :: top_at_1, u2s, partial_block
   logical :: need_dirty_optical_props, need_cloud_optical_props
   logical :: export_clrnoa, export_clrsky, export_allnoa, export_allsky
   logical ::   calc_clrnoa,   calc_clrsky,   calc_allnoa,   calc_allsky
   integer :: ncol, nbnd, ngpt, nmom, nga, icergh
   integer :: b, nBlocks, colS, colE, ncols_subset, &
              partial_blockSize, icol, isub
   integer :: iBeg, iEnd, jBeg, jEnd
   integer :: IM_World, JM_World, dims(3)
   character(len=ESMF_MAXPATHLEN) :: k_dist_file, cloud_optics_file
   character(len=ESMF_MAXSTR)     :: error_msg
   character(len=128)             :: cloud_optics_type, cloud_overlap_type
   real(wp), dimension(LM+1)      :: tlev_wp
   type (ESMF_Time)               :: ReferenceTime
   type (ESMF_TimeInterval)       :: RefreshInterval
!  real :: delTS_r
!  real(wp) delTS

   ! gridcolum presence of liq and ice clouds (ncol,nlay)
   real(wp), dimension(:,:), allocatable :: clwp, ciwp

   ! a column random number generator
#ifdef HAVE_MKL
   type(ty_rng_mklvsl_plus) :: rng
#endif
   integer, dimension(:), allocatable :: seeds

   ! uniform random numbers need by mcICA (ngpt,nlay,cols)
   real(wp), dimension(:,:,:), allocatable :: urand

   ! Cloud mask from overlap scheme (cols,nlay,ngpt)
   logical,  dimension(:,:,:), allocatable :: cld_mask

   ! TEMP ... see below
   real(wp) :: press_ref_min, ptop
   real(wp) ::  temp_ref_min, tmin
   real(wp), parameter :: ptop_increase_OK_fraction = 0.01_wp
   real(wp), parameter :: tmin_increase_OK_Kelvin   = 10.0_wp

   ! block size for column processing
   ! set for efficiency from resource file
   integer :: rrtmgp_blockSize

! For aerosol
   integer                    :: in
   real                       :: xx
   type (ESMF_Time)           :: CURRENTTIME
   real, dimension (LM+1)     :: TLEV
   real, dimension (LM)       :: DP

! pointers to import
!-------------------

   real, pointer, dimension(:    )   :: PREF
   real, pointer, dimension(:,:  )   :: TS
   real, pointer, dimension(:,:  )   :: EMIS
   real, pointer, dimension(:,:,:)   :: PLE, T,  Q,  O3
   real, pointer, dimension(:,:,:)   :: CH4, N2O, CFC11, CFC12, HCFC22
   real, pointer, dimension(:,:,:)   :: QL,  QI, QR, QS
   real, pointer, dimension(:,:,:)   :: RI,  RL, RR, RS, FCLD
   real, pointer, dimension(:,:,:,:) :: RAERO
   real, pointer, dimension(:,:,:)   :: QAERO

! pointers to exports
!--------------------

   real, pointer, dimension(:,:  )   :: CLDPRS
   real, pointer, dimension(:,:  )   :: CLDTMP
   real, pointer, dimension(:,:,:)   :: TAUIR
   real, pointer, dimension(:,:  )   :: CLDTTLW
   real, pointer, dimension(:,:  )   :: CLDHILW
   real, pointer, dimension(:,:  )   :: CLDMDLW
   real, pointer, dimension(:,:  )   :: CLDLOLW
   real, pointer, dimension(:,:  )   :: TSREFF
   real, pointer, dimension(:,:  )   :: SFCEM
   real, pointer, dimension(:,:  )   :: LWS0 
   real, pointer, dimension(:,:  )   :: DSFDTS

   ! for compact multi-export handling
   real, pointer, dimension(:,:  ) :: ptr2d
   real, pointer, dimension(:,:,:) :: ptr3d
   type S_
     character(len=:), allocatable :: str
   end type S_
   type(S_), allocatable :: list(:)

! helper for testing RRTMGP error status on return;
! allows line number reporting cf. original call method
#define TEST_(A) error_msg = A; if (trim(error_msg)/="") then; _ASSERT(.false.,"RRTMGP Error: "//trim(error_msg)); endif

#ifdef _CUDA
! MATMAT CUDA Variables
   type(dim3) :: Grid, Block
   integer :: blocksize
#endif
   integer :: PARTITION_SIZE

!  Begin...
!----------

   IAm = "Lw_Driver"
   call MAPL_TimerOn(MAPL,"--MISC")

! Pointer to Imports used only for full transfer calculation
!-----------------------------------------------------------
   
   call MAPL_GetPointer(IMPORT, PLE,    'PLE',    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT, T,      'T',      RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT, Q,      'QV',     RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT, QL,     'QL',     RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT, QI,     'QI',     RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT, QR,     'QR',     RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT, QS,     'QS',     RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT, RL,     'RL',     RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT, RI,     'RI',     RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT, RR,     'RR',     RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT, RS,     'RS',     RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT, O3,     'O3',     RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT, CH4,    'CH4',    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT, N2O,    'N2O',    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT, CFC11,  'CFC11',  RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT, CFC12,  'CFC12',  RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT, HCFC22, 'HCFC22', RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT, FCLD,   'FCLD',   RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT, EMIS,   'EMIS',   RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT, PREF,   'PREF',   RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT, TS,     'TS',     RC=STATUS); VERIFY_(STATUS)

   PL = 0.5*(PLE(:,:,:UBOUND(PLE,3)-1)+PLE(:,:,LBOUND(PLE,3)+1:))
   RH = Q/GEOS_QSAT(T,PL,PASCALS=.true.)

! Get trace gases concentrations by volume (pppv) from configuration
!-------------------------------------------------------------------

   call MAPL_GetResource( MAPL, CO2, 'CO2:', RC=STATUS)
   VERIFY_(STATUS)

   if(CO2<0.0) then
      call ESMF_ClockGet(CLOCK, currTIME=CURRENTTIME, RC=STATUS)
      VERIFY_(STATUS)
      call ESMF_TimeGet (CURRENTTIME, YY=YY, DayOfYear=DOY, RC=STATUS)
      VERIFY_(STATUS)
      CO2 = GETCO2(YY,DOY)
   endif

   call MAPL_GetResource( MAPL, PRS_LOW_MID,  'PRS_LOW_MID_CLOUDS:',  DEFAULT=70000., &
        RC=STATUS)
   VERIFY_(STATUS)

   call MAPL_GetResource( MAPL, PRS_MID_HIGH, 'PRS_MID_HIGH_CLOUDS:', DEFAULT=40000., &
        RC=STATUS)
   VERIFY_(STATUS)

! Prepare for aerosol optics calculations
!_---------------------------------------

   ! Set the offset for the IRRAD aerosol bands
   if (USE_RRTMGP_SORAD) then
      OFFSET = NB_RRTMGP_SORAD
   else if (USE_RRTMG_SORAD) then
      OFFSET = NB_RRTMG_SORAD 
   else
      OFFSET = NB_CHOU_SORAD
   end if

   ! Set number of IRRAD bands for aerosol optics
   if (USE_RRTMGP) then
      NB_IRRAD = NB_RRTMGP
   else if (USE_RRTMG) then
      NB_IRRAD = NB_RRTMG 
   else
      NB_IRRAD = NB_CHOU
   end if

! Test to see if AGCM.rc is set up correctly for the Radiation selected
!----------------------------------------------------------------------

   call MAPL_GetResource( MAPL, NUM_BANDS ,'NUM_BANDS:', RC=STATUS)
   VERIFY_(STATUS)

   TOTAL_RAD_BANDS = NB_IRRAD
   if (USE_RRTMGP_SORAD) then
     TOTAL_RAD_BANDS = TOTAL_RAD_BANDS + NB_RRTMGP_SORAD
   else if (USE_RRTMG_SORAD) then
     TOTAL_RAD_BANDS = TOTAL_RAD_BANDS + NB_RRTMG_SORAD
   else
     TOTAL_RAD_BANDS = TOTAL_RAD_BANDS + NB_CHOU_SORAD
   end if

   if (NUM_BANDS /= TOTAL_RAD_BANDS) then
      if (MAPL_am_I_Root()) then
         write (*,*) "NUM_BANDS is not set up correctly for the radiation combination selected:"
         write (*,*) "    IRRAD RRTMG: ", USE_RRTMGP      , USE_RRTMG      , USE_CHOU
         write (*,*) "    SORAD RRTMG: ", USE_RRTMGP_SORAD, USE_RRTMG_SORAD, USE_CHOU_SORAD
         write (*,*) "Please check that your optics tables and NUM_BANDS are correct."
      end if
      _ASSERT(.FALSE.,'needs informative message')
   end if

! Compute surface air temperature ("2 m") adiabatically
!------------------------------------------------------ 

   T2M = T(:,:,LM)*(0.5*(1.0 + PLE(:,:,LM-1)/PLE(:,:,LM)))**(-MAPL_KAPPA)

! For now, use the same emissivity for all bands
!-----------------------------------------------

   do K = 1, 10
      EG(:,:,1,K)   = EMIS(:,:)
   end do

! For now, hardwire vegetation and aerosol parameters
!----------------------------------------------------

   FS                  = 1.0
   TG(:,:,1)           = TS
   TV(:,:,1)           = TS
   EV                  = 0.0
   RV                  = 0.0

! Copy cloud constituent properties into contiguous buffers 
!----------------------------------------------------------

   CWC (:,:,:,KICE   ) = QI
   CWC (:,:,:,KLIQUID) = QL
   CWC (:,:,:,KRAIN  ) = QR
   CWC (:,:,:,KSNOW  ) = QS

   REFF(:,:,:,KICE   ) = RI * 1.0e6
   REFF(:,:,:,KLIQUID) = RL * 1.0e6
   REFF(:,:,:,KRAIN  ) = RR * 1.0e6
   REFF(:,:,:,KSNOW  ) = RS * 1.0e6

! Determine the model level seperating high and middle clouds
!------------------------------------------------------------

   LCLDMH = 1
   do K = 1, LM
      if( PREF(K) >= PRS_MID_HIGH ) then
         LCLDMH = K
         exit
      end if
   end do

! Determine the model level seperating low and middle clouds
!-----------------------------------------------------------

   LCLDLM = LM
   do K = 1, LM
      if( PREF(K) >= PRS_LOW_MID  ) then
         LCLDLM = K
         exit
      end if
   end do

   call MAPL_GetPointer(EXPORT, CLDTTLW, 'CLDTTLW', RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT, CLDHILW, 'CLDHILW', RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT, CLDMDLW, 'CLDMDLW', RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT, CLDLOLW, 'CLDLOLW', RC=STATUS); VERIFY_(STATUS)

   ! ------------------
   ! Begin aerosol code
   ! ------------------

   ! Allocate per-band aerosol arrays
   ! --------------------------------

   ALLOCATE(TAUA(IM,JM,LM,NB_IRRAD), STAT = STATUS)
   VERIFY_(STATUS)
   ALLOCATE(SSAA(IM,JM,LM,NB_IRRAD), STAT = STATUS)
   VERIFY_(STATUS)
   ALLOCATE(ASYA(IM,JM,LM,NB_IRRAD), STAT = STATUS)
   VERIFY_(STATUS)

   ! Zero out aerosol arrays. If NA == 0, these zeroes are then used inside IRRAD.
   ! -----------------------------------------------------------------------------
   NA    = 0

   TAUA  = 0.0
   SSAA  = 0.0
   ASYA  = 0.0

   ! If we have aerosols, accumulate the arrays
   ! ------------------------------------------

   call MAPL_TimerOn(MAPL,"---AEROSOLS")

   call ESMF_StateGet(IMPORT, 'AERO', AERO, RC=STATUS)
   VERIFY_(STATUS)

   call ESMF_AttributeGet(aero, name='implements_aerosol_optics_method', &
                                value=implements_aerosol_optics, RC=STATUS)
   VERIFY_(STATUS)

   RADIATIVELY_ACTIVE_AEROSOLS: if (implements_aerosol_optics) then

      ! set RH for aerosol optics
      call ESMF_AttributeGet(AERO, name='relative_humidity_for_aerosol_optics', value=AS_FIELD_NAME, RC=STATUS)
      VERIFY_(STATUS)

      if (AS_FIELD_NAME /= '') then
         call MAPL_GetPointer(AERO, AS_PTR_3D, trim(AS_FIELD_NAME), RC=STATUS)
         VERIFY_(STATUS)

         AS_PTR_3D = RH
      end if

      ! set PLE for aerosol optics
      call ESMF_AttributeGet(AERO, name='air_pressure_for_aerosol_optics', value=AS_FIELD_NAME, RC=STATUS)
      VERIFY_(STATUS)

      if (AS_FIELD_NAME /= '') then
         call MAPL_GetPointer(AERO, AS_PTR_3D, trim(AS_FIELD_NAME), RC=STATUS)
         VERIFY_(STATUS)
           
         AS_PTR_3D = PLE
      end if

      ! allocate memory for total aerosol ext, ssa and asy at all solar bands
      allocate(AEROSOL_EXT(IM,JM,LM,NB_IRRAD),  &
               AEROSOL_SSA(IM,JM,LM,NB_IRRAD),  &
               AEROSOL_ASY(IM,JM,LM,NB_IRRAD),  stat=STATUS)
      VERIFY_(STATUS)

      AEROSOL_EXT = 0.0
      AEROSOL_SSA = 0.0
      AEROSOL_ASY = 0.0

      ! compute aerosol optics at all solar bands
      IR_BANDS: do band = 1, NB_IRRAD
         call ESMF_AttributeSet(AERO, name='band_for_aerosol_optics', value=(OFFSET+band), RC=STATUS)
         VERIFY_(STATUS)

         ! execute the aero provider's optics method
         call ESMF_MethodExecute(AERO, label="aerosol_optics", RC=STATUS)
         VERIFY_(STATUS)

         ! EXT from AERO_PROVIDER
         call ESMF_AttributeGet(AERO, name='extinction_in_air_due_to_ambient_aerosol', value=AS_FIELD_NAME, RC=STATUS)
         VERIFY_(STATUS)

         if (AS_FIELD_NAME /= '') then
            call MAPL_GetPointer(AERO, AS_PTR_3D, trim(AS_FIELD_NAME),  RC=STATUS); VERIFY_(STATUS)

            if (associated(AS_PTR_3D)) then 
               AEROSOL_EXT(:,:,:,band) = AS_PTR_3D
            end if
         end if

         ! SSA from AERO_PROVIDER
         call ESMF_AttributeGet(AERO, name='single_scattering_albedo_of_ambient_aerosol', value=AS_FIELD_NAME, RC=STATUS)
         VERIFY_(STATUS)

         if (AS_FIELD_NAME /= '') then
            call MAPL_GetPointer(AERO, AS_PTR_3D, trim(AS_FIELD_NAME),  RC=STATUS); VERIFY_(STATUS)

            if (associated(AS_PTR_3D)) then 
               AEROSOL_SSA(:,:,:,band) = AS_PTR_3D
            end if
         end if

         ! ASY from AERO_PROVIDER
         call ESMF_AttributeGet(AERO, name='asymmetry_parameter_of_ambient_aerosol', value=AS_FIELD_NAME, RC=STATUS)
         VERIFY_(STATUS)

         if (AS_FIELD_NAME /= '') then
            call MAPL_GetPointer(AERO, AS_PTR_3D, trim(AS_FIELD_NAME),  RC=STATUS)
            VERIFY_(STATUS)

            if (associated(AS_PTR_3D)) then 
               AEROSOL_ASY(:,:,:,band) = AS_PTR_3D
            end if
         end if
      end do IR_BANDS

      NA = 3

      TAUA = AEROSOL_EXT
      SSAA = AEROSOL_SSA
      ASYA = AEROSOL_ASY

      deallocate(AEROSOL_EXT, __STAT__)
      deallocate(AEROSOL_SSA, __STAT__)
      deallocate(AEROSOL_ASY, __STAT__)

   end if RADIATIVELY_ACTIVE_AEROSOLS

   call MAPL_TimerOff(MAPL,"---AEROSOLS")

   call MAPL_TimerOff(MAPL,"--MISC")

   SCHEME: if (USE_CHOU) then

      call MAPL_TimerOn (MAPL,"--IRRAD",RC=STATUS)
      VERIFY_(STATUS)

#ifdef _CUDA

      call MAPL_GetResource(MAPL,BLOCKSIZE,'BLOCKSIZE:',DEFAULT=128,RC=STATUS)
      VERIFY_(STATUS)

      Block = dim3(blocksize,1,1)
      Grid = dim3(ceiling(real(IM*JM)/real(blocksize)),1,1)

      _ASSERT(LM <= GPU_MAXLEVS,'needs informative message') ! If this is tripped, ESMA_arch.mk must be edited.

      _ASSERT(NS == MAXNS,'needs informative message') ! If this is tripped, the local GNUmakefile
                           ! must be edited.

      call MAPL_TimerOn(MAPL,"---IRRAD_ALLOC",RC=STATUS)
      VERIFY_(STATUS)

      ! ----------------------
      ! Allocate device arrays
      ! ----------------------

      ! Inputs
      ! ------

      ALLOCATE(PLE_DEV(IM*JM,LM+1), STAT = STATUS)
      VERIFY_(STATUS)
      ALLOCATE(TA_DEV(IM*JM,LM), STAT = STATUS)
      VERIFY_(STATUS)
      ALLOCATE(WA_DEV(IM*JM,LM), STAT = STATUS)
      VERIFY_(STATUS)
      ALLOCATE(OA_DEV(IM*JM,LM), STAT = STATUS)
      VERIFY_(STATUS)
      ALLOCATE(TB_DEV(IM*JM), STAT = STATUS)
      VERIFY_(STATUS)
      ALLOCATE(N2O_DEV(IM*JM,LM), STAT = STATUS)
      VERIFY_(STATUS)
      ALLOCATE(CH4_DEV(IM*JM,LM), STAT = STATUS)
      VERIFY_(STATUS)
      ALLOCATE(CFC11_DEV(IM*JM,LM), STAT = STATUS)
      VERIFY_(STATUS)
      ALLOCATE(CFC12_DEV(IM*JM,LM), STAT = STATUS)
      VERIFY_(STATUS)
      ALLOCATE(CFC22_DEV(IM*JM,LM), STAT = STATUS)
      VERIFY_(STATUS)
      ALLOCATE(FS_DEV(IM*JM,NS), STAT = STATUS)
      VERIFY_(STATUS)
      ALLOCATE(TG_DEV(IM*JM,NS), STAT = STATUS)
      VERIFY_(STATUS)
      ALLOCATE(EG_DEV(IM*JM,NS,10), STAT = STATUS)
      VERIFY_(STATUS)
      ALLOCATE(TV_DEV(IM*JM,NS), STAT = STATUS)
      VERIFY_(STATUS)
      ALLOCATE(EV_DEV(IM*JM,NS,10), STAT = STATUS)
      VERIFY_(STATUS)
      ALLOCATE(RV_DEV(IM*JM,NS,10), STAT = STATUS)
      VERIFY_(STATUS)
      ALLOCATE(CWC_DEV(IM*JM,LM,4), STAT = STATUS)
      VERIFY_(STATUS)
      ALLOCATE(FCLD_DEV(IM*JM,LM), STAT = STATUS)
      VERIFY_(STATUS)
      ALLOCATE(REFF_DEV(IM*JM,LM,4), STAT=STATUS)
      VERIFY_(STATUS)

      ALLOCATE(TAUA_DEV(IM*JM,LM,NB_CHOU), STAT = STATUS)
      VERIFY_(STATUS)
      ALLOCATE(SSAA_DEV(IM*JM,LM,NB_CHOU), STAT = STATUS)
      VERIFY_(STATUS)
      ALLOCATE(ASYA_DEV(IM*JM,LM,NB_CHOU), STAT = STATUS)
      VERIFY_(STATUS)

      ! Constant arrays in global memory
      ! --------------------------------

      ALLOCATE(C1(NX,NC), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(C2(NX,NC), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(C3(NX,NC), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(OO1(NX,NO), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(OO2(NX,NO), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(OO3(NX,NO), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(H11(NX,NH), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(H12(NX,NH), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(H13(NX,NH), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(H21(NX,NH), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(H22(NX,NH), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(H23(NX,NH), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(H81(NX,NH), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(H82(NX,NH), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(H83(NX,NH), STAT=STATUS)
      VERIFY_(STATUS)

      ! Outputs
      ! -------

      ALLOCATE(FLXU_DEV(IM*JM,LM+1), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(FLXAU_DEV(IM*JM,LM+1), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(FLCU_DEV(IM*JM,LM+1), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(FLAU_DEV(IM*JM,LM+1), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(FLXD_DEV(IM*JM,LM+1), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(FLXAD_DEV(IM*JM,LM+1), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(FLCD_DEV(IM*JM,LM+1), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(FLAD_DEV(IM*JM,LM+1), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(DFDTS_DEV(IM*JM,LM+1), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(SFCEM_DEV(IM*JM), STAT=STATUS)
      VERIFY_(STATUS)
      ALLOCATE(TAUDIAG_DEV(IM*JM,LM,10), STAT=STATUS)
      VERIFY_(STATUS)

      call MAPL_TimerOff(MAPL,"---IRRAD_ALLOC",RC=STATUS)
      VERIFY_(STATUS)

      call MAPL_TimerOn(MAPL,"---IRRAD_DATA",RC=STATUS)
      VERIFY_(STATUS)

      call MAPL_TimerOn(MAPL,"----IRRAD_DATA_DEVICE",RC=STATUS)
      VERIFY_(STATUS)

      ! --------------------------
      ! Copy host arrays to device
      ! --------------------------

      ! Inputs
      ! ------

      STATUS = cudaMemcpy(PLE_DEV,PLE,IM*JM*(LM+1))
      VERIFY_(STATUS)
      STATUS = cudaMemcpy(TA_DEV,T,IM*JM*LM)
      VERIFY_(STATUS)
      STATUS = cudaMemcpy(WA_DEV,Q,IM*JM*LM)
      VERIFY_(STATUS)
      STATUS = cudaMemcpy(OA_DEV,O3,IM*JM*LM)
      VERIFY_(STATUS)
      STATUS = cudaMemcpy(TB_DEV,T2M,IM*JM)
      VERIFY_(STATUS)
      STATUS = cudaMemcpy(N2O_DEV,N2O,IM*JM*LM) 
      VERIFY_(STATUS)
      STATUS = cudaMemcpy(CH4_DEV,CH4,IM*JM*LM) 
      VERIFY_(STATUS)
      STATUS = cudaMemcpy(CFC11_DEV,CFC11,IM*JM*LM) 
      VERIFY_(STATUS)
      STATUS = cudaMemcpy(CFC12_DEV,CFC12,IM*JM*LM) 
      VERIFY_(STATUS)
      STATUS = cudaMemcpy(CFC22_DEV,HCFC22,IM*JM*LM) 
      VERIFY_(STATUS)
      STATUS = cudaMemcpy(FS_DEV,FS,IM*JM*NS) 
      VERIFY_(STATUS)
      STATUS = cudaMemcpy(TG_DEV,TG,IM*JM*NS) 
      VERIFY_(STATUS)
      STATUS = cudaMemcpy(EG_DEV,EG,IM*JM*NS*10) 
      VERIFY_(STATUS)
      STATUS = cudaMemcpy(TV_DEV,TV,IM*JM*NS) 
      VERIFY_(STATUS)
      STATUS = cudaMemcpy(EV_DEV,EV,IM*JM*NS*10) 
      VERIFY_(STATUS)
      STATUS = cudaMemcpy(RV_DEV,RV,IM*JM*NS*10) 
      VERIFY_(STATUS)
      STATUS = cudaMemcpy(CWC_DEV,CWC,IM*JM*LM*4) 
      VERIFY_(STATUS)
      STATUS = cudaMemcpy(FCLD_DEV,FCLD,IM*JM*LM) 
      VERIFY_(STATUS)
      STATUS = cudaMemcpy(REFF_DEV,REFF,IM*JM*LM*4) 
      VERIFY_(STATUS)

      STATUS = cudaMemcpy(TAUA_DEV,TAUA,IM*JM*LM*NB_CHOU)
      VERIFY_(STATUS)
      STATUS = cudaMemcpy(SSAA_DEV,SSAA,IM*JM*LM*NB_CHOU)
      VERIFY_(STATUS)
      STATUS = cudaMemcpy(ASYA_DEV,ASYA,IM*JM*LM*NB_CHOU)
      VERIFY_(STATUS)

      ! ---------------------------------------
      ! Copy Constant Arrays into Global Memory
      ! ---------------------------------------

      C1 = C1_CONST
      C2 = C2_CONST
      C3 = C3_CONST
      OO1 = OO1_CONST
      OO2 = OO2_CONST
      OO3 = OO3_CONST
      H11 = H11_CONST
      H12 = H12_CONST
      H13 = H13_CONST
      H21 = H21_CONST
      H22 = H22_CONST
      H23 = H23_CONST
      H81 = H81_CONST
      H82 = H82_CONST
      H83 = H83_CONST

      call MAPL_TimerOff(MAPL,"----IRRAD_DATA_DEVICE",RC=STATUS)
      VERIFY_(STATUS)

      call MAPL_TimerOn(MAPL,"----IRRAD_DATA_CONST",RC=STATUS)
      VERIFY_(STATUS)

      ! --------------
      ! Copy Constants
      ! --------------

      XKW = XKW_CONST
      XKE = XKE_CONST
      MW = MW_CONST
      AW = AW_CONST
      BW = BW_CONST
      PM = PM_CONST
      FKW = FKW_CONST
      GKW = GKW_CONST
      AIB_IR = AIB_IR_CONST
      AWB_IR = AWB_IR_CONST
      AIW_IR = AIW_IR_CONST
      AWW_IR = AWW_IR_CONST
      AIG_IR = AIG_IR_CONST
      AWG_IR = AWG_IR_CONST
      CB = CB_CONST
      DCB = DCB_CONST
      W11 = W11_CONST
      W12 = W12_CONST
      W13 = W13_CONST
      P11 = P11_CONST
      P12 = P12_CONST
      P13 = P13_CONST
      DWE = DWE_CONST
      DPE = DPE_CONST

      call MAPL_TimerOff(MAPL,"----IRRAD_DATA_CONST",RC=STATUS)
      VERIFY_(STATUS)

      call MAPL_TimerOff(MAPL,"---IRRAD_DATA",RC=STATUS)
      VERIFY_(STATUS)

      ! Do longwave calculations on a list of soundings
      !  This fills the internal state
      !------------------------------------------------
      ! Note: IRRAD wants all species in mole fraction
      ! except O3, which must be in mass mixing ratio.
      !------------------------------------------------

      call MAPL_TimerOn(MAPL,"---IRRAD_RUN",RC=STATUS)
      VERIFY_(STATUS)

      call irrad<<<Grid, Block>>>(IM*JM,LM,CO2,TRACE,NS,NA,NB_CHOU,LCLDMH,LCLDLM)

      STATUS = cudaGetLastError()
      if (STATUS /= 0) then
         write (*,*) "Error code from IRRAD kernel call: ", STATUS
         write (*,*) "Kernel call failed: ", cudaGetErrorString(STATUS)
         _ASSERT(.FALSE.,'needs informative message')
      end if

      call MAPL_TimerOff(MAPL,"---IRRAD_RUN",RC=STATUS)
      VERIFY_(STATUS)

      call MAPL_TimerOn(MAPL,"---IRRAD_DATA",RC=STATUS)
      VERIFY_(STATUS)

      call MAPL_TimerOn(MAPL,"----IRRAD_DATA_DEVICE",RC=STATUS)
      VERIFY_(STATUS)

      ! Copy outputs from device
      ! ------------------------

      STATUS = cudaMemcpy(FLXU_INT,FLXU_DEV,IM*JM*(LM+1))
      VERIFY_(STATUS)
      STATUS = cudaMemcpy(FLXAU_INT,FLXAU_DEV,IM*JM*(LM+1))
      VERIFY_(STATUS)
      STATUS = cudaMemcpy(FLCU_INT,FLCU_DEV,IM*JM*(LM+1))
      VERIFY_(STATUS)
      STATUS = cudaMemcpy(FLAU_INT,FLAU_DEV,IM*JM*(LM+1))
      VERIFY_(STATUS)

      STATUS = cudaMemcpy(FLXD_INT,FLXD_DEV,IM*JM*(LM+1))
      VERIFY_(STATUS)
      STATUS = cudaMemcpy(FLXAD_INT,FLXAD_DEV,IM*JM*(LM+1))
      VERIFY_(STATUS)
      STATUS = cudaMemcpy(FLCD_INT,FLCD_DEV,IM*JM*(LM+1))
      VERIFY_(STATUS)
      STATUS = cudaMemcpy(FLAD_INT,FLAD_DEV,IM*JM*(LM+1))
      VERIFY_(STATUS)

      STATUS = cudaMemcpy(DFDTS,DFDTS_DEV,IM*JM*(LM+1))
      VERIFY_(STATUS)
      STATUS = cudaMemcpy(SFCEM_INT,SFCEM_DEV,IM*JM)
      VERIFY_(STATUS)
      STATUS = cudaMemcpy(TAUDIAG,TAUDIAG_DEV,IM*JM*LM*10)
      VERIFY_(STATUS)

      call MAPL_TimerOff(MAPL,"----IRRAD_DATA_DEVICE",RC=STATUS)
      VERIFY_(STATUS)

      call MAPL_TimerOff(MAPL,"---IRRAD_DATA",RC=STATUS)
      VERIFY_(STATUS)

      call MAPL_TimerOn(MAPL,"---IRRAD_DEALLOC",RC=STATUS)
      VERIFY_(STATUS)

      ! ------------------------
      ! Deallocate device arrays
      ! ------------------------

      ! Inputs
      ! ------

      DEALLOCATE(PLE_DEV)
      DEALLOCATE(TA_DEV)
      DEALLOCATE(WA_DEV)
      DEALLOCATE(OA_DEV)
      DEALLOCATE(TB_DEV)
      DEALLOCATE(N2O_DEV)
      DEALLOCATE(CH4_DEV)
      DEALLOCATE(CFC11_DEV)
      DEALLOCATE(CFC12_DEV)
      DEALLOCATE(CFC22_DEV)
      DEALLOCATE(FS_DEV)
      DEALLOCATE(TG_DEV)
      DEALLOCATE(EG_DEV)
      DEALLOCATE(TV_DEV)
      DEALLOCATE(EV_DEV)
      DEALLOCATE(RV_DEV)
      DEALLOCATE(CWC_DEV)
      DEALLOCATE(FCLD_DEV)
      DEALLOCATE(REFF_DEV)

      DEALLOCATE(TAUA_DEV)
      DEALLOCATE(SSAA_DEV)
      DEALLOCATE(ASYA_DEV)

      ! Constant arrays in global memory
      ! --------------------------------

      DEALLOCATE(C1)
      DEALLOCATE(C2)
      DEALLOCATE(C3)
      DEALLOCATE(OO1)
      DEALLOCATE(OO2)
      DEALLOCATE(OO3)
      DEALLOCATE(H11)
      DEALLOCATE(H12)
      DEALLOCATE(H13)
      DEALLOCATE(H21)
      DEALLOCATE(H22)
      DEALLOCATE(H23)
      DEALLOCATE(H81)
      DEALLOCATE(H82)
      DEALLOCATE(H83)

      ! Outputs
      ! -------

      DEALLOCATE(FLXU_DEV)
      DEALLOCATE(FLXAU_DEV)
      DEALLOCATE(FLCU_DEV)
      DEALLOCATE(FLAU_DEV)
      DEALLOCATE(FLXD_DEV)
      DEALLOCATE(FLXAD_DEV)
      DEALLOCATE(FLCD_DEV)
      DEALLOCATE(FLAD_DEV)
      DEALLOCATE(DFDTS_DEV)
      DEALLOCATE(SFCEM_DEV)
      DEALLOCATE(TAUDIAG_DEV)

      call MAPL_TimerOff(MAPL,"---IRRAD_DEALLOC",RC=STATUS)
      VERIFY_(STATUS)

#else
! Do longwave calculations on a list of soundings
!  This fills the internal state
!------------------------------------------------
! Note: IRRAD wants all species in mole fraction
! except O3, which must be in mass mixing ratio.
!------------------------------------------------

      call MAPL_TimerOn(MAPL,"---IRRAD_RUN",RC=STATUS)
      VERIFY_(STATUS)

      call IRRAD( IM*JM, LM,       PLE,                           &
       T,        Q,      O3,    T2M,    CO2,                      &
       TRACE,    N2O,   CH4,    CFC11,     CFC12, HCFC22,         &
       CWC,    FCLD,  LCLDMH, LCLDLM,    REFF,                    &
       NS,       FS,     TG,    EG,     TV,        EV,    RV,     &
       NA, NB_CHOU, TAUA, SSAA, ASYA,                             &
       FLXU_INT,  FLCU_INT, FLAU_INT, FLXAU_INT,                  &
       FLXD_INT,  FLCD_INT, FLAD_INT, FLXAD_INT,                  &
       DFDTS, SFCEM_INT, TAUDIAG                                  )

      call MAPL_TimerOff(MAPL,"---IRRAD_RUN",RC=STATUS)
      VERIFY_(STATUS)

#endif

      ! pmn:
      ! Chou-Suarez does not provide these derivatives
      ! so clear is set to zero, no-aerosol to aerosol
      DFDTSC = 0.
      DFDTSNA  = DFDTS
      DFDTSCNA = DFDTSC

      call MAPL_TimerOff(MAPL,"--IRRAD",RC=STATUS)
      VERIFY_(STATUS)

   else if (USE_RRTMGP) then

      call MAPL_TimerOn(MAPL,"--RRTMGP",__RC__)

      ! columns are independent so collapse horizontal to 1D
      ncol = IM*JM

      call MAPL_TimerOn(MAPL,"---RRTMGP_SETUP_1",__RC__)

      ! absorbing gas names
      error_msg = gas_concs%init([character(3) :: &
        'h2o','co2','o3','n2o','co','ch4','o2','n2'])
      TEST_(error_msg)

      ! load gas concentrations (volume mixing ratios)
      ! "constant" gases
      TEST_(gas_concs%set_vmr('n2' , real(N2 ,kind=wp)))
      TEST_(gas_concs%set_vmr('o2' , real(O2 ,kind=wp)))
      TEST_(gas_concs%set_vmr('co2', real(CO2,kind=wp)))
      TEST_(gas_concs%set_vmr('co' , real(CO ,kind=wp)))
      ! variable gases
      ! (ozone converted from mass mixing ratio, water vapor from specific humidity)
      TEST_(gas_concs%set_vmr('ch4', real(reshape(CH4                             ,(/ncol,LM/)),kind=wp)))
      TEST_(gas_concs%set_vmr('n2o', real(reshape(N2O                             ,(/ncol,LM/)),kind=wp)))
      TEST_(gas_concs%set_vmr('o3' , real(reshape(O3      *(MAPL_AIRMW/MAPL_O3MW ),(/ncol,LM/)),kind=wp)))
      TEST_(gas_concs%set_vmr('h2o', real(reshape(Q/(1.-Q)*(MAPL_AIRMW/MAPL_H2OMW),(/ncol,LM/)),kind=wp)))

      call MAPL_TimerOff(MAPL,"---RRTMGP_SETUP_1",__RC__)

      call MAPL_TimerOn(MAPL,"---RRTMGP_IO_1",__RC__)

      ! access RRTMGP internal state from the GC
      call ESMF_UserCompGetInternalState(GC, 'RRTMGP_state', wrap, status)
      VERIFY_(status)
      rrtmgp_state => wrap%ptr

      ! initialize k-distribution if not already done
      call MAPL_GetResource( &
        MAPL, k_dist_file, "RRTMGP_DATA_LW:", &
        DEFAULT='rrtmgp-data-lw.nc',__RC__)
      if (.not. rrtmgp_state%initialized) then
        ! gas_concs needed only to access required gas names
        call load_and_init( &
          rrtmgp_state%k_dist, trim(k_dist_file), gas_concs)
        if (.not. rrtmgp_state%k_dist%source_is_internal()) then
          TEST_("RRTMGP-LW: does not seem to be LW")
        endif
        rrtmgp_state%initialized = .true.
      endif

      ! access by shorter name
      k_dist => rrtmgp_state%k_dist

      call MAPL_TimerOff(MAPL,"---RRTMGP_IO_1",__RC__)

      call MAPL_TimerOn(MAPL,"---RRTMGP_SETUP_2",__RC__)

      ! spectral dimensions
      ngpt = k_dist%get_ngpt()
      nbnd = k_dist%get_nband()
      _ASSERT(nbnd == NB_RRTMGP, 'RRTMGP-LW: expected different number of bands')

      ! note: use k_dist%get_band_lims_wavenumber()
      !   to access band limits.
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! pmn: check if same bands as RRTMG --- need code solution to make more flexible !!!!
      ! from RRTMG:
      ! wavenum1(:) = (/ 10., 350., 500., 630., 700., 820.,  980., 1080., 1180., 1390., 1480., 1800., 2080., 2250., 2380., 2600./)
      ! wavenum2(:) = (/350., 500., 630., 700., 820., 980., 1080., 1180., 1390., 1480., 1800., 2080., 2250., 2380., 2600., 3250./)
      ! from RRTMGP:
      !   write(*,*) 'band_lims_wvn(2,nbnd):', k_dist%get_band_lims_wavenumber()
      ! output reordered as above
      !                  10., 250., 500., 630., 700., 820.,  980., 1080., 1180., 1390., 1480., 1800., 2080., 2250., 2390., 2680.
      !                 250., 500., 630., 700., 820., 980., 1080., 1180., 1390., 1480., 1800., 2080., 2250., 2390., 2680., 3250.
      ! clearly there are some differences (250, 2390, 2680) ... have redone aerosol tables
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      ! allocate input arrays
      allocate(t_sfc(ncol), emis_sfc(nbnd,ncol), __STAT__)
      allocate(p_lay(ncol,LM), t_lay(ncol,LM), dp_wp(ncol,LM), cf_wp(ncol,LM), &
               p_lev(ncol,LM+1), t_lev(ncol,LM+1), __STAT__)

      ! load input arrays ...
      ! surface properties
      ! for now, use the same emissivity for all bands
      t_sfc    = real(       reshape(TS  ,(/ncol/))        ,kind=wp)
      emis_sfc = real(spread(reshape(EMIS,(/ncol/)),1,nbnd),kind=wp)

      ! basic profiles
      p_lay = real(reshape(PL  ,(/ncol,LM  /)), kind=wp)
      t_lay = real(reshape(T   ,(/ncol,LM  /)), kind=wp)
      p_lev = real(reshape(PLE ,(/ncol,LM+1/)), kind=wp)
      cf_wp = real(reshape(FCLD,(/ncol,LM  /)), kind=wp)

      ! RRTMGP's rte_lw takes a vertical ordering flag
      ! (no need to flip columns as with RRTMG)
      top_at_1 = p_lay(1, 1) < p_lay(1, LM)
      _ASSERT(top_at_1, 'unexpected vertical ordering')

      ! pmn: pressure KLUGE 
      ! Because currently k_dist%press_ref_min ~ 1.005 > GEOS-5 ptop of 1.0 Pa.
      ! Find better solution, perhaps getting AER to add a higher top.
      press_ref_min = k_dist%get_press_min()
      ptop = minval(p_lev(:,1))
      if (press_ref_min > ptop) then
        ! allow a small increase of ptop
        if (press_ref_min - ptop <= ptop * ptop_increase_OK_fraction) then
          where (p_lev(:,1) < press_ref_min) p_lev(:,1) = press_ref_min
          ! make sure no pressure ordering issues were created
          _ASSERT(all(p_lev(:,1) < p_lay(:,1)), 'pressure kluge causes misordering')
        else
          write(*,*) ' A ', ptop_increase_OK_fraction, &
                       ' fractional increase of ptop was insufficient'
          write(*,*) ' RRTMGP, GEOS-5 top (Pa)', press_ref_min, ptop
          TEST_('Model top too high for RRTMGP')
        endif
      endif

      ! pmn: temperature KLUGE 
      ! Currently k_dist%temp_ref_min = 160K but GEOS-5 has a global minimum
      ! temperature below this occasionally (< 1% of time). (The lowest temp
      ! seen is so far 151K). Consequently we will limit min(t_lay) to 160K.
      ! Find better solution, perhaps getting AER to produce a table with a
      ! lower minimum temperature.
      ! note: add 0.01K to lower limit so that t_lev calculated below will
      !   not fall below k_dist%get_temp_min() due to roundoff issues.
      temp_ref_min = k_dist%get_temp_min() + 0.01_wp
      tmin = minval(t_lay)
      if (temp_ref_min > tmin) then
        ! allow a small increase of tmin
        if (temp_ref_min - tmin <= tmin_increase_OK_Kelvin) then
          where (t_lay < temp_ref_min) t_lay = temp_ref_min
        else
          write(*,*) ' A ', tmin_increase_OK_Kelvin, &
                       'K increase of tmin was insufficient'
          write(*,*) ' RRTMGP, GEOS-5 t_min (K)', temp_ref_min, tmin
          TEST_('Found excessively cold model temperature for RRTMGP')
        endif
      endif

      ! calculate interface temperatures
      ! pmn: note that tlev is an optional argument of gas_optics(), and if not
      !   provided, it will supply its own internally. We could try running
      !   with this latter option to see what difference it makes.
      IJ = 0
      do J=1,JM
        do I=1,IM
          IJ = IJ + 1
          ! pmn: better to use un-kluged top-pressure
          !   for interpolation to tlev and for optical properties later
          ! dp_wp(IJ,1) = p_lev(IJ,2) - p_lev(IJ,1)
          dp_wp(IJ,1) = real(PLE(I,J,1) - PLE(I,J,0), kind=wp)
          do K = 2, LM
            ! interpolation between neighboring t_lay
            dp_wp(IJ,K) = p_lev(IJ,K+1) - p_lev(IJ,K)
            TLEV_wp(K) = (t_lay(IJ,K-1)*dp_wp(IJ,K) + t_lay(IJ,K)*dp_wp(IJ,K-1)) &
                           / (dp_wp(IJ,K-1) + dp_wp(IJ,K))
          enddo
          TLEV_wp(LM+1) = real(T2M(I,J), kind=wp) ! ~surface air temperature
          TLEV_wp(   1) = TLEV_wp(2)              ! assume isotropic at TOA
          t_lev(IJ,:) = TLEV_wp
        enddo
      enddo

      call MAPL_TimerOff(MAPL,"---RRTMGP_SETUP_2",__RC__)

      call MAPL_TimerOn(MAPL,"---RRTMGP_SETUP_3",__RC__)

      ! =================================================================
      ! for efficiency sake, we try to calculate only what we export ...
      ! =================================================================

      ! this line temporarily needed because of compiler bug
      allocate(list(1)); list(1) = S_('dummy')

      ! are clear clean exports requested?
      export_clrnoa = .false.
      list = [S_('FLA'), S_('FLAD'), S_('FLAU')]
      do i = 1, size(list)
        call MAPL_GetPointer(EXPORT, ptr3d, list(i)%str, __RC__)
        export_clrnoa = (export_clrnoa .or. associated(ptr3d))
      end do
      list = [S_('OLA'), S_('FLNSA'), S_('LAS')]
      do i = 1, size(list)
        call MAPL_GetPointer(EXPORT, ptr2d, list(i)%str, __RC__)
        export_clrnoa = (export_clrnoa .or. associated(ptr2d))
      end do

      ! are clear dirty exports requested?
      export_clrsky = .false.
      list = [S_('FLC'), S_('FLCD'), S_('FLCU')]
      do i = 1, size(list)
        call MAPL_GetPointer(EXPORT, ptr3d, list(i)%str, __RC__)
        export_clrsky = (export_clrsky .or. associated(ptr3d))
      end do
      list = [S_('OLC'), S_('OLCC5'), S_('FLNSC'), S_('LCS'), S_('LCSC5')]
      do i = 1, size(list)
        call MAPL_GetPointer(EXPORT, ptr2d, list(i)%str, __RC__)
        export_clrsky = (export_clrsky .or. associated(ptr2d))
      end do

      ! are cloudy clean exports requested?
      export_allnoa = .false.
      list = [S_('FLXA'), S_('FLXAD'), S_('FLXAU')]
      do i = 1, size(list)
        call MAPL_GetPointer(EXPORT, ptr3d, list(i)%str, __RC__)
        export_allnoa = (export_allnoa .or. associated(ptr3d))
      end do
      list = [S_('OLRA'), S_('FLNSNA'), S_('LWSA')]
      do i = 1, size(list)
        call MAPL_GetPointer(EXPORT, ptr2d, list(i)%str, __RC__)
        export_allnoa = (export_allnoa .or. associated(ptr2d))
      end do

      ! are cloudy dirty exports requested?
      export_allsky = .false.
      list = [S_('FLX'), S_('FLXD'), S_('FLXU')]
      do i = 1, size(list)
        call MAPL_GetPointer(EXPORT, ptr3d, list(i)%str, __RC__)
        export_allsky = (export_allsky .or. associated(ptr3d))
      end do
      list = [S_('OLR'), S_('SFCEM'), S_('FLNS'), S_('LWS')]
      do i = 1, size(list)
        call MAPL_GetPointer(EXPORT, ptr2d, list(i)%str, __RC__)
        export_allsky = (export_allsky .or. associated(ptr2d))
      end do
      deallocate(list,__STAT__)

      ! which fluxes to calculate?
      ! the clean fluxes are also used for "dirty" fluxes if no aerosols
      calc_clrnoa = (export_clrnoa .or. (export_clrsky .and. .not.implements_aerosol_optics))
      calc_allnoa = (export_allnoa .or. (export_allsky .and. .not.implements_aerosol_optics))
      calc_clrsky =  export_clrsky
      calc_allsky =  export_allsky

      ! do we actually need dirty optical properties?
      need_dirty_optical_props = &
        (export_clrsky .or. export_allsky) .and. implements_aerosol_optics

      ! do we need cloudy optical properties?
      need_cloud_optical_props = (export_allnoa .or. export_allsky)

      ! allocation of output arrays
      if (calc_clrnoa) then
        allocate(flux_up_clrnoa(ncol,LM+1), &
                 flux_dn_clrnoa(ncol,LM+1), &
                 dfupdts_clrnoa(ncol,LM+1), __STAT__)
      end if
      if (calc_allnoa) then
        allocate(flux_up_allnoa(ncol,LM+1), &
                 flux_dn_allnoa(ncol,LM+1), &
                 dfupdts_allnoa(ncol,LM+1), __STAT__)
      end if
      if (calc_clrsky) then
        allocate(flux_up_clrsky(ncol,LM+1), &
                 flux_dn_clrsky(ncol,LM+1), &
                 dfupdts_clrsky(ncol,LM+1), __STAT__)
      end if
      if (calc_allsky) then
        allocate(flux_up_allsky(ncol,LM+1), &
                 flux_dn_allsky(ncol,LM+1), &
                 dfupdts_allsky(ncol,LM+1), __STAT__)
      end if

      ! =======================================================================================
      ! IMPORTANT: Specify the type (#streams) of the LW RT calculations in clean_optical_props
      ! =======================================================================================
      ! While the aerosol system currently provides two-stream properties, as do the cloud
      ! optics files, we may choose any number of streams for the actual RT calculations by
      ! the appropriate instantiation of clean_optical_props here. The increment() statements
      ! below implicitly convert all component optical properties to this number of streams.
      ! Everything else in the code should adapt polymorphically without modification.
      ! options are: 1scl (no scattering), 2str (2-stream), or nstr (n-stream)
      ! For 1scl, must also specify the number of Gauss angles (nga) below.
      ! For nstr, must also specify the number of phase function moments (nmom) below.
      ! After Feb2020 update:
      ! Even for optical_props_2str, the default rte method is to use rescaled LW transport
      !   to account for scattering (in which case nga is used). To force explicit 2-stream
      !   scattering, must select u2s = .true. and allocate optical_props_2str below.
      ! =======================================================================================

      ! instantiate clean_optical_props with desired streams
      allocate(ty_optical_props_2str::clean_optical_props,__STAT__)

      ! default values
      nga  = 1 ! Used if 1scl or (2str but .not. u2s), in which case must be >= 1
      nmom = 2 ! Used only if nstr, in which case must be >= 2
      u2s = .false. ! forces explicit 2-stream scattering if optical_props_2str

      ! allow user selection of nga and u2s as appropriate
      select type(clean_optical_props)
        class is (ty_optical_props_1scl)
          call MAPL_GetResource( &
            MAPL, nga ,'RRTMGP_LW_N_GAUSS_ANGLES:', DEFAULT=nga, __RC__)
        class is (ty_optical_props_2str)
          call MAPL_GetResource( &
            MAPL, nga ,'RRTMGP_LW_N_GAUSS_ANGLES:', DEFAULT=nga, __RC__)
          call MAPL_GetResource( &
            MAPL, u2s ,'RRTMGP_LW_USE_2STREAM:',    DEFAULT=u2s, __RC__)
          _ASSERT(.not.u2s,'lw_solver_2stream() does not currently support Jacobians')
      end select

      ! the dirty_optical_props have the same number of streams
      if (need_dirty_optical_props) then
        select type(clean_optical_props)
          class is (ty_optical_props_1scl) ! No scattering
            allocate(ty_optical_props_1scl::dirty_optical_props,__STAT__)
          class is (ty_optical_props_2str)
            allocate(ty_optical_props_2str::dirty_optical_props,__STAT__)
          class is (ty_optical_props_nstr)
            allocate(ty_optical_props_nstr::dirty_optical_props,__STAT__)
        end select
      end if

      ! initialize clean_optical_props and sources
      ! dirty_optical_props are copy initialized later if needed
      TEST_(clean_optical_props%init(k_dist))
      TEST_(sources%init(k_dist))

      call MAPL_TimerOff(MAPL,"---RRTMGP_SETUP_3",__RC__)

 ! RRTMGP Jacobian currently uses fixed internal delTS of 1K
 !    ! numerical derivatives wrt surface temperature use a small delta TS
 !    ! (set this to zero exactly to disable linearization for RRTMGP)
 !    call MAPL_GetResource( MAPL, &
 !      delTS_r, "RRTMGP_NUMERICAL_DFDTS_DELTS_IN_K:", DEFAULT=0.1, __RC__)
 !    delTS = real(delTS_r, kind=wp)

      ! get cloud optical properties (band-only)
      if (need_cloud_optical_props) then
        ! pmn: some of this should be done only once per run

        call MAPL_TimerOn(MAPL,"---RRTMGP_IO_2",__RC__)

        ! load and init cloud_optics from file
        call MAPL_GetResource( &
          MAPL, cloud_optics_file, "RRTMGP_CLOUD_OPTICS_COEFFS_LW:", &
          DEFAULT='rrtmgp-cloud-optics-coeffs-lw.nc', __RC__)
        call MAPL_GetResource( &
          MAPL, cloud_optics_type, "RRTMGP_CLOUD_OPTICS_TYPE_LW:", &
          DEFAULT='LUT', __RC__)
        if (trim(cloud_optics_type)=='LUT') then
          call load_cld_lutcoeff (cloud_optics, cloud_optics_file)
        elseif (trim(cloud_optics_type)=='PADE') then
          call load_cld_padecoeff(cloud_optics, cloud_optics_file)
        else
          TEST_('unknown cloud_optics_type: '//trim(cloud_optics_file))
        end if

        call MAPL_TimerOff(MAPL,"---RRTMGP_IO_2",__RC__)

        ! ice surface roughness category for Yang (2013) ice optics
        ! icergh: 1 = none, 2 = medium, 3 = high
        call MAPL_GetResource( &
          MAPL, icergh, "RRTMGP_ICE_ROUGHNESS_LW:", &
          DEFAULT=2, __RC__)
        TEST_(cloud_optics%set_ice_roughness(icergh))

        call MAPL_TimerOn(MAPL,"---RRTMGP_SETUP_4",__RC__)

        ! cloud optics file is currently two-stream
        ! increment() will handle appropriate stream conversions
        allocate(ty_optical_props_2str::cloud_props,__STAT__)
        select type (cloud_props)
          class is (ty_optical_props_2str)
            ! band-only initialization
            TEST_(cloud_props%alloc_2str(ncol, LM, k_dist%get_band_lims_wavenumber()))
            ! allocate a subset for blocking use
            allocate(ty_optical_props_2str::cloud_props_subset,__STAT__)
            ! allocate a gpt version for cloud sampling
            allocate(ty_optical_props_2str::cloud_props_gpt,__STAT__)
          class default
            TEST_('cloud_props hardwired 2-stream for now')
        end select

        call MAPL_TimerOff(MAPL,"---RRTMGP_SETUP_4",__RC__)

        call MAPL_TimerOn(MAPL,"---RRTMGP_CLOUD_OPTICS",__RC__)

        ! make band in-cloud optical properties from cloud_optics
        allocate(clwp(ncol,LM), ciwp(ncol,LM),__STAT__)
        clwp = real(reshape(CWC(:,:,:,KLIQUID),(/ncol,LM/)),kind=wp) * dp_wp * cwp_fac ! in-cloud [g/m2]
        ciwp = real(reshape(CWC(:,:,:,KICE   ),(/ncol,LM/)),kind=wp) * dp_wp * cwp_fac ! in-cloud [g/m2]
        error_msg = cloud_optics%cloud_optics( &
          clwp, ciwp, &
          real(reshape(REFF(:,:,:,KLIQUID),(/ncol,LM/)),kind=wp), &
          real(reshape(REFF(:,:,:,KICE   ),(/ncol,LM/)),kind=wp), &
          cloud_props)
        TEST_(error_msg)
        deallocate(clwp, ciwp, __STAT__)

        call MAPL_TimerOff(MAPL,"---RRTMGP_CLOUD_OPTICS",__RC__)

        ! note: have made cloud_props for all ncol columns
        !   and will subset below into blocks ... we can also
        !   look at option of making cloud_props for each block 
        !   as its needed ... same for aer_props

        ! set desired cloud overlap type
        call MAPL_GetResource( &
          MAPL, cloud_overlap_type, "RRTMGP_CLOUD_OVERLAP_TYPE_LW:", &
          DEFAULT='MAX_RAN_OVERLAP', __RC__)

        call MAPL_TimerOn(MAPL,"---RRTMGP_MCICA",__RC__)

        ! ===============================================================================
        ! Random number setup:
        ! ===============================================================================
        !   We will use the Philox4x32-10 or ARS5 BRNGs from MKL VSL.
        !   Both are keyed families of counter-based PRNGs with a large period 2^130
        ! and a minimal state space (unlike the large Mersenne Twister state).
        !   Philox4x32-10 has a 64-bit key and a 128-bit counter and is very fast on GPUs.
        !   ARS5 has a 128-bit key and a 128-bit counter and is superfast on CPUs for
        ! which AES-NI instructions are hardware implemented.
        !   The SEEDING STRATEGY we will follow is to use a unique key for the gricolumn
        ! location and the simulation time. This gives a repeatable set of random numbers
        ! that remains the same for the members of an ensemble. If a different set is
        ! required for ensemble members, then the model state, such as the fractional
        ! part of the surface pressure, should be incorporated into the key.
        !   To get a different set of random numbers for the SW, for example, either a 
        ! key change or a counter advance will be needed.
        !
        ! Time Component of key:
        ! ~~~~~~~~~~~~~~~~~~~~~~
        ! 1. No need to update more frequently than once per LW refresh.
        ! 2. should reference the number of such intervals since a fixed time, so
        !   that agnostic to stop/restart schedule.
        !
        ! Space component of key:
        ! ~~~~~~~~~~~~~~~~~~~~~~~
        ! 1. should be based on some globally unique index for a gridcolumn, so that
        !   each gridcolumn is independent and so it is agnostic to runs with varying
        !   decompositions among processors.
        ! 2. 2^32 = 4,294,967,296 or about 2.1475e9 positives, which can represent 
        !   globe at over 1/180th degree resolution, so plenty for forseeable
        !   future.
        !
        ! Philox seeding:
        ! ~~~~~~~~~~~~~~~
        !   1. a scalar 32-bit seed sets the lower bits of the key k.
        !   2. a vector of 32-bit seeds of length N is used as follows:
        !      (a) N = 0:      k = c = 0;
        !      (b) N in {1,2}: seeds(1(:2)) set lower (and upper) words of key
        !      (c) N > 2:      ditto plus seeds(3:min(N,6)) set counter c,
        !                        starting from lowest word and working up.
        !
        ! Estimate of maximum LW random numbers needed per gridcolumn:
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! For the current clouds, with homogeneous optical properies in the cloudy part
        ! of each layer, only one random number per layer and gpt is needed. So estimate
        ! LM * ngpt <~ 132 * 256 = 33,792 < 2^16 = 65,536
        ! ===============================================================================

        allocate(seeds(3),__STAT__) ! 2-word key plus word1 of counter

        ! seed(1), the column part (word1) of key is set later
        ! but get required global indicies of local rectangular grid here
        call MAPL_GridGet(ESMFGRID, globalCellCountPerDim=dims, __RC__)
        IM_World = dims(1); JM_World = dims(2)
        call MAPL_GridGetInterior (ESMFGRID,iBeg,iEnd,jBeg,jEnd)

        ! get time part (word2) of key
        call ESMF_ClockGet(CLOCK, currTIME=CurrentTime, __RC__)
        call ESMF_TimeSet (ReferenceTime, yy=2000, mm=1, dd=1, __RC__)
        call ESMF_AlarmGet(ALARM, RINGINTERVAL=RefreshInterval, __RC__)
        seeds(2) = int((CurrentTime - ReferenceTime) / RefreshInterval)

        ! for LW start at counter=0
        seeds(3) = 0

        call MAPL_TimerOn(MAPL,"---RRTMGP_MCICA",__RC__)

      end if ! need_cloud_optical_props

      ! set aerosol optical properties
      if (need_dirty_optical_props) then

        call MAPL_TimerOn(MAPL,"---RRTMGP_AEROSOL_SETUP",__RC__)

        ! aerosol optics system is currently two-stream
        ! increment() will handle appropriate stream conversions
        allocate(ty_optical_props_2str::aer_props,__STAT__)
        select type (aer_props)
          class is (ty_optical_props_2str)

            ! band-only initialization
            TEST_(aer_props%alloc_2str(ncol, LM, k_dist%get_band_lims_wavenumber()))

            ! load unormalized optical properties from aerosol system
            aer_props%tau = real(reshape(TAUA, (/ncol,LM,nbnd/)), kind=wp)
            aer_props%ssa = real(reshape(SSAA, (/ncol,LM,nbnd/)), kind=wp)
            aer_props%g   = real(reshape(ASYA, (/ncol,LM,nbnd/)), kind=wp)

            ! renormalize
            where (aer_props%tau > 0._wp .and. aer_props%ssa > 0._wp )
              aer_props%g   = aer_props%g   / aer_props%ssa
              aer_props%ssa = aer_props%ssa / aer_props%tau
            elsewhere
              aer_props%tau = 0._wp
              aer_props%ssa = 0._wp
              aer_props%g   = 0._wp
            end where

            ! allocate a subset for blocking use
            allocate(ty_optical_props_2str::aer_props_subset,__STAT__)

          class default
            TEST_('aer_props hardwired 2-stream for now')

        end select

        call MAPL_TimerOff(MAPL,"---RRTMGP_AEROSOL_SETUP",__RC__)

      end if

      !--------------------------------------------------!
      ! Loop over subsets (blocks) of blockSize columns  !
      !  - choose rrtmgp_blockSize for efficiency        !
      !  - one possible partial block is done at the end ! 
      !--------------------------------------------------!

      call MAPL_GetResource( MAPL, &
        rrtmgp_blockSize, "RRTMGP_LW_BLOCKSIZE:", DEFAULT=4, __RC__)
      _ASSERT(rrtmgp_blockSize >= 1,'needs informative message')

      ! for random numbers, for efficiency, reserve the maximum possible
      ! subset of subcolumns (rrtmgp_blockSize) since column index is last
      allocate(urand(ngpt,LM,rrtmgp_blocksize), __STAT__)

      ! number of full blocks by integer division
      nBlocks = ncol/rrtmgp_blockSize

      ! allocate intermediate arrays for full blocks
      if (nBlocks > 0) then

        call MAPL_TimerOn(MAPL,"---RRTMGP_SUBSET",__RC__)

        ! block size until possible final partial block
        ncols_subset = rrtmgp_blockSize

        ! allocate clean_optical_props and sources for block
        select type (clean_optical_props)
          class is (ty_optical_props_1scl)
            TEST_(clean_optical_props%alloc_1scl(      ncols_subset, LM))
          class is (ty_optical_props_2str)
            TEST_(clean_optical_props%alloc_2str(      ncols_subset, LM))
          class is (ty_optical_props_nstr)
            TEST_(clean_optical_props%alloc_nstr(nmom, ncols_subset, LM))
        end select
        TEST_(sources%alloc(ncols_subset, LM))
        if (allocated(cld_mask)) then
          deallocate(cld_mask, __STAT__)
        endif
        allocate(cld_mask(ncols_subset, LM, ngpt), __STAT__)

        call MAPL_TimerOff(MAPL,"---RRTMGP_SUBSET",__RC__)

      end if

      ! add final partial block if necessary
      partial_block = mod(ncol, rrtmgp_blockSize) /= 0
      if (partial_block) then
        partial_blockSize = ncol - nBlocks * rrtmgp_blockSize
        nBlocks = nBlocks + 1
      endif

      ! loop over all blocks
      do b = 1, nBlocks

        call MAPL_TimerOn(MAPL,"---RRTMGP_SUBSET",__RC__)

        ! only the final block can be partial
        if (b == nBlocks .and. partial_block) then
          ncols_subset = partial_blockSize
          select type (clean_optical_props)
            class is (ty_optical_props_1scl)
              TEST_(clean_optical_props%alloc_1scl(      ncols_subset, LM))
            class is (ty_optical_props_2str)
              TEST_(clean_optical_props%alloc_2str(      ncols_subset, LM))
            class is (ty_optical_props_nstr)
              TEST_(clean_optical_props%alloc_nstr(nmom, ncols_subset, LM))
          end select
          TEST_(sources%alloc(ncols_subset, LM))
          if (allocated(cld_mask)) then
            deallocate(cld_mask, __STAT__)
          endif
        allocate(cld_mask(ncols_subset, LM, ngpt), __STAT__)
        endif

        ! prepare block
        colS = (b-1) * rrtmgp_blockSize + 1
        colE = colS + ncols_subset - 1
        TEST_(gas_concs%get_subset(colS, ncols_subset, gas_concs_subset))
        if (need_dirty_optical_props) then
          TEST_(aer_props%get_subset(colS, ncols_subset, aer_props_subset))
        end if

        call MAPL_TimerOff(MAPL,"---RRTMGP_SUBSET",__RC__)

        if (need_cloud_optical_props) then

          ! get column subset of the band-space in-cloud optical properties
          call MAPL_TimerOn(MAPL,"---RRTMGP_SUBSET",__RC__)
          TEST_(cloud_props%get_subset(colS, ncols_subset, cloud_props_subset))
          call MAPL_TimerOff(MAPL,"---RRTMGP_SUBSET",__RC__)

          call MAPL_TimerOn(MAPL,"---RRTMGP_MCICA",__RC__)

          ! generate McICA random numbers for subset
          ! Note: really only needed where cloud fraction > 0 (speedup?)
          ! Also, perhaps later this can be parallelized?
#ifdef HAVE_MKL
          do isub = 1, ncols_subset
            ! local 1d column index
            icol = colS + isub - 1
            ! local 2d indicies
            J = (icol-1) / IM + 1
            I = icol - (J-1) * IM
            ! initialize the Philox PRNG
            ! set word1 of key based on GLOBAL location
            ! 32-bits can hold all forseeable resolutions
            seeds(1) = (jBeg + J - 1) * IM_World + (iBeg + I - 1)
            ! instantiate a random number stream for the column
            call rng%init(VSL_BRNG_PHILOX4X32X10,seeds)
            ! draw the random numbers for the column
            urand(:,:,isub) = reshape(rng%get_random(ngpt*LM),(/ngpt,LM/))
            ! free the rng
            call rng%end()
          end do
#endif

          ! cloud sampling to gpoints
          select case (cloud_overlap_type)
            case ("MAX_RAN_OVERLAP")
              error_msg = sampled_mask_max_ran( &
                urand(:,:,1:ncols_subset), cf_wp(colS:colE,:), cld_mask)
              TEST_(error_msg)
            case ("EXP_RAN_OVERLAP")
              TEST_('EXP_RAN_OVERLAP not yet implemnted')
              !error_msg = (sampled_mask_exp_ran())
              !TEST_(error_msg)
            case default
              TEST_('RRTMGP_LW: unknown cloud overlap')
          end select

          ! draw McICA optical property samples (band->gpt)
          select type (cloud_props_gpt)
            class is (ty_optical_props_2str)
              TEST_(cloud_props_gpt%alloc_2str(ncols_subset, LM, k_dist))
            class default
              TEST_('cloud_props_gpt hardwired 2-stream for now')
          end select
          TEST_(draw_samples(cld_mask, cloud_props_subset, cloud_props_gpt))

          call MAPL_TimerOff(MAPL,"---RRTMGP_MCICA",__RC__)

        end if

        call MAPL_TimerOn(MAPL,"---RRTMGP_GAS_OPTICS",__RC__)

        ! get gas optical properties and sources
        error_msg = k_dist%gas_optics( &
          p_lay(colS:colE,:), p_lev(colS:colE,:), t_lay(colS:colE,:), &
          t_sfc(colS:colE), gas_concs_subset, clean_optical_props, sources, &
          tlev = t_lev(colS:colE,:))
        TEST_(error_msg)

        call MAPL_TimerOff(MAPL,"---RRTMGP_GAS_OPTICS",__RC__)

        call MAPL_TimerOn(MAPL,"---RRTMGP_RT",__RC__)

        ! clean clear-sky case
        if (calc_clrnoa) then
          fluxes_clrnoa%flux_up => flux_up_clrnoa(colS:colE,:)
          fluxes_clrnoa%flux_dn => flux_dn_clrnoa(colS:colE,:)
          error_msg = rte_lw( &
            clean_optical_props, &
            top_at_1, sources, emis_sfc(:,colS:colE), &
            fluxes_clrnoa, n_gauss_angles=nga, use_2stream=u2s, &
            flux_up_Jac=dfupdts_clrnoa(colS:colE,:))
          TEST_(error_msg)
        end if

        if (need_dirty_optical_props) then
          ! make copy of clrnoa optical properties as the
          !   starting point for later dirty calculations
          select type (dirty_optical_props)
            class is (ty_optical_props_1scl)
              TEST_(dirty_optical_props%alloc_1scl(ncols_subset, LM, clean_optical_props))
            class is (ty_optical_props_2str)
              TEST_(dirty_optical_props%alloc_2str(ncols_subset, LM, clean_optical_props))
              select type (clean_optical_props)
                class is (ty_optical_props_2str)
                  dirty_optical_props%ssa = clean_optical_props%ssa
                  dirty_optical_props%g   = clean_optical_props%g  
              end select
            class is (ty_optical_props_nstr)
              TEST_(dirty_optical_props%alloc_nstr(nmom, ncols_subset, LM, clean_optical_props))
              select type (clean_optical_props)
                class is (ty_optical_props_nstr)
                  dirty_optical_props%ssa = clean_optical_props%ssa
                  dirty_optical_props%p   = clean_optical_props%p  
              end select
          end select
          ! all streams have tau
          dirty_optical_props%tau = clean_optical_props%tau
        end if
        
        ! clean all-sky case
        if (calc_allnoa) then

          ! add in cloud optical properties
          TEST_(cloud_props_gpt%increment(clean_optical_props))

          ! clean all-sky RT
          fluxes_allnoa%flux_up => flux_up_allnoa(colS:colE,:)
          fluxes_allnoa%flux_dn => flux_dn_allnoa(colS:colE,:)
          error_msg = rte_lw( &
            clean_optical_props, &
            top_at_1, sources, emis_sfc(:,colS:colE), &
            fluxes_allnoa, n_gauss_angles=nga, use_2stream=u2s, &
            flux_up_Jac=dfupdts_allnoa(colS:colE,:))
          TEST_(error_msg)
        end if

        if (export_clrsky .or. export_allsky) then
          if (implements_aerosol_optics) then

            ! dirty flux calculations required ...

            ! "dirty_optical_props" is currently just a copy of the clrnoa optical_props
            !   so must now add in aerosols to make it actually dirty
            TEST_(aer_props_subset%increment(dirty_optical_props))

            ! dirty clear-sky RT
            if (calc_clrsky) then
              fluxes_clrsky%flux_up => flux_up_clrsky(colS:colE,:)
              fluxes_clrsky%flux_dn => flux_dn_clrsky(colS:colE,:)
              error_msg = rte_lw( &
                dirty_optical_props, &
                top_at_1, sources, emis_sfc(:,colS:colE), &
                fluxes_clrsky, n_gauss_angles=nga, use_2stream=u2s, &
                flux_up_Jac=dfupdts_clrsky(colS:colE,:))
              TEST_(error_msg)
            end if

            ! dirty all-sky case
            if (calc_allsky) then

              ! add in cloud optical properties
              TEST_(cloud_props_gpt%increment(dirty_optical_props))

              ! dirty all-sky RT
              fluxes_allsky%flux_up => flux_up_allsky(colS:colE,:)
              fluxes_allsky%flux_dn => flux_dn_allsky(colS:colE,:)
              error_msg = rte_lw( &
                dirty_optical_props, &
                top_at_1, sources, emis_sfc(:,colS:colE), &
                fluxes_allsky, n_gauss_angles=nga, use_2stream=u2s, &
                flux_up_Jac=dfupdts_allsky(colS:colE,:))
              TEST_(error_msg)
            end if

          else

            ! there are no aerosols so we are done because the 
            !   dirty cases are the same as the clean ones
            if (export_clrsky) then
              flux_up_clrsky(colS:colE,:) = flux_up_clrnoa(colS:colE,:)
              flux_dn_clrsky(colS:colE,:) = flux_dn_clrnoa(colS:colE,:)
              dfupdts_clrsky(colS:colE,:) = dfupdts_clrnoa(colS:colE,:)
            end if
            if (export_allsky) then
              flux_up_allsky(colS:colE,:) = flux_up_allnoa(colS:colE,:)
              flux_dn_allsky(colS:colE,:) = flux_dn_allnoa(colS:colE,:)
              dfupdts_allsky(colS:colE,:) = dfupdts_allnoa(colS:colE,:)
            end if
          
          end if ! implements_aerosol_optics
        end if ! export dirty clear-sky or all-sky

        call MAPL_TimerOff(MAPL,"---RRTMGP_RT",__RC__)

      end do ! loop over blocks

      call MAPL_TimerOn(MAPL,"---RRTMGP_POST",__RC__)

      ! load output arrays
      ! note: the DFDTS* are the derivatives of the NEGATED upward fluxes wrt TS
      if (export_clrnoa) then
        FLAU_INT = real(reshape(-flux_up_clrnoa, (/IM,JM,LM+1/)))
        FLAD_INT = real(reshape( flux_dn_clrnoa, (/IM,JM,LM+1/)))
        DFDTSCNA = real(reshape(-dfupdts_clrnoa, (/IM,JM,LM+1/)))
      end if
      if (export_allnoa) then
        FLXAU_INT = real(reshape(-flux_up_allnoa, (/IM,JM,LM+1/)))
        FLXAD_INT = real(reshape( flux_dn_allnoa, (/IM,JM,LM+1/)))
        DFDTSNA   = real(reshape(-dfupdts_allnoa, (/IM,JM,LM+1/)))
      end if
      if (export_clrsky) then
        FLCU_INT = real(reshape(-flux_up_clrsky, (/IM,JM,LM+1/)))
        FLCD_INT = real(reshape( flux_dn_clrsky, (/IM,JM,LM+1/)))
        DFDTSC   = real(reshape(-dfupdts_clrsky, (/IM,JM,LM+1/)))
      end if
      if (export_allsky) then
        FLXU_INT = real(reshape(-flux_up_allsky, (/IM,JM,LM+1/)))
        FLXD_INT = real(reshape( flux_dn_allsky, (/IM,JM,LM+1/)))
        DFDTS    = real(reshape(-dfupdts_allsky, (/IM,JM,LM+1/)))
      end if

      !mjs: Corrected emitted at the surface to remove reflected
      !     from upward. Note that emiss is the same for all bands,
      !     so we use band 1 for the total flux.
      SFCEM_INT = real(reshape( &
        -flux_up_allsky(:,LM+1) + flux_dn_allsky(:,LM+1) * (1._wp - emis_sfc(1,:)), &
        (/IM,JM/)))

      ! clean up
      call sources%finalize()
      call clean_optical_props%finalize()
      if (need_dirty_optical_props) then
        call dirty_optical_props%finalize()
        call aer_props_subset%finalize()
        call aer_props%finalize()
      end if
      if (need_cloud_optical_props) then
        call cloud_props_gpt%finalize()
        call cloud_props_subset%finalize()
        deallocate(seeds,__STAT__)
        deallocate(urand, __STAT__)
        deallocate(cld_mask,__STAT__)
        call cloud_props%finalize()
        call cloud_optics%finalize()
      end if
      if (calc_clrnoa) then
        deallocate(flux_up_clrnoa, flux_dn_clrnoa, dfupdts_clrnoa, __STAT__)
      end if
      if (calc_allnoa) then
        deallocate(flux_up_allnoa, flux_dn_allnoa, dfupdts_allnoa, __STAT__)
      end if
      if (calc_clrsky) then
        deallocate(flux_up_clrsky, flux_dn_clrsky, dfupdts_clrsky, __STAT__)
      end if
      if (calc_allsky) then
        deallocate(flux_up_allsky, flux_dn_allsky, dfupdts_allsky, __STAT__)
      end if
      deallocate(p_lay, t_lay, p_lev, t_lev, dp_wp, cf_wp, __STAT__)
      deallocate(t_sfc, emis_sfc, __STAT__)

      call MAPL_TimerOff(MAPL,"---RRTMGP_POST",__RC__)

      call MAPL_TimerOff(MAPL,"--RRTMGP",__RC__)

   else if (USE_RRTMG) then

      call MAPL_TimerOn(MAPL,"--RRTMG",RC=STATUS)
      VERIFY_(STATUS)
 
      call MAPL_GetResource(MAPL,PARTITION_SIZE,'RRTMGLW_PARTITION_SIZE:',DEFAULT=4,RC=STATUS)
      VERIFY_(STATUS)

      ! reversed profiles for RRTMG (1=bottom layer)
      ! note 0:LM indexing for [PT]LEV_R
      ! but 1:LM+1 for [UD]FLX[C] and DUFLX[C]_DTS
      allocate(FCLD_R(IM*JM,LM),__STAT__)
      allocate(TLEV_R(IM*JM,0:LM),__STAT__)
      allocate(PLE_R(IM*JM,0:LM),__STAT__)
      allocate(ZM_R(IM*JM,LM),__STAT__)
      allocate(EMISS(IM*JM,NB_RRTMG),__STAT__)
      allocate(CLIQWP(IM*JM,LM),__STAT__)
      allocate(CICEWP(IM*JM,LM),__STAT__)
      allocate(RELIQ(IM*JM,LM),__STAT__)
      allocate(REICE(IM*JM,LM),__STAT__)
      allocate(TAUAER(IM*JM,LM,NB_RRTMG),__STAT__)
      allocate(PL_R(IM*JM,LM),__STAT__)
      allocate(T_R(IM*JM,LM),__STAT__)
      allocate(Q_R(IM*JM,LM),__STAT__)
      allocate(O2_R(IM*JM,LM),__STAT__)
      allocate(O3_R(IM*JM,LM),__STAT__)
      allocate(CO2_R(IM*JM,LM),__STAT__)
      allocate(CH4_R(IM*JM,LM),__STAT__)
      allocate(N2O_R(IM*JM,LM),__STAT__)
      allocate(CFC11_R(IM*JM,LM),__STAT__)
      allocate(CFC12_R(IM*JM,LM),__STAT__)
      allocate(CFC22_R(IM*JM,LM),__STAT__)
      allocate(CCL4_R(IM*JM,LM),__STAT__)
      allocate(TSFC(IM*JM),__STAT__)
      allocate(UFLX(IM*JM,LM+1),__STAT__)
      allocate(DFLX(IM*JM,LM+1),__STAT__)
      allocate(UFLXC(IM*JM,LM+1),__STAT__)
      allocate(DFLXC(IM*JM,LM+1),__STAT__)
      allocate(DUFLX_DTS(IM*JM,LM+1),__STAT__)
      allocate(DUFLXC_DTS(IM*JM,LM+1),__STAT__)
      allocate(CLEARCOUNTS(IM*JM,4),__STAT__)
      allocate(ALAT(IM*JM),__STAT__)
      allocate(OLRBRG      (nbndlw,IM*JM),__STAT__)
      allocate(DOLRBRG_DTS (nbndlw,IM*JM),__STAT__)

      ! choices for cloud physical to optical conversion
      ICEFLGLW = 3
      LIQFLGLW = 1

      ! calculate derivatives of upward flux with Tsurf
      Ts_derivs = .true.

      call MAPL_TimerOn(MAPL,"---RRTMG_FLIP",RC=STATUS)
      VERIFY_(STATUS)
 
      ! reverse super-layer interface indicies
      LCLDMH = LM - LCLDMH + 1
      LCLDLM = LM - LCLDLM + 1

      ! collapse horizontal indicies and flip in vertical
      !   (RRTMG indexed bottom to top)
      IJ = 0
      do J = 1,JM
      do I = 1,IM
         IJ = IJ + 1

         TSFC (IJ)   = TS  (I,J)
         EMISS(IJ,:) = EMIS(I,J) ! all bands get same emissivity
         ALAT (IJ)   = LATS(I,J)

         ! calculation of level temperature (still in model ordering)
         ! note: PLE(0:LM) but TLEV(1:LM+1)
         DP(1) = PLE(I,J,1)-PLE(I,J,0)
         do K = 2,LM
            DP(K) = (PLE(I,J,K)-PLE(I,J,K-1) )
            TLEV(K) = (T(I,J,K-1) * DP(K) + T(I,J,K) * DP(K-1)) &
                      / (DP(K-1) + DP(K))
         enddo
         TLEV(LM+1) = T2M(I,J) ! 'surface'
         TLEV(   1) = TLEV(2)  ! model top

         !  Flip in vertical
         do K = 1,LM
            LV = LM-K+1  ! LM --> 1

            ! Convert content [kg/kg] to path [g/m2]
            ! using hydrostatic eqn dp/g ~ rho*dz,
            ! so conversion factor is 1000*dp/g ~ 1.02*100*dp.
            ! pmn: why not use MAPL_GRAV explicitly?
            xx = 1.02*100*DP(LV)
            CLIQWP(IJ,K) = xx*CWC(I,J,LV,KLIQUID)
            CICEWP(IJ,K) = xx*CWC(I,J,LV,KICE)
            RELIQ (IJ,K) =   REFF(I,J,LV,KLIQUID)
            REICE (IJ,K) =   REFF(I,J,LV,KICE   )
               
            ! impose RRTMG re_liq limits
            if    (LIQFLGLW.eq.0) then
               ! pmn: this one not available inside RRTMG_LW
               RELIQ(IJ,K) = min(max(RELIQ(IJ,K),5.0),10.0)
            elseif (LIQFLGLW.eq.1) then
               RELIQ(IJ,K) = min(max(RELIQ(IJ,K),2.5),60.0)
            endif

            ! impose RRTMG re_ice limits
            if     (ICEFLGLW.eq.0) then
               REICE(IJ,K) = min(max(REICE(IJ,K),10.0),30.0)
            elseif (ICEFLGLW.eq.1) then
               REICE(IJ,K) = min(max(REICE(IJ,K),13.0),130.0)
            elseif (ICEFLGLW.eq.2) then
               REICE(IJ,K) = min(max(REICE(IJ,K), 5.0),131.0)
            elseif (ICEFLGLW.eq.3) then
               REICE(IJ,K) = min(max(REICE(IJ,K), 5.0),140.0)
            elseif (ICEFLGLW.eq.4) then
               REICE(IJ,K) = min(max(REICE(IJ,K)*2.,1.0),200.0)
            endif

            ! flipping for LEVEL quantities
            ! PLE_R(0:LM) = PLE(LM:0)
            ! TLEV_R(0:LM) = TLEV(LM+1:1)
            ! top-of-model LEVEL (RRTMG LM) done later
            PLE_R  (IJ,K-1) = PLE(I,J,LV)/100. ! [hPa]
            TLEV_R (IJ,K-1) = TLEV(LV+1)

            ! more flipping for layer quantities
            ! Q [specific humidity] --> Q_R [volume mixing ratio]
            ! O3 [mass mixing ratio] --> O3_R [volume mixing ratio]
            PL_R   (IJ,K) = PL(I,J,LV)/100.  ! [hPa]
            T_R    (IJ,K) = T(I,J,LV)
            Q_R    (IJ,K) = Q(I,J,LV) / (1.-Q(I,J,LV)) * (MAPL_AIRMW/MAPL_H2OMW)
            O3_R   (IJ,K) = O3(I,J,LV) * (MAPL_AIRMW/MAPL_O3MW)
            CH4_R  (IJ,K) = CH4(I,J,LV)
            N2O_R  (IJ,K) = N2O(I,J,LV)
            CO2_R  (IJ,K) = CO2
            O2_R   (IJ,K) = O2
            CCL4_R (IJ,K) = CCL4
            CFC11_R(IJ,K) = CFC11(I,J,LV)
            CFC12_R(IJ,K) = CFC12(I,J,LV)
            CFC22_R(IJ,K) = HCFC22(I,J,LV)
            FCLD_R (IJ,K) = FCLD(I,J,LV)

            ! RRTMG_LW does not scatter, so pass ABSORPTION aerosol
            ! optical thickness to RRTMG. Remember that SSAA is the
            ! aerosol system's *un*-normalized single scattering albedo,
            ! which is actually tau_ext * omega0 = tau_scat, and TAUA
            ! is the aerosol extinction optical thickness.
            TAUAER(IJ,K,:) = TAUA(I,J,LV,:) - SSAA(I,J,LV,:)

         enddo

         ! finish off top-of-model LEVEL
         PLE_R (IJ,LM) = PLE(I,J,0)/100. ! [hPa]
         TLEV_R(IJ,LM) = TLEV(1)

         ! Calculate the LAYER (mid-point) heights.
         ! The interlayer distances are needed for the calculations
         ! of inter-layer correlation for cloud overlapping in RRTMG.
         ! Only *relative* distances matter, so wolog set ZM_R(1) = 0.
         ! pmn: 2021-04-21 this calculation was wrong in earlier revisions.
         ZM_R(IJ,1) = 0.
         do K=2,LM
            ! dz ~ RT/g x dp/p by hysrostatic eqn and ideal gas eqn.
            ! The jump from LAYER k-1 to k is centered on LEVEL k-1
            !   since the RRTMG LEVEL (LE[V]_R) indices are zero-based
            ZM_R(IJ,K) = ZM_R(IJ,K-1) + MAPL_RGAS*TLEV_R(IJ,K-1)/MAPL_GRAV &
                                        * (PL_R(IJ,K-1)-PL_R(IJ,K))/PLE_R(IJ,K-1)
         enddo

      enddo ! IM
      enddo ! JM

      call MAPL_TimerOff(MAPL,"---RRTMG_FLIP",RC=STATUS)
      VERIFY_(STATUS)
 
      call MAPL_TimerOn(MAPL,"---RRTMG_INIT",RC=STATUS)
      VERIFY_(STATUS)
 
! pmn: consider putting futher up calling tree?
! pmn: only needs to be done once per run, but does consume memory
      call RRTMG_LW_INI

      call MAPL_TimerOff(MAPL,"---RRTMG_INIT",RC=STATUS)
      VERIFY_(STATUS)
 
      call MAPL_TimerOn(MAPL,"---RRTMG_RUN",RC=STATUS)
      VERIFY_(STATUS)

      call RRTMG_LW (IM*JM, LM, PARTITION_SIZE, TS_DERIVS, &
              PL_R, PLE_R, T_R, TLEV_R, TSFC, EMISS, &
              Q_R, O3_R, CO2_R, CH4_R, N2O_R, O2_R, &
              CFC11_R, CFC12_R, CFC22_R, CCL4_R, &
              FCLD_R, CICEWP, CLIQWP, REICE, RELIQ, ICEFLGLW, LIQFLGLW, &
              TAUAER, ZM_R, ALAT, DOY, LCLDLM, LCLDMH, CLEARCOUNTS, &
              UFLX, DFLX, UFLXC, DFLXC, DUFLX_DTS, DUFLXC_DTS, &
              BAND_OUTPUT, OLRBRG, DOLRBRG_DTS)

      call MAPL_TimerOff(MAPL,"---RRTMG_RUN",RC=STATUS)
      VERIFY_(STATUS)
 
      call MAPL_TimerOn(MAPL,"---RRTMG_FLIP",RC=STATUS)
      VERIFY_(STATUS)
 
      ! for outputs, unpack flattened horizontal and flip back vertical
      IJ = 0
      do J = 1,JM
      do I = 1,IM
         IJ = IJ + 1

         ! convert super-layer clearCounts to cloud fractions
         if(associated(CLDTTLW)) then
            CLDTTLW(I,J) = 1.0 - CLEARCOUNTS(IJ,1)/float(NGPTLW)
         endif
         if(associated(CLDHILW)) then
            CLDHILW(I,J) = 1.0 - CLEARCOUNTS(IJ,2)/float(NGPTLW)
         endif
         if(associated(CLDMDLW)) then
            CLDMDLW(I,J) = 1.0 - CLEARCOUNTS(IJ,3)/float(NGPTLW)
         endif
         if(associated(CLDLOLW)) then
            CLDLOLW(I,J) = 1.0 - CLEARCOUNTS(IJ,4)/float(NGPTLW)
         endif

         ! upward negative in GEOS-5 convention
         do K = 0,LM
            LV = LM-K+1
            FLXU_INT(I,J,K) =-UFLX      (IJ,LV)
            FLXD_INT(I,J,K) = DFLX      (IJ,LV)
            FLCU_INT(I,J,K) =-UFLXC     (IJ,LV)
            FLCD_INT(I,J,K) = DFLXC     (IJ,LV)
            DFDTS   (I,J,K) =-DUFLX_DTS (IJ,LV)
            DFDTSC  (I,J,K) =-DUFLXC_DTS(IJ,LV)
         enddo

         ! Reflected LW is not counted in surface emitted. Also, for now,
         ! surface emitted is positive downwards consistent with Chou-Suarez.
         ! (Note: All bands use the same emissivity)
         SFCEM_INT(I,J) = -( UFLX(IJ,1) - DFLX(IJ,1)*(1.-EMIS(I,J)) )

      enddo ! IM
      enddo ! JM

      ! band OLR and brightness temperatures
      do ibnd = 1,nbndlw
         if (band_output(ibnd)) then
            write(bb,'(I0.2)') ibnd

            call MAPL_GetPointer(INTERNAL, ptr2d, 'OLRB'//bb//'RG', __RC__)
            ptr2d = reshape(OLRBRG (ibnd,:), [IM,JM])

            call MAPL_GetPointer(INTERNAL, ptr2d, 'DOLRB'//bb//'RGDT', __RC__)
            ptr2d = reshape(DOLRBRG_DTS (ibnd,:), [IM,JM])

         end if
      end do

      call MAPL_TimerOff(MAPL,"---RRTMG_FLIP",RC=STATUS)
      VERIFY_(STATUS)

      ! pmn:
      ! RRTMG does not provide no-aerosol derivatives
      ! so set no-aerosol to aerosol derivatives
      DFDTSNA  = DFDTS
      DFDTSCNA = DFDTSC

      deallocate(FCLD_R,__STAT__)
      deallocate(TLEV_R,__STAT__)
      deallocate(PLE_R,__STAT__)
      deallocate(ZM_R,__STAT__)
      deallocate(EMISS,__STAT__)
      deallocate(CLIQWP,__STAT__)
      deallocate(CICEWP,__STAT__)
      deallocate(RELIQ,__STAT__)
      deallocate(REICE,__STAT__)
      deallocate(TAUAER,__STAT__)
      deallocate(PL_R,__STAT__)
      deallocate(T_R,__STAT__)
      deallocate(Q_R,__STAT__)
      deallocate(O2_R,__STAT__)
      deallocate(O3_R,__STAT__)
      deallocate(CO2_R,__STAT__)
      deallocate(CH4_R,__STAT__)
      deallocate(N2O_R,__STAT__)
      deallocate(CFC11_R,__STAT__)
      deallocate(CFC12_R,__STAT__)
      deallocate(CFC22_R,__STAT__)
      deallocate(CCL4_R,__STAT__)
      deallocate(TSFC,__STAT__)
      deallocate(UFLX,__STAT__)
      deallocate(DFLX,__STAT__)
      deallocate(UFLXC,__STAT__)
      deallocate(DFLXC,__STAT__)
      deallocate(DUFLX_DTS,__STAT__)
      deallocate(DUFLXC_DTS,__STAT__)
      deallocate(CLEARCOUNTS,__STAT__)
      deallocate(ALAT,__STAT__)
      deallocate(OLRBRG,__STAT__)
      deallocate(DOLRBRG_DTS,__STAT__)

      call MAPL_TimerOff(MAPL,"--RRTMG",RC=STATUS)
      VERIFY_(STATUS)

   else

      ! Something is wrong. We've selected neither Chou or RRTMG
      _ASSERT(.false.,'needs informative message')

   end if SCHEME

! Sum up the U and D fluxes to get net downward

   FLX_INT  = FLXD_INT  + FLXU_INT
   FLXA_INT = FLXAD_INT + FLXAU_INT
   FLC_INT  = FLCD_INT  + FLCU_INT
   FLA_INT  = FLAD_INT  + FLAU_INT
   
   ! Revert to SFCEM to a positive quantity.
   ! Earlier surface emitted positive downwards per Chou-Suarez.
   SFCEM_INT = -SFCEM_INT

! Save surface temperature in internal state
!-------------------------------------------

   TS_INT    = TS

! Export some cloud properties in the infrared
!---------------------------------------------

   call MAPL_TimerOn (MAPL,"--MISC")

   call MAPL_GetResource( MAPL, TAUCRIT, 'TAUCRIT:', DEFAULT=0.30, RC=STATUS)
   VERIFY_(STATUS)
   TAUCRIT   = TAUCRIT/2.13

   call MAPL_GetPointer(EXPORT,   CLDPRS,  'CLDPRS'  ,RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,   CLDTMP,  'CLDTMP'  ,RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,    TAUIR,   'TAUIR'  ,RC=STATUS); VERIFY_(STATUS)

   if(associated(TAUIR)) TAUIR = 0.5*(TAUDIAG(:,:,:,3)+TAUDIAG(:,:,:,4))

   if(associated(CLDTMP).or.associated(CLDPRS)) then
      if(associated(CLDTMP)) CLDTMP = MAPL_UNDEF
      if(associated(CLDPRS)) CLDPRS = MAPL_UNDEF
      do j=1,jm
         do i=1,im
            do l=1,lm
               if(0.5*(TAUDIAG(I,J,L,3)+TAUDIAG(I,J,L,4))>TAUCRIT) then
                  if(associated(CLDTMP)) CLDTMP(I,J) = T  (I,J,L)
                  if(associated(CLDPRS)) CLDPRS(I,J) = PLE(I,J,L-1)
                  exit
               end if
            end do
         end do
      end do
   end if

   ! Correcting the timing of the alw and blw (mjs)

   call MAPL_GetPointer(EXPORT,   TSREFF,    'TSREFF' ,RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,   SFCEM,     'SFCEM0' ,RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,   DSFDTS,    'DSFDTS0',RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,   LWS0,      'LWS0'   ,RC=STATUS); VERIFY_(STATUS)

   if(associated(TSREFF)) TSREFF = TS             ! reference TS for linearization
   if(associated(DSFDTS)) DSFDTS =-DFDTS(:,:,LM)  ! d(non-negated upward sfc flux) / dTS
   if(associated(SFCEM )) SFCEM  = SFCEM_INT      ! sfc emitted flux (+ve)
   if(associated(LWS0  )) LWS0   = &              ! absorbed (not reflected)
     FLX_INT(:,:,LM) + SFCEM_INT                  !   downward sfc flux (+ve)

   ! Deallocate per-band aerosol arrays
   ! ----------------------------------

   DEALLOCATE(TAUA)
   DEALLOCATE(SSAA)
   DEALLOCATE(ASYA)

   call MAPL_TimerOff(MAPL,"--MISC")

!  All done
!-----------

   RETURN_(ESMF_SUCCESS)

 end subroutine LW_Driver

!------------------------------------------------
!------------------------------------------------

 subroutine Update_Flx(IM,JM,LM,RC)
   integer,           intent(IN ) :: IM, JM, LM
   integer, optional, intent(OUT) :: RC

!  Locals

   character(len=ESMF_MAXSTR)        :: Iam
   integer                           :: STATUS

   real,          dimension(IM,JM)   :: DELT
   integer                           :: K
   integer                           :: LEV_LOW_MID
   integer                           :: LEV_MID_HIGH
   real                              :: PRS_LOW_MID                   ! pressure separating low and middle clouds
   real                              :: PRS_MID_HIGH                  ! pressure separating low and high   clouds

   ! band wavenumber bounds (m-1)
   real :: wn1, wn2

! pointer to import

   real, pointer, dimension(:,:  )   :: TSINST

! pointers to export

   real, pointer, dimension(:,:,:)   :: FLX
   real, pointer, dimension(:,:,:)   :: FLXA
   real, pointer, dimension(:,:,:)   :: FLC
   real, pointer, dimension(:,:,:)   :: FLA
   real, pointer, dimension(:,:,:)   :: FLXU
   real, pointer, dimension(:,:,:)   :: FLXAU
   real, pointer, dimension(:,:,:)   :: FLCU
   real, pointer, dimension(:,:,:)   :: FLAU
   real, pointer, dimension(:,:,:)   :: FLXD
   real, pointer, dimension(:,:,:)   :: FLXAD
   real, pointer, dimension(:,:,:)   :: FLCD
   real, pointer, dimension(:,:,:)   :: FLAD
   real, pointer, dimension(:,:  )   :: TSREFF
   real, pointer, dimension(:,:  )   :: SFCEM
   real, pointer, dimension(:,:  )   :: DSFDTS
   real, pointer, dimension(:,:  )   :: SFCEM0
   real, pointer, dimension(:,:  )   :: DSFDTS0
   real, pointer, dimension(:,:  )   :: OLR
   real, pointer, dimension(:,:  )   :: OLRA
   real, pointer, dimension(:,:  )   :: OLC
   real, pointer, dimension(:,:  )   :: OLCC5
   real, pointer, dimension(:,:  )   :: OLA
   real, pointer, dimension(:,:  )   :: FLNS
   real, pointer, dimension(:,:  )   :: FLNSNA
   real, pointer, dimension(:,:  )   :: FLNSC
   real, pointer, dimension(:,:  )   :: FLNSA
   real, pointer, dimension(:,:  )   :: LWS
   real, pointer, dimension(:,:  )   :: LWSA
   real, pointer, dimension(:,:  )   :: LCS
   real, pointer, dimension(:,:  )   :: LCSC5
   real, pointer, dimension(:,:  )   :: LAS
   real, pointer, dimension(:,:  )   :: CLDTT
   real, pointer, dimension(:,:  )   :: ptr2d

   real, pointer, dimension(:,:,:)   :: FCLD
   real, pointer, dimension(:    )   :: PREF

   real, allocatable, dimension(:,:) :: DUMTT, OLRB

!  Begin...
!----------

   IAm = "Update_Flx"

! Pointers to Exports
!--------------------

   call MAPL_GetPointer(EXPORT,   FLX   ,    'FLX',   RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,   FLXA  ,    'FLXA',  RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,   FLC   ,    'FLC',   RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,   FLA   ,    'FLA',   RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,   FLXU  ,    'FLXU',  RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,   FLXAU ,    'FLXAU', RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,   FLCU  ,    'FLCU',  RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,   FLAU  ,    'FLAU',  RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,   FLXD  ,    'FLXD',  RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,   FLXAD ,    'FLXAD', RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,   FLCD  ,    'FLCD',  RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,   FLAD  ,    'FLAD',  RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,   TSREFF,    'TSREFF',RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,   SFCEM ,    'SFCEM', RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,   DSFDTS,    'DSFDTS',RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,   SFCEM0,    'SFCEM0',RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,  DSFDTS0,   'DSFDTS0',RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,   OLR   ,    'OLR'   ,RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,   OLRA  ,    'OLRA'  ,RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,   OLC   ,    'OLC'   ,RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,   OLCC5 ,    'OLCC5' ,RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,   OLA   ,    'OLA'   ,RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,   LWS   ,    'LWS'   ,RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,   LWSA  ,    'LWSA'  ,RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,   LCS   ,    'LCS'   ,RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,   LCSC5 ,    'LCSC5' ,RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,   LAS   ,    'LAS'   ,RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,   FLNS  ,    'FLNS'  ,RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,   FLNSNA,    'FLNSNA',RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,   FLNSC ,    'FLNSC' ,RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,   FLNSA ,    'FLNSA' ,RC=STATUS); VERIFY_(STATUS)

   call MAPL_GetPointer(EXPORT,   CLDTT ,  'CLDTT'   ,ALLOC=.TRUE.,RC=STATUS); VERIFY_(STATUS)

! Determine the 2-D Total Cloud Fraction
!---------------------------------------

   call MAPL_GetResource( MAPL, PRS_LOW_MID,    'PRS_LOW_MID_CLOUDS:' ,   DEFAULT=70000.,      RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetResource( MAPL, PRS_MID_HIGH,   'PRS_MID_HIGH_CLOUDS:',   DEFAULT=40000.,      RC=STATUS)
   VERIFY_(STATUS)

   call MAPL_GetPointer( IMPORT, FCLD, 'FCLD', RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer( IMPORT, PREF, 'PREF', RC=STATUS)
   VERIFY_(STATUS)
   
   ALLOCATE( DUMTT(IM,JM), STAT=STATUS)
   VERIFY_(STATUS)

! Determine the model level separating mid and high clouds
!---------------------------------------------------------
   LEV_MID_HIGH = 1
   do K = 1, LM
      if( PREF(K) >= PRS_MID_HIGH ) then
         LEV_MID_HIGH = K
         exit
      end if
   end do

! Determine the model level seperating low and middle clouds
!-----------------------------------------------------------
   LEV_LOW_MID = LM
   do K = 1, LM
      if( PREF(K) >= PRS_LOW_MID  ) then
         LEV_LOW_MID = K
         exit
      end if
   end do

         DUMTT = 0.
         do K=1,LEV_MID_HIGH-1
            DUMTT = max(DUMTT,FCLD(:,:,K))
         end do
         CLDTT = (1-DUMTT)
         DUMTT = 0.
         do K= LEV_MID_HIGH,LEV_LOW_MID-1
            DUMTT = max(DUMTT,FCLD(:,:,K))
         end do
         CLDTT = CLDTT*(1-DUMTT)
         DUMTT = 0.
         do K=LEV_LOW_MID,LM
            DUMTT = max(DUMTT,FCLD(:,:,K))
         end do
         CLDTT = 1.0 - CLDTT*(1-DUMTT)

! Pointers to Imports
!--------------------

   call MAPL_GetPointer(IMPORT,   TSINST, 'TSINST',   RC=STATUS); VERIFY_(STATUS)

! Update fluxes
!--------------

   ! linearization with surface temperature notes:
   ! a. only upward fluxes linearized wrt surface temperature
   ! b. the derivatives DFDTS[C] have the same sign convention as the negated upward fluxes
   !      (i.e., they are the derivatives of negated upward fluxes with surface temperature)

   ! surface temperature change since refresh for linearization
   DELT = TSINST - TS_INT

   if( USE_CHOU .or. USE_RRTMGP ) THEN

       ! fill 3D fluxes
       do K = 0, LM
          ! net downward (downward plus negated upward) fluxes
          if(associated(FLX))    FLX (:,:,K) =   FLX_INT(:,:,K) + DFDTS   (:,:,K) * DELT ! all-sky
          if(associated(FLXA))  FLXA (:,:,K) =  FLXA_INT(:,:,K) + DFDTSNA (:,:,K) * DELT ! all-sky no-aerosol
          if(associated(FLC))    FLC (:,:,K) =   FLC_INT(:,:,K) + DFDTSC  (:,:,K) * DELT ! clr-sky
          if(associated(FLA))    FLA (:,:,K) =   FLA_INT(:,:,K) + DFDTSCNA(:,:,K) * DELT ! clr-sky no-aerosol
          ! negated upward fluxes
          if(associated(FLXU))   FLXU(:,:,K) =  FLXU_INT(:,:,K) + DFDTS   (:,:,K) * DELT
          if(associated(FLXAU)) FLXAU(:,:,K) = FLXAU_INT(:,:,K) + DFDTSNA (:,:,K) * DELT
          if(associated(FLCU))   FLCU(:,:,K) =  FLCU_INT(:,:,K) + DFDTSC  (:,:,K) * DELT
          if(associated(FLAU))   FLAU(:,:,K) =  FLAU_INT(:,:,K) + DFDTSCNA(:,:,K) * DELT
          ! downward fluxes
          if(associated(FLXD))   FLXD(:,:,K) =  FLXD_INT(:,:,K)
          if(associated(FLXAD)) FLXAD(:,:,K) = FLXAD_INT(:,:,K)
          if(associated(FLCD))   FLCD(:,:,K) =  FLCD_INT(:,:,K)
          if(associated(FLAD))   FLAD(:,:,K) =  FLAD_INT(:,:,K)
       end do

       ! fill TOA exports
       ! outgoing longwave radiation
       ! pmn: using FLXU_INT, etc. would be better ... here assuming down at TOA is zero
       if(associated(OLR  )) OLR   = -( FLX_INT(:,:, 0) + DFDTS   (:,:, 0) * DELT)
       if(associated(OLRA )) OLRA  = -(FLXA_INT(:,:, 0) + DFDTSNA (:,:, 0) * DELT)
       if(associated(OLC  )) OLC   = -( FLC_INT(:,:, 0) + DFDTSC  (:,:, 0) * DELT)
       if(associated(OLA  )) OLA   = -( FLA_INT(:,:, 0) + DFDTSCNA(:,:, 0) * DELT)
       if(associated(OLCC5)) then
                          where(CLDTT <= 0.05 )
                             OLCC5 = -( FLC_INT(:,:, 0) + DFDTSC  (:,:, 0) * DELT)
                          elsewhere
                             OLCC5 = MAPL_UNDEF
                          endwhere
       endif

       ! fill surface exports

       ! current surface emitted flux derivative wrt surface temperature (+ve)
       ! pmn: should be deprecated ... same as DSFDTS0
       if(associated(DSFDTS)) DSFDTS = -DFDTS(:,:,LM)

       ! surface emitted flux (+ve)
       if(associated(SFCEM)) SFCEM = SFCEM_INT - DFDTS(:,:,LM) * DELT

       ! absorbed (non-reflected) downward surface fluxes
       ! (remember: downward fluxes are not not linearized)
       if(associated(LWS  )) LWS   =  FLX_INT(:,:,LM) + SFCEM_INT  
       if(associated(LWSA )) LWSA  = FLXA_INT(:,:,LM) + SFCEM_INT  
       if(associated(LCS  )) LCS   =  FLC_INT(:,:,LM) + SFCEM_INT
       if(associated(LAS  )) LAS   =  FLA_INT(:,:,LM) + SFCEM_INT
       if(associated(LCSC5)) then
                          where(CLDTT <= 0.05 )
                             LCSC5 =  FLC_INT(:,:,LM) + SFCEM_INT
                          elsewhere
                             LCSC5 = MAPL_UNDEF
                          endwhere
       endif

       ! surface net downward fluxes
       if(associated(FLNS  )) FLNS   =  FLX_INT(:,:,LM) + DFDTS   (:,:,LM) * DELT
       if(associated(FLNSNA)) FLNSNA = FLXA_INT(:,:,LM) + DFDTSNA (:,:,LM) * DELT
       if(associated(FLNSC )) FLNSC  =  FLC_INT(:,:,LM) + DFDTSC  (:,:,LM) * DELT
       if(associated(FLNSA )) FLNSA  =  FLA_INT(:,:,LM) + DFDTSCNA(:,:,LM) * DELT

   ! RRTMG is a special case because its no-aerosol cases are missing
   else if( USE_RRTMG ) THEN

       ! fill 3D fluxes
       do K = 0, LM
          ! net downward (downward plus negated upward) fluxes
          if(associated(FLX))    FLX (:,:,K) =  FLX_INT(:,:,K) + DFDTS (:,:,K) * DELT ! all-sky
          if(associated(FLXA))  FLXA (:,:,K) = MAPL_UNDEF                             ! all-sky no-aerosol
          if(associated(FLC))    FLC (:,:,K) =  FLC_INT(:,:,K) + DFDTSC(:,:,K) * DELT ! clr-sky
          if(associated(FLA))    FLA (:,:,K) = MAPL_UNDEF                             ! clr-sky no-aerosol
          ! negated upward fluxes
          if(associated(FLXU))   FLXU(:,:,K) = FLXU_INT(:,:,K) + DFDTS (:,:,K) * DELT
          if(associated(FLXAU)) FLXAU(:,:,K) = MAPL_UNDEF
          if(associated(FLCU))   FLCU(:,:,K) = FLCU_INT(:,:,K) + DFDTSC(:,:,K) * DELT
          if(associated(FLAU))   FLAU(:,:,K) = MAPL_UNDEF
          ! downward fluxes
          if(associated(FLXD))   FLXD(:,:,K) = FLXD_INT(:,:,K)
          if(associated(FLXAD)) FLXAD(:,:,K) = MAPL_UNDEF
          if(associated(FLCD))   FLCD(:,:,K) = FLCD_INT(:,:,K)
          if(associated(FLAD))   FLAD(:,:,K) = MAPL_UNDEF
       end do

       ! fill TOA exports
       ! outgoing longwave radiation
       ! pmn: using FLXU_INT, etc. would be better ... here assuming down at TOA is zero
       if(associated(OLR  )) OLR   = -( FLX_INT(:,:, 0) + DFDTS (:,:, 0) * DELT)
       if(associated(OLRA )) OLRA  = MAPL_UNDEF
       if(associated(OLC  )) OLC   = -( FLC_INT(:,:, 0) + DFDTSC(:,:, 0) * DELT)
       if(associated(OLA  )) OLA   = MAPL_UNDEF
       if(associated(OLCC5)) then
                          where(CLDTT <= 0.05 )
                             OLCC5 = -( FLC_INT(:,:, 0) + DFDTSC(:,:, 0) * DELT)
                          elsewhere
                             OLCC5 = MAPL_UNDEF
                          endwhere
       endif

       ! fill surface exports

       ! current surface emitted flux derivative wrt surface temperature (+ve)
       ! pmn: should be deprecated ... same as DSFDTS0
       if(associated(DSFDTS)) DSFDTS = -DFDTS(:,:,LM)

       ! surface emitted flux (+ve)
       if(associated(SFCEM)) SFCEM = SFCEM_INT - DFDTS(:,:,LM) * DELT

       ! absorbed (non-reflected) downward surface fluxes
       ! (remember: downward fluxes are not not linearized)
       if(associated(LWS  )) LWS   = FLX_INT(:,:,LM) + SFCEM_INT  
       if(associated(LWSA )) LWSA  = MAPL_UNDEF
       if(associated(LCS  )) LCS   = FLC_INT(:,:,LM) + SFCEM_INT
       if(associated(LAS  )) LAS   = MAPL_UNDEF
       if(associated(LCSC5)) then
                          where(CLDTT <= 0.05 )
                             LCSC5 = FLC_INT(:,:,LM) + SFCEM_INT
                          elsewhere
                             LCSC5 = MAPL_UNDEF
                          endwhere
       endif

       ! surface net downward fluxes
       if(associated(FLNS  )) FLNS   = FLX_INT(:,:,LM) + DFDTS (:,:,LM) * DELT
       if(associated(FLNSNA)) FLNSNA = MAPL_UNDEF
       if(associated(FLNSC )) FLNSC  = FLC_INT(:,:,LM) + DFDTSC(:,:,LM) * DELT
       if(associated(FLNSA )) FLNSA  = MAPL_UNDEF

       ! band OLR and/or TBR output
       do ibnd = 1,nbndlw
          if (band_output(ibnd)) then

             write(bb,'(I0.2)') ibnd
             allocate(OLRB(IM,JM),__STAT__)

             ! get last full calculation
             call MAPL_GetPointer(INTERNAL, ptr2d, 'OLRB'//bb//'RG', __RC__)
             OLRB = ptr2d

             ! update for surface temperature on heartbeat
             call MAPL_GetPointer(INTERNAL, ptr2d, 'DOLRB'//bb//'RGDT', __RC__)
             OLRB = OLRB + ptr2d * DELT

             ! fill OLRBbbRG if requested
             call MAPL_GetPointer(EXPORT, ptr2d, 'OLRB'//bb//'RG', __RC__)
             if (associated(ptr2d)) ptr2D = OLRB

             ! calculate TBRBbbRG if requested
             call MAPL_GetPointer(EXPORT, ptr2d, 'TBRB'//bb//'RG', __RC__)
             if (associated(ptr2d)) then
                wn1 = wavenum1(ibnd)*100.; wn2 = wavenum2(ibnd)*100.  ! [m-1]
                call Tbr_from_band_flux(IM, JM, OLRB, wn1, wn2, ptr2d, __RC__)
             end if

             deallocate(OLRB)

          end if
       end do

   endif  ! RRTMG

   ! update reference linearization to current temperature
   ! pmn: should be deprecated because its moving along the line passing
   !   through point (TS_INT, SFCEM_INT) with slope -DFDTS (:,:,LM) that
   !   was defined only in the last REFRESH(). Better to just stick with
   !   the exports set in REFRESH() alone.
   if(associated(DSFDTS0)) DSFDTS0 =           - DFDTS(:,:,LM)
   if(associated(SFCEM0 )) SFCEM0  = SFCEM_INT - DFDTS(:,:,LM) * DELT
   if(associated(TSREFF )) TSREFF  = TSINST

!  All done
!-----------
   deallocate( DUMTT )

   RETURN_(ESMF_SUCCESS)

 end subroutine Update_Flx

 ! estimate brightness temperature from a band flux
 subroutine Tbr_from_band_flux(IM, JM, Fband_, wn1, wn2, Tbr_, RC)

   ! input arguments
   integer, intent(in ) :: IM, JM
   real,    intent(in ) :: Fband_(IM,JM) ! band flux [W/m2]
   real,    intent(in ) :: wn1, wn2      ! bounds of band [m-1]

   ! output arguments
   real,    intent(out) :: Tbr_(IM,JM)   ! brightness temp [K]

   ! error code
   integer, optional, intent(out) :: RC

   ! fundamental constants
   double precision, parameter :: h  = 6.626070040d-34  ! Plancks constant         [J.s]
   double precision, parameter :: c  = 2.99792458d8     ! Speed of light in vacuum [m/s]
   double precision, parameter :: kB = 1.38064852d-23   ! Boltzmann constant       [J/K]
   double precision, parameter :: pi = MAPL_PI_R8

   ! other constants
   double precision, parameter :: alT = h * c / kB
   double precision, parameter :: bigS = 2.0d0 * kB**4 * pi / (h**3 * c**2)
   double precision, parameter :: bigC = 2.0d0 * h * c**2

   ! locals
   double precision, dimension(IM,JM) :: Fband, Tbr, Bmean
   real :: wnMid

   if (present(RC)) RC = ESMF_SUCCESS

   ! special case of all zero fluxes before first call to LW_Driver()
   if (all(Fband_ == 0.0)) then
     Tbr_ = MAPL_UNDEF
     return
   end if

   ! calculations done in double precision
   Fband = dble(Fband_)

   ! first guess Tbr from narrow band approximation ...
   ! (1) estimate mean Planck function for a narrow band
   Bmean = Fband / (pi * (wn2 - wn1))
   ! (2) invert Planck function for temp at mid-point wavenumber
   wnMid = (wn1 + wn2) / 2.0d0
   call invert_Planck_for_T(IM, JM, Bmean, wnMid, bigC, alT, Tbr, __RC__)

   ! now refine with a wide band esimate
   ! PMN: Iterative routine not ready for prime time
   !      Produces erroneously large Tbr in cloudy regions
   !call Tbr_wide_band(IM, JM, Fband, wn1, wn2, bigS, alT, Tbr, __RC__)

   ! put output back in real
   Tbr_ = real(Tbr)

 end subroutine Tbr_from_band_flux

 ! invert Planck function for temperature
 subroutine invert_Planck_for_T(IM, JM, Bwn, wn, bigC, alT, T, RC)

   ! input arguments
   integer,          intent(in ) :: IM, JM
   double precision, intent(in ) :: Bwn(IM,JM)  ! PlanckFn(wavenumber)
   real,             intent(in ) :: wn          ! wavenumber [m-1]
   double precision, intent(in ) :: bigC, alT   ! necessary constants

   ! output arguments
   double precision, intent(out) :: T(IM,JM)    ! temperature [K]

   ! error code
   integer, optional, intent(out) :: RC

   ! error checking
   if (present(RC)) RC = ESMF_SUCCESS
   _ASSERT(wn > 0.,'needs informative message')

   ! invert Planck function for temp
   T = alT * wn / log(bigC * wn**3 / Bwn + 1.0d0)

 end subroutine invert_Planck_for_T

 ! Tbr from wide band approximation
 subroutine Tbr_wide_band(IM, JM, Fband, wn1, wn2, bigS, alT, Tbr, RC)

   ! input arguments
   integer,          intent(in   ) :: IM, JM
   double precision, intent(in   ) :: Fband(IM,JM)  ! band flux [W/m2]
   real,             intent(in   ) :: wn1, wn2      ! bounds of band [m-1]
   double precision, intent(in   ) :: bigS, alT     ! necessary constant

   ! Tbr inputs first guess and outputs better estimate
   double precision, intent(inout) :: Tbr(IM,JM)    ! brightness temp [K]

   ! error code
   integer, optional, intent(out) :: RC

   ! number of iterations for wide band estimate (converges slowly)
   integer, parameter :: Nits = 16

   ! locals
   integer :: n
   real    :: alTwn1, alTwn2

   ! error checking
   if (present(RC)) RC = ESMF_SUCCESS
   _ASSERT(wn2 > wn1,'needs informative message')
   _ASSERT(Nits >= 1,'needs informative message')

   ! iterate from first guess Tbr to better estimate
   alTwn1 = alT * wn1
   alTwn2 = alT * wn2
   do n = 1, Nits
    Tbr = ( Fband / (bigS * (Tfunc(alTwn1/Tbr) - Tfunc(alTwn2/Tbr))) ) ** 0.25d0
   end do

 end subroutine Tbr_wide_band

 elemental double precision function Tfunc(x)

   double precision, intent(in) :: x

   ! maximum number of terms in series (converges quickly)
   integer, parameter :: nmax = 4

   ! locals
   integer :: n, n2, n3
   double precision :: emx, cx0, cx1, cx2, cx3, zn

   ! setup
   emx = exp(-x)
   cx0 = 6.0d0
   cx1 = 6.0d0 * x
   cx2 = 3.0d0 * x**2
   cx3 =         x**3

   ! do at least 1st order
   Tfunc = (cx3 + cx2 + cx1 + cx0) * emx
   if (nmax <= 1) return

   ! higher orders
   zn = emx
   do n = 2, nmax
     n2 = n * n
     n3 = n * n2
     zn = zn * emx
     Tfunc = Tfunc + (cx3 + cx2/n + cx1/n2 + cx0/n3) * zn / n
   end do

 end function Tfunc

end subroutine RUN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module GEOS_IrradGridCompMod

