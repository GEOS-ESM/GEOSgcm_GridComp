#include "MAPL_Generic.h"


!=============================================================================
module GEOS_IgniGridCompMod

!BOP

! !MODULE: GEOS_Fires -- child to the "Land" gridded component.  

!DESCRIPTION:
!   {\tt GEOS\_Vegdyn} is a gridded component that performs the
!   necessary interpolation to provide refreshed values of the 
!   dynamic vegetation values prescribed by external data/observations.\\
!
! There are no imports to this routine.
! Exports from this routine are the instaneous values of the
! vegetation parameters on tilespace to be used in other components
! of the land subroutine.  All exports and imports are stored on the
! tile grid inherited from the parent routine.\\
! 
! I. Parameter Class 1: Time AND spatially dependent parameters 
! from a binary data file\\
! 
! Current list: LAI, GRN, NDVI \\
! 
! The gridded component stores the surrounding observations of 
! each parameter in the internal state.  If the run method 
! discovers that the current internal state does not contain the 
! observed values required to interpolate the values at the current 
! time, it performs the required i/o to refresh the values of 
! the internal state.  The first iteration of the run method 
! always has to fill the values.  No restart is required by this 
! gridded component for these parameters.  (A restart *is* now
! required for Vegetation Class 3 \\
!
! INTERNALS: \\
!
! EXPORTS:  \\
!
! !USES:

  use ESMF
  use MAPL

  use GEOS_UtilsMod, only: GEOS_QSAT

  use cffwi, only: fine_fuel_moisture_code, duff_moisture_code,       &
                   drought_code, initial_spread_index, buildup_index, &
                   fire_weather_index, daily_severity_rating,         &
                   FFMC_INIT, DMC_INIT, DC_INIT

  implicit none
  private

! !PUBLIC MEMBER FUNCTIONS:

  public SetServices

!EOP

  integer :: NUM_ENSEMBLE
  integer, parameter :: NTYPS = MAPL_NumVegTypes

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
!
! ErrLog Variables


    character(len=ESMF_MAXSTR)              :: IAm
    integer                                 :: STATUS
    character(len=ESMF_MAXSTR)              :: COMP_NAME

! Local derived type aliases

    type(MAPL_MetaComp),pointer             :: MAPL=>null()

!=============================================================================

! Begin...

!------------------------------------------------------------
! Get my name and set-up traceback handle
!------------------------------------------------------------

    call ESMF_GridCompGet(GC, NAME=COMP_NAME, __RC__)

    Iam = trim(COMP_NAME) // 'SetServices'

! -----------------------------------------------------------
! Set the Run entry point
! -----------------------------------------------------------

    call MAPL_GridCompSetEntryPoint(GC, ESMF_METHOD_RUN, Run, __RC__)

! -----------------------------------------------------------
! Get the configuration
! -----------------------------------------------------------
    call MAPL_GetObjectFromGC(GC, MAPL, __RC__)

! -----------------------------------------------------------
! Get the intervals
! -----------------------------------------------------------
    call MAPL_GetResource ( MAPL, NUM_ENSEMBLE, Label="NUM_LDAS_ENSEMBLE:", DEFAULT=1, __RC__)

    !call MAPL_GetResource ( MAPL,DT, Label="RUN_DT:", RC=STATUS)
    !VERIFY_(STATUS)

    !RUN_DT = nint(DT)

! -----------------------------------------------------------
! At the moment, this will refresh when the land parent 
! needs to refresh.
!
!    call ESMF_ConfigGetFloat ( CF, DT, Label=trim(COMP_NAME)//&
!    "_DT:", default=DT, RC=STATUS)
!     VERIFY_(STATUS)
!
!    MY_STEP = nint(DT)
!
! -----------------------------------------------------------


! -----------------------------------------------------------
! Set the state variable specs.
! -----------------------------------------------------------

!BOS

! -----------------------------------------------------------
!   Import States
! -----------------------------------------------------------
    call MAPL_AddImportSpec(GC,                    &
         SHORT_NAME = 'MOT2M',                     & 
         LONG_NAME  = 'temperature 2m from MO sfc',&
         UNITS      = 'K',                         &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone,          &
         RESTART    = MAPL_RestartSkip, __RC__)

    call MAPL_AddImportSpec(GC,                    &
         SHORT_NAME = 'MOQ2M',                     & 
         LONG_NAME  = 'humidity 2m from MO sfc',   &
         UNITS      = 'kg kg-1',                   &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone,          &
         RESTART    = MAPL_RestartSkip, __RC__)

    call MAPL_AddImportSpec(GC,                    &
         SHORT_NAME = 'MOU10M',                    & 
         LONG_NAME  = 'zonal 10m wind from MO sfc',&
         UNITS      = 'm s-1',                     &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone,          &
         RESTART    = MAPL_RestartSkip, __RC__)

    call MAPL_AddImportSpec(GC,                    &
         SHORT_NAME = 'MOV10M',                    & 
         LONG_NAME  = 'meridional 10m wind from MO sfc', &
         UNITS      = 'm s-1',                     &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone,          &
         RESTART    = MAPL_RestartSkip, __RC__)

#if (1)
    call MAPL_AddImportSpec(GC,                    &
         SHORT_NAME = 'PRLAND',                    & 
         LONG_NAME  = 'total_precipitation_land',  &
         UNITS      = 'kg m-2 s-1',                &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone,          &
         RESTART    = MAPL_RestartOptional, __RC__)
#else
    call MAPL_AddImportSpec(GC,                    &
         LONG_NAME  = 'liquid_water_convective_precipitation', &
         UNITS      = 'kg m-2 s-1',                &
         SHORT_NAME = 'PCU',                       &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone,          &
         RESTART    = MAPL_RestartOptional, __RC__)

    call MAPL_AddImportSpec(GC,                    &
         LONG_NAME  = 'liquid_water_large_scale_precipitation',&
         UNITS      = 'kg m-2 s-1',                &
         SHORT_NAME = 'PLS',                       &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone,           &
         RESTART    = MAPL_RestartOptional, __RC__)

    call MAPL_AddImportSpec(GC,                    &
         LONG_NAME  = 'snowfall',                  &
         UNITS      = 'kg m-2 s-1',                &
         SHORT_NAME = 'SNO',                       &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone,          &
         RESTART    = MAPL_RestartOptional, __RC__) 
#endif

    call MAPL_AddImportSpec(GC,                    &
         LONG_NAME  = 'surface_pressure',          &
         UNITS      = 'Pa',                        &
         SHORT_NAME = 'PS',                        &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone,          &
         RESTART    = MAPL_RestartSkip, __RC__)


! -----------------------------------------------------------
! Internal State 
! -----------------------------------------------------------
    
    call MAPL_AddInternalSpec(GC,                  &
         SHORT_NAME = 'FFMC',                      &
         LONG_NAME  = 'fine_fuel_moisture_code',   &
         UNITS      = '1',                         &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone,          &
         ADD2EXPORT = .true.,                      &
         RESTART    = MAPL_RestartOptional,        &
         DEFAULT    = FFMC_INIT, __RC__)

    call MAPL_AddInternalSpec(GC,                  &
         SHORT_NAME = 'DMC',                       &
         LONG_NAME  = 'duff_moisture_code',        &
         UNITS      = '1',                         &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone,          &
         ADD2EXPORT = .true.,                      &
         RESTART    = MAPL_RestartOptional,        &
         DEFAULT    = DMC_INIT, __RC__)

   call MAPL_AddInternalSpec(GC,                   &
         SHORT_NAME = 'DC',                        &
         LONG_NAME  = 'drought_code',              &
         UNITS      = '1',                         &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone,          &
         ADD2EXPORT = .true.,                      &
         RESTART    = MAPL_RestartOptional,        &
         DEFAULT    = DC_INIT, __RC__) 

    call MAPL_AddInternalSpec(GC,                  &
         SHORT_NAME = 'PRA',                       &
         LONG_NAME  = 'precip_since_sun_noon',     &
         UNITS      = 'kg m-2',                    &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone,          &
         ADD2EXPORT = .true.,                      &
!        PRECISION  = MAPL_R4,                     &
         RESTART    = MAPL_RestartOptional,        &
         DEFAULT    = 0.0, __RC__)

    call MAPL_AddInternalSpec(GC,                  &
         SHORT_NAME = 'PR',                        &
         LONG_NAME  = 'precip_at_sun_noon',        &
         UNITS      = 'kg m-2',                    &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone,          &
         ADD2EXPORT = .true.,                      &
         RESTART    = MAPL_RestartOptional,        &
         DEFAULT    = 0.0, __RC__)


    call MAPL_AddInternalSpec(GC,                  &
         SHORT_NAME = 'T',                         &
         LONG_NAME  = 'temperature_at_sun_noon',   &
         UNITS      = 'K',                         &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone,          &
         ADD2EXPORT = .true.,                      &
         RESTART    = MAPL_RestartOptional,        &
         DEFAULT    = 0.0, __RC__)

    call MAPL_AddInternalSpec(GC,                  &
         SHORT_NAME = 'RH',                        &
         LONG_NAME  = 'humidity_at_sun_noon',      &
         UNITS      = '1',                         &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone,          &
         ADD2EXPORT = .true.,                      &
         RESTART    = MAPL_RestartOptional,        &
         DEFAULT    = 0.50, __RC__)

    call MAPL_AddInternalSpec(GC,                  &
         SHORT_NAME = 'WIND',                      &
         LONG_NAME  = 'wind_at_sun_noon',          &
         UNITS      = 'm s-1',                     &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone,          &
         ADD2EXPORT = .true.,                      &
         RESTART    = MAPL_RestartOptional,        &
         DEFAULT    = 0.0, __RC__)


! -----------------------------------------------------------
! Export Variables
! -----------------------------------------------------------
    call MAPL_AddExportSpec(GC,                    & 
         SHORT_NAME = 'FWI',                       &
         LONG_NAME  = 'fire_weather_index',        &
         UNITS      = '1',                         &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone, __RC__) 

    call MAPL_AddExportSpec(GC,                    & 
         SHORT_NAME = 'BUI',                       &
         LONG_NAME  = 'buildup_index',             &
         UNITS      = '1',                         &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone, __RC__) 

    call MAPL_AddExportSpec(GC,                    &
         SHORT_NAME = 'ISI',                       &
         LONG_NAME  = 'initial_spread_index',      &
         UNITS      = '1',                         &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone, __RC__) 

   call MAPL_AddExportSpec(GC,                     &
         SHORT_NAME = 'DSR',                       &
         LONG_NAME  = 'daily_severity_rating',     &
         UNITS      = '1',                         &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone, __RC__)



    call MAPL_AddExportSpec(GC,                    & 
         SHORT_NAME = 'SNFFMC',                    &
         LONG_NAME  = 'fine_fuel_moisture_code',   &
         UNITS      = '1',                         &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone, __RC__) 

    call MAPL_AddExportSpec(GC,                    & 
         SHORT_NAME = 'SNDMC',                     &
         LONG_NAME  = 'duff_moisture_code',        &
         UNITS      = '1',                         &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone, __RC__) 

    call MAPL_AddExportSpec(GC,                    & 
         SHORT_NAME = 'SNDC',                      &
         LONG_NAME  = 'drought_code',              &
         UNITS      = '1',                         &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone, __RC__) 

    call MAPL_AddExportSpec(GC,                    & 
         SHORT_NAME = 'SNFWI',                     &
         LONG_NAME  = 'fire_weather_index',        &
         UNITS      = '1',                         &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone, __RC__) 

    call MAPL_AddExportSpec(GC,                    & 
         SHORT_NAME = 'SNBUI',                     &
         LONG_NAME  = 'buildup_index',             &
         UNITS      = '1',                         &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone, __RC__) 

    call MAPL_AddExportSpec(GC,                    &
         SHORT_NAME = 'SNISI',                     &
         LONG_NAME  = 'initial_spread_index',      &
         UNITS      = '1',                         &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone, __RC__) 

   call MAPL_AddExportSpec(GC,                     &
         SHORT_NAME = 'SNDSR',                     &
         LONG_NAME  = 'daily_severity_rating',     &
         UNITS      = '1',                         &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone, __RC__)


    call MAPL_AddExportSpec(GC,                    & 
         SHORT_NAME = 'VPD',                       &
         LONG_NAME  = 'vapor_pressure_defficit',   &
         UNITS      = '1',                         &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone, __RC__)


     call MAPL_AddExportSpec(GC,                   & 
         SHORT_NAME = 'DBG1',                      &
         LONG_NAME  = 'DEBUG',                     &
         UNITS      = '1',                         &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone, __RC__)

     call MAPL_AddExportSpec(GC,                   & 
         SHORT_NAME = 'DBG2',                      &
         LONG_NAME  = 'DEBUG',                     &
         UNITS      = '1',                         &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone, __RC__)

     call MAPL_AddExportSpec(GC,                   & 
         SHORT_NAME = 'DBG3',                      &
         LONG_NAME  = 'DEBUG',                     &
         UNITS      = '1',                         &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone, __RC__)

!EOS


! Set the Profiling timers
! ------------------------
    call MAPL_TimerAdd(GC, name='TOTAL'     , __RC__)
!   call MAPL_TimerAdd(GC, name='INITIALIZE', __RC__)
    call MAPL_TimerAdd(GC, name='RUN'       , __RC__)
    call MAPL_TimerAdd(GC, name='-CFFWI'    , __RC__)
!   call MAPL_TimerAdd(GC, name='FINALIZE'  , __RC__)



!------------------------------------------------------------
! Set generic init and final methods
!------------------------------------------------------------

    call MAPL_GenericSetServices(GC, __RC__)

    RETURN_(ESMF_SUCCESS)
  
  end subroutine SetServices

! -----------------------------------------------------------
! RUN -- Run method for the vegdyn component
! -----------------------------------------------------------

  subroutine RUN (GC, IMPORT, EXPORT, CLOCK, RC )

! -----------------------------------------------------------
! !ARGUMENTS:

    type(ESMF_GridComp), intent(inout) :: GC    
    type(ESMF_State),    intent(inout) :: IMPORT
    type(ESMF_State),    intent(inout) :: EXPORT
    type(ESMF_Clock),    intent(inout) :: CLOCK
    integer, optional,   intent(  out) :: RC

!EOP

! ErrLog Variables

    character(len=ESMF_MAXSTR)         :: Iam
    integer                            :: STATUS
    character(len=ESMF_MAXSTR)         :: COMP_NAME

! Locals

    type (MAPL_MetaComp), pointer      :: MAPL => null()
    type (ESMF_State)                  :: INTERNAL

    type(ESMF_TimeInterval)            :: DELT 
    type(MAPL_SunOrbit)                :: ORBIT
    real, allocatable, dimension(:)    :: SLI
    real, allocatable, dimension(:)    :: ZTH
    real, allocatable, dimension(:)    :: ZTH_


! IMPORT pointers
    real, dimension(:), pointer :: T2M
    real, dimension(:), pointer :: Q2M

    real, dimension(:), pointer :: U10M
    real, dimension(:), pointer :: V10M

    real, dimension(:), pointer :: PS
#if (1)
    real, dimension(:), pointer :: PRLAND
#else
    real, dimension(:), pointer :: PCU
    real, dimension(:), pointer :: PLS
    real, dimension(:), pointer :: SNO
#endif

! INTERNAL pointers    
    real, dimension(:), pointer :: FFMC_
    real, dimension(:), pointer :: DMC_
    real, dimension(:), pointer :: DC_

    real(kind=MAPL_R4), dimension(:), pointer :: PRA_
    real, dimension(:), pointer :: PR_
    real, dimension(:), pointer :: RH_
    real, dimension(:), pointer :: T_
    real, dimension(:), pointer :: WIND_
     

! EXPORT pointers
    real, dimension(:), pointer :: FFMC
    real, dimension(:), pointer :: DMC
    real, dimension(:), pointer :: DC
    real, dimension(:), pointer :: ISI
    real, dimension(:), pointer :: FWI
    real, dimension(:), pointer :: BUI
    real, dimension(:), pointer :: DSR

    real, dimension(:), pointer :: SNFFMC
    real, dimension(:), pointer :: SNDMC
    real, dimension(:), pointer :: SNDC
    real, dimension(:), pointer :: SNISI
    real, dimension(:), pointer :: SNFWI
    real, dimension(:), pointer :: SNBUI
    real, dimension(:), pointer :: SNDSR

    real, dimension(:), pointer :: PR

    real, dimension(:), pointer :: VPD

    real, dimension(:), pointer :: DBG1, DBG2, DBG3
  
! Others
    type(ESMF_Time)         :: time
    type(ESMF_Alarm)        :: run_alarm
    type(ESMF_TimeInterval) :: ring_interval
    real(ESMF_KIND_R8)      :: time_step

    integer :: NT

    integer :: year, month, day, hr, mn, sc
    real    :: dt

    real, pointer, dimension(:) :: LATS => null()
    real, pointer, dimension(:) :: LONS => null()

    real, allocatable, dimension(:) :: tmpFFMC
    real, allocatable, dimension(:) :: tmpDMC
    real, allocatable, dimension(:) :: tmpDC
    real, allocatable, dimension(:) :: tmpFWI
    real, allocatable, dimension(:) :: tmpISI
    real, allocatable, dimension(:) :: tmpBUI
    real, allocatable, dimension(:) :: tmpDSR

    real, allocatable, dimension(:) :: solar_time
    logical, allocatable, dimension(:) :: sun_noon


    real, parameter :: f_e = (24*3600)/(2*MAPL_PI)
    real, parameter :: sun_noon_utc = 12*3600.0


    real,    allocatable, dimension(:) :: LSHA0, LSHA1
    logical, allocatable, dimension(:) :: isNoon

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    call ESMF_GridCompGet(GC, name=COMP_NAME, __RC__)
  
    Iam = trim(COMP_NAME) // "Run"

! Get my internal MAPL_Generic state
! -----------------------------------------------------------

    call MAPL_GetObjectFromGC(GC, MAPL, __RC__)

    call MAPL_TimerOn(MAPL, 'TOTAL')
    call MAPL_TimerOn(MAPL, 'RUN'  )

    call MAPL_Get(MAPL, TILELATS=LATS, &
                        TILELONS=LONS, &
                        INTERNAL_ESMF_STATE=INTERNAL, __RC__)

! -----------------------------------------------------------
! Get file names from configuration
! -----------------------------------------------------------

    if(NUM_ENSEMBLE > 1) then
        !comp_name should be vegdynxxxx....
    endif

! get pointers to internal variables
! ----------------------------------
    call MAPL_GetPointer(INTERNAL, FFMC_,  'FFMC',    __RC__)
    call MAPL_GetPointer(INTERNAL, DMC_,   'DMC',     __RC__)
    call MAPL_GetPointer(INTERNAL, DC_,    'DC',      __RC__)
    call MAPL_GetPointer(INTERNAL, PR_,    'PR',      __RC__)
    call MAPL_GetPointer(INTERNAL, PRA_,   'PRA',     __RC__)
    call MAPL_GetPointer(INTERNAL, T_,     'T',       __RC__)
    call MAPL_GetPointer(INTERNAL, RH_,    'RH',      __RC__)
    call MAPL_GetPointer(INTERNAL, WIND_,  'WIND',    __RC__)

! get pointers to internal variables
! ----------------------------------
    call MAPL_GetPointer(IMPORT,   Q2M,    'MOQ2M',   __RC__)
    call MAPL_GetPointer(IMPORT,   T2M,    'MOT2M',   __RC__)

    call MAPL_GetPointer(IMPORT,   U10M,   'MOU10M',  __RC__)
    call MAPL_GetPointer(IMPORT,   V10M,   'MOV10M',  __RC__)

    call MAPL_GetPointer(IMPORT,     PS,   'PS',      __RC__)
#if (1)
    call MAPL_GetPointer(IMPORT, PRLAND,   'PRLAND',  __RC__)
#else
    call MAPL_GetPointer(IMPORT,    PCU,   'PCU',     __RC__)
    call MAPL_GetPointer(IMPORT,    PLS,   'PLS',     __RC__)
    call MAPL_GetPointer(IMPORT,    SNO,   'SNO',     __RC__)
#endif
    

! get pointers to EXPORTS
! -----------------------
    call MAPL_GetPointer(EXPORT,   FFMC,   'FFMC',    __RC__)
    call MAPL_GetPointer(EXPORT,   DMC,    'DMC',     __RC__)
    call MAPL_GetPointer(EXPORT,   DC,     'DC',      __RC__)
    call MAPL_GetPointer(EXPORT,   FWI,    'FWI',     __RC__)
    call MAPL_GetPointer(EXPORT,   ISI,    'ISI',     __RC__)
    call MAPL_GetPointer(EXPORT,   BUI,    'BUI',     __RC__)
    call MAPL_GetPointer(EXPORT,   DSR,    'DSR',     __RC__)

    call MAPL_GetPointer(EXPORT,   SNFFMC, 'SNFFMC',  __RC__)
    call MAPL_GetPointer(EXPORT,   SNDMC,  'SNDMC',   __RC__)
    call MAPL_GetPointer(EXPORT,   SNDC,   'SNDC',    __RC__)
    call MAPL_GetPointer(EXPORT,   SNFWI,  'SNFWI',   __RC__)
    call MAPL_GetPointer(EXPORT,   SNISI,  'SNISI',   __RC__)
    call MAPL_GetPointer(EXPORT,   SNBUI,  'SNBUI',   __RC__)
    call MAPL_GetPointer(EXPORT,   SNDSR,  'SNDSR',   __RC__)


    call MAPL_GetPointer(EXPORT,   PR,    'PR',       __RC__)

    call MAPL_GetPointer(EXPORT,   VPD,   'VPD',      __RC__)

    call MAPL_GetPointer(EXPORT,   DBG1,  'DBG1',  alloc=.true.,   __RC__)
    call MAPL_GetPointer(EXPORT,   DBG2,  'DBG2',  alloc=.true.,   __RC__)
    call MAPL_GetPointer(EXPORT,   DBG3,  'DBG3',  alloc=.true.,   __RC__)


! Fire weather indexes in the CFFWI
! ---------------------------------
    NT = SIZE(LONS)

    call MAPL_TimerOn(MAPL, '-CFFWI')

    if (NT > 0) then

        call MAPL_Get(MAPL, RunAlarm=run_alarm, __RC__)
        call ESMF_AlarmGet(run_alarm, ringInterval=ring_interval, __RC__)

        call ESMF_TimeIntervalGet(ring_interval, s_r8=time_step, __RC__)
        dt = real(time_step)

        allocate(tmpFFMC(NT), tmpDMC(NT), tmpDC(NT), tmpFWI(NT), &
                 tmpISI(NT),  tmpBUI(NT), tmpDSR(NT), __STAT__)

        tmpFFMC = -MAPL_UNDEF
        tmpDMC  = -MAPL_UNDEF
        tmpDC   = -MAPL_UNDEF

        tmpFWI  = -MAPL_UNDEF
        tmpISI  = -MAPL_UNDEF
        tmpBUI  = -MAPL_UNDEF
        tmpDSR  = -MAPL_UNDEF


        ! accumulate precip
#if (1)
        PRA_ = PRA_ + PRLAND*dt
#else
        PRA_ = PRA_ + (PCU+PLS+SNO)*dt
#endif

        call ESMF_ClockGet(CLOCK, currTime=time, __RC__)
        call ESMF_TimeGet(time, yy=year, mm=month, dd=day, h=hr, m=mn, s=sc, __RC__)

        allocate(solar_time(NT), sun_noon(NT), __STAT__)
        solar_time = (hr*3600 + mn*60 + sc) + f_e*LONS
        where (solar_time > 24*3600) solar_time = solar_time - 24*3600

        sun_noon = .false.
        where (((solar_time-sun_noon_utc) > 0) .and. ((solar_time-sun_noon_utc) <= dt))
            sun_noon = .true.
            T_       = T2M
            WIND_    = sqrt(U10M*U10M + V10M*V10M)
            RH_      = min(Q2M / GEOS_QSAT(T2M, PS, PASCALS=.true.), 1.0)
            PR_      = PRA_
            PRA_     = dble(0.0)
        end where

        tmpFFMC = FFMC_
        tmpDMC  = DMC_
        tmpDC   = DC_


        call cffwi_driver(tmpFFMC, tmpDMC, tmpDC, tmpISI, tmpBUI, tmpFWI, tmpDSR, &
                          T_, RH_, WIND_, PR_, month, NT)

        ! update internal state
        where (sun_noon)
            FFMC_ = tmpFFMC
            DMC_  = tmpDMC
            DC_   = tmpDC
        end where

        ! update exports
        if (associated(FFMC)) FFMC = tmpFFMC
        if (associated(DMC))   DMC = tmpDMC
        if (associated(DC))     DC = tmpDC
        if (associated(ISI))   ISI = tmpISI
        if (associated(BUI))   BUI = tmpBUI
        if (associated(FWI))   FWI = tmpFWI
        if (associated(DSR))   DSR = tmpDSR
        if (associated(PR))    PR  = PR_

        where (.not. sun_noon)
            tmpFFMC = MAPL_UNDEF
            tmpDMC  = MAPL_UNDEF
            tmpDC   = MAPL_UNDEF
            tmpISI  = MAPL_UNDEF
            tmpBUI  = MAPL_UNDEF
            tmpFWI  = MAPL_UNDEF
            tmpDSR  = MAPL_UNDEF
        end where

        if (associated(SNFFMC)) SNFFMC = tmpFFMC
        if (associated(SNDMC))   SNDMC = tmpDMC
        if (associated(SNDC))     SNDC = tmpDC
        if (associated(SNISI))   SNISI = tmpISI
        if (associated(SNBUI))   SNBUI = tmpBUI
        if (associated(SNFWI))   SNFWI = tmpFWI
        if (associated(SNDSR))   SNDSR = tmpDSR


        DBG1 = MAPL_UNDEF 
        where (sun_noon) DBG1 = 1.0

        DBG2 = solar_time
 
        deallocate(tmpFFMC, tmpDMC, tmpDC, tmpFWI, tmpISI, tmpBUI, tmpDSR, __STAT__)


!!! test 
        call MAPL_Get(MAPL, ORBIT=ORBIT, __RC__)

        allocate(LSHA0(NT), LSHA1(NT), isNoon(NT), __STAT__)

        call ESMF_ClockGet(CLOCK, CURRTIME=time, __RC__)
        call MAPL_SunGetLocalSolarHourAngle(ORBIT, LONS, LSHA0, TIME=time, __RC__)
        call MAPL_SunGetLocalSolarHourAngle(ORBIT, LONS, LSHA1, TIME=time+ring_interval, __RC__)
        isNoon = (LSHA0 <= 0.0) .and. (LSHA1 > 0.0)

        DBG3 = MAPL_UNDEF 
        where (isNoon) DBG3 = 1.0

        deallocate(LSHA0, LSHA1, isNoon, __STAT__)
!!! test

    end if

    call MAPL_TimerOff(MAPL, '-CFFWI')

!  All done
! ---------
    call MAPL_TimerOff(MAPL, 'RUN'  )
    call MAPL_TimerOff(MAPL, 'TOTAL')

    RETURN_(ESMF_SUCCESS)
  end subroutine RUN



  subroutine cffwi_driver(ffmc, dmc, dc, isi, bui, fwi, dsr, &
                          T, RH, wind, Pr, month, N)
 
    !
    ! Calculates FFMC, DMC, DC, ISI, BUI, FWI and DSR indexes.
    !

    implicit none

    real,    dimension(N), intent(in) :: T, RH, wind, Pr
    integer,               intent(in) :: month
    integer,               intent(in) :: N

    real, dimension(N), intent(inout) :: ffmc, dmc, dc, isi, bui, fwi, dsr
    
    ! local
    integer :: i
    real    :: T_, RH_, Pr_

    do i = 1, N
        T_  = T(i) - 273.15    ! temperature, C
        RH_ = 100 * RH(i)      ! relative humidity, %
        Pr_ = Pr(i)            ! precip, mm

        ffmc(i) = fine_fuel_moisture_code(ffmc(i), T_, RH_, wind(i), Pr_)
        dmc(i)  = duff_moisture_code(dmc(i), T_, RH_, Pr_, month)
        dc(i)   = drought_code(dc(i), T_, Pr_, month)

        isi(i)  = initial_spread_index(ffmc(i), wind(i))
        bui(i)  = buildup_index(dmc(i), dc(i))

        fwi(i)  = fire_weather_index(isi(i), bui(i))
        dsr(i)  = daily_severity_rating(fwi(i))
    end do

  end subroutine cffwi_driver


end module GEOS_IgniGridCompMod
