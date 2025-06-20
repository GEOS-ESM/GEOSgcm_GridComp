#include "MAPL_Generic.h"


!=============================================================================
module GEOS_IgniGridCompMod

!BOP

! !MODULE: GEOS_Igni -- Implements fire weather and fire danger.  

!DESCRIPTION:
!   {\tt GEOS\_Igni} is a gridded component that integrates fire weather 
!   observations into fuel moisture codes and fire behavior indexes 
!   based on the Canadian FWI System. \\
!
!   {\tt GEOS\_Igni} includes hourly and daily variants of the FWI system that
!   predict Fine Fuel Moisture Code (FFMC), Duff Moisture Code (DMC),
!   Drought Code (DC), Initial Spread Index (ISI), Buildup Index (BUI),
!   Fire Weather Index (FWI) and Daily Severity Rating (DSR). These are 
!   calculated on the land tiles. \\
!

! !USES:

  use ESMF
  use MAPL

  use GEOS_UtilsMod, only: GEOS_QSAT

  use cffwi, only: fine_fuel_moisture_code,      &
                   grass_fuel_moisture_code,     & 
                   duff_moisture_code,           &
                   drought_code,                 &
                   initial_spread_index,         &
                   buildup_index,                &
                   fire_weather_index,           &
                   daily_severity_rating,        &
                   FFMC_INIT, DMC_INIT, DC_INIT, &
                   NOMINAL_FINE_FUEL_LOAD

  implicit none
  private

! !PUBLIC MEMBER FUNCTIONS:

  public SetServices


! ! Private state

  integer :: NUM_ENSEMBLE

  integer, parameter :: LOCAL_NOON_SOLAR = 0
  integer, parameter :: LOCAL_NOON_LST   = 1

  type IGNI_State
     private
     integer :: local_noon_method = LOCAL_NOON_SOLAR
  end type IGNI_State

  type wrap_
     type (IGNI_State), pointer :: PTR => null()
  end type wrap_


!EOP

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

! ErrLog Variables

    character(len=ESMF_MAXSTR)   :: Iam
    integer                      :: STATUS
    character(len=ESMF_MAXSTR)   :: COMP_NAME


! Local derived type aliases

    type(MAPL_MetaComp), pointer :: MAPL=>null()

    type (IGNI_State), pointer   :: self
    type (wrap_)                 :: wrap


! Local
    
    real :: run_dt
    real :: dt
    
    character(len=ESMF_MAXSTR) :: resource_file
    type(ESMF_Config)          :: config
    character(len=ESMF_MAXSTR) :: local_noon


!=============================================================================

! Begin...

!------------------------------------------------------------
! Get my name and set-up traceback handle
!------------------------------------------------------------

    call ESMF_GridCompGet(GC, NAME=COMP_NAME, __RC__)

    Iam = trim(COMP_NAME) // 'SetServices'


! -----------------------------------------------------------
! Wrap private internal state for storing in GC
! -----------------------------------------------------------
    allocate (self, __STAT__)
    wrap%ptr => self

! -----------------------------------------------------------
! Set the Run entry point
! -----------------------------------------------------------

    call MAPL_GridCompSetEntryPoint(GC, ESMF_METHOD_RUN, Run1, __RC__)
    call MAPL_GridCompSetEntryPoint(GC, ESMF_METHOD_RUN, Run2, __RC__)

! ----------------------------------
! Store private internal state in GC
! ----------------------------------
    call ESMF_UserCompSetInternalState (GC, 'IGNI_State', wrap, STATUS)
    VERIFY_(STATUS)

! -----------------------------------------------------------
! Get the configuration
! -----------------------------------------------------------
    call MAPL_GetObjectFromGC(GC, MAPL, __RC__)

    call MAPL_GetResource(MAPL, NUM_ENSEMBLE, label='NUM_LDAS_ENSEMBLE:', DEFAULT=1, __RC__)


    ! at the moment, this G.C. will refresh when the land parent refreshes
    call MAPL_GetResource(MAPL, run_dt, label='RUN_DT:', __RC__)
    call MAPL_GetResource(MAPL, dt, label=trim(COMP_NAME)//'_DT:', default=run_dt, __RC__)

    ! set the 'local noon' method property
    call MAPL_GetResource(MAPL, resource_file, label='IGNI_RC:', default='GEOS_IgniGridComp.rc', __RC__)
    
    config = ESMF_ConfigCreate(__RC__)
    call ESMF_ConfigLoadFile(config, resource_file, __RC__)
    call MAPL_GetResource(config, local_noon, label='LOCAL_NOON:', default='SOLAR', __RC__)
    call ESMF_ConfigDestroy(config, __RC__)

    select case (local_noon)
        case ('SOLAR')
            self%local_noon_method = LOCAL_NOON_SOLAR
        case ('LST')
            self%local_noon_method = LOCAL_NOON_LST
        case DEFAULT
            self%local_noon_method = LOCAL_NOON_SOLAR
    end select

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

    call MAPL_AddImportSpec(GC,                    &
         LONG_NAME  = 'surface pressure',          &
         UNITS      = 'Pa',                        &
         SHORT_NAME = 'PS',                        &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone,          &
         RESTART    = MAPL_RestartSkip, __RC__)

#if (1)
    call MAPL_AddImportSpec(GC,                    &
         SHORT_NAME = 'PRLAND',                    & 
         LONG_NAME  = 'total precipitation land',  &
         UNITS      = 'kg m-2 s-1',                &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone,          &
         RESTART    = MAPL_RestartSkip, __RC__)
#else
    call MAPL_AddImportSpec(GC,                    &
         LONG_NAME  = 'liquid water convective precipitation', &
         UNITS      = 'kg m-2 s-1',                &
         SHORT_NAME = 'PCU',                       &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone,          &
         RESTART    = MAPL_RestartOptional, __RC__)

    call MAPL_AddImportSpec(GC,                    &
         LONG_NAME  = 'liquid water large scale precipitation', &
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
         SHORT_NAME = 'ASNOW',                     &
         LONG_NAME  = 'fractional area of land snowcover', &
         UNITS      = '1',                         &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone,          &
         RESTART    = MAPL_RestartSkip, __RC__)

    call MAPL_AddImportSpec(GC,                    &
         SHORT_NAME = 'SNOWDP',                    &
         LONG_NAME  = 'snow depth within snow covered area fraction', &
         UNITS      = 'm',                         &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone,          &
         RESTART    = MAPL_RestartSkip, __RC__)

    call MAPL_AddImportSpec(GC,                    &
         SHORT_NAME = 'SWDOWNLAND',                &
         LONG_NAME  = 'incident shortwave land',   &
         UNITS      = 'W m-2',                     &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone,          &
         RESTART    = MAPL_RestartSkip, __RC__)



! -----------------------------------------------------------
! Internal State 
! -----------------------------------------------------------

    ! hourly

    call MAPL_AddInternalSpec(GC,                  &
         SHORT_NAME = 'FFMC',                      &
         LONG_NAME  = 'fine fuel moisture code',   &
         UNITS      = '1',                         &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone,          &
         RESTART    = MAPL_RestartOptional,        &
         DEFAULT    = FFMC_INIT, __RC__)

    call MAPL_AddInternalSpec(GC,                  &
         SHORT_NAME = 'GFMC',                      &
         LONG_NAME  = 'grass fuel moisture code',  &
         UNITS      = '1',                         &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone,          &
         RESTART    = MAPL_RestartOptional,        &
         DEFAULT    = FFMC_INIT, __RC__)

    ! daily
    call MAPL_AddInternalSpec(GC,                  &
         SHORT_NAME = 'FFMC_DAILY',                &
         LONG_NAME  = 'fine fuel moisture code (daily)',   &
         UNITS      = '1',                         &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone,          &
         RESTART    = MAPL_RestartOptional,        &
         DEFAULT    = FFMC_INIT, __RC__)

    call MAPL_AddInternalSpec(GC,                  &
         SHORT_NAME = 'DMC_DAILY',                 &
         LONG_NAME  = 'duff moisture code (daily)',&
         UNITS      = '1',                         &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone,          &
         RESTART    = MAPL_RestartOptional,        &
         DEFAULT    = DMC_INIT, __RC__)

   call MAPL_AddInternalSpec(GC,                   &
         SHORT_NAME = 'DC_DAILY',                  &
         LONG_NAME  = 'drought code (daily)',      &
         UNITS      = '1',                         &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone,          &
         RESTART    = MAPL_RestartOptional,        &
         DEFAULT    = DC_INIT, __RC__) 

    call MAPL_AddInternalSpec(GC,                  &
         SHORT_NAME = 'DPR_LOCAL_NOON',            &
         LONG_NAME  = 'total precipitation since local noon', &
         UNITS      = 'kg m-2',                    &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone,          &
         ADD2EXPORT = .true.,                      &
!        PRECISION  = MAPL_R4,                     &
         RESTART    = MAPL_RestartOptional,        &
         DEFAULT    = 0.0, __RC__)

    call MAPL_AddInternalSpec(GC,                  &
         SHORT_NAME = 'PR_LOCAL_NOON',             &
         LONG_NAME  = '24-hr total precipitation at local noon', &
         UNITS      = 'kg m-2',                    &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone,          &
         ADD2EXPORT = .true.,                      &
         RESTART    = MAPL_RestartOptional,        &
         DEFAULT    = 0.0, __RC__)


    call MAPL_AddInternalSpec(GC,                  &
         SHORT_NAME = 'T_LOCAL_NOON',              &
         LONG_NAME  = 'temperature at local noon', &
         UNITS      = 'K',                         &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone,          &
         ADD2EXPORT = .true.,                      &
         RESTART    = MAPL_RestartOptional,        &
         DEFAULT    = 0.0, __RC__)

    call MAPL_AddInternalSpec(GC,                  &
         SHORT_NAME = 'RH_LOCAL_NOON',             &
         LONG_NAME  = 'relative humidity at local_noon', &
         UNITS      = '1',                         &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone,          &
         ADD2EXPORT = .true.,                      &
         RESTART    = MAPL_RestartOptional,        &
         DEFAULT    = 0.50, __RC__)

    call MAPL_AddInternalSpec(GC,                  &
         SHORT_NAME = 'WS_LOCAL_NOON',             &
         LONG_NAME  = 'wind speed at local noon',  &
         UNITS      = 'm s-1',                     &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone,          &
         ADD2EXPORT = .true.,                      &
         RESTART    = MAPL_RestartOptional,        &
         DEFAULT    = 3.0, __RC__)


! -----------------------------------------------------------
! Export Variables
! -----------------------------------------------------------

    ! hourly

    call MAPL_AddExportSpec(GC,                    &
         SHORT_NAME = 'FFMC',                      &
         LONG_NAME  = 'fine fuel moisture code',   &
         UNITS      = '1',                         &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone, __RC__)

    call MAPL_AddExportSpec(GC,                    &
         SHORT_NAME = 'GFMC',                      &
         LONG_NAME  = 'grass fuel moisture code',  &
         UNITS      = '1',                         &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone, __RC__)

    call MAPL_AddExportSpec(GC,                    &
         SHORT_NAME = 'DMC',                       &
         LONG_NAME  = 'duff moisture code',        &
         UNITS      = '1',                         &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone, __RC__)

    call MAPL_AddExportSpec(GC,                    &
         SHORT_NAME = 'DC',                        &
         LONG_NAME  = 'drought code',              &
         UNITS      = '1',                         &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone, __RC__)

    call MAPL_AddExportSpec(GC,                    &
         SHORT_NAME = 'FWI',                       &
         LONG_NAME  = 'fire weather index',        &
         UNITS      = '1',                         &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone, __RC__) 

    call MAPL_AddExportSpec(GC,                    &
         SHORT_NAME = 'BUI',                       &
         LONG_NAME  = 'buildup index',             &
         UNITS      = '1',                         &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone, __RC__) 

    call MAPL_AddExportSpec(GC,                    &
         SHORT_NAME = 'ISI',                       &
         LONG_NAME  = 'initial spread index',      &
         UNITS      = '1',                         &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone, __RC__) 

    call MAPL_AddExportSpec(GC,                    &
         SHORT_NAME = 'DSR',                       &
         LONG_NAME  = 'daily severity rating',     &
         UNITS      = '1',                         &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone, __RC__)


    ! daily

    call MAPL_AddExportSpec(GC,                    &
         SHORT_NAME = 'FFMC_DAILY',                &
         LONG_NAME  = 'fine fuel moisture code (daily)', &
         UNITS      = '1',                         &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone, __RC__)

    call MAPL_AddExportSpec(GC,                    &
         SHORT_NAME = 'DMC_DAILY',                 &
         LONG_NAME  = 'duff moisture code (daily)',&
         UNITS      = '1',                         &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone, __RC__)

    call MAPL_AddExportSpec(GC,                    &
         SHORT_NAME = 'DC_DAILY',                  &
         LONG_NAME  = 'drought code (daily)',      &
         UNITS      = '1',                         &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone, __RC__)

    call MAPL_AddExportSpec(GC,                    &
         SHORT_NAME = 'FWI_DAILY',                 &
         LONG_NAME  = 'fire weather index (daily)',&
         UNITS      = '1',                         &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone, __RC__) 

    call MAPL_AddExportSpec(GC,                    &
         SHORT_NAME = 'BUI_DAILY',                 &
         LONG_NAME  = 'buildup index (daily)',     &
         UNITS      = '1',                         &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone, __RC__) 

    call MAPL_AddExportSpec(GC,                    &
         SHORT_NAME = 'ISI_DAILY',                 &
         LONG_NAME  = 'initial spread index (daily)', &
         UNITS      = '1',                         &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone, __RC__) 

    call MAPL_AddExportSpec(GC,                    &
         SHORT_NAME = 'DSR_DAILY',                 &
         LONG_NAME  = 'daily severity rating (daily)', &
         UNITS      = '1',                         &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone, __RC__)


    ! local noon patches

    call MAPL_AddExportSpec(GC,                    &
         SHORT_NAME = 'FFMC_DAILY_',               &
         LONG_NAME  = 'fine fuel moisture code (daily)', &
         UNITS      = '1',                         &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone, __RC__)

    call MAPL_AddExportSpec(GC,                    &
         SHORT_NAME = 'DMC_DAILY_',                &
         LONG_NAME  = 'duff moisture code (daily)', &
         UNITS      = '1',                         &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone, __RC__)

    call MAPL_AddExportSpec(GC,                    &
         SHORT_NAME = 'DC_DAILY_',                 &
         LONG_NAME  = 'drought code (daily)',      &
         UNITS      = '1',                         &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone, __RC__)

    call MAPL_AddExportSpec(GC,                    &
         SHORT_NAME = 'FWI_DAILY_',                &
         LONG_NAME  = 'fire weather index (daily)',&
         UNITS      = '1',                         &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone, __RC__) 

    call MAPL_AddExportSpec(GC,                    &
         SHORT_NAME = 'BUI_DAILY_',                &
         LONG_NAME  = 'buildup index (daily)',     &
         UNITS      = '1',                         &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone, __RC__) 

    call MAPL_AddExportSpec(GC,                    &
         SHORT_NAME = 'ISI_DAILY_',                &
         LONG_NAME  = 'initial spread index(daily)', &
         UNITS      = '1',                         &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone, __RC__) 

    call MAPL_AddExportSpec(GC,                    &
         SHORT_NAME = 'DSR_DAILY_',                &
         LONG_NAME  = 'daily severity rating (daily)', &
         UNITS      = '1',                         &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone, __RC__)


    ! flammability and ignition sources

    call MAPL_AddExportSpec(GC,                    & 
         SHORT_NAME = 'VPD',                       &
         LONG_NAME  = 'vapor pressure deficit',    &
         UNITS      = 'Pa',                        &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone, __RC__)


!EOS


! Set the Profiling timers
! ------------------------
    call MAPL_TimerAdd(GC, name='TOTAL'     , __RC__)
!   call MAPL_TimerAdd(GC, name='INITIALIZE', __RC__)
    call MAPL_TimerAdd(GC, name='RUN'       , __RC__)
    call MAPL_TimerAdd(GC, name='-CFFWI'    , __RC__)
    call MAPL_TimerAdd(GC, name='--daily'   , __RC__)
    call MAPL_TimerAdd(GC, name='--hourly'  , __RC__)
!   call MAPL_TimerAdd(GC, name='FINALIZE'  , __RC__)


!------------------------------------------------------------
! Set generic init and final methods
!------------------------------------------------------------

    call MAPL_GenericSetServices(GC, __RC__)

    RETURN_(ESMF_SUCCESS)
  
  end subroutine SetServices



! -----------------------------------------------------------
! RUN1 -- Run (phase = 1) method of the IGNI component
! -----------------------------------------------------------

  subroutine RUN1 (GC, IMPORT, EXPORT, CLOCK, RC)

! -----------------------------------------------------------
! !ARGUMENTS:

    type(ESMF_GridComp), intent(inout) :: GC    
    type(ESMF_State),    intent(inout) :: IMPORT
    type(ESMF_State),    intent(inout) :: EXPORT
    type(ESMF_Clock),    intent(inout) :: CLOCK
    integer, optional,   intent(  out) :: RC

!EOP

! ErrLog Variables

    character(len=ESMF_MAXSTR) :: Iam
    integer                    :: STATUS
    character(len=ESMF_MAXSTR) :: COMP_NAME

! Locals

    type (MAPL_MetaComp), pointer :: MAPL => null()


! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    call ESMF_GridCompGet(GC, name=COMP_NAME, __RC__)
  
    Iam = trim(COMP_NAME) // 'Run'

! Get my internal MAPL_Generic state
! -----------------------------------------------------------

    call MAPL_GetObjectFromGC(GC, MAPL, __RC__)

    call MAPL_TimerOn(MAPL, 'TOTAL')
    call MAPL_TimerOn(MAPL, 'RUN'  )


! Get file names from configuration
! -----------------------------------------------------------

    if (NUM_ENSEMBLE > 1) then
        !comp_name should be IGNIxxxx...
    end if


!  All done
! ---------
    call MAPL_TimerOff(MAPL, 'RUN'  )
    call MAPL_TimerOff(MAPL, 'TOTAL')

    RETURN_(ESMF_SUCCESS)

  end subroutine RUN1



! -----------------------------------------------------------
! RUN2 -- Run (phase = 2) method of the IGNI component
! -----------------------------------------------------------

  subroutine RUN2 (GC, IMPORT, EXPORT, CLOCK, RC)

! -----------------------------------------------------------
! !ARGUMENTS:

    type(ESMF_GridComp), intent(inout) :: GC    
    type(ESMF_State),    intent(inout) :: IMPORT
    type(ESMF_State),    intent(inout) :: EXPORT
    type(ESMF_Clock),    intent(inout) :: CLOCK
    integer, optional,   intent(  out) :: RC

!EOP

! ErrLog Variables

    character(len=ESMF_MAXSTR) :: Iam
    integer                    :: STATUS
    character(len=ESMF_MAXSTR) :: COMP_NAME

! Locals

    type (MAPL_MetaComp), pointer :: MAPL => null()


! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    call ESMF_GridCompGet(GC, name=COMP_NAME, __RC__)
  
    Iam = trim(COMP_NAME) // 'Run2'

! Get my internal MAPL_Generic state
! -----------------------------------------------------------

    call MAPL_GetObjectFromGC(GC, MAPL, __RC__)

    call MAPL_TimerOn(MAPL, 'TOTAL')
    call MAPL_TimerOn(MAPL, 'RUN'  )


! Get file names from configuration
! -----------------------------------------------------------

    if (NUM_ENSEMBLE > 1) then
        !comp_name should be IGNIxxxx...
    end if


! Fire weather indexes
! --------------------
    call MAPL_TimerOn(MAPL, '-CFFWI')

    call MAPL_TimerOn(MAPL,  '--daily')
    call CFFWI_DAILY (GC, IMPORT, EXPORT, CLOCK, __RC__)
    call MAPL_TimerOff(MAPL, '--daily')

    call MAPL_TimerOn(MAPL,  '--hourly')
    call CFFWI_HOURLY(GC, IMPORT, EXPORT, CLOCK, __RC__)
    call MAPL_TimerOff(MAPL, '--hourly')

    call MAPL_TimerOff(MAPL, '-CFFWI')


! flammability and ignition sources
! ---------------------------------
    call REBURN(GC, IMPORT, EXPORT, CLOCK, __RC__)



!  All done
! ---------
    call MAPL_TimerOff(MAPL, 'RUN'  )
    call MAPL_TimerOff(MAPL, 'TOTAL')

    RETURN_(ESMF_SUCCESS)

  end subroutine RUN2


! -----------------------------------------------------------
! CFFWI_DAILY -- Runs the daily CFFWI
! -----------------------------------------------------------

  subroutine CFFWI_DAILY (GC, IMPORT, EXPORT, CLOCK, RC)

! -----------------------------------------------------------
! !ARGUMENTS:

    type(ESMF_GridComp), intent(inout) :: GC
    type(ESMF_State),    intent(inout) :: IMPORT
    type(ESMF_State),    intent(inout) :: EXPORT
    type(ESMF_Clock),    intent(inout) :: CLOCK
    integer, optional,   intent(  out) :: RC

!EOP

! ErrLog Variables

    character(len=ESMF_MAXSTR)   :: Iam
    character(len=ESMF_MAXSTR)   :: COMP_NAME
    integer                      :: STATUS


! Private internal state
    type (IGNI_State), pointer   :: self
    type (wrap_)                 :: wrap


! IMPORT pointers

    real, dimension(:), pointer :: T2M    => null()
    real, dimension(:), pointer :: Q2M    => null()
    real, dimension(:), pointer :: U10M   => null()
    real, dimension(:), pointer :: V10M   => null()
    real, dimension(:), pointer :: PS     => null()
    real, dimension(:), pointer :: ASNOW  => null()
    real, dimension(:), pointer :: SNOWDP => null()
    real, dimension(:), pointer :: SWDOWN => null()
#if (1)
    real, dimension(:), pointer :: PRLAND => null()
#else
    real, dimension(:), pointer :: PCU    => null()
    real, dimension(:), pointer :: PLS    => null()
    real, dimension(:), pointer :: SNO    => null()
#endif


! INTERNAL pointers

    real, dimension(:), pointer :: FFMC0_daily => null()
    real, dimension(:), pointer :: DMC0_daily  => null()
    real, dimension(:), pointer :: DC0_daily   => null()

    real, dimension(:), pointer :: DPR_noon    => null()
    real, dimension(:), pointer :: PR_noon     => null()
    real, dimension(:), pointer :: RH_noon     => null()
    real, dimension(:), pointer :: T_noon      => null()
    real, dimension(:), pointer :: WS_noon     => null()
     

! EXPORT pointers

    ! daily
    real, dimension(:), pointer :: FFMC_daily  => null()
    real, dimension(:), pointer :: DMC_daily   => null()
    real, dimension(:), pointer :: DC_daily    => null()
    real, dimension(:), pointer :: ISI_daily   => null()
    real, dimension(:), pointer :: FWI_daily   => null()
    real, dimension(:), pointer :: BUI_daily   => null()
    real, dimension(:), pointer :: DSR_daily   => null()

    ! local noon
    real, dimension(:), pointer :: FFMC_daily_ => null()
    real, dimension(:), pointer :: DMC_daily_  => null()
    real, dimension(:), pointer :: DC_daily_   => null()
    real, dimension(:), pointer :: ISI_daily_  => null()
    real, dimension(:), pointer :: FWI_daily_  => null()
    real, dimension(:), pointer :: BUI_daily_  => null()
    real, dimension(:), pointer :: DSR_daily_  => null()


! Misc

    type(MAPL_MetaComp), pointer :: MAPL => null()
    type(ESMF_State)             :: INTERNAL
    type(MAPL_SunOrbit)          :: ORBIT


    type(ESMF_Time)         :: time
    type(ESMF_Alarm)        :: run_alarm
    type(ESMF_TimeInterval) :: ring_interval
    real(ESMF_KIND_R8)      :: time_step

    integer :: NT

    integer :: year, month, day, hr, mn, sc
    real    :: dt

    real, pointer, dimension(:) :: LATS => null()
    real, pointer, dimension(:) :: LONS => null()

    real, allocatable, dimension(:) :: tmpISI
    real, allocatable, dimension(:) :: tmpBUI
    real, allocatable, dimension(:) :: tmpDSR
    real, allocatable, dimension(:) :: tmpFWI

    real, allocatable, dimension(:) :: dt_local_noon

    real, allocatable, dimension(:) :: LSHA0
    real, allocatable, dimension(:) :: LSHA1

    logical, allocatable, dimension(:) :: isNoon


! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    call ESMF_GridCompGet(GC, name=COMP_NAME, __RC__)

    Iam = trim(COMP_NAME) // 'CFFWI_DAILY'


! Get my internal MAPL_Generic state
! -----------------------------------------------------------

    call MAPL_GetObjectFromGC(GC, MAPL, __RC__)


    call MAPL_Get(MAPL, TILELATS=LATS, &
                        TILELONS=LONS, &
                        INTERNAL_ESMF_STATE=INTERNAL, __RC__)

    NT = SIZE(LONS)

    NO_LAND_AREAS: if (NT == 0) then
        RETURN_(ESMF_SUCCESS)
    end if NO_LAND_AREAS


! Get my private internal state
! -----------------------------
    call ESMF_UserCompGetInternalState (GC, 'IGNI_State', wrap, STATUS)
    VERIFY_(STATUS)
    self => wrap%ptr


! Get pointers to internal variables
! ----------------------------------

    call MAPL_GetPointer(INTERNAL, FFMC0_daily, 'FFMC_DAILY',     __RC__)
    call MAPL_GetPointer(INTERNAL, DMC0_daily,  'DMC_DAILY',      __RC__)
    call MAPL_GetPointer(INTERNAL, DC0_daily,   'DC_DAILY',       __RC__)

    call MAPL_GetPointer(INTERNAL, DPR_noon,    'DPR_LOCAL_NOON', __RC__)
    call MAPL_GetPointer(INTERNAL, PR_noon,     'PR_LOCAL_NOON',  __RC__)
    call MAPL_GetPointer(INTERNAL, T_noon,      'T_LOCAL_NOON',   __RC__)
    call MAPL_GetPointer(INTERNAL, RH_noon,     'RH_LOCAL_NOON',  __RC__)
    call MAPL_GetPointer(INTERNAL, WS_noon,     'WS_LOCAL_NOON',  __RC__)


! Get pointers to imports
! -----------------------

    call MAPL_GetPointer(IMPORT, PS,     'PS',     __RC__)
    call MAPL_GetPointer(IMPORT, Q2M,    'MOQ2M',  __RC__)
    call MAPL_GetPointer(IMPORT, T2M,    'MOT2M',  __RC__)
    call MAPL_GetPointer(IMPORT, U10M,   'MOU10M', __RC__)
    call MAPL_GetPointer(IMPORT, V10M,   'MOV10M', __RC__)
    call MAPL_GetPointer(IMPORT, ASNOW,  'ASNOW',  __RC__)
    call MAPL_GetPointer(IMPORT, SNOWDP, 'SNOWDP', __RC__)
#if (1)
    call MAPL_GetPointer(IMPORT, PRLAND, 'PRLAND', __RC__)
#else
    call MAPL_GetPointer(IMPORT, PCU,    'PCU',    __RC__)
    call MAPL_GetPointer(IMPORT, PLS,    'PLS',    __RC__)
    call MAPL_GetPointer(IMPORT, SNO,    'SNO',    __RC__)
#endif


! Get pointers to exports
! -----------------------

    ! global
    call MAPL_GetPointer(EXPORT, FFMC_daily,  'FFMC_DAILY',  __RC__)
    call MAPL_GetPointer(EXPORT, DMC_daily,   'DMC_DAILY',   __RC__)
    call MAPL_GetPointer(EXPORT, DC_daily,    'DC_DAILY',    __RC__)
    call MAPL_GetPointer(EXPORT, FWI_daily,   'FWI_DAILY',   __RC__)
    call MAPL_GetPointer(EXPORT, ISI_daily,   'ISI_DAILY',   __RC__)
    call MAPL_GetPointer(EXPORT, BUI_daily,   'BUI_DAILY',   __RC__)
    call MAPL_GetPointer(EXPORT, DSR_daily,   'DSR_DAILY',   __RC__)
 
    ! local noon patch
    call MAPL_GetPointer(EXPORT, FFMC_daily_, 'FFMC_DAILY_', __RC__)
    call MAPL_GetPointer(EXPORT, DMC_daily_,  'DMC_DAILY_',  __RC__)
    call MAPL_GetPointer(EXPORT, DC_daily_,   'DC_DAILY_',   __RC__)
    call MAPL_GetPointer(EXPORT, FWI_daily_,  'FWI_DAILY_',  __RC__)
    call MAPL_GetPointer(EXPORT, ISI_daily_,  'ISI_DAILY_',  __RC__)
    call MAPL_GetPointer(EXPORT, BUI_daily_,  'BUI_DAILY_',  __RC__)
    call MAPL_GetPointer(EXPORT, DSR_daily_,  'DSR_DAILY_',  __RC__)


! Get the time step
! -----------------

    call MAPL_Get(MAPL, RunAlarm=run_alarm, __RC__)
    call ESMF_AlarmGet(run_alarm, ringInterval=ring_interval, __RC__)

    call ESMF_TimeIntervalGet(ring_interval, s_r8=time_step, __RC__)
    dt = real(time_step)


! Construct local noon mask
! -------------------------

    call ESMF_ClockGet(CLOCK, currTime=time, __RC__)
    call ESMF_TimeGet(time, yy=year, mm=month, dd=day, h=hr, m=mn, s=sc, __RC__)

    allocate(isNoon(NT), __STAT__)

    if (self%local_noon_method == LOCAL_NOON_SOLAR) then
#if (1)
        _ASSERT(.false., 'Precise local solar noon is disabled, please select LST instead.')
#else
        allocate(LSHA0(NT), LSHA1(NT), __STAT__)

        call MAPL_Get(MAPL, ORBIT=ORBIT, __RC__)
        call MAPL_SunGetLocalSolarHourAngle(ORBIT, LONS, LSHA0, TIME=time, __RC__)
        call MAPL_SunGetLocalSolarHourAngle(ORBIT, LONS, LSHA1, TIME=time+ring_interval, __RC__)
    
        isNoon = (LSHA0 <= 0) .and. (LSHA1 > 0)
    
        deallocate(LSHA0, LSHA1, __STAT__)
#endif
    else    
        allocate(dt_local_noon(NT), __STAT__)
        dt_local_noon = ((hr-12)*3600 + mn*60 + sc) + ((24*3600)/(2*MAPL_PI))*LONS

        isNoon = (dt_local_noon >= 0) .and. (dt_local_noon < dt)

        deallocate(dt_local_noon, __STAT__)
    end if


! Accumulate precip
! -----------------
#if (1)
    DPR_NOON = DPR_NOON + PRLAND*dt
#else
    DPR_NOON = DPR_NOON + (PCU + PLS + SNO)*dt
#endif


! Update local noon patches
! -------------------------

    where (isNoon)
        T_noon   = T2M
        WS_noon  = sqrt(U10M*U10M + V10M*V10M)
        RH_noon  = min(Q2M / GEOS_QSAT(T2M, PS, PASCALS=.true.), 1.0)
        PR_noon  = DPR_noon
        DPR_noon = 0.0
    end where


! Run the daily system
! --------------------
    allocate(tmpISI(NT), tmpBUI(NT), tmpDSR(NT), tmpFWI(NT), __STAT__)

    tmpISI = MAPL_UNDEF
    tmpBUI = MAPL_UNDEF
    tmpDSR = MAPL_UNDEF
    tmpFWI = MAPL_UNDEF

    call cffwi_daily_driver(FFMC0_daily, DMC0_daily, DC0_daily, &
                            tmpISI,  tmpBUI, tmpFWI, tmpDSR,    &
                            T_noon, RH_noon, WS_noon, PR_noon,  &
                            ASNOW, SNOWDP, LATS, isNoon, month, NT)

! Update exports
! --------------

    if (associated(FFMC_daily)) FFMC_daily = FFMC0_daily
    if (associated(DMC_daily))   DMC_daily = DMC0_daily
    if (associated(DC_daily))     DC_daily = DC0_daily
    if (associated(ISI_daily))   ISI_daily = tmpISI
    if (associated(BUI_daily))   BUI_daily = tmpBUI
    if (associated(FWI_daily))   FWI_daily = tmpFWI
    if (associated(DSR_daily))   DSR_daily = tmpDSR


    if (associated(FFMC_daily_)) then
        where (isNoon)
            FFMC_daily_ = FFMC0_daily
        elsewhere
            FFMC_daily_ = MAPL_UNDEF
        end where
    end if

    if (associated(DMC_daily_)) then
        where (isNoon)
            DMC_daily_ = DMC0_daily
        elsewhere
            DMC_daily_ = MAPL_UNDEF
        end where
    end if

    if (associated(DC_daily_)) then
        where (isNoon) 
            DC_daily_ = DC0_daily
        elsewhere
            DC_daily_ = MAPL_UNDEF
        end where
    end if 

    if (associated(ISI_daily_)) then
        where (isNoon)
            ISI_daily_ = tmpISI
        elsewhere
            ISI_daily_ = MAPL_UNDEF
        end where 
    end if

    if (associated(BUI_daily_)) then
        where (isNoon)
            BUI_daily_ = tmpBUI
        elsewhere
            BUI_daily_ = MAPL_UNDEF
        end where
    end if

    if (associated(FWI_daily_)) then
        where (isNoon) 
            FWI_daily_ = tmpFWI
        elsewhere
            FWI_daily_ = MAPL_UNDEF
        end where
    end if

    if (associated(DSR_daily_)) then
        where (isNoon)
            DSR_daily_ = tmpDSR
        elsewhere
            DSR_daily_ = MAPL_UNDEF
        end where
    end if


    deallocate(tmpISI, tmpBUI, tmpDSR, tmpFWI, __STAT__)


!  All done
! ---------

    RETURN_(ESMF_SUCCESS)

  end subroutine CFFWI_DAILY


! -----------------------------------------------------------
! CFFWI_HOURLY -- Runs the hourly CFFWI
! -----------------------------------------------------------

  subroutine CFFWI_HOURLY (GC, IMPORT, EXPORT, CLOCK, RC)

! -----------------------------------------------------------
! !ARGUMENTS:

    type(ESMF_GridComp), intent(inout) :: GC
    type(ESMF_State),    intent(inout) :: IMPORT
    type(ESMF_State),    intent(inout) :: EXPORT
    type(ESMF_Clock),    intent(inout) :: CLOCK
    integer, optional,   intent(  out) :: RC

!EOP

! ErrLog Variables

    character(len=ESMF_MAXSTR)   :: Iam
    character(len=ESMF_MAXSTR)   :: COMP_NAME
    integer                      :: STATUS


! IMPORT pointers

    real, dimension(:), pointer :: T2M    => null()
    real, dimension(:), pointer :: Q2M    => null()
    real, dimension(:), pointer :: U10M   => null()
    real, dimension(:), pointer :: V10M   => null()
    real, dimension(:), pointer :: PS     => null()
    real, dimension(:), pointer :: ASNOW  => null()
    real, dimension(:), pointer :: SNOWDP => null()
    real, dimension(:), pointer :: SWDOWN => null()
#if (1)
    real, dimension(:), pointer :: PRLAND => null()
#else
    real, dimension(:), pointer :: PCU    => null()
    real, dimension(:), pointer :: PLS    => null()
    real, dimension(:), pointer :: SNO    => null()
#endif


! INTERNAL pointers

    real, dimension(:), pointer :: FFMC0  => null()
    real, dimension(:), pointer :: GFMC0  => null()
    real, dimension(:), pointer :: DMC0   => null()
    real, dimension(:), pointer :: DC0    => null()


! EXPORT pointers

    ! hourly
    real, dimension(:), pointer :: FFMC => null()
    real, dimension(:), pointer :: GFMC => null()
    real, dimension(:), pointer :: DMC  => null()
    real, dimension(:), pointer :: DC   => null()
    real, dimension(:), pointer :: ISI  => null()
    real, dimension(:), pointer :: FWI  => null()
    real, dimension(:), pointer :: BUI  => null()
    real, dimension(:), pointer :: DSR  => null()


! Misc

    type(MAPL_MetaComp), pointer :: MAPL => null()
    type(ESMF_State)             :: INTERNAL


    type(ESMF_Time)         :: time
    type(ESMF_Alarm)        :: run_alarm
    type(ESMF_TimeInterval) :: ring_interval
    real(ESMF_KIND_R8)      :: time_step

    integer :: NT

    integer :: year, month, day, hr, mn, sc
    real    :: dt

    real, pointer, dimension(:) :: LATS => null()
    real, pointer, dimension(:) :: LONS => null()

    real, allocatable, dimension(:) :: tmpISI
    real, allocatable, dimension(:) :: tmpBUI
    real, allocatable, dimension(:) :: tmpDSR
    real, allocatable, dimension(:) :: tmpFWI


! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    call ESMF_GridCompGet(GC, name=COMP_NAME, __RC__)

    Iam = trim(COMP_NAME) // 'CFFWI_HOURLY'


! Get my internal MAPL_Generic state
! -----------------------------------------------------------

    call MAPL_GetObjectFromGC(GC, MAPL, __RC__)


    call MAPL_Get(MAPL, TILELATS=LATS, &
                        TILELONS=LONS, &
                        INTERNAL_ESMF_STATE=INTERNAL, __RC__)

    NT = SIZE(LONS)

    NO_LAND_AREAS: if (NT == 0) then
        RETURN_(ESMF_SUCCESS)
    end if NO_LAND_AREAS


! Get pointers to internal variables
! ----------------------------------

    call MAPL_GetPointer(INTERNAL, FFMC0, 'FFMC',      __RC__)
    call MAPL_GetPointer(INTERNAL, GFMC0, 'GFMC',      __RC__)
    call MAPL_GetPointer(INTERNAL, DMC0,  'DMC_DAILY', __RC__)  ! requires daily DMC
    call MAPL_GetPointer(INTERNAL, DC0,   'DC_DAILY',  __RC__)  ! requires daily DC 


! Get pointers to imports
! -----------------------

    call MAPL_GetPointer(IMPORT, PS,     'PS',         __RC__)
    call MAPL_GetPointer(IMPORT, Q2M,    'MOQ2M',      __RC__)
    call MAPL_GetPointer(IMPORT, T2M,    'MOT2M',      __RC__)
    call MAPL_GetPointer(IMPORT, U10M,   'MOU10M',     __RC__)
    call MAPL_GetPointer(IMPORT, V10M,   'MOV10M',     __RC__)
    call MAPL_GetPointer(IMPORT, ASNOW,  'ASNOW',      __RC__)
    call MAPL_GetPointer(IMPORT, SNOWDP, 'SNOWDP',     __RC__)
    call MAPL_GetPointer(IMPORT, SWDOWN, 'SWDOWNLAND', __RC__)
#if (1)
    call MAPL_GetPointer(IMPORT, PRLAND, 'PRLAND',     __RC__)
#else
    call MAPL_GetPointer(IMPORT, PCU,    'PCU',        __RC__)
    call MAPL_GetPointer(IMPORT, PLS,    'PLS',        __RC__)
    call MAPL_GetPointer(IMPORT, SNO,    'SNO',        __RC__)
#endif


! Get pointers to exports
! -----------------------

    ! global
    call MAPL_GetPointer(EXPORT, FFMC,  'FFMC',  __RC__)
    call MAPL_GetPointer(EXPORT, GFMC,  'GFMC',  __RC__)
    call MAPL_GetPointer(EXPORT, DMC,   'DMC',   __RC__) ! mostly for symmetry, same as DMC_DAILY
    call MAPL_GetPointer(EXPORT, DC,    'DC',    __RC__) ! mostly for symmetry, same as  DC_DAILY
    call MAPL_GetPointer(EXPORT, FWI,   'FWI',   __RC__)
    call MAPL_GetPointer(EXPORT, ISI,   'ISI',   __RC__)
    call MAPL_GetPointer(EXPORT, BUI,   'BUI',   __RC__)
    call MAPL_GetPointer(EXPORT, DSR,   'DSR',   __RC__)


! Get the time step
! -----------------

    call MAPL_Get(MAPL, RunAlarm=run_alarm, __RC__)
    call ESMF_AlarmGet(run_alarm, ringInterval=ring_interval, __RC__)

    call ESMF_TimeIntervalGet(ring_interval, s_r8=time_step, __RC__)
    dt = real(time_step)


! Get date and time
! -----------------

    call ESMF_ClockGet(CLOCK, currTime=time, __RC__)

    call ESMF_TimeGet(time, yy=year, mm=month, dd=day, h=hr, m=mn, s=sc, __RC__)


! Update local noon patches
! -------------------------

    allocate(tmpISI(NT), tmpBUI(NT), tmpDSR(NT), tmpFWI(NT), __STAT__)

    tmpISI  = MAPL_UNDEF
    tmpBUI  = MAPL_UNDEF
    tmpDSR  = MAPL_UNDEF
    tmpFWI  = MAPL_UNDEF

    call cffwi_hourly_driver(FFMC0, GFMC0, DMC0, DC0, &
                             tmpISI, tmpBUI, tmpFWI, tmpDSR, &
                             T2M, &
                             min(Q2M / GEOS_QSAT(T2M, PS, PASCALS=.true.), 1.0), &
                             sqrt(U10M*U10M + V10M*V10M), &
                             PRLAND*dt, &
                             ASNOW, SNOWDP, &
                             SWDOWN, &
                             month, dt/3600.0, NT)

    ! update exports
    if (associated(FFMC)) FFMC = FFMC0
    if (associated(GFMC)) GFMC = GFMC0
    if (associated(DMC))   DMC = DMC0
    if (associated(DC))     DC = DC0
    if (associated(ISI))   ISI = tmpISI
    if (associated(BUI))   BUI = tmpBUI
    if (associated(FWI))   FWI = tmpFWI
    if (associated(DSR))   DSR = tmpDSR


    deallocate(tmpISI, tmpBUI, tmpDSR, tmpFWI, __STAT__)


!  All done
! ---------

    RETURN_(ESMF_SUCCESS)

  end subroutine CFFWI_HOURLY


! -----------------------------------------------------------
! REBURN -- flammability and ignition sources
! -----------------------------------------------------------

  subroutine REBURN (GC, IMPORT, EXPORT, CLOCK, RC)

! -----------------------------------------------------------
! !ARGUMENTS:

    type(ESMF_GridComp), intent(inout) :: GC
    type(ESMF_State),    intent(inout) :: IMPORT
    type(ESMF_State),    intent(inout) :: EXPORT
    type(ESMF_Clock),    intent(inout) :: CLOCK
    integer, optional,   intent(  out) :: RC

!EOP

! ErrLog Variables

    character(len=ESMF_MAXSTR)   :: Iam
    character(len=ESMF_MAXSTR)   :: COMP_NAME
    integer                      :: STATUS


! IMPORT pointers

    real, dimension(:), pointer :: T2M    => null()
    real, dimension(:), pointer :: Q2M    => null()
    real, dimension(:), pointer :: PS     => null()


! INTERNAL pointers
!   None


! EXPORT pointers

    real, dimension(:), pointer :: VPD => null()
 

! misc

    type(MAPL_MetaComp), pointer :: MAPL => null()
    type(ESMF_State)             :: INTERNAL

    integer :: NT

    real, pointer, dimension(:) :: LATS => null()
    real, pointer, dimension(:) :: LONS => null()


! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    call ESMF_GridCompGet(GC, name=COMP_NAME, __RC__)

    Iam = trim(COMP_NAME) // 'REBURN'


! Get my internal MAPL_Generic state
! -----------------------------------------------------------

    call MAPL_GetObjectFromGC(GC, MAPL, __RC__)


    call MAPL_Get(MAPL, TILELATS=LATS, &
                        TILELONS=LONS, &
                        INTERNAL_ESMF_STATE=INTERNAL, __RC__)

    NT = SIZE(LONS)

    NO_LAND_AREAS: if (NT == 0) then
        RETURN_(ESMF_SUCCESS)
    end if NO_LAND_AREAS



! Get pointers to imports
! -----------------------

    call MAPL_GetPointer(IMPORT, PS,     'PS',         __RC__)
    call MAPL_GetPointer(IMPORT, Q2M,    'MOQ2M',      __RC__)
    call MAPL_GetPointer(IMPORT, T2M,    'MOT2M',      __RC__)


! Get pointers to exports
! -----------------------

    call MAPL_GetPointer(EXPORT, VPD, 'VPD', __RC__)


! Update diagnostics
! -------------------------
    UPDATE_VPD: if (associated(VPD)) then
        ! VPD = e_s - e = e_s * (1 - RH)
        !
        ! e_s = P * Qsat/(MAPL_EPSILON + (1 - MAPL_EPSILON)*Qsat)
        ! MAPL_EQsat(T) is equivalent to the e_s expression

        VPD = MAPL_EQsat(T2M) * (1 - min(Q2M / GEOS_QSAT(T2M, PS, PASCALS=.true.), 1.0))
    end if UPDATE_VPD


!  All done
! ---------

    RETURN_(ESMF_SUCCESS)

  end subroutine REBURN



  subroutine cffwi_daily_driver(ffmc, dmc, dc, isi, bui, fwi, dsr, &
                                T, RH, wind, Pr, &
                                f_snow, snow_depth, &
                                latitude, is_noon, month, N)
 
    !
    ! Calculates daily FFMC, DMC, DC, ISI, BUI, FWI and DSR indexes.
    ! Note that: 
    !     - FFMC, DMC, and DC are updated only in areas where the local 
    !       noon mask is true
    !     - ISI, BUI, FWI and DSR are calculated in all areas using 
    !       the updated FFMC, DMC and DC.
    !

    implicit none

    real,    dimension(N), intent(in) :: T, RH, wind, Pr
    real,    dimension(N), intent(in) :: f_snow, snow_depth
    real,    dimension(N), intent(in) :: latitude
    logical, dimension(N), intent(in) :: is_noon
    integer,               intent(in) :: month
    integer,               intent(in) :: N

    real, dimension(N), intent(inout) :: ffmc, dmc, dc, isi, bui, fwi, dsr
    
    ! local
    integer :: i
    real    :: T_, RH_, Pr_
    real    :: lat_

    real, parameter :: DAILY = 24.0 ! hours

    do i = 1, N
        if (is_noon(i)) then
            T_   = T(i) - 273.15    ! temperature, C
            RH_  = 100 * RH(i)      ! relative humidity, %
            Pr_  = Pr(i)            ! precipitation since local noon, mm
            lat_ = MAPL_RADIANS_TO_DEGREES * latitude(i) 

            ! update FFMC, DMC and DC
            ffmc(i) = fine_fuel_moisture_code(ffmc(i), T_, RH_, wind(i), Pr_, DAILY)
            dmc(i)  = duff_moisture_code(dmc(i), T_, RH_, Pr_, month)
            dc(i)   = drought_code(dc(i), T_, Pr_, lat_, month)
        end if

        ! calculate ISI, BUI, FWI and DSR
        isi(i) = initial_spread_index(ffmc(i), wind(i))

        if (snow_depth(i) > 1e-3) then
            isi(i) = (1 - f_snow(i)) * isi(i)
        end if

        bui(i) = buildup_index(dmc(i), dc(i))

        fwi(i) = fire_weather_index(isi(i), bui(i))
        dsr(i) = daily_severity_rating(fwi(i))
    end do

  end subroutine cffwi_daily_driver



  subroutine cffwi_hourly_driver(ffmc, gfmc, dmc, dc, isi, bui, fwi, dsr, &
                                 T, RH, wind, Pr, &
                                 f_snow, snow_depth, &
                                 swdown, &
                                 month, time_step, N)
 
    !
    ! Calculates daily FFMC, ISI, BUI, FWI and DSR indexes.
    ! Note that DMC and DC are from the daily CFFWI and 
    ! are not modified here. 
    !

    implicit none

    real,    dimension(N), intent(in) :: T, RH, wind, Pr
    real,    dimension(N), intent(in) :: f_snow, snow_depth
    real,    dimension(N), intent(in) :: swdown
    integer,               intent(in) :: month
    real,                  intent(in) :: time_step
    integer,               intent(in) :: N

    real,    dimension(N), intent(in   ) :: dmc, dc
    real,    dimension(N), intent(inout) :: ffmc, gfmc, isi, bui, fwi, dsr
    
    ! local
    integer :: i
    real    :: T_, RH_, Pr_


    do i = 1, N
        T_   = T(i) - 273.15    ! temperature, C
        RH_  = 100 * RH(i)      ! relative humidity, %
        Pr_  = Pr(i)            ! precipitation, mm

        ! update FFMC
        ffmc(i) = fine_fuel_moisture_code(ffmc(i), T_, RH_, wind(i), Pr_, time_step)

        ! update GFMC
        gfmc(i) = grass_fuel_moisture_code(gfmc(i), T_, RH_, wind(i), Pr_, &
                                           swdown(i), NOMINAL_FINE_FUEL_LOAD, time_step)
        
        ! calculate ISI, BUI, FWI and DSR
        isi(i) = initial_spread_index(ffmc(i), wind(i))

        if (snow_depth(1) > 1e-3) then
            isi(i) = (1 - f_snow(i)) * isi(i)
        end if

        bui(i) = buildup_index(dmc(i), dc(i))

        fwi(i) = fire_weather_index(isi(i), bui(i))
        dsr(i) = daily_severity_rating(fwi(i))
    end do

  end subroutine cffwi_hourly_driver


end module GEOS_IgniGridCompMod
