!$Id: GEOS_DataAtmGridComp.F90,v 1.15.6.11.4.18.2.1.2.1.2.5.2.4.4.5.8.4.2.1 2021/05/21 20:59:48 atrayano Exp $

#include "MAPL_Generic.h"

!=============================================================================

module GEOS_DataAtmGridCompMod

!BOP

! !MODULE: GEOS_DataAtm -- A ``fake'' atmospheric component.

! !USES: 

  use ESMF
  use MAPL_Mod
! use ncar_ocean_fluxes_mod
  use GEOS_SurfaceGridCompMod,    only : SurfSetServices      => SetServices
  use GEOS_UtilsMod

  use ice_init,           only: dealloc_column_physics

  implicit none
  private

! !PUBLIC MEMBER FUNCTIONS:

  public SetServices

!=============================================================================

! !DESCRIPTION:
! 
!   {\tt GEOS\_DataAtm  } is a gridded component that ...??

!EOP

  integer            :: SURF

  integer            :: DO_OBIO         ! default (=0) is to run without ocean bio and chem
  integer, parameter :: DO_CO2SC  = 0
  logical, parameter :: DO_GOSWIM = .false.

  integer, parameter :: NUM_DUDP = 5
  integer, parameter :: NUM_DUWT = 5
  integer, parameter :: NUM_DUSD = 5
  integer, parameter :: NUM_BCDP = 2                           ! number of Black Carbon 
  integer, parameter :: NUM_BCWT = 2
  integer, parameter :: NUM_OCDP = 2                           ! number of Organic Carbon 
  integer, parameter :: NUM_OCWT = 2

  integer, parameter :: NB_CHOU_UV  = 5                        ! number of UV bands
  integer, parameter :: NB_CHOU_NIR = 3                        ! number of near-IR bands
  integer, parameter :: NB_CHOU     = NB_CHOU_UV + NB_CHOU_NIR ! total number of bands

  integer, parameter :: ICE   = 1
  integer, parameter :: WATER = 2
  integer, parameter :: OBIO  = 3

  real,    parameter :: KUVR = 0.09
  real,    parameter :: EPS6 = 1.e-7

  type bandptr
   real, pointer, dimension(:)  ::b => null() 
  end type bandptr

  character(len = 2) :: suffix
  character(len = 3) :: label
  integer k 

   contains

!BOP

! !IROUTINE: SetServices -- Sets ESMF services for this component

! !INTERFACE:

  subroutine SetServices ( GC, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    integer, optional                  :: RC  ! return code

!  !DESCRIPTION: This version uses the MAPL_GenericSetServices. This function sets
!                the Initialize and Finalize services, as well as allocating
!   our instance of a generic state and putting it in the 
!   gridded component (GC). Here we only need to set the run method and
!   add the state variable specifications (also generic) to our instance
!   of the generic state. This is the way our true state variables get into
!   the ESMF_State INTERNAL, which is in the MAPL_MetaComp.
!
!EOP

!=============================================================================
!
! ErrLog Variables


    character(len=ESMF_MAXSTR)              :: IAm
    integer                                 :: STATUS
    character(len=ESMF_MAXSTR)              :: COMP_NAME

! Local derived type aliases

    type (ESMF_Config)                      :: CF
    integer                                 :: DO_CICE_THERMO  ! default (=1) is to run with CICE
    type (MAPL_MetaComp),  pointer          :: MAPL

!=============================================================================

! Begin...

! Get my name and set-up traceback handle
! ---------------------------------------

    Iam = "SetServices"
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, CONFIG=CF, __RC__ )
    Iam = trim(COMP_NAME) // Iam

    call MAPL_AddImportSpec(GC,              &
      SHORT_NAME = 'PS',                     &
      LONG_NAME = 'surface_pressure',        &
      UNITS = 'Pa',                          &
      DIMS = MAPL_DimsHorzOnly,              &
      VLOCATION = MAPL_VLocationNone, __RC__)

    call MAPL_AddImportSpec(GC,              &
      SHORT_NAME = 'TA',                     &
      LONG_NAME = 'surface_air_temperature', &
      UNITS = 'K',                           &
      DIMS = MAPL_DimsHorzOnly,              &
      VLOCATION = MAPL_VLocationNone, __RC__)

    call MAPL_AddImportSpec(GC,              &
      SHORT_NAME = 'QA',                     &
      LONG_NAME = 'surface_specific_humidity', &
      UNITS = '1',                           &  ! convert to kg/kg?? or g/kg??
      DIMS = MAPL_DimsHorzOnly,              &
      VLOCATION = MAPL_VLocationNone, __RC__)

    call MAPL_AddImportSpec(GC,              &
      SHORT_NAME = 'UA',                     &
      LONG_NAME = '10-meter_eastward_wind',  &
      UNITS = 'm s-1',                       &
      DIMS = MAPL_DimsHorzOnly,              &
      VLOCATION = MAPL_VLocationNone, __RC__)

    call MAPL_AddImportSpec(GC,              &
      SHORT_NAME = 'VA',                     &
      LONG_NAME = '10-meter_northward_wind', &
      UNITS = 'm s-1',                       &
      DIMS = MAPL_DimsHorzOnly,              &
      VLOCATION = MAPL_VLocationNone, __RC__)

    call MAPL_AddImportSpec(GC,              &
      SHORT_NAME = 'RUNOFF',                 &
      LONG_NAME = 'overland_runoff_including_throughflow', &
      UNITS = 'kg m-2 s-1',                  &
      DIMS = MAPL_DimsHorzOnly,              &
      VLOCATION = MAPL_VLocationNone, __RC__)

    call MAPL_AddImportSpec(GC,              &
      SHORT_NAME = 'PCU',                    &
      LONG_NAME = 'convective_rainfall',     &
      UNITS = 'kg m-2 s-1',                  &
      DIMS = MAPL_DimsHorzOnly,              &
      VLOCATION = MAPL_VLocationNone, __RC__)

    call MAPL_AddImportSpec(GC,              &
      SHORT_NAME = 'PLS',                    &
      LONG_NAME = 'large_scale_rainfall',    &
      UNITS = 'kg m-2 s-1',                  &
      DIMS = MAPL_DimsHorzOnly,              &
      VLOCATION = MAPL_VLocationNone, __RC__)

    call MAPL_AddImportSpec(GC,              &
      SHORT_NAME = 'SNO',                    &
      LONG_NAME = 'snowfall',                &
      UNITS = 'kg m-2 s-1',                  &
      DIMS = MAPL_DimsHorzOnly,              &
      VLOCATION = MAPL_VLocationNone, __RC__)

    call MAPL_AddImportSpec(GC,              &
      SHORT_NAME = 'LWDN',               &
      LONG_NAME = 'open_water_downward_longwave_flux',                &
      UNITS = 'W m-2',                       &
      DIMS = MAPL_DimsHorzOnly,              &
      VLOCATION = MAPL_VLocationNone, __RC__)

    call MAPL_AddImportSpec(GC,              &
      SHORT_NAME = 'SWGDWN',                 &
      LONG_NAME = 'surface_incoming_shortwave_flux',&
      UNITS = 'W m-2',                       &
      DIMS = MAPL_DimsHorzOnly,              &
      VLOCATION = MAPL_VLocationNone, __RC__)

! Internal
    call MAPL_AddInternalSpec(GC,              &
      SHORT_NAME = 'TS',                     &
      LONG_NAME = 'surface_temperature', &
      UNITS = 'K',                           &
      DIMS = MAPL_DimsHorzOnly,              &
      DEFAULT = -1000.0,                     &
      VLOCATION = MAPL_VLocationNone, __RC__)

    call MAPL_AddInternalSpec(GC,              &
         SHORT_NAME         = 'EMIS',          &
         LONG_NAME          = 'surface_emissivity', &
         UNITS              = '1',              &
         DEFAULT            = 1.0,              &
         DIMS = MAPL_DimsHorzOnly,              &
         VLOCATION = MAPL_VLocationNone, __RC__)

! Get my MAPL_Generic state
!--------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, __RC__)

    call MAPL_GetResource ( MAPL, DO_CICE_THERMO, Label="USE_CICE_Thermo:" , DEFAULT=1, __RC__)

! Get constants from CF
! ---------------------

    if (DO_CICE_THERMO == 0) then
       if(MAPL_AM_I_ROOT()) then
          print *, 'Current DATA ATM Gridded Component'
          print *, 'Needs CICE Thermodynamics! You turned it off, so this run will now terminate!'
          ASSERT_(DO_CICE_THERMO == 0)
       endif
    endif

! Ocean biology and chemistry: using OBIO or not?
! ------------------------------------------------

    call MAPL_GetResource ( MAPL, DO_OBIO, Label="USE_OCEANOBIOGEOCHEM:", DEFAULT=0, __RC__)

! Set the Initialize and Run entry points
! ---------------------------------------

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_INITIALIZE, Initialize, __RC__)
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,        Run, __RC__)
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_FINALIZE,  Finalize, __RC__)

! Create children`s gridded components and invoke their SetServices
! -----------------------------------------------------------------
    SURF = MAPL_AddChild(GC, NAME='SURFACE', SS=SurfSetServices, __RC__)

! Set the state variable specs.
! -----------------------------

!BOS

!  !IMPORT STATE:
! none for now

!  !INTERNAL STATE:
! none for now

!  !EXPORT STATE:

!EOS

! Set generic init and final methods
! ----------------------------------

    call MAPL_TimerAdd(GC,    name="INITIALIZE", __RC__)
    call MAPL_TimerAdd(GC,    name="RUN"       , __RC__)
    call MAPL_TimerAdd(GC,    name="FINALIZE"  , __RC__)

    ! This call is needed only when we use ReadForcing.
    ! If we switch to use ExtData, next line has be commented out
    if (DO_CICE_THERMO == 2) then
       call MAPL_TerminateImport    ( GC, SHORT_NAMES=['SURFSTATE'],    &
                                      CHILD_IDS=[SURF],  __RC__  )
    else
       call MAPL_TerminateImport    ( GC, ALL=.true., __RC__ )
    endif

    call MAPL_GenericSetServices    ( GC, __RC__)

    RETURN_(ESMF_SUCCESS)
  
  end subroutine SetServices

!BOP

! !IROUTINE: INITIALIZE -- Initialize stage for the DataAtm component

! !INTERFACE:

subroutine INITIALIZE ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

  type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
  type(ESMF_State),    intent(inout) :: IMPORT ! Import state
  type(ESMF_State),    intent(inout) :: EXPORT ! Export state
  type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
  integer, optional,   intent(  out) :: RC     ! Error code:

! !DESCRIPTION: 

!EOP

! ErrLog Variables

  character(len=ESMF_MAXSTR)          :: IAm
  integer                             :: STATUS
  character(len=ESMF_MAXSTR)          :: COMP_NAME

! Locals

  type (MAPL_MetaComp), pointer   :: MAPL => null()
  type (ESMF_GridComp), pointer   :: GCS(:) => null()
  type (MAPL_LocStream)           :: LOCSTREAM
  type (MAPL_LocStream)           :: EXCH
  real, pointer :: tmp(:,:) ! needed only to force allocation
  type (ESMF_State),         pointer  :: GEX(:)
  type (ESMF_Alarm) :: alarm, solarAlarm

!  Begin...
!----------

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    Iam = "Initialize"
    call ESMF_GridCompGet( GC, name=COMP_NAME, __RC__)
    Iam = trim(COMP_NAME) // Iam

! Get my internal MAPL_Generic state
!----------------------------------

    call MAPL_GetObjectFromGC(GC, MAPL, STATUS)
    VERIFY_(STATUS)

! Start Total timer
!------------------

    call MAPL_TimerOn (MAPL,"TOTAL")
    call MAPL_TimerOn (MAPL,"INITIALIZE"  )

!!! ALT this section below is not needed
!+=======================================+
! Change the location stream to just the ocean part
!--------------------------------------------------

    call MAPL_Get (MAPL, GCS=GCS, __RC__ )

    call MAPL_Get(MAPL, EXCHANGEGRID=EXCH, __RC__ )

    call MAPL_LocStreamCreate(LOCSTREAM, EXCH, NAME='OCEAN', &
                                       MASK=(/MAPL_OCEAN/), __RC__ )

    call MAPL_ExchangeGridSet(GCS(SURF), LOCSTREAM, __RC__)

    call MAPL_TimerOff(MAPL,"INITIALIZE"  )
    call MAPL_TimerOff(MAPL,"TOTAL")

    call MAPL_GenericInitialize ( GC, IMPORT, EXPORT, CLOCK,  __RC__)

!ALT: At this point all the children (i.e. Surface) and grand-children 
! have been initialized
!
! we are now mimicking connections in Physics to force allocation
! for some Surface exports

    call MAPL_Get (MAPL, GEX=GEX, __RC__ )
    call MAPL_GetPointer(GEX(SURF), tmp, 'LWI'  ,  alloc=.true., __RC__)

! also we need to fake the "Solar alarm"
    call MAPL_Get (MAPL, runAlarm=alarm, __RC__ )
    solarAlarm = ESMF_AlarmCreate(alarm, __RC__)
    call ESMF_AlarmSet(solarAlarm, name='SOLAR_Alarm', __RC__)

    RETURN_(ESMF_SUCCESS)

  end subroutine INITIALIZE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!BOP

! ! IROUTINE: RUN -- Run stage for the DataAtm component

! !INTERFACE:

subroutine RUN ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

  type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
  type(ESMF_State),    intent(inout) :: IMPORT ! Import state
  type(ESMF_State),    intent(inout) :: EXPORT ! Export state
  type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
  integer, optional,   intent(  out) :: RC     ! Error code:

! ! DESCRIPTION: Periodically refreshes ...??

!EOP

! ErrLog Variables

  character(len=ESMF_MAXSTR)          :: IAm
  integer                             :: STATUS
  character(len=ESMF_MAXSTR)          :: COMP_NAME

! Locals

  type (MAPL_MetaComp),       pointer :: MAPL       => null()
  type (ESMF_GridComp),       pointer :: GCS(:)     => null()
  type (ESMF_State),          pointer :: GIM(:)     => null()
  type (ESMF_State),          pointer :: GEX(:)     => null()
  character(len=ESMF_MAXSTR), pointer :: GCNAMES(:) => null()
  type (ESMF_State),          pointer :: SurfImport => null()
  type (ESMF_State),          pointer :: SurfExport => null()
  type (ESMF_Time)                    :: currentTime

  integer                                   :: IM, JM
  real, dimension(:,:), allocatable         :: Uskin, Vskin, Qskin
  real, dimension(:,:), allocatable, target :: swrad
  real, dimension(:,:), pointer             :: PS, PSsurf, Tair, Qair, Uair, Vair, DZ
  real, dimension(:,:), pointer             :: ALW, BLW, SPEED, DISCHARGE, rPCU, rPLS, sSNO
  real, dimension(:,:), pointer             :: CT, CQ, CM, SH, EVAP, TAUX, TAUY, Tskin, lwdnsrf
  real, dimension(:,:), pointer             :: DRPARN, DFPARN, DRNIRN, DFNIRN, DRUVRN, DFUVRN
  real, dimension(:,:), pointer             :: DSH, DEVAP
  real, dimension(:,:), pointer             :: EMISSRF

  real, allocatable, dimension(:,:) :: ZTH
  real, allocatable, dimension(:,:) :: SLR
  real,    pointer, dimension(:,:)    :: LATS     => NULL()
  real,    pointer, dimension(:,:)    :: LONS     => NULL()
  real                                :: SC, MG, SB
  logical                             :: USE_NRLSSI2
  character(len=ESMF_MAXPATHLEN) :: SolCycFileName
  logical :: PersistSolar
  type (MAPL_SunOrbit)                :: ORBIT
  type (ESMF_TimeInterval)            :: DELT

! Andrea: do we need these????
  real, parameter :: FRUVR          = 0.07
  real, parameter :: FRPAR          = 0.40
  real, parameter :: FRNIR          = 0.53

  ! partitioning of shortwave radiation for ice 
  real, parameter :: FRVISDIR       = 0.29
  real, parameter :: FRVISDIF       = 0.31
  real, parameter :: FRNIRDIR       = 0.24
  real, parameter :: FRNIRDIF       = 0.16

  real, parameter ::  KUVR          = 0.09

  real, parameter :: alb            = 0.066


! pointers to import - none
! internal pointers to tile variables ???

  real, pointer, dimension(:,:) :: TA       ! => null()
  real, pointer, dimension(:,:) :: QA       ! => null()
  real, pointer, dimension(:,:) :: UA       ! => null()
  real, pointer, dimension(:,:) :: VA       ! => null()
  real, pointer, dimension(:,:) :: RUNOFF   ! => null()
  real, pointer, dimension(:,:) :: PCU      ! => null()
  real, pointer, dimension(:,:) :: PLS      ! => null()
  real, pointer, dimension(:,:) :: SNO      ! => null()
  real, pointer, dimension(:,:) :: LWDN ! => null()
  real, pointer, dimension(:,:) :: SWGDWN   ! => null()

! pointers to export - none???

! pointers to internal
  real, pointer, dimension(:,:) :: TS
  real, pointer, dimension(:,:) :: EMIS
  type(ESMF_STATE) :: internal

! Andrea: OBSERVE????
!  logical, dimension(1) :: OBSERVE
!  real LATSO, LONSO

   real, parameter :: HW_hack = 2.
   logical :: firsttime = .false.

   real :: TAU_TS
   real :: REF_HEIGHT
   real :: DT

   integer :: year, month, day, hr, mn, se

!  Begin...
!----------

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    Iam = "Run"
    call ESMF_GridCompGet( GC, name=COMP_NAME, __RC__ )
    Iam = trim(COMP_NAME) // Iam

! Get my MAPL_Generic (GG) state
!-------------------------------

    call MAPL_GetObjectFromGC(GC, MAPL, status)
    VERIFY_(status)

! Start timers
!-------------

    call MAPL_TimerOn(MAPL,"TOTAL")
    call MAPL_TimerOn(MAPL,"RUN" )

! Get the time step
! -----------------

! Get current time from clock
!----------------------------
    call ESMF_ClockGet(CLOCK, currTime=CurrentTime, __RC__)

    call esmf_timeget(CurrentTime, yy=year, mm=month, dd=day, h=hr, m=mn, s=se, rc = status);

    call MAPL_Get (MAPL, GCS=GCS, GIM=GIM, GEX=GEX, __RC__ )
    call MAPL_Get (MAPL, INTERNAL_ESMF_STATE = INTERNAL, __RC__)

    call MAPL_Get(MAPL, HEARTBEAT = DT, __RC__)
    call MAPL_GetResource ( MAPL, DT, Label="DT:", DEFAULT=DT, __RC__)
    call MAPL_GetResource ( MAPL, TAU_TS, Label="TAU_TS:", DEFAULT=7200.0, __RC__)
    call MAPL_GetResource ( MAPL, REF_HEIGHT, Label="REFERENCE_HEIGHT:", DEFAULT=10.0, __RC__)

! Pointers to Imports
!--------------------

! Pointers to Internals
!----------------------
    call MAPL_GetPointer(Internal, Tskin, 'TS', __RC__)
    call MAPL_GetPointer(Internal, EMIS, 'EMIS', __RC__)

!  Pointers to Exports ????
!---------------------

!new stuff: this is what Surface needs
!======================================

! Get children and their im/ex states from my generic state.
!----------------------------------------------------------
    call MAPL_Get ( MAPL, GCS=GCS, GIM=GIM, GEX=GEX, GCNAMES=GCNAMES, __RC__)
    call MAPL_Get(MAPL, IM=IM, JM=JM, LATS=LATS, LONS=LONS, ORBIT=ORBIT,__RC__)

    allocate( ZTH  (IM,JM), __STAT__)
    allocate( SLR  (IM,JM), __STAT__)

!  !IMPORT STATE:
    SurfImport => GIM(SURF)
    SurfExport => GEX(SURF)

! Read Sea Level Pressure (Pa)
!---------------------------------------------------
    !ALT is this default too low?
!    call ReadForcingData(impName='PS', frcName='SLP', default=90000., __RC__)
    call MAPL_GetPointer(import, PS, 'PS', __RC__)
    call MAPL_GetPointer(SurfImport, PSsurf, 'PS', __RC__)
    where(ps==0.0) ps=100000.
    PSsurf = PS 

! Read 10m temperature (K)
!---------------------------------------------------
!   call ReadForcingData(impName='TA', frcName='T10', default=290., __RC__) ! how about default T10 = 270+30*COS(LAT)
    call MAPL_GetPointer(SurfImport, Tair, 'TA', __RC__)
    call MAPL_GetPointer(import, TA, 'TA', __RC__)
    where(ta==0.0) ta=MAPL_Tice
    Tair = TA

! Read 10m specific humidity (kg kg-1)
!---------------------------------------------------
!   call ReadForcingData(impName='QA', frcName='Q10', default=2.0e-6, __RC__)
    call MAPL_GetPointer(SurfImport, Qair, 'QA', __RC__)
    call MAPL_GetPointer(import, QA, 'QA', __RC__)
!@@    Qair = QA * 0.001 ! SA: convert g/Kg to Kg/Kg
    Qair = QA

! Read 10m zonal wind speed (m s-1)
!---------------------------------------------------
!   call ReadForcingData(impName='UA', frcName='U10', default=0., __RC__)
    call MAPL_GetPointer(SurfImport, Uair, 'UA', __RC__)
    call MAPL_GetPointer(import, UA, 'UA', __RC__)
    Uair = UA
    Uair = merge(tsource = uair, fsource = 0.0, mask = (abs(uair) < 1000.0)); 

! Read 10m meridional wind speed (m s-1)
!---------------------------------------------------
!   call ReadForcingData(impName='VA', frcName='V10', default=0., __RC__)
    call MAPL_GetPointer(SurfImport, Vair, 'VA', __RC__)
    call MAPL_GetPointer(import, VA, 'VA', __RC__)
    Vair = VA
    Vair = merge(tsource = vair, fsource = 0.0, mask = (abs(vair) < 1000.0)); 

    call MAPL_GetPointer(SurfImport, SPEED, 'SPEED', __RC__)
    SPEED = SQRT(Uair**2 + Vair**2)

    IM = size(Uair, 1)
    JM = size(Uair, 2)
    allocate(Uskin(IM,JM), Vskin(IM,JM), Qskin(IM,JM), swrad(IM,JM), __STAT__)

    call MAPL_GetPointer(SurfImport, DZ, 'DZ', __RC__)
    DZ = REF_HEIGHT

! River runoff    
!   call ReadForcingData(impName='DISCHARGE', frcName='RR', default=0., __RC__)
    call MAPL_GetPointer(SurfImport, DISCHARGE, 'DISCHARGE', __RC__)
    call MAPL_GetPointer(import, RUNOFF, 'RUNOFF', __RC__)
    DISCHARGE=RUNOFF

    !ALT: we should read topo, but for now over ocean this is fine
    call SetVarToZero('PHIS', __RC__)

! these are not used, 0 should be OK
    call SetVarToZero('DEWL', __RC__)
    call SetVarToZero('FRSL', __RC__)

    call MAPL_GetPointer(SurfImport, DSH, 'DSH', __RC__)
    call MAPL_GetPointer(SurfImport, DEVAP, 'DEVAP', __RC__)
! these should be set to 0 (for now)
    !call SetVarToZero('DSH', __RC__)
    call SetVarToZero('DFU', __RC__)
    call SetVarToZero('DFV', __RC__)
    !call SetVarToZero('DEVAP', __RC__)
    call SetVarToZero('DDEWL', __RC__)
    call SetVarToZero('DFRSL', __RC__)

! Precipitation
! question for Andrea: the old code reads RAIN. How to partition into
! PCU and PLS? Right now all of RAIN goes into PCU (arbitralily chosen)
! SA: above can be squashed, we can read PCU and PLS and add to get RAIN.
!   call ReadForcingData(impName='PCU', frcName='RAIN', default=0., __RC__)
    call MAPL_GetPointer(SurfImport, rPCU, 'PCU', __RC__)
    call MAPL_GetPointer(import, PCU,    'PCU',    __RC__)
    rPCU=PCU

    call MAPL_GetPointer(SurfImport, rPLS, 'PLS', __RC__)
    call MAPL_GetPointer(import, PLS,    'PLS',    __RC__)
    rPLS=PLS

!   RAIN = PCU + PLS
!   call SetVarToZero('PLS', __RC__)

!   call ReadForcingData(impName='SNO', frcName='SNO', default=0., __RC__)
    call MAPL_GetPointer(SurfImport, sSNO, 'SNO', __RC__)
    call MAPL_GetPointer(import, SNO,    'SNO',    __RC__)
    sSNO=SNO

! Radiation
!   call ReadForcingData(impName='LWDNSRF', frcName='LWRAD', default=100.0, __RC__)
    call MAPL_GetPointer(import, LWDN, 'LWDN', __RC__)
    call MAPL_GetPointer(SurfImport, lwdnsrf, 'LWDNSRF', __RC__)

!   call ReadForcingData(swrad, frcName='SWRAD', default=0.200, __RC__)
    call MAPL_GetPointer(import, SWGDWN, 'SWGDWN', __RC__)
    swrad = SWGDWN

    ! get sw at the TOA (next few lines are copy-and-paste from Surface)

! Get the insolation and zenith angle on grid and tiles
!------------------------------------------------------

    call ESMF_ClockGet(CLOCK,     TIMESTEP=DELT, __RC__)
    ! the line below only works for daily forcing e..g. CORE I
    ! for JRA55-DO or any dataset at higher frequency, this line makes SW much
    ! higher than what data prescribed
    !DELT = DELT * NINT((86400./DT)) ! emulate daily Solar

    call MAPL_SunGetInsolation(LONS, LATS,      &
              ORBIT, ZTH, SLR, &
              INTV  = DELT,    &
              CLOCK = CLOCK,   &
              __RC__ )

    call MAPL_GetResource( MAPL, SC, 'SOLAR_CONSTANT:', __RC__)
    call MAPL_GetResource( MAPL, SolCycFileName, "SOLAR_CYCLE_FILE_NAME:", DEFAULT='/dev/null', __RC__)

    if(SolCycFileName /= '/dev/null') THEN

       call MAPL_GetResource( MAPL, USE_NRLSSI2, "USE_NRLSSI2:", DEFAULT=.TRUE., __RC__)

       if (USE_NRLSSI2) then
          call MAPL_GetResource( MAPL, PersistSolar, "PERSIST_SOLAR:", DEFAULT=.TRUE., __RC__)
          call MAPL_SunGetSolarConstant(CLOCK,trim(SolCycFileName),SC,MG,SB,PersistSolar=PersistSolar,__RC__)
       else
          call MAPL_SunGetSolarConstant(CLOCK,trim(SolCycFileName),SC,__RC__)
       endif
    else if(SC<0.0) then
       call MAPL_SunGetSolarConstant(CURRENTTIME,SC,__RC__)
    end if

!    where (zth <= 0.0 .or. slr == 0.0)
    where (zth <= 1.0e-6)
       swrad = 0.0
    elsewhere
       swrad = swrad / (sc*slr) 
    end where

    call MAPL_GetPointer(SurfImport, DRPARN, 'DRPARN', __RC__)
    call MAPL_GetPointer(SurfImport, DFPARN, 'DFPARN', __RC__)
    call MAPL_GetPointer(SurfImport, DRNIRN, 'DRNIRN', __RC__)
    call MAPL_GetPointer(SurfImport, DFNIRN, 'DFNIRN', __RC__)
    call MAPL_GetPointer(SurfImport, DRUVRN, 'DRUVRN', __RC__)
    call MAPL_GetPointer(SurfImport, DFUVRN, 'DFUVRN', __RC__)

! Andrea: Is partitioning OK here? Or something else?
! Direct vs duffuse, zenith angle, time of the day???
    DRPARN = swrad*FRPAR*0.6
    DFPARN = swrad*FRPAR*0.4
    DRUVRN = swrad*FRUVR*0.6
    DFUVRN = swrad*FRUVR*0.4
    DRNIRN = swrad*FRNIR*0.6
    DFNIRN = swrad*FRNIR*0.4

!@@   print *,'DEBUG:min/max DRPARN',minval(DRPARN),maxval(DRPARN)

! Andrea: Tskin = TS? From previous time step?
! Andrea: do we need radiation vars computed before we call phase 1
!@    call MAPL_GetPointer(SurfExport, Tskin, 'TS_FOUND', __RC__)
!@    call MAPL_GetPointer(SurfExport, Tskin, 'TS', __RC__)

    call MAPL_GetPointer(SurfImport, ALW, 'ALW', __RC__)
    call MAPL_GetPointer(SurfImport, BLW, 'BLW', __RC__)

    if(any(Tskin<0.0)) then !only when DATAATM restart is bootstrapped
       firsttime = .true.
    end if

    if (firsttime) then
       firsttime = .false.
       Tskin = TA
    end if

   BLW = EMIS*4*MAPL_STFBOL * Tskin**3
   ALW = EMIS*MAPL_STFBOL * Tskin**4 - BLW*Tskin

   call MAPL_GetPointer(SurfExport, EMISSRF, 'EMIS', alloc=.true., __RC__)

   lwdnsrf = lwdn
!   lwdnsrf = 0.0
 !  where(LWGNTWTR /=  MAPL_Undef)
 !     lwdnsrf = -(LWGNTWTR - MAPL_STFBOL * Tskin**4)
 !  end where

!    call SetVarToZero('BLW', __RC__)
    
    call SetVarToZero('DTSDT', __RC__)

!!!    if (mapl_am_i_root()) PRINT*, __FILE__, __LINE__
! call Run (or, phase) 1 of Surface
    call ESMF_GridCompRun (GCS(SURF), importState=GIM(SURF), &
         exportState=GEX(SURF), clock=CLOCK, PHASE=1, userRC=status )
    VERIFY_(status)

! we deal with these after Run1 of surf
    call MAPL_GetPointer(SurfImport, Tair, 'TA', __RC__)
    call MAPL_GetPointer(SurfImport, Qair, 'QA', __RC__)
    call MAPL_GetPointer(SurfImport, Uair, 'UA', __RC__)
    call MAPL_GetPointer(SurfImport, Vair, 'VA', __RC__)
    call MAPL_GetPointer(SurfImport, PS,   'PS', __RC__)

    call MAPL_GetPointer(SurfExport, CT, 'CT', __RC__)
    call MAPL_GetPointer(SurfExport, CQ, 'CQ', __RC__)
    call MAPL_GetPointer(SurfExport, CM, 'CM', __RC__)

! Andrea:    -- Tskin can be the SST for this purpose – T foundation   - i think the export version is called HFLX

! Andrea: are you sure about this??? Or is GEOS_QsatLQU(Tskin, PS)

    Qskin = GEOS_Qsat(Tskin, PS, RAMP=0.0, PASCALS=.TRUE.)

    ! here we assume that the wind at the skin is zero
    Uskin = 0.0
    Vskin = 0.0

    call MAPL_GetPointer(SurfImport, SH, 'SH', __RC__)
    call MAPL_GetPointer(SurfImport, EVAP, 'EVAP', __RC__)
    call MAPL_GetPointer(SurfImport, TAUX, 'TAUX', __RC__)
    call MAPL_GetPointer(SurfImport, TAUY, 'TAUY', __RC__)

    SH = CT * MAPL_CP * (Tskin - Tair)
    EVAP = CQ * (Qskin - Qair) 
    TAUX = CM * (Uskin - Uair)
    TAUY = CM * (Vskin - Vair)

    ! these derivatives are important for sea ice
    DSH = CT !* MAPL_CP (MAPL_CP got multiplied in Surf)
    DEVAP = CQ

101 format (A, e20.12, 3I3.2)

!!!    if (mapl_am_i_root()) PRINT*, __FILE__, __LINE__
! call Run (or, phase) 2 of Surface
    call ESMF_GridCompRun (GCS(SURF), importState=GIM(SURF), &
         exportState=GEX(SURF), clock=CLOCK, PHASE=2, userRC=status )
    VERIFY_(status)

! by now Saltwater should be able to provide everything the ocean needs
    call MAPL_GetPointer(SurfExport, TS, 'TS', __RC__)

   WHERE (TS /= MAPL_Undef)   Tskin = Tskin + (TS - Tskin) * DT / (DT + TAU_TS)
   WHERE (EMISSRF /= MAPL_Undef) EMIS = EMISSRF

!-- still left to do
! modify GCM to always get the skin from Saltwater
! modify Surface so we can pass discharge, salinity, etc

! === CHECK_EMAIL 12/22/20

! Andrea: sea surface saloinity???,  river runoff???
! Read sea surface salinity (psu),
!---------------------------------------------------

!  do anything that relates to OBIO via call to obio_extra
!  ONLY if the user wants to run with OBIO, default is not to run it; see ABOVE SetServices
   call obio_extra(DO_OBIO, DO_CO2SC, DO_GOSWIM,    &
                   NUM_DUDP, NUM_DUWT, NUM_DUSD)
! ------------------------------------------------

    deallocate(Uskin, Vskin, Qskin, swrad, __STAT__)

!  All done
!-----------

    call MAPL_TimerOff(MAPL,"RUN"  )
    call MAPL_TimerOff(MAPL,"TOTAL")

    RETURN_(ESMF_SUCCESS)

contains
  subroutine ReadForcingData(dataptr, impName, frcName, default, rc)
    ! note that the variables MAPL, SurfImport and currentTime are "borrowed"
    ! from the parent
    real, pointer,    optional, intent(IN) :: dataptr(:,:)
    character(len=*), optional, intent(IN) :: impName
    character(len=*), optional, intent(IN) :: frcName
    real, intent(IN) :: default
    integer, optional :: rc

    ! local vars
    integer :: status
    real, pointer :: ptr(:,:) => null()
    character(len=ESMF_MAXSTR) :: datafile, frcName_, label_

    if (present(dataptr)) then
       if (mapl_am_i_root()) PRINT*, trim(impName)
       ASSERT_(.not.present(impName))
       ptr => dataPtr
    else
       call MAPL_GetPointer(SurfImport, ptr, impName, __RC__)
    end if
    if (present(frcName)) then
       frcName_ = frcName
    else
       ASSERT_(present(impName))
       frcName_ = impName
    end if
    label_ = trim(frcName_)//'_FILE:'

    call MAPL_GetResource(MAPL, datafile, label=label_, default='none', __RC__)
    if(trim(datafile) == 'none') then
       ptr = default
    else
       call MAPL_ReadForcing(MAPL, frcName_, renamefile(datafile, time = currenttime), currenttime, ptr, __RC__)
    endif 

    RETURN_(ESMF_SUCCESS)
  end subroutine ReadForcingData

  subroutine SetVarToZero(impName, rc)
    ! note that the variable SurfImport is "borrowed" from the parent
    character(len=*), intent(IN) :: impName
    integer, optional :: rc

    ! local vars
    integer :: status
    real, pointer :: ptr(:,:) => null()

    call MAPL_GetPointer(SurfImport, ptr, impName, __RC__)
    ptr = 0.0

    RETURN_(ESMF_SUCCESS)
  end subroutine SetVarToZero

  subroutine obio_extra(DO_OBIO, DO_CO2SC, DO_GOSWIM, &
                        NUM_DUDP, NUM_DUWT, NUM_DUSD)

    integer, intent(IN) :: NUM_DUDP, NUM_DUWT, NUM_DUSD
    integer, intent(IN) :: DO_OBIO, DO_CO2SC
    logical, intent(IN) :: DO_GOSWIM

    integer :: K

! Andrea: dry and wet clay depositions???
#ifdef LEFTOVER_FROM_OLD_CODE

! Andrea: do we need this?
    if((DO_OBIO/=0).OR. (DO_CO2SC /= 0)) then
!            SHORT_NAME         = 'CO2SC',                             &
                 ! ungridded dims NB_CHOU
!            SHORT_NAME         = 'FSWBAND',                           &
!            SHORT_NAME         = 'FSWBANDNA',                         &
    endif

!ALT: this is not done yet.
!SA:  seems we do not need it?!
    if (DO_GOSWIM) then
       ! BUNDLE
!        SHORT_NAME         = 'AERO_DP',                           &
    end if

!   Read Clay-Sized Dry Atmospheric Dust Depositions
!-------------------------------------------------
    do K = 1, NUM_DUDP
       write(label,'(I3.3)') K
       call MAPL_GetResource( MAPL, DATAfile, LABEL='DUDP'//label//'_FILE:', default = 'none', __RC__ )
       if(trim(datafile) == 'none') then; dry_clay = 0.0
       else
          call MAPl_ReadForcing( MAPL, 'DUDP'//label, DATAFILE, CURRENTTIME, dry_clay, __RC__ )
       endif
       if (associated(dry_clayx)) dry_clayx(:,K) = dry_clay
    end do

!   Read Clay-Sized Wet Atmospheric Dust Depositions
!-------------------------------------------------
    do K = 1, NUM_DUWT
       write(label,'(I3.3)') K
       call MAPL_GetResource( MAPL, DATAfile, LABEL='DUWT'//label//'_FILE:', default = 'none', __RC__ )
       if(trim(datafile) == 'none') then; wet_clay = 0.0
       else
          call MAPl_ReadForcing( MAPL, 'DUWT'//label, DATAFILE, CURRENTTIME, wet_clay, __RC__ )
       endif
       if (associated(wet_clayx)) wet_clayx(:,K) = wet_clay
    end do

!   Read Clay-Sized Sedimentary Atmospheric Dust Depositions
!---------------------------------------------------------
    do K = 1, NUM_DUSD
       write(label,'(I3.3)') K
       call MAPL_GetResource( MAPL, DATAfile, LABEL='DUSD'//label//'_FILE:', default = 'none', __RC__ )
       if(trim(datafile) == 'none') then; sed_clay = 0.0
       else
          call MAPl_ReadForcing( MAPL, 'DUSD'//label, DATAFILE, CURRENTTIME, sed_clay, __RC__ )
       endif
       if (associated(sed_clayx)) sed_clayx(:,K) = sed_clay
    end do

!   Read Atmospheric Clouds (Atmospheric Optics)
!---------------------------------------------
    call MAPL_GetResource( MAPL, DATAfile, LABEL='CCOVM_FILE:', default = 'none', __RC__ )
    if(trim(datafile) == 'none') then; ccovm = 0.0; 
    else; call MAPl_ReadForcing( MAPL, 'CCOVM', DATAFILE, CURRENTTIME, ccovm, __RC__ )
    endif; 

    call MAPL_GetResource( MAPL, DATAfile, LABEL='CLDTCM_FILE:', default = 'none', __RC__ )
    if(trim(datafile) == 'none') then; cldtcm = 0.0; 
    else; call MAPl_ReadForcing( MAPL, 'CLDTCM', DATAFILE, CURRENTTIME, cldtcm, __RC__ )
    endif; 

    call MAPL_GetResource( MAPL, DATAfile, LABEL='RLWPM_FILE:', default = 'none', __RC__ )
    if(trim(datafile) == 'none') then; rlwpm = 0.0; 
    else; call MAPl_ReadForcing( MAPL, 'RLWPM', DATAFILE, CURRENTTIME, rlwpm, __RC__ )
    endif; 

    call MAPL_GetResource( MAPL, DATAfile, LABEL='CDREM_FILE:', default = 'none', __RC__ )
    if(trim(datafile) == 'none') then; cdrem = 0.0; 
    else; call MAPl_ReadForcing( MAPL, 'CDREM', DATAFILE, CURRENTTIME, cdrem, __RC__ )
    endif; 

!   Read Atmospheric Properties (Atmospheric Optics)
!-------------------------------------------------
    call MAPL_GetResource( MAPL, DATAfile, LABEL='RH_FILE:', default = 'none', __RC__ )
    if(trim(datafile) == 'none') then; rh = 0.0; 
    else; call MAPl_ReadForcing( MAPL, 'RH', DATAFILE, CURRENTTIME, rh, __RC__ )
    endif; 

    call MAPL_GetResource( MAPL, DATAfile, LABEL='OZ_FILE:', default = 'none', __RC__ )
    if(trim(datafile) == 'none') then; oz = 0.0; 
    else; call MAPl_ReadForcing( MAPL, 'OZ', DATAFILE, CURRENTTIME, oz, __RC__ )
    endif; 

    call MAPL_GetResource( MAPL, DATAfile, LABEL='WV_FILE:', default = 'none', __RC__ )
    if(trim(datafile) == 'none') then; wv = 0.0; 
    else; call MAPl_ReadForcing( MAPL, 'WV', DATAFILE, CURRENTTIME, wv, __RC__ )
    endif; 

!   Read Atmospheric Carbon Dioxide from Carbon Tracker (_2011_OI)
!-----------------------------------------------------
    call MAPL_GetResource( MAPL, DATAfile, LABEL='CO2SC_FILE:', default = 'none', __RC__ )
    if(trim(datafile) == 'none') then; co2sc = 0.0; 
    else; call MAPl_ReadForcing( MAPL, 'CO2SC', DATAFILE, CURRENTTIME, co2sc,  __RC__ )
    endif;
    if ( associated(co2scx) ) co2scx = co2sc


!   Read MODIS Aerosols (Atmospheric Optics)
!-----------------------------------------

    do k=1, 33
     write(unit = suffix, fmt = '(i2.2)') k
     call MAPL_GetResource( MAPL, DATAfile, LABEL='TAUA_FILE:', default = 'none', __RC__ )
     if(trim(datafile) == 'none') then; taua = 0.0
     else; call MAPL_ReadForcing( MAPL, 'TAUA_' // suffix, trim(DATAFILE) // suffix, CURRENTTIME, taua, __RC__)
     endif;
     ataua(k)%b => taua

     call MAPL_GetResource( MAPL, DATAfile, LABEL='ASYMP_FILE:', default = 'none', __RC__ )
     if(trim(datafile) == 'none') then; asymp = 0.0
     else; call MAPL_ReadForcing( MAPL, 'ASYMP_' // suffix, trim(DATAFILE) // suffix, CURRENTTIME, asymp, __RC__)
     endif;  
     aasymp(k)%b => asymp

     call MAPL_GetResource( MAPL, DATAfile, LABEL='SSALB_FILE:', default = 'none', __RC__ )
     if(trim(datafile) == 'none') then; ssalb = 0.0
     else; call MAPL_ReadForcing( MAPL, 'SSALB_' // suffix, trim(DATAFILE) // suffix, CURRENTTIME, ssalb, __RC__)
     endif;
     assalb(k)%b => ssalb
    enddo
#endif

    RETURN_(ESMF_SUCCESS)
  end subroutine obio_extra

end subroutine RUN
!----------------------------------------------------------------------------------------------------------------------------------

    function renamefile(name, time) result(name0)

      character(len = *), intent(in) :: name;
      type(esmf_time), intent(inout) :: time;
      character(len = len(name)) :: name0; 

      integer :: year, month, day, status, i;
          
        name0 = trim(name);
        i = index(string = name, substring = "yyyymmdd");
        if(i == 0) return;

        call esmf_timeget(time, yy = year, mm = month, dd = day, rc = status); 
        write(unit = name0(i:i + 7), fmt = "(i4,i2.2,i2.2)") year, month, day;

    end function;  

!----------------------------------------------------------------------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! !IROUTINE: Finalize        -- Finalize method for CICEThermo wrapper

! !INTERFACE:

  subroutine Finalize ( gc, import, export, clock, rc ) 

! !ARGUMENTS:

  type(ESMF_GridComp), intent(INOUT) :: gc     ! Gridded component 
  type(ESMF_State),    intent(INOUT) :: import ! Import state
  type(ESMF_State),    intent(INOUT) :: export ! Export state
  type(ESMF_Clock),    intent(INOUT) :: clock  ! The supervisor clock
  integer, optional,   intent(  OUT) :: rc     ! Error code:

!EOP

    type (MAPL_MetaComp), pointer:: MAPL 

! ErrLog Variables

    character(len=ESMF_MAXSTR)       :: IAm
    integer                          :: STATUS
    character(len=ESMF_MAXSTR)       :: COMP_NAME

    integer                          :: DO_CICE_THERMO  ! default (=0) is to run without CICE

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    Iam = "Finalize"
    call ESMF_GridCompGet( gc, NAME=comp_name, __RC__)

! Get my internal MAPL_Generic state
!-----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, __RC__)

    call MAPL_GetResource ( MAPL, DO_CICE_THERMO, Label="USE_CICE_Thermo:" , DEFAULT=0, __RC__)

! Profilers
!----------

    call MAPL_TimerOn(MAPL,"TOTAL"   )
    call MAPL_TimerOn(MAPL,"FINALIZE")

    if (DO_CICE_THERMO == 1) call dealloc_column_physics( MAPL_AM_I_Root(), Iam )

    call MAPL_TimerOff(MAPL,"FINALIZE")
    call MAPL_TimerOff(MAPL,"TOTAL"   )

! Generic Finalize
! ------------------
    
    call MAPL_GenericFinalize( GC, IMPORT, EXPORT, CLOCK, __RC__)


! All Done
!---------

    RETURN_(ESMF_SUCCESS)
  end subroutine Finalize

end module GEOS_DataAtmGridCompMod
