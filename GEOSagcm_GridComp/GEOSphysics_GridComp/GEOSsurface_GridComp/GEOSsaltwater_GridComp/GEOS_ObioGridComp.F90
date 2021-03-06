!  $Id$

#include "MAPL_Generic.h"

!=============================================================================
module GEOS_ObioGridCompMod

! !USES:

  use ESMF
  use MAPL

  implicit none
  private

! !PUBLIC MEMBER FUNCTIONS:

  public SetServices

! Following could also be controlled via resource parameter
  integer, parameter :: NUM_DUDP = 5                           ! number of DUst Depositions
  integer, parameter :: NUM_DUWT = 5
  integer, parameter :: NUM_DUSD = 5
  integer, parameter :: NUM_BCDP = 2                           ! number of Black Carbon 
  integer, parameter :: NUM_BCWT = 2
  integer, parameter :: NUM_OCDP = 2                           ! number of Organic Carbon 
  integer, parameter :: NUM_OCWT = 2

  integer, parameter :: NB_CHOU_UV  = 5                        ! number of UV bands
  integer, parameter :: NB_CHOU_NIR = 3                        ! number of near-IR bands
  integer, parameter :: NB_CHOU     = NB_CHOU_UV + NB_CHOU_NIR ! total number of bands

  contains

! !IROUTINE: SetServices -- Sets ESMF services for this component

! !INTERFACE:

  subroutine SetServices ( GC, RC )

    !ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    integer, optional                  :: RC  ! return code

! ErrLog Variables

    character(len=ESMF_MAXSTR)              :: IAm
    integer                                 :: STATUS
    character(len=ESMF_MAXSTR)              :: COMP_NAME

! Local derived type aliases

    type (MAPL_MetaComp),  pointer          :: MAPL
    type (ESMF_Config)                      :: CF

! Get my name and set-up traceback handle
! ---------------------------------------

    call ESMF_GridCompGet( GC, NAME=COMP_NAME, CONFIG=CF, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // 'SetServices'

! Get my MAPL_Generic state
!--------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,  Run1, RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,  Run2, RC=STATUS )
    VERIFY_(STATUS)

! Set the state variable specs.
! -----------------------------

! Following OBIO related exports 
! "passing thru" from atmosphere to ocean, no computation is otherwise done with (on) them.

    call MAPL_AddExportSpec(GC                            ,&
          SHORT_NAME         = 'CO2SC',                     &
          LONG_NAME          = 'CO2 Surface Concentration Bin 001',&
          UNITS              = '1e-6'                      ,&
          DIMS               = MAPL_DimsTileOnly           ,&
          VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
     VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                    &
          SHORT_NAME         = 'DUDP'                      ,&
          LONG_NAME          = 'Dust Dry Deposition'       ,&
          UNITS              = 'kg m-2 s-1'                ,&
          DIMS               = MAPL_DimsTileOnly           ,&
          UNGRIDDED_DIMS     = (/NUM_DUDP/)                ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          RC=STATUS  ) 
     VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                    &
          SHORT_NAME         = 'DUWT'                      ,&
          LONG_NAME          = 'Dust Wet Deposition'       ,&
          UNITS              = 'kg m-2 s-1'                ,&
          DIMS               = MAPL_DimsTileOnly           ,&
          UNGRIDDED_DIMS     = (/NUM_DUWT/)                ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          RC=STATUS  ) 
     VERIFY_(STATUS)
 
    call MAPL_AddExportSpec(GC,                    &
          SHORT_NAME         = 'DUSD'                      ,&
          LONG_NAME          = 'Dust Sedimentation'        ,&
          UNITS              = 'kg m-2 s-1'                ,&
          DIMS               = MAPL_DimsTileOnly           ,&
          UNGRIDDED_DIMS     = (/NUM_DUSD/)                ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          RC=STATUS  ) 
     VERIFY_(STATUS)
     
    call MAPL_AddExportSpec(GC,                    &
          SHORT_NAME         = 'BCDP'                            ,&
          LONG_NAME          = 'Black Carbon Dry Deposition'     ,&
          UNITS              = 'kg m-2 s-1'                      ,&
          DIMS               = MAPL_DimsTileOnly                 ,&
          UNGRIDDED_DIMS     = (/NUM_BCDP/)                      ,&
          VLOCATION          = MAPL_VLocationNone                ,&
          RC=STATUS  ) 
     VERIFY_(STATUS)
     
    call MAPL_AddExportSpec(GC,                    &
          SHORT_NAME         = 'BCWT'                            ,&
          LONG_NAME          = 'Black Carbon Wet Deposition'     ,&
          UNITS              = 'kg m-2 s-1'                      ,&
          DIMS               = MAPL_DimsTileOnly                 ,&
          UNGRIDDED_DIMS     = (/NUM_BCWT/)                      ,&
          VLOCATION          = MAPL_VLocationNone                ,&
          RC=STATUS  ) 
     VERIFY_(STATUS)
     
    call MAPL_AddExportSpec(GC,                    &
          SHORT_NAME         = 'OCDP'                            ,&
          LONG_NAME          = 'Organic Carbon Dry Deposition'   ,&
          UNITS              = 'kg m-2 s-1'                      ,&
          DIMS               = MAPL_DimsTileOnly                 ,&
          UNGRIDDED_DIMS     = (/NUM_OCDP/)                      ,&
          VLOCATION          = MAPL_VLocationNone                ,&
          RC=STATUS  ) 
     VERIFY_(STATUS)
     
    call MAPL_AddExportSpec(GC,                    &
          SHORT_NAME         = 'OCWT'                            ,&
          LONG_NAME          = 'Organic Carbon Wet Deposition'   ,&
          UNITS              = 'kg m-2 s-1'                      ,&
          DIMS               = MAPL_DimsTileOnly                 ,&
          UNGRIDDED_DIMS     = (/NUM_OCWT/)                      ,&
          VLOCATION          = MAPL_VLocationNone                ,&
          RC=STATUS  ) 
     VERIFY_(STATUS)
     
    call MAPL_AddExportSpec(GC,                    &
          SHORT_NAME         = 'FSWBAND'                         ,                   &
          LONG_NAME          = 'net_surface_downward_shortwave_flux_per_band_in_air',&
          UNITS              = 'W m-2'                           ,&
          DIMS               = MAPL_DimsTileOnly                 ,&
          UNGRIDDED_DIMS     = (/NB_CHOU/)                       ,&
          VLOCATION          = MAPL_VLocationNone                ,&
          RC=STATUS  ) 
     VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                    &
          SHORT_NAME         = 'FSWBANDNA'                       ,                                       &
          LONG_NAME          = 'net_surface_downward_shortwave_flux_per_band_in_air_assuming_no_aerosol',&
          UNITS              = 'W m-2'                           ,&
          DIMS               = MAPL_DimsTileOnly                 ,&
          UNGRIDDED_DIMS     = (/NB_CHOU/)                       ,&
          VLOCATION          = MAPL_VLocationNone                ,&
          RC=STATUS  ) 
     VERIFY_(STATUS)

! Following OBIO related imports are
! "passing thru" from atmosphere to ocean, no computation is otherwise done with (on) them.

   call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'CO2SC',                             &
        LONG_NAME          = 'CO2 Surface Concentration Bin 001', &
        UNITS              = '1e-6',                              &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        RESTART            = MAPL_RestartSkip,                    &
                                                       RC=STATUS  )
   VERIFY_(STATUS)

   call MAPL_AddImportSpec(GC,                            &
         SHORT_NAME         = 'DUDP'                      ,&
          LONG_NAME          = 'Dust Dry Deposition'       ,&
          UNITS              = 'kg m-2 s-1'                ,&
          DIMS               = MAPL_DimsTileOnly           ,&
          UNGRIDDED_DIMS     = (/NUM_DUDP/)                ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          RESTART            = MAPL_RestartSkip            ,&
          RC=STATUS  ) 
   VERIFY_(STATUS)

   call MAPL_AddImportSpec(GC,                            &
         SHORT_NAME         = 'DUWT'                      ,&
          LONG_NAME          = 'Dust Wet Deposition'       ,&
          UNITS              = 'kg m-2 s-1'                ,&
          DIMS               = MAPL_DimsTileOnly           ,&
          UNGRIDDED_DIMS     = (/NUM_DUWT/)                ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          RESTART            = MAPL_RestartSkip            ,&
          RC=STATUS  ) 
   VERIFY_(STATUS)

   call MAPL_AddImportSpec(GC,                            &
         SHORT_NAME         = 'DUSD'                      ,&
          LONG_NAME          = 'Dust Sedimentation'        ,&
          UNITS              = 'kg m-2 s-1'                ,&
          DIMS               = MAPL_DimsTileOnly           ,&
          UNGRIDDED_DIMS     = (/NUM_DUSD/)                ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          RESTART            = MAPL_RestartSkip            ,&
          RC=STATUS  ) 
   VERIFY_(STATUS)

   call MAPL_AddImportSpec(GC,                                  &
         SHORT_NAME         = 'BCDP'                            ,&
          LONG_NAME          = 'Black Carbon Dry Deposition'     ,&
          UNITS              = 'kg m-2 s-1'                      ,&
          DIMS               = MAPL_DimsTileOnly                 ,&
          UNGRIDDED_DIMS     = (/NUM_BCDP/)                      ,&
          VLOCATION          = MAPL_VLocationNone                ,&
          RESTART            = MAPL_RestartSkip                  ,&
          RC=STATUS  ) 
   VERIFY_(STATUS)

   call MAPL_AddImportSpec(GC,                                  &
         SHORT_NAME         = 'BCWT'                            ,&
          LONG_NAME          = 'Black Carbon Wet Deposition'     ,&
          UNITS              = 'kg m-2 s-1'                      ,&
          DIMS               = MAPL_DimsTileOnly                 ,&
          UNGRIDDED_DIMS     = (/NUM_BCWT/)                      ,&
          VLOCATION          = MAPL_VLocationNone                ,&
          RESTART            = MAPL_RestartSkip                  ,&
          RC=STATUS  ) 
   VERIFY_(STATUS)

   call MAPL_AddImportSpec(GC,                                  &
         SHORT_NAME         = 'OCDP'                            ,&
          LONG_NAME          = 'Organic Carbon Dry Deposition'   ,&
          UNITS              = 'kg m-2 s-1'                      ,&
          DIMS               = MAPL_DimsTileOnly                 ,&
          UNGRIDDED_DIMS     = (/NUM_OCDP/)                      ,&
          VLOCATION          = MAPL_VLocationNone                ,&
          RESTART            = MAPL_RestartSkip                  ,&
          RC=STATUS  ) 
   VERIFY_(STATUS)

   call MAPL_AddImportSpec(GC,                                  &
         SHORT_NAME         = 'OCWT'                            ,&
          LONG_NAME          = 'Organic Carbon Wet Deposition'   ,&
          UNITS              = 'kg m-2 s-1'                      ,&
          DIMS               = MAPL_DimsTileOnly                 ,&
          UNGRIDDED_DIMS     = (/NUM_OCWT/)                      ,&
          VLOCATION          = MAPL_VLocationNone                ,&
          RESTART            = MAPL_RestartSkip                  ,&
          RC=STATUS  ) 
   VERIFY_(STATUS)

   call MAPL_AddImportSpec(GC,                    &
         SHORT_NAME         = 'FSWBAND'                         ,                   &
          LONG_NAME          = 'net_surface_downward_shortwave_flux_per_band_in_air',&
          UNITS              = 'W m-2'                           ,&
          DIMS               = MAPL_DimsTileOnly                 ,&
          UNGRIDDED_DIMS     = (/NB_CHOU/)                       ,&
          VLOCATION          = MAPL_VLocationNone                ,&
          RESTART            = MAPL_RestartSkip                  ,&
          RC=STATUS  ) 
   VERIFY_(STATUS)

   call MAPL_AddImportSpec(GC,                    &
         SHORT_NAME         = 'FSWBANDNA'                       ,                                       &
          LONG_NAME          = 'net_surface_downward_shortwave_flux_per_band_in_air_assuming_no_aerosol',&
          UNITS              = 'W m-2'                           ,&
          DIMS               = MAPL_DimsTileOnly                 ,&
          UNGRIDDED_DIMS     = (/NB_CHOU/)                       ,&
          VLOCATION          = MAPL_VLocationNone                ,&
          RESTART            = MAPL_RestartSkip                  ,&
          RC=STATUS  ) 
   VERIFY_(STATUS)

! Set the Profiling timers
! ------------------------

    call MAPL_TimerAdd(GC,    name="RUN1"   ,               RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="RUN2"  ,                RC=STATUS)
    VERIFY_(STATUS)

! Set generic init and final methods
! ----------------------------------

    call MAPL_GenericSetServices    ( GC,  RC=STATUS )
    VERIFY_(STATUS)

! Set the Run entry point
! -----------------------

    RETURN_(ESMF_SUCCESS)

  end subroutine SetServices

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! !IROUTINE: RUN1
! !INTERFACE:

subroutine RUN1 ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:
  type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
  type(ESMF_State),    intent(inout) :: IMPORT ! Import state
  type(ESMF_State),    intent(inout) :: EXPORT ! Export state
  type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
  integer, optional,   intent(  out) :: RC     ! Error code:

! !DESCRIPTION: Periodically refreshes the sea-surface conditions

!EOP


! ErrLog Variables

  character(len=ESMF_MAXSTR)      :: IAm
  integer                         :: STATUS
  character(len=ESMF_MAXSTR)      :: COMP_NAME

! Locals

  type (MAPL_MetaComp), pointer   :: MAPL => null()
  type (ESMF_Config)              :: CF
  type (ESMF_State   )            :: INTERNAL

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    Iam = "Run1"
    call ESMF_GridCompGet( GC, name=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // Iam

! Get my internal MAPL_Generic state
!-----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)

! Start Total timer
!------------------

   call MAPL_TimerOn(MAPL,"TOTAL")
   call MAPL_TimerOn(MAPL,"RUN1" )

! Get parameters from generic state.
!-----------------------------------

    call MAPL_Get(MAPL,                          &
         INTERNAL_ESMF_STATE = INTERNAL,         &
                                       RC=STATUS )
    VERIFY_(STATUS)

!  All done
!-----------

   call MAPL_TimerOff(MAPL,"RUN1" )
   call MAPL_TimerOff(MAPL,"TOTAL")

   RETURN_(ESMF_SUCCESS)


 end subroutine RUN1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!BOP
! !IROUTINE: RUN2

! !INTERFACE:

subroutine RUN2 ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

  type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
  type(ESMF_State),    intent(inout) :: IMPORT ! Import state
  type(ESMF_State),    intent(inout) :: EXPORT ! Export state
  type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
  integer, optional,   intent(  out) :: RC     ! Error code:

! !DESCRIPTION: ??

!EOP


! ErrLog Variables

  character(len=ESMF_MAXSTR)          :: IAm
  integer                             :: STATUS
  character(len=ESMF_MAXSTR)          :: COMP_NAME

! Locals

  type (MAPL_MetaComp), pointer       :: MAPL => null()
  type (ESMF_State       )            :: INTERNAL
  type (MAPL_SunOrbit)                :: ORBIT
  type (ESMF_Config      )            :: CF

  real, pointer, dimension(:)         :: LATS => null()
  real, pointer, dimension(:)         :: LONS => null()
  real, pointer, dimension(:)         :: AREA => null()     ! needed to calculate TILEAREA in SaltWaterCore

! Begin... 

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    Iam = "Run2"
    call ESMF_GridCompGet( GC, name=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // Iam

! Get my internal MAPL_Generic state
!-----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)

! Start Total timer
!------------------

   call MAPL_TimerOn(MAPL,"TOTAL")
   call MAPL_TimerOn(MAPL,"RUN2" )

! Get parameters from generic state.
!-----------------------------------

    call MAPL_Get(MAPL,             &
         TILELATS  = LATS ,                      &
         TILELONS  = LONS ,                      &
         TILEAREA  = AREA ,                      &
         ORBIT     = ORBIT,                      &
         INTERNAL_ESMF_STATE = INTERNAL,         &
         CF = CF,                                &
                                       RC=STATUS )
    VERIFY_(STATUS)

! Update the skin variables each step
!------------------------------------

    call OBIOCORE(NT=size(LONS), RC=STATUS )
    VERIFY_(STATUS)

!  All done
!-----------

   call MAPL_TimerOff(MAPL,"RUN2" )
   call MAPL_TimerOff(MAPL,"TOTAL")

   RETURN_(ESMF_SUCCESS)

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine OBIOCORE(NT,RC)
   integer,           intent(IN ) :: NT
   integer, optional, intent(OUT) :: RC

!  Locals

   character(len=ESMF_MAXSTR)     :: IAm
   integer                        :: STATUS

! pointers to export
   real, pointer, dimension(:  )  :: CO2SCEX     => null()
   real, pointer, dimension(:,:)  :: DUDPEX      => null()
   real, pointer, dimension(:,:)  :: DUWTEX      => null()
   real, pointer, dimension(:,:)  :: DUSDEX      => null()
   real, pointer, dimension(:,:)  :: BCDPEX      => null()
   real, pointer, dimension(:,:)  :: BCWTEX      => null()
   real, pointer, dimension(:,:)  :: OCDPEX      => null()
   real, pointer, dimension(:,:)  :: OCWTEX      => null()
   real, pointer, dimension(:,:)  :: FSWBANDEX   => null()
   real, pointer, dimension(:,:)  :: FSWBANDNAEX => null()

! pointers to import
   real, pointer, dimension(:)    :: CO2SC     => null()
   real, pointer, dimension(:,:)  :: DUDP      => null()
   real, pointer, dimension(:,:)  :: DUWT      => null()
   real, pointer, dimension(:,:)  :: DUSD      => null()
   real, pointer, dimension(:,:)  :: BCDP      => null()
   real, pointer, dimension(:,:)  :: BCWT      => null()
   real, pointer, dimension(:,:)  :: OCDP      => null()
   real, pointer, dimension(:,:)  :: OCWT      => null()
   real, pointer, dimension(:,:)  :: FSWBAND   => null()
   real, pointer, dimension(:,:)  :: FSWBANDNA => null()

!  Begin...
!----------

   IAm =  trim(COMP_NAME) // "OBIOCORE"

! Pointers to inputs
!-------------------

    call MAPL_GetPointer(IMPORT,CO2SC  , 'CO2SC'  ,    RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT,DUDP   , 'DUDP'   ,    RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT,DUWT   , 'DUWT'   ,    RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT,DUSD   , 'DUSD'   ,    RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT,BCDP   , 'BCDP'   ,    RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT,BCWT   , 'BCWT'   ,    RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT,OCDP   , 'OCDP'   ,    RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT,OCWT   , 'OCWT'   ,    RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT,FSWBAND ,'FSWBAND' ,   RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT,FSWBANDNA,'FSWBANDNA', RC=STATUS); VERIFY_(STATUS)

! Pointers to outputs
!--------------------

    call MAPL_GetPointer(EXPORT,CO2SCEX,    'CO2SC'   ,    RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT,DUDPEX ,    'DUDP'    ,    RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT,DUWTEX ,    'DUWT'    ,    RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT,DUSDEX ,    'DUSD'    ,    RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT,BCDPEX ,    'BCDP'    ,    RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT,BCWTEX ,    'BCWT'    ,    RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT,OCDPEX ,    'OCDP'    ,    RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT,OCWTEX ,    'OCWT'    ,    RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT,FSWBANDEX,  'FSWBAND',     RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT,FSWBANDNAEX,'FSWBANDNA',   RC=STATUS); VERIFY_(STATUS)


    if  (  associated(CO2SCEX)      )  CO2SCEX      =  CO2SC
    if  (  associated(DUDPEX)       )  DUDPEX       =  DUDP
    if  (  associated(DUWTEX)       )  DUWTEX       =  DUWT
    if  (  associated(DUSDEX)       )  DUSDEX       =  DUSD
    if  (  associated(BCDPEX)       )  BCDPEX       =  BCDP
    if  (  associated(BCWTEX)       )  BCWTEX       =  BCWT
    if  (  associated(OCDPEX)       )  OCDPEX       =  OCDP
    if  (  associated(OCWTEX)       )  OCWTEX       =  OCWT
    if  (  associated(FSWBANDEX)    )  FSWBANDEX    =  FSWBAND
    if  (  associated(FSWBANDNAEX)  )  FSWBANDNAEX  =  FSWBANDNA

!  All done
!-----------

    RETURN_(ESMF_SUCCESS)

  end subroutine OBIOCORE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end subroutine RUN2

end module GEOS_ObioGridCompMod

