!  $Id$

#include "MAPL_Generic.h"

!=============================================================================
!BOP

! !MODULE: GEOS_DataSeaIce -- A {\it fake} seaice model

! !INTERFACE:

module GEOS_DataSeaIceGridCompMod

! !USES:

  use ESMF
  use MAPL

  use ice_state,          only: nt_Tsfc, nt_iage, nt_volpn, init_trcr_depend
  use ice_prescribed_mod

  implicit none
  private


! !PUBLIC MEMBER FUNCTIONS:

  public SetServices

  integer, parameter :: NUM_3D_ICE_TRACERS=3
  integer, parameter :: NUM_SNOW_LAYERS=1
  integer            :: NUM_ICE_CATEGORIES
  integer            :: NUM_ICE_LAYERS
  integer            :: NUM_ICE_LAYERS_ALL
  integer            :: NUM_SNOW_LAYERS_ALL
  integer            :: DO_CICE_THERMO

  logical            :: ocean_extData
  logical            :: seaIceT_extData

! !DESCRIPTION:
!
!   {\tt GEOS\_DataSeaIce} is a gridded component that reads the
!   ocean\_bcs file
!   This module interpolates the sea ice fraction data from
!   either daily or monthly values to the correct time of the simulation.
!   Data are read only if the simulation time is not in the save interval.
!

!EOP

   contains

!BOP

! !IROUTINE: SetServices -- Sets ESMF services for this component

! !INTERFACE:

  subroutine SetServices ( GC, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    integer, optional                  :: RC  ! return code

!  !DESCRIPTION: This version uses the MAPL\_GenericSetServices. This function sets
!                the Initialize and Finalize services, as well as allocating
!   our instance of a generic state and putting it in the
!   gridded component (GC). Here we only need to set the run method and
!   add the state variable specifications (also generic) to our instance
!   of the generic state. This is the way our true state variables get into
!   the ESMF\_State INTERNAL, which is in the MAPL\_MetaComp.
!
!EOP

!=============================================================================
!
! ErrLog Variables


    character(len=ESMF_MAXSTR)              :: IAm
    integer                                 :: STATUS
    character(len=ESMF_MAXSTR)              :: COMP_NAME

! Local derived type aliases

    type (MAPL_MetaComp    ), pointer   :: MAPL => null()
    type (ESMF_Config)                  :: CF

!=============================================================================

! Begin...

! Get my name and set-up traceback handle
! ---------------------------------------

    Iam = "SetServices"
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, CONFIG=CF, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // Iam

! Get my internal MAPL_Generic state
!-----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)

! Get constants from CF
! ---------------------

    call MAPL_GetResource ( MAPL,    DO_CICE_THERMO,     Label="USE_CICE_Thermo:" , DEFAULT=0, RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_GetResource ( MAPL,    ocean_extData, Label="OCEAN_EXT_DATA:",   DEFAULT=.FALSE., __RC__ ) ! .TRUE. or .FALSE.
    call MAPL_GetResource ( MAPL,    seaIceT_extData, Label="SEAICE_THICKNESS_EXT_DATA:",  DEFAULT=.FALSE., _RC ) ! .TRUE. or .FALSE.
    if (seaIceT_extData) then
      if(MAPL_AM_I_ROOT()) print *, 'You have set SEAICE_THICKNESS_EXT_DATA: .TRUE. For this feature, a _newer_: version-2 of data sea ice must be used.'
      _FAIL('Recompile and try.')
    endif

    cice_init_: if (DO_CICE_THERMO /= 0) then
       if(MAPL_AM_I_ROOT()) print *, 'Using Data Sea Ice GC to do CICE Thermo in AMIP mode'
       call ESMF_ConfigGetAttribute(CF, NUM_ICE_CATEGORIES, Label="CICE_N_ICE_CATEGORIES:" , RC=STATUS)
       VERIFY_(STATUS)
       call ESMF_ConfigGetAttribute(CF, NUM_ICE_LAYERS,     Label="CICE_N_ICE_LAYERS:" ,     RC=STATUS)
       VERIFY_(STATUS)
    else
       NUM_ICE_CATEGORIES = 1
       NUM_ICE_LAYERS     = 1
    end if cice_init_

    NUM_ICE_LAYERS_ALL=NUM_ICE_LAYERS*NUM_ICE_CATEGORIES
    NUM_SNOW_LAYERS_ALL=NUM_SNOW_LAYERS*NUM_ICE_CATEGORIES

! Set the Run entry point
! -----------------------

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,  Run, RC=STATUS)
    VERIFY_(STATUS)


! Set the state variable specs.
! -----------------------------

!BOS

! !Import state:

  if (ocean_extData) then
    call MAPL_AddImportSpec(GC,                  &
      SHORT_NAME         = 'DATA_ICE',           &
      LONG_NAME          = 'sea_ice_concentration',        &
      UNITS              = '1',                  &
      DIMS               = MAPL_DimsHorzOnly,    &
      VLOCATION          = MAPL_VLocationNone,   &
      RC=STATUS)
    VERIFY_(status)
  endif

  call MAPL_AddImportSpec(GC,                                 &
    SHORT_NAME         = 'HI',                                &
    LONG_NAME          = 'seaice_skin_layer_depth',           &
    UNITS              = 'm',                                 &
    DIMS               = MAPL_DimsHorzOnly,                   &
    VLOCATION          = MAPL_VLocationNone,                  &
                                                   RC=STATUS  )
  VERIFY_(STATUS)

  if (DO_CICE_THERMO == 0) then
     call MAPL_AddImportSpec(GC,                                 &
          SHORT_NAME         = 'TI',                                &
          LONG_NAME          = 'seaice_skin_temperature',           &
          UNITS              = 'K',                                 &
          DIMS               = MAPL_DimsHorzOnly,                   &
          VLOCATION          = MAPL_VLocationNone,                  &
          RC=STATUS  )
     VERIFY_(STATUS)
  end if

  call MAPL_AddImportSpec(GC,                                 &
    SHORT_NAME         = 'SI',                                &
    LONG_NAME          = 'seaice_skin_salinity',              &
    UNITS              = 'psu',                               &
    DIMS               = MAPL_DimsHorzOnly,                   &
    VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
   VERIFY_(STATUS)

  if (DO_CICE_THERMO /= 0) then
     call MAPL_AddImportSpec(GC,                               &
          SHORT_NAME         = 'FRACICE',                           &
          LONG_NAME          = 'fractional_cover_of_seaice',        &
          UNITS              = '1',                                 &
          DIMS               = MAPL_DimsHorzOnly,                   &
          UNGRIDDED_DIMS     = (/NUM_ICE_CATEGORIES/),              &
          VLOCATION          = MAPL_VLocationNone,                  &
          RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                                 &
          SHORT_NAME         = 'TI',                                &
          LONG_NAME          = 'seaice_skin_temperature',           &
          UNITS              = 'K',                                 &
          DIMS               = MAPL_DimsHorzOnly,                   &
          UNGRIDDED_DIMS     = (/NUM_ICE_CATEGORIES/),              &
          VLOCATION          = MAPL_VLocationNone,                  &
          RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                                 &
          SHORT_NAME         = 'VOLICE',                            &
          LONG_NAME          = 'ice_category_volume_per_unit_area_of_grid_cell',&
          UNITS              = 'm',                                 &
          DIMS               = MAPL_DimsHorzOnly,                   &
          UNGRIDDED_DIMS     = (/NUM_ICE_CATEGORIES/),              &
          VLOCATION          = MAPL_VLocationNone,                  &
          RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                                 &
          SHORT_NAME         = 'VOLSNO',                            &
          LONG_NAME          = 'sno_category_volume_per_unit_area_of_grid_cell',&
          UNITS              = 'm',                                 &
          DIMS               = MAPL_DimsHorzOnly,                   &
          UNGRIDDED_DIMS     = (/NUM_ICE_CATEGORIES/),              &
          VLOCATION          = MAPL_VLocationNone,                  &
          RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                                 &
          SHORT_NAME         = 'ERGICE',                            &
          LONG_NAME          = 'ice_category_layer_internal_energy',&
          UNITS              = 'J m-2',                             &
          DIMS               = MAPL_DimsHorzOnly,                   &
          VLOCATION          = MAPL_VLocationNone,                  &
          UNGRIDDED_DIMS     = (/NUM_ICE_LAYERS_ALL/),              &
          RESTART            = MAPL_RestartSkip,                    &
          RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                                &
          SHORT_NAME         = 'ERGSNO',                            &
          LONG_NAME          = 'snow_category_layer_internal_energy',&
          UNITS              = 'J m-2',                             &
          DIMS               = MAPL_DimsHorzOnly,                   &
          VLOCATION          = MAPL_VLocationNone,                  &
          UNGRIDDED_DIMS     = (/NUM_SNOW_LAYERS_ALL/),             &
          RESTART            = MAPL_RestartSkip,                    &
          RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC                                ,&
          LONG_NAME          = 'melt_pond_volume'                  ,&
          UNITS              = 'm'                                 ,&
          SHORT_NAME         = 'MPOND'                             ,&
          UNGRIDDED_DIMS     = (/NUM_ICE_CATEGORIES/),              &
          DIMS               = MAPL_DimsHorzOnly,                   &
          VLOCATION          = MAPL_VLocationNone,                  &
          RESTART            = MAPL_RestartSkip,                    &
          RC=STATUS                                                 )

     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                                &
          SHORT_NAME         = 'TAUAGE',                            &
          LONG_NAME          = 'volume_weighted_mean_ice_age',      &
          UNITS              = 's',                                 &
          UNGRIDDED_DIMS     = (/NUM_ICE_CATEGORIES/),              &
          DIMS               = MAPL_DimsHorzOnly,                   &
          VLOCATION          = MAPL_VLocationNone,                  &
          RESTART            = MAPL_RestartSkip,                    &
          RC=STATUS  )
     VERIFY_(STATUS)
  end if ! (DO_CICE_THERMO /= 0)

!  !Export state:

  call MAPL_AddExportSpec(GC,                                 &
    SHORT_NAME         = 'UI',                                &
    LONG_NAME          = 'zonal_velocity_of_surface_seaice',  &
    UNITS              = 'm s-1 ',                            &
    DIMS               = MAPL_DimsHorzOnly,                   &
    VLOCATION          = MAPL_VLocationNone,                  &
                                                   RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                                 &
    SHORT_NAME         = 'VI',                                &
    LONG_NAME          = 'meridional_velocity_of_surface_seaice',&
    UNITS              = 'm s-1 ',                            &
    DIMS               = MAPL_DimsHorzOnly,                   &
    VLOCATION          = MAPL_VLocationNone,                  &
                                                   RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                                 &
       SHORT_NAME         = 'FRACICE',                           &
       LONG_NAME          = 'fractional_cover_of_seaice',        &
       UNITS              = '1',                                 &
       DIMS               = MAPL_DimsHorzOnly,                   &
       VLOCATION          = MAPL_VLocationNone,                  &
       RC=STATUS  )
  VERIFY_(STATUS)

  if (DO_CICE_THERMO /= 0) then
     call MAPL_AddExportSpec(GC,                                     &
          SHORT_NAME         = 'TAUYBOT',                           &
          LONG_NAME          = 'northward_stress_at_base_of_ice',   &
          UNITS              = 'N m-2',                             &
          DIMS               = MAPL_DimsHorzOnly,                   &
          VLOCATION          = MAPL_VLocationNone,                  &
          RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                                 &
          SHORT_NAME         = 'HICE',                              &
          LONG_NAME          = 'mean_ice_thickness_of_grid_cell',   &
          UNITS              = 'm',                                 &
          DIMS               = MAPL_DimsHorzOnly,                   &
          VLOCATION          = MAPL_VLocationNone,                  &
          RC=STATUS  )
     VERIFY_(STATUS)


     call MAPL_AddExportSpec(GC,                                 &
          SHORT_NAME         = 'HSNO',                              &
          LONG_NAME          = 'mean_snow_thickness_of_grid_cell',  &
          UNITS              = 'm',                                 &
          DIMS               = MAPL_DimsHorzOnly,                   &
          VLOCATION          = MAPL_VLocationNone,                  &
          RC=STATUS  )
     VERIFY_(STATUS)
  end if
!EOS

  call MAPL_TimerAdd(GC,    name="RUN"     ,RC=STATUS)
  VERIFY_(STATUS)
  call MAPL_TimerAdd(GC,    name="-UPDATE" ,RC=STATUS)
  VERIFY_(STATUS)

! Set generic init and final methods
! ----------------------------------

  call MAPL_GenericSetServices    ( GC, RC=STATUS)
  VERIFY_(STATUS)

  RETURN_(ESMF_SUCCESS)

  end subroutine SetServices



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!BOP

! !IROUTINE: RUN -- Run stage for the DataSeaIce component

! !INTERFACE:

subroutine RUN ( GC, IMPORT, EXPORT, CLOCK, RC )


! !ARGUMENTS:

  type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component
  type(ESMF_State),    intent(inout) :: IMPORT ! Import state
  type(ESMF_State),    intent(inout) :: EXPORT ! Export state
  type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
  integer, optional,   intent(  out) :: RC     ! Error code:

! !DESCRIPTION: Periodically refreshes the SST and Ice information.

!EOP


! ErrLog Variables

  character(len=ESMF_MAXSTR)          :: IAm
  integer                             :: STATUS
  character(len=ESMF_MAXSTR)          :: COMP_NAME

! Locals

  type (MAPL_MetaComp),  pointer      :: MAPL  => null()
  logical                             :: FRIENDLY
  type(ESMF_FIELD)                    :: FIELD
  type (ESMF_Time)                    :: CurrentTime
  character(len=ESMF_MAXSTR)          :: DataFrtFile
  integer                             :: IFCST
  logical                             :: FCST
  real                                :: TAU_SIT
  real                                :: DT
  real                                :: RUN_DT
  real                                :: CTB  ! Ocean-ice turbulent mixing coefficient (m/sec)
  real, parameter                     :: CW=MAPL_CAPWTR
  real                                :: TICE
  integer                             :: I, J, IM, JM, N

! below are for CICE Thermo
  real(kind=ESMF_KIND_R8), allocatable  :: FRCICE(:,:)
  real(kind=ESMF_KIND_R8)               :: FRACICEDB(1)
  real(kind=ESMF_KIND_R8)               :: FRWATERDB(1)
  real(kind=ESMF_KIND_R8)               :: LATSDB(1)
  real(kind=ESMF_KIND_R8)               :: TFDB(1)
  real(kind=ESMF_KIND_R8)               :: FRDB(NUM_ICE_CATEGORIES)
  real(kind=ESMF_KIND_R8)               :: VOLICEDB(NUM_ICE_CATEGORIES)
  real(kind=ESMF_KIND_R8)               :: VOLSNODB(NUM_ICE_CATEGORIES)
  real(kind=ESMF_KIND_R8)               :: ERGICEDB(NUM_ICE_LAYERS_ALL)
  real(kind=ESMF_KIND_R8)               :: ERGSNODB(NUM_SNOW_LAYERS_ALL)

  real(kind=ESMF_KIND_R8), dimension(NUM_3D_ICE_TRACERS, NUM_ICE_CATEGORIES) :: TRACERSDB2
  real,                    dimension(NUM_3D_ICE_TRACERS,NUM_ICE_CATEGORIES)  :: TRACERS
  real                                  :: TNH, TSH
! above were for CICE Thermo

! pointers to export

   real, pointer, dimension(:,:)  :: UI   => null()
   real, pointer, dimension(:,:)  :: VI   => null()
   real, pointer, dimension(:,:)  :: FR   => null()

! pointers to import

   real, pointer, dimension(:,:)  :: TI   => null()
   real, pointer, dimension(:,:)  :: HI   => null()
   real, pointer, dimension(:,:)  :: SI   => null()

! below are for CICE Thermo
   real, pointer, dimension(:,:,:):: FR8     => null()
   real, pointer, dimension(:,:,:):: TI8     => null()

   real, pointer, dimension(:,:,:):: VOLICE  => null()
   real, pointer, dimension(:,:,:):: VOLSNO  => null()
   real, pointer, dimension(:,:,:):: TAUAGE  => null()
   real, pointer, dimension(:,:,:):: MPOND   => null()

   real, pointer, dimension(:,:,:):: ERGICE  => null()
   real, pointer, dimension(:,:,:):: ERGSNO  => null()

   real, pointer, dimension(:,:)  :: LATS    => null()
   real, pointer, dimension(:,:)  :: LONS    => null()

   real, allocatable, dimension(:,:)  :: FRT
   real :: f

   real, pointer :: DATA_ice(:,:) => null()
! above were for CICE Thermo


!  Begin...
!----------

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    Iam = "Run"
    call ESMF_GridCompGet( GC, name=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // Iam

! Get my internal MAPL_Generic state
!----------------------------------

    call MAPL_GetObjectFromGC(GC, MAPL, STATUS)
    VERIFY_(STATUS)

    call MAPL_Get(MAPL,                      &
         LATS  = LATS ,                      &
         LONS  = LONS ,                      &
                                RC=STATUS )
    VERIFY_(STATUS)

! Start Total timer
!------------------

   call MAPL_TimerOn(MAPL,"TOTAL")
   call MAPL_TimerOn(MAPL,"RUN" )


! Pointers to Imports
!--------------------

   if (ocean_extData) then
     call MAPL_GetPointer(IMPORT, DATA_ice     ,  'DATA_ICE', __RC__)
   endif

   if (DO_CICE_THERMO == 0) then
      call MAPL_GetPointer(IMPORT, TI    ,  'TI'   , RC=STATUS)
      VERIFY_(STATUS)
   else
      call MAPL_GetPointer(IMPORT, TI8   ,  'TI'   , RC=STATUS)
      VERIFY_(STATUS)
   end if
   call MAPL_GetPointer(IMPORT, HI      ,  'HI'   , RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT, SI      ,  'SI'   , RC=STATUS)
   VERIFY_(STATUS)

   if (DO_CICE_THERMO /= 0) then
     call MAPL_GetPointer(IMPORT, FR8     ,  'FRACICE', RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(IMPORT, VOLICE  ,  'VOLICE' , RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(IMPORT, ERGICE  ,  'ERGICE' , RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(IMPORT, VOLSNO  ,  'VOLSNO' , RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(IMPORT, ERGSNO  ,  'ERGSNO' , RC=STATUS)
     VERIFY_(STATUS)
   end if

!  Pointers to Exports
!---------------------

    call MAPL_GetPointer(EXPORT,      UI  , 'UI'       , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT,      VI  , 'VI'       , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT,      FR  , 'FRACICE'  , RC=STATUS)
    VERIFY_(STATUS)

! Set current time and calendar
!------------------------------

    call ESMF_ClockGet(CLOCK, currTime=CurrentTime, rc=STATUS)
    VERIFY_(STATUS)

   if (.not. ocean_extData) then
    ! Get the file name from the resource file
    !-----------------------------------------
    call MAPL_GetResource(MAPL,DataFrtFile,LABEL="DATA_FRT_FILE:", RC=STATUS)
    VERIFY_(STATUS)
  endif

! In atmospheric forecast mode we do not have future Sea Ice Conc
!---------------------------------------------------------------

    call MAPL_GetResource(MAPL,IFCST,LABEL="IS_FCST:",default=0,   RC=STATUS)
    VERIFY_(STATUS)

    FCST = IFCST==1

! Get relaxation time
!--------------------

    call MAPL_GetResource(MAPL,TAU_SIT, LABEL="TAU_SIT:", default=86400.0,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource(MAPL,RUN_DT , LABEL="RUN_DT:" ,                 RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource(MAPL,DT     , LABEL="DT:"    , default=RUN_DT, RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource(MAPL,CTB    , LABEL="CTB:"   , default=1.0e-4, RC=STATUS)
    VERIFY_(STATUS)

! If using CICE Thermodynamics in AMIP mode, get prescribed ice thickness
! FOR NOW, following way sets thickness to a constant value in either hemisphere-
! And in future (>04/2016) we will explore "other" ideas. [BZ/SA/MT]
!----------------------------------------------------------------------------------
    if (DO_CICE_THERMO /= 0) then
      call MAPL_GetResource (MAPL, TNH, Label="PRESCRIBED_ICE_NH:" , DEFAULT=1.0, RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetResource (MAPL, TSH, Label="PRESCRIBED_ICE_SH:" , DEFAULT=0.75, RC=STATUS)
      VERIFY_(STATUS)
    end if

!  Update the friendly skin values
!---------------------------------

   call MAPL_TimerOn(MAPL,"-UPDATE" )

   if (DO_CICE_THERMO /= 0) then
     call MAPL_Get(MAPL, IM=IM, JM=JM, RC=STATUS)
     VERIFY_(STATUS)
     allocate(FRT(IM,JM), FRCICE(IM,JM))
   end if

   if (DO_CICE_THERMO == 0) then
     if(associated(FR)) then
       if (ocean_extData) then
         FR = data_ice ! netcdf variable
       else ! binary
         call MAPL_ReadForcing(MAPL,'FRT',DataFrtFile, CURRENTTIME, FR, INIT_ONLY=FCST, __RC__)
       end if

       call MAPL_GetResource ( MAPL, STRICT_ICE_FRACTION, Label="STRICT_ICE_FRACTION:", DEFAULT=.FALSE., __RC__)
       if (STRICT_ICE_FRACTION) then
          if (any(FR < 0.0) .or. any(FR > 1.0)) then
             _ASSERT(.FALSE.,'Error in fraci file. Negative or larger-than-one fraction found')
          endif
       else
          call MAPL_GetResource ( MAPL, ICE_FRACTION_TOLERANCE, Label="ICE_FRACTION_TOLERANCE:", DEFAULT=1.0e-2, __RC__)
          ! Now we can use a tolerance to allow the code to run with ice fraction slight higher than 1.0
          where (FR < (0.0-ICE_FRACTION_TOLERANCE)) FR = 0.0
          where (FR > (1.0+ICE_FRACTION_TOLERANCE)) FR = 1.0
       end if

     end if
   else
     if (ocean_extData) then
       frt = data_ice ! netcdf variable
     else ! binary
       call MAPL_ReadForcing(MAPL,'FRT',DataFrtFile, CURRENTTIME, FRT, INIT_ONLY=FCST, __RC__)
     end if

! Sanity checks
     do I=1, size(FRT,1)
        do J=1, size(FRT,2)
           f=FRT(I,J)
           if (f==MAPL_UNDEF) cycle
           if ((f < 0.0) .or. (f > 1.0)) then
              print *, 'Error in fraci file. Negative or larger-than-one fraction found'
              _ASSERT(.FALSE.,'needs informative message')
           end if
        end do
     end do

     if(associated(FR)) FR = FRT
   end if ! (DO_CICE_THERMO == 0)

   if (DO_CICE_THERMO == 0) then

     where (TI /= MAPL_Undef)
       TI   =(TI + (DT/TAU_SIT)*MAPL_TICE)/(1.+ (DT/TAU_SIT))
     end where

     HI   = HI
   end if

   SI   = 30.0

   call MAPL_TimerOff(MAPL,"-UPDATE" )

!  Update the exports
!--------------------

   if(associated(UI)) UI = 0.0
   if(associated(VI)) VI = 0.0

! Clean-up
!---------

   if (DO_CICE_THERMO /= 0) then
     deallocate(FRT)
     deallocate(FRCICE)
   end if

!  All done
!-----------

   call MAPL_TimerOff(MAPL,"RUN"  )
   call MAPL_TimerOff(MAPL,"TOTAL")

   RETURN_(ESMF_SUCCESS)
end subroutine RUN

end module GEOS_DataSeaIceGridCompMod
