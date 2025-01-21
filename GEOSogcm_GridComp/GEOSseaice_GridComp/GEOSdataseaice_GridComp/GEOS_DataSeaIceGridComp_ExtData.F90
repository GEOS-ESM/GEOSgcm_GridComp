!  $Id$

#include "MAPL_Generic.h"

!=============================================================================
!BOP

! !MODULE: GEOS_DataSeaIce -- A {\it fake} ocean sea ice

! !INTERFACE:

module GEOS_DataSeaIceGridCompMod

! !USES:

  use ESMF
  use MAPL

  implicit none
  private

! !PUBLIC MEMBER FUNCTIONS:

  public SetServices

  integer            :: DO_CICE_THERMO
  logical            :: ocean_extData
  logical            :: seaIceT_extData

  real               :: MAX_SEAICE_THICKNESS
  real               :: MIN_SEAICE_THICKNESS

! !DESCRIPTION:
!
!   {\tt GEOS\_DataSeaIce} is a gridded component that reads the
!   ocean\_bcs file
!   This module interpolates the sea ice fraction and optionally thickness data from
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
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, CONFIG=CF,  _RC )
    Iam = trim(COMP_NAME) // Iam

! Set the Run entry point
! -----------------------

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,  Run,  _RC )

    call MAPL_GetObjectFromGC ( GC, MAPL,  _RC )

    call MAPL_GetResource ( MAPL,    DO_CICE_THERMO,     Label="USE_CICE_Thermo:" , DEFAULT=0,  _RC )

    if (DO_CICE_THERMO /=0) then
      _ASSERT(.FALSE.,'Must not use CICE_Thermo in (Ext) data sea ice. Fix and try.')
    endif

    call MAPL_GetResource ( MAPL,    ocean_extData, Label="OCEAN_EXT_DATA:",    DEFAULT=.FALSE., _RC ) ! .TRUE. or .FALSE.

    ! This feature (ice thickness data) will be with ExtData ONLY.
    call MAPL_GetResource ( MAPL,    seaIceT_extData, Label="SEAICE_THICKNESS_EXT_DATA:",  DEFAULT=.FALSE., _RC ) ! .TRUE. or .FALSE.

    call MAPL_GetResource ( MAPL, MIN_SEAICE_THICKNESS, Label='MIN_SEAICE_THICKNESS:',  DEFAULT=0.0, _RC ) ! .TRUE. or .FALSE.
    call MAPL_GetResource ( MAPL, MAX_SEAICE_THICKNESS, Label='MAX_SEAICE_THICKNESS:',  DEFAULT=2.0,    _RC ) ! .TRUE. or .FALSE.

    if (MAPL_AM_I_ROOT()) then
        print *, '*** WIP: Using V2 of Data Sea Ice GC.'
        print *, '*** WIP: SEAICE_THICKNESS_EXT_DATA: ', seaIceT_extData
        print *, '*** WIP: OCEAN_EXT_DATA:            ', ocean_extData
    end if


! Set the state variable specs.
! -----------------------------

!BOS

! !Import state:

  if (ocean_extData) then
    call MAPL_AddImportSpec(GC,                  &
      SHORT_NAME         = 'DATA_ICEC',          &
      LONG_NAME          = 'sea_ice_concentration',        &
      UNITS              = '1',                  &
      DIMS               = MAPL_DimsHorzOnly,    &
      VLOCATION          = MAPL_VLocationNone,   &
       _RC )
  endif

  if (seaIceT_extData) then
    call MAPL_AddImportSpec(GC,                  &
      SHORT_NAME         = 'DATA_SIT',           &
      LONG_NAME          = 'sea_ice_thickness',  &
      UNITS              = 'm',                  &
      DIMS               = MAPL_DimsHorzOnly,    &
      VLOCATION          = MAPL_VLocationNone,   &
      _RC )
  else
    call MAPL_AddImportSpec(GC,                         &
      SHORT_NAME         = 'HI',                        &
      LONG_NAME          = 'seaice_skin_layer_depth',   &
      UNITS              = 'm',                         &
      DIMS               = MAPL_DimsHorzOnly,           &
      VLOCATION          = MAPL_VLocationNone,          &
      _RC )

    call MAPL_AddImportSpec(GC,                         &
      SHORT_NAME         = 'TI',                        &
      LONG_NAME          = 'seaice_skin_temperature',   &
      UNITS              = 'K',                         &
      DIMS               = MAPL_DimsHorzOnly,           &
      VLOCATION          = MAPL_VLocationNone,          &
      _RC )

    call MAPL_AddImportSpec(GC,                         &
      SHORT_NAME         = 'SI',                        &
      LONG_NAME          = 'seaice_skin_salinity',      &
      UNITS              = 'psu',                       &
      DIMS               = MAPL_DimsHorzOnly,           &
      VLOCATION          = MAPL_VLocationNone,          &
      _RC )
  endif

!  !Export state:

  call MAPL_AddExportSpec(GC,                                 &
    SHORT_NAME         = 'UI',                                &
    LONG_NAME          = 'zonal_velocity_of_surface_seaice',  &
    UNITS              = 'm s-1 ',                            &
    DIMS               = MAPL_DimsHorzOnly,                   &
    VLOCATION          = MAPL_VLocationNone,                  &
    _RC )

  call MAPL_AddExportSpec(GC,                                 &
    SHORT_NAME         = 'VI',                                &
    LONG_NAME          = 'meridional_velocity_of_surface_seaice',&
    UNITS              = 'm s-1 ',                            &
    DIMS               = MAPL_DimsHorzOnly,                   &
    VLOCATION          = MAPL_VLocationNone,                  &
    _RC )

  call MAPL_AddExportSpec(GC,                                 &
    SHORT_NAME         = 'FRACICE',                        &
    LONG_NAME          = 'fractional_cover_of_seaice',     &
    UNITS              = '1',                              &
    DIMS               = MAPL_DimsHorzOnly,                &
    VLOCATION          = MAPL_VLocationNone,               &
    _RC )

  if (seaIceT_extData) then
    call MAPL_AddExportSpec(GC,                            &
      SHORT_NAME         = 'SEAICETHICKNESS',              &
      LONG_NAME          = 'seaice_thickness',             &
      UNITS              = 'm',                            &
      DIMS               = MAPL_DimsHorzOnly,              &
      VLOCATION          = MAPL_VLocationNone,             &
      _RC )
  endif

!EOS

  call MAPL_TimerAdd(GC,    name="RUN"     ,  _RC )
  call MAPL_TimerAdd(GC,    name="-UPDATE" ,  _RC )

! Set generic init and final methods
! ----------------------------------

  call MAPL_GenericSetServices    ( GC,  _RC )

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

! !DESCRIPTION: Periodically refreshes the ocean bcs information.

!EOP

! ErrLog Variables

  character(len=ESMF_MAXSTR)          :: IAm
  integer                             :: STATUS
  character(len=ESMF_MAXSTR)          :: COMP_NAME

! Locals

  type (MAPL_MetaComp),  pointer      :: MAPL  => null()
  type (ESMF_Time)                    :: CurrentTime
  character(len=ESMF_MAXSTR)          :: DataFrtFile
  integer                             :: IFCST
  logical                             :: FCST
  real                                :: TAU_SIT
  real                                :: DT
  real                                :: RUN_DT

! pointers to export

   real, pointer, dimension(:,:)  :: UI   => null()
   real, pointer, dimension(:,:)  :: VI   => null()
   real, pointer, dimension(:,:)  :: FR   => null()
   real, pointer, dimension(:,:)  :: SIT  => null()

! pointers to import

   real, pointer, dimension(:,:)  :: TI   => null()
   real, pointer, dimension(:,:)  :: HI   => null()
   real, pointer, dimension(:,:)  :: SI   => null()
   real, pointer, dimension(:,:)  :: DATA_icec => null()
   real, pointer, dimension(:,:)  :: DATA_sit  => null()

!  Begin...
!----------

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    Iam = "Run"
    call ESMF_GridCompGet( GC, name=COMP_NAME,  _RC )
    Iam = trim(COMP_NAME) // Iam

! Get my internal MAPL_Generic state
!----------------------------------

    call MAPL_GetObjectFromGC(GC, MAPL,  _RC )

! Start Total timer
!------------------

   call MAPL_TimerOn(MAPL,"TOTAL")
   call MAPL_TimerOn(MAPL,"RUN" )

! Pointers to Imports
!--------------------

   if (ocean_extData) then
     call MAPL_GetPointer(IMPORT, DATA_icec    ,  'DATA_ICEC', _RC )
   endif

   if (seaIceT_extData) then
     call MAPL_GetPointer(IMPORT, DATA_sit     ,  'DATA_SIT', _RC )
   else
     call MAPL_GetPointer(IMPORT, TI    ,  'TI'  , _RC )
     call MAPL_GetPointer(IMPORT, HI    ,  'HI'  , _RC )
     call MAPL_GetPointer(IMPORT, SI    ,  'SI'  , _RC )
   endif

!  Pointers to Exports
!---------------------

    call MAPL_GetPointer(EXPORT, UI  , 'UI'       ,  _RC )
    call MAPL_GetPointer(EXPORT, VI  , 'VI'       ,  _RC )
    call MAPL_GetPointer(EXPORT, FR  , 'FRACICE'  ,  _RC )

! Set current time and calendar
!------------------------------

    call ESMF_ClockGet(CLOCK, currTime=CurrentTime,  _RC )

   if (.not. ocean_extData) then
    ! Get the file name from the resource file
    !-----------------------------------------
    call MAPL_GetResource(MAPL,DataFrtFile,LABEL="DATA_FRT_FILE:",  _RC )
  endif

! In atmospheric forecast mode we do not have future Sea Ice Conc
!----------------------------------------------------------------

    call MAPL_GetResource(MAPL,IFCST,LABEL="OGCM_IS_FCST:",default=0, _RC )
    FCST = IFCST==1

! Get relaxation time
!--------------------

   if (.not. seaIceT_extData) then
    call MAPL_GetResource(MAPL,TAU_SIT, LABEL="TAU_SIT:", default=86400.0, _RC )
    call MAPL_GetResource(MAPL,RUN_DT , LABEL="RUN_DT:" ,                  _RC )
    call MAPL_GetResource(MAPL,DT     , LABEL="DT:"    , default=RUN_DT,   _RC )
   endif

!  Update data
!-------------

   call MAPL_TimerOn(MAPL,"-UPDATE" )

   if(associated(FR)) then
     if (ocean_extData) then
       FR = DATA_icec ! netcdf variable
     else ! binary
       call MAPL_ReadForcing(MAPL,'FRT',DataFrtFile, CURRENTTIME, FR, INIT_ONLY=FCST, _RC)
     end if

!    Sanity check
     if (any(FR < 0.0) .or. any(FR > 1.0)) then
        if(MAPL_AM_I_ROOT()) print *, 'Error in sea fraci file: negative or larger-than-one fraction found.'
        _ASSERT(.FALSE.,'Fix sea fraci and try.')
     endif
   end if

   if (seaIceT_extData) then
     call MAPL_GetPointer(EXPORT, SIT , 'SEAICETHICKNESS'  ,  _RC )
     SIT = 0.0

     if (associated(FR)) then
         where ((FR > tiny(FR)) .and. (Data_sit /= MAPL_UNDEF))
             SIT = min(max(MIN_SEAICE_THICKNESS, Data_sit), MAX_SEAICE_THICKNESS)  ! SIT <- DATA_sit | FR * DATA_sit ?
         end where
     else
         where (Data_sit /= MAPL_UNDEF)
             SIT = min(max(MIN_SEAICE_THICKNESS, Data_sit), MAX_SEAICE_THICKNESS)  ! SIT <- DATA_sit | FR * DATA_sit ?
         end where
     end if
   else
     where (TI /= MAPL_Undef)
       TI = (TI + (DT/TAU_SIT)*MAPL_TICE)/(1.+ (DT/TAU_SIT))
     end where
     HI = HI
     SI = 30.0
   end if

   call MAPL_TimerOff(MAPL,"-UPDATE" )

!  Update the exports
!--------------------

   if(associated(UI)) UI = 0.0
   if(associated(VI)) VI = 0.0

!  All done
!-----------

   call MAPL_TimerOff(MAPL,"RUN"  )
   call MAPL_TimerOff(MAPL,"TOTAL")

   RETURN_(ESMF_SUCCESS)
end subroutine RUN

end module GEOS_DataSeaIceGridCompMod
