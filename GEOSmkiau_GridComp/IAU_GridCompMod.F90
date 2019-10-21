!  $Id$

#include "MAPL_Generic.h"

!=============================================================================
!BOP

! !MODULE: IAU -- A Module to handle IAU increments

! !INTERFACE:

module IAU_GridCompMod

! !USES:

  use ESMF
  use MAPL_Mod
  
  implicit none
  private

! !PUBLIC MEMBER FUNCTIONS:

  public SetServices

!=============================================================================

! !DESCRIPTION:
! 
!

!EOP

contains

!BOP

! ! IROUTINE: SetServices -- Sets ESMF services for this component

! ! INTERFACE:

  subroutine SetServices ( GC, RC )

! ! ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    integer, optional                  :: RC  ! return code

! ! DESCRIPTION: This version uses the MAPL_GenericSetServices. This function sets
!                the Initialize and Finalize services, as well as allocating
!   our instance of a generic state and putting it in the 
!   gridded component (GC). Here we only need to set the run method and
!   add the state variable specifications (also generic) to our instance
!   of the generic state. This is the way our true state variables get into
!   the ESMF_State INTERNAL, which is in the MAPL_MetaComp.

!EOP

!=============================================================================
!
! ErrLog Variables

    character(len=ESMF_MAXSTR)              :: IAm
    integer                                 :: STATUS
    character(len=ESMF_MAXSTR)              :: COMP_NAME
    type (MAPL_MetaComp),         pointer   :: MAPL

    integer MKIAU

!=============================================================================

! Begin...

! Get my name and set-up traceback handle
! ---------------------------------------

    Iam = 'SetServices'
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // Iam

! Retrieve the pointer to the state
!----------------------------------

   call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
   VERIFY_(STATUS)

!  MKIAU = MAPL_AddChild(GC, NAME='AIAU', SS=MKIAUSetServices, RC=STATUS)
!  VERIFY_(STATUS)

! Set the Run entry point
! -----------------------

    call MAPL_GridCompSetEntryPoint ( gc, ESMF_METHOD_RUN,  Run1,  &
                                      RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GridCompSetEntryPoint ( gc, ESMF_METHOD_RUN,  Run2,  &
                                      RC=STATUS)
    VERIFY_(STATUS)

! Set the state variable specs.
! -----------------------------

! !IMPORT STATE:

    call MAPL_AddImportSpec ( gc,                                  &
         SHORT_NAME = 'DUDT',                                      &
         LONG_NAME  = 'eastward_wind_analysis_increment',          &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                                  &
         SHORT_NAME = 'DVDT',                                      &
         LONG_NAME  = 'northward_wind_analysis_increment',         &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                                  &
         SHORT_NAME = 'DTDT',                                      &
         LONG_NAME  = 'temperature_analysis_increment',            &
         UNITS      = 'K',                                         &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                                  &
         SHORT_NAME = 'DPEDT',                                     &
         LONG_NAME  = 'edge_pressure_analysis_increment',          &
         UNITS      = 'Pa',                                        &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationEdge,               RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                                  &
         SHORT_NAME = 'DQVDT',                                     &
         LONG_NAME  = 'specific_humidity_analysis_increment',      &
         UNITS      = 'kg kg-1',                                   &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                                  &
         SHORT_NAME = 'DO3DT',                                     &
         LONG_NAME  = 'ozone_analysis_increment',                  &
         UNITS      = 'ppmv',                                      &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                                  &
         SHORT_NAME = 'DTSDT',                                     &
         LONG_NAME  = 'skin_temparature_increment',                &
         UNITS      = 'K',                                         &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
    VERIFY_(STATUS)
    

#ifdef _MAYBE_
    call MAPL_AddImportSpec(GC,                                  &
         SHORT_NAME = 'AGCM_Exports',                            &
         LONG_NAME  = 'export_iau_increments',                   &
         UNITS      = 'X',                                       &
         DATATYPE   = MAPL_BundleItem,                           &
         RC=STATUS  )
    VERIFY_(STATUS)
#endif /* _MAYBE_ */

! !INTERNAL STATE:

    call MAPL_AddInternalSpec ( gc,                                  &
         SHORT_NAME = 'DUDT',                                      &
         LONG_NAME  = 'eastward_wind_analysis_increment',          &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec ( gc,                                  &
         SHORT_NAME = 'DVDT',                                      &
         LONG_NAME  = 'northward_wind_analysis_increment',         &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec ( gc,                                  &
         SHORT_NAME = 'DTDT',                                      &
         LONG_NAME  = 'temperature_analysis_increment',            &
         UNITS      = 'K',                                         &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec ( gc,                                  &
         SHORT_NAME = 'DPEDT',                                     &
         LONG_NAME  = 'edge_pressure_analysis_increment',          &
         UNITS      = 'Pa',                                        &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationEdge,               RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec ( gc,                                  &
         SHORT_NAME = 'DQVDT',                                     &
         LONG_NAME  = 'specific_humidity_analysis_increment',      &
         UNITS      = 'kg kg-1',                                   &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec ( gc,                                  &
         SHORT_NAME = 'DO3DT',                                     &
         LONG_NAME  = 'ozone_analysis_increment',                  &
         UNITS      = 'ppmv',                                      &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec ( gc,                                  &
         SHORT_NAME = 'DTSDT',                                     &
         LONG_NAME  = 'skin_temparature_increment',                &
         UNITS      = 'K',                                         &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
    VERIFY_(STATUS)
    
! !EXPORT STATE:

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'DUDT',                                      &
         LONG_NAME  = 'eastward_wind_analysis_increment',          &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'DVDT',                                      &
         LONG_NAME  = 'northward_wind_analysis_increment',         &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'DTDT',                                      &
         LONG_NAME  = 'temperature_analysis_increment',            &
         UNITS      = 'K',                                         &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'DPEDT',                                     &
         LONG_NAME  = 'edge_pressure_analysis_increment',          &
         UNITS      = 'Pa',                                        &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationEdge,               RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'DQVDT',                                     &
         LONG_NAME  = 'specific_humidity_analysis_increment',      &
         UNITS      = 'kg kg-1',                                   &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'DO3DT',                                     &
         LONG_NAME  = 'ozone_analysis_increment',                  &
         UNITS      = 'ppmv',                                      &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'DTSDT',                                     &
         LONG_NAME  = 'skin_temparature_increment',                &
         UNITS      = 'K',                                         &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
    VERIFY_(STATUS)

! Set generic init and final methods
! ----------------------------------

    call MAPL_GenericSetServices    ( gc, RC=STATUS)
    VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)
  
  end subroutine SetServices

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!BOP

! ! IROUTINE: RUN -- Run method for the MAKEIAU component

! !INTERFACE:

subroutine RUN1 ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

  type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
  type(ESMF_State),    intent(inout) :: IMPORT ! Import state
  type(ESMF_State),    intent(inout) :: EXPORT ! Export state
  type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
  integer, optional,   intent(  out) :: RC     ! Error code:

! ! DESCRIPTION: This version uses the MAPL_GenericSetServices. This function sets
!                the Initialize and Finalize services, as well as allocating

!EOP


! ErrLog Variables

  character(len=ESMF_MAXSTR)          :: IAm
  integer                             :: STATUS
  character(len=ESMF_MAXSTR)          :: COMP_NAME

! Local derived type aliases

  type (MAPL_MetaComp),     pointer   :: MAPL
  type (ESMF_State       )            :: INTERNAL

  integer                             :: IM, JM, LM
  real, pointer, dimension(:,:,:)     ::   du,  dudt
  real, pointer, dimension(:,:,:)     ::   dv,  dvdt
  real, pointer, dimension(:,:,:)     ::   dtv, dtdt
  real, pointer, dimension(:,:,:)     ::   dq,  dqdt
  real, pointer, dimension(:,:,:)     ::   do3, do3dt
  real, pointer, dimension(:,:,:)     ::   dple, dpledt
  real, pointer, dimension(:,:)       ::   dts, dtsdt

  type(ESMF_Grid)                     :: grid
  type(ESMF_FieldBundle)              :: bundle

!=============================================================================

! Begin... 

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

   Iam = "Run1"
   call ESMF_GridCompGet( GC, name=COMP_NAME, RC=STATUS )
   VERIFY_(STATUS)
   Iam = trim(COMP_NAME) // Iam

   if ( MAPL_AM_I_ROOT() ) then
       print *, 'Now running ',trim(Iam)
   endif

! Retrieve the pointer to the state
!----------------------------------

   call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
   VERIFY_(STATUS)

! Local aliases to the state, grid, and configuration
! ---------------------------------------------------

!  call MAPL_TimerOn(MAPL,"TOTAL")

! Get parameters from generic state.
!-----------------------------------

    call MAPL_Get( MAPL, IM=IM, JM=JM, LM=LM,    &
                   INTERNAL_ESMF_STATE=INTERNAL, &
                                       RC=STATUS )
    VERIFY_(STATUS)

    call ESMF_GridCompGet(GC, grid=grid, rc=status)
    VERIFY_(STATUS)


! **********************************************************************
! ****               Get Pointers to BKG Import Data                ****
! **********************************************************************
#if 0
    if ( MAPL_AM_I_ROOT() ) then
       call ESMF_StatePrint(IMPORT)
    end if
#endif

!   Get pointers to internal variables
!   ----------------------------------
    call MAPL_GetPointer(import,   du, 'DUDT',  RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(import,   dv, 'DVDT',  RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(import,  dtv, 'DTDT',  RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(import, dple, 'DPEDT', RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(import,   dq, 'DQVDT', RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(import,  do3, 'DO3DT', RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(import,  dts, 'DTSDT', RC=STATUS)
    VERIFY_(STATUS)
    
!   Extract IAU increments in import and copy then to internal
!   ----------------------------------------------------------
#if _MAYBE_
    call ESMF_StateGet(import, 'AGCM_Exports', bundle, rc=status)
    VERIFY_(STATUS)

    call ESMFL_BundleGetPointerToData(bundle, 'DUDT',    du, RC=STATUS)
    VERIFY_(STATUS)
    call ESMFL_BundleGetPointerToData(bundle, 'DVDT',    dv, RC=STATUS)
    VERIFY_(STATUS)
    call ESMFL_BundleGetPointerToData(bundle, 'DTDT',   dtv, RC=STATUS)
    VERIFY_(STATUS)
    call ESMFL_BundleGetPointerToData(bundle, 'DPEDT', dple, RC=STATUS)
    VERIFY_(STATUS)
    call ESMFL_BundleGetPointerToData(bundle, 'DQVDT',   dq, RC=STATUS)
    VERIFY_(STATUS)
    call ESMFL_BundleGetPointerToData(bundle, 'DO3DT',  do3, RC=STATUS)
    VERIFY_(STATUS)
    call ESMFL_BundleGetPointerToData(bundle, 'DTSDT',  dts, RC=STATUS)
    VERIFY_(STATUS)
#endif /* _MAYBE_ */

!   Get pointers to internal variables
!   ----------------------------------
    call MAPL_GetPointer(internal,   dudt, 'DUDT',  RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(internal,   dvdt, 'DVDT',  RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(internal,   dtdt, 'DTDT',  RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(internal, dpledt, 'DPEDT', RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(internal,   dqdt, 'DQVDT', RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(internal,  do3dt, 'DO3DT', RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(internal,  dtsdt, 'DTSDT', RC=STATUS)
    VERIFY_(STATUS)

    dudt   = du
    dvdt   = dv
    dtdt   = dtv
    dpledt = dple
    dqdt   = dq
    do3dt  = do3
    dtsdt  = dts

    RETURN_(ESMF_SUCCESS)
  end subroutine RUN1

! ! IROUTINE: RUN2 -- Run method for the MAKEIAU component

! !INTERFACE:

subroutine RUN2 ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

  type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
  type(ESMF_State),    intent(inout) :: IMPORT ! Import state
  type(ESMF_State),    intent(inout) :: EXPORT ! Export state
  type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
  integer, optional,   intent(  out) :: RC     ! Error code:

! ! DESCRIPTION: This version uses the MAPL_GenericSetServices. This function sets
!                the Initialize and Finalize services, as well as allocating

!EOP


! ErrLog Variables

  character(len=ESMF_MAXSTR)          :: IAm
  integer                             :: STATUS
  character(len=ESMF_MAXSTR)          :: COMP_NAME

! Local derived type aliases

  type (MAPL_MetaComp),     pointer   :: MAPL
  type (ESMF_State       )            :: INTERNAL

  integer                             :: IM, JM, LM
  real, pointer, dimension(:,:,:)     ::   dudt, ex_du
  real, pointer, dimension(:,:,:)     ::   dvdt, ex_dv
  real, pointer, dimension(:,:,:)     ::   dtdt, ex_dt
  real, pointer, dimension(:,:,:)     ::   dqvdt, ex_dqv
  real, pointer, dimension(:,:,:)     ::   do3dt, ex_do3
  real, pointer, dimension(:,:,:)     ::   dpledt, ex_dple
  real, pointer, dimension(:,:)       ::   dtsdt, ex_dts

  type(ESMF_Grid)                     :: grid
  type(ESMF_FieldBundle)              :: bundle

!=============================================================================

! Begin... 

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

   Iam = "Run2"
   call ESMF_GridCompGet( GC, name=COMP_NAME, RC=STATUS )
   VERIFY_(STATUS)
   Iam = trim(COMP_NAME) // Iam

   if ( MAPL_AM_I_ROOT() ) then
       print *, 'Now running ',trim(Iam)
   endif

! Retrieve the pointer to the state
!----------------------------------

   call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
   VERIFY_(STATUS)

! Local aliases to the state, grid, and configuration
! ---------------------------------------------------

!  call MAPL_TimerOn(MAPL,"TOTAL")

! Get parameters from generic state.
!-----------------------------------

    call MAPL_Get( MAPL, IM=IM, JM=JM, LM=LM,    &
                   INTERNAL_ESMF_STATE=INTERNAL, &
                                       RC=STATUS )
    VERIFY_(STATUS)

    call ESMF_GridCompGet(GC, grid=grid, rc=status)
    VERIFY_(STATUS)


! **********************************************************************
! ****               Get Pointers to BKG Import Data                ****
! **********************************************************************
#if 0
    if ( MAPL_AM_I_ROOT() ) then
       call ESMF_StatePrint(IMPORT)
    end if
#endif

!   Get pointers to internal variables
!   ----------------------------------
    call MAPL_GetPointer(internal,   dudt, 'DUDT',  RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(internal,   dvdt, 'DVDT',  RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(internal,   dtdt, 'DTDT',  RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(internal, dpledt, 'DPEDT', RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(internal,  dqvdt, 'DQVDT', RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(internal,  do3dt, 'DO3DT', RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(internal,  dtsdt, 'DTSDT', RC=STATUS)
    VERIFY_(STATUS)
    
!   Get pointers to export variables
!   --------------------------------
    call MAPL_GetPointer(export,   ex_du, 'DUDT',  RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(export,   ex_dv, 'DVDT',  RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(export,   ex_dt, 'DTDT',  RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(export, ex_dple, 'DPEDT', RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(export,  ex_dqv, 'DQVDT', RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(export,  ex_do3, 'DO3DT', RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(export,  ex_dts, 'DTSDT', RC=STATUS)
    VERIFY_(STATUS)

    if (associated(ex_du)) then
       ex_du=dudt
    endif
    if (associated(ex_dv)) then
       ex_dv=dvdt
    endif
    if (associated(ex_dt)) then
       ex_dt=dtdt
    endif
    if (associated(ex_dple)) then
       ex_dple=dpledt
    endif
    if (associated(ex_dqv)) then
       ex_dqv=dqvdt
    endif
    if (associated(ex_do3)) then
       ex_do3=do3dt
    endif
    if (associated(ex_dts)) then
       ex_dts=dtsdt
    endif

    RETURN_(ESMF_SUCCESS)
  end subroutine RUN2
end module IAU_GridCompMod
