!  $Id$

#include "MAPL_Generic.h"

!=============================================================================

module GEOS_OceanbioGridCompMod

!BOP

! !MODULE: GEOS_OceanbioGridCompMod -- Simulates ocean biology module


! !USES:

  use ESMF
  use MAPL
#ifdef USE_ODAS
      use obio_iodas_iau_mod
#endif
  
  implicit none
  private

! !PUBLIC MEMBER FUNCTIONS:

  public SetServices

!=============================================================================

! !DESCRIPTION:
! 
!   {\tt GEOS\_Obio} is a light-weight gridded component does ocean biology
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

! !DESCRIPTION: This version uses the MAPL\_GenericSetServices, which sets
!                the Initialize and Finalize services, as well as allocating
!   our instance of a generic state and putting it in the 
!   gridded component (GC). Here we only need to set the run method and
!   add the state variable specifications (also generic) to our instance
!   of the generic state. This is the way our true state variables get into
!   the ESMF\_State INTERNAL, which is in the MAPL\_MetaComp. The import
!   and internal variables are allocated and initialized by generic.  Here
!   generic is used for tiles.

!EOP

!=============================================================================

! ErrLog Variables

    character(len=ESMF_MAXSTR)              :: IAm
    integer                                 :: STATUS
    character(len=ESMF_MAXSTR)              :: COMP_NAME

!=============================================================================

! Begin...

! Get my name and set-up traceback handle
! ---------------------------------------

    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // 'SetServices'

! Set the Run entry point
! -----------------------

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_INITIALIZE, Initialize, RC=STATUS )
    VERIFY_(STATUS)
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN, Run, RC=STATUS )
    VERIFY_(STATUS)

! Set the state variable specs.
! -----------------------------

!BOS

!  !Export state:

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'TURB',                              &
        LONG_NAME          = 'water_turbidity',                   &
        UNITS              = 'm-1',                               &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

#ifdef USE_ODAS
    call MAPL_AddInternalSpec(GC,                             &
    SHORT_NAME = 'CHLOROPHYLL',                               &
    LONG_NAME  = 'chlorophyll_concentration',                 &
    UNITS      = 'mg m-3',                                    &
    DIMS       = MAPL_DimsHorzVert,                           &
    VLOCATION  = MAPL_VLocationCenter,                        &
    DEFAULT    = 0.25,                                         &
    FRIENDLYTO = 'ORAD:OANA:OCEAN',                           &
                                                   RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
    SHORT_NAME = 'OCEANCOLOR',                                &
    LONG_NAME  = 'surface_chlorophyll_concentration',         &
    UNITS      = 'mg m-3',                                    &
    DIMS       = MAPL_DimsHorzOnly,                           &
    VLOCATION  = MAPL_VLocationNone,                          &
                                                   RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                               &
    SHORT_NAME = 'CHLOROPHYLLinc',                            &
    LONG_NAME  = 'chlorophyll_concentration_increment',       &
    UNITS      = 'mg m-3',                                    &
    DIMS       = MAPL_DimsHorzVert,                           &
    VLOCATION  = MAPL_VLocationCenter,                        &
    RESTART    = .false.,                                     &
                                                   RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                               &
        SHORT_NAME = 'DH',                                    &
        LONG_NAME  = 'Layer mass',                            &
        UNITS      = 'dyn-m',                                 &
        DIMS       = MAPL_DimsHorzVert,                       &
        VLOCATION  = MAPL_VLocationCenter,                    &
        RESTART    = .false.,                                 &
                                                   RC=STATUS  )
    VERIFY_(STATUS)
#endif

!EOS

! Set the Profiling timers
! ------------------------
    
    call MAPL_TimerAdd(GC,    name="INITIALIZE"   ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="RUN"   ,RC=STATUS)
    VERIFY_(STATUS)
  
! Set generic init and final methods
! ----------------------------------

    call MAPL_GenericSetServices    ( GC, RC=STATUS)
    VERIFY_(STATUS)

!  All done
!-----------

    RETURN_(ESMF_SUCCESS)
  
  end subroutine SetServices

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine INITIALIZE ( GC, IMPORT, EXPORT, CLOCK, RC )

  type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
  type(ESMF_State),    intent(inout) :: IMPORT ! Import state
  type(ESMF_State),    intent(inout) :: EXPORT ! Export state
  type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
  integer, optional,   intent(  out) :: RC     ! Error code:

  character(len=ESMF_MAXSTR)          :: IAm
  integer                             :: STATUS
  character(len=ESMF_MAXSTR)          :: COMP_NAME

  type (mapl_metacomp), pointer   :: mapl
  type (esmf_state)               :: internal

#ifdef USE_ODAS
   real, pointer :: chlorophyll(:, :, :) => null(), ocean_color(:, :) => null()
   integer :: i
#endif


!=============================================================================


! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    call ESMF_GridCompGet( GC, name=COMP_NAME, RC=STATUS)
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // "Initialize"

! Get my internal MAPL_Generic state
!-----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)

! Start Total timer
!------------------

   call MAPL_TimerOn(MAPL,"TOTAL")
   call MAPL_TimerOn(MAPL,"INITIALIZE" )
   call MAPL_GenericInitialize( GC, IMPORT, EXPORT, CLOCK, RC=STATUS)
   VERIFY_(STATUS)

! Pointers to outputs
!--------------------
#ifdef USE_ODAS
   call MAPL_Get(mapl, internal_esmf_state = internal, RC=STATUS )
   VERIFY_(STATUS)
   call MAPL_GetPointer(export, ocean_color, 'OCEANCOLOR', alloc = .true., RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(internal, chlorophyll, 'CHLOROPHYLL', RC=STATUS)
   VERIFY_(STATUS)
   if(associated(chlorophyll)) then
       if(all(mask = (chlorophyll == 0.25))) forall(i = 1:size(array = chlorophyll, dim = 3)) chlorophyll(:, :, i) = 0.25*exp(-0.5*(i - 1.0))  
       if(associated(ocean_color)) ocean_color = maxval(array = chlorophyll, dim = 3)
   endif 
#endif

!  All done
!-----------

   call MAPL_TimerOff(MAPL,"INITIALIZE")

   RETURN_(ESMF_SUCCESS)

 end subroutine INITIALIZE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!BOP

! !IROUTINE: RUN -- First Run stage for the Obio component

! !INTERFACE:

subroutine RUN ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

  type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
  type(ESMF_State),    intent(inout) :: IMPORT ! Import state
  type(ESMF_State),    intent(inout) :: EXPORT ! Export state
  type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
  integer, optional,   intent(  out) :: RC     ! Error code:

! !DESCRIPTION: Periodically refreshes the penetratin radiation

!EOP


! ErrLog Variables

  character(len=ESMF_MAXSTR)          :: IAm
  integer                             :: STATUS
  character(len=ESMF_MAXSTR)          :: COMP_NAME

! Locals

  type (MAPL_MetaComp), pointer   :: MAPL
  type (ESMF_State)               :: INTERNAL

! pointers to export

   real, pointer, dimension(:,:)  :: TURB

   real, parameter  :: PTURB  = 0.1

#ifdef USE_ODAS
   real, pointer :: chlorophyll(:, :, :) => null(), dh(:, :, :) => null(), ocean_color(:, :) => null()
   integer :: i
#endif

!=============================================================================

! Begin... 

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    call ESMF_GridCompGet( GC, name=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // "Run"

! Get my internal MAPL_Generic state
!-----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_Get(MAPL, INTERNAL_ESMF_STATE = INTERNAL, RC=STATUS)
    VERIFY_(STATUS)

! Start Total timer
!------------------

   call MAPL_TimerOn(MAPL,"TOTAL")
   call MAPL_TimerOn(MAPL,"RUN" )

! Pointers to outputs
!--------------------

   call MAPL_GetPointer(EXPORT,TURB  , 'TURB'    ,    RC=STATUS)
   VERIFY_(STATUS)

   if(associated(TURB)) TURB = PTURB

#ifdef USE_ODAS
   call MAPL_GetPointer(INTERNAL, chlorophyll, 'CHLOROPHYLL', RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT, ocean_color, 'OCEANCOLOR', RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT, dh, 'DH', RC=STATUS)
   VERIFY_(STATUS)
   if(associated(chlorophyll)) then  
       chlorophyll = merge(mask = (dh /= mapl_undef), tsource = max(chlorophyll, 0.0), fsource = 0.0)  
       call apply_iau(import, mapl, mask = merge(mask = (dh /= mapl_undef), tsource = 1.0, fsource = 0.0))
       if(associated(ocean_color)) ocean_color = maxval(array = merge(mask = (dh /= mapl_undef), tsource = 1.0, fsource = 0.0)*chlorophyll, dim = 3) 
   endif 
#endif

!  All done
!-----------

   call MAPL_TimerOff(MAPL,"RUN" )
   call MAPL_TimerOff(MAPL,"TOTAL")

   RETURN_(ESMF_SUCCESS)

 end subroutine RUN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module GEOS_OceanbioGridCompMod

