#include "MAPL_Generic.h"

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1     !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: StandAlone_DynAdvCore_GridCompMod
!
! !INTERFACE:

      module StandAlone_DynAdvCore_GridCompMod 
!
! !USES:

      use ESMF
      use MAPL_Mod
      use AdvCore_GridCompMod,    only : AdvCoreSetServices   => SetServices
      use FVdycoreCubed_GridComp, only : DynCoreSetServices   => SetServices

      implicit none
      private

! !PUBLIC MEMBER FUNCTIONS:

      public SetServices
!
      integer :: DynCore  = -1
      integer :: AdvCore  = -1

! !DESCRIPTION: 
!
!  {\tt StandAlone\_DynAdvCore\_GridCompMod} is an ESMF gridded component implementing
!  DynCore and AdvCore.
!
! !REVISION HISTORY:
! 18March2014 Kouatchou First crack.
!
!EOP
!------------------------------------------------------------------------------
contains
!------------------------------------------------------------------------------
!BOP
! !IROUTINE: SetServices - Externally visible registration routine
!
! !INTERFACE:
!
      subroutine SetServices(GC, rc)
!
! !ARGUMENTS:
      type(ESMF_GridComp), intent(inout) :: GC
      integer, optional,   intent(  out) :: RC
!
! !DESCRIPTION:
!
!     User-supplied setservices routine.
!     The register routine sets the subroutines to be called
!     as the init, run, and finalize routines.  Note that those are
!     private to the module.
!
!EOP

      character(len=ESMF_MAXSTR)              :: IAm
      integer                                 :: STATUS
      character(len=ESMF_MAXSTR)              :: COMP_NAME

!=============================================================================

! Begin...

      ! Get my name and set-up traceback handle
      ! ---------------------------------------

      call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
      _VERIFY(STATUS)
      Iam = trim(COMP_NAME) // 'SetServices'

! Register methods with MAPL
! --------------------------

      call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_INITIALIZE,  Initialize, RC=status )
      _VERIFY(STATUS)
      call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,         Run,        RC=status )
      _VERIFY(STATUS)
      call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_FINALIZE,    Finalize,   RC=status )
      _VERIFY(STATUS)

      ! Create childrens gridded components and invoke their SetServices
      !-----------------------------------------------------------------

      DynCore = MAPL_AddChild(GC, NAME='DYN',   SS=DynCoreSetServices, RC=status)
      _VERIFY(STATUS)

      AdvCore = MAPL_AddChild(GC, NAME='ADV',   SS=AdvCoreSetServices, RC=status)
      _VERIFY(STATUS)

      call MAPL_TimerAdd(GC, name="INITIALIZE"    ,RC=STATUS)
      _VERIFY(STATUS)

      call MAPL_TimerAdd(GC, name="RUN"           ,RC=STATUS)
      _VERIFY(STATUS)

      call MAPL_TimerAdd(GC, name="TOTAL"         ,RC=STATUS)
      _VERIFY(STATUS)

      ! AdvCore Imports
      ! ---------------
      CALL MAPL_AddConnectivity ( GC,                                   &
                 SHORT_NAME  = (/'MFX ', 'MFY ', 'CX  ' , 'CY  ', 'PLE0', 'PLE1'/),   &
                 DST_ID      = AdvCore,                                 &
                 SRC_ID      = DynCore,                                 &
                                                             RC=STATUS  )
      _VERIFY(STATUS)


      ! Ending with a Generic SetServices call is a MAPL requirement 
      !-------------------------------------------------------------
      call MAPL_GenericSetServices    ( GC, rc=STATUS)
      _VERIFY(STATUS)

      ! All done
      ! --------
      !if ( MAPL_AM_I_ROOT() ) print *, trim(Iam) // ': done!'

      _RETURN(ESMF_SUCCESS)

      end subroutine SetServices
!EOC
!------------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Initialize_ --- Initialize Example
!
! !INTERFACE:
!
      SUBROUTINE Initialize ( GC, IMPORT, EXPORT, CLOCK, rc )

! !USES:

      implicit NONE

! !INPUT PARAMETERS:

      type(ESMF_Clock),  intent(inout) :: CLOCK     ! The clock

! !OUTPUT PARAMETERS:

      type(ESMF_GridComp), intent(inout)  :: GC     ! Grid Component
      type(ESMF_State), intent(inout) :: IMPORT     ! Import State
      type(ESMF_State), intent(inout) :: EXPORT     ! Export State
      integer, intent(out)            ::  rc        ! Error return code:
                                                 !  0 - all is well
                                                 !  1 - 

! !DESCRIPTION: This is a simple ESMF wrapper.

!EOP

!
! !LOCAL VARIABLES:
      character(len=ESMF_MAXSTR)         :: IAm
      integer                            :: STATUS
      type(ESMF_Grid)               :: GRID        ! Grid
      type(ESMF_Config)             :: CF          ! Universal Config 

      integer                       :: im, jm, lm  ! 3D Dimensions
      real(ESMF_KIND_R4), pointer   :: lons(:,:)   ! Longitudes
      real(ESMF_KIND_R4), pointer   :: lats(:,:)   ! Latitudes

      integer                       :: nymd, nhms  ! date, time
      real                          :: cdt         ! time step in secs

      character(len=ESMF_MAXSTR)    :: comp_name
      type (ESMF_FieldBundle)             :: BUNDLE1, BUNDLE2
      type (ESMF_GridComp),      pointer  :: GCS(:)
      type (ESMF_State),         pointer  :: GIM(:)
      type (ESMF_State),         pointer  :: GEX(:)
      type (MAPL_MetaComp), pointer :: MAPL
      INTEGER :: numTracers, I
      type (ESMF_Field)            :: field
!-------------------------------------------------------------------------
!BOC
      Iam = 'Initialize'

!  Get my name and set-up traceback handle
!  ---------------------------------------
      call ESMF_GridCompGet( GC, name=COMP_NAME, RC=status )
      _VERIFY(STATUS)
      Iam = trim(comp_name) // trim(Iam)

     !if ( MAPL_AM_I_ROOT() ) print *, trim(Iam) // ':  Generic Init...'

      !  Create grid for this GC
      !  ------------------------
      call MAPL_GridCreate  (GC, RC=status )
      _VERIFY(STATUS)

      call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
      _VERIFY(STATUS)

!  Initialize MAPL Generic
!  -----------------------

      call MAPL_GenericInitialize ( gc, IMPORT, EXPORT, clock,  RC=status )
      _VERIFY(STATUS)

      ! Get children and their im/ex states from my generic state.
      !----------------------------------------------------------

      call MAPL_Get ( MAPL, GCS=GCS, GIM=GIM, GEX=GEX, RC=STATUS )
      _VERIFY(STATUS)

      call ESMF_StateGet  (GIM(DynCore), 'TRADV', BUNDLE1, RC=STATUS )
      _VERIFY(STATUS)

      call ESMF_FieldBundleGet(BUNDLE1, fieldCount=numTracers,  rc=STATUS)
      _VERIFY(STATUS)

      IF ( MAPL_AM_I_ROOT() ) PRINT*, 'Number of tracers: ', numTracers

      call ESMF_StateGet  (GIM(AdvCore), 'TRADV', BUNDLE2, RC=STATUS )
      _VERIFY(STATUS)

      if (numTracers > 0) then
         do I=1, numTracers
            call ESMF_FieldBundleGet (BUNDLE1, fieldIndex=I, field=FIELD, RC=STATUS)
            _VERIFY(STATUS)

            call MAPL_FieldBundleAdd ( BUNDLE2, field, rc=STATUS )
            _VERIFY(STATUS)

         end do
      end if


!  All done
!  --------
     !if ( MAPL_AM_I_ROOT() ) print *, trim(Iam) // ': done!'
      _RETURN(ESMF_SUCCESS)

      END SUBROUTINE Initialize
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Run - run routine
!
! !INTERFACE:
!
      subroutine Run(GC, IMPORT, EXPORT, CLOCK, RC)
!
! !INPUT/OUTPUT PARAMETERS:
      type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
      type(ESMF_State),    intent(inout) :: IMPORT ! Import state
      type(ESMF_State),    intent(inout) :: EXPORT ! Export state
      type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
!
! !OUTPUT PARAMETERS:
      integer, optional,   intent(  out) :: RC     ! Error code
!
! !DESCRIPTION:
! 
! The Run method advanced the advection one long time step, as
! specified in the configuration.  This may be broken down int a
! number of internal, small steps, also configurable.
!
!EOP
!=============================================================================
!BOC
! !LOCAL VARIABLES:
      character(len=ESMF_MAXSTR)   :: IAm
      integer                      :: STATUS
      character(len=ESMF_MAXSTR)   :: COMP_NAME

      type(ESMF_Grid)               :: GRID        ! Grid
      type(ESMF_Config)             :: CF          ! Universal Config 
      type (MAPL_MetaComp), pointer :: MAPL

      integer                       :: im, jm, lm  ! 3D Dimensions
      real(ESMF_KIND_R4), pointer   :: lons(:,:)   ! Longitudes
      real(ESMF_KIND_R4), pointer   :: lats(:,:)   ! Latitudes

      integer                       :: nymd, nhms  ! date, time
      real                          :: cdt         ! time step in secs

      type (ESMF_GridComp),      pointer  :: GCS(:)
      type (ESMF_State),         pointer  :: GIM(:)
      type (ESMF_State),         pointer  :: GEX(:)
      type (ESMF_State)                   :: INTERNAL
      character(len=ESMF_MAXSTR),pointer  :: GCNames(:)
      INTEGER :: I

! Get my name and set-up traceback handle
! ---------------------------------------
      Iam = 'Run'
      call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
      _VERIFY(STATUS)
      Iam = trim(COMP_NAME) // trim(Iam)

     !if ( MAPL_AM_I_ROOT() ) print *, trim(Iam) // ':  Generic Run...'

      ! Get my internal MAPL_Generic state
      !-----------------------------------

      call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
      _VERIFY(STATUS)

      call MAPL_TimerOn(MAPL, "TOTAL")
      call MAPL_TimerOn(MAPL, "RUN")

      ! Get the children`s states from the generic state
      !-------------------------------------------------

      call MAPL_Get ( MAPL, &
                      GCS=GCS, GIM=GIM, GEX=GEX,       &
                      GCNames = GCNames, &
                      INTERNAL_ESMF_STATE = INTERNAL,  &
                      RC=STATUS )
      _VERIFY(STATUS)

      ! Call Run Method for Children
      I = DynCore
      call ESMF_GridCompRun(GCS(I), &
                            importState = GIM(I), &
                            exportState = GEX(I), &
                            clock       = CLOCK, &
                            userRC=STATUS)
      _VERIFY(STATUS)

      I = AdvCore
      call ESMF_GridCompRun(GCS(I), &
                            importState = GIM(I), &
                            exportState = GEX(I), &
                            clock       = CLOCK, &
                            userRC=STATUS)
      _VERIFY(STATUS)

      call MAPL_TimerOff(MAPL, "TOTAL")
      call MAPL_TimerOff(MAPL, "RUN")

     !IF ( MAPL_AM_I_ROOT() ) PRINT*, TRIM(Iam) // ': done!'

      _RETURN(ESMF_SUCCESS)

      end subroutine Run
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Finalize - user supplied finalize routine
!
! !INTERFACE:
!
      subroutine Finalize(GC, IMPORT, EXPORT, CLOCK, RC)
!
! !INPUT/OUTPUT PARAMETERS:
      type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
      type(ESMF_State),    intent(inout) :: IMPORT ! Import state
      type(ESMF_State),    intent(inout) :: EXPORT ! Export state
      type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
!
! !OUTPUT PARAMETERS:
      integer, optional,   intent(  out) :: RC     ! Error code
!
! !DESCRIPTION:
!    Finalize merely destroys the FVadv object that was created in Initialize
!    and releases the space for the persistent data .
!
!EOP
!=============================================================================
!BOC
! !LOCAL VARIABLES:

      character(len=ESMF_MAXSTR)    :: IAm
      integer                       :: STATUS
      character(len=ESMF_MAXSTR)    :: COMP_NAME

! Get my name and set-up traceback handle
! ---------------------------------------

      Iam = 'Finalize'
      call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
      _VERIFY(STATUS)
      Iam = trim(COMP_NAME) // TRIM(Iam)

     !if ( MAPL_AM_I_ROOT() ) print *, trim(Iam) // ':  Generic Fin...'

      call MAPL_GenericFinalize(GC, IMPORT, EXPORT, CLOCK, RC)
      _VERIFY(STATUS)

     !IF ( MAPL_AM_I_ROOT() ) PRINT*, TRIM(Iam) // ': done!'

      _RETURN(ESMF_SUCCESS)
      end subroutine Finalize
!EOC
!------------------------------------------------------------------------------




      end module StandAlone_DynAdvCore_GridCompMod 
