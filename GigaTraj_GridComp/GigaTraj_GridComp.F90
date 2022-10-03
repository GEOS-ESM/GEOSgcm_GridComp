#include "MAPL_Generic.h"

!=============================================================================
!BOP

! !MODULE: GigaTraj_GridCompMod -- A Module to run gigatraj

! !INTERFACE:

module GigaTraj_GridCompMod

! !USES:

  use ESMF
  use MAPL

  implicit none
  private

! !PUBLIC MEMBER FUNCTIONS:

  public :: SetServices

contains

!BOP

! !IROUTINE: SetServices -- Sets ESMF services for this component

! !INTERFACE:

    subroutine SetServices ( GC, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    integer, optional                  :: RC  ! return code

    character(len=ESMF_MAXSTR)              :: IAm
    integer                                 :: STATUS
    character(len=ESMF_MAXSTR)              :: COMP_NAME

    type (ESMF_Config)                  :: CF
    type (MAPL_MetaComp),  pointer      :: MAPL

! Begin...

! Get my name and set-up traceback handle
! ---------------------------------------

    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // 'SetServices'

! Register services for this component
! ------------------------------------

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_INITIALIZE, Initialize, RC=STATUS )
    VERIFY_(STATUS)
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,  Run       , RC=STATUS )
    VERIFY_(STATUS)

! Get the configuration from the component
!-----------------------------------------

    call ESMF_GridCompGet( GC, CONFIG = CF, RC=STATUS )
    VERIFY_(STATUS)

! Set the state variable specs.
! -----------------------------
    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)


! Clocks
!-------

    call MAPL_TimerAdd(GC, name="INITIALIZE"    ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC, name="RUN"           ,RC=STATUS)
    VERIFY_(STATUS)

! All done
!---------

    RETURN_(ESMF_SUCCESS)  
  end subroutine SetServices


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!BOP
 
! !IROUTINE: Initialize -- Initialize method for the composite GigaTraj Gridded Component

! !INTERFACE:

  subroutine Initialize ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

  type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
  type(ESMF_State),    intent(inout) :: IMPORT ! Import state
  type(ESMF_State),    intent(inout) :: EXPORT ! Export state
  type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
  integer, optional,   intent(  out) :: RC     ! Error code

! !DESCRIPTION: 
 

!EOP

! ErrLog Variables

  character(len=ESMF_MAXSTR)           :: IAm 
  integer                              :: STATUS
  character(len=ESMF_MAXSTR)           :: COMP_NAME

! Local derived type aliases

   type (MAPL_MetaComp),  pointer  :: STATE
   type (ESMF_State),         pointer  :: GIM(:)
   type (ESMF_State),         pointer  :: GEX(:)
   type (ESMF_Field)                   :: FIELD
   type (ESMF_Time)                    :: CurrTime, RingTime
   type (ESMF_TimeInterval)            :: TIMEINT
   type (ESMF_Alarm)                   :: ALARM
   type (ESMF_Alarm)                   :: ALARM4D
   type (ESMF_Config)                  :: cf
   integer                             :: I, NQ
   real                                :: POFFSET, DT             
   real, pointer, dimension(:,:)       :: PHIS,SGH,VARFLT,PTR
   real, pointer, dimension(:,:,:)     :: TEND!
   character(len=ESMF_MAXSTR)          :: replayMode
   real                                :: RPL_INTERVAL
   real                                :: RPL_SHUTOFF
   real                                :: IAU4dFREQ
   integer                             :: PREDICTOR_DURATION
   integer                             :: MKIAU_FREQUENCY
   character(len=ESMF_MAXSTR), parameter :: INITIALIZED_EXPORTS(3) = &
        (/'PHIS  ', 'SGH   ', 'VARFLT' /)

   logical                             :: DasMode
   character(len=ESMF_MAXSTR)          :: STRING
   character(len=ESMF_MAXSTR)          :: rplMode

! =============================================================================

! Begin... 

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    call ESMF_GridCompGet ( GC, name=COMP_NAME, config=cf, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // "Initialize"

! Get my MAPL_Generic state
!--------------------------

    call MAPL_GetObjectFromGC ( GC, STATE, RC=STATUS)
    VERIFY_(STATUS)


    call MAPL_TimerOn(STATE,"INITIALIZE")

! Call Initialize for every Child

    call MAPL_GenericInitialize ( GC, IMPORT, EXPORT, CLOCK,  RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_TimerOn(STATE,"TOTAL")

    call MAPL_TimerOff(STATE,"TOTAL")
    call MAPL_TimerOff(STATE,"INITIALIZE")


    RETURN_(ESMF_SUCCESS)
 end subroutine Initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!BOP

! !IROUTINE: Run -- Run method for the composite Agcm Gridded Component

! !INTERFACE:

  subroutine Run ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

  type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
  type(ESMF_State),    intent(inout) :: IMPORT ! Import state
  type(ESMF_State),    intent(inout) :: EXPORT ! Export state
  type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
  integer, optional,   intent(  out) :: RC     ! Error code

! !DESCRIPTION: 
 

!EOP

! ErrLog Variables

  character(len=ESMF_MAXSTR)           :: IAm 
  integer                              :: STATUS
  character(len=ESMF_MAXSTR)           :: COMP_NAME


  end subroutine Run

end module GigaTraj_GridCompMod
