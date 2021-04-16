

#include "MAPL_Generic.h"

module GEOSseaice_GridCompMod

!BOP
! !MODULE:  GEOSseaice_GridCompMod -- Implements ESMF wrapper to invoke the DATASEAICE/CICE4/CICE6 seaice models.

! !USES:

  use ESMF
  use MAPL
  use GEOS_CICEDynaGridCompMod,          only : CICE4SeaIceSetServices  => SetServices
  use GEOS_DataSeaIceGridCompMod,        only : DataSeaIceSetServices => SetServices
  use ice_prescribed_mod,                only : ice_nudging

   

  implicit none
  private

! !PUBLIC ROUTINES:

  public SetServices

  character(len=ESMF_MAXSTR)          :: SEAICE_NAME
  integer            :: DO_DATASEAICE

! !DESCRIPTION:
!
!   {\tt GEOSseaice\_GridComp} is a light-weight gridded component that serves an
!   interface to seaice/data\_seaice components.
!
!EOP




  integer ::          ICE
  integer ::          ICEd
  logical ::      DUAL_OCEAN


contains

!BOP

! !IROUTINE: SetServices -- Sets ESMF services for GEOSseaice

! !INTERFACE:

  subroutine SetServices ( GC, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    integer, optional                  :: RC  ! return code

! !DESCRIPTION: This version uses the MAPL\_GenericSetServices,
!       which sets the Run, Initialize, and Finalize services,
!       as well as allocating our instance of a generic state and putting it in the
!   gridded component (GC). Here we override all three methods and declare
!       the specs for the Imports and Export States (no MAPL controlled Internal State).
!
!EOP

!=============================================================================
!

! ErrLog Variables

    character(len=ESMF_MAXSTR)         :: IAm
    integer                            :: STATUS
    character(len=ESMF_MAXSTR)         :: COMP_NAME

! Local vars
    type  (MAPL_MetaComp), pointer     :: MAPL
    type  (ESMF_Config)                :: CF
    integer ::      iDUAL_OCEAN
    character(len=ESMF_MAXSTR)         :: charbuf_
    !character(len=ESMF_MAXSTR)         :: sharedObj

! Begin...

! Get my name and set-up traceback handle
! ---------------------------------------

    Iam = 'SetServices'
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, CONFIG=CF, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // Iam


! Set the state variable specs.
! -----------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)

! Get constants from CF
! ---------------------

    call MAPL_GetResource ( MAPL,       DO_DATASEAICE,     Label="USE_DATASEAICE:" ,       DEFAULT=1, RC=STATUS)
    VERIFY_(STATUS)

    if(DO_DATASEAICE/=0) then
       SEAICE_NAME="DATASEAICE"
       ICE = MAPL_AddChild(GC, NAME=SEAICE_NAME, SS=DataSeaiceSetServices, RC=STATUS)
       VERIFY_(STATUS)
    else
       call MAPL_GetResource ( MAPL, SEAICE_NAME, Label="SEAICE_NAME:", DEFAULT="CICE4", __RC__ )
       select case (trim(SEAICE_NAME))
          case ("CICE4")
             ICE = MAPL_AddChild(GC, NAME=SEAICE_NAME, SS=CICE4SeaIceSetServices, RC=STATUS)
             VERIFY_(STATUS)
          case default
             charbuf_ = "SEAICE_NAME: " // trim(SEAICE_NAME) // " is not implemented, ABORT!"
             call WRITE_PARALLEL(charbuf_)
             VERIFY_(999)
       end select
    endif

    call MAPL_GetResource(MAPL, iDUAL_OCEAN, 'DUAL_OCEAN:', default=0, RC=STATUS )
    VERIFY_(STATUS)
    DUAL_OCEAN = iDUAL_OCEAN /= 0

    ICEd = 0
    if (dual_ocean) then
       ICEd = MAPL_AddChild(GC, NAME=SEAICE_NAME, SS=DataSeaiceSetServices, RC=STATUS) 
       VERIFY_(STATUS)
    endif

! Set the state variable specs.
! -----------------------------
