! $Id: GEOS_AgcmSimpleGridComp.F90,v 1.1 2007/05/16 15:33:09 trayanov Exp $

#include "MAPL_Generic.h"


module GEOS_AgcmSimpleGridCompMod

!=============================================================================
!BOP
! \renewcommand{\comp}{\tt GEOS\_AgcmSimpleGridCompMod}
!
! !MODULE: GEOS_AgcmSimpleGridCompMod
!
! !DESCRIPTION: This gridded component (GC) combines the the GC that
!   implements the Finite-Volume (FV) dynamics, with a simple physics
!   component that implements the Held-Suarez benchmark forcing for 
!   testing dry dynamical cores. 
!
! !USES:

  use ESMF
  use MAPL_Mod
  use GEOS_TopoGetMod

  use GEOS_superdynGridCompMod,  only:  SDYN_SetServices => SetServices
  use GEOS_hsGridCompMod,     only:  PHS_SetServices => SetServices

  implicit none
  private

  integer :: SDYN
  integer :: PHS

! !PUBLIC MEMBER FUNCTIONS:

  public SetServices

 
!EOP

contains

!BOP

! !IROUTINE: SetServices

! !DESCRIPTION:  This is the only method of this component, since all of its
!   registered methods can be defaulted to their MAPL generic versions.
!   SetServices merely creates the children through MAPL and connects
!   their Import-Export states. FV has some Imports that are
!   not used in the Held-Suarez benchmark. Since these have proper defaults
!   in FV, they can simply be terminated here. 


! !INTERFACE:

    subroutine SetServices ( GC, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    integer, optional,   intent(  OUT) :: RC  ! return code


!EOP
!BOC
!=============================================================================
!
! ErrLog Variables

    character(len=ESMF_MAXSTR)              :: IAm
    integer                                 :: STATUS
    character(len=ESMF_MAXSTR)              :: COMP_NAME

! Locals


!=============================================================================

! Begin...

! Get my name and set-up traceback handle
! ---------------------------------------

    Iam = 'SetServices'
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // trim(Iam)

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_INITIALIZE, Initialize, RC=STATUS )
    VERIFY_(STATUS)

! dummy import for testing comcurrent ens
!    call MAPL_AddImportSpec ( gc,                                  &
!         SHORT_NAME = 'DTDT',                                      &
!         LONG_NAME  = 'temperature increment',                     &
!         UNITS      = 'K s-1',                                     &
!         DIMS       = MAPL_DimsHorzVert,                           &
!         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
!    VERIFY_(STATUS)

! these internal spec are "fake" and here only to provide moisture to FV
! this component will not touch them
    call MAPL_AddInternalSpec(GC,                                  &
         SHORT_NAME = 'Q',                                         &
         LONG_NAME  = 'specific_humidity',                         &
         UNITS      = 'kg kg-1',                                   &
         FRIENDLYTO = 'DYNAMICS',    &
         default    = 1.0e-6,                                      &
!         RESTART    = MAPL_RestartRequired,                        &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)                                                                          

    call MAPL_AddInternalSpec(GC,                                        &
         SHORT_NAME = 'QLLS',                                            &
         LONG_NAME  = 'mass_fraction_of_large_scale_cloud_liquid_water', &
         UNITS      = 'kg kg-1',                                         &
         FRIENDLYTO = 'DYNAMICS',                             &
         DIMS       = MAPL_DimsHorzVert,                                 &
         VLOCATION  = MAPL_VLocationCenter,                   RC=STATUS  )  
    VERIFY_(STATUS)                                                                          

    call MAPL_AddInternalSpec(GC,                                       &
         SHORT_NAME = 'QLCN',                                           &
         LONG_NAME  = 'mass_fraction_of_convective_cloud_liquid_water', &
         UNITS      = 'kg kg-1',                                        &
         FRIENDLYTO = 'DYNAMICS',                            &
         DIMS       = MAPL_DimsHorzVert,                                &
         VLOCATION  = MAPL_VLocationCenter,                  RC=STATUS  )  
    VERIFY_(STATUS)                                                                          
    call MAPL_AddInternalSpec(GC,                                     &
         SHORT_NAME = 'QILS',                                         &
         LONG_NAME  = 'mass_fraction_of_large_scale_cloud_ice_water', &
         UNITS      = 'kg kg-1',                                      &
         FRIENDLYTO = 'DYNAMICS',                          &
         DIMS       = MAPL_DimsHorzVert,                              &
         VLOCATION  = MAPL_VLocationCenter,                RC=STATUS  )  
    VERIFY_(STATUS)                                                                          


    call MAPL_AddInternalSpec(GC,                                    &
         SHORT_NAME = 'QICN',                                        &
         LONG_NAME  = 'mass_fraction_of_convective_cloud_ice_water', &
         UNITS      = 'kg kg-1',                                     &
         FRIENDLYTO = 'DYNAMICS',                         &
         DIMS       = MAPL_DimsHorzVert,                             &
         VLOCATION  = MAPL_VLocationCenter,               RC=STATUS  )  
    VERIFY_(STATUS)                                                                          

! Register children with MAPL and go down their SS hierarchy
! ----------------------------------------------------------

    SDYN  = MAPL_AddChild(GC, NAME='SUPERDYNAMICS', SS=SDYN_SetServices, RC=STATUS)
    VERIFY_(STATUS)
    PHS  = MAPL_AddChild(GC, NAME='HSPHYSICS',  SS=PHS_SetServices, RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                          &
         SHORT_NAME = 'T',                                    &
         CHILD_ID   = SDYN,                                    &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                          &
         SHORT_NAME = 'PS',                                   &
         CHILD_ID   = SDYN,                                    &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

! Register connections between children
! -------------------------------------

    call MAPL_AddConnectivity ( GC,                                    &
         SHORT_NAME  = (/ 'DUDT', 'DVDT', 'DTDT' /),                   &
         SRC_ID      =  PHS,                                           &
         DST_ID      =  SDYN,                                           &
                                                            RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddConnectivity ( GC,                                    &
         SRC_NAME = (/ 'U    ', 'V    ', 'T    ', 'PLE  ' /),          &
         DST_NAME = (/ 'U    ', 'V    ', 'TEMP ', 'PLE  ' /),          &
         SRC_ID   =  SDYN,                                              &
         DST_ID   =  PHS,                                              &
                                                            RC=STATUS  )
    VERIFY_(STATUS)


! SetServices clean-up on the way back up through the hierarchy
!--------------------------------------------------------------

    call MAPL_TerminateImport(GC, SHORT_NAME = (/'PHIS ','DPEDT'/), &
         CHILD = SDYN, RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_GenericSetServices( GC, RC=STATUS )
    VERIFY_(STATUS)
   
    RETURN_(ESMF_SUCCESS)
  end subroutine SetServices
!EOC

  subroutine Initialize ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
    type(ESMF_State),    intent(inout) :: IMPORT ! Import state
    type(ESMF_State),    intent(inout) :: EXPORT ! Export state
    type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
    integer, optional,   intent(  out) :: RC     ! Error code

! !DESCRIPTION: The Initialize method of this Gridded Component.

! ErrLog Variables

    character(len=ESMF_MAXSTR)              :: IAm 
    integer                                 :: STATUS
    character(len=ESMF_MAXSTR)              :: COMP_NAME

! local vars
   type (MAPL_MetaComp),      pointer  :: MAPL
   type (ESMF_State),         pointer  :: GIM(:)
   type (ESMF_Field)                   :: FIELD
   type (ESMF_FieldBundle)             :: BUNDLE
   type (ESMF_Config)                  :: cf

    Iam = "Initialize"
    call ESMF_GridCompGet ( GC, name=COMP_NAME, config=cf, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // Iam

    call MAPL_GridCreate(GC, rc=status)
    VERIFY_(STATUS)

    call MAPL_GenericInitialize ( GC, IMPORT, EXPORT, CLOCK,  RC=STATUS)
    VERIFY_(STATUS)

! Get my MAPL_Generic state
!--------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)

! Fill Childrens TOPO variables and Diagnostics
!----------------------------------------------
    call MAPL_Get ( MAPL, GIM=GIM, RC=STATUS )
    VERIFY_(STATUS)

! PHIS ...
!---------
    call ESMF_StateGet( GIM(SDYN), 'PHIS', FIELD, rc=STATUS )
    VERIFY_(STATUS)
    Call GEOS_TopoGet ( cf, MEAN=FIELD, rc=STATUS )
    VERIFY_(STATUS)

! TRADV ...
!----------
    call ESMF_StateGet(GIM(SDYN), 'TRADV', BUNDLE, RC=STATUS )
    VERIFY_(STATUS)
    call MAPL_GridCompGetFriendlies(GC, "DYNAMICS", BUNDLE, RC=STATUS )
    VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)
  end subroutine Initialize

end module GEOS_AgcmSimpleGridCompMod

