! $Id$

! VERIFY_ and RETURN_ macros for error handling.

#include "MAPL_Generic.h"

!=============================================================================
module GEOS_SuperdynGridCompMod

!BOP

! !MODULE: GEOS_SuperdynGridCompMod -- A Module to combine Dynamics and Gravity-Wave-Drag Gridded Components
!
! !USES:

  use ESMF
  use MAPL

  use FVdycore_GridCompMod,     only :    FV_SetServices => SetServices
  use FVdycoreCubed_GridComp,   only :   FV3_SetServices => SetServices
  use ARIESg3_GridCompMod,      only : ARIES_SetServices => SetServices
  use GEOS_DatmoDynGridCompMod, only : DATMO_SetServices => SetServices
  use AdvCore_GridCompMod,      only :   ADV_SetServices => SetServices

  implicit none
  private

! !PUBLIC MEMBER FUNCTIONS:

  public SetServices

! !DESCRIPTION: This gridded component (GC) combines the Dynamics GC and
!   the Gravity Wave Drag GC into a new composite SuperDyn GC.
!
!\vspace{5mm}
!   {\bf Import Couplings}:
!
!   The Import Couplings of the SuperDyn GC are the tendencies of the atmospheric
!   state variables U,V,T,PE (due to external diabatic forcing) in addition to
!   a collection of "Friendly" tracers for advection.  The Friendlies will be searched
!   for moisture for use in virtual effects in both the Dynamics and Gravity Wave Drag
!   parameterization.  If no moisture is found, the SyperDyn will be run dry.
!
!\begin{verbatim}
!       DUDT .... U-Wind                    Tendency (m/s)
!       DVDT .... V-Wind                    Tendency (m/s)
!       DPEDT ... Edge-Pressure             Tendency (Pa/s)
!       DTDT .... Mass-Weighted Temperature Tendency (Pa K/s)
!       TRACER .. Friendly Tracers                   (unknown)
!     If Non-Hydrostatic Dynamics
!       DWDT .... W-Wind                    Tendency (m/s)
!\end{verbatim}
!
!\vspace{5mm}
!   {\bf Run Method}:
!   
!   The run method first calls the Gravity Wave Drag parameterization.
!   The tendencies of the atmospheric state variables created by the GWD are then ADDED to the 
!   SuperDyn Import Couplings (i.e., state variable tendencies due to external diabatic forcing), 
!   which are then used to force the Dynamics GC.
!
!\vspace{5mm}
!   {\bf Export Couplings}:
!
!   The Export Couplings of the SuperDyn GC are the union of the Export
!   Couplings of the individual GCs.  It should be noted that the SuperDyn GC
!   controls the GEOS Topo Utility and produces Topo Variables based on the Grid 
!   defined by the DYN GC.  The Topo Variables are computed during the SuperDyn Initialize 
!   method, and are part of the SuperDyn Export Couplings.  ESMF utilities may be used to regrid 
!   these Topo Variables to other components with differing Grids.
 
!EOP

integer ::          DYN
integer ::          ADV = -1

   contains

!BOP

! !IROUTINE: SetServices -- Sets ESMF services for this component

! !INTERFACE:

    subroutine SetServices ( GC, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT)      :: GC  ! gridded component
    integer, optional,   intent(  OUT)      :: RC  ! return code

! !DESCRIPTION:  The SetServices for the SuperDyn GC needs to register its
!   Initialize, Run, and Finalize methods.  In addition, we need to create the   
!   children GCs (DYN and GWD) and run their respective SetServices.

!EOP

!=============================================================================
!
! ErrLog Variables

    character(len=ESMF_MAXSTR)              :: IAm
    integer                                 :: STATUS
    character(len=ESMF_MAXSTR)              :: COMP_NAME

! Local derived type aliases
    type (MAPL_MetaComp),      pointer      :: MAPL
    character(len=ESMF_MAXSTR)              :: DYCORE

    integer                                 :: SCM_SL

!=============================================================================

! Begin...

! Get my name and set-up traceback handle
! ---------------------------------------

    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // 'SetServices'

    call MAPL_GetObjectFromGC (GC, MAPL,  RC=STATUS )
    VERIFY_(STATUS)

! Register services for this component
! ------------------------------------

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_INITIALIZE,  Initialize, RC=status )
    VERIFY_(STATUS)
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,   Run,        RC=status )
    VERIFY_(STATUS)
    call MAPL_GridCompSetEntryPoint ( gc, ESMF_METHOD_RUN,   RunAddIncs, RC=status)
    VERIFY_(STATUS)


! Create children's gridded components and invoke their SetServices
! -----------------------------------------------------------------

!BOR
    call MAPL_GetResource(MAPL, DYCORE, 'DYCORE:', default="FV3", RC=STATUS )
!EOR
    call MAPL_GetResource(MAPL, SCM_SL, 'SCM_SL:', default=0, RC=STATUS )


    VERIFY_(STATUS)

    if(adjustl(DYCORE)=="FV"   ) then
                  DYN =  MAPL_AddChild(GC, NAME='DYN', SS=   FV_SetServices, RC=STATUS)
                  VERIFY_(STATUS)
    endif
    if(adjustl(DYCORE)=="FV3"  ) then
                  DYN =  MAPL_AddChild(GC, NAME='DYN', SS=  FV3_SetServices, RC=STATUS)
                  VERIFY_(STATUS)
    endif
    if(adjustl(DYCORE)=="FV3+ADV"  ) then
                  DYN =  MAPL_AddChild(GC, NAME='DYN', SS=  FV3_SetServices, RC=STATUS)
                  VERIFY_(STATUS)
                  ADV =  MAPL_AddChild(GC, NAME='ADV', SS=  ADV_SetServices, RC=STATUS)
                  VERIFY_(STATUS)
    endif
    if(adjustl(DYCORE)=="ARIES") then
                  DYN =  MAPL_AddChild(GC, NAME='DYN', SS=ARIES_SetServices, RC=STATUS)
                  VERIFY_(STATUS)
    endif
    if(adjustl(DYCORE)=="DATMO") then
                  DYN =  MAPL_AddChild(GC, NAME='DYN', SS=DATMO_SetServices, RC=STATUS)
                  VERIFY_(STATUS)
    endif

! Set IMPORT specs for the Composite GC. These are the sums of the Phys and Ana tendencies
! ----------------------------------------------------------------------------------------

!BOS
    call MAPL_AddImportSpec ( GC,                                  &
         SHORT_NAME = 'DUDT',                                      &
         LONG_NAME  = 'eastward_wind_tendency',                    &
         UNITS      = 'm s-2',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddImportSpec ( GC,                                  &
         SHORT_NAME = 'DVDT',                                      &
         LONG_NAME  = 'northward_wind_tendency',                   &
         UNITS      = 'm s-2',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddImportSpec ( GC,                                  &
         SHORT_NAME = 'DWDT',                                      &
         LONG_NAME  = 'vertical_velocity_tendency',                &
         UNITS      = 'm s-2',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddImportSpec ( GC,                                  &
         SHORT_NAME = 'DTDT',                                      &
         LONG_NAME  = 'delta-p_weighted_temperature_tendency',     &
         UNITS      = 'Pa K s-1',                                  &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

! Add Exports promoted from child (FV) exports
! --------------------------------------------

    if (SCM_SL /= 0) then
!      print *,'SuperDyn: adding LHOBS and SHOBS exports'
      call MAPL_AddExportSpec(GC,                                  &
           SHORT_NAME='LHOBS',                                     &
           CHILD_ID = DYN,                                         &
                                                         __RC__  )

      call MAPL_AddExportSpec(GC,                                  &
           SHORT_NAME='SHOBS',                                     &
           CHILD_ID = DYN,                                         &
                                                         __RC__  )
    end if

    call MAPL_AddExportSpec ( GC   ,                               &
         SHORT_NAME = 'U',                                         &
         CHILD_ID   = DYN,                                         &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                               &
         SHORT_NAME = 'V',                                         &
         CHILD_ID   = DYN,                                         &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                               &
         SHORT_NAME = 'W',                                         &
         CHILD_ID   = DYN,                                         &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                               &
         SHORT_NAME = 'T',                                         &
         CHILD_ID   = DYN,                                         &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                               &
         SHORT_NAME = 'S',                                         &
         CHILD_ID   = DYN,                                         &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                               &
         SHORT_NAME = 'TH',                                        &
         CHILD_ID   = DYN,                                         &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                               &
         SHORT_NAME = 'PLE',                                       &
         CHILD_ID   = DYN,                                         &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                               &
         SHORT_NAME = 'PL',                                        &
         CHILD_ID   = DYN,                                         &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                               &
         SHORT_NAME = 'ZLE',                                       &
         CHILD_ID   = DYN,                                         &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                               &
         SHORT_NAME = 'PREF',                                      &
         CHILD_ID   = DYN,                                         &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                               &
         SHORT_NAME = 'AK',                                        &
         CHILD_ID   = DYN,                                         &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                               &
         SHORT_NAME = 'BK',                                        &
         CHILD_ID   = DYN,                                         &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                               &
         SHORT_NAME = 'PLK',                                       &
         CHILD_ID   = DYN,                                         &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                               &
         SHORT_NAME = 'PKE',                                       &
         CHILD_ID   = DYN,                                         &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                               &
         SHORT_NAME = 'PS',                                        &
         CHILD_ID   = DYN,                                         &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                               &
         SHORT_NAME = 'DELP',                                      &
         CHILD_ID   = DYN,                                         &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                               &
         SHORT_NAME = 'US',                                        &
         CHILD_ID   = DYN,                                         &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                               &
         SHORT_NAME = 'VS',                                        &
         CHILD_ID   = DYN,                                         &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                               &
         SHORT_NAME = 'TA',                                        &
         CHILD_ID   = DYN,                                         &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                               &
         SHORT_NAME = 'QA',                                        &
         CHILD_ID   = DYN,                                         &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                               &
         SHORT_NAME = 'SPEED',                                     &
         CHILD_ID   = DYN,                                         &
                                                        RC=STATUS  )
    VERIFY_(STATUS)
 
    call MAPL_AddExportSpec ( GC   ,                               &
         SHORT_NAME = 'DZ',                                        &
         CHILD_ID   = DYN,                                         &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                               &
         SHORT_NAME = 'TROPP_BLENDED',                             &
         CHILD_ID   = DYN,                                         &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                               &
         SHORT_NAME = 'TROPK_BLENDED',                             &
         CHILD_ID   = DYN,                                         &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                               &
         SHORT_NAME = 'PV',                                        &
         CHILD_ID   = DYN,                                         &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                               &
         SHORT_NAME = 'TV',                                        &
         CHILD_ID   = DYN,                                         &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                               &
         SHORT_NAME = 'OMEGA',                                     &
         CHILD_ID   = DYN,                                         &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                               &
         SHORT_NAME = 'EPV',                                       &
         CHILD_ID   = DYN,                                         &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                               &
         SHORT_NAME = 'PEANA',                                     &
         CHILD_ID   = DYN,                                         &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                               &
         SHORT_NAME = 'DTHVDTANAINT',                              &
         CHILD_ID   = DYN,                                         &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                               &
         SHORT_NAME = 'PEPHY',                                     &
         CHILD_ID   = DYN,                                         &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                               &
         SHORT_NAME = 'DTHVDTPHYINT',                              &
         CHILD_ID   = DYN,                                         &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                               &
         SHORT_NAME = 'DQVDTANAINT',                               &
         CHILD_ID   = DYN,                                         &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                               &
         SHORT_NAME = 'DQLDTANAINT',                               &
         CHILD_ID   = DYN,                                         &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                               &
         SHORT_NAME = 'DQIDTANAINT',                               &
         CHILD_ID   = DYN,                                         &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                               &
         SHORT_NAME = 'DOXDTANAINT',                               &
         CHILD_ID   = DYN,                                         &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                               &
         SHORT_NAME = 'AREA',                                      &
         CHILD_ID   = DYN,                                         &
                                                        RC=STATUS  )
    VERIFY_(STATUS)
 
    call MAPL_AddExportSpec ( GC   ,                               &
         SHORT_NAME = 'U_DGRID',                                   &
         CHILD_ID   = DYN,                                         &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                               &
         SHORT_NAME = 'V_DGRID',                                   &
         CHILD_ID   = DYN,                                         &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                               &
         SHORT_NAME = 'PT',                                        &
         CHILD_ID   = DYN,                                         &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                               &
         SHORT_NAME = 'PE',                                        &
         CHILD_ID   = DYN,                                         &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                               &
         SHORT_NAME = 'QV_DYN_IN',                                 &
         CHILD_ID   = DYN,                                         &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                               &
         SHORT_NAME = 'T_DYN_IN',                                 &
         CHILD_ID   = DYN,                                         &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                               &
         SHORT_NAME = 'U_DYN_IN',                                 &
         CHILD_ID   = DYN,                                         &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                               &
         SHORT_NAME = 'V_DYN_IN',                                 &
         CHILD_ID   = DYN,                                         &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                               &
         SHORT_NAME = 'PLE_DYN_IN',                                 &
         CHILD_ID   = DYN,                                         &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                               &
         SHORT_NAME = 'DTDTDYN',                                   &
         CHILD_ID   = DYN,                                         &
                                                        RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC   ,                               &
         SHORT_NAME = 'DQVDTDYN',                                   &
         CHILD_ID   = DYN,                                         &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

!needed for NCEP GWD
    call MAPL_AddExportSpec ( GC   ,                               &
         SHORT_NAME = 'DXC',                                       &
         CHILD_ID   = DYN,                                         &
                                                        RC=STATUS  )
    VERIFY_(STATUS)
!EOS

! Create DYN+ADV connectivities
! -----------------------------
    if(adjustl(DYCORE)=="FV3+ADV"  ) then
      CALL MAPL_AddConnectivity ( GC,                               &
                 SHORT_NAME  = (/'MFX ', 'MFY ', 'CX  ' , 'CY  ',   &
                                 'PLE0', 'PLE1', 'QW_BEFORE_DYN',   &
                                 'QW_AFTER_DYN'/),                  &
                 DST_ID      = ADV,                                 &
                 SRC_ID      = DYN,                                 &
                 RC=STATUS  )
      VERIFY_(STATUS)
    endif

! Disable explicitly some connectivities
! the first 3 will be connected explicitly in the run method
! ----------------------------------------------------------

    call MAPL_TerminateImport ( GC,                    &
         SHORT_NAME  = (/'DUDT','DVDT','DTDT','DWDT' /),      &
         CHILD = DYN,                                  &
         RC=STATUS  )
     VERIFY_(STATUS)

! Set the Profiling timers
! ------------------------

    call MAPL_GenericSetServices    ( GC, RC=STATUS)
    VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)
  
  end subroutine SetServices


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!BOP

! !IROUTINE: Initialize -- Initialize method for the composite SuperDyn Gridded Component

! !INTERFACE:

  subroutine Initialize ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

  type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
  type(ESMF_State),    intent(inout) :: IMPORT ! Import state
  type(ESMF_State),    intent(inout) :: EXPORT ! Export state
  type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
  integer, optional,   intent(  out) :: RC     ! Error code

! !DESCRIPTION: The Initialize method of the SuperDyn Composite Gridded Component first 
!   calls the Initialize method of the child Dynamics.  The Dynamics Initialize method will
!   create the ESMF GRID, which will then be used to set the GRID associated with the
!   SuperDyn Composite Component itself.  It should be noted that the 
!   SuperDyn Initialize method also invokes the GEOS Topo Utility which creates all
!   topography related quantities.

!EOP

! ErrLog Variables

  character(len=ESMF_MAXSTR)          :: IAm
  integer                             :: STATUS
  character(len=ESMF_MAXSTR)          :: COMP_NAME

! Local derived type aliases

   type (MAPL_MetaComp),      pointer  :: MAPL
   type (ESMF_Grid )                   :: GRID
   type (ESMF_State),         pointer  :: GIM(:)
   type (ESMF_State),         pointer  :: GEX(:)
   type (ESMF_State)                   :: INTERNAL
  
   type (ESMF_Config)                  :: CF

   real, dimension(:,:,:), pointer     :: PTR3

   integer  unit

!=============================================================================

! Begin... 

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    Iam = "Initialize"
    call ESMF_GridCompGet ( GC, name=COMP_NAME, config=CF, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // trim(Iam)

! Get my internal MAPL_Generic state
!-----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS )
    VERIFY_(STATUS)

! Extract Dynamics ESMF GRID
!---------------------------

    call ESMF_GridCompGet ( GC, GRID=GRID, RC=STATUS )
    VERIFY_(STATUS)

! Call Generic Initialize for SuperDyn GC
!----------------------------------------

    call MAPL_GenericInitialize ( GC, IMPORT, EXPORT, CLOCK, RC=STATUS )
    VERIFY_(STATUS)

! Get parameters from generic state.
!-----------------------------------

    call MAPL_Get ( MAPL, GIM=GIM, GEX=GEX, &
                                INTERNAL_ESMF_STATE=INTERNAL, RC=STATUS )
    VERIFY_(STATUS)

! We will need the Q Export from Dynamics to compute Tv Tendencies from Physics
! -----------------------------------------------------------------------------
    call MAPL_GetPointer(GEX(DYN), PTR3, 'Q', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)

#ifdef PRINT_STATES
    call WRITE_PARALLEL ( trim(Iam)//": IMPORT State" )
    if ( MAPL_am_I_root() ) call ESMF_StatePrint ( IMPORT, rc=STATUS )
    call WRITE_PARALLEL ( trim(Iam)//": EXPORT State" )
    if ( MAPL_am_I_root() ) call ESMF_StatePrint ( EXPORT, rc=STATUS )
#endif

! All done
!---------

   RETURN_(ESMF_SUCCESS)
  end subroutine Initialize


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!BOP

! !IROUTINE: Run -- Run method for the composite SuperDyn Gridded Component

! !INTERFACE:

   subroutine Run ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

  type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
  type(ESMF_State),    intent(inout) :: IMPORT ! Import state
  type(ESMF_State),    intent(inout) :: EXPORT ! Export state
  type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
  integer, optional,   intent(  out) :: RC     ! Error code

! !DESCRIPTION: The run method first calls the Gravity Wave Drag parameterization.
!   The tendencies of the atmospheric state variables created by the GWD are then ADDED to the 
!   SuperDyn Import Couplings (i.e., state variable tendencies due to external diabatic forcing), 
!   which are then used to force the Dynamics GC.
!

!EOP

! ErrLog Variables

   character(len=ESMF_MAXSTR)          :: IAm
   integer                             :: STATUS
   character(len=ESMF_MAXSTR)          :: COMP_NAME

! Local derived type aliases

   type (MAPL_MetaComp),      pointer  :: MAPL
   type (ESMF_State)                   :: INTERNAL

   real, pointer, dimension(:,:,:)     :: DUDT_IMP, DUDT_DYN
   real, pointer, dimension(:,:,:)     :: DVDT_IMP, DVDT_DYN
   real, pointer, dimension(:,:,:)     :: DTDT_IMP, DTDT_DYN
   real, pointer, dimension(:,:,:)     :: DWDT_IMP, DWDT_DYN

   type (ESMF_GridComp),      pointer  :: GCS(:)
   type (ESMF_State),         pointer  :: GIM(:)
   type (ESMF_State),         pointer  :: GEX(:)
  
!=============================================================================

! Begin... 

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    Iam = "Run"
    call ESMF_GridCompGet ( GC, name=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // trim(Iam)

! Get my internal MAPL_Generic state
!-----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_Get ( MAPL, GCS=GCS, GIM=GIM, GEX=GEX, &
         INTERNAL_ESMF_STATE=INTERNAL, RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_TimerOn  (MAPL,"TOTAL")

    call MAPL_GetPointer( IMPORT,   dudt_imp, 'DUDT',  RC=STATUS ); VERIFY_(STATUS)
    call MAPL_GetPointer( IMPORT,   dvdt_imp, 'DVDT',  RC=STATUS ); VERIFY_(STATUS)
    call MAPL_GetPointer( IMPORT,   dtdt_imp, 'DTDT',  RC=STATUS ); VERIFY_(STATUS)
    call MAPL_GetPointer( IMPORT,   dwdt_imp, 'DWDT',  RC=STATUS ); VERIFY_(STATUS)

    call MAPL_GetPointer( GIM(DYN), dudt_dyn, 'DUDT',  RC=STATUS ); VERIFY_(STATUS)
    call MAPL_GetPointer( GIM(DYN), dvdt_dyn, 'DVDT',  RC=STATUS ); VERIFY_(STATUS)
    call MAPL_GetPointer( GIM(DYN), dtdt_dyn, 'DTDT',  RC=STATUS ); VERIFY_(STATUS)
    call MAPL_GetPointer( GIM(DYN), dwdt_dyn, 'DWDT',  RC=STATUS ); VERIFY_(STATUS)

! Add my (SuperDyn) Import tendencies to FV Import tendencies.
!------------------------------------------------------------

    DUDT_DYN = DUDT_IMP
    DVDT_DYN = DVDT_IMP
    DTDT_DYN = DTDT_IMP
    DWDT_DYN = DWDT_IMP

! Call Run Method for Children (FV & ADV)
!----------------------------------
    call ESMF_GridCompRun(GCS(DYN), importState=GIM(DYN), exportState=GEX(DYN), clock=CLOCK, PHASE=1, userRC=STATUS)
    VERIFY_(STATUS)

    if (ADV /= -1) then
      call ESMF_GridCompRun(GCS(ADV), importState=GIM(ADV), exportState=GEX(ADV), clock=CLOCK, userRC=STATUS)
      VERIFY_(STATUS)
    end if

! All done
!---------

    call MAPL_TimerOff (MAPL,"TOTAL")
    RETURN_(ESMF_SUCCESS)
  end subroutine Run


! !IROUTINE: RunAddIncs -- Run method to add increments

! !INTERFACE:

   subroutine RunAddIncs ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

  type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
  type(ESMF_State),    intent(inout) :: IMPORT ! Import state
  type(ESMF_State),    intent(inout) :: EXPORT ! Export state
  type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
  integer, optional,   intent(  out) :: RC     ! Error code

! !DESCRIPTION: The run method to add increments
!

!EOP

! ErrLog Variables

   character(len=ESMF_MAXSTR)          :: IAm
   integer                             :: STATUS
   character(len=ESMF_MAXSTR)          :: COMP_NAME

! Local derived type aliases

   real, pointer, dimension(:,:,:)     :: DUDT_IMP, DUDT_DYN
   real, pointer, dimension(:,:,:)     :: DVDT_IMP, DVDT_DYN
   real, pointer, dimension(:,:,:)     :: DTDT_IMP, DTDT_DYN
   real, pointer, dimension(:,:,:)     :: DWDT_IMP, DWDT_DYN

   type (MAPL_MetaComp),      pointer  :: MAPL
   type (ESMF_GridComp),      pointer  :: GCS(:)
   type (ESMF_State),         pointer  :: GIM(:)
   type (ESMF_State),         pointer  :: GEX(:)
  

!=============================================================================

! Begin... 

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    Iam = "RunAddIncs"
    call ESMF_GridCompGet ( GC, name=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // trim(Iam)

! Get my internal MAPL_Generic state
!-----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_TimerOn  (MAPL,"TOTAL")

! Get the children's states from the generic state
!-------------------------------------------------

    call MAPL_Get ( MAPL, GCS=GCS, GIM=GIM, GEX=GEX, RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_GetPointer( IMPORT,   dudt_imp, 'DUDT',  RC=STATUS ); VERIFY_(STATUS)
    call MAPL_GetPointer( IMPORT,   dvdt_imp, 'DVDT',  RC=STATUS ); VERIFY_(STATUS)
    call MAPL_GetPointer( IMPORT,   dtdt_imp, 'DTDT',  RC=STATUS ); VERIFY_(STATUS)
    call MAPL_GetPointer( IMPORT,   dwdt_imp, 'DWDT',  RC=STATUS ); VERIFY_(STATUS)

    call MAPL_GetPointer( GIM(DYN), dudt_dyn, 'DUDT',  RC=STATUS ); VERIFY_(STATUS)
    call MAPL_GetPointer( GIM(DYN), dvdt_dyn, 'DVDT',  RC=STATUS ); VERIFY_(STATUS)
    call MAPL_GetPointer( GIM(DYN), dtdt_dyn, 'DTDT',  RC=STATUS ); VERIFY_(STATUS)
    call MAPL_GetPointer( GIM(DYN), dwdt_dyn, 'DWDT',  RC=STATUS ); VERIFY_(STATUS)

! Add my (SuperDyn) Import tendencies to FV Import tendencies.
!------------------------------------------------------------

    DUDT_DYN = DUDT_IMP
    DVDT_DYN = DVDT_IMP
    DTDT_DYN = DTDT_IMP
    DWDT_DYN = DWDT_IMP

    call ESMF_GridCompRun(GCS(DYN), importState=GIM(DYN), exportState=GEX(DYN), clock=CLOCK, PHASE=2, userRC=STATUS)
    VERIFY_(STATUS)

    call MAPL_TimerOff  (MAPL,"TOTAL")
    RETURN_(ESMF_SUCCESS)
  end subroutine RunAddIncs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module GEOS_SuperDynGridCompMod
