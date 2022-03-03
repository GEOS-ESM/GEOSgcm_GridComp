!  $Id$

#include "MAPL_Generic.h"

! GEOS   default real kind

#define G5KIND      4
#define REAL_       real(kind=G5KIND)

module CICE_GEOSPlugMod

!BOP
! !MODULE: CICE_GEOSPlugMod -- to couple with CICE6 and later.

!DESCRIPTION:
! A  MAPL/ESMF Gridded Component that acts as a coupler for CICE6 and any future versions.
!

!USES:
  use ESMF
  use MAPL
  use CICE_InitMod                 
  use CICE_FinalMod                 


  implicit none
  private

  !PUBLIC MEMBER FUNCTIONS:
  public :: SetServices



contains

  
!BOP

! !IROUTINE: SetServices -- Sets ESMF services for this component

! !INTERFACE:

    subroutine SetServices ( GC, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    integer, optional                  :: RC  ! return code

! !DESCRIPTION:  The SetServices for the CICE_GEOSplug GC needs to register its
!   Initialize, Run and Finalize.  It uses the MAPL_Generic construct for defining
!   state specs and couplings among its children.  

!EOP

!=============================================================================
!
! ErrLog Variables

    character(len=ESMF_MAXSTR)          :: IAm
    integer                             :: STATUS
    character(len=ESMF_MAXSTR)          :: COMP_NAME

! Locals

!=============================================================================

! Begin...

! Get my name and set-up traceback handle
! ---------------------------------------

    Iam = 'SetServices'
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // Iam

!BOS


  !*CALLBACK*
  !=================================================================================
  ! an ESMF state to pass information b.w. GCs using callback
  ! to be connected to the import in SeaiceInterface
  !
  call MAPL_AddExportSpec(GC                                ,&
          SHORT_NAME         = 'SURFSTATE'                  ,&
          LONG_NAME          = 'surface_state_for_seaice_thermo_coupling',  &
          UNITS              = 'W m-2'                      ,&
          !UNGRIDDED_DIMS     = (/NUM_ICE_CATEGORIES/),      &
          !DIMS               = MAPL_DimsTileOnly,           &
          !VLOCATION          = MAPL_VLocationNone,          &
          DATATYPE           = MAPL_StateItem,               &
                                                       __RC__)

 !=================================================================================



!EOS

! Set the Initialize, Run, Finalize entry points
! ----------------------------------------------

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_INITIALIZE,   Initialize, RC=status)
    VERIFY_(STATUS)
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,          Run,        RC=status)
    VERIFY_(STATUS)
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_FINALIZE,     Finalize,   RC=status)
    VERIFY_(STATUS)
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_WRITERESTART, Record,     RC=status)
    VERIFY_(STATUS)

! Set the Profiling timers
! ------------------------

    call MAPL_TimerAdd(GC,   name="INITIALIZE" ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,   name="RUN"        ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,   name="FINALIZE"   ,RC=STATUS)
    VERIFY_(STATUS)

! Generic SetServices
! -------------------

    call MAPL_GenericSetServices    ( GC, RC=STATUS )
    VERIFY_(STATUS)

! All done
! --------

    RETURN_(ESMF_SUCCESS)

  end subroutine SetServices


!=============================================================================

!BOP

! !IROUTINE: INITIALIZE -- Initialize method for CICE wrapper

! !INTERFACE:

  subroutine Initialize ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

    type(ESMF_GridComp),     intent(INOUT) :: GC     ! Gridded component
    type(ESMF_State),        intent(INOUT) :: IMPORT ! Import state
    type(ESMF_State),        intent(INOUT) :: EXPORT ! Export state
    type(ESMF_Clock),        intent(INOUT) :: CLOCK  ! The clock
    integer, optional,       intent(  OUT) :: RC     ! Error code:

!EOP

! ErrLog Variables

    character(len=ESMF_MAXSTR)         :: IAm
    integer                            :: STATUS
    character(len=ESMF_MAXSTR)         :: COMP_NAME

! Locals

    integer                                :: counts(7)
    integer                                :: Comm
    integer                                :: isc,iec,jsc,jec
    integer                                :: isd,ied,jsd,jed
    integer                                :: IM, JM
    integer                                :: g_isc,g_iec,g_jsc,g_jec
    integer                                :: g_isd,g_ied,g_jsd,g_jed

    integer                                :: YEAR,MONTH,DAY,HR,MN,SC

! Locals with MOM types

    type(time_type)                        :: Time
    type(time_type)                        :: DT

! Locals with ESMF and MAPL types

    type(ESMF_VM)                          :: VM
    type(MAPL_MetaComp), pointer           :: MAPL
    type(ESMF_Grid)                        :: Grid
    type(ESMF_Time)                        :: MyTime
    type(ESMF_TimeInterval)                :: TINT

! Locals

    type(ESMF_State)                       :: SURFST

    type(ice_ocean_boundary_type), pointer :: Boundary                => null()
    type(ocean_public_type),       pointer :: Ocean                   => null()
    type(ocean_state_type),        pointer :: Ocean_State             => null()
    type(MOM_MAPL_Type),           pointer :: MOM_MAPL_internal_state => null()
    type(MOM_MAPLWrap_Type)                :: wrap

    type(ocean_grid_type),         pointer :: Ocean_grid              => null()

    REAL_, pointer                         :: TW  (:,:)        => null()
    REAL_, pointer                         :: SW  (:,:)        => null()
    REAL_, pointer                         :: AREA(:,:)        => null()
    REAL_, pointer                         :: MASK(:,:)        => null()

    integer                                :: DT_SEAICE

! Begin...

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    Iam = "Initialize"
    call ESMF_GridCompGet( gc, NAME=comp_name, RC=status )
    VERIFY_(STATUS)
    Iam = trim(comp_name) // trim(Iam)


! Get my internal MAPL_Generic state
!-----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)

! Profilers
!----------

    call MAPL_TimerOn(MAPL,"TOTAL"     )
    call MAPL_TimerOn(MAPL,"INITIALIZE")

! Get the grid, configuration
!----------------------------

    call ESMF_GridCompGet( GC, grid=Grid,  RC=status )
    VERIFY_(STATUS)

! Get the layout from the grid
!-----------------------------

    call ESMF_VMGetCurrent(VM, rc=STATUS)
    VERIFY_(STATUS)

! Set the time for MOM
!---------------------

    call ESMF_ClockGet(CLOCK, currTIME=MyTime, TimeStep=TINT,  RC=STATUS)
    VERIFY_(STATUS)

    call ESMF_TimeGet (MyTime,                    &
                       YY=YEAR, MM=MONTH, DD=DAY, &
                       H=HR,    M =MN,    S =SC,  &
                                        RC=STATUS )
    VERIFY_(STATUS)

    CALL ESMF_TimeIntervalGet(TINT, S=DT_OCEAN, RC=status)
    VERIFY_(status)

! Allocate this instance of the internal state and wrap
! -----------------------------------------------------

    allocate ( MOM_MAPL_internal_state, stat=status )
    VERIFY_(STATUS)

    wrap%ptr => MOM_MAPL_internal_state

! Save pointer to the wrapped internal state in the GC
! ----------------------------------------------------

    call ESMF_UserCompSetInternalState ( GC, 'MOM_MAPL_state', WRAP, STATUS )
    VERIFY_(STATUS)

    allocate ( Boundary, stat=STATUS); VERIFY_(STATUS)
    allocate ( Ocean,    stat=STATUS); VERIFY_(STATUS)

    MOM_MAPL_internal_state%Ice_ocean_boundary => Boundary
    MOM_MAPL_internal_state%Ocean              => Ocean

! FMS initialization using the communicator from the VM
!------------------------------------------------------

    call ESMF_VMGet(VM, mpiCommunicator=Comm, rc=STATUS)
    VERIFY_(STATUS)


! Init CICE 
!---------------
    ! BZ: need to properly initialze CICE's calendar??
    call CICE_Initialize(Comm)


! Initialize CICE model
!-----------------------

    Ocean%is_ocean_pe = .true.
    call ocean_model_init  (Ocean, Ocean_state, Time, Time, iwind_stagger)
 
    MOM_MAPL_internal_state%Ocean_State => Ocean_State

    call ocean_model_init_sfc(Ocean_state, Ocean)

! Get the ocean grid and sizes of global and computational domains
!-----------------------------------------------------------------

    call get_ocean_grid (Ocean_state, Ocean_grid)
    isc  = Ocean_grid%isc; iec  = Ocean_grid%iec
    isd  = Ocean_grid%isd; ied  = Ocean_grid%ied

    jsc  = Ocean_grid%jsc; jec  = Ocean_grid%jec
    jsd  = Ocean_grid%jsd; jed  = Ocean_grid%jed

! Check local sizes of horizontal dimensions
!--------------------------------------------
    call MAPL_GridGet(GRID, localCellCountPerDim=counts, RC=status)
    VERIFY_(STATUS)

    IM=iec-isc+1
    JM=jec-jsc+1

    ASSERT_(counts(1)==IM)
    ASSERT_(counts(2)==JM)



    !*CALLBACK*
    !=====================================================================================
    call ESMF_StateGet(EXPORT, 'SURFSTATE', SURFST, __RC__)

    !!attach the thermo coupling method
    !
    call ESMF_MethodAdd(SURFST, label='thermo_coupling', userRoutine=thermo_coupling, __RC__)

    ! create and add some fields to the callback state (i.e. SURFSTATE)
    block
      type(ESMF_Field) :: fld
      type (ESMF_Grid) :: grid

      call ESMF_GridCompGet(gc, grid=grid, __RC__)
      fld = MAPL_FieldCreateEmpty('surface_ice_temperature', grid, __RC__)
      call MAPL_FieldAllocCommit(fld, dims=MAPL_DimsHorzOnly,          &
           location=MAPL_VLocationNone, typekind=MAPL_R4, hw=0,        &
           ungrid=[NUM_ICE_CATEGORIES], __RC__)
      call MAPL_StateAdd(SURFST, fld, __RC__)
    end block

    !=====================================================================================


! Profilers
! ---------

    call MAPL_TimerOff(MAPL,"INITIALIZE")
    call MAPL_TimerOff(MAPL,"TOTAL"     )

! Generic initialize
! ------------------

    call MAPL_GenericInitialize( GC, IMPORT, EXPORT, CLOCK, RC=status )
    VERIFY_(STATUS)

! Make sure exports neede by the parent prior to our run call are initialized
!----------------------------------------------------------------------------


! All Done
!---------

    RETURN_(ESMF_SUCCESS)
  end subroutine Initialize


!=================================================================================

!BOP

! !IROUTINE: Run  -- Run method for External Model Plug

! !INTERFACE:

  subroutine Run  ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: GC       ! Gridded component
    type(ESMF_State),    intent(INOUT) :: IMPORT   ! Import state
    type(ESMF_State),    intent(INOUT) :: EXPORT   ! Export state
    type(ESMF_Clock),    intent(INOUT) :: CLOCK    ! The supervisor clock
    integer, optional,   intent(  OUT) :: RC       ! Error code:
    type(ESMF_State)                   :: INTERNAL ! Internal state

!EOP

! ErrLog Variables

    character(len=ESMF_MAXSTR)         :: IAm
    integer                            :: STATUS
    character(len=ESMF_MAXSTR)         :: COMP_NAME

! Locals with ESMF and MAPL types

    type(MAPL_MetaComp),       pointer :: MAPL               => null()
    type(ESMF_Time)                    :: MyTime
    type(ESMF_TimeInterval)            :: TINT

! Locals

    type(ice_ocean_boundary_type), pointer :: Boundary                 => null()
    type(ocean_public_type),       pointer :: Ocean                    => null()
    type(ocean_state_type),        pointer :: Ocean_State              => null()
    type(MOM_MAPL_Type),           pointer :: MOM_MAPL_internal_state  => null()
    type(MOM_MAPLWrap_Type)                :: wrap

!   type(ocean_grid_type),         pointer :: Ocean_grid               => null()

! Required exports

    REAL_, pointer                     :: TW    (:,:)        => null()
    REAL_, pointer                     :: SW    (:,:)        => null()
    REAL_, pointer                     :: UW    (:,:)        => null()
    REAL_, pointer                     :: VW    (:,:)        => null()
    REAL_, pointer                     :: UWB   (:,:)        => null()
    REAL_, pointer                     :: VWB   (:,:)        => null()
    REAL_, pointer                     :: SLV   (:,:)        => null()
    REAL_, pointer                     :: FRAZIL(:,:)        => null()
    REAL_, pointer                     :: MELT_POT(:,:)      => null()
    REAL_, pointer                     :: FRZMLT(:,:)        => null()
    REAL_, pointer                     :: MASK  (:,:)        => null()
    REAL_, pointer                     :: AREA  (:,:)        => null()

! Optional Exports
! none

! Imports
    REAL_, pointer                     :: TAUX(:,:)          => null()
    REAL_, pointer                     :: TAUY(:,:)          => null()
    REAL_, pointer                     :: PS  (:,:)          => null()
    REAL_, pointer                     :: PICE(:,:)          => null()
    REAL_, pointer                     :: LWFLX(:,:)         => null()
    REAL_, pointer                     :: SHFLX(:,:)         => null()
    REAL_, pointer                     :: QFLUX(:,:)         => null()
    REAL_, pointer                     :: RAIN(:,:)          => null()
    REAL_, pointer                     :: SNOW(:,:)          => null()
    REAL_, pointer                     :: SFLX(:,:)          => null()
    REAL_, pointer                     :: PENUVR(:,:)        => null()
    REAL_, pointer                     :: PENPAR(:,:)        => null()
    REAL_, pointer                     :: PENUVF(:,:)        => null()
    REAL_, pointer                     :: PENPAF(:,:)        => null()
    REAL_, pointer                     :: DRNIR(:,:)         => null()
    REAL_, pointer                     :: DFNIR(:,:)         => null()
    REAL_, pointer                     :: DISCHARGE(:,:)     => null()
    REAL_, pointer                     :: AICE(:,:)          => null()
    REAL_, pointer                     :: TAUXBOT(:,:)       => null()
    REAL_, pointer                     :: TAUYBOT(:,:)       => null()

! Temporaries

    real, allocatable                  :: U (:,:),  V(:,:)
    real, allocatable                  :: cos_rot(:,:)
    real, allocatable                  :: sin_rot(:,:)

    integer                            :: IM, JM

    integer                            :: steady_state_ocean = 0       ! SA: Per Atanas T, "name" of this var is misleading
                                                                       ! We run ocean model only when it = 0

    character(len=7)                   :: pres_loading                 ! yes or no

    integer                            :: isc,iec,jsc,jec

    integer                            :: YEAR,MONTH,DAY,HR,MN,SC
    type(time_type)                    :: Time
    type(time_type)                    :: DT

    real                               :: pice_scaling = 1.0
    integer                            :: DT_OCEAN


    REAL_, pointer, dimension(:,:)     :: LATS  => null()
    REAL_, pointer, dimension(:,:)     :: LONS  => null()

! Begin
!------

! Get the component name and set-up traceback handle.
! -----------------------------------------------------
    Iam = "Run"
    call ESMF_GridCompGet( gc, NAME=comp_name, RC=status )
    VERIFY_(status)
    Iam = trim(comp_name) // Iam

! Get my internal MAPL_Generic state
!-----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=status)
    VERIFY_(STATUS)


    call MAPL_Get(MAPL,                      &
         INTERNAL_ESMF_STATE = INTERNAL,     &
         LATS  = LATS ,                      &
         LONS  = LONS ,                      &
                                RC=STATUS )
    VERIFY_(STATUS)

! Profilers
!----------

    call MAPL_TimerOn (MAPL,"TOTAL")
    call MAPL_TimerOn (MAPL,"RUN"  )

! Get the Plug private internal state
!--------------------------------------

    CALL ESMF_UserCompGetInternalState( GC, 'MOM_MAPL_state', WRAP, STATUS )
    VERIFY_(STATUS)

    MOM_MAPL_internal_state => WRAP%PTR

! Aliases to MOM types
!---------------------

    Boundary    => MOM_MAPL_internal_state%Ice_ocean_boundary
    Ocean       => MOM_MAPL_internal_state%Ocean
    Ocean_State => MOM_MAPL_internal_state%Ocean_State

! Get domain size
!----------------

! -------
! do this:
    call mpp_get_compute_domain(Ocean%Domain, isc, iec, jsc, jec)
! instead of:
!   call get_ocean_grid (Ocean_state, Ocean_grid)
!   isc  = Ocean_grid%isc; iec  = Ocean_grid%iec
!   jsc  = Ocean_grid%jsc; jec  = Ocean_grid%jec
! -------

    IM=iec-isc+1
    JM=jec-jsc+1

! Temporaries with MOM default reals
!-----------------------------------

    allocate(U(IM,JM   ),    stat=STATUS); VERIFY_(STATUS)
    allocate(V(IM,JM   ),    stat=STATUS); VERIFY_(STATUS)
    allocate(cos_rot(IM,JM), stat=STATUS); VERIFY_(STATUS)
    allocate(sin_rot(IM,JM), stat=STATUS); VERIFY_(STATUS)

! Get IMPORT pointers
!--------------------

    call MAPL_GetPointer(IMPORT, TAUX,     'TAUX'  ,    RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, TAUY,     'TAUY'  ,    RC=STATUS); VERIFY_(STATUS)


    call MAPL_TimerOff(MAPL,"RUN"   )
    call MAPL_TimerOff(MAPL,"TOTAL" )

! All Done
!---------
    RETURN_(ESMF_SUCCESS)
  end subroutine Run

!BOP

!====================================================================

! !IROUTINE: Finalize        -- Finalize method for CICE wrapper

! !INTERFACE:

  subroutine Finalize ( gc, import, export, clock, rc )

! !ARGUMENTS:

  type(ESMF_GridComp), intent(INOUT) :: GC     ! Gridded component
  type(ESMF_State),    intent(INOUT) :: IMPORT ! Import state
  type(ESMF_State),    intent(INOUT) :: EXPORT ! Export state
  type(ESMF_Clock),    intent(INOUT) :: CLOCK  ! The supervisor clock
  integer, optional,   intent(  OUT) :: RC     ! Error code

!EOP

    type(MAPL_MetaComp),           pointer :: MAPL
    type(ESMF_Time)                        :: MyTime
    type(MOM_MAPL_Type),           pointer :: MOM_MAPL_internal_state => null()
    type(MOM_MAPLWrap_Type)                :: wrap
    type(ocean_public_type),       pointer :: Ocean                   => null()
    type(ocean_state_type),        pointer :: Ocean_State             => null()
    type(ice_ocean_boundary_type), pointer :: Boundary                => null()

! ErrLog Variables

    character(len=ESMF_MAXSTR)       :: IAm
    integer                          :: STATUS
    character(len=ESMF_MAXSTR)       :: COMP_NAME

! Locals with MOM types

    type(time_type)                  :: Time
    integer                          :: YEAR,MONTH,DAY,HR,MN,SC

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    Iam = "Finalize"
    call ESMF_GridCompGet( gc, NAME=comp_name, RC=status )
    VERIFY_(STATUS)
    Iam = trim(comp_name) // Iam

! Get my internal MAPL_Generic state
!-----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=status)
    VERIFY_(STATUS)

! Profilers
!----------

    call MAPL_TimerOn(MAPL,"TOTAL"   )
    call MAPL_TimerOn(MAPL,"FINALIZE")


! Set the times for MOM
!----------------------

    call ESMF_ClockGet( CLOCK, currTime=MyTime, RC=STATUS)
    VERIFY_(status)

    call ESMF_TimeGet (MyTime,      &
         YY=YEAR, MM=MONTH, DD=DAY, &
         H=HR,    M =MN,    S =SC,  &
         RC=STATUS )
    VERIFY_(STATUS)

    call CICE_Finalize !BZ: note save restarts is in ice_step??


    call MAPL_TimerOff(MAPL,"FINALIZE")
    call MAPL_TimerOff(MAPL,"TOTAL"   )

! Generic Finalize
! ------------------

    call MAPL_GenericFinalize( GC, IMPORT, EXPORT, CLOCK, RC=status )
    VERIFY_(STATUS)

! All Done
!---------

    RETURN_(ESMF_SUCCESS)
  end subroutine Finalize

!====================================================================


end module CICE_GEOSPlugMod
