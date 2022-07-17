!  $Id$

#include "MAPL_Generic.h"

! GEOS-5 default real kind

#define G5KIND      4
#define REAL_       real(kind=G5KIND)

#ifndef _RL
#define _RL Real*8
#endif

module MIT_GEOS5PlugMod

!BOP
! !MODULE: MIT_GEOS5PlugMod -- wrapper for MITgcm.

!DESCRIPTION:
! A  MAPL/ESMF Gridded Component that acts as a wrapper for MIT.
! It uses ESMF AND MAPL. It has heavy dependencies on FMS and MIT.
!
! This should be built like MIT, so that its default reals
! are the same as MIT's. It may also be an adequate plug for HIM.
!
! It does not use the configuration.
! Its time step is the clocks time step.
! Each run invocation runs one time step.
!

!USES:
  use ESMF
  use MAPL

  USE MITGCM_STATE_MOD , ONLY :   &
       MITGCM_ISTATE_CONTAINER,       &
       MITGCM_ISTATE,                 &
       MITGCM_ISTATE_WRAP_TYPE,       &
       GETDP

  USE MITGCM_DRIVER_MOD , ONLY :  &
       DRIVER_INIT,                   &
       DRIVER_RUN

  USE STR4C_MOD
  USE DRIVER_SET_IMPORT_STATE_MOD
  USE DRIVER_GET_EXPORT_STATE_MOD

!UDI - What are these? (commente for now)

!  use mpp_parameter_mod,          only: AGRID, SCALAR_PAIR
!  use mpp_io_mod,                 only: MPP_RDONLY, MPP_NETCDF 
!  use mpp_io_mod,                 only: mpp_open, mpp_close
!  use fms_mod,                    only: read_data

! Nothing on the MIT side is visible through this module.

  implicit none
  private

  !PUBLIC MEMBER FUNCTIONS:
  public :: SetServices
!EOP

! These are the MIT-side bulletin boards, where things are in
! MOM's precision and the B grid

  type :: T_PrivateState
     type(MITGCM_ISTATE),  pointer :: istate
     type(ESMF_State)              :: import
     type(ESMF_State)              :: export
  end type T_PrivateState

  type :: T_PrivateState_Wrap
     type(T_PrivateState), pointer :: ptr
  end type T_PrivateState_Wrap

  type(T_PrivateState), pointer :: privateState

  integer            :: NUM_ICE_CATEGORIES
  integer            :: NUM_ICE_LAYERS
  integer, parameter :: NUM_SNOW_LAYERS=1
  integer            :: NUM_ICE_LAYERS_ALL
  integer            :: NUM_SNOW_LAYERS_ALL

  logical :: DUAL_OCEAN

contains


!BOP

! !IROUTINE: SetServices -- Sets ESMF services for this component

! !INTERFACE:

    subroutine SetServices ( GC, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    integer, optional                  :: RC  ! return code

! !DESCRIPTION:  The SetServices for the PhysicsGcm GC needs to register its
!   Initialize and Run.  It uses the MAPL_Generic construct for defining 
!   state specs and couplings among its children.  In addition, it creates the   
!   children GCs (AGCM and OGCM) and runs their
!   respective SetServices.

!EOP

!=============================================================================
!
! ErrLog Variables

    character(len=ESMF_MAXSTR)          :: IAm
    integer                             :: STATUS
    character(len=ESMF_MAXSTR)          :: COMP_NAME

!   Variables for setting MAPL import/export specs
    TYPE MSTATE
     CHARACTER(len=ESMF_MAXSTR)        :: short_name
     CHARACTER(len=ESMF_MAXSTR)        :: long_name
     CHARACTER(len=ESMF_MAXSTR)        :: units
     INTEGER                           :: dims
     INTEGER                           :: vlocation
    END TYPE

    TYPE( MSTATE ), POINTER  :: imports(:)
    TYPE( MSTATE ), POINTER  :: exports(:)
    INTEGER nImports
    INTEGER nExports
    INTEGER I

! Locals
    type (MAPL_MetaComp),  pointer     :: MAPL  
    type  (ESMF_Config)                :: CF
    integer                            :: iDUAL_OCEAN

!=============================================================================

! Begin...

! Get my name and set-up traceback handle
! ---------------------------------------

    Iam = 'SetServices'
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, CONFIG=CF, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // Iam

! Get the MAPL object
! -------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_GetResource(MAPL, iDUAL_OCEAN, 'DUAL_OCEAN:', default=0, RC=STATUS )
    DUAL_OCEAN = iDUAL_OCEAN /= 0


!BOS

!  !IMPORT STATE:

!   Imports and exports specification
!   ---------------------------------

    nimports = 19
    allocate(imports(nimports))
    imports    = (/                                                                                                 &
    mstate('TAUX'     , 'Agrid_eastward_stress_on_skin'           , 'N m-2'     , MAPL_DimsHorzOnly,MAPL_VLocationNone),    &
    mstate('TAUY'     , 'Agrid_northward_stress_on_skin'          , 'N m-2'     , MAPL_DimsHorzOnly,MAPL_VLocationNone),    &
    mstate('LWFLX'    , 'surface_net_downward_longwave_flux'      , 'W m-2'     , MAPL_DimsHorzOnly,MAPL_VLocationNone),    &
    mstate('SHFLX'    , 'upward_sensible_heat_flux'               , 'W m-2'     , MAPL_DimsHorzOnly,MAPL_VLocationNone),    &
    mstate('QFLUX'    , 'evaporation'                             , 'kg m-2 s-1', MAPL_DimsHorzOnly,MAPL_VLocationNone),    & 
    mstate('RAIN'     , 'ocean_rainfall'                          , 'kg m-2 s-1', MAPL_DimsHorzOnly,MAPL_VLocationNone),    &
    mstate('SNOW'     , 'ocean_snowfall'                          , 'kg m-2 s-1', MAPL_DimsHorzOnly,MAPL_VLocationNone),    &
    mstate('DISCHARGE', 'river_discharge_at_ocean_points'         , 'kg m-2 s-1', MAPL_DimsHorzOnly,MAPL_VLocationNone),    &
    mstate('SFLX'     , 'salt_flux_due_to_ice_dynamics'           , 'kg m-2 s-1', MAPL_DimsHorzOnly,MAPL_VLocationNone),    &
    mstate('PS'       , 'Surface Atmospheric Pressure'            , 'Pa'        , MAPL_DimsHorzOnly,MAPL_VLocationNone),    &
    mstate('PENUVR'   ,'net_downward_penetrating_direct_UV_flux'  , 'W m-2'     , MAPL_DimsHorzOnly,MAPL_VLocationNone),    &
    mstate('PENPAR'   ,'net_downward_penetrating_direct_PAR_flux' , 'W m-2'     , MAPL_DimsHorzOnly,MAPL_VLocationNone),    &
    mstate('PENUVF'   ,'net_downward_penetrating_diffuse_UV_flux' , 'W m-2'     , MAPL_DimsHorzOnly,MAPL_VLocationNone),    &
    mstate('PENPAF'   ,'net_downward_penetrating_diffuse_PAR_flux', 'W m-2'     , MAPL_DimsHorzOnly,MAPL_VLocationNone),    &
    mstate('DRNIR'   ,'net_surface_downwelling_nir_beam_flux', 'W m-2'     , MAPL_DimsHorzOnly,MAPL_VLocationNone),    &
    mstate('DFNIR'   ,'net_surface_downwelling_nir_diffuse_flux', 'W m-2'     , MAPL_DimsHorzOnly,MAPL_VLocationNone),    &
    mstate('DEL_TEMP'   ,'temperature correction to top level MIT', 'W m-2'     , MAPL_DimsHorzOnly,MAPL_VLocationNone),    &
    mstate('SWHEAT'   ,'solar_heating_rate'                       , 'W m-2'     , MAPL_DimsHorzVert,MAPL_VLocationCenter),  &
    mstate('WGHT'     , 'weight_for_ocean_grid'                   , '1'         , MAPL_DimsHorzOnly,MAPL_VLocationNone)     &
    /)


    DO I=1,NIMPORTS
     CALL MAPL_AddImportSpec(GC,            &
      SHORT_NAME = imports(i)%short_name,   &
      LONG_NAME  = imports(i)%long_name,    &
      UNITS      = imports(i)%units,        &
      DIMS       = imports(i)%dims,         &
      VLOCATION  = imports(i)%vlocation,    &
      RC         =status); VERIFY_(STATUS)
     call WRITE_PARALLEL("MAPL: adding import "//trim(imports(i)%short_name))

     ! ALT: Mirroring the Imports to Exports for diagnostic purposes
     CALL MAPL_AddExportSpec(GC,            &
      SHORT_NAME = imports(i)%short_name,   &
      LONG_NAME  = imports(i)%long_name,    &
      UNITS      = imports(i)%units,        &
      DIMS       = imports(i)%dims,         &
      VLOCATION  = imports(i)%vlocation,    &
      RC         =status); VERIFY_(STATUS)
     call WRITE_PARALLEL("MAPL: adding export "//trim(imports(i)%short_name))

  ENDDO
    deallocate(imports)

!  !EXPORT STATE:

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME         = 'FRZMLT',                            &
         LONG_NAME          = 'freeze_melt_potential',             &
         UNITS              = 'W m-2',                             &
         DIMS               = MAPL_DimsHorzOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
         RC=STATUS  )
    VERIFY_(STATUS)

! Run1 exports


!   -------------------------

    nexports  = 17
    allocate(exports(nexports))
    exports    = (/ &
    mstate('UW',    'top_layer_Agrid_eastward_velocity',     'm s-1',     MAPL_DimsHorzOnly,MAPL_VLocationNone),    &
    mstate('VW',    'top_layer_Agrid_northward_velocity',    'm s-1',     MAPL_DimsHorzOnly,MAPL_VLocationNone),    &
    mstate('TW',    'top_layer_temperature',                 'K',         MAPL_DimsHorzOnly,MAPL_VLocationNone),    &
    mstate('SW',    'top_layer_salinity',                    'psu',       MAPL_DimsHorzOnly,MAPL_VLocationNone),    &
    mstate('MIT_3D_MASK',  'ocean mask at t-points',         '1',         MAPL_DimsHorzVert,MAPL_VLocationCenter),  &
    mstate('DH',    'layer_thickness',                       'm',         MAPL_DimsHorzVert,MAPL_VLocationCenter),  &
    mstate('SLV',   'sea_level_with_ice_loading',            'm',         MAPL_DimsHorzOnly,MAPL_VLocationNone),    &
    mstate('UWB',   'surface_Bgrid_X_velocity',            'm s-1 ',      MAPL_DimsHorzOnly,MAPL_VLocationNone),    &
    mstate('VWB',   'surface_Bgrid_Y_velocity',            'm s-1 ',      MAPL_DimsHorzOnly,MAPL_VLocationNone),    &
    mstate('SSH',   'sea_level_height',                    'm',           MAPL_DimsHorzOnly,MAPL_VLocationNone),    &
    mstate('T',   'potential_temperature',                 'm',           MAPL_DimsHorzVert,MAPL_VLocationCenter),    &
    mstate('S',   'salinity',                              'psu',         MAPL_DimsHorzVert,MAPL_VLocationCenter),    &
    mstate('DISCHARGEe','river_discharge_at_ocean_points',  'kg m-2 s-1',  MAPL_DimsHorzOnly,MAPL_VLocationNone),    &
    mstate('PBO','pressure_at_bottom_of_ocean',  'N m-2',  MAPL_DimsHorzOnly,MAPL_VLocationNone),    &
    mstate('FRAZIL','heating_from_frazil_formation',  'N m-2',  MAPL_DimsHorzOnly,MAPL_VLocationNone),    &
    mstate('SWFLX'    , 'surface_net_downward_shortwave_flux'     , 'W m-2'     , MAPL_DimsHorzOnly,MAPL_VLocationNone),    &
    mstate('WGHTe', 'weight_for_ocean_grid','1',     MAPL_DimsHorzOnly,MAPL_VLocationNone)    &
     /)

    DO I=1,nexports
     call MAPL_AddExportSpec(GC,                                    &
          SHORT_NAME         = exports(i)%short_name,               &
          LONG_NAME          = exports(i)%long_name,                &
          UNITS              = exports(i)%units,                    &
          DIMS               = exports(i)%dims,                     &
          VLOCATION          = exports(i)%vlocation,                &
          RC=STATUS  )
     VERIFY_(STATUS)
     call WRITE_PARALLEL("MAPL: adding export "//trim(exports(i)%short_name))
    ENDDO

    deallocate(exports)


  call MAPL_AddImportSpec(GC,                            &
    SHORT_NAME         = 'ICESTATES',                      &
    LONG_NAME          = 'container_for_seaice_variables_for_MITgcm', &
    UNITS              = 'N/A',                            &
    DIMS               = MAPL_DimsHorzOnly,                &
    VLOCATION          = MAPL_VLocationNone,               &
    DATATYPE           = MAPL_StateItem,                   &
!    RESTART            = MAPL_RestartSkip,                 &
                                                   RC=STATUS  )
  VERIFY_(STATUS)

!EOS

! Get constants from CF
! ---------------------

    call ESMF_ConfigGetAttribute(CF, NUM_ICE_CATEGORIES, Label="CICE_N_ICE_CATEGORIES:" , RC=STATUS)
    VERIFY_(STATUS)

    call ESMF_ConfigGetAttribute(CF, NUM_ICE_LAYERS,     Label="CICE_N_ICE_LAYERS:" ,     RC=STATUS)
    VERIFY_(STATUS)

    NUM_ICE_LAYERS_ALL  = NUM_ICE_LAYERS  * NUM_ICE_CATEGORIES
    NUM_SNOW_LAYERS_ALL = NUM_SNOW_LAYERS * NUM_ICE_CATEGORIES

!ALT: These were brought from  the original CICEdyna

! Set the Initialize, Run, Finalize entry points
! ----------------------------------------------

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_INITIALIZE,   Initialize, RC=status)
    VERIFY_(STATUS)
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,	    Run,        RC=status)
    VERIFY_(STATUS)
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_FINALIZE,     Finalize,   RC=status)
    VERIFY_(STATUS)
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_WRITERESTART, Record,     RC=status)
    VERIFY_(STATUS)
    if (dual_ocean) then
       call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,	    Run2,        RC=status)
       VERIFY_(STATUS)
    end if

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


!BOP

! !IROUTINE: INITIALIZE -- Initialize method for ExternalOcean wrapper

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

    character(len=ESMF_MAXSTR)		   :: IAm
    integer				   :: STATUS
    character(len=ESMF_MAXSTR)             :: COMP_NAME

! Locals

    integer                                :: Comm

! Locals with ESMF and MAPL types

    type(ESMF_VM)                          :: VM
    type (MAPL_MetaComp), pointer          :: MAPL 

! Locals

!   Variable to hold model state for each instance
    TYPE(MITGCM_ISTATE_CONTAINER) :: mitgcmIState(1)
!   TYPE(MITGCM_ISTATE_WRAP_TYPE) wrap
    TYPE(T_PrivateState_Wrap) wrap
    type(MITGCM_ISTATE),  pointer :: istate

!   Variables for holding and setting run directory
    character(len=ESMF_MAXSTR)            :: ocean_dir
    integer*1, pointer                    :: iarr(:)

!   Local variables used for allocating exports pointers
    REAL_, pointer                         :: TS  (:,:)
    REAL_, pointer                         :: SS  (:,:)
    REAL_, pointer                         :: pMASK(:,:,:)
    REAL_, pointer                         :: DH(:,:,:)
    type(ESMF_State) :: state
    character(len=ESMF_MAXSTR)            :: cname
    integer :: chdir
    external chdir 

! Begin... 

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    Iam = "Initialize"
    call ESMF_GridCompGet( gc, NAME=comp_name, RC=status )
    VERIFY_(STATUS)
    Iam = trim(comp_name) // trim(Iam)

! Allocate the private state...
!------------------------------
    
    allocate( PrivateSTATE , stat=STATUS )
    VERIFY_(STATUS)

    wrap%ptr => PrivateState

! And put it in the GC
!---------------------

    CALL ESMF_UserCompSetInternalState( GC, trim(comp_name)//'_internal_state',&
         WRAP, STATUS )
    VERIFY_(status)

!-------------------
    CALL ESMF_GridCompGet(gc, vm=vm, RC=status); VERIFY_(STATUS)
    CALL ESMF_VMGet(VM, mpiCommunicator=Comm, rc=RC)


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

    call MAPL_GetResource( MAPL, ocean_dir, label='OCEAN_DIR:', rc=status ) ; VERIFY_(STATUS)
    call str4c( iarr, TRIM(ocean_dir) )

! Now do component specific initialization
! ----------------------------------------
!    call system('pushd '//trim(ocean_dir))
!    status = chdir(trim(ocean_dir))
    call mysetdir(iarr)
    call WRITE_PARALLEL("Calling DRIVER_INIT")
    CALL DRIVER_INIT( mitgcmIState=mitgcmIState(1)%p, myComm=Comm)
    call WRITE_PARALLEL("Done DRIVER_INIT")
    deallocate(iarr)
    call popdir
!    status = chdir('..')

    PrivateState%iState => mitgcmIState(1)%p
    istate => PrivateState%iState

    CALL ESMF_UserCompSetInternalState ( GC, 'MITgcm_istate',wrap,status )
    VERIFY_(STATUS)

! Generic initialize
! ------------------

    call MAPL_GenericInitialize( GC, IMPORT, EXPORT, CLOCK, RC=status )
    VERIFY_(STATUS)

!   Force allocation of export arrays for this component
!   ----------------------------------------------------
    call MAPL_GetPointer(EXPORT, DH, 'DH',  alloc=.true., RC=STATUS)
    VERIFY_(STATUS)
!!  DH=1000.
!!  Remove set to 1000 m depths, include info from mitgcm input data file.

!! Hard-wired for now

    DH(:,:, 1) = 50.
    DH(:,:, 2) = 70.
    DH(:,:, 3) = 100.
    DH(:,:, 4) = 140.
    DH(:,:, 5) = 190.
    DH(:,:, 6) = 240.
    DH(:,:, 7) = 290.
    DH(:,:, 8) = 340.
    DH(:,:, 9) = 390.
    DH(:,:,10) = 440.
    DH(:,:,11) = 490.
    DH(:,:,12) = 540.
    DH(:,:,13) = 590.
    DH(:,:,14) = 640.
    DH(:,:,15:) = 690.

    call MAPL_GetPointer(EXPORT, pMASK, trim(COMP_NAME)//'_3D_MASK',  alloc=.true., RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, TS,  'TW',  alloc=.true., RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, SS,  'SW',  alloc=.true., RC=STATUS)
    VERIFY_(STATUS)

    call WRITE_PARALLEL("Calling DRIVER_Get_ExportState")
    CALL DRIVER_GET_EXPORT_STATE(istate, 'MASK', pMASK )
    CALL DRIVER_GET_EXPORT_STATE(istate, 'TS', TS )
    CALL DRIVER_GET_EXPORT_STATE(istate, 'SS', SS )
    call WRITE_PARALLEL("Done DRIVER_Get_ExportState")

    call ESMF_StateGet(IMPORT, 'ICESTATES', state, __RC__)
    call ESMF_AttributeGet(state, name='ICECOMPNAME', value=cname, __RC__)
    call ESMF_StateGet(state, trim(cname)//'_Imports', privateState%import, __RC__)
    call ESMF_StateGet(state, trim(cname)//'_Exports', privateState%export, __RC__)


! Profilers
! ---------

    call MAPL_TimerOff(MAPL,"INITIALIZE")
    call MAPL_TimerOff(MAPL,"TOTAL"     )

! All Done
!---------

    RETURN_(ESMF_SUCCESS)
  end subroutine Initialize


!=================================================================================

!BOP

! !IROUTINE: Run  -- Run method for External Model Plug

! !INTERFACE:

  subroutine Run  ( gc, import, export, clock, rc )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: gc     ! Gridded component 
    type(ESMF_State),    intent(INOUT) :: import ! Import state
    type(ESMF_State),    intent(INOUT) :: export ! Export state
    type(ESMF_Clock),    intent(INOUT) :: clock  ! The supervisor clock
    integer, optional,   intent(  OUT) :: rc     ! Error code:
!    type(ESMF_State)                   :: INTERNAL ! Internal state

!EOP

! ErrLog Variables

    character(len=ESMF_MAXSTR)		   :: IAm
    integer				   :: STATUS
    character(len=ESMF_MAXSTR)             :: COMP_NAME

! Locals

    type (MAPL_MetaComp), pointer          :: MAPL 

    integer :: IM, JM

!   Pointers for passing import state values to component
    REAL_, pointer                         ::   TAUX(:,:  )
    REAL_, pointer                         ::   TAUY(:,:  )
    REAL_, pointer                         ::     PS(:,:  )
    REAL_, pointer                         ::   LWFLX(:,: )
    REAL_, pointer                         ::   SHFLX(:,: )
    REAL_, pointer                         ::   QFLUX(:,: )
    REAL_, pointer                         ::   RAIN(:,: )
    REAL_, pointer                         ::   SNOW(:,: )
    REAL_, pointer                         ::   DISCHARGE(:,:  )
    REAL_, pointer                         ::   SFLX(:,:  )
    REAL_, pointer                         ::   LATS(:,:  )
    REAL_, pointer                         ::   LONS(:,:  )
    REAL_, pointer                         ::   WGHT(:,:  )
    REAL_, pointer                         ::   HFLX(:,:  )
    REAL_, pointer                         ::   QFLX(:,:  )
    REAL_, pointer                         ::   PENUVR(:,:)
    REAL_, pointer                         ::   PENPAR(:,:)
    REAL_, pointer                         ::   PENUVF(:,:)
    REAL_, pointer                         ::   PENPAF(:,:)
    REAL_, pointer                         ::   DRNIR(:,:)
    REAL_, pointer                         ::   DFNIR(:,:)


!   Pointers for mirroring import state values to exports
    REAL_, pointer                         ::   TAUXe(:,:  )
    REAL_, pointer                         ::   TAUYe(:,:  )
    REAL_, pointer                         ::     PSe(:,:  )
    REAL_, pointer                         ::   LWFLXe(:,: )
    REAL_, pointer                         ::   SWFLX(:,: )
    REAL_, pointer                         ::   SHFLXe(:,: )
    REAL_, pointer                         ::   QFLUXe(:,: )
    REAL_, pointer                         ::   RAINe(:,: )
    REAL_, pointer                         ::   SNOWe(:,: )
    REAL_, pointer                         ::   DISCHARGEe(:,:  )
    REAL_, pointer                         ::   SFLXe(:,:  )
    REAL_, pointer                         ::   WGHTe(:,:  )

!   Pointers for fetching export state values from component
    REAL_, pointer                         :: UW  (:,:)
    REAL_, pointer                         :: VW  (:,:)
    REAL_, pointer                         :: TW  (:,:)
    REAL_, pointer                         :: SW  (:,:)
    REAL_, pointer                         :: MASK(:,:,:)

!   Type for getting MITgcm internal state pointer
    TYPE(T_PrivateState_Wrap) wrap
    type(T_PrivateState), pointer :: privateState
    type(MITGCM_ISTATE),  pointer :: istate

!   Variables for holding and setting run directory
    character(len=ESMF_MAXSTR)            :: ocean_dir
    integer*1, pointer                    :: iarr(:)
    integer                               :: active_ocean
    REAL*8                                :: Av

! Sea ice vars
    integer :: i, j
    integer :: C, L, LCI, LCS
    real, parameter :: cutoff = 1.0e-5
    type(ESMF_State) :: importSI, exportSI

    REAL_, pointer                         :: FRACICE(:,:,:)
    REAL_, pointer                         :: VOLICE(:,:,:)
    REAL_, pointer                         :: VOLSNO(:,:,:)
    REAL_, pointer                         :: ERGICE(:,:,:)
    REAL_, pointer                         :: ERGSNO(:,:,:)
    REAL_, pointer                         :: TI(:,:,:)
    REAL_, pointer                         :: SI(:,:)
    REAL_, pointer                         :: HI(:,:)
    REAL_, pointer                         :: MPOND (:,:,:)
    REAL_, pointer                         :: TAUAGE(:,:,:)
    REAL_, pointer                         :: UI(:,:)
    REAL_, pointer                         :: VI(:,:)
    REAL_, pointer                         :: TAUXIe(:,:)
    REAL_, pointer                         :: TAUYIe(:,:)
    REAL_, pointer                         :: TAUXBOT(:,:)
    REAL_, pointer                         :: TAUYBOT(:,:)

    REAL_, pointer                         :: FRACICEe(:,:,:)
    REAL_, pointer                         :: TIe(:,:,:)
    REAL_, pointer                         :: SIe(:,:)
    REAL_, pointer                         :: VOLICEe(:,:,:)
    REAL_, pointer                         :: VOLSNOe(:,:,:)
    REAL_, pointer                         :: ERGICEe(:,:,:)
    REAL_, pointer                         :: ERGSNOe(:,:,:)
    REAL_, pointer                         :: MPONDe(:,:,:)
    REAL_, pointer                         :: TAUAGEe(:,:,:)
    REAL_, pointer                         :: HIe(:,:)

    REAL_, pointer                         :: DEL_FRACICE(:,:,:)
    REAL_, pointer                         :: DEL_TI(:,:,:)
    REAL_, pointer                         :: DEL_SI (:,:)
    REAL_, pointer                         :: DEL_VOLICE(:,:,:)
    REAL_, pointer                         :: DEL_VOLSNO(:,:,:)
    REAL_, pointer                         :: DEL_ERGICE(:,:,:)
    REAL_, pointer                         :: DEL_ERGSNO(:,:,:)
    REAL_, pointer                         :: DEL_MPOND(:,:,:)
    REAL_, pointer                         :: DEL_TAUAGE(:,:,:)
    REAL_, pointer                         :: DEL_HI(:,:)

! Begin
!------


! Get the component's name and set-up traceback handle.
! -----------------------------------------------------
    call WRITE_PARALLEL( ' Starting plug run method ' )
    Iam = "Run"
    call ESMF_GridCompGet( gc, NAME=comp_name, RC=status )
    VERIFY_(status)
    Iam = trim(comp_name) // Iam

! Get the wrapped MIT state
!--------------------------
    call ESMF_UserCompGetInternalState( GC, trim(comp_name)//'_internal_state',&
         WRAP, STATUS )
    VERIFY_(status)

    PrivateState => wrap%ptr
    istate => PrivateState%iState

    importSI = PrivateState%import
    exportSI = PrivateState%export

! Get my internal MAPL_Generic state
!-----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)

! Profilers
!----------

    call MAPL_TimerOn (MAPL,"TOTAL")
    call MAPL_TimerOn (MAPL,"RUN"  )

    Av=2.5E6

! Get IMPORT pointers
!--------------------
    call MAPL_GetPointer(IMPORT,   TAUX,      'TAUX' ,     RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT,   TAUY,      'TAUY' ,     RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT,   PS,        'PS'   ,     RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT,   LWFLX,     'LWFLX',     RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT,   SHFLX,     'SHFLX',     RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT,   QFLUX,     'QFLUX',     RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT,   RAIN,      'RAIN',      RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT,   SNOW,      'SNOW',      RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT,   DISCHARGE, 'DISCHARGE', RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT,   SFLX,      'SFLX',      RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT,   WGHT,      'WGHT',      RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT,   PENUVR,    'PENUVR',    RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT,   PENPAR,    'PENPAR',    RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT,   PENUVF,    'PENUVF',    RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT,   PENPAF,    'PENPAF',    RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT,   DRNIR,     'DRNIR',     RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT,   DFNIR,     'DFNIR',     RC=STATUS); VERIFY_(STATUS)
    call MAPL_Get(MAPL, LATS=LATS, LONS=LONS, RC=status); VERIFY_(STATUS)

! Sea ice vars
    call MAPL_GetPointer(importSI,     HI,     'HI', __RC__)
    call MAPL_GetPointer(importSI,     TI,     'TI', __RC__)
    call MAPL_GetPointer(importSI,     SI,     'SI', __RC__)
    call MAPL_GetPointer(importSI, VOLICE, 'VOLICE', __RC__)
    call MAPL_GetPointer(importSI, VOLSNO, 'VOLSNO', __RC__)
    call MAPL_GetPointer(importSI, ERGICE, 'ERGICE', __RC__)
    call MAPL_GetPointer(importSI, ERGSNO, 'ERGSNO', __RC__)
    call MAPL_GetPointer(importSI, TAUAGE, 'TAUAGE', __RC__)
    call MAPL_GetPointer(importSI,  MPOND,  'MPOND', __RC__)

    call MAPL_GetPointer(importSI, FRACICE,'FRACICE', RC=STATUS)

! Get EXPORT pointers to mirror imports
!--------------------------------------
    call MAPL_GetPointer(EXPORT,   TAUXe,      'TAUX',       RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT,   TAUYe,      'TAUY',       RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT,   PSe,        'PS',         RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT,   LWFLXe,     'LWFLX',      RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT,   SWFLX,      'SWFLX',  alloc=.true.,      RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT,   SHFLXe,     'SHFLX',      RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT,   QFLUXe,     'QFLUX',      RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT,   RAINe,      'RAIN',       RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT,   SNOWe,      'SNOW',       RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT,   WGHTe,      'WGHTe',      RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT,   SFLXe,      'SFLX',      RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT,   DISCHARGEe, 'DISCHARGEe', RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, MASK, trim(COMP_NAME)//'_3D_MASK',  alloc=.true., RC=STATUS); VERIFY_(STATUS)

#ifdef SEAICE_EXPORTS
    CALL MAPL_GetPointer(exportSI, TIe,   'TI', __RC__)
    CALL MAPL_GetPointer(exportSI, SIe,   'SI', __RC__)
    CALL MAPL_GetPointer(exportSI, VOLICEe, 'VOLICE', __RC__)
    CALL MAPL_GetPointer(exportSI, VOLSNOe, 'VOLSNO', __RC__)
    CALL MAPL_GetPointer(exportSI, ERGICEe, 'ERGICE', __RC__)
    CALL MAPL_GetPointer(exportSI, ERGSNOe, 'ERGSNO', __RC__)
    CALL MAPL_GetPointer(exportSI, MPONDe, 'MPOND', __RC__)
    CALL MAPL_GetPointer(exportSI, TAUAGEe, 'TAUAGE', __RC__)
    CALL MAPL_GetPointer(exportSI, HIe, 'HI', __RC__)
#endif

    ! Actual copy (only if needed)
    if (associated(TAUXe)) TAUXe = TAUX
    if (associated(TAUYe)) TAUYe = TAUY
    if (associated(PSe))   PSe = PS
    if (associated(LWFLXe)) LWFLXe = LWFLX
    if (associated(SHFLXe)) SHFLXe = SHFLX
    if (associated(QFLUXe)) QFLUXe = QFLUX
    if (associated(RAINe)) RAINe = RAINe
    if (associated(SNOWe)) SNOWe = SNOWe
    if (associated(SFLXe)) SFLXe = SFLX
    if (associated(DISCHARGEe)) DISCHARGEe = DISCHARGE
    if (associated(WGHTe)) WGHTe = WGHT

    call MAPL_GetResource( MAPL, ocean_dir, label='OCEAN_DIR:', rc=status ) ; VERIFY_(STATUS)
    call str4c( iarr, TRIM(ocean_dir) )

    IM = size(DISCHARGE,1)
    JM = size(DISCHARGE,2)
    allocate(HFLX(IM,JM), STAT=status)
    allocate(QFLX(IM,JM), STAT=status)
    VERIFY_(STATUS)

!US MIT gets net upward heat flux without short-wave radiation
    HFLX=-LWFLX+SHFLX+Av*QFLUX
!US Net fresh water flux, downward positive (flip sign inside the MIT driver)
    QFLX=RAIN+SNOW-QFLUX
!US Net short-wave raditaion, downward positive (flip sign inside the MIT driver)
    SWFLX = PENUVR+PENPAR+PENUVF+PENPAF+DRNIR+DFNIR

!ALT protect agaist "orphan" points
!    where (FRACICE == MAPL_Undef) FRACICE=0.0
!    where (VOLICE == MAPL_Undef) VOLICE=0.0
!    where (VOLSNO == MAPL_Undef) VOLSNO=0.0
!    where (ERGICE == MAPL_Undef) ERGICE=0.0
!    where (ERGSNO == MAPL_Undef) ERGSNO=0.0
!    where (TAUAGE == MAPL_Undef) TAUAGE=0.0
!    where (MPOND == MAPL_Undef) MPOND=0.0
!    where (TI == MAPL_Undef) TI=MAPL_TICE
!    where (SI == MAPL_Undef) SI=30.0
!    where (HI == MAPL_Undef) HI=0.0

    DO J=1,JM
       DO I=1,IM
          if (WGHT(I,J) == 0.0 .and. MASK(I,J,1) /= 0.0) then
             ! WGHT is 0, either because this truely is not a ocean point
             ! or GEOS does not think this is ocean point. In the latter,
             ! the values passed from OGCM are set to MAPL_Undef, 
             ! and we need to protect them. If MITgcm does not think this is
             ! ocean point, the protection is not needed but does not hurt
             ! we are going to be sloppy and overwrite the imports

             ! A related, but somewhat separate question is should be
             ! scale any of these variables by WGHT. If yes, we need to 
             ! un-scale them on the way back

             DISCHARGE(I,J) = 0.0
             TI(I,J,:) = MAPL_TICE
             HI(I,J) = 0.0
             SI(I,J) = 30.0
             FRACICE(I,J,:) = 0.0
             VOLICE(I,J,:) = 0.0
             VOLSNO(I,J,:) = 0.0
             ERGICE(I,J,:) = 0.0
             ERGSNO(I,J,:) = 0.0
             TAUAGE(I,J,:) = 0.0
             MPOND(I,J,:) = 0.0
          end if
       END DO
    END DO
!@@    FRACICE = min(FRACICE, 1.0)
! Put import data into internal state
!------------------------------------
    CALL DRIVER_SET_IMPORT_STATE( istate,   'TAUX',   TAUX )
    CALL DRIVER_SET_IMPORT_STATE( istate,   'TAUY',   TAUY )
    CALL DRIVER_SET_IMPORT_STATE( istate,   'PS',     PS )
    CALL DRIVER_SET_IMPORT_STATE( istate,   'SWHEAT', SWFLX )
    CALL DRIVER_SET_IMPORT_STATE( istate,   'HFLX',   HFLX )
    deallocate(HFLX)
    CALL DRIVER_SET_IMPORT_STATE( istate,   'DISCHARGE',   DISCHARGE )
    CALL DRIVER_SET_IMPORT_STATE( istate,   'QFLX',   QFLX )
    deallocate(QFLX)
    CALL DRIVER_SET_IMPORT_STATE( istate,   'SFLX',   SFLX )
    CALL DRIVER_SET_IMPORT_STATE( istate,   'LATS',   LATS )
    CALL DRIVER_SET_IMPORT_STATE( istate,   'LONS',   LONS )
    CALL DRIVER_SET_IMPORT_STATE( istate,   'WGHT',   WGHT )

    CALL DRIVER_SET_IMPORT_STATE( istate,   'FRACICE', FRACICE )
    CALL DRIVER_SET_IMPORT_STATE( istate,   'VOLICE',  VOLICE )
    CALL DRIVER_SET_IMPORT_STATE( istate,   'VOLSNO',  VOLSNO )
    CALL DRIVER_SET_IMPORT_STATE( istate,   'ERGICE',  ERGICE )
    CALL DRIVER_SET_IMPORT_STATE( istate,   'ERGSNO',  ERGSNO )
    CALL DRIVER_SET_IMPORT_STATE( istate,   'MPOND',   MPOND )
    CALL DRIVER_SET_IMPORT_STATE( istate,   'TAUAGE',  TAUAGE )
    CALL DRIVER_SET_IMPORT_STATE( istate,   'TI',  TI )
    CALL DRIVER_SET_IMPORT_STATE( istate,   'SI',  SI )
    CALL DRIVER_SET_IMPORT_STATE( istate,   'HI',  HI )
    call MAPL_GetResource( MAPL, active_ocean, label='ACTIVE_OCEAN:', &
         default=1, rc=status ) ; VERIFY_(STATUS)

    call mysetdir(iarr)
    if (active_ocean /= 0) CALL DRIVER_RUN( istate, 1 )
    deallocate(iarr)
    call popdir

    CALL MAPL_GetPointer(EXPORT,   UW,   'UW', RC=STATUS); VERIFY_(STATUS)
    CALL MAPL_GetPointer(EXPORT,   VW,   'VW', RC=STATUS); VERIFY_(STATUS)
    CALL MAPL_GetPointer(EXPORT,   TW,   'TW', RC=STATUS); VERIFY_(STATUS)
    CALL MAPL_GetPointer(EXPORT,   SW,   'SW', RC=STATUS); VERIFY_(STATUS)

    CALL MAPL_GetPointer(exportSI, DEL_FRACICE,'DEL_FRACICE', alloc=.true., __RC__)
    CALL MAPL_GetPointer(exportSI, DEL_TI,   'DEL_TI', alloc=.true., __RC__)
    CALL MAPL_GetPointer(exportSI, DEL_SI,   'DEL_SI', alloc=.true., __RC__)
    CALL MAPL_GetPointer(exportSI, DEL_VOLICE, 'DEL_VOLICE', alloc=.true., __RC__)
    CALL MAPL_GetPointer(exportSI, DEL_VOLSNO, 'DEL_VOLSNO', alloc=.true., __RC__)
    CALL MAPL_GetPointer(exportSI, DEL_ERGICE, 'DEL_ERGICE', alloc=.true., __RC__)
    CALL MAPL_GetPointer(exportSI, DEL_ERGSNO, 'DEL_ERGSNO', alloc=.true., __RC__)
    CALL MAPL_GetPointer(exportSI, DEL_MPOND, 'DEL_MPOND', alloc=.true., __RC__)
    CALL MAPL_GetPointer(exportSI, DEL_TAUAGE, 'DEL_TAUAGE', alloc=.true., __RC__)
    CALL MAPL_GetPointer(exportSI, DEL_HI, 'DEL_HI', alloc=.true., __RC__)

    DEL_FRACICE = 0.0
    DEL_TI = 0.0
    DEL_SI = 0.0
    DEL_VOLICE = 0.0
    DEL_VOLSNO = 0.0
    DEL_ERGICE = 0.0
    DEL_ERGSNO = 0.0
    DEL_MPOND = 0.0
    DEL_TAUAGE = 0.0
    DEL_HI = 0.0

    CALL DRIVER_GET_EXPORT_STATE( istate,   'US',   UW )
    CALL DRIVER_GET_EXPORT_STATE( istate,   'VS',   VW )
    CALL DRIVER_GET_EXPORT_STATE( istate,   'TS',   TW )
    CALL DRIVER_GET_EXPORT_STATE( istate,   'SS',   SW )
    CALL DRIVER_GET_EXPORT_STATE( istate, 'MASK', MASK )

    CALL DRIVER_GET_EXPORT_STATE( istate,'DELFRACICE', DEL_FRACICE )
    CALL DRIVER_GET_EXPORT_STATE( istate,'DELTI', DEL_TI )
    CALL DRIVER_GET_EXPORT_STATE( istate,'DELSI', DEL_SI )
    CALL DRIVER_GET_EXPORT_STATE( istate,'DELVOLICE', DEL_VOLICE )
    CALL DRIVER_GET_EXPORT_STATE( istate,'DELVOLSNO', DEL_VOLSNO )
    CALL DRIVER_GET_EXPORT_STATE( istate,'DELERGICE', DEL_ERGICE )
    CALL DRIVER_GET_EXPORT_STATE( istate,'DELERGSNO', DEL_ERGSNO )
    CALL DRIVER_GET_EXPORT_STATE( istate,'DELMPOND', DEL_MPOND )
    CALL DRIVER_GET_EXPORT_STATE( istate,'DELTAUAGE', DEL_TAUAGE )
    CALL DRIVER_GET_EXPORT_STATE( istate,'DELHI', DEL_HI )

    ! ALT: for now, we need to implement ridging
    ! and/or make sure it stays between 0 and 1

    ! ALT: for now we leave FRACICE alone and pass the increment to SALT

    if (associated(FRACICEe)) then
       FRACICEe = FRACICE
    end if

    ! Update the sea-ice fields
 
    DO J = 1, JM
       DO I = 1, IM
          ! US : We decided not to advect seaice in regions which are not 100% ocean
          ! both in GEOS and MITgcm. 0.99 was chosen arbitrarily because ==1 was found 
          ! problematic - many points in the Arctic were ~1e-3 smaller than one
          ! in the c90-llc90_02a experiment and this created many holes of seaice 
          ! advection in there.
          ! In case we would like to revisit this number, advlim in 
          ! import_state_fill_mod.FOR need to match this number.
          if (WGHT(I,J) > 0.99) then

             ! Apply increments
             FRACICE(I,J,:) = FRACICE(I,J,:) + DEL_FRACICE(I,J,:)
             VOLICE(I,J,:) = VOLICE(I,J,:) + DEL_VOLICE(I,J,:)
             VOLSNO(I,J,:) = VOLSNO(I,J,:) + DEL_VOLSNO(I,J,:)
             ERGICE(I,J,:) = ERGICE(I,J,:) + DEL_ERGICE(I,J,:)
             ERGSNO(I,J,:) = ERGSNO(I,J,:) + DEL_ERGSNO(I,J,:)
             TAUAGE(I,J,:) = TAUAGE(I,J,:) + DEL_TAUAGE(I,J,:)
             MPOND(I,J,:) = MPOND(I,J,:) + DEL_MPOND(I,J,:)
             !ALT we do not update skin, and the line below is commented out
             !HI(I,J) = HI(I,J) + DEL_HI(I,J)

             ! Apply a cutoff as additional check to MITgcm regularization
             DO C = 1, NUM_ICE_CATEGORIES
                IF (VOLICE(I,J,C) < cutoff .OR. FRACICE(I,J,C) < cutoff) THEN
                   FRACICE(I,J,C) = 0.0
                   VOLICE(I,J,C) = 0.0
                   TAUAGE(I,J,C) = 0.0
                   MPOND(I,J,C) = 0.0
                   DO L = 1, NUM_ICE_LAYERS
                      LCI = L + (C-1)*NUM_ICE_LAYERS
                      ERGICE(I,J,LCI) = 0.0
                   END DO
                END IF
                IF (VOLSNO(I,J,C) < cutoff .OR. FRACICE(I,J,C) < cutoff) THEN
                   VOLSNO(I,J,C) = 0.0
                   DO L = 1, NUM_SNOW_LAYERS
                      LCS = L + (C-1)*NUM_SNOW_LAYERS
                      ERGSNO(I,J,LCS) = 0.0
                   END DO
                END IF
             END DO
             
             !ALT: workaround to deal with a possible bug in DEL_TI
             where(FRACICE(I,J,:) /= 0.0 )  TI(I,J,:) = TI(I,J,:) + DEL_TI(I,J,:)
             if (any(DEL_VOLICE(I,J,:) /= 0.0)) then ! we advected at least one
                SI(I,J) = SI(I,J) + DEL_SI(I,J)
             end if

             ! Apply ridging algorithm
             if (sum(FRACICE(I,J,:)) > 1) then
                FRACICE(I,J,:) = FRACICE(I,J,:)/sum(FRACICE(I,J,:))
             end if
             
           end if
       END DO
    END DO
  
    call MAPL_TimerOff(MAPL,"RUN"   )
    call MAPL_TimerOff(MAPL,"TOTAL" )

! All Done
!---------
    call WRITE_PARALLEL( ' Finished plug run method ' )
    RETURN_(ESMF_SUCCESS)
  end subroutine Run

!=================================================================================

!BOP

! !IROUTINE: Run2  -- Run2 method, needed only when in dual_ocean mode. Apply correction to top-level MOM temperature, based on DEL_TEMP

! !INTERFACE:

  subroutine Run2  ( gc, import, export, clock, rc )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: gc     ! Gridded component 
    type(ESMF_State),    intent(INOUT) :: import ! Import state
    type(ESMF_State),    intent(INOUT) :: export ! Export state
    type(ESMF_Clock),    intent(INOUT) :: clock  ! The supervisor clock
    integer, optional,   intent(  OUT) :: rc     ! Error code:
!    type(ESMF_State)                   :: INTERNAL ! Internal state

!EOP

! ErrLog Variables

    character(len=ESMF_MAXSTR)		   :: IAm
    integer				   :: STATUS
    character(len=ESMF_MAXSTR)             :: COMP_NAME

! Locals

    integer                                :: counts(7)
    integer                                :: IM, JM

! Imports
    REAL_, pointer                         :: DEL_TEMP(:,:)

! Temporaries

    REAL_, pointer                      :: T(:,:)

! Pointers to export    
    REAL_, pointer                         :: MASK(:,:,:)

    type(MAPL_MetaComp),           pointer :: MAPL 
    TYPE(T_PrivateState_Wrap) wrap
    type(T_PrivateState), pointer :: privateState
    type(MITGCM_ISTATE),  pointer :: istate

    type(ESMF_Grid)                        :: Grid


! Begin
!------


! Get the component's name and set-up traceback handle.
! -----------------------------------------------------
    Iam = "Run2"
    call ESMF_GridCompGet( gc, NAME=comp_name, RC=status )
    VERIFY_(status)
    Iam = trim(comp_name) // Iam

! Get my internal MAPL_Generic state
!-----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=status)
    VERIFY_(STATUS)


! Profilers
!----------

    call MAPL_TimerOn (MAPL,"TOTAL")
    call MAPL_TimerOn (MAPL,"RUN2"  )

! Get the Plug's private internal state
!--------------------------------------

    call ESMF_UserCompGetInternalState( GC, trim(comp_name)//'_internal_state',&
         WRAP, STATUS )
    VERIFY_(status)

    privateState => WRAP%PTR
    istate => PrivateState%iState

! Get the grid, configuration
!----------------------------

    call ESMF_GridCompGet( GC, grid=Grid,  RC=status )
    VERIFY_(STATUS)

! Get IMPORT pointers
!--------------------

    call MAPL_GetPointer(IMPORT, DEL_TEMP, 'DEL_TEMP', RC=STATUS); VERIFY_(STATUS)

! Get EXPORT pointers
!--------------------
    ! by now this should be allocated, so 'alloc=.true.' is not needed
    CALL MAPL_GetPointer(EXPORT, MASK, trim(COMP_NAME)//'_3D_MASK', RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_GridGet(GRID, localCellCountPerDim=counts, RC=status)
    IM=counts(1)
    JM=counts(2)

! Temporaries with MOM default reals
!-----------------------------------

    call MAPL_GetPointer(EXPORT, T,  'TW', RC=STATUS)
    VERIFY_(STATUS)

    where(MASK(:,:,1) > 0.0) ! correct only ocean points
       !ALT: Note that we modify only top level of T
       !     we do not need to worry about temperature units
       !     since we are applying difference

       !     some relaxation ??? here or in guest ???

       T = T + DEL_TEMP

    end where
    CALL DRIVER_SET_IMPORT_STATE( istate,   'TS', T )
    deallocate(T)

! All Done
!---------

    RETURN_(ESMF_SUCCESS)
  end subroutine Run2

!BOP
    
! !IROUTINE: Finalize        -- Finalize method for GuestOcean wrapper

! !INTERFACE:

  subroutine Finalize ( gc, import, export, clock, rc )

! !ARGUMENTS:

  type(ESMF_GridComp), intent(INOUT) :: gc     ! Gridded component 
  type(ESMF_State),    intent(INOUT) :: import ! Import state
  type(ESMF_State),    intent(INOUT) :: export ! Export state
  type(ESMF_Clock),    intent(INOUT) :: clock  ! The supervisor clock
  integer, optional,   intent(  OUT) :: rc     ! Error code:

!EOP

    type (MAPL_MetaComp), pointer:: MAPL 
    type(ESMF_Time)                  :: MyTime

! ErrLog Variables

    character(len=ESMF_MAXSTR)	     :: IAm
    integer			     :: STATUS
    character(len=ESMF_MAXSTR)       :: COMP_NAME
   

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

! Generic Finalize
! ------------------
    
    call MAPL_GenericFinalize( GC, IMPORT, EXPORT, CLOCK, RC=status )
    VERIFY_(STATUS)

! All Done
!---------

    RETURN_(ESMF_SUCCESS)
  end subroutine Finalize



!UDI do we need record?
!====================================================================

! !IROUTINE: Record        -- Record method for GuestOcean wrapper (write intermediate restarts)

! !INTERFACE:

  subroutine Record ( gc, import, export, clock, rc )

! !ARGUMENTS:

  type(ESMF_GridComp), intent(INOUT) :: gc     ! Gridded component 
  type(ESMF_State),    intent(INOUT) :: import ! Import state
  type(ESMF_State),    intent(INOUT) :: export ! Export state
  type(ESMF_Clock),    intent(INOUT) :: clock  ! The supervisor clock
  integer, optional,   intent(  OUT) :: rc     ! Error code:

!EOP

    type (MAPL_MetaComp), pointer    :: MAPL 

! ErrLog Variables

    character(len=ESMF_MAXSTR)	     :: IAm
    integer			     :: STATUS
    character(len=ESMF_MAXSTR)       :: COMP_NAME

! Locals
    logical                          :: doRecord

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    Iam = "Record"
    call ESMF_GridCompGet( gc, NAME=comp_name, RC=status )
    VERIFY_(STATUS)
    Iam = trim(comp_name) // Iam

! Get my internal MAPL_Generic state
!-----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=status)
    VERIFY_(STATUS)

! Profilers
!----------

    call MAPL_TimerOn(MAPL,"TOTAL")

    doRecord = MAPL_RecordAlarmIsRinging(MAPL, RC=status)
    VERIFY_(STATUS)

    call MAPL_TimerOff(MAPL,"TOTAL")
    RETURN_(ESMF_SUCCESS)

  end subroutine Record

!====================================================================

end module MIT_GEOS5PlugMod
  
        

  





