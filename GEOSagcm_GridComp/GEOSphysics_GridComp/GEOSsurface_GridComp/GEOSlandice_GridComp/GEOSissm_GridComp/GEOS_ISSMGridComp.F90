!  $Id$

#include "MAPL_Generic.h"

module GEOS_IssmGridCompMod

!BOP
! !MODULE: GEOS_ISSM --- Calls ISSM (Ice-sheet and Sea-level System Model)

! !DESCRIPTION:
!
!   {\tt GEOS\_ISSM} is a wrapper that calls ISSM (C++) IRF methods
!
! *** NOTES: 
!            (*) currently we run over the Greenland Ice Sheet with all given PETs
!            (*) ISSM mesh is *internal* to ISSM (C++ source) but we create an ESMF_MESH version for regridding     
!            (*) the atmospheric grid is inherited from the parent, so we regrid onto the grid, then regrid onto tiles  
!            (*) for future development, it's important to note that ISSM expects double precision inputs    
!
! next steps:  
!            (*) determine "good" ISSM boundary and initial conditions
!            (*) run basic experiments
!            (*) work on basic coupling w.r.t. dynamic ice-surface elevation   
!            (*) expand to other fields (e.g., iceberg discharge), Antarctica, possibly other glaciers

! !USES:
use iso_fortran_env, only: dp=>real64
use iso_c_binding, only: c_ptr, c_double, c_f_pointer,c_null_char, c_loc, c_int
use ESMF
use MAPL
use GEOS_UtilsMod

implicit none

! declare interface to the ISSM C++ functions
interface
subroutine InitializeISSM(argc, argv, num_elements, num_nodes, comm) bind(C, NAME="InitializeISSM")
    import :: c_ptr, c_int
    integer(c_int), value        :: argc
    type(c_ptr), dimension(argc) :: argv
    integer(c_int)               :: num_elements
    integer(c_int)               :: num_nodes
    integer(c_int)               :: comm
end subroutine InitializeISSM
    
subroutine RunISSM(dt, gcmf, issmouts) bind(C,NAME="RunISSM")
   import :: c_ptr, c_double
   real(c_double),   value   :: dt
   type(c_ptr),      value   :: gcmf
   type(c_ptr),      value   :: issmouts
end subroutine RunISSM

subroutine GetNodesISSM(nodeIds, nodeCoords) bind(C,NAME="GetNodesISSM")
   import :: c_ptr
   type(c_ptr),      value   :: nodeIds
   type(c_ptr),      value   :: nodeCoords ! node coordinates (longitude,latitude)
end subroutine GetNodesISSM

subroutine GetElementsISSM(elementIds, elementConn, elementCoords) bind(C,NAME="GetElementsISSM")
  import :: c_ptr
  type(c_ptr),      value   :: elementIds
  type(c_ptr),      value   :: elementConn
  type(c_ptr),      value   :: elementCoords
end subroutine GetElementsISSM

subroutine FinalizeISSM() bind(C,NAME="FinalizeISSM")
end subroutine FinalizeISSM

end interface


private


! declare any variables here?

public SetServices

! !DESCRIPTION:
! 
!   {\tt GEOS\_ISSM} gridcomp runs NASA's Ice-sheet and Sea-level System Model (ISSM)
!                    currently over the Greenland Ice Sheet. 

type CONNECT_REGRIDHANDLES
  private
  type(ESMF_RouteHandle) :: routehandle_m2g ! routehandle for regridding mesh to grid
  type(ESMF_RouteHandle) :: routehandle_g2m ! routehandle for regridding grid to mesh
end type CONNECT_REGRIDHANDLES

! Wrapper for extracting internal state
! -------------------------------------
type REGRIDHANDLES
  type (CONNECT_REGRIDHANDLES), pointer :: ptr
end type REGRIDHANDLES

contains

!BOP

! !IROUTINE: SetServices -- Sets ESMF services for this component

! !INTERFACE:

subroutine SetServices ( GC, RC )

    ! !ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    integer, optional                  :: RC  ! return code

    ! !DESCRIPTION: 
!                This version uses the MAPL\_GenericSetServices, which sets
!                the Initialize and Finalize services, as well as allocating
!   our instance of a generic state and putting it in the 
!   gridded component (GC). Here we only need to set the run method and
!   add the state variable specifications (also generic) to our instance
!   of the generic state. This is the way our true state variables get into
!   the ESMF\_State INTERNAL, which is in the MAPL\_MetaComp. The import
!   and internal variables are allocated and initialized by generic. 

!EOP

!=============================================================================

! ErrLog Variables

    character(len=ESMF_MAXSTR)              :: IAm
    integer                                 :: STATUS
    character(len=ESMF_MAXSTR)              :: COMP_NAME

!=============================================================================

    type(MAPL_MetaComp), pointer            :: MAPL

    ! Get my internal MAPL_Generic state

    ! Begin...

! Get my name and set-up traceback handle
! ---------------------------------------

    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // 'SetServices'

! Set the Run entry point
! -----------------------

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_INITIALIZE,   Initialize, RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,          Run,        RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_FINALIZE,     Finalize,   RC=STATUS)
    VERIFY_(STATUS)

!-----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)


    ! Set the state variable specs.
! -----------------------------

 !Export states:
    call MAPL_AddExportSpec(GC,                    &
         SHORT_NAME = 'ICEEL',                     &
         LONG_NAME  = 'ice_elevation_tiles',       &
         UNITS      = 'm',                         &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone,          &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                    &
        SHORT_NAME = 'ICESMB',                     &
        LONG_NAME  = 'ice_surface_mass_balance',   &
        UNITS      = 'kg m-2 s-1',                          &
        DIMS       = MAPL_DimsTileOnly,            &
        VLOCATION  = MAPL_VLocationNone,           &
        RC=STATUS  )
    VERIFY_(STATUS)

!  !Import states:
    call MAPL_AddImportSpec(GC,                    &
        SHORT_NAME  = 'ACCUM',                     &
        LONG_NAME   = 'net_ice_accumulation_rate', &
        UNITS       = 'kg m-2 s-1',                &
        DIMS        = MAPL_DimsTileOnly,           &
        VLOCATION   = MAPL_VLocationNone,          &
        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                    &
        SHORT_NAME  = 'RUNOFF',                    &
        LONG_NAME   = 'runoff_total_flux',         &
        UNITS       = 'kg m-2 s-1',                &
        DIMS        = MAPL_DimsTileOnly,           &
        VLOCATION   = MAPL_VLocationNone,          &
        RC=STATUS  ) 
    VERIFY_(STATUS)

! Set the Profiling timers
! ------------------------

    call MAPL_TimerAdd(GC,    name="RUN"   ,RC=STATUS)
    VERIFY_(STATUS)

   
! ----------------------------------
    call MAPL_GenericSetServices    ( GC, RC=STATUS)
    VERIFY_(STATUS)


    RETURN_(ESMF_SUCCESS)
  
  end subroutine SetServices

  ! ! INITIALIZE:
  
  subroutine Initialize ( GC, IMPORT, EXPORT, CLOCK, RC )
    type(ESMF_GridComp),     intent(INOUT) :: GC                  ! Gridded component 
    type(ESMF_State),        intent(INOUT) :: IMPORT              ! Import state
    type(ESMF_State),        intent(INOUT) :: EXPORT              ! Export state
    type(ESMF_Clock),        intent(INOUT) :: CLOCK               ! The clock
    integer, optional,       intent(  OUT) :: RC                  ! Error code
    type(MAPL_MetaComp), pointer           :: MAPL   

    ! ErrLog Variables
    character(len=ESMF_MAXSTR)		   :: IAm
    integer				                   :: STATUS
    character(len=ESMF_MAXSTR)       :: COMP_NAME

    ! virtual machine / mpi comm
    type(ESMF_VM)                  :: vm    
    integer(c_int)                 :: comm                        ! mpi comm to pass to ISSM
    
    ! mesh information
    type(ESMF_Mesh)                :: mesh                        ! ESMF_Mesh representation of ISSM mesh
    integer, allocatable           :: elementTypes(:)             ! element geometry type (triangles)
    integer(c_int)                 :: num_elements                ! number of elements on PET
    integer(c_int)                 :: num_nodes                   ! number of nodes on PET
    integer, pointer, dimension(:) :: elementIds    => null()     ! list of elements local to PET
    integer, pointer, dimension(:) :: elementConn   => null()     ! element connectivity (nodes indices)
    real(dp),pointer, dimension(:) :: elementCoords => null()     ! element centroids
    real(dp),pointer,dimension(:)  :: nodeCoords    => null()     ! node coordinates (longitude,latitude)
    integer, pointer, dimension(:) :: nodeIds       => null()     ! list of nodes local to PET

    ! regridding 
    type(ESMF_Grid)                      :: grid                  ! atmospheric grid
    type(ESMF_RouteHandle)               :: routehandle_m2g       ! routehandle for regridding mesh to grid
    type(ESMF_RouteHandle)               :: routehandle_g2m       ! routehandle for regridding grid to mesh
    type(ESMF_Field)                     :: meshField             ! field on mesh
    type(ESMF_Field)                     :: gridField             ! field on grid
    type(CONNECT_REGRIDHANDLES), pointer :: regrid_handles        ! store the routehandles for access during run
    type(REGRIDHANDLES)                  :: wrap                  ! wrapper for routehandle container

    ! command-line arguments to initialize ISSM 
    integer(c_int)                 :: argc                        ! command line count for ISSM init
    character(len=200), dimension(:), allocatable, target :: argv ! command line args for ISSM init
    type(c_ptr), dimension(:), allocatable :: argv_ptr            ! pointer for passing command line args
    integer :: i                                                  ! loop index


    ! Get the target components name and set-up traceback handle.
    ! -----------------------------------------------------------

    Iam = "Initialize"
    call ESMF_GridCompGet( GC, NAME=comp_name, RC=status )
    VERIFY_(STATUS)
    Iam = trim(comp_name) // trim(Iam)


    ! Get my internal MAPL_Generic state
    !-----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)

    ! generic initialize 
    call MAPL_GenericInitialize( GC, IMPORT, EXPORT, CLOCK, RC=STATUS )
    VERIFY_(STATUS)

    call ESMF_VMGetCurrent(vm, rc=STATUS)
    VERIFY_(STATUS)

    call ESMF_VMGet(vm,mpiCommunicator=comm,rc=STATUS)
    VERIFY_(STATUS)

    ! ****************************************************
    ! call ISSM initialize C++ code so we can set up mesh

    ! Manually set command line argc and argv to initialize ISSM 
    ! note: we will probably want to change the experiment directory argv(3)
    !       to the GEOS experiment directory(?), after copying the input file (GreenlandGEOS.bin)
    !       and toolkits file (GreenlandGEOS.toolkits), TBD... 
    argc = 4  
    allocate(argv(argc))
    argv(1) = "unused"//c_null_char ! this argument is not actually used for anything
    argv(2) = "TransientSolution"//c_null_char
    argv(3) = "/discover/nobackup/agstubbl/ISSM/GEOS-ISSM/ISSM/examples/GEOSInput"//c_null_char
    argv(4) = "GreenlandGEOS"//c_null_char

    ! Convert Fortran strings to C pointers (in argv_ptr)
    allocate(argv_ptr(argc))
    do i = 1, argc
        argv_ptr(i) = c_loc(argv(i))
    end do

    call ESMF_VMBarrier(vm, rc=status)
    
    ! Call the C++ function for initializing ISSM
    ! gets the number of elements and nodes of the mesh
    call InitializeISSM(argc, argv_ptr,num_elements,num_nodes,comm)

    !allocate mesh-related pointers
    allocate(nodeCoords(2*num_nodes))
    allocate(nodeIds(num_nodes))
    allocate(elementTypes(num_elements))
    allocate(elementIds(num_elements))
    allocate(elementConn(3*num_elements))
    allocate(elementCoords(2*num_elements))

    ! create ESMF mesh corresponding to  ISSM mesh 
    ! get information about nodes and elements
    ! node coords and element coords (centroids) are in (lon,lat)
    call GetNodesISSM(c_loc(nodeIds), c_loc(nodeCoords)) 
    call GetElementsISSM(c_loc(elementIds), c_loc(elementConn), c_loc(elementCoords))

    elementTypes(:) = ESMF_MESHELEMTYPE_TRI

    ! create the ESMF mesh from ISSM mesh properties
    mesh = ESMF_MeshCreate(parametricDim=2, spatialDim=2, nodeIds=nodeIds, nodeCoords=nodeCoords, &
            elementIds=elementIds, elementTypes=elementTypes, elementConn=elementConn,& 
            elementCoords=elementCoords,coordSys=ESMF_COORDSYS_SPH_DEG, rc=STATUS)

    VERIFY_(STATUS)

    ! associate ESMF_Mesh representation of ISSM mesh with GC for regridding imports/exports in run method       
    call ESMF_GridCompSet(GC,mesh=mesh,rc=STATUS)
    VERIFY_(STATUS)

    ! set up regridding next
    ! create source field on ISSM mesh
    meshField = ESMF_FieldCreate(mesh=mesh,typekind=ESMF_TYPEKIND_R8,meshloc=ESMF_MESHLOC_ELEMENT,rc=STATUS)
    VERIFY_(STATUS)

    ! get atmospheric grid 
    call ESMF_GridCompGet( GC, GRID=grid, RC=status ); VERIFY_(STATUS)
    
    ! create destination field on atmospheric grid
    gridField = ESMF_FieldCreate(grid=grid,typekind=ESMF_TYPEKIND_R4,rc=STATUS); VERIFY_(STATUS)
    
    ! create routehandle for mesh-to-grid regridding
    call ESMF_FieldRegridStore(srcField=meshField, dstField=gridField,routehandle=routehandle_m2g,unmappedaction=ESMF_UNMAPPEDACTION_IGNORE,rc=STATUS)
    VERIFY_(STATUS)

    ! create routehandle for grid-to-mesh regridding
    call ESMF_FieldRegridStore(srcField=gridField, dstField=meshField,routehandle=routehandle_g2m,unmappedaction=ESMF_UNMAPPEDACTION_IGNORE,rc=STATUS)
    VERIFY_(STATUS)
    
    ! create pointer to routehandle for component's private internal state
    allocate(regrid_handles, stat=status)
    VERIFY_(STATUS)
    regrid_handles%routehandle_m2g  = routehandle_m2g
    regrid_handles%routehandle_g2m  = routehandle_g2m 
    wrap%ptr => regrid_handles
    call ESMF_UserCompSetInternalState ( GC, 'REGRIDHANDLES', wrap, status )
    VERIFY_(STATUS)

    ! deallocate pointers
    deallocate(argv)
    deallocate(argv_ptr)
    deallocate(nodeCoords)
    deallocate(nodeIds)
    deallocate(elementTypes)
    deallocate(elementIds)
    deallocate(elementConn)
    deallocate(elementCoords)

    RETURN_(ESMF_SUCCESS)
  end subroutine Initialize

  !BOP


subroutine RUN ( GC, IMPORT, EXPORT, CLOCK, RC )
! ! Run ISSM ice-sheet model (only over Greenland right now)
! !ARGUMENTS:
  type(ESMF_GridComp), intent(inout)   :: GC                          ! Gridded component 
  type(ESMF_State),    intent(inout)   :: IMPORT                      ! Import state
  type(ESMF_State),    intent(inout)   :: EXPORT                      ! Export state
  type(ESMF_Clock),    intent(inout)   :: CLOCK                       ! The clock
  integer, optional,   intent(  out)   :: RC                          ! Error code
  type (ESMF_Alarm)                    :: ALARM                       ! run alarm for ISSM component 
  
  ! ErrLog Variables
  character(len=ESMF_MAXSTR)           :: IAm
  integer                              :: STATUS
  character(len=ESMF_MAXSTR)           :: COMP_NAME

  ! regridding
  type(ESMF_RouteHandle)               :: routehandle_m2g             ! routehandle for regridding mesh to grid
  type(ESMF_RouteHandle)               :: routehandle_g2m             ! routehandle for regridding grid to mesh
  type(ESMF_Field)                     :: srcField                    ! generic source field for regridding
  type(ESMF_Field)                     :: dstField                    ! generic destination field for regridding
  type(ESMF_Grid)                      :: grid                        ! atmospheric grid
  type(CONNECT_REGRIDHANDLES), pointer :: regrid_handles=>null()      ! store the routehandles for access during run
  type(REGRIDHANDLES)                  :: wrap                        ! wrapper for routehandle container

  ! time stepping information (ISSM_DT set in AGCM.rc)
  real(dp)            :: dt                                           ! time step in seconds
  real(dp)            :: dt_yr                                        ! time step in years (ISSM units)
  real(dp), parameter :: sec_per_year = 31536000.0                    ! conversion factor (ISSM value for consistency)

  type(MAPL_MetaComp), pointer   :: MAPL
  type(ESMF_Mesh)                :: mesh                              ! ESMF version of ISSM mesh
  integer(c_int)                 :: num_elements                      ! number of elements on each PET
  type(ESMF_VM)                  :: vm  

  ! tile information
  type (MAPL_LocStream       )        :: locstream                    ! tile locstream 
  integer, pointer, dimension(:)      :: TILETYPES                    ! types of tiles
  integer                             :: NT                           ! number of tiles

  ! ice elevation on mesh, grid, tile
  real(dp),    pointer, dimension(:)  :: ICEEL_MESH    => null()      ! ice-sheet elevation on mesh elements
  real, pointer, dimension(:,:)       :: ICEEL_GRID    => null()      ! ice-sheet elevation on atmospheric grid 
  real, pointer, dimension(:)         :: ICEEL_TILE(:) => null()      ! ice-sheet elevation on landice tiles
  real, pointer, dimension(:)         :: ICEEL         => null()      ! pointer to ice-sheet elevation export state

  ! need to do the same for SMB
  real(dp),    pointer, dimension(:)  :: ICESMB_MESH   => null()      ! surface mass balce on mesh elements
  real,    pointer, dimension(:,:)    :: ICESMB_GRID   => null()      ! surface mass balance on atmospheric grid
  real, pointer, dimension(:)         :: ICESMB_TILE(:)=> null()      ! surface mass balance on tiles
  real, pointer, dimension(:)         :: ICESMB        => null()      ! pointer to SMB export state (not currently exported by landice)
  real, pointer, dimension(:)         :: ACCUM         => null()      ! accumulation import
  real, pointer, dimension(:)         :: RUNOFF        => null()      ! runoff import

  ! physical parameters
  real(dp), parameter :: rho_ice   = 917.0                            ! pure ice density [kg/m^3]


! Get the target components name, mesh and vm
! -----------------------------------------------------------
  Iam = "Run"
  call ESMF_GridCompGet( GC, name=COMP_NAME,mesh=mesh,vm=vm,RC=STATUS )
  VERIFY_(STATUS)
  Iam = trim(COMP_NAME) // Iam

  ! Get my internal MAPL_Generic state
!----------------------------------

  call MAPL_GetObjectFromGC(GC, MAPL, STATUS)
  VERIFY_(STATUS)

  ! Start Total timer
!------------------
  call MAPL_TimerOn(MAPL,"TOTAL")
  call MAPL_TimerOn(MAPL,"RUN" )

  call MAPL_Get(MAPL, RUNALARM = ALARM, RC=STATUS )
  VERIFY_(STATUS)

  ! run ISSM at specified time steps
  if ( ESMF_AlarmIsRinging (ALARM, RC=STATUS) ) then

    ! *************************************************************************** !
    ! BASIC SETUP
    ! *************************************************************************** !

    ! get timestep for ISSM
    call MAPL_GetResource ( MAPL, dt, Label=trim(COMP_NAME)//"_DT:", RC=STATUS)

    ! convert dt to years for ISSM input
    dt_yr = dt/sec_per_year

    ! get number of mesh elements
    call ESMF_MeshGet(mesh,elementCount=num_elements)

    ! allocate SMB forcing (input to ISSM) and ice-elevation output (export from ISSM)
    allocate(ICESMB_MESH(num_elements))
    allocate(ICEEL_MESH(num_elements))

    call ESMF_VMBarrier(vm, rc=status)
    VERIFY_(STATUS)

    ! set smb and surface to zero (not sure this is needed...)
    ICESMB_MESH(:) = 0.0_dp  
    ICEEL_MESH(:) = 0.0_dp

    ! get LocStream for landice tiles for regridding
    call MAPL_Get(MAPL, LocStream = locstream, TILETYPES = TILETYPES, RC=STATUS)
    VERIFY_(STATUS)
    NT = size(TILETYPES)

    ! get the atmospheric grid 
    call ESMF_GridCompGet( GC, GRID=grid, RC=status )
    VERIFY_(STATUS)
    
    ! get routehandles for regridding
    call ESMF_UserCompGetInternalState(GC, 'REGRIDHANDLES', wrap, status); VERIFY_(STATUS)
    regrid_handles => wrap%ptr
    routehandle_m2g = regrid_handles%routehandle_m2g  ! mesh to grid 
    routehandle_g2m = regrid_handles%routehandle_g2m  ! grid to mesh

    ! *************************************************************************** !
    ! CALCULATE SMB (surface mass balance) & REGRID FROM TILES [to grid] to MESH 
    ! *************************************************************************** !
    call MAPL_GetPointer(EXPORT,ICESMB, 'ICESMB'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT,ACCUM , 'ACCUM'    ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT,RUNOFF, 'RUNOFF'    , RC=STATUS); VERIFY_(STATUS)

    ! allocate tiles for SMB (MKTILE)
    if(associated(ICESMB) .and. .not.associated(ICESMB_TILE)) then
      allocate(ICESMB_TILE(NT), STAT=STATUS)
      VERIFY_(STATUS)
      ICESMB_TILE = MAPL_Undef
    end if

    ! calculate SMB as accumulation minus runoff (check this)
    ICESMB_TILE(:) = ACCUM(:) - RUNOFF(:)
    
    ! pointer to SMB export
    if(associated(ICESMB)) ICESMB(:) = ICESMB_TILE(:)

    ! transform from tile to grid
    call MAPL_LocStreamTransform( LOCSTREAM, ICESMB_GRID,ICESMB_TILE, RC=STATUS)
    VERIFY_(STATUS)

    ! create source field: SMB on grid
    srcField = ESMF_FieldCreate(grid=grid,farrayPtr=ICESMB_GRID, & 
               datacopyflag=ESMF_DATACOPY_VALUE,rc=STATUS)
    VERIFY_(STATUS)

    ! create destination field: SMB on mesh elements
    dstField = ESMF_FieldCreate(mesh=mesh,typekind=ESMF_TYPEKIND_R8,rc=STATUS)
    VERIFY_(STATUS)

    ! regrid SMB from grid to mesh
    call ESMF_FieldRegrid(srcField, dstField, routehandle_g2m, RC=STATUS); VERIFY_(STATUS)

    ! get pointer to SMB on mesh
    call ESMF_FieldGet(dstField,farrayPtr=ICESMB_MESH,RC=STATUS); VERIFY_(STATUS)

    ! convert to SMB to ISSM units [m/yr]
    ICESMB_MESH = ICESMB_MESH*sec_per_year/rho_ice

    ! not sure if these destroy calls are needed because these fields are reused below(?)
    call ESMF_FieldDestroy(srcField,rc=STATUS); VERIFY_(STATUS)
    call ESMF_FieldDestroy(dstField,rc=STATUS); VERIFY_(STATUS)

    ! *************************************************************************** !
    !  RUN ISSM WITH SMB INPUT AND ICE-ELEVATION OUTPUT
    ! *************************************************************************** !


    ! NOTE: do we need the barriers before/after ISSM run?
    call ESMF_VMBarrier(vm, rc=status); VERIFY_(STATUS)

    ! call run method from ISSM library 
    call RunISSM(dt_yr, c_loc(ICESMB_MESH), c_loc(ICEEL_MESH))

    ! *************************************************************************** !
    ! REGRID ICEEL (ice elevation) FROM MESH [to grid] to TILES 
    ! *************************************************************************** !

    ! create source field: ice elevation on mesh elements
    srcField = ESMF_FieldCreate(mesh=mesh,farrayPtr=ICEEL_MESH,meshloc=ESMF_MESHLOC_ELEMENT, & 
    datacopyflag=ESMF_DATACOPY_VALUE,rc=STATUS)
    VERIFY_(STATUS)
    
    ! create destination field: regrid ice elevation onto grid
    dstField = ESMF_FieldCreate(grid=grid,typekind=ESMF_TYPEKIND_R4,rc=STATUS)
    VERIFY_(STATUS)

    ! regrid ice elevation from mesh to grid
    call ESMF_FieldRegrid(srcField, dstField, routehandle_m2g, RC=STATUS); VERIFY_(STATUS)

    ! get pointer to ice elevation on grid
    call ESMF_FieldGet(dstField,farrayPtr=ICEEL_GRID,RC=STATUS); VERIFY_(STATUS)
    
    ! get pointer to export
    call MAPL_GetPointer(EXPORT  , ICEEL , 'ICEEL' ,  RC=STATUS); VERIFY_(STATUS)

    ! allocate tiles for ice elevation (MKTILE)
    if(associated(ICEEL) .and. .not.associated(ICEEL_TILE)) then
      allocate(ICEEL_TILE(NT), STAT=STATUS)
      VERIFY_(STATUS)
      ICEEL_TILE = MAPL_Undef
    end if
   
    ! transform and assign to export pointer
    call MAPL_LocStreamTransform( LOCSTREAM, ICEEL_TILE,ICEEL_GRID, RC=STATUS)
    VERIFY_(STATUS)

    if(associated(ICEEL)) ICEEL(:) = ICEEL_TILE(:)

    call ESMF_FieldDestroy(srcField,rc=STATUS); VERIFY_(STATUS)
    call ESMF_FieldDestroy(dstField,rc=STATUS); VERIFY_(STATUS)

  end if 

  call ESMF_VMBarrier(vm, rc=status)
  VERIFY_(STATUS)

  if(associated(ICESMB_MESH))  deallocate(ICESMB_MESH)
  if(associated(ICEEL_MESH)) deallocate(ICEEL_MESH)
  if(associated(ICEEL_TILE)) deallocate(ICEEL_TILE)

  call MAPL_TimerOff(MAPL,"RUN"  )
  call MAPL_TimerOff(MAPL,"TOTAL")
 
  RETURN_(ESMF_SUCCESS)
 
 end subroutine RUN

 !BOP
    
!IROUTINE: Finalize        -- Finalize method for issm 

!INTERFACE:

 subroutine Finalize ( GC, IMPORT, EXPORT, CLOCK, RC )

    ! !ARGUMENTS:
    
    type(ESMF_GridComp), intent(INOUT) :: GC     ! Gridded component 
    type(ESMF_State),    intent(INOUT) :: IMPORT ! Import state
    type(ESMF_State),    intent(INOUT) :: EXPORT ! Export state
    type(ESMF_Clock),    intent(INOUT) :: CLOCK  ! The supervisor clock
    integer, optional,   intent(  OUT) :: RC     ! Error code:
    
    !EOP
    type(MAPL_MetaComp), pointer       :: MAPL 
    
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

    ! call ISSM finalize (saves binary output .outbin file)
    call FinalizeISSM()

! Generic Finalize
! ------------------
    
    call MAPL_GenericFinalize( GC, IMPORT, EXPORT, CLOCK, RC=status )
    VERIFY_(STATUS)

! All Done
!---------

    RETURN_(ESMF_SUCCESS)
  end subroutine Finalize

end module GEOS_IssmGridCompMod