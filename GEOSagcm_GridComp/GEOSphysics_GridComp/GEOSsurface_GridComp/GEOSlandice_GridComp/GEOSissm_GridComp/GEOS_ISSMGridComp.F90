!  $Id$

#include "MAPL_Generic.h"

module GEOS_IssmGridCompMod

!BOP
! !MODULE: GEOS_ISSM --- Runs ISSM (Ice-sheet and Sea-level System Model)
! 
!
! !DESCRIPTION:
!
!   {\tt GEOS\_ISSM} is a wrapper that calls ISSM (C++) IRF methods
!   Imports: ICESMB
!   Exports: ICEEL
! *** NOTES: 
!            (*) currently we run over the Greenland Ice Sheet with all given PETs
!            (*) ISSM mesh is *internal* to ISSM (C++ source)--we create an ESMF_MESH version for regridding     
!            (*) the attached grid is inherited from parent. we regrid mesh to grid, then transform to tiles  
!            (*) for future development, note that ISSM expects double precision inputs    
!
! next steps:  
!            (*) determine "good" ISSM boundary conditions and initial conditions
!            (*) run basic experiments
!            (*) work on basic coupling w.r.t. dynamic ice-surface elevation   
!            (*) expand to other fields (e.g., ice discharge), Antarctica, possibly other glaciers

! !USES:
use iso_fortran_env, only: dp=>real64
use iso_c_binding, only: c_ptr, c_double, c_f_pointer,c_null_char, c_loc, c_int
use ESMF
use MAPL
use GEOS_UtilsMod

implicit none

! declare interface to the ISSM C++ library (arguments described in Initialize & Run below)
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
   real(c_double),   value       :: dt
   type(c_ptr),      value       :: gcmf
   type(c_ptr),      value       :: issmouts
end subroutine RunISSM

subroutine GetNodesISSM(nodeIds, nodeCoords) bind(C,NAME="GetNodesISSM")
   import :: c_ptr
   type(c_ptr),      value       :: nodeIds
   type(c_ptr),      value       :: nodeCoords !
end subroutine GetNodesISSM

subroutine GetElementsISSM(elementIds, elementConn, elementCoords) bind(C,NAME="GetElementsISSM")
  import :: c_ptr
  type(c_ptr),      value       :: elementIds
  type(c_ptr),      value       :: elementConn
  type(c_ptr),      value       :: elementCoords
end subroutine GetElementsISSM

subroutine FinalizeISSM() bind(C,NAME="FinalizeISSM")
end subroutine FinalizeISSM

end interface


private


public SetServices

! !DESCRIPTION:
! 
!   {\tt GEOS\_ISSM} gridcomp runs NASA's Ice-sheet and Sea-level System Model (ISSM)
!             

! private internal state for regridding 
type T_ISSM_STATE
  private
  type(ESMF_RouteHandle) :: routehandle_m2g ! routehandle for regridding mesh to grid
  type(ESMF_RouteHandle) :: routehandle_g2m ! routehandle for regridding grid to mesh
  type(ESMF_GRID)        :: grid           ! original grid
  type(MAPL_LocStream)   :: locstream   
end type T_ISSM_STATE

! Wrapper for extracting internal state
! -------------------------------------
type ISSM_WRAP
  type (T_ISSM_STATE), pointer :: ptr
end type ISSM_WRAP

contains

!BOP

! !IROUTINE: SetServices -- Sets ESMF services for this component

! !INTERFACE:

subroutine SetServices ( GC, RC )

    ! !ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    integer, optional                  :: RC  ! return code

    ! !DESCRIPTION: 
!   This version uses the MAPL\_GenericSetServices Here we set the initialize method,
!   run method, and finalize method because we are interfacing with the external ISSM
!   library IRF methods. We also add the state variable specifications (also generic) 
!   to our instance of the generic state. This is the way our true state variables get into
!   the ESMF\_State INTERNAL, which is in the MAPL\_MetaComp. The import
!   variables are allocated and initialized by generic. 

!EOP

!=============================================================================

! ErrLog Variables

    character(len=ESMF_MAXSTR)         :: IAm
    integer                            :: STATUS
    character(len=ESMF_MAXSTR)         :: COMP_NAME

!=============================================================================

    type(MAPL_MetaComp), pointer       :: MAPL

    type (ESMF_Config)                 :: CF

    real                               :: dt     ! time step [s] (ISSM_DT set in AGCM.rc)

    ! Get my internal MAPL_Generic state

    ! Begin...

! Get my name and set-up traceback handle
! ---------------------------------------

    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // 'SetServices'

! Set the Initialize, Run, and Finalize entry points
!-----------------------------------

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_INITIALIZE,   Initialize, RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,          Run,        RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_FINALIZE,     Finalize,   RC=STATUS)
    VERIFY_(STATUS)

!-----------------------------------

    call MAPL_GetObjectFromGC (GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)

    call ESMF_GridCompGet(GC, CONFIG = CF, RC=STATUS)
    VERIFY_(STATUS)

    ! get timestep for ISSM
    call ESMF_ConfigGetAttribute(CF, dt, Label=trim(COMP_NAME)//"_DT:", DEFAULT=604800.0, RC=STATUS)
    VERIFY_(STATUS)


    ! Set the state variable specs.
!-----------------------------------

!   Import states:
    call MAPL_AddImportSpec(GC,                    &
        SHORT_NAME = 'ICESMB',                     &
        LONG_NAME  = 'ice_surface_mass_balance',   &
        UNITS      = 'kg m-2 s-1',                 &
        DIMS       = MAPL_DimsTileOnly,            &
        VLOCATION  = MAPL_VLocationNone,           &
        AVERAGING_INTERVAL = nint(dt),             &
        REFRESH_INTERVAL   = nint(dt),             &
        RC=STATUS  )
    VERIFY_(STATUS)

!   Export states:
    call MAPL_AddExportSpec(GC,                    &
         SHORT_NAME = 'ICEEL',                     &
         LONG_NAME  = 'ice_sheet_elevation',       &
         UNITS      = 'm',                         &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone,          &
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
    type(ESMF_GridComp),     intent(INOUT) :: GC                      ! Gridded component 
    type(ESMF_State),        intent(INOUT) :: IMPORT                  ! Import state
    type(ESMF_State),        intent(INOUT) :: EXPORT                  ! Export state
    type(ESMF_Clock),        intent(INOUT) :: CLOCK                   ! The clock
    integer, optional,       intent(OUT)   :: RC                      ! Error code
    type(MAPL_MetaComp), pointer           :: MAPL   

    ! ErrLog Variables
    character(len=ESMF_MAXSTR)             :: IAm
    integer                                :: STATUS
    character(len=ESMF_MAXSTR)             :: COMP_NAME

    ! virtual machine / mpi comm
    type(ESMF_VM)                          :: vm    
    integer(c_int)                         :: comm                    ! mpi comm to pass to ISSM
    
    ! mesh information
    type(ESMF_Mesh)                        :: mesh                    ! ESMF_Mesh representation of ISSM mesh
    integer, allocatable                   :: elementTypes(:)         ! element geometry type (triangles)
    integer(c_int)                         :: num_elements            ! number of elements on PET
    integer(c_int)                         :: num_nodes               ! number of nodes on PET
    integer, pointer, dimension(:)         :: elementIds    => null() ! list of elements local to PET
    integer, pointer, dimension(:)         :: elementConn   => null() ! element connectivity (nodes indices)
    real(dp),pointer, dimension(:)         :: elementCoords => null() ! element centroids
    real(dp),pointer,dimension(:)          :: nodeCoords    => null() ! node coordinates (longitude,latitude)
    integer, pointer, dimension(:)         :: nodeIds       => null() ! list of nodes local to PET

    ! regridding 
    type(ESMF_Grid)                        :: grid                    ! atmospheric grid
    type(ESMF_RouteHandle)                 :: routehandle_m2g         ! routehandle for regridding mesh to grid
    type(ESMF_RouteHandle)                 :: routehandle_g2m         ! routehandle for regridding grid to mesh
    type(ESMF_Field)                       :: meshField               ! field on mesh
    type(ESMF_Field)                       :: gridField               ! field on grid
    type(T_ISSM_STATE), pointer            :: internal_state          ! store the routehandles for access during run
    type(ISSM_WRAP)                        :: wrap                    ! wrapper for routehandle container

    ! command-line arguments to initialize ISSM 
    character(len=ESMF_MAXSTR), dimension(:), allocatable, target :: argv 
    integer(c_int)                         :: argc                    ! command line count for ISSM init   
    type(c_ptr), dimension(:), allocatable :: argv_ptr                ! pointer for passing command line args
    integer                                :: i,j                     ! loop indices
    character(len=ESMF_MAXSTR)             :: ISSM_EXPDIR             ! directory containing ISSM input file
    character(len=ESMF_MAXSTR)             :: ISSM_EXPNAME            ! name of ISSM input file

    ! variables for saving mesh output
    integer                                :: ISSM_SAVEMESH           ! mesh save flag
    type(ESMF_Field)                       :: field                   ! for fieldwrites
    real(dp),    pointer, dimension(:)     :: field_saver => null()
    real(dp),    pointer, dimension(:)     :: nodelons    => null()
    real(dp),    pointer, dimension(:)     :: nodelats    => null()
    integer, pointer, dimension(:)         :: conn_slice  => null()
    character(len=16) :: varname
    character(len=1)  :: istr


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


    call ESMF_VMGetCurrent(vm, rc=STATUS)
    VERIFY_(STATUS)

    call ESMF_VMGet(vm,mpiCommunicator=comm,rc=STATUS)
    VERIFY_(STATUS)

    ! ****************************************************
    ! call ISSM initialize C++ code so we can set up mesh

    ! Manually set command line argc and argv to initialize ISSM 
    call MAPL_GetResource(MAPL, ISSM_EXPDIR, Label=trim(COMP_NAME)//"_EXPDIR:", RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_GetResource(MAPL, ISSM_EXPNAME, Label=trim(COMP_NAME)//"_EXPNAME:", RC=STATUS)
    VERIFY_(STATUS)

    argc = 4  
    allocate(argv(argc))
    argv(1) = "unused"//c_null_char ! executable path: this argument is not used for library calls
    argv(2) = "TransientSolution"//c_null_char
    argv(3) = trim(ISSM_EXPDIR)//c_null_char
    argv(4) = trim(ISSM_EXPNAME)//c_null_char

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

    ! get information about nodes and elements
    ! node coords and element coords (centroids) are in (lon,lat)
    call GetNodesISSM(c_loc(nodeIds), c_loc(nodeCoords)) 
    call GetElementsISSM(c_loc(elementIds), c_loc(elementConn), c_loc(elementCoords))

    elementTypes(:) = ESMF_MESHELEMTYPE_TRI ! triangular elements

    ! create the ESMF mesh from ISSM mesh properties
    mesh = ESMF_MeshCreate(parametricDim=2, spatialDim=2, nodeIds=nodeIds, nodeCoords=nodeCoords, &
            elementIds=elementIds, elementTypes=elementTypes, elementConn=elementConn,& 
            elementCoords=elementCoords,coordSys=ESMF_COORDSYS_SPH_DEG, rc=STATUS)

    VERIFY_(STATUS)

    ! associate ESMF_Mesh representation of ISSM mesh with GC for regridding imports/exports in run method       
    call ESMF_GridCompSet(GC,mesh=mesh,rc=STATUS)
    VERIFY_(STATUS)

    ! Set up regridding next
    !-----------------------------------
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
    allocate(internal_state, stat=status)
    VERIFY_(STATUS)
    internal_state%routehandle_m2g  = routehandle_m2g
    internal_state%routehandle_g2m  = routehandle_g2m 
    internal_state%grid             = grid
    call MAPL_Get(MAPL, LocStream = internal_state%locstream, _RC)

    wrap%ptr => internal_state
    call ESMF_UserCompSetInternalState ( GC, 'ISSM_WRAP', wrap, status )
    VERIFY_(STATUS)

    ! save mesh output if desired
    call MAPL_GetResource(MAPL, ISSM_SAVEMESH, Label=trim(COMP_NAME)//"_SAVEMESH:",default=0,RC=STATUS)
    VERIFY_(STATUS)
    if (ISSM_SAVEMESH/=0) then
        allocate(field_saver(num_elements))
        allocate(nodelons(num_nodes))
        allocate(nodelats(num_nodes))
        allocate(conn_slice(num_elements))

        nodelons = nodeCoords(1::2)
        nodelats = nodeCoords(2::2)
    
        ! initialize ICEEL and ICESMB mesh netcdf files to have time slices
        call ESMF_FieldWrite(meshField, trim(ISSM_EXPDIR)//"/iceel.nc", variableName='ICEEL',timeslice=1,rc=STATUS)
        VERIFY_(STATUS)
        call ESMF_FieldWrite(meshField, trim(ISSM_EXPDIR)//"/icesmb.nc", variableName='ICESMB',timeslice=1,rc=STATUS)
        VERIFY_(STATUS)
        
        ! save nodes for reconstructing mesh output after running
        ! looping over three nodes per triangle
        do i = 1, 3
            conn_slice = [( elementConn(j), j=i, size(elementConn), 3 )]
        
            write(istr, '(I1)') i
        
            ! ---- longitude ----
            varname = "node" // istr // "lon"
            field_saver = nodelons(conn_slice)
            field = ESMF_FieldCreate(mesh=mesh,farrayPtr=field_saver,meshloc=ESMF_MESHLOC_ELEMENT, rc=STATUS)
            VERIFY_(STATUS)
            call ESMF_FieldWrite(field, trim(ISSM_EXPDIR)//"/mesh.nc", variableName=varname, rc=STATUS)
            VERIFY_(STATUS)
            call ESMF_FieldDestroy(field, rc=STATUS)
            VERIFY_(STATUS)
        
            ! ---- latitude ----
            varname = "node" // istr // "lat"
            field_saver = nodelats(conn_slice)
            field = ESMF_FieldCreate(mesh=mesh, farrayPtr=field_saver,meshloc=ESMF_MESHLOC_ELEMENT, rc=STATUS)
            VERIFY_(STATUS)
            call ESMF_FieldWrite(field, trim(ISSM_EXPDIR)//"/mesh.nc", variableName=varname, rc=STATUS)
            VERIFY_(STATUS)
            call ESMF_FieldDestroy(field, rc=STATUS)
            VERIFY_(STATUS)
        
        end do

        ! ---- element IDs ----
        varname = "elementIds"
        field_saver = elementIds(:)
        field = ESMF_FieldCreate(mesh=mesh, farrayPtr=field_saver,meshloc=ESMF_MESHLOC_ELEMENT, rc=STATUS)
        VERIFY_(STATUS)
        call ESMF_FieldWrite(field, trim(ISSM_EXPDIR)//"/mesh.nc", variableName=varname, rc=STATUS)
        VERIFY_(STATUS)
        call ESMF_FieldDestroy(field, rc=STATUS)
        VERIFY_(STATUS)

        deallocate(field_saver)
        deallocate(nodelons)
        deallocate(nodelats)
        deallocate(conn_slice)
    end if 

    ! deallocate pointers
    deallocate(argv)
    deallocate(argv_ptr)
    deallocate(nodeCoords)
    deallocate(nodeIds)
    deallocate(elementTypes)
    deallocate(elementIds)
    deallocate(elementConn)
    deallocate(elementCoords)

    ! Create losctream that match mesh element id, then set it to this GC and MAPL
    

    ! generic initialize 
    call MAPL_GenericInitialize( GC, IMPORT, EXPORT, CLOCK, RC=STATUS )
    VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)
  end subroutine Initialize

  !BOP


subroutine RUN ( GC, IMPORT, EXPORT, CLOCK, RC )
! ! Run ISSM ice-sheet model 
! !ARGUMENTS:
  type(ESMF_GridComp), intent(inout)   :: GC                      ! Gridded component 
  type(ESMF_State),    intent(inout)   :: IMPORT                  ! Import state
  type(ESMF_State),    intent(inout)   :: EXPORT                  ! Export state
  type(ESMF_Clock),    intent(inout)   :: CLOCK                   ! The clock
  integer, optional,   intent(  out)   :: RC                      ! Error code
  type (ESMF_Alarm)                    :: ALARM                   ! run alarm for ISSM component 
  
  ! ErrLog Variables
  character(len=ESMF_MAXSTR)           :: IAm
  integer                              :: STATUS
  character(len=ESMF_MAXSTR)           :: COMP_NAME

  type(MAPL_MetaComp), pointer         :: MAPL
  type(ESMF_VM)                        :: vm  

  ! regridding
  type(ESMF_RouteHandle)               :: routehandle_m2g         ! mesh to grid routehandle
  type(ESMF_RouteHandle)               :: routehandle_g2m         ! grid to mesh routehandle
  type(ESMF_Field)                     :: srcField                ! source field for regridding
  type(ESMF_Field)                     :: dstField                ! destination field for regridding
  type(T_ISSM_STATE), pointer          :: internal_state=>null()  ! stores the routehandles
  type(ISSM_WRAP)                      :: wrap                    ! wrapper for routehandle container
  type(ESMF_Mesh)                      :: mesh                    ! ESMF version of ISSM mesh
  integer(c_int)                       :: num_elements            ! number of elements on each PET
  integer                              :: IM, JM, local_dims(3)   ! local grid size

  ! tile information
  integer                              :: NT                      ! number of tiles

  ! ice elevation on mesh, grid, tile
  real(dp),    pointer, dimension(:)   :: ICEEL_MESH    => null() ! ice-sheet elevation on mesh elements
  real, pointer, dimension(:,:)        :: ICEEL_GRID    => null() ! ice-sheet elevation on grid 
  real, pointer, dimension(:)          :: ICEEL_TILE    => null() ! ice-sheet elevation on landice tiles
  real, pointer, dimension(:)          :: ICEEL         => null() ! pointer to ice-sheet elevation export state

  ! surface mass balance on mesh, grid, tile
  real(dp),    pointer, dimension(:)   :: ICESMB_MESH   => null() ! surface mass balce on mesh elements
  real,    pointer, dimension(:,:)     :: ICESMB_GRID   => null() ! surface mass balance on grid
  real, pointer, dimension(:)          :: ICESMB_TILE   => null() ! surface mass balance on tiles
  real, pointer, dimension(:)          :: ICESMB        => null() ! pointer to SMB import state


  ! physical parameters
  real(dp), parameter                  :: rho_ice   = 917.0       ! pure ice density [kg m^-3]
  real(dp)                             :: dt                      ! time step [s] (ISSM_DT set in AGCM.rc)

  integer                              :: ISSM_SAVEMESH           ! mesh save flag
  character(len=ESMF_MAXSTR)           :: ISSM_EXPDIR             ! directory containing ISSM input file

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

  ! save mesh output if desired
  call MAPL_GetResource(MAPL, ISSM_SAVEMESH, Label=trim(COMP_NAME)//"_SAVEMESH:",default=0,RC=STATUS)
  VERIFY_(STATUS)

  call MAPL_GetResource(MAPL, ISSM_EXPDIR, Label=trim(COMP_NAME)//"_EXPDIR:", RC=STATUS)
  VERIFY_(STATUS)

  ! run ISSM at specified time steps
  if ( ESMF_AlarmIsRinging (ALARM, RC=STATUS) ) then

    ! *************************************************************************** !
    ! BASIC SETUP
    ! *************************************************************************** !

    ! get timestep for ISSM
    call MAPL_GetResource (MAPL, dt, Label=trim(COMP_NAME)//"_DT:",DEFAULT=604800.0, RC=STATUS)
    VERIFY_(STATUS)

    ! get number of mesh elements
    call ESMF_MeshGet(mesh,elementCount=num_elements)

    ! allocate SMB forcing (input to ISSM) and ice-elevation output (export from ISSM)
    allocate(ICEEL_MESH(num_elements))    

    call ESMF_VMBarrier(vm, rc=status)
    VERIFY_(STATUS)

    ! set smb and surface to zero (not sure this is needed...)
    ICEEL_MESH(:) = 0.0_dp

    ! get routehandles for regridding
    call ESMF_UserCompGetInternalState(GC, 'ISSM_WRAP', wrap, status); VERIFY_(STATUS)
    internal_state => wrap%ptr
    routehandle_m2g = internal_state%routehandle_m2g  ! mesh to grid 
    routehandle_g2m = internal_state%routehandle_g2m  ! grid to mesh
    call MAPL_LocStreamGet(internal_state%locstream, NT_LOCAL=NT, _RC)
    call MAPL_GridGet(internal_state%grid, localCellCountPerDim=local_dims, _RC)
    IM = local_dims(1)
    JM = local_dims(2)

    ! *************************************************************************** !
    ! IMPORT SMB (surface mass balance) & REGRID FROM TILES [to grid] to MESH 
    ! *************************************************************************** !
    
    call MAPL_GetPointer(IMPORT,ICESMB, 'ICESMB' , RC=STATUS); VERIFY_(STATUS)

    ! allocate tiles for SMB (MKTILE)
    if(associated(ICESMB) .and. .not.associated(ICESMB_TILE)) then
      allocate(ICESMB_TILE(NT), STAT=STATUS)
      VERIFY_(STATUS)
      ICESMB_TILE = MAPL_Undef
    end if
    
    ! copy import values into tile array (necessary?)
    ICESMB_TILE(:) = ICESMB(:)

    ! allocate SMB on grid
    allocate( ICESMB_GRID(IM,JM), STAT=STATUS ); VERIFY_(STATUS)

    ! transform ICESMB from tile to grid
    ! NOTE: we use the "transpose" option with MAPL_LocStreamTransformG2T 
    ! (rather than MAPL_LocStreamTransformT2G) because the "default" value is zero
    ! (rather than MAPL_UNDEF, which leads to errors when regridding onto mesh)
    call MAPL_LocStreamTransform( internal_state%LOCSTREAM, ICESMB_TILE, ICESMB_GRID, TRANSPOSE=.true., RC=STATUS)
    VERIFY_(STATUS)

    ! create source field: SMB on grid
    srcField = ESMF_FieldCreate(grid=internal_state%grid,farrayPtr=ICESMB_GRID, datacopyflag=ESMF_DATACOPY_VALUE,rc=STATUS)
    VERIFY_(STATUS)

    ! create destination field: SMB on mesh elements
    dstField = ESMF_FieldCreate(mesh=mesh,typekind=ESMF_TYPEKIND_R8,meshloc=ESMF_MESHLOC_ELEMENT, rc=STATUS)
    VERIFY_(STATUS)

    ! regrid SMB from grid to mesh
    call ESMF_FieldRegrid(srcField, dstField, routehandle_g2m, RC=STATUS); VERIFY_(STATUS)

    ! get pointer to SMB on mesh
    call ESMF_FieldGet(dstField,farrayPtr=ICESMB_MESH,RC=STATUS); VERIFY_(STATUS)

    ! save ICESMB on mesh elements 
    if (ISSM_SAVEMESH/=0) then
        call ESMF_FieldWrite(dstField, trim(ISSM_EXPDIR)//"/icesmb.nc", variableName='ICESMB',rc=STATUS)
    end if 
    
    ! convert SMB to units of [m/s] (ice-equivalent) for ISSM
    ICESMB_MESH = ICESMB_MESH/rho_ice

    ! *************************************************************************** !
    !  RUN ISSM WITH SMB INPUT AND ICE-ELEVATION OUTPUT
    ! *************************************************************************** !

    ! NOTE: do we need the barriers before/after ISSM run?
    call ESMF_VMBarrier(vm, rc=status); VERIFY_(STATUS)

    ! call run method from ISSM library 
    call RunISSM(dt, c_loc(ICESMB_MESH), c_loc(ICEEL_MESH))

    call ESMF_VMBarrier(vm, rc=status); VERIFY_(STATUS)

    ! *************************************************************************** !
    ! REGRID ICEEL (ice elevation) FROM MESH [to grid] to TILES 
    ! *************************************************************************** !
     
    ! destroy regridding fields so they can be reused
    call ESMF_FieldDestroy(srcField,rc=STATUS); VERIFY_(STATUS)
    call ESMF_FieldDestroy(dstField,rc=STATUS); VERIFY_(STATUS)

    ! WY note: after getting ICESMB_Mesh, assign it ICEEL_tile directly for the output


    ! create source field: ice elevation on mesh elements
    srcField = ESMF_FieldCreate(mesh=mesh,farrayPtr=ICEEL_MESH,meshloc=ESMF_MESHLOC_ELEMENT, & 
    datacopyflag=ESMF_DATACOPY_VALUE,rc=STATUS)
    VERIFY_(STATUS)

    ! save ICEEL on mesh elements 
    if (ISSM_SAVEMESH/=0) then
        call ESMF_FieldWrite(srcField, trim(ISSM_EXPDIR)//"/iceel.nc", variableName='ICEEL',rc=STATUS)
    end if 
    
    ! create destination field: ice elevation on grid
    dstField = ESMF_FieldCreate(grid=internal_state%grid,typekind=ESMF_TYPEKIND_R4,rc=STATUS)
    VERIFY_(STATUS)

    ! regrid ice elevation from mesh to grid
    call ESMF_FieldRegrid(srcField, dstField, routehandle_m2g, RC=STATUS); VERIFY_(STATUS)

    ! get pointer to ice elevation on grid
    call ESMF_FieldGet(dstField,farrayPtr=ICEEL_GRID,RC=STATUS); VERIFY_(STATUS)
    
    ! get pointer to export
    call MAPL_GetPointer(EXPORT  , ICEEL , 'ICEEL' , alloc=.true. , RC=STATUS); VERIFY_(STATUS)

    ! allocate tiles for ice elevation (MKTILE)
    if(associated(ICEEL) .and. .not.associated(ICEEL_TILE)) then
      allocate(ICEEL_TILE(NT), STAT=STATUS)
      VERIFY_(STATUS)
      ICEEL_TILE = MAPL_Undef
    end if
   
    ! transform from grid to tiles  
    call MAPL_LocStreamTransform( internal_state%LOCSTREAM, ICEEL_TILE,ICEEL_GRID, RC=STATUS)
    VERIFY_(STATUS)

    ! assign to export pointer
    if(associated(ICEEL)) ICEEL(:) = ICEEL_TILE(:)

    call ESMF_FieldDestroy(srcField,rc=STATUS); VERIFY_(STATUS)
    call ESMF_FieldDestroy(dstField,rc=STATUS); VERIFY_(STATUS)

  end if 

  call ESMF_VMBarrier(vm, rc=status)
  VERIFY_(STATUS)

  if(associated(ICEEL_MESH))  deallocate(ICEEL_MESH)
  if(associated(ICEEL_TILE))  deallocate(ICEEL_TILE)
  if(associated(ICESMB_TILE)) deallocate(ICESMB_TILE)
  if(associated(ICESMB_GRID)) deallocate(ICESMB_GRID)

  call MAPL_TimerOff(MAPL,"RUN"  )
  call MAPL_TimerOff(MAPL,"TOTAL")
 
  RETURN_(ESMF_SUCCESS)
 
 end subroutine RUN

 !BOP
    
!IROUTINE: Finalize   -- Finalize method for ISSM 

!INTERFACE:

 subroutine Finalize ( GC, IMPORT, EXPORT, CLOCK, RC )

    !ARGUMENTS:
    type(ESMF_GridComp), intent(INOUT) :: GC     ! Gridded component 
    type(ESMF_State),    intent(INOUT) :: IMPORT ! Import state
    type(ESMF_State),    intent(INOUT) :: EXPORT ! Export state
    type(ESMF_Clock),    intent(INOUT) :: CLOCK  ! The supervisor clock
    integer, optional,   intent(  OUT) :: RC     ! Error code:
    
    !EOP
    type(MAPL_MetaComp), pointer       :: MAPL 
    
    ! ErrLog Variables
    character(len=ESMF_MAXSTR)	       :: IAm
    integer			                   :: STATUS
    character(len=ESMF_MAXSTR)         :: COMP_NAME

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
