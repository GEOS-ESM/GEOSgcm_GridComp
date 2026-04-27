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
!   Imports: ICESMB (defined on landice tiles)
!   Exports: ICESURF, ICETHICK, ICEVEL, EQSMB (defined on mesh)
! *** NOTES: 
!            (*) currently we run over all input files (*.bin) that are found in ISSM_EXPDIR 
!                (e.g., Greenland + Antarctica + any other glaciers that have been configured)
!            (*) ISSM meshes are internal to ISSM (C++ source)--we create an ESMF_MESH version for regridding     
!                imports/exports that is the global combination of all ISSM meshes
!            (*) we transform imports from landice tiles to attached grid, then regrid to the mesh 
!            (*) ISSM outputs are saved with HISTORY via a 'mesh tile space' developed by Weiyuan Jiang (GMAO SI Team)  
!

! !USES:
use iso_fortran_env, only: dp=>real64
use iso_c_binding, only: c_ptr, c_double, c_f_pointer, c_null_char, c_char, c_loc, c_int
use ESMF
use MAPL
use GEOS_UtilsMod

implicit none

! declare interface to the ISSM C++ library (arguments described in Initialize & Run below)
interface
subroutine InitializeISSM(expdir, num_elements, num_nodes, comm) bind(c, name="InitializeISSM")
  import :: c_char, c_int
  character(c_char), dimension(*) :: expdir
  integer(c_int)                  :: num_elements
  integer(c_int)                  :: num_nodes
  integer(c_int)                  :: comm
end subroutine InitializeISSM
    
subroutine RunISSM(dt, gcm_forcings, issm_outputs, elementConn) bind(C,NAME="RunISSM")
   import :: c_ptr, c_double
   real(c_double),   value        :: dt
   type(c_ptr),      value        :: gcm_forcings
   type(c_ptr),      value        :: issm_outputs
   type(c_ptr),      value        :: elementConn
end subroutine RunISSM

subroutine InputFromRestarts(gcm_restarts, elementConn) bind(C,NAME="InputFromRestarts")
  import :: c_ptr
  type(c_ptr),       value        :: gcm_restarts
  type(c_ptr),       value        :: elementConn
end subroutine InputFromRestarts

subroutine GetNodesISSM(nodeIds, nodeCoords) bind(C,NAME="GetNodesISSM")
   import :: c_ptr
   type(c_ptr),      value        :: nodeIds
   type(c_ptr),      value        :: nodeCoords 
end subroutine GetNodesISSM

subroutine GetElementsISSM(elementIds, elementConn, elementCoords, glacIds) bind(C,NAME="GetElementsISSM")
  import :: c_ptr
  type(c_ptr),       value        :: elementIds
  type(c_ptr),       value        :: elementConn
  type(c_ptr),       value        :: elementCoords
  type(c_ptr),       value        :: glacIds
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

public :: T_ISSM_TILE_STATE
public :: ISSM_TILE_WRAP
! define ISSM export as internal variables, will be used by the landice grid comp

type T_ISSM_TILE_STATE
    real, pointer :: ICESURF_TILE(:)
    real, pointer :: ICETHICK_TILE(:)
    real, pointer :: ICEVEL_TILE(:)
    real, pointer :: ICESMB_TILE(:)
end type T_ISSM_TILE_STATE

type ISSM_TILE_WRAP
   type(T_ISSM_TILE_STATE), pointer :: ptr=>null()
end type ISSM_TILE_WRAP

! private internal state for regridding 
type T_ISSM_STATE
  private
  type(ESMF_RouteHandle)        :: routehandle_m2g ! routehandle for regridding mesh to grid
  type(ESMF_RouteHandle)        :: routehandle_g2m ! routehandle for regridding grid to mesh
  type(ESMF_RouteHandle)        :: halohandle      ! routehandle for field halos
  integer, pointer,dimension(:) :: halo_idx        ! indices of halo nodes in arrays
  integer, pointer,dimension(:) :: owned_idx       ! indices of owned nodes in arrays
  integer, pointer,dimension(:) :: halolist        ! list of halo nodeIds
  type(ESMF_DistGrid)           :: nodalDistgrid   ! distgrid (owned nodes)
  type(ESMF_GRID)               :: grid            ! original grid (atmosphere)
  type(MAPL_LocStream)          :: locstream       ! original locstream (landice tiles)
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
    call ESMF_ConfigGetAttribute(CF, dt, Label=trim(COMP_NAME)//"_DT:", DEFAULT=302400.0, RC=STATUS)
    VERIFY_(STATUS)


    ! Set the state variable specs.
!-----------------------------------

!   Import states: ICESMB is imported via the ISSM_TILE internal state

!   Export states:
    call MAPL_AddExportSpec(GC,                    &
         SHORT_NAME = 'ICESURF',                   &
         LONG_NAME  = 'ice_sheet_elevation',       &
         UNITS      = 'm',                         &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone,          &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                    &
         SHORT_NAME = 'ICEVX',                     &
         LONG_NAME  = 'ice_velocity_x_direction',  &
         UNITS      = 'm s-1',                     &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone,          &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                    &
        SHORT_NAME = 'ICEVY',                      &
        LONG_NAME  = 'ice_velocity_y_direction',   &
        UNITS      = 'm s-1',                      &
        DIMS       = MAPL_DimsTileOnly,            &
        VLOCATION  = MAPL_VLocationNone,           &
        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                    &
         SHORT_NAME = 'ICETHICK',                  &
         LONG_NAME  = 'ice_thickness',             &
         UNITS      = 'm',                         &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone,          &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                    &
         SHORT_NAME = 'ICESMB',                    &
         LONG_NAME  = 'ice_surface_mass_balance',  &
         UNITS      = 'kg m-2 s-1',                &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone,          &
         RC=STATUS  )
    VERIFY_(STATUS)

    !   Internal states:
    call MAPL_AddInternalSpec(GC,                  &
         SHORT_NAME = 'ICESURF',                   &
         LONG_NAME  = 'ice_sheet_elevation',       &
         UNITS      = 'm',                         &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone,          &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                  &
         SHORT_NAME = 'ICETHICK',                  &
         LONG_NAME  = 'ice_sheet_thickness',       &
         UNITS      = 'm',                         &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone,          &
         RC=STATUS  )
    VERIFY_(STATUS)
    
    call MAPL_AddInternalSpec(GC,                  &
         SHORT_NAME = 'IMLS',                 &
         LONG_NAME  = 'ice_mask_levelset',         &
         UNITS      = 'none',                      &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone,          &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                  &
         SHORT_NAME = 'OMLS',                 &
         LONG_NAME  = 'ocean_mask_levelset',       &
         UNITS      = 'none',                      &
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
    integer                                :: localPET                ! ~mpi rank
    
    ! mesh information
    type(ESMF_Mesh)                        :: mesh                    ! ESMF_Mesh representation of ISSM mesh
    integer, pointer, dimension(:)         :: elementTypes  => null() ! element geometry type (triangles)
    integer(c_int)                         :: num_elements            ! number of elements on PET
    integer(c_int)                         :: num_nodes               ! number of nodes on PET
    integer(c_int)                         :: num_owned_nodes         ! number of nodes owned by this PET (<=num_nodes)
    integer, pointer, dimension(:)         :: elementIds    => null() ! list of elements local to PET
    integer, pointer, dimension(:)         :: elementConn   => null() ! element connectivity (nodes indices)
    real(dp),pointer, dimension(:)         :: elementCoords => null() ! element centroids
    real(dp),pointer,dimension(:)          :: nodeCoords    => null() ! node coordinates (longitude,latitude)
    integer, pointer, dimension(:)         :: nodeIds       => null() ! Global IDs of nodes local to PET
    integer, pointer, dimension(:)         :: nodeOwners    => null() ! Specify which PET owns each node
    integer, pointer, dimension(:)         :: glacIds       => null() ! glacier ID for each element
    
    ! regridding 
    type(ESMF_Grid)                        :: grid                    ! atmospheric grid
    type(ESMF_RouteHandle)                 :: routehandle_m2g         ! routehandle for regridding mesh to grid
    type(ESMF_RouteHandle)                 :: routehandle_g2m         ! routehandle for regridding grid to mesh
    type(ESMF_Field)                       :: meshField               ! field on mesh
    type(ESMF_Field)                       :: gridField               ! field on grid
    type(T_ISSM_STATE), pointer            :: internal_state          ! store the routehandles for access during run
    type(ISSM_WRAP)                        :: wrap                    ! wrapper for routehandle container

    ! field halo
    integer                                :: num_halo_nodes          ! num_nodes minus num_owned_nodes
    type(ESMF_RouteHandle)                 :: halohandle              ! routehandle for field halos
    integer, pointer, dimension(:)         :: halolist      => null() ! list of halo nodeIds
    integer, pointer, dimension(:)         :: ownedNodeIds  => null() ! nodeIds excluding halolist
    type(ESMF_DistGrid)                    :: nodalDistgrid           ! distgrid (owned nodes)
    type(ESMF_Array)                       :: meshArray               ! array for creating mesh fields
    integer, pointer,dimension(:)          :: halo_idx      => null() ! indices of halo nodes in arrays
    integer, pointer,dimension(:)          :: owned_idx     => null() ! indices of owned nodes in arrays

    ! owned node coordinates (longitude,latitude)
    real(dp),pointer,dimension(:)          :: ownedNodeCoords => null() 
    real, allocatable, dimension(:)        :: ownedNodeLons, ownedNodeLats

    ! command-line arguments to initialize ISSM 
    integer                                :: i,j,k                   ! loop indices
    character(len=ESMF_MAXSTR)             :: ISSM_EXPDIR             ! directory containing ISSM input file
    character(len=ESMF_MAXSTR)             :: EXPDIR                  ! C++ compatible ISSM_EXPDIR string

    ! variables for saving mesh output
    type(ESMF_Field)                       :: field                   ! for fieldwrites
    real(dp),    pointer, dimension(:)     :: field_saver => null()
    real(dp),    pointer, dimension(:)     :: nodelons    => null()
    real(dp),    pointer, dimension(:)     :: nodelats    => null()

    integer, pointer, dimension(:)         :: conn_slice  => null()
    character(len=16) :: varname
    character(len=1)  :: istr

    ! variables for mesh tile space
    type(ESMF_Grid)                        :: mesh_grid               
    type(MAPL_LocStream)                   :: mesh_locstream   

    ! variables for masking the mesh seam (triangles that cross +/-180 longitude)
    ! (needed for elements, but unsure if this is needed for regridding fields defined on nodes)
    real(dp)                               :: dlon,lon1,lon2,lon3
    integer, pointer, dimension(:)         :: elementMask => null()
    integer, pointer, dimension(:)         :: nodeMask => null()
    integer                                :: n1,n2,n3

    ! Get the target components name and set-up traceback handle.
    ! -----------------------------------------------------------

    Iam = "Initialize"
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=status )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // trim(Iam)

    ! Get my internal MAPL_Generic state
    !-----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)

    call ESMF_VMGetCurrent(vm, rc=STATUS)
    VERIFY_(STATUS)

    call ESMF_VMGet(vm,mpiCommunicator=comm,localPet=localPET,rc=STATUS)
    VERIFY_(STATUS)

    ! ****************************************************
    ! call ISSM initialize C++ code so we can set up mesh

    call MAPL_GetResource(MAPL, ISSM_EXPDIR, Label="ISSM_EXPDIR:", RC=STATUS)
    VERIFY_(STATUS)

    EXPDIR = trim(ISSM_EXPDIR)//c_null_char ! create string for C++
    
    ! Call the C++ function for initializing ISSM
    ! gets the number of elements and nodes of the mesh
    call InitializeISSM(EXPDIR, num_elements, num_nodes, comm)

    !allocate mesh-related pointers
    allocate(nodeCoords(2*num_nodes))
    allocate(nodeIds(num_nodes))
    allocate(elementTypes(num_elements))
    allocate(elementIds(num_elements))
    allocate(glacIds(num_elements))
    allocate(elementConn(3*num_elements))
    allocate(elementCoords(2*num_elements))
    allocate(elementMask(num_elements))
    allocate(nodeMask(num_nodes))
    allocate(nodeOwners(num_nodes))

    ! get information about nodes and elements
    ! node coords and element coords (centroids) are in (lon,lat)
    call GetNodesISSM(c_loc(nodeIds), c_loc(nodeCoords)) 
    call GetElementsISSM(c_loc(elementIds), c_loc(elementConn), c_loc(elementCoords),c_loc(glacIds))

    elementTypes(:) = ESMF_MESHELEMTYPE_TRI ! triangular elements

    ! mask triangles that cross the seam (longitude +/- 180)
    ! possibly mask nodes... (note: you don't have to 'activate' these
    ! masks, they can just be associated with the mesh)
    ! note: I think this is only relevant in regridding when fields are
    !       defined on esmf_meshloc_element (rather than esmf_meshloc_node)
    elementMask(:) = 0
    nodeMask(:) = 0
    do j=1,num_elements
      n1 = elementConn(3*(j-1)+1)
      n2 = elementConn(3*(j-1)+2)
      n3 = elementConn(3*(j-1)+3)
      lon1 = nodeCoords(2*n1-1)
      lon2 = nodeCoords(2*n2-1)
      lon3 = nodeCoords(2*n3-1)
      dlon = maxval((/lon1,lon2,lon3/)) - minval((/lon1,lon2,lon3/))
      if ( dlon>180.0 ) then
        elementMask(j) = 1
        nodeMask(n1) = 1
        nodeMask(n2) = 1
        nodeMask(n3) = 1
      end if
    end do

    ! create the ESMF mesh from ISSM mesh properties
    mesh = ESMF_MeshCreate(parametricDim=2, spatialDim=2, nodeIds=nodeIds, nodeCoords=nodeCoords, &
            elementIds=elementIds, elementTypes=elementTypes, elementConn=elementConn,elementMask=elementMask,& 
            nodeMask=nodeMask,elementCoords=elementCoords,coordSys=ESMF_COORDSYS_SPH_DEG, rc=STATUS)
    VERIFY_(STATUS)

    ! associate ESMF_Mesh representation of ISSM mesh with GC for regridding imports/exports in Run method       
    call ESMF_GridCompSet(GC,mesh=mesh,rc=STATUS)
    VERIFY_(STATUS)

    ! set up field halos
    !-----------------------------------
    call ESMF_MeshGet(mesh=mesh,nodeOwners=nodeOwners,numOwnedNodes=num_owned_nodes,nodalDistgrid=nodalDistgrid)

    num_halo_nodes = num_nodes - num_owned_nodes
    allocate(halolist(num_halo_nodes))
    allocate(ownedNodeCoords(2*num_owned_nodes))
    allocate(ownedNodeIds(num_owned_nodes))
    allocate(halo_idx(num_halo_nodes))
    allocate(owned_idx(num_owned_nodes))

    call ESMF_MeshGet(mesh=mesh,ownedNodeCoords=ownedNodeCoords)
    
    ! get list of (global) nodeIds that are halos on this PET
    ! and create a mask to remove these values from arrays
    i=1; k=1
    do j=1,num_nodes
    if (nodeOwners(j)/= localPET) then
      halolist(i) = nodeIds(j)
      halo_idx(i) = j
      i = i+1
    else
      ownedNodeIds(k) = nodeIds(j)
      owned_idx(k) = j
      k = k+1
    end if 
    end do

    ! create array with halo information
    meshArray=ESMF_ArrayCreate(nodalDistgrid,typekind=ESMF_TYPEKIND_R8,haloSeqIndexList=halolist,rc=STATUS)
    VERIFY_(STATUS)

    ! create field on ISSM mesh 
    meshField=ESMF_FieldCreate(mesh, array=meshArray, meshLoc=ESMF_MESHLOC_NODE, rc=STATUS)
    VERIFY_(STATUS)

    ! store the halo operation in a routehandle
    call ESMF_FieldHaloStore(meshField, routehandle=halohandle, rc=STATUS)
    VERIFY_(STATUS)

    ! Set up regridding next
    !-----------------------------------
    ! get atmospheric (attached) grid 
    call ESMF_GridCompGet( GC, GRID=grid, RC=status ); VERIFY_(STATUS)
    
    ! create field on atmospheric grid
    gridField = ESMF_FieldCreate(grid=grid,typekind=ESMF_TYPEKIND_R4,rc=STATUS); VERIFY_(STATUS)
    
    ! create routehandle for mesh-to-grid regridding (set srcMaskValues to 1 if needed... )
    call ESMF_FieldRegridStore(srcField=meshField, dstField=gridField,routehandle=routehandle_m2g,& 
    unmappedaction=ESMF_UNMAPPEDACTION_IGNORE,extrapmethod=ESMF_EXTRAPMETHOD_CREEP,&
    extrapNumLevels=1,rc=STATUS)
    VERIFY_(STATUS)

    ! create routehandle for grid-to-mesh regridding (set dstMaskValues to 1 if needed... )
    call ESMF_FieldRegridStore(srcField=gridField, dstField=meshField,routehandle=routehandle_g2m,& 
    unmappedaction=ESMF_UNMAPPEDACTION_IGNORE,extrapmethod=ESMF_EXTRAPMETHOD_NEAREST_D,rc=STATUS)
    VERIFY_(STATUS)
    
    ! create component's private internal state
    ! stores everything needed for regrid and halo operations during run method
    allocate(internal_state, stat=STATUS)
    VERIFY_(STATUS)
    allocate(internal_state%halo_idx(num_halo_nodes))
    allocate(internal_state%owned_idx(num_owned_nodes))
    allocate(internal_state%halolist(num_halo_nodes))
    internal_state%routehandle_m2g = routehandle_m2g
    internal_state%routehandle_g2m = routehandle_g2m 
    internal_state%halohandle = halohandle 
    internal_state%halo_idx = halo_idx
    internal_state%owned_idx = owned_idx  
    internal_state%grid = grid
    internal_state%halolist = halolist
    internal_state%nodalDistgrid = nodalDistgrid
    call MAPL_Get(MAPL, LocStream = internal_state%locstream, _RC)

    ! wrap the private internal state
    wrap%ptr => internal_state
    call ESMF_UserCompSetInternalState ( GC, 'ISSM_WRAP', wrap, status )
    VERIFY_(STATUS)

  ! ------------------------------ BEGIN SAVE MESH ---------------------------
    allocate(field_saver(num_elements))
    allocate(nodelons(num_nodes))
    allocate(nodelats(num_nodes))
    allocate(conn_slice(num_elements))

    nodelons = nodeCoords(1::2)
    nodelats = nodeCoords(2::2)
    
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
        call ESMF_FieldWrite(field, trim(ISSM_EXPDIR)//"/"//"issm_mesh.nc", variableName=varname, rc=STATUS)
        VERIFY_(STATUS)
        call ESMF_FieldDestroy(field, rc=STATUS)
        VERIFY_(STATUS)
    
        ! ---- latitude ----
        varname = "node" // istr // "lat"
        field_saver = nodelats(conn_slice)
        field = ESMF_FieldCreate(mesh=mesh, farrayPtr=field_saver,meshloc=ESMF_MESHLOC_ELEMENT, rc=STATUS)
        VERIFY_(STATUS)
        call ESMF_FieldWrite(field, trim(ISSM_EXPDIR)//"/"//"issm_mesh.nc", variableName=varname, rc=STATUS)
        VERIFY_(STATUS)
        call ESMF_FieldDestroy(field, rc=STATUS)
        VERIFY_(STATUS)
    end do

    ! ---- save element IDs ----
    varname = "elementIds"
    field_saver = elementIds(:)
    field = ESMF_FieldCreate(mesh=mesh, farrayPtr=field_saver,meshloc=ESMF_MESHLOC_ELEMENT, rc=STATUS)
    VERIFY_(STATUS)
    call ESMF_FieldWrite(field, trim(ISSM_EXPDIR)//"/"//"issm_mesh.nc", variableName=varname, rc=STATUS)
    VERIFY_(STATUS)
    call ESMF_FieldDestroy(field, rc=STATUS)
    VERIFY_(STATUS)
  
    ! ---- save glacier IDs ----
    varname = "glacIds"
    field_saver = glacIds(:)
    field = ESMF_FieldCreate(mesh=mesh, farrayPtr=field_saver,meshloc=ESMF_MESHLOC_ELEMENT, rc=STATUS)
    VERIFY_(STATUS)
    call ESMF_FieldWrite(field, trim(ISSM_EXPDIR)//"/"//"issm_mesh.nc", variableName=varname, rc=STATUS)
    VERIFY_(STATUS)
    call ESMF_FieldDestroy(field, rc=STATUS)
    VERIFY_(STATUS)
  ! ------------------------------ END SAVE MESH ---------------------------

    ! Create losctream that match mesh element id, then set it to this GC and MAPL
    ! note: original attached/atmospheric grid and landice tile locsstream have
    ! been stored in the internal state
    allocate(ownedNodeLons(num_owned_nodes))
    allocate(ownedNodeLats(num_owned_nodes))
    ownedNodeLons = ownedNodeCoords(1::2)*MAPL_DEGREES_TO_RADIANS
    ownedNodeLats = ownedNodeCoords(2::2)*MAPL_DEGREES_TO_RADIANS

    mesh_grid = create_mesh_grid(_RC)
    call MAPL_LocstreamCreate(mesh_locstream, mesh_grid, local_id=ownedNodeIds, &
              tilelons=ownedNodeLons, tilelats=ownedNodeLats,  _RC)
    call MAPL%grid%set(mesh_grid, _RC)
    call ESMF_GridCompSet(gc, grid=mesh_grid, _RC)
    call MAPL_set(MAPL, locstream = mesh_locstream, _RC)

    ! deallocate pointers
    if(associated(field_saver))     deallocate(field_saver)
    if(associated(nodelons))        deallocate(nodelons)
    if(associated(nodelats))        deallocate(nodelats)
    if(associated(conn_slice))      deallocate(conn_slice)
    if(associated(nodeCoords))      deallocate(nodeCoords)
    if(associated(nodeIds))         deallocate(nodeIds)
    if(associated(elementTypes))    deallocate(elementTypes)
    if(associated(elementIds))      deallocate(elementIds)
    if(associated(elementConn))     deallocate(elementConn)
    if(associated(elementCoords))   deallocate(elementCoords)
    if(associated(glacIds))         deallocate(glacIds)
    if(associated(elementMask))     deallocate(elementMask)
    if(associated(nodeMask))        deallocate(nodeMask)
    if(associated(halo_idx))        deallocate(halo_idx)
    if(associated(owned_idx))       deallocate(owned_idx)
    if(associated(halolist))        deallocate(halolist)
    if(associated(ownedNodeCoords)) deallocate(ownedNodeCoords)
    if(associated(ownedNodeIds))    deallocate(ownedNodeIds)
    if(associated(nodeOwners))      deallocate(nodeOwners)

    ! destroy fields and arrays
    call ESMF_FieldDestroy(gridField, rc=STATUS); VERIFY_(STATUS)
    call ESMF_FieldDestroy(meshField, rc=STATUS); VERIFY_(STATUS)
    call ESMF_ArrayDestroy(meshArray, rc=STATUS); VERIFY_(STATUS)

    ! generic initialize 
    call MAPL_GenericInitialize( GC, IMPORT, EXPORT, CLOCK, RC=STATUS )
    VERIFY_(STATUS)
    RETURN_(ESMF_SUCCESS)

    contains

       function create_mesh_grid(rc) result(mesh_grid)
          type (ESMF_Grid) :: mesh_grid
          integer, optional, intent(out) :: rc
          integer :: status, nDEs, num(1)
          real(kind=8), pointer :: centers(:,:)
          integer, allocatable  :: IMs(:)
          

          !comm, VM, num_owned_nodes are from containing subroutine
          call ESMF_VMGet(vm, petcount=nDEs,  _RC) 
          allocate(IMS(nDEs))
          num(1) = num_owned_nodes
          call MAPL_CommsAllGather(vm, num, 1, IMs, 1, _RC) 

          ! create a mesh-grid in 1D
          mesh_grid = ESMF_GridCreate(        &
               name='MESH_GRID',              &
               countsPerDEDim1=IMs,           &
               countsPerDEDim2=[1],           &
               indexFlag=ESMF_INDEX_DELOCAL,  &
               coordDep1 = (/1,2/),           &
               coordDep2 = (/1,2/),           &
               gridEdgeLWidth = (/0,0/),      &
               gridEdgeUWidth = (/0,0/),      &
               _RC)
          ! coord and centers are required for a valid grid,
          ! even if their values don't make sense;
          ! later on, the coord will be set to element's lat lon.
          call ESMF_GridAddCoord(mesh_grid, _RC)
          _VERIFY(status)

          call ESMF_GridGetCoord(mesh_grid, coordDim=1, localDE=0, &
               staggerloc=ESMF_STAGGERLOC_CENTER, &
               farrayPtr=centers, _RC)
          centers = 0 
          call ESMF_GridGetCoord(mesh_grid, coordDim=2, localDE=0, &
               staggerloc=ESMF_STAGGERLOC_CENTER, &
               farrayPtr=centers, _RC)
          centers = 0 

          _RETURN(_SUCCESS)
       end function create_mesh_grid      

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
  type(ESMF_Alarm)                     :: ALARM                   ! run alarm for ISSM component 
  
  ! ErrLog Variables
  character(len=ESMF_MAXSTR)           :: IAm
  integer                              :: STATUS
  character(len=ESMF_MAXSTR)           :: COMP_NAME

  type(MAPL_MetaComp), pointer         :: MAPL
  type(ESMF_State)                     :: INTERNAL
  type(ESMF_VM)                        :: vm  

  ! regridding
  type(ESMF_RouteHandle)               :: routehandle_m2g         ! mesh to grid routehandle
  type(ESMF_RouteHandle)               :: routehandle_g2m         ! grid to mesh routehandle
  type(ESMF_Field)                     :: srcField                ! source field for regridding
  type(ESMF_Field)                     :: dstField                ! destination field for regridding
  type(T_ISSM_STATE), pointer          :: internal_state=>null()  ! stores everything needed for regridding
  type(ISSM_WRAP)                      :: wrap                    ! wrapper for routehandle container
  type(ESMF_Mesh)                      :: mesh                    ! ESMF version of ISSM mesh
  integer                              :: num_elements            ! number of elements on PET
  integer                              :: num_nodes               ! number of nodes on PET
  integer                              :: num_owned_nodes         ! number of nodes owned by PET
  integer                              :: IM, JM, local_dims(3)   ! local grid size

  ! halo information
  integer                              :: num_halo_nodes          ! num_nodes minus num_owned_nodes
  type(ESMF_RouteHandle)               :: halohandle              ! field halo routehandle
  integer, pointer, dimension(:)       :: halolist      => null() ! list of halo nodeIds
  integer, pointer,dimension(:)        :: halo_idx      => null() ! indices of halo nodes in arrays
  integer, pointer,dimension(:)        :: owned_idx     => null() ! indices of owned nodes in arrays
  type(ESMF_DistGrid)                  :: nodalDistgrid           ! distgrid (owned nodes)
  type(ESMF_Array)                     :: meshArray               ! array for creating mesh fields

  ! mesh information
  integer, pointer, dimension(:)       :: elementConn   => null() ! element connectivity (nodes indices)

  ! tile information
  integer                              :: NT                      ! number of landice tiles

  type(T_ISSM_TILE_STATE), pointer     :: issm_tile_state
  type(ISSM_TILE_WRAP)                 :: issm_tile_wrap

  ! surface mass balance on mesh and landice tiles
  real(dp), pointer, dimension(:)      :: ICESMB_MESH   => null() ! surface mass balce on mesh elements
  real, pointer, dimension(:)          :: ICESMB_TILE   => null() ! surface mass balance on landice tiles
  real, pointer, dimension(:)          :: ICESMB_EX     => null() ! pointer to SMB export (mesh)

  ! ISSM Outputs
  integer :: num_outputs = 6                                      ! number of outputs 
  real(dp),    pointer, dimension(:)   :: ISSM_OUTPUTS  => null() ! pointer containing all outputs

  ! ice-surface elevation on mesh and landice tiles
  real(dp),    pointer, dimension(:)   :: ICESURF_MESH  => null() ! ice elevation on mesh
  real, pointer, dimension(:)          :: ICESURF_TILE  => null() ! ice elevation on landice tiles
  real, pointer, dimension(:)          :: ICESURF_EX    => null() ! pointer to export state (mesh tiles)
  real, pointer, dimension(:)          :: ICESURF_IN    => null() ! pointer to internal state (mesh tiles)

  ! ice thickness on mesh and landice tiles
  real(dp),    pointer, dimension(:)   :: ICETHICK_MESH => null() ! ice thickness on mesh
  real, pointer, dimension(:)          :: ICETHICK_TILE => null() ! ice thickness on landice tiles
  real, pointer, dimension(:)          :: ICETHICK_EX   => null() ! pointer to ice thickness export state (mesh tiles)
  real, pointer, dimension(:)          :: ICETHICK_IN   => null() ! pointer to ice thicknesss internal state (mesh tiles)

  ! ice-flow velocity in x direction (in projection coordinates)
  real(dp),    pointer, dimension(:)   :: ICEVX_MESH   => null() ! ice x-velocity on mesh
  real, pointer, dimension(:)          :: ICEVX_EX     => null() ! pointer to export state (mesh tiles)

  ! ice-flow velocity in y direction (in projection coordinates)
  real(dp),    pointer, dimension(:)   :: ICEVY_MESH   => null() ! ice y-velocity on mesh
  real, pointer, dimension(:)          :: ICEVY_EX     => null() ! pointer to export state (mesh tiles)

  ! ice mask level set (tracks glacier terminus)
  real(dp),    pointer, dimension(:)   :: IMLS_MESH   => null() ! ice mask level set
  real, pointer, dimension(:)          :: IMLS_IN     => null() ! pointer to internal state (mesh tiles)

  ! ocean mask level set (tracks grounding line)
  real(dp),    pointer, dimension(:)   :: OMLS_MESH   => null() ! ocean mask level set
  real, pointer, dimension(:)          :: OMLS_IN     => null() ! pointer to internal state (mesh tiles)

  ! ice-flow speed on mesh and landice tiles  
  real(dp),    pointer, dimension(:)   :: ICEVEL_MESH => null() ! ice flow speed on mesh tiles
  real, pointer, dimension(:)          :: ICEVEL_TILE => null() ! ice flow speed on landice tiles
 
  ! physical parameters
  real(dp), parameter                  :: rho_ice   = 917.0     ! pure ice density [kg m-3]
  real(dp)                             :: dt                    ! time step [s] (ISSM_DT set in AGCM.rc)

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

  call MAPL_Get(MAPL,INTERNAL_ESMF_STATE = INTERNAL,RC=STATUS )
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
    call MAPL_GetResource (MAPL, dt, Label=trim(COMP_NAME)//"_DT:",DEFAULT=302400.0, RC=STATUS)
    VERIFY_(STATUS)

    ! get number of mesh elements
    call ESMF_MeshGet(mesh,elementCount=num_elements,nodeCount=num_nodes,numOwnedNodes=num_owned_nodes)
    allocate(elementConn(3*num_elements))
    call ESMF_MeshGet(mesh,elementConn=elementConn)

    ! allocate ice-elevation output (export from ISSM)
    allocate(ISSM_OUTPUTS(num_outputs*num_nodes))  

    ! allocate output arrays defined on mesh nodes
    allocate(ICESURF_MESH(num_nodes))
    allocate(ICETHICK_MESH(num_nodes))
    allocate(ICEVX_MESH(num_nodes))
    allocate(ICEVY_MESH(num_nodes))
    allocate(ICEVEL_MESH(num_nodes))
    allocate(IMLS_MESH(num_nodes))
    allocate(OMLS_MESH(num_nodes))
    
    ! allocate input arrays defined on mesh nodes
    allocate(ICESMB_MESH(num_nodes))  

    ! initialize ISSM outputs to zero 
    ICESURF_MESH(:) = 0.0_dp
    ICETHICK_MESH(:) = 0.0_dp
    ICEVX_MESH(:) = 0.0_dp
    ICEVY_MESH(:) = 0.0_dp
    ICEVEL_MESH(:) = 0.0_dp
    IMLS_MESH(:) = 0.0_dp
    OMLS_MESH(:) = 0.0_dp
    ISSM_OUTPUTS(:) = 0.0_dp

    ! get internal state for regridding
    call ESMF_UserCompGetInternalState(GC, 'ISSM_WRAP', wrap, status); VERIFY_(STATUS)
    internal_state => wrap%ptr
    routehandle_m2g = internal_state%routehandle_m2g  ! mesh to grid 
    routehandle_g2m = internal_state%routehandle_g2m  ! grid to mesh

    ! get field halo information
    num_halo_nodes = num_nodes - num_owned_nodes
    allocate(halo_idx(num_halo_nodes))
    allocate(owned_idx(num_owned_nodes))
    allocate(halolist(num_halo_nodes))
    halo_idx = internal_state%halo_idx
    owned_idx = internal_state%owned_idx
    halohandle = internal_state%halohandle
    halolist = internal_state%halolist
    nodalDistgrid = internal_state%nodalDistgrid

    ! get landice tile dimensions
    call MAPL_LocStreamGet(internal_state%locstream, NT_LOCAL=NT, _RC)

    ! get grid dimensions, used below in regridding subroutines
    call MAPL_GridGet(internal_state%grid, localCellCountPerDim=local_dims, _RC)
    IM = local_dims(1)
    JM = local_dims(2)
     
    call ESMF_UserCompGetInternalState(GC, 'ISSM_TILES', issm_tile_wrap, status); VERIFY_(STATUS)
    issm_tile_state => issm_tile_wrap%ptr
    
    ! *************************************************************************** !
    ! GET ICESMB IMPORT (surface mass balance)
    ! *************************************************************************** !    
    ! allocate tiles for ICESMB 
    if(.not.associated(ICESMB_TILE)) then
      allocate(ICESMB_TILE(NT), STAT=STATUS)
      VERIFY_(STATUS)
      ICESMB_TILE = MAPL_Undef
    end if
    
    ! copy import values into tile array 
    ICESMB_TILE = issm_tile_state%ICESMB_TILE

    ! *************************************************************************** !
    ! TRANSFORM ICESMB FROM TILES TO MESH 
    ! *************************************************************************** !
    call tile_to_mesh(ICESMB_TILE,ICESMB_MESH)

    ! save ICESMB on mesh elements 
    call MAPL_GetPointer(EXPORT  , ICESMB_EX , 'ICESMB' , RC=STATUS); VERIFY_(STATUS)

    if(associated(ICESMB_EX)) ICESMB_EX = ICESMB_MESH(owned_idx)

    ! *************************************************************************** !
    !  RUN ISSM WITH SMB INPUT AND ICE-ELEVATION OUTPUT
    ! *************************************************************************** !
    ! convert SMB to units of [m/s] (ice-equivalent) before passing to ISSM
    ICESMB_MESH = ICESMB_MESH/rho_ice

    call ESMF_VMBarrier(vm, rc=status); VERIFY_(STATUS)

    ! call run method from ISSM library 
    call RunISSM(dt, c_loc(ICESMB_MESH), c_loc(ISSM_OUTPUTS), c_loc(elementConn))

    call ESMF_VMBarrier(vm, rc=status); VERIFY_(STATUS)

    ! *************************************************************************** !
    ! UNPACK AND EXPORT ISSM OUTPUTS ON MESH TILES
    ! *************************************************************************** !
    ! unpack ISSM output pointer
    ICESURF_MESH(:) = ISSM_OUTPUTS(1:num_nodes)
    ICETHICK_MESH(:) = ISSM_OUTPUTS(num_nodes+1:2*num_nodes)
    ICEVX_MESH(:) = ISSM_OUTPUTS(2*num_nodes+1:3*num_nodes)
    ICEVY_MESH(:) = ISSM_OUTPUTS(3*num_nodes+1:4*num_nodes)
    OMLS_MESH(:) = ISSM_OUTPUTS(4*num_nodes+1:5*num_nodes)
    IMLS_MESH(:) = ISSM_OUTPUTS(5*num_nodes+1:6*num_nodes)

    ! calculate ice flow speed
    ICEVEL_MESH = sqrt(ICEVX_MESH**2 + ICEVY_MESH**2)

    ! set pointers to tile-mesh exports
    call MAPL_GetPointer(EXPORT, ICESURF_EX, 'ICESURF', RC=STATUS); VERIFY_(STATUS)
    if(associated(ICESURF_EX)) ICESURF_EX = ICESURF_MESH(owned_idx)

    call MAPL_GetPointer(EXPORT, ICEVX_EX, 'ICEVX', RC=STATUS); VERIFY_(STATUS)
    if(associated(ICEVX_EX)) ICEVX_EX = ICEVX_MESH(owned_idx)

    call MAPL_GetPointer(EXPORT, ICEVY_EX, 'ICEVY', RC=STATUS); VERIFY_(STATUS)
    if(associated(ICEVY_EX)) ICEVY_EX = ICEVY_MESH(owned_idx)

    call MAPL_GetPointer(EXPORT, ICETHICK_EX, 'ICETHICK', RC=STATUS); VERIFY_(STATUS)
    if(associated(ICETHICK_EX)) ICETHICK_EX = ICETHICK_MESH(owned_idx)

    ! set pointers to tile-mesh internals
    call MAPL_GetPointer(INTERNAL, ICESURF_IN, 'ICESURF', RC=STATUS); VERIFY_(STATUS)
    if(associated(ICESURF_IN)) ICESURF_IN = ICESURF_MESH(owned_idx)

    call MAPL_GetPointer(INTERNAL, ICETHICK_IN, 'ICETHICK', RC=STATUS); VERIFY_(STATUS)
    if(associated(ICETHICK_IN)) ICETHICK_IN = ICETHICK_MESH(owned_idx)

    call MAPL_GetPointer(INTERNAL, OMLS_IN, 'OMLS', RC=STATUS); VERIFY_(STATUS)
    if(associated(OMLS_IN)) OMLS_IN = OMLS_MESH(owned_idx)

    call MAPL_GetPointer(INTERNAL, IMLS_IN, 'IMLS', RC=STATUS); VERIFY_(STATUS)
    if(associated(IMLS_IN)) IMLS_IN = IMLS_MESH(owned_idx)

    ! *************************************************************************** !
    ! REGRID MESH FIELDS ONTO LANDICE TILES AND EXPORT VIA INTERNAL STATE
    ! *************************************************************************** !
    ! transform from mesh to tiles
    call mesh_to_tile(ICESURF_MESH,ICESURF_TILE)
    issm_tile_state%ICESURF_TILE = ICESURF_TILE

    call mesh_to_tile(ICETHICK_MESH,ICETHICK_TILE)
    issm_tile_state%ICETHICK_TILE = ICETHICK_TILE

    call mesh_to_tile(ICEVEL_MESH,ICEVEL_TILE)
    issm_tile_state%ICEVEL_TILE = ICEVEL_TILE

  end if 

  
  ! barrier to ensure regridding completes before any deallocates
  call ESMF_VMBarrier(vm, rc=status)
  VERIFY_(STATUS)

  ! deallocates
  if(associated(ICESURF_MESH))  deallocate(ICESURF_MESH)
  if(associated(ICETHICK_MESH)) deallocate(ICETHICK_MESH)
  if(associated(ICEVEL_MESH))   deallocate(ICEVEL_MESH)
  if(associated(ICEVX_MESH))    deallocate(ICEVX_MESH)
  if(associated(ICEVY_MESH))    deallocate(ICEVY_MESH)
  if(associated(IMLS_MESH))     deallocate(IMLS_MESH)
  if(associated(OMLS_MESH))     deallocate(OMLS_MESH)  
  if(associated(ICESMB_MESH))   deallocate(ICESMB_MESH) 
  if(associated(ISSM_OUTPUTS))  deallocate(ISSM_OUTPUTS)
  if(associated(ICESMB_TILE))   deallocate(ICESMB_TILE)
  if(associated(ICESURF_TILE))  deallocate(ICESURF_TILE)
  if(associated(ICETHICK_TILE)) deallocate(ICETHICK_TILE)
  if(associated(ICEVEL_TILE))   deallocate(ICEVEL_TILE)
  if(associated(elementConn))   deallocate(elementConn)
  if(associated(halolist))      deallocate(halolist)
  if(associated(halo_idx))      deallocate(halo_idx)
  if(associated(owned_idx))     deallocate(owned_idx)

  call MAPL_TimerOff(MAPL,"RUN"  )
  call MAPL_TimerOff(MAPL,"TOTAL")
 
  RETURN_(ESMF_SUCCESS)

  contains

  subroutine mesh_to_tile(VAR_MESH,VAR_TILE)
    ! regrid from mesh to grid, then transform from grid to landice tiles
    real(dp),    pointer, dimension(:), intent(inout)   :: VAR_MESH           ! var on mesh nodes
    real, pointer, dimension(:), intent(inout)          :: VAR_TILE           ! var on landice tiles
    real, pointer, dimension(:,:)                       :: VAR_GRID => null() ! var on attached grid
    real(dp),    pointer, dimension(:)                  :: VAR_MESH_OWN       ! var on owned mesh nodes

    allocate(VAR_MESH_OWN(num_owned_nodes))

    VAR_MESH_OWN = VAR_MESH(owned_idx)

    ! allocate tiles 
    if (.not.associated(VAR_TILE)) then
      allocate(VAR_TILE(NT), STAT=STATUS)
      VERIFY_(STATUS)
      VAR_TILE = MAPL_Undef
    end if

    ! create source field: field on mesh elements
    srcField = ESMF_FieldCreate(mesh=mesh,farrayPtr=VAR_MESH_OWN,meshloc=ESMF_MESHLOC_NODE, & 
    datacopyflag=ESMF_DATACOPY_VALUE,rc=STATUS)
    VERIFY_(STATUS)
    
    ! create destination field: field on grid
    dstField = ESMF_FieldCreate(grid=internal_state%grid,typekind=ESMF_TYPEKIND_R4,rc=STATUS)
    VERIFY_(STATUS)

    ! regrid field from mesh to grid
    call ESMF_FieldRegrid(srcField, dstField, routehandle_m2g, RC=STATUS); VERIFY_(STATUS)

    ! get pointer to field on grid
    call ESMF_FieldGet(dstField,farrayPtr=VAR_GRID,RC=STATUS); VERIFY_(STATUS)

    ! transform from grid to tiles  
    call MAPL_LocStreamTransform(internal_state%locstream,VAR_TILE,VAR_GRID, RC=STATUS)
    VERIFY_(STATUS)

    ! destroy regridding fields so they can be reused
    call ESMF_FieldDestroy(srcField,rc=STATUS); VERIFY_(STATUS)
    call ESMF_FieldDestroy(dstField,rc=STATUS); VERIFY_(STATUS)

  end subroutine mesh_to_tile

  subroutine tile_to_mesh(VAR_TILE,VAR_MESH)
    ! transform from landice tile to grid, then regrid onto mesh
    real, pointer, dimension(:), intent(inout)     :: VAR_TILE           ! var on landice tiles
    real(dp), pointer, dimension(:), intent(inout) :: VAR_MESH           ! var on mesh elements
    real, pointer, dimension(:,:)                  :: VAR_GRID => null() ! var on attached grid
    real(dp), pointer, dimension(:)                :: MESH_PTR           ! pointer for ESMF_FieldGet 

    ! allocate pointer on grid for regridding 
    allocate(VAR_GRID(IM,JM), STAT=STATUS); VERIFY_(STATUS)
  
    ! transform from tile to grid
    ! NOTE: we use the "transpose" option with MAPL_LocStreamTransformG2T 
    ! (rather than MAPL_LocStreamTransformT2G) because the "default" value is zero
    ! (rather than MAPL_UNDEF, which leads to errors when regridding onto mesh)
    call MAPL_LocStreamTransform(internal_state%LOCSTREAM, VAR_TILE, VAR_GRID, TRANSPOSE=.true., RC=STATUS)
    VERIFY_(STATUS)

    ! create source field on grid
    srcField = ESMF_FieldCreate(grid=internal_state%grid,farrayPtr=VAR_GRID, datacopyflag=ESMF_DATACOPY_VALUE,rc=STATUS)
    VERIFY_(STATUS)

    ! create destination field on mesh elements
    meshArray=ESMF_ArrayCreate(nodalDistgrid,typekind=ESMF_TYPEKIND_R8,haloSeqIndexList=halolist,rc=STATUS)
    VERIFY_(STATUS)

    ! create field on ISSM mesh
    dstField=ESMF_FieldCreate(mesh, array=meshArray, meshLoc=ESMF_MESHLOC_NODE, rc=STATUS)
    VERIFY_(STATUS)

    ! regrid from grid to mesh
    call ESMF_FieldRegrid(srcField, dstField, routehandle_g2m, RC=STATUS); VERIFY_(STATUS)

    ! append halo values to end of "owned" array
    call ESMF_FieldHalo(dstField, routehandle=halohandle, rc=STATUS); VERIFY_(STATUS)

    ! get pointer to field on mesh
    call ESMF_FieldGet(dstField,farrayPtr=MESH_PTR,RC=STATUS); VERIFY_(STATUS)

    ! copy values into VAR_MESH
    VAR_MESH(owned_idx) = MESH_PTR(1:num_owned_nodes)          ! owned nodes
    VAR_MESH(halo_idx) = MESH_PTR(num_owned_nodes+1:num_nodes) ! halo nodes
    
    ! destroy fields and arrays so they can be reused
    deallocate(VAR_GRID)
    call ESMF_FieldDestroy(srcField,rc=STATUS);  VERIFY_(STATUS)
    call ESMF_FieldDestroy(dstField,rc=STATUS);  VERIFY_(STATUS)
    call ESMF_ArrayDestroy(meshArray,rc=STATUS); VERIFY_(STATUS)

  end subroutine tile_to_mesh  

 
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
    integer			                       :: STATUS
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
