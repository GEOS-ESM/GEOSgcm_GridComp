!  $Id$

#include "MAPL_Generic.h"

module GEOS_IssmGridCompMod

!BOP
! !MODULE: GEOS_ISSM --- Runs ISSM (Ice-sheet and Sea-level System Model)
! 
!
! !DESCRIPTION:
!
!   {\tt GEOS\_ISSM} runs ISSM (Ice-sheet and Sea-level System Model)
!   Imports: ICESMB (defined on landice tiles)  [via private internal state]
!   Exports: ICESURF, ICETHICK, ICESMB_ISSM, ICEVX, ICEVY (defined on mesh) [true export state]
!   Exports: ICESURF, ICETHICK, ICEVEL (defined on landice tiles) [via private internal state]
!   Internals: ICESURF, ICETHICK, IMLS, OMLS, ISSM_NSTEPS (defined on mesh) [true internal state]
! *** NOTES: 
!            (*) currently we run over all input files (*.bin) that are found in ISSM_EXPDIR (scratch directory)
!                (e.g., Greenland + Antarctica + any other glaciers that have been configured)
!            (*) ISSM meshes are internal to ISSM (C++ source)--we create an ESMF_MESH version for regridding     
!                imports/exports that is the global combination of all ISSM meshes
!            (*) we transform imports from landice tiles to attached grid, then regrid to the mesh 
!            (*) ISSM outputs are saved with HISTORY via a 'mesh tile space' developed by Weiyuan Jiang (GMAO SI Team)  
!            (*) ISSM time step is generally larger than LANDICE timestep, or even a job duration. We persist ISSM 
!                variables across job segments through internal state checkpoints (restarts). We make sure that INTERNAL 
!                and EXPORT variables are 'filled in' by Initialize so that LANDICE and HISTORY have access to ISSM 
!                variables before it runs.  
!            (*) Related, we use a custom ISSM run alarm that is keyed to the last time ISSM ran, not the simulation 
!                start time. The number of LANDICE time steps since ISSM last ran is tracked via the internal state.  

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
    
subroutine RunISSM(ISSM_DT, gcm_forcings, issm_outputs) bind(C,NAME="RunISSM")
   import :: c_ptr, c_double
   real(c_double),   value        :: ISSM_DT
   type(c_ptr),      value        :: gcm_forcings
   type(c_ptr),      value        :: issm_outputs
end subroutine RunISSM

subroutine InputFromRestarts(gcm_restarts) bind(C,NAME="InputFromRestarts")
  import :: c_ptr
  type(c_ptr),       value        :: gcm_restarts
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
! produces binary output files from ISSM
! not currently used because we use GEOS restarts, and this
! will throw an error if ISSM has not run during a job segment.
end subroutine FinalizeISSM

end interface

private

public SetServices

! some shared derived types and parameters below:       

public :: T_ISSM_TILE_STATE
public :: ISSM_TILE_WRAP
! define ISSM export as internal variables, will be used by the landice gridcomp

type T_ISSM_TILE_STATE
    real, pointer :: ICESURF_TILE(:)
    real, pointer :: ICETHICK_TILE(:)
    real, pointer :: ICEVEL_TILE(:)
    real, pointer :: ICESMB_ISSM(:)
    integer       :: ISSM_NSTEPS
    real          :: LANDICE_DT
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
  type(ESMF_MESH)               :: mesh            ! ISSM mesh
  type(MAPL_LocStream)          :: locstream       ! original locstream (landice tiles)
end type T_ISSM_STATE

! Wrapper for extracting internal state
! -------------------------------------
type ISSM_WRAP
  type (T_ISSM_STATE), pointer :: ptr
end type ISSM_WRAP

integer                      :: num_outputs = 6          ! number of output fields that ISSM sends to GEOS
logical                      :: ISSM_RST_FOUND = .false. ! restart found flag
type(T_ISSM_STATE), pointer  :: internal_state=>null()   ! internal state for regridding and halo operations

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
!   library IRF methods.

!EOP

!=============================================================================

! ErrLog Variables

    character(len=ESMF_MAXSTR)         :: IAm
    integer                            :: STATUS
    character(len=ESMF_MAXSTR)         :: COMP_NAME

!=============================================================================

    type(MAPL_MetaComp), pointer       :: MAPL

    ! Get my internal MAPL_Generic state

    ! Begin...

! Get my name and set-up traceback handle
! ---------------------------------------

    call ESMF_GridCompGet( GC, NAME=COMP_NAME, _RC )
    Iam = trim(COMP_NAME) // 'SetServices'

! Set the Initialize, Run, and Finalize entry points
!-----------------------------------

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_INITIALIZE,   Initialize, _RC)
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,          Run,        _RC)
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_FINALIZE,     Finalize,   _RC)
    
!-----------------------------------

    call MAPL_GetObjectFromGC (GC, MAPL, _RC)
    
! Set the state variable specs.
!-----------------------------------

!   Import states: ICESMB is imported via the ISSM_TILE private internal state

!   Export states:
    call MAPL_AddExportSpec(GC,                    &
         SHORT_NAME = 'ICESURF',                   &
         LONG_NAME  = 'ice_sheet_elevation',       &
         UNITS      = 'm',                         &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone,          &
         _RC  )
    
    call MAPL_AddExportSpec(GC,                    &
         SHORT_NAME = 'ICEVX',                     &
         LONG_NAME  = 'ice_velocity_x_direction',  &
         UNITS      = 'm s-1',                     &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone,          &
         _RC  )
    
    call MAPL_AddExportSpec(GC,                    &
        SHORT_NAME = 'ICEVY',                      &
        LONG_NAME  = 'ice_velocity_y_direction',   &
        UNITS      = 'm s-1',                      &
        DIMS       = MAPL_DimsTileOnly,            &
        VLOCATION  = MAPL_VLocationNone,           &
        _RC  )
    
    call MAPL_AddExportSpec(GC,                    &
         SHORT_NAME = 'ICETHICK',                  &
         LONG_NAME  = 'ice_thickness',             &
         UNITS      = 'm',                         &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone,          &
         _RC  )
		 
    call MAPL_AddExportSpec(GC,                    &
         SHORT_NAME = 'ICESMB_ISSM',               &
         LONG_NAME  = 'issm_surface_mass_balance', &
         UNITS      = 'kg m-2 s-1',                &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone,          &
         _RC  )	 
    
    !   Internal states:
    call MAPL_AddInternalSpec(GC,                  &
         SHORT_NAME = 'ICESURF',                   &
         LONG_NAME  = 'ice_sheet_elevation',       &
         UNITS      = 'm',                         &
         PRECISION  = ESMF_KIND_R8,                &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone,          &
         RESTART    = MAPL_RestartOptional,        &   
         _RC  )
    
    call MAPL_AddInternalSpec(GC,                  &
         SHORT_NAME = 'ICETHICK',                  &
         LONG_NAME  = 'ice_sheet_thickness',       &
         UNITS      = 'm',                         &
         PRECISION  = ESMF_KIND_R8,                &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone,          &
         RESTART    = MAPL_RestartOptional,        & 
         _RC  )
    
    call MAPL_AddInternalSpec(GC,                  &
         SHORT_NAME = 'IMLS',                      &
         LONG_NAME  = 'ice_mask_levelset',         &
         UNITS      = 'none',                      &
         PRECISION  = ESMF_KIND_R8,                &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone,          &
         RESTART    = MAPL_RestartOptional,        & 
         _RC  )
    
    call MAPL_AddInternalSpec(GC,                  &
         SHORT_NAME = 'OMLS',                      &
         LONG_NAME  = 'ocean_mask_levelset',       &
         UNITS      = 'none',                      &
         PRECISION  = ESMF_KIND_R8,                &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone,          &
         RESTART    = MAPL_RestartOptional,        & 
         _RC  )
    
    call MAPL_AddInternalSpec(GC,                  &
         SHORT_NAME = 'ICEVX',                     &
         LONG_NAME  = 'ice_velocity_x_direction',  &
         UNITS      = 'm s-1',                     &
         PRECISION  = ESMF_KIND_R8,                &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone,          &
         RESTART    = MAPL_RestartOptional,        &
         _RC  )
    
    call MAPL_AddInternalSpec(GC,                  &
         SHORT_NAME = 'ICEVY',                     &
         LONG_NAME  = 'ice_velocity_y_direction',  &
         UNITS      = 'm s-1',                     &
         PRECISION  = ESMF_KIND_R8,                &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone,          &
         RESTART    = MAPL_RestartOptional,        &
         _RC  )
    
    call MAPL_AddInternalSpec(GC,                  &
         SHORT_NAME = 'ISSM_NSTEPS',               &
         LONG_NAME  = 'steps_since_last_issm',     &
         UNITS      = 'none',                      &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone,          &
         RESTART    = MAPL_RestartOptional,        &
         _RC  )
	 call MAPL_AddInternalSpec(GC,                 &
         SHORT_NAME = 'RS_NODEIDS',                &
         LONG_NAME  = 'restart_node_ids',          &
         UNITS      = 'none',                      &
         DIMS       = MAPL_DimsTileOnly,           &
         VLOCATION  = MAPL_VLocationNone,          &
         RESTART    = MAPL_RestartOptional,        &
         _RC  )


! Set the Profiling timers
! ------------------------

    call MAPL_TimerAdd(GC,    name="RUN"   ,_RC)
    call MAPL_TimerAdd(GC,    name="ISSMCore"   ,_RC)
	
       
! ----------------------------------
    call MAPL_GenericSetServices    ( GC, _RC)
    
    _RETURN(_SUCCESS)
  
  end subroutine SetServices

  ! ! INITIALIZE:
  
  subroutine Initialize ( GC, IMPORT, EXPORT, CLOCK, RC )
    type(ESMF_GridComp),     intent(INOUT) :: GC                      ! Gridded component 
    type(ESMF_State),        intent(INOUT) :: IMPORT                  ! Import state
    type(ESMF_State),        intent(INOUT) :: EXPORT                  ! Export state
    type(ESMF_Clock),        intent(INOUT) :: CLOCK                   ! The clock
    integer, optional,       intent(OUT)   :: RC                      ! Error code
    
    type(MAPL_MetaComp), pointer           :: MAPL   
    type(ESMF_State)                       :: INTERNAL                ! internal state

    ! ISSM alarm variables
    type(ESMF_Alarm)                       :: ISSM_ALARM              ! custom ISSM RUNALARM
    integer                                :: sec_to_ring             ! seconds remaining until first ISSM run
    type(ESMF_Time)                        :: startTime               ! initial time
    type(ESMF_TimeInterval)                :: startInterval           ! time interval to first ring
    type(ESMF_Time)                        :: ringTime                ! time of first ring
    type(ESMF_TimeInterval)                :: ringInterval            ! ring time interval (ISSM_DT)
    real                                   :: ISSM_DT                 ! ISSM time step [s] (ISSM_DT set in AGCM.rc)
    real                                   :: LANDICE_DT              ! landice time step [s] 
    integer                                :: NSTEPS_INIT             ! landice timesteps since last ISSM run
    integer                                :: NSTEPS_RING             ! total landice timesteps between ISSM runs
    real, pointer, dimension(:)            :: ISSM_NSTEPS => null()   ! steps since last ISSM run (from internal state)
	
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
    
    ! regridding varibales
    type(ESMF_Grid)                        :: grid                    ! atmospheric grid
    type(ESMF_RouteHandle)                 :: routehandle_m2g         ! routehandle for regridding mesh to grid
    type(ESMF_RouteHandle)                 :: routehandle_g2m         ! routehandle for regridding grid to mesh
    type(ESMF_Field)                       :: meshField               ! field on mesh
    type(ESMF_Field)                       :: gridField               ! field on grid
    type(ISSM_WRAP)                        :: wrap                    ! wrapper for internal state

    ! tile information
    integer                                :: NT                      ! local number of landice tiles
    type(T_ISSM_TILE_STATE), pointer       :: issm_tile_state
    type(ISSM_TILE_WRAP)                   :: issm_tile_wrap

    ! field halo variables
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
    character(len=ESMF_MAXSTR)             :: ISSM_EXPDIR             ! directory containing ISSM input files
    character(len=ESMF_MAXSTR)             :: EXPDIR                  ! C++ compatible ISSM_EXPDIR string

    ! variables for creating mesh tile space
    type(ESMF_Grid)                        :: mesh_grid               
    type(MAPL_LocStream)                   :: mesh_locstream   

    ! variables for masking the mesh seam (triangles that cross +/-180 longitude)
    ! (needed for elements, this is not currently needed for regridding fields defined on nodes)
    real(dp)                               :: dlon,lon1,lon2,lon3
    integer, pointer, dimension(:)         :: elementMask => null()
    integer                                :: n1,n2,n3

    ! pointers to internal state for restarts
    real(dp), pointer, dimension(:)        :: ICESURF_IN    => null() ! ice surface elevation restart
    real(dp), pointer, dimension(:)        :: ICETHICK_IN   => null() ! ice thickness restart
    real(dp), pointer, dimension(:)        :: ICEVX_IN      => null() ! ice velocity (x direction) restart
    real(dp), pointer, dimension(:)        :: ICEVY_IN      => null() ! ice velocity (y direction) restart
    real(dp), pointer, dimension(:)        :: IMLS_IN       => null() ! ice-mask levelset restart
    real(dp), pointer, dimension(:)        :: OMLS_IN       => null() ! ocean-mask levelset restart

    ! restarts with halo points (interleaved), to send to ISSM
    real(dp), pointer, dimension(:)        :: ICESURF_HALO  => null()
    real(dp), pointer, dimension(:)        :: ICETHICK_HALO => null()
    real(dp), pointer, dimension(:)        :: IMLS_HALO     => null()
    real(dp), pointer, dimension(:)        :: OMLS_HALO     => null()
    real(dp), pointer, dimension(:)        :: ICEVX_HALO    => null()
    real(dp), pointer, dimension(:)        :: ICEVY_HALO    => null()
    real(dp), pointer, dimension(:)        :: ICEVEL_HALO   => null()

    real(dp), pointer, dimension(:)        :: GEOS_RESTARTS => null() ! concatenate restart fields
    real(dp), pointer, dimension(:)        :: ZEROS         => null() ! zero input for bootstrapping

    ! export variables on landice tile space
    real, pointer, dimension(:)            :: ICESURF_TILE  => null() ! ice surface elevation on landice tiles
    real, pointer, dimension(:)            :: ICETHICK_TILE => null() ! ice thickness on landice tiles
    real, pointer, dimension(:)            :: ICEVEL_TILE   => null() ! ice flow speed on landice tiles

    ! export variables on mesh tile space
    real, pointer, dimension(:)            :: ICESURF_EX    => null() ! ice surface elevation on mesh tiles
    real, pointer, dimension(:)            :: ICETHICK_EX   => null() ! ice thickness on mesh tiles
    real, pointer, dimension(:)            :: ICEVX_EX      => null() ! ice velocity (x direction) on mesh tiles
    real, pointer, dimension(:)            :: ICEVY_EX      => null() ! ice velocity (y direction) on mesh tiles

	! restart redistribution
    real, pointer, dimension(:)            :: restartNodeIds=> null() ! nodeIds for restart ordering
    type(ESMF_DistGrid)                    :: restartDistgrid         ! distgrid from reading restarts
    logical                                :: distgrid_match          ! check if distgrid from restarts matches nodal disgrid (locally)
    logical                                :: needRedist              ! global check for consistent distgrid across all processes
    integer,allocatable,dimension(:)       :: localFlag, globalFlag   ! arrays for vm operations
    type(ESMF_Array)                       :: restartArray            ! array corresponding to restartDistgrid 
    type(ESMF_Array)                       :: nodalArray              ! array corresponding to nodalDistgrid 
    type(ESMF_RouteHandle)                 :: redisthandle            ! routehandle for redistribution


    ! Get the target components name and set-up traceback handle.
    ! -----------------------------------------------------------

    Iam = "Initialize"
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, _RC )
    
    Iam = trim(COMP_NAME) // trim(Iam)

    ! Get my internal MAPL_Generic state
    !-----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, _RC)

    call ESMF_VMGetCurrent(vm, _RC)
    
    call ESMF_VMGet(vm,mpiCommunicator=comm,localPet=localPET,_RC)
    
    ! ****************************************************
    ! call ISSM initialize C++ code so we can set up mesh

    ! get directory with ISSM binary input files (can modify if needed)
    call GET_ENVIRONMENT_VARIABLE("SCRDIR",ISSM_EXPDIR,STATUS=STATUS); _VERIFY(STATUS)
    
    EXPDIR = trim(ISSM_EXPDIR)//"/"//c_null_char ! create string for C++
    
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
    allocate(nodeOwners(num_nodes))

    ! get information about nodes and elements
    ! node coords and element coords (centroids) are in (lon,lat)
    call GetNodesISSM(c_loc(nodeIds), c_loc(nodeCoords)) 
    call GetElementsISSM(c_loc(elementIds), c_loc(elementConn), c_loc(elementCoords),c_loc(glacIds))

    elementTypes(:) = ESMF_MESHELEMTYPE_TRI ! triangular elements

    ! mask for triangles that cross the seam (longitude +/- 180)
    ! (you don't have to 'activate' this mask, it can just be 'associated' with the mesh)
	!
    ! NOTE: This is only relevant in regridding when fields are defined
    !       on ESMF_MESHLOC_ELEMENT (rather than ESMF_MESHLOC_NODE)
	!       so is NOT CURRENTLY USED, but retained for possible future developments 
    elementMask(:) = 0
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
      end if
    end do

    ! create the ESMF mesh from ISSM mesh properties
    mesh = ESMF_MeshCreate(parametricDim=2, spatialDim=2, nodeIds=nodeIds, nodeCoords=nodeCoords, &
            elementIds=elementIds, elementTypes=elementTypes, elementConn=elementConn,elementMask=elementMask,& 
            elementCoords=elementCoords,coordSys=ESMF_COORDSYS_SPH_DEG, _RC)
    
    ! associate ESMF_Mesh representation of ISSM mesh with GC for regridding imports/exports in Run method       
    call ESMF_GridCompSet(GC,mesh=mesh,_RC)
    
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
    meshArray=ESMF_ArrayCreate(nodalDistgrid,typekind=ESMF_TYPEKIND_R8,haloSeqIndexList=halolist,_RC)
    
    ! create field on ISSM mesh 
    meshField=ESMF_FieldCreate(mesh, array=meshArray, meshLoc=ESMF_MESHLOC_NODE, _RC)
    
    ! store the halo operation in a routehandle
    call ESMF_FieldHaloStore(meshField, routehandle=halohandle, _RC)
    
    ! Set up regridding next
    !-----------------------------------
    ! get atmospheric (attached) grid 
    call ESMF_GridCompGet( GC, GRID=grid, _RC )
    
    ! create field on atmospheric grid
    gridField = ESMF_FieldCreate(grid=grid,typekind=ESMF_TYPEKIND_R4,_RC)
    
    ! create routehandle for mesh-to-grid regridding (set srcMaskValues to 1 if needed... )
    call ESMF_FieldRegridStore(srcField=meshField, dstField=gridField,routehandle=routehandle_m2g,& 
    unmappedaction=ESMF_UNMAPPEDACTION_IGNORE,extrapmethod=ESMF_EXTRAPMETHOD_CREEP,&
    extrapNumLevels=1,_RC)
    
    ! create routehandle for grid-to-mesh regridding (set dstMaskValues to 1 if needed... )
    call ESMF_FieldRegridStore(srcField=gridField, dstField=meshField,routehandle=routehandle_g2m,& 
    unmappedaction=ESMF_UNMAPPEDACTION_IGNORE,extrapmethod=ESMF_EXTRAPMETHOD_NEAREST_D,_RC)
    
    ! create component's private internal state
    ! stores everything needed for regrid and halo operations during run method
    allocate(internal_state, stat=STATUS); _VERIFY(STATUS)
    
    allocate(internal_state%halo_idx(num_halo_nodes))
    allocate(internal_state%owned_idx(num_owned_nodes))
    allocate(internal_state%halolist(num_halo_nodes))
    internal_state%routehandle_m2g = routehandle_m2g
    internal_state%routehandle_g2m = routehandle_g2m 
    internal_state%halohandle = halohandle 
    internal_state%halo_idx = halo_idx
    internal_state%owned_idx = owned_idx  
    internal_state%grid = grid
    internal_state%mesh = mesh
    internal_state%halolist = halolist
    internal_state%nodalDistgrid = nodalDistgrid
    call MAPL_Get(MAPL, LocStream = internal_state%locstream, _RC)

    ! wrap the private internal state
    wrap%ptr => internal_state
    call ESMF_UserCompSetInternalState ( GC, 'ISSM_WRAP', wrap, STATUS ); _VERIFY(STATUS)
    
    ! Create losctream that match mesh element id, then set it to this GC and MAPL
    ! note: original attached/atmospheric grid and landice tile locstream have
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
    call MAPL_Set(MAPL, locstream = mesh_locstream, _RC)

    ! Generic initialize
    !-----------------------------------

    call MAPL_GenericInitialize( GC, IMPORT, EXPORT, CLOCK, _RC )

    ! Get private internal state for sending information to/from LANDICE
	!-----------------------------------
	
    call ESMF_UserCompGetInternalState(GC, 'ISSM_TILES', issm_tile_wrap, status); _VERIFY(STATUS)
    issm_tile_state => issm_tile_wrap%ptr

	! Create Custom ISSM Run Alarm 
    !-----------------------------------

    ! get internal state
    call MAPL_Get(MAPL,INTERNAL_ESMF_STATE = INTERNAL,_RC)

	! get number of time steps since last ISSM run
    call MAPL_GetPointer(INTERNAL, ISSM_NSTEPS, 'ISSM_NSTEPS',_RC)
    NSTEPS_INIT = nint(maxval(ISSM_NSTEPS))

    ! get timestep for landice 
    LANDICE_DT = issm_tile_state%LANDICE_DT
    
    ! get timestep for ISSM
    call MAPL_GetResource(MAPL, ISSM_DT, Label=trim(COMP_NAME)//"_DT:",DEFAULT=302400.0, _RC)
    
    ! total landice time steps between ISSM runs
    NSTEPS_RING = nint(ISSM_DT/LANDICE_DT) 
	
    ! calculate initial ring time from initial time and remaining timesteps
    call ESMF_ClockGet(CLOCK,currTime=startTime)
    sec_to_ring = (NSTEPS_RING-NSTEPS_INIT-1)*nint(LANDICE_DT)
    call ESMF_TimeIntervalSet(startInterval,s = sec_to_ring )
    ringTime = startTime + startInterval

    ! set ring interval to ISSM time step
    call ESMF_TimeIntervalSet(ringInterval,s=nint(ISSM_DT),_RC)
    
    ! create new ISSM_ALARM
    ISSM_ALARM = ESMF_AlarmCreate(CLOCK,ringTime=ringTime,ringInterval=ringInterval,sticky=.false.,_RC)
	  
	! set run alarm
    call MAPL_Set(MAPL, RUNALARM = ISSM_ALARM, _RC)

    ! Next, send GEOS restarts to ISSM
    !-----------------------------------
    ! array holding all restarts to send to/from ISSM
    allocate(GEOS_RESTARTS(num_outputs*num_nodes))
    allocate(ICESURF_HALO(num_nodes))
    allocate(ICETHICK_HALO(num_nodes))
    allocate(ICEVX_HALO(num_nodes))
    allocate(ICEVY_HALO(num_nodes))
    allocate(ICEVEL_HALO(num_nodes))
    allocate(IMLS_HALO(num_nodes))
    allocate(OMLS_HALO(num_nodes))
    allocate(ZEROS(num_nodes))
    
    ! get pointers to restarts
    call MAPL_GetPointer(INTERNAL, ICESURF_IN, 'ICESURF', _RC)
    call MAPL_GetPointer(INTERNAL, ICETHICK_IN, 'ICETHICK',_RC)
    call MAPL_GetPointer(INTERNAL, ICEVX_IN, 'ICEVX',_RC)
    call MAPL_GetPointer(INTERNAL, ICEVY_IN, 'ICEVY',_RC)
    call MAPL_GetPointer(INTERNAL, IMLS_IN, 'IMLS', _RC)
    call MAPL_GetPointer(INTERNAL, OMLS_IN, 'OMLS',_RC)
    call MAPL_GetPointer(INTERNAL, restartNodeIds, 'RS_NODEIDS',_RC)
	
    ! if restart has been read, apply halo operation and send pointers to ISSM
    ! else, ISSM will just use default initial values in ISSM*.bin input files
    if (associated(ICETHICK_IN)) then
        ! simple check for positive ice thickness (initialized to zero if restart not found)
        ! ISSM throws error for zero ice thickness
        ISSM_RST_FOUND = minval(ICETHICK_IN) > epsilon(ICETHICK_IN)
    end if
    
    if (ISSM_RST_FOUND) then
	  ! check if the nodal distgrid created above matches the distgrid read from the restart
	  ! it will only be different if running over a different number of processes than when
	  ! the restart was written. if it is, we redistribute the restart arrays correctly
      allocate(localFlag(1))
      allocate(globalFlag(1))
	  distgrid_match = all(ownedNodeIds==nint(restartNodeIds))
      localFlag(1) = 0
      if (distgrid_match) localFlag(1) = 1
      call ESMF_VMAllReduce(vm, sendData=localFlag, recvData=globalFlag, count=1, reduceflag=ESMF_REDUCE_MIN, _RC)
      needRedist = (globalFlag(1) == 0)
	  
      if (needRedist) then
	  ! create routehandle for redistribution, and redistribute all restarts from the 
	  ! restart distgrid to the current distgrid (nodalDistgrid) 
          restartDistgrid = ESMF_DistGridCreate(arbSeqIndexList=nint(restartNodeIds), _RC)
          restartArray=ESMF_ArrayCreate(distgrid=restartDistgrid,typekind=ESMF_TYPEKIND_R8,_RC)
          nodalArray=ESMF_ArrayCreate(distgrid=nodalDistgrid,typekind=ESMF_TYPEKIND_R8,_RC)
          call ESMF_ArrayRedistStore(srcArray=restartArray, dstArray=nodalArray, routehandle=redisthandle,_RC)

		  call apply_redist(ICESURF_IN,_RC)
		  call apply_redist(ICETHICK_IN,_RC)
		  call apply_redist(ICEVX_IN,_RC)
		  call apply_redist(ICEVY_IN,_RC)
		  call apply_redist(IMLS_IN,_RC)
		  call apply_redist(OMLS_IN,_RC)

          call ESMF_VMBarrier(vm, _RC)   
		  
		  call ESMF_ArrayDestroy(restartArray, _RC)
          call ESMF_ArrayDestroy(nodalArray, _RC)
		  
      end if 
	
      ! apply halo operation to all restart variables
      call apply_halo(ICESURF_IN,ICESURF_HALO,_RC)
      call apply_halo(ICETHICK_IN,ICETHICK_HALO,_RC)
      call apply_halo(ICEVX_IN,ICEVX_HALO,_RC)
      call apply_halo(ICEVY_IN,ICEVY_HALO,_RC)
      call apply_halo(IMLS_IN,IMLS_HALO,_RC)
      call apply_halo(OMLS_IN,OMLS_HALO,_RC)

      ! package restarts into one pointer
      GEOS_RESTARTS(:) = 0.0_dp
      GEOS_RESTARTS(1:num_nodes) = ICESURF_HALO(:) 
      GEOS_RESTARTS(num_nodes+1:2*num_nodes) = ICETHICK_HALO(:) 
      GEOS_RESTARTS(2*num_nodes+1:3*num_nodes) = ICEVX_HALO(:) 
      GEOS_RESTARTS(3*num_nodes+1:4*num_nodes) = ICEVY_HALO(:) 
      GEOS_RESTARTS(4*num_nodes+1:5*num_nodes) = OMLS_HALO(:) 
      GEOS_RESTARTS(5*num_nodes+1:6*num_nodes) = IMLS_HALO(:) 

      ! set restarts on the ISSM side
      call ESMF_VMBarrier(vm, _RC)      
      call InputFromRestarts(c_loc(GEOS_RESTARTS))
      call ESMF_VMBarrier(vm, _RC)

    else
      ! bootstrap restart values from ISSM input files (ISSM*.bin)
      ! by running with 'fake' time step with zero forcing

      call ESMF_VMBarrier(vm, _RC)
      call RunISSM(real(ISSM_DT,kind=dp), c_loc(ZEROS), c_loc(GEOS_RESTARTS))
      call ESMF_VMBarrier(vm, _RC)

      ! Unpack restart array 
      ICESURF_HALO(:)  = GEOS_RESTARTS(1:num_nodes) 
      ICETHICK_HALO(:) = GEOS_RESTARTS(num_nodes+1:2*num_nodes) 
      ICEVX_HALO(:) = GEOS_RESTARTS(2*num_nodes+1:3*num_nodes)
      ICEVY_HALO(:) = GEOS_RESTARTS(3*num_nodes+1:4*num_nodes)
      OMLS_HALO(:)  = GEOS_RESTARTS(4*num_nodes+1:5*num_nodes) 
      IMLS_HALO(:)  = GEOS_RESTARTS(5*num_nodes+1:6*num_nodes) 

      ! filter out halo points (keep the owned indices) for restarts
      if(associated(ICESURF_IN)) ICESURF_IN = ICESURF_HALO(owned_idx)
      if(associated(ICETHICK_IN)) ICETHICK_IN = ICETHICK_HALO(owned_idx)
      if(associated(ICEVX_IN)) ICEVX_IN = ICEVX_HALO(owned_idx)
      if(associated(ICEVY_IN)) ICEVY_IN = ICEVY_HALO(owned_idx)
      if(associated(OMLS_IN)) OMLS_IN = OMLS_HALO(owned_idx)
      if(associated(IMLS_IN)) IMLS_IN = IMLS_HALO(owned_idx)

    end if 

    ! Initialize Export pointers on mesh tile space so history has something to write
    !-----------------------------------

    call MAPL_GetPointer(EXPORT, ICESURF_EX, 'ICESURF',alloc=.true., _RC)
    if(associated(ICESURF_EX)) ICESURF_EX = ICESURF_HALO(owned_idx)

    call MAPL_GetPointer(EXPORT, ICETHICK_EX, 'ICETHICK',alloc=.true., _RC)
    if(associated(ICETHICK_EX)) ICETHICK_EX = ICETHICK_HALO(owned_idx)

    call MAPL_GetPointer(EXPORT, ICEVX_EX, 'ICEVX',alloc=.true.,_RC)
    if(associated(ICEVX_EX)) ICEVX_EX = ICEVX_HALO(owned_idx)

    call MAPL_GetPointer(EXPORT, ICEVY_EX, 'ICEVY',alloc=.true.,_RC)
    if(associated(ICEVY_EX)) ICEVY_EX = ICEVY_HALO(owned_idx)

    ! Finally, set the tile export state so landice can access values before ISSM runs
    !-----------------------------------
	
    ! Regrid from mesh to tile
    call MAPL_LocStreamGet(internal_state%locstream, NT_LOCAL=NT, _RC)

    ! allocate variables on landice tile space
    allocate(ICESURF_TILE(NT))
    allocate(ICETHICK_TILE(NT))
    allocate(ICEVEL_TILE(NT))

    ! calculate ice flow speed
    ICEVEL_HALO = sqrt(ICEVX_HALO**2 + ICEVY_HALO**2)

    call mesh_to_tile(ICESURF_HALO,ICESURF_TILE,_RC)
    issm_tile_state%ICESURF_TILE = ICESURF_TILE

    call mesh_to_tile(ICETHICK_HALO,ICETHICK_TILE,_RC)
    issm_tile_state%ICETHICK_TILE = ICETHICK_TILE

    call mesh_to_tile(ICEVEL_HALO,ICEVEL_TILE,_RC)
    issm_tile_state%ICEVEL_TILE = ICEVEL_TILE

    issm_tile_state%ISSM_NSTEPS = NSTEPS_INIT


	! set nodeIds internal associated with restart
	if(associated(restartNodeIds)) restartNodeIds(:) = ownedNodeIds(:)

    call ESMF_VMBarrier(vm, _RC)

    ! deallocate pointers
    if(associated(nodeCoords))      deallocate(nodeCoords)
    if(associated(nodeIds))         deallocate(nodeIds)
    if(associated(elementTypes))    deallocate(elementTypes)
    if(associated(elementIds))      deallocate(elementIds)
    if(associated(elementConn))     deallocate(elementConn)
    if(associated(elementCoords))   deallocate(elementCoords)
    if(associated(glacIds))         deallocate(glacIds)
    if(associated(elementMask))     deallocate(elementMask)
    if(associated(halo_idx))        deallocate(halo_idx)
    if(associated(owned_idx))       deallocate(owned_idx)
    if(associated(halolist))        deallocate(halolist)
    if(associated(ownedNodeCoords)) deallocate(ownedNodeCoords)
    if(associated(ownedNodeIds))    deallocate(ownedNodeIds)
    if(associated(nodeOwners))      deallocate(nodeOwners)
    if(associated(ICESURF_HALO))    deallocate(ICESURF_HALO)
    if(associated(ICETHICK_HALO))   deallocate(ICETHICK_HALO)
    if(associated(IMLS_HALO))       deallocate(IMLS_HALO)
    if(associated(OMLS_HALO))       deallocate(OMLS_HALO)
    if(associated(ICEVX_HALO))      deallocate(ICEVX_HALO)
    if(associated(ICEVY_HALO))      deallocate(ICEVY_HALO)
    if(associated(GEOS_RESTARTS))   deallocate(GEOS_RESTARTS)
    if(associated(ZEROS))           deallocate(ZEROS)
    if(associated(ICESURF_TILE))    deallocate(ICESURF_TILE)
    if(associated(ICETHICK_TILE))   deallocate(ICETHICK_TILE)
    if(associated(ICEVEL_TILE))     deallocate(ICEVEL_TILE)

    ! destroy fields and arrays
    call ESMF_FieldDestroy(gridField, _RC)
    call ESMF_FieldDestroy(meshField, _RC)
    call ESMF_ArrayDestroy(meshArray, _RC)
    
    _RETURN(_SUCCESS)

    contains
       subroutine apply_halo(VAR_IN,VAR_HALO,RC)
          ! apply halo operation to a restart variable
          ! arguments:
          real(dp), pointer, dimension(:), intent(inout) :: VAR_IN            ! var on owned_nodes
          real(dp), pointer, dimension(:), intent(inout) :: VAR_HALO          ! var on all nodes
          integer, optional, intent(out)                 :: RC

		      ! local variables:
          real(dp), pointer, dimension(:)                :: VAR_DP            ! double version of VAR_IN
          real(dp), pointer, dimension(:)                :: MESH_PTR          ! pointer for ESMF_FieldGet
          real(dp), pointer, dimension(:)                :: ARRAY_PTR         ! pointer for ESMF_ArrayGet 
          type(ESMF_Array)                               :: meshArray         ! array for creating mesh fields
          type(ESMF_Field)                               :: meshField         ! field associated with meshArray

          allocate(VAR_DP(num_nodes))
          VAR_DP(:) = 0.0_dp
          VAR_DP(1:num_owned_nodes) = REAL(VAR_IN,kind=dp)
      
          ! create array with halo information
          meshArray=ESMF_ArrayCreate(nodalDistgrid,typekind=ESMF_TYPEKIND_R8,haloSeqIndexList=halolist,_RC)

          call ESMF_ArrayGet(array=meshArray,farrayPtr=ARRAY_PTR)
          ARRAY_PTR(:) = VAR_DP(:)

          ! create field on ISSM mesh 
          meshField=ESMF_FieldCreate(mesh, array=meshArray, meshLoc=ESMF_MESHLOC_NODE, _RC)
          
          ! append halo values to end of "owned" array
          call ESMF_FieldHalo(meshField, routehandle=halohandle, _RC)
      
          ! get pointer to field on mesh
          call ESMF_FieldGet(meshField,farrayPtr=MESH_PTR,_RC)

          ! copy values into VAR_HALO, interleave according to owned and halo indices
          VAR_HALO(owned_idx) = MESH_PTR(1:num_owned_nodes)          ! owned nodes
          VAR_HALO(halo_idx) = MESH_PTR(num_owned_nodes+1:num_nodes) ! halo nodes

          ! destroy field and array, deallocate pointer
          call ESMF_FieldDestroy(meshField,_RC)  
          call ESMF_ArrayDestroy(meshArray,_RC)  
          deallocate(VAR_DP)

          _RETURN(_SUCCESS)
       end subroutine apply_halo

       subroutine apply_redist(VAR_RS,RC)
	      ! arguments:
          real(dp), pointer, dimension(:), intent(inout) :: VAR_RS       ! var from restsart
          integer, optional, intent(out)                 :: RC

          type(ESMF_Array)                               :: restartArray ! restart array
          type(ESMF_Array)                               :: redistArray  ! redistributed array
          real, pointer, dimension(:)                    :: redistPtr	  

          restartArray=ESMF_ArrayCreate(distgrid=restartDistgrid,farrayPtr=VAR_RS,_RC)
          redistArray=ESMF_ArrayCreate(distgrid=nodalDistgrid,typekind=ESMF_TYPEKIND_R8,_RC)

          ! redistribute the data
          call ESMF_ArrayRedist(srcArray=restartArray, dstArray=redistArray, routehandle=redisthandle,_RC)

          ! get the pointer to the data
		  call ESMF_ArrayGet(redistArray,farrayPtr=redistPtr)

          ! make sure all processes have finished redistribution
          call ESMF_VMBarrier(vm, _RC)
 
		  ! copy values into output
		  VAR_RS(:) = redistPtr(:)

          call ESMF_VMBarrier(vm, _RC)
		  call ESMF_ArrayDestroy(restartArray,_RC)
		  call ESMF_ArrayDestroy(redistArray,_RC)  

          _RETURN(_SUCCESS)
       end subroutine apply_redist

       function create_mesh_grid(rc) result(mesh_grid)
          type (ESMF_Grid) :: mesh_grid
          integer, optional, intent(out) :: RC
          integer :: status, nDEs, num(1)
          real(kind=8), pointer :: centers_lon(:,:)
          real(kind=8), pointer :: centers_lat(:,:)
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
          _VERIFY(STATUS)

          call ESMF_GridGetCoord(mesh_grid, coordDim=1, localDE=0, &
               staggerloc=ESMF_STAGGERLOC_CENTER, &
               farrayPtr=centers_lon, _RC)
          centers_lon(:,1) = ownedNodeLons 
          call ESMF_GridGetCoord(mesh_grid, coordDim=2, localDE=0, &
               staggerloc=ESMF_STAGGERLOC_CENTER, &
               farrayPtr=centers_lat, _RC)
          centers_lat(:,1) = ownedNodeLats 

          _RETURN(_SUCCESS)
       end function create_mesh_grid      

  end subroutine Initialize

  !BOP


  subroutine RUN ( GC, IMPORT, EXPORT, CLOCK, RC )
  ! ! ****** Run ISSM ice-sheet model ******
  ! !  the core C++ solvers and associated pre/post-processing of imports/exports
  ! !  are only performed at ISSM_DT intervals. However, the Run method is engaged
  ! !  at every landice timestep to ensure that ISSM restarts persist  
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

    ! internal state for regridding and halo operations
    type(ESMF_Mesh)                      :: mesh                    ! ESMF version of ISSM mesh
    integer                              :: num_nodes               ! number of nodes on PET

    ! tile information
    integer                              :: NT                      ! number of landice tiles
    type(T_ISSM_TILE_STATE), pointer     :: issm_tile_state
    type(ISSM_TILE_WRAP)                 :: issm_tile_wrap

    ! surface mass balance on mesh and landice tiles
	  ! note: SMB has been time-averaged between ISSM runs
    real(dp), pointer, dimension(:)      :: ICESMB_MESH   => null() ! surface mass balce on mesh elements
    real, pointer, dimension(:)          :: ICESMB_TILE   => null() ! surface mass balance on landice tiles
    real, pointer, dimension(:)          :: ICESMB_EX     => null() ! pointer to export state (mesh tiles)

    ! ISSM Outputs
    real(dp),    pointer, dimension(:)   :: ISSM_OUTPUTS  => null() ! pointer containing all outputs

    ! ice-surface elevation on mesh and landice tiles
    real(dp),    pointer, dimension(:)   :: ICESURF_MESH  => null() ! ice elevation on mesh
    real, pointer, dimension(:)          :: ICESURF_TILE  => null() ! ice elevation on landice tiles
    real, pointer, dimension(:)          :: ICESURF_EX    => null() ! pointer to export state (mesh tiles)
    real(dp), pointer, dimension(:)      :: ICESURF_IN    => null() ! pointer to internal state (mesh tiles)

    ! ice thickness on mesh and landice tiles
    real(dp),    pointer, dimension(:)   :: ICETHICK_MESH => null() ! ice thickness on mesh
    real, pointer, dimension(:)          :: ICETHICK_TILE => null() ! ice thickness on landice tiles
    real, pointer, dimension(:)          :: ICETHICK_EX   => null() ! pointer to ice thickness export state (mesh tiles)
    real(dp), pointer, dimension(:)      :: ICETHICK_IN   => null() ! pointer to ice thicknesss internal state (mesh tiles)

    ! ice-flow velocity in x direction (in projection coordinates)
    real(dp),    pointer, dimension(:)   :: ICEVX_MESH   => null() ! ice x-velocity on mesh
    real, pointer, dimension(:)          :: ICEVX_EX     => null() ! pointer to export state (mesh tiles)
    real(dp), pointer, dimension(:)      :: ICEVX_IN     => null() ! pointer to internal state (mesh tiles)
	

    ! ice-flow velocity in y direction (in projection coordinates)
    real(dp),    pointer, dimension(:)   :: ICEVY_MESH   => null() ! ice y-velocity on mesh
    real, pointer, dimension(:)          :: ICEVY_EX     => null() ! pointer to export state (mesh tiles)
    real(dp), pointer, dimension(:)      :: ICEVY_IN     => null() ! pointer to export state (mesh tiles)

    ! ice mask level set (tracks glacier terminus)
    real(dp),    pointer, dimension(:)   :: IMLS_MESH    => null() ! ice mask level set
    real(dp), pointer, dimension(:)      :: IMLS_IN      => null() ! pointer to internal state (mesh tiles)

    ! ocean mask level set (tracks grounding line)
    real(dp),    pointer, dimension(:)   :: OMLS_MESH    => null() ! ocean mask level set
    real(dp), pointer, dimension(:)      :: OMLS_IN      => null() ! pointer to internal state (mesh tiles)

    ! ice-flow speed on mesh and landice tiles  
    real(dp),    pointer, dimension(:)   :: ICEVEL_MESH  => null() ! ice flow speed on mesh tiles
    real, pointer, dimension(:)          :: ICEVEL_TILE  => null() ! ice flow speed on landice tiles
  
    ! physical parameters
    real(dp), parameter                  :: rho_ice   = 917.0      ! pure ice density [kg m-3]
    real(dp)                             :: ISSM_DT                ! time step [s] (ISSM_DT set in AGCM.rc)

  ! Get the target components name, mesh and vm
  ! -----------------------------------------------------------
    Iam = "Run"
    call ESMF_GridCompGet(GC,name=COMP_NAME,mesh=mesh,vm=vm,_RC)
    
    Iam = trim(COMP_NAME) // Iam

  ! Get my internal MAPL_Generic state
  !----------------------------------
    call MAPL_GetObjectFromGC(GC, MAPL, STATUS)
    _VERIFY(STATUS)

    call MAPL_Get(MAPL,INTERNAL_ESMF_STATE = INTERNAL,_RC )
    

    ! Start Total timer
  !------------------
    call MAPL_TimerOn(MAPL,"TOTAL")
    call MAPL_TimerOn(MAPL,"RUN" )

    call MAPL_Get(MAPL, RUNALARM = ALARM, _RC )
    

    ! run ISSM at specified time steps, 
    ! if bootstrapping restart and issm has run not by final time step, run anyways
    ! with timestep of zero, which just gets restart values
    if (ESMF_AlarmIsRinging (ALARM, RC=STATUS)) then

      ! *************************************************************************** !
      ! BASIC SETUP
      ! *************************************************************************** !

      ! get timestep for ISSM
      call MAPL_GetResource(MAPL, ISSM_DT, Label=trim(COMP_NAME)//"_DT:",DEFAULT=302400.0, _RC)
      
      ! get number of mesh elements
      call ESMF_MeshGet(mesh,nodeCount=num_nodes)

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

      ! get landice tile dimensions
      call MAPL_LocStreamGet(internal_state%locstream, NT_LOCAL=NT, _RC)
      
      call ESMF_UserCompGetInternalState(GC, 'ISSM_TILES', issm_tile_wrap, status); _VERIFY(STATUS)
      issm_tile_state => issm_tile_wrap%ptr
      
      ! *************************************************************************** !
      ! GET ICESMB IMPORT (surface mass balance)
      ! *************************************************************************** ! 
	    ! NOTE: ICESMB (from landice) has been time-averaged between ISSM runs
	    !       hence the name ICESMB_ISSM
	  
      ! allocate tiles for ICESMB 
      if(.not.associated(ICESMB_TILE)) then
        allocate(ICESMB_TILE(NT), STAT=STATUS)
        _VERIFY(STATUS)
        ICESMB_TILE = MAPL_Undef
      end if
      
      ! copy import values into tile array 
      ICESMB_TILE = issm_tile_state%ICESMB_ISSM

	  ! transform ICESMB from landice tiles to mesh 
      call tile_to_mesh(ICESMB_TILE,ICESMB_MESH,_RC)

      ! save ICESMB on mesh elements 
      call MAPL_GetPointer(EXPORT  , ICESMB_EX , 'ICESMB_ISSM' , _RC)

      if(associated(ICESMB_EX)) ICESMB_EX = ICESMB_MESH(internal_state%owned_idx)

      ! *************************************************************************** !
      !  RUN ISSM WITH SMB INPUT AND ICE-ELEVATION OUTPUT
      ! *************************************************************************** !
      ! convert SMB to units of [m/s] (ice-equivalent) before passing to ISSM
      ICESMB_MESH = ICESMB_MESH/rho_ice

      call ESMF_VMBarrier(vm, _RC)
	  call MAPL_TimerOn(MAPL,"ISSMCore" )

      ! call run method from ISSM library 
      call RunISSM(ISSM_DT, c_loc(ICESMB_MESH), c_loc(ISSM_OUTPUTS))

      call ESMF_VMBarrier(vm, _RC)
	  call MAPL_TimerOff(MAPL,"ISSMCore" )

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
      call MAPL_GetPointer(EXPORT, ICESURF_EX, 'ICESURF', _RC)
      if(associated(ICESURF_EX)) ICESURF_EX = ICESURF_MESH(internal_state%owned_idx)

      call MAPL_GetPointer(EXPORT, ICEVX_EX, 'ICEVX', _RC)
      if(associated(ICEVX_EX)) ICEVX_EX = ICEVX_MESH(internal_state%owned_idx)

      call MAPL_GetPointer(EXPORT, ICEVY_EX, 'ICEVY', _RC)
      if(associated(ICEVY_EX)) ICEVY_EX = ICEVY_MESH(internal_state%owned_idx)

      call MAPL_GetPointer(EXPORT, ICETHICK_EX, 'ICETHICK', _RC)
      if(associated(ICETHICK_EX)) ICETHICK_EX = ICETHICK_MESH(internal_state%owned_idx)

      ! set pointers to tile-mesh internals
      call MAPL_GetPointer(INTERNAL, ICESURF_IN, 'ICESURF', _RC)
      if(associated(ICESURF_IN)) ICESURF_IN = ICESURF_MESH(internal_state%owned_idx)

      call MAPL_GetPointer(INTERNAL, ICETHICK_IN, 'ICETHICK', _RC)
      if(associated(ICETHICK_IN)) ICETHICK_IN = ICETHICK_MESH(internal_state%owned_idx)

      call MAPL_GetPointer(INTERNAL, ICEVX_IN, 'ICEVX', _RC)
      if(associated(ICEVX_IN)) ICEVX_IN = ICEVX_MESH(internal_state%owned_idx)

	    call MAPL_GetPointer(INTERNAL, ICEVY_IN, 'ICEVY', _RC)
      if(associated(ICEVY_IN)) ICEVY_IN = ICEVY_MESH(internal_state%owned_idx)

      call MAPL_GetPointer(INTERNAL, OMLS_IN, 'OMLS', _RC)
      if(associated(OMLS_IN)) OMLS_IN = OMLS_MESH(internal_state%owned_idx)

      call MAPL_GetPointer(INTERNAL, IMLS_IN, 'IMLS', _RC)
      if(associated(IMLS_IN)) IMLS_IN = IMLS_MESH(internal_state%owned_idx)

      ! *************************************************************************** !
      ! REGRID MESH FIELDS ONTO LANDICE TILES AND EXPORT VIA INTERNAL STATE
      ! *************************************************************************** !
      ! transform from mesh to tiles
      call mesh_to_tile(ICESURF_MESH,ICESURF_TILE,_RC)
      issm_tile_state%ICESURF_TILE = ICESURF_TILE

      call mesh_to_tile(ICETHICK_MESH,ICETHICK_TILE,_RC)
      issm_tile_state%ICETHICK_TILE = ICETHICK_TILE

      call mesh_to_tile(ICEVEL_MESH,ICEVEL_TILE,_RC)
      issm_tile_state%ICEVEL_TILE = ICEVEL_TILE

    end if 
    
    ! barrier to ensure regridding completes before any deallocates
    call ESMF_VMBarrier(vm,_RC)
    
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

    call MAPL_TimerOff(MAPL,"RUN"  )
    call MAPL_TimerOff(MAPL,"TOTAL")
  
    _RETURN(_SUCCESS)

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

	  type(ESMF_State)                   :: INTERNAL
    
    ! ErrLog Variables
    character(len=ESMF_MAXSTR)	       :: IAm
    integer			                       :: STATUS
    character(len=ESMF_MAXSTR)         :: COMP_NAME


    type(T_ISSM_TILE_STATE), pointer   :: issm_tile_state
    type(ISSM_TILE_WRAP)               :: issm_tile_wrap
	real, pointer, dimension(:)        :: ISSM_NSTEPS

    ! Get the target components name and set-up traceback handle.
    ! -----------------------------------------------------------
    Iam = "Finalize"
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, _RC )
    
    Iam = trim(comp_name) // Iam

    call MAPL_GetObjectFromGC(GC, MAPL, STATUS)
    _VERIFY(STATUS)

	  ! save number of steps since last ISSM run via internal state checkpoints (restarts)
    call ESMF_UserCompGetInternalState(GC, 'ISSM_TILES', issm_tile_wrap, STATUS); _VERIFY(STATUS)
    issm_tile_state => issm_tile_wrap%ptr

    ! get internal state
    call MAPL_Get(MAPL,INTERNAL_ESMF_STATE = INTERNAL,_RC)

    ! get number of time steps since last ISSM run
    call MAPL_GetPointer(INTERNAL, ISSM_NSTEPS, 'ISSM_NSTEPS',_RC)
    ISSM_NSTEPS(:) = real(issm_tile_state%ISSM_NSTEPS)

    ! Generic Finalize
    ! ------------------
    call MAPL_GenericFinalize( GC, IMPORT, EXPORT, CLOCK, _RC )
    
    ! All Done
    ! ------------------

    _RETURN(_SUCCESS)
  end subroutine Finalize


  subroutine mesh_to_tile(VAR_MESH,VAR_TILE,RC)
    ! regrid from mesh to grid, then transform from grid to landice tiles
	! arguments:
    real(dp),    pointer, dimension(:), intent(inout)   :: VAR_MESH           ! var on mesh nodes
    real, pointer, dimension(:), intent(inout)          :: VAR_TILE           ! var on landice tiles
    integer, optional,       intent(OUT)                :: RC                 ! Error code

    ! local variables:
    real, pointer, dimension(:,:)                       :: VAR_GRID => null() ! var on attached grid
    real(dp),    pointer, dimension(:)                  :: VAR_MESH_OWN       ! var on owned mesh nodes
    type(ESMF_Field)                                    :: srcField
    type(ESMF_Field)                                    :: dstField
    integer                                             :: num_owned_nodes
    integer                                             :: NT
    integer                                             :: STATUS
  
    num_owned_nodes = size(internal_state%owned_idx)

    call MAPL_LocStreamGet(internal_state%locstream, NT_LOCAL=NT, _RC)
    
    allocate(VAR_MESH_OWN(num_owned_nodes))

    VAR_MESH_OWN = VAR_MESH(internal_state%owned_idx)

    ! allocate tiles 
    if (.not.associated(VAR_TILE)) then
      allocate(VAR_TILE(NT))
      VAR_TILE = MAPL_Undef
    end if

    ! create source field: field on mesh nodes
    srcField = ESMF_FieldCreate(mesh=internal_state%mesh,farrayPtr=VAR_MESH_OWN,meshloc=ESMF_MESHLOC_NODE, & 
    datacopyflag=ESMF_DATACOPY_VALUE,_RC)
    
    ! create destination field: field on grid
    dstField = ESMF_FieldCreate(grid=internal_state%grid,typekind=ESMF_TYPEKIND_R4,_RC)
    
    ! regrid field from mesh to grid
    call ESMF_FieldRegrid(srcField, dstField, internal_state%routehandle_m2g, _RC)

    ! get pointer to field on grid
    call ESMF_FieldGet(dstField,farrayPtr=VAR_GRID,_RC)

    ! transform from grid to tiles  
    call MAPL_LocStreamTransform(internal_state%locstream,VAR_TILE,VAR_GRID, _RC)
    
    ! destroy regridding fields so they can be reused
    call ESMF_FieldDestroy(srcField,_RC)
    call ESMF_FieldDestroy(dstField,_RC)

    _RETURN(_SUCCESS)

  end subroutine mesh_to_tile

  subroutine tile_to_mesh(VAR_TILE,VAR_MESH,RC)
    ! transform from landice tile to grid, then regrid onto mesh
	! arguments:
    real, pointer, dimension(:), intent(inout)     :: VAR_TILE           ! var on landice tiles
    real(dp), pointer, dimension(:), intent(inout) :: VAR_MESH           ! var on mesh elements
    integer, optional,       intent(OUT)           :: RC                 ! Error code

    ! local variables:
    real, pointer, dimension(:,:)                  :: VAR_GRID => null() ! var on attached grid
    real(dp), pointer, dimension(:)                :: MESH_PTR           ! pointer for ESMF_FieldGet 
    type(ESMF_Field)                               :: srcField
    type(ESMF_Field)                               :: dstField
    type(ESMF_Array)                               :: meshArray
    integer                                        :: num_owned_nodes
    integer                                        :: num_nodes
    integer                                        :: IM, JM, local_dims(3)   
    integer                                        :: STATUS

    ! get number of ndoes
    call ESMF_MeshGet(internal_state%mesh,nodeCount=num_nodes,numOwnedNodes=num_owned_nodes,_RC)
    
    ! get grid dimensions
    call MAPL_GridGet(internal_state%grid, localCellCountPerDim=local_dims, _RC)
    IM = local_dims(1)
    JM = local_dims(2)

    ! allocate pointer on grid for regridding 
    allocate(VAR_GRID(IM,JM))
  
    ! transform from tile to grid
    ! NOTE: we use the "transpose" option with MAPL_LocStreamTransformG2T 
    ! (rather than MAPL_LocStreamTransformT2G) because the "default" value is zero
    ! (rather than MAPL_UNDEF, which leads to errors when regridding onto mesh)
    call MAPL_LocStreamTransform(internal_state%locstream, VAR_TILE, VAR_GRID, TRANSPOSE=.true., _RC)
    
    ! create source field on grid
    srcField = ESMF_FieldCreate(grid=internal_state%grid,farrayPtr=VAR_GRID, datacopyflag=ESMF_DATACOPY_VALUE,_RC)
    
    ! create destination field on mesh elements
    meshArray=ESMF_ArrayCreate(internal_state%nodalDistgrid,typekind=ESMF_TYPEKIND_R8,haloSeqIndexList=internal_state%halolist,_RC)
    
    ! create field on ISSM mesh
    dstField=ESMF_FieldCreate(internal_state%mesh, array=meshArray, meshLoc=ESMF_MESHLOC_NODE, _RC)
    
    ! regrid from grid to mesh
    call ESMF_FieldRegrid(srcField, dstField, internal_state%routehandle_g2m, _RC)

    ! append halo values to end of "owned" array
    call ESMF_FieldHalo(dstField, routehandle=internal_state%halohandle, _RC)

    ! get pointer to field on mesh
    call ESMF_FieldGet(dstField,farrayPtr=MESH_PTR,_RC)

    ! copy values into VAR_MESH
    VAR_MESH(internal_state%owned_idx) = MESH_PTR(1:num_owned_nodes)          ! owned nodes
    VAR_MESH(internal_state%halo_idx) = MESH_PTR(num_owned_nodes+1:num_nodes) ! halo nodes
    
    ! destroy fields and arrays so they can be reused
    deallocate(VAR_GRID)
    call ESMF_FieldDestroy(srcField,_RC)
    call ESMF_FieldDestroy(dstField,_RC)
    call ESMF_ArrayDestroy(meshArray,_RC)

    _RETURN(_SUCCESS)

  end subroutine tile_to_mesh  

end module GEOS_IssmGridCompMod
