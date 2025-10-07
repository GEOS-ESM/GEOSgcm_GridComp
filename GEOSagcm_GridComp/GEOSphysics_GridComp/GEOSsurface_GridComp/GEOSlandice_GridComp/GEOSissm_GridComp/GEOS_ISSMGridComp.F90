!  $Id$

#include "MAPL_Generic.h"

module GEOS_IssmGridCompMod

!BOP
! !MODULE: GEOS_ISSM --- Calls ISSM (Ice-sheet and Sea-level System Model)

! !DESCRIPTION:
!
!   {\tt GEOS\_ISSM} is a wrapper that calls ISSM (C++) IRF methods
!
! *** NOTES: *currently just want to call/"drive" ISSM with no coupling at all
!             to test if the building/linking is working:   
!           
!            *the current grid is just inherited from the parent (landice)
!            *currently just want to run over Greenland with all given PETs
!            *ISSM mesh is *internal* to ISSM but we've created an ESMF_MESH version for future regridding        
!
! next steps:  
!             (1) create GC from ESMF_Mesh corresponding to ISSM's mesh
!                 this should greatly simplify imports/exports. 
!                 then we should be able to call ESMF_GridCompCreate in Initialize
!                 *before* MAPL_GenericInitialize is called
!    
!             (2) make sure to only run ISSM at appropriate timesteps,...
!
!             (3) add import/export states corresponding to surface mass balance and 
!                 ice-surface elevation (regrid from parent grid to mesh?)
! 
!             (4) determine realistic ISSM setup, etc...
! 

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
   type(c_ptr),      value   :: nodeCoords
end subroutine GetNodesISSM

subroutine GetElementsISSM(elementIds, elementConn) bind(C,NAME="GetElementsISSM")
   import :: c_ptr
   type(c_ptr),      value   :: elementIds
   type(c_ptr),      value   :: elementConn
end subroutine GetElementsISSM

subroutine FinalizeISSM() bind(C,NAME="FinalizeISSM")
end subroutine FinalizeISSM

end interface


private


! declare any variables here?

public SetServices

! !DESCRIPTION:
! 
!   {\tt GEOS\_Landice} is a light-weight gridded component that updates
!      the landice tiles
!

!EOP

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

! TODO: placeholders comments for IM/IN/EX states
!  !Export state:
    ! TODO: set export states here, like this: 
!    call MAPL_AddExportSpec(GC,                             &
!    SHORT_NAME         = 'EMIS',                              &
!    LONG_NAME          = 'surface_emissivity',                &
!    UNITS              = '1',                                 &
!    DIMS               = MAPL_DimsTileOnly,                   &
!    VLOCATION          = MAPL_VLocationNone,                  &
!                                                   RC=STATUS  )
! VERIFY_(STATUS)

    !  !Internal state:
    !   ! TODO: add internal states here, like this:

!    call MAPL_AddInternalSpec(GC,                           &
!    SHORT_NAME         = 'TS',                                &
!    LONG_NAME          = 'surface_skin_temperature',          &
!    UNITS              = 'K',                                 &
!    DIMS               = MAPL_DimsTileOnly,                   &
!    UNGRIDDED_DIMS     = (/NUM_SUBTILES/),                    &
!    VLOCATION          = MAPL_VLocationNone,                  &
!    DEFAULT            = 280.0,                               &
!                                                   RC=STATUS  )
! VERIFY_(STATUS)

!  !Import state:
!  ! TODO: add import states here, like this:
!    call MAPL_AddImportSpec(GC,                             &
!    SHORT_NAME         = 'ALW',                               &
!    LONG_NAME          = 'linearization_of_surface_emitted_longwave_flux', &
!    UNITS              = 'W m-2',                             &
!    DIMS               = MAPL_DimsTileOnly,                   &
!    VLOCATION          = MAPL_VLocationNone,                  &
!                                                   RC=STATUS  )
! VERIFY_(STATUS)

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
    type(ESMF_GridComp),     intent(INOUT) :: GC     ! Gridded component 
    type(ESMF_State),        intent(INOUT) :: IMPORT ! Import state
    type(ESMF_State),        intent(INOUT) :: EXPORT ! Export state
    type(ESMF_Clock),        intent(INOUT) :: CLOCK  ! The clock
    integer, optional,       intent(  OUT) :: RC     ! Error code:
    type(MAPL_MetaComp), pointer            :: MAPL

    ! Locals with ESMF and MAPL types
    type(ESMF_VM)                  :: vm    
    type(ESMF_Mesh)                :: mesh
    integer(c_int)                 :: comm
    integer, pointer, dimension(:) :: elementIds
    integer, pointer, dimension(:) :: elementConn    
    integer, allocatable  :: elementTypes(:)

    ! ISSM-related variables
    integer(c_int)                 :: num_elements
    integer(c_int)                 :: num_nodes   
    integer(c_int)                 :: argc
    character(len=200), dimension(:), allocatable, target :: argv
    type(c_ptr), dimension(:), allocatable :: argv_ptr
    integer :: i

    real(dp),    pointer, dimension(:)     :: nodeCoords => null()
    integer,     pointer, dimension(:)     :: nodeIds => null()

    ! ErrLog Variables
    character(len=ESMF_MAXSTR)		   :: IAm
    integer				   :: STATUS
    character(len=ESMF_MAXSTR)             :: COMP_NAME

    ! Begin... 

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

    ! Profilers
    !----------

    ! try generic initialize before doing ISSM stuff ?!?!?
    call MAPL_GenericInitialize( GC, IMPORT, EXPORT, CLOCK, RC=STATUS )
    VERIFY_(STATUS)

    call ESMF_VMGetCurrent(vm, rc=STATUS)
    VERIFY_(STATUS)

    call ESMF_VMGet(vm,mpiCommunicator=comm,rc=STATUS)
    VERIFY_(STATUS)

    ! ****************************************************
    ! call ISSM initial C++ code so we can set up mesh

    ! Manually set command line argc and argv to initialize ISSM 
    argc = 4  
    allocate(argv(argc))
    argv(1) = "/discover/nobackup/projects/gmao/SIteam/ISSM/2025-09-02/ifort_2021.13.0-intelmpi_2021.13.0/ISSM/bin/issm.exe"//c_null_char
    argv(2) = "TransientSolution"//c_null_char
    argv(3) = "/discover/nobackup/agstubbl/ISSM/projs/IRF-ISSM"//c_null_char
    argv(4) = "GreenlandGEOS"//c_null_char

    ! Convert Fortran strings to C pointers (argv)
    allocate(argv_ptr(argc))

    do i = 1, argc
        ! Ensure that we are only getting the memory address once per string
        argv_ptr(i) = c_loc(argv(i))
    end do

    !VERIFY_(STATUS)
    call ESMF_VMBarrier(vm, rc=status)
    
    ! Call the C++ function for initializing ISSM
    ! gets the number of elements and nodes of the mesh
    call InitializeISSM(argc, argv_ptr,num_elements,num_nodes,comm)


    !TO CREATE MESH TO DO THIS:
    !allocate mesh-related pointers
    allocate(nodeCoords(2*num_nodes))
    allocate(nodeIds(num_nodes))
    allocate(elementTypes(num_elements))
    allocate(elementIds(num_elements))
    allocate(elementConn(3*num_elements))

    ! create ESMF mesh corresponding to  ISSM mesh 
    ! get information about nodes and elements
    call GetNodesISSM(c_loc(nodeIds), c_loc(nodeCoords))
    call GetElementsISSM(c_loc(elementIds), c_loc(elementConn))

    elementTypes(:) = ESMF_MESHELEMTYPE_TRI

    ! create the ESMF mesh from ISSM mesh properties
    mesh = ESMF_MeshCreate(parametricDim=2, spatialDim=2, nodeIds=nodeIds, nodeCoords=nodeCoords, &
           elementIds=elementIds, elementTypes=elementTypes, elementConn=elementConn, coordSys=ESMF_COORDSYS_SPH_DEG, rc=rc)

    ! associate ESMF_Mesh representation of ISSM mesh with GC for regridding imports/exports (future work / to-do)       
    call ESMF_GridCompSet(GC,mesh=mesh,rc=STATUS)
    VERIFY_(STATUS)

    ! deallocate pointers
    deallocate(argv)
    deallocate(argv_ptr)
    deallocate(nodeCoords)
    deallocate(nodeIds)
    deallocate(elementTypes)
    deallocate(elementIds)
    deallocate(elementConn)

    RETURN_(ESMF_SUCCESS)
  end subroutine Initialize

  !BOP

! !IROUTINE: RUN -- Run stage for the DataSeaIce component

! !INTERFACE:

subroutine RUN ( GC, IMPORT, EXPORT, CLOCK, RC )
! ! Run ISSM ice-sheet model (only over Greenland right now)
! !ARGUMENTS:
  type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
  type(ESMF_State),    intent(inout) :: IMPORT ! Import state
  type(ESMF_State),    intent(inout) :: EXPORT ! Export state
  type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
  integer, optional,   intent(  out) :: RC     ! Error code:

  type (ESMF_Alarm)                  :: ALARM
  
  ! ErrLog Variables
  character(len=ESMF_MAXSTR)          :: IAm
  integer                             :: STATUS
  character(len=ESMF_MAXSTR)          :: COMP_NAME

  real(dp) :: dt
  real(dp) :: dt_yr
  real(dp) :: sec_per_year
  real(dp),    pointer, dimension(:)     :: SMBToISSM => null()
  real(dp),    pointer, dimension(:)     :: SurfaceToGEOS => null()

  type(MAPL_MetaComp), pointer            :: MAPL

  type(ESMF_Mesh)                :: mesh
  integer(c_int)                 :: num_elements  
  type(ESMF_VM)                  :: vm  

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------
! also get the mesh and vm
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



  ! get timestep for ISSM
  call MAPL_GetResource ( MAPL, dt, Label=trim(COMP_NAME)//"_DT:", RC=STATUS)

  ! convert dt to years for ISSM input
  sec_per_year = 31557600.0
  dt_yr = dt/sec_per_year

  ! get number of mesh elements
  call ESMF_MeshGet(mesh,elementCount=num_elements)

  ! allocate SMB forcing (input to ISSM) and surface output (export from ISSM)
  allocate(SMBToISSM(num_elements))
  allocate(SurfaceToGEOS(num_elements))

  call ESMF_VMBarrier(vm, rc=status)
  VERIFY_(STATUS)

! set smb and surface to zero for simple test 
  SMBToISSM(:) = 0.0_dp     ! placeholder zeros
  SurfaceToGEOS(:) = 0.0_dp ! placeholder zeros
  
  ! NOTE: do we need the barriers before/after ISSM run?
  call ESMF_VMBarrier(vm, rc=status)
  VERIFY_(STATUS)

  ! run ISSM at specified time steps
  if ( ESMF_AlarmIsRinging (ALARM, RC=STATUS) ) then
    call RunISSM(dt_yr, c_loc(SMBToISSM), c_loc(SurfaceToGEOS))
    VERIFY_(STATUS)
  end if 

  call ESMF_VMBarrier(vm, rc=status)
  VERIFY_(STATUS)

  deallocate(SMBToISSM)
  deallocate(SurfaceToGEOS)

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