!  $Id$

#include "MAPL_Generic.h"

module GEOS_IssmGridCompMod

!BOP
! !MODULE: GEOS_ISSM --- Calls ISSM (Ice-sheet and Sea-level System Model)

! !DESCRIPTION:
!
!   {\tt GEOS\_ISSM} is a wrapper that calls ISSM (C++) IRF methods
!
! *** NOTE: currently just want to call/"drive" ISSM with no coupling at all
!           to test if the building/linking is working    
!           ~~the current grid is just inherited from the parent (landice)~~
!           ~~currently only running over Greenland~~
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

! declare interfaces to ISSM HERE ....
    ! Define the interface for the ISSM C++ functions
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

!NOTE: previous version has "save" here?!?!


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

    ! NOTE: there's other stuff here in landice that I don't know if will be needed for our purpose

    ! Set the state variable specs.
! -----------------------------

!BOS
! TODO: lots of placeholders commentary for IM/IN/EX states
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

    ! Set generic init and final methods
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
    character(len=100), dimension(:), allocatable, target :: argv
    type(c_ptr), dimension(:), allocatable :: argv_ptr
    integer :: i

    ! real(dp),    pointer, dimension(:)     :: SMBToISSM => null()
    ! real(dp),    pointer, dimension(:)     :: SurfaceToGEOS => null()
    ! real(dp) :: dt

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
    ! ! not sure if this is needed:
    ! call MAPL_TimerOn(MAPL,"TOTAL"     )
    ! call MAPL_TimerOn(MAPL,"INITIALIZE")


    ! Rule 10: A component’s grid must be fully formed before MAPL_GenericInitialize is invoked
    ! ! NOTE: so currently just getting grid from parent (landice) ***
    ! Generic initialize
    ! ------------------
    
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
    argv(1) = "/discover/nobackup/agstubbl/ISSM/GEOS-ISSM/ISSM/bin/issm.exe"//c_null_char
    argv(2) = "TransientSolution"//c_null_char
    argv(3) = "/discover/nobackup/agstubbl/ISSM/projs/IRF-ISSM"//c_null_char
    argv(4) = "GreenlandGEOS"//c_null_char

    ! Convert Fortran strings to C pointers (argv)
    allocate(argv_ptr(argc))

    do i = 1, argc
        ! Ensure that we are only getting the memory address once per string
        argv_ptr(i) = c_loc(argv(i))
    end do
   
    
    ! ! print the VM information if desired:
    ! call ESMF_VMPrint(vm, rc=rc)
   
    ! Call the C++ function for initializing ISSM
    ! gets the number of elements and nodes of the mesh
    call InitializeISSM(argc, argv_ptr,num_elements,num_nodes,comm)

    !!! TO CREATE MESH TO DO THIS:
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

    ! create the ESMF mesh 
    mesh = ESMF_MeshCreate(parametricDim=2, spatialDim=2, nodeIds=nodeIds, nodeCoords=nodeCoords, &
           elementIds=elementIds, elementTypes=elementTypes, elementConn=elementConn, coordSys=ESMF_COORDSYS_CART, rc=rc)

    call ESMF_GridCompSet(GC,mesh=mesh,rc=STATUS)
    VERIFY_(STATUS)

    ! ! NOTE: How do we set this mesh to be the GC's grid? ^ does that work?
    ! ! Rule 10: A component’s grid must be fully formed before MAPL_GenericInitialize is invoked

    ! ****************************************************

    ! Get the grid, configuration
    !----------------------------

    ! vvvvvvvvv BAD IDEA BUT BE MIGHT INSIGHTFUL vvvvvvvvv
    ! dt = 0.05
    !   ! ! allocate SMB forcing (input to ISSM) and surface output (export from ISSM)
    ! allocate(SMBToISSM(num_elements))
    ! allocate(SurfaceToGEOS(num_elements))

    ! ! set smb and surface for test 
    ! SMBToISSM(:) = 0     ! placeholder zeros
    ! SurfaceToGEOS(:) = 0 ! placeholder zeros
    
    ! ! NOTE: do we need the barriers before/after ISSM run?
    ! call ESMF_VMBarrier(vm, rc=status)
    ! VERIFY_(STATUS)

    ! ! ! call the C++ routine for running a single time step
    ! call RunISSM(dt, c_loc(SMBToISSM), c_loc(SurfaceToGEOS))

    ! call ESMF_VMBarrier(vm, rc=status)
    ! VERIFY_(STATUS)

    ! ^^^^^^^^^^ BAD IDEA BUT BE MIGHT INSIGHTFUL ^^^^^^^^^^

    ! call MAPL_GenericInitialize( GC, IMPORT, EXPORT, CLOCK, RC=STATUS )
    ! VERIFY_(STATUS)

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

  ! ErrLog Variables
  character(len=ESMF_MAXSTR)          :: IAm
  integer                             :: STATUS
  character(len=ESMF_MAXSTR)          :: COMP_NAME

  real(dp) :: dt
  real(dp),    pointer, dimension(:)     :: SMBToISSM => null()
  real(dp),    pointer, dimension(:)     :: SurfaceToGEOS => null()

  type(MAPL_MetaComp), pointer            :: MAPL

  type(ESMF_Mesh)                :: mesh
  integer(c_int)                 :: num_elements  
  type(ESMF_VM)                  :: vm  
  integer(c_int)                 :: comm  

  ! Get the target components name and set-up traceback handle.
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

  ! timestep for issm...
  dt = 0.05   ! timestep in years

  ! ! need to access num_elements 
  call ESMF_MeshGet(mesh,elementCount=num_elements)

  ! allocate SMB forcing (input to ISSM) and surface output (export from ISSM)
  allocate(SMBToISSM(num_elements))
  allocate(SurfaceToGEOS(num_elements))

! set smb and surface for test 
  SMBToISSM(:) = 0     ! placeholder zeros
  SurfaceToGEOS(:) = 0 ! placeholder zeros
  
  ! NOTE: do we need the barriers before/after ISSM run?
  call ESMF_VMBarrier(vm, rc=status)
  VERIFY_(STATUS)

  PRINT *, "Elements of SMBtoISSM:"
  PRINT *, SMBtoISSM

  ! ! call the C++ routine for running a single time step
!   call RunISSM(dt, c_loc(SMBToISSM), c_loc(SurfaceToGEOS))

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

    ! Get my internal MAPL_Generic state
!-----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=status)
    VERIFY_(STATUS)

    ! Profilers
!----------

    call MAPL_TimerOn(MAPL,"TOTAL"   )
    call MAPL_TimerOn(MAPL,"FINALIZE")

    ! call ISSM finalize (saves binary output .outbin file)
    call FinalizeISSM()

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

end module GEOS_IssmGridCompMod