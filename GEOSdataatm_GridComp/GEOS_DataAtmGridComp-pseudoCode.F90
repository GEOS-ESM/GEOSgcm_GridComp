!  $Id$

#include "MAPL_Generic.h"

!=============================================================================

module GEOS_DataAtmGridCompMod_pseudoCode
!BOP

! !MODULE: GEOS_DataAtm -- An ``atmospheric" component for ocean-only simulations, 
!                          therefore will always run with a REALISTIC sea ice (e.g., CICE)

! --------- it does following:
! 1. Get imports from GCM (ocean, sea ice, fields) e.g., UW, VW, UI, VI, KPAR
! 2. Read atmospheric (surface) fields e.g., T10M, Q10M, U10M, etc  <-- could be automated
! 3. Initialze CICE(4) thermodynamics (as in GEOS_CICE4ColumnPhysGridComp.F90 : CICE_PREP_THERMO)
! 4. Compute stress - this invloves Cd, Ch, Cq. 
!    Currently uses ncar_ocean_fluxes, should use helfsurface. 
!    Code from GEOS_SurfaceGridComp.F90 could be used (for grid to tile, and other preparation for helfsurface)
! 5. Compute or update total precipitation, radiation, heat fluxes
! 6. Update water temperature (SST), salinity (SSS)
! 7. Repeat above steps 4- 6 over sea ice, using CICE4 (as in GEOS_CICE4ColumnPhysGridComp.F90)
! 8. Update fr of sea ice and snow
! 9. Fill up exports to "drive" the ocean, sea ice models
! ---------

! !USES: 

! 
! use <modules> = ESMF, MAPL_Mod, ...
! 
  implicit none
  private

! !PUBLIC MEMBER FUNCTIONS:

! public SetServices

!=============================================================================

! !DESCRIPTION:

! 
!   ... will be rewritten ...
!

!EOP

! 
! integer <variables> = WATER, ICE, NUM_SUBTILES, ... from Surface, SaltWater, CICE4ColPhy, OpenWater GCs
!
   contains

!BOP

! !IROUTINE: SetServices -- Sets ESMF services for this component

! !INTERFACE:

! subroutine SetServices ( GC, RC )

! !ARGUMENTS:

!   type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
!   integer, optional                  :: RC  ! return code

!  !DESCRIPTION:
! 
!   ... will be rewritten ...
!
!
!EOP

!=============================================================================
!
! ErrLog Variables


!   character(len=ESMF_MAXSTR)              :: IAm
!   integer                                 :: STATUS
!   character(len=ESMF_MAXSTR)              :: COMP_NAME

! Local derived type aliases

!   type (ESMF_Config)                      :: CF
!   type (MAPL_MetaComp),  pointer          :: MAPL

!=============================================================================

! Begin...

! Get my name and set-up traceback handle
! ---------------------------------------

!   Iam = "SetServices"
!   call ESMF_GridCompGet( GC, NAME=COMP_NAME, CONFIG=CF, RC=STATUS )
!   VERIFY_(STATUS)
!   Iam = trim(COMP_NAME) // Iam

! Get my MAPL_Generic state
!--------------------------

!   call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
!   VERIFY_(STATUS)

! Set the Initialize and Run entry points
! ---------------------------------------

!   call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_INITIALIZE, Initialize, RC=STATUS)
!   VERIFY_(STATUS)
!   call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,        Run, RC=STATUS)
!   VERIFY_(STATUS)
!   call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_FINALIZE,  Finalize, RC=STATUS )
!   VERIFY_(STATUS)

! Set the state variable specs.
! -----------------------------

!BOS

!  !IMPORT STATE:

! 
!  import state variables = ?? they come from the OGCM and Surface
!
!  !INTERNAL STATE:

! 
!  internal state variables = ?? these would be mostly same as those for SaltWater, OpenWater and CICE4ColPhy GCs
!

!  !EXPORT STATE:

! 
!  export state variables = ?? these would be mostly same as those for SaltWater, OpenWater and CICE4ColPhy GCs 
!                              PLUS ocean biology and radiation
!

!EOS

!   call MAPL_TimerAdd(GC,    name="INITIALIZE",RC=STATUS)
!   VERIFY_(STATUS)
!   call MAPL_TimerAdd(GC,    name="RUN"       ,RC=STATUS)
!   VERIFY_(STATUS)
!   call MAPL_TimerAdd(GC,    name="FINALIZE"  ,RC=STATUS)
!   VERIFY_(STATUS)

! Set generic init and final methods
! ----------------------------------

!   call MAPL_GenericSetServices    ( GC, RC=STATUS)
!   VERIFY_(STATUS)

!   RETURN_(ESMF_SUCCESS)
  
! end subroutine SetServices

!BOP

! !IROUTINE: INITIALIZE -- Initialize stage for the DataAtm component

! !INTERFACE:

!subroutine INITIALIZE ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

! type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
! type(ESMF_State),    intent(inout) :: IMPORT ! Import state
! type(ESMF_State),    intent(inout) :: EXPORT ! Export state
! type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
! integer, optional,   intent(  out) :: RC     ! Error code:

! !DESCRIPTION: 

!EOP

! ErrLog Variables

! character(len=ESMF_MAXSTR)      :: IAm
! integer                         :: STATUS
! character(len=ESMF_MAXSTR)      :: COMP_NAME

! Locals

! type (MAPL_MetaComp), pointer   :: STATE => null()
! type (MAPL_LocStream)           :: LOCSTREAM
! type (MAPL_LocStream)           :: EXCH

! 
! integer                         :: ...
! real                            :: ...
! character(len=ESMF_MAXSTR)      :: ...
! logical                         :: ...

!  Begin...
!----------

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

!   Iam = "Initialize"
!   call ESMF_GridCompGet( GC, name=COMP_NAME, RC=STATUS )
!   VERIFY_(STATUS)
!   Iam = trim(COMP_NAME) // Iam

! Get my internal MAPL_Generic state
!----------------------------------

!   call MAPL_GetObjectFromGC(GC, STATE, STATUS)
!   VERIFY_(STATUS)

! Start Total timer
!------------------

!   call MAPL_TimerOn (STATE,"TOTAL")
!   call MAPL_TimerOn (STATE,"INITIALIZE"  )

! Change the location stream to just the ocean part
!--------------------------------------------------

!   call MAPL_Get(STATE, EXCHANGEGRID=EXCH,        RC=STATUS ); VERIFY_(STATUS)
!   call MAPL_LocStreamCreate(LOCSTREAM, EXCH, NAME='OCEAN', MASK=(/MAPL_OCEAN/), RC=STATUS ); VERIFY_(STATUS)
!   call MAPL_Set(STATE, LOCSTREAM=LOCSTREAM,   RC=STATUS ); VERIFY_(STATUS)

!
!   Set resource parameters. Follow 
!
!   https://github.com/GEOS-ESM/GEOSgcm_GridComp/blob/05bf4235504b5b654dc2c388d6b7d5383b61f7fc/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSsurface_GridComp/GEOSsaltwater_GridComp/GEOS_CICE4ColumnPhysGridComp.F90#L2020
!
!   to 
!
!   https://github.com/GEOS-ESM/GEOSgcm_GridComp/blob/05bf4235504b5b654dc2c388d6b7d5383b61f7fc/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSsurface_GridComp/GEOSsaltwater_GridComp/GEOS_CICE4ColumnPhysGridComp.F90#L2081
!

!
!  ! do some necessary CICE initialization stuff. Follow
!
!  https://github.com/GEOS-ESM/GEOSgcm_GridComp/blob/05bf4235504b5b654dc2c388d6b7d5383b61f7fc/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSsurface_GridComp/GEOSsaltwater_GridComp/GEOS_CICE4ColumnPhysGridComp.F90#L2083
!

!   call MAPL_TimerOff(STATE,"INITIALIZE"  )
!   call MAPL_TimerOff(STATE,"TOTAL")

!   call MAPL_GenericInitialize ( GC, IMPORT, EXPORT, CLOCK,  RC=STATUS)
!   VERIFY_(STATUS)

!   RETURN_(ESMF_SUCCESS)

! end subroutine INITIALIZE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!BOP

! !IROUTINE: RUN -- Run stage for the DataAtm component

! !INTERFACE:

!subroutine RUN ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

! type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
! type(ESMF_State),    intent(inout) :: IMPORT ! Import state
! type(ESMF_State),    intent(inout) :: EXPORT ! Export state
! type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
! integer, optional,   intent(  out) :: RC     ! Error code:

! !DESCRIPTION:
! 
!   ... will be rewritten ...
!

!EOP

! ErrLog Variables

! character(len=ESMF_MAXSTR)      :: IAm
! integer                         :: STATUS
! character(len=ESMF_MAXSTR)      :: COMP_NAME

!
! 
! integer                         :: ...
! real                            :: ...
! character(len=ESMF_MAXSTR)      :: ...
! logical                         :: ...


! Locals

! type (MAPL_MetaComp), pointer   :: STATE => null()
! type (MAPL_SunOrbit)            :: ORBIT
! type (ESMF_State)               :: INTERNAL
!!type (ESMF_logical)             :: FRIENDLY     -- NO MORE of SUCH TYPE!
! type (ESMF_FIELD)               :: FIELD
! type (ESMF_Time)                :: CurrentTime
! type (ESMF_GRID)                :: GRID
! type (ESMF_Config  )            :: CF
! character(len=ESMF_MAXSTR)      :: DATAFILE

! integer                         :: NT
! integer                         :: IFCST
! logical                         :: FCST
!
!
! Need a separate interface to handle all the tiled variables
!
! type(MAPL_LocStream)            :: LOCSTREAM
! real, pointer, dimension(  :)   :: TILELONS => null()
! real, pointer, dimension(  :)   :: TILELATS => null()

! type(ESMF_VM)                   :: VM
! integer                         :: MYPE
!
! declare pointers to all the "atmospheric" fields that will be read in from file(s)
! e.g.
!
! real, pointer, dimension(  :)   :: rain => null()
!

! real                             :: DT

!
! pointers to import states
!
!  real, pointer, dimension(:)  :: x1 => null()

!
! pointers to export states
!
!  real, pointer, dimension(:)  :: x3 => null()

!
! internal pointers to tile variables
!
!  real, pointer, dimension(:)  :: x2 => null()

!
! pointers to internal states
!
!  real, pointer, dimension(:)  :: x4 => null()

! local real allocatables, parameters, chars, ints, logicals, ...


!  Begin...
!----------

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

!   Iam = "Run"
!   call ESMF_GridCompGet( GC, NAME=COMP_NAME, GRID=GRID, RC=STATUS )
!   VERIFY_(STATUS)
!   Iam = trim(COMP_NAME) // Iam

! Get my MAPL_Generic (GG) state
!-------------------------------

!   call MAPL_GetObjectFromGC(GC, STATE, STATUS)
!   VERIFY_(STATUS)

!   call MAPL_Get(STATE, LOCSTREAM=LOCSTREAM,   RC=STATUS )
!   VERIFY_(STATUS)

! Start timers
!-------------

!   call MAPL_TimerOn(STATE,"TOTAL")
!   call MAPL_TimerOn(STATE,"RUN" )

! Get info from the GG state
!---------------------------

!   call MAPL_Get(STATE,                        &
!       INTERNAL_ESMF_STATE = INTERNAL,         &
!       TILELONS=TILELONS,                      &
!       TILELATS=TILELATS,                      &
!       ORBIT   = ORBIT,                        &
!                                     RC=STATUS )
!   VERIFY_(STATUS)

! The number of tiles we are working on
!--------------------------------------

!   NT = size(TILELONS)

! Temporary space for reading forcings
!-------------------------------------
!   allocate( rain(NT),STAT=STATUS);  VERIFY_(STATUS)
!

! Get the time step
! -----------------

!   call MAPL_GetResource ( STATE, DT, Label="RUN_DT:"        , RC=STATUS)
!   VERIFY_(STATUS)
!   call MAPL_GetResource ( STATE, DT, Label="DT:", DEFAULT=DT, RC=STATUS)
!   VERIFY_(STATUS)

! Get current time from clock
!----------------------------
!   call ESMF_ClockGet(CLOCK, currTime=CurrentTime, TIMESTEP=DELT, RC=STATUS)
!   VERIFY_(STATUS)

!
! ANY work (getting pointers to import, exports, internals, etc) 
! repetitive should be done in a subroutine. Or else, entire GC is too long and bloated.
!
! Do extra allocation if gridded exports are requested
!-----------------------------------------------------
!   call MK_GRID_OUT(EXPORT, GNAME='SWNg', TNAME='SWN', RC=STATUS)

! Pointers to Imports
!--------------------
!   call GET_POINTER(IMPORT,   KPAR  , 'KPAR'    ,    RC=STATUS); VERIFY_(STATUS)

! Pointers to Internals
!----------------------
!   call GET_POINTER(INTERNAL, TW ,    'TSKINW',    RC=STATUS); VERIFY_(STATUS)

!  Pointers to Exports
!---------------------
!   call GET_POINTER(EXPORT, u10m, 'U10M'   ,    RC=STATUS); VERIFY_(STATUS)

! Read atmospheric (surface) fields. For e.g.,
! XX= t10m, q10m, u10m, v10m, etc; includes aerosols and others for ocean bio application(s)
!-------------------------------------------------------------------------------------------

!   such calls as following >>>>
!
!   call MAPL_GetResource(STATE, DataFile, LABEL = "XX_FILE:", DEFAULT = "none", RC = status)
!   VERIFY_(status)
!   if(trim(DataFile) == 'none') then
!     XX = 290.0
!   else
!    call MAPL_ReadForcing(state, "XX", RenameFile(DataFile, TIME = currenttime), currenttime, XX, rc = status);
!    VERIFY_(status)
!   endif 
!
!   <<<< should be wrapped into a 
!        call read_atm_vars( t10m, q10m, u10m, v10m, ...)
!   it will be "read_atm_vars" responsibility to provide reasonable values (read file names, default values, etc).
!
!   @tclune https://github.com/GEOS-ESM/GEOSgcm_GridComp/issues/140 will be applied there?
!

!   Initialize CICE Thermodynamics, simply by calling CICE_PREP_THERMO(...)
!
!   https://github.com/GEOS-ESM/GEOSgcm_GridComp/blob/05bf4235504b5b654dc2c388d6b7d5383b61f7fc/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSsurface_GridComp/GEOSsaltwater_GridComp/GEOS_CICE4ColumnPhysGridComp.F90#L3698
!
!   thereafter
!   
!   https://github.com/GEOS-ESM/GEOSgcm_GridComp/blob/05bf4235504b5b654dc2c388d6b7d5383b61f7fc/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSsurface_GridComp/GEOSsaltwater_GridComp/GEOS_CICE4ColumnPhysGridComp.F90#L4032
!

! Total precipitation = RAIN + SNOW * (1.-Fr sea ice)
! Shortwave radiation (PRUVF, PRPAF, PRUVR, PRPAR)
! Air density (RHOa) using T10M, Q10M and SLP
! QSAT using RHOa, "water temperature"
!
! Calculate exchange coefficients:
!
! WATER:
!  [cd, ch, ce, ustar, bstar] = ncar_ocean_fluxes(...)

! compute stress over water:   tauxw, tauyw 
! compute stress over sea ice: tauxi, tauyi = atmo_boundary_layer(...)
! total stress over ocean:     tauxo, tauyo = tau{x,y}W*fr W + tau{x,y}I* fr I          

! Update mass of salt and fresh water (kg/m2)
! Relax to SSS (using TAUsss) 
! 

! Update temperature over WATER, after computing
!      1. Latent heat flux   (positive down)
!      2. Absobed SW rad     (positive down)
!      3. Net     LW rad     (positive down) = lwrad - MAPL_STFBOL*(TW)**4
!      4. Latent heat of snow melt
!      5. Sensible heat flux (positive down)

! Update surface salinity over water and sea ice

! Update albedos (water and sea ice)
!
!   water follows: https://github.com/GEOS-ESM/GEOSgcm_GridComp/blob/05bf4235504b5b654dc2c388d6b7d5383b61f7fc/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSsurface_GridComp/GEOSsaltwater_GridComp/GEOS_OpenWaterGridComp.F90#L3034
!
!   sea ice follows: https://github.com/GEOS-ESM/GEOSgcm_GridComp/blob/05bf4235504b5b654dc2c388d6b7d5383b61f7fc/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSsurface_GridComp/GEOSsaltwater_GridComp/GEOS_CICE4ColumnPhysGridComp.F90#L5605
!

! Update temperature over SEA ICE, after computing
!      1. Latent heat flux   (positive down)
!      2. Absobed SW rad     (positive down)
!      3. Net     LW rad     (positive down) = lwrad - MAPL_STFBOL*(TW)**4
!      4. Latent heat of snow melt
!      5. Sensible heat flux (positive down)
!
!   much of this update follows calls to:
!       CICE_THERMO1(...) as in https://github.com/GEOS-ESM/GEOSgcm_GridComp/blob/05bf4235504b5b654dc2c388d6b7d5383b61f7fc/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSsurface_GridComp/GEOSsaltwater_GridComp/GEOS_CICE4ColumnPhysGridComp.F90#L4032
!   
!  and 
!      
!       CICE_THEMO2_STEP1 https://github.com/GEOS-ESM/GEOSgcm_GridComp/blob/05bf4235504b5b654dc2c388d6b7d5383b61f7fc/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSsurface_GridComp/GEOSsaltwater_GridComp/GEOS_CICE4ColumnPhysGridComp.F90#L4317
!
! and   its step2
!    
! Update Fr of sea ice, snow         
! --------------------------

! Fill pointers to gridded exports 
! -----------------------------------
!   call T2G_Regrid(EXPORT, LOCSTREAM, GNAME='SWNg', TNAME='SWN', RC=STATUS)
!


! deallocate arrays

!  All done
!-----------

!   call MAPL_TimerOff(STATE,"RUN"  )
!   call MAPL_TimerOff(STATE,"TOTAL")

!   RETURN_(ESMF_SUCCESS)

!end subroutine RUN

!subroutine T2G_Regrid(STATE, LOCSTREAM, GNAME, TNAME, RC)
!  type(ESMF_State)     :: STATE
!  type(MAPL_LocStream) :: LOCSTREAM
!  character(len=*)     :: GNAME
!  character(len=*)     :: TNAME
!  integer, optional    :: RC
!  
!  real, pointer :: GVAR(:,:) => null()
!  real, pointer :: TVAR(:)   => null()
!  
!  character(len=ESMF_MAXSTR)   :: IAm='T2G_Regrid'
!  integer                      :: STATUS
! 
!  call GET_POINTER(STATE, GVAR, GNAME, RC=STATUS)
!  VERIFY_(STATUS)
!  if(associated(GVAR)) then
!     call GET_POINTER(STATE, TVAR, TNAME, RC=STATUS)
!     VERIFY_(STATUS)
!     call MAPL_LocStreamTransform( LOCSTREAM, GVAR, TVAR, RC=STATUS) 
!     VERIFY_(STATUS)
!  endif
!  RETURN_(ESMF_SUCCESS)
!end subroutine T2G_Regrid

!subroutine MK_GRID_OUT(STATE, GNAME, TNAME, RC)
!  type(ESMF_State) :: STATE
!  character(len=*) :: GNAME
!  character(len=*) :: TNAME
!  integer, optional, intent(OUT) :: RC
  
!  real, pointer                  :: GVAR(:,:) => null()
!!  real, pointer                  :: TVAR(:) => null()
!  character(len=ESMF_MAXSTR)     :: IAm='MK_GRID_OUT'
!  integer                        :: STATUS
!  
!  call GET_POINTER(STATE, GVAR, GNAME, RC=STATUS)
!  VERIFY_(STATUS)
!
!  if(associated(GVAR)) then
!     call GET_POINTER(STATE, TVAR, TNAME, alloc=.true., RC=STATUS)
!     VERIFY_(STATUS)
!     TVAR = MAPL_Undef
!  end if
!  
!  RETURN_(ESMF_SUCCESS)
!end subroutine MK_GRID_OUT

!----------------------------------------------------------------------------------------------------------------------------------

!   function renamefile(name, time) result(name0)

!     CHARACTER(len = *),       intent(in) :: name
!     CHARACTER(len = len(name))           :: name0 
!     TYPE(esmf_time), intent(inout)       :: time

!     integer :: year, month, day, status, i
!         
!       name0 = trim(name)
!       i = index(string = name, substring = "yyyymmdd")
!       if(i == 0) return

!       call ESMF_TIMEGET(TIME, YY = year, MM = month, DD = day, RC = status) 
!       WRITE(unit = name0(i:i + 7), fmt = "(i4,i2.2,i2.2)") year, month, day

!   end function

!----------------------------------------------------------------------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! !IROUTINE: Finalize

! !INTERFACE:

! subroutine Finalize ( GC, IMPORT, EXPORT, CLOCK, RC ) 

! !ARGUMENTS:

! type(ESMF_GridComp), intent(INOUT) :: GC     ! Gridded component 
! type(ESMF_State),    intent(INOUT) :: IMPORT ! Import state
! type(ESMF_State),    intent(INOUT) :: EXPORT ! Export state
! type(ESMF_Clock),    intent(INOUT) :: CLOCK  ! The supervisor clock
! integer, optional,   intent(  OUT) :: RC     ! Error code:

!EOP

!   type (MAPL_MetaComp), pointer:: MAPL 

! ErrLog Variables

!   character(len=ESMF_MAXSTR)       :: IAm
!   integer                          :: STATUS
!   character(len=ESMF_MAXSTR)       :: COMP_NAME

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

!   Iam = "Finalize"
!   call ESMF_GridCompGet( gc, NAME=comp_name, RC=status )
!   VERIFY_(STATUS)
!   Iam = trim(comp_name) // Iam

! Get my internal MAPL_Generic state
!-----------------------------------

!   call MAPL_GetObjectFromGC ( GC, MAPL, RC=status)
!   VERIFY_(STATUS)

! Profilers
!----------

!   call MAPL_TimerOn(MAPL,"TOTAL"   )
!   call MAPL_TimerOn(MAPL,"FINALIZE")


!   call MAPL_TimerOff(MAPL,"FINALIZE")
!   call MAPL_TimerOff(MAPL,"TOTAL"   )

! Generic Finalize
! ------------------
    
!   call MAPL_GenericFinalize( GC, IMPORT, EXPORT, CLOCK, RC=status )
!   VERIFY_(STATUS)


! All Done
!---------

!   RETURN_(ESMF_SUCCESS)
! end subroutine Finalize

end module GEOS_DataAtmGridCompMod_pseudoCode
