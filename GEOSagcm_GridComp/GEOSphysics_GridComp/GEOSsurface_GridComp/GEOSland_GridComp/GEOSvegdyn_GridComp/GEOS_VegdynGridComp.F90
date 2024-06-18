
!  $Id$

#include "MAPL_Generic.h"


!=============================================================================
module GEOS_VegdynGridCompMod

!BOP

! !MODULE: GEOS_Vegdyn -- child to the "Land" gridded component.  

!DESCRIPTION:
!   {\tt GEOS\_Vegdyn} is a gridded component that performs the
!   necessary interpolation to provide refreshed values of the 
!   dynamic vegetation values prescribed by external data/observations.\\
!
! There are no imports to this routine.
! Exports from this routine are the instaneous values of the
! vegetation parameters on tilespace to be used in other components
! of the land subroutine.  All exports and imports are stored on the
! tile grid inherited from the parent routine.\\
! 
! I. Parameter Class 1: Time AND spatially dependent parameters 
! from a binary data file\\
! 
! Current list: LAI, GRN, NDVI \\
! 
! The gridded component stores the surrounding observations of 
! each parameter in the internal state.  If the run method 
! discovers that the current internal state does not contain the 
! observed values required to interpolate the values at the current 
! time, it performs the required i/o to refresh the values of 
! the internal state.  The first iteration of the run method 
! always has to fill the values.  No restart is required by this 
! gridded component for these parameters.  (A restart *is* now
! required for Vegetation Class 3 \\
!
! INTERNALS: ITY, Z2CH, ASCATZ0\\
!
! EXPORTS:  LAI, GRN, NDVI\\
!
! !USES: 

  use ESMF
  use MAPL
  implicit none
  private

! !PUBLIC MEMBER FUNCTIONS:

  public SetServices

!EOP

  integer, parameter		     :: NTYPS = MAPL_NumVegTypes
  real,    dimension(   NTYPS)       :: VGRT
  ! real,    dimension(   NTYPS)       :: VGZ2   

  data VGRT  / 19700., 7000., 9400., 7000., 7000., 14000./
  ! commented out legacy look-up table for veg heights, which are now always from bcs via restarts, - reichle, 17 March 2020
  ! data VGZ2 / 35.0, 20.0, 17.0, 0.6, 0.5, 0.6/ ! Dorman and Sellers (1989)   
  real, pointer :: LAIens0(:),GRNens0(:), NDVIens0(:), ROOTLens0(:)  
contains

!BOP

! !IROUTINE: SetServices -- Sets ESMF services for this component

! !INTERFACE:

  subroutine SetServices ( GC, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    integer, optional                  :: RC  ! return code

! !DESCRIPTION: This version uses the MAPL\_GenericSetServices. This function sets
!                the Initialize and Finalize services, as well as allocating
!   our instance of a generic state and putting it in the 
!   gridded component (GC). Here we only need to set the run method and
!   add the state variable specifications (also generic) to our instance
!   of the generic state. This is the way our true state variables get into
!   the ESMF\_State INTERNAL, which is in the MAPL\_MetaComp.

!EOP

!=============================================================================
!
! ErrLog Variables


    character(len=ESMF_MAXSTR)              :: IAm
    integer                                 :: STATUS
    character(len=ESMF_MAXSTR)              :: COMP_NAME

! Local derived type aliases

    type(MAPL_MetaComp),pointer             :: MAPL=>null()

!=============================================================================

! Begin...

!------------------------------------------------------------
! Get my name and set-up traceback handle
!------------------------------------------------------------

    call ESMF_GridCompGet(GC                             ,&
                          NAME=COMP_NAME                 ,&
                          RC=STATUS )

    VERIFY_(STATUS)

    Iam = trim(COMP_NAME) // 'SetServices'

! -----------------------------------------------------------
! Set the Run entry point
! -----------------------------------------------------------

    call MAPL_GridCompSetEntryPoint (GC, ESMF_METHOD_RUN, Run, RC=STATUS)
    VERIFY_(STATUS)

! -----------------------------------------------------------
! Get the configuration
! -----------------------------------------------------------
    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)

! -----------------------------------------------------------
! Get the intervals
! -----------------------------------------------------------

    !call MAPL_GetResource ( MAPL,DT, Label="RUN_DT:", RC=STATUS)
    !VERIFY_(STATUS)

    !RUN_DT = nint(DT)
    
! -----------------------------------------------------------
! At the moment, this will refresh when the land parent 
! needs to refresh.
!
!    call ESMF_ConfigGetFloat ( CF, DT, Label=trim(COMP_NAME)//&
!    "_DT:", default=DT, RC=STATUS)
!     VERIFY_(STATUS)
!
!    MY_STEP = nint(DT)
!
! -----------------------------------------------------------


! -----------------------------------------------------------
! Set the state variable specs.
! -----------------------------------------------------------

!BOS

! -----------------------------------------------------------
!   Import States
! None at the moment
! -----------------------------------------------------------

! -----------------------------------------------------------
! Internal State 
! -----------------------------------------------------------

    call MAPL_AddInternalSpec(GC                                ,&
         SHORT_NAME = 'ITY'                                     ,&
         LONG_NAME  = 'vegetation_type'			        ,&
         UNITS      = '1'                                       ,&
         DIMS       = MAPL_DimsTileOnly                         ,&
         VLOCATION  = MAPL_VLocationNone                        ,&
         FRIENDLYTO = trim(COMP_NAME)                           ,&
         RC=STATUS  )
    VERIFY_(STATUS)  

    call MAPL_AddInternalSpec(GC                                ,&
         SHORT_NAME = 'Z2CH'                                    ,&
         LONG_NAME  = 'vegetation_height'			,&
         UNITS      = 'm'                                       ,&
         DIMS       = MAPL_DimsTileOnly                         ,&
         VLOCATION  = MAPL_VLocationNone                        ,&
         FRIENDLYTO = trim(COMP_NAME)                           ,&
         RC=STATUS  )
    VERIFY_(STATUS)  

    call MAPL_AddInternalSpec(GC                                ,&
         SHORT_NAME = 'ASCATZ0'                                 ,&
         LONG_NAME  = 'ASCAT_roughness_length'	                ,&
         UNITS      = 'm'                                       ,&
         DIMS       = MAPL_DimsTileOnly                         ,&
         VLOCATION  = MAPL_VLocationNone                        ,&
         FRIENDLYTO = trim(COMP_NAME)                           ,&
         RC=STATUS  )
    VERIFY_(STATUS)  

! -----------------------------------------------------------
! These are variables that are considered time-independent
! and are stored and retrieved as-is
! -----------------------------------------------------------

! -----------------------------------------------------------
! Export Variables
! -----------------------------------------------------------

    call MAPL_AddExportSpec(GC                                ,&
       SHORT_NAME = 'LAI'                                     ,&
       LONG_NAME  = 'leaf_area_index'                         ,&
       UNITS      = '1'                                       ,&
       DIMS       = MAPL_DimsTileOnly                         ,&
       VLOCATION  = MAPL_VLocationNone                        ,&
       RC=STATUS  )

    VERIFY_(STATUS)  

    call MAPL_AddExportSpec(GC                                ,&
       SHORT_NAME = 'GRN'                                     ,&
       LONG_NAME  = 'greeness_fraction'			      ,&
       UNITS      = '1'                                       ,&
       DIMS       = MAPL_DimsTileOnly                         ,&
       VLOCATION  = MAPL_VLocationNone                        ,&
       RC=STATUS  )

    VERIFY_(STATUS)  	 

    call MAPL_AddExportSpec(GC                                ,&
       SHORT_NAME = 'ROOTL'                                   ,&
       LONG_NAME  = 'root_length_density'                     ,&
       UNITS      = 'm+2'                                     ,&
       DIMS       = MAPL_DimsTileOnly                         ,&
       VLOCATION  = MAPL_VLocationNone                        ,&
       RC=STATUS  )

    VERIFY_(STATUS)  	 

    call MAPL_AddExportSpec(GC                                ,&
       SHORT_NAME = 'NDVI'                                    ,&
       LONG_NAME  = 'normalized_difference_vegetation_index'  ,&
       UNITS      = '1'                                       ,&
       DIMS       = MAPL_DimsTileOnly                         ,&
       VLOCATION  = MAPL_VLocationNone                        ,&
       RC=STATUS  )

    VERIFY_(STATUS)  

!EOS

!------------------------------------------------------------
! Set generic init and final methods
!------------------------------------------------------------

    call MAPL_GenericSetServices(GC, RC=STATUS)
    VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)
  
  end subroutine SetServices

! -----------------------------------------------------------
! RUN -- Run method for the vegdyn component
! -----------------------------------------------------------

  subroutine RUN (GC,IMPORT, EXPORT, CLOCK, RC )

! -----------------------------------------------------------
! !ARGUMENTS:

    type(ESMF_GridComp), intent(inout) :: GC    
    type(ESMF_State),    intent(inout) :: IMPORT
    type(ESMF_State),    intent(inout) :: EXPORT
    type(ESMF_Clock),    intent(inout) :: CLOCK
    integer, optional,   intent(  out) :: RC

!EOP

! ErrLog Variables

    character(len=ESMF_MAXSTR)          :: IAm
    integer                             :: STATUS
    character(len=ESMF_MAXSTR)          :: COMP_NAME

! Locals

    type (MAPL_MetaComp),     pointer   :: MAPL=>null()
    type (ESMF_State       )            :: INTERNAL

! INTERNAL pointers
 
    real, dimension(:), pointer :: ITY
    real, dimension(:), pointer :: Z2CH
    real, dimension(:), pointer :: ASCATZ0

! EXPORT pointers 

    real, dimension(:), pointer :: LAI
    real, dimension(:), pointer :: GRN
    real, dimension(:), pointer :: ROOTL
    real, dimension(:), pointer :: NDVI
  
! Time attributes and placeholders

    type(ESMF_Time) :: CURRENT_TIME

! Others

    character(len=ESMF_MAXSTR)         :: LAIFile
    character(len=ESMF_MAXSTR)         :: GRNFile
    character(len=ESMF_MAXSTR)         :: NDVIFile
    character(len=ESMF_MAXSTR)         :: LAItpl
    character(len=ESMF_MAXSTR)         :: GRNtpl
    character(len=ESMF_MAXSTR)         :: NDVItpl
    integer                            :: NUM_LDAS_ENSEMBLE, ens_id_width, ldas_ens_id, ldas_first_ens_id

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    call ESMF_GridCompGet(GC, name=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
  
    Iam = trim(COMP_NAME) // "Run"

! Get my internal MAPL_Generic state
! -----------------------------------------------------------

    call MAPL_GetObjectFromGC(GC, MAPL, STATUS)
    VERIFY_(STATUS)

    call MAPL_TimerOn(MAPL,"TOTAL")

    call MAPL_Get(MAPL, INTERNAL_ESMF_STATE=INTERNAL, RC=STATUS)
    VERIFY_(STATUS) 

    call MAPL_GetResource ( MAPL, NUM_LDAS_ENSEMBLE, Label="NUM_LDAS_ENSEMBLE:", DEFAULT=1, RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource ( MAPL, ens_id_width, Label="ENS_ID_WIDTH:", DEFAULT=0, RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource ( MAPL, ldas_first_ens_id, Label="FIRST_ENS_ID:", DEFAULT=0, RC=STATUS)
    VERIFY_(STATUS)
    ldas_ens_id = ldas_first_ens_id 

! -----------------------------------------------------------
! Get file names from configuration
! -----------------------------------------------------------
    if(NUM_LDAS_ENSEMBLE > 1) then
       !for GEOSldas, the comp_name should be vegdyn_exxxx....
       read(comp_name(9:), *) ldas_ens_id
    endif

    call MAPL_GetResource(MAPL, LAIFILE, label = 'LAI_FILE:', &
        default = 'lai.dat', RC=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, GRNFILE, label = 'GREEN_FILE:', &
         default = 'green.dat', RC=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, NDVIFILE, label = 'NDVI_FILE:', &
        default = 'ndvi.dat', RC=STATUS )
    VERIFY_(STATUS)

! get pointers to internal variables
! ----------------------------------
  
    call MAPL_GetPointer(INTERNAL,      ITY,      'ITY' , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL,      Z2CH,     'Z2CH', RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL,   ASCATZ0,  'ASCATZ0', RC=STATUS)
    VERIFY_(STATUS)

! get pointers to EXPORTS
! -----------------------

    call MAPL_GetPointer(EXPORT, LAI,   'LAI',    RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, GRN,   'GRN',    RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, ROOTL, 'ROOTL',  RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, NDVI,   'NDVI',  RC=STATUS)
    VERIFY_(STATUS)

    if (NUM_LDAS_ENSEMBLE > 1 .and. ldas_ens_id > ldas_first_ens_id) then
       if(associated(LAI))  LAI   = LAIens0  
       if(associated(GRN))  GRN   = GRNens0  
       if(associated(NDVI)) NDVI  = NDVIens0 
       if(associated(ROOTL))ROOTL = ROOTLens0
       call MAPL_TimerOff(MAPL,"TOTAL")
       RETURN_(ESMF_SUCCESS)
    endif

! Do the lai greeness and ndvi interpolation
! ------------------------------------------

    call ESMF_ClockGet  ( CLOCK, currTime=CURRENT_TIME, RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_ReadForcing(MAPL,'LAI',LAIFILE,CURRENT_TIME,LAI,ON_TILES=.true.,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_ReadForcing(MAPL,'GRN',GRNFILE,CURRENT_TIME,GRN,ON_TILES=.true.,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_ReadForcing(MAPL,'NDVI',NDVIFILE,CURRENT_TIME,NDVI,ON_TILES=.true.,RC=STATUS)
    VERIFY_(STATUS)

! Vegetation types used to index into tables
! Root length density no longer depends on time of year
! -----------------------------------------------------

    ROOTL = VGRT(nint(ITY))

    if (NUM_LDAS_ENSEMBLE > 1 .and. ldas_ens_id == ldas_first_ens_id) then
       LAIens0  => LAI
       GRNens0  => GRN
       NDVIens0 => NDVI
       ROOTLens0=> ROOTL
    endif

!  All done
! ---------

    call MAPL_TimerOff(MAPL,"TOTAL")
    RETURN_(ESMF_SUCCESS)

  end subroutine RUN

end module GEOS_VegdynGridCompMod
