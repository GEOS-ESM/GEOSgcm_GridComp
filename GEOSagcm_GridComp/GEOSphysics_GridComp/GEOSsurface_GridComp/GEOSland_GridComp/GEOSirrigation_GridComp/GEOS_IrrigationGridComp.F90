!  $Id$

#include "MAPL_Generic.h"


!=============================================================================
module GEOS_IrrigationGridCompMod

!BOP

! !MODULE: GEOS_Irrigation -- child to the "Land" gridded component.  

!DESCRIPTION:
!   {\tt GEOS\_Irrigation} is a gridded component that performs the
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
  use IRRIGATION_MODULE ONLY: IRRIGATION, NUM_CROPS, NUM_SEASONS
  
  implicit none
  private

! !PUBLIC MEMBER FUNCTIONS:

  public SetServices

!EOP

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
         SHORT_NAME = 'IRRIGFRAC'                               ,&
         LONG_NAME  = 'fraction_of_irrigated_cropland'	        ,&
         UNITS      = '1'                                       ,&
         DIMS       = MAPL_DimsTileOnly                         ,&
         VLOCATION  = MAPL_VLocationNone                        ,&
         FRIENDLYTO = trim(COMP_NAME)                           ,&
         RC=STATUS  )
    VERIFY_(STATUS)  

    call MAPL_AddInternalSpec(GC                                ,&
         SHORT_NAME = 'PADDYFRAC'                               ,&
         LONG_NAME  = 'fraction_of_paddy_cropland'	        ,&
         UNITS      = '1'                                       ,&
         DIMS       = MAPL_DimsTileOnly                         ,&
         VLOCATION  = MAPL_VLocationNone                        ,&
         FRIENDLYTO = trim(COMP_NAME)                           ,&
         RC=STATUS  )
    VERIFY_(STATUS)
        
    call MAPL_AddInternalSpec(GC                                ,&
         SHORT_NAME = 'CROPIRRIGFRAC'                           ,&
         LONG_NAME  = 'Crop_irrigated_fraction'		        ,&
         UNITS      = '1'                                       ,&
         DIMS       = MAPL_DimsTileOnly                         ,&
         VLOCATION  = MAPL_VLocationNone                        ,&
         FRIENDLYTO = trim(COMP_NAME)                           ,&
         UNGRIDDED_DIMS = (/NUM_CROPS/)                         ,&
         RC=STATUS  )
    VERIFY_(STATUS)
    
    call MAPL_AddInternalSpec(GC                                ,&
         SHORT_NAME = 'IRRIGPLANT'                              ,&
         LONG_NAME  = 'DOY_start_planting'			,&
         UNITS      = 'day'                                     ,&
         DIMS       = MAPL_DimsTileOnly                         ,&
         VLOCATION  = MAPL_VLocationNone                        ,&
         FRIENDLYTO = trim(COMP_NAME)                           ,&
         UNGRIDDED_DIMS = (/NUM_SEASONS, NUM_CROPS/)            ,&
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC                                ,&
         SHORT_NAME = 'IRRIGHARVEST'                            ,&
         LONG_NAME  = 'DOY_end_harvesting'			,&
         UNITS      = 'day'                                     ,&
         DIMS       = MAPL_DimsTileOnly                         ,&
         VLOCATION  = MAPL_VLocationNone                        ,&
         FRIENDLYTO = trim(COMP_NAME)                           ,&
         UNGRIDDED_DIMS = (/NUM_SEASONS, NUM_CROPS/)            ,&
         RC=STATUS  )
    VERIFY_(STATUS)
    
    call MAPL_AddInternalSpec(GC                                ,&
         SHORT_NAME = 'IRRIGTYPE'                               ,&
         LONG_NAME  = 'Preferred_Irrig_method=(1)SPRINKLER_(2)DRIP_(3)FLOOD',&
         UNITS      = '1'                                       ,&
         DIMS       = MAPL_DimsTileOnly                         ,&
         VLOCATION  = MAPL_VLocationNone                        ,&
         FRIENDLYTO = trim(COMP_NAME)                           ,&
         UNGRIDDED_DIMS = (/NUM_CROPS/)                         ,&
         RC=STATUS  )
    VERIFY_(STATUS)
    
    call MAPL_AddInternalSpec(GC                                ,&
         SHORT_NAME = 'SPRINKLERFR'                             ,&
         LONG_NAME  = 'fraction_of_sprinkler_irrigation'	,&
         UNITS      = '1'                                       ,&
         DIMS       = MAPL_DimsTileOnly                         ,&
         VLOCATION  = MAPL_VLocationNone                        ,&
         FRIENDLYTO = trim(COMP_NAME)                           ,&
         RC=STATUS  )
    VERIFY_(STATUS)
    
    call MAPL_AddInternalSpec(GC                                ,&
         SHORT_NAME = 'DRIPFR'                                  ,&
         LONG_NAME  = 'fraction_of_drip_irrigation'		,&
         UNITS      = '1'                                       ,&
         DIMS       = MAPL_DimsTileOnly                         ,&
         VLOCATION  = MAPL_VLocationNone                        ,&
         FRIENDLYTO = trim(COMP_NAME)                           ,&
         RC=STATUS  )
    VERIFY_(STATUS)
    
    call MAPL_AddInternalSpec(GC                                ,&
         SHORT_NAME = 'FLOODFR'                                 ,&
         LONG_NAME  = 'fraction_of_flood_irrigation'	        ,&
         UNITS      = '1'                                       ,&
         DIMS       = MAPL_DimsTileOnly                         ,&
         VLOCATION  = MAPL_VLocationNone                        ,&
         FRIENDLYTO = trim(COMP_NAME)                           ,&
         RC=STATUS  )
    VERIFY_(STATUS)
    
    call MAPL_AddInternalSpec(GC                                ,&
         SHORT_NAME = 'LAIMIN'                                  ,&
         LONG_NAME  = 'Minimum_LAI_irrigated_crops'	        ,&
         UNITS      = '1'                                       ,&
         DIMS       = MAPL_DimsTileOnly                         ,&
         VLOCATION  = MAPL_VLocationNone                        ,&
         FRIENDLYTO = trim(COMP_NAME)                           ,&
         RC=STATUS  )
    VERIFY_(STATUS)
    
    call MAPL_AddInternalSpec(GC                                ,&
         SHORT_NAME = 'LAIMAX'                                  ,&
         LONG_NAME  = 'Maximum_LAI_irrigated_crops'             ,&
         UNITS      = '1'                                       ,&
         DIMS       = MAPL_DimsTileOnly                         ,&
         VLOCATION  = MAPL_VLocationNone                        ,&
         FRIENDLYTO = trim(COMP_NAME)                           ,&
         RC=STATUS  )
    VERIFY_(STATUS)  

! -----------------------------------------------------------
! Import Variables
! -----------------------------------------------------------    

    call MAPL_AddImportSpec(GC                                 ,&
         SHORT_NAME = 'POROS'                                  ,&         
         LONG_NAME  = 'soil_porosit'                           ,&
         UNITS      = '1'                                      ,&
         DIMS       = MAPL_DimsTileOnly                        ,&
         VLOCATION  = MAPL_VLocationNone                       ,&
         RC=STATUS  ) 
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC                                ,&
         SHORT_NAME = 'WPWET'                                 ,&         
         LONG_NAME  = 'wetness_at_wilting_point'              ,&
         UNITS      = '1'                                     ,&
         DIMS       = MAPL_DimsTileOnly                       ,&
         VLOCATION  = MAPL_VLocationNone                      ,&
         RC=STATUS  ) 
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC                                ,&
         SHORT_NAME = 'VGWMAX'                                ,&         
         LONG_NAME  = 'max_rootzone_water_content'            ,&
         UNITS      = 'kg m-2'                                ,&
         DIMS       = MAPL_DimsTileOnly                       ,&
         VLOCATION  = MAPL_VLocationNone                      ,&
         RC=STATUS  ) 
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC                                ,&
         SHORT_NAME = 'WCRZ'                                  ,&         
         LONG_NAME  = 'water_root_zone'                       ,&
         UNITS      = 'm3 m-3'                                ,&
         DIMS       = MAPL_DimsTileOnly                       ,&
         VLOCATION  = MAPL_VLocationNone                      ,&
         RC=STATUS  ) 
    VERIFY_(STATUS)


! -----------------------------------------------------------
! Export Variables
! -----------------------------------------------------------

    call MAPL_AddExportSpec(GC                                ,&
         SHORT_NAME = 'SPRINKERRATE'                          ,&
         LONG_NAME  = 'sprinkler_irrigation_rate'             ,&
         UNITS      = 'kg m-2 s-1'                            ,&
         DIMS       = MAPL_DimsTileOnly                       ,&
         VLOCATION  = MAPL_VLocationNone                      ,&
         RC=STATUS  )
    VERIFY_(STATUS)  
    
    call MAPL_AddExportSpec(GC                                ,&
         SHORT_NAME = 'DRIPRATE'                              ,&
         LONG_NAME  = 'drip_irrigation_rate'	              ,&
         UNITS      = 'kg m-2 s-1'                            ,&
         DIMS       = MAPL_DimsTileOnly                       ,&
         VLOCATION  = MAPL_VLocationNone                      ,&
         RC=STATUS  )
    VERIFY_(STATUS)  	 
    
    call MAPL_AddExportSpec(GC                                ,&
         SHORT_NAME = 'FLOODRATE'                             ,&
         LONG_NAME  = 'flood_irrigation_rate'                 ,&
         UNITS      = 'kg m-2 s-1'                            ,&
         DIMS       = MAPL_DimsTileOnly                       ,&
         VLOCATION  = MAPL_VLocationNone                      ,&
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
! RUN -- Run method for the irrigation component
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
       !for GEOSldas, the comp_name should be vegdynxxxx....
       read(comp_name(7:7+ens_id_width-1), *) ldas_ens_id
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

end module GEOS_IrrigationGridCompMod
