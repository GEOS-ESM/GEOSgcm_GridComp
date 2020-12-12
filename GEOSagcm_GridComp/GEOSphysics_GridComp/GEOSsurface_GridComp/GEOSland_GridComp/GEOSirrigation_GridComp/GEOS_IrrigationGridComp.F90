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
! Exports from this routine are the instaneous values of the
! irrigation rates from 3 different irrigation methods on tilespace :
! 1) drip, 2) sprinkler and 3) flood. Land models (CATCH/CATCHCN)  use
! irrigation rates as a water inpput in water budget calculation.
! All exports and imports are stored on the
! tile grid inherited from the parent routine.\\
! 
! I. Parameter Class 1: Time and spatially dependent parameters 
! from a binary data file\\
! 
! EXPORTS:  SPRINKLERRATE, DRIPRATE, FLOODRATE\\
! 
! The gridded component stores the surrounding observations of 
! each parameter in the internal state.  All internals are static parameters.
!
! INTERNALS: IRRIGFRAC, PADDYFRAC, CROPIRRIGFRAC, IRRIGPLANT, IRRIGHARVEST,
!            IRRIGTYPE, SPRINKLERFR, DRIPFR, FLOODFR, LAIMIN, LAIMAX\\
!
! This GC imports soil parameters and root zone soil moisture from land models
!    to compute soil moisture state for IRRIGRATE calculation.
  
! IMPORTS: POROS, WPWET, VGWMAX, WCRZ \\
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
! Import states
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
         SHORT_NAME = 'SPRINKLERRATE'                         ,&
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

    character(len=ESMF_MAXSTR)         :: IAm
    integer                            :: STATUS
    character(len=ESMF_MAXSTR)         :: COMP_NAME

! Locals

    type (MAPL_MetaComp),     pointer  :: MAPL=>null()
    type (ESMF_State       )           :: INTERNAL

! INTERNAL pointers
 
    real, dimension(:),     pointer :: IRRIGFRAC
    real, dimension(:),     pointer :: PADDYFRAC
    real, dimension(:),     pointer :: SPRINKLERFR
    real, dimension(:),     pointer :: DRIPFR
    real, dimension(:),     pointer :: FLOODFR
    real, dimension(:),     pointer :: LAIMIN
    real, dimension(:),     pointer :: LAIMAX
    real, dimension(:,:),   pointer :: CROPIRRIGFRAC
    real, dimension(:,:),   pointer :: IRRIGTYPE
    real, dimension(:,:,:), pointer :: IRRIGPLANT
    real, dimension(:,:,:), pointer :: IRRIGHARVEST
    
! EXPORT pointers 

    real, dimension(:), pointer :: SPRINKLERRATE
    real, dimension(:), pointer :: DRIPRATE
    real, dimension(:), pointer :: FLOODRATE

! IMPORT pointers
    
    real, dimension(:), pointer :: POROS
    real, dimension(:), pointer :: WPWET
    real, dimension(:), pointer :: VGWMAX
    real, dimension(:), pointer :: WCRZ
  
! Time attributes 

    type(ESMF_Time) :: CURRENT_TIME
    logical, save   :: firsttime=.true.
    integer         :: AGCM_YY, AGCM_MM, AGCM_DD, AGCM_MI, AGCM_S, AGCM_HH, dofyr

! Others


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


! get pointers to internal variables
! ----------------------------------
  
    call MAPL_GetPointer(INTERNAL, IRRIGFRAC      ,'IRRIGFRAC',    RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, PADDYFRAC      ,'PADDYFRAC',    RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, CROPIRRIGFRAC  ,'CROPIRRIGFRAC',RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, IRRIGPLANT     ,'IRRIGPLANT',   RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, IRRIGHARVEST   ,'IRRIGHARVEST', RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, IRRIGTYPE      ,'IRRIGTYPE',    RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, SPRINKLERFR    ,'SPRINKLERFR',  RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, DRIPFR         ,'DRIPFR',       RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, FLOODFR        ,'FLOODFR',      RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, LAIMIN         ,'LAIMIN',       RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, LAIMAX         ,'LAIMAX',       RC=STATUS) ; VERIFY_(STATUS)

! get pointers to EXPORTS
! -----------------------

    call MAPL_GetPointer(EXPORT, SPRINKLERRATE, 'SPRINKLERRATE', RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, DRIPRATE,      'DRIPRATE',      RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, FLOODRATE,     'FLOODRATE',     RC=STATUS) ; VERIFY_(STATUS)

! get pointers to IMPORT variables
! --------------------------------

    call MAPL_GetPointer(IMPORT, POROS , 'POROS',  RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, WPWET , 'WPWET',  RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, VGWMAX, 'VGWMAX', RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, WCRZ  , 'WCRZ',   RC=STATUS) ; VERIFY_(STATUS)
    
    if (firsttime) then
       
       if(associated(SPRINKLERRATE)) SPRINKLERRATE = 0.
       if(associated(DRIPRATE))      DRIPRATE      = 0.
       if(associated(FLOODRATE))     FLOODRATE     = 0.
       firsttime = .false.
       call MAPL_TimerOff(MAPL,"TOTAL")
       RETURN_(ESMF_SUCCESS)
       
    endif

! Get time
! --------

    call ESMF_ClockGet  ( CLOCK, currTime=CURRENT_TIME, RC=STATUS )
    VERIFY_(STATUS)

    call ESMF_TimeGet  ( CURRENT_TIME, YY = AGCM_YY,       &
				       MM = AGCM_MM,       &
				       DD = AGCM_DD,       &
                                       H  = AGCM_HH,       &
                                       M  = AGCM_MI,       &
				       S  = AGCM_S ,       &  
				       dayOfYear = dofyr , &
				       rc=status )
    VERIFY_(STATUS)

! call the irrigation model 
! -------------------------
    
    call MAPL_TimerOff(MAPL,"TOTAL")
    RETURN_(ESMF_SUCCESS)

  end subroutine RUN

end module GEOS_IrrigationGridCompMod
