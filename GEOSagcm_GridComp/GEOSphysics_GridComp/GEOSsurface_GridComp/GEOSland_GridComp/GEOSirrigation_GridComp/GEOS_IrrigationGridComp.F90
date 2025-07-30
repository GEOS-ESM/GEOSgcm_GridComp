!  $Id$

#include "MAPL_Generic.h"


!=============================================================================
module GEOS_IrrigationGridCompMod

!BOP

! !MODULE: GEOS_Irrigation -- child to the "Land" gridded component.  

!DESCRIPTION:
!   {\tt GEOS\_Irrigation} is a gridded component that calculates
!   irrigation rates that are (optionally) used in the Catchment[CN] model.\\
!
! Exports from this routine are the instantaneous values of the
! irrigation rates from 4 different irrigation methods on tilespace:
! 1) drip, 2) sprinkler, 3) furrow, and 4) flood. 
! Catchment[CN] (optionally) uses these irrigation rates as inputs to
! the water budget calculations.  All imports and exports are stored 
! in tile space.  Soil parameters and root zone soil moisture are 
! imports from the Catchment[CN] gridded component.\\
! 
! Temporally and spatially varying irrigation model parameters are 
! from a binary data file.\\
!
! All internals are static parameters.\\
!
! IMPORTS:   POROS, WPWET, VGWMAX, WCRZ, LAI\\
!
! EXPORTS:   IRRG_RATE_SPR, IRRG_RATE_DRP, IRRG_RATE_FRW, IRRG_RATE_PDY, IRRG_RATE_TOT\\ 
!  
! INTERNALS: IRRG_IRRIGFRAC, IRRG_PADDYFRAC, IRRG_CROPIRRIGFRAC, IRRG_DOY_PLANT, IRRG_DOY_HARVEST,
!            IRRG_TYPE, IRRG_IRRIGFRAC_SPR, IRRG_IRRIGFRAC_DRP, IRRG_IRRIGFRAC_FRW, IRRG_LAIMIN, IRRG_LAIMAX\\
!
! OPTIONAL INTERNALS:  SRATE, DRATE, FRATE\\
!  
! !USES: 

  use ESMF
  use MAPL
  use IRRIGATION_MODULE
  
  implicit none
  private

! !PUBLIC MEMBER FUNCTIONS:

  public SetServices

  integer    :: IRRG_METHOD, IRRG_TRIGGER
  integer    :: RUN_IRRIG
  
  type IRRIG_WRAP
     type (irrig_params), pointer   :: PTR => null()
  end type IRRIG_WRAP

contains

!BOP

! !IROUTINE: SetServices -- Sets ESMF services for this component

! !INTERFACE:

  subroutine SetServices ( GC, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    integer, optional                  :: RC  ! return code

! !DESCRIPTION: This version uses the MAPL\_GenericSetServices. This function sets
!   the Initialize and Finalize services, as well as allocating
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
    type(ESMF_Config)                       :: SCF
    character(len=ESMF_MAXSTR)              :: SURFRC

    type (irrigation_model), pointer        :: IM => null()
    type (IRRIG_wrap)                       :: wrap

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
! Get the configuration
! -----------------------------------------------------------
    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)

! -----------------------------------------------------------
! Get runtime switches
! -----------------------------------------------------------

    call MAPL_GetResource(MAPL, SURFRC, label='SURFRC:', default='GEOS_SurfaceGridComp.rc', RC=STATUS); VERIFY_(STATUS)

    SCF = ESMF_ConfigCreate(rc=status)                                                                ; VERIFY_(STATUS)

    call ESMF_ConfigLoadFile(    SCF, SURFRC, rc=status)                                              ; VERIFY_(STATUS)

    call ESMF_ConfigGetAttribute(SCF, label='RUN_IRRIG:'   , value=RUN_IRRIG   , DEFAULT=0, __RC__ )
    call ESMF_ConfigGetAttribute(SCF, label='IRRG_TRIGGER:', value=IRRG_TRIGGER, DEFAULT=0, __RC__ )    
    call ESMF_ConfigGetAttribute(SCF, label='IRRG_METHOD:' , value=IRRG_METHOD , DEFAULT=0, __RC__ )
   
    call ESMF_ConfigDestroy     (SCF, __RC__)

    ! Leave GEOSirrigation_GridComp if RUN_IRRIG == 0
    if(RUN_IRRIG == 0) then
       RETURN_(ESMF_SUCCESS)
    endif
    
! -----------------------------------------------------------
! Set the the Initialize and Run entry point
! -----------------------------------------------------------

    call MAPL_GridCompSetEntryPoint (GC, ESMF_METHOD_INITIALIZE, Initialize, RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GridCompSetEntryPoint (GC, ESMF_METHOD_RUN,        Run,        RC=STATUS)
    VERIFY_(STATUS)

! BOS
    
! -----------------------------------------------------------
! Internal State 
! -----------------------------------------------------------

    call MAPL_AddInternalSpec(GC                              ,&
         SHORT_NAME = 'IRRG_IRRIGFRAC'                        ,&
         LONG_NAME  = 'fraction_of_irrigated_cropland'	      ,&
         UNITS      = '1'                                     ,&
         DIMS       = MAPL_DimsTileOnly                       ,&
         VLOCATION  = MAPL_VLocationNone                      ,&
         FRIENDLYTO = trim(COMP_NAME)                         ,&
         RESTART    = MAPL_RestartRequired                    ,&
         RC=STATUS  )
    VERIFY_(STATUS)  
    
    call MAPL_AddInternalSpec(GC                              ,&
         SHORT_NAME = 'IRRG_PADDYFRAC'                        ,&
         LONG_NAME  = 'fraction_of_paddy_cropland'	      ,&
         UNITS      = '1'                                     ,&
         DIMS       = MAPL_DimsTileOnly                       ,&
         VLOCATION  = MAPL_VLocationNone                      ,&
         FRIENDLYTO = trim(COMP_NAME)                         ,&
         RESTART    = MAPL_RestartRequired                    ,&
         RC=STATUS  )
    VERIFY_(STATUS)
       
    call MAPL_AddInternalSpec(GC                              ,&
         SHORT_NAME = 'IRRG_CROPIRRIGFRAC'                    ,&
         LONG_NAME  = 'Crop_irrigated_fraction'		      ,&
         UNITS      = '1'                                     ,&
         DIMS       = MAPL_DimsTileOnly                       ,&
         VLOCATION  = MAPL_VLocationNone                      ,&
         FRIENDLYTO = trim(COMP_NAME)                         ,&
         UNGRIDDED_DIMS = (/IRRG_NCROPS/)                     ,&
         RESTART    = MAPL_RestartRequired                    ,&
         RC=STATUS  )
    VERIFY_(STATUS)
       
    call MAPL_AddInternalSpec(GC                              ,&
         SHORT_NAME = 'IRRG_DOY_PLANT'                        ,&
         LONG_NAME  = 'DOY_start_planting'		      ,&
         UNITS      = 'day'                                   ,&
         DIMS       = MAPL_DimsTileOnly                       ,&
         VLOCATION  = MAPL_VLocationNone                      ,&
         FRIENDLYTO = trim(COMP_NAME)                         ,&
         UNGRIDDED_DIMS = (/IRRG_NSEASONS, IRRG_NCROPS/)      ,&
         RESTART    = MAPL_RestartRequired                    ,&
         RC=STATUS  )
    VERIFY_(STATUS)
    
    call MAPL_AddInternalSpec(GC                              ,&
         SHORT_NAME = 'IRRG_DOY_HARVEST'                      ,&
         LONG_NAME  = 'DOY_end_harvesting'                    ,&
         UNITS      = 'day'                                   ,&
         DIMS       = MAPL_DimsTileOnly                       ,&
         VLOCATION  = MAPL_VLocationNone                      ,&
         FRIENDLYTO = trim(COMP_NAME)                         ,&
         UNGRIDDED_DIMS = (/IRRG_NSEASONS, IRRG_NCROPS/)      ,&
         RESTART    = MAPL_RestartRequired                    ,&
         RC=STATUS  )
    VERIFY_(STATUS)
    
    call MAPL_AddInternalSpec(GC                              ,&
         SHORT_NAME = 'IRRG_TYPE'                             ,&
         LONG_NAME  = 'Preferred_Irrig_method=(0)CONCURRENT_(1)SPRINKLER_(2)DRIP_(3)FLOOD_(negative)AVOID',&
         UNITS      = '1'                                     ,&
         DIMS       = MAPL_DimsTileOnly                       ,&
         VLOCATION  = MAPL_VLocationNone                      ,&
         FRIENDLYTO = trim(COMP_NAME)                         ,&
         UNGRIDDED_DIMS = (/IRRG_NCROPS/)                     ,&
         RESTART    = MAPL_RestartRequired                    ,&
         RC=STATUS  )
    VERIFY_(STATUS)
    
    call MAPL_AddInternalSpec(GC                              ,&
         SHORT_NAME = 'IRRG_IRRIGFRAC_SPR'                    ,&
         LONG_NAME  = 'fraction_of_sprinkler_irrigation'      ,&
         UNITS      = '1'                                     ,&
         DIMS       = MAPL_DimsTileOnly                       ,&
         VLOCATION  = MAPL_VLocationNone                      ,&
         FRIENDLYTO = trim(COMP_NAME)                         ,&
         RC=STATUS  )
    VERIFY_(STATUS)
    
    call MAPL_AddInternalSpec(GC                              ,&
         SHORT_NAME = 'IRRG_IRRIGFRAC_DRP'                    ,&
         LONG_NAME  = 'fraction_of_drip_irrigation'	      ,&
         UNITS      = '1'                                     ,&
         DIMS       = MAPL_DimsTileOnly                       ,&
         VLOCATION  = MAPL_VLocationNone                      ,&
         FRIENDLYTO = trim(COMP_NAME)                         ,&
         RESTART    = MAPL_RestartRequired                    ,&
         RC=STATUS  )
    VERIFY_(STATUS)
    
    call MAPL_AddInternalSpec(GC                              ,&
         SHORT_NAME = 'IRRG_IRRIGFRAC_FRW'                    ,&
         LONG_NAME  = 'fraction_of_flood_irrigation'	      ,&
         UNITS      = '1'                                     ,&
         DIMS       = MAPL_DimsTileOnly                       ,&
         VLOCATION  = MAPL_VLocationNone                      ,&
         FRIENDLYTO = trim(COMP_NAME)                         ,&
         RESTART    = MAPL_RestartRequired                    ,&
         RC=STATUS  )
    VERIFY_(STATUS)
    
    call MAPL_AddInternalSpec(GC                              ,&
         SHORT_NAME = 'IRRG_LAIMIN'                           ,&
         LONG_NAME  = 'Minimum_LAI_irrigated_crops'	      ,&
         UNITS      = '1'                                     ,&
         DIMS       = MAPL_DimsTileOnly                       ,&
         VLOCATION  = MAPL_VLocationNone                      ,&
         FRIENDLYTO = trim(COMP_NAME)                         ,&
         RESTART    = MAPL_RestartRequired                    ,&
         RC=STATUS  )
    VERIFY_(STATUS)
    
    call MAPL_AddInternalSpec(GC                              ,&
         SHORT_NAME = 'IRRG_LAIMAX'                           ,&
         LONG_NAME  = 'Maximum_LAI_irrigated_crops'           ,&
         UNITS      = '1'                                     ,&
         DIMS       = MAPL_DimsTileOnly                       ,&
         VLOCATION  = MAPL_VLocationNone                      ,&
         FRIENDLYTO = trim(COMP_NAME)                         ,&
         RESTART    = MAPL_RestartRequired                    ,&
         RC=STATUS  )
    VERIFY_(STATUS)  

    ! NOTE: UNGRIDDED_DIMS for SRATE, DRATE, and FRATE internal specs depends on IRRG_TRIGGER
    
    if (IRRG_TRIGGER == 0) then

       ! UNGRIDDED_DIMS = 2: irrigated crops (sprinkler/drip/furrow) and paddy (in order)

       call MAPL_AddInternalSpec(GC                              ,&
            SHORT_NAME = 'SRATE'                                 ,&
            LONG_NAME  ='crop_specific_irrigation_flux_sprinkler',&
            UNITS      = 'kg m-2 s-1'                            ,&
            DIMS       = MAPL_DimsTileOnly                       ,&
            VLOCATION  = MAPL_VLocationNone                      ,&
            FRIENDLYTO = trim(COMP_NAME)                         ,&
            RESTART    = MAPL_RestartOptional                    ,&
            UNGRIDDED_DIMS = (/2/)                               ,&
            RC=STATUS  )
       VERIFY_(STATUS)  
       
       call MAPL_AddInternalSpec(GC                              ,&
            SHORT_NAME = 'DRATE'                                 ,&
            LONG_NAME  = 'crop_specific_irrigation_flux_drip'    ,&
            UNITS      = 'kg m-2 s-1'                            ,&
            DIMS       = MAPL_DimsTileOnly                       ,&
            VLOCATION  = MAPL_VLocationNone                      ,&
            FRIENDLYTO = trim(COMP_NAME)                         ,&
            RESTART    = MAPL_RestartOptional                    ,&
            UNGRIDDED_DIMS = (/2/)                               ,&
            RC=STATUS  )
       VERIFY_(STATUS)  	 
       
       call MAPL_AddInternalSpec(GC                              ,&
            SHORT_NAME = 'FRATE'                                 ,&
            LONG_NAME  = 'crop_specific_irrigation_flux_paddy'   ,&
            UNITS      = 'kg m-2 s-1'                            ,&
            DIMS       = MAPL_DimsTileOnly                       ,&
            VLOCATION  = MAPL_VLocationNone                      ,&
            FRIENDLYTO = trim(COMP_NAME)                         ,&
            RESTART    = MAPL_RestartOptional                    ,&
            UNGRIDDED_DIMS = (/2/)                               ,&
            RC=STATUS  )
       VERIFY_(STATUS)

    elseif (IRRG_TRIGGER == 1) then

       ! UNGRIDDED_DIMS = 26 crops of crop calendar
       
       call MAPL_AddInternalSpec(GC                              ,&
            SHORT_NAME = 'SRATE'                                 ,&
            LONG_NAME  ='crop_specific_irrigation_flux_sprinkler',&
            UNITS      = 'kg m-2 s-1'                            ,&
            DIMS       = MAPL_DimsTileOnly                       ,&
            VLOCATION  = MAPL_VLocationNone                      ,&
            FRIENDLYTO = trim(COMP_NAME)                         ,&
            RESTART    = MAPL_RestartOptional                    ,&
            UNGRIDDED_DIMS = (/IRRG_NCROPS/)                     ,&
            RC=STATUS  )
       VERIFY_(STATUS)  
       
       call MAPL_AddInternalSpec(GC                              ,&
            SHORT_NAME = 'DRATE'                                 ,&
            LONG_NAME  = 'crop_specific_irrigation_flux_drip'    ,&
            UNITS      = 'kg m-2 s-1'                            ,&
            DIMS       = MAPL_DimsTileOnly                       ,&
            VLOCATION  = MAPL_VLocationNone                      ,&
            FRIENDLYTO = trim(COMP_NAME)                         ,&
            RESTART    = MAPL_RestartOptional                    ,&
            UNGRIDDED_DIMS = (/IRRG_NCROPS/)                     ,&
            RC=STATUS  )
       VERIFY_(STATUS)  	 
       
       call MAPL_AddInternalSpec(GC                              ,&
            SHORT_NAME = 'FRATE'                                 ,&
            LONG_NAME  = 'crop_specific_irrigation_flux_paddy'   ,&
            UNITS      = 'kg m-2 s-1'                            ,&
            DIMS       = MAPL_DimsTileOnly                       ,&
            VLOCATION  = MAPL_VLocationNone                      ,&
            FRIENDLYTO = trim(COMP_NAME)                         ,&
            RESTART    = MAPL_RestartOptional                    ,&
            UNGRIDDED_DIMS = (/IRRG_NCROPS/)                     ,&
            RC=STATUS  )
       VERIFY_(STATUS)
       
    endif
    
! -----------------------------------------------------------
! Export state
! -----------------------------------------------------------

    call MAPL_AddExportSpec(GC                                ,&
         SHORT_NAME = 'IRRG_RATE_SPR'                         ,&
         LONG_NAME  = 'irrigation_flux_sprinkler'             ,&
         UNITS      = 'kg m-2 s-1'                            ,&
         DIMS       = MAPL_DimsTileOnly                       ,&
         VLOCATION  = MAPL_VLocationNone                      ,&
         RC=STATUS  )
    VERIFY_(STATUS)  
    
    call MAPL_AddExportSpec(GC                                ,&
         SHORT_NAME = 'IRRG_RATE_DRP'                         ,&
         LONG_NAME  = 'irrigation_flux_drip'	              ,&
         UNITS      = 'kg m-2 s-1'                            ,&
         DIMS       = MAPL_DimsTileOnly                       ,&
         VLOCATION  = MAPL_VLocationNone                      ,&
         RC=STATUS  )
    VERIFY_(STATUS)  	 
    
     call MAPL_AddExportSpec(GC                               ,&
         SHORT_NAME = 'IRRG_RATE_FRW'                         ,&
         LONG_NAME  = 'irrigation_flux_furrow'                ,&
         UNITS      = 'kg m-2 s-1'                            ,&
         DIMS       = MAPL_DimsTileOnly                       ,&
         VLOCATION  = MAPL_VLocationNone                      ,&
         RC=STATUS  )
    VERIFY_(STATUS)
  
    call MAPL_AddExportSpec(GC                                ,&
         SHORT_NAME = 'IRRG_RATE_PDY'                         ,&
         LONG_NAME  = 'irrigation_flux_paddy'                 ,&
         UNITS      = 'kg m-2 s-1'                            ,&
         DIMS       = MAPL_DimsTileOnly                       ,&
         VLOCATION  = MAPL_VLocationNone                      ,&
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC                                ,&
         SHORT_NAME = 'IRRG_RATE_TOT'                         ,&
         LONG_NAME  = 'irrigation_flux_total'                 ,&
         UNITS      = 'kg m-2 s-1'                            ,&
         DIMS       = MAPL_DimsTileOnly                       ,&
         VLOCATION  = MAPL_VLocationNone                      ,&
         RC=STATUS  )
    VERIFY_(STATUS)
      
! -----------------------------------------------------------
! Import states
! -----------------------------------------------------------    

    call MAPL_AddImportSpec(GC                                ,&
         SHORT_NAME = 'POROS'                                 ,&         
         LONG_NAME  = 'soil_porosity'                         ,&
         UNITS      = 'm3 m-3'                                ,&
         DIMS       = MAPL_DimsTileOnly                       ,&
         VLOCATION  = MAPL_VLocationNone                      ,&
         RC=STATUS  ) 
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC                                ,&
         SHORT_NAME = 'WPWET'                                 ,&         
         LONG_NAME  = 'soil_wilting_point_in_degree_of_saturation_units'    ,&
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
         LONG_NAME  = 'soil_moisture_rootzone'                ,&
         UNITS      = 'm3 m-3'                                ,&
         DIMS       = MAPL_DimsTileOnly                       ,&
         VLOCATION  = MAPL_VLocationNone                      ,&
         RC=STATUS  ) 
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC                                ,&
         SHORT_NAME = 'LAI'                                   ,&
         LONG_NAME  = 'leaf_area_index'                       ,&
         UNITS      = '1'                                     ,&
         DIMS       = MAPL_DimsTileOnly                       ,&
         VLOCATION  = MAPL_VLocationNone                      ,&
         RC=STATUS  )
    VERIFY_(STATUS)
          
! Allocate this instance of the internal state and put it in wrapper.
! -------------------------------------------------------------------

    allocate(IM,   stat=status )
    VERIFY_(STATUS)
    call IM%init_model (SURFRC)
    wrap%ptr => IM%irrig_params

! Save pointer to the wrapped internal state in the GC
! ----------------------------------------------------

    call ESMF_UserCompSetInternalState ( GC, 'irrigation_state',wrap,status )
    VERIFY_(STATUS)

! Clocks
!-------

    call MAPL_TimerAdd(GC, name="INITIALIZE"    ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC, name="RUN"           ,RC=STATUS)
    VERIFY_(STATUS)

!------------------------------------------------------------
! Set generic init and final methods
!------------------------------------------------------------

    call MAPL_GenericSetServices(GC, RC=STATUS)
    VERIFY_(STATUS)
    
    RETURN_(ESMF_SUCCESS)
  
  end subroutine SetServices

  ! -----------------------------------------------------------
  ! INITIALIZE -- Initialize method for the irrigation component
  ! -----------------------------------------------------------

  subroutine INITIALIZE (GC,IMPORT, EXPORT, CLOCK, RC )
    
    ! ARGUMENTS:
    
    type(ESMF_GridComp), intent(inout) :: GC    
    type(ESMF_State),    intent(inout) :: IMPORT
    type(ESMF_State),    intent(inout) :: EXPORT
    type(ESMF_Clock),    intent(inout) :: CLOCK
    integer, optional,   intent(  out) :: RC
        
    ! ErrLog Variables

    character(len=ESMF_MAXSTR)      :: IAm="Initialize"
    integer                         :: STATUS
    character(len=ESMF_MAXSTR)      :: COMP_NAME

    ! Locals

    type (MAPL_MetaComp),   pointer :: MAPL=>null()
    type (ESMF_State   )            :: INTERNAL

    ! INTERNAL pointers

    real, dimension(:),     pointer :: IRRG_IRRIGFRAC
    real, dimension(:),     pointer :: IRRG_PADDYFRAC
    real, dimension(:,:),   pointer :: IRRG_CROPIRRIGFRAC
    real, dimension(:,:),   pointer :: IRRG_TYPE
    real, dimension(:,:),   pointer :: SRATE
    real, dimension(:,:),   pointer :: DRATE
    real, dimension(:,:),   pointer :: FRATE

! EXPORT pointers
    
    real, dimension(:),     pointer :: IRRG_RATE_SPR
    real, dimension(:),     pointer :: IRRG_RATE_DRP
    real, dimension(:),     pointer :: IRRG_RATE_FRW
    real, dimension(:),     pointer :: IRRG_RATE_PDY
    real, dimension(:),     pointer :: IRRG_RATE_TOT

    type(irrigation_model), pointer :: IM
    type (IRRIG_wrap)               :: wrap
    
! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    call ESMF_GridCompGet(GC, name=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
  
    Iam = trim(COMP_NAME) // "Initialize"

    call MAPL_GenericInitialize ( GC, import, export, clock, rc=status )
    VERIFY_(STATUS)
    call ESMF_UserCompGetInternalState ( GC, 'irrigation_state',wrap,status )
    VERIFY_(STATUS)    
    allocate (IM)
    IM = irrigation_model()
    IM%irrig_params = wrap%ptr

    ! Get my internal MAPL_Generic state
    ! -----------------------------------------------------------

    call MAPL_GetObjectFromGC(GC, MAPL, STATUS)
    VERIFY_(STATUS)
    
    call MAPL_TimerOn(MAPL,"INITIALIZE")

    call MAPL_Get(MAPL, INTERNAL_ESMF_STATE=INTERNAL, RC=STATUS)
    VERIFY_(STATUS)     

    ! get pointers to internal variables
    ! ----------------------------------

    call MAPL_GetPointer(INTERNAL, IRRG_IRRIGFRAC     ,'IRRG_IRRIGFRAC',             RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, IRRG_PADDYFRAC     ,'IRRG_PADDYFRAC',             RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, IRRG_CROPIRRIGFRAC ,'IRRG_CROPIRRIGFRAC',         RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, IRRG_TYPE          ,'IRRG_TYPE',                  RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, SRATE              ,'SRATE',        ALLOC=.true., RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, DRATE              ,'DRATE',        ALLOC=.true., RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, FRATE              ,'FRATE',        ALLOC=.true., RC=STATUS) ; VERIFY_(STATUS)
    
    ! get pointers to EXPORT variable
    ! -------------------------------
    call MAPL_GetPointer(EXPORT, IRRG_RATE_SPR,       'IRRG_RATE_SPR', ALLOC=.true., RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, IRRG_RATE_DRP,       'IRRG_RATE_DRP', ALLOC=.true., RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, IRRG_RATE_FRW,       'IRRG_RATE_FRW', ALLOC=.true., RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, IRRG_RATE_PDY,       'IRRG_RATE_PDY', ALLOC=.true., RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, IRRG_RATE_TOT,       'IRRG_RATE_TOT',               RC=STATUS) ; VERIFY_(STATUS)

    ! Update IRRG_IRRIGFRAC and IRRG_PADDYFRAC for applications that are run on regular tiles in 
    ! which IRRG_IRRIGFRAC and IRRG_PADDYFRAC in BCs are fractions.
    ! The irrigation model would run on tiles with IRRG_IRRIGFRAC + IRRG_PADDYFRAC > IRRG_FRAC_THRES (default is 0.01).

    where (IRRG_IRRIGFRAC + IRRG_PADDYFRAC > IM%irrig_thres)

       ! uncomment the following block to assign the entire cell to the largest fraction:
       
       ! where (IRRG_PADDYFRAC >= IRRG_IRRIGFRAC)
       !    IRRG_PADDYFRAC = 1.
       !    IRRG_IRRIGFRAC = 0.
       ! elsewhere
       !    IRRG_PADDYFRAC = 0.
       !    IRRG_IRRIGFRAC = 1.
       ! endwhere

    elsewhere
       IRRG_PADDYFRAC = 0.
       IRRG_IRRIGFRAC = 0.
    endwhere

    if (IRRG_TRIGGER == 0) then

       ! LAI based trigger: scale soil moisture to LAI seasonal cycle
       ! ============================================================
       
       call IM%update_irates (IRRG_RATE_SPR,IRRG_RATE_DRP,IRRG_RATE_FRW,IRRG_RATE_PDY, & 
            IRRG_IRRIGFRAC,IRRG_PADDYFRAC,SRATE,DRATE,FRATE)
       
    else
       
       ! crop calendar based irrigation
       ! ==============================

       call IM%update_irates (IRRG_RATE_SPR,IRRG_RATE_DRP,IRRG_RATE_FRW,IRRG_RATE_PDY, &
            IRRG_CROPIRRIGFRAC,SRATE,DRATE,FRATE)
       
    endif

    ! Scale computed IRRG_RATE_SPR, IRRG_RATE_DRP, IRRG_RATE_FRW, and IRRG_RATE_PDY to the total
    ! irrigated tile fraction before exporting to Catchment[CN].

    IRRG_RATE_SPR = IRRG_RATE_SPR * IRRG_IRRIGFRAC 
    IRRG_RATE_DRP = IRRG_RATE_DRP * IRRG_IRRIGFRAC 
    IRRG_RATE_FRW = IRRG_RATE_FRW * IRRG_IRRIGFRAC 
    IRRG_RATE_PDY = IRRG_RATE_PDY * IRRG_PADDYFRAC 

    if(associated(IRRG_RATE_TOT))  IRRG_RATE_TOT = IRRG_RATE_SPR + IRRG_RATE_DRP + IRRG_RATE_FRW +IRRG_RATE_PDY
    
    call MAPL_TimerOff(MAPL,"INITIALIZE")
    RETURN_(ESMF_SUCCESS)
    
  end subroutine INITIALIZE
    
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

! ErrLog Variables

    character(len=ESMF_MAXSTR)      :: IAm
    integer                         :: STATUS
    character(len=ESMF_MAXSTR)      :: COMP_NAME

! Locals

    type (MAPL_MetaComp),   pointer :: MAPL=>null()
    type (ESMF_State   )            :: INTERNAL

! INTERNAL pointers

    real, dimension(:),     pointer :: IRRG_IRRIGFRAC
    real, dimension(:),     pointer :: IRRG_PADDYFRAC
    real, dimension(:),     pointer :: IRRG_IRRIGFRAC_SPR
    real, dimension(:),     pointer :: IRRG_IRRIGFRAC_DRP
    real, dimension(:),     pointer :: IRRG_IRRIGFRAC_FRW
    real, dimension(:),     pointer :: IRRG_LAIMIN
    real, dimension(:),     pointer :: IRRG_LAIMAX
    real, dimension(:,:),   pointer :: IRRG_CROPIRRIGFRAC
    real, dimension(:,:),   pointer :: IRRG_TYPE
    real, dimension(:,:,:), pointer :: IRRG_DOY_PLANT
    real, dimension(:,:,:), pointer :: IRRG_DOY_HARVEST
    real, dimension(:,:),   pointer :: SRATE
    real, dimension(:,:),   pointer :: DRATE
    real, dimension(:,:),   pointer :: FRATE

! EXPORT pointers
    
    real, dimension(:),     pointer :: IRRG_RATE_SPR
    real, dimension(:),     pointer :: IRRG_RATE_DRP
    real, dimension(:),     pointer :: IRRG_RATE_FRW    
    real, dimension(:),     pointer :: IRRG_RATE_PDY    
    real, dimension(:),     pointer :: IRRG_RATE_TOT
    
! IMPORT pointers
    
    real, dimension(:),     pointer :: POROS
    real, dimension(:),     pointer :: WPWET
    real, dimension(:),     pointer :: VGWMAX
    real, dimension(:),     pointer :: WCRZ
    real, dimension(:),     pointer :: LAI
  
! Time attributes 

    type(ESMF_Time)                 :: CURRENT_TIME
    integer                         :: AGCM_YY, AGCM_MM, AGCM_DD, AGCM_MI, AGCM_S, AGCM_HH, dofyr

! Others/ Locals

    type(irrigation_model), pointer :: IM
    type (IRRIG_wrap)               :: wrap
    real, dimension(:),     pointer :: lons
    integer                         :: ntiles, n
    real, dimension(:), allocatable :: local_hour, SMWP, SMSAT, SMREF, SMCNT, RZDEF
    real                            :: DT, T1, T2

    ! Get the target components name and set-up traceback handle.
    ! -----------------------------------------------------------

    call ESMF_GridCompGet(GC, name=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
  
    Iam = trim(COMP_NAME) // "Run"

    call ESMF_UserCompGetInternalState ( GC, 'irrigation_state',wrap,status )
    VERIFY_(STATUS)    
    allocate (IM)
    IM = irrigation_model()
    IM%irrig_params = wrap%ptr

    ! Get my internal MAPL_Generic state
    ! -----------------------------------------------------------

    call MAPL_GetObjectFromGC(GC, MAPL, STATUS)
    VERIFY_(STATUS)

    call MAPL_Get(MAPL, HEARTBEAT = DT, RC=STATUS)
    VERIFY_(STATUS)
    
    call MAPL_TimerOn(MAPL,"RUN")

    call MAPL_Get(MAPL, INTERNAL_ESMF_STATE=INTERNAL, RC=STATUS)
    VERIFY_(STATUS) 

    ! get pointers to internal variables
    ! ----------------------------------

    call MAPL_GetPointer(INTERNAL, IRRG_IRRIGFRAC     ,'IRRG_IRRIGFRAC',              RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, IRRG_PADDYFRAC     ,'IRRG_PADDYFRAC',              RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, IRRG_CROPIRRIGFRAC ,'IRRG_CROPIRRIGFRAC',          RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, IRRG_DOY_PLANT     ,'IRRG_DOY_PLANT',              RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, IRRG_DOY_HARVEST   ,'IRRG_DOY_HARVEST',            RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, IRRG_TYPE          ,'IRRG_TYPE',                   RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, IRRG_IRRIGFRAC_SPR ,'IRRG_IRRIGFRAC_SPR',          RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, IRRG_IRRIGFRAC_DRP ,'IRRG_IRRIGFRAC_DRP',          RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, IRRG_IRRIGFRAC_FRW ,'IRRG_IRRIGFRAC_FRW',          RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, IRRG_LAIMIN        ,'IRRG_LAIMIN',                 RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, IRRG_LAIMAX        ,'IRRG_LAIMAX',                 RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, SRATE              ,'SRATE',         ALLOC=.true., RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, DRATE              ,'DRATE',         ALLOC=.true., RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, FRATE              ,'FRATE',         ALLOC=.true., RC=STATUS) ; VERIFY_(STATUS)
    
    ! get pointers to EXPORT variable
    ! -------------------------------
    call MAPL_GetPointer(EXPORT, IRRG_RATE_SPR        ,'IRRG_RATE_SPR', ALLOC=.true., RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, IRRG_RATE_DRP        ,'IRRG_RATE_DRP', ALLOC=.true., RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, IRRG_RATE_FRW        ,'IRRG_RATE_FRW', ALLOC=.true., RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, IRRG_RATE_PDY        ,'IRRG_RATE_PDY', ALLOC=.true., RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, IRRG_RATE_TOT        ,'IRRG_RATE_TOT',               RC=STATUS) ; VERIFY_(STATUS)
  
    
    
    ! get pointers to IMPORT variables
    ! --------------------------------

    call MAPL_GetPointer(IMPORT, POROS                ,'POROS',                       RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, WPWET                ,'WPWET',                       RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, VGWMAX               ,'VGWMAX',                      RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, WCRZ                 ,'WCRZ',                        RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, LAI                  ,'LAI',                         RC=STATUS) ; VERIFY_(STATUS)
        
    ! Get time and parameters from local state
    ! ----------------------------------------

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
    
    call MAPL_Get (MAPL, TILELONS = LONS,                  &                               ! longitude in [radians]
         INTERNAL_ESMF_STATE = INTERNAL, RC=STATUS )
    VERIFY_(STATUS)
    

    ! call irrigation model 
    ! ---------------------

    NTILES = SIZE (LONS)
    
    allocate (local_hour (1:NTILES))
    allocate (SMWP       (1:NTILES))
    allocate (SMSAT      (1:NTILES))
    allocate (SMREF      (1:NTILES))
    allocate (SMCNT      (1:NTILES))
    allocate (RZDEF      (1:NTILES))

    ! soil moisture state
    SMWP  =  VGWMAX        * WPWET   ! RZ soil moisture content at wilting point [mm]
    SMSAT =  VGWMAX                  ! RZ soil moisture at saturation  [mm]
    SMCNT = (VGWMAX/POROS) * WCRZ    ! actual RZ soil moisture content [mm]

    DO N = 1, NTILES

       ! local time [hour]

       local_hour(n) = AGCM_HH + AGCM_MI / 60. + AGCM_S / 3600. + 12.* (lons(n)/MAPL_PI)   ! longitude in [radians]
       IF (local_hour(n) >= 24.) local_hour(n) = local_hour(n) - 24.
       IF (local_hour(n) <   0.) local_hour(n) = local_hour(n) + 24.
       T1 = CEILING (local_hour(n))     - DT/3600.
       T2 = FLOOR   (local_hour(n) + 1) + DT/3600.
       if((local_hour(n) >= T1).and.(local_hour(n) < T2))then
          local_hour(n) = real(NINT(local_hour(n)))
       end if

       ! The reference soil moisture content is set to lower tercile of the root zone soil 
       ! moisture range [mm] to be consistent with ASTRFR = 0.333 used in Catchment[CN].
       !
       ! Note on choice of SMREF:
       !    Perhaps, soil field capacity (FIELDCAP) is the desired parameter here - the upper limit
       !    of water content that the soil can hold for plants after excess water drained off downward 
       !    quickly.  In future, could switch to FIELDCAP, which has already been derived
       !    on tiles and is available in the irrigation_IMxJM_DL.dat file.
       
       SMREF (n) = VGWMAX (n) * (wpwet (n) + (1. - wpwet (n))/ 3.)

       ! rootzone moisture deficit to reach complete soil saturation for paddy [mm]

       RZDEF (n) = MAX(SMSAT(n) - SMCNT(n), 0.)  
                                                 
    END DO
        
    if (IRRG_TRIGGER == 0) then

       ! LAI based trigger: scale soil moisture to LAI seasonal cycle
       ! ============================================================
                    
       call IM%run_model(IRRG_METHOD, local_hour,                                                       &
            IRRG_IRRIGFRAC, IRRG_PADDYFRAC, IRRG_IRRIGFRAC_SPR, IRRG_IRRIGFRAC_DRP, IRRG_IRRIGFRAC_FRW, &           
            SMWP,SMSAT,SMREF,SMCNT, LAI, IRRG_LAIMIN, IRRG_LAIMAX, RZDEF,                               &
            IRRG_RATE_SPR, IRRG_RATE_DRP, IRRG_RATE_FRW, IRRG_RATE_PDY,                                 &
            SRATE, DRATE, FRATE) 
       
    else
       
       ! crop calendar based irrigation
       ! ==============================

       call IM%run_model (dofyr,local_hour,                                                             &
            IRRG_IRRIGFRAC_SPR, IRRG_IRRIGFRAC_DRP, IRRG_IRRIGFRAC_FRW,                                 &
            IRRG_CROPIRRIGFRAC,IRRG_DOY_PLANT,IRRG_DOY_HARVEST,IRRG_TYPE ,                              &
            SMWP,SMSAT,SMREF,SMCNT, RZDEF,                                                              & 
            IRRG_RATE_SPR, IRRG_RATE_DRP, IRRG_RATE_FRW, IRRG_RATE_PDY,                                 &
            SRATE, DRATE, FRATE) 

    endif

    ! Scale computed IRRG_RATE_SPR, IRRG_RATE_DRP, IRRG_RATE_FRW, and IRRG_RATE_PDY to the total
    ! irrigated tile fraction before exporting to Catchment[CN].
 
    IRRG_RATE_SPR = IRRG_RATE_SPR * IRRG_IRRIGFRAC 
    IRRG_RATE_DRP = IRRG_RATE_DRP * IRRG_IRRIGFRAC 
    IRRG_RATE_FRW = IRRG_RATE_FRW * IRRG_IRRIGFRAC 
    IRRG_RATE_PDY = IRRG_RATE_PDY * IRRG_PADDYFRAC 
    
    if(associated(IRRG_RATE_TOT))  IRRG_RATE_TOT = IRRG_RATE_SPR + IRRG_RATE_DRP + IRRG_RATE_FRW +IRRG_RATE_PDY
    
    deallocate (local_hour, SMWP, SMSAT, SMREF, SMCNT, RZDEF, IM)

    call MAPL_TimerOff(MAPL,"RUN")
    RETURN_(ESMF_SUCCESS)
    
  end subroutine RUN

end module GEOS_IrrigationGridCompMod


! ========================= EOF ==================================================================================
