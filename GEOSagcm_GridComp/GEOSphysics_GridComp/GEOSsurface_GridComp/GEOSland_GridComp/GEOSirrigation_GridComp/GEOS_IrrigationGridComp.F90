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
  
! IMPORTS: POROS, WPWET, VGWMAX, WCRZ, LAI \\
!
! !USES: 

  use ESMF
  use MAPL
  use IRRIGATION_MODULE ONLY: NUM_CROPS, NUM_SEASONS, IRRIGATION_MODEL
  
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
    type(ESMF_Config)                       :: SCF
    character(len=ESMF_MAXSTR)              :: SURFRC
    integer                                 :: RUN_IRRIG, IRRIG_METHOD, IRRIG_TRIGGER
    
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
! Get runtime switches
! -----------------------------------------------------------

    call MAPL_GetResource (MAPL, SURFRC, label = 'SURFRC:', default = 'GEOS_SurfaceGridComp.rc', RC=STATUS) ; VERIFY_(STATUS)
    SCF = ESMF_ConfigCreate(rc=status) ; VERIFY_(STATUS)
    call ESMF_ConfigLoadFile(SCF,SURFRC,rc=status) ; VERIFY_(STATUS)
    call ESMF_ConfigGetAttribute (SCF, label='RUN_IRRIG:'    , value=RUN_IRRIG   , DEFAULT=0, __RC__ )
    call ESMF_ConfigGetAttribute (SCF, label='IRRIG_METHOD:' , value=IRRIG_METHOD, DEFAULT=0, __RC__ )
    call ESMF_ConfigGetAttribute (SCF, label='IRRIG_TRIGGER:', value=RUN_IRRIG   , DEFAULT=0, __RC__ )    
    call ESMF_ConfigDestroy      (SCF, __RC__)
    
!BOS

! -----------------------------------------------------------
! Internal State 
! -----------------------------------------------------------

    if(RUN_IRRIG == 1) then
       call MAPL_AddInternalSpec(GC                                ,&
            SHORT_NAME = 'IRRIGFRAC'                               ,&
            LONG_NAME  = 'fraction_of_irrigated_cropland'	   ,&
            UNITS      = '1'                                       ,&
            DIMS       = MAPL_DimsTileOnly                         ,&
            VLOCATION  = MAPL_VLocationNone                        ,&
            FRIENDLYTO = trim(COMP_NAME)                           ,&
            RESTART    = MAPL_RestartRequired                      ,&
            RC=STATUS  )
       VERIFY_(STATUS)  
       
       call MAPL_AddInternalSpec(GC                                ,&
            SHORT_NAME = 'PADDYFRAC'                               ,&
            LONG_NAME  = 'fraction_of_paddy_cropland'	           ,&
            UNITS      = '1'                                       ,&
            DIMS       = MAPL_DimsTileOnly                         ,&
            VLOCATION  = MAPL_VLocationNone                        ,&
            FRIENDLYTO = trim(COMP_NAME)                           ,&
            RESTART    = MAPL_RestartRequired                      ,&
            RC=STATUS  )
       VERIFY_(STATUS)
       
       call MAPL_AddInternalSpec(GC                                ,&
            SHORT_NAME = 'CROPIRRIGFRAC'                           ,&
            LONG_NAME  = 'Crop_irrigated_fraction'		   ,&
            UNITS      = '1'                                       ,&
            DIMS       = MAPL_DimsTileOnly                         ,&
            VLOCATION  = MAPL_VLocationNone                        ,&
            FRIENDLYTO = trim(COMP_NAME)                           ,&
            UNGRIDDED_DIMS = (/NUM_CROPS/)                         ,&
            RESTART    = MAPL_RestartRequired                      ,&
            RC=STATUS  )
       VERIFY_(STATUS)
       
       call MAPL_AddInternalSpec(GC                                ,&
            SHORT_NAME = 'IRRIGPLANT'                              ,&
            LONG_NAME  = 'DOY_start_planting'			   ,&
            UNITS      = 'day'                                     ,&
            DIMS       = MAPL_DimsTileOnly                         ,&
            VLOCATION  = MAPL_VLocationNone                        ,&
            FRIENDLYTO = trim(COMP_NAME)                           ,&
            UNGRIDDED_DIMS = (/NUM_SEASONS, NUM_CROPS/)            ,&
            RESTART    = MAPL_RestartRequired                      ,&
            RC=STATUS  )
       VERIFY_(STATUS)
       
       call MAPL_AddInternalSpec(GC                                ,&
            SHORT_NAME = 'IRRIGHARVEST'                            ,&
            LONG_NAME  = 'DOY_end_harvesting'		           ,&
            UNITS      = 'day'                                     ,&
            DIMS       = MAPL_DimsTileOnly                         ,&
            VLOCATION  = MAPL_VLocationNone                        ,&
            FRIENDLYTO = trim(COMP_NAME)                           ,&
            UNGRIDDED_DIMS = (/NUM_SEASONS, NUM_CROPS/)            ,&
            RESTART    = MAPL_RestartRequired                      ,&
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
            RESTART    = MAPL_RestartRequired                      ,&
            RC=STATUS  )
       VERIFY_(STATUS)
       
       call MAPL_AddInternalSpec(GC                                ,&
            SHORT_NAME = 'SPRINKLERFR'                             ,&
            LONG_NAME  = 'fraction_of_sprinkler_irrigation'	   ,&
            UNITS      = '1'                                       ,&
            DIMS       = MAPL_DimsTileOnly                         ,&
            VLOCATION  = MAPL_VLocationNone                        ,&
            FRIENDLYTO = trim(COMP_NAME)                           ,&
            RC=STATUS  )
       VERIFY_(STATUS)
       
       call MAPL_AddInternalSpec(GC                                ,&
            SHORT_NAME = 'DRIPFR'                                  ,&
            LONG_NAME  = 'fraction_of_drip_irrigation'		   ,&
            UNITS      = '1'                                       ,&
            DIMS       = MAPL_DimsTileOnly                         ,&
            VLOCATION  = MAPL_VLocationNone                        ,&
            FRIENDLYTO = trim(COMP_NAME)                           ,&
            RESTART    = MAPL_RestartRequired                      ,&
            RC=STATUS  )
       VERIFY_(STATUS)
       
       call MAPL_AddInternalSpec(GC                                ,&
            SHORT_NAME = 'FLOODFR'                                 ,&
            LONG_NAME  = 'fraction_of_flood_irrigation'	           ,&
            UNITS      = '1'                                       ,&
            DIMS       = MAPL_DimsTileOnly                         ,&
            VLOCATION  = MAPL_VLocationNone                        ,&
            FRIENDLYTO = trim(COMP_NAME)                           ,&
            RESTART    = MAPL_RestartRequired                      ,&
         RC=STATUS  )
       VERIFY_(STATUS)
    
       call MAPL_AddInternalSpec(GC                                ,&
            SHORT_NAME = 'LAIMIN'                                  ,&
            LONG_NAME  = 'Minimum_LAI_irrigated_crops'	           ,&
            UNITS      = '1'                                       ,&
            DIMS       = MAPL_DimsTileOnly                         ,&
            VLOCATION  = MAPL_VLocationNone                        ,&
            FRIENDLYTO = trim(COMP_NAME)                           ,&
            RESTART    = MAPL_RestartRequired                      ,&
            RC=STATUS  )
       VERIFY_(STATUS)
    
       call MAPL_AddInternalSpec(GC                                ,&
            SHORT_NAME = 'LAIMAX'                                  ,&
            LONG_NAME  = 'Maximum_LAI_irrigated_crops'             ,&
            UNITS      = '1'                                       ,&
            DIMS       = MAPL_DimsTileOnly                         ,&
            VLOCATION  = MAPL_VLocationNone                        ,&
            FRIENDLYTO = trim(COMP_NAME)                           ,&
            RESTART    = MAPL_RestartRequired                      ,&
            RC=STATUS  )
       VERIFY_(STATUS)  
    endif

    call MAPL_AddInternalSpec(GC                              ,&
         SHORT_NAME = 'SPRINKLERRATE'                         ,&
         LONG_NAME  = 'sprinkler_irrigation_rate'             ,&
         UNITS      = 'kg m-2 s-1'                            ,&
         DIMS       = MAPL_DimsTileOnly                       ,&
         VLOCATION  = MAPL_VLocationNone                      ,&
         FRIENDLYTO = trim(COMP_NAME)                         ,&
         RESTART    = RESTART                                 ,&
         RC=STATUS  )
    VERIFY_(STATUS)  
    
    call MAPL_AddInternalSpec(GC                              ,&
         SHORT_NAME = 'DRIPRATE'                              ,&
         LONG_NAME  = 'drip_irrigation_rate'	              ,&
         UNITS      = 'kg m-2 s-1'                            ,&
         DIMS       = MAPL_DimsTileOnly                       ,&
         VLOCATION  = MAPL_VLocationNone                      ,&
         FRIENDLYTO = trim(COMP_NAME)                         ,&
         RESTART    = RESTART                                 ,&
         RC=STATUS  )
    VERIFY_(STATUS)  	 
    
    call MAPL_AddInternalSpec(GC                              ,&
         SHORT_NAME = 'FLOODRATE'                             ,&
         LONG_NAME  = 'flood_irrigation_rate'                 ,&
         UNITS      = 'kg m-2 s-1'                            ,&
         DIMS       = MAPL_DimsTileOnly                       ,&
         VLOCATION  = MAPL_VLocationNone                      ,&
         FRIENDLYTO = trim(COMP_NAME)                         ,&
         RESTART    = RESTART                                 ,&
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

    call MAPL_AddImportSpec(GC                                ,&
       SHORT_NAME = 'LAI'                                     ,&
       LONG_NAME  = 'leaf_area_index'                         ,&
       UNITS      = '1'                                       ,&
       DIMS       = MAPL_DimsTileOnly                         ,&
       VLOCATION  = MAPL_VLocationNone                        ,&
       RC=STATUS  )
    VERIFY_(STATUS)

! -----------------------------------------------------------
! Export Variables
! -----------------------------------------------------------


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
    real, dimension(:), pointer :: LAI
  
! Time attributes 

    type(ESMF_Time)             :: CURRENT_TIME
    logical, save               :: firsttime=.true.
    integer                     :: AGCM_YY, AGCM_MM, AGCM_DD, AGCM_MI, AGCM_S, AGCM_HH, dofyr

! Others/ Locals

    type(irrigation_model),save :: IM
    real,pointer,dimension(:)   :: lons
    integer                     :: ntiles, n
    real, dimension (:)         :: local_hour, SMWP, SMSAT, SMREF, SMCNT, LAIFAC, LAITHRES


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

! get pointers to EXPORTS
! -----------------------

    call MAPL_GetPointer(EXPORT, SPRINKLERRATE, 'SPRINKLERRATE',ALLOC=.true., RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, DRIPRATE,      'DRIPRATE',     ALLOC=.true., RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, FLOODRATE,     'FLOODRATE',    ALLOC=.true., RC=STATUS) ; VERIFY_(STATUS)

    if (RUN_IRRI G ==0) then
       if(associated(SPRINKLERRATE)) SPRINKLERRATE = 0.
       if(associated(DRIPRATE))      DRIPRATE      = 0.
       if(associated(FLOODRATE))     FLOODRATE     = 0.
       call MAPL_TimerOff(MAPL,"TOTAL")
       RETURN_(ESMF_SUCCESS)
    endif

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

! get pointers to IMPORT variables
! --------------------------------

    call MAPL_GetPointer(IMPORT, POROS , 'POROS',  RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, WPWET , 'WPWET',  RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, VGWMAX, 'VGWMAX', RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, WCRZ  , 'WCRZ',   RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, LAI   , 'LAI',    RC=STATUS) ; VERIFY_(STATUS)
    
    if (firsttime) then
       
       if(associated(SPRINKLERRATE)) SPRINKLERRATE = 0.
       if(associated(DRIPRATE))      DRIPRATE      = 0.
       if(associated(FLOODRATE))     FLOODRATE     = 0.
       call IM%init_model (SURFRC)
       firsttime = .false.
       call MAPL_TimerOff(MAPL,"TOTAL")
       RETURN_(ESMF_SUCCESS)
       
    endif

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

    call MAPL_Get (MAPL, TILELONS = LONS,                  &
         INTERNAL_ESMF_STATE = INTERNAL, RC=STATUS )
    VERIFY_(STATUS)
    

! call the irrigation model 
! -------------------------

    NTILES = SIZE (LONS)
    allocate (local_hour (1:NTILES))
    allocate (SMWP       (1:NTILES))
    allocate (SMSAT      (1:NTILES))
    allocate (SMREF      (1:NTILES))
    allocate (SMCNT      (1:NTILES))

    ! soil moisture state
    SMWP  = VGWMAX * WPWET           ! RZ soil moisture content at wilting point [mm]
    SMSAT = VGWMAX                   ! RZ soil moisture at saturation  [mm]
    SMCNT = (VGWMAX/POROS) * WCRZ    ! actual RZ soil moisture content [mm]
       
    DO N = 1, NTILES
       ! local times
       local_hour(n) = AGCM_HH + AGCM_MI / 60. + AGCM_S / 3600. + 12.* (lons(n)/MAPL_PI)
       IF (local_hour(n) >= 24.) local_hour(n) = local_hour(n) - 24.
       IF (local_hour(n) <   0.) local_hour(n) = local_hour(n) + 24.

       ! reference soil moisture content is at lower tercile of RZ soil moisture
       SMREF (n) = VGWMAX (n) * poros(n) * (wpwet (n) + (1. - wpwet (n))/3.) 
       
    END DO

    if (IRRIG_TRIGGER == 0) then

       ! LAI based trigger: scale soil moisture to LAI seasonal cycle
       ! ============================================================
       
       allocate (LAIFAC   (1:NTILES))
       allocate (LAITHRES (1:NTILES))
       
       DO N = 1, NTILES
          LAITHRES (N)   = LAIMIN (N) + 0.60 * (LAIMAX (N) - LAIMIN (N))
          IF(LAIMAX (N) > LAIMIN (N)) THEN
             LAIFAC (N) = (LAI(N) - LAIMIN (N)) / (LAIMAX(N) - LAIMIN(N))
          ELSE
             LAIFAC (N) = 0.
          ENDIF
       END DO
      
       SMWP  = LAIFAC * SMWP  
       SMSAT = LAIFAC * SMSAT
       SMREF = LAIFAC * SMREF
       SMCNT = LAIFAC * SMCNT
             
       call IM%run_model(IRRIG_METHOD, local_hour,                 &
            IRRIGFRAC, PADDYFRAC, SPRINKLERFR, DRIPFR, FLOODFR,    &           
            SMWP, SMSAT, SMREF, SMCNT, LAI, LAITHRES,              &
            SPRINKLERRATE, DRIPRATE, FLOODRATE) 

       deallocate (LAIFAC, LAITHRES)
       
    else
       
       ! crop calendar based irrigation
       ! ==============================
       call IM%run_model

    endif
    
    deallocate (local_hour, SMWP, SMSAT, SMREF, SMCNT)
    
    call MAPL_TimerOff(MAPL,"TOTAL")
    RETURN_(ESMF_SUCCESS)
    
  end subroutine RUN

end module GEOS_IrrigationGridCompMod
