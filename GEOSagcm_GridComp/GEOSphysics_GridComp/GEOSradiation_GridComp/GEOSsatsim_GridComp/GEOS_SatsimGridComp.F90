
! VERIFY_ and RETURN_ macros for error handling.
! test may 7
#include "MAPL_Generic.h"
#undef RTTOV

module GEOS_SatsimGridCompMod

!BOP
! !MODULE: GEOS_Satsim  -- A Module to drive satellite simulators using grid mean cloud parameters

! !DESCRIPTION:
! 
!   {\tt GEOS\_MoistGridCompMod} implements mosit processes in GEOS-5. These
!   include all processes that involve phase changes in the atmosphere, such
!   as large-scale condensation, convective clouds, and all rain and cloud
!   formation. It's state consists of water vapor, various types of condensate,
!   and fractions of various cloud types.
!
!
! !USES:
#define USE_MAPL_UNDEF

  use ESMF
  use MAPL
  use GEOS_UtilsMod

  use gettau

  use MOD_COSP_TYPES
  use MOD_COSP, only: COSP
!  use MOD_COSP_Modis_Simulator, only: cosp_modis
  use MOD_COSP_Modis_Simulator

  implicit none
  private

! !PUBLIC MEMBER FUNCTIONS:

  public SetServices

! Private State
  TYPE SatSim_State
     PRIVATE

     integer    :: nmask_vars ! number of masked variables
     character(len=ESMF_MAXSTR), pointer :: export_name(:) => null()
     character(len=ESMF_MAXSTR), pointer :: mask_name(:)   => null()
     character(len=ESMF_MAXSTR), pointer :: newvar_name(:) => null()
     logical, pointer :: newvar(:) => null()

  END TYPE SatSim_State

! Hook for ESMF
! -------------

  TYPE SatSim_Wrap
     TYPE(SatSim_State), pointer :: PTR => null()
  END TYPE SatSim_Wrap

!EOP

contains

!BOP
! !IROUTINE: SetServices -- Sets ESMF services for this component

! !INTERFACE:
    subroutine SetServices ( GC, RC )

! !ARGUMENTS:
    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    integer, optional                  :: RC  ! return code
    
! !DESCRIPTION:  {\tt GEOS\_MoistGridCompMod} uses the default Initialize and Finalize 
!                services, but registers its own Run method.

!EOP

!=============================================================================
!
! ErrLog Variables

    character(len=ESMF_MAXSTR)              :: IAm
    integer                                 :: STATUS
    character(len=ESMF_MAXSTR)              :: COMP_NAME

! Local derived type aliases

    type (MAPL_MetaComp    ), pointer   :: STATE 
    type (ESMF_Config          )            :: CF

    integer      :: DIMS(3)

!   Local derived type aliases
!   --------------------------
    type (SatSim_State), pointer  :: self => null()   ! internal, that is
    type (SatSim_wrap)            :: wrap

    integer :: nLines,nCols,m,vindex
    character(len=ESMF_MAXSTR)    :: tmpname
    character(len=ESMF_MAXSTR)    :: long_name
    character(len=ESMF_MAXSTR)    :: units
    integer                       :: mapl_dims,vlocation
    integer, pointer              :: ungridded_dims(:) => null()
    character(len=ESMF_MAXSTR)    :: exportName(1)
    type(MAPL_VarSpec), pointer :: ExportSpec(:) => null()
    logical :: found
   
!=============================================================================

! Begin...

! Get my name and set-up traceback handle
! ---------------------------------------

    Iam = 'SetServices'
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // Iam

!   Wrap internal state for storing in GC; rename legacyState
!   -------------------------------------
    allocate ( self, stat=STATUS )
    VERIFY_(STATUS)
    wrap%ptr => self

! Set the Run entry point
! -----------------------

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,  Run,  &
                                      RC=STATUS)
    VERIFY_(STATUS)
    

! Get the configuration from the component
!-----------------------------------------

    call ESMF_GridCompGet( GC, CONFIG = CF, RC=STATUS )
    VERIFY_(STATUS)

!BOS
! !IMPORT STATE:

    call MAPL_AddImportSpec(GC,                           &
        SHORT_NAME ='T',                                          &
        LONG_NAME  ='air_temperature',                          &
        UNITS      ='K',                                          &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
     VERIFY_(STATUS)                                                                          

     call MAPL_AddImportSpec(GC,                           &
        SHORT_NAME ='QV',                                          &
        LONG_NAME  ='specific_humidity',                          &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
     VERIFY_(STATUS)                                                                          

     call MAPL_AddImportSpec(GC,                           &
        SHORT_NAME ='FCLD',                                       &
        LONG_NAME  ='cloud_area_fraction',             &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                               &
         SHORT_NAME='RL',                                          & 
         LONG_NAME ='liquid_cloud_particle_effective_radius',      &
         UNITS     ='m',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                               &
         SHORT_NAME='RI',                                          & 
         LONG_NAME ='ice_phase_cloud_particle_effective_radius',   &
         UNITS     ='m',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                               &
         SHORT_NAME='RS',                                          & 
         LONG_NAME ='snow_particle_effective_radius',      &
         UNITS     ='m',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                               &
         SHORT_NAME='RR',                                          & 
         LONG_NAME ='rain_cloud_particle_effective_radius',   &
         UNITS     ='m',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                           &
        SHORT_NAME ='PLE',                                       &
        LONG_NAME  ='Edge_pressures',             &
        UNITS      ='Pa',                                          &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationEdge,             RC=STATUS  )  
     VERIFY_(STATUS)                                                                          

     call MAPL_AddImportSpec(GC,                           &
        SHORT_NAME ='ZLE',                                       &
        LONG_NAME  ='Edge_heights',             &
        UNITS      ='m',                                          &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationEdge,             RC=STATUS  )  
     VERIFY_(STATUS)                                                                          

    call MAPL_AddImportSpec(GC,                           &
        SHORT_NAME ='QLLS',                                       &
        LONG_NAME  ='mass_fraction_of_large_scale_cloud_liquid_water', &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                           &
        SHORT_NAME ='QLCN',                                       &
        LONG_NAME  ='mass_fraction_of_convective_cloud_liquid_water', &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                           &
        SHORT_NAME ='QILS',                                       &
        LONG_NAME  ='mass_fraction_of_large_scale_cloud_ice_water', &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                           &
        SHORT_NAME ='QICN',                                       &
        LONG_NAME  ='mass_fraction_of_convective_cloud_ice_water', &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME ='QRTOT',                                      &
        LONG_NAME  ='mass_fraction_of_falling_rain',              &
        UNITS      ='kg kg-1',                                    &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME ='QSTOT',                                      &
        LONG_NAME  ='mass_fraction_of_falling_snow',              &
        UNITS      ='kg kg-1',                                    &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME ='TS',                                         &
        LONG_NAME  ='skin_temperature'                          , &
        UNITS      ='K',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                           &
        SHORT_NAME ='MCOSZ',                                      &
        LONG_NAME  = 'mean_cosine_of_the_solar_zenith_angle',     &
        UNITS      = '1',                                         &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                           &
        SHORT_NAME ='FRLAND',                                     &
        LONG_NAME  = 'fraction_of_land',                          &
        UNITS      = '1',                                         &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                           &
        SHORT_NAME ='FROCEAN',                                    &
        LONG_NAME  = 'fraction_of_ocean',                         &
        UNITS      = '1',                                         &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)




! !EXPORT STATE:


!! Fractions in 49 ISCCP Cloud categories
!!    9 super categories
!///////////////////////////////////////////
!     0.   |    .    |    .    |    .       c
! 1        |    .    |    .    |    .       l
!   180.  .|...VII...|..VIII...|....IX....  o 
! 2        |    .    |    .    |    .       u
!   310.  .|.........|.........|..........  d
! 3        |    .    |    .    |    .       
!   440. ---------------------------------  t
! 4        |    .    |    .    |    .       o
!   560.  .|....IV...|....V....|....VI....  p
! 5        |    .    |    .    |    .       
!   680. ---------------------------------  p
! 6        | oa . ob |    .    |    .       r
!   800.  .|....I....|...II....|...III....  e
! 7        | ua . ub |    .    |    .       s
!   sfc    |    .    |    .    |    .       s
!      ====================================
!      0.  x   1.3  3.6  9.4  23.  60.   -> 
!        1   2     3    4    5    6    7
!            optical depth increasing ->
!
!   9 Supercategories:
!     I=cumulus(Cu), II=stratocumulus(StCu), III=stratus(St)
!    IV=altocumulus(ACu), V=altostratus(ASt), VI=nimbostratus(NSt)
!    VII=cirrus(Ci), VIII=cirrostratus(CiSt), IX=Deep convection(Cb)
!  
!   Supercategories further subvided into high/over(O),{middle(M)}, or 
!   low/under(U) and thin (A) and thick (B)
!
!   In addition 7 categories of subvisual clouds (Sub) defined by 
!   P_top with Sub_1 above 180 hPa ... etc
!
!   In ISCCP simulator 7x7 frequency for these types is arranged thus: 
!
!         fq_isccp(  itau(1-7),ipres(1-7) )
!   
!   So frequencies for 7 subvisible cloud types are:
!
!          fq_isccp(  1 , 1 ) = subvis above 180.
!          fq_isccp(  1 , 2 ) = subvis between 180. and 310.
!               etc. ..
!


!
! 4 Cumulus (CU) subcategories
!    fq_isccp(2:3,6:7)  
!-------------------------------------------------
!    fq_isccp(2,7)
     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='ISCCP_CU_UA',                                   &
        LONG_NAME  ='isccp_fraction_of_thin_lower_cumulus',         &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

!    fq_isccp(3,7)
     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='ISCCP_CU_UB',                                   &
        LONG_NAME  ='isccp_fraction_of_thick_lower_cumulus',         &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

!    fq_isccp(2,6)
     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='ISCCP_CU_OA',                                   &
        LONG_NAME  ='isccp_fraction_of_thin_higher_cumulus',         &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

!    fq_isccp(3,6)
     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='ISCCP_CU_OB',                                   &
        LONG_NAME  ='isccp_fraction_of_thick_higher_cumulus',         &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

! 4 Stratocumulus STCU subcategories
!    fq_isccp(4:5,6:7)  
!-------------------------------------------------
!    fq_isccp(4,7)
     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='ISCCP_STCU_UA',                                 &
        LONG_NAME  ='isccp_fraction_of_thin_lower_stratocumulus',       &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

!    fq_isccp(5,7)
     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='ISCCP_STCU_UB',                                 &
        LONG_NAME  ='isccp_fraction_of_thick_lower_stratocumulus',      &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

!    fq_isccp(4,6)
     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='ISCCP_STCU_OA',                                 &
        LONG_NAME  ='isccp_fraction_of_thin_higher_stratocumulus',      &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

!    fq_isccp(5,6)
     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='ISCCP_STCU_OB',                                 &
        LONG_NAME  ='isccp_fraction_of_thick_higher_stratocumulus',     &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)


! 4 Stratus ST subcategories
!    fq_isccp(6:7,6:7)  
!-------------------------------------------------
!    fq_isccp(6,7)
     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='ISCCP_ST_UA',                                 &
        LONG_NAME  ='isccp_fraction_of_thin_lower_stratus',       &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

!    fq_isccp(7,7)
     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='ISCCP_ST_UB',                                 &
        LONG_NAME  ='isccp_fraction_of_thick_lower_stratus',      &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

!    fq_isccp(6,6)
     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='ISCCP_ST_OA',                                 &
        LONG_NAME  ='isccp_fraction_of_thin_higher_stratus',      &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

!    fq_isccp(7,6)
     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='ISCCP_ST_OB',                                 &
        LONG_NAME  ='isccp_fraction_of_thick_higher_stratus',     &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

! 4 Altocumulus ACU subcategories
!    fq_isccp(2:3,4:5)  
!-------------------------------------------------
!    fq_isccp(2,5)
     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='ISCCP_ACU_UA',                                 &
        LONG_NAME  ='isccp_fraction_of_thin_lower_altocumulus',       &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

!    fq_isccp(3,5)
     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='ISCCP_ACU_UB',                                 &
        LONG_NAME  ='isccp_fraction_of_thick_lower_altocumulus',      &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

!    fq_isccp(2,4)
     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='ISCCP_ACU_OA',                                 &
        LONG_NAME  ='isccp_fraction_of_thin_higher_altocumulus',      &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

!    fq_isccp(3,4)
     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='ISCCP_ACU_OB',                                 &
        LONG_NAME  ='isccp_fraction_of_thick_higher_altocumulus',     &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)


! 4 Altostratus AST subcategories
!    fq_isccp(4:5,4:5)  
!-------------------------------------------------
!    fq_isccp(4,5)
     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='ISCCP_AST_UA',                                 &
        LONG_NAME  ='isccp_fraction_of_thin_lower_altostratus',       &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

!    fq_isccp(5,5)
     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='ISCCP_AST_UB',                                 &
        LONG_NAME  ='isccp_fraction_of_thick_lower_altostratus',      &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

!    fq_isccp(4,4)
     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='ISCCP_AST_OA',                                 &
        LONG_NAME  ='isccp_fraction_of_thin_higher_altostratus',      &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

!    fq_isccp(5,4)
     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='ISCCP_AST_OB',                                 &
        LONG_NAME  ='isccp_fraction_of_thick_higher_altostratus',     &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)



! 4 Nimbostratus NST subcategories
!    fq_isccp(6:7,4:5)  
!-------------------------------------------------
!    fq_isccp(6,5)
     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='ISCCP_NST_UA',                                 &
        LONG_NAME  ='isccp_fraction_of_thin_lower_nimbostratus',       &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

!    fq_isccp(7,5)
     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='ISCCP_NST_UB',                                 &
        LONG_NAME  ='isccp_fraction_of_thick_lower_nimbostratus',      &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

!    fq_isccp(6,4)
     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='ISCCP_NST_OA',                                 &
        LONG_NAME  ='isccp_fraction_of_thin_higher_nimbostratus',      &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

!    fq_isccp(7,4)
     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='ISCCP_NST_OB',                                 &
        LONG_NAME  ='isccp_fraction_of_thick_higher_nimbostratus',     &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)



! 6 Cirrus Ci subcategories
!    fq_isccp(2:3,1:3)  
!-------------------------------------------------
!    fq_isccp(2,3)
     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='ISCCP_CI_UA',                                 &
        LONG_NAME  ='isccp_fraction_of_thin_lower_cirrus',       &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

!    fq_isccp(3,3)
     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='ISCCP_CI_UB',                                 &
        LONG_NAME  ='isccp_fraction_of_thick_lower_cirrus',      &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

!    fq_isccp(2,2)
     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='ISCCP_CI_MA',                                 &
        LONG_NAME  ='isccp_fraction_of_thin_middle_cirrus',       &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

!    fq_isccp(3,2)
     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='ISCCP_CI_MB',                                 &
        LONG_NAME  ='isccp_fraction_of_thick_middle_cirrus',      &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

!    fq_isccp(2,1)
     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='ISCCP_CI_OA',                                 &
        LONG_NAME  ='isccp_fraction_of_thin_higher_cirrus',      &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

!    fq_isccp(3,1)
     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='ISCCP_CI_OB',                                 &
        LONG_NAME  ='isccp_fraction_of_thick_higher_cirrus',     &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)


! 6 Cirrostratus CIST subcategories
!    fq_isccp(4:5,1:3)  
!-------------------------------------------------
!    fq_isccp(4,3)
     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='ISCCP_CIST_UA',                                 &
        LONG_NAME  ='isccp_fraction_of_thin_lower_cirrostratus',       &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

!    fq_isccp(5,3)
     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='ISCCP_CIST_UB',                                 &
        LONG_NAME  ='isccp_fraction_of_thick_lower_cirrostratus',      &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

!    fq_isccp(4,2)
     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='ISCCP_CIST_MA',                                 &
        LONG_NAME  ='isccp_fraction_of_thin_middle_cirrostratus',       &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

!    fq_isccp(5,2)
     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='ISCCP_CIST_MB',                                 &
        LONG_NAME  ='isccp_fraction_of_thick_middle_cirrostratus',      &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

!    fq_isccp(4,1)
     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='ISCCP_CIST_OA',                                 &
        LONG_NAME  ='isccp_fraction_of_thin_higher_cirrostratus',      &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

!    fq_isccp(5,1)
     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='ISCCP_CIST_OB',                                 &
        LONG_NAME  ='isccp_fraction_of_thick_higher_cirrostratus',     &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)


! 6 Cumulonimbus CB subcategories
!    fq_isccp(6:7,1:3)  
!-------------------------------------------------
!    fq_isccp(6,3)  
     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='ISCCP_CB_UA',                                 &
        LONG_NAME  ='isccp_fraction_of_thin_lower_cumulonimbus',       &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

!    fq_isccp(7,3)  
     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='ISCCP_CB_UB',                                 &
        LONG_NAME  ='isccp_fraction_of_thick_lower_cumulonimbus',      &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

!    fq_isccp(6,2)  
     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='ISCCP_CB_MA',                                 &
        LONG_NAME  ='isccp_fraction_of_thin_middle_cumulonimbus',       &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

!    fq_isccp(7,2)  
     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='ISCCP_CB_MB',                                 &
        LONG_NAME  ='isccp_fraction_of_thick_middle_cumulonimbus',      &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

!    fq_isccp(6,1)  
     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='ISCCP_CB_OA',                                 &
        LONG_NAME  ='isccp_fraction_of_thin_higher_cumulonimbus',      &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

!    fq_isccp(7,1)  
     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='ISCCP_CB_OB',                                 &
        LONG_NAME  ='isccp_fraction_of_thick_higher_cumulonimbus',     &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)


! 7 subvisible/subdetection (SUBV) subcategories
!    fq_isccp(1,1:7)  
!-------------------------------------------------
!    fq_isccp(1,1)  
     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='ISCCP_SUBV1'  ,                                 &
        LONG_NAME  ='isccp_fraction_of_subvisible_cloud_0_180_hPa',     &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

!    fq_isccp(1,2)  
     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='ISCCP_SUBV2'  ,                                 &
        LONG_NAME  ='isccp_fraction_of_subvisible_cloud_180_310_hPa',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

!    fq_isccp(1,3)  
     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='ISCCP_SUBV3'  ,                                 &
        LONG_NAME  ='isccp_fraction_of_subvisible_cloud_310_440_hPa',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

!    fq_isccp(1,4)  
     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='ISCCP_SUBV4'  ,                                 &
        LONG_NAME  ='isccp_fraction_of_subvisible_cloud_440_560_hPa',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

!    fq_isccp(1,5)  
     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='ISCCP_SUBV5'  ,                                 &
        LONG_NAME  ='isccp_fraction_of_subvisible_cloud_560_680_hPa',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

!    fq_isccp(1,6)  
     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='ISCCP_SUBV6'  ,                                 &
        LONG_NAME  ='isccp_fraction_of_subvisible_cloud_680_800_hPa',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

!    fq_isccp(1,7)  
     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='ISCCP_SUBV7'  ,                                 &
        LONG_NAME  ='isccp_fraction_of_subvisible_cloud_800_SFC_hPa',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

! 7 exports fields corresponding to pressure bands, with 7 thickess
! classes in each -- this is to partially comply with CFMIP data specs
!    fq_isccp(:,1:7)  
!-------------------------------------------------
!    fq_isccp(:,1)  
     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='CLISCCP1'  ,                                   &
        LONG_NAME  ='isccp_cloud_area_fraction_0_180_hPa',                     &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        UNGRIDDED_DIMS = (/ 7 /),                                 &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

!    fq_isccp(:,2)  
     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='CLISCCP2'  ,                                 &
        LONG_NAME  ='isccp_cloud_area_fraction_180_310_hPa',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        UNGRIDDED_DIMS = (/ 7 /),                                 &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

!    fq_isccp(:,3)  
     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='CLISCCP3'  ,                                 &
        LONG_NAME  ='isccp_cloud_area_fraction_310_440_hPa',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        UNGRIDDED_DIMS = (/ 7 /),                                 &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

!    fq_isccp(:,4)  
     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='CLISCCP4'  ,                                 &
        LONG_NAME  ='isccp_cloud_area_fraction_440_560_hPa',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        UNGRIDDED_DIMS = (/ 7 /),                                 &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

!    fq_isccp(:,5)  
     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='CLISCCP5'  ,                                 &
        LONG_NAME  ='isccp_cloud_area_fraction_560_680_hPa',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        UNGRIDDED_DIMS = (/ 7 /),                                 &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

!    fq_isccp(:,6)  
     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='CLISCCP6'  ,                                 &
        LONG_NAME  ='isccp_cloud_area_fraction_680_800_hPa',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        UNGRIDDED_DIMS = (/ 7 /),                                 &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

!    fq_isccp(:,7)  
     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='CLISCCP7'  ,                                 &
        LONG_NAME  ='isccp_cloud_area_fraction_800_SFC_hPa',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        UNGRIDDED_DIMS = (/ 7 /),                                 &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

! The next 7 exports (ISCCP[1-7]) are deprecated
!    fq_isccp(:,1)  
     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='ISCCP1'  ,                                   &
        LONG_NAME  ='isccp_cloud_area_fraction_0_180_hPa',                     &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        UNGRIDDED_DIMS = (/ 7 /),                                 &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

!    fq_isccp(:,2)  
     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='ISCCP2'  ,                                 &
        LONG_NAME  ='isccp_cloud_area_fraction_180_310_hPa',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        UNGRIDDED_DIMS = (/ 7 /),                                 &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

!    fq_isccp(:,3)  
     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='ISCCP3'  ,                                 &
        LONG_NAME  ='isccp_cloud_area_fraction_310_440_hPa',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        UNGRIDDED_DIMS = (/ 7 /),                                 &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

!    fq_isccp(:,4)  
     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='ISCCP4'  ,                                 &
        LONG_NAME  ='isccp_cloud_area_fraction_440_560_hPa',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        UNGRIDDED_DIMS = (/ 7 /),                                 &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

!    fq_isccp(:,5)  
     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='ISCCP5'  ,                                 &
        LONG_NAME  ='isccp_cloud_area_fraction_560_680_hPa',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        UNGRIDDED_DIMS = (/ 7 /),                                 &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

!    fq_isccp(:,6)  
     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='ISCCP6'  ,                                 &
        LONG_NAME  ='isccp_cloud_area_fraction_680_800_hPa',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        UNGRIDDED_DIMS = (/ 7 /),                                 &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

!    fq_isccp(:,7)  
     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='ISCCP7'  ,                                 &
        LONG_NAME  ='isccp_cloud_area_fraction_800_SFC_hPa',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        UNGRIDDED_DIMS = (/ 7 /),                                 &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)


! Other ISCCP ouputs
!---------------------------------------------------


    call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='CLTISCCP'  ,                                 &
        LONG_NAME  ='isccp_cloud_area_fraction',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

! same as CLTISCCP, deprecated in favor of CFMIP short name
    call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='TCLISCCP'  ,                                 &
        LONG_NAME  ='isccp_cloud_area_fraction',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='PCTISCCP'  ,                                 &
        LONG_NAME  ='isccp_air_pressure_at_cloud_top',   &
        UNITS      ='Pa',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

! same as PCTISCCP, deprecated in favor of CFMIP short name
    call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='CTPISCCP'  ,                                 &
        LONG_NAME  ='isccp_air_pressure_at_cloud_top',   &
        UNITS      ='Pa',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='ALBISCCP'  ,                                 &
        LONG_NAME  ='isccp_cloud_albedo',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='TBISCCP'  ,                                 &
        LONG_NAME  ='isccp_mean_all_sky_10.5_micron_brightness_temp',   &
        UNITS      ='K',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

! cloud fraction diagnostics from SCOPS 

    call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='SGFCLD'  ,                                 &
        LONG_NAME  ='summed_subgrid_cloud_fraction_from_scops',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

! lidar_simulator output - not used by CFMIP
! MAPL_VLocationCenter is a guess - afe
! ---------------------------------------------
    
     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='LIDARPMOL',                                       &
        LONG_NAME  ='calipso_molecular_attenuated_backscatter_signal_power',   &
        UNITS      ='m-1 sr-1',                                          &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)
    
     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='LIDARPTOT',                                       &
        LONG_NAME  ='calipso_total_attenuated_backscatter_signal_power',   &
        UNITS      ='m-1 sr-1',                                          &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='LIDARTAUTOT',                                     &
        LONG_NAME  ='calipso_optical_thickess_integrated_from_top_to_level_z',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='RADARZETOT',                                     &
        LONG_NAME  ='cloudsat_total_reflectivity',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='CLCALIPSO',                                       &
        LONG_NAME  ='calipso_total_cloud_fraction',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='CLLCALIPSO'  ,                                 &
        LONG_NAME  ='calipso_low_level_cloud_fraction',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='CLMCALIPSO'  ,                                 &
        LONG_NAME  ='calipso_mid_level_cloud_fraction',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='CLHCALIPSO'  ,                                 &
        LONG_NAME  ='calipso_high_level_cloud_fraction',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)
 
     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='CLTCALIPSO'  ,                                 &
        LONG_NAME  ='calipso_total_cloud_fraction',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

      call MAPL_AddExportSpec(GC, &
        SHORT_NAME ='PARASOLREFL0'  , &
        LONG_NAME  ='parasol_reflectance', &
        UNITS      ='1', &
        DIMS       = MAPL_DimsHorzOnly, & 
        UNGRIDDED_DIMS = (/ PARASOL_NREFL /), &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )

     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='PARASOLREFL1'  ,                                 &
        LONG_NAME  ='parasol_reflectance_1',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

      call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='PARASOLREFL2'  ,                                 &
        LONG_NAME  ='parasol_reflectance_2',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)
 
     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='PARASOLREFL3'  ,                                 &
        LONG_NAME  ='parasol_reflectance_3',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)
 
     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='PARASOLREFL4'  ,                                 &
        LONG_NAME  ='parasol_reflectance_4',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)
 
     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='PARASOLREFL5'  ,                                 &
        LONG_NAME  ='parasol_reflectance_5',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)
 

      call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='RADARLTCC'  ,                                 &
        LONG_NAME  ='cloudsat_calipso_total_cloud_amount',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)
 
     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='CLCALIPSO2'  ,                                 &
        LONG_NAME  ='calipso_no_cloudsat_cloud_fraction',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

      call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='CFADLIDARSR532_01'  ,                                 &
        LONG_NAME  ='calipso_scattering_ratio_cfad_01',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)
 
       call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='CFADLIDARSR532_02'  ,                                 &
        LONG_NAME  ='calipso_scattering_ratio_cfad_02',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)
 
       call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='CFADLIDARSR532_03'  ,                                 &
        LONG_NAME  ='calipso_scattering_ratio_cfad_03',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)
 
       call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='CFADLIDARSR532_04'  ,                                 &
        LONG_NAME  ='calipso_scattering_ratio_cfad_04',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)
 
       call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='CFADLIDARSR532_05'  ,                                 &
        LONG_NAME  ='calipso_scattering_ratio_cfad_05',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)
 
       call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='CFADLIDARSR532_06'  ,                                 &
        LONG_NAME  ='calipso_scattering_ratio_cfad_06',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)
 
       call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='CFADLIDARSR532_07'  ,                                 &
        LONG_NAME  ='calipso_scattering_ratio_cfad_07',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)
 
      call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='CFADLIDARSR532_08'  ,                                 &
        LONG_NAME  ='calipso_scattering_ratio_cfad_08',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)
 
       call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='CFADLIDARSR532_09'  ,                                 &
        LONG_NAME  ='calipso_scattering_ratio_cfad_09',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)
 
       call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='CFADLIDARSR532_10'  ,                                 &
        LONG_NAME  ='calipso_scattering_ratio_cfad_10',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)
 
       call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='CFADLIDARSR532_11'  ,                                 &
        LONG_NAME  ='calipso_scattering_ratio_cfad_11',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)
 
       call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='CFADLIDARSR532_12'  ,                                 &
        LONG_NAME  ='calipso_scattering_ratio_cfad_12',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)
 
       call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='CFADLIDARSR532_13'  ,                                 &
        LONG_NAME  ='calipso_scattering_ratio_cfad_13',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)
 
       call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='CFADLIDARSR532_14'  ,                                 &
        LONG_NAME  ='calipso_scattering_ratio_cfad_14',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)
 
       call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='CFADLIDARSR532_15'  ,                                 &
        LONG_NAME  ='calipso_scattering_ratio_cfad_15',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)
 
 
! Radar simulator exports
!


      call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='CLOUDSATCFAD01'  ,                                 &
        LONG_NAME  ='cloudsat_radar_reflectivity_cfad',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)
 
       call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='CLOUDSATCFAD02'  ,                                 &
        LONG_NAME  ='cloudsat_radar_reflectivity_cfad',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)
 
       call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='CLOUDSATCFAD03'  ,                                 &
        LONG_NAME  ='cloudsat_radar_reflectivity_cfad',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)
 
       call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='CLOUDSATCFAD04'  ,                                 &
        LONG_NAME  ='cloudsat_radar_reflectivity_cfad',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)
 
       call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='CLOUDSATCFAD05'  ,                                 &
        LONG_NAME  ='cloudsat_radar_reflectivity_cfad',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)
 
       call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='CLOUDSATCFAD06'  ,                                 &
        LONG_NAME  ='cloudsat_radar_reflectivity_cfad',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)
 
       call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='CLOUDSATCFAD07'  ,                                 &
        LONG_NAME  ='cloudsat_radar_reflectivity_cfad',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)
 
      call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='CLOUDSATCFAD08'  ,                                 &
        LONG_NAME  ='cloudsat_radar_reflectivity_cfad',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)
 
       call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='CLOUDSATCFAD09'  ,                                 &
        LONG_NAME  ='cloudsat_radar_reflectivity_cfad',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)
 
       call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='CLOUDSATCFAD10'  ,                                 &
        LONG_NAME  ='cloudsat_radar_reflectivity_cfad',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)
 
       call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='CLOUDSATCFAD11'  ,                                 &
        LONG_NAME  ='cloudsat_radar_reflectivity_cfad',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)
 
       call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='CLOUDSATCFAD12'  ,                                 &
        LONG_NAME  ='cloudsat_radar_reflectivity_cfad',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)
 
       call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='CLOUDSATCFAD13'  ,                                 &
        LONG_NAME  ='cloudsat_radar_reflectivity_cfad',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)
 
       call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='CLOUDSATCFAD14'  ,                                 &
        LONG_NAME  ='cloudsat_radar_reflectivity_cfad',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)
 
       call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='CLOUDSATCFAD15'  ,                                 &
        LONG_NAME  ='cloudsat_radar_reflectivity_cfad',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

!
! MODIS simulator output
!-------------------------------------------------------------------- 
      call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='MDSCLDFRCTTL'  ,                                 &
        LONG_NAME  ='modis_cloud_fraction_total_mean',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

      call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='MDSCLDFRCWTR'  ,                                 &
        LONG_NAME  ='modis_cloud_fraction_water_mean',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)
 
! the following is deprecated, use MDSCLDFRCWTR
     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='MDSCLDFRCH2O'  ,                                 &
        LONG_NAME  ='modis_cloud_fraction_water_mean',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)
 
         call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='MDSCLDFRCICE'  ,                                 &
        LONG_NAME  ='modis_cloud_fraction_ice_mean',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='MDSCLDFRCHI'  ,                                 &
        LONG_NAME  ='modis_cloud_fraction_high_mean',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='MDSCLDFRCMID'  ,                                 &
        LONG_NAME  ='modis_cloud_fraction_mid_mean',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='MDSCLDFRCLO'  ,                                 &
        LONG_NAME  ='modis_cloud_fraction_low_mean',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='MDSOPTHCKTTL'  ,                                 &
        LONG_NAME  ='modis_optical_thickness_total_mean',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='MDSOPTHCKWTR'  ,                                 &
        LONG_NAME  ='modis_optical_thickness_water_mean',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

! the following is deprecated, use MDSOPTHCKWTR
    call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='MDSOPTHCKH2O'  ,                                 &
        LONG_NAME  ='modis_optical_thickness_water_mean',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='MDSOPTHCKICE'  ,                                 &
        LONG_NAME  ='modis_optical_thickness_ice_mean',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='MDSOPTHCKTTLLG'  ,                                 &
        LONG_NAME  ='modis_optical_thickness_total_logmean',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='MDSOPTHCKWTRLG'  ,                                 &
        LONG_NAME  ='modis_optical_thickness_water_logmean',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

! the following is deprecated, use MDSOPTHCKWTRLG
     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='MDSOPTHCKH2OLG'  ,                                 &
        LONG_NAME  ='modis_optical_thickness_water_logmean',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='MDSOPTHCKICELG'  ,                                 &
        LONG_NAME  ='modis_optical_thickness_ice_logmean',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='MDSCLDSZWTR'  ,                                 &
        LONG_NAME  ='modis_cloud_particle_size_water_mean',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

! the following is deprecated, use MDSCLDSZWTR
     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='MDSCLDSZH20'  ,                                 &
        LONG_NAME  ='modis_cloud_particle_size_water_mean',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='MDSCLDSZICE'  ,                                 &
        LONG_NAME  ='modis_cloud_particle_size_ice_mean',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='MDSCLDTOPPS'  ,                                 &
        LONG_NAME  ='modis_cloud_top_pressure_total_mean',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='MDSWTRPATH'  ,                                 &
        LONG_NAME  ='modis_liquid_water_path_mean',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
    VERIFY_(STATUS)

! the following is deprecated, use MDSWTRPATH
    call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='MDSH2OPATH'  ,                                 &
        LONG_NAME  ='modis_liquid_water_path_mean',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='MDSICEPATH'  ,                                 &
        LONG_NAME  ='modis_ice_water_path_mean',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='MDSTAUPRSHIST11'  ,                                 &
        LONG_NAME  ='modis_tau_pressure_histogram_bin_1_1',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)


     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='MDSTAUPRSHIST12'  ,                                 &
        LONG_NAME  ='modis_tau_pressure_histogram_bin_1_2',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)


     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='MDSTAUPRSHIST13'  ,                                 &
        LONG_NAME  ='modis_tau_pressure_histogram_bin_1_3',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)


     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='MDSTAUPRSHIST14'  ,                                 &
        LONG_NAME  ='modis_tau_pressure_histogram_bin_1_4',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)


     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='MDSTAUPRSHIST15'  ,                                 &
        LONG_NAME  ='modis_tau_pressure_histogram_bin_1_5',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)


     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='MDSTAUPRSHIST16'  ,                                 &
        LONG_NAME  ='modis_tau_pressure_histogram_bin_1_6',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)


    call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='MDSTAUPRSHIST17'  ,                                 &
        LONG_NAME  ='modis_tau_pressure_histogram_bin_1_7',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='MDSTAUPRSHIST21'  ,                                 &
        LONG_NAME  ='modis_tau_pressure_histogram_bin_2_1',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)


     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='MDSTAUPRSHIST22'  ,                                 &
        LONG_NAME  ='modis_tau_pressure_histogram_bin_2_2',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)


     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='MDSTAUPRSHIST23'  ,                                 &
        LONG_NAME  ='modis_tau_pressure_histogram_bin_2_3',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)


     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='MDSTAUPRSHIST24'  ,                                 &
        LONG_NAME  ='modis_tau_pressure_histogram_bin_2_4',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)


     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='MDSTAUPRSHIST25'  ,                                 &
        LONG_NAME  ='modis_tau_pressure_histogram_bin_2_5',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)


     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='MDSTAUPRSHIST26'  ,                                 &
        LONG_NAME  ='modis_tau_pressure_histogram_bin_2_6',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)


    call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='MDSTAUPRSHIST27'  ,                                 &
        LONG_NAME  ='modis_tau_pressure_histogram_bin_2_7',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='MDSTAUPRSHIST31'  ,                                 &
        LONG_NAME  ='modis_tau_pressure_histogram_bin_3_1',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)


     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='MDSTAUPRSHIST32'  ,                                 &
        LONG_NAME  ='modis_tau_pressure_histogram_bin_3_2',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)


     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='MDSTAUPRSHIST33'  ,                                 &
        LONG_NAME  ='modis_tau_pressure_histogram_bin_3_3',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)


     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='MDSTAUPRSHIST34'  ,                                 &
        LONG_NAME  ='modis_tau_pressure_histogram_bin_3_4',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)


     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='MDSTAUPRSHIST35'  ,                                 &
        LONG_NAME  ='modis_tau_pressure_histogram_bin_3_5',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)


     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='MDSTAUPRSHIST36'  ,                                 &
        LONG_NAME  ='modis_tau_pressure_histogram_bin_3_6',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)


    call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='MDSTAUPRSHIST37'  ,                                 &
        LONG_NAME  ='modis_tau_pressure_histogram_bin_3_7',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='MDSTAUPRSHIST41'  ,                                 &
        LONG_NAME  ='modis_tau_pressure_histogram_bin_4_1',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)


     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='MDSTAUPRSHIST42'  ,                                 &
        LONG_NAME  ='modis_tau_pressure_histogram_bin_4_2',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)


     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='MDSTAUPRSHIST43'  ,                                 &
        LONG_NAME  ='modis_tau_pressure_histogram_bin_4_3',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)


     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='MDSTAUPRSHIST44'  ,                                 &
        LONG_NAME  ='modis_tau_pressure_histogram_bin_4_4',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)


     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='MDSTAUPRSHIST45'  ,                                 &
        LONG_NAME  ='modis_tau_pressure_histogram_bin_4_5',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)


     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='MDSTAUPRSHIST46'  ,                                 &
        LONG_NAME  ='modis_tau_pressure_histogram_bin_4_6',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)


    call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='MDSTAUPRSHIST47'  ,                                 &
        LONG_NAME  ='modis_tau_pressure_histogram_bin_4_7',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='MDSTAUPRSHIST51'  ,                                 &
        LONG_NAME  ='modis_tau_pressure_histogram_bin_5_1',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)


     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='MDSTAUPRSHIST52'  ,                                 &
        LONG_NAME  ='modis_tau_pressure_histogram_bin_5_2',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)


     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='MDSTAUPRSHIST53'  ,                                 &
        LONG_NAME  ='modis_tau_pressure_histogram_bin_5_3',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)


     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='MDSTAUPRSHIST54'  ,                                 &
        LONG_NAME  ='modis_tau_pressure_histogram_bin_5_4',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)


     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='MDSTAUPRSHIST55'  ,                                 &
        LONG_NAME  ='modis_tau_pressure_histogram_bin_5_5',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)


     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='MDSTAUPRSHIST56'  ,                                 &
        LONG_NAME  ='modis_tau_pressure_histogram_bin_5_6',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)


    call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='MDSTAUPRSHIST57'  ,                                 &
        LONG_NAME  ='modis_tau_pressure_histogram_bin_5_7',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='MDSTAUPRSHIST61'  ,                                 &
        LONG_NAME  ='modis_tau_pressure_histogram_bin_6_1',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)


     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='MDSTAUPRSHIST62'  ,                                 &
        LONG_NAME  ='modis_tau_pressure_histogram_bin_6_2',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)


     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='MDSTAUPRSHIST63'  ,                                 &
        LONG_NAME  ='modis_tau_pressure_histogram_bin_6_3',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)


     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='MDSTAUPRSHIST64'  ,                                 &
        LONG_NAME  ='modis_tau_pressure_histogram_bin_6_4',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)


     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='MDSTAUPRSHIST65'  ,                                 &
        LONG_NAME  ='modis_tau_pressure_histogram_bin_6_5',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)


     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='MDSTAUPRSHIST66'  ,                                 &
        LONG_NAME  ='modis_tau_pressure_histogram_bin_6_6',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)


    call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='MDSTAUPRSHIST67'  ,                                 &
        LONG_NAME  ='modis_tau_pressure_histogram_bin_6_7',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='MDSTAUPRSHIST71'  ,                                 &
        LONG_NAME  ='modis_tau_pressure_histogram_bin_7_1',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)


     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='MDSTAUPRSHIST72'  ,                                 &
        LONG_NAME  ='modis_tau_pressure_histogram_bin_7_2',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)


     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='MDSTAUPRSHIST73'  ,                                 &
        LONG_NAME  ='modis_tau_pressure_histogram_bin_7_3',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)


     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='MDSTAUPRSHIST74'  ,                                 &
        LONG_NAME  ='modis_tau_pressure_histogram_bin_7_4',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)


     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='MDSTAUPRSHIST75'  ,                                 &
        LONG_NAME  ='modis_tau_pressure_histogram_bin_7_5',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)


     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='MDSTAUPRSHIST76'  ,                                 &
        LONG_NAME  ='modis_tau_pressure_histogram_bin_7_6',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)


    call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME ='MDSTAUPRSHIST77'  ,                                 &
        LONG_NAME  ='modis_tau_pressure_histogram_bin_7_7',   &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)









!
! MISR simulator output
!-------------------------------------------------------------------- 
      call MAPL_AddExportSpec(GC,                                 &
        SHORT_NAME ='MISRMNCLDTP'  ,                              &
        LONG_NAME  ='MISR_mead_cloud_top_height',                 &
        UNITS      ='m',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

      call MAPL_AddExportSpec(GC,                                 &
        SHORT_NAME ='MISRCLDAREA'  ,                              &
        LONG_NAME  ='MISR_cloud_area',                            &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

      call MAPL_AddExportSpec(GC,                                 &
        SHORT_NAME ='MISRLYRTP0'  ,                              &
        LONG_NAME  ='MISR_layer_top',                            &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                                 &
        SHORT_NAME ='MISRLYRTP250'  ,                            &
        LONG_NAME  ='MISR_layer_top',                            &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)
 
     call MAPL_AddExportSpec(GC,                                 &
        SHORT_NAME ='MISRLYRTP750'  ,                            &
        LONG_NAME  ='MISR_layer_top',                            &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)
 
     call MAPL_AddExportSpec(GC,                                 &
        SHORT_NAME ='MISRLYRTP1250'  ,                           &
        LONG_NAME  ='MISR_layer_top',                            &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)
 
     call MAPL_AddExportSpec(GC,                                 &
        SHORT_NAME ='MISRLYRTP1750',                              &
        LONG_NAME  ='MISR_layer_top',                            &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)
 
     call MAPL_AddExportSpec(GC,                                 &
        SHORT_NAME ='MISRLYRTP2250',                             &
        LONG_NAME  ='MISR_layer_top',                            &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)
 
     call MAPL_AddExportSpec(GC,                                 &
        SHORT_NAME ='MISRLYRTP2750',                             &
        LONG_NAME  ='MISR_layer_top',                            &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)
 
     call MAPL_AddExportSpec(GC,                                 &
        SHORT_NAME ='MISRLYRTP3500',                             &
        LONG_NAME  ='MISR_layer_top',                            &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)
 
     call MAPL_AddExportSpec(GC,                                 &
        SHORT_NAME ='MISRLYRTP4500',                             &
        LONG_NAME  ='MISR_layer_top',                            &
        UNITS      ='1',                                         &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)
 
     call MAPL_AddExportSpec(GC,                                 &
        SHORT_NAME ='MISRLYRTP6000' ,                              &
        LONG_NAME  ='MISR_layer_top',                            &
        UNITS      ='1',                                         &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)
 
     call MAPL_AddExportSpec(GC,                                 &
        SHORT_NAME ='MISRLYRTP8000'  ,                              &
        LONG_NAME  ='MISR_layer_top',                            &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)
 
     call MAPL_AddExportSpec(GC,                                 &
        SHORT_NAME ='MISRLYRTP10000'  ,                              &
        LONG_NAME  ='MISR_layer_top',                            &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)
 
     call MAPL_AddExportSpec(GC,                                 &
        SHORT_NAME ='MISRLYRTP12000'  ,                              &
        LONG_NAME  ='MISR_layer_top',                            &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)
 
     call MAPL_AddExportSpec(GC,                                 &
        SHORT_NAME ='MISRLYRTP14000'  ,                              &
        LONG_NAME  ='MISR_layer_top',                            &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)
 
     call MAPL_AddExportSpec(GC,                                 &
        SHORT_NAME ='MISRLYRTP16000'  ,                              &
        LONG_NAME  ='MISR_layer_top',                            &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)
 
     call MAPL_AddExportSpec(GC,                                 &
        SHORT_NAME ='MISRLYRTP18000'  ,                              &
        LONG_NAME  ='MISR_layer_top',                            &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                                 &
        SHORT_NAME ='MISRFQ0'  ,                              &
        LONG_NAME  ='MISR_cloud_area',                            &
        UNITS      ='1',                                          &
        UNGRIDDED_DIMS = (/ 7 /),                                 &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                                 &
        SHORT_NAME ='MISRFQ250'  ,                            &
        LONG_NAME  ='MISR_cloud_area',                            &
        UNITS      ='1',                                          &
        UNGRIDDED_DIMS = (/ 7 /),                                 &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)
 
     call MAPL_AddExportSpec(GC,                                 &
        SHORT_NAME ='MISRFQ750'  ,                            &
        LONG_NAME  ='MISR_cloud_area',                            &
        UNITS      ='1',                                          &
        UNGRIDDED_DIMS = (/ 7 /),                                 &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)
 
     call MAPL_AddExportSpec(GC,                                 &
        SHORT_NAME ='MISRFQ1250'  ,                           &
        LONG_NAME  ='MISR_cloud_area',                            &
        UNITS      ='1',                                          &
        UNGRIDDED_DIMS = (/ 7 /),                                 &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)
 
     call MAPL_AddExportSpec(GC,                                 &
        SHORT_NAME ='MISRFQ1750',                              &
        LONG_NAME  ='MISR_cloud_area',                            &
        UNITS      ='1',                                          &
        UNGRIDDED_DIMS = (/ 7 /),                                 &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)
 
     call MAPL_AddExportSpec(GC,                                 &
        SHORT_NAME ='MISRFQ2250',                             &
        LONG_NAME  ='MISR_cloud_area',                            &
        UNITS      ='1',                                          &
        UNGRIDDED_DIMS = (/ 7 /),                                 &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)
 
     call MAPL_AddExportSpec(GC,                                 &
        SHORT_NAME ='MISRFQ2750',                             &
        LONG_NAME  ='MISR_cloud_area',                            &
        UNITS      ='1',                                          &
        UNGRIDDED_DIMS = (/ 7 /),                                 &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)
 
     call MAPL_AddExportSpec(GC,                                 &
        SHORT_NAME ='MISRFQ3500',                             &
        LONG_NAME  ='MISR_cloud_area',                            &
        UNITS      ='1',                                          &
        UNGRIDDED_DIMS = (/ 7 /),                                 &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)
 
     call MAPL_AddExportSpec(GC,                                 &
        SHORT_NAME ='MISRFQ4500',                             &
        LONG_NAME  ='MISR_cloud_area',                            &
        UNITS      ='1',                                         &
        UNGRIDDED_DIMS = (/ 7 /),                                 &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)
 
     call MAPL_AddExportSpec(GC,                                 &
        SHORT_NAME ='MISRFQ6000' ,                              &
        LONG_NAME  ='MISR_cloud_area',                            &
        UNITS      ='1',                                         &
        UNGRIDDED_DIMS = (/ 7 /),                                 &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)
 
     call MAPL_AddExportSpec(GC,                                 &
        SHORT_NAME ='MISRFQ8000'  ,                              &
        LONG_NAME  ='MISR_cloud_area',                            &
        UNITS      ='1',                                          &
        UNGRIDDED_DIMS = (/ 7 /),                                 &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)
 
     call MAPL_AddExportSpec(GC,                                 &
        SHORT_NAME ='MISRFQ10000'  ,                              &
        LONG_NAME  ='MISR_cloud_area',                            &
        UNITS      ='1',                                          &
        UNGRIDDED_DIMS = (/ 7 /),                                 &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)
 
     call MAPL_AddExportSpec(GC,                                 &
        SHORT_NAME ='MISRFQ12000'  ,                              &
        LONG_NAME  ='MISR_cloud_area',                            &
        UNITS      ='1',                                          &
        UNGRIDDED_DIMS = (/ 7 /),                                 &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)
 
     call MAPL_AddExportSpec(GC,                                 &
        SHORT_NAME ='MISRFQ14000'  ,                              &
        LONG_NAME  ='MISR_cloud_area',                            &
        UNITS      ='1',                                          &
        UNGRIDDED_DIMS = (/ 7 /),                                 &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)
 
     call MAPL_AddExportSpec(GC,                                 &
        SHORT_NAME ='MISRFQ16000'  ,                              &
        LONG_NAME  ='MISR_cloud_area',                            &
        UNITS      ='1',                                          &
        UNGRIDDED_DIMS = (/ 7 /),                                 &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)
 
     call MAPL_AddExportSpec(GC,                                 &
        SHORT_NAME ='MISRFQ18000'  ,                              &
        LONG_NAME  ='MISR_cloud_area',                            &
        UNITS      ='1',                                          &
        UNGRIDDED_DIMS = (/ 7 /),                                 &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,             RC=STATUS  )
     VERIFY_(STATUS)
 
     call MAPL_AddImportSpec(GC, &
        SHORT_NAME     = 'SATORB', &
        LONG_NAME      = 'Satellite_orbits', &
        UNITS          = 'days' , &
        DIMS           = MAPL_DimsHorzOnly , &
        DATATYPE       = MAPL_BundleItem , &
                RC = STATUS )
    VERIFY_(STATUS)
!EOS
!
! Set the Profiling timers
! ------------------------

    call MAPL_TimerAdd(GC,name="DRIVER" ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,name="-MISC"    ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,name="-COSP"   ,RC=STATUS)
    VERIFY_(STATUS)

! parse resource file for masks
    CF = ESMF_ConfigCreate(rc=status)
    VERIFY_(STATUS)
    inquire(file='SatSim.rc',exist=found)
    if (found) then
       call ESMF_ConfigLoadFile ( CF,'SatSim.rc',rc=status)
       VERIFY_(STATUS) 
       call ESMF_ConfigGetDim(CF, nLines,nCols, LABEL='Masked_Exports::',RC=STATUS)
       if (status== ESMF_SUCCESS) then
          self%nmask_vars = nLines
          allocate(self%export_name(nLines),self%mask_name(nLines),self%newvar_name(nLines),self%newvar(nLines),stat=status)
          VERIFY_(STATUS)
          call ESMF_ConfigFindLabel(CF, 'Masked_Exports::',__RC__)
          VERIFY_(STATUS)
          do m=1,nLines
             call ESMF_ConfigNextLine(CF,RC=STATUS)
             VERIFY_(STATUS)
             call ESMF_ConfigGetAttribute(CF,self%export_name(m),RC=STATUS)
             VERIFY_(STATUS)
             call ESMF_ConfigGetAttribute(CF,self%mask_name(m),RC=STATUS)
             VERIFY_(STATUS)
             call ESMF_ConfigGetAttribute(CF,tmpname,RC=STATUS)
             if (status == ESMF_SUCCESS) then
                self%newvar_name(m)=tmpname
                self%newvar(m)=.true.
             else
                self%newvar_name(m)=self%export_name(m)
                self%newvar(m)=.false.
             end if
          end do
       end if
    else
       self%nmask_vars = 0
    end if
    if (MAPL_AM_I_Root()) then
       write(*,*)'Parsing of satsim.rc'
       do m=1,self%nmask_vars
          write(*,*)m,self%newvar(m)
          write(*,*)trim(self%export_name(m)),' ',trim(self%mask_name(m)),' ',trim(self%newvar_name(m))
       end do
    end if

    call MAPL_GetObjectFromGC ( GC, STATE, RC=STATUS)
    VERIFY_(STATUS)

    do m=1,self%nmask_vars
       if (self%newvar(m)) then
          call MAPL_Get(STATE,ExportSpec=ExportSpec,RC=STATUS)
          VERIFY_(STATUS)
          vindex = MAPL_VarSpecGetIndex(ExportSpec,self%export_name(m),RC=STATUS)
          VERIFY_(STATUS)
          call MAPL_VarSpecGet(ExportSpec(vindex),LONG_NAME=long_name, &
             UNITS=units,Dims=MAPL_dims,VLocation=VLocation,ungridded_dims=Ungridded_dims,RC=STATUS)
          VERIFY_(STATUS)
          if (associated(Ungridded_dims)) then
             call MAPL_AddExportSpec(GC,                       &
             SHORT_NAME =self%newvar_name(m)  ,                &
             LONG_NAME  =long_name,                            &
             UNITS      =units,                                &
             UNGRIDDED_DIMS = Ungridded_dims,                  &
             DIMS       = MAPL_dims,                           &
             VLOCATION  = VLocation,             RC=STATUS  )
             VERIFY_(STATUS)
             nullify(ungridded_dims)
          else
             call MAPL_AddExportSpec(GC,                       &
             SHORT_NAME =self%newvar_name(m)  ,                &
             LONG_NAME  =long_name,                            &
             UNITS      =units,                                &
             DIMS       = MAPL_dims,                           &
             VLOCATION  = VLocation,             RC=STATUS  )
             VERIFY_(STATUS)
          end if
          exportName(1) = self%newvar_name(m)
          call MAPL_DoNotDeferExport(GC,exportName,RC=STATUS)
          VERIFY_(STATUS)
       end if
       call MAPL_DoNotDeferExport(GC,self%export_name,RC=STATUS)
       VERIFY_(STATUS)
    end do

!   Store internal state in GC
!   --------------------------
    call ESMF_UserCompSetInternalState ( GC, 'SatSim_State', wrap, STATUS )
    VERIFY_(STATUS)

! Set generic init and final methods
! ----------------------------------

    call MAPL_GenericSetServices    ( GC, RC=STATUS)
    VERIFY_(STATUS)
     
    RETURN_(ESMF_SUCCESS)
     
  end subroutine SetServices

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!BOP
! !IROUTINE: RUN -- Run method for the SATSIM component

! !INTERFACE:
  subroutine RUN ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:
    type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
    type(ESMF_State),    intent(inout) :: IMPORT ! Import state
    type(ESMF_State),    intent(inout) :: EXPORT ! Export state
    type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
    integer, optional,   intent(  out) :: RC     ! Error code:
    
! !DESCRIPTION: This version uses the MAPL\_GenericSetServices. This function sets
!                the Initialize and Finalize services, as well as allocating

!EOP


! ErrLog Variables

    character(len=ESMF_MAXSTR)          :: IAm
    integer                             :: STATUS
    character(len=ESMF_MAXSTR)          :: COMP_NAME

! Local derived type aliases

    type (MAPL_MetaComp), pointer   :: STATE
    type (ESMF_Config      )            :: CF
    type (ESMF_State       )            :: INTERNAL
    type (ESMF_Alarm)                   :: ALARM

! Local variables

    integer                             :: IM,JM,LM
    real, pointer, dimension(:,:)       :: LONS
    real, pointer, dimension(:,:)       :: LATS


! pointers for Export outputs

! ISCCP/Icarus

!CU family
    real, pointer, dimension(:,:)       :: ISCCP_CU_OA,  ISCCP_CU_OB 
    real, pointer, dimension(:,:)       :: ISCCP_CU_UA,  ISCCP_CU_UB 

!STCU family
    real, pointer, dimension(:,:)       :: ISCCP_STCU_OA,ISCCP_STCU_OB 
    real, pointer, dimension(:,:)       :: ISCCP_STCU_UA,ISCCP_STCU_UB 
                                           
!ST family
    real, pointer, dimension(:,:)       :: ISCCP_ST_OA,  ISCCP_ST_OB 
    real, pointer, dimension(:,:)       :: ISCCP_ST_UA,  ISCCP_ST_UB 
                                           
!ACU family
    real, pointer, dimension(:,:)       :: ISCCP_ACU_OA,   ISCCP_ACU_OB 
    real, pointer, dimension(:,:)       :: ISCCP_ACU_UA,   ISCCP_ACU_UB 

!AST family
    real, pointer, dimension(:,:)       :: ISCCP_AST_OA,   ISCCP_AST_OB 
    real, pointer, dimension(:,:)       :: ISCCP_AST_UA,   ISCCP_AST_UB 
                                           
!NST family
    real, pointer, dimension(:,:)       :: ISCCP_NST_OA,   ISCCP_NST_OB 
    real, pointer, dimension(:,:)       :: ISCCP_NST_UA,   ISCCP_NST_UB 
                                           
!CI family
    real, pointer, dimension(:,:)       :: ISCCP_CI_OA,    ISCCP_CI_OB 
    real, pointer, dimension(:,:)       :: ISCCP_CI_MA,    ISCCP_CI_MB 
    real, pointer, dimension(:,:)       :: ISCCP_CI_UA,    ISCCP_CI_UB 

!CIST family
    real, pointer, dimension(:,:)       :: ISCCP_CIST_OA,  ISCCP_CIST_OB 
    real, pointer, dimension(:,:)       :: ISCCP_CIST_MA,  ISCCP_CIST_MB 
    real, pointer, dimension(:,:)       :: ISCCP_CIST_UA,  ISCCP_CIST_UB 

!Cb family
    real, pointer, dimension(:,:)       :: ISCCP_CB_OA,    ISCCP_CB_OB 
    real, pointer, dimension(:,:)       :: ISCCP_CB_MA,    ISCCP_CB_MB 
    real, pointer, dimension(:,:)       :: ISCCP_CB_UA,    ISCCP_CB_UB 

!SubVisible Family
    real, pointer, dimension(:,:)       :: ISCCP_SUBV1,    ISCCP_SUBV2
    real, pointer, dimension(:,:)       :: ISCCP_SUBV3,    ISCCP_SUBV4
    real, pointer, dimension(:,:)       :: ISCCP_SUBV5,    ISCCP_SUBV6
    real, pointer, dimension(:,:)       :: ISCCP_SUBV7

! "stacked" fields (im*jm*7)
    real, pointer, dimension(:,:,:)       :: CLISCCP1,    CLISCCP2
    real, pointer, dimension(:,:,:)       :: CLISCCP3,    CLISCCP4
    real, pointer, dimension(:,:,:)       :: CLISCCP5,    CLISCCP6
    real, pointer, dimension(:,:,:)       :: CLISCCP7

! "stacked" fields (im*jm*7)
    real, pointer, dimension(:,:,:)       :: ISCCP1,    ISCCP2
    real, pointer, dimension(:,:,:)       :: ISCCP3,    ISCCP4
    real, pointer, dimension(:,:,:)       :: ISCCP5,    ISCCP6
    real, pointer, dimension(:,:,:)       :: ISCCP7

!other 2Ds
    real, pointer, dimension(:,:)       :: TCLISCCP, CLTISCCP
    real, pointer, dimension(:,:)       :: CTPISCCP, PCTISCCP
    real, pointer, dimension(:,:)       :: ALBISCCP
    real, pointer, dimension(:,:)       :: TBISCCP


! LIDAR/Calipso exports

    real, pointer, dimension(:,:)       :: RADARLTCC
    real, pointer, dimension(:,:)       :: CLLCALIPSO
    real, pointer, dimension(:,:)       :: CLMCALIPSO
    real, pointer, dimension(:,:)       :: CLHCALIPSO
    real, pointer, dimension(:,:)       :: CLTCALIPSO

    real, pointer, dimension(:,:,:)       :: PARASOLREFL0
    real, pointer, dimension(:,:)       :: PARASOLREFL1
    real, pointer, dimension(:,:)       :: PARASOLREFL2
    real, pointer, dimension(:,:)       :: PARASOLREFL3
    real, pointer, dimension(:,:)       :: PARASOLREFL4
    real, pointer, dimension(:,:)       :: PARASOLREFL5

    real, pointer, dimension(:,:,:)     :: SGFCLD
    real, pointer, dimension(:,:,:)     :: LIDARPMOL,   LIDARPTOT
    real, pointer, dimension(:,:,:)     :: LIDARTAUTOT 
    real, pointer, dimension(:,:,:)     :: RADARZETOT
    real, pointer, dimension(:,:,:)     :: CLCALIPSO2
    real, pointer, dimension(:,:,:)     :: CLCALIPSO

    real, pointer, dimension(:,:,:)     :: CFADLIDARSR532_01
    real, pointer, dimension(:,:,:)     :: CFADLIDARSR532_02
    real, pointer, dimension(:,:,:)     :: CFADLIDARSR532_03
    real, pointer, dimension(:,:,:)     :: CFADLIDARSR532_04
    real, pointer, dimension(:,:,:)     :: CFADLIDARSR532_05
    real, pointer, dimension(:,:,:)     :: CFADLIDARSR532_06
    real, pointer, dimension(:,:,:)     :: CFADLIDARSR532_07
    real, pointer, dimension(:,:,:)     :: CFADLIDARSR532_08
    real, pointer, dimension(:,:,:)     :: CFADLIDARSR532_09
    real, pointer, dimension(:,:,:)     :: CFADLIDARSR532_10
    real, pointer, dimension(:,:,:)     :: CFADLIDARSR532_11
    real, pointer, dimension(:,:,:)     :: CFADLIDARSR532_12
    real, pointer, dimension(:,:,:)     :: CFADLIDARSR532_13
    real, pointer, dimension(:,:,:)     :: CFADLIDARSR532_14
    real, pointer, dimension(:,:,:)     :: CFADLIDARSR532_15

! Radar/cloudsat exports

    real, pointer, dimension(:,:,:)     :: CLOUDSATCFAD01
    real, pointer, dimension(:,:,:)     :: CLOUDSATCFAD02
    real, pointer, dimension(:,:,:)     :: CLOUDSATCFAD03
    real, pointer, dimension(:,:,:)     :: CLOUDSATCFAD04
    real, pointer, dimension(:,:,:)     :: CLOUDSATCFAD05
    real, pointer, dimension(:,:,:)     :: CLOUDSATCFAD06
    real, pointer, dimension(:,:,:)     :: CLOUDSATCFAD07
    real, pointer, dimension(:,:,:)     :: CLOUDSATCFAD08
    real, pointer, dimension(:,:,:)     :: CLOUDSATCFAD09
    real, pointer, dimension(:,:,:)     :: CLOUDSATCFAD10
    real, pointer, dimension(:,:,:)     :: CLOUDSATCFAD11
    real, pointer, dimension(:,:,:)     :: CLOUDSATCFAD12
    real, pointer, dimension(:,:,:)     :: CLOUDSATCFAD13
    real, pointer, dimension(:,:,:)     :: CLOUDSATCFAD14
    real, pointer, dimension(:,:,:)     :: CLOUDSATCFAD15

! MODIS simulator exports

    real, pointer, dimension(:,:)       :: MDSCLDFRCTTL
    real, pointer, dimension(:,:)       :: MDSCLDFRCWTR
    real, pointer, dimension(:,:)       :: MDSCLDFRCH2O
    real, pointer, dimension(:,:)       :: MDSCLDFRCICE
    real, pointer, dimension(:,:)       :: MDSCLDFRCHI
    real, pointer, dimension(:,:)       :: MDSCLDFRCMID
    real, pointer, dimension(:,:)       :: MDSCLDFRCLO
    real, pointer, dimension(:,:)       :: MDSOPTHCKTTL
    real, pointer, dimension(:,:)       :: MDSOPTHCKWTR
    real, pointer, dimension(:,:)       :: MDSOPTHCKH2O
    real, pointer, dimension(:,:)       :: MDSOPTHCKICE
    real, pointer, dimension(:,:)       :: MDSOPTHCKTTLLG
    real, pointer, dimension(:,:)       :: MDSOPTHCKWTRLG
    real, pointer, dimension(:,:)       :: MDSOPTHCKH2OLG
    real, pointer, dimension(:,:)       :: MDSOPTHCKICELG
    real, pointer, dimension(:,:)       :: MDSCLDSZWTR
    real, pointer, dimension(:,:)       :: MDSCLDSZH20
    real, pointer, dimension(:,:)       :: MDSCLDSZICE
    real, pointer, dimension(:,:)       :: MDSCLDTOPPS
    real, pointer, dimension(:,:)       :: MDSWTRPATH
    real, pointer, dimension(:,:)       :: MDSH2OPATH
    real, pointer, dimension(:,:)       :: MDSICEPATH

    real, pointer, dimension(:,:)       :: MDSTAUPRSHIST11
    real, pointer, dimension(:,:)       :: MDSTAUPRSHIST12
    real, pointer, dimension(:,:)       :: MDSTAUPRSHIST13
    real, pointer, dimension(:,:)       :: MDSTAUPRSHIST14
    real, pointer, dimension(:,:)       :: MDSTAUPRSHIST15
    real, pointer, dimension(:,:)       :: MDSTAUPRSHIST16
    real, pointer, dimension(:,:)       :: MDSTAUPRSHIST17

    real, pointer, dimension(:,:)       :: MDSTAUPRSHIST21
    real, pointer, dimension(:,:)       :: MDSTAUPRSHIST22
    real, pointer, dimension(:,:)       :: MDSTAUPRSHIST23
    real, pointer, dimension(:,:)       :: MDSTAUPRSHIST24
    real, pointer, dimension(:,:)       :: MDSTAUPRSHIST25
    real, pointer, dimension(:,:)       :: MDSTAUPRSHIST26
    real, pointer, dimension(:,:)       :: MDSTAUPRSHIST27

    real, pointer, dimension(:,:)       :: MDSTAUPRSHIST31
    real, pointer, dimension(:,:)       :: MDSTAUPRSHIST32
    real, pointer, dimension(:,:)       :: MDSTAUPRSHIST33
    real, pointer, dimension(:,:)       :: MDSTAUPRSHIST34
    real, pointer, dimension(:,:)       :: MDSTAUPRSHIST35
    real, pointer, dimension(:,:)       :: MDSTAUPRSHIST36
    real, pointer, dimension(:,:)       :: MDSTAUPRSHIST37

    real, pointer, dimension(:,:)       :: MDSTAUPRSHIST41
    real, pointer, dimension(:,:)       :: MDSTAUPRSHIST42
    real, pointer, dimension(:,:)       :: MDSTAUPRSHIST43
    real, pointer, dimension(:,:)       :: MDSTAUPRSHIST44
    real, pointer, dimension(:,:)       :: MDSTAUPRSHIST45
    real, pointer, dimension(:,:)       :: MDSTAUPRSHIST46
    real, pointer, dimension(:,:)       :: MDSTAUPRSHIST47

    real, pointer, dimension(:,:)       :: MDSTAUPRSHIST51
    real, pointer, dimension(:,:)       :: MDSTAUPRSHIST52
    real, pointer, dimension(:,:)       :: MDSTAUPRSHIST53
    real, pointer, dimension(:,:)       :: MDSTAUPRSHIST54
    real, pointer, dimension(:,:)       :: MDSTAUPRSHIST55
    real, pointer, dimension(:,:)       :: MDSTAUPRSHIST56
    real, pointer, dimension(:,:)       :: MDSTAUPRSHIST57

    real, pointer, dimension(:,:)       :: MDSTAUPRSHIST61
    real, pointer, dimension(:,:)       :: MDSTAUPRSHIST62
    real, pointer, dimension(:,:)       :: MDSTAUPRSHIST63
    real, pointer, dimension(:,:)       :: MDSTAUPRSHIST64
    real, pointer, dimension(:,:)       :: MDSTAUPRSHIST65
    real, pointer, dimension(:,:)       :: MDSTAUPRSHIST66
    real, pointer, dimension(:,:)       :: MDSTAUPRSHIST67

    real, pointer, dimension(:,:)       :: MDSTAUPRSHIST71
    real, pointer, dimension(:,:)       :: MDSTAUPRSHIST72
    real, pointer, dimension(:,:)       :: MDSTAUPRSHIST73
    real, pointer, dimension(:,:)       :: MDSTAUPRSHIST74
    real, pointer, dimension(:,:)       :: MDSTAUPRSHIST75
    real, pointer, dimension(:,:)       :: MDSTAUPRSHIST76
    real, pointer, dimension(:,:)       :: MDSTAUPRSHIST77









! MISR simulator exports

    real, pointer, dimension(:,:)       :: MISRMNCLDTP
    real, pointer, dimension(:,:)       :: MISRCLDAREA
    real, pointer, dimension(:,:)       :: MISRLYRTP0
    real, pointer, dimension(:,:)       :: MISRLYRTP250
    real, pointer, dimension(:,:)       :: MISRLYRTP750
    real, pointer, dimension(:,:)       :: MISRLYRTP1250
    real, pointer, dimension(:,:)       :: MISRLYRTP1750
    real, pointer, dimension(:,:)       :: MISRLYRTP2250
    real, pointer, dimension(:,:)       :: MISRLYRTP2750
    real, pointer, dimension(:,:)       :: MISRLYRTP3500
    real, pointer, dimension(:,:)       :: MISRLYRTP4500
    real, pointer, dimension(:,:)       :: MISRLYRTP6000
    real, pointer, dimension(:,:)       :: MISRLYRTP8000
    real, pointer, dimension(:,:)       :: MISRLYRTP10000
    real, pointer, dimension(:,:)       :: MISRLYRTP12000
    real, pointer, dimension(:,:)       :: MISRLYRTP14000
    real, pointer, dimension(:,:)       :: MISRLYRTP16000
    real, pointer, dimension(:,:)       :: MISRLYRTP18000
    real, pointer, dimension(:,:,:)       :: MISRFQ0
    real, pointer, dimension(:,:,:)       :: MISRFQ250
    real, pointer, dimension(:,:,:)       :: MISRFQ750
    real, pointer, dimension(:,:,:)       :: MISRFQ1250
    real, pointer, dimension(:,:,:)       :: MISRFQ1750
    real, pointer, dimension(:,:,:)       :: MISRFQ2250
    real, pointer, dimension(:,:,:)       :: MISRFQ2750
    real, pointer, dimension(:,:,:)       :: MISRFQ3500
    real, pointer, dimension(:,:,:)       :: MISRFQ4500
    real, pointer, dimension(:,:,:)       :: MISRFQ6000
    real, pointer, dimension(:,:,:)       :: MISRFQ8000
    real, pointer, dimension(:,:,:)       :: MISRFQ10000
    real, pointer, dimension(:,:,:)       :: MISRFQ12000
    real, pointer, dimension(:,:,:)       :: MISRFQ14000
    real, pointer, dimension(:,:,:)       :: MISRFQ16000
    real, pointer, dimension(:,:,:)       :: MISRFQ18000
 
!=============================================================================

! Begin... 

! Get my name and set-up traceback handle
! ---------------------------------------

    Iam = 'Run'
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // Iam


! Get my internal MAPL_Generic state
!-----------------------------------

    call MAPL_GetObjectFromGC ( GC, STATE, RC=STATUS)
    VERIFY_(STATUS)


! Get parameters from generic state.
!-----------------------------------

    call MAPL_Get( STATE, IM=IM, JM=JM, LM=LM,   &
                               CF       = CF,                &
                               LONS     = LONS,              &
                               LATS     = LATS,              &
                               INTERNAL_ESMF_STATE=INTERNAL, &
                               RUNALARM = ALARM,             &
                                                   RC=STATUS )
    VERIFY_(STATUS)

    if( .not. ESMF_AlarmIsRinging( ALARM)) then
       RETURN_(ESMF_SUCCESS)
    end if

    call MAPL_TimerOn(STATE,"DRIVER")

    call SIM_DRIVER(IM,JM,LM, RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_TimerOff(STATE,"DRIVER")

    RETURN_(ESMF_SUCCESS)

  contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine SIM_DRIVER(IM,JM,LM, RC)
      integer,           intent(IN ) :: IM, JM, LM
      integer, optional, intent(OUT) :: RC


!  Locals

      character(len=ESMF_MAXSTR)      :: IAm
      character(len=ESMF_MAXSTR)      :: NAME
      integer                         :: STATUS

! Import arrays

! These are the gen-u-ine import arrays in GEOS-5 dimensions
! PLE and ZLE are indexed 0:LM 
      real, pointer, dimension(:,:,:) :: T, PLE, QV, FCLD, ZLE
      real, pointer, dimension(:,:,:) :: RDFL, RDFI, QLLS, QILS, QLCN, QICN
      real, pointer, dimension(:,:,:) :: RDFR, RDFS
      real, pointer, dimension(:,:,:) :: QRTOT, QSTOT
      real, dimension(IM,JM,LM) :: QLTOT, QITOT
      real, pointer, dimension(:,:) :: MCOSZ,FRLAND,TS,FROCEAN

! These are the same, in logical 2 dimensions for icarus's sake
      real, dimension(IM*JM,0:LM) :: ZLE2D
      real, dimension(IM*JM) :: MCOSZCOSP,FRLANDCOSP,TSCOSP,FROCEANCOSP

! These are the same, converted to upside-down and 2d 
      real, dimension(IM*JM,LM) :: TCOSP, PLOCOSP, QVCOSP, FCLDCOSP
      real, dimension(IM*JM,LM) :: RDFLCOSP, RDFICOSP, QLLSCOSP, QILSCOSP, QLCNCOSP, QICNCOSP
      real, dimension(IM*JM,LM) :: QLTOTCOSP, QITOTCOSP
      real, dimension(IM*JM,LM) :: QRTOTCOSP, QSTOTCOSP
      real, dimension(IM*JM,LM) :: RDFRCOSP, RDFSCOSP
      real, dimension(IM*JM,0:LM) :: PLECOSP, PLE2DTOA0INV, ZLECOSP 

! Local variables needed for all simulators
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    real, dimension(IM,JM,LM) :: PLO, RH, QST
    real, dimension(IM*JM,LM) :: RHCOSP
    integer :: Npoints 
    real, dimension(:,:,:),target,allocatable :: frac_out
    !real, dimension(:,:,:),target,allocatable :: frac_outinv
    real, dimension(IM*JM,LM) :: frac_ls
    integer  :: isccp_overlap
    logical                         :: DEBUG_GC
    integer  :: i,j,k,l,icb,ict
    integer  :: scops_debug = 0
    integer  :: idebug      = 0
    integer  :: idebugcol   = 0
    real,dimension(:,:),pointer :: column_frac_out

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!  The follwing assignments are for variables needed for COSP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    integer  :: melt_lay        ! melting layer model off=0, on=1
    integer  :: Naero    ! Number of aerosol species
    integer  :: surface_radar   ! surface=1,spaceborne=0
    ! Time [days] - need to do something about this - afe
    double precision :: time=1.D0
    double precision :: time_bands(2,1)=1.D0


! COSP types

      type(cosp_config)               :: cfg     ! Configuration options
      type(cosp_gridbox)              :: gbx     ! Gridbox information. Input for COSP
      type(cosp_subgrid)              :: sgx     ! Subgrid outputs
      type(cosp_sgradar)              :: sgradar ! Output from radar simulator
      type(cosp_sghydro)              :: sghydro ! Input to radar simulator
      type(cosp_sglidar)              :: sglidar ! Output from lidar simulator
      type(cosp_isccp)                :: isccp   ! Output from ISCCP simulator
      type(cosp_vgrid)                :: vgrid   ! Information on vertical grid of stats
      type(cosp_radarstats)           :: stradar ! Summary statistics from radar simulator
      type(cosp_lidarstats)           :: stlidar ! Summary statistics from lidar simulator
      type(cosp_modis)   :: modis   ! Output from MODIS simulator
      type(cosp_misr)    :: misr    ! Output from MISR simulator

! RRTOV paramters needed only by COSP, otherwise unused
! -------------------------------------------------------

    integer :: Npoints_it
    logical :: use_precipitation_fluxes,use_reff
    integer :: Plat
    integer :: Sat
    integer :: Inst
    integer :: Nchan
    integer :: Ichan(0)
    real    :: SurfEm(0)
    real    :: ZenAng
    real    :: co2,ch4,n2o,co 
!


!  The follwing definitions are for variables needed for COSP and isccp/icarus
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real :: CCA(IM,JM,LM) ! convective clouds fraction, goes into scops
      real :: CCACOSP(IM*JM,LM) ! convective clouds fraction, goes into scops
      real :: DTAU_S(IM,JM,LM)
      real :: dtau_sCOSP(IM*JM,LM)
      real, dimension(IM,JM,LM,4) :: CWC, REFF
      real, dimension(IM,JM,LM,4) :: tausw,taulw
      real, dimension(      LM,4) :: dumtaubeam
      real, dimension(      LM  ) :: dumasycl
      real, dimension(    0:LM  ) :: dumtcldlyr
      real, dimension(    0:LM  ) :: dumenn
      real                        :: taucir
      real, dimension(IM,JM,LM,10) :: TAUDIAG
      real :: sunlit(IM*JM)
      real, dimension(IM,JM,LM) :: DELP, TAUSWICE, TAUSWLIQ, EMISS
      real, dimension(IM*JM,LM) :: EMISS2D, EMISSCOSP
    real    :: isccp_emsfc_lw     
   integer  :: isccp_top_height
    integer  :: isccp_top_height_direction

! icarus/isccp outputs
      real, dimension(IM*JM      ) :: ISCCP_totalcldarea
      real, dimension(IM*JM      ) :: ISCCP_meanptop
      real, dimension(IM*JM,7,7  ) :: fq_isccp
      real, dimension(IM,JM,7,7  ) :: fq_isccp3D
      real, dimension(IM*JM      ) :: ISCCP_meanalbedocld
      real, dimension(IM*JM      ) :: ISCCP_meanallskybrighttemp
 

! variables for lidar simulator
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real, dimension(IM*JM,LM   ) :: lidar_pmol 
      real, dimension(:,:,:), allocatable :: lidar_beta_tot,lidar_tau_tot
      real, dimension(IM*JM,PARASOL_NREFL) :: REFL

 
! variables for lidar and radar simulators
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    real, dimension(IM*JM,LM) :: ZLOCOSP
    real, dimension(IM*JM,LM) :: ZLO2D
     real                              :: K2
    integer  :: use_mie_tables  ! use a precomputed lookup table? yes=1,no=0,2=use first 
                                !  column everywhere
    integer  :: use_gas_abs     ! include gaseous absorption? yes=1,no=0
    integer  :: do_ray          ! calculate/output Rayleigh refl=1, not=0
    integer  :: Nprmts_max_hydro    ! Max number of parameters for hydrometeor size distributions
    integer  :: Nprmts_max_aero    ! Max number of parameters for aerosol size distributions
    real                              :: radar_freq
   integer  :: lidar_ice_type !Ice particle shape in lidar calculations 
                               ! (0=ice-spheres;1=ice-non-spherical)



! radar and lidar simulator output -- the admittedly confusing lidar_ and radar_
! suffixes indicate the cosp return structures that the variables are copied out of
    real,    dimension(:,:,:), allocatable :: radar_ze_tot
    real,    dimension(:,:,:), allocatable :: radar_ze_tot_tmp
    logical, dimension(:,:,:), allocatable :: radar_ze_tot_mask
    integer, dimension(IM*JM,LM) :: ncvalid
    real,    dimension(IM*JM,LM) :: radar_ze_tot_max
    real,    dimension(IM*JM,LM) :: radar_ze_tot_mean

    real, dimension (IM*JM,dBZe_bins,LM) :: radar_cfad_ze
    real, dimension(IM*JM,LM) :: radar_lidar_only_freq_cloud
    real,dimension(IM*JM) :: radar_lidar_tcc
    real,dimension(IM*JM) :: land_mask  

    real, dimension(IM*JM,SR_BINS,LM) :: lidar_cfad_sr        ! CFAD of scattering ratio
    real, dimension(IM*JM,LM) :: lidar_lidarcld               ! 3D "lidar" cloud fraction 
    real, dimension(IM*JM,LIDAR_NCAT) :: lidar_cldlayer       ! low, mid, high-level lidar cloud cover
    real, dimension(IM*JM,PARASOL_NREFL) :: lidar_parasolrefl ! mean parasol reflectance

    real, dimension(IM*JM,7,MISR_N_CTH) :: fq_MISR
    real, dimension(IM*JM,MISR_N_CTH) :: MISR_dist_model_layertops
    real, dimension(IM*JM):: MISR_meanztop
    real, dimension(IM*JM):: MISR_cldarea

    real, dimension(IM*JM):: MODIS_Cloud_Fraction_Total_Mean
    real, dimension(IM*JM):: MODIS_Cloud_Fraction_Water_Mean
    real, dimension(IM*JM):: MODIS_Cloud_Fraction_Ice_Mean
    real, dimension(IM*JM):: MODIS_Cloud_Fraction_High_Mean
    real, dimension(IM*JM):: MODIS_Cloud_Fraction_Mid_Mean
    real, dimension(IM*JM):: MODIS_Cloud_Fraction_Low_Mean
    real, dimension(IM*JM):: MODIS_Optical_Thickness_Total_Mean
    real, dimension(IM*JM):: MODIS_Optical_Thickness_Water_Mean
    real, dimension(IM*JM):: MODIS_Optical_Thickness_Ice_Mean
    real, dimension(IM*JM):: MODIS_Optical_Thickness_Total_LogMean
    real, dimension(IM*JM):: MODIS_Optical_Thickness_Water_LogMean
    real, dimension(IM*JM):: MODIS_Optical_Thickness_Ice_LogMean
    real, dimension(IM*JM):: MODIS_Cloud_Particle_Size_Water_Mean
    real, dimension(IM*JM):: MODIS_Cloud_Particle_Size_Ice_Mean
    real, dimension(IM*JM):: MODIS_Cloud_Top_Pressure_Total_Mean
    real, dimension(IM*JM):: MODIS_Liquid_Water_Path_Mean
    real, dimension(IM*JM):: MODIS_Ice_Water_Path_Mean
    real, dimension(IM*JM,7,7):: MODIS_Optical_Thickness_vs_Cloud_Top_Pressure


    integer :: BEGSEG, ENDSEG

    type(SatSim_state), pointer     :: self => null()
    type(SatSim_Wrap)               :: wrap
    integer                       :: mapl_dims, vindex
    integer, pointer              :: ungridded_dims(:)
    type(MAPL_VarSpec), pointer :: ExportSpec(:)
    real, pointer, dimension(:,:)   :: ptr2d, ptr2d_new
    real, pointer, dimension(:,:)   :: ptr_mask
    real, pointer, dimension(:,:,:) :: ptr3d, ptr3d_new
    type(ESMF_FieldBundle) :: bundle
    type (MAPL_MetaComp),      pointer      :: MAPL


    integer  :: use_satsim, use_satsim_isccp, use_satsim_modis, use_satsim_lidar, use_satsim_radar, use_satsim_misr
    integer  :: ncolumns

    character(len=ESMF_MAXSTR) :: GRIDNAME
    character(len=5)           :: imchar
    character(len=2)           :: dateline
    integer                    :: imsize,nn

!  Begin...
!----------
      Iam = trim(COMP_NAME) // 'Sim_Driver'
      DEBUG_GC=.FALSE.



! Get my MAPL_Generic state
!--------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_GetResource(MAPL,USE_SATSIM,LABEL="USE_SATSIM:",default=0,   RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_GetResource(MAPL,USE_SATSIM_ISCCP,LABEL="USE_SATSIM_ISCCP:",default=0,   RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_GetResource(MAPL,USE_SATSIM_MODIS,LABEL="USE_SATSIM_MODIS:",default=0,   RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_GetResource(MAPL,USE_SATSIM_RADAR,LABEL="USE_SATSIM_RADAR:",default=0,   RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_GetResource(MAPL,USE_SATSIM_LIDAR,LABEL="USE_SATSIM_LIDAR:",default=0,   RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_GetResource(MAPL,USE_SATSIM_MISR,LABEL="USE_SATSIM_MISR:",default=0,   RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_GetResource(MAPL,GRIDNAME,'AGCM_GRIDNAME:', RC=STATUS)
    VERIFY_(STATUS)
    GRIDNAME =  AdjustL(GRIDNAME)
    nn = len_trim(GRIDNAME)
    dateline = GRIDNAME(nn-1:nn)
    imchar = GRIDNAME(3:index(GRIDNAME,'x')-1)
    read(imchar,*) imsize
    if(dateline.eq.'CF') imsize = imsize*4
    associate (default => MIN(30,MAX(1,INT(4*5760*4/imsize))) )
      call MAPL_GetResource(MAPL,Ncolumns,LABEL="SATSIM_NCOLUMNS:",default=default,  RC=STATUS)
      VERIFY_(STATUS)
    end associate

    call MAPL_GetResource(MAPL,Npoints_it,LABEL="SATSIM_POINTS_PER_ITERATION:",default=-999,   RC=STATUS)
    VERIFY_(STATUS)

    allocate(         frac_out(IM*JM,NCOLUMNS,LM), __STAT__)
    !allocate(      frac_outinv(IM*JM,NCOLUMNS,LM), __STAT__)
    allocate(   lidar_beta_tot(IM*JM,NCOLUMNS,LM), __STAT__)
    allocate(    lidar_tau_tot(IM*JM,NCOLUMNS,LM), __STAT__)
    allocate(     radar_ze_tot(IM*JM,NCOLUMNS,LM), __STAT__)
    allocate( radar_ze_tot_tmp(IM*JM,NCOLUMNS,LM), __STAT__)
    allocate(radar_ze_tot_mask(IM*JM,NCOLUMNS,LM), __STAT__)


! Pointers to imports
!--------------------
      call MAPL_GetPointer(IMPORT, PLE,  'PLE'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT, ZLE,  'ZLE'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT, T,    'T'     , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT, QV,   'QV'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT, RDFL, 'RL'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT, RDFI, 'RI'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT, RDFR, 'RR'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT, RDFS, 'RS'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT, FCLD, 'FCLD'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT, QLLS, 'QLLS', RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT, QILS, 'QILS', RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT, QLCN, 'QLCN', RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT, QICN, 'QICN', RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT, QRTOT , 'QRTOT', RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT, QSTOT , 'QSTOT', RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT, MCOSZ, 'MCOSZ', RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT, FRLAND,'FRLAND', RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT, FROCEAN,'FROCEAN', RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT, TS    ,'TS', RC=STATUS); VERIFY_(STATUS)

! Pointers to Exports
!--------------------

! ISCCP/Icarus

      call MAPL_GetPointer(EXPORT, ISCCP_CU_OA,  'ISCCP_CU_OA'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, ISCCP_CU_UA,  'ISCCP_CU_UA'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, ISCCP_CU_OB,  'ISCCP_CU_OB'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, ISCCP_CU_UB,  'ISCCP_CU_UB'    , RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(EXPORT, ISCCP_STCU_OA, 'ISCCP_STCU_OA' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, ISCCP_STCU_UA, 'ISCCP_STCU_UA' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, ISCCP_STCU_OB, 'ISCCP_STCU_OB' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, ISCCP_STCU_UB, 'ISCCP_STCU_UB' , RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(EXPORT, ISCCP_ST_OA,   'ISCCP_ST_OA' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, ISCCP_ST_UA,   'ISCCP_ST_UA' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, ISCCP_ST_OB,   'ISCCP_ST_OB' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, ISCCP_ST_UB,   'ISCCP_ST_UB' , RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(EXPORT, ISCCP_ACU_OA,  'ISCCP_ACU_OA' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, ISCCP_ACU_UA,  'ISCCP_ACU_UA' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, ISCCP_ACU_OB,  'ISCCP_ACU_OB' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, ISCCP_ACU_UB,  'ISCCP_ACU_UB' , RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(EXPORT, ISCCP_AST_OA,  'ISCCP_AST_OA' , RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(EXPORT, ISCCP_AST_UA,  'ISCCP_AST_UA' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, ISCCP_AST_OB,  'ISCCP_AST_OB' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, ISCCP_AST_UB,  'ISCCP_AST_UB' , RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(EXPORT, ISCCP_NST_OA,  'ISCCP_NST_OA' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, ISCCP_NST_UA,  'ISCCP_NST_UA' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, ISCCP_NST_OB,  'ISCCP_NST_OB' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, ISCCP_NST_UB,  'ISCCP_NST_UB' , RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(EXPORT, ISCCP_CI_OA,   'ISCCP_CI_OA' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, ISCCP_CI_MA,   'ISCCP_CI_MA' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, ISCCP_CI_UA,   'ISCCP_CI_UA' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, ISCCP_CI_OB,   'ISCCP_CI_OB' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, ISCCP_CI_MB,   'ISCCP_CI_MB' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, ISCCP_CI_UB,   'ISCCP_CI_UB' , RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(EXPORT, ISCCP_CIST_OA, 'ISCCP_CIST_OA' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, ISCCP_CIST_MA, 'ISCCP_CIST_MA' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, ISCCP_CIST_UA, 'ISCCP_CIST_UA' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, ISCCP_CIST_OB, 'ISCCP_CIST_OB' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, ISCCP_CIST_MB, 'ISCCP_CIST_MB' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, ISCCP_CIST_UB, 'ISCCP_CIST_UB' , RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(EXPORT, ISCCP_CB_OA,   'ISCCP_CB_OA' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, ISCCP_CB_MA,   'ISCCP_CB_MA' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, ISCCP_CB_UA,   'ISCCP_CB_UA' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, ISCCP_CB_OB,   'ISCCP_CB_OB' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, ISCCP_CB_MB,   'ISCCP_CB_MB' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, ISCCP_CB_UB,   'ISCCP_CB_UB' , RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(EXPORT, ISCCP_SUBV1,   'ISCCP_SUBV1' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, ISCCP_SUBV2,   'ISCCP_SUBV2' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, ISCCP_SUBV3,   'ISCCP_SUBV3' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, ISCCP_SUBV4,   'ISCCP_SUBV4' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, ISCCP_SUBV5,   'ISCCP_SUBV5' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, ISCCP_SUBV6,   'ISCCP_SUBV6' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, ISCCP_SUBV7,   'ISCCP_SUBV7' , RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(EXPORT, CLISCCP1,   'CLISCCP1' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CLISCCP2,   'CLISCCP2' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CLISCCP3,   'CLISCCP3' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CLISCCP4,   'CLISCCP4' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CLISCCP5,   'CLISCCP5' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CLISCCP6,   'CLISCCP6' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CLISCCP7,   'CLISCCP7' , RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(EXPORT, ISCCP1,   'ISCCP1' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, ISCCP2,   'ISCCP2' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, ISCCP3,   'ISCCP3' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, ISCCP4,   'ISCCP4' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, ISCCP5,   'ISCCP5' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, ISCCP6,   'ISCCP6' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, ISCCP7,   'ISCCP7' , RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(EXPORT, TCLISCCP, 'TCLISCCP' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CLTISCCP, 'CLTISCCP' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, PCTISCCP, 'PCTISCCP' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CTPISCCP, 'CTPISCCP' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, ALBISCCP, 'ALBISCCP' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, TBISCCP, 'TBISCCP' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, SGFCLD, 'SGFCLD' , RC=STATUS); VERIFY_(STATUS)

! LIDAR/CALIPSO

      call MAPL_GetPointer(EXPORT, RADARLTCC, 'RADARLTCC' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CLCALIPSO2, 'CLCALIPSO2' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, LIDARPMOL, 'LIDARPMOL' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, LIDARPTOT, 'LIDARPTOT' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, LIDARTAUTOT, 'LIDARTAUTOT' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, RADARZETOT, 'RADARZETOT' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CLCALIPSO, 'CLCALIPSO' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CLLCALIPSO, 'CLLCALIPSO' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CLMCALIPSO, 'CLMCALIPSO' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CLHCALIPSO, 'CLHCALIPSO' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CLTCALIPSO, 'CLTCALIPSO' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, PARASOLREFL0, 'PARASOLREFL0' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, PARASOLREFL1, 'PARASOLREFL1' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, PARASOLREFL2, 'PARASOLREFL2' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, PARASOLREFL3, 'PARASOLREFL3' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, PARASOLREFL4, 'PARASOLREFL4' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, PARASOLREFL5, 'PARASOLREFL5' , RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(EXPORT, CFADLIDARSR532_01,'CFADLIDARSR532_01' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CFADLIDARSR532_02,'CFADLIDARSR532_02' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CFADLIDARSR532_03,'CFADLIDARSR532_03' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CFADLIDARSR532_04,'CFADLIDARSR532_04' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CFADLIDARSR532_05,'CFADLIDARSR532_05' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CFADLIDARSR532_06,'CFADLIDARSR532_06' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CFADLIDARSR532_07,'CFADLIDARSR532_07' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CFADLIDARSR532_08,'CFADLIDARSR532_08' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CFADLIDARSR532_09,'CFADLIDARSR532_09' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CFADLIDARSR532_10,'CFADLIDARSR532_10' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CFADLIDARSR532_11,'CFADLIDARSR532_11' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CFADLIDARSR532_12,'CFADLIDARSR532_12' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CFADLIDARSR532_13,'CFADLIDARSR532_13' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CFADLIDARSR532_14,'CFADLIDARSR532_14' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CFADLIDARSR532_15,'CFADLIDARSR532_15' , RC=STATUS); VERIFY_(STATUS)

! RADAR/Cloudsat

      call MAPL_GetPointer(EXPORT, CLOUDSATCFAD01,'CLOUDSATCFAD01' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CLOUDSATCFAD02,'CLOUDSATCFAD02' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CLOUDSATCFAD03,'CLOUDSATCFAD03' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CLOUDSATCFAD04,'CLOUDSATCFAD04' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CLOUDSATCFAD05,'CLOUDSATCFAD05' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CLOUDSATCFAD06,'CLOUDSATCFAD06' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CLOUDSATCFAD07,'CLOUDSATCFAD07' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CLOUDSATCFAD08,'CLOUDSATCFAD08' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CLOUDSATCFAD09,'CLOUDSATCFAD09' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CLOUDSATCFAD10,'CLOUDSATCFAD10' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CLOUDSATCFAD11,'CLOUDSATCFAD11' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CLOUDSATCFAD12,'CLOUDSATCFAD12' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CLOUDSATCFAD13,'CLOUDSATCFAD13' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CLOUDSATCFAD14,'CLOUDSATCFAD14' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CLOUDSATCFAD15,'CLOUDSATCFAD15' , RC=STATUS); VERIFY_(STATUS)

! MODIS

      call MAPL_GetPointer(EXPORT, MDSCLDFRCTTL,'MDSCLDFRCTTL' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MDSCLDFRCWTR,'MDSCLDFRCWTR' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MDSCLDFRCH2O,'MDSCLDFRCH2O' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MDSCLDFRCICE,'MDSCLDFRCICE' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MDSCLDFRCHI,'MDSCLDFRCHI' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MDSCLDFRCMID,'MDSCLDFRCMID' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MDSCLDFRCLO,'MDSCLDFRCLO' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MDSOPTHCKTTL,'MDSOPTHCKTTL' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MDSOPTHCKWTR,'MDSOPTHCKWTR' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MDSOPTHCKH2O,'MDSOPTHCKH2O' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MDSOPTHCKICE,'MDSOPTHCKICE' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MDSOPTHCKTTLLG,'MDSOPTHCKTTLLG' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MDSOPTHCKWTRLG,'MDSOPTHCKWTRLG' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MDSOPTHCKH2OLG,'MDSOPTHCKH2OLG' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MDSOPTHCKICELG,'MDSOPTHCKICELG' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MDSCLDSZWTR,'MDSCLDSZWTR' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MDSCLDSZH20,'MDSCLDSZH20' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MDSCLDSZICE,'MDSCLDSZICE' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MDSCLDTOPPS,'MDSCLDTOPPS' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MDSWTRPATH,'MDSWTRPATH' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MDSH2OPATH,'MDSH2OPATH' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MDSICEPATH,'MDSICEPATH' , RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST11,'MDSTAUPRSHIST11' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST12,'MDSTAUPRSHIST12' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST13,'MDSTAUPRSHIST13' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST14,'MDSTAUPRSHIST14' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST15,'MDSTAUPRSHIST15' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST16,'MDSTAUPRSHIST16' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST17,'MDSTAUPRSHIST17' , RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST21,'MDSTAUPRSHIST21' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST22,'MDSTAUPRSHIST22' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST23,'MDSTAUPRSHIST23' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST24,'MDSTAUPRSHIST24' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST25,'MDSTAUPRSHIST25' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST26,'MDSTAUPRSHIST26' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST27,'MDSTAUPRSHIST27' , RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST31,'MDSTAUPRSHIST31' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST32,'MDSTAUPRSHIST32' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST33,'MDSTAUPRSHIST33' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST34,'MDSTAUPRSHIST34' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST35,'MDSTAUPRSHIST35' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST36,'MDSTAUPRSHIST36' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST37,'MDSTAUPRSHIST37' , RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST41,'MDSTAUPRSHIST41' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST42,'MDSTAUPRSHIST42' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST43,'MDSTAUPRSHIST43' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST44,'MDSTAUPRSHIST44' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST45,'MDSTAUPRSHIST45' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST46,'MDSTAUPRSHIST46' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST47,'MDSTAUPRSHIST47' , RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST51,'MDSTAUPRSHIST51' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST52,'MDSTAUPRSHIST52' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST53,'MDSTAUPRSHIST53' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST54,'MDSTAUPRSHIST54' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST55,'MDSTAUPRSHIST55' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST56,'MDSTAUPRSHIST56' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST57,'MDSTAUPRSHIST57' , RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST61,'MDSTAUPRSHIST61' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST62,'MDSTAUPRSHIST62' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST63,'MDSTAUPRSHIST63' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST64,'MDSTAUPRSHIST64' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST65,'MDSTAUPRSHIST65' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST66,'MDSTAUPRSHIST66' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST67,'MDSTAUPRSHIST67' , RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST71,'MDSTAUPRSHIST71' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST72,'MDSTAUPRSHIST72' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST73,'MDSTAUPRSHIST73' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST74,'MDSTAUPRSHIST74' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST75,'MDSTAUPRSHIST75' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST76,'MDSTAUPRSHIST76' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST77,'MDSTAUPRSHIST77' , RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(EXPORT, MISRMNCLDTP,'MISRMNCLDTP' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MISRCLDAREA,'MISRCLDAREA' , RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(EXPORT, MISRLYRTP0,'MISRLYRTP0' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MISRLYRTP250,'MISRLYRTP250' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MISRLYRTP750,'MISRLYRTP750' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MISRLYRTP1250,'MISRLYRTP1250' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MISRLYRTP1750,'MISRLYRTP1750' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MISRLYRTP2250,'MISRLYRTP2250' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MISRLYRTP2750,'MISRLYRTP2750' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MISRLYRTP3500,'MISRLYRTP3500' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MISRLYRTP4500,'MISRLYRTP4500' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MISRLYRTP6000,'MISRLYRTP6000' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MISRLYRTP8000,'MISRLYRTP8000' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MISRLYRTP10000,'MISRLYRTP10000' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MISRLYRTP12000,'MISRLYRTP12000' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MISRLYRTP14000,'MISRLYRTP14000' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MISRLYRTP16000,'MISRLYRTP16000' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MISRLYRTP18000,'MISRLYRTP18000' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MISRFQ0,'MISRFQ0' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MISRFQ250,'MISRFQ250' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MISRFQ750,'MISRFQ750' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MISRFQ1250,'MISRFQ1250' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MISRFQ1750,'MISRFQ1750' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MISRFQ2250,'MISRFQ2250' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MISRFQ2750,'MISRFQ2750' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MISRFQ3500,'MISRFQ3500' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MISRFQ4500,'MISRFQ4500' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MISRFQ6000,'MISRFQ6000' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MISRFQ8000,'MISRFQ8000' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MISRFQ10000,'MISRFQ10000' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MISRFQ12000,'MISRFQ12000' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MISRFQ14000,'MISRFQ14000' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MISRFQ16000,'MISRFQ16000' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MISRFQ18000,'MISRFQ18000' , RC=STATUS); VERIFY_(STATUS)



!  The follwing assignments are for variables needed for COSP and isccp/icarus
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ISCCP_TOP_HEIGHT=1           !  1 = adjust top height using both a computed
   ISCCP_TOP_HEIGHT_DIRECTION=2 ! direction for finding atmosphere pressure level
  isccp_emsfc_lw=0.999              
   sunlit        = 0.0
   ICB=58  
   ICT=48

! anything that needs to be calculated in 3D is done so, to
! be converted to 2D and inverted later
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

     PLO  = 0.5*( PLE(:,:,1:LM)+PLE(:,:,0:lm-1))

! The following are convective terms, which we don't want now
      CCA           	= 0.0

      QLTOT = QLLS +  QLCN
      QITOT = QILS +  QICN

! the following is in lieu of 
!      where(FCLD.gt.0.)
! to avoid stupid high cwc, which makes for stupid high dtau_s, which
! crashes icarus
      where(FCLD.gt.0.01)
      cwc(:,:,:,1) = MAX( QITOT/FCLD , 1.0e-12 )
      cwc(:,:,:,2) = MAX( QLTOT/FCLD , 1.0e-12 )
      cwc(:,:,:,3) = MAX( QRTOT/FCLD , 1.0e-12 )
      cwc(:,:,:,4) = MAX( QSTOT/FCLD , 1.0e-12 )
      elsewhere
      cwc(:,:,:,1) = 0.
      cwc(:,:,:,2) = 0.
      cwc(:,:,:,3) = 0.
      cwc(:,:,:,4) = 0.
      endwhere

! delp in Pascals, reff in microns

      reff(:,:,:,1) = RDFI *1.e6
      reff(:,:,:,2) = RDFL *1.e6
      reff(:,:,:,3) = RDFR *1.e6
      reff(:,:,:,4) = RDFS *1.e6

      DELP = PLE(:,:,1:LM) - PLE(:,:,0:lm-1)


      do i = 1, im
         do j = 1, jm

            ! NOTE: mcosz is being passed in but not used since we aren't doing cloud scaling 
            !       as the 0's at the end of inputs determine. If scaling is needed, mcosz might
            !       be the wrong cosine to use (i.e., not the mean). Also, dummies will need to
            !       be replaced if needed.
            call getvistau(lm,mcosz(i,j),delp(i,j,:),FCLD(i,j,:),reff(i,j,:,:),cwc(i,j,:,:),0,0,&
                           dumtaubeam(:,:),tausw(i,j,:,:),dumasycl(:))

            ! IR Taus: Here we calculate band 4 only. Again, dummies for unneeded arrays.
            call getirtau(4,lm,delp(i,j,:),FCLD(i,j,:),reff(i,j,:,:),cwc(i,j,:,:),&
                          taulw(i,j,:,:),dumtcldlyr(:),dumenn(:))

            do k = 1, lm

               ! NOTE: If you wish to include both falling rain and falling snow to 
               !       DTAU_S and EMISS, use the lines below with all four
               !       tausw and taulw returns

               ! VIS
               ! ---

               !DTAU_S(i,j,k) = tausw(i,j,k,1)+tausw(i,j,k,2)+tausw(i,j,k,3)+tausw(i,j,k,4)
               DTAU_S(i,j,k) = tausw(i,j,k,1)+tausw(i,j,k,2)

               ! IR
               ! --

               !taucir = taulw(i,j,k,1) + taulw(i,j,k,2) + taulw(i,j,k,3) + taulw(i,j,k,4)
               taucir = taulw(i,j,k,1) + taulw(i,j,k,2)

               EMISS(i,j,k) = 1 - EXP( -1.0 * taucir )

            end do

         end do
      end do

!  The follwing assignments are for variables needed for all simulators
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   NPOINTS   = IM*JM 
   ISCCP_OVERLAP  =3            !  overlap type: 1=max, 2=rand, 3=max/rand
				! used by SCOPS, COSP and icarus


!  The follwing assignments are for variables needed for COSP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   Naero=1              ! Number of aerosol species (Not used)
   if (Npoints_it == -999)  Npoints_it=NPOINTS         ! Number of gridpoints processed in one iteration (probably not used)
   use_precipitation_fluxes=.false.    
   use_reff=.true.                    
   SURFACE_RADAR=0   ! surface=1, spaceborne=0
   melt_lay=0        ! melting layer model off=0, on=1
   time=0			! Time since start of run [days] - something should be done about this

! Switches in cfg need to be toggled to make things run right

! remember these are integers so + is logical or, not and

    if ( USE_SATSIM + USE_SATSIM_ISCCP + USE_SATSIM_MODIS > 0 ) then ! isccp needed for modis
      cfg%Lisccp_sim = .true.
    else 
      cfg%Lisccp_sim = .false.
    endif 

    if ( USE_SATSIM + USE_SATSIM_MODIS > 0 ) then
      cfg%Lmodis_sim= .true.
    else
      cfg%Lmodis_sim= .false.
    endif 

    if ( USE_SATSIM + USE_SATSIM_RADAR > 0 ) then
      cfg%Lradar_sim= .true.
    else
      cfg%Lradar_sim= .false.
    endif

    if ( USE_SATSIM + USE_SATSIM_LIDAR > 0 ) then
      cfg%Llidar_sim = .true.
      cfg%LCFADLIDARSR532 = .true.
    else
      cfg%Llidar_sim = .false.
      cfg%LCFADLIDARSR532 = .false.
    endif

    if ( USE_SATSIM + USE_SATSIM_MISR > 0 ) then
      cfg%Lmisr_sim= .true.
    else
      cfg%Lmisr_sim= .false.
    endif

    if ( USE_SATSIM /= 0 .or. ( USE_SATSIM_LIDAR + USE_SATSIM_RADAR .eq. 2 ) ) then
      cfg%Lstats= .true.
    else
      cfg%Lstats= .false.
    endif 

! RRTOV paramters (not used)
   plat=0                        
   sat=0                         
   inst=0                        
   nchan=0                       
   ZenAng=0.0                      
   Ichan=0              
   surfem=0.0             
   co2=0.0                         
   ch4=0.0                         
   n2o=0.0                         
   co=0.0                          



!  The follwing assignments are for variables needed for COSP and lidar_simulator
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   lidar_ice_type=0     ! Ice particle shape in lidar calculations (0=ice-spheres ; 1=ice-non-spherical)

!  The follwing assignments are for variables needed for COSP and radar_simulator
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   Nprmts_max_aero=1    ! Max number of parameters for aerosol size distributions (Not used)
   Nprmts_max_hydro=12  ! Max number of parameters for hydrometeor size distributions
   RADAR_FREQ=94.0   ! CloudSat radar frequency (GHz)
   k2=-1.            ! |K|^2, -1=use frequency dependent default
   do_ray=0          ! calculate/output Rayleigh refl=1, not=0
   use_gas_abs=1     ! include gaseous absorption? yes=1,no=0
   use_mie_tables=0  ! use a precomputed lookup table? yes=1,no=0
   QST  = GEOS_QSAT( T ,  PLO/100. )
   RH      = ( QV / QST ) * 100.0 ! rel. humidity (in percent)
   land_mask=0.0 !  needed for stats -- should be changed [0 - Ocean, 1 - Land]  

      RHCOSP = reshape( RH(:,:,LM:1:-1) , (/ IM*JM , LM /) )
      PLOCOSP  = reshape( PLO(:,:,LM:1:-1), (/ IM*JM , LM /) )
      QVCOSP  = reshape( QV(:,:,LM:1:-1), (/ IM*JM , LM /) )
      CCACOSP = reshape( CCA(:,:,LM:1:-1), (/ IM*JM , LM /) )
      DTAU_SCOSP = reshape( DTAU_S(:,:,LM:1:-1), (/ IM*JM , LM /) )
      TCOSP = reshape( T(:,:,LM:1:-1), (/ IM*JM , LM /) )
      EMISSCOSP = reshape( EMISS(:,:,LM:1:-1), (/ IM*JM , LM /) )
      FCLDCOSP  = reshape( FCLD(:,:,LM:1:-1), (/ IM*JM , LM /) )
      RDFLCOSP  = reshape( RDFL(:,:,LM:1:-1), (/ IM*JM , LM /) )
      RDFICOSP  = reshape( RDFI(:,:,LM:1:-1), (/ IM*JM , LM /) )
      RDFRCOSP  = reshape( RDFR(:,:,LM:1:-1), (/ IM*JM , LM /) )
      RDFSCOSP  = reshape( RDFS(:,:,LM:1:-1), (/ IM*JM , LM /) )
      QLTOTCOSP  = reshape( QLTOT(:,:,LM:1:-1), (/ IM*JM , LM /) )
      QITOTCOSP  = reshape( QITOT(:,:,LM:1:-1), (/ IM*JM , LM /) )
      QRTOTCOSP  = reshape( QRTOT(:,:,LM:1:-1), (/ IM*JM , LM /) )
      QSTOTCOSP  = reshape( QSTOT(:,:,LM:1:-1), (/ IM*JM , LM /) )
      QLCNCOSP  = reshape( QLCN(:,:,LM:1:-1), (/ IM*JM , LM /) )
      QICNCOSP  = reshape( QICN(:,:,LM:1:-1), (/ IM*JM , LM /) )
      PLECOSP  = reshape( PLE(:,:,LM:0:-1), (/ IM*JM , LM+1 /) )
      ZLECOSP = ZLE2D(:,LM:0:-1)
      MCOSZCOSP = reshape( MCOSZ , (/ IM*JM /) )
      FRLANDCOSP = reshape( FRLAND , (/ IM*JM /) )
      FROCEANCOSP = reshape( FROCEAN , (/ IM*JM /) )
      TSCOSP = reshape( TS , (/ IM*JM /) )

      ZLE2D  = reshape( ZLE, (/ IM*JM , LM+1 /) )
      zlo2d =  0.5*( zle2d(:,0:LM-1)+zle2d(:,1:lm))
      ZLOCOSP = ZLO2D(:,LM:1:-1)


   where ( MCOSZCOSP > 0.0 ) SUNLIT = 1.0
   where ( FROCEANCOSP < 0.5 ) LAND_MASK = 1.0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! zero out exports
      ISCCP_totalcldarea = 0.0
      ISCCP_meanptop     = 0.0
      ISCCP_meanalbedocld= 0.0
      ISCCP_meanallskybrighttemp= 0.0
      frac_out     = 0.0
      fq_isccp     = 0.0
      lidar_lidarcld     = 0.0
      lidar_cldlayer     = 0.0
      lidar_parasolrefl  = 0.0
      lidar_cfad_sr      = 0.0
      radar_lidar_tcc = 0.0
      radar_lidar_only_freq_cloud = 0.0
      radar_cfad_ze = 0.0 
      radar_ze_tot = 0.0 

   lidar_pmol=0.0
   REFL=0.0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 
  call construct_cosp_gridbox( time                        &	
                             , time_bands                  &
                             , radar_freq                  &
                             , surface_radar               &
                             , use_mie_tables              &
                             , use_gas_abs                 &
                             , do_ray                      &
                             , melt_lay                    &
                             , k2                          &
                             , Npoints                     & 
                             , LM                          &
                             , Ncolumns                    &
                             , N_HYDRO                     &
                             , Nprmts_max_hydro            &
                             , Naero                       &
                             , Nprmts_max_aero             &
                             , Npoints_it                  &
                             , lidar_ice_type              &
                             , isccp_top_height            &
                             , isccp_top_height_direction  &
                             , isccp_overlap               &
                             , isccp_emsfc_lw              &
                             , use_precipitation_fluxes    &
                             , use_reff                    &
                             , plat                        &
                             , sat                         &
                             , inst                        &
                             , nchan                       &
                             , ZenAng                      &
                             , Ichan                       &
                             , surfem                      &
                             , co2                         &
                             , ch4                         &
                             , n2o                         &
                             , co                          &
                             , gbx                          )

   call construct_cosp_subgrid(Npoints, Ncolumns, LM, sgx)
   if (cfg%Lisccp_sim) call construct_cosp_isccp(cfg,Npoints,Ncolumns,LM,isccp)
   !call construct_cosp_sghydro(Npoints,Ncolumns,LM,N_hydro,sghydro)
   if (cfg%Lradar_sim) call construct_cosp_sgradar(cfg,Npoints,Ncolumns,LM,N_hydro,sgradar)
   if (cfg%Llidar_sim) call construct_cosp_sglidar(cfg,Npoints,Ncolumns,LM,N_hydro,PARASOL_NREFL,sglidar)
   if (cfg%Lmodis_sim) call construct_cosp_modis(cfg, Npoints, modis )
   if (cfg%Lmisr_sim ) call construct_cosp_misr(cfg, Npoints, misr )

   !   vgrid is use for lidar and radar simulation stats, but only if the stats are on 
   ! a different grid from the model.  So mostly dummy params get passed`
   call construct_cosp_vgrid(gbx,0,.false.,.false.,vgrid)

   if (cfg%Lradar_sim) call construct_cosp_radarstats(cfg,Npoints,Ncolumns,LM,N_hydro,stradar)
   if (cfg%Llidar_sim) call construct_cosp_lidarstats(cfg,Npoints,Ncolumns,LM,N_hydro,PARASOL_NREFL,stlidar)
      
   BEGSEG = 1
   ENDSEG = NPOINTS

   gbx%T   =  TCOSP(BEGSEG:ENDSEG,:)
   gbx%SH  =  QVCOSP(BEGSEG:ENDSEG,:)
   gbx%Q   =  RHCOSP(BEGSEG:ENDSEG,:)
   gbx%CCA =   CCACOSP(BEGSEG:ENDSEG,:)
   gbx%TCA =  FCLDCOSP(BEGSEG:ENDSEG,:)
! in cosp_types, gbx%ph is supposed to be half pressure levels, with the bottom as surface
!   pressure -- see psfc in cosp_test  
   gbx%ph =  PLECOSP(BEGSEG:ENDSEG,1:LM)
   gbx%p  =  PLOCOSP(BEGSEG:ENDSEG,:)
   gbx%psfc = PLECOSP(BEGSEG:ENDSEG,0) ! surface pressure
   gbx%dtau_s   = dtau_sCOSP(BEGSEG:ENDSEG,:)
   gbx%dtau_c   = 0.
   gbx%dem_s    = emissCOSP(BEGSEG:ENDSEG,:)
   gbx%dem_c    = 0. 
   gbx%sunlit=sunlit(BEGSEG:ENDSEG)
   gbx%land=LAND_MASK(BEGSEG:ENDSEG)
   gbx%skt = TSCOSP(BEGSEG:ENDSEG)
   gbx%mr_hydro(:,:,I_LSCLIQ) =  QLTOTCOSP(BEGSEG:ENDSEG,:)
   gbx%mr_hydro(:,:,I_LSCICE) =  QITOTCOSP(BEGSEG:ENDSEG,:)
   gbx%mr_hydro(:,:,I_LSRAIN) =  QRTOTCOSP(BEGSEG:ENDSEG,:)
   gbx%mr_hydro(:,:,I_LSSNOW) =  QSTOTCOSP(BEGSEG:ENDSEG,:)
   gbx%mr_hydro(:,:,I_CVCLIQ) = 0.0
   gbx%mr_hydro(:,:,I_CVCICE) = 0.0
   gbx%reff(:,:,I_LSCLIQ) =   RDFLCOSP(BEGSEG:ENDSEG,:)
   gbx%reff(:,:,I_LSCICE) =   RDFICOSP(BEGSEG:ENDSEG,:)
   gbx%reff(:,:,I_LSRAIN) =   RDFRCOSP(BEGSEG:ENDSEG,:)
   gbx%reff(:,:,I_LSSNOW) =   RDFSCOSP(BEGSEG:ENDSEG,:)
   gbx%zlev =  ZLOCOSP(BEGSEG:ENDSEG,:)
   gbx%zlev_half =  ZLECOSP(BEGSEG:ENDSEG,1:LM)


   where (gbx%mr_hydro < 0.0) gbx%mr_hydro=0.0
   
if ( MAPL_Am_I_Root() .and. DEBUG_GC ) then

  write(*,*) 'ncolumns: ', ncolumns
  write(*,*) 'lm: ', lm
  write(*,*) 'npoints: ', npoints
  write(*,*) 'npoints_it: ', npoints_it
  
  !write(*,*) 'FCLD: ', FCLD

  !write(*,*) 'gbx%dtau_s: ', gbx%dtau_s
  !write(*,*) 'gbx%dem_s: ', gbx%dem_s

  !write(*,*) 'gbx%tca: ', gbx%tca
  !write(*,*) 'sgx%frac_out: ', sgx%frac_out

endif


    call MAPL_TimerOn(STATE,"-COSP")

   call COSP(isccp_overlap,Ncolumns,cfg,vgrid,gbx,sgx,sgradar,sglidar,isccp,misr,modis,stradar,stlidar)

    call MAPL_TimerOff(STATE,"-COSP")

! Dims of sgx%frac_out (Npoints,Ncolumns,Nlevels)
   if ( associated( SGFCLD ) ) then
      sgfcld = reshape(                      &
                   sum( sgx%frac_out(:,:,lm:1:-1) , 2)   &
               ,      (/ IM, JM, LM /)   )   &
                            / (1.*NCOLUMNS)
   endif


! more outputs are available in the structures -- the current ones are those
! needed for CFMIP

   if ( USE_SATSIM_LIDAR + USE_SATSIM > 0 ) then

      lidar_parasolrefl(BEGSEG:ENDSEG,:)=stlidar%parasolrefl
      lidar_cldlayer(BEGSEG:ENDSEG,:)=stlidar%cldlayer
      lidar_lidarcld(BEGSEG:ENDSEG,:)=stlidar%lidarcld
      lidar_cfad_sr(BEGSEG:ENDSEG,:,:)=stlidar%cfad_sr

! used for diagnostics
      lidar_pmol(BEGSEG:ENDSEG,:) = sglidar%beta_mol
      lidar_beta_tot(BEGSEG:ENDSEG,:,:) = sglidar%beta_tot
      lidar_tau_tot(BEGSEG:ENDSEG,:,:) = sglidar%tau_tot

   endif

   if ( USE_SATSIM_RADAR + USE_SATSIM > 0 ) then
      radar_ze_tot(BEGSEG:ENDSEG,:,:)  = sgradar%ze_tot
      radar_lidar_only_freq_cloud(BEGSEG:ENDSEG,:)=stradar%lidar_only_freq_cloud
      radar_lidar_tcc(BEGSEG:ENDSEG)=stradar%radar_lidar_tcc
      radar_cfad_ze(BEGSEG:ENDSEG,:,:)=stradar%cfad_ze
   endif


   if ( USE_SATSIM + USE_SATSIM_ISCCP > 0 ) then 
! cosp reverses the pressure levls of the isccp matrix for some reason
      fq_isccp(BEGSEG:ENDSEG,:,:) = isccp%fq_isccp(:,:,7:1:-1) 
      ISCCP_totalcldarea(BEGSEG:ENDSEG) = isccp%totalcldarea
      ISCCP_meanptop(BEGSEG:ENDSEG) = isccp%meanptop 
      ISCCP_meanalbedocld(BEGSEG:ENDSEG) = isccp%meanalbedocld
      ISCCP_meanallskybrighttemp(BEGSEG:ENDSEG) = isccp%meantb
      fq_isccp3D=reshape( fq_isccp, (/ IM, JM, 7 , 7 /) )

#ifdef USE_MAPL_UNDEF

      where (fq_isccp3D < -10.0) fq_isccp3D=MAPL_UNDEF
      where (fq_isccp3D > 10.0) fq_isccp3D=MAPL_UNDEF
      where (ISCCP_totalcldarea < 0.0) ISCCP_totalcldarea=MAPL_UNDEF
      where (ISCCP_totalcldarea > 1.0) ISCCP_totalcldarea=MAPL_UNDEF
      where (ISCCP_meanptop .le. 0.0 ) ISCCP_meanptop=MAPL_UNDEF
      where (ISCCP_meanalbedocld < 0.0) ISCCP_meanalbedocld=MAPL_UNDEF
      where (ISCCP_meanalbedocld > 1.0) ISCCP_meanalbedocld=MAPL_UNDEF
      where (ISCCP_meanallskybrighttemp < 0.0) ISCCP_meanallskybrighttemp=MAPL_UNDEF

#endif 

   endif  ! ( USE_SATSIM + USE_SATSIM_ISCCP > 0 )

   if ( USE_SATSIM_MISR + USE_SATSIM > 0 ) then
      MISR_meanztop(BEGSEG:ENDSEG) = misr%MISR_meanztop
      fq_MISR(BEGSEG:ENDSEG,:,:) = misr%fq_MISR
      MISR_cldarea(BEGSEG:ENDSEG) = misr%MISR_cldarea
      MISR_dist_model_layertops(BEGSEG:ENDSEG,:) = misr%MISR_dist_model_layertops
   endif

   if ( USE_SATSIM_MODIS + USE_SATSIM > 0 ) then
      MODIS_Cloud_Fraction_Total_Mean(BEGSEG:ENDSEG) = modis%Cloud_Fraction_Total_Mean
      MODIS_Cloud_Fraction_Water_Mean(BEGSEG:ENDSEG) = modis%Cloud_Fraction_Water_Mean
      MODIS_Cloud_Fraction_Ice_Mean(BEGSEG:ENDSEG) = modis%Cloud_Fraction_Ice_Mean 
      MODIS_Cloud_Fraction_High_Mean(BEGSEG:ENDSEG) = modis%Cloud_Fraction_High_Mean 
      MODIS_Cloud_Fraction_Mid_Mean(BEGSEG:ENDSEG) = modis%Cloud_Fraction_Mid_Mean 
      MODIS_Cloud_Fraction_Low_Mean(BEGSEG:ENDSEG) = modis%Cloud_Fraction_Low_Mean 
      MODIS_Optical_Thickness_Total_Mean(BEGSEG:ENDSEG) = modis%Optical_Thickness_Total_Mean 
      MODIS_Optical_Thickness_Water_Mean(BEGSEG:ENDSEG) = modis%Optical_Thickness_Water_Mean
      MODIS_Optical_Thickness_Ice_Mean(BEGSEG:ENDSEG) = modis%Optical_Thickness_Ice_Mean
      MODIS_Optical_Thickness_Total_LogMean(BEGSEG:ENDSEG) = modis%Optical_Thickness_Total_LogMean
      MODIS_Optical_Thickness_Water_LogMean(BEGSEG:ENDSEG) = modis%Optical_Thickness_Water_LogMean
      MODIS_Optical_Thickness_Ice_LogMean(BEGSEG:ENDSEG) = modis%Optical_Thickness_Ice_LogMean
      MODIS_Cloud_Particle_Size_Water_Mean(BEGSEG:ENDSEG) = modis%Cloud_Particle_Size_Water_Mean
      MODIS_Cloud_Particle_Size_Ice_Mean(BEGSEG:ENDSEG) = modis%Cloud_Particle_Size_Ice_Mean
      MODIS_Cloud_Top_Pressure_Total_Mean(BEGSEG:ENDSEG) = modis%Cloud_Top_Pressure_Total_Mean
      MODIS_Liquid_Water_Path_Mean(BEGSEG:ENDSEG) = modis%Liquid_Water_Path_Mean
      MODIS_Ice_Water_Path_Mean(BEGSEG:ENDSEG) = modis%Ice_Water_Path_Mean
      MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(BEGSEG:ENDSEG,:,:) = modis%Optical_Thickness_vs_Cloud_Top_Pressure
   endif 
   
#ifdef USE_MAPL_UNDEF

   if ( USE_SATSIM_MODIS + USE_SATSIM > 0 ) then
! setting to undef since 0 implies no clouds
      where (MODIS_Optical_Thickness_Total_Mean    .le. 0.0 )  MODIS_Optical_Thickness_Total_Mean    = MAPL_UNDEF
      where (modis_Optical_Thickness_Water_Mean    .le. 0.0 )  modis_Optical_Thickness_Water_Mean    = MAPL_UNDEF
      where (modis_Optical_Thickness_Ice_Mean      .le. 0.0 )  modis_Optical_Thickness_Ice_Mean      = MAPL_UNDEF
      where (modis_Optical_Thickness_Total_LogMean .le. 0.0 )  modis_Optical_Thickness_Total_LogMean = MAPL_UNDEF
      where (modis_Optical_Thickness_Water_LogMean .le. 0.0 )  modis_Optical_Thickness_Water_LogMean = MAPL_UNDEF
      where (modis_Optical_Thickness_Ice_LogMean   .le. 0.0 )  modis_Optical_Thickness_Ice_LogMean   = MAPL_UNDEF
      where (modis_Cloud_Particle_Size_Water_Mean  .le. 0.0 )  modis_Cloud_Particle_Size_Water_Mean  = MAPL_UNDEF
      where (modis_Cloud_Particle_Size_Ice_Mean    .le. 0.0 )  modis_Cloud_Particle_Size_Ice_Mean    = MAPL_UNDEF

      where (MODIS_Cloud_Fraction_Total_Mean == R_UNDEF )  MODIS_Cloud_Fraction_Total_Mean = MAPL_UNDEF
      where (MODIS_Cloud_Fraction_Water_Mean == R_UNDEF )  MODIS_Cloud_Fraction_Water_Mean = MAPL_UNDEF
      where (MODIS_Cloud_Fraction_Ice_Mean == R_UNDEF )  MODIS_Cloud_Fraction_Ice_Mean  = MAPL_UNDEF
      where (MODIS_Cloud_Fraction_High_Mean == R_UNDEF )  MODIS_Cloud_Fraction_High_Mean  = MAPL_UNDEF
      where (MODIS_Cloud_Fraction_Mid_Mean == R_UNDEF )  MODIS_Cloud_Fraction_Mid_Mean  = MAPL_UNDEF
      where (MODIS_Cloud_Fraction_Low_Mean == R_UNDEF )  MODIS_Cloud_Fraction_Low_Mean  = MAPL_UNDEF
      where (MODIS_Optical_Thickness_Total_Mean == R_UNDEF )  MODIS_Optical_Thickness_Total_Mean  = MAPL_UNDEF
      where (MODIS_Optical_Thickness_Water_Mean == R_UNDEF )  MODIS_Optical_Thickness_Water_Mean = MAPL_UNDEF
      where (MODIS_Optical_Thickness_Ice_Mean == R_UNDEF )  MODIS_Optical_Thickness_Ice_Mean = MAPL_UNDEF
      where (MODIS_Optical_Thickness_Total_LogMean == R_UNDEF )  MODIS_Optical_Thickness_Total_LogMean = MAPL_UNDEF
      where (MODIS_Optical_Thickness_Water_LogMean == R_UNDEF )  MODIS_Optical_Thickness_Water_LogMean = MAPL_UNDEF
      where (MODIS_Optical_Thickness_Ice_LogMean == R_UNDEF )  MODIS_Optical_Thickness_Ice_LogMean = MAPL_UNDEF
      where (MODIS_Cloud_Particle_Size_Water_Mean == R_UNDEF )  MODIS_Cloud_Particle_Size_Water_Mean = MAPL_UNDEF
      where (MODIS_Cloud_Particle_Size_Ice_Mean == R_UNDEF )  MODIS_Cloud_Particle_Size_Ice_Mean = MAPL_UNDEF
      where (MODIS_Cloud_Top_Pressure_Total_Mean == R_UNDEF )  MODIS_Cloud_Top_Pressure_Total_Mean = MAPL_UNDEF
      where (MODIS_Cloud_Top_Pressure_Total_Mean .le. 0.0 )  MODIS_Cloud_Top_Pressure_Total_Mean = MAPL_UNDEF
      where (MODIS_Liquid_Water_Path_Mean == R_UNDEF )  MODIS_Liquid_Water_Path_Mean = MAPL_UNDEF
      where (MODIS_Ice_Water_Path_Mean == R_UNDEF )  MODIS_Ice_Water_Path_Mean = MAPL_UNDEF
      where (MODIS_Optical_Thickness_vs_Cloud_Top_Pressure == R_UNDEF )  MODIS_Optical_Thickness_vs_Cloud_Top_Pressure = MAPL_UNDEF
   endif

   if ( USE_SATSIM_LIDAR + USE_SATSIM > 0 ) then
      where (lidar_cfad_sr   == R_UNDEF) lidar_cfad_sr   = MAPL_UNDEF
      where (lidar_lidarcld  .gt. 1.0 ) lidar_lidarcld  = MAPL_UNDEF 
      where (lidar_lidarcld  .lt. 0.0 ) lidar_lidarcld  = MAPL_UNDEF 
      where (lidar_cldlayer  == R_UNDEF) lidar_cldlayer  = MAPL_UNDEF
      where (lidar_parasolrefl == R_UNDEF) lidar_parasolrefl = MAPL_UNDEF 
   endif

   if ( USE_SATSIM_RADAR + USE_SATSIM > 0 ) then
      where (radar_ze_tot == R_UNDEF) radar_ze_tot = MAPL_UNDEF 
      where ( radar_lidar_only_freq_cloud == R_UNDEF) radar_lidar_only_freq_cloud = MAPL_UNDEF 
      where ( radar_lidar_tcc == R_UNDEF) radar_lidar_tcc = MAPL_UNDEF 
      where ( radar_cfad_ze == R_UNDEF) radar_cfad_ze = MAPL_UNDEF 
   endif

   if ( USE_SATSIM_MISR + USE_SATSIM > 0 ) then
      where (MISR_meanztop == R_UNDEF )  MISR_meanztop = MAPL_UNDEF
      where (MISR_meanztop .le. 0.0 )  MISR_meanztop = MAPL_UNDEF
      where (fq_MISR == R_UNDEF )  fq_MISR = MAPL_UNDEF
      where (MISR_cldarea == R_UNDEF ) MISR_cldarea = MAPL_UNDEF
      where (MISR_dist_model_layertops == R_UNDEF ) MISR_dist_model_layertops = MAPL_UNDEF
   endif
#endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! ISCCP exports
!-------------------------------------------------
!-------------------------------------------------


! 4 Cumulus (CU) subcategories
!    fq_isccp3D(2:3,6:7)  
!-------------------------------------------------

    if (associated(ISCCP_CU_OA)) then
       ISCCP_CU_OA  =  fq_isccp3D(:,:,2,6)
    endif

    if (associated(ISCCP_CU_OB)) then
       ISCCP_CU_OB  =  fq_isccp3D(:,:,3,6)
    endif

    if (associated(ISCCP_CU_UA)) then
       ISCCP_CU_UA  =  fq_isccp3D(:,:,2,7)
    endif

    if (associated(ISCCP_CU_UB)) then
       ISCCP_CU_UB  =  fq_isccp3D(:,:,3,7)
   endif


! 4 Stratocumulus STCU subcategories
!    ISCCP_isccp3D(4:5,6:7)  
!-------------------------------------------------

    if (associated(ISCCP_STCU_OA)) then
       ISCCP_STCU_OA  =  fq_isccp3D(:,:,4,6)
    endif

    if (associated(ISCCP_STCU_OB)) then
       ISCCP_STCU_OB  =  fq_isccp3D(:,:,5,6)
    endif

    if (associated(ISCCP_STCU_UA)) then
       ISCCP_STCU_UA  =  fq_isccp3D(:,:,4,7)
    endif

    if (associated(ISCCP_STCU_UB)) then
       ISCCP_STCU_UB  =  fq_isccp3D(:,:,5,7)
    endif

! 4 Stratus ST subcategories
!    fq_isccp3D(6:7,6:7)  
!-------------------------------------------------

    if (associated(ISCCP_ST_OA)) then
       ISCCP_ST_OA  =  fq_isccp3D(:,:,6,6)
    endif

    if (associated(ISCCP_ST_OB)) then
       ISCCP_ST_OB  =  fq_isccp3D(:,:,7,6)
    endif

    if (associated(ISCCP_ST_UA)) then
       ISCCP_ST_UA  =  fq_isccp3D(:,:,6,7)
    endif

    if (associated(ISCCP_ST_UB)) then
       ISCCP_ST_UB  =  fq_isccp3D(:,:,7,7)
    endif

! 4 Altocumulus ACU subcategories
!    fq_isccp3D(2:3,4:5)  
!-------------------------------------------------

    if (associated(ISCCP_ACU_OA)) then
       ISCCP_ACU_OA  =  fq_isccp3D(:,:,2,4)
    endif

    if (associated(ISCCP_ACU_OB)) then
       ISCCP_ACU_OB  =  fq_isccp3D(:,:,3,4)
    endif

    if (associated(ISCCP_ACU_UA)) then
       ISCCP_ACU_UA  =  fq_isccp3D(:,:,2,5)
    endif

    if (associated(ISCCP_ACU_UB)) then
       ISCCP_ACU_UB  =  fq_isccp3D(:,:,3,5)
    endif

! 4 Altostratus AST subcategories
!    fq_isccp3D(4:5,4:5)  
!-------------------------------------------------

    if (associated(ISCCP_AST_OA)) then
       ISCCP_AST_OA  =  fq_isccp3D(:,:,4,4)
    endif

    if (associated(ISCCP_AST_OB)) then
       ISCCP_AST_OB  =  fq_isccp3D(:,:,5,4)
    endif

    if (associated(ISCCP_AST_UA)) then
       ISCCP_AST_UA  =  fq_isccp3D(:,:,4,5)
    endif

    if (associated(ISCCP_AST_UB)) then
       ISCCP_AST_UB  =  fq_isccp3D(:,:,5,5)
    endif

! 4 Nimbostratus NST subcategories
!    fq_isccp3D(6:7,4:5)  
!-------------------------------------------------

    if (associated(ISCCP_NST_OA)) then
       ISCCP_NST_OA  =  fq_isccp3D(:,:,6,4)
    endif

    if (associated(ISCCP_NST_OB)) then
       ISCCP_NST_OB  =  fq_isccp3D(:,:,7,4)
    endif

    if (associated(ISCCP_NST_UA)) then
       ISCCP_NST_UA  =  fq_isccp3D(:,:,6,5)
    endif

    if (associated(ISCCP_NST_UB)) then
       ISCCP_NST_UB  =  fq_isccp3D(:,:,7,5)
    endif

! 6 Cirrus CI subcategories
!    fq_isccp3D(2:3,1:3)  
!-------------------------------------------------

    if (associated(ISCCP_CI_OA)) then
       ISCCP_CI_OA  =  fq_isccp3D(:,:,2,1)
    endif

    if (associated(ISCCP_CI_OB)) then
       ISCCP_CI_OB  =  fq_isccp3D(:,:,3,1)
    endif

    if (associated(ISCCP_CI_MA)) then
       ISCCP_CI_MA  =  fq_isccp3D(:,:,2,2)
    endif

    if (associated(ISCCP_CI_MB)) then
       ISCCP_CI_MB  =  fq_isccp3D(:,:,3,2)
    endif

    if (associated(ISCCP_CI_UA)) then
       ISCCP_CI_UA  =  fq_isccp3D(:,:,2,3)
    endif

    if (associated(ISCCP_CI_UB)) then
       ISCCP_CI_UB  =  fq_isccp3D(:,:,3,3)
    endif

! 6 Cirrostratus CIST subcategories
!    fq_isccp3D(4:5,1:3)  
!-------------------------------------------------

    if (associated(ISCCP_CIST_OA)) then
       ISCCP_CIST_OA  =  fq_isccp3D(:,:,4,1)
    endif

    if (associated(ISCCP_CIST_OB)) then
       ISCCP_CIST_OB  =  fq_isccp3D(:,:,5,1)
    endif

    if (associated(ISCCP_CIST_MA)) then
       ISCCP_CIST_MA  =  fq_isccp3D(:,:,4,2)
    endif

    if (associated(ISCCP_CIST_MB)) then
       ISCCP_CIST_MB  =  fq_isccp3D(:,:,5,2)
    endif

    if (associated(ISCCP_CIST_UA)) then
       ISCCP_CIST_UA  =  fq_isccp3D(:,:,4,3)
    endif

    if (associated(ISCCP_CIST_UB)) then
       ISCCP_CIST_UB  =  fq_isccp3D(:,:,5,3)
    endif

! 6 Cumulonimbus CB subcategories
!    fq_isccp3D(6:7,1:3)  
!-------------------------------------------------

    if (associated(ISCCP_CB_OA)) then
       ISCCP_CB_OA  =  fq_isccp3D(:,:,6,1)
    endif

    if (associated(ISCCP_CB_OB)) then
       ISCCP_CB_OB  =  fq_isccp3D(:,:,7,1)
    endif

    if (associated(ISCCP_CB_MA)) then
       ISCCP_CB_MA  =  fq_isccp3D(:,:,6,2)
    endif

    if (associated(ISCCP_CB_MB)) then
       ISCCP_CB_MB  =  fq_isccp3D(:,:,7,2)
    endif

    if (associated(ISCCP_CB_UA)) then
       ISCCP_CB_UA  =  fq_isccp3D(:,:,6,3)
    endif

    if (associated(ISCCP_CB_UB)) then
       ISCCP_CB_UB  =  fq_isccp3D(:,:,7,3)
    endif

! 7 Subvisible/Subdetection (SUBV) subcategories  
!    ISCCP_isccp3D(1,1:7)  
!-------------------------------------------------

    if (associated(ISCCP_SUBV1)) then
       ISCCP_SUBV1  =  fq_isccp3D(:,:,1,1)
    endif

    if (associated(ISCCP_SUBV2)) then
       ISCCP_SUBV2  =  fq_isccp3D(:,:,1,2)
    endif

    if (associated(ISCCP_SUBV3)) then
       ISCCP_SUBV3  =  fq_isccp3D(:,:,1,3)
    endif

    if (associated(ISCCP_SUBV4)) then
       ISCCP_SUBV4  =  fq_isccp3D(:,:,1,4)
    endif

    if (associated(ISCCP_SUBV5)) then
       ISCCP_SUBV5  =  fq_isccp3D(:,:,1,5)
    endif

    if (associated(ISCCP_SUBV6)) then
       ISCCP_SUBV6  =  fq_isccp3D(:,:,1,6)
    endif

    if (associated(ISCCP_SUBV7)) then
       ISCCP_SUBV7  =  fq_isccp3D(:,:,1,7)
    endif


! 7 pressure levels of all thicknesses  
!    ISCCP_isccp3D(:,:,:,1:7)  
!-------------------------------------------------

    if (associated(CLISCCP1)) then
       CLISCCP1  =  fq_isccp3D(:,:,:,1)
    endif

    if (associated(CLISCCP2)) then
       CLISCCP2  =  fq_isccp3D(:,:,:,2)
    endif

    if (associated(CLISCCP3)) then
       CLISCCP3  =  fq_isccp3D(:,:,:,3)
    endif

    if (associated(CLISCCP4)) then
       CLISCCP4  =  fq_isccp3D(:,:,:,4)
    endif

    if (associated(CLISCCP5)) then
       CLISCCP5  =  fq_isccp3D(:,:,:,5)
    endif

    if (associated(CLISCCP6)) then
       CLISCCP6  =  fq_isccp3D(:,:,:,6)
    endif

    if (associated(CLISCCP7)) then
       CLISCCP7  =  fq_isccp3D(:,:,:,7)
    endif



    if (associated(ISCCP1)) then
       ISCCP1  =  fq_isccp3D(:,:,:,1)
    endif

    if (associated(ISCCP2)) then
       ISCCP2  =  fq_isccp3D(:,:,:,2)
    endif

    if (associated(ISCCP3)) then
       ISCCP3  =  fq_isccp3D(:,:,:,3)
    endif

    if (associated(ISCCP4)) then
       ISCCP4  =  fq_isccp3D(:,:,:,4)
    endif

    if (associated(ISCCP5)) then
       ISCCP5  =  fq_isccp3D(:,:,:,5)
    endif

    if (associated(ISCCP6)) then
       ISCCP6  =  fq_isccp3D(:,:,:,6)
    endif

    if (associated(ISCCP7)) then
       ISCCP7  =  fq_isccp3D(:,:,:,7)
    endif


    if (associated(CLTISCCP)) then
       CLTISCCP  = reshape ( ISCCP_totalcldarea, (/ IM , JM /) ) 
    endif

   if (associated(TCLISCCP)) then
       TCLISCCP  = reshape ( ISCCP_totalcldarea, (/ IM , JM /) ) 
    endif

   if (associated(PCTISCCP)) then
      PCTISCCP  = reshape ( ISCCP_meanptop, (/ IM , JM /) )
    endif

   if (associated(CTPISCCP)) then
      CTPISCCP  = reshape ( ISCCP_meanptop, (/ IM , JM /) )
    endif

   if (associated(ALBISCCP)) then
       ALBISCCP  = reshape ( ISCCP_meanalbedocld, (/ IM , JM /) )
    endif

   if (associated(TBISCCP)) then
       TBISCCP  = reshape ( ISCCP_meanallskybrighttemp, (/ IM , JM /) )
    endif


! CALIPSO/lidar exports
!-------------------------------------------------
!-------------------------------------------------

    if (associated(CLCALIPSO)) then
       CLCALIPSO = reshape( lidar_lidarcld(:,LM:1:-1), (/IM , JM , LM /) )
    endif

    if (associated(CLLCALIPSO)) then
       CLLCALIPSO = reshape( lidar_cldlayer(:,1), (/IM , JM /) ) ! low
    endif

    if (associated(CLMCALIPSO)) then
       CLMCALIPSO = reshape( lidar_cldlayer(:,2), (/IM , JM /) ) ! middle
    endif

    if (associated(CLHCALIPSO)) then
       CLHCALIPSO = reshape( lidar_cldlayer(:,3), (/IM , JM /) ) ! high
    endif

    if (associated(CLTCALIPSO)) then
       CLTCALIPSO = reshape( lidar_cldlayer(:,4), (/IM , JM /) ) ! total
    endif

    if (associated(PARASOLREFL0)) then
       PARASOLREFL0 = reshape( lidar_parasolrefl, (/IM , JM, PARASOL_NREFL /) )
    endif

    if (associated(PARASOLREFL1)) then
       PARASOLREFL1 = reshape( lidar_parasolrefl(:,1), (/IM , JM /) ) 
    endif

    if (associated(PARASOLREFL2)) then
       PARASOLREFL2 = reshape( lidar_parasolrefl(:,2), (/IM , JM /) ) 
    endif

    if (associated(PARASOLREFL3)) then
       PARASOLREFL3 = reshape( lidar_parasolrefl(:,3), (/IM , JM /) ) 
    endif

    if (associated(PARASOLREFL4)) then
       PARASOLREFL4 = reshape( lidar_parasolrefl(:,4), (/IM , JM /) ) 
    endif

   if (associated(PARASOLREFL5)) then
      PARASOLREFL5 = reshape( lidar_parasolrefl(:,5), (/IM , JM /) ) 
   endif

    if (associated(CFADLIDARSR532_01)) then
       CFADLIDARSR532_01 = reshape( lidar_cfad_sr(:,1,LM:1:-1) , (/IM , JM , LM /) )
    endif

    if (associated(CFADLIDARSR532_02)) then
       CFADLIDARSR532_02 = reshape( lidar_cfad_sr(:,2,LM:1:-1) , (/IM , JM , LM /) )
    endif


    if (associated(CFADLIDARSR532_03)) then
       CFADLIDARSR532_03 = reshape( lidar_cfad_sr(:,3,LM:1:-1) , (/IM , JM , LM /) )
    endif


    if (associated(CFADLIDARSR532_04)) then
       CFADLIDARSR532_04 = reshape( lidar_cfad_sr(:,4,LM:1:-1) , (/IM , JM , LM /) )
    endif


    if (associated(CFADLIDARSR532_05)) then
       CFADLIDARSR532_05= reshape( lidar_cfad_sr (:,5,LM:1:-1) , (/IM , JM , LM /) )
    endif


    if (associated(CFADLIDARSR532_06)) then
       CFADLIDARSR532_06 = reshape( lidar_cfad_sr (:,6,LM:1:-1) , (/IM , JM , LM /) )
    endif


    if (associated(CFADLIDARSR532_07)) then
       CFADLIDARSR532_07 = reshape( lidar_cfad_sr(:,7,LM:1:-1) , (/IM , JM , LM /) ) 
    endif


    if (associated(CFADLIDARSR532_08)) then
       CFADLIDARSR532_08 = reshape( lidar_cfad_sr (:,8,LM:1:-1) , (/IM , JM , LM /) )
    endif


    if (associated(CFADLIDARSR532_09)) then
       CFADLIDARSR532_09 = reshape( lidar_cfad_sr(:,9,LM:1:-1) , (/IM , JM , LM /) ) 
    endif


    if (associated(CFADLIDARSR532_10)) then
       CFADLIDARSR532_10 = reshape( lidar_cfad_sr(:,10,LM:1:-1) , (/IM , JM , LM /) ) 
    endif


    if (associated(CFADLIDARSR532_11)) then
       CFADLIDARSR532_11 = reshape( lidar_cfad_sr(:,11,LM:1:-1) , (/IM , JM , LM /) ) 
    endif


    if (associated(CFADLIDARSR532_12)) then
       CFADLIDARSR532_12 = reshape( lidar_cfad_sr(:,12,LM:1:-1) , (/IM , JM , LM /) ) 
    endif


    if (associated(CFADLIDARSR532_13)) then
       CFADLIDARSR532_13 = reshape( lidar_cfad_sr(:,13,LM:1:-1) , (/IM , JM , LM /) ) 
    endif


    if (associated(CFADLIDARSR532_14)) then
       CFADLIDARSR532_14 = reshape( lidar_cfad_sr(:,14,LM:1:-1) , (/IM , JM , LM /) ) 
    endif


    if (associated(CFADLIDARSR532_15)) then
       CFADLIDARSR532_15 = reshape( lidar_cfad_sr(:,15,LM:1:-1) , (/IM , JM , LM /) ) 
    endif

    if (associated(LIDARPMOL)) then
      LIDARPMOL = reshape( lidar_pmol(:,LM:1:-1), (/IM , JM , LM /) )
    endif

    if (associated(LIDARPTOT)) then
     LIDARPTOT =reshape(  sum( lidar_beta_tot(:,:,LM:1:-1) , 2 )/(1.*NCOLumns) , (/IM , JM , LM /) )
    endif

    if (associated(LIDARTAUTOT)) then
     LIDARTAUTOT =reshape(  sum( lidar_tau_tot(:,:,LM:1:-1) , 2 )/(1.*NCOLumns) , (/IM , JM , LM /) )
    endif



! CLOUDSAT/radar exports
!-------------------------------------------------
!-------------------------------------------------

#define SATSIMTEMP
#ifdef SATSIMTEMP 
       if (associated(RADARZETOT)) then 
       radar_ze_tot_mask = (radar_ze_tot /= MAPL_UNDEF)                            ! mask out undefined columns
       ncvalid = count(radar_ze_tot_mask, dim=2)                                   ! number of defined columns
       radar_ze_tot_max = maxval(radar_ze_tot,dim=2,mask=radar_ze_tot_mask)        ! maximum among columns
       radar_ze_tot_tmp = spread(radar_ze_tot_max,dim=2,ncopies=NCOLUMNS)          ! max copied to all columns
       where (radar_ze_tot_mask) &
         radar_ze_tot_tmp = 10.0**((radar_ze_tot-radar_ze_tot_tmp)/10.0)           ! convert to rel reflectivity
       where (ncvalid > 0)
         radar_ze_tot_mean = &
           sum(radar_ze_tot_tmp,dim=2,mask=radar_ze_tot_mask)/real(ncvalid)        ! mean rel reflectivity
         radar_ze_tot_mean = radar_ze_tot_max + 10.0*log10(radar_ze_tot_mean)      ! convert back to dBZ
       elsewhere
         radar_ze_tot_mean = MAPL_UNDEF
       end where
       RADARZETOT=reshape(radar_ze_tot_mean(:,LM:1:-1),(/IM,JM,LM/))
    end if 
   !if (associated(RADARZETOT)) then 
   ! radar_ze_tot_mask = (radar_ze_tot .ne. MAPL_UNDEF) 
   ! ncvalid = count( radar_ze_tot_mask, 2 ) 
   ! where (ncvalid > 0) 
   !    radar_ze_tot_max = maxval( radar_ze_tot, 2, radar_ze_tot_mask )
   !    radar_ze_tot_mean = radar_ze_tot_max + 10.0*log10( sum( 10.0**((radar_ze_tot - spread(radar_ze_tot_max, 2, NCOLUMNS))/10.0), 2, radar_ze_tot_mask ) / ncvalid )
   ! elsewhere 
   !    radar_ze_tot_mean = MAPL_UNDEF 
   ! end where 
   ! RADARZETOT=reshape( radar_ze_tot_mean(:,LM:1:-1), (/IM , JM , LM /) ) 
   !end if 
#else 
    if (associated(RADARZETOT)) then
     RADARZETOT=reshape(  sum( radar_ze_tot(:,:,LM:1:-1) , 2 )/(1.*NCOLumns) , (/IM , JM , LM /) )
    endif
#endif ! SATSIMTEMP


    if (associated(RADARLTCC)) then
       RADARLTCC= reshape( radar_lidar_tcc, (/IM , JM /) )
    endif

    if (associated(CLCALIPSO2)) then
       CLCALIPSO2= reshape( radar_lidar_only_freq_cloud(:,LM:1:-1), (/IM , JM , LM /) )
    endif

    if (associated(CLOUDSATCFAD01)) then
       CLOUDSATCFAD01 = reshape( radar_cfad_ze(:,1,LM:1:-1) , (/IM , JM , LM /) )
    endif

    if (associated(CLOUDSATCFAD02)) then
       CLOUDSATCFAD02 = reshape( radar_cfad_ze(:,2,LM:1:-1) , (/IM , JM , LM /) )
    endif


    if (associated(CLOUDSATCFAD03)) then
       CLOUDSATCFAD03 = reshape( radar_cfad_ze(:,3,LM:1:-1) , (/IM , JM , LM /) )
    endif


    if (associated(CLOUDSATCFAD04)) then
       CLOUDSATCFAD04 = reshape( radar_cfad_ze(:,4,LM:1:-1) , (/IM , JM , LM /) )
    endif


    if (associated(CLOUDSATCFAD05)) then
       CLOUDSATCFAD05= reshape( radar_cfad_ze (:,5,LM:1:-1) , (/IM , JM , LM /) )
    endif


    if (associated(CLOUDSATCFAD06)) then
       CLOUDSATCFAD06 = reshape( radar_cfad_ze (:,6,LM:1:-1) , (/IM , JM , LM /) )
    endif


    if (associated(CLOUDSATCFAD07)) then
       CLOUDSATCFAD07 = reshape( radar_cfad_ze(:,7,LM:1:-1) , (/IM , JM , LM /) ) 
    endif


    if (associated(CLOUDSATCFAD08)) then
       CLOUDSATCFAD08 = reshape( radar_cfad_ze (:,8,LM:1:-1) , (/IM , JM , LM /) )
    endif


    if (associated(CLOUDSATCFAD09)) then
       CLOUDSATCFAD09 = reshape( radar_cfad_ze(:,9,LM:1:-1) , (/IM , JM , LM /) ) 
    endif


    if (associated(CLOUDSATCFAD10)) then
       CLOUDSATCFAD10 = reshape( radar_cfad_ze(:,10,LM:1:-1) , (/IM , JM , LM /) ) 
    endif


    if (associated(CLOUDSATCFAD11)) then
       CLOUDSATCFAD11 = reshape( radar_cfad_ze(:,11,LM:1:-1) , (/IM , JM , LM /) ) 
    endif


    if (associated(CLOUDSATCFAD12)) then
       CLOUDSATCFAD12 = reshape( radar_cfad_ze(:,12,LM:1:-1) , (/IM , JM , LM /) ) 
    endif


    if (associated(CLOUDSATCFAD13)) then
       CLOUDSATCFAD13 = reshape( radar_cfad_ze(:,13,LM:1:-1) , (/IM , JM , LM /) ) 
    endif

    if (associated(CLOUDSATCFAD14)) then
       CLOUDSATCFAD14 = reshape( radar_cfad_ze(:,14,LM:1:-1) , (/IM , JM , LM /) ) 
    endif


    if (associated(CLOUDSATCFAD15)) then
       CLOUDSATCFAD15 = reshape( radar_cfad_ze(:,15,LM:1:-1) , (/IM , JM , LM /) ) 
    endif



! MODIS exports -- currently available only from COSP

    if (associated(MDSCLDFRCTTL)) then
      MDSCLDFRCTTL = reshape( MODIS_Cloud_Fraction_Total_Mean, (/IM , JM /) ) 
    endif

    if (associated(MDSCLDFRCWTR)) then
      MDSCLDFRCWTR = reshape( MODIS_Cloud_Fraction_Water_Mean, (/IM , JM /) ) 
    endif

    if (associated(MDSCLDFRCH2O)) then
      MDSCLDFRCH2O = reshape( MODIS_Cloud_Fraction_Water_Mean, (/IM , JM /) ) 
    endif

    if (associated(MDSCLDFRCICE)) then
      MDSCLDFRCICE = reshape( MODIS_Cloud_Fraction_Ice_Mean, (/IM , JM /) ) 
    endif

    if (associated(MDSCLDFRCHI)) then
      MDSCLDFRCHI = reshape( MODIS_Cloud_Fraction_High_Mean, (/IM , JM /) ) 
    endif

    if (associated(MDSCLDFRCMID)) then
      MDSCLDFRCMID = reshape( MODIS_Cloud_Fraction_Mid_Mean, (/IM , JM /) ) 
    endif

    if (associated(MDSCLDFRCLO)) then
      MDSCLDFRCLO = reshape( MODIS_Cloud_Fraction_Low_Mean, (/IM , JM /) ) 
    endif

    if (associated(MDSOPTHCKTTL)) then
      MDSOPTHCKTTL = reshape( MODIS_Optical_Thickness_Total_Mean, (/IM , JM /) ) 
    endif

    if (associated(MDSOPTHCKWTR)) then
      MDSOPTHCKWTR = reshape( MODIS_Optical_Thickness_Water_Mean, (/IM , JM /) ) 
    endif

    if (associated(MDSOPTHCKH2O)) then
      MDSOPTHCKH2O = reshape( MODIS_Optical_Thickness_Water_Mean, (/IM , JM /) ) 
    endif

    if (associated(MDSOPTHCKICE)) then
      MDSOPTHCKICE = reshape( MODIS_Optical_Thickness_Ice_Mean, (/IM , JM /) ) 
    endif

    if (associated(MDSOPTHCKTTLLG)) then
      MDSOPTHCKTTLLG = reshape( MODIS_Optical_Thickness_Total_LogMean, (/IM , JM /) ) 
    endif

    if (associated(MDSOPTHCKWTRLG)) then
      MDSOPTHCKWTRLG = reshape( MODIS_Optical_Thickness_Water_LogMean, (/IM , JM /) ) 
    endif

    if (associated(MDSOPTHCKH2OLG)) then
      MDSOPTHCKH2OLG = reshape( MODIS_Optical_Thickness_Water_LogMean, (/IM , JM /) ) 
    endif

    if (associated(MDSOPTHCKICELG)) then
      MDSOPTHCKICELG = reshape( MODIS_Optical_Thickness_Ice_LogMean, (/IM , JM /) ) 
    endif

    if (associated(MDSCLDSZWTR)) then
      MDSCLDSZWTR = reshape( MODIS_Cloud_Particle_Size_Water_Mean, (/IM , JM /) ) 
    endif

    if (associated(MDSCLDSZH20)) then
      MDSCLDSZH20 = reshape( MODIS_Cloud_Particle_Size_Water_Mean, (/IM , JM /) ) 
    endif

    if (associated(MDSCLDSZICE)) then
      MDSCLDSZICE = reshape( MODIS_Cloud_Particle_Size_Ice_Mean, (/IM , JM /) ) 
    endif

    if (associated(MDSCLDTOPPS)) then
      MDSCLDTOPPS = reshape( MODIS_Cloud_Top_Pressure_Total_Mean, (/IM , JM /) ) 
    endif

    if (associated(MDSWTRPATH)) then
      MDSWTRPATH = reshape( MODIS_Liquid_Water_Path_Mean, (/IM , JM /) ) 
    endif

    if (associated(MDSH2OPATH)) then
      MDSH2OPATH = reshape( MODIS_Liquid_Water_Path_Mean, (/IM , JM /) ) 
    endif

    if (associated(MDSICEPATH)) then
      MDSICEPATH = reshape( MODIS_Ice_Water_Path_Mean, (/IM , JM /) ) 
    endif

    if (associated(MDSTAUPRSHIST11)) then
      MDSTAUPRSHIST11 = reshape( MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:,1,1),(/IM , JM /))
    endif

    if (associated(MDSTAUPRSHIST12)) then
      MDSTAUPRSHIST12 = reshape( MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:,1,2),(/IM , JM /))
    endif

    if (associated(MDSTAUPRSHIST13)) then
      MDSTAUPRSHIST13 = reshape( MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:,1,3),(/IM , JM /))
    endif

    if (associated(MDSTAUPRSHIST14)) then
      MDSTAUPRSHIST14 = reshape( MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:,1,4),(/IM , JM /))
    endif

   if (associated(MDSTAUPRSHIST15)) then
      MDSTAUPRSHIST15 = reshape( MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:,1,5),(/IM , JM /))
    endif

    if (associated(MDSTAUPRSHIST16)) then
      MDSTAUPRSHIST16 = reshape( MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:,1,6),(/IM , JM /))
    endif

   if (associated(MDSTAUPRSHIST17)) then
      MDSTAUPRSHIST17 = reshape( MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:,1,7),(/IM , JM /))
    endif

    if (associated(MDSTAUPRSHIST21)) then
      MDSTAUPRSHIST21 = reshape( MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:,2,1),(/IM , JM /))
    endif

    if (associated(MDSTAUPRSHIST22)) then
      MDSTAUPRSHIST22 = reshape( MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:,2,2),(/IM , JM /))
    endif

    if (associated(MDSTAUPRSHIST23)) then
      MDSTAUPRSHIST23 = reshape( MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:,2,3),(/IM , JM /))
    endif

    if (associated(MDSTAUPRSHIST24)) then
      MDSTAUPRSHIST24 = reshape( MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:,2,4),(/IM , JM /))
    endif

   if (associated(MDSTAUPRSHIST25)) then
      MDSTAUPRSHIST25 = reshape( MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:,2,5),(/IM , JM /))
    endif

    if (associated(MDSTAUPRSHIST26)) then
      MDSTAUPRSHIST26 = reshape( MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:,2,6),(/IM , JM /))
    endif

   if (associated(MDSTAUPRSHIST27)) then
      MDSTAUPRSHIST27 = reshape( MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:,2,7),(/IM , JM /))
    endif

    if (associated(MDSTAUPRSHIST31)) then
      MDSTAUPRSHIST31 = reshape( MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:,3,1),(/IM , JM /))
    endif

    if (associated(MDSTAUPRSHIST32)) then
      MDSTAUPRSHIST32 = reshape( MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:,3,2),(/IM , JM /))
    endif

    if (associated(MDSTAUPRSHIST33)) then
      MDSTAUPRSHIST33 = reshape( MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:,3,3),(/IM , JM /))
    endif

    if (associated(MDSTAUPRSHIST34)) then
      MDSTAUPRSHIST34 = reshape( MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:,3,4),(/IM , JM /))
    endif

   if (associated(MDSTAUPRSHIST35)) then
      MDSTAUPRSHIST35 = reshape( MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:,3,5),(/IM , JM /))
    endif

    if (associated(MDSTAUPRSHIST36)) then
      MDSTAUPRSHIST36 = reshape( MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:,3,6),(/IM , JM /))
    endif

   if (associated(MDSTAUPRSHIST37)) then
      MDSTAUPRSHIST37 = reshape( MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:,3,7),(/IM , JM /))
    endif

    if (associated(MDSTAUPRSHIST41)) then
      MDSTAUPRSHIST41 = reshape( MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:,4,1),(/IM , JM /))
    endif

    if (associated(MDSTAUPRSHIST42)) then
      MDSTAUPRSHIST42 = reshape( MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:,4,2),(/IM , JM /))
    endif

    if (associated(MDSTAUPRSHIST43)) then
      MDSTAUPRSHIST43 = reshape( MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:,4,3),(/IM , JM /))
    endif

    if (associated(MDSTAUPRSHIST44)) then
      MDSTAUPRSHIST44 = reshape( MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:,4,4),(/IM , JM /))
    endif

   if (associated(MDSTAUPRSHIST45)) then
      MDSTAUPRSHIST45 = reshape( MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:,4,5),(/IM , JM /))
    endif

    if (associated(MDSTAUPRSHIST46)) then
      MDSTAUPRSHIST46 = reshape( MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:,4,6),(/IM , JM /))
    endif

   if (associated(MDSTAUPRSHIST47)) then
      MDSTAUPRSHIST47 = reshape( MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:,4,7),(/IM , JM /))
    endif

    if (associated(MDSTAUPRSHIST51)) then
      MDSTAUPRSHIST51 = reshape( MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:,5,1),(/IM , JM /))
    endif

    if (associated(MDSTAUPRSHIST52)) then
      MDSTAUPRSHIST52 = reshape( MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:,5,2),(/IM , JM /))
    endif

    if (associated(MDSTAUPRSHIST53)) then
      MDSTAUPRSHIST53 = reshape( MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:,5,3),(/IM , JM /))
    endif

    if (associated(MDSTAUPRSHIST54)) then
      MDSTAUPRSHIST54 = reshape( MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:,5,4),(/IM , JM /))
    endif

   if (associated(MDSTAUPRSHIST55)) then
      MDSTAUPRSHIST55 = reshape( MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:,5,5),(/IM , JM /))
    endif

    if (associated(MDSTAUPRSHIST56)) then
      MDSTAUPRSHIST56 = reshape( MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:,5,6),(/IM , JM /))
    endif

   if (associated(MDSTAUPRSHIST57)) then
      MDSTAUPRSHIST57 = reshape( MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:,5,7),(/IM , JM /))
    endif

    if (associated(MDSTAUPRSHIST61)) then
      MDSTAUPRSHIST61 = reshape( MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:,6,1),(/IM , JM /))
    endif

    if (associated(MDSTAUPRSHIST62)) then
      MDSTAUPRSHIST62 = reshape( MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:,6,2),(/IM , JM /))
    endif

    if (associated(MDSTAUPRSHIST63)) then
      MDSTAUPRSHIST63 = reshape( MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:,6,3),(/IM , JM /))
    endif

    if (associated(MDSTAUPRSHIST64)) then
      MDSTAUPRSHIST64 = reshape( MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:,6,4),(/IM , JM /))
    endif

   if (associated(MDSTAUPRSHIST65)) then
      MDSTAUPRSHIST65 = reshape( MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:,6,5),(/IM , JM /))
    endif

    if (associated(MDSTAUPRSHIST66)) then
      MDSTAUPRSHIST66 = reshape( MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:,6,6),(/IM , JM /))
    endif

   if (associated(MDSTAUPRSHIST67)) then
      MDSTAUPRSHIST67 = reshape( MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:,6,7),(/IM , JM /))
    endif

    if (associated(MDSTAUPRSHIST71)) then
      MDSTAUPRSHIST71 = reshape( MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:,7,1),(/IM , JM /))
    endif

    if (associated(MDSTAUPRSHIST72)) then
      MDSTAUPRSHIST72 = reshape( MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:,7,2),(/IM , JM /))
    endif

    if (associated(MDSTAUPRSHIST73)) then
      MDSTAUPRSHIST73 = reshape( MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:,7,3),(/IM , JM /))
    endif

    if (associated(MDSTAUPRSHIST74)) then
      MDSTAUPRSHIST74 = reshape( MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:,7,4),(/IM , JM /))
    endif

   if (associated(MDSTAUPRSHIST75)) then
      MDSTAUPRSHIST75 = reshape( MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:,7,5),(/IM , JM /))
    endif

    if (associated(MDSTAUPRSHIST76)) then
      MDSTAUPRSHIST76 = reshape( MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:,7,6),(/IM , JM /))
    endif

   if (associated(MDSTAUPRSHIST77)) then
      MDSTAUPRSHIST77 = reshape( MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:,7,7),(/IM , JM /))
    endif

    if (associated(MISRMNCLDTP)) then
      MISRMNCLDTP = reshape( MISR_meanztop, (/IM , JM /) ) 
    endif

    if (associated(MISRCLDAREA)) then
      MISRCLDAREA = reshape( MISR_cldarea, (/IM , JM /) ) 
    endif

    if (associated(MISRLYRTP0)) then
      MISRLYRTP0 = reshape( MISR_dist_model_layertops(:,1), (/IM , JM /) )
    endif

    if (associated(MISRLYRTP250)) then
      MISRLYRTP250 = reshape( MISR_dist_model_layertops(:,2), (/IM , JM /) )
    endif

    if (associated(MISRLYRTP750)) then
      MISRLYRTP750 = reshape( MISR_dist_model_layertops(:,3), (/IM , JM /) )
    endif

    if (associated(MISRLYRTP1250)) then
      MISRLYRTP1250 = reshape( MISR_dist_model_layertops(:,4), (/IM , JM /) )
    endif

    if (associated(MISRLYRTP1750)) then
      MISRLYRTP1750 = reshape( MISR_dist_model_layertops(:,5), (/IM , JM /) )
    endif

    if (associated(MISRLYRTP2250)) then
      MISRLYRTP2250 = reshape( MISR_dist_model_layertops(:,6), (/IM , JM /) )
    endif

    if (associated(MISRLYRTP2750)) then
      MISRLYRTP2750 = reshape( MISR_dist_model_layertops(:,7), (/IM , JM /) )
    endif

    if (associated(MISRLYRTP3500)) then
      MISRLYRTP3500 = reshape( MISR_dist_model_layertops(:,8), (/IM , JM /) )
    endif

    if (associated(MISRLYRTP4500)) then
      MISRLYRTP4500 = reshape( MISR_dist_model_layertops(:,9), (/IM , JM /) )
    endif

    if (associated(MISRLYRTP6000)) then
      MISRLYRTP6000 = reshape( MISR_dist_model_layertops(:,10), (/IM , JM /) )
    endif

    if (associated(MISRLYRTP8000)) then
      MISRLYRTP8000 = reshape( MISR_dist_model_layertops(:,11), (/IM , JM /) )
    endif

    if (associated(MISRLYRTP10000)) then
      MISRLYRTP10000 = reshape( MISR_dist_model_layertops(:,12), (/IM , JM /) )
    endif

    if (associated(MISRLYRTP12000)) then
      MISRLYRTP12000 = reshape( MISR_dist_model_layertops(:,13), (/IM , JM /) )
    endif

    if (associated(MISRLYRTP14000)) then
      MISRLYRTP14000 = reshape( MISR_dist_model_layertops(:,14), (/IM , JM /) )
    endif

    if (associated(MISRLYRTP16000)) then
      MISRLYRTP16000 = reshape( MISR_dist_model_layertops(:,15), (/IM , JM /) )
    endif

    if (associated(MISRLYRTP18000)) then
      MISRLYRTP18000 = reshape( MISR_dist_model_layertops(:,16), (/IM , JM /) )
    endif


    if (associated(MISRFQ0)) then
      MISRFQ0 = reshape( fq_MISR(:,:,1), (/IM , JM , 7 /) )
    endif

    if (associated(MISRFQ250)) then
      MISRFQ250 = reshape( fq_MISR(:,:,2), (/IM , JM , 7 /) )
    endif

    if (associated(MISRFQ750)) then
      MISRFQ750 = reshape( fq_MISR(:,:,3), (/IM , JM , 7 /) )
    endif

    if (associated(MISRFQ1250)) then
      MISRFQ1250 = reshape( fq_MISR(:,:,4), (/IM , JM , 7 /) )
    endif

    if (associated(MISRFQ1750)) then
      MISRFQ1750 = reshape( fq_MISR(:,:,5), (/IM , JM , 7 /) )
    endif

    if (associated(MISRFQ2250)) then
      MISRFQ2250 = reshape( fq_MISR(:,:,6), (/IM , JM , 7 /) )
    endif

    if (associated(MISRFQ2750)) then
      MISRFQ2750 = reshape( fq_MISR(:,:,7), (/IM , JM , 7 /) )
    endif

    if (associated(MISRFQ3500)) then
      MISRFQ3500 = reshape( fq_MISR(:,:,8), (/IM , JM , 7 /) )
    endif

    if (associated(MISRFQ4500)) then
      MISRFQ4500 = reshape( fq_MISR(:,:,9), (/IM , JM , 7 /) )
    endif

    if (associated(MISRFQ6000)) then
      MISRFQ6000 = reshape( fq_MISR(:,:,10), (/IM , JM , 7 /) )
    endif

    if (associated(MISRFQ8000)) then
      MISRFQ8000 = reshape( fq_MISR(:,:,11), (/IM , JM , 7 /) )
    endif

    if (associated(MISRFQ10000)) then
      MISRFQ10000 = reshape( fq_MISR(:,:,12), (/IM , JM , 7 /) )
    endif

    if (associated(MISRFQ12000)) then
      MISRFQ12000 = reshape( fq_MISR(:,:,13), (/IM , JM , 7 /) )
    endif

    if (associated(MISRFQ14000)) then
      MISRFQ14000 = reshape( fq_MISR(:,:,14), (/IM , JM , 7 /) )
    endif

    if (associated(MISRFQ16000)) then
      MISRFQ16000 = reshape( fq_MISR(:,:,15), (/IM , JM , 7 /) )
    endif

    if (associated(MISRFQ18000)) then
      MISRFQ18000 = reshape( fq_MISR(:,:,16), (/IM , JM , 7 /) )
    endif

    call ESMF_UserCompGetInternalState(gc, 'SatSim_State', WRAP, STATUS)
    VERIFY_(STATUS)
    self => wrap%ptr

    if (self%nmask_vars > 0) then
       ! get orb bundle
       call ESMF_StateGet(IMPORT,'SATORB',BUNDLE,RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_Get(STATE,ExportSpec=ExportSpec,RC=STATUS)
       VERIFY_(STATUS)
       do i=1,self%nmask_vars
          vindex = MAPL_VarSpecGetIndex(ExportSpec,self%export_name(i),RC=STATUS)
          VERIFY_(STATUS)
          call MAPL_VarSpecGet(ExportSpec(vindex),Dims=MAPL_dims,ungridded_dims=Ungridded_dims,RC=STATUS)
          VERIFY_(STATUS)
          ! mask from bundle
          call ESMFL_BundleGetPointerToData(BUNDLE,self%mask_name(i),ptr_mask,rc=status)
          VERIFY_(STATUS)
          if (self%newvar(i)) then

             ! determine what dimension pointer we need
             if (MAPL_dims == MAPL_DimsHorzOnly .and. (.not.associated(ungridded_dims))) then
                call MAPL_GetPointer(EXPORT, ptr2d_new,self%newvar_name(i), RC=STATUS)
                VERIFY_(STATUS)
                if (associated(ptr2d_new)) then
                   call MAPL_GetPointer(EXPORT,ptr2d,self%export_name(i),RC=STATUS)
                   VERIFY_(STATUS)
                   ptr2d_new=ptr2d
                   where(ptr_mask == MAPL_UNDEF)
                      ptr2d_new = MAPL_UNDEF
                   end where
                   nullify(ptr2d_new)
                   nullify(ptr2d)
                end if
             end if
             if (MAPL_dims == MAPL_DimsHorzOnly .and. associated(ungridded_dims)) then
                call MAPL_GetPointer(EXPORT, ptr3d_new,self%newvar_name(i), RC=STATUS)
                VERIFY_(STATUS)
                if (associated(ptr3d_new)) then
                   call MAPL_GetPointer(EXPORT,ptr3d,self%export_name(i),RC=STATUS)
                   VERIFY_(STATUS)
                   ptr3d_new=ptr3d
                   do j=1,ungridded_dims(1)
                      where(ptr_mask == MAPL_UNDEF)
                         ptr3d_new(:,:,j) = MAPL_UNDEF
                      end where
                   nullify(ptr3d_new)
                   nullify(ptr3d)
                   end do
                end if
             end if
             if (MAPL_dims == MAPL_DimsHorzVert) then
                call MAPL_GetPointer(EXPORT, ptr3d_new,self%newvar_name(i), RC=STATUS)
                VERIFY_(STATUS)
                if (associated(ptr3d_new)) then
                   call MAPL_GetPointer(EXPORT,ptr3d,self%export_name(i),RC=STATUS)
                   VERIFY_(STATUS)
                   ptr3d_new=ptr3d
                   do j=1,lm
                      where(ptr_mask == MAPL_UNDEF)
                         ptr3d_new(:,:,j) = MAPL_UNDEF
                      end where
                   end do
                   nullify(ptr3d_new)
                   nullify(ptr3d)
                end if
             end if

          else
             ! determine what dimension pointer we need
             if (MAPL_dims == MAPL_DimsHorzOnly .and. (.not.associated(ungridded_dims))) then
                call MAPL_GetPointer(EXPORT, ptr2d,self%newvar_name(i), RC=STATUS)
                VERIFY_(STATUS)
                if (associated(ptr2d)) then
                   where(ptr_mask == MAPL_UNDEF)
                      ptr2d = MAPL_UNDEF
                   end where
                   nullify(ptr2d)
                end if
             end if
             if (MAPL_dims == MAPL_DimsHorzOnly .and. associated(ungridded_dims)) then
                call MAPL_GetPointer(EXPORT, ptr3d,self%newvar_name(i), RC=STATUS)
                VERIFY_(STATUS)
                if (associated(ptr3d)) then
                   do j=1,ungridded_dims(1)
                      where(ptr_mask == MAPL_UNDEF)
                         ptr3d(:,:,j) = MAPL_UNDEF
                      end where
                   end do
                   nullify(ptr3d)
                end if
             end if
             if (MAPL_dims == MAPL_DimsHorzVert) then
                call MAPL_GetPointer(EXPORT, ptr3d,self%newvar_name(i), RC=STATUS)
                VERIFY_(STATUS)
                if (associated(ptr3d)) then
                   do j=1,lm
                      where(ptr_mask == MAPL_UNDEF)
                         ptr3d(:,:,j) = MAPL_UNDEF
                      end where
                   end do
                   nullify(ptr3d)
                end if
             end if

          end if

          nullify(ptr_mask)

       end do
    end if

     ! deallocate stuff
     !-----------------
     deallocate(         frac_out, __STAT__)
     !deallocate(      frac_outinv, __STAT__)
     deallocate(   lidar_beta_tot, __STAT__)
     deallocate(    lidar_tau_tot, __STAT__)
     deallocate(     radar_ze_tot, __STAT__)
     deallocate( radar_ze_tot_tmp, __STAT__)
     deallocate(radar_ze_tot_mask, __STAT__)

     call FREE_COSP_GRIDBOX(gbx)
     call FREE_COSP_SUBGRID(sgx)
     if (cfg%Lisccp_sim) call FREE_COSP_ISCCP(isccp)
     !call FREE_COSP_SGHYDRO(sghydro)
     if (cfg%Lradar_sim) call FREE_COSP_SGRADAR(sgradar)
     if (cfg%Llidar_sim) call FREE_COSP_SGLIDAR(sglidar)
     call FREE_COSP_VGRID(vgrid)
     if (cfg%Llidar_sim) call FREE_COSP_LIDARSTATS(stlidar)
     if (cfg%Lradar_sim) call FREE_COSP_RADARSTATS(stradar)
     if (cfg%Lmodis_sim) call FREE_COSP_MODIS(modis)
     if (cfg%Lmisr_sim ) call FREE_COSP_MISR(misr)

      RETURN_(ESMF_SUCCESS)

    end subroutine SIM_DRIVER

  end subroutine RUN


end module GEOS_SatsimGridCompMod

