! $Id$

! VERIFY_ and RETURN_ macros for error handling.

!#define UWDIAG 1
!#define PDFDIAG 1

#include "MAPL_Generic.h"

!=============================================================================
!BOP

! !MODULE: GEOS_Moist -- A Module to compute moist processes, including convection,
!   large-scale condensation and precipitation and cloud parameters.

! !INTERFACE:

module GEOS_MoistGridCompMod

  ! !USES:

  use RAS       ! using module that contains ras code

  use UWSHCU   ! using module that contains uwshcu code

  use gfdl2_cloud_microphys_mod

#ifndef _CUDA
  use CLOUDNEW, only: PROGNO_CLOUD, ICE_FRACTION, T_CLOUD_CTL, pdfcondensate, pdffrac, RADCOUPLE
#else
  use CLOUDNEW, only: &
       ! Subroutines
       PROGNO_CLOUD, ICE_FRACTION, &
       ! Derived Data Types
       T_CLOUD_CTL, &
       ! Inputs
       PP_DEV, EXNP_DEV, PPE_DEV, KH_DEV, DTS_DEV, SNOMAS_DEV, FRLAND_DEV, FRLANDICE_DEV, &
       RMFDTR_DEV, QLWDTR_DEV, U_DEV, V_DEV, QST3_DEV, &
       DZET_DEV, QDDF3_DEV, TEMPOR_DEV, CNV_FRACTION_DEV, TROPP_DEV, &
       ! Inoutputs
       TH_DEV, Q_DEV, QRN_CU_DEV, CNV_UPDFRC_DEV, QLW_LS_DEV, &
       QLW_AN_DEV, QIW_LS_DEV, QIW_AN_DEV, ANVFRC_DEV, CLDFRC_DEV, &
       QRN_SC_DEV, QSN_SC_DEV, SC_UPDFRC_DEV, &
       ! Outputs
       RAD_CLDFRC_DEV, RAD_QL_DEV, RAD_QI_DEV, RAD_QR_DEV, RAD_QS_DEV, RAD_QG_DEV, QPLS_DEV, &
       CLDREFFL_DEV, CLDREFFI_DEV, PRELS_DEV, PRECU_DEV, PREAN_DEV, &
       LSARF_DEV, CUARF_DEV, ANARF_DEV, SNRLS_DEV, SNRCU_DEV, SNRAN_DEV, SNRSC_DEV, &
       ! Working arrays 
       PFL_CN_DEV, PFI_CN_DEV, PFL_AN_DEV, PFL_SC_DEV, &
       PFI_AN_DEV, PFL_LS_DEV, PFI_LS_DEV, PFI_SC_DEV, &
       ! Diagnostics
       RHX_DEV, &
       REV_LS_DEV, REV_AN_DEV, REV_CN_DEV, REV_SC_DEV, &
       RSU_LS_DEV, RSU_AN_DEV, RSU_CN_DEV, RSU_SC_DEV, &
       ACLL_CN_DEV, ACIL_CN_DEV, ACLL_AN_DEV, ACLL_SC_DEV, &
       ACIL_AN_DEV, ACLL_LS_DEV, ACIL_LS_DEV, ACIL_SC_DEV, &
       PDFL_DEV, PDFI_DEV, FIXL_DEV, FIXI_DEV, &
       AUT_DEV, EVAPC_DEV, SDM_DEV, &
       SUBLC_DEV, FRZ_TT_DEV, DCNVL_DEV, DCNVI_DEV, &
       ALPHT_DEV, ALPH1_DEV, ALPH2_DEV, &
       CFPDF_DEV, RHCLR_DEV, DQRL_DEV, FRZ_PP_DEV, &
       VFALLICE_AN_DEV, VFALLICE_LS_DEV, &
       VFALLWAT_AN_DEV, VFALLWAT_LS_DEV, &
       VFALLSN_AN_DEV, VFALLSN_LS_DEV, &
       VFALLSN_CN_DEV, VFALLRN_AN_DEV, &
       VFALLRN_LS_DEV, VFALLRN_CN_DEV, &
       VFALLSN_SC_DEV, VFALLRN_SC_DEV, &
!!$       LIQANMOVE_DEV,  ICEANMOVE_DEV, &
!!$       DANCLD_DEV,     DLSCLD_DEV, &
!!$       CURAINMOVE_DEV, CUSNOWMOVE_DEV, &
       !CFPDFX is no longer calculated in PROGNO_CLOUD, but still an export
       !CFPDFX_DEV, &
       
       ! Constants
       ! CLDPARAMS Constants are loaded into constant memory
       CNV_BETA, ANV_BETA, LS_BETA, RH00, C_00, LWCRIT, C_ACC, &
       C_EV_R, C_EV_S, CLDVOL2FRC, RHSUP_ICE, SHR_EVAP_FAC, MIN_CLD_WATER, &
       CLD_EVP_EFF, NSMAX, LS_SDQV2, LS_SDQV3, LS_SDQVT1, ANV_SDQV2, &
       ANV_SDQV3, ANV_SDQVT1, ANV_TO_LS, N_WARM, N_ICE, N_ANVIL, &
       N_PBL, DISABLE_RAD, ANV_ICEFALL_C, LS_ICEFALL_C, REVAP_OFF_P, CNVENVFC, ANVENVFC, SCENVFC, LSENVFC, &
       WRHODEP, T_ICE_ALL, CNVICEPARAM, ICEFRPWR, CNVDDRFC, ANVDDRFC, &
       LSDDRFC, TANHRHCRIT, MINRHCRIT, MAXRHCRIT, TURNRHCRIT, MAXRHCRITLAND, &
       FR_LS_WAT, FR_LS_ICE, FR_AN_WAT, FR_AN_ICE, MIN_RL, MIN_RI, MAX_RL, &
       MAX_RI, FAC_RL, FAC_RI, PDFFLAG, ICE_SETTLE
  use cudafor
#endif

  use DDF

  use ESMF
  use MAPL, r8 => MAPL_R8
  use GEOS_UtilsMod
  use cldwat2m_micro
  use cldmacro
  use aer_cloud
  use RASPARAMS
  use CLDPARAMS
  use SHLWPARAMS
  use micro_mg3_0
 
!-srf-gf-scheme
  USE ConvPar_GF_GEOS5, only: gf_geos5_interface &
      ,maxiens, icumulus_gf, closure_choice, deep, shal, mid&
      ,DEBUG_GF,USE_SCALE_DEP,DICYCLE,Hcts&
      ,USE_TRACER_TRANSP,USE_TRACER_SCAVEN, TAU_DEEP, TAU_MID&
      ,USE_FLUX_FORM, USE_FCT, USE_TRACER_EVAP,ALP1
!-srf-gf-scheme

!ALT-protection for GF
! USE ConvectionMod, only: Disable_Convection
!ALT-protection for GF

!--kml--------------
 USE Aer_Actv_Single_Moment,only: Aer_Actv_1M_interface,INT_USE_AEROSOL_NN,USE_AEROSOL_NN
!--kml--------------
!  use m_zeit

 ! External lightning module
 USE Lightning_mod, only: HEMCO_FlashRate

  implicit none

!-srf-gf-scheme
  character(LEN=ESMF_MAXSTR):: CONVPAR_OPTION  ! GF, RAS, BOTH, NONE
  character(LEN=ESMF_MAXSTR):: AERO_PROVIDER  
!-srf-gf-scheme

  character(len=ESMF_MAXSTR) :: HYDROSTATIC ! TRUE or FALSE
  logical                    :: LHYDROSTATIC

  character(len=ESMF_MAXSTR) :: PHYS_HYDROSTATIC ! TRUE or FALSE
  logical                    :: LPHYS_HYDROSTATIC

  character(len=ESMF_MAXSTR) :: DIAGNOSE_PRECIP_TYPE ! TRUE or FALSE
  logical                    :: LDIAGNOSE_PRECIP_TYPE

  integer :: DOSHLW
  real    :: MGVERSION
  integer :: DOGRAUPEL

  private

  ! !PUBLIC MEMBER FUNCTIONS:

  public SetServices

  ! !DESCRIPTION:
  ! 
  !   {\tt GEOS\_MoistGridCompMod} implements moist processes in GEOS-5. These
  !   include all processes that involve phase changes in the atmosphere, such
  !   as large-scale condensation, convective clouds, and all rain and cloud
  !   formation. Its state consists of water vapor, various types of condensate,
  !   and fractions of various cloud types. 
  !   two moment cloud microphysics (Barahona et al., GMD, 2014.) can be run by setting CLDMACRO==2MOMENT. 
  !   When using 2-moment microphysics the number concentration of ice crystals and cloud droplets 
  !   are part of the state.
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
    type (ESMF_Config      )            :: CF

    integer      :: RFRSHINT
    integer      :: AVRGNINT
    integer      :: IQVAINC
    real         :: DT
    
    character(len=ESMF_MAXSTR) :: FRIENDLIES_NCPL , FRIENDLIES_NCPI , &
                                  FRIENDLIES_NRAIN, FRIENDLIES_NSNOW, FRIENDLIES_NGRAUPEL
    character(len=ESMF_MAXSTR) :: FRIENDLIES_QRAIN, FRIENDLIES_QSNOW, FRIENDLIES_QGRAUPEL

    !=============================================================================

    ! Begin...

    ! Get my name and set-up traceback handle
    ! ---------------------------------------

    Iam = 'SetServices'
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // Iam


    ! Set the Init entry point
    call MAPL_GridCompSetEntryPoint ( GC,  ESMF_METHOD_INITIALIZE,  Initialize, RC=status )
    VERIFY_(STATUS)

    ! Get the configuration from the component
    !-----------------------------------------

    call ESMF_GridCompGet( GC, CONFIG = CF, RC=STATUS )
    VERIFY_(STATUS)

    call ESMF_ConfigGetAttribute( CF, IQVAINC, Label='ALLOW_MOIST_AINC_UPDATE:',   default=0,        RC=STATUS)
    VERIFY_(STATUS)


    call ESMF_ConfigGetAttribute( CF, CONVPAR_OPTION, Label='CONVPAR_OPTION:', RC=STATUS) ! Note: Default set in GEOS_GcmGridComp.F90
    VERIFY_(STATUS)
    !~ call ESMF_ConfigGetAttribute( CF, DEBUG_GF, Label='DEBUG_GF:', default=0,RC=STATUS)
    !~ VERIFY_(STATUS)
    !ALT

    call ESMF_ConfigGetAttribute( CF, AERO_PROVIDER, Label='AERO_PROVIDER:', RC=STATUS) ! Note: Default set in GEOS_GcmGridComp.F90
    VERIFY_(STATUS)

    ! Set the Run entry point
    ! -----------------------

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,  Run,  &
         RC=STATUS)
    VERIFY_(STATUS)
    if ( iqvainc/=0 ) then
       call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,  AINC_UPDATE, &
            RC=STATUS)
       VERIFY_(STATUS)
    endif


    ! Set the state variable specs.
    ! -----------------------------

    call ESMF_ConfigGetAttribute ( CF, DT, Label="RUN_DT:",                                   RC=STATUS)
    VERIFY_(STATUS)

    call ESMF_ConfigGetAttribute( CF, RFRSHINT, Label="REFRESH_INTERVAL:",  default=nint(DT), RC=STATUS)
    VERIFY_(STATUS)

    call ESMF_ConfigGetAttribute( CF, AVRGNINT, Label='AVERAGING_INTERVAL:',default=RFRSHINT, RC=STATUS)
    VERIFY_(STATUS)

    ! Inititialize cloud microphysics (Options: 1MOMENT, 2MOMENT or GFDL)
    !--------------------------------------------------------------
    call ESMF_ConfigGetAttribute( CF, CLDMICRO, Label="CLDMICRO:",  default="1MOMENT", RC=STATUS)
    VERIFY_(STATUS)
    LCLDMICRO = adjustl(CLDMICRO)=="1MOMENT" .or. &
                adjustl(CLDMICRO)=="2MOMENT" .or. &
                adjustl(CLDMICRO)=="GFDL"
    _ASSERT( LCLDMICRO, 'needs informative message' )
    if (adjustl(CLDMICRO)=="2MOMENT") then
      call ESMF_ConfigGetAttribute( CF, MGVERSION, Label="MGVERSION:",  default=0.0, RC=STATUS)
    endif
    call ESMF_ConfigGetAttribute( CF, DOSHLW, Label="DOSHLW:",  default=0, RC=STATUS)

    call ESMF_ConfigGetAttribute( CF, HYDROSTATIC, Label="HYDROSTATIC:",  default="TRUE", RC=STATUS)
    VERIFY_(STATUS)
    if (adjustl(HYDROSTATIC)=="TRUE" ) LHYDROSTATIC=.true.
    if (adjustl(HYDROSTATIC)=="FALSE") LHYDROSTATIC=.false.

    call ESMF_ConfigGetAttribute( CF, PHYS_HYDROSTATIC, Label="PHYS_HYDROSTATIC:",  default="TRUE", RC=STATUS)
    VERIFY_(STATUS)
    if (adjustl(PHYS_HYDROSTATIC)=="TRUE" ) LPHYS_HYDROSTATIC=.true.
    if (adjustl(PHYS_HYDROSTATIC)=="FALSE") LPHYS_HYDROSTATIC=.false.

    call ESMF_ConfigGetAttribute( CF, DIAGNOSE_PRECIP_TYPE, Label="DIAGNOSE_PRECIP_TYPE:",  default="TRUE", RC=STATUS)
    VERIFY_(STATUS)
    if (adjustl(DIAGNOSE_PRECIP_TYPE)=="TRUE" ) LDIAGNOSE_PRECIP_TYPE=.true.
    if (adjustl(DIAGNOSE_PRECIP_TYPE)=="FALSE") LDIAGNOSE_PRECIP_TYPE=.false.

    FRIENDLIES_NCPI     = trim(COMP_NAME)
    FRIENDLIES_NCPL     = trim(COMP_NAME)
    FRIENDLIES_NRAIN    = trim(COMP_NAME)    
    FRIENDLIES_NSNOW    = trim(COMP_NAME)
    FRIENDLIES_NGRAUPEL = trim(COMP_NAME)
    FRIENDLIES_QRAIN    = trim(COMP_NAME)
    FRIENDLIES_QSNOW    = trim(COMP_NAME)
    FRIENDLIES_QGRAUPEL = trim(COMP_NAME)
   
    if(adjustl(CLDMICRO)=="2MOMENT") then
      if (MGVERSION==0) then    
        FRIENDLIES_NCPI = 'DYNAMICS:TURBULENCE'      
        FRIENDLIES_NCPL = 'DYNAMICS:TURBULENCE'
      endif
      if(MGVERSION==2) then
        call ESMF_ConfigGetAttribute( CF, DOGRAUPEL, Label="DOGRAUPEL:",  default=0, RC=STATUS)
        if (DOGRAUPEL == 0) then
          FRIENDLIES_NCPI = 'DYNAMICS:TURBULENCE'
          FRIENDLIES_NCPL = 'DYNAMICS:TURBULENCE'
          FRIENDLIES_NRAIN = 'DYNAMICS:TURBULENCE'
          FRIENDLIES_QRAIN = 'DYNAMICS:TURBULENCE'
          FRIENDLIES_NSNOW = 'DYNAMICS:TURBULENCE'
          FRIENDLIES_QSNOW = 'DYNAMICS:TURBULENCE'
        else
          FRIENDLIES_NCPI = 'DYNAMICS:TURBULENCE'
          FRIENDLIES_NCPL = 'DYNAMICS:TURBULENCE'
          FRIENDLIES_NRAIN = 'DYNAMICS:TURBULENCE'
          FRIENDLIES_QRAIN = 'DYNAMICS:TURBULENCE'
          FRIENDLIES_NSNOW = 'DYNAMICS:TURBULENCE'
          FRIENDLIES_QSNOW = 'DYNAMICS:TURBULENCE'
          FRIENDLIES_NGRAUPEL = 'DYNAMICS:TURBULENCE'
          FRIENDLIES_QGRAUPEL = 'DYNAMICS:TURBULENCE'
        endif
      endif
    endif
    if(adjustl(CLDMICRO)=="GFDL") then
      FRIENDLIES_QRAIN = 'DYNAMICS:TURBULENCE'
      FRIENDLIES_QSNOW = 'DYNAMICS:TURBULENCE'
      FRIENDLIES_QGRAUPEL = 'DYNAMICS:TURBULENCE'
    endif

    ! !INTERNAL STATE:

    !BOS

    call MAPL_AddInternalSpec(GC,                                  &
         SHORT_NAME = 'Q',                                         &
         LONG_NAME  = 'specific_humidity',                         &
         UNITS      = 'kg kg-1',                                   &
         FRIENDLYTO = 'DYNAMICS:TURBULENCE:CHEMISTRY:ANALYSIS',    &
         default    = 1.0e-6,                                      &
         RESTART    = MAPL_RestartRequired,                        &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)                                                                          

    call MAPL_AddInternalSpec(GC,                                        &
         SHORT_NAME = 'QLLS',                                            &
         LONG_NAME  = 'mass_fraction_of_large_scale_cloud_liquid_water', &
         UNITS      = 'kg kg-1',                                         &
         FRIENDLYTO = 'DYNAMICS:TURBULENCE',                             &
         DIMS       = MAPL_DimsHorzVert,                                 &
         VLOCATION  = MAPL_VLocationCenter,                   RC=STATUS  )  
    VERIFY_(STATUS)                                                                          

    call MAPL_AddInternalSpec(GC,                                       &
         SHORT_NAME = 'QLCN',                                           &
         LONG_NAME  = 'mass_fraction_of_convective_cloud_liquid_water', &
         UNITS      = 'kg kg-1',                                        &
         FRIENDLYTO = 'DYNAMICS:TURBULENCE',                            &
         DIMS       = MAPL_DimsHorzVert,                                &
         VLOCATION  = MAPL_VLocationCenter,                  RC=STATUS  )  
    VERIFY_(STATUS)                                                                          


    call MAPL_AddInternalSpec(GC,                                  &
         SHORT_NAME = 'CLLS',                                      &
         LONG_NAME  = 'large_scale_cloud_area_fraction',           &
         UNITS      = '1',                                         &
         FRIENDLYTO = 'DYNAMICS',                                  &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)                                                                          

    call MAPL_AddInternalSpec(GC,                                  &
         SHORT_NAME = 'CLCN',                                      &
         LONG_NAME  = 'convective_cloud_area_fraction',            &
         UNITS      = '1',                                         &
         FRIENDLYTO = 'DYNAMICS',                                  &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)                                                                          

    call MAPL_AddInternalSpec(GC,                                     &
         SHORT_NAME = 'QILS',                                         &
         LONG_NAME  = 'mass_fraction_of_large_scale_cloud_ice_water', &
         UNITS      = 'kg kg-1',                                      &
         FRIENDLYTO = 'DYNAMICS:TURBULENCE',                          &
         DIMS       = MAPL_DimsHorzVert,                              &
         VLOCATION  = MAPL_VLocationCenter,                RC=STATUS  )  
    VERIFY_(STATUS)                                                                          


    call MAPL_AddInternalSpec(GC,                                    &
         SHORT_NAME = 'QICN',                                        &
         LONG_NAME  = 'mass_fraction_of_convective_cloud_ice_water', &
         UNITS      = 'kg kg-1',                                     &
         FRIENDLYTO = 'DYNAMICS:TURBULENCE',                         &
         DIMS       = MAPL_DimsHorzVert,                             &
         VLOCATION  = MAPL_VLocationCenter,               RC=STATUS  )  
    VERIFY_(STATUS)                                                                          

 

    call MAPL_AddInternalSpec(GC,                                  &
         SHORT_NAME ='NCPL',                                       &
         LONG_NAME  ='particle_number_for_liquid_cloud',           &
         UNITS      ='kg-1',                                       &
         FRIENDLYTO = trim(FRIENDLIES_NCPL),                       &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,          &
         DEFAULT = 50.0e6 ,   RC=STATUS  )  
     
         VERIFY_(STATUS)                                                          

    call MAPL_AddInternalSpec(GC,                                  &
         SHORT_NAME ='NCPI',                                       &
         LONG_NAME  ='particle_number_for_ice_cloud',              &
         UNITS      ='kg-1',                                       &
         FRIENDLYTO = trim(FRIENDLIES_NCPI),                       &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,        &
         DEFAULT = 1.0e3,   RC=STATUS  )  
    
        VERIFY_(STATUS)                                                   


    call MAPL_AddInternalSpec(GC,                                  &
         SHORT_NAME ='NRAIN',                                       &
         LONG_NAME  ='particle_number_for_rain',           &
         UNITS      ='kg-1',                                       &
         FRIENDLYTO = trim(FRIENDLIES_NRAIN),                       &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,          &
         DEFAULT = 0.0 ,   RC=STATUS  )  
     
         VERIFY_(STATUS)                                                          

    call MAPL_AddInternalSpec(GC,                                  &
         SHORT_NAME ='NSNOW',                                       &
         LONG_NAME  ='particle_number_for_snow',              &
         UNITS      ='kg-1',                                       &
         FRIENDLYTO = trim(FRIENDLIES_NSNOW),                       &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,        &
         DEFAULT = 0.0,   RC=STATUS  )  
    
        VERIFY_(STATUS)                                                   

       call MAPL_AddInternalSpec(GC,                                  &
         SHORT_NAME ='NGRAUPEL',                                       &
         LONG_NAME  ='particle_number_for_graupel',              &
         UNITS      ='kg-1',                                       &
         FRIENDLYTO = trim(FRIENDLIES_NGRAUPEL),                       &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,        &
         DEFAULT = 0.0,   RC=STATUS  )  
    
        VERIFY_(STATUS)


    call MAPL_AddInternalSpec(GC,                                  &
         SHORT_NAME = 'QRAIN',                                     &
         LONG_NAME  = 'mass_fraction_of_rain',                     & 
         UNITS      = 'kg kg-1',                                   &
         FRIENDLYTO = trim(FRIENDLIES_QRAIN),                       &
         default    = 0.0,                                         &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                                  &
         SHORT_NAME = 'QSNOW',                                     &
         LONG_NAME  = 'mass_fraction_of_snow',                     &
         UNITS      = 'kg kg-1',                                   &
         FRIENDLYTO = trim(FRIENDLIES_QSNOW),                       &
         default    = 0.0,                                         &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                                  &
         SHORT_NAME = 'QGRAUPEL',                                  &
         LONG_NAME  = 'mass_fraction_of_graupel',                  &
         UNITS      = 'kg kg-1',                                   &
         FRIENDLYTO = trim(FRIENDLIES_QGRAUPEL),                       &
         default    = 0.0,                                         &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                                  &
         SHORT_NAME = 'QW',                                        &
         LONG_NAME  = 'mass_fraction_of_wet_air',                  &
         UNITS      = 'kg kg-1',                                   &
         RESTART    = MAPL_RestartSkip,                            &
         FRIENDLYTO = 'TURBULENCE:' // trim(COMP_NAME),            &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)                                                                          

    call MAPL_AddInternalSpec(GC,                                  &
         SHORT_NAME ='PDF_A',                                       &
          LONG_NAME = 'SHOC_PDF_relative_area_fraction',            &
         UNITS      ='1',                                           &       
         DIMS      = MAPL_DimsHorzVert,                             &
         VLOCATION = MAPL_VLocationCenter,                          &
         DEFAULT= 0.5,                                              &
         RC=STATUS  )  
    VERIFY_(STATUS)                         

     
    
  if (DOSHLW /= 0) then
          
    call MAPL_AddInternalSpec(GC,                                  &
         SHORT_NAME ='CUSH',                                       &
          LONG_NAME = 'Cumulus_scale_height_from_UW_shlw_convection',&
         UNITS      ='m',                                           &       
         DIMS      = MAPL_DimsHorzOnly,                             &
         VLOCATION = MAPL_VLocationNone,                            &
         DEFAULT= 1000.0,                                           &
         RC=STATUS  )  
         VERIFY_(STATUS)                         
    
  end if 
  
    ! !IMPORT STATE:

    call MAPL_AddImportSpec(GC,                              &
         SHORT_NAME = 'PLE',                                         &
         LONG_NAME  = 'air_pressure',                                &
         UNITS      = 'Pa',                                          &
         DIMS       = MAPL_DimsHorzVert,                            &
         VLOCATION  = MAPL_VLocationEdge,                           &
         AVERAGING_INTERVAL = AVRGNINT,                             &
         REFRESH_INTERVAL   = RFRSHINT,                             &
         RC=STATUS  )
    VERIFY_(STATUS)                                                                          

    call MAPL_AddImportSpec(GC,                              &
         SHORT_NAME = 'PKE',                                         &
         LONG_NAME  = 'edge_p$^\kappa$',                      &
         UNITS      = 'Pa$^\kappa$',                               &
         DIMS       = MAPL_DimsHorzVert,                            &
         VLOCATION  = MAPL_VLocationEdge,                           &
         AVERAGING_INTERVAL = AVRGNINT,                             &
         REFRESH_INTERVAL   = RFRSHINT,                             &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                              &
         SHORT_NAME = 'PLK',                                       &
         LONG_NAME  = 'mid-layer_p$^\kappa$',                      &
         UNITS      = 'Pa$^\kappa$',                               &
         DIMS       = MAPL_DimsHorzVert,                            &
         VLOCATION  = MAPL_VLocationCenter,                           &
         AVERAGING_INTERVAL = AVRGNINT,                             &
         REFRESH_INTERVAL   = RFRSHINT,                             &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                              &
         SHORT_NAME = 'PREF',                                       &
         LONG_NAME  = 'reference_air_pressure',                     &
         UNITS      = 'Pa',                                         &
         DIMS       = MAPL_DimsVertOnly,                            &
         VLOCATION  = MAPL_VLocationEdge,                           &
         AVERAGING_INTERVAL = AVRGNINT,                             &
         REFRESH_INTERVAL   = RFRSHINT,                             &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                              &
         SHORT_NAME = 'KH',                                         &
         LONG_NAME  = 'scalar_diffusivity',                         &
         UNITS      = 'm+2 s-1',                                    &
         DIMS       = MAPL_DimsHorzVert,                            &
         VLOCATION  = MAPL_VLocationEdge,                           &
         AVERAGING_INTERVAL = AVRGNINT,                             &
         REFRESH_INTERVAL   = RFRSHINT,                             &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                              &
         SHORT_NAME = 'TKE',                                        &
         LONG_NAME  = 'turbulent_kinetic_energy',                   &
         UNITS      = 'm+2 s-2',                                    &
         DIMS       = MAPL_DimsHorzVert,                            &
         VLOCATION  = MAPL_VLocationEdge,                           &
         AVERAGING_INTERVAL = AVRGNINT,                             &
         REFRESH_INTERVAL   = RFRSHINT,                             &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                                              &
         SHORT_NAME = 'WQT',                                                   &
         LONG_NAME  = 'Total_water_flux',                                      &
         UNITS      = 'kg kg-1 m s-1',                                               &
         DIMS       = MAPL_DimsHorzVert,                                       &
         VLOCATION  = MAPL_VLocationCenter,                                    &
         AVERAGING_INTERVAL = AVRGNINT,                             &
         REFRESH_INTERVAL   = RFRSHINT,                             &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                                              &
         SHORT_NAME = 'WHL',                                                   &
         LONG_NAME  = 'Liquid_water_static_energy_flux',                       &
         UNITS      = 'K m s-1',                                               &
         DIMS       = MAPL_DimsHorzVert,                                       &
         VLOCATION  = MAPL_VLocationCenter,                                    &
         AVERAGING_INTERVAL = AVRGNINT,                             &
         REFRESH_INTERVAL   = RFRSHINT,                             &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

       call MAPL_AddImportSpec(GC,                                  &
            SHORT_NAME = 'W2',                                       &
            LONG_NAME  = 'variance_of_vertical_velocity',             &
            UNITS      = 'm2 s-2',                                       &
            DIMS       = MAPL_DimsHorzVert,                           &
            VLOCATION  = MAPL_VLocationCenter,                          &
            AVERAGING_INTERVAL = AVRGNINT,                            &
            REFRESH_INTERVAL   = RFRSHINT,                            &
            RC=STATUS )
       VERIFY_(STATUS)

       call MAPL_AddImportSpec(GC,                                  &
            SHORT_NAME = 'W3',                                       &
            LONG_NAME  = 'third_moment_of_vertical_velocity',         &
            UNITS      = 'm3 s-3',                                     &
            DIMS       = MAPL_DimsHorzVert,                           &
            VLOCATION  = MAPL_VLocationCenter,                          &
            AVERAGING_INTERVAL = AVRGNINT,                            &
            REFRESH_INTERVAL   = RFRSHINT,                            &
            RC=STATUS )
       VERIFY_(STATUS)

       call MAPL_AddImportSpec(GC,                                  &
            SHORT_NAME = 'HL3',                                       &
            LONG_NAME  = 'third_moment_of_liquid_water_static_energy',    &
            UNITS      = 'K+3',                                       &
            DIMS       = MAPL_DimsHorzVert,                           &
            VLOCATION  = MAPL_VLocationCenter,                          &
            AVERAGING_INTERVAL = AVRGNINT,                            &
            REFRESH_INTERVAL   = RFRSHINT,                            &
            RC=STATUS )
       VERIFY_(STATUS)

       call MAPL_AddImportSpec(GC,                                  &
            SHORT_NAME = 'EDMF_FRC',                                       &
            LONG_NAME  = 'Mass_Flux_Fractional_Area',                 &
            UNITS      = '1',                                         &
            DIMS       = MAPL_DimsHorzVert,                           &
            VLOCATION  = MAPL_VLocationCenter,                          &
            AVERAGING_INTERVAL = AVRGNINT,                            &
            REFRESH_INTERVAL   = RFRSHINT,                            &
            RC=STATUS )
       VERIFY_(STATUS)

       call MAPL_AddImportSpec(GC,                                  &
            SHORT_NAME = 'HL2',                                       &
            LONG_NAME  = 'variance_of_liquid_water_static_energy',    &
            UNITS      = 'K+2',                                       &
            DIMS       = MAPL_DimsHorzVert,                           &
            VLOCATION  = MAPL_VLocationCenter,                          &
            AVERAGING_INTERVAL = AVRGNINT,                            &
            REFRESH_INTERVAL   = RFRSHINT,                            &
            RC=STATUS )
       VERIFY_(STATUS)

       call MAPL_AddImportSpec(GC,                                  &
            SHORT_NAME = 'QT2',                                       &
            LONG_NAME  = 'variance_of_total_water_specific_humidity', &
            UNITS      = '1',                                         &
            DIMS       = MAPL_DimsHorzVert,                           &
            VLOCATION  = MAPL_VLocationCenter,                          &
            AVERAGING_INTERVAL = AVRGNINT,                            &
            REFRESH_INTERVAL   = RFRSHINT,                            &
            RC=STATUS )
       VERIFY_(STATUS)

       call MAPL_AddImportSpec(GC,                                  &
            SHORT_NAME = 'QT3',                                       &
            LONG_NAME  = 'third_moment_of_total_water_specific_humidity', &
            UNITS      = '1',                                         &
            DIMS       = MAPL_DimsHorzVert,                           &
            VLOCATION  = MAPL_VLocationCenter,                          &
            AVERAGING_INTERVAL = AVRGNINT,                            &
            REFRESH_INTERVAL   = RFRSHINT,                            &
            RC=STATUS )
       VERIFY_(STATUS)

       call MAPL_AddImportSpec(GC,                                  &
            SHORT_NAME = 'HLQT',                                      &
            LONG_NAME  = 'covariance_of_liquid_water_static_energy_and_total_water_specific_humidity', &
            UNITS      = 'K',                                         &
            DIMS       = MAPL_DimsHorzVert,                           &
            VLOCATION  = MAPL_VLocationCenter,                          &
            AVERAGING_INTERVAL = AVRGNINT,                            &
            REFRESH_INTERVAL   = RFRSHINT,                            &
            RC=STATUS )
       VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                             &
         SHORT_NAME = 'TH',                                        &
         LONG_NAME  = 'potential_temperature',                     &
         UNITS      = 'K',                                         &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,                        &
         AVERAGING_INTERVAL = AVRGNINT,                            &
         REFRESH_INTERVAL   = RFRSHINT,                            &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                             &
         SHORT_NAME = 'U',                                         &
         LONG_NAME  = 'eastward_wind',                             &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,                        &
         AVERAGING_INTERVAL = AVRGNINT,                            &
         REFRESH_INTERVAL   = RFRSHINT,                            &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                             &
         SHORT_NAME = 'V',                                         &
         LONG_NAME  = 'northward_wind',                            &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,                        &
         AVERAGING_INTERVAL = AVRGNINT,                            &
         REFRESH_INTERVAL   = RFRSHINT,                            &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                             &
         SHORT_NAME = 'W',                                         &
         LONG_NAME  = 'vertical_velocity',                         &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,                        &
         AVERAGING_INTERVAL = AVRGNINT,                            &
         REFRESH_INTERVAL   = RFRSHINT,                            &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                             &
         SHORT_NAME = 'TS',                                        &
         LONG_NAME  = 'surface temperature',                       &
         UNITS      = 'K',                                         &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,                          &
         AVERAGING_INTERVAL = AVRGNINT,                            &
         REFRESH_INTERVAL   = RFRSHINT,                            &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                             &
         SHORT_NAME = 'SNOMAS',                                    &
         LONG_NAME  = 'snow_mass',                       &
         UNITS      = 'kg/m2',                                     &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,                        &
         AVERAGING_INTERVAL = AVRGNINT,                            &
         REFRESH_INTERVAL   = RFRSHINT,                            &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                             &
         SHORT_NAME = 'FRLANDICE',                                    &
         LONG_NAME  = 'areal_landice_fraction',                       &
         UNITS      = '1',                                     &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,                        &
         AVERAGING_INTERVAL = AVRGNINT,                            &
         REFRESH_INTERVAL   = RFRSHINT,                            &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                             &
         SHORT_NAME = 'FRLAND',                                    &
         LONG_NAME  = 'areal_land_fraction',                       &
         UNITS      = '1',                                     &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,                        &
         AVERAGING_INTERVAL = AVRGNINT,                            &
         REFRESH_INTERVAL   = RFRSHINT,                            &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                                    &
         SHORT_NAME = 'FRACI',                                     &
         LONG_NAME  = 'ice_covered_fraction_of_tile',              &
         UNITS      = '1',                                         &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,                          &
         AVERAGING_INTERVAL = AVRGNINT,                            &
         REFRESH_INTERVAL   = RFRSHINT,                            &
         RC=STATUS  )
    VERIFY_(STATUS)

    ! These bundles should be changed when we merge w/ the head.

    call MAPL_AddImportSpec(GC,                             &
         SHORT_NAME = 'MTR',                                        &
         LONG_NAME  = 'tracers_for_moist',                          &
         UNITS      = 'X',                                          &
         DATATYPE   = MAPL_BundleItem,                             &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,                        &
         RESTART    = MAPL_RestartSkip,                            &
         RC=STATUS )
    VERIFY_(STATUS)                                                                           

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME = 'KPBL',                                       &
        LONG_NAME  = 'planetary_boundary_layer_level',             &
        UNITS      = '1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                            &
        VLOCATION  = MAPL_VLocationNone,                           &
        AVERAGING_INTERVAL = AVRGNINT,                             &
        REFRESH_INTERVAL   = RFRSHINT,                             &
                                                        RC=STATUS  )
     VERIFY_(STATUS)      

     call MAPL_AddImportSpec(GC,                                   &
        SHORT_NAME = 'KPBL_SC',                                    &
        LONG_NAME  = 'boundary_layer_level_for_UW_shlw',           &
        UNITS      = '1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                            &
        VLOCATION  = MAPL_VLocationNone,                           &
        AVERAGING_INTERVAL = AVRGNINT,                             &
        REFRESH_INTERVAL   = RFRSHINT,                             &
        DEFAULT    = 72.,                                          &
                                                        RC=STATUS  )
     VERIFY_(STATUS)      

    call MAPL_AddImportSpec(GC,                                     &
         LONG_NAME  = 'aerosols',                                   &
         UNITS      = '1',                                          &
         SHORT_NAME = 'AERO',                                   &
         DIMS       = MAPL_DimsHorzVert,                            &
         VLOCATION  = MAPL_VLocationCenter,                         &
         DATATYPE   = MAPL_StateItem,                               &
         RC=STATUS  )
    VERIFY_(STATUS)      


    !new imports required for Aer-Cloud Interactions

    call MAPL_AddImportSpec(GC,                             &
         SHORT_NAME = 'TAUGWX',                                    &
         LONG_NAME  = 'surface_eastward_gravity_wave_stress',      &
         UNITS      = 'N m-2',                                     &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
    VERIFY_(STATUS)   

    call MAPL_AddImportSpec(GC,                             &
         SHORT_NAME = 'TAUGWY',                                    &
         LONG_NAME  = 'surface_northward_gravity_wave_stress',     &
         UNITS      = 'N m-2',                                     &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
    VERIFY_(STATUS)   


    call MAPL_AddImportSpec(GC,                             &
         LONG_NAME          = 'eastward_surface_stress_on_air',    &
         UNITS              = 'N m-2',                             &
         SHORT_NAME         = 'TAUX',                              &
         DIMS               = MAPL_DimsHorzOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                             &
         LONG_NAME          = 'northward_surface_stress_on_air',   &
         UNITS              = 'N m-2',                             &
         SHORT_NAME         = 'TAUY',                              &
         DIMS               = MAPL_DimsHorzOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
         RC=STATUS  )

    call MAPL_AddImportSpec(GC,                             &
         SHORT_NAME = 'TAUOROX',                                   &
         LONG_NAME  = 'surface_eastward_orographic_gravity_wave_stress',      &
         UNITS      = 'N m-2',                                     &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                             &
         SHORT_NAME = 'TAUOROY',                                   &
         LONG_NAME  = 'surface_northward_orographic_gravity_wave_stress',     &
         UNITS      = 'N m-2',                                     &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
    VERIFY_(STATUS)


    call MAPL_AddImportSpec ( gc,                                  &
         SHORT_NAME = 'OMEGA',                                     &
         LONG_NAME  = 'vertical_pressure_velocity',                &
         UNITS      = 'Pa s-1',                                  &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                                              &
         LONG_NAME  = 'Blackadar_length_scale_for_scalars',                    &
         UNITS      = 'm',                                                     &
         SHORT_NAME = 'ALH',                                                   &
         DIMS       = MAPL_DimsHorzVert,                                       &
         VLOCATION  = MAPL_VLocationEdge,                                      &
         RC=STATUS  )                                                   
    VERIFY_(STATUS)

    call MAPL_AddImportSpec ( GC,                                   &
         SHORT_NAME = 'RADLW',                                           &
         LONG_NAME  = 'air_temperature_tendency_due_to_longwave',        &
         UNITS      = 'K s-1',                                           &
         DIMS       = MAPL_DimsHorzVert,                                 &
         VLOCATION  = MAPL_VLocationCenter,                              &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec ( GC,                                   &
         SHORT_NAME = 'RADSW',                                           &
         LONG_NAME  = 'air_temperature_tendency_due_to_longwave',        &
         UNITS      = 'K s-1',                                           &
         DIMS       = MAPL_DimsHorzVert,                                 &
         VLOCATION  = MAPL_VLocationCenter,                              &
         RC=STATUS  )
    VERIFY_(STATUS)

    !============================================

   call MAPL_AddImportSpec ( GC,                                          & !USe the nature run to force cirrus
         SHORT_NAME = 'WSUB_NATURE',                                 &
         LONG_NAME  =  'variance in wsub from the nature run',     &
         UNITS      = 'm2 s-2',                                    &
         RESTART    = MAPL_RestartSkip,                            &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                             &
         SHORT_NAME = 'TROPP',                                          &
         LONG_NAME  = 'tropopause_pressure_based_on_blended_estimate',  &
         UNITS      = 'Pa',                                             &
         DIMS       = MAPL_DimsHorzOnly,                                &
         VLOCATION  = MAPL_VLocationNone,                               &
         AVERAGING_INTERVAL = AVRGNINT,                                 &
         REFRESH_INTERVAL   = RFRSHINT,                                 &
         RC=STATUS  )
    VERIFY_(STATUS)

    if ( IQVAINC /=0 ) then
       ! The following import is only for offline purposes
       ! NOTE: This is only used in offline applications so when adding new 
       !       fields to IMPORT state the suggestion is to add them BEFORE
       !       this state - unlike the usual procedure of always appending
       !       to the end of the state.
       call MAPL_AddImportSpec(GC,                                    &
            SHORT_NAME = 'QVAINC',                                    &
            LONG_NAME  = 'specific_humidity_analysis_increment',      &
            UNITS      = 'kg kg-1',                                   &
            default    = 0.0,                                         &
            DIMS       = MAPL_DimsHorzVert,                           &
            VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
       VERIFY_(STATUS)                                                                          
    endif

   call MAPL_AddImportSpec(GC,                                    &
         SHORT_NAME = 'T',                                         &
         LONG_NAME  = 'air_temperature',                           &
         UNITS      = 'K',                                         &
         DIMS       =  MAPL_DimsHorzVert,                          &
         VLOCATION  =  MAPL_VLocationCenter,                       &
                                                        RC=STATUS  )
   VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                                    &
         SHORT_NAME = 'DTDTDYN',                                     &
         LONG_NAME  = 'tendency_of_air_temperature_due_to_dynamics', &
         UNITS      = 'K s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                             &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS    )
     VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                                      &
         SHORT_NAME = 'DQVDTDYN',                                      &
         LONG_NAME  = 'tendency_of_specific_humidity_due_to_dynamics', &
         UNITS      = 'kg/kg/s',                                     &
         DIMS       = MAPL_DimsHorzVert,                               &
         VLOCATION  = MAPL_VLocationCenter,                 RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                                 &
         SHORT_NAME = 'QV_DYN_IN',                                 &
         LONG_NAME  = 'spec_humidity_at_begin_of_time_step',       &
         UNITS      = 'kg kg-1',                                   &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                                 &
         SHORT_NAME = 'T_DYN_IN',                                 &
         LONG_NAME  = 'temperature_at_begin_of_time_step',       &
         UNITS      = 'K',                                   &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                                 &
         SHORT_NAME = 'U_DYN_IN',                                 &
         LONG_NAME  = 'u_wind_at_begin_of_time_step',       &
         UNITS      = 'm s-1',                                   &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                                 &
         SHORT_NAME = 'V_DYN_IN',                                 &
         LONG_NAME  = 'v_wind_at_begin_of_time_step',       &
         UNITS      = 'm s-1',                                   &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                                 &
         SHORT_NAME = 'PLE_DYN_IN',                                 &
         LONG_NAME  = 'edge_pressure_at_begin_of_time_step',       &
         UNITS      = 'Pa',                                   &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationEdge,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                                   &
        SHORT_NAME         = 'AREA',                              &
        LONG_NAME          = 'agrid_cell_area',            &
        UNITS              = 'm+2',                             &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                                   &
        SHORT_NAME         = 'USTAR',                             &
        LONG_NAME          = 'surface_velocity_scale',            &
        UNITS              = 'm s-1',                             &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'TSTAR',                             &
        LONG_NAME          = 'surface_temperature_scale',         &
        UNITS              = 'K',                                 &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'QSTAR',                             &
        LONG_NAME          = 'surface_moisture_scale',            &
        UNITS              = 'kg kg-1',                           &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                                         &
       LONG_NAME  = '2-meter_air_temperature',                               &
       UNITS      = 'K',                                                     &
       SHORT_NAME = 'T2M',                                                   &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                                         &
       LONG_NAME  = '2-meter_specific_humidity',                             &
       UNITS      = 'kg kg-1',                                               &
       SHORT_NAME = 'Q2M',                                                   &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                                              &
       LONG_NAME  = 'surface_air_temperature',                               &
       UNITS      = 'K',                                                     &
       SHORT_NAME = 'TA',                                                    &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                                              &
       LONG_NAME  = 'surface_air_specific_humidity',                         &
       UNITS      = 'kg kg-1',                                               &
       SHORT_NAME = 'QA',                                                    &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)
    
    call MAPL_AddImportSpec(GC,                                   &
        LONG_NAME          = 'sensible_heat_flux',                &
        UNITS              = 'W m-2',                             &
        SHORT_NAME         = 'SH',                                &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'evaporation',                       &
        UNITS              = 'kg m-2 s-1',                        &
        SHORT_NAME         = 'EVAP',                              &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'surface geopotential height',       &
        UNITS              = 'm+2 s-2',                           &
        SHORT_NAME         = 'PHIS',                              &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
    VERIFY_(STATUS)
!-srf-gf-scheme

    ! !EXPORT STATE:

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'SHOC_PDF_relative_area_fraction',                       &
       UNITS      = '1',                                                     &
       SHORT_NAME = 'PDF_A',                                                 &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

#ifdef PDFDIAG
    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'SHOC_PDF_vertical_velocity_standard_deviation_first_plume', &
       UNITS      = '1',                                                     &
       SHORT_NAME = 'PDF_SIGW1',                                             &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'SHOC_PDF_vertical_velocity_standard_deviation_second_plume', &
       UNITS      = '1',                                                     &
       SHORT_NAME = 'PDF_SIGW2',                                             &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'SHOC_PDF_avg_vertical_velocity_of_first_plume',         &
       UNITS      = 'm s-1',                                                 &
       SHORT_NAME = 'PDF_W1',                                                &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'SHOC_PDF_avg_vertical_velocity_of_second_plume',        &
       UNITS      = 'm s-1',                                                 &
       SHORT_NAME = 'PDF_W2',                                                &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'SHOC_PDF_stddev_liq_wat_pot_temp_of_first_plume',       &
       UNITS      = 'K',                                                     &
       SHORT_NAME = 'PDF_SIGTH1',                                            &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'SHOC_PDF_stddev_liq_wat_pot_temp_of_second_plume',      &
       UNITS      = 'K',                                                     &
       SHORT_NAME = 'PDF_SIGTH2',                                            &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'SHOC_PDF_avg_liq_wat_pot_temp_of_first_plume',          &
       UNITS      = 'K',                                                     &
       SHORT_NAME = 'PDF_TH1',                                               &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'SHOC_PDF_avg_liq_wat_pot_temp_of_second_plume',         &
       UNITS      = 'K',                                                     &
       SHORT_NAME = 'PDF_TH2',                                               &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'SHOC_PDF_stddev_total_water_of_first_plume',            &
       UNITS      = 'kg kg-1',                                               &
       SHORT_NAME = 'PDF_SIGQT1',                                            &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'SHOC_PDF_stddev_total_water_of_second_plume',           &
       UNITS      = 'kg kg-1',                                               &
       SHORT_NAME = 'PDF_SIGQT2',                                            &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'SHOC_PDF_avg_total_water_of_first_plume',               &
       UNITS      = 'kg kg-1',                                               &
       SHORT_NAME = 'PDF_QT1',                                               &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'SHOC_PDF_avg_total_water_of_second_plume',              &
       UNITS      = 'kg kg-1',                                               &
       SHORT_NAME = 'PDF_QT2',                                               &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'SHOC_PDF_corr_total_water_liq_wat_pot_temp',            &
       UNITS      = '1',                                                     &
       SHORT_NAME = 'PDF_RQTTH',                                               &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'SHOC_PDF_corr_vertical_velocity_liq_wat_pot_temp',      &
       UNITS      = '1',                                                     &
       SHORT_NAME = 'PDF_RWTH',                                              &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'SHOC_PDF_corr_vertical_velocity_total_water',           &
       UNITS      = '1',                                                     &
       SHORT_NAME = 'PDF_RWQT',                                              &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)
#endif

    call MAPL_AddExportSpec(GC,                                              &
       SHORT_NAME = 'WTHV2',                                                 &
       LONG_NAME  = 'Buoyancy_flux_for_SHOC',                                &
       UNITS      = '1',                                                     &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)


    call MAPL_AddExportSpec(GC,                                              &
       SHORT_NAME = 'WQL',                                                   &
       LONG_NAME  = 'Liquid_water_flux',                                     &
       UNITS      = 'kg kg-1 m s-1',                                         &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)


    call MAPL_AddExportSpec(GC,                             &
         SHORT_NAME = 'QCTOT',                                      &
         LONG_NAME  = 'mass_fraction_of_total_cloud_water',         &
         UNITS      = 'kg kg-1',                                    &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                             &
         SHORT_NAME = 'QLTOT',                                      &
         LONG_NAME  = 'grid_box_mass_fraction_of_cloud_liquid_water',        &
         UNITS      = 'kg kg-1',                                    &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                             &
         SHORT_NAME = 'QITOT',                                      &
         LONG_NAME  = 'grid_box_mass_fraction_of_cloud_ice_water',           &
         UNITS      = 'kg kg-1',                                    &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                             &
         SHORT_NAME = 'QRTOT',                                      &
         LONG_NAME  = 'mass_fraction_of_falling_rain',              &
         UNITS      = 'kg kg-1',                                    &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                             &
         SHORT_NAME = 'QSTOT',                                      &
         LONG_NAME  = 'mass_fraction_of_falling_snow',              &
         UNITS      = 'kg kg-1',                                    &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                             &
         SHORT_NAME = 'QPTOTLS',                                    &
         LONG_NAME  = 'mass_fraction_of_large_scale_falling_precip', & 
         UNITS      = 'kg kg-1',                                    &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                              &
         SHORT_NAME = 'DTHDT ',                                     &
         LONG_NAME = 'pressure_weighted_potential_temperature_tendency_due_to_moist',&
         UNITS     = 'Pa K s-1',                                   &
         DIMS      = MAPL_DimsHorzVert,                           &
         VLOCATION = MAPL_VLocationCenter,                        &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                              &
         SHORT_NAME = 'DTDTFRIC',                                   &
         LONG_NAME = 'pressure_weighted_temperature_tendency_due_to_moist_friction',&
         UNITS     = 'Pa K s-1',                                   &
         DIMS      = MAPL_DimsHorzVert,                           &
         VLOCATION = MAPL_VLocationCenter,                        &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                              &
         SHORT_NAME = 'DQDT  ',                                     &
         LONG_NAME = 'specific_humidity_tendency_due_to_moist',    &
         UNITS     = 'kg kg-1 s-1',                                &
         DIMS      = MAPL_DimsHorzVert,                           &
         VLOCATION = MAPL_VLocationCenter,                        &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'DUDT  ',                                      &
         LONG_NAME = 'zonal_wind_tendency_due_to_moist',            &
         UNITS     = 'm s-2',                                       &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &                  
         SHORT_NAME = 'DVDT  ',                                      &
         LONG_NAME = 'meridional_wind_tendency_due_to_moist',       &
         UNITS     = 'm s-2',                                       &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'DWDT  ',                                      &
         LONG_NAME = 'vertical_velocity_tendency_due_to_moist',       &
         UNITS     = 'm s-2',                                       &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'DTHDTCN',                                     &
         LONG_NAME = 'potential_temperature_tendency_due_to_convection',&
         UNITS     = 'K s-1',                                       &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'DQDTCN',                                      &
         LONG_NAME = 'specific_humidity_tendency_due_to_convection',&
         UNITS     = 'kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'DQLDT ',                                      &
         LONG_NAME = 'total_liq_water_tendency_due_to_moist',       &
         UNITS     = 'kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME= 'DQIDT ',                                      &
         LONG_NAME = 'total_ice_water_tendency_due_to_moist',       &
         UNITS     = 'kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'DQCDTCN ',                                    &
         LONG_NAME = 'condensate_tendency_due_to_convection',       &
         UNITS     = 'kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME= 'CNV_DQLDT ',                                  &
         LONG_NAME = 'convective_condensate_source',                &
         UNITS     = 'kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME= 'CNV_PRC3 ',                                   &
         LONG_NAME = 'convective_precipitation_from_RAS',           &
         UNITS     = 'kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'DQRL   ',                                     &
         LONG_NAME = 'large_scale_rainwater_source',                &
         UNITS     = 'kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME= 'DQRC   ',                                     &
         LONG_NAME = 'net_convective_rainwater_source',             &
         UNITS     = 'kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
         RC=STATUS  )
    VERIFY_(STATUS)


    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME= 'CNV_MF0',                                     &
         LONG_NAME = 'cloud_base_mass_flux',                        &
         UNITS     = 'kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               & 
         SHORT_NAME = 'CNV_MFD',                                     & 
         LONG_NAME = 'detraining_mass_flux',                        &
         UNITS     = 'kg m-2 s-1',                                  &    
         DIMS      = MAPL_DimsHorzVert,                            &  
         VLOCATION = MAPL_VLocationCenter,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &                  
         SHORT_NAME = 'CNV_MFC',                                     & 
         LONG_NAME = 'cumulative_mass_flux',                        &
         UNITS     = 'kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationEdge,                           &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME = 'CNV_FREQ',                                    & 
         LONG_NAME = 'convective_frequency',                        &
         UNITS     = 'fraction',                                    &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                           &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME = 'CNV_BASEP',                                   & 
         LONG_NAME = 'pressure_at_convective_cloud_base',           &
         UNITS     = 'Pa',                                          &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                           &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME = 'CNV_TOPP',                                    & 
         LONG_NAME = 'pressure_at_convective_cloud_top',            &
         UNITS     = 'Pa',                                          &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                           &
         RC=STATUS  )
    VERIFY_(STATUS)



    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'CNV_UPDF',                                    &
         LONG_NAME = 'updraft_areal_fraction',                      &
         UNITS     = '1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
         RC=STATUS  )

    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'CNV_CVW',                                     &
         LONG_NAME = 'updraft_vertical_velocity',                   &
         UNITS     = 'hPa s-1',                                     &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
         RC=STATUS  )

    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME ='CNV_QC',                                      &
         LONG_NAME ='grid_mean_convective_condensate',             &
         UNITS     ='kg kg-1',                                     &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'CLDBASEHGT',                                &
         LONG_NAME = 'Height_of_cloud_base',                       &
         UNITS     = 'm',                                          &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                           &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME= 'SHLW_PRC3 ',                                   &
         LONG_NAME = 'shallow_convective_rain',           &
         UNITS     = 'kg kg-1 s-1',                                  &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME= 'SHLW_SNO3 ',                                   &
         LONG_NAME = 'shallow_convective_snow',           &
         UNITS     = 'kg kg-1 s-1',                                  &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'CUFRC_SC',                                  &
         LONG_NAME  = 'shallow_cumulus_cloud_fraction',            &
         UNITS      = 'fraction',                                  &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,                        &
         RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'DQRDT_SC',                                  &
         LONG_NAME  = 'shallow_cumulus_precipitating_condensate',            &
         UNITS      = 'kg kg-1 s-1',                                  &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,                        &
         RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'DQCDT_SC',                             &
         LONG_NAME  = 'shallow_cumulus_condensate_tendency',            &
         UNITS      = 'kg kg-1 s-1',                               &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,                        &
         RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'DQSDT_SC',                                  &
         LONG_NAME  = 'shallow_cumulus_precipitating_frozen_condensate',            &
         UNITS      = 'kg kg-1 s-1',                                  &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,                        &
         RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'DTHDT_SC',                                  &
         LONG_NAME  = 'Potential_temperature_tendency_from_shallow_convection',            &
         UNITS      = 'K s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,                        &
         RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'DQIDT_SC',                                  &
         LONG_NAME  = 'Ice_tendency_from_shallow_convection',      &
         UNITS      = 'kg ks-1 s-1',                               &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,                        &
         RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'DQLDT_SC',                                      &
         LONG_NAME = 'Liquid_water_tendency_from_shallow_convection',&
         UNITS     = 'kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'DQVDT_SC',                                      &
         LONG_NAME = 'Specific_humidity_tendency_from_shallow_convection',&
         UNITS     = 'kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'DUDT_SC',                                      &
         LONG_NAME = 'Zonal_wind_tendency_from_shallow_convection',&
         UNITS     = 'm s-2',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'DVDT_SC',                                      &
         LONG_NAME = 'Meridional_wind_tendency_from_shallow_convection',&
         UNITS     = 'm s-2',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'WLCL_SC',                                      &
         LONG_NAME = 'Vertical_velocity_at_LCL_for_shallow_convection',&
         UNITS     = 'm s-1',                                 &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'QTSRC_SC',                                      &
         LONG_NAME = 'Total_water_in_source_air_for_shallow_convection',&
         UNITS     = 'kg kg-1',                                 &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'THLSRC_SC',                                      &
         LONG_NAME = 'Liquid_potential_temperature_of_source_air_for_shallow_convection',&
         UNITS     = 'K',                                 &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'THVLSRC_SC',                                      &
         LONG_NAME = 'Liquid_virtual_potential_temperature_source_air_shallow_convection',&
         UNITS     = 'K',                                 &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'TKEAVG_SC',                                      &
         LONG_NAME = 'Average_boundary_layer_TKE_used_for_shallow_convection',&
         UNITS     = 'm2 s-2',                                 &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'CLDTOP_SC',                                      &
         LONG_NAME = 'Cloud_top_height_from_shallow_convection',&
         UNITS     = 'm',                                 &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                         &
         RC=STATUS  )
    VERIFY_(STATUS)


    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'WUP_SC',                                      &
         LONG_NAME = 'Vertical_velocity_in_shallow_convection_updraft',&
         UNITS     = 'm s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationEdge,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'QTUP_SC',                                      &
         LONG_NAME = 'Total_water_in_shallow_convection_updraft',&
         UNITS     = 'kg kg-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationEdge,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'THLUP_SC',                                      &
         LONG_NAME = 'Liquid_potential_temperature_in_shallow_convection_updraft',&
         UNITS     = 'K',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationEdge,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'THVUP_SC',                                      &
         LONG_NAME = 'Virtual_potential_temperature_in_shallow_convection_updraft',&
         UNITS     = 'K',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationEdge,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'UUP_SC',                                      &
         LONG_NAME = 'Zonal_wind_in_shallow_convection_updraft',&
         UNITS     = 'm s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationEdge,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'VUP_SC',                                      &
         LONG_NAME = 'Meridional_wind_in_shallow_convection_updraft',&
         UNITS     = 'm s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationEdge,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'ENTR_SC',                                      &
         LONG_NAME = 'Lateral_entrainment_rate_in_shallow_convection_updraft',&
         UNITS     = 'Pa-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'DETR_SC',                                      &
         LONG_NAME = 'Lateral_detrainment_rate_in_shallow_convection_updraft',&
         UNITS     = 'Pa-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'XC_SC',                                      &
         LONG_NAME = 'Critical_mixing_fraction_in_shallow_updraft',&
         UNITS     = 'Pa-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'QLDET_SC',                                      &
         LONG_NAME = 'Detrained_liquid_water_from_shallow_convection',&
         UNITS     = 'kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'QIDET_SC',                                      &
         LONG_NAME = 'Detrained_ice_water_from_shallow_convection',&
         UNITS     = 'kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'QLENT_SC',                                      &
         LONG_NAME = 'Sink_from_entrained_liquid_water_by_shallow_convection',&
         UNITS     = 'kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'QIENT_SC',                                      &
         LONG_NAME = 'Sink_from_entrained_ice_water_by_shallow_convection',&
         UNITS     = 'kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'QLSUB_SC',                                      &
         LONG_NAME = 'Shallow_convective_subsidence_liquid_water_tendency',&
         UNITS     = 'kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'QISUB_SC',                                      &
         LONG_NAME = 'Shallow_convective_subsidence_ice_water_tendency',&
         UNITS     = 'kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'CNT_SC',                                      &
         LONG_NAME = 'Shallow_cloud_top_interface',            &
         UNITS     = 'index',                                 &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'CNB_SC',                                      &
         LONG_NAME = 'Shallow_cloud_bottom_interface',            &
         UNITS     = 'index',                                 &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'QCU_SC',                                      &
         LONG_NAME = 'Shallow_updraft_condensate',            &
         UNITS     = 'kg kg-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'QLU_SC',                                      &
         LONG_NAME = 'Shallow_updraft_liquid_condensate',            &
         UNITS     = 'kg kg-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'QIU_SC',                                      &
         LONG_NAME = 'Shallow_updraft_frozen_condensate',            &
         UNITS     = 'kg kg-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'UMF_SC',                                      &
         LONG_NAME = 'Shallow_updraft_mass_flux_at_interfaces', &
         UNITS     = 'kg m-2 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationEdge,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'MFD_SC',                                      &
         LONG_NAME = 'Shallow_updraft_detrained_mass_flux', &
         UNITS     = 'kg m-2 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'CIN_SC',                                      &
         LONG_NAME = 'Convective_inhibition_for_shallow_convection', &
         UNITS     = 'J kg-1',                                 &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'PLCL_SC',                                      &
         LONG_NAME = 'Pressure_at_lift_condensation_level',   &
         UNITS     = 'Pa',                                 &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'PLFC_SC',                                      &
         LONG_NAME = 'Pressure_at_level_of_free_convection',   &
         UNITS     = 'Pa',                                 &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'PINV_SC',                                      &
         LONG_NAME = 'Pressure_at_inversion_level',   &
         UNITS     = 'Pa',                                 &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'PREL_SC',                                      &
         LONG_NAME = 'Pressure_at_release_level',   &
         UNITS     = 'Pa',                                 &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'PBUP_SC',                                      &
         LONG_NAME = 'Pressure_at_level_neutral_buoyancy',   &
         UNITS     = 'Pa',                                 &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'CBMF_SC',                                      &
         LONG_NAME = 'cloud_base_mass_flux_due_to_shallow_convection',&
         UNITS     = 'kg s-1 m-2',                                 &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'CLDBASEHGT',                                &
         LONG_NAME = 'Height_of_cloud_base',                       &
         UNITS     = 'm',                                          &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                           &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME ='RL',                                          & 
         LONG_NAME ='liquid_cloud_particle_effective_radius',      &
         UNITS     ='m',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'RI',                                          & 
         LONG_NAME = 'ice_phase_cloud_particle_effective_radius',   &
         UNITS     = 'm',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'RR',                                          & 
         LONG_NAME = 'falling_rain_particle_effective_radius',      &
         UNITS     = 'm',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'RS',                                          & 
         LONG_NAME  = 'falling_snow_particle_effective_radius',       &
         UNITS     = 'm',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

 call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'RG',                                          & 
         LONG_NAME  = 'falling_graupel_particle_effective_radius',       &
         UNITS     = 'm',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)
    
    
    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME ='CLDNCCN',                                     & 
         LONG_NAME ='number_concentration_of_cloud_particles',     &
         UNITS     ='m-3',                                         &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'QSATI'  ,                                     & 
         LONG_NAME = 'saturation_spec_hum_over_ice',                &
         UNITS     = 'kg kg-1',                                     &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'QSATL'  ,                                     & 
         LONG_NAME = 'saturation_spec_hum_over_liquid',             &
         UNITS     = 'kg kg-1',                                     &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'ALPHT'  ,                                     & 
         LONG_NAME = 'pdf_spread_for_condensation_over_qsat_total', &
         UNITS     = '1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='ALPH1'  ,                                     & 
         LONG_NAME ='pdf_spread_for_condensation_over_qsat_term1', &
         UNITS     ='1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'ALPH2'  ,                                     & 
         LONG_NAME = 'pdf_spread_for_condensation_over_qsat_term2', &
         UNITS     = '1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'CFPDFX'  ,                                    & 
         LONG_NAME = 'cloud_fraction_internal_in_PDF_scheme',       &
         UNITS     = '1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'RHCLR'  ,                                     & 
         LONG_NAME = 'RH_clear_sky',                                &
         UNITS     = '1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'CFPDF'  ,                                     & 
         LONG_NAME = 'cloud_fraction_after_PDF',                    &
         UNITS     = '1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'FCLD'  ,                                      & 
         LONG_NAME = 'cloud_fraction_for_radiation',                &
         UNITS     = '1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME ='QV',                                          & 
         LONG_NAME ='water_vapor_for_radiation',                   &
         UNITS     ='kg kg-1',                                     &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'QL',                                          & 
         LONG_NAME = 'in_cloud_cloud_liquid_for_radiation',                  &
         UNITS     = 'kg kg-1',                                     &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'QI',                                          & 
         LONG_NAME = 'in_cloud_cloud_ice_for_radiation',                     &
         UNITS     = 'kg kg-1',                                     &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'QR',                                          & 
         LONG_NAME = 'Falling_rain_for_radiation',                  &
         UNITS     = 'kg kg-1',                                     &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'QS',                                          & 
         LONG_NAME = 'Falling_snow_for_radiation',                  &
         UNITS     = 'kg kg-1',                                     &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'QG',                                          &
         LONG_NAME = 'Falling_graupel_for_radiation',                  &
         UNITS     = 'kg kg-1',                                     &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'DBZ',                                          &
         LONG_NAME = 'Simulated_radar_reflectivity',                  &
         UNITS     = 'dBZ',                                     &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'DBZ_MAX',                                          &
         LONG_NAME = 'Maximum_simulated_radar_reflectivity',                  &
         UNITS     = 'dBZ',                                     &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME ='PRCP_RAIN',                                     &
         LONG_NAME ='falling_rain_precipitation_at_surface',          &
         UNITS     ='kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME ='PRCP_SNOW',                                     &
         LONG_NAME ='falling_snow_precipitation_at_surface',          &
         UNITS     ='kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME ='PRCP_ICE',                                     &
         LONG_NAME ='falling_ice_precipitation_at_surface',          &
         UNITS     ='kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME ='PRCP_GRAUPEL',                                     &
         LONG_NAME ='falling_graupel_precipitation_at_surface',          &
         UNITS     ='kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME ='LS_PRCP',                                     & 
         LONG_NAME ='nonanvil_large_scale_precipitation',          &
         UNITS     ='kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'AN_PRCP',                                     & 
         LONG_NAME = 'anvil_precipitation',                         &
         UNITS     = 'kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME ='CN_PRCP',                                     & 
         LONG_NAME ='convective_precipitation',                    &
         UNITS     ='kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME ='SC_PRCP',                                    & 
         LONG_NAME ='shallow_convective_precipitation',            &
         UNITS     ='kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'ER_PRCP',                                     & 
         LONG_NAME = 'spurious_rain_from_RH_cleanup',          &
         UNITS     = 'kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'FILLNQV',                                     & 
         LONG_NAME = 'filling_of_negative_Q',          &
         UNITS     = 'kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'PGENTOT',                                     & 
         LONG_NAME = 'Total_column_production_of_precipitation',    &
         UNITS     = 'kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'PREVTOT',                                     & 
         LONG_NAME = 'Total_column_re-evap/subl_of_precipitation',    &
         UNITS     = 'kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'LS_ARF',                                      & 
         LONG_NAME = 'areal_fraction_of_nonanvil_large_scale_showers',&
         UNITS     = '1',                                           &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'AN_ARF',                                      & 
         LONG_NAME = 'areal_fraction_of_anvil_showers',             &
         UNITS     = '1',                                           &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'CN_ARF',                                      & 
         LONG_NAME = 'areal_fraction_of_convective_showers',        &
         UNITS     = '1',                                           &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)


    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'SC_ARF',                                      & 
         LONG_NAME = 'areal_fraction_of_shallow_showers',        &
         UNITS     = '1',                                           &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'PTYPE',                                         &
         LONG_NAME = 'surface_precipitation_type',                 &
         UNITS     = '1',                                  &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'FRZR',                                         &
         LONG_NAME = 'freezing_rain_fall',                                    &
         UNITS     = 'kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'ICE',                                         &
         LONG_NAME = 'icefall',                                    &
         UNITS     = 'kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'SNO',                                         & 
         LONG_NAME = 'snowfall',                                    &
         UNITS     = 'kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'PCU',                                         & 
         LONG_NAME = 'convective_rainfall',                         &
         UNITS     = 'kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'PLS',                                         & 
         LONG_NAME = 'large_scale_rainfall',                        &
         UNITS     = 'kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='TPREC',                                       & 
         LONG_NAME ='total_precipitation',                         &
         UNITS     ='kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='HOURNORAIN',                                  & 
         LONG_NAME ='time-during_an_hour_with_no_precipitation',   &
         UNITS     ='s',                                           &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='TPW',                                         & 
         LONG_NAME ='total_precipitable_water',                    &
         UNITS     ='kg m-2'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='CCWP',                                        & 
         LONG_NAME ='grid_mean_conv_cond_water_path_diagnostic',   &
         UNITS     ='kg m-2'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='CWP',                                         & 
         LONG_NAME ='condensed_water_path',                        &
         UNITS     ='kg m-2'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='LWP',                                         & 
         LONG_NAME ='liquid_water_path',                           &
         UNITS     ='kg m-2'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='IWP',                                         & 
         LONG_NAME ='ice_water_path',                              &
         UNITS     ='kg m-2'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='BYNCY',                                       & 
         LONG_NAME ='buoyancy_of surface_parcel',                  &
         UNITS     ='m s-2',                                       &
         DIMS      = MAPL_DimsHorzVert,                            & 
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='CAPE',                                        & 
         LONG_NAME ='cape_for_surface_parcel',                     &
         UNITS     ='J kg-1',                                      &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='INHB',                                        & 
         LONG_NAME ='inhibition_for_surface_parcel',               &
         UNITS     ='J kg-1',                                      &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='TVQ0',                                        & 
         LONG_NAME ='Total_Water_Substance_Before',                &
         UNITS     ='kg m-2'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='TVQ1',                                        & 
         LONG_NAME ='Total_Water_Substance_After',                 &
         UNITS     ='kg m-2'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='DCPTE',                                        & 
         LONG_NAME ='Total_VI_DcpT',                         &
         UNITS     ='J m-2'  ,                                     &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='TVE0',                                        & 
         LONG_NAME ='Total_VI_MSE_Before',                         &
         UNITS     ='J m-2'  ,                                     &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='TVE1',                                        & 
         LONG_NAME ='Total_VI_MSE_After',                          &
         UNITS     ='J m-2'  ,                                     &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='TVEX',                                        & 
         LONG_NAME ='Total_VI_MSE_Somewhere',                      &
         UNITS     ='J m-2'  ,                                     &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='ZPBLCN',                                      & 
         LONG_NAME ='boundary_layer_depth',                        &
         UNITS     ='m'   ,                                        &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='ZLCL',                                        & 
         LONG_NAME ='lifting_condensation_level',                  &
         UNITS     ='m'  ,                                         &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='ZLFC',                                        & 
         LONG_NAME ='level_of_free_convection',                    &
         UNITS     ='m'  ,                                         &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='ZCBL',                                        & 
         LONG_NAME ='height_of_cloud_base_layer',                  &
         UNITS     ='m'  ,                                         &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='MXDIAM',                                      & 
         LONG_NAME ='diameter_of_largest_RAS_plume',               &
         UNITS     ='m'  ,                                         &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='RH600',                                      &
         LONG_NAME ='rh_at_600mb',               &
         UNITS     ='kg/kg'  ,                                         &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='Q600',                                      &
         LONG_NAME ='q_at_600mb',               &
         UNITS     ='kg/kg'  ,                                         &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='QCBL',                                      &
         LONG_NAME ='q_at_cloud_base_level',               &
         UNITS     ='kg/kg'  ,                                         &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='CNV_FRC',                                      &
         LONG_NAME ='convective_fraction',               &
         UNITS     =''  ,                                         &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='STOCH_CNV',                                     &
         LONG_NAME ='stochastic_factor_for_convection',               &
         UNITS     =''  ,                                         &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='SIGMA_DEEP',                                     &
         LONG_NAME ='sigma_for_deep_in_convection',               &
         UNITS     =''  ,                                         &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='SIGMA_MID',                                     &
         LONG_NAME ='sigma_for_mid_in_convection',               &
         UNITS     =''  ,                                         &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='RAS_TIME',                                     & 
         LONG_NAME ='timescale_for_RAS_plumes',               &
         UNITS     ='s'  ,                                         &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='RAS_TRG',                                     &
         LONG_NAME ='rh_trigger_for_RAS_plumes',               &
         UNITS     =''  ,                                         &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='RAS_TOKI',                                     &
         LONG_NAME ='tokioka_factor_for_RAS_plumes',               &
         UNITS     =''  ,                                         &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='RAS_PBL',                                     &
         LONG_NAME ='pbl_fraction_for_RAS_plumes',               &
         UNITS     =''  ,                                         &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='RAS_WFN',                                     &
         LONG_NAME ='RAS_work_function_before_scaling',               &
         UNITS     =''  ,                                         &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='ENTLAM',                                      &
         LONG_NAME ='entrainment parameter',                       &
         UNITS     ='m-1',                                  &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

#if 0 
    ! taken out since they are now friendly to dynamics
    call MAPL_AddExportSpec(GC,                             &
         SHORT_NAME ='QLCN',                                       &
         LONG_NAME  ='mass_fraction_of_convective_cloud_liquid_water', &
         UNITS      ='1',                                          &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)                                                                          

    call MAPL_AddExportSpec(GC,                             &
         SHORT_NAME ='QICN',                                       &
         LONG_NAME  ='mass_fraction_of_convective_cloud_ice_water', &
         UNITS      ='1',                                          &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)                                                                          

    call MAPL_AddExportSpec(GC,                             &
         SHORT_NAME ='CLLS',                                       &
         LONG_NAME  ='large_scale_cloud_volume_fraction',            &
         UNITS      ='1',                                          &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)                                                                          

    call MAPL_AddExportSpec(GC,                             &
         SHORT_NAME ='CLCN',                                       &
         LONG_NAME  ='convective_cloud_volume_fraction',             &
         UNITS      ='1',                                          &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)                                                                          
#endif


    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='RH1',                                         & 
         LONG_NAME ='relative_humidity_before_moist',              &
         UNITS     ='1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='RH2',                                         & 
         LONG_NAME ='relative_humidity_after_moist',               &
         UNITS     ='1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    !Outputs to give model trajectory in the moist TLM/ADJ
    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='TH_moist',                                    & 
         LONG_NAME ='potential_temp_before_moist',                 &
         UNITS     ='K',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='Q_moist',                                     & 
         LONG_NAME ='specific_humidity_before_moist',              &
         UNITS     ='kg kg-1',                                     &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='SL',                                          & 
         LONG_NAME ='liquid_water_static_energy',                  &
         UNITS     ='J kg-1',                                      &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='QT',                                          & 
         LONG_NAME ='total_water_specific_humidity',               &
         UNITS     ='kg kg-1',                                     &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='T',                                           & 
         LONG_NAME ='temperature',                                 &
         UNITS     ='K',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='KCBL_moist',                                  & 
         LONG_NAME ='KCBL_before_moist',                           &
         UNITS     ='1',                                           &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='ctop_moist',                                  & 
         LONG_NAME ='ctop_after_ras',                              &
         UNITS     ='1',                                           &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='TS_moist',                                    & 
         LONG_NAME ='surface_temp_before_moist',                   &
         UNITS     ='K',                                           &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='KHu_moist',                                   &
         LONG_NAME ='upper_index_where_Kh_greater_than_2',         &
         UNITS     ='1',                                           &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='KHl_moist',                                   &
         LONG_NAME ='lower_index_where_Kh_greater_than_2',         &
         UNITS     ='1',                                           &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
    VERIFY_(STATUS)

    !End outputs for trajectory of moist TLM/ADJ

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='KHX',                                         &
         LONG_NAME ='scalar_diffusivity_for_pbl_clouds',               &
         UNITS     ='m+2 s-1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='DTSX',                                   &
         LONG_NAME ='dts_for_pbl_clouds',         &
         UNITS     ='1',                                           &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='RHX',                                         & 
         LONG_NAME ='relative_humidity_after_PDF',               &
         UNITS     ='1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='REVSU_CN',                                    & 
         LONG_NAME ='evap_subl_of_convective_precipitation',       &
         UNITS     ='kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='REVSU_LSAN',                                    & 
         LONG_NAME ='evap_subl_of_non_convective_precipitation',       &
         UNITS     ='kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='REV_CN',                                      & 
         LONG_NAME ='diagnostic_evaporation_of_convective_precipitation',     &
         UNITS     ='kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='REV_CN_GF',                                   & 
         LONG_NAME ='diagnostic_evaporation_of_convective_precipitation',     &
         UNITS     ='kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='REV_SC',                                      & 
         LONG_NAME ='evaporation_of_shallow_precipitation',     &
         UNITS     ='kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='REV_AN',                                      & 
         LONG_NAME ='evaporation_of_anvil_precipitation',          &
         UNITS     ='kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='REV_LS',                                      & 
         LONG_NAME ='evaporation_of_nonanvil_large_scale_precipitation',&
         UNITS     ='kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='RSU_CN',                                      & 
         LONG_NAME ='diagnostic_sublimation_of_convective_precipitation',     &
         UNITS     ='kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='RSU_CN_GF',                                   & 
         LONG_NAME ='diagnostic_sublimation_of_convective_precipitation',     &
         UNITS     ='kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='RSU_SC',                                      & 
         LONG_NAME ='sublimation_of_shallow_precipitation',     &
         UNITS     ='kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='RSU_AN',                                      & 
         LONG_NAME ='sublimation_of_anvil_precipitation',          &
         UNITS     ='kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='RSU_LS',                                      & 
         LONG_NAME ='sublimation_of_nonanvil_large_scale_precipitation',&
         UNITS     ='kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='ACR_TOT',                                      & 
         LONG_NAME ='total_accretion_of__precipitation',     &
         UNITS     ='kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='ACRLL_CN',                                      & 
         LONG_NAME ='liq_liq_accretion_of_convective_precipitation',     &
         UNITS     ='kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='ACRLL_SC',                                      & 
         LONG_NAME ='liq_liq_accretion_of_shallow_precipitation',     &
         UNITS     ='kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='ACRLL_AN',                                      & 
         LONG_NAME ='liq_liq_accretion_of_anvil_precipitation',          &
         UNITS     ='kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='ACRLL_LS',                                      & 
         LONG_NAME ='liq_liq_accretion_of_nonanvil_large_scale_precipitation',&
         UNITS     ='kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='ACRIL_CN',                                      & 
         LONG_NAME ='ice_liq_accretion_of_convective_precipitation',     &
         UNITS     ='kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='ACRIL_SC',                                      & 
         LONG_NAME ='ice_liq_accretion_of_shallow_convective_precipitation',     &
         UNITS     ='kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='ACRIL_AN',                                      & 
         LONG_NAME ='ice_liq_accretion_of_anvil_precipitation',          &
         UNITS     ='kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='ACRIL_LS',                                      & 
         LONG_NAME ='ice_liq_accretion_of_nonanvil_large_scale_precipitation',&
         UNITS     ='kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='PFI_CN',                                      & 
         LONG_NAME ='3D_flux_of_ice_convective_precipitation',     &
         UNITS     ='kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationEdge,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='PFI_CN_GF',                                   & 
         LONG_NAME ='3D_flux_of_ice_convective_precipitation',     &
         UNITS     ='kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationEdge,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='PFI_SC',                                      & 
         LONG_NAME ='3D_flux_of_ice_shallow_convective_precipitation',     &
         UNITS     ='kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationEdge,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='PFI_AN',                                      & 
         LONG_NAME ='3D_flux_of_ice_anvil_precipitation',          &
         UNITS     ='kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationEdge,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='PFI_LS',                                      & 
         LONG_NAME ='3D_flux_of_ice_nonanvil_large_scale_precipitation',&
         UNITS     ='kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationEdge,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='PFI_LSAN',                                      & 
         LONG_NAME ='3D_flux_of_ice_nonconvective_precipitation'  ,&
         UNITS     ='kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationEdge,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='PFL_CN',                                      & 
         LONG_NAME ='3D_flux_of_liquid_convective_precipitation',     &
         UNITS     ='kg m-2 s-1',                                   &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationEdge,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='PFL_CN_GF',                                    & 
         LONG_NAME ='3D_flux_of_liquid_convective_precipitation',   &
         UNITS     ='kg m-2 s-1',                                   &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationEdge,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='PFL_SC',                                      & 
         LONG_NAME ='3D_flux_of_liquid_shallow_convective_precipitation',     &
         UNITS     ='kg m-2 s-1',                                   &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationEdge,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='PFL_AN',                                      & 
         LONG_NAME ='3D_flux_of_liquid_anvil_precipitation',          &
         UNITS     ='kg m-2 s-1',                                         &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationEdge,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='PFL_LS',                                      & 
         LONG_NAME ='3D_flux_of_liquid_nonanvil_large_scale_precipitation',&
         UNITS     ='kg m-2 s-1',                                   &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationEdge,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='PFL_LSAN',                                      & 
         LONG_NAME ='3D_flux_of_liquid_nonconvective_precipitation'  ,&
         UNITS     ='kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationEdge,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                     &
         SHORT_NAME='DPDTMST',                                      & 
         LONG_NAME ='layer_pressure_thickness_tendency_from_moist', &
         UNITS     ='Pa s-1',                                       &
         DIMS      = MAPL_DimsHorzVert,                             &
         VLOCATION = MAPL_VLocationCenter,               RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='DCNVL',                                       & 
         LONG_NAME ='convective_source_of_cloud_liq',   &
         UNITS     ='kg kg-1 s-1',                                   &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='DCNVI',                                       & 
         LONG_NAME ='convective_source_of_cloud_ice',        &
         UNITS     ='kg kg-1 s-1',                                   &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='DLPDF',                                    & 
         LONG_NAME ='pdf_source_sink_of_cloud_liq',    &
         UNITS     ='kg kg-1 s-1',                                   &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='DIPDF',                                    & 
         LONG_NAME ='pdf_source_sink_of_cloud_ice',    &
         UNITS     ='kg kg-1 s-1',                                   &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='DLFIX',                                    & 
         LONG_NAME ='fix_source_sink_of_cloud_liq',          &
         UNITS     ='kg kg-1 s-1',                                   &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='DIFIX',                                    & 
         LONG_NAME ='fix_source_sink_of_cloud_ice',          &
         UNITS     ='kg kg-1 s-1',                                   &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='AUT',                                      & 
         LONG_NAME ='autoconv_sink_of_cloud_liq',               &
         UNITS     ='kg kg-1 s-1',                                   &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='EVAPC',                                  & 
         LONG_NAME ='evaporation_of_cloud_liq',               &
         UNITS     ='kg kg-1 s-1',                            &
         DIMS      = MAPL_DimsHorzVert,                       &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='SDM',                                    & 
         LONG_NAME ='sedimentation_sink_of_cloud_ice',        &
         UNITS     ='kg kg-1 s-1',                            &
         DIMS      = MAPL_DimsHorzVert,                       &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                 &
         SHORT_NAME='VFALLICE_AN',                              & 
         LONG_NAME ='autoconversion_fall_velocity_of_anvil_snow',       &
         UNITS     ='m s-1',                                    &
         DIMS      = MAPL_DimsHorzVert,                         &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                 &
         SHORT_NAME='VFALLICE_LS',                              & 
         LONG_NAME ='autoconversion_fall_velocity_of_largescale_snow',  &
         UNITS     ='m s-1',                                    &
         DIMS      = MAPL_DimsHorzVert,                         &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                 &
         SHORT_NAME='VFALLWAT_AN',                              & 
         LONG_NAME ='autoconversion_fall_velocity_of_anvil_rain',     &
         UNITS     ='m s-1',                                    &
         DIMS      = MAPL_DimsHorzVert,                         &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                 &
         SHORT_NAME='VFALLWAT_LS',                              & 
         LONG_NAME ='autoconversion_fall_velocity_of_largescale_rain',&
         UNITS     ='m s-1',                                    &
         DIMS      = MAPL_DimsHorzVert,                         &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                 &
         SHORT_NAME='VFALLRN_AN',                               & 
         LONG_NAME ='reevaporation_fall_velocity_of_anvil_rain',              &
         UNITS     ='m s-1',                                    &
         DIMS      = MAPL_DimsHorzVert,                         &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                 &
         SHORT_NAME='VFALLRN_LS',                               & 
         LONG_NAME ='reevaporation_fall_velocity_of_largescale_rain',         &
         UNITS     ='m s-1',                                    &
         DIMS      = MAPL_DimsHorzVert,                         &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                 &
         SHORT_NAME='VFALLRN_CN',                               & 
         LONG_NAME ='reevaporation_fall_velocity_of_convective_rain',         &
         UNITS     ='m s-1',                                    &
         DIMS      = MAPL_DimsHorzVert,                         &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                 &
         SHORT_NAME='VFALLRN_SC',                               & 
         LONG_NAME ='reevaporation_fall_velocity_of_shallow_rain',         &
         UNITS     ='m s-1',                                    &
         DIMS      = MAPL_DimsHorzVert,                         &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                 &
         SHORT_NAME='VFALLSN_AN',                               & 
         LONG_NAME ='reevaporation_fall_velocity_of_anvil_snow',              &
         UNITS     ='m s-1',                                    &
         DIMS      = MAPL_DimsHorzVert,                         &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                 &
         SHORT_NAME='VFALLSN_LS',                               & 
         LONG_NAME ='reevaporation_fall_velocity_of_largescale_snow',         &
         UNITS     ='m s-1',                                    &
         DIMS      = MAPL_DimsHorzVert,                         &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                 &
         SHORT_NAME='VFALLSN_CN',                               & 
         LONG_NAME ='reevaporation_fall_velocity_of_convective_snow',         &
         UNITS     ='m s-1',                                    &
         DIMS      = MAPL_DimsHorzVert,                         &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                 &
         SHORT_NAME='VFALLSN_SC',                               & 
         LONG_NAME ='reevaporation_fall_velocity_of_shallow_snow',         &
         UNITS     ='m s-1',                                    &
         DIMS      = MAPL_DimsHorzVert,                         &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='SUBLC',                                  & 
         LONG_NAME ='sublimation_of_cloud_ice',               &
         UNITS     ='kg kg-1 s-1',                            &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='FRZ_TT',                                 & 
         LONG_NAME ='freezing_of_cloud_condensate',           &
         UNITS     ='kg kg-1 s-1',                            &
         DIMS      = MAPL_DimsHorzVert,                       &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='FRZ_PP',                                 & 
         LONG_NAME ='freezing_of_precip_condensate',          &
         UNITS     ='kg kg-1 s-1',                            &
         DIMS      = MAPL_DimsHorzVert,                       &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='PFRZ',                                 & 
         LONG_NAME ='Probability_of_freezing_of_aerosol_part',          &
         UNITS     ='1',                            &
         DIMS      = MAPL_DimsHorzVert,                       &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

!!$    call MAPL_AddExportSpec(GC,                               &
!!$         SHORT_NAME='LIQANMOVE',                              &
!!$         LONG_NAME ='move2anv_source_of_anvil_liq',           &
!!$         UNITS     ='kg kg-1 s-1',                            &
!!$         DIMS      = MAPL_DimsHorzVert,                       &
!!$         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
!!$    VERIFY_(STATUS)
!!$
!!$    call MAPL_AddExportSpec(GC,                               &
!!$         SHORT_NAME='ICEANMOVE',                              & 
!!$         LONG_NAME ='move2anv_source_of_anvil_ice',           &
!!$         UNITS     ='kg kg-1 s-1',                            &
!!$         DIMS      = MAPL_DimsHorzVert,                       &
!!$         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
!!$    VERIFY_(STATUS)
!!$
!!$    call MAPL_AddExportSpec(GC,                                   &
!!$         SHORT_NAME = 'DLSCLD'  ,                                 & 
!!$         LONG_NAME = 'move2anv_change_in_large_scale_cloud_fraction',   &
!!$         UNITS     = '1',                                         &
!!$         DIMS      = MAPL_DimsHorzVert,                           &
!!$         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
!!$    VERIFY_(STATUS)
!!$
!!$    call MAPL_AddExportSpec(GC,                                   &
!!$         SHORT_NAME = 'DANCLD'  ,                                 & 
!!$         LONG_NAME = 'move2anv_change_in_anvil_cloud_fraction',   &
!!$         UNITS     = '1',                                         &
!!$         DIMS      = MAPL_DimsHorzVert,                           &
!!$         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
!!$    VERIFY_(STATUS)
!!$
!!$    call MAPL_AddExportSpec(GC,                               &
!!$         SHORT_NAME='CURAINMOVE',                             &
!!$         LONG_NAME ='movels2conv_source_of_cnv_rain',         &
!!$         UNITS     ='kg kg-1 s-1',                            &
!!$         DIMS      = MAPL_DimsHorzVert,                       &
!!$         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
!!$    VERIFY_(STATUS)
!!$
!!$    call MAPL_AddExportSpec(GC,                               &
!!$         SHORT_NAME='CUSNOWMOVE',                             & 
!!$         LONG_NAME ='movels2conv_source_of_cnv_snow',         &
!!$         UNITS     ='kg kg-1 s-1',                            &
!!$         DIMS      = MAPL_DimsHorzVert,                       &
!!$         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
!!$    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='PFLCNMOVE',                              &
         LONG_NAME ='moved_source_of_cnv_rain',               &
         UNITS     ='kg kg-1 s-1',                            &
         DIMS      = MAPL_DimsHorzVert,                       &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='PFICNMOVE',                              & 
         LONG_NAME ='moved_source_of_cnv_snow',               &
         UNITS     ='kg kg-1 s-1',                            &
         DIMS      = MAPL_DimsHorzVert,                       &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='CU2DRAINMOVE',                           &
         LONG_NAME ='moved_2d_source_of_cnv_rain',            &
         UNITS     ='kg kg-1 s-1',                            &
         DIMS      = MAPL_DimsHorzOnly,                       &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='CU2DSNOWMOVE',                           & 
         LONG_NAME ='moved_2d_source_of_cnv_snow',            &
         UNITS     ='kg kg-1 s-1',                            &
         DIMS      = MAPL_DimsHorzOnly,                       &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    ! Vertically integrated water substance conversions


    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='PDFLZ',                                        & 
         LONG_NAME ='statistical_source_of_cloud_water',          &
         UNITS     ='kg m-2 s-1'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='PDFIZ',                                        & 
         LONG_NAME ='statistical_source_of_cloud_ice',          &
         UNITS     ='kg m-2 s-1'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='CNVRNZ',                                        & 
         LONG_NAME ='convective_production_of_rain_water',          &
         UNITS     ='kg m-2 s-1'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='CNVLZ',                                        & 
         LONG_NAME ='convective_source_of_cloud_water',          &
         UNITS     ='kg m-2 s-1'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='CNVIZ',                                        & 
         LONG_NAME ='convective_source_of_cloud_ice',          &
         UNITS     ='kg m-2 s-1'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='EVPCZ',                                        & 
         LONG_NAME ='evaporation_loss_of_cloud_water',          &
         UNITS     ='kg m-2 s-1'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='SUBCZ',                                        & 
         LONG_NAME ='sublimation_loss_of_cloud_ice',          &
         UNITS     ='kg m-2 s-1'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='EVPPZ',                                        & 
         LONG_NAME ='evaporation_loss_of_precip_water',          &
         UNITS     ='kg m-2 s-1'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='SUBPZ',                                        & 
         LONG_NAME ='sublimation_loss_of_precip_ice',          &
         UNITS     ='kg m-2 s-1'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='AUTZ',                                        & 
         LONG_NAME ='autoconversion_loss_of_cloud_water',          &
         UNITS     ='kg m-2 s-1'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='SDMZ',                                        & 
         LONG_NAME ='sedimentation_loss_of_cloud_ice',          &
         UNITS     ='kg m-2 s-1'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='COLLLZ',                                        & 
         LONG_NAME ='accretion_loss_of_cloud_water_to_rain',          &
         UNITS     ='kg m-2 s-1'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='COLLIZ',                                        & 
         LONG_NAME ='accretion_loss_of_cloud_water_to_snow',          &
         UNITS     ='kg m-2 s-1'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='FRZCZ',                                        & 
         LONG_NAME ='net_freezing_of_cloud_condensate',          &
         UNITS     ='kg m-2 s-1'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='FRZPZ',                                        & 
         LONG_NAME ='net_freezing_of_precip_condensate',          &
         UNITS     ='kg m-2 s-1'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)


    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='RCCODE',                                      & 
         LONG_NAME ='Convection_return_codes',                     &
         UNITS     ='codes',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='TRIEDLV',                                     & 
         LONG_NAME ='Tested_for_convection_at_this_level',         &
         UNITS     ='0 or 1',                                      &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    ! MATMAT Exports for after-RAS inoutputs

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='QVRAS',                                       & 
         LONG_NAME ='water_vapor_after_ras',                       &
         UNITS     ='kg kg-1',                                     &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='THRAS',                                       & 
         LONG_NAME ='potential_temperature_after_ras',             &
         UNITS     ='K',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='URAS',                                        & 
         LONG_NAME ='eastward_wind_after_ras',                     &
         UNITS     ='m s-1',                                       &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='VRAS',                                        & 
         LONG_NAME ='northward_wind_after_ras',                    &
         UNITS     ='m s-1',                                       &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    ! MATMAT Exports for before-RAS inputs for RAStest

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='THOI',                                        & 
         LONG_NAME ='potential_temperature_before_ras',            &
         UNITS     ='K',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='QHOI',                                        & 
         LONG_NAME ='specific_humidity_before_ras',                &
         UNITS     ='kg kg-1',                                     &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='QSSI',                                        & 
         LONG_NAME ='saturation_specific_humidity_before_ras',     &
         UNITS     ='kg kg-1',                                     &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='DQSI',                                        & 
         LONG_NAME ='deriv_sat_specific_humidity_wrt_t_before_ras',&
         UNITS     ='kg kg-1 K-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='PLEI',                                        & 
         LONG_NAME ='air_pressure_before_ras',                     &
         UNITS     ='Pa',                                          &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationEdge,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='TPERTI',                                      & 
         LONG_NAME ='temperature_perturbation_before_ras',         &
         UNITS     ='K',                                           &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='KCBLI',                                       & 
         LONG_NAME ='cloud_base_layer_before_ras',                 &
         UNITS     ='1',                                           &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)


    ! !RECORD IMPORTS AND INTERNALS AT TOP OF MOIST:

    call MAPL_AddExportSpec(GC,                           &
         SHORT_NAME ='QX0',                                          &
         LONG_NAME  ='specific_humidity',                          &
         UNITS      ='kg kg-1',                                    &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)                                                                          

    call MAPL_AddExportSpec(GC,                           &
         SHORT_NAME ='QX1',                                        &
         LONG_NAME  ='specific_humidity_after_moist_physics',      &
         UNITS      ='kg kg-1',                                    &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)                                                                          

    call MAPL_AddExportSpec(GC,                           &
         SHORT_NAME ='QLLSX0',                                       &
         LONG_NAME  ='initial_mass_fraction_of_large_scale_cloud_liquid_water', &
         UNITS      ='kg kg-1',                                    &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)                                                                          

    call MAPL_AddExportSpec(GC,                           &
         SHORT_NAME ='QLLSX1',                                       &
         LONG_NAME  ='final_mass_fraction_of_large_scale_cloud_liquid_water', &
         UNITS      ='kg kg-1',                                    &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)                                                                          

    call MAPL_AddExportSpec(GC,                           &
         SHORT_NAME ='QLCNX0',                                       &
         LONG_NAME  ='initial_mass_fraction_of_convective_cloud_liquid_water', &
         UNITS      ='kg kg-1',                                    &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)                                                                          

    call MAPL_AddExportSpec(GC,                           &
         SHORT_NAME ='QLCNX1',                                       &
         LONG_NAME  ='final_mass_fraction_of_convective_cloud_liquid_water', &
         UNITS      ='kg kg-1',                                    &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)                                                                          

    call MAPL_AddExportSpec(GC,                           &
         SHORT_NAME ='CLLSX0',                                       &
         LONG_NAME  ='large_scale_cloud_area_fraction',            &
         UNITS      ='1',                                          &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)                                                                          

    call MAPL_AddExportSpec(GC,                           &
         SHORT_NAME ='CLCNX0',                                       &
         LONG_NAME  ='convective_cloud_area_fraction',             &
         UNITS      ='1',                                          &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)                                                                          

    call MAPL_AddExportSpec(GC,                           &
         SHORT_NAME ='QILSX0',                                       &
         LONG_NAME  ='initial_mass_fraction_of_large_scale_cloud_ice_water', &
         UNITS      ='kg kg-1',                                    &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)                                                                          

    call MAPL_AddExportSpec(GC,                           &
         SHORT_NAME ='QILSX1',                                       &
         LONG_NAME  ='final_mass_fraction_of_large_scale_cloud_ice_water', &
         UNITS      ='kg kg-1',                                    &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)                                                                          


    call MAPL_AddExportSpec(GC,                           &
         SHORT_NAME ='QICNX0',                                       &
         LONG_NAME  ='initial_mass_fraction_of_convective_cloud_ice_water', &
         UNITS      ='kg kg-1',                                    &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)                                                                          

    call MAPL_AddExportSpec(GC,                           &
         SHORT_NAME ='QICNX1',                                       &
         LONG_NAME  ='final_mass_fraction_of_convective_cloud_ice_water', &
         UNITS      ='kg kg-1',                                    &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)                                                                          

    call MAPL_AddExportSpec(GC,                           &
         SHORT_NAME ='QCCNX0',                                       &
         LONG_NAME  ='initial_mass_fraction_of_convective_cloud_condensate', &
         UNITS      ='kg kg-1',                                    &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)                                                                          

    call MAPL_AddExportSpec(GC,                           &
         SHORT_NAME ='QCLSX0',                                       &
         LONG_NAME  ='initial_mass_fraction_of_large_scale_cloud_condensate', &
         UNITS      ='kg kg-1',                                    &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)                                                                          


    call MAPL_AddExportSpec(GC,                              &
         SHORT_NAME = 'W_DIAG',                                        &
         LONG_NAME  = 'Diagnostic_vertical_velocity',                  &
         UNITS      = 'm s-1',                                        &
         DIMS       = MAPL_DimsHorzVert,                              &
         VLOCATION  = MAPL_VLocationCenter,                           &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                              &
         SHORT_NAME = 'KHX0',                                         &
         LONG_NAME  = 'scalar_diffusivity',                         &
         UNITS      = 'm+2 s-1',                                    &
         DIMS       = MAPL_DimsHorzVert,                            &
         VLOCATION  = MAPL_VLocationEdge,                           &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                             &
         SHORT_NAME = 'THX0',                                        &
         LONG_NAME  = 'potential_temperature',                     &
         UNITS      = 'K',                                         &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,                        &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                             &
         SHORT_NAME = 'UX0',                                         &
         LONG_NAME  = 'eastward_wind',                             &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,                        &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                             &
         SHORT_NAME = 'VX0',                                         &
         LONG_NAME  = 'northward_wind',                            &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,                        &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                             &
         SHORT_NAME = 'TSX0',                                        &
         LONG_NAME  = 'surface temperature',                       &
         UNITS      = 'K',                                         &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,                          &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                             &
         SHORT_NAME = 'FRLANDX0',                                    &
         LONG_NAME  = 'areal_land_fraction',                       &
         UNITS      = '1',                                     &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,                        &
         RC=STATUS  )
    VERIFY_(STATUS)


!!! downdraft diags
    call MAPL_AddExportSpec(GC,                             &
         SHORT_NAME = 'DDF_MFC',                                   &
         LONG_NAME  = 'Downdraft_mass_flux',                       &
         UNITS      = 'kg m-2 s-1',                                &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationEdge,                          &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                             &
         SHORT_NAME = 'DDF_RH1',                                  &
         LONG_NAME  = 'Downdraft_in_cloud_RH_before',       &
         UNITS      = '1',                                       &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,                        &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                             &
         SHORT_NAME = 'DDF_RH2',                                  &
         LONG_NAME  = 'Downdraft_in_cloud_RH_after',       &
         UNITS      = '1',                                       &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,                        &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                             &
         SHORT_NAME = 'DDF_TC',                                    &
         LONG_NAME  = 'Temperature_excess_in_DDF',                 &
         UNITS      = 'K',                                         &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,                        &
         RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC,                             &
         SHORT_NAME = 'DDF_QVC',                                   &
         LONG_NAME  = 'Spec_hum_excess_in_DDF',                    &
         UNITS      = 'kg kg-1',                                   &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,                        &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                             &
         SHORT_NAME = 'DDF_BYNC',                                  &
         LONG_NAME  = 'Buoyancy_of_DDF',                           &
         UNITS      = 'm s-2',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,                        &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                             &
         SHORT_NAME = 'DDF_MUPH',                                  &
         LONG_NAME  = 'Downdraft_moistening_from_evap_subl',       &
         UNITS      = 'kg kg-1 s-1',                               &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,                        &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                             &
         SHORT_NAME = 'DDF_DQDT',                                  &
         LONG_NAME  = 'Total_Downdraft_moistening',                      &
         UNITS      = 'kg kg-1 s-1',                               &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,                        &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                             &
         SHORT_NAME = 'DDF_DTDT',                                  &
         LONG_NAME  = 'Total_Downdraft_heating',                         &
         UNITS      = 'K s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,                        &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                             &
         SHORT_NAME = 'DDF_ZSCALE',                                &
         LONG_NAME  = 'vertical_scale_for_downdraft',              &
         UNITS      = 'm',                                         &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,                        &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                            &
         SHORT_NAME = 'KEDISS',                                    &
         LONG_NAME  = 'kinetic_energy_diss_in_RAS',       &
         UNITS      = 'W m-2',                                    &
         DIMS       = MAPL_DimsHorzVert,                          &
         VLOCATION  = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                            &
         SHORT_NAME = 'KEMST',                                    &
         LONG_NAME  = 'vertically_integrated_kinetic_energy_tendency_across_moist',       &
         UNITS      = 'W m-2',                                    &
         DIMS       = MAPL_DimsHorzOnly,                          &
         VLOCATION  = MAPL_VLocationNone,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                            &
         SHORT_NAME = 'KEMST2',                                    &
         LONG_NAME  = 'vertically_integrated_KE_dissipation_in_RAS',       &
         UNITS      = 'W m-2',                                    &
         DIMS       = MAPL_DimsHorzOnly,                          &
         VLOCATION  = MAPL_VLocationNone,              RC=STATUS  )
    VERIFY_(STATUS)

    ! CAR 12/5/08
    ! Aerosol Scavenging diagnostics/export states
    ! ------------------------------
    call MAPL_AddExportSpec(GC,                                       &
         SHORT_NAME='DDU2gDT ',                                         &
         LONG_NAME ='dust_tendency_due_to_conv_scav',                 &
         UNITS     ='kg m-2 s-1',                                     &
         DIMS      = MAPL_DimsHorzOnly,                               &
         __RC__  )

    call MAPL_AddExportSpec(GC,                                       &
         SHORT_NAME='DSS2gDT ',                                         &
         LONG_NAME ='sea_salt_tendency_due_to_conv_scav',             &
         UNITS     ='kg m-2 s-1',                                     &
         DIMS      = MAPL_DimsHorzOnly,                               &
         __RC__  )

    call MAPL_AddExportSpec(GC,                                       &
         SHORT_NAME='DOC2gDT ',                                         &
         LONG_NAME ='organic_carbon_tendency_due_to_conv_scav',       &
         UNITS     ='kg m-2 s-1',                                     &
         DIMS      = MAPL_DimsHorzOnly,                               &
         __RC__  )

    call MAPL_AddExportSpec(GC,                                       &
         SHORT_NAME='DBC2gDT ',                                         &
         LONG_NAME ='black_carbon_tendency_due_to_conv_scav',         &
         UNITS     ='kg m-2 s-1',                                     &
         DIMS      = MAPL_DimsHorzOnly,                               &
         __RC__  )

    call MAPL_AddExportSpec(GC,                                       &
         SHORT_NAME='DSU2gDT ',                                         &
         LONG_NAME ='sulfate_tendency_due_to_conv_scav',              &
         UNITS     ='kg m-2 s-1',                                     &
         DIMS      = MAPL_DimsHorzOnly,                               &
         __RC__  )

    call MAPL_AddExportSpec(GC,                                       &
         SHORT_NAME='DNI2gDT',                                          &
         LONG_NAME ='nitrate_tendency_due_to_conv_scav',              &
         UNITS     ='kg m-2 s-1',                                     &
         DIMS      = MAPL_DimsHorzOnly,                               &
         __RC__  )

    call MAPL_AddExportSpec(GC,                                       &
         SHORT_NAME='DNH4A2gDT',                                        &
         LONG_NAME ='ammonium_aerosol_tendency_due_to_conv_scav',     &
         UNITS     ='kg m-2 s-1',                                     &
         DIMS      = MAPL_DimsHorzOnly,                               &
         __RC__  )

    call MAPL_AddExportSpec(GC,                                       &
         SHORT_NAME='DNH32gDT',                                         &
         LONG_NAME ='ammonia_tendency_due_to_conv_scav',              &
         UNITS     ='kg m-2 s-1',                                     &
         DIMS      = MAPL_DimsHorzOnly,                               &
         __RC__  )

    call MAPL_AddExportSpec(GC,                                       &
         SHORT_NAME='DBRC2gDT',                                          &
         LONG_NAME ='brown_carbon_tendency_due_to_conv_scav',              &
         UNITS     ='kg m-2 s-1',                                     &
         DIMS      = MAPL_DimsHorzOnly,                               &
         __RC__  )





    call MAPL_AddExportSpec(GC,                                       &
         SHORT_NAME='DDUDT ',                                         & 
         LONG_NAME ='dust_tendency_due_to_conv_scav',                 &
         UNITS     ='kg m-2 s-1',                                     &
         DIMS      = MAPL_DimsHorzOnly,                               &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                       &
         SHORT_NAME='DSSDT ',                                         &
         LONG_NAME ='sea_salt_tendency_due_to_conv_scav',             &
         UNITS     ='kg m-2 s-1',                                     &
         DIMS      = MAPL_DimsHorzOnly,                               &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                       &
         SHORT_NAME='DOCDT ',                                         &
         LONG_NAME ='organic_carbon_tendency_due_to_conv_scav',       &
         UNITS     ='kg m-2 s-1',                                     &
         DIMS      = MAPL_DimsHorzOnly,                               &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                       &
         SHORT_NAME='DBCDT ',                                         &
         LONG_NAME ='black_carbon_tendency_due_to_conv_scav',         &
         UNITS     ='kg m-2 s-1',                                     &
         DIMS      = MAPL_DimsHorzOnly,                               &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                       &
         SHORT_NAME='DSUDT ',                                         &
         LONG_NAME ='sulfate_tendency_due_to_conv_scav',              &
         UNITS     ='kg m-2 s-1',                                     &
         DIMS      = MAPL_DimsHorzOnly,                               &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                       &
         SHORT_NAME='DNIDT',                                          &
         LONG_NAME ='nitrate_tendency_due_to_conv_scav',              &
         UNITS     ='kg m-2 s-1',                                     &
         DIMS      = MAPL_DimsHorzOnly,                               &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                       &
         SHORT_NAME='DNH4ADT',                                        &
         LONG_NAME ='ammonium_aerosol_tendency_due_to_conv_scav',     &
         UNITS     ='kg m-2 s-1',                                     &
         DIMS      = MAPL_DimsHorzOnly,                               &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                       &
         SHORT_NAME='DNH3DT',                                         &
         LONG_NAME ='ammonia_tendency_due_to_conv_scav',              &
         UNITS     ='kg m-2 s-1',                                     &
         DIMS      = MAPL_DimsHorzOnly,                               &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                       &
         SHORT_NAME='DBRCDT',                                          &
         LONG_NAME ='brown_carbon_tendency_due_to_conv_scav',              &
         UNITS     ='kg m-2 s-1',                                     &
         DIMS      = MAPL_DimsHorzOnly,                               &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                       &
         SHORT_NAME='DDUDTcarma ',                                    &
         LONG_NAME ='carma_dust_tendency_due_to_conv_scav',           &
         UNITS     ='kg m-2 s-1',                                     &
         DIMS      = MAPL_DimsHorzOnly,                               &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                       &
         SHORT_NAME='DSSDTcarma ',                                    &
         LONG_NAME ='carma_seasalt_tendency_due_to_conv_scav',        &
         UNITS     ='kg m-2 s-1',                                     &
         DIMS      = MAPL_DimsHorzOnly,                               &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                       &
         SHORT_NAME='LFR_GCC',                                        &
         LONG_NAME ='lightning_flash_rate_for_GEOSCHEMchem',          &
         UNITS     ='km-2 s-1',                                       &
         DIMS      = MAPL_DimsHorzOnly,                               &
         VLOCATION = MAPL_VLocationNone,                              &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='THMOIST',                                     & 
         LONG_NAME ='potential_temperature_after_all_of_moist',   &
         UNITS     ='K',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='SMOIST',                                      & 
         LONG_NAME ='dry_static_energy_after_all_of_moist',        &
         UNITS     ='m+2 s-2',                                     &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)


    !-------Aerosol Cloud Interactions Diagnostics  

    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME='SMAX_LIQ',                                          & 
         LONG_NAME ='Maximum incloud supersaturation for liquid',        &
         UNITS     ='%',                                                 &
         DIMS      = MAPL_DimsHorzVert,                                  &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME='WSUB',                                          & 
         LONG_NAME ='Subgrid Scale in-cloud vertical velocity',          &
         UNITS     ='m s-1',                                             &
         DIMS      = MAPL_DimsHorzVert,                                  &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                        &
         SHORT_NAME='CCN01',                                             & 
         LONG_NAME ='CCN conc at 0.1 % supersaturation (grid_avg)',                 &
         UNITS     ='m-3',                                             &
         DIMS      = MAPL_DimsHorzVert,                                  &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME='CCN04',                                             & 
         LONG_NAME ='CCN conc at 0.4 % supersaturation (grid_avg)',                 &
         UNITS     ='m-3',                                             &
         DIMS      = MAPL_DimsHorzVert,                                  &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME='CCN1',                                              & 
         LONG_NAME ='CCN conc at 1.0 % supersaturation (grid_avg)',                 &
         UNITS     ='m-3',                                             &
         DIMS      = MAPL_DimsHorzVert,                                  &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME='SMAX_ICE',                                          & 
         LONG_NAME ='Maximum incloud supersaturation for ice',        &
         UNITS     ='%',                                                 &
         DIMS      = MAPL_DimsHorzVert,                                  &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)  

    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME='CDNC_NUC',                                          & 
         LONG_NAME ='Nucleated cloud droplet concentration (grid_avg)',        &
         UNITS     ='m-3',                                                 &
         DIMS      = MAPL_DimsHorzVert,                                  &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS) 

    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME='INC_NUC',                                          & 
         LONG_NAME ='Nucleated ice crystal concentration (grid_avg)',        &
         UNITS     ='m-3',                                                 &
         DIMS      = MAPL_DimsHorzVert,                                  &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS) 

    call MAPL_AddExportSpec(GC,                           &
         SHORT_NAME ='NCPL_VOL',                                       &
         LONG_NAME  ='particle_number_for_liquid_cloud', &
         UNITS      ='m-3',                                          &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)                                                                          

    call MAPL_AddExportSpec(GC,                           &
         SHORT_NAME ='NCPI_VOL',                                       &
         LONG_NAME  ='particle_number_for_ice_cloud', &
         UNITS      ='m-3',                                          &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS) 

    call MAPL_AddExportSpec(GC,                           &
         SHORT_NAME ='SO4',                                       &
         LONG_NAME  ='Sulfate number conc.', &
         UNITS      ='m-3',                                          &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS) 

    call MAPL_AddExportSpec(GC,                           &
         SHORT_NAME ='ORG',                                       &
         LONG_NAME  ='Organic number conc. (hydrophilic)',        &
         UNITS      ='m-3',                                          &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS) 

    call MAPL_AddExportSpec(GC,                           &
         SHORT_NAME ='BCARBON',                                       &
         LONG_NAME  ='Black carbon number conc. (hydrophilic)',        &
         UNITS      ='m-3',                                          &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS) 

    call MAPL_AddExportSpec(GC,                           &
         SHORT_NAME ='DUST',                                       &
         LONG_NAME  ='Total dust number conc', &
         UNITS      ='m-3',                                          &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS) 

    call MAPL_AddExportSpec(GC,                           &
         SHORT_NAME ='SEASALT',                                       &
         LONG_NAME  ='Total sea number conc', &
         UNITS      ='m-3',                                          &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS) 

    call MAPL_AddExportSpec(GC,                           &
         SHORT_NAME ='NHET_NUC',                                       &
         LONG_NAME  ='Nucleated ice crystal concentration by het freezing (grid_avg)', &
         UNITS      ='m-3',                                          &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS) 

    call MAPL_AddExportSpec(GC,                           &
         SHORT_NAME ='NLIM_NUC',                                       &
         LONG_NAME  ='Limiting IN concentration allowing hom freezing', &
         UNITS      ='m-3',                                          &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='SAT_RAT',                                         & 
         LONG_NAME ='saturation_ratio_after_moist',               &
         UNITS     ='1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS) 



    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='DQVDT_micro',                                         & 
         LONG_NAME ='Q tendency due to microphysics ',               &
         UNITS     ='kg kg-1 s-1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)  

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='DQLDT_micro',                                         & 
         LONG_NAME ='QL tendency due to microphysics ',               &
         UNITS     ='kg kg-1 s-1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)  

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='DQIDT_micro',                                         & 
         LONG_NAME ='QI tendency due to microphysics ',               &
         UNITS     ='kg kg-1 s-1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)  

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='DQRDT_micro',                                         &
         LONG_NAME ='QR tendency due to microphysics ',               &
         UNITS     ='kg kg-1 s-1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='DQSDT_micro',                                         &
         LONG_NAME ='QS tendency due to microphysics ',               &
         UNITS     ='kg kg-1 s-1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='DQGDT_micro',                                         &
         LONG_NAME ='QG tendency due to microphysics ',               &
         UNITS     ='kg kg-1 s-1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='DQADT_micro',                                         &
         LONG_NAME ='QA tendency due to microphysics ',               &
         UNITS     ='kg kg-1 s-1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)


    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='DTDT_micro',                                         & 
         LONG_NAME ='T tendency due to microphysics ',               &
         UNITS     ='K s-1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)  

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='DUDT_micro',                                         &
         LONG_NAME ='U tendency due to microphysics ',               &
         UNITS     ='m s-2',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='DVDT_micro',                                         &
         LONG_NAME ='V tendency due to microphysics ',               &
         UNITS     ='m s-2',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)


    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='TVQX1',                                        & 
         LONG_NAME ='Total_Water_Substance_bef_macro',                 &
         UNITS     ='kg m-2'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='TVQX2',                                        & 
         LONG_NAME ='Total_Water_Substance_bef_micro',                 &
         UNITS     ='kg m-2'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)


    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='RL_MASK',                                          & 
         LONG_NAME ='volumetric_liquid_cloud_particle_volume_radius',      &
         UNITS     ='m',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='RI_MASK',                                          & 
         LONG_NAME ='volumetric_ice_cloud_particle_volume_radius',      &
         UNITS     ='m',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                   &
         SHORT_NAME='KAPPA',                                       & 
         LONG_NAME ='kappa parameter for activation',              &
         UNITS     ='1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                             &
         SHORT_NAME ='CFLIQ',                                       &
         LONG_NAME  ='liquid_cloud_area_fraction (LS)',            &
         UNITS      ='1',                                          &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)      

    call MAPL_AddExportSpec(GC,                             &
         SHORT_NAME ='CFICE',                                       &
         LONG_NAME  ='ice_cloud_area_fraction (LS)',            &
         UNITS      ='1',                                          &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS) 

    call MAPL_AddExportSpec(GC,                             &
         SHORT_NAME ='RHICE',                                       &
         LONG_NAME  ='Relative humidity wrt ice',            &
         UNITS      ='1',                                          &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                             &
         SHORT_NAME ='RHLIQ',                                       &
         LONG_NAME  ='Relative humidity wrt liquid',            &
         UNITS      ='1',                                          &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                             &
         SHORT_NAME ='SC_ICE',                                       &
         LONG_NAME  ='Effective Freezing RHi',            &
         UNITS      ='1',                                          &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)


    call MAPL_AddExportSpec(GC,                           &
         SHORT_NAME ='NHET_IMM',                                       &
         LONG_NAME  ='Immersion_IN', &
         UNITS      ='m-3',                                          &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                           &
         SHORT_NAME ='NHET_DEP',                                       &
         LONG_NAME  ='Deposition_IN', &
         UNITS      ='m-3',                                          &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)


    call MAPL_AddExportSpec(GC,                           &
         SHORT_NAME ='DUST_IMM',                                       &
         LONG_NAME  ='Immersion IN from dust', &
         UNITS      ='m-3',                                          &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                           &
         SHORT_NAME ='DUST_DEP',                                       &
         LONG_NAME  ='deposition IN from dust', &
         UNITS      ='m-3',                                          &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                           &
         SHORT_NAME ='SCF',                                       &
         LONG_NAME  ='Supercooled cloud fraction', &
         UNITS      ='1',                                          &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                           &
         SHORT_NAME ='SCF_ALL',                                       &
         LONG_NAME  ='Supercooled cloud fraction including snow and rain', &
         UNITS      ='1',                                          &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME='SIGW_GW',                                          & 
         LONG_NAME ='Subgrid Scale vertical velocity variance from GW',  &
         UNITS     ='m2 s-2',                                             &
         DIMS      = MAPL_DimsHorzVert,                                  &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME='SIGW_CNV',                                          & 
         LONG_NAME ='Subgrid Scale vertical velocity variance from convection',  &
         UNITS     ='m2 s-2',                                             &
         DIMS      = MAPL_DimsHorzVert,                                  &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)


    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME='SIGW_TURB',                                          & 
         LONG_NAME ='Subgrid Scale vertical velocity variance from turbulence', &
         UNITS     ='m2 s-2',                                             &
         DIMS      = MAPL_DimsHorzVert,                                  &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME='SIGW_RC',                                          & 
         LONG_NAME ='Mean subgrid Scale vertical velocity from rad cooling', &
         UNITS     ='m s-1',                                             &
         DIMS      = MAPL_DimsHorzVert,                                  &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME='CNV_FICE',                                          & 
         LONG_NAME ='Ice fraction in convective tower', &
         UNITS     ='1',                                             &
         DIMS      = MAPL_DimsHorzVert,                                  &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME='CNV_NDROP',                                          & 
         LONG_NAME ='Droplet number conc. in conv. detrainment', &
         UNITS     ='m-3',                                             &
         DIMS      = MAPL_DimsHorzVert,                                  &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME='CNV_NICE',                                          & 
         LONG_NAME ='Ice crystal number conc. in conv. detrainment', &
         UNITS     ='m-3',                                             &
         DIMS      = MAPL_DimsHorzVert,                                  &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME='SC_NDROP',                                          & 
         LONG_NAME ='Droplet number conc. in shallow detrainment', &
         UNITS     ='m-3',                                             &
         DIMS      = MAPL_DimsHorzVert,                                  &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME='SC_NICE',                                          & 
         LONG_NAME ='Ice crystal number conc. in shallow detrainment', &
         UNITS     ='m-3',                                             &
         DIMS      = MAPL_DimsHorzVert,                                  &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME='RHCmicro',                                          & 
         LONG_NAME ='Corrected RHc after micro', &
         UNITS     ='1',                                             &
         DIMS      = MAPL_DimsHorzVert,                                  &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                   &
         SHORT_NAME='CCNCOLUMN',                                   & 
         LONG_NAME ='Vertically integrated CCN at 1% ssat',        &
         UNITS     ='m-2'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                   &
         SHORT_NAME='NDCOLUMN',                                   & 
         LONG_NAME ='Vertically integrated NCPL',        &
         UNITS     ='m-2'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                   &
         SHORT_NAME='NCCOLUMN',                                   & 
         LONG_NAME ='Vertically integrated NCPI',        &
         UNITS     ='m-2'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME ='BERG',                                                 &
         LONG_NAME  ='ice mixing ration tendency due to Bergeron process ',  &
         UNITS      ='Kg Kg-1 s-1',                                            &
         DIMS       = MAPL_DimsHorzVert,                                     &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME ='BERGS',                                                 &
         LONG_NAME  ='Snow mixing ration tendency due to Bergeron process ',  &
         UNITS      ='Kg Kg-1 s-1',                                            &
         DIMS       = MAPL_DimsHorzVert,                                     &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME ='MELT',                                                 &
         LONG_NAME  ='Melting of cloud ice',  &
         UNITS      ='Kg Kg-1 s-1',                                            &
         DIMS       = MAPL_DimsHorzVert,                                     &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)


    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME ='QCRES',                                                 &
         LONG_NAME  ='Residual cloud tendency in micro',  &
         UNITS      ='Kg Kg-1 s-1',                                            &
         DIMS       = MAPL_DimsHorzVert,                                     &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME ='QIRES',                                                 &
         LONG_NAME  ='Residual ice cloud tendency in micro',  &
         UNITS      ='Kg Kg-1 s-1',                                            &
         DIMS       = MAPL_DimsHorzVert,                                     &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='AUTICE',                                      & 
         LONG_NAME ='autoconv_sink_of_cloud_ice',               &
         UNITS     ='kg kg-1 s-1',                                   &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)


    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME ='DT_RASP',                                                 &
         LONG_NAME  ='T tendency from ras precip',  &
         UNITS      ='K s-1',                                            &
         DIMS       = MAPL_DimsHorzVert,                                     &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)


    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME ='DNHET_IMM',                                                 &
         LONG_NAME  ='Ice number tendency due to immersion freezing',  &
         UNITS      ='Kg-1 s-1',                                            &
         DIMS       = MAPL_DimsHorzVert,                                     &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)


    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME ='DNHET_CT',                                                 &
         LONG_NAME  ='Ice number tendency due to contact freezing',  &
         UNITS      ='Kg-1 s-1',                                            &
         DIMS       = MAPL_DimsHorzVert,                                     &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='DTDT_macro',                                         & 
         LONG_NAME ='T tendency due to macrophysics ',               &
         UNITS     ='K s-1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)  

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='DQVDT_macro',                                         &
         LONG_NAME ='QV tendency due to macrophysics ',               &
         UNITS     ='kg kg-1 s-1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='DQLDT_macro',                                         &
         LONG_NAME ='QL tendency due to macrophysics ',               &
         UNITS     ='kg kg-1 s-1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='DQIDT_macro',                                         &
         LONG_NAME ='QI tendency due to macrophysics ',               &
         UNITS     ='kg kg-1 s-1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='DQADT_macro',                                         &
         LONG_NAME ='QA tendency due to macrophysics ',               &
         UNITS     ='kg kg-1 s-1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='DTDT_moist',                                         & 
         LONG_NAME ='T tendency due to moist',               &
         UNITS     ='K s-1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)  

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='DTDTCN',                                         & 
         LONG_NAME ='T tendency due to convection',               &
         UNITS     ='K s-1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)  

    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME ='FRZPP_LS',                                                 &
         LONG_NAME  ='Tendency in T from micro precip freezing',  &
         UNITS      ='K s-1',                                            &
         DIMS       = MAPL_DimsHorzVert,                                     &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME ='SNOWMELT_LS',                                                 &
         LONG_NAME  ='Tendency in T from micro snow melting',  &
         UNITS      ='K s-1',                                            &
         DIMS       = MAPL_DimsHorzVert,                                     &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)


    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME ='DNCNUC',                                                 &
         LONG_NAME  ='Ice number tendency due to nucleation on aerosol',  &
         UNITS      ='m-3 s-1',                                            &
         DIMS       = MAPL_DimsHorzVert,                                     &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)


    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME ='DNCHMSPLIT',                                                 &
         LONG_NAME  ='Ice number tendency due to H-M splittering',  &
         UNITS      ='m-3 s-1',                                            &
         DIMS       = MAPL_DimsHorzVert,                                     &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME ='DNCSUBL',                                                 &
         LONG_NAME  ='Ice number tendency due to sublimation',  &
         UNITS      ='m-3 s-1',                                            &
         DIMS       = MAPL_DimsHorzVert,                                     &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME ='DNCAUTICE',                                                 &
         LONG_NAME  ='Ice number tendency due to autoconversion to snow',  &
         UNITS      ='m-3 s-1',                                            &
         DIMS       = MAPL_DimsHorzVert,                                     &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME ='DNCACRIS',                                                 &
         LONG_NAME  ='Ice number tendency due to accretion by snow',  &
         UNITS      ='m-3 s-1',                                            &
         DIMS       = MAPL_DimsHorzVert,                                     &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)


    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME ='DNCCNV',                                                 &
         LONG_NAME  ='Ice crystal number tendency from convective detrainment',  &
         UNITS      ='m-3 s-1',                                            &
         DIMS       = MAPL_DimsHorzVert,                                     &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)


    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME ='DNDCCN',                                                 &
         LONG_NAME  ='Cloud droplet number tendency due to activation',  &
         UNITS      ='m-3 s-1',                                            &
         DIMS       = MAPL_DimsHorzVert,                                     &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME ='DNDACRLS',                                                 &
         LONG_NAME  ='Cloud droplet number tendency due to accretion by snow',  &
         UNITS      ='m-3 s-1',                                            &
         DIMS       = MAPL_DimsHorzVert,                                     &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)


    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME ='DNDEVAPC',                                                 &
         LONG_NAME  ='Cloud droplet number tendency due to evaporation',  &
         UNITS      ='m-3 s-1',                                            &
         DIMS       = MAPL_DimsHorzVert,                                     &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME ='DNDACRLR',                                                 &
         LONG_NAME  ='Cloud droplet number tendency due to accretion by rain',  &
         UNITS      ='m-3 s-1',                                            &
         DIMS       = MAPL_DimsHorzVert,                                     &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME ='DNDAUTLIQ',                                                 &
         LONG_NAME  ='Cloud droplet number tendency due to autoconversion',  &
         UNITS      ='m-3 s-1',                                            &
         DIMS       = MAPL_DimsHorzVert,                                     &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)


    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME ='DNDCNV',                                                 &
         LONG_NAME  ='Cloud droplet number tendency from convective detrainment',  &
         UNITS      ='m-3 s-1',                                            &
         DIMS       = MAPL_DimsHorzVert,                                     &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)


    call MAPL_AddExportSpec(GC,                                   &
         SHORT_NAME='CLDREFFL_TOP',                                   & 
         LONG_NAME ='Droplet effective radius at cloud top',        &
         UNITS     ='m'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                   &
         SHORT_NAME='CLDREFFI_TOP',                                   & 
         LONG_NAME ='ice crystal effective radius at cloud top',        &
         UNITS     ='m'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                   &
         SHORT_NAME='NCPL_TOP',                                   & 
         LONG_NAME ='Grid-averaged NCPL at cloud top',        &
         UNITS     ='m-3'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                   &
         SHORT_NAME='NCPI_TOP',                                   & 
         LONG_NAME ='Grid-averaged NCPI at cloud top',        &
         UNITS     ='m-3'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                   &
         SHORT_NAME='NCPL_CLDBASE',                                   & 
         LONG_NAME ='IN-CLOUD NCPL at cloud base',        &
         UNITS     ='m-3'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)
    


    call MAPL_AddExportSpec(GC,                             &
         SHORT_NAME = 'LWC',                                      &
         LONG_NAME  = 'liquid water content',        &
         UNITS      = 'kg m-3',                                    &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                             &
         SHORT_NAME = 'IWC',                                      &
         LONG_NAME  = 'ice water content',        &
         UNITS      = 'kg m-3',                                    &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                             &
         SHORT_NAME='QCVAR_EXP',                                   & 
         LONG_NAME ='inverse relative variance of cloud water',        &
         UNITS      = '1',                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)
        
                
    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME='LTS',                                          & 
         LONG_NAME ='Lower tropospheric stability', &
         UNITS     ='K',                                             &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
     VERIFY_(STATUS)
        
      call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME='EIS',                                          & 
         LONG_NAME ='Estimated Inversion Strength', &
         UNITS     ='K',                                             &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
     VERIFY_(STATUS)
    
    
     call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME='RAS_ALPHA',                                          & 
         LONG_NAME ='RAS relaxation parameter', &
         UNITS     ='1',                                             &
         DIMS      = MAPL_DimsHorzVert,                                  &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                                        &
         SHORT_NAME='RAS_TAU',                                          & 
         LONG_NAME ='RAS total relaxation timescale',                   &
         UNITS     ='1',                                                &
         DIMS      = MAPL_DimsHorzVert,                                 &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)
    
    call MAPL_AddImportSpec(GC,                                        &
         SHORT_NAME = 'DTDT_BL',                                       &
         LONG_NAME  = 'tendency_of_air_temperature_due_to_bound_layer', &
         UNITS      = 'K s-1',                                         &
         DIMS       = MAPL_DimsHorzVert,                               &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)
    
    call MAPL_AddImportSpec(GC,                                       &
         SHORT_NAME = 'DQDT_BL',                                      &
         LONG_NAME  = 'tendency_of_spec_humidity_due_to_bound_layer',  &
         UNITS      = 'kg kg-1 s-1',                                  &
         DIMS       = MAPL_DimsHorzVert,                              &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'DNDDT ',                                      &
         LONG_NAME = 'total_liq_droplet_number_tendency_due_to_moist',       &
         UNITS     = 'm-3 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME= 'DNCDT ',                                      &
         LONG_NAME = 'total_ice_crystal_number_tendency_due_to_moist',       &
         UNITS     = 'm-3 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
         RC=STATUS  )
    VERIFY_(STATUS)
    
!-srf-gf-scheme
    IF(DEBUG_GF==1 .and. (CONVPAR_OPTION=='GF' .or. CONVPAR_OPTION=='BOTH')) THEN

       call MAPL_AddExportSpec(GC,                                     &
         SHORT_NAME = 'DTRDT_GF',                                    &
         LONG_NAME  = 'tendency_of_tracer_due_GF',        &
         UNITS      = 'kg kg-1 s-1',                                  &
         DIMS       = MAPL_DimsHorzVert,                              &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
       VERIFY_(STATUS)    
       call MAPL_AddExportSpec(GC,                                     &
         SHORT_NAME = 'DQDT_GF',                                    &
         LONG_NAME  = 'tendency_of_spec_humidity_due_GF',        &
         UNITS      = 'kg kg-1 s-1',                                  &
         DIMS       = MAPL_DimsHorzVert,                              &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec(GC,                                     &
         SHORT_NAME = 'DTDT_GF',                                    &
         LONG_NAME  = 'tendency_of_temp_due_GF',        &
         UNITS      = 'K s-1',                                  &
         DIMS       = MAPL_DimsHorzVert,                              &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec(GC,                                     &
         SHORT_NAME = 'MUPDP',                                    &
         LONG_NAME  = 'Mass_flux_deep_GF',        &
         UNITS      = 'kg m-2 s-1',                                  &
         DIMS       = MAPL_DimsHorzVert,                              &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec(GC,                                     &
         SHORT_NAME = 'MUPSH',                                    &
         LONG_NAME  = 'Mass_flux_shallow_GF',        &
         UNITS      = 'kg m-2 s-1',                                  &
         DIMS       = MAPL_DimsHorzVert,                              &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec(GC,                                     &
         SHORT_NAME = 'MUPMD',                                    &
         LONG_NAME  = 'Mass_flux_congestus_GF',        &
         UNITS      = 'kg m-2 s-1',                                  &
         DIMS       = MAPL_DimsHorzVert,                              &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
       VERIFY_(STATUS)
       !-2d
       call MAPL_AddExportSpec(GC,                                &
        SHORT_NAME         = 'MFDP',                              &
        LONG_NAME          = 'mass_flux_cloud_base_deep_GF',      &
        UNITS              = 'kg m-2 s-1',                        &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec(GC,                                &
        SHORT_NAME         = 'MFSH',                              &
        LONG_NAME          = 'mass_flux_cloud_base_shal_GF',      &
        UNITS              = 'kg m-2 s-1',                        &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec(GC,                                &
        SHORT_NAME         = 'MFMD',                              &
        LONG_NAME          = 'mass_flux_cloud_base_mid_GF',       &
        UNITS              = 'kg m-2 s-1',                        &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec(GC,                                &
        SHORT_NAME         = 'ERRDP',                             &
        LONG_NAME          = 'convection_code_deep_GF',                   & 
        UNITS              = '1',                                &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec(GC,                                &
        SHORT_NAME         = 'ERRSH',                             &
        LONG_NAME          = 'convection_code_shallow_GF',                   & 
        UNITS              = '1',                                &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
        VERIFY_(STATUS)
        call MAPL_AddExportSpec(GC,                               &
        SHORT_NAME         = 'ERRMD',                             &
        LONG_NAME          = 'convection_code_mid_GF',            & 
        UNITS              = '1',                                &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
        VERIFY_(STATUS)
        call MAPL_AddExportSpec(GC,                               &
        SHORT_NAME         = 'AA0',                               &
        LONG_NAME          = 'cloud work function 0',             & 
        UNITS              = 'J kg-1',                              &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
        VERIFY_(STATUS)

        call MAPL_AddExportSpec(GC,                               &
        SHORT_NAME         = 'AA1',                               &
        LONG_NAME          = 'cloud work function 1',             & 
        UNITS              = 'J kg-1',                              &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
        VERIFY_(STATUS)

        call MAPL_AddExportSpec(GC,                               &
        SHORT_NAME         = 'AA2',                               &
        LONG_NAME          = 'cloud work function 2',             & 
        UNITS              = 'J kg-1',                              &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
        VERIFY_(STATUS)

        call MAPL_AddExportSpec(GC,                               &
        SHORT_NAME         = 'AA3',                               &
        LONG_NAME          = 'cloud work function 3',             & 
        UNITS              = 'J kg-1',                              &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
        VERIFY_(STATUS)
	
        call MAPL_AddExportSpec(GC,                               &
        SHORT_NAME         = 'AA1_CIN',                           &
        LONG_NAME          = 'cloud work function CIN',           & 
        UNITS              = 'J kg-1',                              &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
        VERIFY_(STATUS)

        call MAPL_AddExportSpec(GC,                               &
        SHORT_NAME         = 'AA1_BL',                            &
        LONG_NAME          = 'Bound layer AA1',                   & 
        UNITS              = 'J kg-1 s-1',                            &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
        VERIFY_(STATUS)
	
        call MAPL_AddExportSpec(GC,                               &
        SHORT_NAME         = 'TAU_BL',                             &
        LONG_NAME          = 'Bound layer time scale',            & 
        UNITS              = 's',                                &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
        VERIFY_(STATUS)
	
        call MAPL_AddExportSpec(GC,                               &
        SHORT_NAME         = 'TAU_EC',                             &
        LONG_NAME          = 'cape removal time scale',           & 
        UNITS              = 's',                                &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
        VERIFY_(STATUS)

    ENDIF
!-srf-gf-scheme
!========================================================================================
!--kml--- activation for single-moment uphysics
        call MAPL_AddExportSpec(GC,                               &
        SHORT_NAME         = 'NACTL',                             &
        LONG_NAME          = 'activ aero # conc liq phase for 1-mom',              & 
        UNITS              = 'm-3',                               &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,     RC=STATUS  )
        VERIFY_(STATUS)

        call MAPL_AddExportSpec(GC,                               &
        SHORT_NAME         = 'NACTI',                             &
        LONG_NAME          = 'activ aero # conv ice phase for 1-mom',               & 
        UNITS              = 'm-3',                               &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,     RC=STATUS  )
        VERIFY_(STATUS)
!--kml--- activation for single-moment uphysics

!========================================================================================             


    !EOS


    ! Set the Profiling timers
    ! ------------------------

    call MAPL_TimerAdd(GC,name="DRIVER" ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,name="-PRE_RAS"    ,RC=STATUS)
    VERIFY_(STATUS)
!-srf-gf-scheme
    call MAPL_TimerAdd(GC,name="-GF"    ,RC=STATUS)
    VERIFY_(STATUS)
!-srf-gf-scheme    
    call MAPL_TimerAdd(GC,name="-RAS"    ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,name="--RAS_RUN"    ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,name="-POST_RAS"    ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,name="--UWSHCU"    ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,name="--UWSHCU_LOOP"    ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,name="--UWSHCU_LOOP1"    ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,name="-CLOUD"    ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,name="--CLOUD_RUN"    ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,name="---ACTIV"    ,RC=STATUS)
    VERIFY_(STATUS)    
    call MAPL_TimerAdd(GC,name="---CLDMACRO"    ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,name="---GFDL_CLDMICRO" ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,name="---MGMICRO"    ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,name="--CLOUD_DATA"    ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,name="---CLOUD_DATA_DEVICE"    ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,name="---CLOUD_DATA_CONST"    ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,name="--CLOUD_ALLOC"    ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,name="--CLOUD_DEALLOC"    ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,name="--USE_AEROSOL_NN1"   ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,name="--USE_AEROSOL_NN2"   ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,name="--USE_AEROSOL_NN3"   ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,name="-MISC1"    ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,name="-MISC2"    ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,name="-MISC3"    ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,name="--FLASH"    ,RC=STATUS)
    VERIFY_(STATUS)


    ! Set generic init and final methods
    ! ----------------------------------

    call MAPL_GenericSetServices    ( GC, RC=STATUS)
    VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)

  end subroutine SetServices

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!-!-!-!!!!!!!!!!!!!!!

  !EOP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !BOP

  ! !IROUTINE: Initialize -- Initialize method for the composite Moist Gridded Component

  ! !INTERFACE:

  subroutine Initialize ( GC, IMPORT, EXPORT, CLOCK, RC )

    ! !ARGUMENTS:

    type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
    type(ESMF_State),    intent(inout) :: IMPORT ! Import state
    type(ESMF_State),    intent(inout) :: EXPORT ! Export state
    type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
    integer, optional,   intent(  out) :: RC     ! Error code

    ! !DESCRIPTION: The Initialize method of the Moist Physics Gridded Component first 
    !   calls the Initialize method of the child Dynamics.  The Dynamics Initialize method will
    !   create the ESMF GRID, which will then be used to set the GRID associated with the
    !   SuperDyn Composite Component itself.  It should be noted that the 
    !   SuperDyn Initialize method also invokes the GEOS Topo Utility which creates all
    !   topography related quantities.

    !EOP


    ! ErrLog Variables

    character(len=ESMF_MAXSTR)          :: IAm
    integer                             :: STATUS
    character(len=ESMF_MAXSTR)          :: COMP_NAME

    ! Local derived type aliases

    type (MAPL_MetaComp),      pointer  :: MAPL
    type (ESMF_Grid )                   :: GRID
    type (ESMF_State),         pointer  :: GIM(:)
    type (ESMF_State),         pointer  :: GEX(:)
    type (ESMF_State)                   :: INTERNAL

    type (ESMF_Config)                  :: CF

    real, pointer, dimension(:,:,:)     :: Q, QLLS, QLCN, QILS, QICN, QRAIN, QSNOW, QGRAUPEL, QW
    real, dimension(:,:,:), pointer     :: PTR3

    integer  unit

    real DCS, QCVAR_, WBFFACTOR, NC_CST, NI_CST, NG_CST
    logical  :: nccons, nicons, ngcons, do_graupel
    integer  :: LM
 
    real(ESMF_KIND_R8)  Dcsr8, qcvarr8,  micro_mg_berg_eff_factor_in, ncnstr8, ninstr8, ngnstr8
    !=============================================================================

    ! Begin... 

    ! Get the target components name and set-up traceback handle.
    ! -----------------------------------------------------------
    Iam = "Initialize"
    call ESMF_GridCompGet ( GC, name=COMP_NAME, config=CF, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // trim(Iam)

    ! Get my internal MAPL_Generic state
    !-----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS )
    VERIFY_(STATUS)

    ! Call Generic Initialize for MOIST GC
    !----------------------------------------

    call MAPL_GenericInitialize ( GC, IMPORT, EXPORT, CLOCK, RC=STATUS )
    VERIFY_(STATUS)

    ! Get parameters from generic state.
    !-----------------------------------

    call MAPL_Get ( MAPL, LM=LM, GIM=GIM, GEX=GEX, INTERNAL_ESMF_STATE=INTERNAL, RC=STATUS )
    VERIFY_(STATUS)

    ! Inititialize cloud microphysics (Options: 1MOMENT, 2MOMENT or GFDL)
    !--------------------------------------------------------------
    call MAPL_GetResource( MAPL, CLDMICRO, Label="CLDMICRO:",  default="1MOMENT", RC=STATUS)
    VERIFY_(STATUS)
    LCLDMICRO = adjustl(CLDMICRO)=="1MOMENT" .or. &
                adjustl(CLDMICRO)=="2MOMENT" .or. &
                adjustl(CLDMICRO)=="GFDL"
    _ASSERT( LCLDMICRO, 'needs informative message' )
    if (adjustl(CLDMICRO)=="2MOMENT") then
      call MAPL_GetResource( MAPL, MGVERSION, Label="MGVERSION:",  default=0.0, RC=STATUS)
    endif
    call MAPL_GetResource( MAPL, DOSHLW, Label="DOSHLW:",  default=0, RC=STATUS)
 
    ! Inititialize QW Passive Tracer
    !-------------------------------

    call MAPL_GetPointer(INTERNAL, Q,        'Q'       , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QRAIN,    'QRAIN'   , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QSNOW,    'QSNOW'   , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QGRAUPEL, 'QGRAUPEL', RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QLLS,     'QLLS'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QLCN,     'QLCN'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QILS,     'QILS'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QICN,     'QICN'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QW,       'QW'      , RC=STATUS); VERIFY_(STATUS)

    QW = Q+QLLS+QLCN+QILS+QICN+QRAIN+QSNOW+QGRAUPEL

    if(adjustl(CLDMICRO)=="GFDL") then
       call gfdl_cloud_microphys_init()
       call WRITE_PARALLEL ("INITIALIZED GFDL microphysics in non-generic GC INIT")
    end if

    if(adjustl(CLDMICRO)=="2MOMENT") then

       call MAPL_GetResource( MAPL, DOGRAUPEL, Label="DOGRAUPEL:",  default=0, RC=STATUS)
                         do_graupel = .false.
       if (DOGRAUPEL/=0) do_graupel = .true.

       call MAPL_GetResource(MAPL, DCS, 'DCS:', default=350.0e-6, RC=STATUS )
       VERIFY_(STATUS)    
       Dcsr8 = DCS    

       call MAPL_GetResource(MAPL, QCVAR_,   'QCVAR:', DEFAULT= 2.0 ,RC=STATUS) !variance of the QL distribution     
       qcvarr8=QCVAR_

       call MAPL_GetResource(MAPL, WBFFACTOR,   'WBFFACTOR:', DEFAULT= 1.0 ,RC=STATUS) !variance of the QL distribution     
       micro_mg_berg_eff_factor_in = WBFFACTOR
       
       call MAPL_GetResource(MAPL, NC_CST,  'NC_CST:', DEFAULT= 0.0 ,RC=STATUS) !constant nd (set if greather than zero)     
       call MAPL_GetResource(MAPL, NI_CST,  'NI_CST:', DEFAULT= 0.0 ,RC=STATUS) !constant nd (set if greather than zero) 
       call MAPL_GetResource(MAPL, NG_CST,  'NG_CST:', DEFAULT= 0.0 ,RC=STATUS) !constant ng (set if greather than zero) 
       
       ncnstr8 = NC_CST
       if  (NC_CST .gt. 0.0)  nccons =.true.       
       ninstr8 = NC_CST
       if  (NI_CST .gt. 0.0)  nicons =.true.
       ngnstr8 = NC_CST
       if  (NG_CST .gt. 0.0)  ngcons =.true.
         
       if  (MGVERSION .gt. 1.0) then 
          call micro_mg_init(Dcsr8, do_graupel,  micro_mg_berg_eff_factor_in, &
                         nccons, nicons, ncnstr8, ninstr8, ngcons, ngnstr8)       
       else     
           call ini_micro(Dcsr8, micro_mg_berg_eff_factor_in, &
                          nccons, nicons, ncnstr8, ninstr8, qcvarr8)
       end if 
                         
       call aer_cloud_init()
       call WRITE_PARALLEL ("INITIALIZED MG in non-generic GC INIT")
    end if
    
!--kml
    call MAPL_GetResource(MAPL,INT_USE_AEROSOL_NN,'USE_AEROSOL_NN:',default=1, RC=STATUS )
    VERIFY_(STATUS)
    if (INT_USE_AEROSOL_NN == 0) then
       USE_AEROSOL_NN = .FALSE.
    else
       USE_AEROSOL_NN = .TRUE.
    end if
 
    if( (adjustl(CLDMICRO)=="1MOMENT" .or. adjustl(CLDMICRO)=="GFDL") .and. USE_AEROSOL_NN) then
       call aer_cloud_init()
       call WRITE_PARALLEL ("INITIALIZED aer_cloud_init for 1moment")
    end if
!--kml

!-srf-gf-scheme
    ! Inititialize configuration parameters for GF convection scheme 
    IF(ADJUSTL(CONVPAR_OPTION) == 'GF') THEN
        ! Inititialize parameters of convection scheme GF
        call MAPL_GetResource(MAPL, icumulus_gf(deep), 'DEEP:'     , default=1, RC=STATUS )
        VERIFY_(STATUS)
        call MAPL_GetResource(MAPL, icumulus_gf(shal), 'SHALLOW:'  , default=1, RC=STATUS )
        VERIFY_(STATUS) 
        call MAPL_GetResource(MAPL, icumulus_gf(mid) , 'CONGESTUS:', default=1, RC=STATUS )
        VERIFY_(STATUS)
        call MAPL_GetResource(MAPL, closure_choice(deep), 'CLOSURE_DEEP:'     , default= 0, RC=STATUS )
        VERIFY_(STATUS)
        call MAPL_GetResource(MAPL, closure_choice(shal), 'CLOSURE_SHALLOW:'  , default= 7, RC=STATUS )
        VERIFY_(STATUS) 
        call MAPL_GetResource(MAPL, closure_choice(mid) , 'CLOSURE_CONGESTUS:', default= 3, RC=STATUS )
        VERIFY_(STATUS)
        call MAPL_GetResource(MAPL,USE_TRACER_TRANSP   ,'USE_TRACER_TRANSP:',default= 1, RC=STATUS )
        VERIFY_(STATUS)
        call MAPL_GetResource(MAPL,USE_TRACER_SCAVEN   ,'USE_TRACER_SCAVEN:',default= 2, RC=STATUS )
        VERIFY_(STATUS)
        call MAPL_GetResource(MAPL,USE_SCALE_DEP       ,'USE_SCALE_DEP:'    ,default= 1, RC=STATUS )
        VERIFY_(STATUS)
        call MAPL_GetResource(MAPL,DICYCLE             ,'DICYCLE:'          ,default= 1, RC=STATUS )
        VERIFY_(STATUS)
        call MAPL_GetResource(MAPL,TAU_DEEP   ,'TAU_DEEP:',default= 5400.0, RC=STATUS )
        VERIFY_(STATUS)
        call MAPL_GetResource(MAPL,TAU_MID    ,'TAU_MID:' ,default= 3600.0, RC=STATUS )
        VERIFY_(STATUS)

        call MAPL_GetResource(MAPL, USE_FLUX_FORM   ,'USE_FLUX_FORM:'   ,default= 1,  RC=STATUS );VERIFY_(STATUS)
        call MAPL_GetResource(MAPL, USE_FCT         ,'USE_FCT:'	        ,default= 0,  RC=STATUS );VERIFY_(STATUS)
        call MAPL_GetResource(MAPL, USE_TRACER_EVAP ,'USE_TRACER_EVAP:' ,default= 1,  RC=STATUS );VERIFY_(STATUS)
        call MAPL_GetResource(MAPL, ALP1            ,'ALP1:'            ,default= 1., RC=STATUS );VERIFY_(STATUS)

       ! IF(ADJUSTL(AERO_PROVIDER) == 'GOCART.data' .AND. USE_TRACER_TRANSP == 1) THEN
       !    call WRITE_PARALLEL ("AERO_PROVIDER: GOCART.data detected, disabling tracer transport for GF")
       !    USE_TRACER_TRANSP = 0
       ! END IF

       ! IF(ADJUSTL(AERO_PROVIDER) == 'GOCART.data' .AND. USE_TRACER_SCAVEN == 1) THEN
       !    call WRITE_PARALLEL ("AERO_PROVIDER: GOCART.data detected, disabling scavenging for GF")
       !    USE_TRACER_SCAVEN = 0
       ! END IF

       !IF(USE_TRACER_TRANSP == 1) THEN
       !    call WRITE_PARALLEL ("GEOS_MoistGridCompMod: GF tracer transport detected, disabling transport in GOCART")
       !    call Disable_Convection
       ! ELSE
       !    call WRITE_PARALLEL ("GEOS_MoistGridCompMod: Using GOCART for tracer transport")
       !END IF

    ENDIF
!-srf-gf-scheme
    call MAPL_GetResource(MAPL,USE_TRACER_TRANSP_UW,'USE_TRACER_TRANSP_UW:',default= 1, RC=STATUS )
    VERIFY_(STATUS)
    !IF(ADJUSTL(AERO_PROVIDER) == 'GOCART.data' .AND. USE_TRACER_TRANSP_UW == 1) THEN
    !       call WRITE_PARALLEL ("AERO_PROVIDER: GOCART.data detected, disabling tracer transport for UW")
    !       USE_TRACER_TRANSP_UW = 0
    !END IF

    ! All done
    !---------

    RETURN_(ESMF_SUCCESS)
  end subroutine Initialize


  !===================================================================================

  !BOP

  ! !IROUTINE: RUN -- Run method for the CONVECT component

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

    character(len=ESMF_MAXSTR)      :: IAm
    integer                         :: STATUS
    character(len=ESMF_MAXSTR)      :: COMP_NAME

    type( ESMF_VM )                 :: VMG

    ! Local derived type aliases

    type (MAPL_MetaComp), pointer   :: STATE
    type (ESMF_Config  )            :: CF
    type (ESMF_State   )            :: INTERNAL
    type (ESMF_Alarm   )            :: ALARM

    type (RASPARAM_TYPE)            :: RASPARAMS
    type (CLDPARAM_TYPE)            :: CLDPARAMS
    type (SHLWPARAM_TYPE)            :: SHLWPARAMS

    ! Local variables

    integer                         :: IM,JM,LM
    real, pointer, dimension(:,:)   :: LONS
    real, pointer, dimension(:,:)   :: LATS

    !=============================================================================

    ! Begin... 

    ! Get my name and set-up traceback handle
    ! ---------------------------------------

    Iam = 'Run'
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, CONFIG=CF, VM=VMG, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // Iam

    ! Get my internal MAPL_Generic state
    !-----------------------------------

    call MAPL_GetObjectFromGC ( GC, STATE, RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_TimerOn (STATE,"TOTAL")

    ! Get parameters from generic state.
    !-----------------------------------

    call MAPL_Get( STATE, IM=IM, JM=JM, LM=LM,   &
         RUNALARM = ALARM,             &
         CF       = CF,                &
         LONS     = LONS,              &
         LATS     = LATS,              &
         INTERNAL_ESMF_STATE=INTERNAL, &
         RC=STATUS )
    VERIFY_(STATUS)




    ! If its time, calculate convective tendencies
    ! --------------------------------------------

    if ( ESMF_AlarmIsRinging( ALARM, RC=status) ) then
       VERIFY_(STATUS)
       call ESMF_AlarmRingerOff(ALARM, RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_TimerOn(STATE,"DRIVER")
       call MOIST_DRIVER(IM,JM,LM, RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_TimerOff(STATE,"DRIVER")
    endif
    VERIFY_(STATUS)

    call MAPL_TimerOff(STATE,"TOTAL")

    RETURN_(ESMF_SUCCESS)

  contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine MOIST_DRIVER(IM,JM,LM, RC)
      integer,           intent(IN ) :: IM, JM, LM
      integer, optional, intent(OUT) :: RC

      type RAS_Tracer_T
         real, pointer :: Q(:,:,:) => null()
      end type RAS_Tracer_T

      !  Locals

      character(len=ESMF_MAXSTR)      :: IAm
      character(len=ESMF_MAXSTR)      :: NAME
      integer                         :: STATUS

      type (ESMF_TimeInterval)        :: TINT
      type (RAS_Tracer_T ), pointer   :: TRPtrs (:)

      !  LOCAL COPY OF VARIABLES

      real, pointer, dimension(:,:,:) :: DQDTCN, DTHDTCN,DQCDTCN,DTDTFRIC

      integer :: DOCLDMACRO
      real, pointer, dimension(:,:,:) :: UMF_SC, MFD_SC, WUP_SC, QTUP_SC, &
                                         THLUP_SC, THVUP_SC, UUP_SC, VUP_SC 
      real, pointer, dimension(:,:,:) :: QCU_SC, QLU_SC, QIU_SC
      real, pointer, dimension(:,:,:) :: DTHDT_SC, DQVDT_SC, DQRDT_SC, DQSDT_SC, &
                                         DQIDT_SC, DQLDT_SC, DQCDT_SC, CUFRC_SC
      real, pointer, dimension(:,:,:) :: DUDT_SC, DVDT_SC, &
                                         ENTR_SC, DETR_SC, XC_SC,QLDET_SC, &
                                         QIDET_SC, QLENT_SC, QIENT_SC, &
                                         QLSUB_SC, QISUB_SC, SC_NDROP, SC_NICE
      real, pointer, dimension(:,:  ) :: CBMF_SC
      real, pointer, dimension(:,:  ) :: CIN_SC, PINV_SC, PLCL_SC, PLFC_SC, &
                                         PREL_SC, PBUP_SC
      real, pointer, dimension(:,:  ) :: WLCL_SC, QTSRC_SC, THLSRC_SC, &
                                         THVLSRC_SC, TKEAVG_SC, CLDTOP_SC, CUSH
      real, pointer, dimension(:,:  ) :: CNT_SC, CNB_SC, CLDBASEHGT

! Diagnostic output from ADG PDF
    real, dimension(:,:,:),pointer     :: PDF_A,      &
                                          PDF_AX
#ifdef PDFDIAG
    real, dimension(:,:,:),pointer     :: PDF_SIGW1,  &
                                          PDF_SIGW2,  &
                                          PDF_W1,     &
                                          PDF_W2,     &
                                          PDF_SIGTH1, &
                                          PDF_SIGTH2, &
                                          PDF_TH1,    &
                                          PDF_TH2,    &
                                          PDF_SIGQT1, &
                                          PDF_SIGQT2, &
                                          PDF_QT1,    &
                                          PDF_QT2,    &
                                          PDF_RQTTH,  &
                                          PDF_RWTH,   &
                                          PDF_RWQT
#endif

! Inputs for ADG PDF
    real, dimension(:,:,:),pointer     :: HL2,       &
                                          HL3,       &
                                          QT2,       &
                                          QT3,       &
                                          W2,        &
                                          W3,        &
                                          HLQT,      &
                                          WQT,       &
                                          WQL,       &
                                          WHL,       &
                                          EDMF_FRC

    real, dimension(:,:,:),pointer     :: WTHV2,WTHV2_RAD


      real, pointer, dimension(:,:,:) :: CNV_DQLDT            , &
           CNV_MF0              , &
           CNV_MFD              , &
           CNV_MFC              , &
           CNV_UPDF             , &
           CNV_CVW              , &
           CNV_QC               , &
           RAD_CF               , &
           RAD_QL               , &
           RAD_QI               , &
           RAD_QR               , &
           RAD_QS               , &
           RAD_QG               , &
           RAD_QV               , &
           CLDNCCN              , &                
           CLDREFFR             , &                
           CLDREFFS             , &
           CLDREFFG             , &                                
           CLDREFFL             , &                
           CLDREFFI                             


      real, pointer, dimension(:,:,:) :: T, PLE, U, V, W, TH
      real, pointer, dimension(:,:)   :: TROPP
      real, pointer, dimension(:,:,:) :: DQDT, UI, VI, WI, TI, KH, TKE
      real, pointer, dimension(    :) :: PREF
      real, pointer, dimension(:,:,:) :: Q, QRAIN, QSNOW, QGRAUPEL, QLLS, QLCN, CLLS, CLCN, BYNCY, QILS, QICN, QCTOT,QITOT,QLTOT
      real, pointer, dimension(:,:,:) :: QPTOTLS, QRTOT, QSTOT,  CFLIQ, CFICE !DONIF

      real, pointer, dimension(:,:,:) :: NCPL,NCPI, NRAIN, NSNOW, NGRAUPEL
 
 
      real, pointer, dimension(:,:,:) :: DBZ
      real, pointer, dimension(:,:  ) :: DBZ_MAX

      real, pointer, dimension(:,:  ) :: PRCP_RAIN, PRCP_SNOW, PRCP_ICE, PRCP_GRAUPEL
      real, pointer, dimension(:,:  ) :: LS_PRCP,CN_PRCP,AN_PRCP,SC_PRCP,TT_PRCP,ER_PRCP,FILLNQV
      real, pointer, dimension(:,:  ) :: HOURNORAIN
      integer                         :: YEAR, MONTH, DAY, HR, SE, MN
      type (ESMF_Time)                :: CurrentTime
      real, pointer, dimension(:,:  ) :: LS_ARF, CN_ARF, AN_ARF, SC_ARF
      real, pointer, dimension(:,:  ) :: PTYPE,FRZR,ICE,SNR,PRECU,PRELS,TS,SNOMAS,FRLANDICE,FRLAND
      real, pointer, dimension(:,:  ) :: IWP,LWP,CWP,TPW,CAPE,ZPBLCN,INHB,ZLCL,ZLFC,ZCBL,CCWP , KPBLIN, KPBLSC
      real, pointer, dimension(:,:  ) :: TVQ0,TVQ1,TVE0,TVE1,TVEX,DCPTE, TVQX2, TVQX1, CCNCOLUMN, NDCOLUMN, NCCOLUMN  !DONIF
      real, pointer, dimension(:,:,:,:) :: XHO
      real, pointer, dimension(:,:  ) ::  MXDIAM, RH600, Q600, QCBL, QRATIO, CNV_FRC

      real, pointer, dimension(:,:  ) :: RAS_TIME, RAS_TRG, RAS_TOKI, RAS_PBL, RAS_WFN 
      real, pointer, dimension(:,:,:) :: RAS_ALPHA, RAS_TAU
      real, pointer, dimension(:,:  ) :: STOCH_CNV, SIGMA_DEEP, SIGMA_MID

!!$      real, pointer, dimension(:,:,:) :: LIQANMOVE, ICEANMOVE, DANCLD, DLSCLD
!!$      real, pointer, dimension(:,:,:) :: CURAINMOVE, CUSNOWMOVE
      real, pointer, dimension(:,:,:) :: PFLCNMOVE, PFICNMOVE
      real, pointer, dimension(:,:)   :: CU2DRAINMOVE, CU2DSNOWMOVE

      real, pointer, dimension(:,:  ) :: KEMST,KEMST2
      real, pointer, dimension(:,:,:) :: QSN, QRN, QPLS

      logical, pointer                :: IS_FRIENDLY(:)
      logical, pointer                :: IS_WEIGHTED(:)

      !     Tracer scavenging coefficients
      !     FSCAV is the Fraction of tracer scavanged per km (=0 no scavenging, =1 full scavenging)
      real, pointer, dimension(:    ) :: FSCAV_, &  ! holding array for all tracers
                                         FSCAV      ! container for friendly to moist tracers
      real, dimension(4) :: Vect_Hcts

      ! Aerosol convective scavenging internal pointers (2D column-averages);  must deallocate!!!
      ! CAR 
      real, pointer, dimension(:,: )                           :: DDUDT, &
           DSSDT, DOCDT, DBCDT, DSUDT,  DNIDT, DNH4ADT, DNH3DT, DBRCDT, DDUDTcarma, DSSDTcarma
      real, pointer, dimension(:,: )                           :: DDU2gDT, &
           DSS2gDT, DOC2gDT, DBC2gDT, DSU2gDT,  DNI2gDT, DNH4A2gDT, DNH32gDT, DBRC2gDT
      character(len=ESMF_MAXSTR)                               :: QNAME,  CNAME, ENAME
      character(len=ESMF_MAXSTR), pointer, dimension(:)        :: QNAMES, CNAMES
      integer                                                  :: ind

      integer                                                  :: i_src_mode
      integer                                                  :: i_dst_mode

      real, pointer, dimension(:,:)   :: PGENTOT , PREVTOT

      !Record vars at top pf moist
      real, pointer, dimension(:,:,:) :: Ux0, Vx0, THx0, KHx0, W_DIAG
      real, pointer, dimension(:,:)   :: TSx0, FRLANDx0
      real, pointer, dimension(:,:,:) :: Qx0, Qx1, QLLSx0, QLCNx0, CLLSx0, CLCNx0, QILSx0, QICNx0, QCLSX0, QCCNX0

      ! MATMAT Exports for pre-ras inputs for RAStest
      real, pointer, dimension(:,:,:) :: THOI,QHOI,QSSI,DQSI
      real, pointer, dimension(:,:,:) :: PLEI
      real, pointer, dimension(:,:  ) :: TPERTI,KCBLI

      !Extra outputs

      real, pointer, dimension(:,:,:) :: DQRL, DQRC, DQLDT, DQIDT, KEDISS
      real, pointer, dimension(:,:,:) :: RH1, RHX, RH2, XQLLS, XQLCN, XCLLS, XCLCN, XQILS, XQICN
      real, pointer, dimension(:,:,:) :: REV_CN, REV_AN, REV_LS, REV_SC, RSU_SC
      real, pointer, dimension(:,:,:) :: RSU_CN, RSU_AN, RSU_LS,    ALPHT, ALPH1, ALPH2
      real, pointer, dimension(:,:,:) :: ENTLAM
      real, pointer, dimension(:,:,:) :: KHX
      real, pointer, dimension(:,:,:) :: SLX, TX, QTX
      real, pointer, dimension(:,:  ) :: DTSX

      real, pointer, dimension(:,:,:) :: REVSU_CN, REVSU_LSAN
      real, pointer, dimension(:,:  ) :: CNV_FREQ, CNV_BASEP, CNV_TOPP


      real, pointer, dimension(:,:,:) :: ACLL_CN, ACLL_AN, ACLL_LS, ACLL_SC
      real, pointer, dimension(:,:,:) :: ACIL_CN, ACIL_AN, ACIL_LS, ACIL_SC
      real, pointer, dimension(:,:,:) :: ACR_TOT

      real, pointer, dimension(:,:,:) :: PFL_CN, PFL_AN, PFL_LS, PFL_SC
      real, pointer, dimension(:,:,:) :: PFI_CN, PFI_AN, PFI_LS, PFI_SC
      real, pointer, dimension(:,:,:) :: PFI_LSAN, PFL_LSAN
      real, pointer, dimension(:,:,:) :: DPDTMST

      real, pointer, dimension(:,:,:) :: DlPDF, DiPDF, DlFIX, DiFIX, FRZ_PP
      real, pointer, dimension(:,:,:) :: AUT, EVAPC, SDM, SUBLC, CFPDF, CFPDFX, RHCLR
      real, pointer, dimension(:,:,:) :: VFALLICE_AN, VFALLICE_LS
      real, pointer, dimension(:,:,:) :: VFALLWAT_AN,VFALLWAT_LS
      real, pointer, dimension(:,:,:) :: VFALLSN_AN,VFALLSN_LS,VFALLSN_CN,VFALLSN_SC
      real, pointer, dimension(:,:,:) :: VFALLRN_AN,VFALLRN_LS,VFALLRN_CN,VFALLRN_SC
      real, pointer, dimension(:,:,:) :: FRZ_TT, DCNVL, DCNVi,QSATi,QSATl,RCCODE,TRIEDLV,QVRAS

      real, pointer, dimension(:,:  ) :: LFR_GCC

      !Whether to guard against negatives
      logical                         :: RAS_NO_NEG

      !Trajectory for Moist TLM/ADJ
      real, pointer, dimension(:,:,:) :: TH_moist, Q_moist
      real, pointer, dimension(:,:  ) :: KCBL_moist, TS_moist, ctop_moist, KHu_moist, KHl_moist


      ! MATMAT Additional after-RAS exports
      real, pointer, dimension(:,:,:) :: THRAS,URAS,VRAS
      !AMM sync t,q extra diagnostic
      real, pointer, dimension(:,:,:) :: THMOIST, SMOIST
      real, dimension(IM,JM,0:LM)     :: geopenew

      ! vapor-to-liquid  [ *ALHL ]
      real, pointer, dimension(:,:) :: PDFLZ,CNVLZ,CNVRNZ

      ! liquid-to-vapor  [ *ALHL ]
      real, pointer, dimension(:,:) :: EVPCZ, EVPPZ

      ! vapor-to-ice  [ *ALHS ]
      real, pointer, dimension(:,:) :: PDFIZ,CNVIZ 

      ! ice-to-vapor  [ *ALHS ]
      real, pointer, dimension(:,:) :: SUBCZ, SUBPZ

      ! liquid-to-ice  [ *ALHF ]
      real, pointer, dimension(:,:) :: FRZCZ, COLLIZ, FRZPZ

      ! liquid-to-liquid  [ * 0 ]
      real, pointer, dimension(:,:) :: AUTZ, COLLLZ

      ! ice-to-ice  [ * 0 ]
      real, pointer, dimension(:,:) :: SDMZ

!-srf-gf-scheme    
      !~ character(LEN=ESMF_MAXSTR) :: CONVPAR_OPTION
      real, pointer, dimension(:,:  )       :: AREA
      real, pointer, dimension(:,:,:)       :: DQVDTDYN
      real, pointer, dimension(:,:,:)       :: DTDTDYN  
      real, pointer, dimension(:,:,:)       :: QV_DYN_IN
      real, pointer, dimension(:,:,:)       :: U_DYN_IN
      real, pointer, dimension(:,:,:)       :: V_DYN_IN
      real, pointer, dimension(:,:,:)       :: T_DYN_IN
      real, pointer, dimension(:,:,:)       :: PLE_DYN_IN
      real, pointer, dimension(:,:,:)       :: DTDT_BL
      real, pointer, dimension(:,:,:)       :: DQDT_BL
      real, pointer, dimension(:,:,:)       :: DQDT_GF,DTDT_GF,MUPDP,MUPSH,MUPMD,DTRDT_GF
      real, pointer, dimension(:,:  )       :: USTAR,TSTAR,QSTAR,T2M,Q2M,TA,QA,SH,EVAP,PHIS
      real, pointer, dimension(:,:  )       :: MFDP,MFSH,MFMD,ERRDP,ERRSH,ERRMD
      real, pointer, dimension(:,:  )       :: AA0,AA1,AA2,AA3,AA1_BL,AA1_CIN,TAU_BL,TAU_EC
      real, pointer, dimension(:,:,:)       :: RSU_CN_GF,REV_CN_GF,PFL_CN_GF,PFI_CN_GF
!-srf-gf-scheme    
!--kml--- activation for single-moment uphysics
      real, pointer, dimension(:,:,:)       :: NACTL,NACTI
!--kml--- activation for single-moment uphysics

      ! Aerosol-Cloud interactions

      real, pointer, dimension(:,:,:) :: SMAXL, WSUB, CCN01, CCN04, CCN1, SMAXI, & !DONIF
           CDNC_NUC, INC_NUC, NCPL_VOL, NCPI_VOL, SO4, &
           ORG, DUST, SEASALT, BCARBON, NHET_NUC, NLIM_NUC, SAT_RAT, &
           DQVDT_micro, DQIDT_micro, DQLDT_micro, &
           DQRDT_micro, DQSDT_micro, DQGDT_micro, DQADT_micro, &
           DUDT_micro, DVDT_micro, DTDT_micro, RL_MASK, RI_MASK, KAPPA, RHICE, RHLIQ, SC_ICE, &
           NHET_IMM, NHET_DEP,  &
           DUST_IMM, DUST_DEP, SCF, SCF_ALL, SIGW_GW, SIGW_CNV, & 
           SIGW_TURB, SIGW_RC, CNV_FICE, &
           CNV_NDROP, CNV_NICE, RHCmicro, DNHET_IMM,  &
            BERG, BERGS, MELT,  DNHET_CT, QCRES, DT_RASP, &
           DTDT_macro, DQVDT_macro, DQLDT_macro, DQIDT_macro, DQADT_macro, &
           FRZPP_LS, SNOWMELT_LS, QIRES, AUTICE, PFRZ, DNCNUC, DNCHMSPLIT, DNCSUBL, &
           DNCAUTICE, DNCACRIS, DNDCCN, DNDACRLS, DNDEVAPC, DNDACRLR, DNDAUTLIQ, & 
           DNDCNV, DNCCNV, DTDT_moist, DTDTCN, ALPHA_RAS, DNDDT, DNCDT



      real, dimension(IM,JM) :: OMEGA500_X

      real, pointer, dimension(:,:,:) :: LWC, IWC, CIRRUS_FREQ_HOM, CIRRUS_FREQ_HET
      real, pointer, dimension(:,:) ::  CLDREFFI_TOP, CLDREFFL_TOP, NCPL_TOP, NCPI_TOP, QCVAR_EXP, &
                                        EIS, LTS, NCPL_CLDBASE


      real, allocatable, dimension(:,:,:)  ::  FQAl, FQAI, QCNTOT, FQA
      real, dimension(IM,JM,LM)  ::  QTOT, QL_TOT, QI_TOT, CFX, QRAIN_CN, QSNOW_CN,  &
                                                 CNV_NDROP_X, CNV_NICE_X, CNV_MFD_X, CNV_DQLDT_X, &
                                                 CNV_PRC3_X ,  CNV_UPDF_X  
       
      real, dimension(IM,JM) :: CLDBASEx
      real, dimension(IM,JM)  :: ZPBL
      integer, dimension(IM,JM)  :: KMIN_TROP
      real, parameter :: r_air = 3.47d-3 !m3 Pa kg-1K-1
      real, parameter :: pmin_trop = 10.0 !mbar minimum pressure to do cloud microphysics

      integer,  parameter :: ncolmicro = 1
      integer,  parameter :: nmods_gocart = 15

!!!!! 2-moment muphys declarations !!!!!!!!!!!!!!!!!!!!!!!!
!!! mmicro_pcond real*8 scalars and arrays 
      
      REAL, dimension(1,1:LM) :: SCICE_tmp, FQA_tmp, ALPH_tmp
      real(ESMF_KIND_R8), dimension(1,1:LM)   ::  so4x, seasaltx, dustx, orgx, bcx 
        real(ESMF_KIND_R8), dimension(1,1:LM)  ::ttendr8, qtendr8, cwtendr8, &
           cldor8,  rpdelr8, zmr8, omegr8, rhdfdar8, rhu00r8, ficer8 , &
            ndropr8, nimmr8
      
         real(ESMF_KIND_R8), dimension(1,LM+1)  :: pintr8, kkvhr8
 
       real(r8), dimension(1,1:LM)  ::                                                     &  
                             ter8,                          qvr8,                              &
                             qcr8,                          qir8,                          &
                             ncr8,                          nir8,                          &
                             qrr8,                          qsr8,                          &
                             nrr8,                          nsr8,                          &
                             qgr8,                          ngr8,                         &
                             relvarr8,                      accre_enhanr8,                  &
                             plevr8,                       pdelr8,                         &
                             cldfr8,               liqcldfr8,            icecldfr8,  qsatfacr8,          &
                             qcsinksum_rate1ordr8,                                         &
                             naair8,                         npccninr8,                        &
                             tlatr8,                         qvlatr8,                        &
                             qctendr8,                       qitendr8,                       &
                             nctendr8,                       nitendr8,                       &
                             qrtendr8,                       qstendr8,   qgtendr8,                     &
                             nrtendr8,                       nstendr8,   ngtendr8,                   &
                             effcr8,               effc_fnr8,            effir8,               &
                             sadicer8,                       sadsnowr8,                      &                            
                             nevaprr8,                       evapsnowr8,                     &
                             am_evp_str8,                                                  &
                             prainr8,                        prodsnowr8,                     &
                             cmeoutr8,                       deffir8,                        &
                             pgamradr8,                      lamcradr8,                      &
                             qsoutr8,                        dsoutr8,                        &
                             qgoutr8,     ngoutr8,           dgoutr8,                        &
                             qroutr8,          &
                             reff_rainr8,                    reff_snowr8, reff_graur8,        &
                             qcsevapr8,            qisevapr8,            qvresr8,              &
                             cmeioutr8,            vtrmcr8,              vtrmir8,              &
                             umrr8,                          umsr8,                          &
                             umgr8,                          qgsedtendr8,                    &    
                             qcsedtenr8,                     qisedtenr8,                     &
                             qrsedtenr8,                     qssedtenr8,                     &
                             praor8,                       prcor8,                       &
                             mnucccor8,          mnucctor8,          msacwior8,          &
                             psacwsor8,          bergsor8,           bergor8,            &
                             meltor8,                      homoor8,                      &
                             qcresor8,           prcior8,            praior8,            &
                             qirestotr8,           mnuccrtotr8,          mnuccritotr8, pracstotr8,           &                           
                             meltsdtr8,         frzrdtr8,          mnuccdor8,          &
                             pracgtotr8,           psacwgtotr8,          pgsacwtotr8,          &
                             pgracstotr8,          prdgtotr8,           &
                             qmultgtotr8,          qmultrgtotr8,         psacrtotr8,           &
                             npracgtotr8,          nscngtotr8,           ngracstotr8,          &
                             nmultgtotr8,          nmultrgtotr8,         npsacwgtotr8,         & 
                             nroutr8,                            nsoutr8,                        &
                             reflr8,               areflr8,              areflzr8,             &
                             freflr8,              csrflr8,              acsrflr8,             &
                             fcsrflr8,                       rercldr8,                       &
                             ncair8,                         ncalr8,                         &
                             qrout2r8,                       qsout2r8,                       &
                             nrout2r8,                       nsout2r8,                       &
                             drout2r8,                       dsout2r8,                       &
                             qgout2r8,     ngout2r8,           dgout2r8,  freqgr8,   &
                             freqsr8,                        freqrr8,                        &
                             nficer8,                        qcratr8,                        &
!                             errstring, & ! Below arguments are "optional" (pass null pointers to omit).
                             tnd_qsnow,          tnd_nsnow,          re_ice,    &
                             prer_evap, &
                             frzimmr8,             frzcntr8,              frzdepr8,  & ! contact is not passed since it depends on the droplet size dist
                             nsootr8, rnsootr8,  & ! soot for contact IN
                             npccnor8, npsacwsor8,npraor8,nsubcor8, nprc1or8, &  ! Number tendencies for liquid
                             npraior8, nnucctor8, nnucccor8, nnuccdor8, nsubior8, nprcior8, nsacwior8, &
                             
                             mnuccror8,pracsor8, qiresor8, rate1ord_cw2pr, & !only MG1
                             
                             
                             sc_icer8, nhet_immr8, dnhet_immr8, nhet_depr8,  &  ! activation
                              dust_immr8, dust_depr8,  dpre8, npre8
     
    
      real(ESMF_KIND_R8), dimension(1,1:LM,10)  ::    rndstr8,naconr8  !Assume maximum 5 dust bins
      real(ESMF_KIND_R8), dimension(1,LM+1)  :: rflxr8, sflxr8, lflxr8, iflxr8, gflxr8                
      real(ESMF_KIND_R8), dimension(1)       :: prectr8, precir8
      real(ESMF_KIND_R8)                     :: disp_liu, ui_scale, & 
           dcrit, tfreez, qcvar8, ts_autice, dcsr8, qcvarr8, scale_ri, mtimesc, urscale
          

       integer :: num_steps_micro,  pcnst, n_modes, kbmin, kcldtop, kcldbot , NAUX, kcldtopcvn, nbincontactdust, index    


      ! Aerosol cloud interactions internal variables 
      type(ESMF_State)                :: aero_aci
      character(len=ESMF_MAXSTR)      :: aci_field_name
      real, pointer, dimension(:,:)   :: aci_ptr_2d
      real, pointer, dimension(:,:,:) :: aci_ptr_3d
      real, pointer, dimension(:,:,:) :: aci_num
      real, pointer, dimension(:,:,:) :: aci_dgn
      real, pointer, dimension(:,:,:) :: aci_sigma
      real, pointer, dimension(:,:,:) :: aci_density
      real, pointer, dimension(:,:,:) :: aci_hygroscopicity
      real, pointer, dimension(:,:,:) :: aci_f_dust
      real, pointer, dimension(:,:,:) :: aci_f_soot
      real, pointer, dimension(:,:,:) :: aci_f_organic
      logical                         :: implements_aerosol_activation_properties
      character(len=ESMF_MAXSTR), allocatable, dimension(:) :: aero_aci_modes
      integer                         :: ACI_STATUS


      type  (AerProps), dimension (IM, JM, LM) :: AeroProps !Storages aerosol properties for activation 
      type  (AerProps) :: AeroAux, AeroAux_b


      real, pointer, dimension(:,:,:)    :: AERTEMP
      real, pointer, dimension(:,:)      :: TAUGWX, TAUGWY, TAUX, TAUY, TAUOROX, TAUOROY    
      

      real, pointer, dimension(:,:, :)   :: OMEGA, RADLW, RADSW, WSUB_NATURE
      real, pointer, dimension(:,:, :)      :: ALH   

      type (ESMF_Field)                  :: AERFIELD
      character(len=ESMF_MAXSTR)         :: AERNAME

      real                               ::  T_ICE_ALL , T_ICE_MAX, USURF, USURF2
      real (ESMF_KIND_R8), dimension(:), allocatable ::  aer_aux
      real (ESMF_KIND_R8), dimension(1, 1:LM) :: wparc, smaxliq, atot, &
           smaxicer8, nheticer8, incr8, swparc, &
           nhetr8, nlimicer8, qilsr8, wparc_gw, wparc_ls, &
           wparc_turb, wparc_cnv, lc_turb, rad_cooling, wparc_rc, &
	       uwind_gw, wparc_cgw, pfrz_inc_r8


      real (ESMF_KIND_R8), dimension(3)       :: ccn_diag
      integer, dimension(IM, JM)              :: KCT
            
      ! Subgrid velocity parameterization 

      real, dimension (1, 1:LM) :: tm_gw, pm_gw, nm_gw, rho_gw, theta_tr, khaux, qcaux, dummyW , &
            c2_gw, fcn, cfaux    
	   
      real, dimension (1, 0:LM) :: pi_gw, rhoi_gw, ni_gw, ti_gw
      real                      :: maxkhpbl, tausurf_gw, overscale, fracover, cfc_aux 
      real (ESMF_KIND_R8)       :: tauxr8, fsoot_drop, fdust_drop, sigma_nuc_r8, rh1_r8, frachet_dust, &
                                   frachet_bc, frachet_org, frachet_ss
      logical                   :: ismarine, is_stable, use_average_v                  
      real                      :: Nct, Wct, DX, ksa1, Xscale



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !==================================================



      !Outputs for downdraft
      real, pointer, dimension(:,:,:) :: DDF_DQDT , DDF_DTDT , DDF_MFC , DDF_MUPH , DDF_RH1 , DDF_RH2 , DDF_BYNC
      real, pointer, dimension(:,:,:) :: DDF_QVC ,  DDF_TC
      real, pointer, dimension(:,:)   :: DDF_ZSCALE
      real, pointer, dimension(:)     :: DDF_DQDTz, DDF_DTDTz, DDF_MFCz, DDF_MUPHz, DDF_RH1z, DDF_RH2z, DDF_BYNCz
      real, pointer, dimension(:)     :: DDF_QVCz,  DDF_TCz
      real                            :: DDF_ZSCALEz

      type (ESMF_Field)               :: FIELD
      type (ESMF_FieldBundle)         :: TR
      type (ESMF_FieldBundle)         :: TRI

      real                            :: DT, DT_MOIST
      real(ESMF_KIND_R8)              :: DT_R8

      real                            :: PMIN_DET,AUTOC_CN_LAND,AUTOC_CN_OCN, LCCIRRUS, & 
           UISCALE, AUTO_CNV, SS_SCALE, REEVAP_MICRO, LIU_MU, TFRZ, & 
           NPRE_FRAC, QCVAR, ZPBLMAXLL, TMAXLL,  &
           LTS_LOW, LTS_UP, MIN_EXP, BKGTAU, DCRIT_, USE_AV_V, TS_AUTO_ICE, CCN_PARAM, IN_PARAM, &
	   FDROP_DUST, FDROP_SOOT,  USE_NATURE_WSUB, SIGMA_NUC,  MIN_ALH, DCS, HMOIST_950, & 
       HSMOIST_500, SINST, MAX_EXP, MAX_CAPE, MIN_CAPE, DUST_INFAC, ORG_INFAC, BC_INFAC, SS_INFAC, &
        MAPL, RRTMG_IRRAD, RRTMG_SORAD, SCWST, MTIME, SWCIRRUS, MINCDNC, TMAXBASELQ, TMAXCFCORR, Immersion_param, &
        DT_MICRO, DT_AUX, UR_SCALE    
        

      real                            :: THLSRC_PERT, QTSRC_PERT
      real                            :: UWTOLS
      real                            :: PMIN_CBL

      real                            :: CBL_TPERTi, CBL_TPERT, CBL_QPERT, RASAL1, RASAL2
      real                            :: RHEXCESS,CBL_TPERT_MXOCN,CBL_TPERT_MXLND
      real                            :: HEARTBEAT

      integer                         :: L, K, I, J, KM, KK , Kii, NAER, N !DONIF
     

      integer                         :: IDIM, IRUN, ICMIN, K0
      integer                         :: ITRCR,KSTRAP,CBL_METHOD,KCBLMIN,CLEANUP_RH

      real,    dimension(IM,JM,0:LM)  :: PKE
      real,    dimension(IM,JM,  LM)  :: DQS, QSS, PLO, ZLO, TEMP, PK, DP, DQSDT, DBZ3D
      real,    dimension(IM,JM,  LM)  :: KEX, DKEX
      real,    dimension(IM,JM,  LM)  :: Q1, W1, U1, V1, TH1, CNV_PRC3,fQi,CFPBL,CNV_HAIL

      integer                         :: SHLWDIAG
      real,    dimension(IM,JM,  LM)  :: SHLW_PRC3,SHLW_SNO3,UFRC_SC
      real,    dimension(IM,JM,0:LM)  :: CNV_PLE,ZLE
      real,    dimension(      0:LM)  :: SIGE
      real,    dimension(IM,JM,  LM)  :: GZLO, HHO,HSO
      real,    dimension(IM,JM,0:LM)  :: GZLE
      real,    dimension(IM,JM)       :: RASPRCP
      real,    dimension(IM,JM)       :: CO_AUTO, CCNSCALE
      integer, dimension(IM,JM,2)     :: SEEDRAS
      real,    dimension(IM,JM)       :: SEEDCNV
      integer, dimension(IM,JM)       :: KLCL, KLFC, KPBL, KPBL_SC

      real,    dimension(IM,JM)       :: LS_ARFX,CN_ARFX,AN_ARFX,SC_ARFX,QSSFC,IKEX,IKEX2
      real,    dimension(IM,JM)       :: LS_SNR, CN_SNR, AN_SNR, SC_SNR, ZCBLx, FILLQ, MXDIAMx
      real,    dimension(IM,JM)       :: LS_PRC2,CN_PRC2,AN_PRC2,SC_PRC2,ER_PRC2, TPREC, TVQX
      real,    dimension(IM,JM)       :: TPERT, QPERT, TSFCAIR, DTS, RASAL2_2d, NPRE_FRAC_2d
      integer, dimension(IM,JM)       :: IRAS, JRAS, KCBL
      real,    dimension(IM,JM,LM)    :: WGT0, WGT1
      real,    dimension(IM,JM,LM)    :: TRDLX
      integer, dimension(IM,JM,LM)    :: irccode
      real,    dimension(IM,JM)       :: ZWS
      REAL                            :: aux1,aux2,aux3,hfs,hfl,IFRC,TEND

      type (T_CLOUD_CTL)              :: CLOUD_CTL
      type (T_DDF_CTL)                :: DDF_CTL

      logical                         :: isPresent

      real, dimension(:,:,:,:,:), allocatable :: buffer
      logical :: USE_MOIST_BUFFER

      !!real,    dimension(IM,JM,  LM)  :: QILS, QICN ! Soon to be moved into internal state

      logical ALLOC_PTYPE

      logical ALLOC_CNV_DQLDT 
      logical ALLOC_CNV_MF0   
      logical ALLOC_CNV_MFD   
      logical ALLOC_CNV_MFC
      logical ALLOC_CNV_TOPP   
      logical ALLOC_CNV_UPDF  
      logical ALLOC_CNV_CVW  
      logical ALLOC_CNV_QC 
      logical ALLOC_RAD_CF    
      logical ALLOC_RAD_QV    
      logical ALLOC_RAD_QL    
      logical ALLOC_RAD_QI    
      logical ALLOC_RAD_QR    
      logical ALLOC_RAD_QS  
      logical ALLOC_RAD_QG  
      logical ALLOC_CLDREFFL  
      logical ALLOC_CLDREFFI  
      logical ALLOC_CLDREFFR  
      logical ALLOC_CLDREFFS  
      logical ALLOC_CLDREFFG  
      logical ALLOC_CLDNCCN   
      logical ALLOC_BYNCY
      logical ALLOC_CAPE 
      logical ALLOC_INHB 

      logical ALLOC_DQRC

      logical ALLOC_ENTLAM

      logical ALLOC_SMAXL !DONIF
      logical ALLOC_SMAXI
      logical ALLOC_WSUB
      logical ALLOC_CCN01
      logical ALLOC_CCN04
      logical ALLOC_CCN1
      logical ALLOC_CDNC_NUC
      logical ALLOC_INC_NUC
      logical ALLOC_NCPL_VOL
      logical ALLOC_NCPI_VOL
      logical ALLOC_SO4
      logical ALLOC_ORG
      logical ALLOC_DUST
      logical ALLOC_SEASALT
      logical ALLOC_BCARBON
      logical ALLOC_NHET_NUC
      logical ALLOC_NLIM_NUC
      logical ALLOC_SAT_RAT
      logical ALLOC_QSTOT
      logical ALLOC_QRTOT
      logical ALLOC_QPTOTLS

      logical ALLOC_PRCP_RAIN
      logical ALLOC_PRCP_SNOW
      logical ALLOC_PRCP_ICE
      logical ALLOC_PRCP_GRAUPEL
     
      logical ALLOC_DQVDT_micro      
      logical ALLOC_DQIDT_micro
      logical ALLOC_DQLDT_micro
      logical ALLOC_DQRDT_micro
      logical ALLOC_DQSDT_micro
      logical ALLOC_DQGDT_micro
      logical ALLOC_DQADT_micro
      logical ALLOC_DUDT_micro
      logical ALLOC_DVDT_micro
      logical ALLOC_DTDT_micro 
      logical ALLOC_DTDT_macro   
      logical ALLOC_DQVDT_macro
      logical ALLOC_DQLDT_macro
      logical ALLOC_DQIDT_macro
      logical ALLOC_DQADT_macro
      logical ALLOC_DTDT_moist  
      logical ALLOC_DTDTCN 
      logical ALLOC_RL_MASK  
      logical ALLOC_RI_MASK  
      logical ALLOC_KAPPA  
      logical ALLOC_SC_ICE  
      logical ALLOC_CFICE  
      logical ALLOC_CFLIQ  
      logical ALLOC_RHICE
      logical ALLOC_RHLIQ
      logical ALLOC_NHET_IMM
      logical ALLOC_DNHET_IMM
      logical ALLOC_NHET_DEP
      logical ALLOC_DUST_IMM
      logical ALLOC_DUST_DEP
      logical ALLOC_SCF
      logical ALLOC_SCF_ALL
      logical ALLOC_SIGW_GW
      logical ALLOC_SIGW_CNV
      logical ALLOC_SIGW_TURB
      logical ALLOC_SIGW_RC
      logical ALLOC_CNV_FICE       
      logical ALLOC_CNV_NDROP
      logical ALLOC_CNV_NICE
      logical ALLOC_RHCmicro
      logical ALLOC_DZover
      logical ALLOC_BERG
      logical ALLOC_BERGS
      logical ALLOC_MELT
      logical ALLOC_DNHET_CT
      logical ALLOC_PFRZ

      logical ALLOC_QCRES
      logical ALLOC_QIRES
      logical ALLOC_AUTICE

      logical ALLOC_DT_RASP
      logical ALLOC_FRZPP_LS
      logical ALLOC_SNOWMELT_LS


      logical :: ALLOC_DNCNUC, ALLOC_DNCHMSPLIT, ALLOC_DNCSUBL, &
           ALLOC_DNCAUTICE, ALLOC_DNCACRIS, ALLOC_DNDCCN, ALLOC_DNDACRLS, &
           ALLOC_DNDEVAPC, ALLOC_DNDACRLR,ALLOC_DNDAUTLIQ, ALLOC_DNDCNV, ALLOC_DNCCNV

      !------ logicals for allocating shallow vars -------

      logical :: ALLOC_UMF
      logical :: ALLOC_MFD
      logical :: ALLOC_DQVDT
      logical :: ALLOC_DQLDT
      logical :: ALLOC_DQIDT
      logical :: ALLOC_DTHDT
      logical :: ALLOC_DUDT
      logical :: ALLOC_DVDT
      logical :: ALLOC_DQRDT
      logical :: ALLOC_DQSDT
      logical :: ALLOC_CUFRC
      logical :: ALLOC_QCU
      logical :: ALLOC_QLU
      logical :: ALLOC_QIU
      logical :: ALLOC_CBMF
      logical :: ALLOC_DQCDT
      logical :: ALLOC_CNT
      logical :: ALLOC_CNB
      logical :: ALLOC_CIN
      logical :: ALLOC_PLCL
      logical :: ALLOC_PLFC
      logical :: ALLOC_PINV
      logical :: ALLOC_PREL
      logical :: ALLOC_PBUP
      logical :: ALLOC_WLCL
      logical :: ALLOC_QTSRC
      logical :: ALLOC_THLSRC
      logical :: ALLOC_THVLSRC
      logical :: ALLOC_TKEAVG
      logical :: ALLOC_CLDTOP
      logical :: ALLOC_WUP
      logical :: ALLOC_QTUP
      logical :: ALLOC_THLUP
      logical :: ALLOC_THVUP
      logical :: ALLOC_UUP
      logical :: ALLOC_VUP
      logical :: ALLOC_ENTR
      logical :: ALLOC_DETR
      logical :: ALLOC_XC
      logical :: ALLOC_QLDET
      logical :: ALLOC_QIDET
      logical :: ALLOC_QLENT
      logical :: ALLOC_QIENT
      logical :: ALLOC_QLSUB
      logical :: ALLOC_QISUB
      logical :: ALLOC_NDROP
      logical :: ALLOC_NICE

      !---------------------------------------------------

      character(len=ESMF_MAXSTR) :: GRIDNAME
      character(len=4)           :: imchar
      character(len=2)           :: dateline
      integer                    :: imsize,nn

      ! Needed for PTYPE diagnostic
      integer                    :: KTOP
      real                       :: PA2, PA, NA, TH_TOP, TH_BOT, TL_MEAN, Z_LAYER, ZTHICK

      integer                    :: levs925, levs600
      real, dimension(IM,JM   )  :: tempor2d, GF_AREA

      ! Manage diagnostic outputs for re-evaporation
      !---------------------------------------------------
      real, dimension(IM,JM,LM) :: RSU_CN_X
      real, dimension(IM,JM,LM) :: RSU_AN_X
      real, dimension(IM,JM,LM) :: RSU_LS_X
      real, dimension(IM,JM,LM) :: RSU_SC_X
      real, dimension(IM,JM,LM) :: RSU_MC_X

      real, dimension(IM,JM,LM) :: REV_CN_X
      real, dimension(IM,JM,LM) :: REV_AN_X
      real, dimension(IM,JM,LM) :: REV_LS_X
      real, dimension(IM,JM,LM) :: REV_SC_X
      real, dimension(IM,JM,LM) :: REV_MC_X

      ! Manage diagnostic outputs for accretion 
      !---------------------------------------------------
      real, dimension(IM,JM,LM) :: ACIL_CN_X
      real, dimension(IM,JM,LM) :: ACIL_AN_X
      real, dimension(IM,JM,LM) :: ACIL_LS_X
      real, dimension(IM,JM,LM) :: ACIL_SC_X

      real, dimension(IM,JM,LM) :: ACLL_CN_X
      real, dimension(IM,JM,LM) :: ACLL_AN_X
      real, dimension(IM,JM,LM) :: ACLL_LS_X
      real, dimension(IM,JM,LM) :: ACLL_SC_X

      ! Manage diagnostic outputs for 3D precip fluxes
      !---------------------------------------------------
      real, dimension(IM,JM,0:LM) :: PFI_CN_X
      real, dimension(IM,JM,0:LM) :: PFI_AN_X
      real, dimension(IM,JM,0:LM) :: PFI_LS_X
      real, dimension(IM,JM,0:LM) :: PFI_SC_X

      real, dimension(IM,JM,0:LM) :: PFL_CN_X
      real, dimension(IM,JM,0:LM) :: PFL_AN_X
      real, dimension(IM,JM,0:LM) :: PFL_LS_X
      real, dimension(IM,JM,0:LM) :: PFL_SC_X

      real, dimension(IM,JM,LM) :: DLPDF_X
      real, dimension(IM,JM,LM) :: DIPDF_X
      real, dimension(IM,JM,LM) :: DLFIX_X
      real, dimension(IM,JM,LM) :: DIFIX_X

      real, dimension(IM,JM,LM) :: RHX_X
      real, dimension(IM,JM,LM) :: AUT_X
      real, dimension(IM,JM,LM) :: EVAPC_X
      real, dimension(IM,JM,LM) :: SDM_X
      real, dimension(IM,JM,LM) :: SUBLC_X
      real, dimension(IM,JM,LM) :: FRZ_TT_X
      real, dimension(IM,JM,LM) :: FRZ_PP_X
      real, dimension(IM,JM,LM) :: DCNVL_X
      real, dimension(IM,JM,LM) :: DCNVI_X

      real, dimension(IM,JM,LM) :: DNDCNV_X
      real, dimension(IM,JM,LM) :: DNCCNV_X

      real, dimension(IM,JM,LM) :: LWC_X, IWC_X, CIRRUS_FREQ_HOM_X, CIRRUS_FREQ_HET_X
      real, dimension(IM,JM) ::  CLDREFFI_TOP_X, CLDREFFL_TOP_X,  NCPL_TOP_X, NCPI_TOP_X, NCPL_CLDBASEX


      real, dimension(IM,JM,LM) :: ALPHT_X
      real, dimension(IM,JM,LM) :: ALPH1_X
      real, dimension(IM,JM,LM) :: ALPH2_X

      real, dimension(IM,JM,LM) :: CFPDF_X
      real, dimension(IM,JM,LM) :: RHCLR_X
      real, dimension(IM,JM,LM) :: DQRL_X

      real, dimension(IM,JM,LM) :: VFALLICE_AN_X
      real, dimension(IM,JM,LM) :: VFALLICE_LS_X
      real, dimension(IM,JM,LM) :: VFALLWAT_AN_X
      real, dimension(IM,JM,LM) :: VFALLWAT_LS_X

      real, dimension(IM,JM,LM) :: VFALLRN_AN_X
      real, dimension(IM,JM,LM) :: VFALLRN_LS_X
      real, dimension(IM,JM,LM) :: VFALLRN_CN_X
      real, dimension(IM,JM,LM) :: VFALLRN_SC_X
      real, dimension(IM,JM,LM) :: VFALLSN_AN_X
      real, dimension(IM,JM,LM) :: VFALLSN_LS_X
      real, dimension(IM,JM,LM) :: VFALLSN_CN_X
      real, dimension(IM,JM,LM) :: VFALLSN_SC_X

!!$      real, dimension(IM,JM,LM) :: LIQANMOVE_X 
!!$      real, dimension(IM,JM,LM) :: ICEANMOVE_X 
!!$      real, dimension(IM,JM,LM) :: DANCLD_X 
!!$      real, dimension(IM,JM,LM) :: DLSCLD_X 
!!$      real, dimension(IM,JM,LM) :: CURAINMOVE_X 
!!$      real, dimension(IM,JM,LM) :: CUSNOWMOVE_X 

      ! MATMAT Additional variables for GPUization
      real, dimension(IM,JM,LM)   :: DQST3,QST3,TEMP_0
      real, dimension(IM,JM,LM)   :: DZ, DZET, iMASS, MASS, QDDF3
      real, dimension(IM,JM)      :: VMIP
      real, dimension(IM,JM,LM+1) :: ZET

      ! Convective Fraction
      integer, dimension(IM,JM)           :: L600
      real   , dimension(IM,JM)           :: RHat600, QV600, QVCBL
      real   , dimension(IM,JM)           :: CNV_FRACTION
      real                                :: CNV_FRACTION_MIN
      real                                :: CNV_FRACTION_MAX
      real                                :: CNV_FRACTION_EXP
      real                                :: GF_MIN_AREA
      ! Control for stochasticity in CNV
      integer                             :: STOCHASTIC_CNV

      real :: icefall
      real :: cNN, cNN_OCEAN, cNN_LAND, CONVERT

      real   , dimension(IM,JM)           :: CMDU, CMSS, CMOC, CMBC, CMSU, CMNI, CMNH3, CMNH4A, CMBRC
      real   , dimension(IM,JM)           :: CMDU2g, CMSS2g, CMOC2g, CMBC2g, CMSU2g, CMNI2g, CMNH32g, CMNH4A2g, CMBRC2g
      real   , dimension(IM,JM)           :: CMDUcarma, CMSScarma
       
      real :: sigmaqt, qcn, cfn, qsatn, dqlls, dqils, qt

      ! MATMAT CUDA Variables
#ifdef _CUDA
      type(dim3) :: Grid, Block
      integer :: blocksize

      logical :: COPY_RSU_CN
      logical :: COPY_RSU_AN
      logical :: COPY_RSU_LS
      logical :: COPY_RSU_SC

      logical :: COPY_REV_CN
      logical :: COPY_REV_AN
      logical :: COPY_REV_LS
      logical :: COPY_REV_SC

      logical :: COPY_ACIL_CN
      logical :: COPY_ACIL_AN
      logical :: COPY_ACIL_LS
      logical :: COPY_ACIL_SC

      logical :: COPY_ACLL_CN
      logical :: COPY_ACLL_AN
      logical :: COPY_ACLL_LS
      logical :: COPY_ACLL_SC

      logical :: COPY_DLPDF
      logical :: COPY_DIPDF
      logical :: COPY_DLFIX
      logical :: COPY_DIFIX

      logical :: COPY_PFI_CN
      logical :: COPY_PFI_AN
      logical :: COPY_PFI_LS
      logical :: COPY_PFI_SC

      logical :: COPY_PFL_CN
      logical :: COPY_PFL_AN
      logical :: COPY_PFL_LS
      logical :: COPY_PFL_SC

      logical :: COPY_RHX
      logical :: COPY_AUT     
      logical :: COPY_EVAPC   
      logical :: COPY_SDM    
      logical :: COPY_SUBLC    
      logical :: COPY_FRZ_TT 
      logical :: COPY_FRZ_PP   
      logical :: COPY_DCNVL    
      logical :: COPY_DCNVI   

      logical :: COPY_ALPHT
      logical :: COPY_ALPH1
      logical :: COPY_ALPH2

      logical :: COPY_CFPDF
      logical :: COPY_RHCLR
      logical :: COPY_DQRL    

      logical :: COPY_VFALLICE_AN
      logical :: COPY_VFALLICE_LS
      logical :: COPY_VFALLWAT_AN
      logical :: COPY_VFALLWAT_LS

      logical :: COPY_VFALLRN_AN
      logical :: COPY_VFALLRN_LS
      logical :: COPY_VFALLRN_CN
      logical :: COPY_VFALLRN_SC
      logical :: COPY_VFALLSN_AN
      logical :: COPY_VFALLSN_LS
      logical :: COPY_VFALLSN_CN
      logical :: COPY_VFALLSN_SC

!!$      logical :: COPY_LIQANMOVE
!!$      logical :: COPY_ICEANMOVE
!!$      logical :: COPY_DANCLD
!!$      logical :: COPY_DLSCLD
!!$      logical :: COPY_CURAINMOVE
!!$      logical :: COPY_CUSNOWMOVE

#endif

      !  Begin...
      !----------
      Iam = trim(COMP_NAME) // 'Convect_Driver'

      call MAPL_TimerOn(STATE,"-MISC1")

      ! Get parameters from configuration
      !----------------------------------
      call ESMF_ConfigGetAttribute (CF, HEARTBEAT, Label="RUN_DT:", RC=STATUS)
      VERIFY_(STATUS)

      call ESMF_ConfigGetAttribute( CF, RAS_NO_NEG, Label='RAS_NO_NEG:', default=.FALSE. , RC=STATUS)
      VERIFY_(STATUS)

      call MAPL_GetResource(STATE, CLEANUP_RH,               'CLEANUP_RH:',     DEFAULT= 1,     RC=STATUS)
      call MAPL_GetResource(STATE, RASPARAMS%CUFRICFAC,      'CUFRICFAC:',      DEFAULT= 1.000, RC=STATUS)
      call MAPL_GetResource(STATE, RASPARAMS%SHR_LAMBDA_FAC, 'SHR_LAMBDA_FAC:', DEFAULT= 0.05,  RC=STATUS)

      call MAPL_GetResource(STATE, AUTOC_CN_OCN, 'AUTOC_CN:',      DEFAULT= 2.5e-3,RC=STATUS)


! gf convective param
    call MAPL_GetResource(STATE, CONVPAR_OPTION, Label='CONVPAR_OPTION:', RC=STATUS)
    VERIFY_(STATUS)
!


      !======================================================================
      !2-M microphysics "tuning nobs"
      
      call MAPL_GetResource(STATE, AUTO_CNV,       'AUTO_CNV:',       DEFAULT= 1.0,    RC=STATUS) !scaling factor for ras liquid precip critical size
      call MAPL_GetResource(STATE, LCCIRRUS,       'LCCIRRUS:',       DEFAULT= 500.0,  RC=STATUS) !Characteristic Length (m) of high freq gravity waves
      call MAPL_GetResource(STATE, UISCALE,        'UISCALE:',        DEFAULT= 1.0,    RC=STATUS) !Scaling factor for sed vel of ice      
      call MAPL_GetResource(STATE, LIU_MU,         'LIU_MU:',         DEFAULT= 2.0,    RC=STATUS) !Liu autoconversion parameter
      call MAPL_GetResource(STATE, NPRE_FRAC,      'NPRE_FRAC:',      DEFAULT= -1.0,   RC=STATUS) !Fraction of preexisting ice affecting ice nucleationn            
      call MAPL_GetResource(STATE, CLDPARAMS%MIN_LTS,'LTS_LOW:',      DEFAULT= 20.0,   RC=STATUS) !lower LTS for morphology correction
      call MAPL_GetResource(STATE, LTS_UP,         'LTS_UP:',         DEFAULT= 22.0,   RC=STATUS) !Upper LTS for morphology correction
      call MAPL_GetResource(STATE, MIN_EXP,        'MIN_EXP:',        DEFAULT= 0.5,   RC=STATUS) !Exponent of the relation CFA=CFV^n
       call MAPL_GetResource(STATE, MAX_EXP,        'MAX_EXP:',       DEFAULT= 1.0,   RC=STATUS) !Exponent of the relation CFA=CFV^n
      call MAPL_GetResource(STATE, DCRIT_,         'DCRIT_DRIZZLE:',  DEFAULT= 0.6,    RC=STATUS) !scale factor for critical size for drizzle
      call MAPL_GetResource(STATE, USE_AV_V,       'USE_AV_V:',       DEFAULT= 1.0,    RC=STATUS) !Set to > 0 to use an average velocity for activation
      call MAPL_GetResource(STATE, TS_AUTO_ICE,    'TS_AUTO_ICE:',    DEFAULT= 4.0, RC=STATUS) !Ice autoconversion time scale
      call MAPL_GetResource(STATE, TMAXLL,         'TMAXLL:',         DEFAULT= 250.0,  RC=STATUS) !Liquid clouds min T
      call MAPL_GetResource(STATE, CCN_PARAM,      'CCNPARAM:',       DEFAULT= 2.0,    RC=STATUS) !CCN activation param
      call MAPL_GetResource(STATE, IN_PARAM,       'INPARAM:',        DEFAULT= 6.0,    RC=STATUS) !IN param
      call MAPL_GetResource(STATE, Immersion_param,'ImmersionPARAM:', DEFAULT= 6.0,    RC=STATUS) !Immersion param
      
      call MAPL_GetResource(STATE, FDROP_DUST,     'FDROP_DUST:',     DEFAULT= 0.04,    RC=STATUS) !Fraction of dust within droplets for immersion freezing
      call MAPL_GetResource(STATE, FDROP_SOOT,     'FDROP_SOOT:',     DEFAULT= 0.01,   RC=STATUS) !Fraction of soot within droplets for immersion freezing	
      call MAPL_GetResource(STATE, SIGMA_NUC,      'SIGMA_NUC:',      DEFAULT= 1.0,   RC=STATUS) !Widht of the in-cloud distribution of relative humidity in cirrus
      call MAPL_GetResource(STATE, MIN_ALH,            'MIN_ALH:',      DEFAULT= 5.0,  RC=STATUS) !scale factor for vertical velocity in sttratocumulus
      call MAPL_GetResource(STATE, SCWST,            'SCWST:',      DEFAULT= 5.0,  RC=STATUS) !scale factor for vertical velocity in sttratocumulus
     
      call MAPL_GetResource(STATE, MINCDNC,          'MINCDNC:',      DEFAULT= 0.0,  RC=STATUS) !min nucleated droplet conc. cm-3

      call MAPL_GetResource(STATE, TMAXCFCORR,         'TMAXCFCORR:',     DEFAULT= 285.0,  RC=STATUS) !Minimum T for CF correction
      call MAPL_GetResource(STATE, MTIME,         'MTIME:',  DEFAULT= -1.0,    RC=STATUS) !Mixing time scale for aerosol within the cloud. Default is time step
      call MAPL_GetResource(STATE, SWCIRRUS, 'SWCIRRUS:', DEFAULT= 1.0, RC=STATUS) !Tunes vertical velocity in cirrus
      
      call MAPL_GetResource(STATE, DUST_INFAC,    'DUST_INFAC:',        DEFAULT= 0.5,   RC=STATUS)  !work on this
      call MAPL_GetResource(STATE, BC_INFAC,        'BC_INFAC:',        DEFAULT= 0.1,   RC=STATUS) 
      call MAPL_GetResource(STATE, ORG_INFAC,     'ORG_INFAC:',        DEFAULT= 1.0,   RC=STATUS)   
	 call MAPL_GetResource(STATE, SS_INFAC,          'SS_INFAC:',        DEFAULT= 1.0,   RC=STATUS)   
     	  
      call MAPL_GetResource(STATE, DT_MICRO,          'DT_MICRO:',        DEFAULT= HEARTBEAT,   RC=STATUS)    ! time step of the microphysics substepping (s) (MG2) (5 min)
      call MAPL_GetResource(STATE, UR_SCALE,        'URSCALE:',        DEFAULT= 1.0,    RC=STATUS) !Scaling factor for sed vel of rain    
          
      call MAPL_GetResource(STATE, USE_NATURE_WSUB,     'USE_NAT_WSUB:',     DEFAULT= 1.0  ,RC=STATUS) !greater than zero reads wsub from nature run	             
      call MAPL_GetResource(STATE, DCS, 'DCS:', default=350.0e-6, RC=STATUS )
      call MAPL_GetResource(STATE, CLDPARAMS%DISP_FACTOR_LIQ,         'DISP_FACTOR_LIQ:',     DEFAULT= 20.0,   RC=STATUS) ! Scales the droplet/ice crystal number in convective detrainment 
      call MAPL_GetResource(STATE, CLDPARAMS%DISP_FACTOR_ICE,         'DISP_FACTOR_ICE:',     DEFAULT= 10.0,   RC=STATUS) ! Scales the droplet/ice crystal number in convective detrainment 
      
      call MAPL_GetResource( STATE, RRTMG_IRRAD ,'USE_RRTMG_IRRAD:', DEFAULT=0.0, RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetResource( STATE, RRTMG_SORAD ,'USE_RRTMG_SORAD:', DEFAULT=0.0, RC=STATUS)
      VERIFY_(STATUS)
    
      
!!!!!!!!!!!!!!!!!
      call MAPL_GetResource(STATE,GRIDNAME,'AGCM_GRIDNAME:', RC=STATUS)
      VERIFY_(STATUS)
      GRIDNAME =  AdjustL(GRIDNAME)
      nn = len_trim(GRIDNAME)
      dateline = GRIDNAME(nn-1:nn)
      imchar = GRIDNAME(3:index(GRIDNAME,'x')-1)
      read(imchar,*) imsize
      if(dateline.eq.'CF') imsize = imsize*4

      call MAPL_GetResource(STATE, RASPARAMS%QC_CRIT_CN,   'QC_CRIT_CN:',    DEFAULT= 8.0e-4,  RC=STATUS)

      call MAPL_GetResource(STATE, RASPARAMS%RASAL1,       'RASAL1:',        DEFAULT=   1800.0,RC=STATUS)
      call MAPL_GetResource(STATE, RASPARAMS%RASAL2,       'RASAL2:',        DEFAULT= -86400.0,RC=STATUS)

      call MAPL_GetResource(STATE, RASPARAMS%RASNCL,       'RASNCL:',        DEFAULT= -300.,   RC=STATUS)
      call MAPL_GetResource(STATE, RASPARAMS%LAMBDA_FAC,   'LAMBDA_FAC:',    DEFAULT= 4.0 ,    RC=STATUS)
      call MAPL_GetResource(STATE, RASPARAMS%LAMBMX_FAC,   'LAMBMX_FAC:',    DEFAULT= 0.0 ,    RC=STATUS)
      call MAPL_GetResource(STATE, RASPARAMS%MIN_DIAMETER, 'MIN_DIAMETER:',  DEFAULT= 200.,    RC=STATUS)
      call MAPL_GetResource(STATE, RASPARAMS%CUFRICLAMBDA, 'CUFRICLAMBDA:',  DEFAULT= 7.5e-4,  RC=STATUS)
      call MAPL_GetResource(STATE, RASPARAMS%RDTLEXPON,    'RDTLEXPON:',     DEFAULT= 1.0 ,    RC=STATUS)
      call MAPL_GetResource(STATE, RASPARAMS%STRAPPING,    'STRAPPING:',     DEFAULT=-1.0 ,    RC=STATUS)
      call MAPL_GetResource(STATE, RASPARAMS%SDQV2,        'SDQV2:',         DEFAULT= 1.3 ,    RC=STATUS)
      call MAPL_GetResource(STATE, RASPARAMS%SDQV3,        'SDQV3:',         DEFAULT= 1.3 ,    RC=STATUS)
      call MAPL_GetResource(STATE, RASPARAMS%SDQVT1,       'SDQVT1:',        DEFAULT= 263.,    RC=STATUS)
      call MAPL_GetResource(STATE, RASPARAMS%ACRITFAC,     'ACRITFAC:',      DEFAULT= 0.5 ,    RC=STATUS)
      call MAPL_GetResource(STATE, RASPARAMS%HMINTRIGGER,  'HMINTRIGGER:',   DEFAULT= 1.0 ,    RC=STATUS)
      call MAPL_GetResource(STATE, RASPARAMS%LLDISAGGXP,   'LLDISAGGXP:',    DEFAULT= 0.0 ,    RC=STATUS)
      call MAPL_GetResource(STATE, RASPARAMS%PBLFRAC,      'PBLFRAC:',       DEFAULT= 0.1 ,    RC=STATUS)
      call MAPL_GetResource(STATE, RASPARAMS%RASAUTORAMPB, 'RASAUTORAMPB:',  DEFAULT= 0.8 ,    RC=STATUS)
      call MAPL_GetResource(STATE, PMIN_DET,               'PMIN_DET',       DEFAULT= 3000.0,  RC=STATUS)
      call MAPL_GetResource(STATE, PMIN_CBL,               'PMIN_CBL',       DEFAULT= 50000.0, RC=STATUS)

      call MAPL_GetResource(STATE,           AUTOC_CN_LAND,'AUTOC_CN_LAND:', DEFAULT= AUTOC_CN_OCN, RC=STATUS)
      call MAPL_GetResource(STATE, RASPARAMS%AUTOC_CN_ZDEP,'AUTOC_CN_ZDEP:', DEFAULT= 1.0          ,RC=STATUS)


      if( imsize.le.200       ) call MAPL_GetResource(STATE, RASPARAMS%MAXDALLOWED_S, 'MAXDALLOWED_S:', DEFAULT= 4000.0 ,RC=STATUS)
      if( imsize.gt.200 .and. &
          imsize.le.400       ) call MAPL_GetResource(STATE, RASPARAMS%MAXDALLOWED_S, 'MAXDALLOWED_S:', DEFAULT= 2000.0 ,RC=STATUS)
      if( imsize.gt.400 .and. &
          imsize.le.800       ) call MAPL_GetResource(STATE, RASPARAMS%MAXDALLOWED_S, 'MAXDALLOWED_S:', DEFAULT=  700.0 ,RC=STATUS)
      if( imsize.gt.800 .and. &
          imsize.le.1600      ) call MAPL_GetResource(STATE, RASPARAMS%MAXDALLOWED_S, 'MAXDALLOWED_S:', DEFAULT=  450.0 ,RC=STATUS)
      if( imsize.gt.1600      ) call MAPL_GetResource(STATE, RASPARAMS%MAXDALLOWED_S, 'MAXDALLOWED_S:', DEFAULT=  450.0 ,RC=STATUS)
                                call MAPL_GetResource(STATE, RASPARAMS%MAXDALLOWED_D, 'MAXDALLOWED_D:', DEFAULT= RASPARAMS%MAXDALLOWED_S ,RC=STATUS)
                                call MAPL_GetResource(STATE, RASPARAMS%MAXDALLOWED_E, 'MAXDALLOWED_E:', DEFAULT= -0.500 ,RC=STATUS)

      call MAPL_GetResource(STATE, RASPARAMS%RASAL_SLOPE , 'RASAL_SLOPE:',  DEFAULT= 8000.0  ,  RC=STATUS)

      call MAPL_GetResource(STATE, RASPARAMS%RAS_RHMIN , 'RAS_RHMIN:',  DEFAULT= 0.5  ,  RC=STATUS)
      call MAPL_GetResource(STATE, RASPARAMS%RAS_RHFULL, 'RAS_RHFULL:', DEFAULT= 0.65 ,  RC=STATUS)

      call MAPL_GetResource(STATE,CBL_METHOD,   'CBL_METHOD:',    DEFAULT= 6     , RC=STATUS)
      call MAPL_GetResource(STATE,CBL_QPERT,    'CBL_QPERT:',     DEFAULT= 0.0   , RC=STATUS)
      call MAPL_GetResource(STATE,CBL_TPERT,    'CBL_TPERT:',     DEFAULT=-1.0   , RC=STATUS)

      call MAPL_GetResource(STATE,CBL_TPERT_MXOCN, 'CBL_TPERT_MXOCN:',     DEFAULT= 2.0   , RC=STATUS)
      call MAPL_GetResource(STATE,CBL_TPERT_MXLND, 'CBL_TPERT_MXLND:',     DEFAULT= 0.0   , RC=STATUS)

     !!! MODIFIED by npa: remove when done testing shallow
      call MAPL_GetResource(STATE,THLSRC_PERT, 'THLSRC_PERT:',     DEFAULT= 0.0   , RC=STATUS)     
      call MAPL_GetResource(STATE,QTSRC_PERT, 'QTSRC_PERT:',     DEFAULT= 1.0   , RC=STATUS)     
      call MAPL_GetResource(STATE,UWTOLS, 'UWTOLS:',     DEFAULT= 0.0   , RC=STATUS)     

      KSTRAP = INT( RASPARAMS%STRAPPING )


     ! Get parameters for shallow convection
     !----------------------------------------------
      if (DOSHLW /= 0) then
         call MAPL_GetResource(STATE, RASPARAMS%MIN_DIAMETER, 'MIN_DIAMETER:',  DEFAULT= 400.,    RC=STATUS)
      else
         call MAPL_GetResource(STATE, RASPARAMS%MIN_DIAMETER, 'MIN_DIAMETER:',  DEFAULT= 200.,    RC=STATUS)
      end if

      call MAPL_GetResource(STATE, SHLWPARAMS%NITER_XC,         'NITER_XC:'       ,DEFAULT=2, RC=STATUS)
      call MAPL_GetResource(STATE, SHLWPARAMS%ITER_CIN,         'ITER_CIN:'       ,DEFAULT=2, RC=STATUS)
      call MAPL_GetResource(STATE, SHLWPARAMS%USE_CINCIN,       'USE_CINCIN:'     ,DEFAULT=1, RC=STATUS)
      call MAPL_GetResource(STATE, SHLWPARAMS%USE_SELF_DETRAIN, 'USE_SELF_DETRAIN:',DEFAULT=0, RC=STATUS)
      call MAPL_GetResource(STATE, SHLWPARAMS%USE_MOMENFLX,     'USE_MOMENFLX:'   ,DEFAULT=1, RC=STATUS)
      call MAPL_GetResource(STATE, SHLWPARAMS%USE_CUMPENENT,    'USE_CUMPENENT:'   ,DEFAULT=1, RC=STATUS)
      call MAPL_GetResource(STATE, SHLWPARAMS%SCVERBOSE,        'SCVERBOSE:',      DEFAULT=0, RC=STATUS)
      call MAPL_GetResource(STATE, SHLWPARAMS%RPEN,    'RPEN:'    ,DEFAULT=3.0, RC=STATUS)
      call MAPL_GetResource(STATE, SHLWPARAMS%RLE,     'RLE:'     ,DEFAULT=0.1, RC=STATUS)
      call MAPL_GetResource(STATE, SHLWPARAMS%RKM,     'RKM:'     ,DEFAULT=12.0, RC=STATUS)
      call MAPL_GetResource(STATE, SHLWPARAMS%RKFRE,   'RKFRE:'   ,DEFAULT=1.0, RC=STATUS)
      call MAPL_GetResource(STATE, SHLWPARAMS%RMAXFRAC,'RMAXFRAC:',DEFAULT=0.1,  RC=STATUS)
      call MAPL_GetResource(STATE, SHLWPARAMS%MUMIN1,  'MUMIN1:'  ,DEFAULT=0.906, RC=STATUS)
      call MAPL_GetResource(STATE, SHLWPARAMS%RBUOY,   'RBUOY:'   ,DEFAULT=1.0, RC=STATUS)
      call MAPL_GetResource(STATE, SHLWPARAMS%RDRAG,   'RDRAG:'   ,DEFAULT=1.0, RC=STATUS)
      call MAPL_GetResource(STATE, SHLWPARAMS%EPSVARW, 'EPSVARW:' ,DEFAULT=5.e-4, RC=STATUS)
      call MAPL_GetResource(STATE, SHLWPARAMS%PGFC,    'PGFC:'    ,DEFAULT=0.7, RC=STATUS)
      call MAPL_GetResource(STATE, SHLWPARAMS%CRIQC,   'CRIQC:'   ,DEFAULT=1.0e-3, RC=STATUS)
      call MAPL_GetResource(STATE, SHLWPARAMS%KEVP,    'KEVP:'    ,DEFAULT=2.e-6,    RC=STATUS)
      call MAPL_GetResource(STATE, SHLWPARAMS%RDROP,   'SHLW_RDROP:',DEFAULT=8.e-6,    RC=STATUS)

      if(adjustl(CLDMICRO)=="GFDL") then
        call MAPL_GetResource(STATE, DOCLDMACRO,         'DOCLDMACRO:' ,DEFAULT=0  , RC=STATUS)
        call MAPL_GetResource(STATE, SHLWPARAMS%FRC_RASN,'FRC_RASN:'   ,DEFAULT=1.0, RC=STATUS)
      else
        call MAPL_GetResource(STATE, DOCLDMACRO,         'DOCLDMACRO:' ,DEFAULT=1  , RC=STATUS)
        call MAPL_GetResource(STATE, SHLWPARAMS%FRC_RASN,'FRC_RASN:'   ,DEFAULT=0.0, RC=STATUS)
      endif

      ! Get the time step from the alarm
      !---------------------------------

      call ESMF_AlarmGet( ALARM, RingInterval=TINT, RC=STATUS)
      VERIFY_(STATUS)
      call ESMF_TimeIntervalGet(TINT, S_R8=DT_R8, RC=STATUS)
      VERIFY_(STATUS)

      DT_MOIST = DT_R8

      ! Pointers to internals
      !----------------------

      call MAPL_GetPointer(INTERNAL, Q,        'Q'       , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(INTERNAL, QRAIN,    'QRAIN'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(INTERNAL, QSNOW,    'QSNOW'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(INTERNAL, QGRAUPEL, 'QGRAUPEL', RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(INTERNAL, QLLS,     'QLLS'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(INTERNAL, QLCN,     'QLCN'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(INTERNAL, CLCN,     'CLCN'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(INTERNAL, CLLS,     'CLLS'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(INTERNAL, QILS,     'QILS'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(INTERNAL, QICN,     'QICN'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(INTERNAL, NCPL,     'NCPL'    , RC=STATUS); VERIFY_(STATUS)  !DONIF
      call MAPL_GetPointer(INTERNAL, NCPI,     'NCPI'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(INTERNAL, NRAIN,    'NRAIN'    , RC=STATUS); VERIFY_(STATUS)  
      call MAPL_GetPointer(INTERNAL, NSNOW,    'NSNOW'    , RC=STATUS); VERIFY_(STATUS)      
      call MAPL_GetPointer(INTERNAL, NGRAUPEL, 'NGRAUPEL'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(INTERNAL, PDF_A,    'PDF_A'   , RC=STATUS); VERIFY_(STATUS) 
       
      if (DOSHLW /= 0) then
       call MAPL_GetPointer(INTERNAL, CUSH,  'CUSH'    , RC=STATUS); VERIFY_(STATUS)  !DONIF
      end if
      ! Pointers to imports
      !--------------------
      !==AER_CLOUD===

      call MAPL_GetPointer(IMPORT, TAUGWX, 'TAUGWX'     , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT, TAUGWY, 'TAUGWY'     , RC=STATUS); VERIFY_(STATUS)     
      call MAPL_GetPointer(IMPORT, TAUX,   'TAUX'     , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT, TAUY,   'TAUY'     , RC=STATUS); VERIFY_(STATUS)            
      call MAPL_GetPointer(IMPORT, OMEGA,  'OMEGA'     , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT, ALH,    'ALH'     , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT, RADLW,  'RADLW'     , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT, RADSW,  'RADSW'     , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT, TAUOROX, 'TAUOROX'     , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT, TAUOROY, 'TAUOROY'     , RC=STATUS); VERIFY_(STATUS)  
      call MAPL_GetPointer(IMPORT, WSUB_NATURE,  'WSUB_NATURE'     , RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(IMPORT, AREA,    'AREA'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT, PLE,     'PLE'     , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT, PREF,    'PREF'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT, KH,      'KH'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT, TKE,     'TKE'     , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT, TH,      'TH'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT, U,       'U'       , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT, V,       'V'       , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT, W,       'W'       , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT, SNOMAS,  'SNOMAS'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT, FRLANDICE,  'FRLANDICE'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT, FRLAND,  'FRLAND'  , RC=STATUS); VERIFY_(STATUS)
      ! frland =0.0
      call MAPL_GetPointer(IMPORT, TS,      'TS'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT, TROPP,   'TROPP'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT, KPBLIN,  'KPBL'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT, KPBLSC,  'KPBL_SC' , RC=STATUS); VERIFY_(STATUS)

      call ESMF_StateGet (IMPORT, 'MTR' ,   TR,        RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(IMPORT, SH       ,'SH   '   ,RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT, EVAP     ,'EVAP '   ,RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT, T        ,'T'       ,RC=STATUS); VERIFY_(STATUS)

      ! Get pointers to inputs for ADG PDF
      call MAPL_GetPointer(IMPORT, EDMF_FRC,'EDMF_FRC',RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT, W2 ,    'W2' ,    RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT, W3 ,    'W3' ,    RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT, WQT ,   'WQT',    RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT, WHL ,   'WHL',    RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT, HL2 ,   'HL2',    RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT, HL3 ,   'HL3',    RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT, QT2  ,  'QT2',    RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT, QT3  ,  'QT3',    RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT, HLQT,   'HLQT',   RC=STATUS); VERIFY_(STATUS)


      ! Pointers to exports
      !--------------------

      call MAPL_GetPointer(EXPORT, TI,       'DTHDT'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, UI,       'DUDT'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, VI,       'DVDT'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, WI,       'DWDT'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DQDT,     'DQDT'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, PTYPE,    'PTYPE'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, FRZR,     'FRZR'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, ICE,      'ICE'     , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, SNR,      'SNO'     , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, PRELS,    'PLS'     , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, PRECU,    'PCU'     , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, RH1   ,   'RH1'     , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, RH2   ,   'RH2'     , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, RHX   ,   'RHX'     , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, REV_CN,   'REV_CN'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, REV_SC,   'REV_SC'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, REV_AN,   'REV_AN'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, REV_LS,   'REV_LS'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, RSU_CN,   'RSU_CN'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, RSU_SC,   'RSU_SC'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, RSU_AN,   'RSU_AN'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, RSU_LS,   'RSU_LS'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DQRL,     'DQRL'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DQRC,     'DQRC'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DQDTCN,   'DQDTCN'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DQLDT,    'DQLDT'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DQIDT,    'DQIDT'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DQCDTCN,  'DQCDTCN' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DTHDTCN,  'DTHDTCN' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CCWP,     'CCWP'    , RC=STATUS); VERIFY_(STATUS) 
      call MAPL_GetPointer(EXPORT, CWP,      'CWP'     , RC=STATUS); VERIFY_(STATUS)  
      call MAPL_GetPointer(EXPORT, IWP,      'IWP'     , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, LWP,      'LWP'     , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, TPW,      'TPW'     , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, BYNCY,    'BYNCY'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CAPE,     'CAPE'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, INHB,     'INHB'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, TVQ0,     'TVQ0'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, TVQ1,     'TVQ1'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DCPTE,    'DCPTE'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, TVE0,     'TVE0'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, TVE1,     'TVE1'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, TVEX,     'TVEX'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, ZPBLCN,   'ZPBLCN'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, ZLCL,     'ZLCL'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, ZLFC,     'ZLFC'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, ZCBL,     'ZCBL'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, TT_PRCP,  'TPREC'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, HOURNORAIN,  'HOURNORAIN'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, TX    ,   'T'       , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, SLX   ,   'SL'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, QTX   ,   'QT'      , RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(EXPORT, PRCP_RAIN,    'PRCP_RAIN' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, PRCP_SNOW,    'PRCP_SNOW'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, PRCP_ICE,     'PRCP_ICE'     , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, PRCP_GRAUPEL, 'PRCP_GRAUPEL' , RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(EXPORT, LS_PRCP,  'LS_PRCP' , RC=STATUS); VERIFY_(STATUS) 
      call MAPL_GetPointer(EXPORT, CN_PRCP,  'CN_PRCP' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, AN_PRCP,  'AN_PRCP' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, ER_PRCP,  'ER_PRCP' , RC=STATUS); VERIFY_(STATUS) 
      call MAPL_GetPointer(EXPORT, FILLNQV,  'FILLNQV' , RC=STATUS); VERIFY_(STATUS) 
      call MAPL_GetPointer(EXPORT, XQLLS,    'QLLSX1'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, XQLCN,    'QLCNX1'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, XQILS,    'QILSX1'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, XQICN,    'QICNX1'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, XCLCN,    'CLCN'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, XCLLS,    'CLLS'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, QX1,      'QX1'     , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, QCTOT,    'QCTOT'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, QITOT,    'QITOT'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, QRTOT,    'QRTOT'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, QSTOT,    'QSTOT'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, QPTOTLS,  'QPTOTLS' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, QLTOT,    'QLTOT'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, LS_ARF,   'LS_ARF'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CN_ARF,   'CN_ARF'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, SC_ARF,   'SC_ARF'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, AN_ARF,   'AN_ARF'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CNV_DQLDT,'CNV_DQLDT',RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CNV_MF0,  'CNV_MF0' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CNV_MFD,  'CNV_MFD' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CNV_MFC,  'CNV_MFC' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CNV_FREQ, 'CNV_FREQ', RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CNV_BASEP,'CNV_BASEP', RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CNV_TOPP, 'CNV_TOPP', RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CNV_UPDF, 'CNV_UPDF', RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CNV_CVW,  'CNV_CVW' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CNV_QC,   'CNV_QC'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, RAD_CF,   'FCLD'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, RAD_QV,   'QV'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, RAD_QL,   'QL'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, RAD_QI,   'QI'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, RAD_QR,   'QR'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, RAD_QS,   'QS'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, RAD_QG,   'QG'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DBZ,      'DBZ'     , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DBZ_MAX,  'DBZ_MAX' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CLDREFFL, 'RL'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CLDREFFI, 'RI'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CLDREFFR, 'RR'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CLDREFFS, 'RS'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CLDREFFG, 'RG'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CLDNCCN,  'CLDNCCN' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT,DTDTFRIC, 'DTDTFRIC' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CLDBASEHGT,'CLDBASEHGT', RC=STATUS); VERIFY_(STATUS)


     call MAPL_GetPointer(EXPORT, PDF_AX,     'PDF_A',   RC=STATUS)
     VERIFY_(STATUS) 
     call MAPL_GetPointer(EXPORT,  WTHV2,      'WTHV2',     ALLOC=.TRUE., RC=STATUS)
     VERIFY_(STATUS) 
     call MAPL_GetPointer(EXPORT, WQL,        'WQL',        ALLOC=.TRUE., RC=STATUS)
     VERIFY_(STATUS)

!============ DG PDF diagnostics ==============
#ifdef PDFDIAG
     call MAPL_GetPointer(EXPORT, PDF_SIGW1,  'PDF_SIGW1', ALLOC=.TRUE.,  RC=STATUS)
     VERIFY_(STATUS) 
     call MAPL_GetPointer(EXPORT, PDF_SIGW2,  'PDF_SIGW2', ALLOC=.TRUE.,  RC=STATUS)
     VERIFY_(STATUS) 
     call MAPL_GetPointer(EXPORT, PDF_W1,     'PDF_W1', ALLOC=.TRUE.,    RC=STATUS)
     VERIFY_(STATUS) 
     call MAPL_GetPointer(EXPORT, PDF_W2,     'PDF_W2', ALLOC=.TRUE.,    RC=STATUS)
     VERIFY_(STATUS) 
     call MAPL_GetPointer(EXPORT, PDF_SIGTH1, 'PDF_SIGTH1', ALLOC=.TRUE., RC=STATUS)
     VERIFY_(STATUS) 
     call MAPL_GetPointer(EXPORT, PDF_SIGTH2, 'PDF_SIGTH2', ALLOC=.TRUE., RC=STATUS)
     VERIFY_(STATUS) 
     call MAPL_GetPointer(EXPORT, PDF_TH1,    'PDF_TH1',    ALLOC=.TRUE., RC=STATUS)
     VERIFY_(STATUS) 
     call MAPL_GetPointer(EXPORT, PDF_TH2,    'PDF_TH2',    ALLOC=.TRUE., RC=STATUS)
     VERIFY_(STATUS) 
     call MAPL_GetPointer(EXPORT, PDF_SIGQT1, 'PDF_SIGQT1', ALLOC=.TRUE., RC=STATUS)
     VERIFY_(STATUS) 
     call MAPL_GetPointer(EXPORT, PDF_SIGQT2, 'PDF_SIGQT2', ALLOC=.TRUE., RC=STATUS)
     VERIFY_(STATUS) 
     call MAPL_GetPointer(EXPORT, PDF_QT1,    'PDF_QT1',    ALLOC=.TRUE., RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT, PDF_QT2,    'PDF_QT2',    ALLOC=.TRUE., RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(EXPORT, PDF_RQTTH,  'PDF_RQTTH',  ALLOC=.TRUE., RC=STATUS)
     VERIFY_(STATUS) 
     call MAPL_GetPointer(EXPORT, PDF_RWTH,   'PDF_RWTH',   ALLOC=.TRUE., RC=STATUS)
     VERIFY_(STATUS) 
     call MAPL_GetPointer(EXPORT, PDF_RWQT,   'PDF_RWQT',   ALLOC=.TRUE., RC=STATUS)
     VERIFY_(STATUS) 
#endif

!!! shallow vars
      call MAPL_GetPointer(EXPORT, CBMF_SC,  'CBMF_SC' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, UMF_SC,  'UMF_SC' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MFD_SC,  'MFD_SC' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CUFRC_SC,'CUFRC_SC' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, PLCL_SC,'PLCL_SC' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, PLFC_SC,'PLFC_SC' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, PINV_SC,'PINV_SC' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, PREL_SC,'PREL_SC' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, PBUP_SC,'PBUP_SC' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CIN_SC,  'CIN_SC' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CNT_SC,  'CNT_SC' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CNB_SC,  'CNB_SC' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CNB_SC,  'CNB_SC' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DQRDT_SC,  'DQRDT_SC' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DQSDT_SC,  'DQSDT_SC' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DTHDT_SC,  'DTHDT_SC'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DQVDT_SC,  'DQVDT_SC'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DQIDT_SC,  'DQIDT_SC'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DQLDT_SC,  'DQLDT_SC'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DQCDT_SC,  'DQCDT_SC'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DUDT_SC,   'DUDT_SC'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DVDT_SC,   'DVDT_SC'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, WLCL_SC,   'WLCL_SC'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, QTSRC_SC,  'QTSRC_SC'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, THLSRC_SC, 'THLSRC_SC' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, THVLSRC_SC,'THVLSRC_SC', RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, TKEAVG_SC, 'TKEAVG_SC' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CLDTOP_SC, 'CLDTOP_SC' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, WUP_SC,    'WUP_SC'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, QTUP_SC,   'QTUP_SC'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, THLUP_SC,  'THLUP_SC'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, THVUP_SC,  'THVUP_SC'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, UUP_SC,    'UUP_SC'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, VUP_SC,    'VUP_SC'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, QCU_SC,    'QCU_SC'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, QLU_SC,    'QLU_SC'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, QIU_SC,    'QIU_SC'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, ENTR_SC,   'ENTR_SC'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DETR_SC,   'DETR_SC'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, XC_SC,     'XC_SC'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, QLDET_SC,  'QLDET_SC'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, QIDET_SC,  'QIDET_SC'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, QLENT_SC,  'QLENT_SC'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, QIENT_SC,  'QIENT_SC'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, QLSUB_SC,  'QLSUB_SC'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, QISUB_SC,  'QISUB_SC'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, SC_NDROP, 'SC_NDROP' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, SC_NICE,  'SC_NICE' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, SC_PRCP,  'SC_PRCP'  , RC=STATUS); VERIFY_(STATUS)


!      call MAPL_GetPointer(EXPORT, LIQANMOVE, 'LIQANMOVE' , RC=STATUS); VERIFY_(STATUS)
!      call MAPL_GetPointer(EXPORT, ICEANMOVE, 'ICEANMOVE' , RC=STATUS); VERIFY_(STATUS)
!      call MAPL_GetPointer(EXPORT, DANCLD,    'DANCLD'    , RC=STATUS); VERIFY_(STATUS)
!      call MAPL_GetPointer(EXPORT, DLSCLD,    'DLSCLD'    , RC=STATUS); VERIFY_(STATUS)
!      call MAPL_GetPointer(EXPORT, CURAINMOVE,'CURAINMOVE', RC=STATUS); VERIFY_(STATUS)
!      call MAPL_GetPointer(EXPORT, CUSNOWMOVE,'CUSNOWMOVE', RC=STATUS); VERIFY_(STATUS)
!      call MAPL_GetPointer(EXPORT, PFLCNMOVE, 'PFLCNMOVE',  RC=STATUS); VERIFY_(STATUS)
!      call MAPL_GetPointer(EXPORT, PFICNMOVE, 'PFICNMOVE',  RC=STATUS); VERIFY_(STATUS)

      !Trajectory for the moist TLM/ADJ
      call MAPL_GetPointer(EXPORT, TH_moist,   'TH_moist'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, Q_moist,    'Q_moist'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, KCBL_moist, 'KCBL_moist', RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, TS_moist,   'TS_moist'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, ctop_moist, 'ctop_moist', RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, KHu_moist,  'KHu_moist',  RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, KHl_moist,  'KHl_moist',  RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(EXPORT, REVSU_LSAN, 'REVSU_LSAN'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, REVSU_CN,   'REVSU_CN'    , RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(EXPORT, ACLL_CN,   'ACRLL_CN'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, ACLL_SC,   'ACRLL_SC'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, ACLL_AN,   'ACRLL_AN'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, ACLL_LS,   'ACRLL_LS'  , RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(EXPORT, ACIL_CN,   'ACRIL_CN'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, ACIL_SC,   'ACRIL_SC'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, ACIL_AN,   'ACRIL_AN'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, ACIL_LS,   'ACRIL_LS'  , RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(EXPORT, ACR_TOT,   'ACR_TOT'  , RC=STATUS); VERIFY_(STATUS)


      call MAPL_GetPointer(EXPORT, PFL_CN,   'PFL_CN'   , ALLOC=.TRUE. , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, PFL_SC,   'PFL_SC'   , ALLOC=.TRUE. , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, PFL_AN,   'PFL_AN'   , ALLOC=.TRUE. , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, PFL_LS,   'PFL_LS'   , ALLOC=.TRUE. , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, PFL_LSAN, 'PFL_LSAN' , ALLOC=.TRUE. , RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(EXPORT, PFI_CN,   'PFI_CN'   , ALLOC=.TRUE. , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, PFI_SC,   'PFI_SC'   , ALLOC=.TRUE. , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, PFI_AN,   'PFI_AN'   , ALLOC=.TRUE. , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, PFI_LS,   'PFI_LS'   , ALLOC=.TRUE. , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, PFI_LSAN, 'PFI_LSAN' , ALLOC=.TRUE. , RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(EXPORT, DPDTMST, 'DPDTMST' , RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(EXPORT, DlPDF,  'DLPDF'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DiPDF,  'DIPDF'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DlFIX,  'DLFIX'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DiFIX,  'DIFIX'  , RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(EXPORT, AUT,    'AUT'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, EVAPC,  'EVAPC'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, SDM,    'SDM'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, VFALLICE_AN, 'VFALLICE_AN' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, VFALLICE_LS, 'VFALLICE_LS' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, VFALLWAT_AN, 'VFALLWAT_AN' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, VFALLWAT_LS, 'VFALLWAT_LS' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, VFALLRN_AN, 'VFALLRN_AN' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, VFALLRN_LS, 'VFALLRN_LS' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, VFALLRN_CN, 'VFALLRN_CN' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, VFALLRN_SC, 'VFALLRN_SC' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, VFALLSN_AN, 'VFALLSN_AN' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, VFALLSN_LS, 'VFALLSN_LS' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, VFALLSN_CN, 'VFALLSN_CN' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, VFALLSN_SC, 'VFALLSN_SC' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, SUBLC,  'SUBLC'    , RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(EXPORT, FRZ_TT,    'FRZ_TT'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, FRZ_PP,    'FRZ_PP'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DCNVL,     'DCNVL'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DCNVI,     'DCNVI'    , RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(EXPORT, QSATl,     'QSATL'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, QSATi,     'QSATI'    , RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(EXPORT, ALPHT,     'ALPHT'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, ALPH1,     'ALPH1'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, ALPH2,     'ALPH2'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CFPDF,     'CFPDF'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CFPDFX,    'CFPDFX'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, RHCLR,     'RHCLR'    , RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(EXPORT, RCCODE,    'RCCODE'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, TRIEDLV,   'TRIEDLV'   , RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(EXPORT, QVRAS,     'QVRAS'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, THRAS,     'THRAS'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, URAS,      'URAS '   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, VRAS,      'VRAS '   , RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(EXPORT, SIGMA_DEEP,  'SIGMA_DEEP' , ALLOC=.TRUE., RC=STATUS); _VERIFY(STATUS)
      call MAPL_GetPointer(EXPORT, SIGMA_MID,  'SIGMA_MID' , ALLOC=.TRUE., RC=STATUS); _VERIFY(STATUS)

      call MAPL_GetPointer(EXPORT, RAS_TIME,  'RAS_TIME' , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, RAS_TRG,   'RAS_TRG'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, RAS_TOKI,  'RAS_TOKI' , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, RAS_PBL,   'RAS_PBL'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, RAS_WFN,   'RAS_WFN'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(EXPORT, MXDIAM,    'MXDIAM'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, RH600,     'RH600'     , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, Q600,      'Q600'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, QCBL,      'QCBL'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CNV_FRC,   'CNV_FRC'   , RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(EXPORT, PREVTOT,   'PREVTOT'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, PGENTOT,   'PGENTOT'   , RC=STATUS); VERIFY_(STATUS)


      call MAPL_GetPointer(EXPORT, FRZCZ,     'FRZCZ'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, FRZPZ,     'FRZPZ'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, COLLIZ,    'COLLIZ'  , RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(EXPORT, CNVLZ,     'CNVLZ'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, PDFLZ,     'PDFLZ'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CNVRNZ,    'CNVRNZ'   , RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(EXPORT, EVPCZ,     'EVPCZ'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, EVPPZ,     'EVPPZ'   , RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(EXPORT, CNVIZ,     'CNVIZ'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, PDFIZ,     'PDFIZ'   , RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(EXPORT, SUBCZ,     'SUBCZ'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, SUBPZ,     'SUBPZ'   , RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(EXPORT, AUTZ,      'AUTZ'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, COLLLZ,    'COLLLZ'  , RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(EXPORT, SDMZ,      'SDMZ'    , RC=STATUS); VERIFY_(STATUS)


      call MAPL_GetPointer(EXPORT, DDF_RH1,    'DDF_RH1'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DDF_RH2,    'DDF_RH2'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DDF_MUPH,   'DDF_MUPH'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DDF_BYNC,   'DDF_BYNC'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DDF_TC,     'DDF_TC'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DDF_QVC,    'DDF_QVC'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DDF_DQDT,   'DDF_DQDT'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DDF_DTDT,   'DDF_DTDT'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DDF_MFC,    'DDF_MFC'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DDF_ZSCALE, 'DDF_ZSCALE', RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, KEMST,      'KEMST'     , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, KEMST2,     'KEMST2'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, KEDISS,     'KEDISS'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, ENTLAM,     'ENTLAM'    , RC=STATUS); VERIFY_(STATUS)
      ! Aerosol Scavenging
      ! CAR 12/5/08
      call MAPL_GetPointer(EXPORT, DDU2gDT,     'DDU2gDT'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DSS2gDT,     'DSS2gDT'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DBC2gDT,     'DBC2gDT'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DOC2gDT,     'DOC2gDT'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DSU2gDT,     'DSU2gDT'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DNI2gDT,     'DNI2gDT'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DNH4A2gDT,   'DNH4A2gDT'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DNH32gDT,    'DNH32gDT'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DBRC2gDT,    'DBRC2gDT'   , RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(EXPORT, DDUDT,     'DDUDT'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DSSDT,     'DSSDT'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DBCDT,     'DBCDT'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DOCDT,     'DOCDT'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DSUDT,     'DSUDT'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DNIDT,     'DNIDT'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DNH4ADT,   'DNH4ADT'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DNH3DT,    'DNH3DT'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DBRCDT,    'DBRCDT'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DDUDTcarma,  'DDUDTcarma'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DSSDTcarma,  'DSSDTcarma'    , RC=STATUS); VERIFY_(STATUS)

      !AMM sync t,q extra diagnostic
      call MAPL_GetPointer(EXPORT, THMOIST,  'THMOIST' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT,  SMOIST,  'SMOIST'  , RC=STATUS); VERIFY_(STATUS)

      !--------Aerosol-Cloud interactions 

      call MAPL_GetPointer(EXPORT, SMAXL, 'SMAX_LIQ'      , RC=STATUS); VERIFY_(STATUS)    
      call MAPL_GetPointer(EXPORT, WSUB,  'WSUB'      , RC=STATUS); VERIFY_(STATUS)  
      call MAPL_GetPointer(EXPORT, CCN01, 'CCN01'         , RC=STATUS); VERIFY_(STATUS) 
      call MAPL_GetPointer(EXPORT, CCN04, 'CCN04'         , RC=STATUS); VERIFY_(STATUS) 
      call MAPL_GetPointer(EXPORT, CCN1,  'CCN1'          , RC=STATUS); VERIFY_(STATUS) 
      call MAPL_GetPointer(EXPORT, SMAXI, 'SMAX_ICE'      , RC=STATUS); VERIFY_(STATUS) 
      call MAPL_GetPointer(EXPORT, CDNC_NUC, 'CDNC_NUC'      , RC=STATUS); VERIFY_(STATUS) 
      call MAPL_GetPointer(EXPORT, INC_NUC, 'INC_NUC'      , RC=STATUS); VERIFY_(STATUS) 
      call MAPL_GetPointer(EXPORT, NCPL_VOL, 'NCPL_VOL'      , RC=STATUS); VERIFY_(STATUS)  
      call MAPL_GetPointer(EXPORT, NCPI_VOL, 'NCPI_VOL'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, SO4, 'SO4'      , RC=STATUS); VERIFY_(STATUS)  
      call MAPL_GetPointer(EXPORT, ORG, 'ORG'      , RC=STATUS); VERIFY_(STATUS)  
      call MAPL_GetPointer(EXPORT, DUST, 'DUST'      , RC=STATUS); VERIFY_(STATUS)  
      call MAPL_GetPointer(EXPORT, SEASALT, 'SEASALT'      , RC=STATUS); VERIFY_(STATUS)      
      call MAPL_GetPointer(EXPORT, BCARBON, 'BCARBON'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, NHET_NUC, 'NHET_NUC'      , RC=STATUS); VERIFY_(STATUS) 
      call MAPL_GetPointer(EXPORT, NLIM_NUC, 'NLIM_NUC'      , RC=STATUS); VERIFY_(STATUS) 
      call MAPL_GetPointer(EXPORT, SAT_RAT, 'SAT_RAT'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DQVDT_micro, 'DQVDT_micro'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DQIDT_micro, 'DQIDT_micro'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DQLDT_micro, 'DQLDT_micro'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DQRDT_micro, 'DQRDT_micro'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DQSDT_micro, 'DQSDT_micro'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DQGDT_micro, 'DQGDT_micro'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DQADT_micro, 'DQADT_micro'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DUDT_micro, 'DUDT_micro'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DVDT_micro, 'DVDT_micro'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DTDT_micro, 'DTDT_micro'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DTDT_macro, 'DTDT_macro'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DQVDT_macro, 'DQVDT_macro'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DQLDT_macro, 'DQLDT_macro'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DQIDT_macro, 'DQIDT_macro'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DQADT_macro, 'DQADT_macro'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DTDT_moist, 'DTDT_moist'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DTDTCN, 'DTDTCN'      , RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(EXPORT, TVQX1,     'TVQX1'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, TVQX2,     'TVQX2'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, RL_MASK, 'RL_MASK'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, RI_MASK, 'RI_MASK'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, KAPPA,   'KAPPA'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, SC_ICE,  'SC_ICE'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CFICE,   'CFICE'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CFLIQ,   'CFLIQ'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, RHICE,   'RHICE'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, RHLIQ,   'RHLIQ'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, NHET_IMM, 'NHET_IMM'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, NHET_DEP, 'NHET_DEP'      , RC=STATUS); VERIFY_(STATUS)      
      call MAPL_GetPointer(EXPORT, DUST_IMM, 'DUST_IMM'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DUST_DEP, 'DUST_DEP'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, SCF     , 'SCF'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, SCF_ALL,  'SCF_ALL'      , RC=STATUS); VERIFY_(STATUS)      
      call MAPL_GetPointer(EXPORT, SIGW_GW,   'SIGW_GW'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, SIGW_TURB, 'SIGW_TURB'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, SIGW_CNV,   'SIGW_CNV'      , RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(EXPORT, SIGW_RC,   'SIGW_RC'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CNV_FICE,  'CNV_FICE' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CNV_NDROP, 'CNV_NDROP' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CNV_NICE,  'CNV_NICE' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, RHCmicro,  'RHCmicro' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CCNCOLUMN,      'CCNCOLUMN'     , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, NDCOLUMN,      'NDCOLUMN'     , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, NCCOLUMN,      'NCCOLUMN'     , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, BERG, 'BERG'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, BERGS, 'BERGS'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, MELT, 'MELT'      , RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(EXPORT, QCRES,   'QCRES'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, QIRES,   'QIRES'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, AUTICE,  'AUTICE'      , RC=STATUS); VERIFY_(STATUS)      
      call MAPL_GetPointer(EXPORT, DT_RASP, 'DT_RASP'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, FRZPP_LS,   'FRZPP_LS'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, SNOWMELT_LS, 'SNOWMELT_LS'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DNHET_IMM, 'DNHET_IMM'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DNHET_CT, 'DNHET_CT'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, PFRZ, 'PFRZ'      , RC=STATUS); VERIFY_(STATUS)  
      call MAPL_GetPointer(EXPORT, DNCNUC, 'DNCNUC'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DNCHMSPLIT, 'DNCHMSPLIT'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DNCAUTICE, 'DNCAUTICE'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DNCSUBL, 'DNCSUBL'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DNCACRIS, 'DNCACRIS'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DNDCCN, 'DNDCCN'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DNDACRLS, 'DNDACRLS'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DNDEVAPC, 'DNDEVAPC'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DNDACRLR, 'DNDACRLR'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DNDAUTLIQ, 'DNDAUTLIQ'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DNDCNV, 'DNDCNV'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DNCCNV, 'DNCCNV'      , RC=STATUS); VERIFY_(STATUS)


      call MAPL_GetPointer(EXPORT, CLDREFFL_TOP,      'CLDREFFL_TOP'     , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CLDREFFI_TOP,      'CLDREFFI_TOP'     , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, NCPL_TOP,      'NCPL_TOP'     , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, NCPI_TOP,      'NCPI_TOP'     , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(EXPORT, NCPL_CLDBASE,      'NCPL_CLDBASE'     , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, LWC,      'LWC'     , RC=STATUS); VERIFY_(STATUS)      
      call MAPL_GetPointer(EXPORT, IWC,      'IWC'     , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, QCVAR_EXP,      'QCVAR_EXP',  ALLOC=.TRUE.  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, LTS,      'LTS',  ALLOC=.TRUE.  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, EIS,      'EIS',  ALLOC=.TRUE.  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, RAS_ALPHA, 'RAS_ALPHA' , ALLOC=.TRUE.  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, RAS_TAU, 'RAS_TAU' , ALLOC=.TRUE.  , RC=STATUS); VERIFY_(STATUS)

       call MAPL_GetPointer(EXPORT, DNDDT,    'DNDDT'   , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(EXPORT, DNCDT,    'DNCDT'   , RC=STATUS); VERIFY_(STATUS)
      
      call MAPL_GetPointer(EXPORT, DTSX  ,   'DTSX'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, KHX   ,   'KHX'     , RC=STATUS); VERIFY_(STATUS)
      
      !---------

!--kml--- activation for single-moment uphysics      
!!!      if(adjustl(CLDMICRO)=="1MOMENT" .and. USE_AEROSOL_NN) then 
        call MAPL_GetPointer(EXPORT, NACTL,'NACTL' ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS);NACTL=0.0
        call MAPL_GetPointer(EXPORT, NACTI,'NACTI' ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS);NACTI=0.0
!!!      endif
!--kml-----------------------------------------------------------------------------

      ! Count the fields in TR...
      !--------------------------

      call ESMF_FieldBundleGet(TR,FieldCount=KM , RC=STATUS)
      VERIFY_(STATUS)


      ! Allocate tracer stuff
      !----------------------

      allocate(IS_FRIENDLY(KM),stat=STATUS)
      VERIFY_(STATUS)
      allocate(IS_WEIGHTED(KM),stat=STATUS)
      VERIFY_(STATUS)
      allocate(TRPtrs     (KM),stat=STATUS)
      VERIFY_(STATUS)
      allocate(FSCAV_     (KM),stat=STATUS)
      VERIFY_(STATUS)
      allocate(QNAMES(KM), CNAMES(KM), stat=STATUS)
      VERIFY_(STATUS)
      IF(ADJUSTL(CONVPAR_OPTION) == 'GF') THEN
        ALLOCATE(Hcts(KM), stat=STATUS); VERIFY_(STATUS)
	IF(STATUS==0) THeN
              Hcts(1:KM)%hstar = -1.
              Hcts(1:KM)%dhr   = -1.
              Hcts(1:KM)%ak0   = -1.
              Hcts(1:KM)%dak   = -1.
        ENDIF
      ENDIF

      QNAMES = " "
      CNAMES = " "
      QNAME = " "
      CNAME = " "

      !CAR get the names of things
      call ESMF_FieldBundleGet(TR, fieldNameList=QNAMES, rc=STATUS )
      VERIFY_(STATUS) 

      ! Loop over all quantities to be diffused.
      !----------------------------------------

      ITRCR = 0

      do K=1,KM

         ! Get the Kth Field from tracer bundle
         !-------------------------------------

         NAME = trim(QNAMES(K))

         call ESMF_FieldBundleGet(TR, fieldName=trim(NAME), Field=FIELD, RC=STATUS)
         VERIFY_(STATUS)

         ! Get items friendly status (default is not friendly)
         !-----------------------------------------------------

         call ESMF_AttributeGet  (FIELD, NAME="FriendlyToMOIST",isPresent=isPresent, RC=STATUS)
         VERIFY_(STATUS)
         if(isPresent) then
            call ESMF_AttributeGet  (FIELD, NAME="FriendlyToMOIST",VALUE=IS_FRIENDLY(K), RC=STATUS)
            VERIFY_(STATUS)
         else
            IS_FRIENDLY(K) = .false.
         end if

         ! Get items weighting (default is unweighted tendencies)
         !--------------------------------------------------------

         call ESMF_AttributeGet  (FIELD, NAME="WeightedTendency",isPresent=isPresent, RC=STATUS)
         VERIFY_(STATUS)
         if(isPresent) then
            call ESMF_AttributeGet  (FIELD, NAME="WeightedTendency",VALUE=IS_WEIGHTED(K), RC=STATUS)
            VERIFY_(STATUS)
         else
            IS_WEIGHTED(K) = .false.
         end if

         ! Get items scavenging fraction
         !-------------------------------

         call ESMF_AttributeGet  (FIELD, NAME="ScavengingFractionPerKm",isPresent=isPresent, RC=STATUS)
         VERIFY_(STATUS)
         if(isPresent) then
            call ESMF_AttributeGet  (FIELD, NAME="ScavengingFractionPerKm",VALUE=FSCAV_(K), RC=STATUS)
            VERIFY_(STATUS)
         else
            FSCAV_(K) = 0.0 ! no scavenging
         end if

         ! Get items for the wet removal parameterization for gases based on the Henry's Law
         !-------------------------------
         IF(ADJUSTL(CONVPAR_OPTION) == 'GF') THEN
           Vect_Hcts(:)=-99.
           call ESMF_AttributeGet  (FIELD,"SetofHenryLawCts",isPresent=isPresent,  RC=STATUS)
           VERIFY_(STATUS)
           if (isPresent) then
              call ESMF_AttributeGet  (FIELD,"SetofHenryLawCts",Vect_Hcts,  RC=STATUS)
              !.. if (MAPL_AM_I_ROOT()) then
                 !.. PRINT*,"spcname=",k,trim(QNAMES(K)),FSCAV_(K)
                 !.. print*,"Vect_Hcts=",Vect_Hcts(:),statUS
                 !.. flush(6)
              !.. ENDIF
              Hcts(k)%hstar = Vect_Hcts(1)
              Hcts(k)%dhr   = Vect_Hcts(2)
              Hcts(k)%ak0   = Vect_Hcts(3)
              Hcts(k)%dak   = Vect_Hcts(4)
           ENDIF
         ENDIF

         ! Check aerosol names
         ! CAR 12/5/08

         !PRINT *, "*********NAME CHECKING:*******"
         !PRINT *, trim(QNAMES(K))

         ! Remove qualifier from variable name:  GOCART::du001->du001
         ! CAR 12/5/08
         !-------------------------------

         QNAME = QNAMES(K)
         ind= index(QNAME, '::')
         if (ind> 0) then
            CNAMES(K) = trim(QNAME(1:ind-1))  ! Component name (e.g., GOCART, CARMA)
            QNAMES(K) = trim(QNAME(ind+2:))
         end if
         !PRINT *, "******CROPPED NAME CHECKING*******"
         !PRINT *, trim(QNAMES(K)), FSCAV_(K)


         ! Get pointer to the quantity, its tendency, its surface value,
         !   the surface flux, and the sensitivity of the surface flux.
         ! -------------------------------------------------------------

         if (IS_FRIENDLY(K)) then
            ITRCR = ITRCR + 1
            call ESMFL_BundleGetPointerToData(TR    , trim(NAME),        TRPtrs (K)%Q, RC=STATUS)
            VERIFY_(STATUS)
         else
            ASSERT_(.not.associated(TRPtrs (K)%Q))
            TRPtrs(K)%Q  => null()
         end if

      end do
      ! CAR END LOOP OVER K

      ! Allocate space for tracers
      !---------------------------

      allocate(XHO (IM,JM,LM,ITRCR),stat=STATUS)
      VERIFY_(STATUS)
      !     FSCAV changes dimensions of FSCAV_
      allocate(FSCAV(ITRCR),stat=STATUS)
      VERIFY_(STATUS)

      ALLOC_PTYPE     = .not.associated(PTYPE)
      if(ALLOC_PTYPE) allocate(PTYPE(IM,JM))

      ALLOC_CNV_DQLDT = .not.associated(CNV_DQLDT )
      ALLOC_CNV_MF0   = .not.associated(CNV_MF0   )
      ALLOC_CNV_MFD   = .not.associated(CNV_MFD   )
      ALLOC_CNV_MFC   = .not.associated(CNV_MFC   )
      ALLOC_CNV_TOPP  = .not.associated(CNV_TOPP  )
      ALLOC_CNV_UPDF  = .not.associated(CNV_UPDF  )
      ALLOC_CNV_CVW   = .not.associated(CNV_CVW   )
      ALLOC_CNV_QC    = .not.associated(CNV_QC    )
      ALLOC_RAD_CF    = .not.associated(RAD_CF    )
      ALLOC_RAD_QV    = .not.associated(RAD_QV    )
      ALLOC_RAD_QL    = .not.associated(RAD_QL    )
      ALLOC_RAD_QI    = .not.associated(RAD_QI    )
      ALLOC_RAD_QR    = .not.associated(RAD_QR    )
      ALLOC_RAD_QS    = .not.associated(RAD_QS    )
      ALLOC_RAD_QG    = .not.associated(RAD_QG    )
      ALLOC_CLDREFFL  = .not.associated(CLDREFFL  )
      ALLOC_CLDREFFI  = .not.associated(CLDREFFI  )
      ALLOC_CLDREFFR  = .not.associated(CLDREFFR  )
      ALLOC_CLDREFFS  = .not.associated(CLDREFFS  )
      ALLOC_CLDREFFG  = .not.associated(CLDREFFG  )
      ALLOC_CLDNCCN   = .not.associated(CLDNCCN   )
      ALLOC_ENTLAM    = .not.associated(ENTLAM    )

      if(ALLOC_CNV_DQLDT) allocate(CNV_DQLDT(IM,JM,LM))
      if(ALLOC_CNV_MF0  ) allocate(CNV_MF0  (IM,JM,LM))
      if(ALLOC_CNV_MFD  ) allocate(CNV_MFD  (IM,JM,LM))
      if(ALLOC_CNV_MFC  ) allocate(CNV_MFC  (IM,JM,0:LM))
      if(ALLOC_CNV_TOPP ) allocate(CNV_TOPP (IM,JM))
      if(ALLOC_CNV_UPDF ) allocate(CNV_UPDF (IM,JM,LM))
      if(ALLOC_CNV_CVW  ) allocate(CNV_CVW  (IM,JM,LM))
      if(ALLOC_CNV_QC   ) allocate(CNV_QC   (IM,JM,LM))
      if(ALLOC_RAD_CF   ) allocate(RAD_CF   (IM,JM,LM))
      if(ALLOC_RAD_QV   ) allocate(RAD_QV   (IM,JM,LM))
      if(ALLOC_RAD_QL   ) allocate(RAD_QL   (IM,JM,LM))
      if(ALLOC_RAD_QI   ) allocate(RAD_QI   (IM,JM,LM))
      if(ALLOC_RAD_QR   ) allocate(RAD_QR   (IM,JM,LM))
      if(ALLOC_RAD_QS   ) allocate(RAD_QS   (IM,JM,LM))
      if(ALLOC_RAD_QG   ) allocate(RAD_QG   (IM,JM,LM))
      if(ALLOC_CLDREFFL ) allocate(CLDREFFL (IM,JM,LM))
      if(ALLOC_CLDREFFI ) allocate(CLDREFFI (IM,JM,LM))
      if(ALLOC_CLDREFFR ) allocate(CLDREFFR (IM,JM,LM))
      if(ALLOC_CLDREFFS ) allocate(CLDREFFS (IM,JM,LM))
      if(ALLOC_CLDREFFG ) allocate(CLDREFFG (IM,JM,LM))
      if(ALLOC_CLDNCCN  ) allocate(CLDNCCN  (IM,JM,LM))
      if(ALLOC_ENTLAM   ) allocate(ENTLAM   (IM,JM,LM))

      ALLOC_DQRC     = .not. associated( DQRC  )

      if(ALLOC_DQRC  ) allocate ( DQRC     (IM,JM,LM))

      ! MAT These are 2MOMENT exports that are now necessary
      !     inside of RASE
      ALLOC_CNV_FICE     = .not.associated(CNV_FICE)
      ALLOC_CNV_NDROP    = .not.associated(CNV_NDROP)
      ALLOC_CNV_NICE     = .not.associated(CNV_NICE)           

      if(ALLOC_CNV_FICE  ) allocate (CNV_FICE (IM,JM,LM))
      if(ALLOC_CNV_NDROP  ) allocate (CNV_NDROP (IM,JM,LM))
      if(ALLOC_CNV_NICE  ) allocate (CNV_NICE (IM,JM,LM))              

      ALLOC_CFICE          = .not.associated(CFICE)
      ALLOC_CFLIQ          = .not.associated(CFLIQ)

      if(ALLOC_CFICE     ) allocate (CFICE  (IM,JM,LM)) 
      if(ALLOC_CFLIQ     ) allocate (CFLIQ   (IM,JM,LM))   

      if(adjustl(CLDMICRO)=="GFDL") then
         ALLOC_PRCP_RAIN    = .not.associated(PRCP_RAIN)
         ALLOC_PRCP_SNOW    = .not.associated(PRCP_SNOW)
         ALLOC_PRCP_ICE     = .not.associated(PRCP_ICE)
         ALLOC_PRCP_GRAUPEL = .not.associated(PRCP_GRAUPEL)
         if(ALLOC_PRCP_RAIN     ) allocate(PRCP_RAIN    (IM,JM))
         if(ALLOC_PRCP_SNOW     ) allocate(PRCP_SNOW    (IM,JM))
         if(ALLOC_PRCP_ICE      ) allocate(PRCP_ICE     (IM,JM))
         if(ALLOC_PRCP_GRAUPEL  ) allocate(PRCP_GRAUPEL (IM,JM))
      endif

      if(adjustl(CLDMICRO)=="2MOMENT" .or. adjustl(CLDMICRO)=="GFDL") then
         ALLOC_DQVDT_micro     = .not.associated(DQVDT_micro)
         ALLOC_DQIDT_micro     = .not.associated(DQIDT_micro)
         ALLOC_DQLDT_micro     = .not.associated(DQLDT_micro)
         ALLOC_DQRDT_micro     = .not.associated(DQRDT_micro)
         ALLOC_DQSDT_micro     = .not.associated(DQSDT_micro)
         ALLOC_DQGDT_micro     = .not.associated(DQGDT_micro)
         ALLOC_DQADT_micro     = .not.associated(DQADT_micro)
         ALLOC_DUDT_micro      = .not.associated(DUDT_micro)
         ALLOC_DVDT_micro      = .not.associated(DVDT_micro)
         ALLOC_DTDT_micro      = .not.associated(DTDT_micro)
         if(ALLOC_DQVDT_micro   ) allocate(DQVDT_micro   (IM,JM,LM))
         if(ALLOC_DQIDT_micro   ) allocate(DQIDT_micro   (IM,JM,LM))
         if(ALLOC_DQLDT_micro   ) allocate(DQLDT_micro   (IM,JM,LM))
         if(ALLOC_DQRDT_micro   ) allocate(DQRDT_micro   (IM,JM,LM))
         if(ALLOC_DQSDT_micro   ) allocate(DQSDT_micro   (IM,JM,LM))
         if(ALLOC_DQGDT_micro   ) allocate(DQGDT_micro   (IM,JM,LM))
         if(ALLOC_DQADT_micro   ) allocate(DQADT_micro   (IM,JM,LM))
         if(ALLOC_DUDT_micro    ) allocate(DUDT_micro    (IM,JM,LM))
         if(ALLOC_DVDT_micro    ) allocate(DVDT_micro    (IM,JM,LM))
         if(ALLOC_DTDT_micro    ) allocate(DTDT_micro    (IM,JM,LM))
      ! CLDMACRO stuff
         ALLOC_PFRZ         = .not.associated(PFRZ)
         ALLOC_DTDT_macro   = .not.associated(DTDT_macro)
         ALLOC_DQVDT_macro  = .not.associated(DQVDT_macro)
         ALLOC_DQLDT_macro  = .not.associated(DQLDT_macro)
         ALLOC_DQIDT_macro  = .not.associated(DQIDT_macro)
         ALLOC_DQADT_macro  = .not.associated(DQADT_macro)
         ALLOC_SC_ICE       = .not.associated(SC_ICE)
         ALLOC_DT_RASP      = .not.associated(DT_RASP)
         if(ALLOC_PFRZ )       allocate(PFRZ      (IM,JM,LM))
         if(ALLOC_DTDT_macro ) allocate(DTDT_macro(IM,JM,LM))
         if(ALLOC_DQVDT_macro) allocate(DQVDT_macro(IM,JM,LM))
         if(ALLOC_DQLDT_macro) allocate(DQLDT_macro(IM,JM,LM))
         if(ALLOC_DQIDT_macro) allocate(DQIDT_macro(IM,JM,LM))
         if(ALLOC_DQADT_macro) allocate(DQADT_macro(IM,JM,LM))
         if(ALLOC_SC_ICE     ) allocate(SC_ICE    (IM,JM,LM))
         if(ALLOC_DT_RASP    ) allocate(DT_RASP   (IM,JM,LM))
      endif

      allocate(FQAl(IM,JM,LM), stat=STATUS)
      VERIFY_(STATUS)
      allocate(FQAI(IM,JM,LM), stat=STATUS)
      VERIFY_(STATUS)
      allocate(FQA(IM,JM,LM), stat=STATUS)
      VERIFY_(STATUS)

      if(adjustl(CLDMICRO)=="2MOMENT") then

         allocate(QCNTOT(IM,JM,LM), stat=STATUS)
         VERIFY_(STATUS)

         !Aerosol-Cloud Interactions ---------------------
         ALLOC_SMAXL    = .not.associated(SMAXL    ) !DONIF
         ALLOC_SMAXI     = .not.associated(SMAXI    )
         ALLOC_WSUB     = .not.associated(WSUB     )
         ALLOC_CCN01    = .not.associated(CCN01    )
         ALLOC_CCN04    = .not.associated(CCN04    )
         ALLOC_CCN1     = .not.associated(CCN1     )
         ALLOC_NHET_NUC    = .not.associated(NHET_NUC   )
         ALLOC_NLIM_NUC     = .not.associated(NLIM_NUC   )
         ALLOC_SO4        = .not.associated(SO4   )
         ALLOC_ORG        = .not.associated(ORG   )
         ALLOC_BCARBON    = .not.associated(BCARBON  )
         ALLOC_DUST       = .not.associated(DUST   )
         ALLOC_SEASALT    = .not.associated(SEASALT   )
         ALLOC_NCPL_VOL   = .not.associated(NCPL_VOL  )
         ALLOC_NCPI_VOL   = .not.associated(NCPI_VOL   )

         ALLOC_CDNC_NUC   = .not.associated(CDNC_NUC )
         ALLOC_INC_NUC    = .not.associated(INC_NUC   )
         ALLOC_SAT_RAT    = .not.associated(SAT_RAT )
         ALLOC_QSTOT    = .not.associated(QSTOT)
         ALLOC_QRTOT    = .not.associated(QRTOT)
         ALLOC_QPTOTLS    = .not.associated(QPTOTLS)

         ALLOC_DTDT_moist     = .not.associated(DTDT_moist)
         ALLOC_DTDTCN     = .not.associated(DTDTCN)    

         ALLOC_RL_MASK        = .not.associated(RL_MASK)
         ALLOC_RI_MASK        = .not.associated(RI_MASK)
         ALLOC_KAPPA          = .not.associated(KAPPA)
         ALLOC_RHICE          = .not.associated(RHICE)
         ALLOC_RHLIQ          = .not.associated(RHLIQ)

         ALLOC_NHET_IMM     = .not.associated(NHET_IMM)    
         ALLOC_NHET_DEP     = .not.associated(NHET_DEP)
         ALLOC_DUST_IMM     = .not.associated(DUST_IMM)
         ALLOC_DUST_DEP     = .not.associated(DUST_DEP)

         ALLOC_SCF          = .not.associated(SCF)
         ALLOC_SCF_ALL      = .not.associated(SCF_ALL)     
         ALLOC_SIGW_GW      = .not.associated(SIGW_GW)
         ALLOC_SIGW_CNV      = .not.associated(SIGW_CNV)

         ALLOC_SIGW_TURB    = .not.associated(SIGW_TURB)
         ALLOC_SIGW_RC      = .not.associated(SIGW_RC)
         ALLOC_RHCmicro     = .not.associated(RHCmicro) 
         ALLOC_DNHET_IMM     = .not.associated(DNHET_IMM)
         ALLOC_BERG            = .not.associated(BERG)
         ALLOC_BERGS           = .not.associated(BERGS)
         ALLOC_MELT            = .not.associated(MELT)
         ALLOC_DNHET_CT     = .not.associated(DNHET_CT)

         ALLOC_QCRES       = .not.associated(QCRES)
         ALLOC_QIRES   = .not.associated(QIRES)
         ALLOC_AUTICE   = .not.associated(AUTICE)
         ALLOC_FRZPP_LS        = .not.associated(FRZPP_LS)
         ALLOC_SNOWMELT_LS     = .not.associated(SNOWMELT_LS)

         ALLOC_DNCNUC        = .not.associated(DNCNUC)
         ALLOC_DNCSUBL       = .not.associated(DNCSUBL)
         ALLOC_DNCHMSPLIT    = .not.associated(DNCHMSPLIT)      
         ALLOC_DNCAUTICE     = .not.associated(DNCAUTICE)
         ALLOC_DNCACRIS      = .not.associated(DNCACRIS)      
         ALLOC_DNDCCN        = .not.associated(DNDCCN)
         ALLOC_DNDACRLS        = .not.associated(DNDACRLS)
         ALLOC_DNDACRLR        = .not.associated(DNDACRLR)
         ALLOC_DNDEVAPC        = .not.associated(DNDEVAPC)
         ALLOC_DNDAUTLIQ        = .not.associated(DNDAUTLIQ)      
         ALLOC_DNDCNV       = .not.associated(DNDCNV)
         ALLOC_DNCCNV        = .not.associated(DNCCNV)





         if(ALLOC_SMAXL   ) allocate(SMAXL   (IM,JM,LM))
         if(ALLOC_SMAXI   ) allocate(SMAXI   (IM,JM,LM))
         if(ALLOC_WSUB    ) allocate(WSUB    (IM,JM,LM)) !DONIF
         if(ALLOC_CCN01   ) allocate(CCN01   (IM,JM,LM))
         if(ALLOC_CCN04   ) allocate(CCN04   (IM,JM,LM))
         if(ALLOC_CCN1    ) allocate(CCN1    (IM,JM,LM))
         if(ALLOC_NHET_NUC   ) allocate(NHET_NUC   (IM,JM,LM))
         if(ALLOC_NLIM_NUC   ) allocate(NLIM_NUC   (IM,JM,LM))
         if(ALLOC_SO4        ) allocate(SO4        (IM,JM,LM))
         if(ALLOC_ORG        ) allocate(ORG        (IM,JM,LM))
         if(ALLOC_BCARBON    ) allocate(BCARBON    (IM,JM,LM))
         if(ALLOC_DUST       ) allocate(DUST       (IM,JM,LM))
         if(ALLOC_SEASALT    ) allocate(SEASALT    (IM,JM,LM))
         if(ALLOC_NCPI_VOL   ) allocate(NCPI_VOL   (IM,JM,LM))
         if(ALLOC_NCPL_VOL   ) allocate(NCPL_VOL   (IM,JM,LM))

         if(ALLOC_CDNC_NUC   ) allocate(CDNC_NUC   (IM,JM,LM))
         if(ALLOC_INC_NUC    ) allocate(INC_NUC    (IM,JM,LM))
         if(ALLOC_SAT_RAT    ) allocate(SAT_RAT    (IM,JM,LM))
         if(ALLOC_QSTOT    ) allocate(QSTOT    (IM,JM,LM))
         if(ALLOC_QRTOT    ) allocate(QRTOT    (IM,JM,LM))
         if(ALLOC_QPTOTLS    ) allocate(QPTOTLS    (IM,JM,LM))

         if(ALLOC_DTDT_moist    ) allocate(DTDT_moist    (IM,JM,LM))    
         if(ALLOC_DTDTCN) allocate(DTDTCN  (IM,JM,LM))    

         if(ALLOC_RL_MASK   ) allocate(RL_MASK (IM,JM,LM))  
         if(ALLOC_RI_MASK   ) allocate(RI_MASK (IM,JM,LM))  
         if(ALLOC_KAPPA     ) allocate(KAPPA   (IM,JM,LM))  
         if(ALLOC_RHICE     ) allocate(RHICE   (IM,JM,LM)) 
         if(ALLOC_RHLIQ     ) allocate(RHLIQ   (IM,JM,LM)) 

         if(ALLOC_NHET_IMM  ) allocate (NHET_IMM   (IM,JM,LM))     
         if(ALLOC_NHET_DEP  ) allocate (NHET_DEP   (IM,JM,LM))
         if(ALLOC_DUST_IMM  ) allocate (DUST_IMM   (IM,JM,LM))
         if(ALLOC_DUST_DEP  ) allocate (DUST_DEP   (IM,JM,LM))


         if(ALLOC_SCF )      allocate (SCF     (IM,JM,LM))
         if(ALLOC_SCF_ALL  ) allocate (SCF_ALL (IM,JM,LM))      
         if(ALLOC_SIGW_GW  ) allocate (SIGW_GW (IM,JM,LM))
         if(ALLOC_SIGW_CNV  ) allocate (SIGW_CNV (IM,JM,LM))

         if(ALLOC_SIGW_TURB  ) allocate (SIGW_TURB (IM,JM,LM))
         if(ALLOC_SIGW_RC  ) allocate (SIGW_RC (IM,JM,LM))
         if(ALLOC_RHCmicro  ) allocate (RHCmicro (IM,JM,LM)) 
         if(ALLOC_DNHET_IMM  ) allocate (DNHET_IMM   (IM,JM,LM))
         if(ALLOC_BERG  )  allocate  (BERG   (IM,JM,LM))
         if(ALLOC_BERGS )  allocate  (BERGS   (IM,JM,LM))
         if(ALLOC_MELT  )  allocate  (MELT   (IM,JM,LM))
         if(ALLOC_DNHET_CT  )  allocate  (DNHET_CT   (IM,JM,LM))      


         if(ALLOC_QCRES  )    allocate  (QCRES   (IM,JM,LM))
         if(ALLOC_QIRES  )    allocate  (QIRES   (IM,JM,LM))
         if(ALLOC_AUTICE  )    allocate  (AUTICE   (IM,JM,LM))      
         if(ALLOC_FRZPP_LS  )    allocate  (FRZPP_LS   (IM,JM,LM))
         if(ALLOC_SNOWMELT_LS )  allocate  (SNOWMELT_LS   (IM,JM,LM))   

         if(ALLOC_DNCNUC  ) allocate (DNCNUC   (IM,JM,LM))   
         if(ALLOC_DNCSUBL ) allocate (DNCSUBL   (IM,JM,LM))
         if(ALLOC_DNCHMSPLIT  ) allocate (DNCHMSPLIT   (IM,JM,LM))
         if(ALLOC_DNCACRIS  ) allocate (DNCACRIS   (IM,JM,LM))
         if(ALLOC_DNCAUTICE ) allocate (DNCAUTICE   (IM,JM,LM))      
         if(ALLOC_DNDCCN ) allocate (DNDCCN   (IM,JM,LM))
         if(ALLOC_DNDACRLS ) allocate (DNDACRLS   (IM,JM,LM))
         if(ALLOC_DNDACRLR ) allocate (DNDACRLR   (IM,JM,LM))
         if(ALLOC_DNDEVAPC ) allocate (DNDEVAPC   (IM,JM,LM))
         if(ALLOC_DNDAUTLIQ ) allocate (DNDAUTLIQ   (IM,JM,LM))      
         if(ALLOC_DNDCNV ) allocate (DNDCNV   (IM,JM,LM))
         if(ALLOC_DNCCNV ) allocate (DNCCNV   (IM,JM,LM))




        ! if (associated(CLDREFFI_TOP)) allocate (CLDREFFI_TOP(IM, JM)) 
        ! if (associated(CLDREFFL_TOP)) allocate (CLDREFFL_TOP(IM, JM))
        ! if (associated(NCPI_TOP)) allocate (NCPI_TOP(IM, JM)) 
        ! if (associated(NCPL_TOP)) allocate (NCPL_TOP(IM, JM)) 
        ! if (associated(QCVAR_EXP)) allocate (QCVAR_EXP(IM, JM))  
                  



             
	   

      end if ! 2MOMENT


      !-----------------------------------------------------------------

      !  Allocate local space for some diagnostics that have to be in arg list to progno_cloud
      allocate(QSN(IM,JM,LM))
      allocate(QRN(IM,JM,LM))
      allocate(QPLS(IM,JM,LM))
      QSN = 0.
      QRN = 0.
      QPLS  = 0.

      ! Recording of import/internal vars into export if desired
      !---------------------------------------------------------

      call MAPL_GetPointer(EXPORT, Qx0,     'QX0'       , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, QLLSx0,  'QLLSX0'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, QLCNx0,  'QLCNX0'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CLCNx0,  'CLCNX0'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CLLSx0,  'CLLSX0'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, QILSx0,  'QILSX0'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, QICNx0,  'QICNX0'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, QCCNx0,  'QCCNX0'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, QCLSx0,  'QCLSX0'    , RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(EXPORT, KHx0,    'KHX0'     , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, THx0,    'THX0'     , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, Ux0,     'UX0'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, Vx0,     'VX0'      , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, TSx0,    'TSX0'     , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, FRLANDx0,'FRLANDX0' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, W_DIAG,  'W_DIAG'   , RC=STATUS); VERIFY_(STATUS)


      if(associated(Qx0       )) Qx0        = Q
      if(associated(QLLSx0    )) QLLSx0     = QLLS
      if(associated(QILSx0    )) QILSx0     = QILS
      if(associated(QLCNx0    )) QLCNx0     = QLCN
      if(associated(QICNx0    )) QICNx0     = QICN
      if(associated(CLLSx0    )) CLLSx0     = CLLS
      if(associated(CLCNx0    )) CLCNx0     = CLCN
      if(associated(QCCNx0    )) QCCNx0     = QICN+QLCN
      if(associated(QCLSx0    )) QCLSx0     = QILS+QLLS

      if(associated(KHx0      )) KHx0       = KH
      if(associated(THx0      )) THx0       = TH
      if(associated(Ux0       )) Ux0        = U
      if(associated(Vx0       )) Vx0        = V
      if(associated(TSx0      )) TSx0       = TS
      if(associated(FRLANDx0  )) FRLANDx0   = FRLAND

      !------ allocate UW vars ------

      ALLOC_UMF    = .not.associated(UMF_SC    )
      ALLOC_MFD    = .not.associated(MFD_SC    )
      ALLOC_DQVDT  = .not.associated(DQVDT_SC  )
      ALLOC_DQLDT  = .not.associated(DQLDT_SC  )
      ALLOC_DQIDT  = .not.associated(DQIDT_SC  )
      ALLOC_DTHDT  = .not.associated(DTHDT_SC  )
      ALLOC_DUDT   = .not.associated(DUDT_SC   )
      ALLOC_DVDT   = .not.associated(DVDT_SC   )
      ALLOC_DQRDT  = .not.associated(DQRDT_SC  )
      ALLOC_DQSDT  = .not.associated(DQSDT_SC  )
      ALLOC_CUFRC  = .not.associated(CUFRC_SC  )
      ALLOC_QCU    = .not.associated(QCU_SC    )
      ALLOC_QLU    = .not.associated(QLU_SC    )
      ALLOC_QIU    = .not.associated(QIU_SC    )
      ALLOC_CBMF   = .not.associated(CBMF_SC   )
      ALLOC_DQCDT  = .not.associated(DQCDT_SC  )
      ALLOC_CNT    = .not.associated(CNT_SC    )
      ALLOC_CNB    = .not.associated(CNB_SC    )
      ALLOC_CIN    = .not.associated(CIN_SC    )
      ALLOC_PLCL   = .not.associated(PLCL_SC    )
      ALLOC_PLFC   = .not.associated(PLFC_SC    )
      ALLOC_PINV   = .not.associated(PINV_SC    )
      ALLOC_PREL   = .not.associated(PREL_SC    )
      ALLOC_PBUP   = .not.associated(PBUP_SC    )
      ALLOC_WLCL   = .not.associated(WLCL_SC    )
      ALLOC_QTSRC  = .not.associated(QTSRC_SC   )
      ALLOC_THLSRC = .not.associated(THLSRC_SC  )
      ALLOC_THVLSRC= .not.associated(THVLSRC_SC )
      ALLOC_TKEAVG = .not.associated(TKEAVG_SC  )
      ALLOC_CLDTOP = .not.associated(CLDTOP_SC  )
      ALLOC_WUP    = .not.associated(WUP_SC     )
      ALLOC_QTUP   = .not.associated(QTUP_SC    )
      ALLOC_THLUP  = .not.associated(THLUP_SC   )
      ALLOC_THVUP  = .not.associated(THVUP_SC   )
      ALLOC_UUP    = .not.associated(UUP_SC     )
      ALLOC_VUP    = .not.associated(VUP_SC     )
      ALLOC_ENTR   = .not.associated(ENTR_SC    )
      ALLOC_DETR   = .not.associated(DETR_SC    )
      ALLOC_XC     = .not.associated(XC_SC      )
      ALLOC_QLDET  = .not.associated(QLDET_SC   )
      ALLOC_QIDET  = .not.associated(QIDET_SC   )
      ALLOC_QLENT  = .not.associated(QLENT_SC   )
      ALLOC_QIENT  = .not.associated(QIENT_SC   )
      ALLOC_QLSUB  = .not.associated(QLSUB_SC   )
      ALLOC_QISUB  = .not.associated(QISUB_SC   )
      ALLOC_NDROP  = .not.associated(SC_NDROP   )
      ALLOC_NICE   = .not.associated(SC_NICE    )

      if(ALLOC_UMF)    allocate( UMF_SC(IM,JM,0:LM) )
      if(ALLOC_MFD)    allocate( MFD_SC(IM,JM,LM) )
      if(ALLOC_DQVDT)  allocate( DQVDT_SC(IM,JM,LM) )
      if(ALLOC_DQLDT)  allocate( DQLDT_SC(IM,JM,LM) )
      if(ALLOC_DQIDT)  allocate( DQIDT_SC(IM,JM,LM) )
      if(ALLOC_DTHDT)  allocate( DTHDT_SC(IM,JM,LM) )
      if(ALLOC_DUDT)   allocate( DUDT_SC(IM,JM,LM)  )
      if(ALLOC_DVDT)   allocate( DVDT_SC(IM,JM,LM)  )
      if(ALLOC_DQRDT)  allocate( DQRDT_SC(IM,JM,LM) )
      if(ALLOC_DQSDT)  allocate( DQSDT_SC(IM,JM,LM) )
      if(ALLOC_CUFRC)  allocate( CUFRC_SC(IM,JM,LM) )
      if(ALLOC_QCU)    allocate( QCU_SC(IM,JM,LM)   )
      if(ALLOC_QLU)    allocate( QLU_SC(IM,JM,LM)   )
      if(ALLOC_QIU)    allocate( QIU_SC(IM,JM,LM)   )
      if(ALLOC_CBMF)   allocate( CBMF_SC(IM,JM)     )
      if(ALLOC_DQCDT)  allocate( DQCDT_SC(IM,JM,LM) )
      if(ALLOC_CNT)    allocate( CNT_SC(IM,JM)      )
      if(ALLOC_CNB)    allocate( CNB_SC(IM,JM)      )
      if(ALLOC_CIN)    allocate( CIN_SC(IM,JM)      )
      if(ALLOC_PLCL)   allocate( PLCL_SC(IM,JM)     )
      if(ALLOC_PLFC)   allocate( PLFC_SC(IM,JM)     )
      if(ALLOC_PINV)   allocate( PINV_SC(IM,JM)     )
      if(ALLOC_PREL)   allocate( PREL_SC(IM,JM)     )
      if(ALLOC_PBUP)   allocate( PBUP_SC(IM,JM)     )
      if(ALLOC_WLCL)   allocate( WLCL_SC(IM,JM)     )
      if(ALLOC_QTSRC)  allocate( QTSRC_SC(IM,JM)    )
      if(ALLOC_THLSRC) allocate( THLSRC_SC(IM,JM)   )
      if(ALLOC_THVLSRC)allocate( THVLSRC_SC(IM,JM)  )
      if(ALLOC_TKEAVG) allocate( TKEAVG_SC(IM,JM)   )
      if(ALLOC_CLDTOP) allocate( CLDTOP_SC(IM,JM)   )
      if(ALLOC_WUP)    allocate( WUP_SC(IM,JM,0:LM) )
      if(ALLOC_QTUP)   allocate( QTUP_SC(IM,JM,0:LM))
      if(ALLOC_THLUP)  allocate( THLUP_SC(IM,JM,0:LM))
      if(ALLOC_THVUP ) allocate( THVUP_SC(IM,JM,0:LM))
      if(ALLOC_UUP)    allocate( UUP_SC(IM,JM,0:LM) )
      if(ALLOC_VUP)    allocate( VUP_SC(IM,JM,0:LM) )
      if(ALLOC_ENTR)   allocate( ENTR_SC(IM,JM,LM)  )
      if(ALLOC_DETR)   allocate( DETR_SC(IM,JM,LM)  )
      if(ALLOC_XC)     allocate( XC_SC(IM,JM,LM)    )
      if(ALLOC_QLDET)  allocate( QLDET_SC(IM,JM,LM) )
      if(ALLOC_QIDET)  allocate( QIDET_SC(IM,JM,LM) )
      if(ALLOC_QLENT)  allocate( QLENT_SC(IM,JM,LM) )
      if(ALLOC_QIENT)  allocate( QIENT_SC(IM,JM,LM) )
      if(ALLOC_QLSUB)  allocate( QLSUB_SC(IM,JM,LM) )
      if(ALLOC_QISUB)  allocate( QISUB_SC(IM,JM,LM) )
      if(ALLOC_NDROP)  allocate( SC_NDROP(IM,JM,LM) )
      if(ALLOC_NICE)   allocate( SC_NICE(IM,JM,LM)  )

      IDIM = IM*JM
      IRUN = IM*JM


    ! define some default effective radii
      CLDREFFR = 10.0e-6
      CLDREFFS = 90.0e-6
      CLDREFFG = 90.0e-6
      CLDREFFI = 25.0e-6
      CLDREFFL = 10.0e-6


      !  Copy incoming state vars to local arrays that will be adjusted
      !  by physics.  Untouched state vars will later be used for 
      !  post facto tendency calculations.
      !----------------------------------------------------------------

      TH1      = TH
      Q1       = Q
      U1       = U
      V1       = V
       
      CNV_HAIL = 0.0
      CNV_PLE  = PLE*.01
      PLO      = 0.5*(CNV_PLE(:,:,0:LM-1) +  CNV_PLE(:,:,1:LM  ) )
      PKE      = (      PLE/MAPL_P00)**(MAPL_KAPPA)
      PK       = (100.0*PLO/MAPL_P00)**(MAPL_KAPPA)
      DP       = ( PLE(:,:,1:LM)-PLE(:,:,0:LM-1) )
      MASS     = DP/MAPL_GRAV
      iMASS    = 1.0/MASS

      TEMP     = TH1*PK

      DQS      = GEOS_DQSAT(TEMP, PLO, qsat=QSS)

      QSSFC    = GEOS_QSAT( TS , CNV_PLE(:,:,LM) )

!--srf-gf-scheme
      if(associated(DTDTCN)) DTDTCN=TEMP
!--srf-gf-scheme

!-srf - placed here before the cumulus parameterizations are called.
      if(associated(CNV_DQLDT)) CNV_DQLDT =  0.0
      ZLE(:,:,LM) = 0.
      do L=LM,1,-1
         ZLE(:,:,L-1) = TH (:,:,L) * (1.+MAPL_VIREPS*Q(:,:,L))    ! This term is really THV
         ZLO(:,:,L  ) = ZLE(:,:,L) + (MAPL_CP/MAPL_GRAV)*( PKE(:,:,L)-PK (:,:,L  ) ) * ZLE(:,:,L-1)
         ZLE(:,:,L-1) = ZLO(:,:,L) + (MAPL_CP/MAPL_GRAV)*( PK (:,:,L)-PKE(:,:,L-1) ) * ZLE(:,:,L-1)
         DZET(:,:,L ) = ZLE(:,:,L-1) - ZLE(:,:,L)
      end do

      GZLE  = MAPL_GRAV * ZLE
      GZLO  = MAPL_GRAV * ZLO
      RASPRCP  = 0.

      TSFCAIR  = TEMP(:,:,LM) + (MAPL_GRAV/MAPL_CP)*ZLO(:,:,LM)
      DTS      = TS - TSFCAIR

      IF (ASSOCIATED(DTSX)) DTSX = DTS
      IF (ASSOCIATED(KHX )) KHX  = KH(:,:,0:LM-1)


      if (associated(W_DIAG)) W_DIAG = OMEGA*(ZLE(:,:,1:LM)-ZLE(:,:,0:LM-1))/(PLE(:,:,1:LM)-PLE(:,:,0:LM-1))

      IRAS       = nint(LONS*100)
      JRAS       = nint(LATS*100)

      RASPARAMS%CLDMICRO = 0.0

      ! Some export work prior to doing calculations
      !---------------------------------------------
      if(associated(CNV_MFC)) CNV_MFC(:,:,LM) = 0.
      if(associated(RH1    )) RH1     = Q1/QSS
      if(associated(TVQ0   )) TVQ0    = SUM( (  Q +  QLLS + QLCN + QILS + QICN + QSNOW + QRAIN + QGRAUPEL )*MASS , 3 )
      if(associated(TVE0   )) TVE0    = SUM( (  MAPL_CP*TEMP + MAPL_ALHL*Q           & 
           -  MAPL_ALHF*(QILS+QICN) )*MASS , 3 )
      if(associated(DCPTE  )) DCPTE   = SUM( MAPL_CP*TEMP*MASS , 3 )
      if(associated(DQDTCN )) DQDTCN  = Q1
      if(associated(DTHDTCN)) DTHDTCN = TH1

      if(associated(DQLDT  )) DQLDT   = QLLS + QLCN
      if(associated(DQIDT  )) DQIDT   = QILS + QICN
      
      if(associated(DNDDT  )) DNDDT   = NCPL
      if(associated(DNCDT  )) DNCDT   = NCPI
      

      if(associated(DTDT_MOIST)) DTDT_MOIST   = TEMP

      if(associated(DDF_RH1   )) allocate( DDF_RH1z(LM) )
      if(associated(DDF_RH2   )) allocate( DDF_RH2z(LM) )
      if(associated(DDF_MUPH  )) allocate( DDF_MUPHz(LM) )
      if(associated(DDF_BYNC  )) allocate( DDF_BYNCz(LM) )
      if(associated(DDF_QVC   )) allocate( DDF_QVCz(LM) )
      if(associated(DDF_TC    )) allocate( DDF_TCz(LM) )
      if(associated(DDF_DQDT  )) allocate( DDF_DQDTz(LM) )
      if(associated(DDF_DTDT  )) allocate( DDF_DTDTz(LM) )
      if(associated(DDF_MFC   )) allocate( DDF_MFCz(0:LM) )

      !Trajectory for Moist TLM/ADJ
      if(associated(TH_moist)) TH_moist = TH1 
      if(associated(Q_moist))  Q_moist  = Q1

      ALLOC_BYNCY     = .not.associated(BYNCY     )
      ALLOC_CAPE      = .not.associated(CAPE      )
      ALLOC_INHB      = .not.associated(INHB      )

      if(ALLOC_BYNCY    ) allocate(BYNCY    (IM,JM,LM))
      if(ALLOC_CAPE     ) allocate(CAPE     (IM,JM   ))
      if(ALLOC_INHB     ) allocate(INHB     (IM,JM   ))

      call BUOYANCY( TEMP, Q1, QSS, DQS, DZET, ZLO, BYNCY, CAPE, INHB)

    ! Determine Deep Convective Fraction for RAS and PROGNO_CLOUD
    !    CNV_FRACTION = 0.0  --->  Large-scale
    !    CNV_FRACTION = 1.0  --->  Deep-Convection
         CNV_FRACTION = 0.0 

     ! Find RH at 600mb level
      call VertInterp(RHat600,Q1/QSS,log(PLE),log(60000.),STATUS)
      VERIFY_(STATUS)
     ! Fill undefs (600mb below the surface) with surface RH values L=LM
         levs600  = max(1,count(PREF < 60000.))
      WHERE (RHat600 == MAPL_UNDEF)
         RHat600 = Q1(:,:,levs600)/QSS(:,:,levs600)
      END WHERE

     ! Find QV at 600mb level
      call VertInterp(QV600,Q,log(PLE),log(60000.),STATUS)
      VERIFY_(STATUS)
     ! Fill undefs (600mb below the surface) with surface QV values L=LM
         levs600  = max(1,count(PREF < 60000.))
      WHERE (QV600 == MAPL_UNDEF)
         QV600 = Q(:,:,levs600)
      END WHERE

    ! Compute the deep convective fraction based on mid-tropospheric moisture (QV at 600mb)
    !     mid-tropospheric moisture is used as an indicator of vertical motion
    !     associated with active deep convection lifting moisture to the mid-troposphere
    ! QV at 600mb Criteria 
    ! call MAPL_GetResource(STATE,CNV_FRACTION_MIN, 'CNV_FRACTION_MIN:', DEFAULT= 0.00250, RC=STATUS)
    ! VERIFY_(STATUS)
    ! call MAPL_GetResource(STATE,CNV_FRACTION_MAX, 'CNV_FRACTION_MAX:', DEFAULT= 0.00600, RC=STATUS)
    ! VERIFY_(STATUS)

    ! CAPE Criteria
      if( LM .ne. 72 ) then
        call MAPL_GetResource(STATE,CNV_FRACTION_MIN, 'CNV_FRACTION_MIN:', DEFAULT=    0.0, RC=STATUS)
        VERIFY_(STATUS)
        call MAPL_GetResource(STATE,CNV_FRACTION_MAX, 'CNV_FRACTION_MAX:', DEFAULT= 1500.0, RC=STATUS)
        VERIFY_(STATUS)
        call MAPL_GetResource(STATE,GF_MIN_AREA, 'GF_MIN_AREA:', DEFAULT= 1.e6, RC=STATUS)
        VERIFY_(STATUS)
        call MAPL_GetResource(STATE,STOCHASTIC_CNV, 'STOCHASTIC_CNV:', DEFAULT= 1, RC=STATUS)
        VERIFY_(STATUS)
      else
        call MAPL_GetResource(STATE,CNV_FRACTION_MIN, 'CNV_FRACTION_MIN:', DEFAULT=  500.0, RC=STATUS)
        VERIFY_(STATUS)
        call MAPL_GetResource(STATE,CNV_FRACTION_MAX, 'CNV_FRACTION_MAX:', DEFAULT= 1500.0, RC=STATUS)
        VERIFY_(STATUS)
        call MAPL_GetResource(STATE,GF_MIN_AREA, 'GF_MIN_AREA:', DEFAULT= 1.e6, RC=STATUS)
        VERIFY_(STATUS)
        call MAPL_GetResource(STATE,STOCHASTIC_CNV, 'STOCHASTIC_CNV:', DEFAULT= 0, RC=STATUS)
        VERIFY_(STATUS)
      endif

      if( CNV_FRACTION_MAX > CNV_FRACTION_MIN ) then
         if (CNV_FRACTION_MAX < 1.0) then
          ! QV at 600mb
           DO J=1, JM
             DO I=1, IM
               CNV_FRACTION(I,J) =MAX(0.0,MIN(1.0,(QV600(I,J)-CNV_FRACTION_MIN)/(CNV_FRACTION_MAX-CNV_FRACTION_MIN)))
             END DO
           END DO
         else
          ! CAPE
           WHERE (CAPE .ne. MAPL_UNDEF)   
              CNV_FRACTION =MAX(0.0,MIN(1.0,(CAPE-CNV_FRACTION_MIN)/(CNV_FRACTION_MAX-CNV_FRACTION_MIN)))
           END WHERE
         endif
      endif

      if(associated(CNV_FRC )) CNV_FRC  = CNV_FRACTION
      if(associated(Q600    )) Q600     = QV600
      if(associated(RH600   )) RH600    = RHat600

      K0 = LM
      ICMIN    = max(1,count(PREF < PMIN_DET))
      KCBLMIN  =       count(PREF < PMIN_CBL)
      levs925  = max(1,count(PREF < 92500.))

      ! MAT This line kills zero-diff with Icarus-3_2_p9 when
      !     bootstrapping turbulence (at least). Either GF or
      !     shallow must need this, but this might not be the
      !     place for it
      !
      !     Also, KPBLIN is always associated so that
      !     is unnessary.
      !if(associated(KPBLIN)) KPBLIN = MAX(KPBLIN, 1.000001)
      DO J=1, JM
         DO I=1, IM
            if (KPBLIN(I,J) == 0) KPBLIN(I,J) = LM-1
         END DO
      END DO

      ! MAT The code below should use nint(KPBLIN)
      DO J=1, JM
         DO I=1, IM
            ZPBL(I, J) =ZLE(I, J,  NINT(KPBLIN(I, J )))
         END DO
      END DO

#ifdef DONT_SKIP_cloud_ptr_stubs
    if (.false.) then 
       call cloud_ptr_stubs (SMAXL, SMAXI, WSUB, CCN01, CCN04, CCN1, &
            NHET_NUC, NLIM_NUC, SO4, ORG, BCARBON, &
            DUST, SEASALT, NCPL_VOL, NCPI_VOL, NRAIN, NSNOW, &
            CDNC_NUC, INC_NUC, SAT_RAT, QSTOT, QRTOT, CLDREFFS, CLDREFFR, & 
            DQVDT_micro,DQIDT_micro, DQLDT_micro, DTDT_micro, RL_MASK, RI_MASK, &
            KAPPA, SC_ICE, CFICE, CFLIQ, RHICE, RHLIQ,  &
            RAD_CF, RAD_QL, RAD_QI, RAD_QS, RAD_QR, RAD_QV, &
            CLDREFFI, CLDREFFL, NHET_IMM, NHET_DEP, & 
            DUST_IMM, DUST_DEP,  SCF, SCF_ALL, &
            SIGW_GW, SIGW_CNV, &
            SIGW_TURB, SIGW_RC, RHCmicro, DNHET_IMM, &
            BERG, BERGS, MELT, DNHET_CT, DTDT_macro, QCRES, DT_RASP, FRZPP_LS, &
            SNOWMELT_LS, QIRES, AUTICE, PFRZ,  DNCNUC, DNCHMSPLIT, DNCSUBL, &
            DNCAUTICE, DNCACRIS, DNDCCN, DNDACRLS, DNDEVAPC, DNDACRLR, DNDAUTLIQ, &
            DNDCNV, DNCCNV) 
    end if 
#endif

!--kml
      if(adjustl(CLDMICRO)=="2MOMENT" .or. USE_AEROSOL_NN) then
!--kml

      call MAPL_TimerOn(STATE,"--USE_AEROSOL_NN1")

      do k = 1, LM
          do j = 1, JM
              do i = 1, IM
                  AeroProps(i,j,k)%num = 0.0
              end do
          end do
      end do    


      call ESMF_StateGet(IMPORT, 'AERO', aero_aci, __RC__)

!      call ESMF_AttributeGet(aero_aci, name='implements_aerosol_activation_properties_method', &
!                                       value=implements_aerosol_activation_properties, __RC__)

!      if (implements_aerosol_activation_properties) then

          call ESMF_AttributeGet(aero_aci, name='number_of_aerosol_modes', value=n_modes, __RC__)

          if (n_modes > 0) then

              allocate(aero_aci_modes(n_modes), __STAT__)
              call ESMF_AttributeGet(aero_aci, name='aerosol_modes', itemcount=n_modes, valuelist=aero_aci_modes, __RC__)

              call ESMF_AttributeGet(aero_aci, name='air_pressure_for_aerosol_optics', value=aci_field_name, __RC__)
              if (aci_field_name /= '') then
                  call MAPL_GetPointer(aero_aci, aci_ptr_3d, trim(aci_field_name), __RC__)
                  aci_ptr_3d = PLE
              end if

              call ESMF_AttributeGet(aero_aci, name='air_temperature', value=aci_field_name, __RC__)
              if (aci_field_name /= '') then
                  call MAPL_GetPointer(aero_aci, aci_ptr_3d, trim(aci_field_name), __RC__)
                  aci_ptr_3d = TEMP
              end if    

              call ESMF_AttributeGet(aero_aci, name='fraction_of_land_type', value=aci_field_name, __RC__)
              if (aci_field_name /= '') then
                  call MAPL_GetPointer(aero_aci, aci_ptr_2d, trim(aci_field_name), __RC__)
                  aci_ptr_2d = FRLAND
              end if


              call MAPL_GetResource(STATE,USE_MOIST_BUFFER, 'USE_MOIST_BUFFER:', DEFAULT=.TRUE., RC=STATUS)
              if (USE_MOIST_BUFFER) then
                 allocate(buffer(im,jm,lm,n_modes,8), __STAT__)
              end if

              ACTIVATION_PROPERTIES: do n = 1, n_modes
                 call ESMF_AttributeSet(aero_aci, name='aerosol_mode', value=trim(aero_aci_modes(n)), __RC__)
                 
                 ! execute the aerosol activation properties method 
                 call ESMF_MethodExecute(aero_aci, label='aerosol_activation_properties', userRC=ACI_STATUS, RC=STATUS)
                 VERIFY_(ACI_STATUS)
                 VERIFY_(STATUS)

                 ! copy out aerosol activation properties
                 call ESMF_AttributeGet(aero_aci, name='aerosol_number_concentration', value=aci_field_name, __RC__)
                 call MAPL_GetPointer(aero_aci, aci_num, trim(aci_field_name), __RC__)

                 call ESMF_AttributeGet(aero_aci, name='aerosol_dry_size', value=aci_field_name, __RC__)
                 call MAPL_GetPointer(aero_aci, aci_dgn, trim(aci_field_name), __RC__)

                 call ESMF_AttributeGet(aero_aci, name='width_of_aerosol_mode', value=aci_field_name, __RC__)
                 call MAPL_GetPointer(aero_aci, aci_sigma, trim(aci_field_name), __RC__)

                 call ESMF_AttributeGet(aero_aci, name='aerosol_density', value=aci_field_name, __RC__)
                 call MAPL_GetPointer(aero_aci, aci_density, trim(aci_field_name), __RC__)

                 call ESMF_AttributeGet(aero_aci, name='aerosol_hygroscopicity', value=aci_field_name, __RC__)
                 call MAPL_GetPointer(aero_aci, aci_hygroscopicity, trim(aci_field_name), __RC__)

                 call ESMF_AttributeGet(aero_aci, name='fraction_of_dust_aerosol', value=aci_field_name, __RC__)
                 call MAPL_GetPointer(aero_aci, aci_f_dust, trim(aci_field_name), __RC__)

                 call ESMF_AttributeGet(aero_aci, name='fraction_of_soot_aerosol', value=aci_field_name, __RC__)
                 call MAPL_GetPointer(aero_aci, aci_f_soot, trim(aci_field_name), __RC__)

                 call ESMF_AttributeGet(aero_aci, name='fraction_of_organic_aerosol', value=aci_field_name, __RC__)
                 call MAPL_GetPointer(aero_aci, aci_f_organic, trim(aci_field_name), __RC__)

#if (0)
                 if (MAPL_AM_I_ROOT()) then
                    print *
                    print *, 'AERO::' // trim(aero_aci_modes(n))

                    print *, 'num            : ', aci_num(1,1,LM)
                    print *, 'dgn            : ', aci_dgn(1,1,LM)
                    print *, 'sigma          : ', aci_sigma(1,1,LM)
                    print *, 'hygroscopicity : ', aci_hygroscopicity(1,1,LM)
                    print *, 'density        : ', aci_density(1,1,LM)
                    print *, 'f_dust         : ', aci_f_dust(1,1,LM)
                    print *, 'f_soot         : ', aci_f_soot(1,1,LM)
                    print *, 'f_organic      : ', aci_f_organic(1,1,LM)
                    print *
                 END IF
#endif

                 if (USE_MOIST_BUFFER) then
                    buffer(:,:,:,n,1) = aci_num
                    buffer(:,:,:,n,2) = aci_dgn
                    buffer(:,:,:,n,3) = aci_sigma
                    buffer(:,:,:,n,4) = aci_hygroscopicity
                    buffer(:,:,:,n,5) = aci_density
                    buffer(:,:,:,n,6) = aci_f_dust
                    buffer(:,:,:,n,7) = aci_f_soot
                    buffer(:,:,:,n,8) = aci_f_organic
                 else
                    AeroProps(:,:,:)%num(n)   = aci_num
                    AeroProps(:,:,:)%dpg(n)   = aci_dgn
                    AeroProps(:,:,:)%sig(n)   = aci_sigma
                    AeroProps(:,:,:)%kap(n)   = aci_hygroscopicity
                    AeroProps(:,:,:)%den(n)   = aci_density
                    AeroProps(:,:,:)%fdust(n) = aci_f_dust
                    AeroProps(:,:,:)%fsoot(n) = aci_f_soot
                    AeroProps(:,:,:)%forg(n)  = aci_f_organic
                    AeroProps(:,:,:)%nmods    = n_modes                 ! no need of a 3D field: aero provider specific
                 end if

              end do ACTIVATION_PROPERTIES

              if (USE_MOIST_BUFFER) then
                 do k = 1, LM
                    do j = 1, JM
                       do i = 1, IM
                          do n = 1, n_modes
                             AeroProps(i,j,k)%num(n)   = buffer(i,j,k,n,1)
                             AeroProps(i,j,k)%dpg(n)   = buffer(i,j,k,n,2)
                             AeroProps(i,j,k)%sig(n)   = buffer(i,j,k,n,3)
                             AeroProps(i,j,k)%kap(n)   = buffer(i,j,k,n,4)
                             AeroProps(i,j,k)%den(n)   = buffer(i,j,k,n,5)
                             AeroProps(i,j,k)%fdust(n) = buffer(i,j,k,n,6)
                             AeroProps(i,j,k)%fsoot(n) = buffer(i,j,k,n,7)
                             AeroProps(i,j,k)%forg(n)  = buffer(i,j,k,n,8)
                          end do
                          AeroProps(i,j,k)%nmods       = n_modes                 ! no need of a 3D field: aero provider specific
                       end do
                    end do
                 end do

                 deallocate(buffer, __STAT__)
              end if

              deallocate(aero_aci_modes, __STAT__)
          end if

      call MAPL_TimerOff(STATE,"--USE_AEROSOL_NN1")

!      else
          ! options: 
          !     *) set aerosol concentrations to 0.0, i.e., no aerosol
          !     *) raise an exception if aerosol is required!
!      end if


      call init_Aer(AeroAux)
      call init_Aer(AeroAux_b)

!--kml
      if((adjustl(CLDMICRO)=="1MOMENT" .or. adjustl(CLDMICRO)=="GFDL") .and. USE_AEROSOL_NN) GOTO 31416
!--kml

      call MAPL_TimerOn(STATE,"--USE_AEROSOL_NN2")

         !initialize some values 
         CNV_PRC3 =0.0    
         CNV_UPDF =0.0    
         CNV_FICE =0.0    
         CNV_NICE =0.0
         CNV_NDROP=0.0 !(1e-6 cm-3)
         RHCmicro = 0.0            
         AUTOC_CN_OCN  = AUTO_CNV
         AUTOC_CN_LAND = AUTO_CNV
         CNV_CVW  = 0.0
         ENTLAM   = 0.0
         RASPARAMS%CLDMICRO = 1.0
	 RASPARAMS%FDROP_DUST = FDROP_DUST
	 RASPARAMS%FDROP_SOOT = FDROP_SOOT
     RASPARAMS%FDROP_SEASALT = SS_INFAC
         SHLW_PRC3 = 0.0
         SHLW_SNO3 = 0.0
         UFRC_SC   = 0.0	 

         wparc_rc = 0.0

         ! find the minimun level for cloud micro calculations
	 call find_l(KMIN_TROP, PLO, pmin_trop, IM, JM, LM, 10, LM-2) 
	 

      call MAPL_TimerOff(STATE,"--USE_AEROSOL_NN2")

      end if          !CLDMICRO
      !===================================================================================
!--kml
31416 CONTINUE
!--kml

      call MAPL_TimerOff(STATE,"-MISC1")

      ! Do convection on IDIM soundings
      !-------------------------------------------------

      call MAPL_TimerOn (STATE,"-PRE_RAS")

      ! Find Cloud base level (LCB), strap if desired
      ! and reset K0 accordingly
      !-----------------------------------------------



      ! Copy tracers to local array
      !----------------------------

      KK=0
      do K=1,KM
         if(IS_FRIENDLY(K)) then
            KK = KK+1
            !PRINT *, "*******TESTING: QNAME, FSCAV_, FSCAV********"
            FSCAV(KK) = FSCAV_(K)
            !PRINT *, QNAMES(K), FSCAV_(K), FSCAV(KK)
            XHO(:,:,:,KK) = TRPtrs(K)%Q(:,:,:)
         end if
      end do

      ! Determine how to do cloud base 
      !-------------------------------

      KLCL = FINDLCL( TH1, Q1, PLO, PK, IM, JM, LM )
      KLFC = FINDLFC( BYNCY, IM, JM, LM )
      !!    KPBL = FINDPBL( KH, IM, JM, LM )
      !! Set subcloud layer height to one level below PBL height level
      !!   make sure subcloud layer is at least 2 levels thick
      do j = 1,jm
         do i = 1,im
            if(nint(kpblin(i,j)).ne.0) then
               kpbl(i,j) = max(min(nint(kpblin(i,j))+1,LM-1), 1)
            else
               kpbl(i,j) = LM-1
            endif
         enddo
      enddo

      SIGE = PREF/PREF(LM) ! this should eventually change

      do J=1,JM
         do I=1,IM

            SELECT CASE( CBL_METHOD )

            CASE( 1 )
               KCBL(I,J)   =  K0 - KSTRAP
               KCBL(I,J)   =  MAX( KCBL(I,J), KCBLMIN )
               WGT0(I,J,:)            = 0.
               WGT0(I,J,KCBL(I,J):K0) = 1.0 
               WGT1(I,J,:)            = 0.
               WGT1(I,J,KCBL(I,J):K0) = 1.0 

            CASE( 2 )
               KCBL(I,J)   = KLCL(I,J)
               KCBL(I,J)   =  MAX( KCBL(I,J), KCBLMIN )
               WGT0(I,J,:)            = 0.
               WGT0(I,J,KCBL(I,J):K0) = 1.0 
               WGT1(I,J,:)            = 0.
               WGT1(I,J,KCBL(I,J):K0) = 1.0 

            CASE ( 3 )
               KCBL(I,J)   = KPBL(I,J)
               KCBL(I,J)   =  MAX( KCBL(I,J), KCBLMIN )
               WGT0(I,J,:)            = 0.
               WGT0(I,J,KCBL(I,J)   ) = 1.0 
               WGT1(I,J,:)            = 0.
               WGT1(I,J,KCBL(I,J)   ) = 1.0 

            CASE ( 4 )
               KCBL(I,J)   = KLCL(I,J)
               KCBL(I,J)   =  MAX( KCBL(I,J), KCBLMIN )
               WGT0(I,J,:)            = 0.
               WGT0(I,J,KCBL(I,J)   ) = 1.0 
               WGT1(I,J,:)            = 0.
               WGT1(I,J,KCBL(I,J)   ) = 1.0 
               do kii = kcbl(i,j)+1,k0
                  WGT0(I,J,kii  )      = wgt1(I,J,kii-1)*0.6
                  WGT1(I,J,kii  )      = wgt1(I,J,kii-1)*0.6  
               end do

            CASE ( 5 )
               KCBL(I,J)   = KPBL(I,J)
               KCBL(I,J)   =  MAX( KCBL(I,J), KCBLMIN )
               WGT0(I,J,:)            = 0.
               WGT0(I,J,KCBL(I,J)   ) = 1.0 
               WGT1(I,J,:)            = 0.
               WGT1(I,J,KCBL(I,J)   ) = 1.0 
               do kii = kcbl(i,j)+1,k0
                  WGT0(i,j,kii  )      = wgt1(i,j,kii-1)*0.6
                  WGT1(i,j,kii  )      = wgt1(i,j,kii-1)*0.6  
               end do

            CASE( 6 )
               KCBL(I,J)   = KPBL(I,J)
               KCBL(I,J)   =  MAX( KCBL(I,J), KCBLMIN )
               WGT0(I,J,:)            = 0.
               WGT0(I,J,KCBL(I,J):K0) = 1.0 
               WGT1(I,J,:)            = 0.
               WGT1(I,J,KCBL(I,J):K0) = 1.0 

            END SELECT

            ZCBLx(I,J) = ZLE( I, J, KCBL(I,J)-1 )

            if(associated(ZLCL)) then
               if(KLCL(I,J)>0) then
                  ZLCL(I,J) = ZLE(I,J,KLCL(I,J)-1)
               else
                  ZLCL(I,J) = MAPL_UNDEF
               end if
            end if

            if(associated(ZLFC)) then
               if(KLFC(I,J)>0) then
                  ZLFC(I,J) = ZLE(I,J,KLFC(I,J)-1)
               else
                  ZLFC(I,J) = MAPL_UNDEF
               end if
            end if


         end do
      end do
     
      if(associated(QCBL)) then
        do J=1,JM
           do I=1,IM
              QVCBL(I,J) = Q(I,J,KCBL(I,J)-1)
           end do
        end do
        QCBL = QVCBL
      end if
     !Option to Use Q at CBL to adjust convective intensities based Q at CBL
     !if( CNV_FRACTION_MAX > CNV_FRACTION_MIN ) then
     !   if (CNV_FRACTION_MAX < 1.0) then ! QV at CBL
     !     DO J=1, JM
     !       DO I=1, IM
     !         CNV_FRACTION(I,J) = MAX(0.0,MIN(1.0,(QVCBL(I,J)-CNV_FRACTION_MIN)/(CNV_FRACTION_MAX-CNV_FRACTION_MIN)))
     !       END DO
     !     END DO
     !     if(associated(CNV_FRC )) CNV_FRC = CNV_FRACTION
     !   endif
     !endif

      if (ADJUSTL(CONVPAR_OPTION) == "RAS") then 
       RASAL1 = RASPARAMS%RASAL1
       RASAL2 = RASPARAMS%RASAL2

       if (RASAL2 > 0.0) then
         RASAL2_2d(:,:) = RASAL2
       else
         ! include CNV dependence
         DO J=1, JM
            DO I=1, IM
            RASAL2_2d(I,J) = CNV_FRACTION(I,J)*ABS(RASAL2) + (1-CNV_FRACTION(I,J))*RASAL1
            END DO
         END DO
       endif

       ! Cheat by adding a kick to CB temp and q
       ! ---------------------------------------
       TPERT  = ABS(CBL_TPERT) * ( TS - ( TEMP(:,:,LM)+ MAPL_GRAV*ZLO(:,:,LM)/MAPL_CP )  ) 
       if (CBL_TPERT < 0) then
        ! Make TPERT 0 in areas of deep convection
         TPERT = TPERT*(1.0-CNV_FRACTION)
       endif
       QPERT  = CBL_QPERT * ( QSSFC - Q(:,:,LM) )          
       TPERT  = MAX( TPERT , 0.0 )
       QPERT  = MAX( QPERT , 0.0 )

       where (FRLAND<0.1) 
          TPERT = MIN( TPERT , CBL_TPERT_MXOCN ) ! ocean
       elsewhere
          TPERT = MIN( TPERT , CBL_TPERT_MXLND ) ! land
       end where
      end if

       ! Myong-In I just   
       ! put these 100s    
       ! in. Think about  
       ! whether you want           |                         |
       ! to keep them               V                         V
       SEEDRAS(:,:,1) = 1000000 * ( 100*TEMP(:,:,LM)   - INT( 100*TEMP(:,:,LM) ) )
       SEEDRAS(:,:,2) = 1000000 * ( 100*TEMP(:,:,LM-1) - INT( 100*TEMP(:,:,LM-1) ) )

      if (STOCHASTIC_CNV /= 0) then
      ! Create bit-processor-reproducible random white noise for convection [0:1]
       SEEDCNV(:,:)   = SEEDRAS(:,:,1)/1000000.0
       where (SEEDCNV > 1.0)
          SEEDCNV = 1.0
       end where
       where (SEEDCNV < 0.0)
          SEEDCNV = 0.0
       end where 
       SEEDCNV = SEEDCNV*(1.875-0.5)+0.5
      else
       SEEDCNV(:,:) = 1.0
      endif

       CALL MAPL_GetPointer(EXPORT, STOCH_CNV,  'STOCH_CNV', RC=STATUS)
       VERIFY_(STATUS)
       if (associated(STOCH_CNV)) STOCH_CNV = SEEDCNV

      if(adjustl(CLDMICRO)=="2MOMENT") then
       if (NPRE_FRAC > 0.0) then
         NPRE_FRAC_2d(:,:) = NPRE_FRAC
       else
         ! include CNV_FRACTION dependence
         DO J=1, JM
            DO I=1, IM
            NPRE_FRAC_2d(I,J) = CNV_FRACTION(I,J)*ABS(NPRE_FRAC) + (1-CNV_FRACTION(I,J))*0.05
            END DO
         END DO
       endif
      end if

      ! Compute initial mass loading for aerosols; CAR 12/19/08
      ! -------------------------------------------------------
      !! First initialize everything to zero
      if(associated(DDU2gDT))   DDU2gDT =  0.0
      if(associated(DSS2gDT))   DSS2gDT =  0.0
      if(associated(DBC2gDT))   DBC2gDT =  0.0
      if(associated(DOC2gDT))   DOC2gDT =  0.0
      if(associated(DSU2gDT))   DSU2gDT =  0.0
      if(associated(DNI2gDT))   DNI2gDT =  0.0
      if(associated(DNH4A2gDT)) DNH4A2gDT =  0.0
      if(associated(DNH32gDT))  DNH32gDT =  0.0
      if(associated(DBRC2gDT))  DBRC2gDT=  0.0

      if(associated(DDUDT))   DDUDT =  0.0
      if(associated(DSSDT))   DSSDT =  0.0
      if(associated(DBCDT))   DBCDT =  0.0
      if(associated(DOCDT))   DOCDT =  0.0
      if(associated(DSUDT))   DSUDT =  0.0
      if(associated(DNIDT))   DNIDT =  0.0
      if(associated(DNH4ADT)) DNH4ADT =  0.0
      if(associated(DNH3DT))  DNH3DT =  0.0
      if(associated(DBRCDT))  DBRCDT=  0.0
      if(associated(DDUDTcarma)) DDUDTcarma =  0.0
      if(associated(DSSDTcarma)) DSSDTcarma =  0.0

      CMDU2g   = 0.0
      CMSS2g   = 0.0
      CMOC2g   = 0.0
      CMBC2g   = 0.0
      CMSU2g   = 0.0
      CMNI2g   = 0.0
      CMNH4A2g = 0.0
      CMNH32g  = 0.0
      CMBRC2g  = 0.0

      CMDU   = 0.0
      CMSS   = 0.0
      CMOC   = 0.0
      CMBC   = 0.0
      CMSU   = 0.0
      CMNI   = 0.0
      CMNH4A = 0.0
      CMNH3  = 0.0
      CMBRC  = 0.0
      CMDUcarma = 0.0
      CMSScarma = 0.0

      !! Now loop over tracers and accumulate initial column loading
      !! tendency  kg/m2/s CAR
!if(mapl_am_i_root()) print*,'MOIST CNAMES = ',CNAMES
      KK=0
      do K=1,KM
         if(IS_FRIENDLY(K)) then
            KK = KK + 1
            QNAME = trim(QNAMES(K))
            CNAME = trim(CNAMES(K))
            if((CNAME == 'DU') .or. (CNAME == 'SS') .or. (CNAME == 'NI') .or. (CNAME == 'SU') .or. &
               (CNAME == 'CA.oc') .or. (CNAME == 'CA.bc') .or. (CNAME == 'CA.br')) then   ! Diagnostics for GOCART2G tracers
               SELECT CASE (QNAME(1:3))
               CASE ('DU0')
                  if(associated(DDU2gDT)) then
                     CMDU2g = CMDU2g + sum(XHO(:,:,:,KK)*DP(:,:,:),dim=3)
                  end if
               CASE ('SS0')
                  if(associated(DSS2gDT)) then
                     CMSS2g = CMSS2g + sum(XHO(:,:,:,KK)*DP(:,:,:),dim=3)
                  end if
               CASE ('SO4')
                  if(associated(DSU2gDT)) then
                     CMSU2g = CMSU2g + sum(XHO(:,:,:,KK)*DP(:,:,:),dim=3)
                  end if
               CASE ('NO3')
                  if(associated(DNI2gDT)) then
                     CMNI2g = CMNI2g + sum(XHO(:,:,:,KK)*DP(:,:,:),dim=3)
                  end if
               CASE ('NH3')
                  if(associated(DNH32gDT)) then
                     CMNH32g = CMNH32g + sum(XHO(:,:,:,KK)*DP(:,:,:),dim=3)
                  end if
               CASE ('NH4')
                  if(associated(DNH4A2gDT)) then
                     CMNH4A2g = CMNH4A2g + sum(XHO(:,:,:,KK)*DP(:,:,:),dim=3)
                  end if
               END SELECT

               SELECT CASE (QNAME(1:13))
               CASE ('CA.bcphilic')
                  if(associated(DBC2gDT)) then
                     CMBC2g = CMBC2g + sum(XHO(:,:,:,KK)*DP(:,:,:),dim=3)
                  end if
!               CASE ('OCp')
               CASE ('CA.ocphilic')
                  if(associated(DOC2gDT)) then
                     CMOC2g = CMOC2g + sum(XHO(:,:,:,KK)*DP(:,:,:),dim=3)
                  end if
               CASE ('CA.brphilic')
                  if(associated(DBRC2gDT)) then
                     CMBRC2g = CMBRC2g + sum(XHO(:,:,:,KK)*DP(:,:,:),dim=3)
                  end if
               END SELECT
            endif

            if(CNAME == 'GOCART') then   ! Diagnostics for GOCART tracers
               SELECT CASE (QNAME(1:3))
               CASE ('du0')
                  if(associated(DDUDT)) then
                     CMDU = CMDU + sum(XHO(:,:,:,KK)*DP(:,:,:),dim=3)
                  end if
               CASE ('ss0')
                  if(associated(DSSDT)) then
                     CMSS = CMSS + sum(XHO(:,:,:,KK)*DP(:,:,:),dim=3)
                  end if
               CASE ('BCp')
                  if(associated(DBCDT)) then
                     CMBC = CMBC + sum(XHO(:,:,:,KK)*DP(:,:,:),dim=3)
                  end if
               CASE ('OCp')
                  if(associated(DOCDT)) then
                     CMOC = CMOC + sum(XHO(:,:,:,KK)*DP(:,:,:),dim=3)
                  end if
               CASE ('SO4')
                  if(associated(DSUDT)) then
                     CMSU = CMSU + sum(XHO(:,:,:,KK)*DP(:,:,:),dim=3)
                  end if
               CASE ('NO3')
                  if(associated(DNIDT)) then
                     CMNI = CMNI + sum(XHO(:,:,:,KK)*DP(:,:,:),dim=3)
                  end if
               CASE ('NH3')
                  if(associated(DNH3DT)) then
                     CMNH3 = CMNH3 + sum(XHO(:,:,:,KK)*DP(:,:,:),dim=3)
                  end if
               CASE ('NH4')
                  if(associated(DNH4ADT)) then
                     CMNH4A = CMNH4A + sum(XHO(:,:,:,KK)*DP(:,:,:),dim=3)
                  end if
               CASE ('BRC')
                  if(associated(DBRCDT)) then
                     CMBRC = CMBRC + sum(XHO(:,:,:,KK)*DP(:,:,:),dim=3)
                  end if
               END SELECT
            endif

            if(CNAME == 'CARMA') then   ! Diagnostics for CARMA tracers
               ! Check name to see if it is a "pc" element
               ENAME = ''
               ind= index(QNAME, '::')
               if (ind> 0) then
                  ENAME = trim(QNAME(ind+2:ind+3))  ! Component name (e.g., GOCART, CARMA)
                  if(ENAME == 'pc') then
                     SELECT CASE (QNAME(1:4))
                     CASE ('dust') ! CARMA DUST
                        if(associated(DDUDTcarma)) then
                           CMDUcarma = CMDUcarma + sum(XHO(:,:,:,KK)*DP(:,:,:),dim=3)
                        end if
                     CASE ('seas') ! CARMA SEASALT
                        if(associated(DSSDTcarma)) then
                           CMSScarma = CMSScarma + sum(XHO(:,:,:,KK)*DP(:,:,:),dim=3) 
                        end if
                     END SELECT
                  endif
               endif
            endif
         end if
      end do
      
      call MAPL_TimerOff(STATE,"-PRE_RAS")

!=====================================================================================
!-srf-gf-scheme
      IF(ADJUSTL(CONVPAR_OPTION) == 'GF' .or. ADJUSTL(CONVPAR_OPTION) == 'BOTH') THEN
     
         call MAPL_TimerOn (STATE,"-GF")
	 !- block shallow plume if UWShallow scheme is being applied
	 if(DOSHLW==1 .and.  icumulus_gf(shal)==1) print*,"Moist is running with UWshallow and GFshallow"
         !
         !-initialize/reset output arrays 
         CNV_MF0  =0.0 ! 'cloud_base_mass_flux'              - 'kg m-2 s-1'
         CNV_MFD  =0.0 ! 'detraining_mass_flux',             - 'kg m-2 s-1'    
         CNV_MFC  =0.0 ! 'cumulative_mass_flux',             - 'kg m-2 s-1'   
         CNV_CVW  =0.0 ! 'updraft_vertical_velocity',        - 'hPa s-1', 
         CNV_QC   =0.0 ! 'grid_mean_convective_condensate',  - 'kg kg-1',   
         CNV_UPDF =0.0 ! 'updraft_areal_fraction',           - '1', 
         ENTLAM   =0.0 ! 'entrainment parameter',            - 'm-1',  
         CNV_PRC3 =0.0 ! 'convective_precipitation           - 'kg m-2 s-1'
    
         call MAPL_GetPointer(IMPORT, USTAR    ,'USTAR'   ,RC=STATUS); VERIFY_(STATUS)
         call MAPL_GetPointer(IMPORT, TSTAR    ,'TSTAR'   ,RC=STATUS); VERIFY_(STATUS)
         call MAPL_GetPointer(IMPORT, QSTAR    ,'QSTAR'   ,RC=STATUS); VERIFY_(STATUS)
         call MAPL_GetPointer(IMPORT, T2M      ,'T2M  '   ,RC=STATUS); VERIFY_(STATUS)
         call MAPL_GetPointer(IMPORT, Q2M      ,'Q2M  '   ,RC=STATUS); VERIFY_(STATUS)
         call MAPL_GetPointer(IMPORT, TA       ,'TA   '   ,RC=STATUS); VERIFY_(STATUS) 
         call MAPL_GetPointer(IMPORT, QA       ,'QA   '   ,RC=STATUS); VERIFY_(STATUS) 
         call MAPL_GetPointer(IMPORT, PHIS     ,'PHIS '   ,RC=STATUS); VERIFY_(STATUS)
         call MAPL_GetPointer(IMPORT, DQVDTDYN ,'DQVDTDYN' ,RC=STATUS); VERIFY_(STATUS)
         call MAPL_GetPointer(IMPORT, DTDTDYN  ,'DTDTDYN'  ,RC=STATUS); VERIFY_(STATUS)
         call MAPL_GetPointer(IMPORT, QV_DYN_IN,'QV_DYN_IN',RC=STATUS); VERIFY_(STATUS)
         call MAPL_GetPointer(IMPORT, U_DYN_IN ,'U_DYN_IN' ,RC=STATUS); VERIFY_(STATUS)
         call MAPL_GetPointer(IMPORT, V_DYN_IN ,'V_DYN_IN' ,RC=STATUS); VERIFY_(STATUS)
         call MAPL_GetPointer(IMPORT, T_DYN_IN ,'T_DYN_IN' ,RC=STATUS); VERIFY_(STATUS)
         call MAPL_GetPointer(IMPORT, PLE_DYN_IN,'PLE_DYN_IN',RC=STATUS); VERIFY_(STATUS)
         call MAPL_GetPointer(IMPORT, DTDT_BL  ,'DTDT_BL'  ,RC=STATUS); VERIFY_(STATUS)
         call MAPL_GetPointer(IMPORT, DQDT_BL  ,'DQDT_BL'  ,RC=STATUS); VERIFY_(STATUS)
         
         IF(DEBUG_GF==1) THEN
          call MAPL_GetPointer(EXPORT, DTRDT_GF,'DTRDT_GF' ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS);DTRDT_GF=0.0
          call MAPL_GetPointer(EXPORT, DQDT_GF,'DQDT_GF' ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS);DQDT_GF=0.0
          call MAPL_GetPointer(EXPORT, DTDT_GF,'DTDT_GF' ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS);DTDT_GF=0.0
          call MAPL_GetPointer(EXPORT, MUPDP  ,'MUPDP'   ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS);MUPDP=0.0
          call MAPL_GetPointer(EXPORT, MUPSH  ,'MUPSH'   ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS);MUPSH=0.0
          call MAPL_GetPointer(EXPORT, MUPMD  ,'MUPMD'   ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS);MUPMD=0.0
          call MAPL_GetPointer(EXPORT, MFDP   ,'MFDP'    ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS);MFDP=0.0
          call MAPL_GetPointer(EXPORT, MFSH   ,'MFSH'    ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS);MFSH=0.0
          call MAPL_GetPointer(EXPORT, MFMD   ,'MFMD'    ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS);MFMD=0.0
          call MAPL_GetPointer(EXPORT, ERRDP  ,'ERRDP'   ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS);ERRDP=0.0
          call MAPL_GetPointer(EXPORT, ERRSH  ,'ERRSH'   ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS);ERRSH=0.0
          call MAPL_GetPointer(EXPORT, ERRMD  ,'ERRMD'   ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS);ERRMD=0.0
    
          call MAPL_GetPointer(EXPORT, RSU_CN_GF  ,'RSU_CN_GF'   ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS);RSU_CN_GF=0.0
          call MAPL_GetPointer(EXPORT, REV_CN_GF  ,'REV_CN_GF'   ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS);REV_CN_GF=0.0
          call MAPL_GetPointer(EXPORT, PFI_CN_GF  ,'PFI_CN_GF'   ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS);PFI_CN_GF=0.0
          call MAPL_GetPointer(EXPORT, PFL_CN_GF  ,'PFL_CN_GF'   ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS);PFL_CN_GF=0.0

          call MAPL_GetPointer(EXPORT, AA0      ,'AA0'     ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS);AA0=0.0
          call MAPL_GetPointer(EXPORT, AA1      ,'AA1'     ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS);AA1=0.0
          call MAPL_GetPointer(EXPORT, AA2      ,'AA2'     ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS);AA2=0.0
          call MAPL_GetPointer(EXPORT, AA3      ,'AA3'     ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS);AA3=0.0
          call MAPL_GetPointer(EXPORT, AA1_BL   ,'AA1_BL'  ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS);AA1_BL =0.0
          call MAPL_GetPointer(EXPORT, AA1_CIN  ,'AA1_CIN' ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS);AA1_CIN=0.0
          call MAPL_GetPointer(EXPORT, TAU_BL   ,'TAU_BL'  ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS);TAU_BL =0.0
          call MAPL_GetPointer(EXPORT, TAU_EC   ,'TAU_EC'  ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS);TAU_EC =0.0
          !print*,"sizes=",size(ERRMD),size(MUPMD);flush(6)
         ENDIF

! WMP
! Modify AREA (m^2) here so GF scale dependence has a CNV_FRACTION dependence
         if (GF_MIN_AREA > 0) then
           GF_AREA = GF_MIN_AREA*CNV_FRACTION + AREA*(1.0-CNV_FRACTION)
         else
           GF_AREA = AREA
         endif
! WMP
         
         !- call GF/GEOS5 interface routine
         call GF_GEOS5_Interface( IM,JM,LM,KM,ITRCR,LONS,LATS,DT_MOIST                       &
                                 ,T, PLE, PLO, ZLE, ZLO, PK, U, V, OMEGA            & 
                                 ,TH1, Q1, U1, V1, QLCN, QICN,QLLS,QILS, RASPRCP    &
                                 ,CNV_MF0, CNV_PRC3, CNV_MFD, CNV_DQLDT,ENTLAM      &
                                 ,CNV_MFC, CNV_UPDF, CNV_CVW, CNV_QC , CLCN         &                           
                                 ,QV_DYN_IN,PLE_DYN_IN,U_DYN_IN,V_DYN_IN,T_DYN_IN   &
                                 ,RADSW   ,RADLW  ,DQDT_BL  ,DTDT_BL                &
                                 ,FRLAND, GF_AREA,USTAR,TSTAR,QSTAR,T2M             &
                                 ,Q2M ,TA ,QA ,SH ,EVAP ,PHIS                       &
                                 ,KPBLIN    &
                                 ,MAPL_GRAV &
                                 ,SEEDCNV, SIGMA_DEEP, SIGMA_MID                    &
                                 ,DQDT_GF,DTDT_GF,MUPDP,MUPSH,MUPMD                 &
                                 ,MFDP,MFSH,MFMD,ERRDP,ERRSH,ERRMD                  &
                                 ,AA0,AA1,AA2,AA3,AA1_BL,AA1_CIN,TAU_BL,TAU_EC      &
                                 ,DTDTDYN,DQVDTDYN                                  &
                                 ,NCPL, NCPI, CNV_NICE, CNV_NDROP, CNV_FICE, CLDMICRO &
                                 ,RASPARAMS%QC_CRIT_CN, AUTOC_CN_OCN                &
                                 ,XHO,FSCAV,CNAMES,QNAMES,DTRDT_GF                  &
				 ,RSU_CN_GF,REV_CN_GF, PFI_CN_GF, PFL_CN_GF)
                                                                   
         HHO      =  0.0
         HSO      =  0.0    
         irccode  = -99
         trdlx    =  1.0   
         KEX      =  1.E-04

         ! adjust units to be [kg kg-1 s-1]
         REV_CN_GF = REV_CN_GF*iMASS
         RSU_CN_GF = RSU_CN_GF*iMASS

         call MAPL_TimerOff(STATE,"-GF")
      ENDIF
!-srf-gf-scheme 
!=====================================================================================


     IF(ADJUSTL(CONVPAR_OPTION) == 'RAS' .or. ADJUSTL(CONVPAR_OPTION) == 'BOTH') THEN

      call MAPL_TimerOn (STATE,"-RAS")

      ! MATMAT Export out the inputs before RAS needed for RAStest
      call MAPL_GetPointer(EXPORT, THOI,   'THOI'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, QHOI,   'QHOI'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, QSSI,   'QSSI'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DQSI,   'DQSI'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, PLEI,   'PLEI'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, TPERTI, 'TPERTI', RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, KCBLI,  'KCBLI' , RC=STATUS); VERIFY_(STATUS)

      if(associated(THOI  )) THOI   = TH1
      if(associated(QHOI  )) QHOI   = Q1
      if(associated(QSSI  )) QSSI   = QSS
      if(associated(DQSI  )) DQSI   = DQS
      if(associated(PLEI  )) PLEI   = CNV_PLE
      if(associated(TPERTI)) TPERTI = TPERT
      if(associated(KCBLI )) KCBLI  = KCBL

      ! temp kluge to differentiate ocean,land convective autoc (jtb 6/29/05)
      ! ---------------------------------------------------------------------
      where (FRLAND<0.1) 
         CO_AUTO = AUTOC_CN_OCN   ! ocean value
      elsewhere
         CO_AUTO = AUTOC_CN_LAND  ! land value
      end where

     call MAPL_TimerOn (STATE,"--RAS_RUN",RC=STATUS)
     VERIFY_(STATUS)

      RAS_ALPHA   = MAPL_UNDEF
      RAS_TAU     = MAPL_UNDEF
      RAS_TIME    = MAPL_UNDEF
      RAS_TRG     = MAPL_UNDEF 
      RAS_TOKI    = MAPL_UNDEF 
      RAS_PBL     = MAPL_UNDEF 
      RAS_WFN     = MAPL_UNDEF 

      call RASE(                        &                    
           IDIM                 , &
           IRUN                 , &                 
           K0                   , &
           ICMIN                , &
           DT_MOIST             , &  !!where? see below.
           MAPL_CP              , &
           MAPL_ALHL            , &
           MAPL_ALHS            , &
           MAPL_TICE            , &
           MAPL_GRAV            , &
           SEEDRAS              , &
           IRAS                 , &
           JRAS                 , &
           SIGE                 , &
           ! inputs for CBL
           KCBL                 , &
           WGT0                 , &
           WGT1                 , &
           ZCBLx                , &
           MXDIAMx              , &
           TPERT                , &
           QPERT                , &
           ! inputs
           TH1                  , &
           Q1                   , & 
           U1                   , &
           V1                   , &
           QSS                  , & 
           DQS                  , &
           CNV_FRACTION         , &
           RASAL2_2d            , &
           ! Pass in CO_AUTO
           CO_AUTO              , &
           ! - new for ras 2
           PK                   , &
           PLO                  , &
           GZLO                 , &   
           GZLE                 , &   
           QLCN                 , &   
           QICN                 , &   
           !
           CNV_PLE              , &
           PKE                  , &
           ! outputs
           CNV_DQLDT            , &   ! -> progno_cloud
           CNV_MF0              , &   ! -> diag
           CNV_MFD              , &   ! -> progno_cloud
           CNV_MFC              , &   ! -> diag
           CNV_PRC3             , &   ! -> progno_cloud 
           CNV_UPDF             , &   ! -> progno_cloud
           CNV_CVW              , &   ! 
           CNV_QC               , &   ! -> progno_cloud ???
           ENTLAM               , &   ! -> diag
           CLCN                 , &   ! -> upd if RAS-2 
           HHO                  , &
           HSO                  , &   ! -> diag
           RASPRCP              , &
           
           RASPARAMS            , & ! params
           RAS_NO_NEG           , &
           RAS_TIME, RAS_TRG, RAS_TOKI, RAS_PBL, RAS_WFN, &
           RAS_TAU        , &

!!!=======AER_CLOUD=======
           
           AeroProps    , &  !-> Aerosol properties
           CNV_FICE       , & !-> Fraction of ice in detrainment 
           CNV_NICE       , & !-> Detrained ice crystal concentration
           CNV_NDROP      , & !-> Detrained cloud droplet concentration
           RAS_ALPHA      , &

           
!!!==============              
           ITRCR                , &
           irccode              , &
           XHO = XHO            , &
           TRIEDLEV_DIAG = trdlx, &
           FSCAV  = FSCAV       , &
           DISSKE = KEX           )

      call MAPL_TimerOff(STATE,"--RAS_RUN",RC=STATUS)
      VERIFY_(STATUS)

      if(associated(QVRAS  )) QVRAS    = Q1
      if(associated(THRAS  )) THRAS    = TH1
      if(associated(URAS   )) URAS     = U1
      if(associated(VRAS   )) VRAS     = V1

      if(associated(MXDIAM )) MXDIAM  = MXDIAMx
      if(associated(RCCODE )) RCCODE  = 1.0*IRCCODE
      if(associated(TRIEDLV)) TRIEDLV = TRDLX

      if(associated(TVEX   )) TVEX     = SUM( (MAPL_CP*TEMP + MAPL_ALHL*Q)*MASS, 3 )

      call MAPL_TimerOff(STATE,"-RAS")
     ENDIF

     IF(ADJUSTL(CONVPAR_OPTION) == 'NONE') THEN
        !-zero output convpar arrays 
         CNV_MF0  =0.0 ! 'cloud_base_mass_flux'              - 'kg m-2 s-1'
         CNV_MFD  =0.0 ! 'detraining_mass_flux',             - 'kg m-2 s-1'    
         CNV_MFC  =0.0 ! 'cumulative_mass_flux',             - 'kg m-2 s-1'   
         CNV_CVW  =0.0 ! 'updraft_vertical_velocity',        - 'hPa s-1', 
         CNV_QC   =0.0 ! 'grid_mean_convective_condensate',  - 'kg kg-1',   
         CNV_UPDF =0.0 ! 'updraft_areal_fraction',           - '1', 
         ENTLAM   =0.0 ! 'entrainment parameter',            - 'm-1',  
         CNV_PRC3 =0.0 ! 'convective_precipitation           - 'kg m-2 s-1'
        ! Move any convective condensate/cloud to large-scale
         QLLS = QLLS+QLCN
         QILS = QILS+QICN
         CLLS = MIN(CLLS+CLCN,1.0)
        ! Zero out convective condensate/cloud
         QLCN = 0.0
         QICN = 0.0
         CLCN = 0.0
     ENDIF

      call MAPL_TimerOn (STATE,"-POST_RAS")

      CNV_NICE_X = CNV_NICE*CNV_MFD
      CNV_NDROP_X = CNV_NDROP*CNV_MFD
        
      where (CNV_MFD   .le. 0.0)
          CNV_FICE = 0.0
      end where
        
      if(associated(ZCBL   )) ZCBL     = ZCBLx
      if(associated(DQRC   )) DQRC     = CNV_PRC3 / DT_MOIST
      if(associated(DQDTCN )) DQDTCN   = ( Q1  -  DQDTCN  ) / DT_MOIST
      if(associated(DTHDTCN)) DTHDTCN  = ( TH1 -  DTHDTCN ) / DT_MOIST
      if(associated(DQCDTCN)) DQCDTCN  = CNV_DQLDT * MAPL_GRAV / ( PLE(:,:,1:LM)-PLE(:,:,0:LM-1) )
      if(associated(DTDTCN )) DTDTCN   = ( TH1*PK -  DTDTCN ) / DT_MOIST

#ifdef USEDDF
!-srf-gf-scheme
   IF(ADJUSTL(CONVPAR_OPTION) == 'RAS' .or. ADJUSTL(CONVPAR_OPTION) == 'BOTH') THEN
!-srf-gf-scheme

      DDF_CTL%ALPHA_DDF_UDF  = 0.5
      DDF_CTL%AREAL_FRACTION = 0.1

      DDF_CTL%PARTITION%TYPE  = "HEIGHT_DEP"
      DDF_CTL%PARTITION%MAX_RATIO = 0.9
      DDF_CTL%PARTITION%MIN_RATIO = 0.1
      DDF_CTL%THETA_IS_THVAR      = .TRUE.

      do I = 1, IM
         do J = 1, JM

            ! Call to ddf v 1.10.2.4 (tag=merra_dndrft_03)
            ! needs to change if more recent ddf used
            !----------------------------------------------
            call DDF1( DDF_CTL , DT_MOIST , TH1(I,J,:) , & 
                 Q1(I,J,:), QLCN(I,J,:), QICN(I,J,:), & 
                 CNV_PRC3(I,j,:), CNV_HAIL(I,j,:), &
                 ZLO(I,J,:), ZLE(I,J,:), PLO(I,J,:),  & 
                 CNV_PLE(I,J,:),  PK(I,J,:), CNV_MFC(I,J,:) ,  & 
                 DDF_ZSCALEz, DDF_DQDTz , DDF_DTDTz, DDF_MFCz )

            if(associated(DDF_DQDT   ))  DDF_DQDT   (I,J,:)  = DDF_DQDTz
            if(associated(DDF_DTDT   ))  DDF_DTDT   (I,J,:)  = DDF_DTDTz
            if(associated(DDF_MFC    ))  DDF_MFC    (I,J,:)  = DDF_MFCz
            if(associated(DDF_ZSCALE ))  DDF_ZSCALE (I,J)    = DDF_ZSCALEz

         end do
      end do
!-srf-gf-scheme
   ENDIF
!-srf-gf-scheme
#endif

      call MAPL_TimerOff(STATE,"-POST_RAS")

      call MAPL_TimerOn (STATE,"--UWSHCU")
      if (DOSHLW /= 0) then

      !  Call UW shallow convection
      !----------------------------------------------------------------

      call compute_uwshcu_inv(IDIM, K0, ITRCR, DT_MOIST,  & ! IN
            PLO*100., ZLO, PK, PLE, ZLE, PKE, DP,         &
            U1, V1, Q1, QLLS, QILS, TH1, TKE, KPBLSC,     &
            THLSRC_PERT,                                  &
            CUSH, XHO,                                    & ! INOUT
            UMF_SC, DQVDT_SC, DQLDT_SC, DQIDT_SC,         & ! OUT
            DTHDT_SC, DUDT_SC, DVDT_SC, DQRDT_SC,         &
            DQSDT_SC, CUFRC_SC, ENTR_SC, DETR_SC,         &
            QLDET_SC, QIDET_SC, QLSUB_SC, QISUB_SC,       &
            SC_NDROP, SC_NICE,                            &
#ifdef UWDIAG 
            QCU_SC, QLU_SC,                               & ! DIAG ONLY 
            QIU_SC, CBMF_SC, DQCDT_SC, CNT_SC, CNB_SC,    &
            CIN_SC, PLCL_SC, PLFC_SC, PINV_SC, PREL_SC,   &
            PBUP_SC, WLCL_SC, QTSRC_SC, THLSRC_SC,        &
            THVLSRC_SC, TKEAVG_SC, CLDTOP_SC, WUP_SC,     &
            QTUP_SC, THLUP_SC, THVUP_SC, UUP_SC, VUP_SC,  &
            XC_SC,                                        &
#endif 
            USE_TRACER_TRANSP_UW, SHLWPARAMS )


      !  Apply tendencies
      !--------------------------------------------------------------
      Q1  = Q1  + DQVDT_SC * DT_MOIST    ! note this adds to the convective
      TH1 = TH1 + DTHDT_SC * DT_MOIST    !  tendencies calculated below
      U1  = U1  + DUDT_SC * DT_MOIST
      V1  = V1  + DVDT_SC * DT_MOIST

      !  Calculate detrained mass flux
      !--------------------------------------------------------------
      where (DETR_SC.ne.MAPL_UNDEF)
        MFD_SC = 0.5*(UMF_SC(:,:,1:LM)+UMF_SC(:,:,0:LM-1))*DETR_SC*DP
      elsewhere
        MFD_SC = 0.0
      end where
 
      !  Convert detrained water units before passing to cloud
      !---------------------------------------------------------------
        QLENT_SC = 0.
        QIENT_SC = 0.
        WHERE (QLDET_SC.lt.0.) 
          QLENT_SC = QLDET_SC
          QLDET_SC = 0.
        END WHERE
        WHERE (QIDET_SC.lt.0.) 
          QIENT_SC = QIDET_SC
          QIDET_SC = 0.
        END WHERE
        QLDET_SC = QLDET_SC*MASS
        QIDET_SC = QIDET_SC*MASS
      !  Apply condensate tendency from subsidence, and sink from
      !  condensate entrained into shallow updraft. 
      !  Detrained condensate added in microphysics below.
      !-------------------------------------------------------------
        QLLS = QLLS + (QLSUB_SC+QLENT_SC)*DT_MOIST
        QILS = QILS + (QISUB_SC+QIENT_SC)*DT_MOIST

      !  Calculate updraft core fraction from cumulus fraction.
      !  CUFRC is assumed in compute_uwshcu to be twice updraft frac
      !--------------------------------------------------------------
      UFRC_SC = 0.5 * CUFRC_SC

      !  Number concentrations for 2-moment microphysics
      !--------------------------------------------------------------
      SC_NDROP = SC_NDROP*MASS
      SC_NICE = SC_NICE*MASS

      !  Precipitation
      !--------------------------------------------------------------
      SHLW_PRC3 = DQRDT_SC    ! [kg/kg/s]
      SHLW_SNO3 = DQSDT_SC    ! [kg/kg/s]

      else   ! if UW shallow scheme not called

        UMF_SC    = 0.
        MFD_SC    = 0.
        SHLW_PRC3 = 0.
        SHLW_SNO3 = 0.      
        UFRC_SC   = 0.
        SC_NDROP  = 0.
        SC_NICE   = 0.
        QLDET_SC  = 0.
        QIDET_SC  = 0.
        QLSUB_SC  = 0.
        QISUB_SC  = 0.

      end if

      
      call MAPL_TimerOff (STATE,"--UWSHCU")


      call MAPL_TimerOn(STATE,"-POST_RAS")

      !     Compute new mass loading for aerosols; CAR 12/19/08
      !     -----------------------------------------------------
      !     Loop over tracers, accumulate, subtract new from initial, divide
      !     time to get kg/m2/s

      KK=0
      do K=1,KM
         if(IS_FRIENDLY(K)) then
            KK = KK + 1
            QNAME = trim(QNAMES(K))
            CNAME = trim(CNAMES(K))
            if((CNAME == 'DU') .or. (CNAME == 'SS') .or. (CNAME == 'NI') .or. (CNAME == 'SU') .or. &
               (CNAME == 'CA.oc') .or. (CNAME == 'CA.bc') .or. (CNAME == 'CA.br')) then   ! Diagnostics for GOCART2G tracers
               SELECT CASE (QNAME(1:3))
               CASE ('DU0')
                  if(associated(DDU2gDT)) then
                     DDU2gDT = DDU2gDT + sum(XHO(:,:,:,KK)*DP(:,:,:),dim=3)
                  end if
              CASE ('SS0')
                  if(associated(DSS2gDT)) then
                     DSS2gDT = DSS2gDT + sum(XHO(:,:,:,KK)*DP(:,:,:),dim=3)
                  end if
               CASE ('SO4')
                  if(associated(DSU2gDT)) then
                     DSU2gDT = DSU2gDT + sum(XHO(:,:,:,KK)*DP(:,:,:),dim=3)
                  end if
               CASE ('NO3')
                  if(associated(DNI2gDT)) then
                     DNI2gDT = DNI2gDT + sum(XHO(:,:,:,KK)*DP(:,:,:),dim=3)
                  end if
               CASE ('NH3')
                  if(associated(DNH32gDT)) then
                     DNH32gDT = DNH32gDT + sum(XHO(:,:,:,KK)*DP(:,:,:),dim=3)
                  end if
               CASE ('NH4')
                  if(associated(DNH4A2gDT)) then
                     DNH4A2gDT = DNH4A2gDT + sum(XHO(:,:,:,KK)*DP(:,:,:),dim=3)
                  end if
               END SELECT

               SELECT CASE (QNAME(1:13))
               CASE ('CA.bcphilic')
                  if(associated(DBC2gDT)) then
                     DBC2gDT = DBC2gDT + sum(XHO(:,:,:,KK)*DP(:,:,:),dim=3)
                  end if
!               CASE ('OCp')
               CASE ('CA.ocphilic')
                  if(associated(DOC2gDT)) then
                     DOC2gDT = DOC2gDT + sum(XHO(:,:,:,KK)*DP(:,:,:),dim=3)
                  end if
               CASE ('CA.brphilic')
                  if(associated(DBRC2gDT)) then
                     DBRC2gDT = DBRC2gDT + sum(XHO(:,:,:,KK)*DP(:,:,:),dim=3)
                  end if
               END SELECT
            endif

            if(CNAME == 'GOCART') then   ! Diagnostics for GOCART tracers
               SELECT CASE (QNAME(1:3))
               CASE ('du0')
                  if(associated(DDUDT)) then
                     DDUDT = DDUDT + sum(XHO(:,:,:,KK)*DP(:,:,:),dim=3)
                  end if
               CASE ('ss0')
                  if(associated(DSSDT)) then
                     DSSDT = DSSDT + sum(XHO(:,:,:,KK)*DP(:,:,:),dim=3)
                  end if
               CASE ('BCp')
                  if(associated(DBCDT)) then
                     DBCDT = DBCDT + sum(XHO(:,:,:,KK)*DP(:,:,:),dim=3)
                  end if
               CASE ('OCp')
                  if(associated(DOCDT)) then
                     DOCDT = DOCDT + sum(XHO(:,:,:,KK)*DP(:,:,:),dim=3)
                  end if
               CASE ('SO4')
                  if(associated(DSUDT)) then
                     DSUDT = DSUDT + sum(XHO(:,:,:,KK)*DP(:,:,:),dim=3)
                  end if
               CASE ('NO3')
                  if(associated(DNIDT)) then
                     DNIDT = DNIDT + sum(XHO(:,:,:,KK)*DP(:,:,:),dim=3)
                  end if
               CASE ('NH3')
                  if(associated(DNH3DT)) then
                     DNH3DT = DNH3DT + sum(XHO(:,:,:,KK)*DP(:,:,:),dim=3)
                  end if
               CASE ('NH4')
                  if(associated(DNH4ADT)) then
                     DNH4ADT = DNH4ADT + sum(XHO(:,:,:,KK)*DP(:,:,:),dim=3)
                  end if
               CASE ('BRC')
                  if(associated(DBRCDT)) then
                     DBRCDT = DBRCDT + sum(XHO(:,:,:,KK)*DP(:,:,:),dim=3)
                  end if 
               END SELECT
            endif
         end if
         if(CNAME == 'CARMA') then   ! Diagnostics for CARMA tracers
            ! Check name to see if it is a "pc" element
            ENAME = ''
            ind= index(QNAME, '::')
            if (ind> 0) then
               ENAME = trim(QNAME(ind+2:ind+3))  ! Component name (e.g., GOCART, CARMA)
               if(ENAME == 'pc') then
                  SELECT CASE (QNAME(1:4))
                  CASE ('dust') ! CARMA DUST
                     if(associated(DDUDTcarma)) then
                        DDUDTcarma = DDUDTcarma + sum(XHO(:,:,:,KK)*DP(:,:,:),dim=3) 
                     end if
                  CASE ('seas') ! CARMA SEASALT
                     if(associated(DSSDTcarma)) then
                        DSSDTcarma = DSSDTcarma + sum(XHO(:,:,:,KK)*DP(:,:,:),dim=3) 
                     end if
                  END SELECT
               endif
            endif
         endif
      end do

      if (associated(DDU2gDT))   DDU2gDT = (DDU2gDT - CMDU2g) / (MAPL_GRAV*DT_MOIST)
      if (associated(DSS2gDT))   DSS2gDT = (DSS2gDT - CMSS2g) / (MAPL_GRAV*DT_MOIST)
      if (associated(DBC2gDT))   DBC2gDT = (DBC2gDT - CMBC2g) / (MAPL_GRAV*DT_MOIST)
      if (associated(DOC2gDT))   DOC2gDT = (DOC2gDT - CMOC2g) / (MAPL_GRAV*DT_MOIST)
      if (associated(DSU2gDT))   DSU2gDT = (DSU2gDT - CMSU2g) / (MAPL_GRAV*DT_MOIST)
      if (associated(DNI2gDT))   DNI2gDT = (DNI2gDT - CMNI2g) / (MAPL_GRAV*DT_MOIST)
      if (associated(DNH32gDT))  DNH32gDT = (DNH32gDT - CMNH32g) / (MAPL_GRAV*DT_MOIST)
      if (associated(DNH4A2gDT)) DNH4A2gDT = (DNH4A2gDT - CMNH4A2g) / (MAPL_GRAV*DT_MOIST)
      if (associated(DBRC2gDT))  DBRC2gDT= (DBRC2gDT- CMBRC2g)/ (MAPL_GRAV*DT_MOIST)

      if (associated(DDUDT))   DDUDT = (DDUDT - CMDU) / (MAPL_GRAV*DT_MOIST)
      if (associated(DSSDT))   DSSDT = (DSSDT - CMSS) / (MAPL_GRAV*DT_MOIST)
      if (associated(DBCDT))   DBCDT = (DBCDT - CMBC) / (MAPL_GRAV*DT_MOIST)
      if (associated(DOCDT))   DOCDT = (DOCDT - CMOC) / (MAPL_GRAV*DT_MOIST)
      if (associated(DSUDT))   DSUDT = (DSUDT - CMSU) / (MAPL_GRAV*DT_MOIST)
      if (associated(DNIDT))   DNIDT = (DNIDT - CMNI) / (MAPL_GRAV*DT_MOIST)
      if (associated(DNH3DT))  DNH3DT = (DNH3DT - CMNH3) / (MAPL_GRAV*DT_MOIST)
      if (associated(DNH4ADT)) DNH4ADT = (DNH4ADT - CMNH4A) / (MAPL_GRAV*DT_MOIST)
      if (associated(DBRCDT))  DBRCDT= (DBRCDT- CMBRC)/ (MAPL_GRAV*DT_MOIST)

      if (associated(DDUDTcarma))  DDUDTcarma = (DDUDTcarma - CMDUcarma) / (MAPL_GRAV*DT_MOIST)
      if (associated(DSSDTcarma))  DSSDTcarma = (DSSDTcarma - CMSScarma) / (MAPL_GRAV*DT_MOIST)


      ! Fill in tracer tendencies
      !--------------------------

      KK=0
      do K=1,KM
         if(IS_FRIENDLY(K)) then
            KK = KK+1
            TRPtrs(K)%Q(:,:,:) =  XHO(:,:,:,KK)
         end if
      end do

      call MAPL_TimerOff(STATE,"-POST_RAS")


      call MAPL_GetResource( STATE, CLDPARAMS%PDFSHAPE,  'PDFSHAPE:',   DEFAULT= 1.0    )

      call MAPL_GetResource( STATE, CLDPARAMS%TURNRHCRIT_UP, 'TURNRHCRIT_UP:', DEFAULT= 300.0  )
      call MAPL_GetResource( STATE, CLDPARAMS%SLOPERHCRIT, 'SLOPERHCRIT:', DEFAULT= 20.0  )

      ! Horizontal resolution dependant defaults for minimum RH crit
      if( imsize.le.200       ) call MAPL_GetResource( STATE, CLDPARAMS%MINRHCRIT, 'MINRHCRIT:', DEFAULT=0.80, RC=STATUS)
      if( imsize.gt.200 .and. &
          imsize.le.400       ) call MAPL_GetResource( STATE, CLDPARAMS%MINRHCRIT, 'MINRHCRIT:', DEFAULT=0.90, RC=STATUS)
      if( imsize.gt.400 .and. &
          imsize.le.800       ) call MAPL_GetResource( STATE, CLDPARAMS%MINRHCRIT, 'MINRHCRIT:', DEFAULT=0.93, RC=STATUS)
      if( imsize.gt.800 .and. &
          imsize.le.1600      ) call MAPL_GetResource( STATE, CLDPARAMS%MINRHCRIT, 'MINRHCRIT:', DEFAULT=0.95, RC=STATUS)
      if( imsize.gt.1600 .and. &
          imsize.le.3200      ) call MAPL_GetResource( STATE, CLDPARAMS%MINRHCRIT, 'MINRHCRIT:', DEFAULT=0.97 ,RC=STATUS)
      if( imsize.gt.3200 .and. &
          imsize.le.6400      ) call MAPL_GetResource( STATE, CLDPARAMS%MINRHCRIT, 'MINRHCRIT:', DEFAULT=0.98 ,RC=STATUS)
      if( imsize.gt.6400 .and. &
          imsize.le.12800     ) call MAPL_GetResource( STATE, CLDPARAMS%MINRHCRIT, 'MINRHCRIT:', DEFAULT=0.99 ,RC=STATUS)
      if( imsize.gt.12800     ) call MAPL_GetResource( STATE, CLDPARAMS%MINRHCRIT, 'MINRHCRIT:', DEFAULT=0.99 ,RC=STATUS)

      call MAPL_GetResource( STATE, CLDPARAMS%MAXRHCRIT    , 'MAXRHCRIT:'    , DEFAULT= 1.0 )
      call MAPL_GetResource( STATE, CLDPARAMS%MAXRHCRITLAND, 'MAXRHCRITLAND:', DEFAULT= 1.0 )

      if (DOCLDMACRO==0) then
        call MAPL_TimerOn(STATE,"---CLDMACRO")
        TEMP = TH1*PK
        DTDT_macro=TEMP
        DQVDT_macro=Q1
        DQLDT_macro=QLCN+QLLS
        DQIDT_macro=QICN+QILS
        DQADT_macro=CLCN+CLLS
       ! add DeepCu QL/QI/CL to Convective
        do K=1,LM
          do J=1,JM
           do I=1,IM
            IFRC = ICE_FRACTION( TEMP(I,J,K), CNV_FRACTION(I,J), &
                                 SNOMAS(I,J), FRLANDICE(I,J), FRLAND(I,J) )
            TEND = CNV_DQLDT(I,J,K)*iMASS(I,J,K)*DT_MOIST
            QLCN(I,J,K) = QLCN(I,J,K) + (1.0-IFRC)*TEND
            QICN(I,J,K) = QICN(I,J,K) +      IFRC *TEND
           ! dont forget that conv cond has never frozen !!!!
            IF(ADJUSTL(CONVPAR_OPTION) /= 'GF') THEN
              TEMP(I,J,K) = TEMP(I,J,K) + (MAPL_ALHS-MAPL_ALHL) * IFRC * TEND / MAPL_CP
            ENDIF
           enddo
          enddo
        enddo
       ! add DeepCu Clouds to Convective
        CLCN = CLCN + CNV_MFD*iMASS*DT_MOIST
        if (UWTOLS/=0) then
       ! add ShallowCu CL/QL/QI tendencies to Large-Scale
          CLLS = CLLS +   MFD_SC*iMASS*DT_MOIST
          QLLS = QLLS + QLDET_SC*iMASS*DT_MOIST
          QILS = QILS + QIDET_SC*iMASS*DT_MOIST
        else
          CLCN = CLCN +   MFD_SC*iMASS*DT_MOIST
!          CLCN = CLCN +   DCM_SC*iMASS*DT_MOIST
          QLCN = QLCN + QLDET_SC*iMASS*DT_MOIST
          QICN = QICN + QIDET_SC*iMASS*DT_MOIST

          CLCN = max(min(CLCN,1.0),0.0)

          do K=1,LM
            do J=1,JM
              do I=1,IM

                if (CLCN(i,j,k).lt.0.99) then

                  QT = Q1(i,j,k) + (QLLS(i,j,k)+QILS(i,j,k))/(1.-CLCN(i,j,k))   ! QT in non-convective area
                 
                  do n = 1,5

                    qsatn = GEOS_QSAT( TEMP(i,j,k), PLO(i,j,k) )

                    sigmaqt = CLDPARAMS%MINRHCRIT + (CLDPARAMS%MAXRHCRIT-CLDPARAMS%MINRHCRIT)/(19.) * &
                        ((atan( (2.*(PLO(i,j,k)-CNV_PLE(i,j,LM)+260.)/(260.)-1.) * &
                        tan(20.*MAPL_PI/21.-0.5*MAPL_PI) ) + 0.5*MAPL_PI) * 21./MAPL_PI - 1.)

                    sigmaqt = (1.0-sigmaqt)*qsatn

                    call pdffrac(INT(CLDPARAMS%PDFSHAPE),QT,sigmaqt,sigmaqt,qsatn,cfn)
                    call pdfcondensate(INT(CLDPARAMS%PDFSHAPE),QT,sigmaqt,sigmaqt,qsatn,qcn)

                    IFRC = ICE_FRACTION( TEMP(I,J,K), 0.0, SNOMAS(I,J), FRLANDICE(I,J), FRLAND(I,J) )
                   
                    ! calculate change in grid mean QLLS and QILS
                    ! delta = new QXLS - old QXLS
                    dqils = 0.75*qcn*IFRC*(1.-CLCN(i,j,k)) - QILS(i,j,k)
                    dqlls = 0.75*qcn*(1.-IFRC)*(1.-CLCN(i,j,k)) - QLLS(i,j,k)
                   
                    TEMP(i,j,k) = TEMP(i,j,k) + (MAPL_ALHL*(dqils+dqlls)+MAPL_ALHF*dqils)/ MAPL_CP
                    Q1(i,j,k) = Q1(i,j,k) - dqils - dqlls
                    QLLS(i,j,k) = QLLS(i,j,k) + dqlls
                    QILS(i,j,k) = QILS(i,j,k) + dqils

                  end do ! n convergence loop

                  ! cfn is fraction in non-convective area. convert to grid area.
                  CLLS(i,j,k) = cfn*(1.-CLCN(i,j,k))

                end if ! if clcn<0.99

              end do ! IM loop
            end do ! JM loop
          end do ! LM loop
        
        endif
       ! add ShallowCu rain/snow tendencies
        QRAIN = QRAIN + SHLW_PRC3*DT_MOIST
        QSNOW = QSNOW + SHLW_SNO3*DT_MOIST
       ! add DeepCu rain/snow
        if (ADJUSTL(CONVPAR_OPTION) /= 'GF') then
          where (TEMP < MAPL_TICE) !SNOW
            QSNOW = QSNOW + CNV_PRC3*iMASS*DT_MOIST
            TEMP  = TEMP  + CNV_PRC3*iMASS*(MAPL_ALHS-MAPL_ALHL) / MAPL_CP
          else where ! RAIN
            QRAIN = QRAIN + CNV_PRC3*iMASS*DT_MOIST
          end where
        end if
     ! Clean up clouds before microphysics
        CALL FIX_UP_CLOUDS( TEMP, Q1, QLLS, QILS, CLLS, QLCN, QICN, CLCN )
     ! Clean up any negative specific humidity
        CALL FILLQ2ZERO2( Q1, MASS, FILLQ  )
        DTDT_macro=  (TEMP-DTDT_macro)/DT_MOIST
        DQVDT_macro=(Q1-DQVDT_macro)/DT_MOIST
        DQLDT_macro=((QLCN+QLLS)-DQLDT_macro)/DT_MOIST
        DQIDT_macro=((QICN+QILS)-DQIDT_macro)/DT_MOIST
        DQADT_macro=((CLCN+CLLS)-DQADT_macro)/DT_MOIST
        TH1 = TEMP/PK
     ! Zero-out 3D CNV/ANV/SHL CLDMACRO Precipitation & Fluxes
        PFI_CN_X = 0.0
        PFL_CN_X = 0.0
        PFI_AN_X = 0.0
        PFL_AN_X = 0.0
        PFI_SC_X = 0.0
        PFL_SC_X = 0.0
        AN_PRC2  = 0.0
        CN_PRC2  = 0.0
        SC_PRC2  = 0.0
        AN_SNR   = 0.0
        CN_SNR   = 0.0
        SC_SNR   = 0.0
        call MAPL_TimerOff(STATE,"---CLDMACRO")
      endif

      call MAPL_TimerOn (STATE,"-MISC2")

      !Trajectory for Moist TLM/ADJ
      if(associated(TS_moist)) TS_moist = TS
      if(associated(KCBL_moist)) KCBL_moist = KCBL
      if(associated(KHu_moist)) then
         KHu_moist = -1.0
         DO i = 1,IM
            DO j = 1,JM
               DO l = 0,LM
                  if (KH(i,j,l) .gt. 2.0) then
                     KHu_moist(i,j) = l * 1.0
                     exit
                  endif
               endDO
            endDO
         endDO
      endif
      if(associated(KHl_moist)) then
         KHl_moist = -1.0
         DO i = 1,IM
            DO j = 1,JM
               DO l = LM,0,-1
                  if (KH(i,j,l) .gt. 2.0) then
                     KHl_moist(i,j) = l * 1.0
                     exit
                  endif
               endDO
            endDO
         endDO
      endif
      if(associated(ctop_moist)) then
         ctop_moist(:,:) = 0.0
         do j=1,jm
            do i=1,im
               do l=1,lm
                  if( ( TH1(i,j,l) - TH(i,j,l) )/DT_moist .gt. 10e-6 )   then
                     ctop_moist(i,j) = l
                     exit
                  endif
               enddo
            enddo
         enddo
      endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!AER_CLOOUD

      if(adjustl(CLDMICRO)=="2MOMENT") then

        KCT = 20 !default upper limit. Less than 20 makes no difference
         ! Find Convective Cloud Top
	   call find_l(KCT, CNV_DQLDT, 1.0e-9, IM, JM, LM, 20, LM-2) 
	 
      end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!AER_CLOUD


      !-- intermediate water column (all phases and reservoirs)------
!      if(associated(TVQX1))  TVQX1     = SUM( (  Q1 +  QLLS + QLCN + QILS + QICN + CNV_PRC3)*MASS + CNV_DQLDT*DT_MOIST , 3 )

      deallocate(IS_FRIENDLY,stat=STATUS)
      VERIFY_(STATUS)
      deallocate(IS_WEIGHTED,stat=STATUS)
      VERIFY_(STATUS)
      deallocate(TRPtrs     ,stat=STATUS)
      VERIFY_(STATUS)
      deallocate(XHO        ,stat=STATUS)
      VERIFY_(STATUS)
      ! CAR 12/5/08
      deallocate(FSCAV      ,stat=STATUS)
      VERIFY_(STATUS)
      deallocate(FSCAV_     ,stat=STATUS)
      VERIFY_(STATUS)
      deallocate(QNAMES, CNAMES ,stat=STATUS)
      VERIFY_(STATUS)
      IF(ADJUSTL(CONVPAR_OPTION) == 'GF') THEN
         DEALLOCATE(Hcts , stat=STATUS); VERIFY_(STATUS)
      ENDIF

      call MAPL_GetResource( STATE, CLDPARAMS%CNV_BETA,       'CNV_BETA:',       DEFAULT= 10.0    )
      call MAPL_GetResource( STATE, CLDPARAMS%ANV_BETA,       'ANV_BETA:',       DEFAULT= 4.0     )
      call MAPL_GetResource( STATE, CLDPARAMS%LS_BETA,        'LS_BETA:',        DEFAULT= 4.0     )
      call MAPL_GetResource( STATE, CLDPARAMS%RH_CRIT,        'RH_CRIT:',        DEFAULT= 1.0     )
      call MAPL_GetResource( STATE, CLDPARAMS%QC_CRIT_LS,     'QC_CRIT_LS:',     DEFAULT= 8.0e-4  )
      call MAPL_GetResource( STATE, CLDPARAMS%ACCRETION,      'ACCRETION:',      DEFAULT= 2.0     )

      call MAPL_GetResource(STATE, CLDPARAMS%RAIN_REVAP_FAC, 'RAIN_REVAP_FAC:', DEFAULT= 1.0 ,RC=STATUS)

      call MAPL_GetResource( STATE, CLDPARAMS%VOL_TO_FRAC,    'VOL_TO_FRAC:',    DEFAULT= -1.0    )
      call MAPL_GetResource( STATE, CLDPARAMS%SUPERSAT,       'SUPERSAT:',       DEFAULT= 0.0     )
      call MAPL_GetResource( STATE, CLDPARAMS%SHEAR_EVAP_FAC, 'SHEAR_EVAP_FAC:', DEFAULT= 1.3     )
      call MAPL_GetResource( STATE, CLDPARAMS%MIN_ALLOW_CCW,  'MIN_ALLOW_CCW:',  DEFAULT= 1.0e-9  )
      
      IF(ADJUSTL(CONVPAR_OPTION) == 'GF' .and. icumulus_gf(shal) == 1) THEN
        call MAPL_GetResource( STATE, CLDPARAMS%CCW_EVAP_EFF,   'CCW_EVAP_EFF:',   DEFAULT= 5.0e-4  )
        call MAPL_GetResource( STATE, CLDPARAMS%AUTOC_LS,       'AUTOC_LS:',       DEFAULT= 1.5e-3  )
        call MAPL_GetResource( STATE, CLDPARAMS%AUTOC_ANV,      'AUTOC_ANV:',      DEFAULT= 1.5e-3  )
      ELSE
        call MAPL_GetResource( STATE, CLDPARAMS%CCW_EVAP_EFF,   'CCW_EVAP_EFF:',   DEFAULT= 4.0e-3  )
        call MAPL_GetResource( STATE, CLDPARAMS%AUTOC_LS,       'AUTOC_LS:',       DEFAULT= 1.0e-3  )
        call MAPL_GetResource( STATE, CLDPARAMS%AUTOC_ANV,      'AUTOC_ANV:',      DEFAULT= 1.0e-3  )
      ENDIF
      
      call MAPL_GetResource( STATE, CLDPARAMS%NSUB_AUTOCONV,  'NSUB_AUTOCONV:',  DEFAULT= 20.     )
      call MAPL_GetResource( STATE, CLDPARAMS%LS_SUND_INTER,  'LS_SUND_INTER:',  DEFAULT= 1.0     )
      call MAPL_GetResource( STATE, CLDPARAMS%LS_SUND_COLD,   'LS_SUND_COLD:',   DEFAULT= 1.0     )
      call MAPL_GetResource( STATE, CLDPARAMS%LS_SUND_TEMP1,  'LS_SUND_TEMP1:',  DEFAULT= 230.    )
      call MAPL_GetResource( STATE, CLDPARAMS%ANV_SUND_INTER, 'ANV_SUND_INTER:', DEFAULT= 1.0     )
      call MAPL_GetResource( STATE, CLDPARAMS%ANV_SUND_COLD,  'ANV_SUND_COLD:',  DEFAULT= 1.0     )
      call MAPL_GetResource( STATE, CLDPARAMS%ANV_SUND_TEMP1, 'ANV_SUND_TEMP1:', DEFAULT= 230.    )
      call MAPL_GetResource( STATE, CLDPARAMS%ANV_TO_LS_TIME, 'ANV_TO_LS_TIME:', DEFAULT= 14400.  )
      call MAPL_GetResource( STATE, CLDPARAMS%NCCN_WARM,      'NCCN_WARM:',      DEFAULT= 50.     )
      call MAPL_GetResource( STATE, CLDPARAMS%NCCN_ICE,       'NCCN_ICE:',       DEFAULT= 0.01    )
      call MAPL_GetResource( STATE, CLDPARAMS%NCCN_ANVIL,     'NCCN_ANVIL:',     DEFAULT= 0.1     )
      call MAPL_GetResource( STATE, CLDPARAMS%NCCN_PBL,       'NCCN_PBL:',       DEFAULT= 200.    )
      call MAPL_GetResource( STATE, CLDPARAMS%DISABLE_RAD,    'DISABLE_RAD:',    DEFAULT= 0.      )
      call MAPL_GetResource( STATE, CLDPARAMS%REVAP_OFF_P,    'REVAP_OFF_P:',    DEFAULT= 2000.   )
      call MAPL_GetResource( STATE, CLDPARAMS%ICE_RAMP,       'ICE_RAMP:',       DEFAULT= -27.0   )
      call MAPL_GetResource( STATE, CLDPARAMS%CNV_ICEPARAM,   'CNV_ICEPARAM:',   DEFAULT= 1.0     )
      call MAPL_GetResource( STATE, CLDPARAMS%CNV_ICEFRPWR,   'CNV_ICEFRPWR:',   DEFAULT= 4.0     )
      call MAPL_GetResource( STATE, CLDPARAMS%CNV_DDRF,       'CNV_DDRF:',       DEFAULT= 0.0     )
      call MAPL_GetResource( STATE, CLDPARAMS%ANV_DDRF,       'ANV_DDRF:',       DEFAULT= 0.0     )
      call MAPL_GetResource( STATE, CLDPARAMS%LS_DDRF,        'LS_DDRF:',        DEFAULT= 0.0     )
      call MAPL_GetResource( STATE, CLDPARAMS%QC_CRIT_ANV,    'QC_CRIT_ANV:',    DEFAULT= 8.0e-4  )
      call MAPL_GetResource( STATE, CLDPARAMS%TANHRHCRIT,     'TANHRHCRIT:',     DEFAULT= 1.0     )

      if( LM .eq. 72 ) then
        call MAPL_GetResource( STATE, CLDPARAMS%ICE_SETTLE,     'ICE_SETTLE:',     DEFAULT= 1.      )
        call MAPL_GetResource( STATE, CLDPARAMS%ANV_ICEFALL,    'ANV_ICEFALL:',    DEFAULT= 1.0     )
        call MAPL_GetResource( STATE, CLDPARAMS%LS_ICEFALL,     'LS_ICEFALL:',     DEFAULT= 1.0     )
        call MAPL_GetResource( STATE, CLDPARAMS%WRHODEP,        'WRHODEP:',        DEFAULT= 0.5     )
      else
        call MAPL_GetResource( STATE, CLDPARAMS%ICE_SETTLE,     'ICE_SETTLE:',     DEFAULT= 1.      )
        call MAPL_GetResource( STATE, CLDPARAMS%ANV_ICEFALL,    'ANV_ICEFALL:',    DEFAULT= 0.2     )
        call MAPL_GetResource( STATE, CLDPARAMS%LS_ICEFALL,     'LS_ICEFALL:',     DEFAULT= 0.2     )
        call MAPL_GetResource( STATE, CLDPARAMS%WRHODEP,        'WRHODEP:',        DEFAULT= 0.0     )
      endif



      if(adjustl(CLDMICRO)=="2MOMENT") then
         call MAPL_GetResource( STATE, CLDPARAMS%FAC_RI,         'FAC_RI:',         DEFAULT= 1.0     )
         call MAPL_GetResource( STATE, CLDPARAMS%MIN_RI,         'MIN_RI:',         DEFAULT= 15.e-6  )
         call MAPL_GetResource( STATE, CLDPARAMS%MAX_RI,         'MAX_RI:',         DEFAULT= 150.e-6 )
         call MAPL_GetResource( STATE, CLDPARAMS%FAC_RL,         'FAC_RL:',         DEFAULT= 1.0     )
         call MAPL_GetResource( STATE, CLDPARAMS%MIN_RL,         'MIN_RL:',         DEFAULT= 5.e-6   )
         call MAPL_GetResource( STATE, CLDPARAMS%MAX_RL,         'MAX_RL:',         DEFAULT= 21.e-6  )
         call MAPL_GetResource( STATE, CLDPARAMS%PRECIPRAD,      'PRECIPRAD:',      DEFAULT= 1.0     )
         call MAPL_GetResource( STATE, CLDPARAMS%SNOW_REVAP_FAC, 'SNOW_REVAP_FAC:', DEFAULT= 0.5     )
         call MAPL_GetResource( STATE, CLDPARAMS%TURNRHCRIT,     'TURNRHCRIT:',     DEFAULT= 884.0   )
      elseif (adjustl(CLDMICRO) =="GFDL") then
         call MAPL_GetResource( STATE, CLDPARAMS%FAC_RI,         'FAC_RI:',         DEFAULT= 0.5     )
         call MAPL_GetResource( STATE, CLDPARAMS%MIN_RI,         'MIN_RI:',         DEFAULT= 15.e-6  )
         call MAPL_GetResource( STATE, CLDPARAMS%MAX_RI,         'MAX_RI:',         DEFAULT= 150.e-6 )
         call MAPL_GetResource( STATE, CLDPARAMS%FAC_RL,         'FAC_RL:',         DEFAULT= 1.0     )
         call MAPL_GetResource( STATE, CLDPARAMS%MIN_RL,         'MIN_RL:',         DEFAULT= 5.e-6   )
         call MAPL_GetResource( STATE, CLDPARAMS%MAX_RL,         'MAX_RL:',         DEFAULT= 21.e-6  )
         call MAPL_GetResource( STATE, CLDPARAMS%PRECIPRAD,      'PRECIPRAD:',      DEFAULT= 0.0     )
         call MAPL_GetResource( STATE, CLDPARAMS%SNOW_REVAP_FAC, 'SNOW_REVAP_FAC:', DEFAULT= 1.0     ) ! irrelevant
         call MAPL_GetResource( STATE, CLDPARAMS%TURNRHCRIT,     'TURNRHCRIT:',     DEFAULT= 884.0   ) ! irrelevant
      else
         call MAPL_GetResource( STATE, CLDPARAMS%FAC_RI,         'FAC_RI:',         DEFAULT= 1.0     )
         call MAPL_GetResource( STATE, CLDPARAMS%MIN_RI,         'MIN_RI:',         DEFAULT= 15.e-6  )
         call MAPL_GetResource( STATE, CLDPARAMS%MAX_RI,         'MAX_RI:',         DEFAULT= 150.e-6 )
         call MAPL_GetResource( STATE, CLDPARAMS%FAC_RL,         'FAC_RL:',         DEFAULT= 1.0     )
         call MAPL_GetResource( STATE, CLDPARAMS%MIN_RL,         'MIN_RL:',         DEFAULT= 5.e-6   )
         call MAPL_GetResource( STATE, CLDPARAMS%MAX_RL,         'MAX_RL:',         DEFAULT= 21.e-6  )
         call MAPL_GetResource( STATE, CLDPARAMS%PRECIPRAD,      'PRECIPRAD:',      DEFAULT= 0.0     )
         call MAPL_GetResource( STATE, CLDPARAMS%SNOW_REVAP_FAC, 'SNOW_REVAP_FAC:', DEFAULT= 1.0     )
         call MAPL_GetResource( STATE, CLDPARAMS%TURNRHCRIT,     'TURNRHCRIT:',     DEFAULT= 750.0   )
      end if

     call MAPL_GetResource( STATE, CLOUD_CTL%SCLMFDFR,       'SCLMFDFR:',       DEFAULT= 1.0   )
     call MAPL_GetResource( STATE, CLDPARAMS%SCLM_SHW,       'SCLM_SHW:',       DEFAULT= 1.0   )
     
     
      call MAPL_GetResource( STATE, CLDPARAMS%CNV_ENVF,  'CNV_ENVF:',   DEFAULT= 1.0    )
      call MAPL_GetResource( STATE, CLDPARAMS%ANV_ENVF,  'ANV_ENVF:',   DEFAULT= 1.0    )
      call MAPL_GetResource( STATE, CLDPARAMS%SC_ENVF,    'SC_ENVF:',   DEFAULT= 1.0    )
      call MAPL_GetResource( STATE, CLDPARAMS%LS_ENVF,    'LS_ENVF:',   DEFAULT= 1.0     )

      call MAPL_GetResource( STATE, CLDPARAMS%FR_LS_WAT, 'FR_LS_WAT:',  DEFAULT= 1.0    )
      call MAPL_GetResource( STATE, CLDPARAMS%FR_AN_WAT, 'FR_AN_WAT:',  DEFAULT= 1.0    )
      call MAPL_GetResource( STATE, CLDPARAMS%FR_LS_ICE, 'FR_LS_ICE:',  DEFAULT= 0.0    )
      call MAPL_GetResource( STATE, CLDPARAMS%FR_AN_ICE, 'FR_AN_ICE:',  DEFAULT= 0.0    )

   
      call MAPL_GetResource( STATE, CLDPARAMS%CFPBL_EXP,      'CFPBL_EXP:',      DEFAULT= 1 )
      
      call MAPL_GetResource( STATE, CLOUD_CTL%RSUB_RADIUS, 'RSUB_RADIUS:',   DEFAULT= 1.00e-3  )

      call MAPL_TimerOff(STATE,"-MISC2")

      call MAPL_TimerOn (STATE,"-CLOUD")

      ! -----------------------------------------
      ! Preliminary calculations for progno_cloud
      ! -----------------------------------------

      ! GPU These are done here due to the backwards loop for ZET
      !     Done here, the progno_cloud code is completely local.

      ! Calculate QST3 and pass in to progno_cloud
      ! ------------------------------------------

      TEMP_0 = TH1 * PK
      DQST3 = GEOS_DQSAT(TEMP_0, PLO, QSAT=QST3)

      ! Calculate QDDF3 and pass in to progno_cloud=>precip3
      ! ----------------------------------------------------

      DZET(:,:,1:LM) = TH1(:,:,1:LM) * (PKE(:,:,1:LM) - PKE(:,:,0:LM-1)) * MAPL_CP/MAPL_GRAV

      ZET(:,:,LM+1) = 0.0
      DO K = LM, 1, -1
         ZET(:,:,K) = ZET(:,:,K+1)+DZET(:,:,K)
      END DO

      WHERE ( ZET(:,:,1:LM) < 3000. )
         QDDF3 = -( ZET(:,:,1:LM)-3000. ) * ZET(:,:,1:LM) * MASS
      ELSEWHERE
         QDDF3 = 0.
      END WHERE

      VMIP = SUM(QDDF3, 3)
      DO K = 1,LM
         QDDF3(:,:,K) = QDDF3(:,:,K) / VMIP
      END DO

      ! Calculate tempor2d based off of levs925 => tempor in RADCOUPLE
      ! --------------------------------------------------------------

      tempor2d = 0.
      do l = levs925, lm
         where (u1(:,:,l).gt.4.) tempor2d(:,:) = 1.
      end do
   
   
    
         do j = 1, JM
            do i = 1, IM
	     aux1=PLE(i,j,LM)/(287.04*(T(i,j,LM)*(1.+0.608*Q1(i,j,LM)))) ! air_dens (kg m^-3)
	     hfs = -SH  (i,j) ! W m^-2
	     hfl = -EVAP(i,j) ! kg m^-2 s^-1
             aux2= (hfs/MAPL_CP + 0.608*T(i,j,LM)*hfl)/aux1 ! buoyancy flux (h+le)
             aux3= ZLE(I, J,  nint(KPBLIN(I, J)))           ! pbl height (m)
             !-convective velocity scale W* (m/s)
             zws(i,j) = max(0.,0.001-1.5*0.41*MAPL_GRAV*aux2*aux3/T(i,j,LM))
             zws(i,j) = 1.2*zws(i,j)**0.3333 ! m/s	     
   	 enddo; enddo
     
     
      !==========AER_CLOUD===================microphysics "if"===========================
      if(adjustl(CLDMICRO) /="2MOMENT") then

        if(USE_AEROSOL_NN) then
!--kml--- aerosol activation (single-moment uphysics)      
         call MAPL_TimerOn(STATE,"--USE_AEROSOL_NN3")
         call Aer_Actv_1M_interface(IM,JM,LM,N_MODES,T, PLO, ZLO, ZLE, QLCN, QICN, QLLS, QILS &
                                   ,KPBLIN,ZWS,OMEGA, FRLAND ,AeroProps, NACTL,NACTI)
         call MAPL_TimerOff(STATE,"--USE_AEROSOL_NN3")
        ENDIF

       if(adjustl(CLDMICRO) == "GFDL") then

         call MAPL_TimerOn (STATE,"--CLOUD_RUN",RC=STATUS)
         VERIFY_(STATUS)

        ! Temperature (K)
         TEMP = TH1*PK
        ! Delta-Z layer thickness (gfdl expects this to be negative)
       ! DZ = TH1 * (PKE(:,:,0:LM-1) - PKE(:,:,1:LM)) * MAPL_CP/MAPL_GRAV
       ! DZ = ( ZLE(:,:,1:LM)-ZLE(:,:,0:LM-1) )
         DZ = -1.0*DZET
        ! Get cloud nuclei particle numbers
         if (USE_AEROSOL_NN) then
           CFX =100.*PLO*r_air/TEMP !density times conversion factor
           NCPL = NACTL/CFX ! kg-1
           NCPI = NACTI/CFX ! kg-1
         else
           NCPL = 0.
           NCPI = 0.
         endif

         call MAPL_TimerOn(STATE,"---CLDMACRO")
         if (DOCLDMACRO==1) then
        ! Fill LTS and EIS
         call FIND_EIS(TH1, QSS, TEMP, ZLE, PLO, KLCL, IM, JM, LM, LTS, EIS)
        ! Clean up any negative specific humidity before macro+microphysics
         call FILLQ2ZERO2( Q1, MASS, FILLQ)
         REV_CN_X  = 0.0
         REV_AN_X  = 0.0 
         REV_LS_X  = 0.0
         REV_SC_X  = 0.0
         RSU_CN_X  = 0.0
         RSU_AN_X  = 0.0
         RSU_LS_X  = 0.0     
         RSU_SC_X  = 0.0     
         CN_PRC2   = 0.0 
         LS_PRC2   = 0.0
         AN_PRC2   = 0.0 
         SC_PRC2   = 0.0 
         CN_SNR    = 0.0 
         LS_SNR    = 0.0  
         AN_SNR    = 0.0
         SC_SNR    = 0.0
         DTDT_macro=TEMP
         DQVDT_macro=Q1
         DQLDT_macro=QLCN+QLLS
         DQIDT_macro=QICN+QILS
         DQADT_macro=CLCN+CLLS
         PFI_CN_X  = 0.0
         PFI_AN_X  = 0.0
         PFI_LS_X  = 0.0
         PFI_SC_X  = 0.0
         PFL_CN_X  = 0.0
         PFL_AN_X  = 0.0
         PFL_LS_X  = 0.0    
         PFL_SC_X  = 0.0    
         SC_ICE    = 1.0
         QSNOW_CN  = 0.0
         QRAIN_CN  = 0.0
         PFRZ      = 0.0
         CNV_MFD_X   = CNV_MFD     ! needed for cloud fraction
         CNV_DQLDT_X = CNV_DQLDT
         CNV_PRC3_X  = CNV_PRC3
         CNV_UPDF_X  = CNV_UPDF 
        ! Fill CNV_FICE,CNV_NICE,CNV_NDROP assume still 1-moment
         do K=1,LM
          do J=1,JM
           do I=1,IM
             CNV_FICE(I,J,K) = ICE_FRACTION( TEMP(I,J,K), CNV_FRACTION(I,J), SNOMAS(I,J), FRLANDICE(I,J), FRLAND(I,J) )
           enddo
          enddo
         enddo
         IF(ADJUSTL(CONVPAR_OPTION) == 'GF') THEN
        ! GF scheme handles its own conv precipitation
              CNV_PRC3_X  = 0.0
              CNV_NDROP_X = 0.0
              CNV_NICE_X  = 0.0
         END IF
         call  macro_cloud (      &
              IM*JM, LM         , &
              DT_MOIST          , &
              PLO               , &
              CNV_PLE           , &
              PK                , &
              SNOMAS            , &   ! <- surf
              FRLANDICE         , &   ! <- surf
              FRLAND            , &   ! <- surf
              KH                , &   ! <- turb
              CNV_MFD_X         , &   ! <- ras/gf
              CNV_DQLDT_X       , &   ! <- ras/gf
              CNV_PRC3_X        , &   ! <- ras/gf
              CNV_UPDF_X        , &   ! <- ras/gf
              MFD_SC            , &   ! <- shcu   
              QLDET_SC          , &   ! <- shcu   
              QIDET_SC          , &   ! <- shcu   
              SHLW_PRC3         , &   ! <- shcu   
              SHLW_SNO3         , &   ! <- shcu   
              UFRC_SC           , &   ! <- shcu 
              U1                , &
              V1                , & 
              TH1               , &              
              Q1                , &
              QLLS              , &
              QLCN              , &
              QILS              , &
              QICN              , &
              CLCN              , &
              CLLS              , &           
              CN_PRC2           , &            
              CN_ARFX           , &
              CN_SNR            , &
              CLDPARAMS         , &
              CLOUD_CTL%SCLMFDFR, &
              QST3              , &
              DZET              , &
              CNV_FRACTION      ,  &
              QDDF3             , &
                                ! Diagnostics
              RHX_X             , &
              REV_CN_X          , &
              RSU_CN_X          , &
              ACLL_CN_X,ACIL_CN_X   , &
              PFL_CN_X,PFI_CN_X     , &
              DLPDF_X,DIPDF_X,DLFIX_X,DIFIX_X,    &
              DCNVL_X, DCNVi_X,       &
              ALPHT_X, CFPDF_X, &
              DQRL_X            , &
              VFALLSN_CN_X      , &
              VFALLRN_CN_X      , &
              EVAPC_X , SUBLC_X ,  &
                                ! End diagnostics
             !!====2-Moment============
              CNV_FICE          , &
              CNV_NDROP_X       , &
              CNV_NICE_X        , &
              SC_NDROP          , &
              SC_NICE           , &
              SC_ICE            , &
              NCPL              , &
              NCPI              , &
              PFRZ              , &
              DNDCNV_X          , &
              DNCCNV_X          , &
              DT_RASP           , &
              QRAIN_CN          , & !grid av
              QSNOW_CN          , &
              KCBL, LTS,  CONVPAR_OPTION)              
      

        ! Fill DTDT_MACRO diagnostic
         TEMP    = TH1*PK
         DTDT_macro=  (TEMP-DTDT_macro)/DT_MOIST
         DQVDT_macro=(Q1-DQVDT_macro)/DT_MOIST
         DQLDT_macro=((QLCN+QLLS)-DQLDT_macro)/DT_MOIST
         DQIDT_macro=((QICN+QILS)-DQIDT_macro)/DT_MOIST
         DQADT_macro=((CLCN+CLLS)-DQADT_macro)/DT_MOIST
         else
         REV_CN_X  = 0.0
         REV_AN_X  = 0.0
         REV_LS_X  = 0.0
         REV_SC_X  = 0.0
         RSU_CN_X  = 0.0
         RSU_AN_X  = 0.0
         RSU_LS_X  = 0.0
         RSU_SC_X  = 0.0
         PFI_CN_X  = 0.0
         PFI_AN_X  = 0.0
         PFI_LS_X  = 0.0
         PFI_SC_X  = 0.0
         PFL_CN_X  = 0.0
         PFL_AN_X  = 0.0
         PFL_LS_X  = 0.0
         PFL_SC_X  = 0.0
         EVAPC_X = 0.0
         SUBLC_X = 0
         endif
         call MAPL_TimerOff(STATE,"---CLDMACRO")

         call MAPL_TimerOn(STATE,"---GFDL_CLDMICRO")
        ! Cloud
         FQA= 0.0
         RAD_CF = MIN(CLCN+CLLS,1.0)
         FQA =  MIN(1.0,MAX(CLCN/MAX(RAD_CF,1.e-5),0.0))
        ! Liquid
         FQAl = 0.0
         RAD_QL = QLCN+QLLS
         FQAl =  MIN(1.0,MAX(QLCN/MAX(RAD_QL,1.E-8),0.0))
        ! Ice
         FQAi = 0.0
         RAD_QI = QICN+QILS
         FQAi =  MIN(1.0,MAX(QICN/MAX(RAD_QI,1.E-8),0.0))
        ! VAPOR
         RAD_QV = Q1
        ! RAIN
         RAD_QR = QRAIN
        ! SNOW
         RAD_QS = QSNOW
        ! GRAUPEL
         RAD_QG = QGRAUPEL
        ! Vertical velocity
         W1 = W
      !  W1 = -W*(1.0+MAPL_VIREPS*RAD_QV) * TEMP / PLO * (MAPL_RDRY/MAPL_GRAV)

        ! Zero-out microphysics tendencies
         DQVDT_micro = 0.
         DQLDT_micro = 0.
         DQRDT_micro = 0.
         DQIDT_micro = 0.
         DQSDT_micro = 0.
         DQGDT_micro = 0.
         DQADT_micro = 0.
          DUDT_micro = 0.
          DVDT_micro = 0.
          DTDT_micro = 0.
     ! Zero-out 3D Precipitation Fluxes 
        ! Ice
            PFI_LS_X = 0.
        ! Liquid
            PFL_LS_X = 0.

        ! Execute GFDL microphysics
         call gfdl_cloud_microphys_driver( &
                             ! Input water/cloud species and liquid+ice CCN [NACTL+NACTI]
                               RAD_QV, RAD_QL, RAD_QR, RAD_QI, RAD_QS, RAD_QG, RAD_CF, (NACTL+NACTI)/1.e6, &
                             ! Output tendencies
                               DQVDT_micro, DQLDT_micro, DQRDT_micro, DQIDT_micro, &
                               DQSDT_micro, DQGDT_micro, DQADT_micro, DTDT_micro, &
                             ! Input fields
                               TEMP, W1, U1, V1, DUDT_micro, DVDT_micro, DZ, DP, &
                             ! constant inputs
                               AREA, DT_MOIST, FRLAND, CNV_FRACTION, &
                             ! Output rain re-evaporation and sublimation
                               REV_MC_X, RSU_MC_X, & 
                             ! Output precipitates
                               PRCP_RAIN, PRCP_SNOW, PRCP_ICE, PRCP_GRAUPEL, &
                             ! Output mass flux during sedimentation (Pa kg/kg)
                               PFL_LS_X(:,:,1:LM), PFI_LS_X(:,:,1:LM), &
                             ! constant grid/time information
                               LHYDROSTATIC, LPHYS_HYDROSTATIC, &
                               1,IM, 1,JM, 1,LM, 1, LM)

     ! Convert precip diagnostics from mm/day to kg m-2 s-1
         PRCP_RAIN    = PRCP_RAIN    / 86400.
         PRCP_SNOW    = PRCP_SNOW    / 86400.
         PRCP_ICE     = PRCP_ICE     / 86400.
         PRCP_GRAUPEL = PRCP_GRAUPEL / 86400.
     ! Convert precipitation fluxes from (Pa kg/kg) to (kg m-2 s-1)
         PFL_LS_X = PFL_LS_X/(MAPL_GRAV*DT_MOIST)
         PFI_LS_X = PFI_LS_X/(MAPL_GRAV*DT_MOIST)
     ! Redistribute precipitation fluxes for chemistry
         PFL_AN_X(:,:,1:LM) = PFL_LS_X(:,:,1:LM) * FQAl
         PFL_LS_X(:,:,1:LM) = PFL_LS_X(:,:,1:LM) - PFL_AN_X(:,:,1:LM)
         PFI_AN_X(:,:,1:LM) = PFI_LS_X(:,:,1:LM) * FQAi
         PFI_LS_X(:,:,1:LM) = PFI_LS_X(:,:,1:LM) - PFI_AN_X(:,:,1:LM)
     ! Fill precip diagnostics
         LS_PRC2 = PRCP_RAIN
         LS_SNR  = PRCP_SNOW + PRCP_ICE + PRCP_GRAUPEL
     ! Apply tendencies
         TEMP = TEMP + DTDT_micro  * DT_MOIST
         U1   = U1   + DUDT_micro  * DT_MOIST
         V1   = V1   + DVDT_micro  * DT_MOIST
     ! W1 was updated in gfdl_microphsycis, fill WI tendency export
         if (associated(WI)) WI =  (W1 - W)/DT_MOIST
     ! Apply moist/cloud species tendencies
         RAD_QV = RAD_QV + DQVDT_micro * DT_MOIST
         RAD_QL = RAD_QL + DQLDT_micro * DT_MOIST
         RAD_QR = RAD_QR + DQRDT_micro * DT_MOIST
         RAD_QI = RAD_QI + DQIDT_micro * DT_MOIST
         RAD_QS = RAD_QS + DQSDT_micro * DT_MOIST
         RAD_QG = RAD_QG + DQGDT_micro * DT_MOIST
         RAD_CF = RAD_CF + DQADT_micro * DT_MOIST
     ! when do_qa=.true. in GFDL_MP RAD_CF is update internally and DQADT_micro is zero
     ! so lets be sure we get the real cloud tendency from micro here
         DQADT_micro = ( RAD_CF - CLCN - CLLS ) / DT_MOIST
     ! Cloud liquid & Ice tendencies (these exports are confusing, for now keep them zeros)
         REV_LS_X = REV_LS_X + REV_MC_X
         RSU_LS_X = RSU_LS_X + RSU_MC_X
    !    EVAPC_X = ( EVAPC_X - RAD_QL ) / DT_MOIST
    !    SUBLC_X = ( SUBLC_X - RAD_QI ) / DT_MOIST
     ! Fill vapor/rain/snow/graupel state
         Q1       = RAD_QV
         QRAIN    = RAD_QR
         QSNOW    = RAD_QS
         QGRAUPEL = RAD_QG
         if (associated(QRTOT)) QRTOT = QRAIN
         if (associated(QSTOT)) QSTOT = QSNOW
     ! Redistribute cloud/liquid/ice species...
         CLCN = RAD_CF*     FQA
         CLLS = RAD_CF*(1.0-FQA)
         QLCN = RAD_QL*     FQAl
         QLLS = RAD_QL*(1.0-FQAl)
         QICN = RAD_QI*     FQAi
         QILS = RAD_QI*(1.0-FQAi)
     ! Clean up clouds after microphysics
         CALL FIX_UP_CLOUDS( TEMP, Q1, QLLS, QILS, CLLS, QLCN, QICN, CLCN )
     ! Clean up any negative specific humidity
         CALL FILLQ2ZERO2( Q1, MASS, FILLQ  )
     ! Convert back to PT
         TH1 = TEMP/PK
     ! Radiation Coupling
      if (CLDPARAMS%DISABLE_RAD==1) then
               RAD_QL     = 0.
               RAD_QI     = 0.
               RAD_QR     = 0.
               RAD_QS     = 0.
               RAD_QG     = 0.
               RAD_CF     = 0.
               CLDREFFL   = 0.
               CLDREFFI   = 0.
      else
         do K = 1, LM
           do J = 1, JM
             do I = 1, IM
               RHX_X(I,J,K) = Q1(I,J,K)/GEOS_QSAT( TEMP(I,J,K), PLO(I,J,K) )
               call RADCOUPLE ( TEMP(I,J,K), PLO(I,J,K), CLLS(I,J,K), CLCN(I,J,K), &
                     Q1(I,J,K), QLLS(I,J,K), QILS(I,J,K), QLCN(I,J,K), QICN(I,J,K), QRAIN(I,J,K), QSNOW(I,J,K), NACTL(I,J,K), NACTI(I,J,K), &
                     RAD_QV(I,J,K), RAD_QL(I,J,K), RAD_QI(I,J,K), RAD_QR(I,J,K), RAD_QS(I,J,K), RAD_CF(I,J,K), &
                     CLDREFFL(I,J,K), CLDREFFI(I,J,K), FRLAND(I,J), CNV_FRACTION(I,J), INT(CLDPARAMS%FR_AN_WAT), & 
                     CLDPARAMS%FAC_RL, CLDPARAMS%MIN_RL, CLDPARAMS%MAX_RL, CLDPARAMS%FAC_RI, CLDPARAMS%MIN_RI, CLDPARAMS%MAX_RI, &
                     RHX(I,J,K) )
            enddo
          enddo
        enddo
        RAD_QG = 0.0
      endif

         if (USE_AEROSOL_NN) then
           CFX =100.*PLO*r_air/TEMP !density times conversion factor
           NCPL = NACTL/CFX ! kg-1
           NCPI = NACTI/CFX ! kg-1
         else
           NCPL = 0.
           NCPI = 0.
         endif
       ! Exports required
         CFLIQ  = 0.0
         CFICE  = 0.0
         if(CLDPARAMS%PRECIPRAD.eq.0.) then
           QTOT   = QICN+QILS+QLCN+QLLS
           QL_TOT = QLCN+QLLS
           QI_TOT = QICN+QILS
         else
           QTOT   = QICN+QILS+QLCN+QLLS+QRAIN+QSNOW+QGRAUPEL
           QL_TOT = QLCN+QLLS+QRAIN
           QI_TOT = QICN+QILS+QSNOW+QGRAUPEL
         endif
         WHERE (QTOT .gt. 1.0e-12)
            CFLIQ=RAD_CF*QL_TOT/QTOT
            CFICE=RAD_CF*QI_TOT/QTOT
         END WHERE
         CFLIQ=MAX(MIN(CFLIQ, 1.0), 0.0)
         CFICE=MAX(MIN(CFICE, 1.0), 0.0)
         where (QI_TOT .le. 0.0)
            CFICE =0.0
            NCPI=0.0
            CLDREFFI = CLDPARAMS%MIN_RI
         end where

         where (QL_TOT .le. 0.0)
            CFLIQ =0.0
            NCPL  =0.0
            CLDREFFL = CLDPARAMS%MIN_RL
         end where

         call MAPL_TimerOff(STATE,"---GFDL_CLDMICRO",RC=STATUS)
         VERIFY_(STATUS)
         call MAPL_TimerOff(STATE,"--CLOUD_RUN",RC=STATUS)
         VERIFY_(STATUS)

     ! endif
       else !===== 1-moment microphysics
       
         if (associated(DQVDT_micro)) DQVDT_micro = Q1
         if (associated(DQLDT_micro)) DQLDT_micro = QLLS + QLCN
         if (associated(DQIDT_micro)) DQIDT_micro = QILS + QICN
         if (associated(DQRDT_micro)) DQRDT_micro = 0.0
         if (associated(DQSDT_micro)) DQSDT_micro = 0.0
         if (associated(DQGDT_micro)) DQGDT_micro = 0.0
         if (associated(DQADT_micro)) DQADT_micro = CLLS + CLCN
         if (associated(DUDT_micro) ) DUDT_micro  = U1
         if (associated(DVDT_micro) ) DVDT_micro  = V1
         if (associated(DTDT_micro) ) DTDT_micro  = TH1*PK
         if (associated(DTDT_macro) ) DTDT_macro  = 0.0
         if (associated(DQVDT_macro)) DQVDT_macro  = 0.0
         if (associated(DQLDT_macro)) DQLDT_macro  = 0.0
         if (associated(DQIDT_macro)) DQIDT_macro  = 0.0
         if (associated(DQADT_macro)) DQADT_macro  = 0.0
 
#ifdef _CUDA

         call MAPL_GetResource(STATE,BLOCKSIZE,'BLOCKSIZE:',DEFAULT=128,RC=STATUS)
         VERIFY_(STATUS)

         Block = dim3(blocksize,1,1)
         Grid = dim3(ceiling(real(IM*JM)/real(blocksize)),1,1)

         ! GPU In order to save time on the GPU, we need to know the logic of all the 
         ! diagnostics that rely on these. When we integrate the code after PROGNO_CLOUD this will be redundant.
         COPY_RSU_CN = ASSOCIATED(RSU_CN) .OR. ASSOCIATED(REVSU_CN)   .OR. ASSOCIATED(PREVTOT) .OR. ASSOCIATED(SUBPZ)
         COPY_RSU_AN = ASSOCIATED(RSU_AN) .OR. ASSOCIATED(REVSU_LSAN) .OR. ASSOCIATED(PREVTOT) .OR. ASSOCIATED(SUBPZ)
         COPY_RSU_LS = ASSOCIATED(RSU_LS) .OR. ASSOCIATED(REVSU_LSAN) .OR. ASSOCIATED(PREVTOT) .OR. ASSOCIATED(SUBPZ)
         COPY_RSU_SC = ASSOCIATED(RSU_SC) .OR. ASSOCIATED(REVSU_SC) .OR. ASSOCIATED(PREVTOT) .OR. ASSOCIATED(SUBPZ)

         COPY_REV_CN = ASSOCIATED(REV_CN) .OR. ASSOCIATED(REVSU_CN)   .OR. ASSOCIATED(PREVTOT) .OR. ASSOCIATED(EVPPZ)
         COPY_REV_SC = ASSOCIATED(REV_SC) .OR. ASSOCIATED(REVSU_SC)   .OR. ASSOCIATED(PREVTOT) .OR. ASSOCIATED(EVPPZ)
         COPY_REV_AN = ASSOCIATED(REV_AN) .OR. ASSOCIATED(REVSU_LSAN) .OR. ASSOCIATED(PREVTOT) .OR. ASSOCIATED(EVPPZ)
         COPY_REV_LS = ASSOCIATED(REV_LS) .OR. ASSOCIATED(REVSU_LSAN) .OR. ASSOCIATED(PREVTOT) .OR. ASSOCIATED(EVPPZ)

         COPY_ACIL_CN = ASSOCIATED(ACIL_CN) .OR. ASSOCIATED(ACR_TOT) .OR. ASSOCIATED(PGENTOT) .OR. ASSOCIATED(COLLIZ)
         COPY_ACIL_SC = ASSOCIATED(ACIL_SC) .OR. ASSOCIATED(ACR_TOT) .OR. ASSOCIATED(PGENTOT) .OR. ASSOCIATED(COLLIZ)
         COPY_ACIL_AN = ASSOCIATED(ACIL_AN) .OR. ASSOCIATED(ACR_TOT) .OR. ASSOCIATED(PGENTOT) .OR. ASSOCIATED(COLLIZ)
         COPY_ACIL_LS = ASSOCIATED(ACIL_LS) .OR. ASSOCIATED(ACR_TOT) .OR. ASSOCIATED(PGENTOT) .OR. ASSOCIATED(COLLIZ)

         COPY_ACLL_CN = ASSOCIATED(ACLL_CN) .OR. ASSOCIATED(ACR_TOT) .OR. ASSOCIATED(PGENTOT) .OR. ASSOCIATED(COLLLZ)
         COPY_ACLL_SC = ASSOCIATED(ACLL_SC) .OR. ASSOCIATED(ACR_TOT) .OR. ASSOCIATED(PGENTOT) .OR. ASSOCIATED(COLLLZ)
         COPY_ACLL_AN = ASSOCIATED(ACLL_AN) .OR. ASSOCIATED(ACR_TOT) .OR. ASSOCIATED(PGENTOT) .OR. ASSOCIATED(COLLLZ)
         COPY_ACLL_LS = ASSOCIATED(ACLL_LS) .OR. ASSOCIATED(ACR_TOT) .OR. ASSOCIATED(PGENTOT) .OR. ASSOCIATED(COLLLZ)

         COPY_PFI_CN = ASSOCIATED(PFI_CN)
         COPY_PFI_SC = ASSOCIATED(PFI_SC)
         COPY_PFI_AN = ASSOCIATED(PFI_AN) .OR. ASSOCIATED(PFI_LSAN)
         COPY_PFI_LS = ASSOCIATED(PFI_LS) .OR. ASSOCIATED(PFI_LSAN)

         COPY_PFL_CN = ASSOCIATED(PFL_CN)
         COPY_PFL_SC = ASSOCIATED(PFL_SC)
         COPY_PFL_AN = ASSOCIATED(PFL_AN) .OR. ASSOCIATED(PFL_LSAN)
         COPY_PFL_LS = ASSOCIATED(PFL_LS) .OR. ASSOCIATED(PFL_LSAN)

         COPY_DLPDF = ASSOCIATED(DLPDF) .OR. ASSOCIATED(PDFLZ)
         COPY_DIPDF = ASSOCIATED(DIPDF) .OR. ASSOCIATED(PDFIZ)
         COPY_DLFIX = ASSOCIATED(DLFIX) .OR. ASSOCIATED(EVPCZ)
         COPY_DIFIX = ASSOCIATED(DIFIX) .OR. ASSOCIATED(SUBCZ)

         COPY_RHX    = ASSOCIATED(RHX)
         COPY_AUT    = ASSOCIATED(AUT)    .OR. ASSOCIATED(AUTZ)
         COPY_EVAPC  = ASSOCIATED(EVAPC)  .OR. ASSOCIATED(EVPCZ)
         COPY_SDM    = ASSOCIATED(SDM)    .OR. ASSOCIATED(SDMZ)
         COPY_SUBLC  = ASSOCIATED(SUBLC)  .OR. ASSOCIATED(SUBCZ)
         COPY_FRZ_TT = ASSOCIATED(FRZ_TT) .OR. ASSOCIATED(FRZCZ)
         COPY_FRZ_PP = ASSOCIATED(FRZ_PP) .OR. ASSOCIATED(FRZPZ)
         COPY_DCNVL  = ASSOCIATED(DCNVL)  .OR. ASSOCIATED(CNVLZ)
         COPY_DCNVI  = ASSOCIATED(DCNVI)  .OR. ASSOCIATED(CNVIZ)

         COPY_ALPHT = ASSOCIATED(ALPHT)
         COPY_ALPH1 = ASSOCIATED(ALPH1)
         COPY_ALPH2 = ASSOCIATED(ALPH2)

         COPY_CFPDF = ASSOCIATED(CFPDF)
         COPY_RHCLR = ASSOCIATED(RHCLR)
         COPY_DQRL  = ASSOCIATED(DQRL) .OR. ASSOCIATED(PGENTOT)

         COPY_VFALLICE_AN = ASSOCIATED(VFALLICE_AN)
         COPY_VFALLICE_LS = ASSOCIATED(VFALLICE_LS)
         COPY_VFALLWAT_AN = ASSOCIATED(VFALLWAT_AN)
         COPY_VFALLWAT_LS = ASSOCIATED(VFALLWAT_LS)

         COPY_VFALLRN_AN = ASSOCIATED(VFALLRN_AN)
         COPY_VFALLRN_LS = ASSOCIATED(VFALLRN_LS)
         COPY_VFALLRN_CN = ASSOCIATED(VFALLRN_CN)
         COPY_VFALLRN_SC = ASSOCIATED(VFALLRN_SC)
         COPY_VFALLSN_AN = ASSOCIATED(VFALLSN_AN)
         COPY_VFALLSN_LS = ASSOCIATED(VFALLSN_LS)
         COPY_VFALLSN_CN = ASSOCIATED(VFALLSN_CN)
         COPY_VFALLSN_SC = ASSOCIATED(VFALLSN_SC)

!!$         COPY_LIQANMOVE = ASSOCIATED(LIQANMOVE)
!!$         COPY_ICEANMOVE = ASSOCIATED(ICEANMOVE)
!!$         COPY_DANCLD = ASSOCIATED(DANCLD)
!!$         COPY_DLSCLD = ASSOCIATED(DLSCLD)
!!$         COPY_CURAINMOVE = ASSOCIATED(CURAINMOVE)
!!$         COPY_CUSNOWMOVE = ASSOCIATED(CUSNOWMOVE)

         call MAPL_TimerOn (STATE,"--CLOUD_ALLOC",RC=STATUS)
         VERIFY_(STATUS)

         ! ----------------------
         ! Allocate device arrays
         ! ----------------------

         ! Inputs
         ! ------

         ALLOCATE(PP_DEV(IM*JM,LM), STAT=STATUS)
         VERIFY_(STATUS)
         ALLOCATE(EXNP_DEV(IM*JM,LM), STAT=STATUS)
         VERIFY_(STATUS)
         ALLOCATE(PPE_DEV(IM*JM,0:LM), STAT=STATUS)
         VERIFY_(STATUS)
         ALLOCATE(KH_DEV(IM*JM,0:LM), STAT=STATUS)
         VERIFY_(STATUS)
         ALLOCATE(SNOMAS_DEV(IM*JM), STAT=STATUS)
         VERIFY_(STATUS)
         ALLOCATE(FRLAND_DEV(IM*JM), STAT=STATUS)
         VERIFY_(STATUS)
         ALLOCATE(FRLANDICE_DEV(IM*JM), STAT=STATUS)
         VERIFY_(STATUS)
         ALLOCATE(RMFDTR_DEV(IM*JM,LM), STAT=STATUS)
         VERIFY_(STATUS)
         ALLOCATE(SMFDTR_DEV(IM*JM,LM), STAT=STATUS)
         VERIFY_(STATUS)
         ALLOCATE(QLWDTR_DEV(IM*JM,LM), STAT=STATUS)
         VERIFY_(STATUS)
         ALLOCATE(SQLWDTR_DEV(IM*JM,LM), STAT=STATUS)
         VERIFY_(STATUS)
         ALLOCATE(SQIWDTR_DEV(IM*JM,LM), STAT=STATUS)
         VERIFY_(STATUS)
         ALLOCATE(U_DEV(IM*JM,LM), STAT=STATUS)
         VERIFY_(STATUS)
         ALLOCATE(V_DEV(IM*JM,LM), STAT=STATUS)
         VERIFY_(STATUS)
         ALLOCATE(QST3_DEV(IM*JM,LM), STAT=STATUS)
         VERIFY_(STATUS)
         ALLOCATE(DZET_DEV(IM*JM,LM), STAT=STATUS)
         VERIFY_(STATUS)
         ALLOCATE(QDDF3_DEV(IM*JM,LM), STAT=STATUS)
         VERIFY_(STATUS)
         ALLOCATE(TEMPOR_DEV(IM*JM), STAT=STATUS)
         VERIFY_(STATUS)
         ALLOCATE(CNV_FRACTION_DEV(IM*JM), STAT=STATUS)
         VERIFY_(STATUS) 
         ALLOCATE(TROPP_DEV(IM*JM), STAT=STATUS)
         VERIFY_(STATUS)

         ! Inoutputs
         ! ---------

         ALLOCATE(TH_DEV(IM*JM,LM), STAT=STATUS)
         VERIFY_(STATUS)
         ALLOCATE(Q_DEV(IM*JM,LM), STAT=STATUS)
         VERIFY_(STATUS)
         ALLOCATE(QRN_CU_DEV(IM*JM,LM), STAT=STATUS)
         VERIFY_(STATUS)
         ALLOCATE(QRN_SC_DEV(IM*JM,LM), STAT=STATUS)
         VERIFY_(STATUS)
         ALLOCATE(QSN_SC_DEV(IM*JM,LM), STAT=STATUS)
         VERIFY_(STATUS)
         ALLOCATE(CNV_UPDFRC_DEV(IM*JM,LM), STAT=STATUS)
         VERIFY_(STATUS)
         ALLOCATE(SC_UPDFRC_DEV(IM*JM,LM), STAT=STATUS)
         VERIFY_(STATUS)
         ALLOCATE(QLW_LS_DEV(IM*JM,LM), STAT=STATUS)
         VERIFY_(STATUS)
         ALLOCATE(QLW_AN_DEV(IM*JM,LM), STAT=STATUS)
         VERIFY_(STATUS)
         ALLOCATE(QIW_LS_DEV(IM*JM,LM), STAT=STATUS)
         VERIFY_(STATUS)
         ALLOCATE(QIW_AN_DEV(IM*JM,LM), STAT=STATUS)
         VERIFY_(STATUS)
         ALLOCATE(ANVFRC_DEV(IM*JM,LM), STAT=STATUS)
         VERIFY_(STATUS)
         ALLOCATE(CLDFRC_DEV(IM*JM,LM), STAT=STATUS)
         VERIFY_(STATUS)

         ! Outputs
         ! -------

         ALLOCATE(RAD_CLDFRC_DEV(IM*JM,LM), STAT=STATUS)
         VERIFY_(STATUS)
         ALLOCATE(RAD_QL_DEV(IM*JM,LM), STAT=STATUS)
         VERIFY_(STATUS)
         ALLOCATE(RAD_QI_DEV(IM*JM,LM), STAT=STATUS)
         VERIFY_(STATUS)
         ALLOCATE(RAD_QR_DEV(IM*JM,LM), STAT=STATUS)
         VERIFY_(STATUS)
         ALLOCATE(RAD_QS_DEV(IM*JM,LM), STAT=STATUS)
         VERIFY_(STATUS)
         ALLOCATE(CLDREFFL_DEV(IM*JM,LM), STAT=STATUS)
         VERIFY_(STATUS)
         ALLOCATE(CLDREFFI_DEV(IM*JM,LM), STAT=STATUS)
         VERIFY_(STATUS)
         ALLOCATE(QPLS_DEV(IM*JM,LM), STAT=STATUS)
         VERIFY_(STATUS)
         ALLOCATE(PRELS_DEV(IM*JM), STAT=STATUS)
         VERIFY_(STATUS)
         ALLOCATE(PRECU_DEV(IM*JM), STAT=STATUS)
         VERIFY_(STATUS)
         ALLOCATE(PREAN_DEV(IM*JM), STAT=STATUS)
         VERIFY_(STATUS)
         ALLOCATE(LSARF_DEV(IM*JM), STAT=STATUS)
         VERIFY_(STATUS)
         ALLOCATE(CUARF_DEV(IM*JM), STAT=STATUS)
         VERIFY_(STATUS)
         ALLOCATE(ANARF_DEV(IM*JM), STAT=STATUS)
         VERIFY_(STATUS)
         ALLOCATE(SNRLS_DEV(IM*JM), STAT=STATUS)
         VERIFY_(STATUS)
         ALLOCATE(SNRCU_DEV(IM*JM), STAT=STATUS)
         VERIFY_(STATUS)
         ALLOCATE(SNRAN_DEV(IM*JM), STAT=STATUS)
         VERIFY_(STATUS)
         ALLOCATE(SNRSC_DEV(IM*JM), STAT=STATUS)
         VERIFY_(STATUS)

         ! Because of use_autoconv_timescale in PROGNO_CLOUD, the PFL and PFI
         ! arrays can "silently" be used as working arrays
         ! They must always be allocated
         ALLOCATE(PFL_CN_DEV(IM*JM,0:LM), STAT=STATUS)
         VERIFY_(STATUS)
         ALLOCATE(PFI_CN_DEV(IM*JM,0:LM), STAT=STATUS)
         VERIFY_(STATUS)
         ALLOCATE(PFL_AN_DEV(IM*JM,0:LM), STAT=STATUS)
         VERIFY_(STATUS)
         ALLOCATE(PFI_AN_DEV(IM*JM,0:LM), STAT=STATUS)
         VERIFY_(STATUS)
         ALLOCATE(PFL_LS_DEV(IM*JM,0:LM), STAT=STATUS)
         VERIFY_(STATUS)
         ALLOCATE(PFI_LS_DEV(IM*JM,0:LM), STAT=STATUS)
         VERIFY_(STATUS)
         ALLOCATE(PFL_SC_DEV(IM*JM,0:LM), STAT=STATUS)
         VERIFY_(STATUS)
         ALLOCATE(PFI_SC_DEV(IM*JM,0:LM), STAT=STATUS)
         VERIFY_(STATUS)

         ! Diagnostics
         ! -----------

         ALLOCATE(RHX_DEV(IM*JM,LM), STAT=STATUS)
         ALLOCATE(REV_LS_DEV(IM*JM,LM), STAT=STATUS)
         ALLOCATE(REV_AN_DEV(IM*JM,LM), STAT=STATUS)
         ALLOCATE(REV_CN_DEV(IM*JM,LM), STAT=STATUS)
         ALLOCATE(REV_SC_DEV(IM*JM,LM), STAT=STATUS)
         ALLOCATE(RSU_LS_DEV(IM*JM,LM), STAT=STATUS)
         ALLOCATE(RSU_AN_DEV(IM*JM,LM), STAT=STATUS)
         ALLOCATE(RSU_CN_DEV(IM*JM,LM), STAT=STATUS)
         ALLOCATE(RSU_SC_DEV(IM*JM,LM), STAT=STATUS)
         ALLOCATE(ACLL_CN_DEV(IM*JM,LM), STAT=STATUS)
         ALLOCATE(ACIL_CN_DEV(IM*JM,LM), STAT=STATUS)
         ALLOCATE(ACLL_SC_DEV(IM*JM,LM), STAT=STATUS)
         ALLOCATE(ACIL_SC_DEV(IM*JM,LM), STAT=STATUS)
         ALLOCATE(ACLL_AN_DEV(IM*JM,LM), STAT=STATUS)
         ALLOCATE(ACIL_AN_DEV(IM*JM,LM), STAT=STATUS)
         ALLOCATE(ACLL_LS_DEV(IM*JM,LM), STAT=STATUS)
         ALLOCATE(ACIL_LS_DEV(IM*JM,LM), STAT=STATUS)
         ALLOCATE(PDFL_DEV(IM*JM,LM), STAT=STATUS)
         ALLOCATE(PDFI_DEV(IM*JM,LM), STAT=STATUS)
         ALLOCATE(FIXL_DEV(IM*JM,LM), STAT=STATUS)
         ALLOCATE(FIXI_DEV(IM*JM,LM), STAT=STATUS)
         ALLOCATE(AUT_DEV(IM*JM,LM), STAT=STATUS)
         ALLOCATE(EVAPC_DEV(IM*JM,LM), STAT=STATUS)
         ALLOCATE(SDM_DEV(IM*JM,LM), STAT=STATUS)
         ALLOCATE(SUBLC_DEV (IM*JM,LM), STAT=STATUS)
         ALLOCATE(FRZ_TT_DEV(IM*JM,LM), STAT=STATUS)
         ALLOCATE(FRZ_PP_DEV(IM*JM,LM), STAT=STATUS)
         ALLOCATE(DCNVL_DEV(IM*JM,LM), STAT=STATUS)
         ALLOCATE(DCNVI_DEV(IM*JM,LM), STAT=STATUS)
         ALLOCATE(ALPHT_DEV(IM*JM,LM), STAT=STATUS)
         ALLOCATE(ALPH1_DEV(IM*JM,LM), STAT=STATUS)
         ALLOCATE(ALPH2_DEV(IM*JM,LM), STAT=STATUS)
         ALLOCATE(CFPDF_DEV(IM*JM,LM), STAT=STATUS)
         ALLOCATE(RHCLR_DEV(IM*JM,LM), STAT=STATUS)
         ALLOCATE(DQRL_DEV(IM*JM,LM), STAT=STATUS)
         ALLOCATE(VFALLICE_AN_DEV(IM*JM,LM), STAT=STATUS)
         ALLOCATE(VFALLICE_LS_DEV(IM*JM,LM), STAT=STATUS)
         ALLOCATE(VFALLWAT_AN_DEV(IM*JM,LM), STAT=STATUS)
         ALLOCATE(VFALLWAT_LS_DEV(IM*JM,LM), STAT=STATUS)
         ALLOCATE(VFALLSN_AN_DEV(IM*JM,LM), STAT=STATUS)
         ALLOCATE(VFALLSN_LS_DEV(IM*JM,LM), STAT=STATUS)
         ALLOCATE(VFALLSN_CN_DEV(IM*JM,LM), STAT=STATUS)
         ALLOCATE(VFALLSN_SC_DEV(IM*JM,LM), STAT=STATUS)
         ALLOCATE(VFALLRN_AN_DEV(IM*JM,LM), STAT=STATUS)
         ALLOCATE(VFALLRN_LS_DEV(IM*JM,LM), STAT=STATUS)
         ALLOCATE(VFALLRN_CN_DEV(IM*JM,LM), STAT=STATUS)
!!$         ALLOCATE(LIQANMOVE_DEV(IM*JM,LM), STAT=STATUS)
!!$         ALLOCATE(ICEANMOVE_DEV(IM*JM,LM), STAT=STATUS)
!!$         ALLOCATE(DANCLD_DEV(IM*JM,LM), STAT=STATUS)
!!$         ALLOCATE(DLSCLD_DEV(IM*JM,LM), STAT=STATUS)
!!$         ALLOCATE(CURAINMOVE_DEV(IM*JM,LM), STAT=STATUS)
!!$         ALLOCATE(CUSNOWMOVE_DEV(IM*JM,LM), STAT=STATUS)

         call MAPL_TimerOff(STATE,"--CLOUD_ALLOC",RC=STATUS)
         VERIFY_(STATUS)

         ! Copy host inputs to device
         ! --------------------------

         call MAPL_TimerOn (STATE,"--CLOUD_DATA",RC=STATUS)
         VERIFY_(STATUS)

         call MAPL_TimerOn (STATE,"---CLOUD_DATA_DEVICE",RC=STATUS)
         VERIFY_(STATUS)

         STATUS = cudaMemcpy(PP_DEV,PLO,IM*JM*LM) 
         STATUS = cudaMemcpy(PPE_DEV,CNV_PLE,IM*JM*(LM+1))
         STATUS = cudaMemcpy(EXNP_DEV,PK,IM*JM*LM) 
         STATUS = cudaMemcpy(SNOMAS_DEV,SNOMAS,IM*JM)
         STATUS = cudaMemcpy(FRLANDICE_DEV,FRLANDICE,IM*JM)
         STATUS = cudaMemcpy(FRLAND_DEV,FRLAND,IM*JM) 
         STATUS = cudaMemcpy(KH_DEV,KH,IM*JM*(LM+1))
         STATUS = cudaMemcpy(RMFDTR_DEV,CNV_MFD,IM*JM*LM) 
         STATUS = cudaMemcpy(SMFDTR_DEV,MFD_SC,IM*JM*LM) 
         STATUS = cudaMemcpy(QLWDTR_DEV,CNV_DQLDT,IM*JM*LM) 
         STATUS = cudaMemcpy(SQLWDTR_DEV,QLDET_SC,IM*JM*LM) 
         STATUS = cudaMemcpy(QQIWDTR_DEV,QIDET_SC,IM*JM*LM) 
         STATUS = cudaMemcpy(U_DEV,U1,IM*JM*LM) 
         STATUS = cudaMemcpy(V_DEV,V1,IM*JM*LM) 
         STATUS = cudaMemcpy(QST3_DEV,QST3,IM*JM*LM) 
         STATUS = cudaMemcpy(DZET_DEV,DZET,IM*JM*LM) 
         STATUS = cudaMemcpy(QDDF3_DEV,QDDF3,IM*JM*LM) 
         STATUS = cudaMemcpy(TEMPOR_DEV,TEMPOR2D,IM*JM) 
         STATUS = cudaMemcpy(CNV_FRACTION_DEV,CNV_FRACTION,IM*JM)  
         STATUS = cudaMemcpy(TROPP_DEV,TROPP,IM*JM) 

         ! Inoutputs
         ! ---------

         STATUS = cudaMemcpy(QRN_CU_DEV,CNV_PRC3,IM*JM*LM) 
         STATUS = cudaMemcpy(QRN_SC_DEV,SHLW_PRC3,IM*JM*LM) 
         STATUS = cudaMemcpy(QSN_SC_DEV,SHLW_SNO3,IM*JM*LM) 
         STATUS = cudaMemcpy(CNV_UPDFRC_DEV,CNV_UPDF,IM*JM*LM) 
         STATUS = cudaMemcpy(SC_UPDFRC_DEV,UFRC_SC,IM*JM*LM) 
         STATUS = cudaMemcpy(TH_DEV,TH1,IM*JM*LM) 
         STATUS = cudaMemcpy(Q_DEV,Q1,IM*JM*LM) 
         STATUS = cudaMemcpy(QLW_LS_DEV,QLLS,IM*JM*LM) 
         STATUS = cudaMemcpy(QLW_AN_DEV,QLCN,IM*JM*LM) 
         STATUS = cudaMemcpy(QIW_LS_DEV,QILS,IM*JM*LM) 
         STATUS = cudaMemcpy(QIW_AN_DEV,QICN,IM*JM*LM) 
         STATUS = cudaMemcpy(ANVFRC_DEV,CLCN,IM*JM*LM) 
         STATUS = cudaMemcpy(CLDFRC_DEV,CLLS,IM*JM*LM) 

         call MAPL_TimerOff(STATE,"---CLOUD_DATA_DEVICE",RC=STATUS)
         VERIFY_(STATUS)

         call MAPL_TimerOn (STATE,"---CLOUD_DATA_CONST",RC=STATUS)
         VERIFY_(STATUS)

         ! ------------------------
         ! Load constants in device
         ! ------------------------

         ! CLDPARAMS Constants are loaded into constant memory
         ! ---------------------------------------------------

         CNV_BETA      = CLDPARAMS%CNV_BETA  ! Area factor for convective rain showers (non-dim)
         ANV_BETA      = CLDPARAMS%ANV_BETA  ! Area factor for anvil rain showers (non-dim)
         LS_BETA       = CLDPARAMS%LS_BETA   ! Area factor for Large Scale rain showers (non-dim)
         RH00          = CLDPARAMS%RH_CRIT   ! Critical relative humidity
         C_00          = CLDPARAMS%AUTOC_LS
         LWCRIT        = CLDPARAMS%QC_CRIT_LS
         C_ACC         = CLDPARAMS%ACCRETION
         C_EV_R        = CLDPARAMS%RAIN_REVAP_FAC
         C_EV_S        = CLDPARAMS%SNOW_REVAP_FAC
         CLDVOL2FRC    = CLDPARAMS%VOL_TO_FRAC
         RHSUP_ICE     = CLDPARAMS%SUPERSAT
         SHR_EVAP_FAC  = CLDPARAMS%SHEAR_EVAP_FAC
         MIN_CLD_WATER = CLDPARAMS%MIN_ALLOW_CCW
         CLD_EVP_EFF   = CLDPARAMS%CCW_EVAP_EFF
         NSMAX         = INT( CLDPARAMS%NSUB_AUTOCONV  )
         LS_SDQV2      = CLDPARAMS%LS_SUND_INTER
         LS_SDQV3      = CLDPARAMS%LS_SUND_COLD
         LS_SDQVT1     = CLDPARAMS%LS_SUND_TEMP1
         ANV_SDQV2     = CLDPARAMS%ANV_SUND_INTER
         ANV_SDQV3     = CLDPARAMS%ANV_SUND_COLD
         ANV_SDQVT1    = CLDPARAMS%ANV_SUND_TEMP1
         ANV_TO_LS     = CLDPARAMS%ANV_TO_LS_TIME
         N_WARM        = CLDPARAMS%NCCN_WARM
         N_ICE         = CLDPARAMS%NCCN_ICE
         N_ANVIL       = CLDPARAMS%NCCN_ANVIL
         N_PBL         = CLDPARAMS%NCCN_PBL
         DISABLE_RAD   = INT( CLDPARAMS%DISABLE_RAD )
         ICE_SETTLE    = NINT( CLDPARAMS%ICE_SETTLE )
         ANV_ICEFALL_C = CLDPARAMS%ANV_ICEFALL
         LS_ICEFALL_C  = CLDPARAMS%LS_ICEFALL
         REVAP_OFF_P   = CLDPARAMS%REVAP_OFF_P
         CNVENVFC      = CLDPARAMS%CNV_ENVF
         ANVENVFC      = CLDPARAMS%ANV_ENVF
         SCENVFC       = CLDPARAMS%SC_ENVF
         LSENVFC       = CLDPARAMS%LS_ENVF
         WRHODEP       = CLDPARAMS%WRHODEP 
         T_ICE_ALL     = CLDPARAMS%ICE_RAMP + MAPL_TICE
         CNVICEPARAM   = CLDPARAMS%CNV_ICEPARAM
         ICEFRPWR      = INT( CLDPARAMS%CNV_ICEFRPWR + .001 )
         CNVDDRFC      = CLDPARAMS%CNV_DDRF
         ANVDDRFC      = CLDPARAMS%ANV_DDRF
         LSDDRFC       = CLDPARAMS%LS_DDRF
         TANHRHCRIT    = INT( CLDPARAMS%TANHRHCRIT )
         MINRHCRIT     = CLDPARAMS%MINRHCRIT
         MAXRHCRIT     = CLDPARAMS%MAXRHCRIT
         TURNRHCRIT    = CLDPARAMS%TURNRHCRIT
         MAXRHCRITLAND = CLDPARAMS%MAXRHCRITLAND
         FR_LS_WAT     = INT( CLDPARAMS%FR_LS_WAT )
         FR_LS_ICE     = INT( CLDPARAMS%FR_LS_ICE )
         FR_AN_WAT     = INT( CLDPARAMS%FR_AN_WAT )
         FR_AN_ICE     = INT( CLDPARAMS%FR_AN_ICE )
         MIN_RL        = CLDPARAMS%MIN_RL
         MIN_RI        = CLDPARAMS%MIN_RI
         MAX_RL        = CLDPARAMS%MAX_RL
         MAX_RI        = CLDPARAMS%MAX_RI
         FAC_RL        = CLDPARAMS%FAC_RL
         FAC_RI        = CLDPARAMS%FAC_RI
         PDFFLAG       = INT(CLDPARAMS%PDFSHAPE)

         call MAPL_TimerOff(STATE,"---CLOUD_DATA_CONST",RC=STATUS)
         VERIFY_(STATUS)

         call MAPL_TimerOff(STATE,"--CLOUD_DATA",RC=STATUS)
         VERIFY_(STATUS)

         call MAPL_TimerOn (STATE,"--CLOUD_RUN",RC=STATUS)
         VERIFY_(STATUS)

         call PROGNO_CLOUD<<<Grid, Block>>>(IM*JM,LM,DT_MOIST,CLOUD_CTL%SCLMFDFR)

         STATUS = cudaGetLastError()
         if (STATUS /= 0) then 
            write (*,*) "Error code from PROGNO_CLOUD kernel call: ", STATUS
            write (*,*) "Kernel call failed: ", cudaGetErrorString(STATUS)
            _ASSERT(.FALSE.,'needs informative message')
         end if

         call MAPL_TimerOff(STATE,"--CLOUD_RUN",RC=STATUS)
         VERIFY_(STATUS)

         call MAPL_TimerOn (STATE,"--CLOUD_DATA",RC=STATUS)
         VERIFY_(STATUS)

         call MAPL_TimerOn (STATE,"---CLOUD_DATA_DEVICE",RC=STATUS)
         VERIFY_(STATUS)

         ! Move device outputs to host
         ! ---------------------------

         ! Outputs
         ! -------

         STATUS = cudaMemcpy(RAD_CF,RAD_CLDFRC_DEV,IM*JM*LM)
         STATUS = cudaMemcpy(RAD_QL,RAD_QL_DEV,IM*JM*LM)
         STATUS = cudaMemcpy(RAD_QI,RAD_QI_DEV,IM*JM*LM)
         STATUS = cudaMemcpy(QRN,RAD_QR_DEV,IM*JM*LM)
         STATUS = cudaMemcpy(QSN,RAD_QS_DEV,IM*JM*LM)
         STATUS = cudaMemcpy(CLDREFFL,CLDREFFL_DEV,IM*JM*LM)
         STATUS = cudaMemcpy(CLDREFFI,CLDREFFI_DEV,IM*JM*LM)
         STATUS = cudaMemcpy(QPLS,QPLS_DEV,IM*JM*LM)
         STATUS = cudaMemcpy(LS_PRC2,PRELS_DEV,IM*JM)
         STATUS = cudaMemcpy(CN_PRC2,PRECU_DEV,IM*JM)
         STATUS = cudaMemcpy(AN_PRC2,PREAN_DEV,IM*JM)
         STATUS = cudaMemcpy(SC_PRC2,PRESC_DEV,IM*JM)
         STATUS = cudaMemcpy(LS_ARFX,LSARF_DEV,IM*JM)
         STATUS = cudaMemcpy(CN_ARFX,CUARF_DEV,IM*JM)
         STATUS = cudaMemcpy(AN_ARFX,ANARF_DEV,IM*JM)
         STATUS = cudaMemcpy(SC_ARFX,SCARF_DEV,IM*JM)
         STATUS = cudaMemcpy(LS_SNR,SNRLS_DEV,IM*JM)
         STATUS = cudaMemcpy(CN_SNR,SNRCU_DEV,IM*JM)
         STATUS = cudaMemcpy(AN_SNR,SNRAN_DEV,IM*JM)
         STATUS = cudaMemcpy(SC_SNR,SNRSC_DEV,IM*JM)

         ! Inoutputs
         ! ---------

         STATUS = cudaMemcpy(CNV_PRC3,QRN_CU_DEV,IM*JM*LM) 
         STATUS = cudaMemcpy(SHLW_PRC3,QRN_SC_DEV,IM*JM*LM) 
         STATUS = cudaMemcpy(SHLW_SNO3,QSN_SC_DEV,IM*JM*LM) 
         STATUS = cudaMemcpy(CNV_UPDF,CNV_UPDFRC_DEV,IM*JM*LM) 
         STATUS = cudaMemcpy(UPDF_SC,SC_UPDFRC_DEV,IM*JM*LM) 
         STATUS = cudaMemcpy(TH1,TH_DEV,IM*JM*LM) 
         STATUS = cudaMemcpy(Q1,Q_DEV,IM*JM*LM) 
         STATUS = cudaMemcpy(QLLS,QLW_LS_DEV,IM*JM*LM) 
         STATUS = cudaMemcpy(QLCN,QLW_AN_DEV,IM*JM*LM) 
         STATUS = cudaMemcpy(QILS,QIW_LS_DEV,IM*JM*LM) 
         STATUS = cudaMemcpy(QICN,QIW_AN_DEV,IM*JM*LM) 
         STATUS = cudaMemcpy(CLCN,ANVFRC_DEV,IM*JM*LM) 
         STATUS = cudaMemcpy(CLLS,CLDFRC_DEV,IM*JM*LM) 

         ! Diagnostics
         ! -----------

         IF (COPY_RHX)     STATUS = cudaMemcpy(RHX_X,RHX_DEV,SIZE(RHX_X))
         IF (COPY_REV_LS)  STATUS = cudaMemcpy(REV_LS_X,REV_LS_DEV,SIZE(REV_LS_X))
         IF (COPY_REV_AN)  STATUS = cudaMemcpy(REV_AN_X,REV_AN_DEV,SIZE(REV_AN_X))
         IF (COPY_REV_CN)  STATUS = cudaMemcpy(REV_CN_X,REV_CN_DEV,SIZE(REV_CN_X))
         IF (COPY_REV_SC)  STATUS = cudaMemcpy(REV_SC_X,REV_SC_DEV,SIZE(REV_SC_X))
         IF (COPY_RSU_LS)  STATUS = cudaMemcpy(RSU_LS_X,RSU_LS_DEV,SIZE(RSU_LS_X))
         IF (COPY_RSU_AN)  STATUS = cudaMemcpy(RSU_AN_X,RSU_AN_DEV,SIZE(RSU_AN_X))
         IF (COPY_RSU_CN)  STATUS = cudaMemcpy(RSU_CN_X,RSU_CN_DEV,SIZE(RSU_CN_X))
         IF (COPY_RSU_SC)  STATUS = cudaMemcpy(RSU_SC_X,RSU_SC_DEV,SIZE(RSU_SC_X))
         IF (COPY_ACLL_CN) STATUS = cudaMemcpy(ACLL_CN_X,ACLL_CN_DEV,SIZE(ACLL_CN_X))
         IF (COPY_ACIL_CN) STATUS = cudaMemcpy(ACIL_CN_X,ACIL_CN_DEV,SIZE(ACIL_CN_X))
         IF (COPY_ACLL_AN) STATUS = cudaMemcpy(ACLL_AN_X,ACLL_AN_DEV,SIZE(ACLL_AN_X))
         IF (COPY_ACIL_AN) STATUS = cudaMemcpy(ACIL_AN_X,ACIL_AN_DEV,SIZE(ACIL_AN_X))
         IF (COPY_ACLL_LS) STATUS = cudaMemcpy(ACLL_LS_X,ACLL_LS_DEV,SIZE(ACLL_LS_X))
         IF (COPY_ACIL_LS) STATUS = cudaMemcpy(ACIL_LS_X,ACIL_LS_DEV,SIZE(ACIL_LS_X))
         IF (COPY_ACLL_SC) STATUS = cudaMemcpy(ACLL_SC_X,ACLL_SC_DEV,SIZE(ACLL_SC_X))
         IF (COPY_ACIL_SC) STATUS = cudaMemcpy(ACIL_SC_X,ACIL_SC_DEV,SIZE(ACIL_SC_X))
         IF (COPY_PFL_CN)  STATUS = cudaMemcpy(PFL_CN_X,PFL_CN_DEV,SIZE(PFL_CN_X))
         IF (COPY_PFI_CN)  STATUS = cudaMemcpy(PFI_CN_X,PFI_CN_DEV,SIZE(PFI_CN_X))
         IF (COPY_PFL_AN)  STATUS = cudaMemcpy(PFL_AN_X,PFL_AN_DEV,SIZE(PFL_AN_X))
         IF (COPY_PFI_AN)  STATUS = cudaMemcpy(PFI_AN_X,PFI_AN_DEV,SIZE(PFI_AN_X))
         IF (COPY_PFL_LS)  STATUS = cudaMemcpy(PFL_LS_X,PFL_LS_DEV,SIZE(PFL_LS_X))
         IF (COPY_PFI_LS)  STATUS = cudaMemcpy(PFI_LS_X,PFI_LS_DEV,SIZE(PFI_LS_X))
         IF (COPY_PFL_SC)  STATUS = cudaMemcpy(PFL_SC_X,PFL_SC_DEV,SIZE(PFL_SC_X))
         IF (COPY_PFI_SC)  STATUS = cudaMemcpy(PFI_SC_X,PFI_SC_DEV,SIZE(PFI_SC_X))
         IF (COPY_DLPDF)   STATUS = cudaMemcpy(DLPDF_X,PDFL_DEV,SIZE(DLPDF_X))
         IF (COPY_DIPDF)   STATUS = cudaMemcpy(DIPDF_X,PDFI_DEV,SIZE(DIPDF_X))
         IF (COPY_DLFIX)   STATUS = cudaMemcpy(DLFIX_X,FIXL_DEV,SIZE(DLFIX_X))
         IF (COPY_DIFIX)   STATUS = cudaMemcpy(DIFIX_X,FIXI_DEV,SIZE(DIFIX_X))
         IF (COPY_AUT)     STATUS = cudaMemcpy(AUT_X,AUT_DEV,SIZE(AUT_X))
         IF (COPY_EVAPC)   STATUS = cudaMemcpy(EVAPC_X,EVAPC_DEV,SIZE(EVAPC_X))
         IF (COPY_SDM)     STATUS = cudaMemcpy(SDM_X,SDM_DEV,SIZE(SDM_X))
         IF (COPY_SUBLC)   STATUS = cudaMemcpy(SUBLC_X,SUBLC_DEV,SIZE(SUBLC_X))
         IF (COPY_FRZ_TT)  STATUS = cudaMemcpy(FRZ_TT_X,FRZ_TT_DEV,SIZE(FRZ_TT_X))
         IF (COPY_FRZ_PP)  STATUS = cudaMemcpy(FRZ_PP_X,FRZ_PP_DEV,SIZE(FRZ_PP_X))
         IF (COPY_DCNVL)   STATUS = cudaMemcpy(DCNVL_X,DCNVL_DEV,SIZE(DCNVL_X))
         IF (COPY_DCNVi)   STATUS = cudaMemcpy(DCNVi_X,DCNVI_DEV,SIZE(DCNVi_X))
         IF (COPY_ALPHT)   STATUS = cudaMemcpy(ALPHT_X,ALPHT_DEV,SIZE(ALPHT_X))
         IF (COPY_ALPH1)   STATUS = cudaMemcpy(ALPH1_X,ALPH1_DEV,SIZE(ALPH1_X))
         IF (COPY_ALPH2)   STATUS = cudaMemcpy(ALPH2_X,ALPH2_DEV,SIZE(ALPH2_X))
         IF (COPY_CFPDF)   STATUS = cudaMemcpy(CFPDF_X,CFPDF_DEV,SIZE(CFPDF_X))
         IF (COPY_RHCLR)   STATUS = cudaMemcpy(RHCLR_X,RHCLR_DEV,SIZE(RHCLR_X))
         IF (COPY_DQRL)    STATUS = cudaMemcpy(DQRL_X,DQRL_DEV,SIZE(DQRL_X))
         IF (COPY_VFALLICE_AN) STATUS = cudaMemcpy(VFALLICE_AN_X,VFALLICE_AN_DEV,SIZE(VFALLICE_AN_X))
         IF (COPY_VFALLICE_LS) STATUS = cudaMemcpy(VFALLICE_LS_X,VFALLICE_LS_DEV,SIZE(VFALLICE_LS_X))
         IF (COPY_VFALLWAT_AN) STATUS = cudaMemcpy(VFALLWAT_AN_X,VFALLWAT_AN_DEV,SIZE(VFALLWAT_AN_X))
         IF (COPY_VFALLWAT_LS) STATUS = cudaMemcpy(VFALLWAT_LS_X,VFALLWAT_LS_DEV,SIZE(VFALLWAT_LS_X))
         IF (COPY_VFALLSN_AN)  STATUS = cudaMemcpy(VFALLSN_AN_X,VFALLSN_AN_DEV,SIZE(VFALLSN_AN_X))
         IF (COPY_VFALLSN_LS)  STATUS = cudaMemcpy(VFALLSN_LS_X,VFALLSN_LS_DEV,SIZE(VFALLSN_LS_X))
         IF (COPY_VFALLSN_CN)  STATUS = cudaMemcpy(VFALLSN_CN_X,VFALLSN_CN_DEV,SIZE(VFALLSN_CN_X))
         IF (COPY_VFALLSN_SC)  STATUS = cudaMemcpy(VFALLSN_SC_X,VFALLSN_SC_DEV,SIZE(VFALLSN_SC_X))
         IF (COPY_VFALLRN_AN)  STATUS = cudaMemcpy(VFALLRN_AN_X,VFALLRN_AN_DEV,SIZE(VFALLRN_AN_X))
         IF (COPY_VFALLRN_LS)  STATUS = cudaMemcpy(VFALLRN_LS_X,VFALLRN_LS_DEV,SIZE(VFALLRN_LS_X))
         IF (COPY_VFALLRN_CN)  STATUS = cudaMemcpy(VFALLRN_CN_X,VFALLRN_CN_DEV,SIZE(VFALLRN_CN_X))
!!$         IF (COPY_LIQANMOVE)  STATUS = cudaMemcpy(LIQANMOVE_X,LIQANMOVE_DEV,SIZE(LIQANMOVE_X))
!!$         IF (COPY_ICEANMOVE)  STATUS = cudaMemcpy(ICEANMOVE_X,ICEANMOVE_DEV,SIZE(ICEANMOVE_X))
!!$         IF (COPY_DANCLD)  STATUS = cudaMemcpy(DANCLD_X,DANCLD_DEV,SIZE(DANCLD_X))
!!$         IF (COPY_DLSCLD)  STATUS = cudaMemcpy(DLSCLD_X,DLSCLD_DEV,SIZE(DLSCLD_X))
!!$         IF (COPY_CURAINMOVE)  STATUS = cudaMemcpy(CURAINMOVE_X,CURAINMOVE_DEV,SIZE(CURAINMOVE_X))
!!$         IF (COPY_CUSNOWMOVE)  STATUS = cudaMemcpy(CUSNOWMOVE_X,CUSNOWMOVE_DEV,SIZE(CUSNOWMOVE_X))

         call MAPL_TimerOff(STATE,"---CLOUD_DATA_DEVICE",RC=STATUS)
         VERIFY_(STATUS)

         call MAPL_TimerOff(STATE,"--CLOUD_DATA",RC=STATUS)
         VERIFY_(STATUS)

         call MAPL_TimerOn (STATE,"--CLOUD_DEALLOC",RC=STATUS)
         VERIFY_(STATUS)

         ! ------------------------
         ! Deallocate device arrays
         ! ------------------------

         ! Inputs
         ! ------

         DEALLOCATE(PP_DEV)
         DEALLOCATE(EXNP_DEV)
         DEALLOCATE(PPE_DEV)
         DEALLOCATE(KH_DEV)
         DEALLOCATE(SNOMAS_DEV)
         DEALLOCATE(FRLANDICE_DEV)
         DEALLOCATE(FRLAND_DEV)
         DEALLOCATE(RMFDTR_DEV)
         DEALLOCATE(SMFDTR_DEV)
         DEALLOCATE(QLWDTR_DEV)
         DEALLOCATE(SQLWDTR_DEV)
         DEALLOCATE(SQIWDTR_DEV)
         DEALLOCATE(U_DEV)
         DEALLOCATE(V_DEV)
         DEALLOCATE(QST3_DEV)
         DEALLOCATE(DZET_DEV)
         DEALLOCATE(QDDF3_DEV)
         DEALLOCATE(TEMPOR_DEV)
         DEALLOCATE(CNV_FRACTION_DEV) 
         DEALLOCATE(TROPP_DEV)

         ! Inoutputs
         ! ---------

         DEALLOCATE(TH_DEV)
         DEALLOCATE(Q_DEV)
         DEALLOCATE(QRN_CU_DEV)
         DEALLOCATE(CNV_UPDFRC_DEV)
         DEALLOCATE(SC_UPDFRC_DEV)
         DEALLOCATE(QLW_LS_DEV)
         DEALLOCATE(QLW_AN_DEV)
         DEALLOCATE(QIW_LS_DEV)
         DEALLOCATE(QIW_AN_DEV)
         DEALLOCATE(ANVFRC_DEV)
         DEALLOCATE(CLDFRC_DEV)

         ! Outputs
         ! -------

         DEALLOCATE(RAD_CLDFRC_DEV)
         DEALLOCATE(RAD_QL_DEV)
         DEALLOCATE(RAD_QI_DEV)
         DEALLOCATE(RAD_QR_DEV)
         DEALLOCATE(RAD_QS_DEV)
         DEALLOCATE(CLDREFFL_DEV)
         DEALLOCATE(CLDREFFI_DEV)
         DEALLOCATE(QPLS_DEV)
         DEALLOCATE(PRELS_DEV)
         DEALLOCATE(PRECU_DEV)
         DEALLOCATE(PREAN_DEV)
         DEALLOCATE(PRESC_DEV)
         DEALLOCATE(LSARF_DEV)
         DEALLOCATE(CUARF_DEV)
         DEALLOCATE(ANARF_DEV)
         DEALLOCATE(SNRLS_DEV)
         DEALLOCATE(SNRCU_DEV)
         DEALLOCATE(SNRSC_DEV)
         DEALLOCATE(SNRAN_DEV)

         ! Because of use_autoconv_timescale in PROGNO_CLOUD, the PFL and PFI
         ! arrays can "silently" be used as working arrays
         ! They must always be allocated
         DEALLOCATE(PFL_CN_DEV)
         DEALLOCATE(PFI_CN_DEV)
         DEALLOCATE(PFL_AN_DEV)
         DEALLOCATE(PFI_AN_DEV)
         DEALLOCATE(PFL_LS_DEV)
         DEALLOCATE(PFI_LS_DEV)
         DEALLOCATE(PFL_SC_DEV)
         DEALLOCATE(PFI_SC_DEV)

         ! Diagnostics
         ! -----------


         DEALLOCATE(RHX_DEV)
         DEALLOCATE(REV_LS_DEV)
         DEALLOCATE(REV_AN_DEV)
         DEALLOCATE(REV_CN_DEV)
         DEALLOCATE(REV_SC_DEV)
         DEALLOCATE(RSU_LS_DEV)
         DEALLOCATE(RSU_AN_DEV)
         DEALLOCATE(RSU_CN_DEV)
         DEALLOCATE(RSU_SC_DEV)
         DEALLOCATE(ACLL_CN_DEV)
         DEALLOCATE(ACIL_CN_DEV)
         DEALLOCATE(ACLL_AN_DEV)
         DEALLOCATE(ACIL_AN_DEV)
         DEALLOCATE(ACLL_LS_DEV)
         DEALLOCATE(ACIL_LS_DEV)
         DEALLOCATE(ACLL_SC_DEV)
         DEALLOCATE(ACIL_SC_DEV)
         DEALLOCATE(PDFL_DEV)
         DEALLOCATE(PDFI_DEV)
         DEALLOCATE(FIXL_DEV)
         DEALLOCATE(FIXI_DEV)
         DEALLOCATE(AUT_DEV)
         DEALLOCATE(EVAPC_DEV)
         DEALLOCATE(SDM_DEV)
         DEALLOCATE(SUBLC_DEV)
         DEALLOCATE(FRZ_TT_DEV)
         DEALLOCATE(FRZ_PP_DEV)
         DEALLOCATE(DCNVL_DEV)
         DEALLOCATE(DCNVI_DEV)
         DEALLOCATE(ALPHT_DEV)
         DEALLOCATE(ALPH1_DEV)
         DEALLOCATE(ALPH2_DEV)
         DEALLOCATE(CFPDF_DEV)
         DEALLOCATE(RHCLR_DEV)
         DEALLOCATE(DQRL_DEV)
         DEALLOCATE(VFALLICE_AN_DEV)
         DEALLOCATE(VFALLICE_LS_DEV)
         DEALLOCATE(VFALLWAT_AN_DEV)
         DEALLOCATE(VFALLWAT_LS_DEV)
         DEALLOCATE(VFALLSN_AN_DEV)
         DEALLOCATE(VFALLSN_LS_DEV)
         DEALLOCATE(VFALLSN_CN_DEV)
         DEALLOCATE(VFALLSN_SC_DEV)
         DEALLOCATE(VFALLRN_AN_DEV)
         DEALLOCATE(VFALLRN_LS_DEV)
         DEALLOCATE(VFALLRN_CN_DEV)
!!$         DEALLOCATE(LIQANMOVE_DEV)
!!$         DEALLOCATE(ICEANMOVE_DEV)
!!$         DEALLOCATE(DANCLD_DEV)
!!$         DEALLOCATE(DLSCLD_DEV)
!!$         DEALLOCATE(CURAINMOVE_DEV)
!!$         DEALLOCATE(CUSNOWMOVE_DEV)

         STATUS = cudaDeviceSynchronize()

         call MAPL_TimerOff(STATE,"--CLOUD_DEALLOC",RC=STATUS)
         VERIFY_(STATUS)

#else 
 !CUDA DEF

         call MAPL_TimerOn (STATE,"--CLOUD_RUN",RC=STATUS)
         VERIFY_(STATUS)

         SC_SNR = 0.0
         SC_PRC2 = 0.0
         SC_ARFX = 0.0
         SC_SNR  = 0.0
         REV_SC_X = 0.0
         RSU_SC_X = 0.0
         PFL_SC_X = 0.0
         PFI_SC_X = 0.0

!         print *,"QLDET_SC=",QLDET_SC
!         print *,"CNV_DQLDT=",CNV_DQLDT

         call  PROGNO_CLOUD (                    &
              IM*JM, LM         , &
              DT_MOIST          , &
              LATS              , &
              PLO               , &
              ZLO               , &
              CNV_PLE           , &
              PK                , &
              SNOMAS            , &   ! <- surf
              FRLANDICE         , &   ! <- surf
              FRLAND            , &   ! <- surf
              KH                , &   ! <- turb
              EDMF_FRC          , &   ! <- turb
              WQT               , &   ! <- turb
              WHL               , &   ! <- turb
              QT2               , &   ! <- turb
              HL2               , &   ! <- turb
              HLQT              , &   ! <- turb
              W2                , &   ! <- turb
              W3                , &   ! <- turb
              QT3               , &
              HL3               , &
              DTS               , &
              CNV_MFD           , &   ! <- ras
              CNV_DQLDT         , &   ! <- ras              
              CNV_PRC3          , &   ! <- ras   
              CNV_UPDF          , &   ! <- ras
              MFD_SC            , &   ! <- shlw   
              QLDET_SC          , &   ! <- shlw   
              QIDET_SC          , &   ! <- shlw   
              SHLW_PRC3         , &   ! <- shlw   
              SHLW_SNO3         , &   ! <- shlw   
              UFRC_SC           , &   ! <- shlw   
              U1                , &
              V1                , & 
              TH1               , &              
              Q1                , &
              QLLS              , &
              QLCN              , &
              QILS              , &
              QICN              , &
              CLCN              , &
              CLLS              , &
              RAD_CF            , &
              RAD_QV            , &
              RAD_QL            , &
              RAD_QI            , &
              QRN               , &
              QSN               , &
              QPLS              , &
              CLDREFFL          , &
              CLDREFFI          , &
              LS_PRC2           , &
              CN_PRC2           , &
              AN_PRC2           , &
              SC_PRC2           , &
              LS_ARFX           , &
              CN_ARFX           , &
              AN_ARFX           , &
              SC_ARFX           , &
              LS_SNR            , &
              CN_SNR            , &
              AN_SNR            , &
              SC_SNR            , &
              CLDPARAMS         , &
              CLOUD_CTL%SCLMFDFR, &
              QST3              , &
              DZET              , &
              QDDF3             , &
              CNV_FRACTION      , &
              TROPP             , &
                                ! Diagnostics
              RHX_X             , &
              REV_LS_X          , &
              REV_AN_X          , &
              REV_CN_X          , &
              REV_SC_X          , &
              RSU_LS_X          , &
              RSU_AN_X          , &
              RSU_CN_X          , &
              RSU_SC_X          , &
              ACLL_CN_X,ACIL_CN_X   , &
              ACLL_AN_X,ACIL_AN_X   , &
              ACLL_LS_X,ACIL_LS_X   , &
              ACLL_SC_X,ACIL_SC_X   , &
              PFL_CN_X,PFI_CN_X     , &
              PFL_AN_X,PFI_AN_X     , &
              PFL_LS_X,PFI_LS_X     , &
              PFL_SC_X,PFI_SC_X     , &
              DLPDF_X,DIPDF_X,DLFIX_X,DIFIX_X,    &
              AUT_X, EVAPC_X , SDM_X , SUBLC_X ,  &
              FRZ_TT_X, DCNVL_X, DCNVi_X,       &
              ALPHT_X, ALPH1_X, ALPH2_X, CFPDF_X, &
              RHCLR_X,                      &
              DQRL_X, FRZ_PP_X,               &
              VFALLICE_AN_X,VFALLICE_LS_X,    &
              VFALLWAT_AN_X,VFALLWAT_LS_X,    &
              VFALLSN_AN_X,VFALLSN_LS_X,VFALLSN_CN_X,VFALLSN_SC_X,  &
              VFALLRN_AN_X,VFALLRN_LS_X,VFALLRN_CN_X,VFALLRN_SC_X,  &
              PDF_A, &
#ifdef PDFDIAG
              PDF_SIGW1, PDF_SIGW2, PDF_W1, PDF_W2, & 
              PDF_SIGTH1, PDF_SIGTH2, PDF_TH1, PDF_TH2, &
              PDF_SIGQT1, PDF_SIGQT2, PDF_QT1, PDF_QT2, &
              PDF_RQTTH, PDF_RWTH, PDF_RWQT,            &
#endif
              WTHV2, WQL, &
              TEMPOR2D, &
              DOSHLW,   &
              NACTL,    &
              NACTI,    &
              CONVPAR_OPTION )

         VERIFY_(STATUS)

         if (associated(PDF_AX)) PDF_AX = PDF_A

         call MAPL_TimerOff(STATE,"--CLOUD_RUN",RC=STATUS)
         VERIFY_(STATUS)

         if (USE_AEROSOL_NN) then
           CFX =100.*PLO*r_air/TEMP !density times conversion factor
           NCPL = NACTL/CFX ! kg-1
           NCPI = NACTI/CFX ! kg-1
         else
           NCPL = 0.
           NCPI = 0.
         endif
#endif
 !CUDA DEF


         !-----------------------------------------
         !! kluge in some skewness for PBL clouds 
         !! where DTS:=TS-TSFCAIR > 0. i.e., unstable
         RAD_QV   = max( Q1 , 0. )

         if (CLDPARAMS%CFPBL_EXP > 1.0) then
          CFPBL = RAD_CF
          do L = 1,LM
             where( (KH(:,:,L-1) > 5.) .AND. (DTS > 0.) )
                CFPBL(:,:,L) = RAD_CF(:,:,L)**CLDPARAMS%CFPBL_EXP
             endwhere
          enddo
          where( CFPBL > 0.0 )
             RAD_QL = RAD_QL*(RAD_CF/CFPBL)
             RAD_QI = RAD_QI*(RAD_CF/CFPBL)
             QRN    = QRN   *(RAD_CF/CFPBL)
             QSN    = QSN   *(RAD_CF/CFPBL)
          endwhere
          RAD_CF = CFPBL
         endif

         RAD_QL = MIN( RAD_QL , 0.001 )  ! Still a ridiculously large
         RAD_QI = MIN( RAD_QI , 0.001 )  ! value.
         QRN    = MIN( QRN   , 0.01 )  ! value.
         QSN    = MIN( QSN   , 0.01 )  ! value.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         ! Set rain water for radiation to 0 if preciprad flag is off (set to 0)
         if(CLDPARAMS%PRECIPRAD.eq.0.) then
            RAD_QR = 0.
            RAD_QS = 0.
            RAD_QG = 0.
         else
            RAD_QR = QRN
            RAD_QS = QSN
         endif

         if (associated(QRTOT)) QRTOT = QRN*RAD_CF
         if (associated(QSTOT)) QSTOT = QSN*RAD_CF

         if (associated(QPTOTLS)) QPTOTLS = QPLS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         CLDREFFR = 100.e-6
         CLDREFFS = 140.e-6

         !Calculate CFICE and CFLIQ 
         CFLIQ=0.0
         CFICE=0.0
         QTOT= QICN+QILS+QLCN+QLLS
         QL_TOT = QLCN+QLLS
         QI_TOT = QICN+QILS
         WHERE (QTOT .gt. 1.0e-12)
            CFLIQ=RAD_CF*QL_TOT/QTOT
            CFICE=RAD_CF*QI_TOT/QTOT
         END WHERE
         CFLIQ=MAX(MIN(CFLIQ, 1.0), 0.0)
         CFICE=MAX(MIN(CFICE, 1.0), 0.0)

         !======================================================================================================================
         !===========================Clean stuff and send it to radiation ======================================================
         !======================================================================================================================

         where (QI_TOT .le. 0.0)
            CFICE =0.0
            NCPI=0.0
            CLDREFFI = 36.0e-6
         end where

         where (QL_TOT .le. 0.0)
            CFLIQ =0.0
            NCPL  =0.0
            CLDREFFL = 14.0e-6
         end where

         if (associated(DQVDT_micro)) DQVDT_micro = (Q1 - DQVDT_micro         ) / DT_MOIST
         if (associated(DQLDT_micro)) DQLDT_micro = ((QLLS+QLCN) - DQLDT_micro) / DT_MOIST
         if (associated(DQIDT_micro)) DQIDT_micro = ((QILS+QICN) - DQIDT_micro) / DT_MOIST
         if (associated(DQRDT_micro)) DQRDT_micro = 0.0
         if (associated(DQSDT_micro)) DQSDT_micro = 0.0
         if (associated(DQGDT_micro)) DQGDT_micro = 0.0
         if (associated(DQADT_micro)) DQADT_micro = ((CLLS+CLCN) - DQADT_micro) / DT_MOIST
         if (associated(DUDT_micro) ) DUDT_micro  = (U1 - DUDT_micro - U1) / DT_MOIST
         if (associated(DVDT_micro) ) DVDT_micro  = (V1 - DVDT_micro - V1) / DT_MOIST
         if (associated(DTDT_micro) ) DTDT_micro  = (TH1*PK - DTDT_micro) / DT_MOIST
         if (associated(DTDT_macro) ) DTDT_macro  = 0.0
         if (associated(DQVDT_macro) ) DQVDT_macro  = 0.0
         if (associated(DQLDT_macro) ) DQLDT_macro  = 0.0
         if (associated(DQIDT_macro) ) DQIDT_macro  = 0.0
         if (associated(DQADT_macro) ) DQADT_macro  = 0.0

       endif !====1-moment Microphysics=

         !==================================
      else !====2-moment Microphysics=
         !==================================

         !============================================= Start Stratiform cloud processes==========================================
         !set up initial values

         call MAPL_TimerOn (STATE,"--CLOUD_RUN",RC=STATUS)
         VERIFY_(STATUS)

         TEMP    = TH1*PK   
         DTDT_macro=TEMP
         DQVDT_macro=Q1
         DQLDT_macro=QLCN+QLLS
         DQIDT_macro=QICN+QILS
         DQADT_macro=CLCN+CLLS
 
         SC_ICE=1.0
         NCPL=MAX( NCPL , 0. )
         NCPI=MAX( NCPI , 0. )
         CLDREFFR = 10.0e-6 
         CLDREFFS = 90.0e-6
         CLDREFFG = 90.0e-6
         CLDREFFI = 25.0e-6
         CLDREFFL = 10.0e-6
         RAD_CF   = min(CLLS+CLCN, 1.0)
         RAD_QL   = 0.0
         RAD_QI   = 0.0
         RAD_QR   = 0.0
         RAD_QS   = 0.0
         RAD_QV   = Q1
         CDNC_NUC = 0.0
         INC_NUC  = 0.0
         T_ICE_MAX = MAPL_TICE     ! -7.0+TICE DONIFF
         T_ICE_ALL = CLDPARAMS%ICE_RAMP + MAPL_TICE 
	 QSNOW_CN = 0.0
	 QRAIN_CN = 0.0
	 PFRZ= 0.0

      ! Find estimated inversion strength (DONIF)


       call    FIND_EIS(TH1, QSS, TEMP, ZLE, PLO, KLCL, IM, JM, LM, LTS, EIS)


         ! Clean up any negative specific humidity before the microphysics scheme
         !-----------------------------------------
         call FILLQ2ZERO2( Q1, MASS, FILLQ)  


          PFL_AN_X   = 0.0 
	  PFL_LS_X  = 0.0


         !=======================================================================================================================
         !=======================================================================================================================
         !===================================Nucleation of cloud droplets and ice crystals ======================================
         ! Aerosol cloud interactions. Calculate maxCCN tendency using Fountoukis and nenes (2005) or Abdul Razzak and Ghan (2002)
         ! liquid Activation Parameterization
         ! Ice activation follows the Barahona & Nenes ice activation scheme, ACP, (2008, 2009). 
         ! Written by Donifan Barahona and described in Barahona et al. (2013)
         !=======================================================================================================================
         !=======================================================================================================================
         !=======================================================================================================================

         call MAPL_TimerOn(STATE,"---ACTIV") !Activation timer


         use_average_v = .false.  
         if (USE_AV_V .gt. 0.0) then   
            use_average_v = .true.
         end if
        fdust_drop =  FDROP_DUST
	fsoot_drop =  FDROP_SOOT
	sigma_nuc_r8  = SIGMA_NUC
     frachet_org = 	ORG_INFAC
     frachet_dust = 	DUST_INFAC
     frachet_bc = 	BC_INFAC
     frachet_ss = 	SS_INFAC

	CFX=0.0
	  where (QSS > 0.0) 
	    CFX =Q1/(QSS)
	  end where 
 !recalculate bkgtau: scaling of W variance with respect to Nature run
	 if (USE_NATURE_WSUB .gt. 0.) then 
            xscale = (72000.0/imsize)            
            BKGTAU=  1.472/sqrt(1.0+ (xscale/6.0)) 
            BKGTAU = max((1.71 - BKGTAU), 0.0)*SWCIRRUS
             
      end if 
	
	 do J=1,JM
            do I=1,IM

               	     ccn_diag  = 0.0
                     smaxliq   = 0.0
                     smaxicer8 = 0.0
                     nheticer8 = 0.0
                     sc_icer8  = 1.0
                     naair8    = 0.0
                     npccninr8 = 0.0
                     nlimicer8 = 0.0    
                     nhet_immr8 = 0.0
                     dnhet_immr8 = 0.0
                     nhet_depr8 = 0.0
                     dust_immr8 = 0.0
                     dust_depr8 = 0.0
                      so4x = 0.0
   	                  dustx = 0.0
                      bcx= 0.0
	                  orgx=0.0
	                  seasaltx=0.0
                     
                     
                     
                     uwind_gw(1,1:LM)           = min(0.5*SQRT( U1(I,J,1:LM)**2+  V1(I,J,1:LM)**2), 50.0)
                     tausurf_gw   = min(0.5*SQRT(TAUOROX(I , J)**2+TAUOROY(I , J)**2), 10.0) !limit to a very high value     
                     if (USE_NATURE_WSUB .le. 0.) then 
                     tausurf_gw   =tausurf_gw  + min(0.5*SQRT(TAUX(I , J)**2+TAUY(I , J)**2), 5.0)*BKGTAU !adds a minimum value from unresolved sources (rewritten 04/01/15)
                     end if 
                     
                     pi_gw(1, 0:LM)        = 100.0*CNV_PLE(I,J,0:LM)                     
                     theta_tr(1, 1:LM)     = TH1(I,J,1:LM)
		             rhoi_gw = 0.0  
                     ni_gw   = 0.0 
                     ti_gw   = 0.0                                         
                     ter8(1,1:LM)      = TEMP(I,J,1:LM)   
		             pi_gw(1, 0:LM)   = 100.0*CNV_PLE(I,J,0:LM) 
                     plevr8(1,1:LM)    = 100.*PLO(I,J,:)
                     ndropr8(1, 1:LM) = NCPL(I, J, 1:LM)
		            qir8(1, 1:LM)     =  QILS(I, J,1:LM)+QICN(I, J,1:LM)
                    qcr8(1, 1:LM)     =  QLLS(I, J,1:LM)+QLCN(I, J,1:LM)
                    
                  ! where (RAD_CF(I, J, 1:LM) .gt. 0.01)
                  !  npre8(1,1:LM)     = NPRE_FRAC_2d(I,J)*NCPI(I,J,1:LM)/RAD_CF(I, J, 1:LM)
                  ! elsewhere 
                    npre8(1,1:LM)     = NPRE_FRAC_2d(I,J)*NCPI(I,J,1:LM)
                  ! end where 
                    
                    
                    omegr8(1,1:LM)    = OMEGA(I,J,1:LM)  
                    
                    rad_cooling(1,1:LM) = RADLW(I,J,1:LM)+RADSW(I,J,1:LM)
		    wparc_ls = 0.0
		    wparc_gw = 0.0
		    wparc_cgw= 0.0
		    wparc_turb = 0.0
		    swparc=0.0
		    tm_gw =ter8
		    pm_gw =plevr8
		    pfrz_inc_r8 = 0.0
		    Ksa1= 1.0
                    
		    if (FRLAND(I, J) .lt. 0.1) then 
 		        lc_turb(1,1:LM)   =  max(ALH(I,J,1:LM), MIN_ALH) 
                     else
                        lc_turb(1,1:LM)   =  max(ALH(I,J,1:LM), 50.0)
	            end if 
			 
                    where ((npre8 .gt. 0.0)   .and. (qir8 .gt. 0.0))
                         dpre8    = ( qir8/(5400.0*npre8*MAPL_PI))**(0.33) !Assume exponential distribution
                    elsewhere
                        dpre8=1.0e-9
                   end where

                   call   gw_prof (1, LM, 1, tm_gw, pm_gw, pi_gw, &
                                 rhoi_gw, ni_gw, ti_gw, nm_gw) !get Brunt_Vaisala Frequency and midpoint densities 

                    kcldtopcvn=KCT(I, J)
		            Nct =nm_gw(1, kcldtopcvn)      !BV frequency ar cloud top
                    Wct = max(CNV_CVW(I, J, kcldtopcvn), 0.0)
                    fcn = maxval(CNV_UPDF(I, J, kcldtopcvn:LM))    

		            kbmin= min(KCBL(I, J), LM -1)-2
                    maxkhpbl=maxval(KH(I, J, kbmin:LM-1))    
		    
		    !ksa1= max(real(kbmin - kcldtopcvn), 1.0)

               ! ==========================================================================================    
               ! ========================Activate the aerosols ============================================ 
	   
	       
	       
                do K = KMIN_TROP(I, J), LM-1 !limit to troposphere and no activation at the surface
    
                ! find vertical velocity variance 
!		       call zeit_ci("MOIST::aero_vvar")


		       call vertical_vel_variance(omegr8(1, K), lc_turb(1, K), ter8(1, K), plevr8(1, K), rad_cooling(1,K),  uwind_gw(1,K), &
		                                         tausurf_gw, nm_gw(1, K), LCCIRRUS, Nct, Wct, &
                                                         ksa1, fcn(1, K), KH(I, J, K), FRLAND(I, J), ZPBL(I, J), ZLE(I, J, k), maxkhpbl, &
						            wparc_ls(1, K), wparc_gw(1, K), wparc_cgw(1, K), wparc_turb(1, K), EIS(I, J), TKE(I, J, K))
		     
!			       call zeit_co("MOIST::aero_vvar")
					
                   if (FRLAND(I, J) .lt. 0.1) then 
                    if (LTS(I, J) .gt. CLDPARAMS%MIN_LTS) then                     
	                   if (K .ge. kbmin-2) wparc_ls(1, K)=max(wparc_ls(1,K)+ zws(i, j), 0.00)*SCWST ! add convective velocity within the PBL
                    end if 
                   else
                      if (K .ge. kbmin-2) wparc_ls(1, K)=max(wparc_ls(1,K)+ zws(i, j), 0.00) 
                   end if 

                     if (K .ge. kbmin-2) wparc_turb(1, K)=max(wparc_turb(1,K), 0.04)    !minimum velocity within the PBL (not resolved by RAS)
                     			          
                     if (K .ge.  kcldtopcvn) wparc_cgw(1, K) = 0.0                     
                    

			
			 if (USE_NATURE_WSUB .gt. 0.) then !use climatology from the Nature run (only for cirrus)
			 	
				  !wparc_cgw(1, k)= max(WSUB_NATURE(I, J, K)+BKGTAU*BKGTAU, 0.0)!BKG accounts for unresolved vertical velocity at 7 km				  
				   wparc_cgw(1, k)= max(WSUB_NATURE(I, J, K)*BKGTAU*BKGTAU, 0.0)!BKG accounts for unresolved vertical velocity at 7 km				  
				   wparc_gw(1, k) = 0.0

                        end if 

                             swparc(1, K)=sqrt(wparc_gw(1, K)+wparc_turb(1, K)+ wparc_cgw(1, K))
 
		        
				
                       !Supersaturations to calculate CCN diagnostics
                        ccn_diag(1)=0.001
                        ccn_diag(2)=0.004
                        ccn_diag(3)=0.01

                         AeroAux%nmods = 0
                         AeroAux%num   = 0.0
                         do i_src_mode = 1, AeroProps(I,J,K)%nmods
                            if (AeroProps(I,J,K)%num(i_src_mode) > 0.1) then
                               AeroAux%nmods = AeroAux%nmods + 1
                               i_dst_mode = AeroAux%nmods

                               AeroAux%num(i_dst_mode)   = AeroProps(I,J,K)%num(i_src_mode)
                               AeroAux%dpg(i_dst_mode)   = AeroProps(I,J,K)%dpg(i_src_mode)
                               AeroAux%sig(i_dst_mode)   = AeroProps(I,J,K)%sig(i_src_mode)
                               AeroAux%den(i_dst_mode)   = AeroProps(I,J,K)%den(i_src_mode)
                               AeroAux%kap(i_dst_mode)   = AeroProps(I,J,K)%kap(i_src_mode)
                               AeroAux%fdust(i_dst_mode) = AeroProps(I,J,K)%fdust(i_src_mode)
                               AeroAux%fsoot(i_dst_mode) = AeroProps(I,J,K)%fsoot(i_src_mode)
                               AeroAux%forg(i_dst_mode)  = AeroProps(I,J,K)%forg(i_src_mode)
                             end if
                        end do

		       rh1_r8=CFX(I, J, K)
		       
		       !if ((K.gt. 2) .and. (K .lt. LM-1)) then 
		       !    tauxr8 =  (ter8(1, K+1) +  ter8(1, K) + ter8(1, K-1))/3.0
			   !else 
			     tauxr8 = ter8(1, K)
             ! end if 
                     	
!                        call zeit_co("MOIST::aero_unpack")
                                   
		     !!Subroutine aerosol_activate contains the CCN activation and ice nucleation parameterizations. Lives in aer_cloud.F90.

                     call   aerosol_activate(tauxr8, plevr8(1, K), swparc(1, K), wparc_ls(1, K),  AeroAux, &
                          npre8(1, k), dpre8(1, k), ccn_diag, ndropr8(1, k), qcr8(1, K), &
                          npccninr8(1, K), smaxliq(1, K), naair8(1, K), smaxicer8(1, K), nheticer8(1, K), &
                          nhet_immr8(1, K), dnhet_immr8(1, K), nhet_depr8(1, k), sc_icer8(1, k), &
                          dust_immr8(1, K), dust_depr8(1, k), nlimicer8(1, k), use_average_v, int(CCN_PARAM), int(IN_PARAM),  &
                          so4x(1, k), seasaltx(1, k), dustx(1, k), orgx(1, K), bcx(1, k), &                          
			              fdust_drop, fsoot_drop, pfrz_inc_r8(1, K), rh1_r8, frachet_dust, frachet_bc, frachet_org, frachet_ss, int(Immersion_PARAM))

                      CCN01(I, J, K) = max(ccn_diag(1), 0.0)
                      CCN04(I, J, K) = max(ccn_diag(2), 0.0)
                      CCN1 (I, J, K) = max(ccn_diag(3), 0.0)
                      if (K .ge. kbmin) npccninr8(1, K) =max(npccninr8(1, K), MINCDNC*1.e6)
                       
               end do

               SMAXL(I, J, 1:LM) = smaxliq(1, 1:LM)*100.0         
               SMAXI(I, J, 1:LM) = smaxicer8(1, 1:LM)*100.0
               NHET_NUC(I, J, 1:LM)  = nheticer8(1, 1:LM)
               NLIM_NUC(I, J, 1:LM) =  nlimicer8(1, 1:LM)            
               SC_ICE(I, J, 1:LM) = sc_icer8(1, 1:LM)                  
               CDNC_NUC(I,J,1:LM)    = npccninr8(1, 1:LM)
               INC_NUC (I,J,1:LM)    = naair8(1, 1:LM)         
               NHET_IMM(I, J, 1:LM)  = max(nhet_immr8(1, 1:LM), 0.0)
               DNHET_IMM(I, J, 1:LM)  = max(dnhet_immr8(1, 1:LM), 0.0)
               NHET_DEP(I, J, 1:LM)  = nhet_depr8(1, 1:LM)
               DUST_IMM(I, J, 1:LM)  = max(dust_immr8(1, 1:LM), 0.0)
               DUST_DEP(I, J, 1:LM)  = max(dust_depr8(1, 1:LM), 0.0)
               WSUB (I, J, 1:LM) =  wparc_ls(1, 1:LM)+swparc(1, 1:LM)*0.8        
               SIGW_GW (I, J, 1:LM)   =  wparc_gw(1, 1:LM)
               SIGW_CNV (I, J, 1:LM)   =  wparc_cgw(1, 1:LM)
               SIGW_TURB (I, J, 1:LM) = wparc_turb(1, 1:LM)
               SIGW_RC (I, J, 1:LM)   =  wparc_ls(1, 1:LM)
	           PFRZ (I, J, 1:LM)   =  pfrz_inc_r8(1, 1:LM)
                
               SO4(I, J, 1:LM)=so4x(1, 1:LM)
   	           DUST(I, J, 1:LM)=dustx(1, 1:LM)
               BCARBON(I, J, 1:LM)=bcx(1, 1:LM)
	           ORG(I, J, 1:LM)=orgx(1, 1:LM)
	           SEASALT(I, J, 1:LM)=seasaltx(1, 1:LM)


	       
            enddo
         enddo

    
	
         !=============================================End cloud particle nucleation=====================================
         !===============================================================================================================

         call MAPL_TimerOff(STATE,"---ACTIV", RC=STATUS)
         call MAPL_TimerOn(STATE,"---CLDMACRO")

         !==========================================================================================================
         !===================================Cloud Macrophysics ====================================================
         !==========================================================================================================

         REV_CN_X  = 0.0
         REV_AN_X  = 0.0 
         REV_LS_X  = 0.0
         REV_SC_X  = 0.0
         RSU_CN_X  = 0.0
         RSU_AN_X  = 0.0
         RSU_LS_X  = 0.0     
         RSU_SC_X  = 0.0     

         CFX=INC_NUC + NHET_IMM 
      
         CN_PRC2 = 0.0 
         LS_PRC2= 0.0
         AN_PRC2 = 0.0 
         SC_PRC2 = 0.0 
         CN_SNR  =0.0 
         LS_SNR = 0.0  
         AN_SNR = 0.0
         SC_SNR = 0.0
         DTDT_macro=TEMP     
         DQVDT_macro=Q1
         DQLDT_macro=QLCN+QLLS
         DQIDT_macro=QICN+QILS
         DQADT_macro=CLCN+CLLS
         PFI_CN_X= 0.0
         PFI_AN_X=0.0
         PFI_LS_X=0.0
         PFI_SC_X=0.0
         PFL_CN_X= 0.0
         PFL_AN_X=0.0
         PFL_LS_X=0.0    
         PFL_SC_X=0.0    
      

    ! if(associated(TVQX1))  TVQX1     = SUM( (  Q1 +  QLLS + QLCN + QILS + QICN + CNV_PRC3)*DM + CNV_DQLDT*DT_MOIST , 3 )
    
         
              CNV_MFD_X    =  CNV_MFD     ! needed for cloud fraction
              CNV_DQLDT_X    =  CNV_DQLDT
              CNV_PRC3_X     = CNV_PRC3
              CNV_UPDF_X     = CNV_UPDF 
              CNV_NICE_X =    CNV_NICE_X +  SC_NICE
              CNV_NDROP_X =   CNV_NDROP_X +  SC_NDROP
  
    IF(ADJUSTL(CONVPAR_OPTION) == 'GF') THEN    
              CNV_PRC3_X     = 0.0
     END IF 
        
        
        if (MGVERSION .gt. 1.0) then 
        
        
              if(associated(TVQX1))  TVQX1     =  SUM( (  Q1 +  QLLS + QLCN + QILS + QICN + QRAIN + QSNOW + QGRAUPEL + SHLW_PRC3 + SHLW_SNO3)*MASS &
        
      
                 + (CNV_DQLDT)*DT_MOIST &
                 
                 + (QLDET_SC  + QIDET_SC)*DT_MOIST &
                  
                 , 3 ) + RASPRCP*DT_MOIST - TVQ0 ! up to here water is conserved Donif 01/2020
                 
        else
        
          if(associated(TVQX1))  TVQX1     =  SUM( (  Q1 +  QLLS + QLCN + QILS + QICN +  SHLW_PRC3 + SHLW_SNO3)*MASS &
        
      
                 + (CNV_DQLDT)*DT_MOIST &
                 
                 + (QLDET_SC  + QIDET_SC)*DT_MOIST &
                  
                 , 3 ) + RASPRCP*DT_MOIST - TVQ0 ! up to here water is conserved Donif 01/2020
 
         end if     
 
                       
  if (DOCLDMACRO/=0) then   
  call  macro_cloud (                    &
              IM*JM, LM         , &
              DT_MOIST          , &
              PLO               , &
              CNV_PLE           , &
              PK                , &
              SNOMAS            , &   ! <- surf
              FRLANDICE         , &   ! <- surf
              FRLAND            , &   ! <- surf
              KH                , &   ! <- turb
              CNV_MFD_X           , &   ! <- ras
              CNV_DQLDT_X         , &   ! <- ras              
              CNV_PRC3_X          , &   ! <- ras   
              CNV_UPDF_X          , &   ! <- ras
              MFD_SC            , &   ! <- shcu   
              QLDET_SC          , &   ! <- shcu   
              QIDET_SC          , &   ! <- shcu   
              SHLW_PRC3         , &   ! <- shcu   
              SHLW_SNO3         , &   ! <- shcu   
              UFRC_SC           , &   ! <- shcu 
              U1                , &
              V1                , & 
              TH1               , &              
              Q1                , &
              QLLS              , &
              QLCN              , &
              QILS              , &
              QICN              , &
              CLCN              , &
              CLLS              , &           
              CN_PRC2           , &            
              CN_ARFX           , &
              CN_SNR            , &
              CLDPARAMS         , &
              CLOUD_CTL%SCLMFDFR, &
              QST3              , &
              DZET              , &
              CNV_FRACTION      ,  &
              QDDF3             , &
                                ! Diagnostics
              RHX_X             , &
              REV_CN_X          , &
              RSU_CN_X          , &
              ACLL_CN_X,ACIL_CN_X   , &
              PFL_CN_X,PFI_CN_X     , &
              DLPDF_X,DIPDF_X,DLFIX_X,DIFIX_X,    &
              DCNVL_X, DCNVi_X,       &
              ALPHT_X, CFPDF_X, &
              DQRL_X,               &
              VFALLSN_CN_X,  &
              VFALLRN_CN_X,  &
              EVAPC_X , SUBLC_X ,  &
                                ! End diagnostics
            !!====2-Moment============
              CNV_FICE,  &
              CNV_NDROP_X, &
              CNV_NICE_X,  &
              SC_NDROP,  &
              SC_NICE,   &
              SC_ICE,    &
              NCPL, &
              NCPI, &
              PFRZ, &
              DNDCNV_X, &
              DNCCNV_X, &
              DT_RASP , &
              QRAIN_CN, & !grid av
              QSNOW_CN, &
              KCBL, LTS,  CONVPAR_OPTION)
       endif ! DOCLDMACRO


       IF(ADJUSTL(CONVPAR_OPTION) == 'GF') THEN   
            where (CNV_MFD .gt. 0.0) 
              CNV_NICE = CNV_NICE_X/CNV_MFD
              CNV_NDROP = CNV_NDROP_X/CNV_MFD
            end where
       END IF     

       TPREC = CN_PRC2 + LS_PRC2 + AN_PRC2 + SC_PRC2 + &
              CN_SNR  + LS_SNR  + AN_SNR + SC_SNR
      
if (MGVERSION .gt. 1.0) then 
      if(associated(TVQX2)) TVQX2    = SUM( ( Q1 +  QLLS + QLCN + QILS + QICN +  QRAIN +  QSNOW + QGRAUPEL)*MASS , 3 )  + TPREC*DT_MOIST -TVQ0
else

      if(associated(TVQX2)) TVQX2    = SUM( ( Q1 +  QLLS + QLCN + QILS + QICN)*MASS , 3 )  + TPREC*DT_MOIST -TVQ0
end if 


         TEMP    = TH1*PK

         DTDT_macro=  (TEMP-DTDT_macro)/DT_MOIST
         DQVDT_macro=(Q1-DQVDT_macro)/DT_MOIST
         DQLDT_macro=((QLCN+QLLS)-DQLDT_macro)/DT_MOIST
         DQIDT_macro=((QICN+QILS)-DQIDT_macro)/DT_MOIST
         DQADT_macro=((CLCN+CLLS)-DQADT_macro)/DT_MOIST

      
         !make sure QI , NI stay within T limits 
         call meltfrz_inst  (     &
              IM,JM,LM    , &
              TEMP              , &
              QLLS          , &
              QLCN         , &
              QILS           , &
              QICN          , &               
              NCPL         , &
              NCPI          )


         !=============================================End cloud macrophysics=====================================
         !======================================================================================================================
         !


         FQAI = 0.0
         FQAL = 0.0
         FQA  = 0.0 
         QCNTOT = QLCN+QICN
         QTOT   =  QCNTOT+ QLLS+QILS
         QL_TOT  = QLCN+QLLS 
         QI_TOT  = QICN+QILS   

         where (QTOT .gt. 0.0)
            FQA= min(max(QCNTOT/QTOT, 0.0), 1.0)    
         end where
	 
	 CFLIQ=0.0
         CFICE=0.0
       
         RAD_CF   = min(CLLS+CLCN, 1.0)

         WHERE (QTOT .gt. 0.0) 
            CFLIQ=RAD_CF*QL_TOT/QTOT
            CFICE=RAD_CF*QI_TOT/QTOT
         END WHERE


      
      
            INC_NUC = INC_NUC*PFRZ!!Accounts for ice crystal dilution after nucleation. 
            NHET_NUC = NHET_NUC*PFRZ!
      

         !==================================================================================================================
         !===============================================Two-moment stratiform microphysics ================================
         !================This is the implementation of the Morrison and Gettelman (2008) microphysics =====================
         !==================================================================================================================




         rhdfdar8   = 1.e-8_r8
         rhu00r8    = 0.95_r8
         ttendr8=0._r8
         qtendr8=0._r8
         cwtendr8=0._r8
         naair8=0.
         rndstr8 = 2.0e-7
         npccninr8 = 0.
         naconr8   = 0
         scale_ri =  1.3 ! scaling factor to account for the different definition of Ri in Chao and Suarez

        if ((RRTMG_SORAD .gt. 0.0) .or. (RRTMG_IRRAD .gt. 0.0)) then 
        scale_ri =  1.0
        end if 


         call MAPL_TimerOff(STATE,"---CLDMACRO", RC=STATUS)
         call MAPL_TimerOn(STATE,"---MGMICRO")
	 
	     !initialize MG variables
         nimmr8 = 0.0_r8
         cldfr8 = 0.0_r8 
         prectr8 = 0.0_r8 
	     precir8 = 0.0_r8
	     qctendr8 = 0.0_r8
	     qitendr8 = 0.0_r8
	     qvlatr8 = 0.0_r8
	     tlatr8 = 0.0_r8
	     nctendr8 = 0.0_r8
	     nitendr8 = 0.0_r8
	     effcr8 = 0.0_r8
	     effir8 = 0.0_r8
	     drout2r8 =0.0_r8
	     dsout2r8 = 0.0_r8
         dgout2r8 = 0.0_r8
	     qrout2r8 = 0.0_r8
	     qsout2r8 =0.0_r8
         qgout2r8 =0.0_r8
	     nrout2r8 = 0.0_r8
	     nsout2r8 =0.0_r8
         ngout2r8 =0.0_r8         
	     evapsnowr8 =0.0_r8
	     nevaprr8 =0.0_r8
	     cmeioutr8 =0.0_r8
	     bergsor8 =0.0_r8
	     mnucccor8 =0.0_r8
	     mnucctor8 =0.0_r8
	     homoor8 = 0.0_r8
	     mnuccror8 = 0.0_r8
	     pracsor8 = 0.0_r8
	     meltor8 =0.0_r8
	     qisedtenr8 =0.0_r8
	     bergor8 =0.0_r8
	     psacwsor8 = 0.0_r8
	     qcresor8 =0.0_r8
	     qiresor8 = 0.0_r8
	     praor8 =0.0_r8
	     prcor8 = 0.0_r8
	     prcior8 =0.0_r8
	     praior8 = 0.0_r8
	     msacwior8 =0.0_r8
	     frzrdtr8 =0.0_r8
	     meltsdtr8 = 0.0_r8
	     nnucctor8 =0.0_r8
	     nnucccor8 = 0.0_r8
	     nnuccdor8 =0.0_r8
	     nsacwior8 =0.0_r8
	     nsubior8 = 0.0_r8
	     npraior8 =0.0_r8
	     nprcior8 =0.0_r8
	     npccnor8 = 0.0_r8
	     npsacwsor8 =0.0_r8
	     npraor8 =0.0_r8
	     nsubcor8 =0.0_r8
	     nprc1or8 =0.0_r8
         rndstr8 = 2.0e-7
         naconr8   = 0.
     
         lflxr8 = 0.0_r8             
         iflxr8 = 0.0_r8
         rflxr8 = 0.0_r8
         sflxr8 = 0.0_r8
         gflxr8 = 0.0_r8    

         frzcntr8 =0.0_r8 
         qrtendr8 =  0.0_r8
         nrtendr8 =  0.0_r8
         qstendr8 =  0.0_r8
         nstendr8 =  0.0_r8
     
         qgtendr8 =  0.0_r8
         ngtendr8 =  0.0_r8
     
         accre_enhanr8= 1.0_r8 
         AN_PRC2     = 0. !prectr8(1)
         AN_SNR      = 0. !precir8(1)
         AN_ARFX     = 0. !maxval( cldfr8(1,1:LM) )    
         PFL_LS_X = 0.0
         PFI_LS_X= 0.0
         QCVAR_EXP = 2.0
         do J=1,JM
            do I=1,IM

              
	           kbmin =1 	   
               npccninr8  = 0.0
               naair8     = 0.0
               omegr8     = 0.0
               rndstr8 = 2.0e-7
               naconr8   = 0.

               cldfr8(1,1:LM)      = RAD_CF(I,J,1:LM) !Assume minimum overlap 
             
              ! liqcldfr8(1, 1:LM)  = cldfr8(1,1:LM) 
              ! icecldfr8(1, 1:LM)  = cldfr8(1,1:LM) 


                  
             
               liqcldfr8(1, 1:LM)  = cldfr8(1,1:LM) 
               icecldfr8(1, 1:LM)  = cldfr8(1,1:LM)  ! this is better to avoid removing liq clouds in the high lats
             !  liqcldfr8(1, 1:LM)  = CFLIQ(I, J,1:LM) 
             !  icecldfr8(1, 1:LM)  = CFICE(I, J,1:LM) 
	           cldor8           = cldfr8  
               ter8(1,1:LM)        = TEMP(I,J,1:LM)
               qvr8(1,1:LM)        = Q1(I,J,1:LM)

               qcr8(1,1:LM)        = QL_TOT(I,J,1:LM)
               qir8(1,1:LM)        = QI_TOT(I,J,1:LM)
               ncr8(1,1:LM)        = MAX(NCPL(I,J,1:LM), 0.0) 
               nir8(1,1:LM)        = MAX(NCPI(I,J,1:LM), 0.0) 

               ! Nucleation variables 
               naair8(1, 1:LM)     = INC_NUC(I, J, 1:LM)
               npccninr8(1, 1:LM)  = CDNC_NUC(I, J, 1:LM)

               where  ((naair8-ncr8  .gt. 1.0e3)) ! add cloud fraction if nucleation is happening 2018
                   icecldfr8 = max(0.05,  icecldfr8)
               end where 
             

               where (cldfr8(1,:) .ge. 0.001) 
                  nimmr8(1, 1:LM)     = MIN(DNHET_IMM(I, J, 1:LM), ncr8(1,1:LM)/cldfr8(1,1:LM)/DT_MOIST) !tendency    
               elsewhere 
                  nimmr8(1, 1:LM)   = 0.0 
               end where
               
               nhet_depr8(1, 1:LM) = NHET_DEP(I, J, 1:LM)/DT_MOIST !becomes a tendency (could be done a bit better)
               nbincontactdust = 1


               DO K=kbmin, LM
                  AeroAux = AeroProps(I, J, K)
                  ! Get dust properties for contact ice nucleation 
                  call getINsubset(1,  AeroAux, AeroAux_b) 
                  naux = AeroAux_b%nmods 
                  if (nbincontactdust  .lt. naux) then 
                     nbincontactdust =  naux 
                  end if
                  naconr8(1, K, 1:naux) =  AeroAux_b%num(1:naux)
                  rndstr8( 1, K, 1:naux)=AeroAux_b%dpg(1:naux)/2.0

                  ! Get black carbon properties for contact ice nucleation      
                  call init_Aer(AeroAux_b)   
                  call getINsubset(2,  AeroAux, AeroAux_b)          
                  nsootr8  (1, K)  = sum(AeroAux_b%num) !
                  naux = AeroAux_b%nmods 
                  rnsootr8 (1, K)  = sum(AeroAux_b%dpg(1:naux))/naux
               END DO


               pdelr8(1,1:LM)  = PLE(I,J,1:LM) - PLE(I,J,0:LM-1)  
               rpdelr8      = 1./pdelr8 
               pintr8(1,1:LM+1) = PLE(I,J,0:LM)  
               plevr8(1,1:LM)      = 100.*PLO(I,J,1:LM)
               zmr8(1,1:LM)        = ZLO(I,J,1:LM)     
               kkvhr8(1,1:LM+1) = KH(I,J,0:LM)  
               ficer8 = qir8 /( qcr8+qir8 + 1.e-10 )  
               omegr8(1,1:LM)=WSUB(I, J, 1:LM)
               
	       
               !Tuning factors
               disp_liu = LIU_MU
               ui_scale = UISCALE
               urscale  = URSCALE
               ts_autice = DT_R8*TS_AUTO_ICE 
               dcrit = DCRIT_
               if (MTIME .le. 0.0) then 
                   mtimesc  = DT_MOIST
               else               
                  mtimesc=MTIME
               end if 
  
  !!!!================Estimate qcvar following Xie and Zhang, JGR, 2015
                HMOIST_950 = 0.0
                HSMOIST_500 = 0.0
                              
                 IF (PLO(I, J, LM) .le. 500.0) then                                        
                      qcvarr8  = 2.0
                 ELSEIF (PLO(I, J, LM) .lt. 950.0) then 
                   
                    DO K=LM, 1, -1       
                         if (PLO(I,J,K) .lt. 500.0) exit  
                         HSMOIST_500 = MAPL_CP*TEMP(I, J, K) + GZLO(I, J, K) + QST3(I, J, K)*MAPL_ALHL
                    END DO 
                            
                      HMOIST_950 = MAPL_CP*TEMP(I, J, LM) + GZLO(I, J, LM) + Q1(I, J, LM)*MAPL_ALHL               
                      SINST = (HMOIST_950 -  HSMOIST_500)/(PLO(I,J,LM)*100.0- 50000.0)                   
                  else
                     DO K=LM, 1, -1       
                         if (PLO(I,J,K) .lt. 500.0) exit  
                         HSMOIST_500 = MAPL_CP*TEMP(I, J, K) + GZLO(I, J, K) + QST3(I, J, K)*MAPL_ALHL
                    END DO 

                    DO K=LM, 1, -1       
                     if (PLO(I,J,K) .lt. 950.0) exit  
                     HMOIST_950 = MAPL_CP*TEMP(I, J, K) + GZLO(I, J, K) + Q1(I, J, K)*MAPL_ALHL
                    END DO                                          
                     SINST = (HMOIST_950 -  HSMOIST_500)/45000.0                  
               
                   end if  
               
                  xscale = (36000.0/imsize)**(-0.666)
                  qcvarr8 =  0.67 -0.38*SINST +  4.96*xscale - 8.32*SINST*xscale  
                  qcvarr8 = min(max(qcvarr8, 0.5), 50.0)
                  if (associated(QCVAR_EXP)) QCVAR_EXP(I, J) = real(qcvarr8)
                  relvarr8 = qcvarr8
                  
                
               ! for MG23 (initial values)     
                        frzimmr8 =  nimmr8
                        frzcntr8 = nimmr8*0.0  
                        frzdepr8 = nhet_depr8
                        qrr8(1, 1:LM)     =  QRAIN(I, J,1:LM)
                        qsr8(1, 1:LM)     =  QSNOW(I, J,1:LM)
                        qgr8(1, 1:LM)     =  QGRAUPEL(I, J,1:LM)                        
                        nrr8(1, 1:LM)     =  NRAIN(I, J,1:LM)
                        nsr8(1, 1:LM)     =  NSNOW(I, J,1:LM)
                        ngr8(1, 1:LM)     =  NGRAUPEL(I, J,1:LM)                         
                        qsatfacr8 = 1.0                        
                        SCICE_tmp(1, 1:LM)  =  SC_ICE(I, J, 1:LM)
                        FQA_tmp(1, 1:LM)  = FQA(I, J, 1:LM) 
                        ALPH_tmp(1, 1:LM)  = ALPHT_X(I, J, 1:LM)
                        
    !                     if (0) then 
   ! print *, '=========before mG=========='
  
                      DO NAUX = 1, LM
                       
                         if (TEMP(I,J,NAUX) .lt. 150.0) then 
                          print *, '========beforemg========'
                          print *,  I, J, NAUX, TEMP(I,J,NAUX)
                        end if 
                        
                        if (isnan(TEMP (I,J,NAUX))) then 
                          print *, '========beforemg========'
                          print *,  I, J, NAUX, 'tnan'
                        end if 
                        
                     end do 


 
     
     
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
  !CALLS to MG versions
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  
!!!Call to MG microphysics. Lives in cldwat2m_micro.F90

               
   if (MGVERSION < 2.0)  then          
               
               call set_qcvar (qcvarr8)
               
               call mmicro_pcond (                           &
                    ncolmicro, ncolmicro, dt_r8, ter8, ttendr8,                   &
                    ncolmicro, LM ,                                      &
                    qvr8, qtendr8, cwtendr8, qcr8, qir8,          &
                    ncr8, nir8, plevr8, pdelr8, cldfr8,           &
                    liqcldfr8, icecldfr8,                         &
                    cldor8, pintr8, rpdelr8, zmr8, omegr8,        &
                    rate1ord_cw2pr,                               &  ! <= actually an output
                    naair8, npccninr8, rndstr8,naconr8,           &
                    rhdfdar8, rhu00r8, ficer8,                    & 
                    !                          ** outputs **
                    tlatr8, qvlatr8,        &
                    qctendr8, qitendr8, nctendr8, nitendr8, effcr8,    &
                    effc_fnr8, effir8, prectr8, precir8,             &  
                    nevaprr8, evapsnowr8,      &
                    prainr8, prodsnowr8,       &
                    cmeoutr8, deffir8, pgamradr8, &
                    lamcradr8,qsout2r8,dsout2r8, qrout2r8,drout2r8, &
                    qcsevapr8,qisevapr8,   &
                    qvresr8,cmeioutr8, &
                    vtrmcr8,vtrmir8,   &
                    qcsedtenr8,qisedtenr8, &
                    praor8,prcor8,mnucccor8, &
                    mnucctor8,msacwior8,psacwsor8,&
                    bergsor8,bergor8,meltor8, &
                    homoor8,qcresor8,prcior8, &
                    praior8,qiresor8,  &
                    mnuccror8,pracsor8, &
                    meltsdtr8,frzrdtr8, ncalr8, ncair8, mnuccdor8, nnucctor8, &
                    nsoutr8, nroutr8, nimmr8, disp_liu, &
                    nsootr8, rnsootr8, ui_scale, dcrit, mtimesc, &
                    nnuccdor8, nnucccor8, nsacwior8, nsubior8, nprcior8, &
                    npraior8, npccnor8, npsacwsor8, nsubcor8, npraor8, nprc1or8,  nbincontactdust, &
                    ts_autice, rflxr8, sflxr8)

    else ! MG2/3
        
         call  micro_mg_tend_interface ( DT_MICRO, INT(CLDPARAMS%PDFSHAPE), ALPH_tmp, SCICE_tmp, FQA_tmp, &
                                        CNV_FRACTION(I, J), SNOMAS(I, J), FRLANDICE(I, J), FRLAND(I, J), & 
                             ncolmicro,             LM,               dt_r8,       & 
                             ter8,                            qvr8,                              &
                             qcr8,                          qir8,                          &
                             ncr8,                          nir8,                          &
                             qrr8,                          qsr8,                          &
                             nrr8,                          nsr8,                          &
                             qgr8,                          ngr8,                         &
                             relvarr8,                     accre_enhanr8,                  &
                             plevr8,                       pdelr8,                         &
                             cldfr8,               liqcldfr8,            icecldfr8,  qsatfacr8,          &
                             qcsinksum_rate1ordr8,                                         &
                             naair8,                         npccninr8,                        &
                             rndstr8,                        naconr8,                        &
                             tlatr8,                         qvlatr8,                        &
                             qctendr8,                       qitendr8,                       &
                             nctendr8,                       nitendr8,                       &
                             qrtendr8,                       qstendr8,   qgtendr8,                     &
                             nrtendr8,                       nstendr8,   ngtendr8,                   &
                             effcr8,               effc_fnr8,            effir8,               &
                             sadicer8,                       sadsnowr8,                      &
                             prectr8,                        precir8,                        &
                             nevaprr8,                       evapsnowr8,                     &
                             am_evp_str8,                                                  &
                             prainr8,                        prodsnowr8,                     &
                             cmeoutr8,                       deffir8,                        &
                             pgamradr8,                      lamcradr8,                      &
                             qsoutr8,                        dsoutr8,                        &
                             qgoutr8,     ngoutr8,           dgoutr8,                        &
                             lflxr8,               iflxr8,   gflxr8,                           &
                             rflxr8,               sflxr8,    qroutr8,          &
                             reff_rainr8,                    reff_snowr8, reff_graur8,        &
                             qcsevapr8,            qisevapr8,            qvresr8,              &
                             cmeioutr8,            vtrmcr8,              vtrmir8,              &
                             umrr8,                          umsr8,                          &
                             umgr8,                          qgsedtendr8,                    &    
                             qcsedtenr8,                     qisedtenr8,                     &
                             qrsedtenr8,                     qssedtenr8,                     &
                             praor8,                       prcor8,                       &
                             mnucccor8,          mnucctor8,          msacwior8,          &
                             psacwsor8,          bergsor8,           bergor8,            &
                             meltor8,                      homoor8,                      &
                             qcresor8,           prcior8,            praior8,            &
                             qirestotr8,           mnuccrtotr8,          mnuccritotr8, pracstotr8,           &                           
                             meltsdtr8,         frzrdtr8,          mnuccdor8,          &
                             pracgtotr8,           psacwgtotr8,          pgsacwtotr8,          &
                             pgracstotr8,          prdgtotr8,           &
                             qmultgtotr8,          qmultrgtotr8,         psacrtotr8,           &
                             npracgtotr8,          nscngtotr8,           ngracstotr8,          &
                             nmultgtotr8,          nmultrgtotr8,         npsacwgtotr8,         & 
                             nroutr8,                            nsoutr8,                        &
                             reflr8,               areflr8,              areflzr8,             &
                             freflr8,              csrflr8,              acsrflr8,             &
                             fcsrflr8,                       rercldr8,                       &
                             ncair8,                         ncalr8,                         &
                             qrout2r8,                       qsout2r8,                       &
                             nrout2r8,                       nsout2r8,                       &
                             drout2r8,                       dsout2r8,                       &
                             qgout2r8,     ngout2r8,         dgout2r8,   freqgr8,                     &
                             freqsr8,                        freqrr8,                        &
                             nficer8,                        qcratr8,                        &
!                             errstring, & ! Below arguments are "optional" (pass null pointers to omit).
                             tnd_qsnow,          tnd_nsnow,          re_ice,    &
                             prer_evap, &
                             frzimmr8,             frzcntr8,              frzdepr8,  & ! contact is not passed since it depends on the droplet size dist
                             nsootr8, rnsootr8,  & ! soot for contact IN
                             npccnor8, npsacwsor8,npraor8,nsubcor8, nprc1or8, &  ! Number tendencies for liquid
                             npraior8, nnucctor8, nnucccor8, nnuccdor8, nsubior8, nprcior8, nsacwior8,  &  ! Number tendencies for ice
                             ts_autice, ui_scale, dcrit, disp_liu, nbincontactdust, urscale)



    end if 

        IF (MGVERSION > 1.0) then 

#ifdef FAILS 
                  QRAIN(I,J,1:LM)  = max(QRAIN(I,J,1:LM) + REAL(qrtendr8(1, 1:LM)*DT_R8), 0.0) ! grid average 
                  QSNOW(I,J,1:LM)  = max(QSNOW(I,J,1:LM) + REAL(qstendr8(1, 1:LM)*DT_R8), 0.0) ! grid average                     
                  NRAIN(I,J,1:LM)  = max(NRAIN(I,J,1:LM) + REAL(nrtendr8(1, 1:LM)*DT_R8), 0.0)
                  NSNOW(I,J,1:LM)  = max(NSNOW(I,J,1:LM) + REAL(nstendr8(1, 1:LM)*DT_R8), 0.0)                  
                  CLDREFFR(I,J,1:LM) = REAL(reff_rainr8(1, 1:LM))        
                  CLDREFFS(I,J,1:LM) = REAL(reff_snowr8(1, 1:LM))/scale_ri 
                  CLDREFFG(I,J,1:LM) = REAL(reff_graur8(1, 1:LM))/scale_ri  
                  DQRL_X(I,J,1:LM)   = REAL(   qrtendr8(1, 1:LM)) !rain mixing ratio tendency from micro
                  
              if (adjustl(CLDMICRO)=="MG3") then                   
                  QGRAUPEL(I,J,1:LM)  = max(QGRAUPEL(I,J,1:LM) + REAL(qgtendr8(1, 1:LM)*DT_R8), 0.0) ! grid average 
                  NGRAUPEL(I,J,1:LM)  = max(NGRAUPEL(I,J,1:LM) + REAL(ngtendr8(1, 1:LM)*DT_R8), 0.0)
               else
                  QGRAUPEL(I,J,1:LM)  = qgout2r8(1, 1:LM) ! grid average                    
                  NGRAUPEL(I,J,1:LM)  = ngout2r8(1, 1:LM) ! grid average 
               end if                         
#else
                   QRAIN   (I,J,1:LM)  = max(REAL(qroutr8(1, 1:LM)), 0.0)
                   QSNOW   (I,J,1:LM)  = max(REAL(qsoutr8(1, 1:LM)), 0.0)
                   QGRAUPEL(I,J,1:LM)  = max(REAL(qgoutr8(1, 1:LM)), 0.0)
                   NRAIN   (I,J,1:LM)  = max(REAL(nroutr8(1, 1:LM)), 0.0)
                   NSNOW   (I,J,1:LM)  = max(REAL(nsoutr8(1, 1:LM)), 0.0)
                   NGRAUPEL(I,J,1:LM)  = max(REAL(ngoutr8(1, 1:LM)), 0.0)
                   CLDREFFR(I,J,1:LM)  = REAL(reff_rainr8(1, 1:LM))
                   CLDREFFS(I,J,1:LM)  = REAL(reff_snowr8(1, 1:LM))/scale_ri
                   CLDREFFG(I,J,1:LM)  = REAL(reff_graur8(1, 1:LM))/scale_ri
                   DQRL_X(I,J,1:LM)    = REAL(qroutr8(1, 1:LM)/DT_R8) !rain mixing ratio tendency from micro
#endif
            
        else
                    
                   QRAIN(I,J,1:LM)  = max(REAL(qrout2r8(1, 1:LM)), 0.0) ! grid average 
                   QSNOW(I,J,1:LM)  = max(REAL(qsout2r8(1, 1:LM)), 0.0)                      
                   NRAIN(I,J,1:LM)  = max(REAL(nrout2r8(1, 1:LM)), 0.0)
                   NSNOW(I,J,1:LM)  = max(REAL(nsout2r8(1, 1:LM)), 0.0)
                   CLDREFFR(I,J,1:LM) = REAL(drout2r8(1, 1:LM))/2.0        
                   CLDREFFS(I,J,1:LM) = REAL(dsout2r8(1, 1:LM))/2.0/scale_ri
                   DQRL_X(I,J,1:LM)   = REAL(qrout2r8(1, 1:LM)/DT_R8) !rain mixing ratio tendency from micro
                 
         end if          
         
         
  
               PFL_LS_X(I, J, 1:LM) = rflxr8(1, 1:LM) !+ lflxr8(1, 1:LM)
               PFI_LS_X(I, J, 1:LM) = sflxr8(1, 1:LM) !+ gflxr8(1, 1:LM) +  iflxr8(1, 1:LM)
              
               !Update state after microphysisc
               LS_PRC2(I,J)     = max(1000.*REAL((prectr8(1)-precir8(1))), 0.0)
               LS_SNR(I,J)      = max(1000.*REAL(precir8(1)), 0.0)          
               QL_TOT(I,J,1:LM) = max(QL_TOT(I,J,1:LM)   + REAL(qctendr8(1,1:LM)) * DT_R8, 0.0)
               QI_TOT(I,J,1:LM) = max(QI_TOT(I,J,1:LM)   + REAL(qitendr8(1,1:LM)) * DT_R8, 0.0)    
               Q1(I,J,1:LM)   = MAX(Q1(I,J,1:LM)     + REAL(qvlatr8(1,1:LM)) * DT_R8, 0.0)
               TEMP(I,J,1:LM) = TEMP(I,J,1:LM)   + REAL(tlatr8(1,1:LM)) * DT_R8 / (MAPL_CP)  
               NCPL(I,J,1:LM) = MAX(NCPL(I,J,1:LM)   + REAL(nctendr8(1,1:LM)) * DT_R8, 0.0) 
               NCPI(I,J,1:LM) = MAX(NCPI(I,J,1:LM)   + REAL(nitendr8(1,1:LM)) * DT_R8, 0.0)  
	                   
                       
                        DO NAUX = 1, LM
                       
                        
                         if (TEMP(I,J,NAUX) .lt. 150.0) then 
                          print *, '========aftermg========'
                          print *,  I, J, NAUX, TEMP(I,J,NAUX)
                        end if 
                        
                        
                        if (isnan(TEMP (I,J,NAUX))) then 
                          print *, '========aftermg========'
                          print *,  I, J, NAUX, 'tnan'
                        end if 
                        
                        if (isnan(Q1 (I,J,NAUX))) then 
                          print *, '========aftermg========'
                          print *,  I, J, NAUX, 'qnan'
                        end if 
                        
                         
                        if (isnan(QLLS (I,J,NAUX))) then 
                          print *, '========aftermg========'
                          print *,  I, J, NAUX, 'qnan'
                        end if 
                         
                        if (isnan(QLCN (I,J,NAUX))) then 
                          print *, '========aftermg========'
                          print *,  I, J, NAUX, 'qnan'
                        end if 
                         
                        if (isnan(CLLS (I,J,NAUX))) then 
                          print *, '========aftermg========'
                          print *,  I, J, NAUX, 'qnan'
                        end if 
                        
                          if (isnan(CLCN (I,J,NAUX))) then 
                          print *, '========aftermg========'
                          print *,  I, J, NAUX, 'qnan'
                        end if 
                        
                          if (isnan(QRAIN (I,J,NAUX))) then 
                          print *, '========aftermg========'
                          print *,  I, J, NAUX, 'qnan'
                        end if 
                          if (isnan(QSNOW (I,J,NAUX))) then 
                          print *, '========aftermg========'
                          print *,  I, J, NAUX, 'qnan'
                        end if 
                        
                          if (isnan(QGRAUPEL (I,J,NAUX))) then 
                          print *, '========aftermg========'
                          print *,  I, J, NAUX, 'qnan'
                        end if 
                        
                         
                        if (isnan(NCPL (I,J,NAUX))) then 
                          print *, '========aftermg========'
                          print *,  I, J, NAUX, 'qnan'
                        end if 
                        
                          if (isnan(NCPI (I,J,NAUX))) then 
                          print *, '========aftermg========'
                          print *,  I, J, NAUX, 'qnan'
                        end if 
                        
                          if (isnan(NRAIN (I,J,NAUX))) then 
                          print *, '========aftermg========'
                          print *,  I, J, NAUX, 'qnan'
                        end if 
                          if (isnan(NSNOW (I,J,NAUX))) then 
                          print *, '========aftermg========'
                          print *,  I, J, NAUX, 'qnan'
                        end if 
                        
                          if (isnan(NGRAUPEL (I,J,NAUX))) then 
                          print *, '========aftermg========'
                          print *,  I, J, NAUX, 'qnan'
                        end if 
                        
                        
                        
                     end do 
                       

               LS_ARFX(I,J)     = maxval( REAL(cldfr8(1,1:LM)) )
                            
               CLDREFFL(I,J,1:LM) = max(REAL(effcr8(1,1:LM))*1.0e-6, 1.0e-6)             
               CLDREFFI(I,J,1:LM) = max(REAL(effir8(1,1:LM))*1.0e-6, 1.0e-6)/scale_ri !scale to match the Dge definition of Fu 1996                    

               ! diagnostics from the microphysics********************

               RSU_LS_X(I,J,1:LM) = REAL(evapsnowr8(1, 1:LM))                    !Snow evap   
               REV_LS_X(I,J,1:LM) = REAL(nevaprr8(1, 1:LM) )                        !rain evap
               SUBLC_X(I,J,1:LM) = REAL(cmeioutr8(1, 1:LM)) + SUBLC_X(I,J,1:LM)     ! Ice subl already grid -ave
               BERGS(I,J,1:LM) = REAL(bergsor8(1,1:LM))                               ! Snow Bergeron
               FRZ_TT_X(I,J,1:LM) = REAL(mnucccor8(1,1:LM)+ mnucctor8(1,1:LM) + homoor8(1,1:LM))!ice mixing ratio from nucleation (hom+het)
               FRZ_PP_X(I,J,1:LM) = REAL(mnuccror8(1, 1:LM) + pracsor8(1, 1:LM)) !freezing of rain (hom+ het freezing and accretion by snow) 
               MELT(I,J,1:LM) = REAL(meltor8(1,1:LM))                            !melting of cloud ice  and snow
               SDM_X(I,J,1:LM) = REAL(qisedtenr8(1, 1:LM))                    ! ice sed   
               EVAPC_X(I,J,1:LM) = REAL(qcsevapr8(1,1:LM) )  +  EVAPC_X(I,J,1:LM)   ! cloud evap
               BERG(I,J,1:LM)  =  REAL(bergor8(1,1:LM))                          ! Bergeron process  
               ACIL_LS_X(I,J,1:LM) =REAL(psacwsor8(1, 1:LM) + msacwior8(1, 1:LM))   !Acc + collection of cloud  by snow
               QCRES(I,J,1:LM) =REAL(qcresor8(1, 1:LM) )             !residual liquid condensation
               QIRES(I,J,1:LM) =REAL(qiresor8(1, 1:LM))                  !residual ice condensation   

               ACLL_LS_X(I,J,1:LM) =REAL(praor8(1, 1:LM) )          ! Acc cloud by rain    
               AUT_X(I,J,1:LM) = REAL(prcor8(1, 1:LM))                  ! Aut liquid 
               AUTICE(I,J,1:LM) = REAL(prcior8(1, 1:LM))                !Aut  ice    
               ACIL_AN_X(I,J,1:LM) = REAL(praior8(1, 1:LM))           !Acc ice  by snow
               ACLL_AN_X(I,J,1:LM) = REAL(msacwior8(1, 1:LM))   !HM process

               FRZPP_LS(I,J,1:LM) = REAL(frzrdtr8(1, 1:LM)) / MAPL_CP !precip freezing latent heat rate
               SNOWMELT_LS(I,J,1:LM) =REAL(meltsdtr8(1, 1:LM))/ MAPL_CP !melting of snow latent heat rate

               !diagnostics for number concentration budget (all grid-average, Kg-1 s-1)
              

               DNDCCN(I,J,1:LM)       = REAL(npccnor8(1,1:LM))   !droplet number tendency from CCN activation
               DNDACRLS(I,J,1:LM)     = REAL(npsacwsor8(1,1:LM)) !droplet number tendency from accretion by snow
               DNDACRLR(I,J,1:LM)     = REAL(npraor8(1,1:LM))    !droplet number tendency from accretion by rain
               DNDEVAPC(I,J,1:LM)     = REAL(nsubcor8(1,1:LM))   !droplet number tendency from evaporation 
               DNDAUTLIQ(I,J,1:LM)    = REAL(nprc1or8(1,1:LM))   !droplet number tendency from autoconversion
	       
	       DNCACRIS (I,J,1:LM)     = REAL(npraior8(1,1:LM))  !ice number tendency from accretion by snow
	       DNHET_CT(I,J,1:LM)      = REAL(nnucctor8(1,1:LM)) !ice number tendency from contact IN  
               DNHET_IMM(I,J,1:LM)     = REAL(nnucccor8(1,1:LM)) !ice number tendency from immersion IN    
               DNCNUC(I,J,1:LM)        = REAL(nnuccdor8(1,1:LM)) !ice number tendency from nucleation on aerosol             
               DNCSUBL (I,J,1:LM)      = REAL(nsubior8(1,1:LM))  !ice number tendency from sublimation               
               DNCAUTICE (I,J,1:LM)    = REAL(nprcior8(1,1:LM))  !ice number tendency from autoconversion
	       DNCHMSPLIT(I,J,1:LM)    = REAL(nsacwior8(1,1:LM)) !ice number tendency from H-M process


               ! Total tendencies

               
               DQVDT_micro(I,J,1:LM)   = REAL(qvlatr8(1,1:LM))  
               DQIDT_micro(I,J,1:LM)   = REAL(qitendr8(1,1:LM))   
               DQLDT_micro(I,J,1:LM)   = REAL(qctendr8(1,1:LM) )    
               DTDT_micro(I,J,1:LM)    = REAL((tlatr8(1,1:LM) / MAPL_CP))
               DNDCCN(I,J,1:LM)       = REAL(npccnor8(1,1:LM))   !droplet number tendency from CCN activation
               DNDACRLS(I,J,1:LM)     = REAL(npsacwsor8(1,1:LM) )!droplet number tendency from accretion by snow
               DNDACRLR(I,J,1:LM)     = REAL(npraor8(1,1:LM))    !droplet number tendency from accretion by rain
               DNDEVAPC(I,J,1:LM)     = REAL(nsubcor8(1,1:LM))   !droplet number tendency from evaporation 
               DNDAUTLIQ(I,J,1:LM)    = REAL(nprc1or8(1,1:LM))   !droplet number tendency from autoconversion


            enddo !I
         enddo !J
         !============================================Finish 2-moment micro implementation===========================


    !update water tracers
         QLCN=QL_TOT*FQA
         QLLS=QL_TOT-QLCN
         QICN=QI_TOT*FQA
         QILS=QI_TOT-QICN
         PFL_AN_X(:,:,1:LM) = PFL_LS_X(:,:,1:LM) * FQA
         PFL_LS_X(:,:,1:LM) = PFL_LS_X(:,:,1:LM) - PFL_AN_X(:,:,1:LM)
         PFI_AN_X(:,:,1:LM) = PFI_LS_X(:,:,1:LM) * FQA
         PFI_LS_X(:,:,1:LM) = PFI_LS_X(:,:,1:LM) - PFI_AN_X(:,:,1:LM)
         QTOT= QICN+QILS+QLCN+QLLS

       TPREC = CN_PRC2 + LS_PRC2 + AN_PRC2 + SC_PRC2 + &
              CN_SNR  + LS_SNR  + AN_SNR + SC_SNR


   
         !============ Recalculate cloud fraction back in contact with the PDF and create new condensate if neccesary (Barahona et al., GMD, 2014)============
   !IF (MGVERSION <= 1.0) then

         !============ Put cloud fraction back in contact with the PDF and create new condensate if neccesary (Barahona et al., GMD, 2014)============


 DLPDF_X=  QLLS +QLCN
 DIPDF_X =  QILS +QICN

do K= 1, LM
   do J=1,JM
            do I=1,IM

		  call update_cld( &
        		 DT_MOIST                , &
        		 ALPHT_X(I, J, K)        , &
        		 INT(CLDPARAMS%PDFSHAPE) , &
        		 PLO(I, J, K)            , &
        		 Q1 (I, J, K)            , &
        		 QLLS(I, J, K)           , &
        		 QLCN(I, J, K)           , &
        		 QILS(I, J, K)           , &
        		 QICN(I, J, K)           , &
        		 TEMP(I, J, K)           , &
        		 CLLS(I, J, K)           , &
        		 CLCN(I, J, K)           , &
			 SC_ICE(I, J, K)         , &
			 NCPI(I, J, K)           , &
			 NCPL(I, J, K)           , &
			 INC_NUC(I, J, K)        , &
			 RHCmicro(I, J, K), &
             CNV_FRACTION(I, J), SNOMAS(I, J), FRLANDICE(I, J), FRLAND(I, J))
	      
           end do 
	end do
  end do 
  
  
   DLPDF_X=((QLLS+QLCN)  - DLPDF_X)/DT_MOIST
   DIPDF_X=((QILS+QICN)  - DIPDF_X)/DT_MOIST

         ! Make sure ice and liquid stay within T limits  


  call meltfrz_inst  (     &
              IM,JM,LM    , &
              TEMP              , &
              QLLS          , &
              QLCN         , &
              QILS           , &
              QICN          , &               
              NCPL         , &
              NCPI          )
	      	    
         RAD_CF =min(CLLS+CLCN, 1.0)


         TPREC = CN_PRC2 + LS_PRC2 + AN_PRC2 + SC_PRC2 + &
              CN_SNR  + LS_SNR  + AN_SNR + SC_SNR + RASPRCP

  
 !if(associated(TVQX2)) TVQX2    = SUM( ( Q1 +  QLLS + QLCN + QILS + QICN )*MASS , 3 )  + TPREC*DT_MOIST

     
         !=============================================End Stratiform cloud processes==========================================
         !======================================================================================================================
         call MAPL_TimerOff(STATE,"---MGMICRO", RC=STATUS)
         VERIFY_(STATUS)
         call MAPL_TimerOff(STATE,"--CLOUD_RUN",RC=STATUS)
         VERIFY_(STATUS)


         !Calculate CFICE and CFLIQ 

         CFLIQ=0.0
         CFICE=0.0
         QTOT= QICN+QILS+QLCN+QLLS
         QL_TOT = QLCN+QLLS
         QI_TOT = QICN+QILS

         WHERE (QTOT .gt. 1.0e-12) 
            CFLIQ=RAD_CF*QL_TOT/QTOT
            CFICE=RAD_CF*QI_TOT/QTOT
         END WHERE

         CFLIQ=MAX(MIN(CFLIQ, 1.0), 0.0)
         CFICE=MAX(MIN(CFICE, 1.0), 0.0)


         !======================================================================================================================
         !===========================Clean stuff and send it to radiation ======================================================
         !======================================================================================================================

         where (QI_TOT .le. 0.0)
            CFICE =0.0
            NCPI=0.0
            CLDREFFI = 36.0e-6
         end where

         where (QL_TOT .le. 0.0)
            CFLIQ =0.0
            NCPL  =0.0
            CLDREFFL = 14.0e-6
         end where

         WHERE  (RAD_CF > 1e-4)
            RAD_QL = min((QLLS+QLCN)/RAD_CF, 1.0e-3)
            RAD_QI = min((QILS+QICN)/RAD_CF, 1.0e-3) !
            RAD_QR =  QRAIN/RAD_CF  
            RAD_QS =  QSNOW/RAD_CF
            RAD_QG =  QGRAUPEL/RAD_CF
         ELSEWHERE 
            RAD_QL = 0.0         
            RAD_QI = 0.0
            RAD_QR = 0.0
            RAD_QS = 0.0
            RAD_QG = 0.0
         end where


         
          
       
         !Everything in-cloud for radiation============== 

         RAD_QV = MAX( Q1 , 0. )
         RAD_QL = MAX(MIN( RAD_QL , 0.001 ), 0.0)  ! Still a ridiculously large
         RAD_QI = MAX(MIN( RAD_QI , 0.001 ), 0.0)  ! value.
         RAD_QR = MAX(MIN( RAD_QR , 0.01 ), 0.0)  ! value.
         RAD_QS = MAX(MIN( RAD_QS , 0.01 ), 0.0)  ! value
         RAD_QG = MAX(MIN( RAD_QG , 0.01 ), 0.0)  ! value
         
         !=================================================================================
         !    Units conversion for diagnostics


         CFX =100.*PLO*r_air/TEMP !density times conversion factor
 
  
         !to m-3
         NCPL_VOL=NCPL*CFX !
         NCPI_VOL=NCPI*CFX
         CDNC_NUC=CDNC_NUC*CFX 
         INC_NUC =INC_NUC*CFX             
         CNV_NICE=CNV_NICE*CFX 
         CNV_NDROP=CNV_NDROP*CFX         
	 
        !to m-3 s-1

         DNHET_CT    = DNHET_CT*CFX 
         DNHET_IMM   = DNHET_IMM*CFX 
         DNCNUC      = DNCNUC*CFX 
         DNCHMSPLIT  = DNCHMSPLIT*CFX
         DNCSUBL     = DNCSUBL*CFX 
         DNCACRIS    = DNCACRIS*CFX 
         DNCAUTICE   = DNCAUTICE*CFX 
         DNCCNV      = DNCCNV_X*CFX

         DNDCCN       = DNDCCN*CFX   
         DNDACRLS     = DNDACRLS*CFX
         DNDACRLR     = DNDACRLR*CFX    
         DNDEVAPC     = DNDEVAPC*CFX   
         DNDAUTLIQ    = DNDAUTLIQ*CFX 
         DNDCNV       = DNDCNV_X*CFX

         !make grid averages

         !NHET_IMM=NHET_IMM*CFICE
         !INC_NUC = INC_NUC*CFICE
         !CDNC_NUC = CDNC_NUC*CFLIQ
         !NHET_NUC = NHET_NUC*CFICE
         !DUST_DEP = DUST_DEP*CFICE
         !DUST_IMM = DUST_IMM*CFICE


         !Grid average  volumetric  radius for comparison against field data
         WHERE  ((CFICE > 0.001) .and. (NCPI .gt. 1.0))
            RI_MASK =  CFICE*((QILS+QICN)/(800.0*1.333*MAPL_PI*NCPI))**0.333
         elsewhere
            RI_MASK=0.0
         end  where

         WHERE  ((CFLIQ > 0.001) .and. (NCPL .gt. 1.0)) 
            RL_MASK =  CFLIQ*((QLLS+QLCN)/(900.0*1.333*MAPL_PI*NCPL))**0.333
         ELSEWHERE
            RL_MASK = 0.0 
         END WHERE



         TH1 = TEMP / PK

  
         !Set rain water for radiation to 0 if preciprad flag is off (set to 0)
         if(CLDPARAMS%PRECIPRAD .eq. 0.) then
            RAD_QR = 0.
            RAD_QS = 0.
            RAD_QG = 0.      
         endif

         if (associated(QRTOT)) QRTOT = QRAIN
         if (associated(QSTOT)) QSTOT = QSNOW


         CLDREFFL = MAX(4.1e-6, CLDREFFL) !DONIF Limits according to MG2008-I 
         CLDREFFL = MIN(29.e-6, CLDREFFL)
         CLDREFFI = MAX(6.e-6, CLDREFFI)   
         CLDREFFI = MIN(89.e-6, CLDREFFI)  !maximum number for the correlation and modis sim 
  
         CLDREFFR = MAX(4.1e-6, CLDREFFR) 
         CLDREFFR = MIN(29.e-6, CLDREFFR)
         CLDREFFS = MAX(6.e-6, CLDREFFS)   
         CLDREFFS = MIN(89.e-6, CLDREFFS)  !maximum number for the correlation and modis sim   
         CLDREFFG = MAX(6.e-6, CLDREFFG)   
         CLDREFFG = MIN(89.e-6, CLDREFFG)  !maximum number for the correlation and modis sim 
          

         !===========================

         ! Diagnostic cloud top/base properties
         
         CLDREFFL_TOP_X = MAPL_UNDEF
         CLDREFFI_TOP_X =MAPL_UNDEF
         NCPL_TOP_X = MAPL_UNDEF
         NCPL_CLDBASEX = MAPL_UNDEF
         NCPI_TOP_X = MAPL_UNDEF
         kbmin = LM
         LWC_X = 0.0
         IWC_X = 0.0

         CFX =100.*PLO*r_air/TEMP
         DO I=1, IM
            DO J= 1 , JM

               cfaux(1, 1:LM) =  CFLIQ(I, J, 1:LM)
               call find_cldtop(1, LM, cfaux, kbmin)

               if (kbmin .ge. LM-1) then 
                  CLDREFFL_TOP_X (I, J)  = 8.0e-6 
                  NCPL_TOP_X (I, J)  = 0.0
               else

                  CLDREFFL_TOP_X (I, J)  = CLDREFFL(I, J,  kbmin) 
                  NCPL_TOP_X (I, J)  = NCPL_VOL(I, J,  kbmin)*CFX(I, J, kbmin) 
               end if


                call find_cldbase(1, LM, cfaux, kbmin)
                if (kbmin .gt. 10) then 
                  NCPL_CLDBASEX (I, J)  = NCPL_VOL(I, J,  kbmin)/max(cfaux(1, kbmin), 0.01) 
                end if
                
               cfaux(1, 1:LM) =CFICE(I, J, 1:LM)   
               call find_cldtop(1, LM, cfaux, kbmin)

               if (kbmin .ge. LM-1) then 
                  CLDREFFI_TOP_X (I, J)=20.0E-6
                  NCPI_TOP_X (I, J)  = 0.0
               else       
                  CLDREFFI_TOP_X (I, J)  = CLDREFFI(I, J,  kbmin) 
                  NCPI_TOP_X (I, J)  = NCPI_VOL(I, J,  kbmin)*CFX(I, J, kbmin)
               end if

            END DO
         END DO

    
      !=====Tune area cloud fraction of extended PBL clouds. Area cloud fraction may be different from Volume cloud fraction 
         FQA= 0.0
         where (RAD_CF .gt. 0.0)
              FQA =  CLCN/RAD_CF
         end where
         
       
 
         DO I=1, IM
            DO J = 1, JM    


               if (FRLAND(I, J) .lt. 0.1) then 

                  
                  cfc_aux =(TEMP(I, J, LM) - TMAXCFCORR)/2.0
                  cfc_aux =  min(max(cfc_aux,-20.0), 20.0)
                  cfc_aux=   1.0/(1.0+exp(-cfc_aux))
                  
                  DO K=LM-1, 2, -1
                     if ((RAD_CF(I, J, K) .gt. 0.01) .and. (RAD_CF(I, J, K) .lt. 0.99)) then  

                           USURF=1.0                        
                           USURF= (LTS_UP-LTS(I, J))/(LTS_UP -  CLDPARAMS%MIN_LTS) 
                           USURF=min(max(USURF, MIN_EXP), MAX_EXP)                            
                           
                           fracover=min(max((TEMP(I, J, K) -TMAXLL)/2.0, -20.0), 20.0)
                           fracover = 1.0/(1.0+exp(-fracover))
                          
                           USURF = USURF*fracover + 1.0-fracover   !only near the surface                                  
                           USURF = USURF*cfc_aux + 1.0-cfc_aux !only for the subtropics                          
                           USURF =  usurf+ (1.0-USURF)*CNV_FRACTION(I, J) !only non-convective
                                            
                           RAD_CF(I, J, K)=RAD_CF(I, J, K)**USURF								
                     END IF

                  END DO
               END IF

            end do
         end do
         
    
           CLCN   =  FQA*RAD_CF
           CLLS =  (1.0-FQA)*RAD_CF   

       WHERE (QTOT .gt. 1.0e-12) 
            CFLIQ=RAD_CF*QL_TOT/QTOT
            CFICE=RAD_CF*QI_TOT/QTOT
        END WHERE
         
      if(associated(CNV_FRC )) CNV_FRC  = CNV_FRACTION

         !=======================================================================

      end if 
      !=================End microphysics conditional===========================
      !=======================================================================


      !----------Other diagnostics

      QTOT = QLLS+QILS+QLCN+QICN
        
      if (associated(SCF)) then 
         WHERE (QTOT .gt. 1.0e-15)
            SCF=min(max((QLCN+QLLS)/QTOT, 0.0), 1.0)
         ELSEWHERE 
           SCF=MAPL_UNDEF
         END WHERE
      endif

      if (associated(SCF_ALL)) then
         WHERE (QRAIN+QTOT+QSNOW .gt. 1.0e-15)
            SCF_ALL=min(max((QRAIN+QLCN+QLLS)/(QRAIN+QSNOW+QTOT), 0.0), 1.0)
         ELSEWHERE 
            SCF_ALL=MAPL_UNDEF
         END WHERE
      endif

      if (associated(CLDBASEHGT)) then
         CLDBASEx = MAPL_UNDEF
         do i = 1,IM
           do j = 1,JM
             do k =  LM-1, 1, -1
               if (ZLE(i,j,k).gt.20000.) exit
               if ( ( RAD_CF(i,j,k) .ge. 1e-2 ) .and. ( QTOT(i,j,k) .ge. 1e-6 ) ) then
                 CLDBASEx(i,j)  = ZLE(i,j,k)
                 exit
               end if
             end do
           end do
         end do
         CLDBASEHGT = CLDBASEx
      end if


      ! Compute DBZ radar reflectivity
      if (associated(DBZ) .OR. associated(DBZ_MAX)) then
         if(adjustl(CLDMICRO)=="1MOMENT") then
           call CALCDBZ(DBZ3D,100*PLO,TEMP,Q1,QRN*RAD_CF,QSN*RAD_CF,QSN*RAD_CF,IM,JM,LM,1,0,0)
         else
           call CALCDBZ(DBZ3D,100*PLO,TEMP,Q1,QRAIN,QSNOW,QGRAUPEL,IM,JM,LM,1,0,0)
         endif
         if (associated(DBZ)) DBZ = DBZ3D
         if (associated(DBZ_MAX)) then
            DBZ_MAX=-9999.0
            DO K=1,LM
            DO J=1,JM
            DO I=1,IM
               DBZ_MAX(I,J) = MAX(DBZ_MAX(I,J),DBZ3D(I,J,K))
            END DO
            END DO
            END DO
         endif
      endif

      !ice water and liquid water content 

      CFX =100.*PLO*r_air/TEMP
      LWC_X = (QLLS+QLCN)*CFX
      IWC_X = (QILS+QICN)*CFX    

     IF(ADJUSTL(CONVPAR_OPTION) == 'GF') THEN
         REV_CN_X = REV_CN_GF
	 RSU_CN_X = RSU_CN_GF 
         ACLL_CN_X = 0.
         ACIL_CN_X = 0.
     ENDIF 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
      if (associated(LWC    ))   LWC      = LWC_X
      if (associated(IWC    ))   IWC      = IWC_X
      if (associated(CLDREFFI_TOP))   CLDREFFI_TOP    = CLDREFFI_TOP_X
      if (associated(CLDREFFL_TOP))   CLDREFFL_TOP    = CLDREFFL_TOP_X
      if (associated(NCPL_TOP))   NCPL_TOP    = NCPL_TOP_X
      if (associated(NCPI_TOP))   NCPI_TOP    = NCPI_TOP_X
      if (associated(NCPL_CLDBASE))   NCPL_CLDBASE    = NCPL_CLDBASEX

      IF (ASSOCIATED(RSU_CN)) RSU_CN = RSU_CN_X
      IF (ASSOCIATED(RSU_AN)) RSU_AN = RSU_AN_X
      IF (ASSOCIATED(RSU_LS)) RSU_LS = RSU_LS_X
      IF (ASSOCIATED(RSU_SC)) RSU_SC = RSU_SC_X

      IF (ASSOCIATED(REV_CN)) REV_CN = REV_CN_X
      IF (ASSOCIATED(REV_AN)) REV_AN = REV_AN_X
      IF (ASSOCIATED(REV_LS)) REV_LS = REV_LS_X
      IF (ASSOCIATED(REV_SC)) REV_SC = REV_SC_X

      IF (ASSOCIATED(ACIL_CN)) ACIL_CN = ACIL_CN_X
      IF (ASSOCIATED(ACIL_AN)) ACIL_AN = ACIL_AN_X
      IF (ASSOCIATED(ACIL_LS)) ACIL_LS = ACIL_LS_X
      IF (ASSOCIATED(ACIL_SC)) ACIL_SC = ACIL_SC_X

      IF (ASSOCIATED(ACLL_CN)) ACLL_CN = ACLL_CN_X
      IF (ASSOCIATED(ACLL_AN)) ACLL_AN = ACLL_AN_X
      IF (ASSOCIATED(ACLL_LS)) ACLL_LS = ACLL_LS_X
      IF (ASSOCIATED(ACLL_SC)) ACLL_SC = ACLL_SC_X

      IF (ASSOCIATED(PFI_CN)) PFI_CN = PFI_CN_X
      IF (ASSOCIATED(PFI_AN)) PFI_AN = PFI_AN_X
      IF (ASSOCIATED(PFI_LS)) PFI_LS = PFI_LS_X
      IF (ASSOCIATED(PFI_SC)) PFI_SC = PFI_SC_X

      IF (ASSOCIATED(PFL_CN)) PFL_CN = PFL_CN_X
      IF (ASSOCIATED(PFL_AN)) PFL_AN = PFL_AN_X
      IF (ASSOCIATED(PFL_LS)) PFL_LS = PFL_LS_X
      IF (ASSOCIATED(PFL_SC)) PFL_SC = PFL_SC_X

!!$      IF (ASSOCIATED(LIQANMOVE)) LIQANMOVE = LIQANMOVE_X
!!$      IF (ASSOCIATED(ICEANMOVE)) ICEANMOVE = ICEANMOVE_X
!!$      IF (ASSOCIATED(DANCLD)) DANCLD = DANCLD_X
!!$      IF (ASSOCIATED(DLSCLD)) DLSCLD = DLSCLD_X
!!$      IF (ASSOCIATED(CURAINMOVE)) CURAINMOVE = CURAINMOVE_X
!!$      IF (ASSOCIATED(CUSNOWMOVE)) CUSNOWMOVE = CUSNOWMOVE_X

      if (associated(dpdtmst)) dpdtmst = MAPL_GRAV * ( &
           (PFI_CN_X(:,:,0:lm-1)-PFI_CN_X(:,:,1:lm)) + &
           (PFI_AN_X(:,:,0:lm-1)-PFI_AN_X(:,:,1:lm)) + &
           (PFI_LS_X(:,:,0:lm-1)-PFI_LS_X(:,:,1:lm)) + &
           (PFI_SC_X(:,:,0:lm-1)-PFI_SC_X(:,:,1:lm)) + &
           (PFL_CN_X(:,:,0:lm-1)-PFL_CN_X(:,:,1:lm)) + &
           (PFL_AN_X(:,:,0:lm-1)-PFL_AN_X(:,:,1:lm)) + &
           (PFL_SC_X(:,:,0:lm-1)-PFL_SC_X(:,:,1:lm)) + &
           (PFL_LS_X(:,:,0:lm-1)-PFL_LS_X(:,:,1:lm)) )

      IF (ASSOCIATED(DLPDF)) DLPDF = DLPDF_X
      IF (ASSOCIATED(DIPDF)) DIPDF = DIPDF_X
      IF (ASSOCIATED(DLFIX)) DLFIX = DLFIX_X
      IF (ASSOCIATED(DIFIX)) DIFIX = DIFIX_X

      IF (ASSOCIATED(RHX))       RHX = RHX_X
      IF (ASSOCIATED(AUT))       AUT = AUT_X
      IF (ASSOCIATED(EVAPC))   EVAPC = EVAPC_X
      IF (ASSOCIATED(SDM))       SDM = SDM_X
      IF (ASSOCIATED(SUBLC))   SUBLC = SUBLC_X
      IF (ASSOCIATED(FRZ_TT)) FRZ_TT = FRZ_TT_X
      IF (ASSOCIATED(FRZ_PP)) FRZ_PP = FRZ_PP_X
      IF (ASSOCIATED(DCNVL))   DCNVL = DCNVL_X
      IF (ASSOCIATED(DCNVI))   DCNVI = DCNVI_X

      IF (ASSOCIATED(DNDCNV))   DNDCNV = DNDCNV_X
      IF (ASSOCIATED(DNCCNV))   DNCCNV = DNCCNV_X



      IF (ASSOCIATED(ALPHT)) ALPHT = ALPHT_X
      IF (ASSOCIATED(ALPH1)) ALPH1 = ALPH1_X
      IF (ASSOCIATED(ALPH2)) ALPH2 = ALPH2_X

      IF (ASSOCIATED(CFPDF)) CFPDF = CFPDF_X
      IF (ASSOCIATED(RHCLR)) RHCLR = RHCLR_X
      IF (ASSOCIATED(DQRL))   DQRL = DQRL_X

      IF (ASSOCIATED(VFALLICE_AN)) VFALLICE_AN = VFALLICE_AN_X
      IF (ASSOCIATED(VFALLICE_LS)) VFALLICE_LS = VFALLICE_LS_X
      IF (ASSOCIATED(VFALLWAT_AN)) VFALLWAT_AN = VFALLWAT_AN_X
      IF (ASSOCIATED(VFALLWAT_LS)) VFALLWAT_LS = VFALLWAT_LS_X

      IF (ASSOCIATED(VFALLRN_AN)) VFALLRN_AN = VFALLRN_AN_X
      IF (ASSOCIATED(VFALLRN_LS)) VFALLRN_LS = VFALLRN_LS_X
      IF (ASSOCIATED(VFALLRN_CN)) VFALLRN_CN = VFALLRN_CN_X
      IF (ASSOCIATED(VFALLRN_SC)) VFALLRN_SC = VFALLRN_SC_X
      IF (ASSOCIATED(VFALLSN_AN)) VFALLSN_AN = VFALLSN_AN_X
      IF (ASSOCIATED(VFALLSN_LS)) VFALLSN_LS = VFALLSN_LS_X
      IF (ASSOCIATED(VFALLSN_CN)) VFALLSN_CN = VFALLSN_CN_X
      IF (ASSOCIATED(VFALLSN_SC)) VFALLSN_SC = VFALLSN_SC_X
      

      call MAPL_TimerOff(STATE,"-CLOUD")

      call MAPL_TimerOn (STATE,"-MISC3")

      RAD_QV   = max( Q1 , 0. )

      
      IF ( INT(CLDPARAMS%DISABLE_RAD)==1 ) THEN
         RAD_QL     = 0.
         RAD_QI     = 0.
         RAD_CF     = 0.
         RAD_QR     = 0.
         RAD_QS     = 0.  
         CLDREFFL    = 10.e-6
         CLDREFFI    = 25.0e-6
      ENDIF


      CN_PRC2 = CN_PRC2 + RASPRCP
      CN_PRC2 = max(CN_PRC2 , 0.)
      LS_PRC2 = max(LS_PRC2 , 0.)
      AN_PRC2 = max(AN_PRC2 , 0.)
      SC_PRC2 = max(SC_PRC2 , 0.)
      CN_SNR  = max(CN_SNR  , 0.)
      LS_SNR  = max(LS_SNR  , 0.)
      AN_SNR  = max(AN_SNR  , 0.)
      SC_SNR  = max(SC_SNR  , 0.)

      TEMP    = TH1*PK

      ! Clean up Relative Humidity where RH > 110%
      !---------------------------------------------
      if (CLEANUP_RH == 1)  then

         RHEXCESS = 1.1
         

         !-----If using 2-moment microphysics, allow for supersaturation w.r.t ice ---
         if(adjustl(CLDMICRO)=="2MOMENT") then

            QSS=GEOS_QsatICE (TH1*PK, PLO*100.0)
            where (CFICE .lt. 0.99 .and. QSS .gt. 1.0e-20)            
               CFX = max(( Q1 - QSS*CFICE), 0.0)/(1.0-CFICE)     
               SAT_RAT = min(CFX/QSS, 2.0) ! this is just diagnostic   
            elsewhere
               SAT_RAT = 1.0
            end where

           if (associated(RHICE    )) RHICE  =  Q1/QSS
              where (TEMP>273.15)
                 RHICE=0.0
              end where
          
      
            QSS  = GEOS_QsatLQU(         &
                 TEMP   , &
                 PLO*100.0 , DQ=DQSDT     ) !clean up only with respect to liquid water  DONIF
          if (associated(RHLIQ    ))   RHLIQ     =  Q1/QSS  !DONIF
                 
         else
             DQSDT = GEOS_DQSAT (TEMP , PLO, QSAT=QSS ) 
                 
         end if
         !------

         where ( Q1 > RHEXCESS*QSS )
            DQS = (Q1 - RHEXCESS*QSS)/( 1.0 + RHEXCESS*DQSDT*MAPL_ALHL/MAPL_CP )
         elsewhere
            DQS = 0.0
         endwhere


         Q1      =  Q1    - DQS
         TEMP    =  TEMP  + (MAPL_ALHL/MAPL_CP)*DQS
         TH1     =  TEMP  / PK

         if(associated(dpdtmst)) dpdtmst = dpdtmst - (dqs*dp)/DT_MOIST


         ER_PRC2 =  SUM( DQS * MASS, 3)/DT_MOIST
         LS_PRC2 =  LS_PRC2 + ER_PRC2
      ELSE
         ER_PRC2 =  0.0
      ENDIF

      TPREC = CN_PRC2 + LS_PRC2 + AN_PRC2 + SC_PRC2 + &
           CN_SNR  + LS_SNR  + AN_SNR + SC_SNR


      ! Clean up any negative specific humidity
      !-----------------------------------------

      if(adjustl(CLDMICRO)=="1MOMENT") then
         call FILLQ2ZERO( Q1, MASS, FILLQ )
      else
         call FILLQ2ZERO2( Q1, MASS, FILLQ ) !Slightly different formulation
      end if

      ! Outputs for Spec
      !----------------------------------------------


!      if(associated(TVQX2)) TVQX2    = SUM( ( Q1 +  QLLS + QLCN + QILS + QICN )*MASS , 3 )  + TPREC*DT_MOIST
	 
      if (associated(ACR_TOT   ))   ACR_TOT    = ACLL_CN_X + ACIL_CN_X + ACLL_AN_X + & 
                                    ACIL_AN_X + ACLL_LS_X + ACIL_LS_X + ACLL_SC_X + ACIL_SC_X
      if (associated(REVSU_CN  ))   REVSU_CN   = REV_CN_X  + RSU_CN_X 
      if (associated(REVSU_LSAN))   REVSU_LSAN = REV_AN_X  + RSU_AN_X  + REV_LS_X  + RSU_LS_X 
      if (associated(PFL_LSAN  ))   then
         PFL_LSAN         = PFL_AN_X  + PFL_LS_X  
         PFL_LSAN(:,:,LM) = PFL_LSAN(:,:,LM) + ER_PRC2
      endif
      if (associated(PFI_LSAN  ))   PFI_LSAN   = PFI_AN_X  + PFI_LS_X  

      if (associated(PGENTOT   ))   PGENTOT    =  SUM( ( DQRL_X    + DQRC          &  
           + ACLL_CN_X + ACIL_CN_X       & 
           + ACLL_AN_X + ACIL_AN_X       & 
           + ACLL_SC_X + ACIL_SC_X       & 
           + ACLL_LS_X + ACIL_LS_X ) * MASS , 3 ) + ER_PRC2

      if (associated(PREVTOT   ))   PREVTOT    =  SUM( ( REV_AN_X  + RSU_AN_X        & 
           + REV_SC_X  + RSU_SC_X        & 
           + REV_LS_X  + RSU_LS_X ) * MASS , 3 )

      if(associated(cnv_mfc)) then
         if(associated(cnv_freq)) then
            cnv_freq(:,:) = 0.0
            do j=1,jm
               do i=1,im
                  do l=lm,1,-1
                     if(cnv_mfc(i,j,l)/=0.0) then
                        cnv_freq(i,j) = 1.0
                        exit
                     endif
                  enddo
               enddo
            enddo
         endif

         if(associated(cnv_basep)) then
            cnv_basep(:,:) = MAPL_UNDEF
            do j=1,jm
               do i=1,im
                  do l=lm,1,-1
                     if(cnv_mfc(i,j,l)/=0.0) then
                        cnv_basep(i,j) = ple(i,j,l-1)
                        exit
                     endif
                  enddo
               enddo
            enddo
         endif

         if(associated(cnv_topp)) then
            cnv_topp(:,:) = MAPL_UNDEF
            do j=1,jm
               do i=1,im
                  do l=1,lm
                     if(cnv_mfc(i,j,l)/=0.0) then
                        cnv_topp(i,j) = ple(i,j,l)
                        exit
                     endif
                  enddo
               enddo
            enddo
         endif
      endif

      if (associated(LS_ARF ))   LS_ARF  = LS_ARFX
      if (associated(AN_ARF ))   AN_ARF  = AN_ARFX
      if (associated(CN_ARF ))   CN_ARF  = CN_ARFX
      if (associated(SC_ARF ))   SC_ARF  = SC_ARFX
      if (associated(QX1    ))   QX1     = Q1
      if (associated(XQLLS  ))   XQLLS   = QLLS
      if (associated(XQILS  ))   XQILS   = QILS
      if (associated(XCLLS  ))   XCLLS   = CLLS
      if (associated(XQLCN  ))   XQLCN   = QLCN
      if (associated(XQICN  ))   XQICN   = QICN
      if (associated(XCLCN  ))   XCLCN   = CLCN
      if (associated(QITOT  ))   QITOT   = QICN + QILS + QSNOW + QGRAUPEL
      if (associated(QLTOT  ))   QLTOT   = QLCN + QLLS + QRAIN
      if (associated(QCTOT  ))   QCTOT   = QLCN + QLLS + QICN + QILS + QRAIN + QSNOW + QGRAUPEL
      if (associated(TVQ1   ))   TVQ1    = SUM( ( Q1 +  QLLS + QLCN + QILS + QICN + QRAIN + QSNOW + QGRAUPEL )*MASS , 3 ) & 
           +  TPREC*DT_MOIST
      if (associated(TVE1   ))   TVE1    = SUM( (  MAPL_CP*TEMP + MAPL_ALHL*Q1             & 
           -  MAPL_ALHF*(QILS+QICN) )*MASS , 3 )        &
           -  MAPL_ALHF*( CN_SNR  + LS_SNR  + AN_SNR + SC_SNR )*DT_MOIST

      if (associated(DCPTE  ))   DCPTE   = (  SUM(  MAPL_CP*TEMP *MASS , 3) - DCPTE )/DT_MOIST 
      if (associated(CWP    ))   CWP     = SUM( ( QLCN+QLLS+QICN+QILS+QRAIN+QSNOW+QGRAUPEL )*MASS , 3 )
      if (associated(LWP    ))   LWP     = SUM( ( QLCN+QLLS+QRAIN ) *MASS , 3 )
      if (associated(IWP    ))   IWP     = SUM( ( QICN+QILS+QSNOW+QGRAUPEL ) *MASS , 3 )
      if (associated(CCWP   ))   CCWP    = SUM(   CNV_QC *MASS , 3 )
      if (associated(TPW    ))   TPW     = SUM(   Q1         *MASS , 3 )
      if (associated(RH2    ))   RH2     = max(MIN( Q1/GEOS_QSAT (TH1*PK, PLO) , 1.02 ),0.0)
      if (associated(SLX    ))   SLX     = GZLO + MAPL_CP*TEMP - MAPL_ALHL*(QLLS+QLCN) - MAPL_ALHS*(QILS+QICN)
      if (associated(QTX    ))   QTX     = Q1 + QLCN + QLLS + QICN + QILS
      if (associated(TX     ))   TX      = TEMP

      if(adjustl(CLDMICRO)=="2MOMENT") then
         if (associated(CCNCOLUMN    ))   CCNCOLUMN      = SUM(CCN1*MASS/(100.*PLO*r_air/TEMP) , 3)
         if (associated(NDCOLUMN    ))    NDCOLUMN      =  SUM(NCPL_VOL*MASS/(100.*PLO*r_air/TEMP) , 3)
         if (associated(NCCOLUMN    ))    NCCOLUMN      =SUM(NCPI_VOL*MASS/(100.*PLO*r_air/TEMP) , 3)
      end if

      if (associated(PRECU  ))   PRECU   = CN_PRC2 + SC_PRC2
      if (associated(PRELS  ))   PRELS   = LS_PRC2 + AN_PRC2
      if (associated(TT_PRCP))   TT_PRCP = TPREC

   !  if(adjustl(CLDMICRO)=="GFDL") then
   !   ! Separate PRCP_ICE from LS_SNR which is the sum of PRCP_ICE + PRCP_SNOW + PRCP_GRAUPEL
   !    if (associated(SNR    ))   SNR     = PRCP_SNOW + PRCP_GRAUPEL + AN_SNR + CN_SNR + SC_SNR
   !    if (associated(ICE    ))   ICE     = PRCP_ICE
   !    if (associated(FRZR   ))   FRZR    = 0.0
   !  else
   !   ! Other microphysics for now have just snow/rain unless diagnosed later...
        if (associated(SNR    ))   SNR     = LS_SNR  + AN_SNR + CN_SNR + SC_SNR
        if (associated(ICE    ))   ICE     = 0.0
        if (associated(FRZR   ))   FRZR    = 0.0
   !  endif

      if (associated(HOURNORAIN))   then
         call ESMF_ClockGet(CLOCK, currTime=CurrentTime, rc=STATUS)
         call ESMF_TimeGet (CurrentTime, YY=YEAR, MM=MONTH, DD=DAY, H=HR, M=MN, S=SE, RC=STATUS )
         if(SE==0 .and. MN==0) then
            HOURNORAIN  = 0.
         else
            where(TPREC.EQ.0.)
               HOURNORAIN  = HOURNORAIN + DT_MOIST
            endwhere
         endif
      endif


      !--------------------------------------------------------------
      !  add ShallowCu contribution to total/detraining mass flux exports
      !--------------------------------------------------------------
      CNV_MFC = CNV_MFC + UMF_SC
      CNV_MFD = CNV_MFD + MFD_SC
      !--------------------------------------------------------------

! For 2 moment, move some LS precip/flux into the CN precip/flux category for use by chemistry
! --------------------------------------------------------------------------------------------
      if(adjustl(CLDMICRO)=="2MOMENT") then

      call MAPL_GetPointer(EXPORT, CU2DRAINMOVE,'CU2DRAINMOVE', RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CU2DSNOWMOVE,'CU2DSNOWMOVE', RC=STATUS); VERIFY_(STATUS)

      if(associated(CU2DRAINMOVE)) cu2drainmove = cn_prc2
      if(associated(CU2DSNOWMOVE)) cu2dsnowmove = cn_snr

       CN_PRC2 = CN_PRC2 + LS_PRC2*cnv_fraction
       LS_PRC2 = LS_PRC2 - LS_PRC2*cnv_fraction

       CN_SNR = CN_SNR + LS_SNR*cnv_fraction
       LS_SNR = LS_SNR - LS_SNR*cnv_fraction

      if(associated(CU2DRAINMOVE)) cu2drainmove = cn_prc2 - cu2drainmove
      if(associated(CU2DSNOWMOVE)) cu2dsnowmove = cn_snr - cu2dsnowmove


      call MAPL_GetPointer(EXPORT, PFLCNMOVE,'PFLCNMOVE', RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, PFICNMOVE,'PFICNMOVE', RC=STATUS); VERIFY_(STATUS)

      if(associated(PFLCNMOVE)) pflcnmove = pfl_cn
      if(associated(PFICNMOVE)) pficnmove = pfi_cn

       do l=1,lm

       pfl_cn  (:,:,L) = pfl_cn  (:,:,L) + pfl_ls(:,:,L)*cnv_fraction
       pfl_lsan(:,:,L) = pfl_lsan(:,:,L) - pfl_ls(:,:,L)*cnv_fraction
       pfl_ls  (:,:,L) = pfl_ls  (:,:,L) - pfl_ls(:,:,L)*cnv_fraction

       pfi_cn  (:,:,L) = pfi_cn  (:,:,L) + pfi_ls(:,:,L)*cnv_fraction
       pfi_lsan(:,:,L) = pfi_lsan(:,:,L) - pfi_ls(:,:,L)*cnv_fraction
       pfi_ls  (:,:,L) = pfi_ls  (:,:,L) - pfi_ls(:,:,L)*cnv_fraction

       enddo

      if(associated(PFLCNMOVE)) pflcnmove = pfl_cn - pflcnmove
      if(associated(PFICNMOVE)) pficnmove = pfi_cn - pficnmove

      endif
! --------------------------------------------------------------------------------------------

      DIAGNOSE_PTYPE: if (LDIAGNOSE_PRECIP_TYPE) then
         PTYPE(:,:) = 0 ! default PTYPE to rain
         ! Surface Precip Type diagnostic
         !
         !   PTYPE = 0  => Rain
         !   PTYPE = 1  => Freezing Rain
         !   PTYPE = 2  => Ice Pellets (sleet)
         !   PTYPE = 3  => Rain mixed with Snow
         !   PTYPE = 4  => Snow
         !
         ! Based on Pierre Bourgouin, 1999, "A Method to Determine Precipitation Types"
         !       in WEATHER AND FORECASTING Vol 15 pp 583-592
         do J=1,JM 
         do I=1,IM
! WMP: 2019-04-01 For now keep all frozen precip from microphysics frozen
         if (SNR(I,J) > 0.0) then ! Diagnostic PTYPE when GEOS thinks there is frozen precip...
           PTYPE(I,J) = 4 ! Start as snow
           PA2 = -999
           ! Sweep down the column from ~300mb looking for freezing/melting layers
           KTOP = max(1,count(PREF < 30000.))
           do while (KTOP < LM)
            NA = 0.0
            PA = 0.0
            TH_TOP  = TH1(I,J,KTOP) ! Potential Temp at top of layer
            TL_MEAN = 0.0 ! Layer mean absolute temperature
            Z_LAYER = 0.0 ! Layer thickness
            if (TEMP(I,J,KTOP) > MAPL_TICE) then
               do L=KTOP,LM-1
                  KTOP = L
                  if (TEMP(I,J,L) > MAPL_TICE) then
                      ZTHICK = TH1(I,J,L) * (PKE(I,J,L) - PKE(I,J,L-1)) * MAPL_CP/MAPL_GRAV
                      TL_MEAN = TL_MEAN + TEMP(I,J,L)*ZTHICK 
                      Z_LAYER = Z_LAYER + ZTHICK
                   else
                      if (Z_LAYER > 0) then
                         TL_MEAN = TL_MEAN/Z_LAYER
                         TH_BOT = TH1(I,J,L)
                        ! Determine depth of the warm layer [Positive Area (PA)]
                         PA = MAPL_CP*TL_MEAN*log( TH_TOP/TH_BOT )
                      endif
                      EXIT
                   endif
               enddo
            else
               do L=KTOP,LM-1
                  KTOP = L
                  if (TEMP(I,J,L) <= MAPL_TICE) then
                      ZTHICK = TH1(I,J,L) * (PKE(I,J,L) - PKE(I,J,L-1)) * MAPL_CP/MAPL_GRAV
                      TL_MEAN = TL_MEAN + TH1(I,J,L)*ZTHICK
                      Z_LAYER = Z_LAYER + ZTHICK
                  else
                     if (Z_LAYER > 0) then
                        TL_MEAN = TL_MEAN/Z_LAYER
                        TH_BOT = TH1(I,J,L)
                       ! Determine depth of the freezing layer [Negative Area (NA)]
                        NA = MAPL_CP*TL_MEAN*log( TH_TOP/TH_BOT )
                     endif
                     EXIT
                  endif
               enddo
            endif
            if (KTOP == LM-1) then
               TH_BOT = (TEMP(I,J,LM) + (MAPL_GRAV/MAPL_CP)*ZLO(I,J,LM)) / PKE(I,J,LM)
               ZTHICK = TH1(I,J,LM) * (PKE(I,J,LM) - PKE(I,J,LM-1)) * MAPL_CP/MAPL_GRAV
               TL_MEAN = TL_MEAN + TEMP(I,J,LM)*ZTHICK
               Z_LAYER = Z_LAYER + ZTHICK
               if (TEMP(I,J,LM) > MAPL_TICE) then
                  if (Z_LAYER > 0) then
                     TL_MEAN = TL_MEAN/Z_LAYER
                     ! Determine depth of the warm layer [Positive Area (PA)]
                     PA = MAPL_CP*TL_MEAN*log( TH_TOP/TH_BOT )
                  endif
               else
                  if (Z_LAYER > 0) then
                     TL_MEAN = TL_MEAN/Z_LAYER
                     ! Determine depth of the freezing layer [Negative Area (NA)]
                     NA = MAPL_CP*TL_MEAN*log( TH_TOP/TH_BOT )
                  endif
               endif
               if (SNR(I,J) > 0.0) then ! Diagnostic PTYPE when GEOS thinks there is frozen precip...
                 if (PTYPE(I,J) == 4) then ! No Warm layer found above the surface yet
                   if (PA <  5.6) PTYPE(I,J) = 4 ! SNOW
                   if (PA >= 5.6) PTYPE(I,J) = 3 ! Mix of Snow Ice and Rain
! WMP: 2019-04-01 For now keep all frozen precip from microphysics frozen
! WMP: 2019-04-01  if (PA > 13.2) PTYPE(I,J) = 0 ! RAIN
                 else
                   if ( NA <  (46.0 + 0.66*PA2) ) PTYPE(I,J) = 1   ! Freezing Rain
                   if ( NA >= (46.0 + 0.66*PA2) ) PTYPE(I,J) = 1.5 ! Freezing Rain & Ice Pellets (sleet)
                   if ( NA >  (66.0 + 0.66*PA2) ) PTYPE(I,J) = 2   ! Ice Pellets (sleet)
                 endif
               else
                                 PTYPE(I,J) = 0 ! Rain
                 if (NA > 50.0 ) PTYPE(I,J) = 1 ! Freezing Rain
                 if (NA > 200.0) PTYPE(I,J) = 2 ! Ice Pellets (sleet)
               endif
               KTOP = LM
            else
               if (PA > 2.0) then ! Found a warm layer...
                  if (PA <  5.6) PTYPE(I,J) = 4 ! SNOW
                  if (PA >= 5.6) PTYPE(I,J) = 3 ! Mix of Snow and Rain
! WMP: 2019-04-01 For now keep all frozen precip from microphysics frozen
! WMP: 2019-04-01 if (PA > 13.2) PTYPE(I,J) = 0 ! RAIN
                  PA2 = PA
               else ! Found a freezing layer
                  PA2 = 0
                  if ( (PTYPE(I,J) <= 3) ) then
                     if (NA > 50.0 ) PTYPE(I,J) = 1 ! Freezing Rain
                     if (NA > 200.0) PTYPE(I,J) = 2 ! Ice Pellets (sleet)
                  endif      
               endif
            endif
           enddo
         endif
         enddo
         enddo

         ! Zero out snowfall exports when using PTYPE diagnostics if not SNOW/MIX
         if (associated(SNR) .AND. associated(PTYPE)) then
         do J=1,JM
            do I=1,IM
               if (PTYPE(I,J) <= 2) SNR(I,J) = 0.0
            enddo
         enddo
         endif

         ! Fill ICE fall diagnostic
         if (associated(ICE) .AND. associated(PTYPE)) then
            ICE = 0.0
            do J=1,JM
               do I=1,IM
                  if (PTYPE(I,J) == 2) ICE(I,J) = TPREC(I,J) 
               enddo
            enddo
         else
            if (associated(ICE))   ICE = 0.0
         endif

         ! Fill FRZR fall diagnostic
         if (associated(FRZR) .AND. associated(PTYPE)) then
            FRZR = 0.0
            do J=1,JM
               do I=1,IM
                  if (PTYPE(I,J) == 1  ) FRZR(I,J) = TPREC(I,J)
                  if (PTYPE(I,J) == 1.5) FRZR(I,J) = TPREC(I,J)
               enddo
            enddo
         else
            if (associated(FRZR ))   FRZR = 0.0
         endif

! WMP: 2019-04-01
! For now keep all frozen precip from microphysics frozen
!        ! Move convective frozen precip at surface to liquid based on PTYPE
!        if (associated(PRECU) .AND. associated(PTYPE)) then
!         do J=1,JM
!           do I=1,IM
!              if (PTYPE(I,J) <= 2) then
!                 CN_PRC2(I,J) = CN_PRC2(I,J) + CN_SNR(I,J)
!                 SC_PRC2(I,J) = SC_PRC2(I,J) + SC_SNR(I,J)
!                 CN_SNR(I,J) = 0.0
!                 SC_SNR(I,J) = 0.0
!              endif
!           enddo
!         enddo
!        endif
!
!        ! Move large-scale frozen precip at surface to liquid based on PTYPE
!        if (associated(PRECU) .AND. associated(PTYPE)) then
!         do J=1,JM
!           do I=1,IM
!              if (PTYPE(I,J) <= 2) then
!                 LS_PRC2(I,J) = LS_PRC2(I,J) + LS_SNR(I,J)
!                 AN_PRC2(I,J) = AN_PRC2(I,J) + AN_SNR(I,J)
!                 LS_SNR(I,J) = 0.0
!                 AN_SNR(I,J) = 0.0
!              endif
!           enddo
!         enddo
!        endif
! WMP: 2019-04-01

      endif DIAGNOSE_PTYPE

      if (associated(CN_PRCP))   CN_PRCP = CN_PRC2 + CN_SNR
      if (associated(SC_PRCP))   SC_PRCP = SC_PRC2 + SC_SNR
      if (associated(AN_PRCP))   AN_PRCP = AN_PRC2 + AN_SNR
      if (associated(LS_PRCP))   LS_PRCP = LS_PRC2 + LS_SNR
      if (associated(ER_PRCP))   ER_PRCP = ER_PRC2 
      if (associated(FILLNQV))   FILLNQV = FILLQ / DT_MOIST 
      if (associated(DQDT   ))   DQDT    = (Q1  - Q )/DT_MOIST
      if (associated(DTDT_MOIST)) DTDT_MOIST = (TEMP - DTDT_MOIST)/DT_MOIST
      !!if(associated(QSATI   ))   DQS      = GEOS_DQSAT(TEMP, PLO, qsat=QSATi, OVER_ICE=.TRUE. )
      !!if(associated(QSATl   ))   DQS      = GEOS_DQSAT(TEMP, PLO, qsat=QSATl, OVER_LIQUID=.TRUE. )
      if (associated(TI     ))   TI      = (TH1 - TH)*(PLE(:,:,1:LM)-PLE(:,:,0:LM-1))/DT_MOIST

      if (associated(UI     ))   UI      = (U1  - U )/DT_MOIST
      if (associated(VI     ))   VI      = (V1  - V )/DT_MOIST
      if(associated(DQLDT   ))   then 
         DQLDT   = (QLLS + QLCN - DQLDT)
         where( ABS(DQLDT) < 1.e-20 )
            DQLDT=0.0
         endwhere
         DQLDT   = DQLDT / DT_MOIST
      endif
      if(associated(DQIDT   ))   then 
         DQIDT   = (QILS + QICN - DQIDT)
         where( ABS(DQIDT) < 1.e-20 )
            DQIDT=0.0
         endwhere
         DQIDT   = DQIDT / DT_MOIST
      endif
         
      if (associated(CLDBASEHGT)) then
         CLDBASEx = MAPL_UNDEF
         do i = 1,IM
           do j = 1,JM
             do k =  LM, 1, -1
               if (ZLE(i,j,k).gt.20000.) exit
               if ( ( RAD_CF(i,j,k) .ge. 1e-2 ) ) then
                 CLDBASEx(i,j)  = ZLE(i,j,k)
                 exit
               end if
             end do
           end do
         end do
         CLDBASEHGT = CLDBASEx
      end if
         
      CFX =100.*PLO*r_air/TEMP

      if(associated(DNDDT   ))   then 
         DNDDT   = (NCPL - DNDDT)
         where( ABS(DNDDT) < 1.e-20 )
            DNDDT=0.0
         endwhere
         DNDDT   = CFX*DNDDT / DT_MOIST
      endif
    if(associated(DNCDT   ))   then 
         DNCDT   = (NCPI - DNCDT)
         where( ABS(DNCDT) < 1.e-20 )
            DNCDT=0.0
         endwhere
         DNCDT   = CFX*DNCDT / DT_MOIST
      endif
     
      temp = (0.5/DT_MOIST)*((V1**2+U1**2)  - (V**2+U**2))*MASS
      !temp = (1.0/DT_MOIST)*( U *( U1 - U )  + V * ( V1 - V ) ) *MASS

      IKEX  = SUM( TEMP , 3 )     
      IKEX2 = MAX(  SUM( KEX * MASS , 3 ) ,  1.0e-6 ) ! floor at 1e-6 W m-2   

      if (associated(KEMST   ))  KEMST    = IKEX
      if (associated(KEMST2  ))  KEMST2   = IKEX2

      !scaled 3D kinetic energy dissipation
      if (associated(KEDISS  ))  then 
         do L=1,LM
            KEDISS(:,:,L)   = (IKEX/IKEX2) * KEX(:,:,L)
         enddo
      end if


      if(associated(DTDTFRIC)) then 
         do L=1,LM
            DTDTFRIC(:,:,l) = -(1./MAPL_CP)*(IKEX/IKEX2) * KEX(:,:,L) * (PLE(:,:,L)-PLE(:,:,L-1))
         end do
      end if

      ! Compute aerosol tendncies due to moist
      ! --------------------------------------
      ! if ( associated(ddudt) ) then
      !    do k = 1, km
      !       if it is dust: (based on qnames)
      !          ddudt = sum_k q(k) * delp(k) - ddudt
      !    else if is sea salt ...
      !
      ! end do


      if (associated(CNVRNZ     ))   CNVRNZ    =  SUM( DQRC  * MASS , 3 )
      if (associated(PDFLZ      ))   PDFLZ     =  SUM( DLPDF_X * MASS , 3 )
      if (associated(CNVLZ      ))   CNVLZ     =  SUM( (DCNVL_X + QLENT_SC + QLSUB_SC) * MASS , 3 )

      if (associated(PDFIZ      ))   PDFIZ     =  SUM( DIPDF_X * MASS , 3 )
      if (associated(CNVIZ      ))   CNVIZ     =  SUM( DCNVI_X * MASS , 3 )

      if (associated(EVPCZ      ))   EVPCZ     =  SUM( ( EVAPC_X  + DLFIX_X ) * MASS , 3 )
      if (associated(EVPPZ      ))   EVPPZ     =  SUM( ( REV_AN_X + REV_LS_X + REV_SC_X ) * MASS , 3 )

      if (associated(SUBCZ      ))   SUBCZ     =  SUM( ( SUBLC_X  + DIFIX_X ) * MASS , 3 )
      if (associated(SUBPZ      ))   SUBPZ     =  SUM( ( RSU_AN_X + RSU_LS_X + RSU_SC_X ) * MASS , 3 )

      if (associated(FRZCZ      ))   FRZCZ     =  SUM( FRZ_TT_X * MASS , 3 )
      if (associated(FRZPZ      ))   FRZPZ     =  SUM( FRZ_PP_X * MASS , 3 )
      if (associated(COLLIZ     ))   COLLIZ    =  SUM( ( ACIL_CN_X + ACIL_AN_X + ACIL_LS_X + ACIL_SC_X ) * MASS , 3 )

      if (associated(COLLLZ     ))   COLLLZ    =  SUM( ( ACLL_CN_X + ACLL_AN_X + ACLL_LS_X + ACLL_SC_X ) * MASS , 3 )
      if (associated(AUTZ       ))   AUTZ      =  SUM( AUT_X * MASS , 3 )

      if (associated(SDMZ       ))   SDMZ      =  SUM( SDM_X * MASS , 3 )


      ! Replace the modified humidity
      !------------------------------

      Q = Q1
      !AMM to sync up T and Q also update to the modified TH
      if(associated(THMOIST  )) THMOIST    = TH1
      ! edge geopotential - array defined from 0 to LM (ie., bottom edge), pke defined that way too
      geopenew(:,:,LM) = 0.
      do k=lm-1,0,-1
         geopenew(:,:,k) = geopenew(:,:,k+1) + mapl_cp*(th1(:,:,k+1)*(1.+mapl_vireps*q1(:,:,k+1)))*( pke(:,:,k+1)-pke(:,:,k) )
      enddo
      !!   if(associated(SMOIST)) SMOIST = mapl_cp*th1*pk + (0.5*( geopenew(:,:,0:lm-1)+geopenew(:,:,1:lm) ))
      if(associated(SMOIST)) SMOIST = mapl_cp*th1*pk + gzlo


      ! Calculate flash rate following Murray et al. (2012), as used by GEOS-Chem 
      !-------------------------------------------------------------------------------------
      CALL MAPL_GetPointer( EXPORT, LFR_GCC, 'LFR_GCC', NotFoundOk=.TRUE., __RC__ )
      IF ( ASSOCIATED( LFR_GCC ) ) THEN
         CALL MAPL_TimerOn (STATE,"--FLASH", __RC__ )
         CALL Get_hemcoFlashrate ( STATE,       &
                                   IMPORT,      &
                                   IM, JM, LM,  &
                                   T,           &
                                   PLE,         &
                                   ZLE,         &
                                   CNV_MFC,     &
                                   AREA,        &
                                   TS,          &
                                   LFR_GCC,     &
                                          __RC__ )
         CALL MAPL_TimerOff(STATE,"--FLASH", __RC__ )
      ENDIF

      ! Deallocate temp space if necessary
      !-----------------------------------

      if(ALLOC_PTYPE) deallocate(PTYPE)

      if(associated(DDF_BYNC  )) deallocate( DDF_BYNCz )
      if(associated(DDF_MUPH  )) deallocate( DDF_MUPHz )
      if(associated(DDF_RH1   )) deallocate( DDF_RH1z )
      if(associated(DDF_RH2   )) deallocate( DDF_RH2z )
      if(associated(DDF_DQDT  )) deallocate( DDF_DQDTz )
      if(associated(DDF_DTDT  )) deallocate( DDF_DTDTz )
      if(associated(DDF_TC    )) deallocate( DDF_TCz  )
      if(associated(DDF_QVC   )) deallocate( DDF_QVCz )
      if(associated(DDF_MFC   )) deallocate( DDF_MFCz )

      if(ALLOC_BYNCY    ) deallocate(BYNCY    )
      if(ALLOC_CAPE     ) deallocate(CAPE     )
      if(ALLOC_INHB     ) deallocate(INHB     )
      if(ALLOC_CNV_DQLDT) deallocate(CNV_DQLDT)
      if(ALLOC_CNV_MF0  ) deallocate(CNV_MF0  )
      if(ALLOC_CNV_MFD  ) deallocate(CNV_MFD  )
      if(ALLOC_CNV_MFC  ) deallocate(CNV_MFC  )
      if(ALLOC_CNV_TOPP ) deallocate(CNV_TOPP )
      if(ALLOC_CNV_UPDF ) deallocate(CNV_UPDF )
      if(ALLOC_CNV_CVW  ) deallocate(CNV_CVW  )
      if(ALLOC_CNV_QC )   deallocate(CNV_QC   )
      if(ALLOC_RAD_CF   ) deallocate(RAD_CF   )
      if(ALLOC_RAD_QV   ) deallocate(RAD_QV   )
      if(ALLOC_RAD_QL   ) deallocate(RAD_QL   )
      if(ALLOC_RAD_QI   ) deallocate(RAD_QI   )
      if(ALLOC_RAD_QR   ) deallocate(RAD_QR   )
      if(ALLOC_RAD_QS   ) deallocate(RAD_QS   )
      if(ALLOC_RAD_QG   ) deallocate(RAD_QG   )

      if(ALLOC_CLDREFFL ) deallocate(CLDREFFL )
      if(ALLOC_CLDREFFI ) deallocate(CLDREFFI )
      if(ALLOC_CLDREFFR ) deallocate(CLDREFFR )
      if(ALLOC_CLDREFFS ) deallocate(CLDREFFS )
      if(ALLOC_CLDREFFG ) deallocate(CLDREFFG )
 
      if(ALLOC_CLDNCCN  ) deallocate(CLDNCCN  )

      if(ALLOC_DQRC )      deallocate ( DQRC  )

      if(ALLOC_ENTLAM )    deallocate ( ENTLAM )

      deallocate ( QRN )
      deallocate ( QSN )
      deallocate ( QPLS )

      if(ALLOC_CNV_FICE)  deallocate(CNV_FICE)
      if(ALLOC_CNV_NDROP)  deallocate(CNV_NDROP)
      if(ALLOC_CNV_NICE)  deallocate(CNV_NICE)     

      if(ALLOC_CFICE  )        deallocate(CFICE )
      if(ALLOC_CFLIQ  )        deallocate(CFLIQ )


      !!! Deallocate UW shallow vars

      if(ALLOC_UMF)    deallocate( UMF_SC )
      if(ALLOC_MFD)    deallocate( MFD_SC )
      if(ALLOC_DQVDT)  deallocate( DQVDT_SC )
      if(ALLOC_DQLDT)  deallocate( DQLDT_SC )
      if(ALLOC_DQIDT)  deallocate( DQIDT_SC )
      if(ALLOC_DTHDT)  deallocate( DTHDT_SC )
      if(ALLOC_DUDT)   deallocate( DUDT_SC )
      if(ALLOC_DVDT)   deallocate( DVDT_SC )
      if(ALLOC_DQRDT)  deallocate( DQRDT_SC )
      if(ALLOC_DQSDT)  deallocate( DQSDT_SC )
      if(ALLOC_CUFRC)  deallocate( CUFRC_SC )
      if(ALLOC_QCU)    deallocate( QCU_SC )
      if(ALLOC_QLU)    deallocate( QLU_SC )
      if(ALLOC_QIU)    deallocate( QIU_SC )
      if(ALLOC_CBMF)   deallocate( CBMF_SC )
      if(ALLOC_DQCDT)  deallocate( DQCDT_SC )
      if(ALLOC_CNT)    deallocate( CNT_SC )
      if(ALLOC_CNB)    deallocate( CNB_SC )
      if(ALLOC_CIN)    deallocate( CIN_SC )
      if(ALLOC_PLCL)   deallocate( PLCL_SC )
      if(ALLOC_PLFC)   deallocate( PLFC_SC )
      if(ALLOC_PINV)   deallocate( PINV_SC )
      if(ALLOC_PREL)   deallocate( PREL_SC )
      if(ALLOC_PBUP)   deallocate( PBUP_SC )
      if(ALLOC_WLCL)   deallocate( WLCL_SC )
      if(ALLOC_QTSRC)  deallocate( QTSRC_SC )
      if(ALLOC_THLSRC) deallocate( THLSRC_SC )
      if(ALLOC_THVLSRC)deallocate( THVLSRC_SC )
      if(ALLOC_TKEAVG) deallocate( TKEAVG_SC )
      if(ALLOC_CLDTOP) deallocate( CLDTOP_SC )

      if(ALLOC_WUP)    deallocate( WUP_SC )
      if(ALLOC_QTUP)   deallocate( QTUP_SC )
      if(ALLOC_THLUP)  deallocate( THLUP_SC )
      if(ALLOC_THVUP ) deallocate( THVUP_SC )
      if(ALLOC_UUP)    deallocate( UUP_SC )
      if(ALLOC_VUP)    deallocate( VUP_SC )
      if(ALLOC_ENTR)   deallocate( ENTR_SC )
      if(ALLOC_DETR)   deallocate( DETR_SC )
      if(ALLOC_XC)     deallocate( XC_SC )
      if(ALLOC_QLDET)  deallocate( QLDET_SC )
      if(ALLOC_QIDET)  deallocate( QIDET_SC )
      if(ALLOC_QLENT)  deallocate( QLENT_SC )
      if(ALLOC_QIENT)  deallocate( QIENT_SC )
      if(ALLOC_QLSUB)  deallocate( QLSUB_SC )
      if(ALLOC_QISUB)  deallocate( QISUB_SC )
      if(ALLOC_NDROP)  deallocate( SC_NDROP )
      if(ALLOC_NICE)   deallocate( SC_NICE )

      if(adjustl(CLDMICRO)=="GFDL") then
         if(ALLOC_PRCP_RAIN  )  deallocate(PRCP_RAIN  )
         if(ALLOC_PRCP_SNOW  )  deallocate(PRCP_SNOW  )
         if(ALLOC_PRCP_ICE  )  deallocate(PRCP_ICE  )
         if(ALLOC_PRCP_GRAUPEL  )  deallocate(PRCP_GRAUPEL  )
      endif

      if(adjustl(CLDMICRO)=="2MOMENT" .or. adjustl(CLDMICRO)=="GFDL") then
         if(ALLOC_DQVDT_micro  )  deallocate(DQVDT_micro )
         if(ALLOC_DQIDT_micro  )  deallocate(DQIDT_micro )
         if(ALLOC_DQLDT_micro  )  deallocate(DQLDT_micro )
         if(ALLOC_DQRDT_micro  )  deallocate(DQRDT_micro )
         if(ALLOC_DQSDT_micro  )  deallocate(DQSDT_micro )
         if(ALLOC_DQGDT_micro  )  deallocate(DQGDT_micro )
         if(ALLOC_DQADT_micro  )  deallocate(DQADT_micro )
         if(ALLOC_DUDT_micro   )  deallocate(DUDT_micro  )
         if(ALLOC_DVDT_micro   )  deallocate(DVDT_micro  )
         if(ALLOC_DTDT_micro   )  deallocate(DTDT_micro  )
         if(ALLOC_DTDT_macro   )  deallocate(DTDT_macro  )
         if(ALLOC_DQVDT_macro  )  deallocate(DQVDT_macro )
         if(ALLOC_DQLDT_macro  )  deallocate(DQLDT_macro )
         if(ALLOC_DQIDT_macro  )  deallocate(DQIDT_macro )
         if(ALLOC_DQADT_macro  )  deallocate(DQADT_macro )
         if(ALLOC_PFRZ         )  deallocate(PFRZ        )
         if(ALLOC_SC_ICE       )  deallocate(SC_ICE      )
         if(ALLOC_DT_RASP      )  deallocate(DT_RASP     )
      endif

      deallocate(FQAl, stat=STATUS)
      VERIFY_(STATUS)
      deallocate(FQAI, stat=STATUS)
      VERIFY_(STATUS)
      deallocate(FQA, stat=STATUS)
      VERIFY_(STATUS)

      if(adjustl(CLDMICRO)=="2MOMENT") then

         deallocate(QCNTOT, stat=STATUS)
         VERIFY_(STATUS)

         !Aer_cloud----------------------

         if(ALLOC_SMAXL  ) deallocate(SMAXL  )
         if(ALLOC_WSUB   ) deallocate(WSUB   )
         if(ALLOC_CCN01  ) deallocate(CCN01  )
         if(ALLOC_CCN04  ) deallocate(CCN04  )
         if(ALLOC_CCN1   ) deallocate(CCN1  )
         if(ALLOC_SMAXI  ) deallocate(SMAXI  )
         if(ALLOC_CDNC_NUC  ) deallocate(CDNC_NUC  )
         if(ALLOC_INC_NUC  )  deallocate(INC_NUC  )
         if(ALLOC_NCPL_VOL  ) deallocate(NCPL_VOL  )
         if(ALLOC_NCPI_VOL  ) deallocate(NCPI_VOL  )

         if(ALLOC_SO4  )      deallocate(SO4  )
         if(ALLOC_ORG  )      deallocate(ORG  )
         if(ALLOC_DUST  )     deallocate(DUST  )
         if(ALLOC_SEASALT )  deallocate(SEASALT  )
         if(ALLOC_BCARBON )  deallocate(BCARBON  )
         if(ALLOC_NHET_NUC  )  deallocate(NHET_NUC  )
         if(ALLOC_NLIM_NUC  )  deallocate(NLIM_NUC  )
         if(ALLOC_SAT_RAT  )  deallocate(SAT_RAT  )
         if(ALLOC_QSTOT )  deallocate(QSTOT)
         if(ALLOC_QRTOT )  deallocate(QRTOT)
         if(ALLOC_QPTOTLS )  deallocate(QPTOTLS)

         if(ALLOC_DTDT_moist  )   deallocate(DTDT_moist )            
         if(ALLOC_DTDTCN  )   deallocate(DTDTCN )            

         if(ALLOC_RL_MASK  )      deallocate(RL_MASK )
         if(ALLOC_RI_MASK  )      deallocate(RI_MASK )
         if(ALLOC_KAPPA  )        deallocate(KAPPA )
         if(ALLOC_RHICE  )        deallocate(RHICE)
         if(ALLOC_RHLIQ  )        deallocate(RHLIQ )


         if(ALLOC_NHET_IMM) deallocate(NHET_IMM)    
         if(ALLOC_NHET_DEP) deallocate(NHET_DEP)
         if(ALLOC_DUST_IMM) deallocate(DUST_IMM)
         if(ALLOC_DUST_DEP) deallocate(DUST_DEP)
      
         if(ALLOC_SCF)      deallocate(SCF)
         if(ALLOC_SCF_ALL)  deallocate(SCF_ALL)
         if(ALLOC_SIGW_GW)  deallocate(SIGW_GW)
         if(ALLOC_SIGW_CNV)  deallocate(SIGW_CNV)
         if(ALLOC_SIGW_TURB)  deallocate(SIGW_TURB)
         if(ALLOC_SIGW_RC)  deallocate(SIGW_RC)
         if(ALLOC_RHCmicro)  deallocate(RHCmicro)
         if(ALLOC_DNHET_IMM) deallocate(DNHET_IMM)
         if(ALLOC_BERG) deallocate(BERG)
         if(ALLOC_BERGS) deallocate(BERGS)
         if(ALLOC_MELT) deallocate(MELT)
         if(ALLOC_DNHET_CT) deallocate(DNHET_CT)
         if(ALLOC_QCRES) deallocate(QCRES)
         if(ALLOC_QIRES) deallocate(QIRES)
         if(ALLOC_AUTICE  )   deallocate  (AUTICE)       
         if(ALLOC_FRZPP_LS) deallocate(FRZPP_LS)
         if(ALLOC_SNOWMELT_LS) deallocate(SNOWMELT_LS)

         if(ALLOC_DNCNUC) deallocate(DNCNUC)
         if(ALLOC_DNCSUBL) deallocate(DNCSUBL)
         if(ALLOC_DNCHMSPLIT) deallocate(DNCHMSPLIT)    
         if(ALLOC_DNCAUTICE) deallocate(DNCAUTICE)
         if(ALLOC_DNCACRIS) deallocate(DNCACRIS)       
         if(ALLOC_DNDCCN) deallocate(DNDCCN)
         if(ALLOC_DNDACRLS) deallocate(DNDACRLS)
         if(ALLOC_DNDACRLR) deallocate(DNDACRLR)
         if(ALLOC_DNDEVAPC) deallocate(DNDEVAPC)
         if(ALLOC_DNDAUTLIQ) deallocate(DNDAUTLIQ)      
         if(ALLOC_DNDCNV ) deallocate (DNDCNV)
         if(ALLOC_DNCCNV ) deallocate (DNCCNV) 
         



         !if (associated(CLDREFFI_TOP)) deallocate (CLDREFFI_TOP) 
         !if (associated(CLDREFFL_TOP)) deallocate (CLDREFFL_TOP)
         !if (associated(NCPI_TOP)) deallocate (NCPI_TOP) 
         !if (associated(NCPL_TOP)) deallocate (NCPL_TOP) 
         !if (associated(QCVAR_EXP)) deallocate (QCVAR_EXP)  



      end if

      call MAPL_TimerOff(STATE,"-MISC3")



      !  All done
      !-----------

      RETURN_(ESMF_SUCCESS)

    end subroutine MOIST_DRIVER

!!!!!!!!-!-!-!!!!!!!!
  end subroutine RUN

!BOP

! !IROUTINE: AINC_UPDATE -- Update Q with analysis increment

! !INTERFACE:

  subroutine AINC_UPDATE ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

  type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
  type(ESMF_State),    intent(inout) :: IMPORT ! Import state
  type(ESMF_State),    intent(inout) :: EXPORT ! Export state
  type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
  integer, optional,   intent(  out) :: RC     ! Error code:

! ! DESCRIPTION: This RUN method simply updates Q with the analysis increment.

!EOP

! ErrLog Variables

  character(len=ESMF_MAXSTR)          :: IAm
  integer                             :: STATUS
  character(len=ESMF_MAXSTR)          :: COMP_NAME

! Local derived type aliases

  type (MAPL_MetaComp),     pointer   :: MAPL
  type (ESMF_State       )            :: INTERNAL

  integer                             :: IM, JM, LM
  real, pointer, dimension(:,:,:)     :: dqv
  real, pointer, dimension(:,:,:)     :: qv

  type(ESMF_Grid)                     :: grid

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

   Iam = "AINC_UPDATE"
   call ESMF_GridCompGet( GC, name=COMP_NAME, RC=STATUS )
   VERIFY_(STATUS)
   Iam = trim(COMP_NAME) // Iam

   if ( MAPL_AM_I_ROOT() ) then
       print *, 'Now running ',trim(Iam)
   endif

! Retrieve the pointer to the state
!----------------------------------

   call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
   VERIFY_(STATUS)

! Local aliases to the state, grid, and configuration
! ---------------------------------------------------

!  call MAPL_TimerOn(MAPL,"TOTAL")

! Get parameters from generic state.
!-----------------------------------

    call MAPL_Get( MAPL, IM=IM, JM=JM, LM=LM,    &
                   INTERNAL_ESMF_STATE=INTERNAL, &
                                       RC=STATUS )
    VERIFY_(STATUS)

    call ESMF_GridCompGet(GC, grid=grid, rc=status)
    VERIFY_(STATUS)

! **********************************************************************
! ****               Get Pointers to BKG Import Data                ****
! **********************************************************************
#if 0
    if ( MAPL_AM_I_ROOT() ) then
       call ESMF_StatePrint(IMPORT)
    end if
#endif

!   Get pointers to import variables
!   --------------------------------
    call MAPL_GetPointer(import,   dqv, 'QVAINC',  RC=STATUS)
    VERIFY_(STATUS)

!   Get pointers to internal variables
!   ----------------------------------
    call MAPL_GetPointer(internal,   qv, 'Q',  RC=STATUS)
    VERIFY_(STATUS)

    if (associated(qv) .and. associated(dqv) ) then
       qv = max(0.0,qv + dqv) ! this update needs to be revised in view of Takacs, Suarez & Todling (2016)
    endif


  end subroutine AINC_UPDATE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine FILLQ2ZERO( Q, MASS, FILLQ  )

    real, dimension(:,:,:),   intent(inout)  :: Q
    real, dimension(:,:,:),   intent(in)     :: MASS
    real, dimension(:,:),     intent(  out)  :: FILLQ
    integer                                  :: IM,JM,LM
    integer                                  :: I,J,K,L

    real                                     :: TPW, NEGTPW

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Fills in negative q values in a mass conserving way.
    ! Conservation of TPW was checked.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    IM = SIZE( Q, 1 )
    JM = SIZE( Q, 2 )
    LM = SIZE( Q, 3 )


    do j=1,JM
       do i=1,IM

          TPW = SUM( Q(i,j,:)*MASS(i,j,:) )

          NEGTPW = 0.
          do l=1,LM
             if ( Q(i,j,l) < 0.0 ) then 
                NEGTPW   = NEGTPW + ( Q(i,j,l)*MASS( i,j,l ) )
                Q(i,j,l) = 0.0
             endif
          enddo

          do l=1,LM
             if ( Q(i,j,l) >= 0.0 ) then 
                Q(i,j,l) = Q(i,j,l)*( 1.0+NEGTPW/(TPW-NEGTPW) )
             endif
          enddo

          FILLQ(i,j) = -NEGTPW

       end do
    end do



  end subroutine FILLQ2ZERO



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine FILLQ2ZERO2( Q, MASS, FILLQ  )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! New algorithm to fill the negative q values in a mass conserving way.
    ! Conservation of TPW was checked. Donifan Barahona
    ! Updated from FILLQ2ZERO, avoids the usage of scalars
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    real, dimension(:,:,:),   intent(inout)  :: Q
    real, dimension(:,:,:),   intent(in)     :: MASS
    real, dimension(:,:),     intent(  out)  :: FILLQ
    real, dimension(:,:), allocatable        :: TPW1, TPW2, TPWC
    integer                                  :: IM,JM,LM, l

    IM = SIZE( Q, 1 )
    JM = SIZE( Q, 2 )
    LM = SIZE( Q, 3 )  

    ALLOCATE(TPW1(IM, JM))
    ALLOCATE(TPW2(IM, JM))
    ALLOCATE(TPWC(IM, JM))

    TPW2 =0.0
    TPWC= 0.0
    TPW1 = SUM( Q*MASS, 3 ) 

    WHERE (Q < 0.0)  
       Q=0.0
    END WHERE

    TPW2 = SUM( Q*MASS, 3 )

    WHERE (TPW2 > 0.0)       
       TPWC=(TPW2-TPW1)/TPW2
    END WHERE

    do l=1,LM       
       Q(:, :, l)= Q(:, :, l)*(1.0-TPWC) !reduce Q proportionally to the increase in TPW
    end do

    FILLQ = TPW2-TPW1 

    DEALLOCATE(TPW1)
    DEALLOCATE(TPW2)
    DEALLOCATE(TPWC)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  end subroutine FILLQ2ZERO2



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine BUOYANCY( T, Q, QS, DQS, DZ, ZLO, BUOY, CAPE, INHB)


    ! !DESCRIPTION: Computes the buoyancy $ g \frac{T_c-T_e}{T_e} $ at each level
    !  for a parcel raised from the surface. $T_c$ is the virtual temperature of
    !  the parcel and $T_e$ is the virtual temperature of the environment.

    real, dimension(:,:,:),   intent(in)  :: T, Q, QS, DQS, DZ, ZLO
    real, dimension(:,:,:),   intent(out) :: BUOY
    real, dimension(:,:),     intent(out) :: CAPE, INHB

    integer :: L, LM

    LM = size(T,3)

    BUOY(:,:,LM) =  T(:,:,LM) + (MAPL_GRAV/MAPL_CP)*ZLO(:,:,LM) + (MAPL_ALHL/MAPL_CP)*Q(:,:,LM)

    do L=LM-1,1,-1
       BUOY(:,:,L) = BUOY(:,:,LM) - (T(:,:,L) + (MAPL_GRAV/MAPL_CP)*ZLO(:,:,L) + (MAPL_ALHL/MAPL_CP)*QS(:,:,L))
       BUOY(:,:,L) = MAPL_GRAV*BUOY(:,:,L) / ( (1.+ (MAPL_ALHL/MAPL_CP)*DQS(:,:,L))*T(:,:,L) )
    enddo

    BUOY(:,:,LM) = 0.0

    CAPE = 0.
    INHB = 0.

    do L=1,LM-1
       where(BUOY(:,:,L)>0.) 
          CAPE = CAPE + BUOY(:,:,L)*DZ(:,:,L)
       end where
       where(BUOY(:,:,L)<0.) 
          INHB = INHB - BUOY(:,:,L)*DZ(:,:,L)
       end where
    end do

    where(CAPE <= 0.0) 
       CAPE=MAPL_UNDEF
       INHB=MAPL_UNDEF
    end where

  end subroutine BUOYANCY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  function FINDLFC( BUOY, IM, JM, LM ) result( KLFC )
    ! !DESCRIPTION: 
    integer,                   intent(in)  :: IM, JM, LM
    real, dimension(IM,JM,LM), intent(in)  :: BUOY

    integer, dimension(IM,JM)              :: KLFC

    integer                                :: I, J, L

    do I = 1, IM
       do J = 1, JM

          KLFC(I,J) = 0    
          do L = LM,1,-1
             IF( BUOY(I,J,L) > 0. ) THEN 
                KLFC(I,J) = L
                EXIT
             ENDIF
          enddo

       end do
    end do
  end function FINDLFC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine FIND_EIS(TH1, QSAT, TEMP, ZET, PLO, KLCL, IM, JM, LM, LTS, EIS)
    ! !DESCRIPTION:  Returns ESrtimated Inversion Strength (K) according to Wood and Betherton, J.Climate, 2006    
   ! Written by Donifan Barahona
    
    integer                    , intent(in) :: IM,JM,LM
    real, dimension(IM,JM,LM), intent(in) :: TH1, QSAT, TEMP, PLO    
    real, dimension(IM,JM,0:LM), intent(in) :: ZET   
        
    integer, dimension(IM,JM), intent(in)     :: KLCL
    
    real, dimension(IM,JM), intent(out)     :: EIS, LTS
    real, dimension(IM,JM)        ::  Z700, ZLCL, QS700, QSLCL, T700, TLCL, GAMMA700, GAMMALCL

    integer                                 :: I, J, K


    do I = 1, IM
       do J = 1, JM

           LTS(I, J) =  0.0 
	       DO K = LM-1, 2, -1 
                  If (PLO(I, J, K) .lt. 700.0) then
                     LTS(I, J) =  TH1(I, J, K + 1)
                     Z700(I, J) = ZET(I, J, K + 1)
                     QS700(I, J) =  QSAT(I, J, K + 1) 
                     T700(I, J) =  TEMP(I, J, K + 1) 
                     exit
                  end if
             END DO
             
              LTS(I, J)  =  LTS(I, J)-TH1(I, J, LM)   
              
              ZLCL(I, J) =  ZET(I, J, KLCL(I, J)-1)             
              QSLCL(I, J) =  QSAT(I, J, KLCL(I, J)-1)
              TLCL(I, J) =  TEMP(I, J, KLCL(I, J)-1)                       
       end do
    end do
    
   
   
    GAMMA700 =  (1.0+(MAPL_ALHL*QS700/MAPL_RGAS/T700))/(1.0 + (MAPL_ALHL*MAPL_ALHL*QS700/MAPL_RVAP/T700/T700))        
    GAMMA700 =  (MAPL_GRAV/MAPL_CP)*(1.0-GAMMA700)
    GAMMALCL =  (1.0+(MAPL_ALHL*QSLCL/MAPL_RGAS/TLCL))/(1.0 + (MAPL_ALHL*MAPL_ALHL*QSLCL/MAPL_RVAP/TLCL/TLCL))        
    GAMMALCL =  (MAPL_GRAV/MAPL_CP)*(1.0-GAMMALCL)
    
    EIS =  LTS -  GAMMA700*Z700 + GAMMALCL*ZLCL
    
    end subroutine find_eis

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function FINDPBL( KH, IM, JM, LM ) result( KPBL )
    ! !DESCRIPTION: 
    integer                    , intent(in) :: IM,JM,LM
    real, dimension(IM,JM,0:LM), intent(in) :: KH

    integer, dimension(IM,JM)               :: KPBL

    integer                                 :: I, J, L
    real                                    :: KHCRIT

    KHCRIT = 2.0  ! m+2 s-1

    do I = 1, IM
       do J = 1, JM

          KPBL(I,J) = LM    
          do L = LM-1,1,-1
             IF( ( KH(I,J,L) >= KHCRIT ).AND.( KH(I,J,L-1) < KHCRIT ) ) THEN ! "top" is between L and L-1
                KPBL(I,J) = L+1   ! returned index for CBL q,t etc. is just below PBL top
                EXIT
             ENDIF
          enddo

          KPBL(I,J)=MIN( LM-1, KPBL(I,J) )      

       end do
    end do


  end function FINDPBL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  function FINDLCL( THM, QM, PL, PK, IM, JM, LM ) result( KLCL )
    ! !DESCRIPTION: 
    integer,                      intent(in) :: IM, JM, LM
    real,    dimension(IM,JM,LM), intent(in) :: THM, QM
    real,    dimension(IM,JM,LM), intent(in) :: PL, PK

    integer, dimension(IM,JM)             :: KLCL

    real, dimension(LM) :: TPCL, QSPCL
    integer             :: I, J, L, KOFFSET

    do I = 1, IM
       do J = 1, JM

          TPCL  = THM(I,J,LM) * PK(I,J,:)
          QSPCL = GEOS_QSAT(TPCL, PL(I,J,:) )

          KLCL(I,J) = 0

          do L = LM,1,-1
             if( QM(I,J,LM) >= QSPCL(L) ) then
                KLCL(I,J) = L
                exit
             endif
          enddo


          !! ------------------------------------
          !!   Disabled for Daedalus (e0203) merge
          !! ------------------------------------
          !!AMM      KOFFSET   = INT ( (LM - KLCL(I,J))/2 )   !! disable for Gan4_0
          KOFFSET   = 0
          KLCL(I,J) = MIN ( LM-1,  KLCL(I,J)+KOFFSET )
          KLCL(I,J) = MAX (    2,  KLCL(I,J)+KOFFSET )

       end do
    end do

  end function FINDLCL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE Get_hemcoFlashrate ( STATE, IMPORT, IM, JM, LM, T, PLE, ZLE, CNV_MFC, &
                                  AREA,  TS, LFR, RC )

    !=====================================================================================
    !BOP
    ! !DESCRIPTION:
    !  Wrapper routine to call HEMCO flashrate, which computes the lightning flash rate
    !  following Murray et al. 2012. 
    !
    !EOP
    ! !REVISION HISTORY
    ! 15 Jan 2020 - christoph.a.keller@nasa.gov - Initial version 
    !=====================================================================================

    ! Args
    TYPE(MAPL_MetaComp),          POINTER       :: STATE       ! Internal MAPL_Generic state
    type(ESMF_State),             INTENT(INOUT) :: IMPORT      ! Import state
    INTEGER,                      INTENT(IN)    :: IM, JM, LM
    REAL, DIMENSION(IM,JM,LM),    INTENT(IN)    :: T
    REAL, DIMENSION(IM,JM,0:LM),  INTENT(IN)    :: PLE
    REAL, DIMENSION(IM,JM,0:LM),  INTENT(IN)    :: ZLE
    REAL, DIMENSION(IM,JM,0:LM),  INTENT(IN)    :: CNV_MFC
    REAL, DIMENSION(IM,JM),       INTENT(IN)    :: AREA
    REAL, DIMENSION(IM,JM),       INTENT(IN)    :: TS
    REAL, DIMENSION(:,:),         POINTER       :: LFR
    INTEGER,                      INTENT(INOUT) :: RC

    ! Local variables
    INTEGER                          :: AGCM_IM
    REAL, ALLOCATABLE                :: LONS(:,:)
    REAL, ALLOCATABLE                :: LATS(:,:)
    REAL, ALLOCATABLE                :: LWI(:,:)
    real, pointer, dimension(:,:)    :: LONS_RAD
    real, pointer, dimension(:,:)    :: LATS_RAD
    real, pointer, dimension(:,:)    :: FRLAND 
    real, pointer, dimension(:,:)    :: FRACI
    REAL, SAVE                       :: OTDLISSCAL = -1.0

!---Initialize
    __Iam__('Get_hemcoFlashrate')
    LFR = 0.0

!---Calculate LWI, make water default value 
    CALL MAPL_GetPointer(IMPORT, FRLAND,  'FRLAND'  , __RC__ ) 
    CALL MAPL_GetPointer(IMPORT, FRACI,   'FRACI'   , __RC__ ) 
    ALLOCATE(LWI(IM,JM),STAT=RC)
    ASSERT_(RC==0)
                                       LWI = 0.0  ! Water 
    where ( FRLAND > 0.4 )             LWI = 1.0  ! Land
    where ( LWI==0.0 .and. FRACI>0.5 ) LWI = 2.0  ! Ice
    where ( LWI==0.0 .and. TS<271.40 ) LWI = 2.0  ! Ice

!---Get lat/lon in degrees
    ALLOCATE(LONS(IM,JM),LATS(IM,JM),STAT=RC)
    ASSERT_(RC==0)
    CALL MAPL_Get( STATE, lons=LONS_RAD, lats=LATS_RAD, __RC__ )
    LONS = LONS_RAD * MAPL_RADIANS_TO_DEGREES
    LATS = LATS_RAD * MAPL_RADIANS_TO_DEGREES

!---Scale factor
    IF ( OTDLISSCAL < 0.0 ) THEN
       CALL MAPL_GetResource(STATE,OTDLISSCAL,'LFR_GCC_OTDLISSCAL:',DEFAULT=-1.0, __RC__ )
       ! Estimate it from grid resolution
       IF ( OTDLISSCAL < 0.0 ) THEN
          CALL MAPL_GetResource(STATE,AGCM_IM,'AGCM_IM:', __RC__ )
          SELECT CASE ( AGCM_IM )
             CASE ( 48 )
                OTDLISSCAL = 0.355
             CASE ( 90 )
                OTDLISSCAL = 0.1
             CASE ( 180 )
                OTDLISSCAL = 1.527e-2
             CASE ( 360 )
                OTDLISSCAL = 6.32e-3
             CASE ( 720 )
                OTDLISSCAL = 1.4152e-3
             CASE DEFAULT
                OTDLISSCAL = -999.0
          END SELECT
       ENDIF
       IF ( OTDLISSCAL < 0.0 ) THEN
          WRITE(*,*) 'Invalid OTDLISSCAL scaling factor, needed to compute flash rate for GEOSCHEMchem'
          WRITE(*,*) 'Please specify parameter "LFR_GCC_OTD_LISSCAL" or make sure that there is a default'
          WRITE(*,*) 'values for this grid resolution in GEOS_MoistGridComp.F90.'
          FLUSH(6)
          ASSERT_(OTDLISSCAL>0.0)
       ENDIF
       IF ( MAPL_AM_I_ROOT() ) THEN
          WRITE(*,*) TRIM(Iam),': OTD-LIS scale factor for GCC LFR = ',OTDLISSCAL
       ENDIF 
    ENDIF

!---Get flashrate
    CALL HEMCO_FlashRate(cellArea=AREA,          &
                         lwi=LWI,                &
                         lonslocal=LONS,         &
                         latslocal=LATS,         &
                         airTemp=T,              &
                         ple=PLE,                &
                         geoPotHeight=ZLE,       &
                         cnvMfc=CNV_MFC,         &
                         otdLisScale=OTDLISSCAL, &
                         flashRate=LFR,          &
                                           __RC__ )

!---Cleanup
    IF(ALLOCATED(LONS)) DEALLOCATE(LONS)
    IF(ALLOCATED(LATS)) DEALLOCATE(LATS)
    IF(ALLOCATED(LWI))  DEALLOCATE(LWI)
    RETURN_(ESMF_SUCCESS)

END SUBROUTINE Get_hemcoFlashrate


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !DONIF Calculate the Brunt_Vaisala frequency

  !===============================================================================
  subroutine gw_prof (pcols, pver, ncol, t, pm, pi, rhoi, ni, ti, nm)
    !-----------------------------------------------------------------------
    ! Compute profiles of background state quantities for the multiple
    ! gravity wave drag parameterization.
    ! 
    ! The parameterization is assumed to operate only where water vapor 
    ! concentrations are negligible in determining the density.
    !-----------------------------------------------------------------------
    !------------------------------Arguments--------------------------------
    integer, intent(in)  :: ncol               ! number of atmospheric columns
    integer, intent(in)  :: pcols              ! number of atmospheric columns
    integer, intent(in)  :: pver               ! number of vertical layers

    !real,    intent(in)  :: u(pcols,pver)      ! midpoint zonal wind
    !real,    intent(in)  :: v(pcols,pver)      ! midpoint meridional wind
    real,    intent(in)  :: t(pcols,pver)      ! midpoint temperatures
    real,    intent(in)  :: pm(pcols,pver)     ! midpoint pressures
    real,    intent(in)  :: pi(pcols,0:pver)   ! interface pressures

    real,    intent(out) :: rhoi(pcols,0:pver) ! interface density
    real,    intent(out) :: ni(pcols,0:pver)   ! interface Brunt-Vaisalla frequency
    real,    intent(out) :: ti(pcols,0:pver)   ! interface temperature
    real,    intent(out) :: nm(pcols,pver)     ! midpoint Brunt-Vaisalla frequency

    !---------------------------Local storage-------------------------------
    integer :: ix,kx                            ! loop indexes

    real    :: dtdp
    real    :: n2, cpair, r,g                              ! Brunt-Vaisalla frequency squared
    real :: n2min   = 1.e-8 
    r=MAPL_RGAS
    cpair=MAPL_CP
    g=MAPL_GRAV

    !-----------------------------------------------------------------------------
    ! Determine the interface densities and Brunt-Vaisala frequencies.
    !-----------------------------------------------------------------------------

    ! The top interface values are calculated assuming an isothermal atmosphere 
    ! above the top level.
    kx = 0
    do ix = 1, ncol
       ti(ix,kx) = t(ix,kx+1)
       rhoi(ix,kx) = pi(ix,kx) / (r*ti(ix,kx))
       ni(ix,kx) = sqrt (g*g / (cpair*ti(ix,kx)))
    end do

    ! Interior points use centered differences
    do kx = 1, pver-1
       do ix = 1, ncol
          ti(ix,kx) = 0.5 * (t(ix,kx) + t(ix,kx+1))
          rhoi(ix,kx) = pi(ix,kx) / (r*ti(ix,kx))
          dtdp = (t(ix,kx+1)-t(ix,kx)) / (pm(ix,kx+1)-pm(ix,kx))
          n2 = g*g/ti(ix,kx) * (1./cpair - rhoi(ix,kx)*dtdp)
          ni(ix,kx) = sqrt (max (n2min, n2))
       end do
    end do

    ! Bottom interface uses bottom level temperature, density; next interface
    ! B-V frequency.
    kx = pver
    do ix = 1, ncol
       ti(ix,kx) = t(ix,kx)
       rhoi(ix,kx) = pi(ix,kx) / (r*ti(ix,kx))
       ni(ix,kx) = ni(ix,kx-1)
    end do

    !-----------------------------------------------------------------------------
    ! Determine the midpoint Brunt-Vaisala frequencies.
    !-----------------------------------------------------------------------------
    do kx=1,pver
       do ix=1,ncol
          nm(ix,kx) = 0.5 * (ni(ix,kx-1) + ni(ix,kx))
       end do
    end do

    return
  end subroutine gw_prof
!************************************************
!finds the level closets to X2=criteria 
  subroutine find_l(kp, X, crit, im, jm, lm, kmin, kmax)
  
      real, intent(in):: crit, X(im, jm, lm)
      integer, intent (in) :: im, jm, lm, kmin, kmax
      integer, intent (out) :: kp(im, jm)
      integer :: i, j , k
       
        DO j = 1, jm
	     DO i = 1, im
	           DO k = lm-1, kmin, -1
		     if ((X(i, j, k) .lt. crit) .and.  (X(i, j, k+1) .gt. crit))then      
	               kp(i, j) =  max(min(k  + 1, kmax), kmin)		       
		       exit
		     end if 
                    end do 
               end do
       end do
  end subroutine find_l
!************************************************

  !Find cloud top based on cloud fraction 

  subroutine find_cldtop(ncol, pver, cf, kcldtop)

    integer, intent(in)  :: pver , ncol ! number of vertical layers
    real,    intent(in)  :: cf(ncol,pver)     ! midpoint potential temperatures
    integer, intent(out) ::  kcldtop
    integer              ::  kuppest, ibot, k
    real                 ::  stab,  cfcrit, cf00, cfp1


    ibot = pver-1
    kcldtop  = ibot+1
    kuppest = 20   
    cfcrit = 1.0e-2


    do k =  kuppest , ibot
       cfp1 = cf(ncol, k+1)  ! qc one level down           

       if ( ( cfp1  .ge. cfcrit ) ) then
          kcldtop  = k +1 
          exit
       end if
    end do

    if (kcldtop .ge. ibot) then 
       kcldtop = pver
       return
    endif


  end subroutine find_cldtop


!Find cloud base  for a given cloud fraction 

  subroutine find_cldbase(ncol, pver, cf, kcldbase)

    integer, intent(in)  :: pver , ncol ! number of vertical layers
    real,    intent(in)  :: cf(ncol,pver)     ! midpoint potential temperatures
    integer, intent(out) ::  kcldbase
    integer              ::  kuppest, ibot, k
    real                 ::  stab,  cfcrit, cf00, cfp1


    ibot = pver-1
    kcldbase  = 20
    kuppest = 20   
    cfcrit = 1.0e-3


    do k =  ibot, kuppest, -1
       cfp1 = cf(ncol, k)  !           

       if ( ( cfp1  .ge. cfcrit ) ) then
          kcldbase  = k
          exit
       end if
    end do

    if (kcldbase .le. kuppest) then 
       kcldbase = 1
       return
    endif


  end subroutine find_cldbase
  

 
 
  subroutine VertInterp(v2,v3,ple,pp,rc)

    real    , intent(OUT) :: v2(:,:)
    real    , intent(IN ) :: v3(:,:,:)
    real    , intent(IN ) :: ple(:,:,:)
    real    , intent(IN ) :: pp
    integer, optional, intent(OUT) :: rc

    real, dimension(size(v2,1),size(v2,2)) :: al,PT,PB
    integer k,km
    logical edge

    character*(10) :: Iam='VertInterp'

    km   = size(ple,3)-1
    edge = size(v3,3)==km+1

    _ASSERT(edge .or. size(v3,3)==km,'needs informative message')

    v2   = MAPL_UNDEF

    if(EDGE) then
       pb   = ple(:,:,km+1)
       do k=km,1,-1
          pt = ple(:,:,k)
          if(all(pb<pp)) exit
          where(pp>pt .and. pp<=pb)
             al = (pb-pp)/(pb-pt)
             v2 = v3(:,:,k)*al + v3(:,:,k+1)*(1.0-al)
          end where
          pb = pt
       end do
    else
       pb = 0.5*(ple(:,:,km)+ple(:,:,km+1))
       do k=km,2,-1
          pt = 0.5*(ple(:,:,k-1)+ple(:,:,k))
          if(all(pb<pp)) exit
          where( (pp>pt.and.pp<=pb) )
             al = (pb-pp)/(pb-pt)
             v2 = v3(:,:,k-1)*al + v3(:,:,k)*(1.0-al)
          end where
          pb = pt
       end do
       pt = 0.5*(ple(:,:,km)+ple(:,:,km-1))
       pb = 0.5*(ple(:,:,km)+ple(:,:,km+1))
          where( (pp>pb.and.pp<=ple(:,:,km+1)) )
             v2 = v3(:,:,km)
          end where
    end if

    RETURN_(ESMF_SUCCESS)
  end subroutine VertInterp

   subroutine meltfrz_all  ( IM, JM, LM, Dt, TE, QV, QL, QI, &
                             CNV_FRACTION, SNOMAS, FRLANDICE, FRLAND)

      integer,                   intent(in)    :: IM, JM, LM
      real, dimension(IM,JM,LM), intent(inout) :: TE,QV,QL,QI
      real,                      intent(in)    :: Dt
      real, dimension(IM,JM),    intent(in)    :: CNV_FRACTION, SNOMAS, FRLANDICE, FRLAND
      real                   :: fQi,dQil,DQmax, QLTOT, QITOT, FQA
      integer                :: i, j, k, n
      integer, parameter     :: MaxIterations=1
      logical                :: converged
      real                   :: taufrz=450 ! timescale for freezing (seconds)

      do k=1,LM
        do j=1,JM
          do i=1,IM

      convergence: do n=1,MaxIterations

      ! melt ice using ICE_FRACTION
      fQi = ice_fraction( TE(i,j,k), CNV_FRACTION(i,j), SNOMAS(i,j), FRLANDICE(i,j), FRLAND(i,j) )
      if ( fQi < 1.0 ) then
         DQmax = (TE(i,j,k)-MAPL_TICE)*MAPL_CP/(MAPL_ALHS-MAPL_ALHL)
         dQil  = QI(i,j,k)*(1.0-fQi)
         dQil  = min(dQil, DQmax)
         dQil  = max(  0., dQil )
         QL(i,j,k) = max(QL(i,j,k)+dQil, 0.)
         QI(i,j,k) = max(QI(i,j,k)-dQil, 0.)
         TE(i,j,k) = TE(i,j,k) - (MAPL_ALHS-MAPL_ALHL)*dQil/MAPL_CP
      end if

      ! freeze liquid using ICE_FRACTION 
      fQi = ice_fraction( TE(i,j,k), CNV_FRACTION(i,j), SNOMAS(i,j), FRLANDICE(i,j), FRLAND(i,j) )
      if ( fQi > 0.0 ) then
         DQmax = (MAPL_TICE-TE(i,j,k))*MAPL_CP/(MAPL_ALHS-MAPL_ALHL)
         dQil  = QL(i,j,k)*(1.0 - EXP( -Dt * fQi / taufrz ) )
         dQil  = min(dQil, DQmax)
         dQil  = max(  0., dQil )
         QLTOT = max(QL(i,j,k)-dQil, 0.)
         QITOT = max(QI(i,j,k)+dQil, 0.)
         TE(i,j,k) = TE(i,j,k) + (MAPL_ALHS-MAPL_ALHL)*dQil/MAPL_CP
      end if

      converged = (ABS(dQil/QL(i,j,k)) <= 0.01) ! Converge at <=1% of liquid changing phase 
      if (converged) exit convergence     

      end do convergence

          enddo
        enddo
      enddo

   end subroutine meltfrz_all

   function LDRADIUS(PL,TE,QC,NN,RHX,NNL,NNI,ITYPE) RESULT(RADIUS)
   
       real, parameter :: betai     = -2.262e+3
       real, parameter :: gamai     =  5.113e+6
       real, parameter :: deltai    =  2.809e+3
       real, parameter :: densic    =  500.e0   !Ice crystal density in kgm-3

       real, parameter :: abeta = 0.07
       real, parameter :: bbeta = -0.14
       real, parameter :: bx = 100.* (3./(4.*MAPL_PI))**(1./3.)  
       real, parameter :: r13 = 1./3.  
       real, parameter :: r13bbeta = 1./3. - 0.14
       real, parameter :: premit= 750.e2            ! top pressure bound for mid level cloud
       real, parameter :: smnum = 1.e-1            !
       real, parameter :: rhminh= 0.80              ! minimum rh for high stable clouds      
       
       REAL   , INTENT(IN) :: TE,PL,NN,QC,RHX,NNL,NNI
       INTEGER, INTENT(IN) :: ITYPE
       REAL  :: RADIUS
       INTEGER, PARAMETER  :: CLOUD=1, ICE=2
       
       REAL :: NNX,RHO,IWL,BB, RIV,CLOUD_HIGH,RHDIF,WC
 
       IF(ITYPE == CLOUD) THEN        
       !- liquid cloud effective radius ----- 
          !- [liu&daum, 2000 and 2005. liu et al 2008]
          !- air density (kg/m^3)
          RHO = 100.*PL / (MAPL_RGAS*TE )
          IF(USE_AEROSOL_NN) THEN
            !- cloud drop number concentration 
            !- from the aerosol model + ....
            NNX = NNL * 1.e-6  !#/cm3
          ELSE
            !- cloud drop number concentration :u[NNX]= #/cm^3
            NNX = NN * 1.e-6  !#/cm3
          ENDIF
          !- liquid water content : u[lwl] = g/m3
          WC = RHO*QC* 1.e+3  !g/m3
          !- radius in micrometers
          RADIUS= bx *  ( WC /NNX)**r13bbeta*abeta*6.92 !6.92=(1.e-6)**bbeta	      
          !- RADIUS is limited between 2.5 and 60 micrometers as 
          !- required by rrtm parameterization
          RADIUS = max(2.5, min( 60.0, RADIUS ) )
          !- convert to meter
          RADIUS = RADIUS*1.e-6
  
       ELSEIF(ITYPE == ICE) THEN

          !------ice cloud effective radius ----- [klaus wyser, 1998]
          !- air density (kg/m^3)
          RHO = 100.*PL / (MAPL_RGAS*TE )
          !- ice water content
          iwl = RHO*QC* 1.e+3  !g/m3
          if(iwl<1.0e-6 .or. TE>273.0) then
             RADIUS = 5.0*1.e-6
          else
             BB     = -2. + log10(iwl/50.)*(1.e-3*(273.15-max(210.15,TE))**1.5)
             RADIUS =377.4 + 203.3 * bb+ 37.91 * bb **2 + 2.3696 * bb **3
             RADIUS =RADIUS * 1.e-6 !- convert to meter
          endif
 
       ELSE
          STOP "WRONG HYDROMETEOR type: CLOUD = 1 OR ICE = 2"
       ENDIF

   end function LDRADIUS

   subroutine FIX_UP_CLOUDS( TE, QV, QLLS, QILS, CLLS, QLCN, QICN, CLCN )
      real, dimension(:,:,:), intent(inout) :: TE, QV, QLLS, QILS, CLLS, QLCN, QICN, CLCN 
   ! Clouds too small
      where (CLCN < 1.E-5)
         QV   = QV + QLCN + QICN
         TE   = TE - (MAPL_ALHL/MAPL_CP)*QLCN - (MAPL_ALHS/MAPL_CP)*QICN
         CLCN = 0.
         QLCN = 0.
         QICN = 0.
      end where
      where (CLLS < 1.E-5) 
         QV   = QV + QLLS + QILS
         TE   = TE - (MAPL_ALHL/MAPL_CP)*QLLS - (MAPL_ALHS/MAPL_CP)*QILS
         CLLS = 0.
         QLLS = 0.
         QILS = 0.
      end where
   ! Liquid too small
      where ( QLCN < 1.E-8 )
         QV   = QV + QLCN
         TE   = TE - (MAPL_ALHL/MAPL_CP)*QLCN
         QLCN = 0.
      end where
      where ( QLLS < 1.E-8 )
         QV   = QV + QLLS
         TE   = TE - (MAPL_ALHL/MAPL_CP)*QLLS
         QLLS = 0.
      end where
   ! Ice too small
      where ( QICN < 1.E-8 )
         QV   = QV + QICN
         TE   = TE - (MAPL_ALHS/MAPL_CP)*QICN
         QICN = 0.
      end where
      where ( QILS < 1.E-8 )
         QV   = QV + QILS
         TE   = TE - (MAPL_ALHS/MAPL_CP)*QILS
         QILS = 0.
      end where
   ! Fix ALL clouds if cloud LIQUID+ICE too small
      where ( ( QLCN + QICN ) < 1.E-8 )
         QV   = QV + QLCN + QICN
         TE   = TE - (MAPL_ALHL/MAPL_CP)*QLCN - (MAPL_ALHS/MAPL_CP)*QICN
         CLCN = 0.
         QLCN = 0.
         QICN = 0.
      end where
      where ( ( QLLS + QILS ) < 1.E-8 )
         QV   = QV + QLLS + QILS
         TE   = TE - (MAPL_ALHL/MAPL_CP)*QLLS - (MAPL_ALHS/MAPL_CP)*QILS
         CLLS = 0.
         QLLS = 0.
         QILS = 0.
      end where

   end subroutine FIX_UP_CLOUDS

end module GEOS_MoistGridCompMod

