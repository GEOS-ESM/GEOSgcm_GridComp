! $Id$

#include "MAPL_Generic.h"

!=============================================================================
!BOP

! !MODULE: GEOS_RAS_InterfaceMod -- A Module to interface with the
!   RAS convection

module GEOS_RAS_InterfaceMod

  use ESMF
  use MAPL, r8 => MAPL_R8
  use aer_cloud
  use GEOS_UtilsMod
  use RAS       ! using module that contains ras code
  use RASPARAMS
  implicit none

  private

  character(len=ESMF_MAXSTR)              :: IAm
  integer                                 :: STATUS

  type   (RASPARAM_TYPE) :: RASPARAMS

  public :: RAS_Setup, RAS_Initialize, RAS_Run

contains

subroutine RAS_Setup (GC, CF, RC)
    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    type(ESMF_Config),   intent(inout) :: CF
    integer, optional                  :: RC  ! return code

    !EOP

    !=============================================================================
    !
    ! ErrLog Variables

    character(len=ESMF_MAXSTR)              :: IAm
    integer                                 :: STATUS
    character(len=ESMF_MAXSTR)              :: COMP_NAME

    integer      :: RFRSHINT
    integer      :: AVRGNINT

    IAm = "GEOS_RAS_InterfaceMod"
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // Iam

    call MAPL_TimerAdd(GC, name="--RAS", RC=STATUS)
    VERIFY_(STATUS)

    ! !IMPORT STATE:

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
         SHORT_NAME = 'PLE',                                         &
         LONG_NAME  = 'air_pressure',                                &
         UNITS      = 'Pa',                                          &
         DIMS       = MAPL_DimsHorzVert,                            &
         VLOCATION  = MAPL_VLocationEdge,                           &
         AVERAGING_INTERVAL = AVRGNINT,                             &
         REFRESH_INTERVAL   = RFRSHINT,                             &
         RC=STATUS  )
    VERIFY_(STATUS)                                                                          

    call MAPL_AddImportSpec(GC,                                    &
         SHORT_NAME = 'ZLE',                                       &
         LONG_NAME  = 'geopotential_height',                       &
         UNITS      = 'm',                                         &
         DIMS       =  MAPL_DimsHorzVert,                          &
         VLOCATION  =  MAPL_VLocationEdge,                         &
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
         SHORT_NAME = 'U',                                         &
         LONG_NAME  = 'eastward_wind',                             &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,                        &
         AVERAGING_INTERVAL = AVRGNINT,                            &
         REFRESH_INTERVAL   = RFRSHINT,                            &
         RC=STATUS  )
    VERIFY_(STATUS)

    ! INTERNAL STATE:
    call MAPL_AddInternalSpec(GC,                                       &
         SHORT_NAME = 'QLCN',                                           &
         LONG_NAME  = 'mass_fraction_of_convective_cloud_liquid_water', &
         UNITS      = 'kg kg-1',                                        &
         FRIENDLYTO = 'DYNAMICS:TURBULENCE',                            &
         DIMS       = MAPL_DimsHorzVert,                                &
         VLOCATION  = MAPL_VLocationCenter,                  RC=STATUS  )  
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
         SHORT_NAME = 'CLCN',                                      &
         LONG_NAME  = 'convective_cloud_area_fraction',            &
         UNITS      = '1',                                         &
         FRIENDLYTO = 'DYNAMICS',                                  &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)                                                                          

    ! EXPORT STATE:
    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME= 'CNV_DQLDT ',                                  &
         LONG_NAME = 'convective_condensate_source',                &
         UNITS     = 'kg m-2 s-1',                                  &
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
         SHORT_NAME= 'CNV_PRC3 ',                                   &
         LONG_NAME = 'convective_precipitation_from_RAS',           &
         UNITS     = 'kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
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
         SHORT_NAME='ENTLAM',                                      &
         LONG_NAME ='entrainment parameter',                       &
         UNITS     ='m-1',                                  &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
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

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='BYNCY',                                       & 
         LONG_NAME ='buoyancy_of surface_parcel',                  &
         UNITS     ='m s-2',                                       &
         DIMS      = MAPL_DimsHorzVert,                            & 
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

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

end subroutine RAS_Setup

subroutine RAS_Initialize (MAPL, RC)
    type (MAPL_MetaComp), intent(inout) :: MAPL
    integer, optional                   :: RC  ! return code

    character(len=ESMF_MAXSTR) :: GRIDNAME
    character(len=4)           :: imchar
    character(len=2)           :: dateline
    integer                    :: imsize,nn

!#ifdef NODISABLE

      call MAPL_GetResource(MAPL,GRIDNAME,'AGCM_GRIDNAME:', RC=STATUS)
      VERIFY_(STATUS)
      GRIDNAME =  AdjustL(GRIDNAME)
      nn = len_trim(GRIDNAME)
      dateline = GRIDNAME(nn-1:nn)
      imchar = GRIDNAME(3:index(GRIDNAME,'x')-1)
      read(imchar,*) imsize
      if(dateline.eq.'CF') imsize = imsize*4

      call MAPL_GetResource(MAPL, RASPARAMS%CUFRICFAC,     'CUFRICFAC:',      DEFAULT= 1.000, RC=STATUS)
      call MAPL_GetResource(MAPL, RASPARAMS%SHR_LAMBDA_FAC,'SHR_LAMBDA_FAC:', DEFAULT= 0.05,  RC=STATUS)
      call MAPL_GetResource(MAPL, RASPARAMS%QC_CRIT_CN,    'QC_CRIT_CN:',    DEFAULT= 8.0e-4,  RC=STATUS)
      call MAPL_GetResource(MAPL, RASPARAMS%RASAL1,        'RASAL1:',        DEFAULT=   1800.0,RC=STATUS)
      call MAPL_GetResource(MAPL, RASPARAMS%RASAL2,        'RASAL2:',        DEFAULT= -86400.0,RC=STATUS)
      call MAPL_GetResource(MAPL, RASPARAMS%RASNCL,        'RASNCL:',        DEFAULT= -300.,   RC=STATUS)
      call MAPL_GetResource(MAPL, RASPARAMS%LAMBDA_FAC,    'LAMBDA_FAC:',    DEFAULT= 4.0 ,    RC=STATUS)
      call MAPL_GetResource(MAPL, RASPARAMS%LAMBMX_FAC,    'LAMBMX_FAC:',    DEFAULT= 0.0 ,    RC=STATUS)
      call MAPL_GetResource(MAPL, RASPARAMS%STRAPPING,     'STRAPPING:',     DEFAULT=-1, RC=STATUS)
      call MAPL_GetResource(MAPL, RASPARAMS%MIN_DIAMETER,  'MIN_DIAMETER:',  DEFAULT= 400.,    RC=STATUS)
      call MAPL_GetResource(MAPL, RASPARAMS%CUFRICLAMBDA,  'CUFRICLAMBDA:',  DEFAULT= 7.5e-4,  RC=STATUS)
      call MAPL_GetResource(MAPL, RASPARAMS%RDTLEXPON,     'RDTLEXPON:',     DEFAULT= 1.0 ,    RC=STATUS)
      call MAPL_GetResource(MAPL, RASPARAMS%SDQV2,         'SDQV2:',         DEFAULT= 1.3 ,    RC=STATUS)
      call MAPL_GetResource(MAPL, RASPARAMS%SDQV3,         'SDQV3:',         DEFAULT= 1.3 ,    RC=STATUS)
      call MAPL_GetResource(MAPL, RASPARAMS%SDQVT1,        'SDQVT1:',        DEFAULT= 263.,    RC=STATUS)
      call MAPL_GetResource(MAPL, RASPARAMS%ACRITFAC,      'ACRITFAC:',      DEFAULT= 0.5 ,    RC=STATUS)
      call MAPL_GetResource(MAPL, RASPARAMS%HMINTRIGGER,   'HMINTRIGGER:',   DEFAULT= 1.0 ,    RC=STATUS)
      call MAPL_GetResource(MAPL, RASPARAMS%LLDISAGGXP,    'LLDISAGGXP:',    DEFAULT= 0.0 ,    RC=STATUS)
      call MAPL_GetResource(MAPL, RASPARAMS%PBLFRAC,       'PBLFRAC:',       DEFAULT= 0.1 ,    RC=STATUS)
      call MAPL_GetResource(MAPL, RASPARAMS%RASAUTORAMPB,  'RASAUTORAMPB:',  DEFAULT= 0.8 ,    RC=STATUS)
      call MAPL_GetResource(MAPL, RASPARAMS%AUTOC_CN_ZDEP, 'AUTOC_CN_ZDEP:', DEFAULT= 1.0          ,RC=STATUS)
      if( imsize.le.200       ) call MAPL_GetResource(MAPL, RASPARAMS%MAXDALLOWED_S, 'MAXDALLOWED_S:', DEFAULT= 4000.0 ,RC=STATUS)
      if( imsize.gt.200 .and. &
          imsize.le.400       ) call MAPL_GetResource(MAPL, RASPARAMS%MAXDALLOWED_S, 'MAXDALLOWED_S:', DEFAULT= 2000.0 ,RC=STATUS)
      if( imsize.gt.400 .and. &
          imsize.le.800       ) call MAPL_GetResource(MAPL, RASPARAMS%MAXDALLOWED_S, 'MAXDALLOWED_S:', DEFAULT=  700.0 ,RC=STATUS)
      if( imsize.gt.800 .and. &
          imsize.le.1600      ) call MAPL_GetResource(MAPL, RASPARAMS%MAXDALLOWED_S, 'MAXDALLOWED_S:', DEFAULT=  450.0 ,RC=STATUS)
      if( imsize.gt.1600      ) call MAPL_GetResource(MAPL, RASPARAMS%MAXDALLOWED_S, 'MAXDALLOWED_S:', DEFAULT=  450.0 ,RC=STATUS)
                                call MAPL_GetResource(MAPL, RASPARAMS%MAXDALLOWED_D, 'MAXDALLOWED_D:', DEFAULT= RASPARAMS%MAXDALLOWED_S ,RC=STATUS)
                                call MAPL_GetResource(MAPL, RASPARAMS%MAXDALLOWED_E, 'MAXDALLOWED_E:', DEFAULT= -0.500 ,RC=STATUS)
      call MAPL_GetResource(MAPL, RASPARAMS%RASAL_SLOPE ,  'RASAL_SLOPE:',     DEFAULT= 8000.0  ,  RC=STATUS)
      call MAPL_GetResource(MAPL, RASPARAMS%RAS_RHMIN ,    'RAS_RHMIN:',       DEFAULT= 0.5  ,  RC=STATUS)
      call MAPL_GetResource(MAPL, RASPARAMS%RAS_RHFULL,    'RAS_RHFULL:',      DEFAULT= 0.65 ,  RC=STATUS)
!#endif

end subroutine RAS_Initialize

subroutine RAS_Run (GC, IMPORT, EXPORT, CLOCK, RC)
    type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
    type(ESMF_State),    intent(inout) :: IMPORT ! Import state
    type(ESMF_State),    intent(inout) :: EXPORT ! Export state
    type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
    integer, optional,   intent(  out) :: RC     ! Error code:

    ! Local derived type aliases

    type (MAPL_MetaComp), pointer   :: MAPL
    type (ESMF_Config  )            :: CF
    type (ESMF_State   )            :: INTERNAL
    type (ESMF_Alarm   )            :: ALARM
    type (ESMF_TimeInterval)        :: TINT
    real(ESMF_KIND_R8)              :: DT_R8
  
    ! Local variables

    integer                         :: IM, JM, LM, L, I, J, kii
    integer                         :: IDIM, IRUN, K0, ICMIN
    integer                         :: ITRCR,KSTRAP,CBL_METHOD,KCBLMIN,CLEANUP_RH
    real                            :: PMIN_DET
    real                            :: PMIN_CBL
    real                            :: DT_MOIST
    real                            :: CBL_TPERT, CBL_QPERT
    real                            :: CBL_TPERT_MXOCN,CBL_TPERT_MXLND
    real                            :: AUTOC_CN_OCN, AUTOC_CN_LAND
    real                            :: FDROP_DUST

    !Whether to guard against negatives
    logical                         :: RAS_NO_NEG

    ! pointers
    real, pointer, dimension(:)     :: PREF
    real, pointer, dimension(:    ) :: FSCAV      ! container for friendly to moist tracers
    real, pointer, dimension(:,:)     :: LONS
    real, pointer, dimension(:,:)     :: LATS
    real, pointer, dimension(:,:)     :: TS,SNOMAS,FRLAND
    real, pointer, dimension(:,:)     :: RAS_TIME, RAS_TRG, RAS_TOKI, RAS_PBL, RAS_WFN 
    real, pointer, dimension(:,:)     :: ZLCL, ZLFC
    real, pointer, dimension(:,:)     :: TPERTI,KCBLI
    real, pointer, dimension(:,:)     :: TVEX 
    real, pointer, dimension(:,:)     :: MXDIAM
    real, pointer, dimension(:,:,:)   :: TH, PLE
    real, pointer, dimension(:,:,:)   :: Q, QLCN, CLCN, BYNCY, QICN
    real, pointer, dimension(:,:,:)   :: CNV_DQLDT, CNV_MF0, CNV_MFD , CNV_MFC
    real, pointer, dimension(:,:,:)   :: CNV_UPDF, CNV_CVW, CNV_QC
    real, pointer, dimension(:,:,:)   :: ENTLAM
    real, pointer, dimension(:,:,:)   :: RAS_ALPHA, RAS_TAU
    real, pointer, dimension(:,:,:)   :: CNV_FICE, CNV_NDROP, CNV_NICE
    real, pointer, dimension(:,:,:)   :: KH
    real, pointer, dimension(:,:,:)   :: ZL0, V, U
    real, pointer, dimension(:,:,:)   :: THOI,QHOI,QSSI,DQSI
    real, pointer, dimension(:,:,:)   :: PLEI
    real, pointer, dimension(:,:,:)   :: THRAS,URAS,VRAS
    real, pointer, dimension(:,:)     ::  TRIEDLV
    real, pointer, dimension(:,:,:)   :: RCCODE,QVRAS

    ! allocatable arrays
    real, allocatable, dimension(:)     :: SIGE
    real, allocatable, dimension(:,:)   :: RASPRCP
    real, allocatable, dimension(:,:)   :: ZCBLx, MXDIAMx
    real, allocatable, dimension(:,:)   :: TPERT, QPERT, RASAL2_2d
    real, allocatable, dimension(:,:)   :: CNV_FRACTION
    real, allocatable, dimension(:,:)   :: CO_AUTO
    real, allocatable, dimension(:,:)   :: QSSFC
    real, allocatable, dimension(:,:,:) :: DQS, QSS, PLO, ZLO, TEMP, PK
    real, allocatable, dimension(:,:,:) :: Q1, W1, U1, V1, TH1, CNV_PRC3,fQi,CFPBL,CNV_HAIL
    real, allocatable, dimension(:,:,:) :: WGT0, WGT1
    real, allocatable, dimension(:,:,:) :: GZLO, HHO,HSO
    real, allocatable, dimension(:,:,:) :: GZLE
    real, allocatable, dimension(:,:,:) :: PKE
    real, allocatable, dimension(:,:,:) :: CNV_PLE,ZLE
    real, allocatable, dimension(:,:)   :: KEX
    real, allocatable, dimension(:,:,:) :: MASS
    real, allocatable, dimension(:,:,:) :: XHO
    real, allocatable, dimension(:,:) :: TRDLX

    integer, allocatable, dimension(:,:) :: IRAS, JRAS, KCBL ! Need IM and JM to initialize
    integer, allocatable, dimension(:,:) :: KLCL, KLFC, KPBL, KPBL_SC
    integer, allocatable, dimension(:,:,:) :: SEEDINI ! Need IM and JM to initialize
    integer, allocatable, dimension(:,:,:) :: IRCCODE

    type  (AerProps), allocatable, dimension(:,:,:) :: AeroProps !Storages aerosol properties for activation 

    call ESMF_GridCompGet( GC, CONFIG=CF, RC=STATUS ) 
    VERIFY_(STATUS)

    ! Get my internal MAPL_Generic state
    !-----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_TimerOn (MAPL,"--RAS")

    ! Get parameters from generic state.
    !-----------------------------------
  
    call MAPL_Get( MAPL, IM=IM, JM=JM, LM=LM,   &
         RUNALARM = ALARM,             &
         CF       = CF,                &
         LONS     = LONS,              &
         LATS     = LATS,              &
         INTERNAL_ESMF_STATE=INTERNAL, &
         RC=STATUS )
    VERIFY_(STATUS)
    
    call ESMF_AlarmGet(ALARM, RingInterval=TINT, RC=STATUS); VERIFY_(STATUS)
    call ESMF_TimeIntervalGet(TINT,   S_R8=DT_R8,RC=STATUS); VERIFY_(STATUS)

!#ifdef NODISABLE

    ! get resources

    call MAPL_GetResource(MAPL, RAS_NO_NEG, 'RAS_NO_NEG:', default=.FALSE. , RC=STATUS)
    call MAPL_GetResource(MAPL, PMIN_DET, 'PMIN_DET', DEFAULT= 3000.0, RC=STATUS)
    call MAPL_GetResource(MAPL, CBL_METHOD, 'CBL_METHOD:', DEFAULT= 6, RC=STATUS)
    call MAPL_GetResource(MAPL, PMIN_CBL, 'PMIN_CBL', DEFAULT= 50000.0, RC=STATUS)
    call MAPL_GetResource(MAPL, CBL_QPERT, 'CBL_QPERT:', DEFAULT= 0.0, RC=STATUS)
    call MAPL_GetResource(MAPL, CBL_TPERT, 'CBL_TPERT:',       DEFAULT=-1.0   , RC=STATUS)
    call MAPL_GetResource(MAPL, CBL_TPERT_MXOCN, 'CBL_TPERT_MXOCN:', DEFAULT= 2.0   , RC=STATUS)
    call MAPL_GetResource(MAPL, CBL_TPERT_MXLND, 'CBL_TPERT_MXLND:', DEFAULT= 0.0   , RC=STATUS)
    call MAPL_GetResource(MAPL, CBL_QPERT, 'CBL_QPERT:', DEFAULT= 0.0   , RC=STATUS)
    call MAPL_GetResource(MAPL, AUTOC_CN_OCN, 'AUTOC_CN:', DEFAULT= 2.5e-3, RC=STATUS)
    call MAPL_GetResource(MAPL, AUTOC_CN_LAND, 'AUTOC_CN_LAND:', DEFAULT= AUTOC_CN_OCN, RC=STATUS)
    call MAPL_GetResource(MAPL, FDROP_DUST,     'FDROP_DUST:',     DEFAULT= 0.5,    RC=STATUS) !Fraction of dust within droplets for immersion freezing
    
    ! get pointers

    call MAPL_GetPointer(IMPORT, TH, 'TH', RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, PREF, 'PREF', RC=STATUS); VERIFY_(STATUS) 
    call MAPL_GetPointer(IMPORT, PLE, 'PLE' , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, BYNCY, 'BYNCY', RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, KH, 'KH', RC=STATUS); VERIFY_(STATUS)
    ! frland =0.0
    call MAPL_GetPointer(IMPORT, TS, 'TS', RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QLCN, 'QLCN', RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QICN, 'QICN', RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, CNV_DQLDT,'CNV_DQLDT',RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, CNV_MF0,  'CNV_MF0' , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, CNV_MFD,  'CNV_MFD' , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, CNV_MFC,  'CNV_MFC' , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, CNV_MFD,  'CNV_MFD' , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, CNV_UPDF, 'CNV_UPDF', RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, CNV_CVW,  'CNV_CVW' , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, CNV_QC,   'CNV_QC'  , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, ENTLAM,     'ENTLAM'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, CLCN,     'CLCN'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, RAS_TIME,  'RAS_TIME' , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, RAS_TRG,   'RAS_TRG'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, RAS_TOKI,  'RAS_TOKI' , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, RAS_PBL,   'RAS_PBL'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, RAS_WFN,   'RAS_WFN'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, RAS_ALPHA, 'RAS_ALPHA' , ALLOC=.TRUE.  , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, RAS_TAU, 'RAS_TAU' , ALLOC=.TRUE.  , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, CNV_FICE,  'CNV_FICE' , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, CNV_NDROP, 'CNV_NDROP' , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, CNV_NICE,  'CNV_NICE' , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, ZL0,     'ZLE'     , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, V,       'V'       , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, U,       'U'       , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, ZLCL,     'ZLCL'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, ZLFC,     'ZLFC'    , RC=STATUS); VERIFY_(STATUS)

    ! allocate arrays
    if(.not.allocated(SEEDINI)) allocate(SEEDINI(IM,JM,2)); VERIFY_(STATUS)
    if(.not.allocated(IRAS)) allocate(IRAS(IM,JM)); VERIFY_(STATUS)
    if(.not.allocated(JRAS)) allocate(JRAS(IM,JM)); VERIFY_(STATUS)
    if(.not.allocated(SIGE)) allocate(SIGE(0:LM)); VERIFY_(STATUS)
    if(.not.allocated(KCBL)) allocate(KCBL(IM,JM)); VERIFY_(STATUS)
    if(.not.allocated(WGT0)) allocate(WGT0(IM,JM,LM)); VERIFY_(STATUS)
    if(.not.allocated(WGT1)) allocate(WGT1(IM,JM,LM)); VERIFY_(STATUS)
    if(.not.allocated(TH1)) allocate(TH1(IM,JM,LM)); VERIFY_(STATUS)
    if(.not.allocated(Q1)) allocate(Q1(IM,JM,LM)); VERIFY_(STATUS)
    if(.not.allocated(U1)) allocate(U1(IM,JM,LM)); VERIFY_(STATUS)
    if(.not.allocated(V1)) allocate(V1(IM,JM,LM)); VERIFY_(STATUS)
    if(.not.allocated(DQS)) allocate(DQS(IM,JM,LM)); VERIFY_(STATUS)
    if(.not.allocated(QSS)) allocate(QSS(IM,JM,LM)); VERIFY_(STATUS)
    if(.not.allocated(PLO)) allocate(PLO(IM,JM,LM)); VERIFY_(STATUS)
    if(.not.allocated(ZLO)) allocate(ZLO(IM,JM,LM)); VERIFY_(STATUS)
    if(.not.allocated(TEMP)) allocate(TEMP(IM,JM,LM)); VERIFY_(STATUS)
    if(.not.allocated(PK)) allocate(PK(IM,JM,LM)); VERIFY_(STATUS)
    if(.not.allocated(KLCL)) allocate(KLCL(IM,JM)); VERIFY_(STATUS) 
    if(.not.allocated(KLFC)) allocate(KLFC(IM,JM)); VERIFY_(STATUS) 
    if(.not.allocated(KPBL)) allocate(KPBL(IM,JM)); VERIFY_(STATUS) 
    if(.not.allocated(KPBL_SC)) allocate(KPBL_SC(IM,JM)); VERIFY_(STATUS) 
    if(.not.allocated(CNV_PLE)) allocate(CNV_PLE(IM,JM,0:LM)); VERIFY_(STATUS)
    if(.not.allocated(ZLE)) allocate(ZLE(IM,JM,0:LM)); VERIFY_(STATUS)
    if(.not.allocated(ZCBLx)) allocate(ZCBLx(IM,JM)); VERIFY_(STATUS) 
    if(.not.allocated(MXDIAMx)) allocate(MXDIAMx(IM,JM)); VERIFY_(STATUS) 
    if(.not.allocated(TPERT)) allocate(TPERT(IM,JM)); VERIFY_(STATUS) 
    if(.not.allocated(QPERT)) allocate(QPERT(IM,JM)); VERIFY_(STATUS) 
    if(.not.allocated(CNV_FRACTION)) allocate(CNV_FRACTION(IM,JM)); VERIFY_(STATUS) 
    if(.not.allocated(RASAL2_2d)) allocate(RASAL2_2d(IM,JM)); VERIFY_(STATUS) 
    if(.not.allocated(CO_AUTO)) allocate(CO_AUTO(IM,JM)); VERIFY_(STATUS) 
    if(.not.allocated(GZLO)) allocate(GZLO(IM,JM,LM)); VERIFY_(STATUS)
    if(.not.allocated(GZLE)) allocate(GZLE(IM,JM,0:LM)); VERIFY_(STATUS) 
    if(.not.allocated(HHO)) allocate(HHO(IM,JM,LM)); VERIFY_(STATUS)
    if(.not.allocated(HSO)) allocate(HSO(IM,JM,LM)); VERIFY_(STATUS)
    if(.not.allocated(PKE)) allocate(PKE(IM,JM,0:LM)); VERIFY_(STATUS) 
    if(.not.allocated(CNV_PRC3)) allocate(CNV_PRC3(IM,JM,LM)); VERIFY_(STATUS) 
    if(.not.allocated(RASPRCP)) allocate(RASPRCP(IM,JM)); VERIFY_(STATUS)
    if(.not.allocated(AeroProps)) allocate(AeroProps(IM,JM,LM)); VERIFY_(STATUS)
    if(.not.allocated(IRCCODE)) allocate(IRCCODE(IM,JM,LM)); VERIFY_(STATUS)
    if(.not.allocated(KEX)) allocate(KEX(IM,JM)); VERIFY_(STATUS) 
    if(.not.allocated(RASAL2_2d)) allocate(RASAL2_2d(IM,JM)); VERIFY_(STATUS)
    if(.not.allocated(QSSFC)) allocate(QSSFC(IM,JM)); VERIFY_(STATUS)
    if(.not.allocated(XHO)) allocate(XHO (IM,JM,ITRCR),stat=STATUS); VERIFY_(STATUS)

   ! initialize local variables

    ITRCR       = 0 ! This meeds to be before XHO is allocated.
    
    if(.not.associated(FSCAV)) allocate(FSCAV(ITRCR),stat=STATUS)

    !  Copy incoming state vars to local arrays that will be adjusted
    !  by physics.  Untouched state vars will later be used for 
    !  post facto tendency calculations.
    !----------------------------------------------------------------

    if(associated(CNV_DQLDT)) CNV_DQLDT =  0.0
    if(associated(TH)) TH1 = TH 
    if(associated(Q)) Q1 = Q 
    if(associated(V)) V1 = V 
    if(associated(U)) U1 = U 

    DT_MOIST = DT_R8
    IDIM = IM * JM 
    IRUN = IM * JM 
    K0 = LM 
    KCBLMIN  = count(PREF < PMIN_CBL)
    CNV_PLE  = PLE*.01
    PLO      = 0.5*(CNV_PLE(:,:,0:LM-1) +  CNV_PLE(:,:,1:LM  ) )
    PKE      = (      PLE/MAPL_P00)**(MAPL_KAPPA)
    PK       = (100.0*PLO/MAPL_P00)**(MAPL_KAPPA)
    CNV_FRACTION = 0.0 
    CNV_PRC3 = 0.0    
    ENTLAM   = 0.0
    !-srf - placed here before the cumulus parameterizations are called.
    do L=LM,0,-1
       ZLE(:,:,L) = ZL0(:,:,L) - ZL0(:,:,LM)
    end do
    do L=LM,1,-1
       ZLO(:,:,L) = 0.5*(ZLE(:,:,L-1) + ZLE(:,:,L))
    end do

    GZLE  = MAPL_GRAV * ZLE
    GZLO  = MAPL_GRAV * ZLO

    CNV_UPDF = 0.0    
    CNV_CVW  = 0.0
    RASPRCP  = 0.0

    HHO      =  0.0
    HSO      =  0.0    
    IRCCODE  = -99
    TRDLX    = 1.0

    ! Some export work prior to doing calculations
    !---------------------------------------------
    if(associated(CNV_MFC)) CNV_MFC(:,:,LM) = 0.

    ! Determine how to do cloud base 
    !-------------------------------

    KLCL = FINDLCL( TH1, Q1, PLO, PK, IM, JM, LM )
    KLFC = FINDLFC( BYNCY, IM, JM, LM )
    KPBL = FINDPBL( KH, IM, JM, LM )

    IRAS       = nint(LONS*100)
    JRAS       = nint(LATS*100)
    RASPARAMS%CLDMICRO = 0.0
!    if(adjustl(CLDMICRO)=="MGB2_2M") then !WDB Maybe, this block is unnecessary.
!       RASPARAMS%CLDMICRO = 1.0
!       RASPARAMS%FDROP_DUST = FDROP_DUST !WDB Note multiple statements.
!       RASPARAMS%FDROP_SOOT = FDROP_SOOT
!       RASPARAMS%FDROP_SEASALT = SS_INFAC
!    endif
    SEEDINI(:,:,1) = 1000000 * ( 100*TEMP(:,:,LM)   - INT( 100*TEMP(:,:,LM) ) )
    SEEDINI(:,:,2) = 1000000 * ( 100*TEMP(:,:,LM-1) - INT( 100*TEMP(:,:,LM-1) ) )

    KSTRAP = INT( RASPARAMS%STRAPPING )

    if (RASPARAMS%RASAL2 > 0.0) then
       RASAL2_2d(:,:) = RASPARAMS%RASAL2
    else
       ! include CNV dependence
       DO J=1, JM
          DO I=1, IM
          RASAL2_2d(I,J) = CNV_FRACTION(I,J)*ABS(RASPARAMS%RASAL2) + (1-CNV_FRACTION(I,J))*RASPARAMS%RASAL1
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

      ICMIN    = max(1,count(PREF < PMIN_DET))
      ! temp kluge to differentiate ocean,land convective autoc (jtb 6/29/05)
      ! ---------------------------------------------------------------------
      where (FRLAND<0.1)
         CO_AUTO = AUTOC_CN_OCN   ! ocean value
      elsewhere
         CO_AUTO = AUTOC_CN_LAND  ! land value
      end where

      RAS_ALPHA   = MAPL_UNDEF
      RAS_TAU     = MAPL_UNDEF
      RAS_TIME    = MAPL_UNDEF
      RAS_TRG     = MAPL_UNDEF
      RAS_TOKI    = MAPL_UNDEF
      RAS_PBL     = MAPL_UNDEF
      RAS_WFN     = MAPL_UNDEF
      CNV_FICE    = 0.0    
      CNV_NICE    = 0.0
      CNV_NDROP   = 0.0 !(1e-6 cm-3)

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
           SEEDINI              , &
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
           IRCCODE              , &
           XHO = XHO            , &
           TRIEDLEV_DIAG = TRDLX, &
           FSCAV  = FSCAV       , &
           DISSKE = KEX           )

      if(associated(QVRAS  )) QVRAS    = Q1
      if(associated(THRAS  )) THRAS    = TH1
      if(associated(URAS   )) URAS     = U1
      if(associated(VRAS   )) VRAS     = V1
      if(associated(MXDIAM )) MXDIAM  = MXDIAMx
      if(associated(RCCODE )) RCCODE  = 1.0*IRCCODE
      if(associated(TRIEDLV)) TRIEDLV = TRDLX
      if(associated(TVEX   )) TVEX     = SUM( (MAPL_CP*TEMP + MAPL_ALHL*Q)*MASS, 3 )
!#endif
      if(associated(FSCAV))     deallocate(FSCAV, stat=STATUS); VERIFY_(STATUS)

    call MAPL_TimerOff(MAPL,"--RAS")

end subroutine RAS_Run

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
end module GEOS_RAS_InterfaceMod
