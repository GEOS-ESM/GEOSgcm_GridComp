!  $Id: GEOS_OceanbiogeochemGridComp.F90,v 1.5.2.12.8.9.2.8.6.14.2.2.8.5.8.2.4.4 2018/10/23 13:44:33 croussea Exp $

#include "MAPL_Generic.h"
!=============================================================================

module GEOS_OceanbiogeochemGridCompMod

!BOP

! !MODULE: GEOS_OceanbiogeochemGridCompMod -- Implements ocean biology

! !USES:
  use ESMF
  use MAPL

  implicit none
  private

! !PUBLIC ROUTINES:
  public SetServices

  integer            :: DO_DATA_ATM4OCN
  integer            :: NUM_ICE_CATEGORIES
  integer, parameter :: NUM_DUDP           = 5
  integer, parameter :: NUM_DUWT           = 5
  integer, parameter :: NUM_DUSD           = 5
  integer, parameter :: NUM_BCDP           = 2
  integer, parameter :: NUM_BCWT           = 2
  integer, parameter :: NUM_OCDP           = 2
  integer, parameter :: NUM_OCWT           = 2
  real, parameter    :: MOLWGHT_FE = 55.8452 ! g/mol

! !DEBUG-ing flags:
!------------------------------------------

!  integer, save :: DEBUG

! Private state containing static constants
!------------------------------------------
#include "definebio.h"
    
  type T_OBIO_STATE
    private
    real :: rmumax(nchl)     !max phyto growth rate at 20oC, /d
    real :: cnratio,cfratio,drate,Rm,rlamz,greff,dratez1,dratez2
    real :: regen,solFe
    real :: remin(ndet)      !detrital remineralization rate /s
    real :: rkn(nchl),rks(nchl) !half-sat constants for N and S (uM)
    real :: rkf(nchl)        !half-sat constant for Fe (nM)
    real :: rik(3,nchl)      !light saturation parameter umol quanta/m2/s
    real :: cchl(3)          !C:chl ratio g:g for 3 light states
!    real :: fescavrate(2)    !Fe scavenging rate /s
    real :: wss(ntyp)        !particle sinking rate m/s
    real :: axs(ntyp,2)      !detrital sinking a and b coefficients
    real :: Pdeep(ntyp)      !deep BC's (not needed for MOM)
  end type T_OBIO_STATE

   type OBIO_WRAP
      type (T_OBIO_STATE), pointer :: PTR => null()
   end type OBIO_WRAP

!=============================================================================

! !DESCRIPTION:
! 
!   {\tt GEOS\_Obio} is a light-weight gridded component does ocean biology
!

!EOP
   contains
!BOP

! ! IROUTINE: SetServices -- Sets ESMF services for this component

! ! INTERFACE:
  subroutine SetServices ( GC, RC )

! ! ARGUMENTS:
    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    integer, optional                  :: RC  ! return code

! ! DESCRIPTION: This version uses the MAPL_GenericSetServices, which sets
!                the Initialize and Finalize services, as well as allocating
!   our instance of a generic state and putting it in the 
!   gridded component (GC). Here we only need to set the run method and
!   add the state variable specifications (also generic) to our instance
!   of the generic state. This is the way our true state variables get into
!   the ESMF_State INTERNAL, which is in the MAPL_MetaComp. The import
!   and internal variables are allocated and initialized by generic.  Here
!   generic is used for tiles.

!EOP

!=============================================================================
! ErrLog Variables

    character(len=ESMF_MAXSTR)              :: IAm
    integer                                 :: STATUS
    character(len=ESMF_MAXSTR)              :: COMP_NAME

    type (MAPL_MetaComp),  pointer          :: MAPL
    type  (ESMF_Config)                     :: CF
    integer                                 :: DO_CICE_THERMO

!=============================================================================

! Begin...

! Get my name and set-up traceback handle
! ---------------------------------------

    call ESMF_GridCompGet( GC, NAME=COMP_NAME, CONFIG=CF, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // 'SetServices'

  
! Set the state variable specs.
! -----------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_GetResource ( MAPL,       DO_CICE_THERMO,     Label="USE_CICE_Thermo:" ,       DEFAULT=0, RC=STATUS)
    VERIFY_(STATUS)

! Get constants from CF
! ---------------------

    if (DO_CICE_THERMO /= 0) then     ! Before merging CICEthermo with SaltWater, following was set via makefile.
       call ESMF_ConfigGetAttribute(CF, NUM_ICE_CATEGORIES, Label="CICE_N_ICE_CATEGORIES:" , RC=STATUS)
       VERIFY_(STATUS)
    else
       NUM_ICE_CATEGORIES = 1
    endif

    call MAPL_GetResource ( MAPL, DO_DATA_ATM4OCN, Label="USE_DATA_ATM4OCN:" , DEFAULT=0, _RC)

! Set the Run entry point
! -----------------------

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_INITIALIZE, Initialize, RC=STATUS )
    VERIFY_(STATUS)
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,  Run, RC=STATUS )
    VERIFY_(STATUS)

! Set the state variable specs.
! -----------------------------

!BOC

!  !EXPORT STATE:

    call MAPL_AddExportSpec(GC,                               &
    SHORT_NAME = 'PCO2',                                      &
    LONG_NAME  = 'partial_pressure_of_co2',                   &
    UNITS      = 'uatm',                                      &
    DIMS       = MAPL_DimsHorzOnly,                           &
    VLOCATION  = MAPL_VLocationNone,                          &
    RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
    SHORT_NAME = 'FCO2',                                      &
    LONG_NAME  = 'co2 flux atm-to-ocean',                     &
    UNITS      = 'ugC m-2 s-1',                               &
    DIMS       = MAPL_DimsHorzOnly,                           &
    VLOCATION  = MAPL_VLocationNone,                          &
    RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
    SHORT_NAME = 'BIO',                                       &
    LONG_NAME  = 'bio variables',                             &
    UNITS      = 'variable',                                  &
    DIMS       = MAPL_DimsHorzOnly,                           &
    VLOCATION  = MAPL_VLocationCenter,                        &
    DATATYPE   = MAPL_BundleItem,                             &
    RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
    SHORT_NAME = 'PPDIATOM',                                  &
    LONG_NAME  = 'primary_production_diatoms',                &
    UNITS      = 'mg C m-2 d-1',                              &
    DIMS       = MAPL_DimsHorzOnly,                           &
    VLOCATION  = MAPL_VLocationNone,                          &
    RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
    SHORT_NAME = 'PPCHLORO',                                  &
    LONG_NAME  = 'primary_production_chlorophytes',           &
    UNITS      = 'mg C m-2 d-1',                              &
    DIMS       = MAPL_DimsHorzOnly,                           &
    VLOCATION  = MAPL_VLocationNone,                          &
    RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
    SHORT_NAME = 'PPCYANO',                                   &
    LONG_NAME  = 'primary_production_cyanobacteria',          &
    UNITS      = 'mg C m-2 d-1',                              &
    DIMS       = MAPL_DimsHorzOnly,                           &
    VLOCATION  = MAPL_VLocationNone,                          &
    RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
    SHORT_NAME = 'PPCOCCO',                                   &
    LONG_NAME  = 'primary_production_coccolithophores',       &
    UNITS      = 'mg C m-2 d-1',                              &
    DIMS       = MAPL_DimsHorzOnly,                           &
    VLOCATION  = MAPL_VLocationNone,                          &
    RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
    SHORT_NAME = 'PPDINO',                                    &
    LONG_NAME  = 'primary_production_dinoflagellates',        &
    UNITS      = 'mg C m-2 d-1',                              &
    DIMS       = MAPL_DimsHorzOnly,                           &
    VLOCATION  = MAPL_VLocationNone,                          &
    RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
    SHORT_NAME = 'PPPHAEO',                                   &
    LONG_NAME  = 'primary_production_phaeocystis',            &
    UNITS      = 'mg C m-2 d-1',                              &
    DIMS       = MAPL_DimsHorzOnly,                           &
    VLOCATION  = MAPL_VLocationNone,                          &
    RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
    SHORT_NAME = 'DISCHARGE',                                   &
    LONG_NAME  = 'river_discharge_at_ocean_points',            &
    UNITS      = 'kg m-2 s-1',                              &
    DIMS       = MAPL_DimsHorzOnly,                           &
    VLOCATION  = MAPL_VLocationNone,                          &
    RC=STATUS  )
    VERIFY_(STATUS)


!  !IMPORT STATE:

    call MAPL_AddImportSpec(GC,                               &
    SHORT_NAME = 'TIRRQ',                                     &
    LONG_NAME  = 'total_irradiance in quanta',                &
    UNITS      = 'umol quanta m-2 s-1',                       &
    DIMS       = MAPL_DimsHorzVert,                           &
    VLOCATION  = MAPL_VLocationCenter,                        &
    RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                               &
    SHORT_NAME              = 'CDOMABSQ',                     &
    LONG_NAME               = 'quanta_absorbed_by_CDOM',      &
    UNITS                   = 'umol quanta k-2 s-1',          &
    DIMS                    = MAPL_DimsHorzVert,              &
    VLOCATION               = MAPL_VLocationCenter,           &
    RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                               &
    SHORT_NAME = 'PS',                                        &
    LONG_NAME  = 'surface_pressure',                          &
    UNITS      = 'Pa',                                        &
    DIMS       = MAPL_DimsHorzOnly,                           &
    VLOCATION  = MAPL_VLocationNone,                          &
    RESTART    = MAPL_RestartSkip,                            &
    RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                               &
    SHORT_NAME = 'UU',                                        &
    LONG_NAME  = 'surface_wind_speed',                        &
    UNITS      = 'm s-1',                                     &
    DIMS       = MAPL_DimsHorzOnly,                           &
    VLOCATION  = MAPL_VLocationNone,                          &
    RESTART    = MAPL_RestartSkip,                            &
    RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                               &
    SHORT_NAME = 'DH',                                        &
    LONG_NAME  = 'layer_thickness',                           &
    UNITS      = 'm',                                         &
    DIMS       = MAPL_DimsHorzVert,                           &
    VLOCATION  = MAPL_VLocationCenter,                        &
!    RESTART    = MAPL_RestartSkip,                            &
    RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                               &
    SHORT_NAME = 'T',                                         &
    LONG_NAME  = 'potential_temperature',                     &
    UNITS      = 'K',                                         &
    DIMS       = MAPL_DimsHorzVert,                           &
    VLOCATION  = MAPL_VLocationCenter,                        &
!    RESTART    = MAPL_RestartSkip,                            &
    RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                               &
    SHORT_NAME = 'S',                                         &
    LONG_NAME  = 'Salinity',                                  &
    UNITS      = 'psu',                                       &
    DIMS       = MAPL_DimsHorzVert,                           &
    VLOCATION  = MAPL_VLocationCenter,                        &
!    RESTART    = MAPL_RestartSkip,                            &
    RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                               &
    SHORT_NAME = 'FRACICE',                                   &
    LONG_NAME  = 'fractional_cover_of_seaice',                &
    UNITS      = '1',                                         &
    DIMS       = MAPL_DimsHorzOnly,                           &
!    UNGRIDDED_DIMS = (/NUM_ICE_CATEGORIES/),                  &
    VLOCATION  = MAPL_VLocationNone,                          &
!    RESTART    = MAPL_RestartSkip,                            &
    RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                               &
    SHORT_NAME = 'FROCEAN',                                   &
    LONG_NAME  = 'fraction_of_gridbox_covered_by_skin',                &
    UNITS      = '1',                                         &
    DIMS       = MAPL_DimsHorzOnly,                           &
    VLOCATION  = MAPL_VLocationNone,                          &
    RESTART    = MAPL_RestartSkip,                            &
    RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                               &
    SHORT_NAME = 'MASKO',                                   &
    LONG_NAME  = 'ocean_mask',                &
    UNITS      = '1',                                         &
    DIMS       = MAPL_DimsHorzOnly,                           &
    VLOCATION  = MAPL_VLocationNone,                          &
!    RESTART    = MAPL_RestartSkip,                            &
    RC=STATUS  )
    VERIFY_(STATUS)


    call MAPL_AddImportSpec(GC,                               &
    SHORT_NAME              = 'DUDP',                         &
    LONG_NAME               = 'dry dust deposition',          &
    UNITS                   = 'kg m-2 s-1',                   &
    DIMS                    = MAPL_DimsHorzOnly,              &
    UNGRIDDED_DIMS          = (/NUM_DUDP/),                   &
    VLOCATION               = MAPL_VLocationNone,             &
    RESTART                 = MAPL_RestartSkip,               &
    RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                               &
    LONG_NAME               = 'Wet Dust Deposition',          &
    UNITS                   = 'kg m-2 s-1',                   &
    SHORT_NAME              = 'DUWT',                         &
    DIMS                    = MAPL_DimsHorzOnly,              &
    UNGRIDDED_DIMS          = (/NUM_DUWT/),                   &
    VLOCATION               = MAPL_VLocationNone,             &
    RESTART                 = MAPL_RestartSkip,               &
    RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                               &
    LONG_NAME               = 'Dust Sedimentation',           &
    UNITS                   = 'kg m-2 s-1',                   &
    SHORT_NAME              = 'DUSD',                         &
    DIMS                    = MAPL_DimsHorzOnly,              &
    UNGRIDDED_DIMS          = (/NUM_DUSD/),                   &
    VLOCATION               = MAPL_VLocationNone,             &
    RESTART                 = MAPL_RestartSkip,               &
    RC=STATUS  )
    VERIFY_(STATUS)

    if(DO_DATA_ATM4OCN==0) then
      call MAPL_AddImportSpec(GC,                               &
      LONG_NAME               = 'Black Carbon Dry Deposition',  &
      UNITS                   = 'kg m-2 s-1',                   &
      SHORT_NAME              = 'BCDP',                         &
      DIMS                    = MAPL_DimsHorzOnly,              &
      UNGRIDDED_DIMS          = (/NUM_BCDP/),                   &
      VLOCATION               = MAPL_VLocationNone,             &
      RESTART                 = MAPL_RestartSkip,               &
      RC=STATUS  )
      VERIFY_(STATUS)

      call MAPL_AddImportSpec(GC,                               &
      LONG_NAME               = 'Black Carbon Wet Deposition',  &
      UNITS                   = 'kg m-2 s-1',                   &
      SHORT_NAME              = 'BCWT',                         &
      DIMS                    = MAPL_DimsHorzOnly,              &
      UNGRIDDED_DIMS          = (/NUM_BCWT/),                   &
      VLOCATION               = MAPL_VLocationNone,             &
      RESTART                 = MAPL_RestartSkip,               &
      RC=STATUS  )
      VERIFY_(STATUS)

      call MAPL_AddImportSpec(GC,                               &
      LONG_NAME               = 'Organic Carbon Dry Deposition',&
      UNITS                   = 'kg m-2 s-1',                   &
      SHORT_NAME              = 'OCDP',                         &
      DIMS                    = MAPL_DimsHorzOnly,              &
      UNGRIDDED_DIMS          = (/NUM_OCDP/),                   &
      VLOCATION               = MAPL_VLocationNone,             &
      RESTART                 = MAPL_RestartSkip,               &
      RC=STATUS  )
      VERIFY_(STATUS)

      call MAPL_AddImportSpec(GC,                               &
      LONG_NAME               = 'Organic Carbon Wet Deposition',&
      UNITS                   = 'kg m-2 s-1',                   &
      SHORT_NAME              = 'OCWT',                         &
      DIMS                    = MAPL_DimsHorzOnly,              &
      UNGRIDDED_DIMS          = (/NUM_OCWT/),                   &
      VLOCATION               = MAPL_VLocationNone,             &
      RESTART                 = MAPL_RestartSkip,               &
      RC=STATUS  )
      VERIFY_(STATUS)
    endif

    call MAPL_AddImportSpec(GC,                               &
    SHORT_NAME              = 'CO2SC',                        &
    LONG_NAME               = 'atmospheric co2',              &
    UNITS                   = '1e-6',                         &
    DIMS                    = MAPL_DimsHorzOnly,              &
    VLOCATION               = MAPL_VLocationNone,             &
    RESTART                 = MAPL_RestartSkip,               &
    RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                               &
    SHORT_NAME              = 'DISCHARGE',                     &
    LONG_NAME               = 'river_discharge_at_ocean_points',      &
    UNITS                   = 'kg m-2 s-1',          &
    DIMS                    = MAPL_DimsHorzOnly,              &
    VLOCATION               = MAPL_VLocationCenter,           &
    RESTART                 = MAPL_RestartSkip,               &
    RC=STATUS )
    VERIFY_(STATUS)

! !INTERNAL STATE:

    call MAPL_AddInternalSpec(GC,                             &
    SHORT_NAME = 'NITRATE',                                   &
    LONG_NAME  = 'nitrate_concentration',                     &
    UNITS      = 'uM',                                        &
    DIMS       = MAPL_DimsHorzVert,                           &
    VLOCATION  = MAPL_VLocationCenter,                        &
    DEFAULT    = 5.0,                                         &
    FRIENDLYTO = 'OCEAN:OANA',                                     &
    RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                             &
    SHORT_NAME = 'AMMON',                                     &
    LONG_NAME  = 'ammonium_concentration',                    &
    UNITS      = 'uM',                                        &
    DIMS       = MAPL_DimsHorzVert,                           &
    VLOCATION  = MAPL_VLocationCenter,                        &
    DEFAULT    = 0.05,                                        &
    FRIENDLYTO = 'OCEAN:OANA',                                     &
    RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                             &
    SHORT_NAME = 'SILICA',                                    &
    LONG_NAME  = 'silica_concentration',                      &
    UNITS      = 'uM',                                        &
    DIMS       = MAPL_DimsHorzVert,                           &
    VLOCATION  = MAPL_VLocationCenter,                        &
    DEFAULT    = 5.0,                                         &
    FRIENDLYTO = 'OCEAN:OANA',                                     &
    RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                             &
    SHORT_NAME = 'IRON',                                      &
    LONG_NAME  = 'dissolved_iron_concentration',              &
    UNITS      = 'nM',                                        &
    DIMS       = MAPL_DimsHorzVert,                           &
    VLOCATION  = MAPL_VLocationCenter,                        &
    DEFAULT    = 0.2,                                         &
    FRIENDLYTO = 'OCEAN:OANA',                                     &
    RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                             &
    SHORT_NAME = 'DIATOM',                                    &
    LONG_NAME  = 'diatom_concentration',                      &
    UNITS      = 'mg m-3',                                    &
    DIMS       = MAPL_DimsHorzVert,                           &
    VLOCATION  = MAPL_VLocationCenter,                        &
    DEFAULT    = 0.05,                                        &
    FRIENDLYTO = 'OCEAN',                                     &
    RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                             &
    SHORT_NAME = 'CHLORO',                                    &
    LONG_NAME  = 'chlorophyte_concentration',                 &
    UNITS      = 'mg m-3',                                    &
    DIMS       = MAPL_DimsHorzVert,                           &
    VLOCATION  = MAPL_VLocationCenter,                        &
    DEFAULT    = 0.05,                                        &
    FRIENDLYTO = 'OCEAN',                                     &
    RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                             &
    SHORT_NAME = 'CYANO',                                     &
    LONG_NAME  = 'cyanobacteria_concentration',               &
    UNITS      = 'mg m-3',                                    &
    DIMS       = MAPL_DimsHorzVert,                           &
    VLOCATION  = MAPL_VLocationCenter,                        &
    DEFAULT    = 0.05,                                        &
    FRIENDLYTO = 'OCEAN',                                     &
    RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                             &
    SHORT_NAME = 'COCCO',                                     &
    LONG_NAME  = 'coccolithophore_concentration',             &
    UNITS      = 'mg m-3',                                    &
    DIMS       = MAPL_DimsHorzVert,                           &
    VLOCATION  = MAPL_VLocationCenter,                        &
    DEFAULT    = 0.05,                                        &
    FRIENDLYTO = 'OCEAN',                                     &
    RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                             &
    SHORT_NAME = 'DINO',                                      &
    LONG_NAME  = 'dinoflagellate_concentration',              &
    UNITS      = 'mg m-3',                                    &
    DIMS       = MAPL_DimsHorzVert,                           &
    VLOCATION  = MAPL_VLocationCenter,                        &
    DEFAULT    = 0.05,                                        &
    FRIENDLYTO = 'OCEAN',                                     &
    RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                             &
    SHORT_NAME = 'PHAEO',                                     &
    LONG_NAME  = 'phaeocystis_concentration',                 &
    UNITS      = 'mg m-3',                                    &
    DIMS       = MAPL_DimsHorzVert,                           &
    VLOCATION  = MAPL_VLocationCenter,                        &
    DEFAULT    = 0.05,                                        &
    FRIENDLYTO = 'OCEAN',                                     &
    RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                             &
    SHORT_NAME = 'HERB',                                      &
    LONG_NAME  = 'herbivore_concentration',                   &
    UNITS      = 'mg m-3',                                    &
    DIMS       = MAPL_DimsHorzVert,                           &
    VLOCATION  = MAPL_VLocationCenter,                        &
    DEFAULT    = 0.05,                                        &
    FRIENDLYTO = 'OCEAN',                                     &
    RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                             &
    SHORT_NAME = 'CDET',                                      &
    LONG_NAME  = 'carbon/nitrogen_detritus_concentration',    &
    UNITS      = 'ugC l-1',                                   &
    DIMS       = MAPL_DimsHorzVert,                           &
    VLOCATION  = MAPL_VLocationCenter,                        &
    DEFAULT    = 0.0,                                         &
    FRIENDLYTO = 'OCEAN',                                     &
    RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                             &
    SHORT_NAME = 'SDET',                                      &
    LONG_NAME  = 'silica_detritus_concentration',             &
    UNITS      = 'uM',                                        &
    DIMS       = MAPL_DimsHorzVert,                           &
    VLOCATION  = MAPL_VLocationCenter,                        &
    DEFAULT    = 0.0,                                         &
    FRIENDLYTO = 'OCEAN',                                     &
    RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                             &
    SHORT_NAME = 'FDET',                                      &
    LONG_NAME  = 'iron_detritus_concentration',               &
    UNITS      = 'nM',                                        &
    DIMS       = MAPL_DimsHorzVert,                           &
    VLOCATION  = MAPL_VLocationCenter,                        &
    DEFAULT    = 0.0,                                         &
    FRIENDLYTO = 'OCEAN',                                     &
    RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                             &
    SHORT_NAME = 'DOC',                                       &
    LONG_NAME  = 'dissolved_organic_carbon_concentration',    &
    UNITS      = 'uM',                                        &
    DIMS       = MAPL_DimsHorzVert,                           &
    VLOCATION  = MAPL_VLocationCenter,                        &
    DEFAULT    = 0.0,                                         &
    FRIENDLYTO = 'OCEAN',                                     &
    RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                             &
    SHORT_NAME = 'DIC',                                       &
    LONG_NAME  = 'dissolved_inorganic_carbon_concentration',  &
    UNITS      = 'uM',                                        &
    DIMS       = MAPL_DimsHorzVert,                           &
    VLOCATION  = MAPL_VLocationCenter,                        &
    DEFAULT    = 2055.0,                                      &
    FRIENDLYTO = 'OCEAN',                                     &
    RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                             &
    SHORT_NAME = 'ALK',                                       &
    LONG_NAME  = 'Alkalinity',                                &
    UNITS      = 'uEq l-1',                                   &
    DIMS       = MAPL_DimsHorzVert,                           &
    VLOCATION  = MAPL_VLocationCenter,                        &
    DEFAULT    = 2055.0,                                      &
    FRIENDLYTO = 'OCEAN',                                     &
    RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                             &
    SHORT_NAME = 'PIC',                                       &
    LONG_NAME  = 'Particulate inorganic carbon',              &
    UNITS      = 'ug l-1',                                    &
    DIMS       = MAPL_DimsHorzVert,                           &
    VLOCATION  = MAPL_VLocationCenter,                        &
    DEFAULT    = 0.0,                                         &
    FRIENDLYTO = 'OCEAN',                                     &
    RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                             &
    SHORT_NAME = 'CDC',                                       &
    LONG_NAME  = 'Colored dissolved organic carbon',          &
    UNITS      = 'uM',                                        &
    DIMS       = MAPL_DimsHorzVert,                           &
    VLOCATION  = MAPL_VLocationCenter,                        &
    DEFAULT    = 0.0,                                         &
    FRIENDLYTO = 'OCEAN',                                     &
    RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                             &
    SHORT_NAME = 'COCCOGRO',                                  &
    LONG_NAME  = 'maximum_growth_rate_of_coccos',             &
    UNITS      = 's-1',                                       &
    DIMS       = MAPL_DimsHorzVert,                           &
    VLOCATION  = MAPL_VLocationCenter,                        &
    DEFAULT    = 0.00000579,                                  &
    RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                             &
    SHORT_NAME = 'COCCOSIN',                                  &
    LONG_NAME  = 'maximum_sinking_rate_of_coccos',            &
    UNITS      = 'm s-1',                                     &
    DIMS       = MAPL_DimsHorzVert,                           &
    VLOCATION  = MAPL_VLocationCenter,                        &
    DEFAULT    = 0.00000347,                                  &
    RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                             &
    SHORT_NAME = 'AVGQ',                                      &
    LONG_NAME  = 'average_quantum_irradiance',                &
    UNITS      = 'umol m-2 s-1',                              &
    DIMS       = MAPL_DimsHorzVert,                           &
    VLOCATION  = MAPL_VLocationCenter,                        &
    DEFAULT    = 100.0,                                       &
    FRIENDLYTO = 'ORAD',                                      &
    RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                             &
    SHORT_NAME = 'CCHLRATIO',                                 &
    LONG_NAME  = 'carbon-to-chlorophyll ratio',               &
    UNITS      = 'g to g-1',                                  &
    DIMS       = MAPL_DimsHorzVert,                           &
    VLOCATION  = MAPL_VLocationCenter,                        &
    DEFAULT    = 50.0,                                        &
    ADD2EXPORT = .TRUE.,                                      &
    RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                             &
    SHORT_NAME = 'RIKDIA',                                    &
    LONG_NAME  = 'daily light saturation parameter diatoms',  &
    UNITS      = 'umol m-2 s-1',                              &
    DIMS       = MAPL_DimsHorzVert,                           &
    VLOCATION  = MAPL_VLocationCenter,                        &
    DEFAULT    = 25.0,                                        &
    RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                             &
    SHORT_NAME = 'RIKCHL',                                    &
    LONG_NAME  = 'daily light saturation parameter chloro',   &
    UNITS      = 'umol m-2 s-1',                              &
    DIMS       = MAPL_DimsHorzVert,                           &
    VLOCATION  = MAPL_VLocationCenter,                        &
    DEFAULT    = 25.0,                                        &
    RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                             &
    SHORT_NAME = 'RIKCYA',                                    &
    LONG_NAME  = 'daily light saturation parameter cyano',    &
    UNITS      = 'umol m-2 s-1',                              &
    DIMS       = MAPL_DimsHorzVert,                           &
    VLOCATION  = MAPL_VLocationCenter,                        &
    DEFAULT    = 25.0,                                        &
    RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                             &
    SHORT_NAME = 'RIKCOC',                                    &
    LONG_NAME  = 'daily light saturation parameter cocco',    &
    UNITS      = 'umol m-2 s-1',                              &
    DIMS       = MAPL_DimsHorzVert,                           &
    VLOCATION  = MAPL_VLocationCenter,                        &
    DEFAULT    = 25.0,                                        &
    RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                             &
    SHORT_NAME = 'RIKDIN',                                    &
    LONG_NAME  = 'daily light saturation parameter dino',     &
    UNITS      = 'umol m-2 s-1',                              &
    DIMS       = MAPL_DimsHorzVert,                           &
    VLOCATION  = MAPL_VLocationCenter,                        &
    DEFAULT    = 25.0,                                        &
    RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                             &
    SHORT_NAME = 'RIKPHA',                                    &
    LONG_NAME  = 'daily light saturation parameter phaeo',    &
    UNITS      = 'umol m-2 s-1',                              &
    DIMS       = MAPL_DimsHorzVert,                           &
    VLOCATION  = MAPL_VLocationCenter,                        &
    DEFAULT    = 25.0,                                        &
    RC=STATUS  )
    VERIFY_(STATUS)



!!EOC

! Set the Profiling timers
! ------------------------

    call MAPL_TimerAdd(GC, NAME="INITIALIZE", RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC, NAME="RUN",        RC=STATUS)
    VERIFY_(STATUS)
  
! Set generic init and final methods
! ----------------------------------

    call MAPL_GenericSetServices    ( GC, RC=STATUS)
    VERIFY_(STATUS)

!  All done
!-----------

    RETURN_(ESMF_SUCCESS)
  
  end subroutine SetServices

!-------------------------------------------------------------------
!BOPI

! !IROUTINE: Initialize -- Initialize method

! !INTERFACE:

  subroutine Initialize ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:
    type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component
    type(ESMF_State),    intent(inout) :: IMPORT ! Import state
    type(ESMF_State),    intent(inout) :: EXPORT ! Export state
    type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
    integer, optional,   intent(  out) :: RC     ! Error code:

! ! DESCRIPTION: This function initializes the OBIO gridded component.

!EOPI
!-------------------------------------------------------------------------

! Local Variables

    character(len=ESMF_MAXSTR)          :: IAm
    integer                             :: STATUS
    character(len=ESMF_MAXSTR)          :: COMP_NAME

    type (MAPL_MetaComp),   pointer     :: MAPL => null()
    type (T_OBIO_STATE),    pointer     :: State => null()
    type (OBIO_wrap)                    :: WRAP

! Begin
!------

    Iam = "Initialize"
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
    Iam = trim(COMP_NAME) // Iam

! Do Generic Initialize first
!----------------------------

    call MAPL_GenericInitialize( GC, IMPORT, EXPORT, CLOCK, RC=STATUS )

! Get my internal MAPL object
!----------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS )

! Start Total timer
!------------------

    call MAPL_TimerOn(MAPL,"TOTAL")
    call MAPL_TimerOn(MAPL,"INITIALIZE" )

! Allocate this instance of the internal state and put it in wrapper.
! -------------------------------------------------------------------

    allocate( State, __STAT__ )
    WRAP%PTR => State

! Save pointer to the wrapped internal state in the GC
! ----------------------------------------------------

    call ESMF_UserCompSetInternalState ( GC, 'OBIO_state', WRAP, STATUS ); VERIFY_(STATUS)

! Initialize biology parameters in private state
!-----------------------------------------------
#include "definetab.h"
    call setbio(State%rmumax,State%cnratio,State%cfratio,State%solFe,  &
                State%drate, State%Rm,     State%rlamz,  State%greff,  &
                State%dratez1,State%dratez2,State%regen, State%remin,  &
                State%rkn,   State%rks,    State%rkf,    State%rik,    &
                State%cchl,  State%wss,    State%axs,                  & 
                State%Pdeep)

! Stop Total timer
!-----------------

    call MAPL_TimerOff(MAPL,"INITIALIZE" )
    call MAPL_TimerOff(MAPL,"TOTAL")

!   end if

! All Done
!---------

    RETURN_(ESMF_SUCCESS)

  end subroutine Initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!BOP
! ! IROUTINE: RUN -- First Run stage for the Obio component

! !INTERFACE:

  subroutine RUN ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
    type(ESMF_State),    intent(inout) :: IMPORT ! Import state
    type(ESMF_State),    intent(inout) :: EXPORT ! Export state
    type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
    integer, optional,   intent(  out) :: RC     ! Error code:

! ! DESCRIPTION: Ocean biology

!EOP


! ErrLog Variables

    character(len=ESMF_MAXSTR)       :: IAm
    integer                          :: STATUS
    character(len=ESMF_MAXSTR)       :: COMP_NAME

! Locals

    type (MAPL_MetaComp), pointer    :: MAPL => null()
    type (T_OBIO_STATE),  pointer    :: State => null()
    type (OBIO_wrap)                 :: WRAP
!    type (ESMF_Bundle)              :: INTERNAL
    type (ESMF_State)                :: INTERNAL
    type (ESMF_Time)                 :: CURRENTTIME
    character(len=ESMF_MAXSTR)       :: DATAFILE
    type (ESMF_TimeInterval)         :: DTint
    type (MAPL_SunOrbit)             :: ORBIT

    real, dimension(:,:), allocatable :: BIO
    real, dimension(:),   allocatable :: AVGQ1
    real, dimension(:),   allocatable :: CCHLRATIO1
    real, dimension(:),   allocatable :: GCMAX
    real, dimension(:),   allocatable :: WSSC
    real, dimension(:),   allocatable :: PPZ
    real, dimension(:,:), allocatable :: RIKD
    real, dimension(:,:), allocatable :: COSZ
    real, dimension(:,:), allocatable :: SLR
    real, dimension(:,:), allocatable :: WGHT
    real, dimension(:,:), allocatable :: DISCHRG

    real, pointer, dimension(:,:) :: LONS => null()
    real, pointer, dimension(:,:) :: LATS => null()

    logical :: IS_MIDNIGHT
    real    :: PCO
    real    :: CFLX
    real    :: DT
    real*8  :: DT8
#include "definebio.h"
    integer :: IM,JM,LM,i,j,k,n
    integer :: YY,DOY,MM,DD,hr,mn,sec
    real    :: slp,wspd,Fedep,Fe,fnoice
    real    :: aco2
    real    :: discharg
    real, save :: CO2_0
    data    CO2_0 /0.0/


! pointers to import

    real, pointer, dimension(:,:)   :: UU => null()
    real, pointer, dimension(:,:)   :: PS => null()
    real, pointer, dimension(:,:,:) :: DH => null()
    real, pointer, dimension(:,:,:) :: T => null()
    real, pointer, dimension(:,:,:) :: S => null()
    real, pointer, dimension(:,:)   :: MASKO => null()
    real, pointer, dimension(:,:,:) :: TIRRQ => null()
    real, pointer, dimension(:,:,:) :: CDOMABSQ => null()
!    real, pointer, dimension(:,:,:) :: FRACICE => null()
    real, pointer, dimension(:,:)   :: FR => null()
    real, pointer, dimension(:,:,:)   :: DRY_DUST => null()
    real, pointer, dimension(:,:,:)   :: WET_DUST => null()
    real, pointer, dimension(:,:,:)   :: SED_DUST => null()
    real, pointer, dimension(:,:)   :: CO2SC => null()
    real, pointer, dimension(:,:)   :: DISCHARGE => null()


! pointers to export

    real, pointer, dimension(:,:)   :: PCO2 => null()
    real, pointer, dimension(:,:)   :: FCO2 => null()
    real, pointer, dimension(:,:)   :: ppDIATOM => null()
    real, pointer, dimension(:,:)   :: ppCHLORO => null()
    real, pointer, dimension(:,:)   :: ppCYANO => null()
    real, pointer, dimension(:,:)   :: ppCOCCO => null()
    real, pointer, dimension(:,:)   :: ppDINO => null()
    real, pointer, dimension(:,:)   :: ppPHAEO => null()
    real, pointer, dimension(:,:)   :: dout => null()


! pointers to internal

    real, pointer, dimension(:,:,:) :: NITRATE => null()
    real, pointer, dimension(:,:,:) :: AMMON => null()
    real, pointer, dimension(:,:,:) :: SILICA => null()
    real, pointer, dimension(:,:,:) :: IRON => null()
    real, pointer, dimension(:,:,:) :: DIATOM => null()
    real, pointer, dimension(:,:,:) :: CHLORO => null()
    real, pointer, dimension(:,:,:) :: CYANO => null()
    real, pointer, dimension(:,:,:) :: COCCO => null()
    real, pointer, dimension(:,:,:) :: DINO => null()
    real, pointer, dimension(:,:,:) :: PHAEO => null()
    real, pointer, dimension(:,:,:) :: HERB => null()
    real, pointer, dimension(:,:,:) :: CDET => null()
    real, pointer, dimension(:,:,:) :: SDET => null()
    real, pointer, dimension(:,:,:) :: FDET => null()
    real, pointer, dimension(:,:,:) :: DOC => null()
    real, pointer, dimension(:,:,:) :: DIC => null()
    real, pointer, dimension(:,:,:) :: ALK => null()
    real, pointer, dimension(:,:,:) :: PIC => null()
    real, pointer, dimension(:,:,:) :: CDC => null()
    real, pointer, dimension(:,:,:) :: COCCOGRO => null()
    real, pointer, dimension(:,:,:) :: COCCOSIN => null()
    real, pointer, dimension(:,:,:) :: AVGQ => null()
    real, pointer, dimension(:,:,:) :: CCHLRATIO => null()
    real, pointer, dimension(:,:,:) :: RIKDIA => null()
    real, pointer, dimension(:,:,:) :: RIKCHL => null()
    real, pointer, dimension(:,:,:) :: RIKCYA => null()
    real, pointer, dimension(:,:,:) :: RIKCOC => null()
    real, pointer, dimension(:,:,:) :: RIKDIN => null()
    real, pointer, dimension(:,:,:) :: RIKPHA => null()

    real, pointer, dimension(:,:)   :: FRICE => null()

!=============================================================================

! Begin... 

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // "Run"

! Get my internal MAPL_Generic state
!-----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)

! Get parameters from MAPL object
!--------------------------------

    call MAPL_Get(MAPL,     &
         LATS      = LATS,  &
         LONS      = LONS,  &
         ORBIT     = ORBIT, &
         RC=STATUS )

! Start Total timer
!------------------

    call MAPL_TimerOn(MAPL,"TOTAL")
    call MAPL_TimerOn(MAPL,"RUN" )

! Get parameters from MAPL object
!--------------------------------

    call MAPL_Get(MAPL,                  &
         IM                  = IM,       &
         JM                  = JM,       &
         LM                  = LM,       &
         INTERNAL_ESMF_STATE = INTERNAL, &
         RC=STATUS )
    VERIFY_(STATUS)

! Get my internal private state
!------------------------------

    call ESMF_UserCompGetInternalState(GC, 'OBIO_state', WRAP, STATUS); VERIFY_(STATUS)
    State => WRAP%PTR

! Get Time step
!--------------

    call ESMF_ClockGet       (CLOCK, TimeStep=DTint, RC=STATUS)
    VERIFY_(STATUS)
    call ESMF_TimeIntervalGet(DTint, s_r8=DT8,       RC=STATUS)
    VERIFY_(STATUS)
    DT = DT8

! Pointers to Exports
!--------------------

    call MAPL_GetPointer(EXPORT, PCO2,     'PCO2',     RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, FCO2,     'FCO2',     RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, ppDIATOM, 'PPDIATOM', RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, ppCHLORO, 'PPCHLORO', RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, ppCYANO,  'PPCYANO',  RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, ppCOCCO,  'PPCOCCO',  RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, ppDINO,  'PPDINO',  RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, ppPHAEO,  'PPPHAEO',  RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, DOUT,  'DISCHARGE',  RC=STATUS)
    VERIFY_(STATUS)


! Pointers to Imports
!--------------------

    call MAPL_GetPointer(IMPORT, TIRRQ,   'TIRRQ',   RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, CDOMABSQ, 'CDOMABSQ',   RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, UU,      'UU',      RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, PS,      'PS',      RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, DH,      'DH',      RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, T,       'T',       RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, S,       'S',       RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, MASKO,   'MASKO',   RC=STATUS)
    VERIFY_(STATUS)
!    call MAPL_GetPointer(IMPORT, FRACICE, 'FRACICE', RC=STATUS)
!    VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, FRICE, 'FRACICE', RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, FR, 'FROCEAN', RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, DRY_DUST,'DUDP',    RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, WET_DUST,'DUWT',    RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, SED_DUST,'DUSD',    RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, CO2SC,   'CO2SC',   RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, DISCHARGE,   'DISCHARGE',   RC=STATUS)
    VERIFY_(STATUS)


! Pointers to Internals
!----------------------

    call MAPL_GetPointer(INTERNAL, NITRATE,  'NITRATE',  RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, AMMON,    'AMMON',    RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, SILICA,   'SILICA',   RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, IRON,     'IRON',     RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, DIATOM,   'DIATOM',   RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, CHLORO,   'CHLORO',   RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, CYANO,    'CYANO',    RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, COCCO,    'COCCO',    RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, DINO,     'DINO',     RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, PHAEO,    'PHAEO',    RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, HERB,     'HERB',     RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, CDET,     'CDET',     RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, SDET,     'SDET',     RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, FDET,     'FDET',     RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, DOC,      'DOC',      RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, DIC,      'DIC',      RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, ALK,      'ALK',      RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, PIC,      'PIC',      RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, CDC,      'CDC',      RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, COCCOGRO, 'COCCOGRO', RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, COCCOSIN, 'COCCOSIN', RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, AVGQ,     'AVGQ',     RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, CCHLRATIO, 'CCHLRATIO', RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, RIKDIA,   'RIKDIA',   RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, RIKCHL,   'RIKCHL',   RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, RIKCYA,   'RIKCYA',   RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, RIKCOC,   'RIKCOC',   RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, RIKDIN,   'RIKDIN',   RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, RIKPHA,   'RIKPHA',   RC=STATUS)
    VERIFY_(STATUS)

! Weight for ocean grid
!----------------------
       allocate(WGHT(IM,JM), __STAT__)
       allocate(DISCHRG(IM,JM), __STAT__)
       where(MASKO>0)
          WGHT = FR/MASKO
       elsewhere
          WGHT = 0.0
       end where

       where(WGHT > 0.0)
          dischrg=discharge
       elsewhere
          dischrg=0.0
       end where

       if(associated(dout)) then
          where(masko > 0.0)
             dout=dischrg
          elsewhere
             dout=mapl_undef
          end where
       end if

! Get solar zenith angle
!-----------------------

    allocate(COSZ(IM,JM),__STAT__)
    allocate( SLR(IM,JM),__STAT__)

    call MAPL_SunGetInsolation(LONS, LATS, &
               ORBIT, COSZ, SLR,           &
               CLOCK = CLOCK,              &
               RC=STATUS )
    VERIFY_(STATUS)

! Get parameters from configuration
!----------------------------------

    call ESMF_ClockGet(CLOCK,       currTime=CURRENTTIME, RC=STATUS)
    VERIFY_(STATUS)

    call ESMF_TimeGet (CURRENTTIME, YY=YY, DayOfYear=DOY,   &
          mm=MM, dd=DD, h=hr, m=mn, s=sec, RC=STATUS)
    VERIFY_(STATUS)

! Alarm for midnight
!-------------------

    if (hr==0 .and. mn==0 .and. sec==0) then
      IS_MIDNIGHT = .true.
    else
      IS_MIDNIGHT = .false.
    endif

!    aco2=368.6
    aco2 = CO2SC(im/2,jm/2)*1.0E6
    discharg = DISCHARGE(im,jm)
    if ( CO2_0 /= aco2 ) then
     CO2_0 = aco2
     if(mapl_am_i_root()) write(6,'(A,I4,A,I3,A,E12.5)') &
     " Updated CO2 in OBIO for year/day ", YY, "/", DOY, " is ", aco2
    endif

! Get total seaice fraction
!--------------------------
!    allocate(FRICE(IM,JM), __STAT__)
!    FRICE = SUM(FRACICE, DIM=3)

!  FRACICE from Import state is already in 2d. These two lines of code must have been 
!  before FRACICE was a 2d variable. EMS


! Allocate BIO and other variables
!---------------------------------

    allocate(BIO(LM,ntyp),                   __STAT__)
    allocate(GCMAX(LM), WSSC(LM), AVGQ1(LM), __STAT__)
    allocate(CCHLRATIO1(LM),                 __STAT__)
    allocate(RIKD(LM,nchl),                  __STAT__)
    allocate(PPZ(nchl),                      __STAT__)

    do j = 1, JM
     do i = 1, IM
      if (MASKO(i,j) > 0.0 ) then
       slp = PS(i,j)*0.01       ! convert from Pa to mbar
       wspd = UU(i,j)
       fnoice = 1.0 - FRICE(i,j)
       aco2 = CO2SC(i,j)
!       aco2 = aco2*1.0E6
       discharg = DISCHARGE(i,j)

!    Dust is assumed to have an iron fraction of 5% for clay and 1.2%
!    for silt (Fung et al. 2000)
       Fedep = 0.05*(DRY_DUST(i,j,1)+WET_DUST(i,j,1)+SED_DUST(i,j,1))  &
              + 0.012*(DRY_DUST(i,j,2)+WET_DUST(i,j,2)+SED_DUST(i,j,2) &
              + DRY_DUST(i,j,2)+WET_DUST(i,j,2)+SED_DUST(i,j,2)        &
              + DRY_DUST(i,j,3)+WET_DUST(i,j,3)+SED_DUST(i,j,3)        &
              + DRY_DUST(i,j,4)+WET_DUST(i,j,4)+SED_DUST(i,j,4)        &
              + DRY_DUST(i,j,5)+WET_DUST(i,j,5)+SED_DUST(i,j,5))
!    Convert Fe from kg m-2 s-1 to nmol Fe m-2 s-1
       Fe = (1.0E12/MOLWGHT_FE)*Fedep
!#ifdef DATAATM
!!    For data atmosphere, use pre-calculated fractions and content
!!    and use only the first fraction
!      Fe = DRY_CLAY(i,j)/3600.0   ! convert from nmFe m-2 h-1 to s-1
!#endif

!
!  There are mismatches between the ocean, land and atmosphere in GEOS-5
!  and live ocean points for MOM that do not have a corresponding 
!  atmosphere.  These are called "grottoes" because they are assumed 
!  to have ocean underneath with land overhead.  Use defaults for 
!  atmosphere that represent this condition.
       if (slp < 100.0 .or. slp > 1.0E10) then
        slp = 1013.25
        wspd = 0.0
!        Fedep = 0.05*5.0E-13 + 0.012*2.2E-13 ! Reference Fe is mean
        Fe = 0.0
        aco2 = 368.6
       endif
       if (fnoice < 0.0 .or. fnoice > 1.0)fnoice = 1.0

!if (aco2 > 500.0)then
!write(6,*)'i,j,co2sc,aco2 = ',i,j,CO2SC(i,j),aco2
!endif

!if (discharg > 0.0)then
!write(6,*)'i,j,discharg = ',i,j,discharg
!endif

       BIO(:,1)   = NITRATE(i,j,:)
       BIO(:,2)   = AMMON(i,j,:)
       BIO(:,3)   = SILICA(i,j,:)
       BIO(:,4)   = IRON(i,j,:)
       BIO(:,5)   = DIATOM(i,j,:)
       BIO(:,6)   = CHLORO(i,j,:)
       BIO(:,7)   = CYANO(i,j,:)
       BIO(:,8)   = COCCO(i,j,:)
       BIO(:,9)   = DINO(i,j,:)
       BIO(:,10)  = PHAEO(i,j,:)
       BIO(:,11)  = HERB(i,j,:)
       BIO(:,12)  = CDET(i,j,:)
       BIO(:,13)  = SDET(i,j,:)
       BIO(:,14)  = FDET(i,j,:)
       BIO(:,15)  = DOC(i,j,:)
       BIO(:,16)  = DIC(i,j,:)
       BIO(:,17)  = ALK(i,j,:)
       BIO(:,18)  = PIC(i,j,:)
       BIO(:,19)  = CDC(i,j,:)
       GCMAX(:)   = COCCOGRO(i,j,:)
       WSSC(:)    = COCCOSIN(i,j,:)
       AVGQ1(:)   = AVGQ(i,j,:)
       CCHLRATIO1(:) = CCHLRATIO(i,j,:)
       RIKD(:,1)  = RIKDIA(i,j,:)
       RIKD(:,2)  = RIKCHL(i,j,:)
       RIKD(:,3)  = RIKCYA(i,j,:)
       RIKD(:,4)  = RIKCOC(i,j,:)
       RIKD(:,5)  = RIKDIN(i,j,:)
       RIKD(:,6)  = RIKPHA(i,j,:)
       if ( IS_MIDNIGHT ) then
         call daysetbio( LM, State%rik, State%cchl, CCHLRATIO1,     &
                         AVGQ1, GCMAX, BIO, DH(i,j,:), RIKD, WSSC)
       endif

       call kloop ( LM, DT,        State%solFe,      State%remin,   &
                    State%rkn,     State%rks,        State%rkf,     &
                    State%cnratio, State%cfratio,    CCHLRATIO1,    &
                    State%wss,     WSSC,             State%drate,   &
                    State%Rm,      State%rlamz,      State%greff,   &
                    State%dratez1, State%dratez2,    State%regen,   &
                    State%axs,     State%rmumax,                    &
                    BIO, DH(i,j,:), Fe, fnoice, RIKD, GCMAX,        &
                    TIRRQ(i,j,:),CDOMABSQ(i,j,:), aco2, wspd, slp,  &
                    T(i,j,:)-MAPL_TICE, S(i,j,:), PCO, CFLX, PPZ)


BIO = max(BIO,0.0)    !reduce MOM4 propensity for negative values
       NITRATE(i,j,:)  = BIO(:,1)
       AMMON(i,j,:)    = BIO(:,2) 
       SILICA(i,j,:)   = BIO(:,3) 
       IRON(i,j,:)     = BIO(:,4) 
       DIATOM(i,j,:)   = BIO(:,5) 
       CHLORO(i,j,:)   = BIO(:,6) 
       CYANO(i,j,:)    = BIO(:,7) 
       COCCO(i,j,:)    = BIO(:,8) 
       DINO(i,j,:)     = BIO(:,9)
       PHAEO(i,j,:)    = BIO(:,10)
       HERB(i,j,:)     = BIO(:,11)
       CDET(i,j,:)     = BIO(:,12)
       SDET(i,j,:)     = BIO(:,13)
       FDET(i,j,:)     = BIO(:,14)
       DOC(i,j,:)      = BIO(:,15)
       DIC(i,j,:)      = BIO(:,16)
       ALK(i,j,:)      = BIO(:,17)
       PIC(i,j,:)      = BIO(:,18)
       CDC(i,j,:)      = BIO(:,19)
       COCCOGRO(i,j,:) = GCMAX(:)
       COCCOSIN(i,j,:) = WSSC(:)
       AVGQ(i,j,:)     = AVGQ1(:)
       CCHLRATIO(i,j,:) = CCHLRATIO1(:)
       RIKDIA(i,j,:)   = RIKD(:,1)
       RIKCHL(i,j,:)   = RIKD(:,2)
       RIKCYA(i,j,:)   = RIKD(:,3)
       RIKCOC(i,j,:)   = RIKD(:,4)
       RIKDIN(i,j,:)   = RIKD(:,5)
       RIKPHA(i,j,:)   = RIKD(:,6)
       if ( associated(ppDIATOM) ) ppDIATOM(i,j) = PPZ(1)
       if ( associated(ppCHLORO) ) ppCHLORO(i,j) = PPZ(2)
       if ( associated(ppCYANO) )  ppCYANO(i,j)  = PPZ(3)
       if ( associated(ppCOCCO) )  ppCOCCO(i,j)  = PPZ(4)
       if ( associated(ppDINO) )   ppDINO(i,j)   = PPZ(5)
       if ( associated(ppPHAEO) )  ppPHAEO(i,j)  = PPZ(6)
       if ( associated(PCO2) ) PCO2(i,j) = PCO
       if ( associated(FCO2) ) FCO2(i,j) = CFLX
      endif
     enddo
    enddo

    if ( associated(ppDIATOM) ) &
      where ( MASKO == 0.0 ) ppDIATOM = MAPL_UNDEF
    if ( associated(ppCHLORO) ) &
      where ( MASKO == 0.0 ) ppCHLORO = MAPL_UNDEF
    if ( associated(ppCYANO) ) &
      where ( MASKO == 0.0 ) ppCYANO  = MAPL_UNDEF
    if ( associated(ppCOCCO) ) &
      where ( MASKO == 0.0 ) ppCOCCO  = MAPL_UNDEF
    if ( associated(ppDINO) ) & 
      where ( MASKO == 0.0 ) ppDINO   = MAPL_UNDEF
    if ( associated(ppPHAEO) ) & 
      where ( MASKO == 0.0 ) ppPHAEO  = MAPL_UNDEF
    if ( associated(PCO2) ) &
      where ( MASKO == 0.0 ) PCO2 = MAPL_UNDEF
    if ( associated(FCO2) ) &
      where ( MASKO == 0.0 ) FCO2 = MAPL_UNDEF

    deallocate(COSZ  )
    deallocate(SLR   )
!    deallocate(FRICE)
    deallocate(BIO   )
    deallocate(GCMAX ) ; deallocate(WSSC  ) ; deallocate(AVGQ1 )
    deallocate(CCHLRATIO1 )
    deallocate(RIKD  )
    deallocate(PPZ   )
    deallocate(WGHT)
    deallocate(dischrg)

!  All done
!-----------
    call MAPL_TimerOff(MAPL,"RUN" )
    call MAPL_TimerOff(MAPL,"TOTAL")

    RETURN_(ESMF_SUCCESS)

  end subroutine RUN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module GEOS_OceanbiogeochemGridCompMod
