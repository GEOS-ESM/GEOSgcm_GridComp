
!  $Id$

#include "MAPL_Generic.h"

module GEOS_RadiationGridCompMod

!BOP

! !MODULE: GEOS_RadiationGridCompMod--Container for atmospheric radiation calculations

! !DESCRIPTION: A Composite MAPL/ESMF gridded component (GC) containing the
!    longwave and shortwave radiation GCs. It is intended as a container
!    for the ESMF/MAPL Solar and Irrad gridded components used in GEOS-5.
!    In SetServices, it creates one instance of each of its two children
!    (SOLAR and IRRAD) Its Run method combines results from the children to produce
!    total radiative exports. \newline
!
!  It follows the standard rules for composite ESMF/MAPL GCs.
!  It passes the ESMF grid that appears in the
!  gridded component to both children and all their Imports
!  and Exports are assumed to be on this grid. The grid must
!  be present in the GC and properly initialized before Initialize
!  is called, since {\tt GEOS\_RadiationGridCompMod} has its own 
!  explicitly declared Import state variables. \newline
!
!  The restrictions on the grid are that it be 3-dimensional
!  with two horizontal and one vertical dimension and
!  with only the horizontal dimensions decomposed. The vertical dimension
!  is also assumed to the the third dimension of the Fortran arrays and
!  is indexed from the top down. No particular vertical coordinate is assumed,
!  rather the 3-dimensional field of air pressure at the layer interfaces is 
!  a required Import. \newline
!
!  This module contains only SetServices, Initialize, and Run methods. 
!  The Finalize method
!  being defaulted to the MAPL\_Generic version. 
!  The SetServices method is the only public
!  entity. There are no public types or data. \newline
!  \newline
!
! !USES:

  use ESMF
  use MAPL

  use GEOS_SolarGridCompMod,  only : solarSetServices  => SetServices
  use GEOS_IrradGridCompMod,  only : irradSetServices  => SetServices
  use GEOS_SatsimGridCompMod, only : satsimSetServices => SetServices

  implicit none
  private

! !PUBLIC MEMBER FUNCTIONS:

  public SetServices

!EOP

  integer ::          SOL
  integer ::          IRR
  integer ::          STM

   contains

!BOP

! !IROUTINE: SetServices -- Sets ESMF services for this component

! !INTERFACE:

  subroutine SetServices ( GC, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    integer, optional                  :: RC  ! return code

! !DESCRIPTION:  
!
!  {\tt GEOS\_RadiationGridCompMod} registers Initialize and  Run methods with ESMF. 
!  The Finalize method
!  being defaulted to the MAPL\_Generic version and automatically
!  registered by MAPL\_GenericSetServices. 
!  SetServices registers the two children (SOLAR and IRRAD) with MAPL.
! \newline
!
!  Fields in its Import and Export states that are needed or produced by
!  the component are explicitly
!  registered in SetServices and are described in tables in this documentation.
!  The Import and Export States may also contain fields implicitly
!  inherited from the from the children.
!  Finally, a number of child exports are ``promoted'' to be
!  explicit Exports of {\tt GEOS\_RadiationGridCompMod}.  \newline
!
!  {\tt GEOS\_RadiationGridCompMod} has no Internal state.  \newline
!
! 

!EOP

! ErrLog Variables

    character(len=ESMF_MAXSTR)              :: IAm
    integer                                 :: STATUS
    character(len=ESMF_MAXSTR)              :: COMP_NAME
    integer                                 :: USE_SATSIM
    integer                                 :: USE_SATSIM_ISCCP
    integer                                 :: USE_SATSIM_MODIS
    integer                                 :: USE_SATSIM_RADAR
    integer                                 :: USE_SATSIM_LIDAR
    integer                                 :: USE_SATSIM_MISR

! Locals

    integer                                 :: I
    type (MAPL_MetaComp),      pointer      :: MAPL


! Begin...

! Get my name and set-up traceback handle
! ---------------------------------------

    call ESMF_GridCompGet ( GC, name=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // 'SetServices'


! Get my MAPL_Generic state
!--------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)

! Set the Run entry point
! -----------------------

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_INITIALIZE,  Initialize, RC=status )
    VERIFY_(STATUS)
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,  Run, RC=STATUS)
    VERIFY_(STATUS)

    SOL = MAPL_AddChild(GC, NAME='SOLAR', SS=solarSetServices, RC=STATUS)
    VERIFY_(STATUS)
    IRR = MAPL_AddChild(GC, NAME='IRRAD', SS=irradSetServices, RC=STATUS)
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
   
! set use_satsim if anything is toggled
 
   USE_SATSIM = USE_SATSIM + USE_SATSIM_ISCCP + USE_SATSIM_MODIS + USE_SATSIM_RADAR + USE_SATSIM_LIDAR + USE_SATSIM_MISR 

   if (USE_SATSIM > 0 ) then
       STM = MAPL_AddChild(GC, NAME='SATSIM', SS=satsimSetServices, RC=STATUS)
       VERIFY_(STATUS)
   end if    

! Set the state variable specs.
! -----------------------------

!BOS

! !IMPORT STATE:

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'PLEINST',                           &
        LONG_NAME          = 'air_pressure',                      &
        UNITS              = 'Pa',                                &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationEdge,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

! !EXPORT STATE:

    call MAPL_AddExportSpec ( GC,                                   &
         SHORT_NAME = 'DTDT',                                            &
         LONG_NAME  = 'pressure_weighted_air_temperature_tendency_due_to_radiation',&
         UNITS      = 'Pa K s-1',                                        &
         DIMS       = MAPL_DimsHorzVert,                                 &
         VLOCATION  = MAPL_VLocationCenter,                              &
                                                              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC,                                   &
         SHORT_NAME = 'RADLW',                                           &
         LONG_NAME  = 'air_temperature_tendency_due_to_longwave',        &
         UNITS      = 'K s-1',                                           &
         DIMS       = MAPL_DimsHorzVert,                                 &
         VLOCATION  = MAPL_VLocationCenter,                              &
                                                              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC,                                   &
         SHORT_NAME = 'RADSW',                                           &
         LONG_NAME  = 'air_temperature_tendency_due_to_shortwave',       &
         UNITS      = 'K s-1',                                           &
         DIMS       = MAPL_DimsHorzVert,                                 &
         VLOCATION  = MAPL_VLocationCenter,                              &
                                                              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC,                                   &
         SHORT_NAME = 'RADLWC',                                          &
         LONG_NAME  = 'air_temperature_tendency_due_to_longwave_for_clear_skies',&
         UNITS      = 'K s-1',                                           &
         DIMS       = MAPL_DimsHorzVert,                                 &
         VLOCATION  = MAPL_VLocationCenter,                              &
                                                              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC,                                   &
         SHORT_NAME = 'RADSWC',                                          &
         LONG_NAME  = 'air_temperature_tendency_due_to_shortwave_for_clear_skies',&
         UNITS      = 'K s-1',                                           &
         DIMS       = MAPL_DimsHorzVert,                                 &
         VLOCATION  = MAPL_VLocationCenter,                              &
                                                              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC,                                   &
         SHORT_NAME = 'RADSWNA',                                         &
         LONG_NAME  = 'air_temperature_tendency_due_to_shortwave_no_aerosol', &
         UNITS      = 'K s-1',                                           &
         DIMS       = MAPL_DimsHorzVert,                                 &
         VLOCATION  = MAPL_VLocationCenter,                              &
                                                              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC,                                   &
         SHORT_NAME = 'RADLWCNA',                                        &
         LONG_NAME  = 'air_temperature_tendency_due_to_longwave_for_clear_skies_no_aerosol',&
         UNITS      = 'K s-1',                                           &
         DIMS       = MAPL_DimsHorzVert,                                 &
         VLOCATION  = MAPL_VLocationCenter,                              &
                                                              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC,                                   &
         SHORT_NAME = 'RADSWCNA',                                        &
         LONG_NAME  = 'air_temperature_tendency_due_to_shortwave_for_clear_skies_no_aerosol',&
         UNITS      = 'K s-1',                                           &
         DIMS       = MAPL_DimsHorzVert,                                 &
         VLOCATION  = MAPL_VLocationCenter,                              &
                                                              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC,                                   &
         SHORT_NAME = 'RADSRF',                                          &
         LONG_NAME  = 'net_downwelling_radiation_at_surface',            &
         UNITS      = 'W m-2',                                           &
         DIMS       = MAPL_DimsHorzOnly,                                 &
         VLOCATION  = MAPL_VLocationNone,                                &
                                                              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC,                                   &
         SHORT_NAME = 'ALW',                                             &
         LONG_NAME  = 'linearization_of_surface_upwelling_longwave_flux',&
         UNITS      = 'W m-2',                                           &
         DIMS       = MAPL_DimsHorzOnly,                                 &
         VLOCATION  = MAPL_VLocationNone,                                &
                                                              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC,                                   &
         SHORT_NAME = 'BLW',                                             &
         LONG_NAME  = 'linearization_of_surface_upwelling_longwave_flux',&
         UNITS      = 'W m-2 K-1',                                       &
         DIMS       = MAPL_DimsHorzOnly,                                 &
         VLOCATION  = MAPL_VLocationNone,                                &
                                                              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                                &
         SHORT_NAME = 'DRPAR',                                           &
         CHILD_ID = SOL,                                                 &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                                &
         SHORT_NAME = 'DFPAR',                                           &
         CHILD_ID = SOL,                                                 &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                                &
         SHORT_NAME = 'DRNIR',                                           &
         CHILD_ID = SOL,                                                 &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                                &
         SHORT_NAME = 'DFNIR',                                           &
         CHILD_ID = SOL,                                                 &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                                &
         SHORT_NAME = 'DRUVR',                                           &
         CHILD_ID = SOL,                                                 &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                                &
         SHORT_NAME = 'DFUVR',                                           &
         CHILD_ID = SOL,                                                 &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                                &
         SHORT_NAME = 'DRPARN',                                          &
         CHILD_ID = SOL,                                                 &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                                &
         SHORT_NAME = 'DFPARN',                                          &
         CHILD_ID = SOL,                                                 &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                                &
         SHORT_NAME = 'DRNIRN',                                          &
         CHILD_ID = SOL,                                                 &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                                &
         SHORT_NAME = 'DFNIRN',                                          &
         CHILD_ID = SOL,                                                 &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                                &
         SHORT_NAME = 'DRUVRN',                                          &
         CHILD_ID = SOL,                                                 &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                                &
         SHORT_NAME = 'DFUVRN',                                          &
         CHILD_ID = SOL,                                                 &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                                &
         SHORT_NAME = 'FCLD',                                            &
         CHILD_ID = SOL,                                                 &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                            &
         SHORT_NAME =  'TAUCLI',                                         &
         CHILD_ID = SOL,                                                 &
         RC=STATUS  )
    VERIFY_(STATUS)

        call MAPL_AddExportSpec ( GC   ,                            &
         SHORT_NAME =  'TAUCLW',                                         &
         CHILD_ID = SOL,                                                 &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                                &
         SHORT_NAME = 'LWS',                                             &
         CHILD_ID = IRR,                                                 &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                                &
         SHORT_NAME = 'LWS0',                                            &
         CHILD_ID = IRR,                                                 &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                                &
         SHORT_NAME = 'CLDTT',                                           &
         CHILD_ID = SOL,                                                 &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                                &
         SHORT_NAME = 'ALBEDO',                                          &
         CHILD_ID = SOL,                                                 &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                                &
         SHORT_NAME = 'FSWBAND',                                         &
         CHILD_ID = SOL,                                                 &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                                &
         SHORT_NAME = 'FSWBANDNA',                                       &
         CHILD_ID = SOL,                                                 &
         RC=STATUS  )
    VERIFY_(STATUS)

!EOS

! Set internal connections between the children`s IMPORTS and EXPORTS
! -------------------------------------------------------------------

!!BOP

! !CONNECTIONS:

! Solar imports
!-------------------

   if (USE_SATSIM > 0 ) then
    call MAPL_AddConnectivity ( GC,                                &
         SHORT_NAME  = (/ 'MCOSZ' /),                              &
         DST_ID      =  STM,                                       &
         SRC_ID      =  SOL,                                       &
                                                        RC=STATUS  )
    VERIFY_(STATUS)
    endif

!!EOP


! Set generic init and final methods for us and for our children
! --------------------------------------------------------------

    call MAPL_GenericSetServices    ( gc, RC=STATUS)
    VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)

  end subroutine SetServices


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !IROUTINE: Initialize -- Initialize method for the composite SuperDyn Gridded Component

! !INTERFACE:

  subroutine Initialize ( GC, IMPORT, EXPORT, CLOCK, RC )

! !USES:
  use cloud_condensate_inhomogeneity, only: initialize_inhomogeneity
  use cloud_subcol_gen, only : initialize_cloud_subcol_gen, &
    def_aam1, def_aam2, def_aam30, def_aam4, &
    def_ram1, def_ram2, def_ram30, def_ram4

! !ARGUMENTS:

  type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
  type(ESMF_State),    intent(inout) :: IMPORT ! Import state
  type(ESMF_State),    intent(inout) :: EXPORT ! Export state
  type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
  integer, optional,   intent(  out) :: RC     ! Error code

! !DESCRIPTION: 
!
!  A component specific Initialize is present to deal with Aerosols.
!  This should be done away with when a better treatment of 
!  aerosol optical properties.
!  The Aerosol optical tables are currently in a global module
!  shared by the two children; these are initialized here. 

!EOP

! ErrLog Variables

  character(len=ESMF_MAXSTR)          :: IAm
  integer                             :: STATUS
  character(len=ESMF_MAXSTR)          :: COMP_NAME

! Local derived type aliases

  type (MAPL_MetaComp),      pointer  :: MAPL
! type (ESMF_State),         pointer  :: GIM(:)
! type (ESMF_Config)                  :: CF

! Correlation length parameters from resource file

  real :: aam1, aam2, aam30, aam4
  real :: ram1, ram2, ram30, ram4

!=============================================================================

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    Iam = "Initialize"
    call ESMF_GridCompGet (GC, NAME=COMP_NAME, __RC__)  ! CONFIG=CF, __RC__)
    Iam = trim(COMP_NAME) // trim(Iam)

! Generic initialize
!-------------------

    call MAPL_GenericInitialize (GC, IMPORT, EXPORT, CLOCK, __RC__)

! Get my MAPL object
!-------------------

    call MAPL_GetObjectFromGC (GC, MAPL, __RC__)

! Start Total timer after generic initialize
!-------------------------------------------

    call MAPL_TimerOn (MAPL,"TOTAL")

! Get parameters from MAPL
!-------------------------

!   call MAPL_Get (MAPL, GIM=GIM, __RC__)

! Initialize module-level cloud generator details.
! Currently these are only used by RRTMG SW and LW, so should probably
! make conditional on either RRTMG SW and LW being used. Dont bother
! for now. Also, may later use the cloud generator for RRTMGP as well.
!---------------------------------------------------------------------

! Set up RRTMG condensate inhomogeneity tables

    call initialize_inhomogeneity(1)

! Set RRTMG cloud subcolumn generator correlation length parameters to non-default values
! from MAPL resource parameters. Comment out to just use defaults in module cloud_subcol_gen.

    call MAPL_GetResource(MAPL,aam1 ,LABEL="ADL_AM1:" ,default=def_aam1 ,__RC__)
    call MAPL_GetResource(MAPL,aam2 ,LABEL="ADL_AM2:" ,default=def_aam2 ,__RC__)
    call MAPL_GetResource(MAPL,aam30,LABEL="ADL_AM30:",default=def_aam30,__RC__)
    call MAPL_GetResource(MAPL,aam4 ,LABEL="ADL_AM4:" ,default=def_aam4 ,__RC__)
    call MAPL_GetResource(MAPL,ram1 ,LABEL="RDL_AM1:" ,default=def_ram1 ,__RC__)
    call MAPL_GetResource(MAPL,ram2 ,LABEL="RDL_AM2:" ,default=def_ram2 ,__RC__)
    call MAPL_GetResource(MAPL,ram30,LABEL="RDL_AM30:",default=def_ram30,__RC__)
    call MAPL_GetResource(MAPL,ram4 ,LABEL="RDL_AM4:" ,default=def_ram4 ,__RC__)
    call initialize_cloud_subcol_gen( &
      adl_am1=aam1, adl_am2=aam2, adl_am30=aam30, adl_am4=aam4, &
      rdl_am1=ram1, rdl_am2=ram2, rdl_am30=ram30, rdl_am4=ram4)

! All done
!---------

    call MAPL_TimerOff (MAPL,"TOTAL")
    RETURN_(ESMF_SUCCESS)

  end subroutine Initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!BOP

! !IROUTINE: RUN -- Run method for the composite radiation component

! !INTERFACE:

subroutine RUN ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

  type(ESMF_GridComp),  intent(inout) :: GC     ! Gridded component 
  type(ESMF_State),     intent(inout) :: IMPORT ! Import state
  type(ESMF_State),     intent(inout) :: EXPORT ! Export state
  type(ESMF_Clock),     intent(inout) :: CLOCK  ! The clock
  integer, optional,    intent(  out) :: RC     ! Error code:

! !DESCRIPTION: Calls the run methods of solar and irrad and combines
!     their fluxes into a single pressure-weighted temperature tendency.
!  \newline
!
!  The Children's Run method is called with the Clock that appears in 
!  {\tt GEOS\_RadiationGridCompMod}'s Run method.
!  It assumes that, on every call to
!  their Run method,  both children update their Exports
!  to values consisten with the current time on that Clock. \newline
!
!  The code assumes a single Aerosol bundle will go to both
!  children. The bundle is passed to the children with a MAPL call. \newline
!
!EOP

! ErrLog Variables

  character(len=ESMF_MAXSTR)          :: IAm
  integer                             :: STATUS
  character(len=ESMF_MAXSTR)          :: COMP_NAME

! Local derived type aliases

  type (MAPL_MetaComp),  pointer  :: MAPL
  type (ESMF_State       ),  pointer  :: GEX(:  )
  
  integer                             :: I, L

!  Pointers to IR exports

  real, pointer, dimension(:,:  )     :: SFCEM
  real, pointer, dimension(:,:  )     :: DSFDTS
  real, pointer, dimension(:,:  )     :: TRD
  real, pointer, dimension(:,:,:)     :: FLW
  real, pointer, dimension(:,:,:)     :: FLWCLR
  real, pointer, dimension(:,:,:)     :: FLA

!  Pointers to SOLAR exports

  real, pointer, dimension(:,:,:)     :: FSW
  real, pointer, dimension(:,:,:)     :: FSWCLR
  real, pointer, dimension(:,:,:)     :: FSWNA
  real, pointer, dimension(:,:,:)     :: FSCNA  
 
!  Pointers to imports

  real, pointer, dimension(:,:,:)     :: PLE

!  Pointers to exports

  real, pointer, dimension(:,:,:)     :: DTDT
  real, pointer, dimension(:,:,:)     :: RADLW
  real, pointer, dimension(:,:,:)     :: RADSW
  real, pointer, dimension(:,:,:)     :: RADLWC
  real, pointer, dimension(:,:,:)     :: RADSWC
  real, pointer, dimension(:,:,:)     :: RADSWNA
  real, pointer, dimension(:,:,:)     :: RADLWCNA
  real, pointer, dimension(:,:,:)     :: RADSWCNA
  real, pointer, dimension(:,:  )     :: ALW
  real, pointer, dimension(:,:  )     :: BLW
  real, pointer, dimension(:,:  )     :: RADSRF

! Locals

  real, pointer, dimension(:,:,:)     :: DMI
  integer                             :: IM
  integer                             :: JM
  integer                             :: LM

!=============================================================================

! Begin... 

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    call ESMF_GridCompGet( GC, name=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // "Run"

! Get my internal MAPL_Generic state
!-----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)

! Start Total timer
!------------------

    call MAPL_TimerOn(MAPL,"TOTAL")

! Get parameters from generic state.
!-----------------------------------

    call MAPL_Get(MAPL,  &
         GEX=GEX,        &
         IM=IM, JM=JM, LM=LM,         &
                            RC=STATUS )
    VERIFY_(STATUS)

! Get pointers to exports
!------------------------

    call MAPL_GetPointer ( IMPORT, PLE    , 'PLEINST'  ,  RC=STATUS )
    VERIFY_(STATUS)

! Get pointers to exports
!------------------------

    call MAPL_GetPointer ( EXPORT, DTDT    , 'DTDT'    ,  RC=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetPointer ( EXPORT, ALW     , 'ALW'     ,  RC=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetPointer ( EXPORT, BLW     , 'BLW'     ,  RC=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetPointer ( EXPORT, RADSRF  , 'RADSRF'  ,  RC=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetPointer ( EXPORT, RADLW   , 'RADLW'   ,  RC=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetPointer ( EXPORT, RADSW   , 'RADSW'   ,  RC=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetPointer ( EXPORT, RADLWC  , 'RADLWC'  ,  RC=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetPointer ( EXPORT, RADSWC  , 'RADSWC'  ,  RC=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetPointer ( EXPORT, RADSWNA , 'RADSWNA' ,  RC=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetPointer ( EXPORT, RADLWCNA, 'RADLWCNA',  RC=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetPointer ( EXPORT, RADSWCNA, 'RADSWCNA',  RC=STATUS )
    VERIFY_(STATUS)

! Allocate children's exports that we need
!-----------------------------------------

    if (associated(RADSW  ) .or. associated(DTDT   ) .or. associated(RADSRF) ) then
       call MAPL_GetPointer ( GEX(SOL), FSW   , 'FSW'    ,  alloc=.TRUE.,RC=STATUS )
       VERIFY_(STATUS)
    end if

    if (associated(RADSWC)                            ) then
       call MAPL_GetPointer ( GEX(SOL), FSWCLR, 'FSC'    ,  alloc=.TRUE.,RC=STATUS )
       VERIFY_(STATUS)
    end if

    if (associated(RADLW) .or. associated(DTDT   )   .or. associated(RADSRF)) then
       call MAPL_GetPointer ( GEX(IRR), FLW   , 'FLX'    ,  alloc=.TRUE.,RC=STATUS )
       VERIFY_(STATUS)
    end if

    if (associated(RADLWC)                            ) then
       call MAPL_GetPointer ( GEX(IRR), FLWCLR, 'FLC'    ,  alloc=.TRUE.,RC=STATUS )
       VERIFY_(STATUS)
    end if

    if (associated(RADSWNA) ) then
       call MAPL_GetPointer ( GEX(SOL), FSWNA , 'FSWNA'  ,  alloc=.TRUE.,RC=STATUS )
       VERIFY_(STATUS)
    end if

    if (associated(RADSWCNA)                          ) then
       call MAPL_GetPointer ( GEX(SOL), FSCNA,  'FSCNA'  ,  alloc=.TRUE.,RC=STATUS )
       VERIFY_(STATUS)
    end if

    if (associated(RADLWCNA)                          ) then
       call MAPL_GetPointer ( GEX(IRR), FLA   , 'FLA'    ,  alloc=.TRUE.,RC=STATUS )
       VERIFY_(STATUS)
    end if

    if (  associated(ALW   )  .or. associated(BLW     ) ) then
       call MAPL_GetPointer ( GEX(IRR), DSFDTS, 'DSFDTS0',  alloc=.TRUE.,RC=STATUS )
       VERIFY_(STATUS)
       call MAPL_GetPointer ( GEX(IRR), SFCEM , 'SFCEM0' ,  alloc=.TRUE.,RC=STATUS )
       VERIFY_(STATUS)
       call MAPL_GetPointer ( GEX(IRR), TRD   , 'TSREFF' ,  alloc=.TRUE.,RC=STATUS )
       VERIFY_(STATUS)
    end if

! Run the child components and their couplers
!--------------------------------------------

    call MAPL_TimerOff(MAPL,"TOTAL")
    call MAPL_GenericRunChildren (GC, IMPORT, EXPORT, CLOCK, RC=STATUS )
    VERIFY_(STATUS)
    call MAPL_TimerOn (MAPL,"TOTAL")

! Prepare exports
!----------------

    if( associated (BLW     ) ) BLW      =  DSFDTS
    if( associated (ALW     ) ) ALW      =  SFCEM - DSFDTS*TRD
    if( associated (RADSRF  ) ) RADSRF   =   (FSW(:,:,  LM  ) + FLW(:,:,  LM))
    if( associated (DTDT    ) ) DTDT     = ( (FLW(:,:,0:LM-1) - FLW(:,:,1:LM)) + &
                                             (FSW(:,:,0:LM-1) - FSW(:,:,1:LM)) ) * (MAPL_GRAV/MAPL_CP)

    if( associated (RADLW ) .or. associated (RADSW )   .or. &
        associated (RADLWC) .or. associated (RADSWC)   .or. &
        associated (RADSWNA).or. associated (RADSWCNA) .or. &
        associated (RADLWCNA)                             ) then

       allocate(DMI(IM,JM,LM),stat=STATUS)
       VERIFY_(STATUS)
       DMI = MAPL_GRAV/(MAPL_CP*(PLE(:,:,1:LM)-PLE(:,:,0:LM-1)))

       if( associated (RADLW   ) ) RADLW    = (FLW   (:,:,0:LM-1) - FLW   (:,:,1:LM))*DMI
       if( associated (RADSW   ) ) RADSW    = (FSW   (:,:,0:LM-1) - FSW   (:,:,1:LM))*DMI
       if( associated (RADLWC  ) ) RADLWC   = (FLWCLR(:,:,0:LM-1) - FLWCLR(:,:,1:LM))*DMI
       if( associated (RADSWC  ) ) RADSWC   = (FSWCLR(:,:,0:LM-1) - FSWCLR(:,:,1:LM))*DMI
       if( associated (RADSWNA ) ) RADSWNA  = (FSWNA (:,:,0:LM-1) - FSWNA (:,:,1:LM))*DMI
       if( associated (RADLWCNA) ) RADLWCNA = (FLA   (:,:,0:LM-1) - FLA   (:,:,1:LM))*DMI
       if( associated (RADSWCNA) ) RADSWCNA = (FSCNA (:,:,0:LM-1) - FSCNA (:,:,1:LM))*DMI

       deallocate(DMI,stat=STATUS)
       VERIFY_(STATUS)

    end if

    call MAPL_TimerOff(MAPL,"TOTAL")

    RETURN_(ESMF_SUCCESS)

  end subroutine RUN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module GEOS_RadiationGridCompMod

