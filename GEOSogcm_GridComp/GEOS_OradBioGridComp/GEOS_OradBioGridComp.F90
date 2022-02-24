!  $Id: GEOS_OradBioGridComp.F90,v 1.8.2.2.4.1.2.5.6.2.4.2 2018/10/15 12:27:11 croussea Exp $

#include "MAPL_Generic.h"
!=============================================================================

module GEOS_OradBioGridCompMod

!BOP

! !MODULE: GEOS_OradBioGridCompMod -- Solar radiation on the ocean


! !USES:

  use ESMF
  use MAPL
  
  implicit none
  private

! !PUBLIC ROUTINES:
  public SetServices

  integer, parameter :: NB_CHOU_UV  = 5 ! Number of UV bands
  integer, parameter :: NB_CHOU_NIR = 3 ! Number of near-IR bands
  integer, parameter :: NB_CHOU     = NB_CHOU_UV + NB_CHOU_NIR ! Total number of bands
  
  logical            :: ocean_extData

#include "definebio.h"
#include "comlte.h"

  type T_ORADBIO_STATE
    private
    integer :: lam(nlt)
    real    :: aw(nlt),bw(nlt),excdom(nlt),exdet(nlt),WtoQ(nlt),wfac(nlt)
    real    :: ac(nchl,nlt),bc(nchl,nlt)
    real    :: bpic(nlt)
    real    :: rad,pi2
    real    :: Fobar(nlt),thray(nlt),oza(nlt),awv(nlt),ao(nlt),aco2(nlt)
    real    :: am,Vi
    real   :: asl(ncld),bsl(ncld),csl(ncld),dsl(ncld),esl(ncld),fsl(ncld)
    integer :: ica(nlt)
    type(ESMF_Alarm) :: DayAlarm
  end type T_ORADBIO_STATE

   type ORADBIO_WRAP
      type (T_ORADBIO_STATE), pointer :: PTR => null()
   end type ORADBIO_WRAP

!  Structure to create arrays for aerosol variables
  type bandptr
    real, pointer, dimension(:,:) :: b => null()
  end type bandptr

  integer :: nl
  character(len = 2) :: suffix

!========================================================================

! !DESCRIPTION:
! 
!   {\tt GEOS\_Orad} is a light-weight gridded component that updates
!      the solar radiation penetrating the ocean. It implements 
!      absorption and scattering of solar radiation in the ocean using 
!      the OASIM (Gregg and Casey, 2009, JMS) radiation package. It 
!      interfaces with the biology package, producing radiation in 33 
!      spectral bands and with the GEOS AGCM radiation package by 
!      accepting solar irradiance at the surface in 8 bands and 
!      returning a flux difference for each layer
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

! ! DESCRIPTION: This version uses the GEOS_GenericSetServices, which 
!   sets the Initialize and Finalize services, as well as allocating
!   our instance of a generic state and putting it in the 
!   gridded component (GC). Here we only need to set the run method and
!   add the state variable specifications (also generic) to our instance
!   of the generic state. This is the way our true state variables get 
!   into the ESMF_State INTERNAL, which is in the GEOS_GenericState. The 
!   import and internal variables are allocated and initialized by 
!   generic.  Here generic is used for the ocean grid.

!EOP

!=============================================================================

! ErrLog Variables

    character(len=ESMF_MAXSTR)              :: IAm
    integer                                 :: STATUS
    character(len=ESMF_MAXSTR)              :: COMP_NAME

! Local derived type aliases

    type (MAPL_MetaComp    ), pointer       :: MAPL => null()
!=============================================================================

! Begin...

! Get my name and set-up traceback handle
! ---------------------------------------

    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // 'SetServices'

! Set the Run entry point
! -----------------------

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_INITIALIZE, Initialize, RC=STATUS )
    VERIFY_(STATUS)
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,  Run,        RC=STATUS )
    VERIFY_(STATUS)

! Set the state variable specs.
! -----------------------------

   call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
   VERIFY_(STATUS)

   call MAPL_GetResource (MAPL,   ocean_extData, Label="OCEAN_EXT_DATA:",   DEFAULT=.FALSE., __RC__ ) ! .TRUE. or .FALSE.

!BOC

!  !IMPORT STATE:

    call MAPL_AddImportSpec(GC,                               &
    LONG_NAME  = 'cosine_of_the_solar_zenith_angle',          &
    UNITS      = '1',                                         &
    SHORT_NAME = 'COSZ',                                      &
    DIMS       = MAPL_DimsHorzOnly,                           &
    VLOCATION  = MAPL_VLocationNone,                          &
    RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                               &
    SHORT_NAME = 'PENUVR',                                    &
    LONG_NAME  = 'net_downward_penetrating_direct_UV_flux',   &
    UNITS      = 'W m-2',                                     &
    DIMS       = MAPL_DimsHorzOnly,                           &
    VLOCATION  = MAPL_VLocationNone,                          &
    RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                               &
    SHORT_NAME = 'PENUVF',                                    &
    LONG_NAME  = 'net_downward_penetrating_diffuse_UV_flux',  &
    UNITS      = 'W m-2',                                     &
    DIMS       = MAPL_DimsHorzOnly,                           &
    VLOCATION  = MAPL_VLocationNone,                          &
    RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                               &
    SHORT_NAME = 'PENPAR',                                    &
    LONG_NAME  = 'net_downward_penetrating_direct_PAR_flux',  &
    UNITS      = 'W m-2',                                     &
    DIMS       = MAPL_DimsHorzOnly,                           &
    VLOCATION  = MAPL_VLocationNone,                          &
    RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                               &
    SHORT_NAME = 'PENPAF',                                    &
    LONG_NAME  = 'net_downward_penetrating_diffuse_PAR_flux', &
    UNITS      = 'W m-2',                                     &
    DIMS       = MAPL_DimsHorzOnly,                           &
    VLOCATION  = MAPL_VLocationNone,                          &
    RC=STATUS  )
    VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC                         ,&
          LONG_NAME          = 'net_surface_downwelling_nir_beam_flux',&
          UNITS              = 'W m-2'                       ,&
          SHORT_NAME         = 'DRNIR'                       ,&
          DIMS               = MAPL_DimsHorzOnly             ,&
          VLOCATION          = MAPL_VLocationNone            ,&
          RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC                         ,&
          LONG_NAME          = 'net_surface_downwelling_nir_diffuse_flux',&
          UNITS              = 'W m-2'                       ,&
          SHORT_NAME         = 'DFNIR'                       ,&
          DIMS               = MAPL_DimsHorzOnly             ,&
          VLOCATION          = MAPL_VLocationNone            ,&
          RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                               &
    SHORT_NAME = 'FROCEAN',                                   &
    LONG_NAME  = 'ocean_fraction_of_grid_cell',               &
    UNITS      = '1',                                         &
    DIMS       = MAPL_DimsHorzOnly,                           &
    VLOCATION  = MAPL_VLocationNone,                          &
    RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                               &
    SHORT_NAME = 'DH',                                        &
    LONG_NAME  = 'layer_thickness',                           &
    UNITS      = 'm',                                         &
    DIMS       = MAPL_DimsHorzVert,                           &
    VLOCATION  = MAPL_VLocationCenter,                        &
    RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                               &
    SHORT_NAME = 'DIATOM',                                    &
    LONG_NAME  = 'diatom_concentration',                      &
    UNITS      = 'mg m-3',                                    &
    DIMS       = MAPL_DimsHorzVert,                           &
    VLOCATION  = MAPL_VLocationCenter,                        &
    RESTART    = MAPL_RestartSkip,                            &
    RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                               &
    SHORT_NAME = 'CHLORO',                                    &
    LONG_NAME  = 'chlorophyte_concentration',                 &
    UNITS      = 'mg m-3',                                    &
    DIMS       = MAPL_DimsHorzVert,                           &
    VLOCATION  = MAPL_VLocationCenter,                        &
    RESTART    = MAPL_RestartSkip,                            &
    RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                               &
    SHORT_NAME = 'CYANO',                                     &
    LONG_NAME  = 'cyano-bacteria_concentration',              &
    UNITS      = 'mg m-3',                                    &
    DIMS       = MAPL_DimsHorzVert,                           &
    VLOCATION  = MAPL_VLocationCenter,                        &
    RESTART    = MAPL_RestartSkip,                            &
    RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                               &
    SHORT_NAME = 'COCCO',                                     &
    LONG_NAME  = 'coccolithophore_concentration',             &
    UNITS      = 'mg m-3',                                    &
    DIMS       = MAPL_DimsHorzVert,                           &
    VLOCATION  = MAPL_VLocationCenter,                        &
    RESTART    = MAPL_RestartSkip,                            &
    RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                               &
    SHORT_NAME = 'DINO',                                      &
    LONG_NAME  = 'dinoflagellate_concentration',              &
    UNITS      = 'mg m-3',                                    &
    DIMS       = MAPL_DimsHorzVert,                           &
    VLOCATION  = MAPL_VLocationCenter,                        &
    RESTART    = MAPL_RestartSkip,                            &
    RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                               &
    SHORT_NAME = 'PHAEO',                                     &
    LONG_NAME  = 'phaeocystis_concentration',                 &
    UNITS      = 'mg m-3',                                    &
    DIMS       = MAPL_DimsHorzVert,                           &
    VLOCATION  = MAPL_VLocationCenter,                        &
    RESTART    = MAPL_RestartSkip,                            &
    RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                               &
    SHORT_NAME = 'CDET',                                      &
    LONG_NAME  = 'carbon/nitrogen_detritus_concentration',    &
    UNITS      = 'ugC l-1',                                   &
    DIMS       = MAPL_DimsHorzVert,                           &
    VLOCATION  = MAPL_VLocationCenter,                        &
    RESTART    = MAPL_RestartSkip,                            &
    RC=STATUS  )
    VERIFY_(STATUS)
    
    call MAPL_AddImportSpec(GC,                               &
    SHORT_NAME = 'PIC',                                       &
    LONG_NAME  = 'particulate inorganic carbon',              &
    UNITS      = 'ugC l-1',                                   &
    DIMS       = MAPL_DimsHorzVert,                           &
    VLOCATION  = MAPL_VLocationCenter,                        &
    RESTART    = MAPL_RestartSkip,                            &
    RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                               &
    SHORT_NAME = 'CDC',                                       &
    LONG_NAME  = 'Colored dissolved organic carbon',          &
    UNITS      = 'uM',                                        &
    DIMS       = MAPL_DimsHorzVert,                           &
    VLOCATION  = MAPL_VLocationCenter,                        &
    RESTART    = MAPL_RestartSkip,                            &
    RC=STATUS  )
    VERIFY_(STATUS)
    
    call MAPL_AddImportSpec(GC,                               &
    SHORT_NAME = 'AVGQ',                                      &
    LONG_NAME  = 'average quantum irradiance',                &
    UNITS      = 'umol m-2 s-1',                              &
    DIMS       = MAPL_DimsHorzVert,                           &
    VLOCATION  = MAPL_VLocationCenter,                        &
    RESTART    = MAPL_RestartSkip,                            &
    RC=STATUS  )
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
    RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                               &
    SHORT_NAME = 'OZ',                                        &
    LONG_NAME  = 'ozone thickness',                           &
    UNITS      = 'Dobson units',                              &
    DIMS       = MAPL_DimsHorzOnly,                           &
    VLOCATION  = MAPL_VLocationNone,                          &
    RESTART    = MAPL_RestartSkip,                            &
    RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                               &
    SHORT_NAME = 'WV',                                        &
    LONG_NAME  = 'water vapor',                               &
    UNITS      = 'cm',                                        &
    DIMS       = MAPL_DimsHorzOnly,                           &
    VLOCATION  = MAPL_VLocationNone,                          &
    RESTART    = MAPL_RestartSkip,                            &
    RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                               &
    SHORT_NAME = 'RH',                                        &
    LONG_NAME  = 'relative humidity',                         &
    UNITS      = 'percent',                                   &
    DIMS       = MAPL_DimsHorzOnly,                           &
    VLOCATION  = MAPL_VLocationNone,                          &
    RESTART    = MAPL_RestartSkip,                            &
    RC=STATUS  )
    VERIFY_(STATUS)

    do nl = 1,nlt
     write(unit = suffix, fmt = '(i2.2)') nl
     call MAPL_AddImportSpec(GC,                              &
     SHORT_NAME = 'TAUA_'//suffix,                            &
     LONG_NAME  = 'aerosol optical thickness',                &
     UNITS      = '',                                         &
     DIMS       = MAPL_DimsHorzOnly,                          &
     VLOCATION  = MAPL_VLocationNone,                         &
     RESTART    = MAPL_RestartSkip,                           &
     RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                              &
     SHORT_NAME = 'SSALB_'//suffix,                           &
     LONG_NAME  = 'single scattering albedo',                 &
     UNITS      = '',                                         &
     DIMS       = MAPL_DimsHorzOnly,                          &
     VLOCATION  = MAPL_VLocationNone,                         &
     RESTART    = MAPL_RestartSkip,                           &
     RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                              &
     SHORT_NAME = 'ASYMP_'//suffix,                           &
     LONG_NAME  = 'asymmetry parameter',                      &
     UNITS      = '',                                         &
     DIMS       = MAPL_DimsHorzOnly,                          &
     VLOCATION  = MAPL_VLocationNone,                         &
     RESTART    = MAPL_RestartSkip,                           &
     RC=STATUS  )
     VERIFY_(STATUS)
    enddo

    call MAPL_AddImportSpec(GC,                               &
    SHORT_NAME = 'CCOVM',                                     &
    LONG_NAME  = 'cloud cover',                               &
    UNITS      = 'percent',                                   &
    DIMS       = MAPL_DimsHorzOnly,                           &
    VLOCATION  = MAPL_VLocationNone,                          &
    RESTART    = MAPL_RestartSkip,                            &
    RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                               &
    SHORT_NAME = 'CLDTCM',                                    &
    LONG_NAME  = 'cloud optical thickness',                   &
    UNITS      = 'dimensionless',                             &
    DIMS       = MAPL_DimsHorzOnly,                           &
    VLOCATION  = MAPL_VLocationNone,                          &
    RESTART    = MAPL_RestartSkip,                            &
    RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                               &
    SHORT_NAME = 'RLWPM',                                     &
    LONG_NAME  = 'cloud liquid water path',                   &
    UNITS      = 'dimensionless',                             &
    DIMS       = MAPL_DimsHorzOnly,                           &
    VLOCATION  = MAPL_VLocationNone,                          &
    RESTART    = MAPL_RestartSkip,                            &
    RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                               &
    SHORT_NAME = 'CDREM',                                     &
    LONG_NAME  = 'cloud droplet effective radius',            &
    UNITS      = 'dimensionless',                             &
    DIMS       = MAPL_DimsHorzOnly,                           &
    VLOCATION  = MAPL_VLocationNone,                          &
    RESTART    = MAPL_RestartSkip,                            &
    RC=STATUS  )
    VERIFY_(STATUS)

    if (ocean_extData) then
      call MAPL_AddImportSpec(GC,                             &
           SHORT_NAME = 'DATA_KPAR',                          &
           LONG_NAME  = 'PAR_extinction_coefficient',         &
           UNITS      = 'm-1',                                &
           DIMS       = MAPL_DimsHorzOnly,                    &
           VLOCATION  = MAPL_VLocationNone,                   &
           RC=STATUS  )
      VERIFY_(STATUS)
    endif

!  !EXPORT STATE:

    call MAPL_AddExportSpec(GC,                               &
    SHORT_NAME = 'SWHEAT',                                    &
    LONG_NAME  = 'delta_irradiance_at_each_depth',            &
    UNITS      = 'W m-2',                                     &
    DIMS       = MAPL_DimsHorzVert,                           &
    VLOCATION  = MAPL_VLocationCenter,                        &
    RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
    SHORT_NAME = 'KPAR',                                      &
    LONG_NAME  = 'PAR_extinction_coefficient',                &
    UNITS      = 'm-1',                                       &
    DIMS       = MAPL_DimsHorzOnly,                           &
    VLOCATION  = MAPL_VLocationNone,                          &
    RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
    SHORT_NAME = 'TIRRQ',                                     &
    LONG_NAME  = 'total_irradiance',                          &
    UNITS      = 'umol quanta m-2 s-1',                       &
    DIMS       = MAPL_DimsHorzVert,                           &
    VLOCATION  = MAPL_VLocationCenter,                        &
    RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
    SHORT_NAME = 'CDOMABSQ',                                  &
    LONG_NAME  = 'absorption of quanta by CDOM',              &
    UNITS      = 'umol quanta m-2 s-1',                       &
    DIMS       = MAPL_DimsHorzVert,                           &
    VLOCATION  = MAPL_VLocationCenter,                        &
    RC=STATUS  )
    VERIFY_(STATUS)
   
 !EOC
 
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

! ! DESCRIPTION: This function initializes the ORADBIO gridded component.

!EOPI

!------------------------------------------------------------------------

! Local Variables

    character(len=ESMF_MAXSTR)          :: IAm
    integer                             :: STATUS
    character(len=ESMF_MAXSTR)          :: COMP_NAME

    type (MAPL_MetaComp),   pointer     :: MAPL => null()
    type (T_ORADBIO_STATE), pointer     :: State => null()
    type (ORADBIO_wrap)                 :: WRAP
    type (ESMF_TimeInterval)            :: RingInterval
!    type (ESMF_Bundle)                  :: BUNDLE

! Begin
!------

    Iam = "Initialize"
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
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

    call ESMF_UserCompSetInternalState ( GC, 'ORADBIO_state', WRAP, STATUS ); VERIFY_(STATUS)

! Initialize OASIM parameters in private state
!-----------------------------------------------

    call setlte(State%lam,   State%aw,   State%bw,                    &
                State%ac,    State%bc,   State%bpic,                  &
                State%excdom,State%exdet,State%WtoQ,State%wfac)
    call setsfclte(State%rad,State%pi2,State%lam,State%Fobar,         &
                   State%thray,State%oza,State%awv,State%ao,          &
                   State%aco2, State%am, State%Vi,                    &
                   State%asl,State%bsl,State%csl,State%dsl,State%esl, &
                   State%fsl,State%ica)

! Set Daily alarm
!----------------

    call ESMF_TimeIntervalSet(RingInterval,  D=1, RC=STATUS)
    VERIFY_(STATUS)

    State%DayAlarm = ESMF_AlarmCreate(NAME="OradAlarm", &
                     CLOCK        = CLOCK,         &
                     RingInterval = RingInterval,  &
                     Sticky       = .FALSE.,       &
                     RC=STATUS  )
    VERIFY_(STATUS)

! Stop Total timer
!-----------------

    call MAPL_TimerOff(MAPL,"INITIALIZE" )
    call MAPL_TimerOff(MAPL,"TOTAL")

! All Done
!---------

    RETURN_(ESMF_SUCCESS)

  end subroutine Initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!BOP

! ! IROUTINE: RUN -- First Run stage for the Orad component

! !INTERFACE:

  subroutine RUN ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
    type(ESMF_State),    intent(inout) :: IMPORT ! Import state
    type(ESMF_State),    intent(inout) :: EXPORT ! Export state
    type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
    integer, optional,   intent(  out) :: RC     ! Error code:

! ! DESCRIPTION: full radiative transfer calculation (OASIM).

!EOP

! ErrLog Variables

    character(len=ESMF_MAXSTR)          :: IAm
    integer                             :: STATUS
    character(len=ESMF_MAXSTR)          :: COMP_NAME

! Locals

    type (MAPL_MetaComp),   pointer :: MAPL => null()
    type (T_ORADBIO_STATE), pointer :: State => null()
    type (ORADBIO_wrap)             :: WRAP
    type (ESMF_State)               :: INTERNAL
    type (ESMF_Time)                :: currTime
    type (ESMF_TimeInterval)        :: DTint
    type (MAPL_SunOrbit)            :: ORBIT

    real, dimension(:,:),   allocatable :: COSZ
    real, dimension(:,:),   allocatable :: SLR
    real, dimension(:,:),   allocatable :: PHYTO
    real, dimension(:,:),   pointer     :: LATS => null()
    real, dimension(:,:),   pointer     :: LONS => null()
    real, dimension(:,:,:), pointer     :: DIATOM => null()
    real, dimension(:,:,:), pointer     :: CHLORO => null()
    real, dimension(:,:,:), pointer     :: CYANO => null()
    real, dimension(:,:,:), pointer     :: COCCO => null()
    real, dimension(:,:,:), pointer     :: DINO => null()
    real, dimension(:,:,:), pointer     :: PHAEO => null()
    real, dimension(:,:,:), pointer     :: AVGQ => null()
    real, dimension(:,:,:), pointer     :: TIRRQ => null()
    real, dimension(:,:,:), pointer     :: CDOMABSQ => null()
    real, dimension(:,:,:), pointer     :: CDET => null()
    real, dimension(:,:,:), pointer     :: PIC => null()
    real, dimension(:,:,:), pointer     :: CDC => null()
    real, dimension(:,:,:), pointer     :: DH => null()

    logical :: IS_MIDNIGHT
    integer :: YY,DOY
    real    :: DT
    real*8  :: DT8
    integer :: NumDBands, NumSBands
    integer :: i,j,n,nl
    integer :: IM, JM, LM
    real    :: hr,rday,daycor,sunz
    real    :: slp,wspd,ozone,wvapor,relhum
    real    :: cov,cldtau,clwp,cldre
    real, dimension(nlt) :: Ed,Es
    real, dimension(nlt) :: ta,wa,asym
    real, dimension(nlt) :: rod,ros

    real, parameter                 :: KUVR = 0.09
    integer                         :: L
    type (ESMF_Time)                :: CurrentTime
    character(len=ESMF_MAXSTR)      :: DATAFILE
    real, pointer, dimension(:,:)   :: FR => null()
    real, pointer, dimension(:,:)   :: PRUVR => null()
    real, pointer, dimension(:,:)   :: PRPAR => null()
    real, pointer, dimension(:,:)   :: PRUVF => null()
    real, pointer, dimension(:,:)   :: PRPAF => null()
    real, pointer, dimension(:,:,:) :: QSW => null()
    real, pointer, dimension(:,:)   :: KPAR => null()
    real, pointer, dimension(:,:)   :: KPARX => null()
    real, pointer, dimension(:,:)   :: Z => null()
    real, pointer, dimension(:,:)   :: UVR => null()
    real, pointer, dimension(:,:)   :: PAR => null()
    real, pointer, dimension(:,:)   :: PS => null()
    real, pointer, dimension(:,:)   :: WSM => null()
    real, pointer, dimension(:,:)   :: OZ => null()
    real, pointer, dimension(:,:)   :: WV => null()
    real, pointer, dimension(:,:)   :: RH => null()
    real, pointer, dimension(:,:)   :: TAUA_T => null()
    real, pointer, dimension(:,:)   :: SSALB_T => null()
    real, pointer, dimension(:,:)   :: ASYMP_T => null()
    type(bandptr), dimension(nlt) :: TAUA
    type(bandptr), dimension(nlt) :: SSALB
    type(bandptr), dimension(nlt) :: ASYMP
    real, pointer, dimension(:,:)   :: CCOV => null()
    real, pointer, dimension(:,:)   :: CLDTC => null()
    real, pointer, dimension(:,:)   :: RLWP => null()
    real, pointer, dimension(:,:)   :: CDRE => null()
    real, pointer, dimension(:,:)   :: data_kpar => null()

!=============================================================================

! Begin... 

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
    Iam = trim(COMP_NAME) // "Run"

! Get the MAPL object from the GC
!--------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS )

! Get parameters from MAPL object
!--------------------------------

    call MAPL_Get(MAPL,     &
         LATS      = LATS,  &
         LONS      = LONS,  &
         ORBIT     = ORBIT, &
         RC=STATUS )
    VERIFY_(STATUS)

! Start Total timer
!------------------

    call MAPL_TimerOn(MAPL,"TOTAL")
    call MAPL_TimerOn(MAPL,"RUN")

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

    call ESMF_UserCompGetInternalState(GC, 'ORADBIO_state', WRAP, STATUS); VERIFY_(STATUS)
    State => WRAP%PTR

! Get Time step
!--------------

    call ESMF_ClockGet       (Clock, TimeStep=DTint, RC=STATUS)
    call ESMF_TimeIntervalGet(DTint, s_r8=DT8,       RC=STATUS)

    DT = DT8

!==== this section was copied from GEOS_OradGridComp.F90 =====
!=============================================================

    call MAPL_GetPointer(IMPORT, FR,    'FROCEAN', RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, PRUVF, 'PENUVF',  RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, PRPAF, 'PENPAF',  RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, PRUVR, 'PENUVR',  RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, PRPAR, 'PENPAR',  RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, DH,    'DH',      RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_GetPointer(EXPORT, QSW,   'SWHEAT', RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, KPARX, 'KPAR',   RC=STATUS)
    VERIFY_(STATUS)

    if( associated(KPARX) ) then
      KPAR => KPARX
    else
      allocate(KPAR(IM,JM), __STAT__)
    end if
    allocate(   Z(IM,JM), __STAT__)
    allocate(UVR (IM,JM), __STAT__)
    allocate(PAR (IM,JM), __STAT__)

! Get current time from clock
!----------------------------

    call ESMF_ClockGet(CLOCK, currTime=CurrentTime, RC=STATUS)
    call ESMF_TimeGet (CurrentTime, YY=YY, DayOfYear=DOY, RC=STATUS)

! Get KPAR from data file
!------------------------

    if (.not. ocean_extData) then
      call MAPL_GetResource(MAPL,DATAFILE,LABEL="KPAR_FILE:"     , RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_ReadForcing(MAPL,'KPAR',DATAFILE,CURRENTTIME,KPAR, RC=STATUS)
      VERIFY_(STATUS)
    else
      call MAPL_GetPointer(import, data_kpar, 'DATA_KPAR', __RC__)
      KPAR = data_kpar
    end if

! Use Beer'S Law to compute flux divergence
!------------------------------------------
   
    UVR = (PRUVF+PRUVR)*FR
    PAR = (PRPAF+PRPAR)*FR
    Z   = 0.0

    if ( associated(QSW) ) then
      QSW(:,:,1) = UVR + PAR

      do L=2,LM   
        Z = Z + DH(:,:,L-1)
        QSW(:,:,L  ) = UVR*exp(-KUVR*Z) + PAR*exp(-KPAR*Z)
        QSW(:,:,L-1) = QSW(:,:,L-1) - QSW(:,:,L)
      enddo
      Z = Z + DH(:,:,LM)
      QSW(:,:,LM) = QSW(:,:,LM) - (UVR*exp(-KUVR*Z) + PAR*exp(-KPAR*Z))
    end if

    if( .not. associated(KPARX) ) deallocate(KPAR)
    deallocate(Z   )
    deallocate(PAR )
    deallocate(UVR )

! ==== end of copied section from Orad ======
!============================================

! Pointers to exports
!--------------------

    call MAPL_GetPointer(EXPORT,TIRRQ , 'TIRRQ' , RC=STATUS)
    VERIFY_(STATUS)
    ASSERT_(associated(TIRRQ ))
    call MAPL_GetPointer(EXPORT,CDOMABSQ , 'CDOMABSQ' , RC=STATUS)
    VERIFY_(STATUS)
    ASSERT_(associated(CDOMABSQ ))

! Pointers to imports
!--------------------
   
    call MAPL_GetPointer(IMPORT, DIATOM,    'DIATOM',    RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, CHLORO,    'CHLORO',    RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, CYANO,     'CYANO',     RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, COCCO,     'COCCO',     RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, DINO,      'DINO',      RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, PHAEO,     'PHAEO',     RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, AVGQ,      'AVGQ',      RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, CDET,      'CDET',      RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, PIC,       'PIC',       RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, CDC,       'CDC',       RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, PS,        'PS',        RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, WSM,       'UU',        RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, OZ,        'OZ',        RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, WV,        'WV',        RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, RH,        'RH',        RC=STATUS)
    VERIFY_(STATUS)
    do nl = 1,nlt
      write(unit = suffix, fmt = '(i2.2)') nl
      call MAPL_GetPointer(IMPORT, TAUA_T, 'TAUA_'//suffix,  RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT, SSALB_T,'SSALB_'//suffix, RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT, ASYMP_T,'ASYMP_'//suffix, RC=STATUS)
      VERIFY_(STATUS)
      TAUA(nl)%b => TAUA_T
      SSALB(nl)%b => SSALB_T
      ASYMP(nl)%b => ASYMP_T
    enddo

    call MAPL_GetPointer(IMPORT, CCOV,      'CCOVM',     RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, CLDTC,     'CLDTCM',    RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, RLWP,      'RLWPM',     RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, CDRE,      'CDREM',     RC=STATUS)
    VERIFY_(STATUS)
!
! Get solar zenith angle
!-----------------------

    allocate(COSZ(IM,JM),__STAT__)
    allocate( SLR(IM,JM),__STAT__)

   call MAPL_SunGetInsolation(LONS, LATS, &
              ORBIT,COSZ, SLR,            &
              CLOCK = CLOCK,              &
              RC=STATUS )
    VERIFY_(STATUS)

! Alarm for midnight
!-------------------

   IS_MIDNIGHT = ESMF_AlarmIsRinging(State%DayAlarm, RC=STATUS)
    VERIFY_(STATUS)

! Allocate 4D array for phytoplankton
!----------------------------------------------

    allocate(PHYTO(LM,ntyp), __STAT__)

! Invoke radiation code
!----------------------
   
    if ( associated(TIRRQ) ) TIRRQ = 0.0
    if ( associated(CDOMABSQ) ) CDOMABSQ = 0.0
hr=1.0
! Obtain Earth-Sun distance 
    rday = float(DOY) + hr/24.0
    daycor = (1.0+1.67E-2*cos(State%pi2*(rday-3.0)/365.0))**2
!   if (Is_Leap)daycor = (1.0+1.67E-2*cos(State%pi2*(rday-3.0)/366.0))**2
!
    do j = 1, JM
     do i = 1, IM
      if (DH(i,j,1) < 1.0E10 .and. COSZ(i,j) > 0.0)then
       slp = PS(i,j)*0.01  ! convert from Pa to mbar
       wspd = WSM(i,j)
       ozone = OZ(i,j)
       wvapor = WV(i,j)
       relhum = RH(i,j)
       do nl = 1,nlt
        ta(nl) = TAUA(nl)%b(i,j)
        wa(nl) = SSALB(nl)%b(i,j)
        asym(nl) = ASYMP(nl)%b(i,j)
       enddo
       cov = CCOV(i,j)
       cldtau = CLDTC(i,j)
       clwp = RLWP(i,j)
       cldre = CDRE(i,j)

!  There are mismatches between the ocean, land and atmosphere in GEOS-5
!  and live ocean points for MOM that do not have a corresponding
!  atmosphere.  These are called "grottoes" because they are assumed
!  to have ocean underneath with land overhead.  Set irradiance to 0
!  to represent this condition.
!       if (slp < 0.0 .or. slp >1.0E10)then ! original line, EMS
       if ((slp < 0.0 .or. slp >1.0E10) .or. (wspd == MAPL_UNDEF) .or. (ozone == MAPL_UNDEF) &
           .or. (wvapor == MAPL_UNDEF) .or. (relhum == MAPL_UNDEF)) then ! extended to include additional variables, EMS
        Ed = 0.0
        Es = 0.0
        TIRRQ(i,j,:) = 0.0
        CDOMABSQ(i,j,:) = 0.0
        AVGQ(i,j,:) = 0.0
       else
!   Spectral irradiance just above surface
        call sfcirr(State%lam,State%Fobar,State%thray,State%oza,       &
                    State%awv,State%ao,State%aco2,                     &
                    State%asl,State%bsl,State%csl,State%dsl,           &
                    State%esl,State%fsl,State%ica,                     &
                    daycor,COSZ(i,j),                                  &
                    slp,wspd,ozone,wvapor,relhum,                      &
                    ta,wa,asym,State%am,State%Vi,                      &
                    cov,cldtau,clwp,cldre,Ed,Es)

!   Spectral irradiance just below surface
        sunz = acos(COSZ(i,j))*State%rad
        call ocalbedo(State%rad,State%lam,State%aw,State%bw,State%wfac, &
         sunz,wspd,rod,ros)
        do nl = 1,nlt
         Ed(nl) = Ed(nl)*(1.0-rod(nl))
         Es(nl) = Es(nl)*(1.0-ros(nl))
        enddo

        PHYTO(:,1) = DIATOM(i,j,:)
        PHYTO(:,2) = CHLORO(i,j,:)
        PHYTO(:,3) = CYANO(i,j,:)
        PHYTO(:,4) = COCCO(i,j,:)
        PHYTO(:,5) = DINO(i,j,:)
        PHYTO(:,6) = PHAEO(i,j,:)
!   Spectral irradiance in the water column
        call glight(LM, IS_MIDNIGHT, COSZ(i,j),                        &
             State%lam,State%aw,State%bw,State%ac,State%bc,State%bpic, &
             State%excdom,State%exdet,State%WtoQ,                      &
             Ed,Es,DH(i,j,:),PHYTO,CDET(i,j,:),PIC(i,j,:),CDC(i,j,:),  &
             TIRRQ(i,j,:),CDOMABSQ(i,j,:),AVGQ(i,j,:),DT)
!TIRRQ(i,j,:) = 1.0
!CDOMABSQ(i,j,:) = 0.0
       endif
      endif
     enddo
    enddo

!  Clean-up
!----------

    if( .not. associated(KPARX) ) deallocate(KPAR)

    deallocate(PHYTO, __STAT__)
    deallocate(COSZ,  __STAT__)
    deallocate(SLR,   __STAT__)

!  All done
!-----------

    call MAPL_TimerOff(MAPL,"RUN" )
    call MAPL_TimerOff(MAPL,"TOTAL")

    RETURN_(ESMF_SUCCESS)

  end subroutine RUN

end module GEOS_OradBioGridCompMod
