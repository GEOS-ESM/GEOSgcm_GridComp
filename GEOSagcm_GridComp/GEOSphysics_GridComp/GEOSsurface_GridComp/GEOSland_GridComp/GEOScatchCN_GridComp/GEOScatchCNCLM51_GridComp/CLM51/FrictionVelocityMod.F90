module FrictionVelocityMod

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculation of the friction velocity, relation for potential
  ! temperature and humidity profiles of surface boundary layer.
  !
  ! !USES:
#include "shr_assert.h"
  use shr_kind_mod            , only : r8 => shr_kind_r8
  use nanMod                  , only : nan
  use shr_log_mod             , only : errMsg => shr_log_errMsg
  use shr_const_mod           , only : SHR_CONST_PI
  use decompMod               , only : bounds_type
  use clm_varcon              , only : spval
  use clm_varctl              , only : use_cn, use_luna
  use LandunitType            , only : lun
  use ColumnType              , only : col
  use PatchType               , only : patch
  use landunit_varcon         , only : istsoil, istcrop, istice_mec, istwet
  use ncdio_pio               , only : file_desc_t
  use paramUtilMod            , only : readNcdioScalar
  use atm2lndType             , only : atm2lnd_type
  use WaterDiagnosticBulkType , only : waterdiagnosticbulk_type
  use CanopyStateType         , only : canopystate_type
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  save

! !PUBLIC MEMBER FUNCTIONS:

  type, public :: frictionvel_type
     private

     ! Scalar parameters
     real(r8), public :: zetamaxstable = -999._r8  ! Max value zeta ("height" used in Monin-Obukhov theory) can go to under stable conditions
     real(r8) :: zsno = -999._r8  ! Momentum roughness length for snow (m)
     real(r8) :: zlnd = -999._r8  ! Momentum roughness length for soil, glacier, wetland (m)

     ! Roughness length/resistance for friction velocity calculation

     real(r8), pointer, public :: forc_hgt_u_patch (:)   ! patch wind forcing height (10m+z0m+d) (m)
     real(r8), pointer, public :: forc_hgt_t_patch (:)   ! patch temperature forcing height (10m+z0m+d) (m)
     real(r8), pointer, public :: forc_hgt_q_patch (:)   ! patch specific humidity forcing height (10m+z0m+d) (m)
     real(r8), pointer, public :: u10_patch        (:)   ! patch 10-m wind (m/s) (for dust model)
     real(r8), pointer, public :: u10_clm_patch    (:)   ! patch 10-m wind (m/s) (for clm_map2gcell)
     real(r8), pointer, public :: va_patch         (:)   ! patch atmospheric wind speed plus convective velocity (m/s)
     real(r8), pointer, public :: vds_patch        (:)   ! patch deposition velocity term (m/s) (for dry dep SO4, NH4NO3)
     real(r8), pointer, public :: fv_patch         (:)   ! patch friction velocity (m/s) (for dust model)
     real(r8), pointer, public :: rb1_patch        (:)   ! patch aerodynamical resistance (s/m) (for dry deposition of chemical tracers)
     real(r8), pointer, public :: rb10_patch       (:)   ! 10-day mean patch aerodynamical resistance (s/m) (for LUNA model)
     real(r8), pointer, public :: ram1_patch       (:)   ! patch aerodynamical resistance (s/m)
     real(r8), pointer, public :: z0mv_patch       (:)   ! patch roughness length over vegetation, momentum [m]
     real(r8), pointer, public :: z0hv_patch       (:)   ! patch roughness length over vegetation, sensible heat [m]
     real(r8), pointer, public :: z0qv_patch       (:)   ! patch roughness length over vegetation, latent heat [m]
     real(r8), pointer, public :: z0mg_col         (:)   ! col roughness length over ground, momentum  [m] 
     real(r8), pointer, public :: z0hg_col         (:)   ! col roughness length over ground, sensible heat [m]
     real(r8), pointer, public :: z0qg_col         (:)   ! col roughness length over ground, latent heat [m]
     ! variables to add history output from CanopyFluxesMod
     real(r8), pointer, public :: rah1_patch       (:)   ! patch sensible heat flux resistance [s/m]
     real(r8), pointer, public :: rah2_patch       (:)   ! patch below-canopy sensible heat flux resistance [s/m]
     real(r8), pointer, public :: raw1_patch       (:)   ! patch moisture flux resistance [s/m]
     real(r8), pointer, public :: raw2_patch       (:)   ! patch below-canopy moisture flux resistance [s/m]
     real(r8), pointer, public :: ustar_patch      (:)   ! patch friction velocity [m/s]
     real(r8), pointer, public :: um_patch         (:)   ! patch wind speed including the stablity effect [m/s]
     real(r8), pointer, public :: uaf_patch        (:)   ! patch canopy air speed [m/s]
     real(r8), pointer, public :: taf_patch        (:)   ! patch canopy air temperature [K]
     real(r8), pointer, public :: qaf_patch        (:)   ! patch canopy humidity [kg/kg]
     real(r8), pointer, public :: obu_patch        (:)   ! patch Monin-Obukhov length [m]
     real(r8), pointer, public :: zeta_patch       (:)   ! patch dimensionless stability parameter
     real(r8), pointer, public :: vpd_patch        (:)   ! patch vapor pressure deficit [Pa]
     real(r8), pointer, public :: num_iter_patch   (:)   ! patch number of iterations
     real(r8), pointer, public :: z0m_actual_patch (:)   ! patch roughness length actually used in flux calculations, momentum [m]

   contains

     procedure , public :: Init

  end type frictionvel_type
  type(frictionvel_type), public, target, save :: frictionvel_inst

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !------------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init( this, bounds)

  !  use shr_infnan_mod , only : nan => shr_infnan_nan
    
    type(bounds_type), intent(in) :: bounds
    class(frictionvel_type)       :: this
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc

    allocate(this%forc_hgt_u_patch (begp:endp)) ; this%forc_hgt_u_patch (:)   = nan
    allocate(this%forc_hgt_t_patch (begp:endp)) ; this%forc_hgt_t_patch (:)   = nan
    allocate(this%forc_hgt_q_patch (begp:endp)) ; this%forc_hgt_q_patch (:)   = nan
    allocate(this%u10_patch        (begp:endp)) ; this%u10_patch        (:)   = nan
    allocate(this%u10_clm_patch    (begp:endp)) ; this%u10_clm_patch    (:)   = nan
    allocate(this%va_patch         (begp:endp)) ; this%va_patch         (:)   = nan
    allocate(this%vds_patch        (begp:endp)) ; this%vds_patch        (:)   = nan
    allocate(this%fv_patch         (begp:endp)) ; this%fv_patch         (:)   = nan
    allocate(this%rb1_patch        (begp:endp)) ; this%rb1_patch        (:)   = nan
    allocate(this%rb10_patch       (begp:endp)) ; this%rb10_patch       (:)   = spval
    allocate(this%ram1_patch       (begp:endp)) ; this%ram1_patch       (:)   = nan
    allocate(this%z0mv_patch       (begp:endp)) ; this%z0mv_patch       (:)   = nan
    allocate(this%z0hv_patch       (begp:endp)) ; this%z0hv_patch       (:)   = nan
    allocate(this%z0qv_patch       (begp:endp)) ; this%z0qv_patch       (:)   = nan
    allocate(this%z0mg_col         (begc:endc)) ; this%z0mg_col         (:)   = nan
    allocate(this%z0qg_col         (begc:endc)) ; this%z0qg_col         (:)   = nan
    allocate(this%z0hg_col         (begc:endc)) ; this%z0hg_col         (:)   = nan
    allocate(this%rah1_patch       (begp:endp)) ; this%rah1_patch       (:)   = nan
    allocate(this%rah2_patch       (begp:endp)) ; this%rah2_patch       (:)   = nan
    allocate(this%raw1_patch       (begp:endp)) ; this%raw1_patch       (:)   = nan
    allocate(this%raw2_patch       (begp:endp)) ; this%raw2_patch       (:)   = nan
    allocate(this%um_patch         (begp:endp)) ; this%um_patch         (:)   = nan
    allocate(this%uaf_patch        (begp:endp)) ; this%uaf_patch        (:)   = nan
    allocate(this%taf_patch        (begp:endp)) ; this%taf_patch        (:)   = nan
    allocate(this%qaf_patch        (begp:endp)) ; this%qaf_patch        (:)   = nan
    allocate(this%ustar_patch      (begp:endp)) ; this%ustar_patch      (:)   = nan
    allocate(this%obu_patch        (begp:endp)) ; this%obu_patch        (:)   = nan
    allocate(this%zeta_patch       (begp:endp)) ; this%zeta_patch       (:)   = nan
    allocate(this%vpd_patch        (begp:endp)) ; this%vpd_patch        (:)   = nan
    allocate(this%num_iter_patch   (begp:endp)) ; this%num_iter_patch   (:)   = nan
    allocate(this%z0m_actual_patch (begp:endp)) ; this%z0m_actual_patch (:)   = nan

  end subroutine Init


end module FrictionVelocityMod
