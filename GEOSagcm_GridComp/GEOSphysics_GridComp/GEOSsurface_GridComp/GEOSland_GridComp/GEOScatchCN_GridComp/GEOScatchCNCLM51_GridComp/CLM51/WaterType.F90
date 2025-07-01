module WaterType

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Container for derived types relating to water, both for bulk water and for isotopes
  ! and other tracers.
  !
  ! Variables pertaining to bulk water can be accessed in two ways:
  !
  !    (1) Using water_inst%water*bulk_inst
  !
  !    (2) As one of the indices in water_inst%bulk_and_tracers(:)%water*_inst
  !
  !    Method (1) is greatly preferable when you are just operating on bulk water. Method
  !    (2) is just meant to be used when you are doing the same operation on bulk water
  !    and all water tracers.
  !
  ! To loop through bulk and all tracers, use code like this:
  !    do i = water_inst%bulk_and_tracers_beg, water_inst%bulk_and_tracers_end
  !       associate( &
  !            waterflux_inst => water_inst%bulk_and_tracers(i)%waterflux_inst, &
  !            [and other associations, as necessary])
  !       [Do calculations involving waterflux_inst, etc.]
  !       end associate
  !    end do
  !
  ! To loop through all tracers (not bulk), use code like this:
  !    do i = water_inst%tracers_beg, water_inst%tracers_end
  !       associate( &
  !            waterflux_inst => water_inst%bulk_and_tracers(i)%waterflux_inst, &
  !            [and other associations, as necessary])
  !       [Do calculations involving waterflux_inst, etc.]
  !       end associate
  !    end do
  !
  ! To loop through all isotopes (not bulk or other water tracers), use code like this:
  !    type(water_info_isotope_type), pointer :: iso_info
  !
  !    do i = water_inst%tracers_beg, water_inst%tracers_end
  !       if (water_inst%IsIsotope(i)) then
  !          call water_inst%GetIsotopeInfo(i, iso_info)
  !          associate( &
  !               waterflux_inst => water_inst%bulk_and_tracers(i)%waterflux_inst, &
  !               [and other associations, as necessary])
  !          [Do calculations involving iso_info, waterflux_inst, etc.]
  !          end associate
  !       end if
  !    end do
  !
  ! The associate statements given above aren't crucial. If the block of code refers to
  ! multiple instances (waterstate, waterflux, etc.), but only refers to each one once or
  ! twice, it can be best to just have:
  !    associate(bulk_or_tracer => water_inst%bulk_and_tracers(i))
  !
  ! !USES:
#include "shr_assert.h"
  use shr_kind_mod             , only : r8 => shr_kind_r8
  use shr_log_mod              , only : errMsg => shr_log_errMsg
  use abortutils               , only : endrun
  use decompMod                , only : bounds_type
  use clm_varctl               , only : iulog
  use clm_varpar               , only : nlevsno
  use ncdio_pio                , only : file_desc_t
  use WaterFluxBulkType        , only : waterfluxbulk_type
  use WaterFluxType            , only : waterflux_type
  use WaterStateBulkType       , only : waterstatebulk_type
  use WaterStateType           , only : waterstate_type
  use WaterDiagnosticType      , only : waterdiagnostic_type
  use WaterDiagnosticBulkType  , only : waterdiagnosticbulk_type
!  use WaterBalanceType         , only : waterbalance_type
!  use WaterInfoBaseType        , only : water_info_base_type
!  use WaterInfoBulkType        , only : water_info_bulk_type
!  use WaterInfoTracerType      , only : water_info_tracer_type
!  use WaterInfoIsotopeType     , only : water_info_isotope_type
!  use Waterlnd2atmType         , only : waterlnd2atm_type
!  use Waterlnd2atmBulkType     , only : waterlnd2atmbulk_type
  use Wateratm2lndType         , only : wateratm2lnd_type
  use Wateratm2lndBulkType     , only : wateratm2lndbulk_type
 ! use WaterTracerContainerType , only : water_tracer_container_type
 ! use WaterTracerUtils         , only : CompareBulkToTracer, SetTracerToBulkTimesRatio

  implicit none
  private

  !
  ! !PRIVATE TYPES:

  ! This type holds instances needed for bulk water or for a single tracer
  type, private :: bulk_or_tracer_type
     private

     ! ------------------------------------------------------------------------
     ! Public data members
     ! ------------------------------------------------------------------------

     class(waterflux_type)       , pointer, public :: waterflux_inst
     class(waterstate_type)      , pointer, public :: waterstate_inst
     class(waterdiagnostic_type) , pointer, public :: waterdiagnostic_inst
    ! class(waterbalance_type)    , pointer, public :: waterbalance_inst
    ! class(waterlnd2atm_type)    , pointer, public :: waterlnd2atm_inst
     class(wateratm2lnd_type)    , pointer, public :: wateratm2lnd_inst

     ! ------------------------------------------------------------------------
     ! Private data members
     ! ------------------------------------------------------------------------

  !   logical :: is_isotope = .false.
  !   class(water_info_base_type) , pointer :: info
  !   type(water_tracer_container_type) :: vars

  end type bulk_or_tracer_type

  !
  ! !PUBLIC TYPES:

!  ! water_params_type is public for the sake of unit tests
!  type, public :: water_params_type
!     private
!
!     ! Whether we add tracers that are used for the tracer consistency checks
!     logical :: enable_consistency_checks
!
!     ! Whether we add tracers that are used for isotopes
!     logical :: enable_isotopes
!  end type water_params_type

  type, public :: water_type
     private

     ! ------------------------------------------------------------------------
     ! Public data members
     ! ------------------------------------------------------------------------

     ! indices into the bulk_and_tracers array
     integer, public :: bulk_and_tracers_beg  ! first index when looping over bulk & tracers
     integer, public :: bulk_and_tracers_end  ! last index when looping over bulk & tracers
     integer, public :: tracers_beg           ! first index when looping over just tracers
     integer, public :: tracers_end           ! last index when looping over just tracers
     integer, public :: i_bulk                ! index of bulk in arrays that contain both bulk and tracers

     type(waterfluxbulk_type), pointer, public       :: waterfluxbulk_inst
     type(waterstatebulk_type), pointer, public      :: waterstatebulk_inst
     type(waterdiagnosticbulk_type), pointer, public :: waterdiagnosticbulk_inst
    ! type(waterbalance_type), pointer, public        :: waterbalancebulk_inst
    ! type(waterlnd2atmbulk_type), pointer, public    :: waterlnd2atmbulk_inst
     type(wateratm2lndbulk_type), pointer, public    :: wateratm2lndbulk_inst

     type(bulk_or_tracer_type),  public  :: bulk_and_tracers(1)

     ! ------------------------------------------------------------------------
     ! Private data members
     ! ------------------------------------------------------------------------

  !   type(water_params_type) :: params
     integer :: bulk_tracer_index  ! index of the tracer that replicates bulk water (-1 if it doesn't exist)

   contains
     ! Public routines
     procedure, public :: Init

  end type water_type
  type(water_type), public, target, save :: water_inst

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains


  !-----------------------------------------------------------------------
  subroutine Init(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize all water variables
    !
    ! !ARGUMENTS:
    class(water_type), intent(inout) :: this
    type(bounds_type), intent(in)    :: bounds
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'Init'
    !-----------------------------------------------------------------------

    allocate(this%waterfluxbulk_inst)
    this%bulk_and_tracers(1)%waterflux_inst => this%waterfluxbulk_inst

    allocate(this%waterstatebulk_inst)
    this%bulk_and_tracers(1)%waterstate_inst => this%waterstatebulk_inst

    allocate(this%waterdiagnosticbulk_inst)
    this%bulk_and_tracers(1)%waterdiagnostic_inst => this%waterdiagnosticbulk_inst

    allocate(this%wateratm2lndbulk_inst)
    this%bulk_and_tracers(1)%wateratm2lnd_inst => this%wateratm2lndbulk_inst

    allocate(waterflux_type :: this%bulk_and_tracers(1)%waterflux_inst)
    allocate(waterdiagnostic_type :: this%bulk_and_tracers(1)%waterdiagnostic_inst)
    allocate(waterstate_type :: this%bulk_and_tracers(1)%waterstate_inst)
    allocate(wateratm2lnd_type :: this%bulk_and_tracers(1)%wateratm2lnd_inst)

    call this%bulk_and_tracers(1)%waterflux_inst%Init            (bounds)
    call this%bulk_and_tracers(1)%wateratm2lnd_inst%Init         (bounds)
    call this%bulk_and_tracers(1)%waterstate_inst%Init           (bounds)
    call this%waterfluxbulk_inst%InitBulk                        (bounds)
    call this%waterdiagnosticbulk_inst%InitBulk                  (bounds)
    call this%wateratm2lndbulk_inst%InitBulk                     (bounds)
    call this%waterstatebulk_inst%InitBulk                           (bounds)


  end subroutine Init


end module WaterType
