module dynSubgridControlMod

  !---------------------------------------------------------------------------
  !
  ! !DESCRIPTION:
  ! Defines a class for storing and querying control flags related to dynamic subgrid
  ! operation.
  !
  ! Note that this is implemented (essentially) as a singleton, so the only instance of
  ! this class is stored in this module. This is done for convenience, to avoid having to
  ! pass around the single instance just to query these control flags.
  !
  ! !USES:
#include "shr_assert.h"
  use shr_log_mod        , only : errMsg => shr_log_errMsg
  use abortutils         , only : endrun
  !
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: dynSubgridControl_init
  public :: get_do_transient_pfts   ! return the value of the do_transient_pfts control flag
  public :: get_do_transient_crops  ! return the value of the do_transient_crops control flag
  public :: get_do_harvest          ! return the value of the do_harvest control flag
  public :: run_has_transient_landcover ! returns true if any aspects of prescribed transient landcover are enabled
  !
  ! !PRIVATE TYPES:
  type dyn_subgrid_control_type
     private
     logical :: do_transient_pfts  = .false. ! whether to apply transient natural PFTs from dataset
     logical :: do_transient_crops = .false. ! whether to apply transient crops from dataset
     logical :: do_transient_lakes = .false. ! whether to apply transient lakes from dataset 
     logical :: do_harvest         = .false. ! whether to apply harvest from dataset

     logical :: reset_dynbal_baselines = .false. ! whether to reset baseline values of total column water and energy in the first step of the run

     ! The following is only meant for testing: Whether area changes are allowed at times
     ! other than the year boundary. This should only arise in some test configurations
     ! where we artificially create changes more frequently so that we can run short
     ! tests. This flag is only used for error-checking, not controlling any model
     ! behavior.
     logical :: for_testing_allow_non_annual_changes = .false.

     ! The following is only meant for testing: If .true., set the dynbal water and
     ! energy fluxes to zero. This is needed in some tests where we have daily rather
     ! than annual glacier dynamics: if we allow the true dynbal adjustment fluxes in
     ! those tests, we end up with sensible heat fluxes of thousands of W m-2 or more,
     ! which causes CAM to blow up. However, note that setting it to true will break
     ! water and energy conservation!
     logical :: for_testing_zero_dynbal_fluxes = .false.

     logical :: initialized        = .false. ! whether this object has been initialized
  end type dyn_subgrid_control_type
  
  type(dyn_subgrid_control_type) :: dyn_subgrid_control_inst

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains

    !-----------------------------------------------------------------------
  subroutine dynSubgridControl_init(  )
    !
    ! !DESCRIPTION:
    ! Initialize the dyn_subgrid_control settings.
    !
    ! !USES:
    use spmdMod           , only : masterproc
    !
    ! !ARGUMENTS:

    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'dynSubgridControl_init'
    !-----------------------------------------------------------------------

    dyn_subgrid_control_inst%initialized = .true.

  end subroutine dynSubgridControl_init

  !-----------------------------------------------------------------------
  logical function get_do_transient_pfts()
    ! !DESCRIPTION:
    ! Return the value of the do_transient_pfts control flag
    !-----------------------------------------------------------------------
    
    SHR_ASSERT_FL(dyn_subgrid_control_inst%initialized, sourcefile, __LINE__)

    get_do_transient_pfts = dyn_subgrid_control_inst%do_transient_pfts

  end function get_do_transient_pfts

  !-----------------------------------------------------------------------
  logical function get_do_transient_crops()
    ! !DESCRIPTION:
    ! Return the value of the do_transient_crops control flag
    !-----------------------------------------------------------------------
    
    SHR_ASSERT_FL(dyn_subgrid_control_inst%initialized, sourcefile, __LINE__)

    get_do_transient_crops = dyn_subgrid_control_inst%do_transient_crops

  end function get_do_transient_crops
  
  !-----------------------------------------------------------------------
  logical function run_has_transient_landcover()
    ! !DESCRIPTION:
    ! Returns true if any aspects of prescribed transient landcover are enabled
    !-----------------------------------------------------------------------

    run_has_transient_landcover = &
         (get_do_transient_pfts() .or. &
         get_do_transient_crops())
  end function run_has_transient_landcover

  !-----------------------------------------------------------------------

  logical function get_do_harvest()
    ! !DESCRIPTION:
    ! Return the value of the do_harvest control flag
    !-----------------------------------------------------------------------

    SHR_ASSERT_FL(dyn_subgrid_control_inst%initialized, sourcefile, __LINE__)

    get_do_harvest = dyn_subgrid_control_inst%do_harvest

  end function get_do_harvest

end module dynSubgridControlMod
