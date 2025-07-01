module FireDataBaseType

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! module for handling of fire data
  !
  ! !USES:
  use shr_kind_mod                       , only : r8 => shr_kind_r8, CL => shr_kind_CL
  use shr_log_mod                        , only : errMsg => shr_log_errMsg
  use clm_varctl                         , only : iulog
  use fileutils                          , only : getavu, relavu
  use abortutils                         , only : endrun
  use decompMod                          , only : bounds_type
  use FireMethodType                     , only : fire_method_type
  !
  implicit none
  private
  !
  ! !PUBLIC TYPES:
  public :: fire_base_type

  !
  type, abstract, extends(fire_method_type) :: fire_base_type
    private
      ! !PRIVATE MEMBER DATA:

      real(r8), public, pointer :: forc_lnfm(:)        ! Lightning frequency
      real(r8), public, pointer :: forc_hdm(:)         ! Human population density

    contains
      !
      ! !PUBLIC MEMBER FUNCTIONS:
      procedure, public :: FireInit => BaseFireInit     ! Initialization of Fire
      procedure, public :: BaseFireInit                 ! Initialization of Fire
      procedure(FireReadNML_interface), public, deferred :: FireReadNML       ! Read in namelist for Fire
      procedure(need_lightning_and_popdens_interface), public, deferred :: &
           need_lightning_and_popdens ! Returns true if need lightning & popdens
      !
  end type fire_base_type

  !-------------------------------------------------------------------------

  abstract interface
     !-----------------------------------------------------------------------
     function need_lightning_and_popdens_interface(this) result(need_lightning_and_popdens)
       !
       ! !DESCRIPTION:
       ! Returns true if need lightning and popdens, false otherwise
       !
       ! USES
       import :: fire_base_type
       !
       ! !ARGUMENTS:
       class(fire_base_type), intent(in) :: this
       logical :: need_lightning_and_popdens  ! function result
       !-----------------------------------------------------------------------
     end function need_lightning_and_popdens_interface
  end interface

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains

  !-----------------------------------------------------------------------
  subroutine FireReadNML_interface( this, NLFilename )
    !
    ! !DESCRIPTION:
    ! Read the namelist for Fire
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(fire_base_type) :: this
    character(len=*), intent(in) :: NLFilename ! Namelist filename
  end subroutine FireReadNML_interface

  !-----------------------------------------------------------------------
  subroutine BaseFireInit( this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize CN Fire module
    ! !USES:
      use nanMod      , only : nan
    !
    ! !ARGUMENTS:
    class(fire_base_type) :: this
    type(bounds_type), intent(in) :: bounds
    !character(len=*),  intent(in) :: NLFilename
    !-----------------------------------------------------------------------

    if ( this%need_lightning_and_popdens() ) then
       ! Allocate lightning forcing data
       allocate( this%forc_lnfm(bounds%begg:bounds%endg) )
       this%forc_lnfm(bounds%begg:) = nan
       ! Allocate pop dens forcing data
       allocate( this%forc_hdm(bounds%begg:bounds%endg) )
       this%forc_hdm(bounds%begg:) = nan
    end if

  end subroutine BaseFireInit


end module FireDataBaseType
