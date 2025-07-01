module initVerticalMod

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Initialize vertical components of column datatype
  !
  use shr_kind_mod      , only : r8 => shr_kind_r8
  use shr_infnan_mod    , only : nan => shr_infnan_nan
  use shr_log_mod       , only : errMsg => shr_log_errMsg
  use shr_sys_mod       , only : shr_sys_abort
  use decompMod         , only : bounds_type
  use spmdMod           , only : masterproc
  use clm_varpar        , only : nlevsno, nlevgrnd, nlevlak
  use clm_varpar        , only : nlevsoi, nlevurb, nlevmaxurbgrnd
  use clm_varctl        , only : iulog
  use clm_varctl        , only : use_vertsoilc
  use clm_varctl        , only : use_fates
  use clm_varcon        , only : zsoi, dzsoi, zisoi, dzsoi_decomp, spval, ispval
  use fileutils         , only : getfil
  use LandunitType      , only : lun                
  use GridcellType      , only : grc                
  use ColumnType        , only : col                          
  use abortUtils        , only : endrun    
  use ncdio_pio
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
 ! public :: initVertical
  public :: find_soil_layer_containing_depth

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

  !
  !------------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine find_soil_layer_containing_depth(depth, layer)
    !
    ! !DESCRIPTION:
    ! Find the soil layer that contains the given depth
    !
    ! Aborts if the given depth doesn't exist in the soil profile
    !
    ! We consider the interface between two layers to belong to the layer *above* that
    ! interface. This implies that the top interface (at exactly 0 m) is not considered
    ! to be part of the soil profile.
    !
    ! !ARGUMENTS:
    real(r8), intent(in)  :: depth  ! target depth, m
    integer , intent(out) :: layer  ! layer containing target depth
    !
    ! !LOCAL VARIABLES:
    logical :: found
    integer :: i

    character(len=*), parameter :: subname = 'find_soil_layer_containing_depth'
    !-----------------------------------------------------------------------

    if (depth <= zisoi(0)) then
       write(iulog,*) subname, ': ERROR: depth above top of soil'
       write(iulog,*) 'depth = ', depth
       write(iulog,*) 'zisoi = ', zisoi
       call endrun(msg=subname//': depth above top of soil')
    end if

    found = .false.
    do i = 1, nlevgrnd
       if (depth <= zisoi(i)) then
          layer = i
          found = .true.
          exit
       end if
    end do

    if (.not. found) then
       write(iulog,*) subname, ': ERROR: depth below bottom of soil'
       write(iulog,*) 'depth = ', depth
       write(iulog,*) 'zisoi = ', zisoi
       call endrun(msg=subname//': depth below bottom of soil')
    end if

  end subroutine find_soil_layer_containing_depth

end module initVerticalMod
