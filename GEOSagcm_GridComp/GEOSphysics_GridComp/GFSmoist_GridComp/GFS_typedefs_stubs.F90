#undef MULTI_GASES

module GFS_typedefs
       
       implicit none

!----------------------------------------------------------------------------------
! GFS_control_type
!   model control parameters input from a namelist and/or derived from others
!   list of those that can be modified during the run are at the bottom of the list 
!----------------------------------------------------------------------------------
!! \section arg_table_GFS_control_type
!! \htmlinclude GFS_control_type.html
!!
  type GFS_control_type
  end type GFS_control_type

!---------------------------------------------------------------------
! GFS_interstitial_type
!   fields required for interstitial code in CCPP schemes, previously
!   in GFS_{physics,radiation}_driver.F90
!---------------------------------------------------------------------
!! \section arg_table_GFS_interstitial_type
!! \htmlinclude GFS_interstitial_type.html
!!
  type GFS_interstitial_type
    contains
      procedure :: rad_reset   => interstitial_rad_reset  !<   reset array data for radiation
      procedure :: phys_reset  => interstitial_phys_reset !<   reset array data for physics
  end type GFS_interstitial_type

!----------------
! PUBLIC ENTITIES
!----------------
  public GFS_control_type
  public GFS_interstitial_type

!*******************************************************************************************
  CONTAINS

  subroutine interstitial_rad_reset (Interstitial, Model)
    !
    implicit none
    !
    class(GFS_interstitial_type) :: Interstitial
    type(GFS_control_type), intent(in) :: Model
  end subroutine interstitial_rad_reset

  subroutine interstitial_phys_reset (Interstitial, Model)
    !
    implicit none
    !
    class(GFS_interstitial_type) :: Interstitial
    type(GFS_control_type), intent(in) :: Model

  end subroutine interstitial_phys_reset

end module GFS_typedefs
