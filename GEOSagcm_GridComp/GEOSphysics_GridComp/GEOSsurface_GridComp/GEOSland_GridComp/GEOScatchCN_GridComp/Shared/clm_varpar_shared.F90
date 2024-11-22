module clm_varpar_shared

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: clm_varpar_shared
!
! !DESCRIPTION:
! Module containing CLM parameters
!
! !USES:
!
! !PUBLIC TYPES:
  implicit none
  save
!
! Define number of levels

  integer, parameter :: numpft_CN         = 19   ! actual # of pfts (without bare) for Catchment-CN4.0 and Catchment-CN4.5
  integer, parameter :: numpft_CN51       = 15   ! actual # of pfts (without bare) for Catchment-CN5.1

  integer, parameter, PUBLIC :: NUM_ZON_CN=3   ! number of CN hydrology zones per tile
  integer, parameter, PUBLIC :: NUM_VEG_CN=4   ! number of CN PFTs per zone for Catchment-CN4.0 and Catchment-CN4.5
  integer, parameter, PUBLIC :: NUM_VEG_CN51 = 2 ! number of CN PFTs per zone for Catchment-CN5.1
  integer, parameter, PUBLIC :: VAR_COL_40=40 ! number of CN column restart variables
  integer, parameter, PUBLIC :: VAR_PFT_40=74 ! number of CN PFT variables per column 
  integer, parameter, PUBLIC :: VAR_COL_45=35  ! number of CN column restart variables
  integer, parameter, PUBLIC :: VAR_PFT_45=75  ! number of CN PFT variables per column  
  integer, parameter, PUBLIC :: VAR_COL_51=37  ! number of CN column restart variables
  integer, parameter, PUBLIC :: VAR_PFT_51=83  ! number of CN PFT restart variables  


!------------------------------------------------------------------------------
end module clm_varpar_shared
