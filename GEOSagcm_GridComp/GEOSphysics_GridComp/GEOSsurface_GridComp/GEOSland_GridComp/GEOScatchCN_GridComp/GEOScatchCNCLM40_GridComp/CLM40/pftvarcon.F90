module pftvarcon

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: pftvarcon
!
! !DESCRIPTION:
! Module containing vegetation constants and method to
! read and initialize vegetation (PFT) constants.
!
! !USES:
!
! !PUBLIC TYPES:
  implicit none
  save
!
! Vegetation type constants
!
! gkw: these are CLM types
!
  integer, parameter :: noveg                =  0  !value for not vegetated 
  integer, parameter :: ndllf_evr_tmp_tree   =  1  !value for Needleleaf evergreen temperate tree
  integer, parameter :: ndllf_evr_brl_tree   =  2  !value for Needleleaf evergreen boreal tree
  integer, parameter :: ndllf_dcd_brl_tree   =  3  !value for Needleleaf deciduous boreal tree
  integer, parameter :: nbrdlf_evr_trp_tree  =  4  !value for Broadleaf evergreen tropical tree
  integer, parameter :: nbrdlf_evr_tmp_tree  =  5  !value for Broadleaf evergreen temperate tree
  integer, parameter :: nbrdlf_dcd_trp_tree  =  6  !value for Broadleaf deciduous tropical tree
  integer, parameter :: nbrdlf_dcd_tmp_tree  =  7  !value for Broadleaf deciduous temperate tree
  integer, parameter :: nbrdlf_dcd_brl_tree  =  8  !value for Broadleaf deciduous boreal tree
  integer, parameter :: ntree                =  nbrdlf_dcd_brl_tree  !value for last type of tree
  integer, parameter :: nbrdlf_evr_shrub     =  9  !value for Broadleaf evergreen shrub
  integer, parameter :: nbrdlf_dcd_tmp_shrub = 10  !value for Broadleaf deciduous temperate shrub
  integer, parameter :: nbrdlf_dcd_tmp_shrub2= 11  !value for Broadleaf deciduous temperate shrub
  integer, parameter :: nbrdlf_dcd_brl_shrub = 12  !value for Broadleaf deciduous boreal shrub
  integer, parameter :: nc3_arctic_grass     = 13  !value for C3 arctic grass
  integer, parameter :: nc3_nonarctic_grass  = 14  !value for C3 non-arctic grass
  integer, parameter :: nc3_nonarctic_grass2 = 15  !value for C3 non-arctic grass
  integer, parameter :: nc4_grass            = 16  !value for C4 grass
  integer, parameter :: nc4_grass2           = 17  !value for C4 grass
  integer, parameter :: ncrop                = 18  !value for crop
  integer, parameter :: ncrop2               = 19  !value for crop

end module pftvarcon

