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
! fzeng: these are CLM types
!
  integer, parameter :: noveg                 =  0   ! Bare                                 
  integer, parameter :: ndllf_evr_tmp_tree    =  1   ! Needleleaf evergreen temperate tree           
  integer, parameter :: ndllf_evr_brl_tree    =  2   ! Needleleaf evergreen boreal tree             
  integer, parameter :: ndllf_dcd_brl_tree    =  3   ! Needleleaf deciduous boreal tree              
  integer, parameter :: nbrdlf_evr_trp_tree   =  4   ! Broadleaf evergreen tropical tree            
  integer, parameter :: nbrdlf_evr_tmp_tree   =  5   ! Broadleaf evergreen temperate tree            
  integer, parameter :: nbrdlf_dcd_trp_tree   =  6   ! Broadleaf deciduous tropical tree             
  integer, parameter :: nbrdlf_dcd_tmp_tree   =  7   ! Broadleaf deciduous temperate tree            
  integer, parameter :: nbrdlf_dcd_brl_tree   =  8   ! Broadleaf deciduous boreal tree               
  integer, parameter :: nbrdlf_evr_shrub      =  9   ! Broadleaf evergreen temperate shrub                     
  integer, parameter :: nbrdlf_dcd_tmp_shrub  = 10   ! Broadleaf deciduous temperate shrub [moisture + deciduous]
  integer, parameter :: nbrdlf_dcd_tmp_shrub2 = 11   ! Broadleaf deciduous temperate shrub [moisture stress only]          
  integer, parameter :: nbrdlf_dcd_brl_shrub  = 12   ! Broadleaf deciduous boreal shrub            
  integer, parameter :: nc3_arctic_grass      = 13   ! Arctic c3 grass                              
  integer, parameter :: nc3_nonarctic_grass   = 14   ! Cool c3 grass [moisture + deciduous]
  integer, parameter :: nc3_nonarctic_grass2  = 15   ! Cool c3 grass [moisture stress only]                       
  integer, parameter :: nc4_grass             = 16   ! Warm c4 grass [moisture + deciduous]
  integer, parameter :: nc4_grass2            = 17   ! Warm c4 grass [moisture stress only]                                  
  integer, parameter :: nc3crop               = 18   ! C3_crop [moisture + deciduous]                                      
  integer, parameter :: nc3crop2              = 19   ! C3_crop [moisture stress only]                                  
! Although we don't use prognostic crops, we need to keep these below because the CN code needs npcropmin, fzeng, 10 May 2018
  integer, parameter :: ncorn                 = 20   ! Corn                                          
  integer, parameter :: ncornirrig            = 21   ! Irrigated corn                                
  integer, parameter :: nscereal              = 22   ! Spring temperate cereal                       
  integer, parameter :: nscerealirrig         = 23   ! Irrigated spring temperate cereal            
  integer, parameter :: nwcereal              = 24   ! winter temperate cereal                       
  integer, parameter :: nwcerealirrig         = 25   ! Irrigated winter temperate cereal            
  integer, parameter :: nsoybean              = 26   ! Soybean                                       
  integer, parameter :: nsoybeanirrig         = 27   ! Irrigated Soybean                             
  
  integer, parameter :: ntree                 = nbrdlf_dcd_brl_tree  !value for last type of tree
  integer, parameter :: npcropmin             = ncorn                ! first prognostic crop
  integer, parameter :: npcropmax             = nsoybeanirrig        ! last prognostic crop in list

end module pftvarcon

