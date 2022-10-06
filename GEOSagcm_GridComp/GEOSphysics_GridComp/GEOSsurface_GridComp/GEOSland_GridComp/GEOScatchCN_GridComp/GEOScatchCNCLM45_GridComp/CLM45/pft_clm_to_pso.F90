module compute_rc_functions
    
  implicit none

  private
  public :: pft_clm_to_pso

  contains
  
    subroutine pft_clm_to_pso(clm_pft,pso_pft)
  
  
      !  map catchment type into PFT
      !  ---------------------------						! MAPPED TO G1 EF PFTS! (by twr)
      !PFT 	Description							! --------------------------- 
      ! 0 	bare 	                                                    ! Bare
      ! 1 	needleleaf evergreen temperate tree				            ! Forest
      ! 2 	needleleaf evergreen boreal tree				            ! Forest
      ! 3 	needleleaf deciduous boreal tree				            ! Forest
      ! 4 	broadleaf evergreen tropical tree				            ! Forest
      ! 5 	broadleaf evergreen temperate tree				            ! Forest
      ! 6 	broadleaf deciduous tropical tree				            ! Forest
      ! 7 	broadleaf deciduous temperate tree				            ! Forest
      ! 8 	broadleaf deciduous boreal tree					            ! Forest
      ! 9     broadleaf evergreen temperate shrub                         ! Shrublands
      ! 10 	broadleaf deciduous temperate shrub [moisture + deciduous]	! Shrublands
      ! 11 	broadleaf deciduous temperate shrub [moisture stress only]	! Shrublands
      ! 12 	broadleaf deciduous boreal shrub				            ! Shrublands
      ! 13 	arctic c3 grass							                    ! Grasslands
      ! 14 	cool c3 grass [moisture + deciduous]				        ! Grasslands
      ! 15 	cool c3 grass [moisture stress only]   				        ! Grasslands
      ! 16 	warm c4 grass [moisture + deciduous]				        ! Grasslands
      ! 17 	warm c4 grass [moisture stress only]				        ! Grasslands
      ! 18 	crop          [moisture + deciduous]				        ! Croplands 
      ! 19 	crop          [moisture stress only]				        ! Croplands
      !
      ! ------------------------------------------------------------------------------------
      !
      ! PSO .csv PFT #s
      !
      ! 1 Forests
      ! 2 Shrublands
      ! 3 Savannas
      ! 4 Grasslands
      ! 5 Croplands
      ! 6 Bare
  
  
  
      integer, intent(in) :: clm_pft
      integer, intent(out) :: pso_pft
  
      if (clm_pft == 0) then
        pso_pft = 6
      elseif (clm_pft == 1) then
        pso_pft = 1
      elseif (clm_pft == 2) then
        pso_pft = 1
      elseif (clm_pft == 3) then
        pso_pft = 1
      elseif (clm_pft == 4) then
        pso_pft = 1
      elseif (clm_pft == 5) then
        pso_pft = 1
      elseif (clm_pft == 6) then
        pso_pft = 1
      elseif (clm_pft == 7) then
        pso_pft = 1
      elseif (clm_pft == 8) then
        pso_pft = 1
      elseif (clm_pft == 9) then
        pso_pft = 2
      elseif (clm_pft == 10) then
        pso_pft = 2
      elseif (clm_pft == 11) then
        pso_pft = 2
      elseif (clm_pft == 12) then
        pso_pft = 2
      elseif (clm_pft == 13) then
        pso_pft = 4
      elseif (clm_pft == 14) then
        pso_pft = 4
      elseif (clm_pft == 15) then
        pso_pft = 4
      elseif (clm_pft == 16) then
        pso_pft = 4
      elseif (clm_pft == 17) then
        pso_pft = 4
      elseif (clm_pft == 18) then
        pso_pft = 5
      elseif (clm_pft == 19) then
        pso_pft = 5
      endif
  
    end subroutine pft_clm_to_pso
  
end module compute_rc_functions 
