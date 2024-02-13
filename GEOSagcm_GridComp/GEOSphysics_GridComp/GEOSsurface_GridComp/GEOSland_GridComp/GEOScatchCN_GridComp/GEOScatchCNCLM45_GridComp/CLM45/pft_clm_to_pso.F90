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
      ! 1 	needleleaf evergreen temperate tree				            ! needleleaf tree
      ! 2 	needleleaf evergreen boreal tree				            ! needleleaf tree
      ! 3 	needleleaf deciduous boreal tree				            ! needleleaf tree
      ! 4 	broadleaf evergreen tropical tree				            ! broadleaf evergreen tree
      ! 5 	broadleaf evergreen temperate tree				            ! broadleaf evergreen tree
      ! 6 	broadleaf deciduous tropical tree				            ! broadleaf deciduous tree
      ! 7 	broadleaf deciduous temperate tree				            ! broadleaf deciduous
      ! 8 	broadleaf deciduous boreal tree					            ! broadleaf deciduous
      ! 9     broadleaf evergreen temperate shrub                       ! shrub
      ! 10 	broadleaf deciduous temperate shrub [moisture + deciduous]	! shrub
      ! 11 	broadleaf deciduous temperate shrub [moisture stress only]	! shrub
      ! 12 	broadleaf deciduous boreal shrub				            ! shrub
      ! 13 	arctic c3 grass							                    ! arctic c3 grass
      ! 14 	cool c3 grass [moisture + deciduous]				        ! c3 grass
      ! 15 	cool c3 grass [moisture stress only]   				        ! c3 grass
      ! 16 	warm c4 grass [moisture + deciduous]				        ! c4 grass
      ! 17 	warm c4 grass [moisture stress only]				        ! c4 grass
      ! 18 	crop          [moisture + deciduous]				        ! crop
      ! 19 	crop          [moisture stress only]				        ! crop
      !
      ! ------------------------------------------------------------------------------------
      !
      ! PSO PFT #s
      !
      ! 1 needleleaf tree
      ! 2 broadleaf evergreen tree
      ! 3 broadleaf deciduous tree
      ! 4 shrub
      ! 5 arctic c3 grass
      ! 6 c3 grass
      ! 7 c4 grass
      ! 8 crop
  
  
  
      integer, intent(in) :: clm_pft
      integer, intent(out) :: pso_pft
  
      if (clm_pft == 0) then
        pso_pft = 8
      elseif (clm_pft == 1) then
        pso_pft = 1
      elseif (clm_pft == 2) then
        pso_pft = 1
      elseif (clm_pft == 3) then
        pso_pft = 1
      elseif (clm_pft == 4) then
        pso_pft = 2
      elseif (clm_pft == 5) then
        pso_pft = 2
      elseif (clm_pft == 6) then
        pso_pft = 3
      elseif (clm_pft == 7) then
        pso_pft = 3
      elseif (clm_pft == 8) then
        pso_pft = 3
      elseif (clm_pft == 9) then
        pso_pft = 4
      elseif (clm_pft == 10) then
        pso_pft = 4
      elseif (clm_pft == 11) then
        pso_pft = 4
      elseif (clm_pft == 12) then
        pso_pft = 4
      elseif (clm_pft == 13) then
        pso_pft = 5
      elseif (clm_pft == 14) then
        pso_pft = 6
      elseif (clm_pft == 15) then
        pso_pft = 6
      elseif (clm_pft == 16) then
        pso_pft = 7
      elseif (clm_pft == 17) then
        pso_pft = 7
      elseif (clm_pft == 18) then
        pso_pft = 8
      elseif (clm_pft == 19) then
        pso_pft = 8
      endif
  
    end subroutine pft_clm_to_pso
  
end module compute_rc_functions 
