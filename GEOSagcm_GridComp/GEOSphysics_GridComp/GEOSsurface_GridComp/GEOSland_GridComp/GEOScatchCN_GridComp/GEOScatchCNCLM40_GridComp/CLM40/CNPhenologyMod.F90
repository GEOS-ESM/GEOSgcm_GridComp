module CNPhenologyMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: CNPhenologyMod
!
! !DESCRIPTION:
! Module holding routines used in phenology model for coupled carbon
! nitrogen code.
!
! !USES:
  use clmtype
  implicit none
  save
  private

! local variables to the whole module

! !PUBLIC MEMBER FUNCTIONS:
  public :: CNPhenology
!
! !REVISION HISTORY:
! 8/1/03: Created by Peter Thornton
! 10/23/03, Peter Thornton: migrated all routines to vector data structures
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNPhenology
!
! !INTERFACE:
subroutine CNPhenology (num_soilc, filter_soilc, num_soilp, filter_soilp)
!
! !DESCRIPTION:
! Dynamic phenology routine for coupled carbon-nitrogen code (CN)
! 1. grass phenology
!
! !USES:
!
! !ARGUMENTS:
   integer, intent(in) :: num_soilc       ! number of soil columns in filter
   integer, intent(in) :: filter_soilc(:) ! filter for soil columns
   integer, intent(in) :: num_soilp       ! number of soil pfts in filter
   integer, intent(in) :: filter_soilp(:) ! filter for soil pfts
!
! !CALLED FROM:
! subroutine CNEcosystemDyn in module CNEcosystemDynMod.F90
!
! !REVISION HISTORY:
! 7/28/03: Created by Peter Thornton
! 9/05/03, Peter Thornton: moved from call with (p) to call with (c)
! 10/3/03, Peter Thornton: added subroutine calls for different phenology types
! 11/7/03, Peter Thornton: moved phenology type tests into phenology type
!    routines, and moved onset, offset, background litfall routines into
!    main phenology call.
! !LOCAL VARIABLES:
! local pointers to implicit in arrays
!
! local pointers to implicit in/out scalars
!
! local pointers to implicit out scalars
!
! !OTHER LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------

   ! each of the following phenology type routines includes a filter
   ! to operate only on the relevant pfts

   call CNPhenologyClimate(num_soilp, filter_soilp)
   
   call CNEvergreenPhenology(num_soilp, filter_soilp)

   call CNSeasonDecidPhenology(num_soilp, filter_soilp)

   call CNStressDecidPhenology(num_soilp, filter_soilp)

   ! the same onset and offset routines are called regardless of
   ! phenology type - they depend only on onset_flag, offset_flag, bglfr, and bgtr

   call CNOnsetGrowth(num_soilp, filter_soilp)

   call CNOffsetLitterfall(num_soilp, filter_soilp)

   call CNBackgroundLitterfall(num_soilp, filter_soilp)

   call CNLivewoodTurnover(num_soilp, filter_soilp)

   ! gather all pft-level litterfall fluxes to the column
   ! for litter C and N inputs

   call CNLitterToColumn(num_soilc, filter_soilc)

end subroutine CNPhenology
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNPhenologyClimate
!
! !INTERFACE:
subroutine CNPhenologyClimate (num_soilp, filter_soilp)
!
! !DESCRIPTION:
! For coupled carbon-nitrogen code (CN).
!
! !USES:
   use clm_time_manager, only: get_step_size
!
! !ARGUMENTS:
   integer, intent(in) :: num_soilp       ! number of soil pfts in filter
   integer, intent(in) :: filter_soilp(:) ! filter for soil pfts
!
! !CALLED FROM:
! subroutine CNPhenology
!
! !REVISION HISTORY:
! 3/13/07: Created by Peter Thornton
!
! !LOCAL VARIABLES:
! local pointers to implicit in scalars
!
   integer , pointer :: ivt(:)       ! pft vegetation type
   ! ecophysiological constants
   real, pointer :: t_ref2m(:)            ! 2m air temperature (K)
   real, pointer :: tempavg_t2m(:)     ! temp. avg 2m air temperature (K)
!
! local pointers to implicit in/out scalars
!
!
! local pointers to implicit out scalars
!
! !OTHER LOCAL VARIABLES:
   integer :: p                      ! indices
   integer :: fp             !lake filter pft index
   real:: dt             !radiation time step delta t (seconds)
   real:: fracday        !dtime as a fraction of day
!EOP
!-----------------------------------------------------------------------

   ! assign local pointers to derived type arrays
   ivt       => clm3%g%l%c%p%itype
   t_ref2m                       => clm3%g%l%c%p%pes%t_ref2m
   tempavg_t2m                   => clm3%g%l%c%p%pepv%tempavg_t2m

   ! set time steps
   dt = real( get_step_size() )
   fracday = dt/86400.0

   do fp = 1,num_soilp
      p = filter_soilp(fp)
	  tempavg_t2m(p) = tempavg_t2m(p) + t_ref2m(p) * (fracday/365.)
   end do

end subroutine CNPhenologyClimate
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNEvergreenPhenology
!
! !INTERFACE:
subroutine CNEvergreenPhenology (num_soilp, filter_soilp)
!
! !DESCRIPTION:
! For coupled carbon-nitrogen code (CN).
!
! !USES:
!
! !ARGUMENTS:
   integer, intent(in) :: num_soilp       ! number of soil pfts in filter
   integer, intent(in) :: filter_soilp(:) ! filter for soil pfts
!
! !CALLED FROM:
! subroutine CNPhenology
!
! !REVISION HISTORY:
! 10/2/03: Created by Peter Thornton
!
! !LOCAL VARIABLES:
! local pointers to implicit in scalars
!
   integer , pointer :: ivt(:)       ! pft vegetation type
   ! ecophysiological constants
   real, pointer :: evergreen(:) ! binary flag for evergreen leaf habit (0 or 1)
   real, pointer :: leaf_long(:) ! leaf longevity (yrs)
!
! local pointers to implicit in/out scalars
!
   real, pointer :: bglfr(:)     ! background litterfall rate (1/s)
   real, pointer :: bgtr(:)      ! background transfer growth rate (1/s)
   real, pointer :: lgsf(:)      ! long growing season factor [0-1]
!
! local pointers to implicit out scalars
!
! !OTHER LOCAL VARIABLES:
   integer :: p                      ! indices
   integer :: fp                     ! lake filter pft index
!EOP
!-----------------------------------------------------------------------

   ! assign local pointers to derived type arrays
   ivt       => clm3%g%l%c%p%itype
   evergreen => pftcon%evergreen
   leaf_long => pftcon%leaf_long
   bglfr     => clm3%g%l%c%p%pepv%bglfr
   bgtr      => clm3%g%l%c%p%pepv%bgtr
   lgsf      => clm3%g%l%c%p%pepv%lgsf

   do fp = 1,num_soilp
      p = filter_soilp(fp)
      if (evergreen(ivt(p)) == 1.) then
          bglfr(p) = 1./(leaf_long(ivt(p))*365.*86400.)
          bgtr(p)  = 0.
          lgsf(p)  = 0.
      end if
   end do

end subroutine CNEvergreenPhenology
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNSeasonDecidPhenology
!
! !INTERFACE:
subroutine CNSeasonDecidPhenology (num_soilp, filter_soilp)
!
! !DESCRIPTION:
! For coupled carbon-nitrogen code (CN).
! This routine handles the seasonal deciduous phenology code (temperate
! deciduous vegetation that has only one growing season per year).
!
! !USES:
   use clm_time_manager, only: get_step_size
   use shr_const_mod, only: SHR_CONST_TKFRZ
!
! !ARGUMENTS:
   integer, intent(in) :: num_soilp       ! number of soil pfts in filter
   integer, intent(in) :: filter_soilp(:) ! filter for soil pfts
!
! !CALLED FROM:
! subroutine CNPhenology
!
! !REVISION HISTORY:
! 10/6/03: Created by Peter Thornton
! 10/24/03, Peter Thornton: migrated to vector data structures
!
! !LOCAL VARIABLES:
! local pointers to implicit in scalars
   integer , pointer :: ivt(:)                ! pft vegetation type
   integer , pointer :: pcolumn(:)            ! pft's column index
   integer , pointer :: pgridcell(:)          ! pft's gridcell index
   real, pointer :: t_soisno(:,:)         ! soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)
   real, pointer :: soilpsi(:,:)          ! soil water potential in each soil layer (MPa)
   real, pointer :: leafc_storage(:)      ! (kgC/m2) leaf C storage
   real, pointer :: frootc_storage(:)     ! (kgC/m2) fine root C storage
   real, pointer :: livestemc_storage(:)  ! (kgC/m2) live stem C storage
   real, pointer :: deadstemc_storage(:)  ! (kgC/m2) dead stem C storage
   real, pointer :: livecrootc_storage(:) ! (kgC/m2) live coarse root C storage
   real, pointer :: deadcrootc_storage(:) ! (kgC/m2) dead coarse root C storage
   real, pointer :: gresp_storage(:)      ! (kgC/m2) growth respiration storage
   real, pointer :: leafn_storage(:)      ! (kgN/m2) leaf N storage
   real, pointer :: frootn_storage(:)     ! (kgN/m2) fine root N storage
   real, pointer :: livestemn_storage(:)  ! (kgN/m2) live stem N storage
   real, pointer :: deadstemn_storage(:)  ! (kgN/m2) dead stem N storage
   real, pointer :: livecrootn_storage(:) ! (kgN/m2) live coarse root N storage
   real, pointer :: deadcrootn_storage(:) ! (kgN/m2) dead coarse root N storage
   real, pointer :: t_grnd(:)             ! ground temperature (Kelvin)
   ! ecophysiological constants
   real, pointer :: season_decid(:) ! binary flag for seasonal-deciduous leaf habit (0 or 1)
   real, pointer :: woody(:)        ! binary flag for woody lifeform (1=woody, 0=not woody)
!
! local pointers to implicit in/out scalars
   real, pointer :: dormant_flag(:)    ! dormancy flag
   real, pointer :: days_active(:)     ! number of days since last dormancy
   real, pointer :: onset_flag(:)      ! onset flag
   real, pointer :: onset_counter(:)   ! onset counter (seconds)
   real, pointer :: onset_gddflag(:)   ! onset freeze flag
   real, pointer :: onset_gdd(:)       ! onset growing degree days
   real, pointer :: offset_flag(:)     ! offset flag
   real, pointer :: offset_counter(:)  ! offset counter (seconds)
   real, pointer :: dayl(:)            ! daylength (seconds)
   real, pointer :: prev_dayl(:)       ! daylength from previous albedo timestep (seconds)
   real, pointer :: annavg_t2m(:)      ! annual average 2m air temperature (K)
   real, pointer :: prev_leafc_to_litter(:)  ! previous timestep leaf C litterfall flux (gC/m2/s)
   real, pointer :: prev_frootc_to_litter(:) ! previous timestep froot C litterfall flux (gC/m2/s)
   real, pointer :: lgsf(:)            ! long growing season factor [0-1]
   real, pointer :: bglfr(:)           ! background litterfall rate (1/s)
   real, pointer :: bgtr(:)            ! background transfer growth rate (1/s)
   real, pointer :: leafc_xfer_to_leafc(:)
   real, pointer :: frootc_xfer_to_frootc(:)
   real, pointer :: livestemc_xfer_to_livestemc(:)
   real, pointer :: deadstemc_xfer_to_deadstemc(:)
   real, pointer :: livecrootc_xfer_to_livecrootc(:)
   real, pointer :: deadcrootc_xfer_to_deadcrootc(:)
   real, pointer :: leafn_xfer_to_leafn(:)
   real, pointer :: frootn_xfer_to_frootn(:)
   real, pointer :: livestemn_xfer_to_livestemn(:)
   real, pointer :: deadstemn_xfer_to_deadstemn(:)
   real, pointer :: livecrootn_xfer_to_livecrootn(:)
   real, pointer :: deadcrootn_xfer_to_deadcrootn(:)
   real, pointer :: leafc_xfer(:)      ! (kgC/m2) leaf C transfer
   real, pointer :: frootc_xfer(:)     ! (kgC/m2) fine root C transfer
   real, pointer :: livestemc_xfer(:)  ! (kgC/m2) live stem C transfer
   real, pointer :: deadstemc_xfer(:)  ! (kgC/m2) dead stem C transfer
   real, pointer :: livecrootc_xfer(:) ! (kgC/m2) live coarse root C transfer
   real, pointer :: deadcrootc_xfer(:) ! (kgC/m2) dead coarse root C transfer
   real, pointer :: leafn_xfer(:)      ! (kgN/m2) leaf N transfer
   real, pointer :: frootn_xfer(:)     ! (kgN/m2) fine root N transfer
   real, pointer :: livestemn_xfer(:)  ! (kgN/m2) live stem N transfer
   real, pointer :: deadstemn_xfer(:)  ! (kgN/m2) dead stem N transfer
   real, pointer :: livecrootn_xfer(:) ! (kgN/m2) live coarse root N transfer
   real, pointer :: deadcrootn_xfer(:) ! (kgN/m2) dead coarse root N transfer
   real, pointer :: leafc_storage_to_xfer(:)
   real, pointer :: frootc_storage_to_xfer(:)
   real, pointer :: livestemc_storage_to_xfer(:)
   real, pointer :: deadstemc_storage_to_xfer(:)
   real, pointer :: livecrootc_storage_to_xfer(:)
   real, pointer :: deadcrootc_storage_to_xfer(:)
   real, pointer :: gresp_storage_to_xfer(:)
   real, pointer :: leafn_storage_to_xfer(:)
   real, pointer :: frootn_storage_to_xfer(:)
   real, pointer :: livestemn_storage_to_xfer(:)
   real, pointer :: deadstemn_storage_to_xfer(:)
   real, pointer :: livecrootn_storage_to_xfer(:)
   real, pointer :: deadcrootn_storage_to_xfer(:)
!
! local pointers to implicit out scalars
!
! !OTHER LOCAL VARIABLES:
   integer :: c,p            !indices
   integer :: fp             !lake filter pft index
   real:: dt             !radiation time step delta t (seconds)
   real:: fracday        !dtime as a fraction of day
   real:: crit_dayl      !critical daylength for offset (seconds)
   real:: ws_flag        !winter-summer solstice flag (0 or 1)
   real:: crit_onset_gdd !critical onset growing degree-day sum
   real:: ndays_on       !number of days to complete onset
   real:: ndays_off      !number of days to complete offset
   real:: soilt
   real:: fstor2tran     !fraction of storage to move to transfer on each onset

!EOP
!-----------------------------------------------------------------------
   ! Assign local pointers to derived type arrays (in)
   ivt                           => clm3%g%l%c%p%itype
   pcolumn                       => clm3%g%l%c%p%column
   pgridcell                     => clm3%g%l%c%p%gridcell
   t_soisno                      => clm3%g%l%c%ces%t_soisno
   t_grnd                        => clm3%g%l%c%ces%t_grnd
   leafc_storage                 => clm3%g%l%c%p%pcs%leafc_storage
   frootc_storage                => clm3%g%l%c%p%pcs%frootc_storage
   livestemc_storage             => clm3%g%l%c%p%pcs%livestemc_storage
   deadstemc_storage             => clm3%g%l%c%p%pcs%deadstemc_storage
   livecrootc_storage            => clm3%g%l%c%p%pcs%livecrootc_storage
   deadcrootc_storage            => clm3%g%l%c%p%pcs%deadcrootc_storage
   gresp_storage                 => clm3%g%l%c%p%pcs%gresp_storage
   leafn_storage                 => clm3%g%l%c%p%pns%leafn_storage
   frootn_storage                => clm3%g%l%c%p%pns%frootn_storage
   livestemn_storage             => clm3%g%l%c%p%pns%livestemn_storage
   deadstemn_storage             => clm3%g%l%c%p%pns%deadstemn_storage
   livecrootn_storage            => clm3%g%l%c%p%pns%livecrootn_storage
   deadcrootn_storage            => clm3%g%l%c%p%pns%deadcrootn_storage
   season_decid                  => pftcon%season_decid
   woody                         => pftcon%woody

   ! Assign local pointers to derived type arrays (out)
   dormant_flag                  => clm3%g%l%c%p%pepv%dormant_flag
   days_active                   => clm3%g%l%c%p%pepv%days_active
   onset_flag                    => clm3%g%l%c%p%pepv%onset_flag
   onset_counter                 => clm3%g%l%c%p%pepv%onset_counter
   onset_gddflag                 => clm3%g%l%c%p%pepv%onset_gddflag
   onset_gdd                     => clm3%g%l%c%p%pepv%onset_gdd
   offset_flag                   => clm3%g%l%c%p%pepv%offset_flag
   offset_counter                => clm3%g%l%c%p%pepv%offset_counter
   dayl                          => clm3%g%l%c%p%pepv%dayl
   prev_dayl                     => clm3%g%l%c%p%pepv%prev_dayl
   annavg_t2m                    => clm3%g%l%c%p%pepv%annavg_t2m
   prev_leafc_to_litter          => clm3%g%l%c%p%pepv%prev_leafc_to_litter
   prev_frootc_to_litter         => clm3%g%l%c%p%pepv%prev_frootc_to_litter
   bglfr                         => clm3%g%l%c%p%pepv%bglfr
   bgtr                          => clm3%g%l%c%p%pepv%bgtr
   lgsf                          => clm3%g%l%c%p%pepv%lgsf
   leafc_xfer_to_leafc           => clm3%g%l%c%p%pcf%leafc_xfer_to_leafc
   frootc_xfer_to_frootc         => clm3%g%l%c%p%pcf%frootc_xfer_to_frootc
   livestemc_xfer_to_livestemc   => clm3%g%l%c%p%pcf%livestemc_xfer_to_livestemc
   deadstemc_xfer_to_deadstemc   => clm3%g%l%c%p%pcf%deadstemc_xfer_to_deadstemc
   livecrootc_xfer_to_livecrootc => clm3%g%l%c%p%pcf%livecrootc_xfer_to_livecrootc
   deadcrootc_xfer_to_deadcrootc => clm3%g%l%c%p%pcf%deadcrootc_xfer_to_deadcrootc
   leafn_xfer_to_leafn           => clm3%g%l%c%p%pnf%leafn_xfer_to_leafn
   frootn_xfer_to_frootn         => clm3%g%l%c%p%pnf%frootn_xfer_to_frootn
   livestemn_xfer_to_livestemn   => clm3%g%l%c%p%pnf%livestemn_xfer_to_livestemn
   deadstemn_xfer_to_deadstemn   => clm3%g%l%c%p%pnf%deadstemn_xfer_to_deadstemn
   livecrootn_xfer_to_livecrootn => clm3%g%l%c%p%pnf%livecrootn_xfer_to_livecrootn
   deadcrootn_xfer_to_deadcrootn => clm3%g%l%c%p%pnf%deadcrootn_xfer_to_deadcrootn
   leafc_xfer                    => clm3%g%l%c%p%pcs%leafc_xfer
   frootc_xfer                   => clm3%g%l%c%p%pcs%frootc_xfer
   livestemc_xfer                => clm3%g%l%c%p%pcs%livestemc_xfer
   deadstemc_xfer                => clm3%g%l%c%p%pcs%deadstemc_xfer
   livecrootc_xfer               => clm3%g%l%c%p%pcs%livecrootc_xfer
   deadcrootc_xfer               => clm3%g%l%c%p%pcs%deadcrootc_xfer
   leafn_xfer                    => clm3%g%l%c%p%pns%leafn_xfer
   frootn_xfer                   => clm3%g%l%c%p%pns%frootn_xfer
   livestemn_xfer                => clm3%g%l%c%p%pns%livestemn_xfer
   deadstemn_xfer                => clm3%g%l%c%p%pns%deadstemn_xfer
   livecrootn_xfer               => clm3%g%l%c%p%pns%livecrootn_xfer
   deadcrootn_xfer               => clm3%g%l%c%p%pns%deadcrootn_xfer
   leafc_storage_to_xfer         => clm3%g%l%c%p%pcf%leafc_storage_to_xfer
   frootc_storage_to_xfer        => clm3%g%l%c%p%pcf%frootc_storage_to_xfer
   livestemc_storage_to_xfer     => clm3%g%l%c%p%pcf%livestemc_storage_to_xfer
   deadstemc_storage_to_xfer     => clm3%g%l%c%p%pcf%deadstemc_storage_to_xfer
   livecrootc_storage_to_xfer    => clm3%g%l%c%p%pcf%livecrootc_storage_to_xfer
   deadcrootc_storage_to_xfer    => clm3%g%l%c%p%pcf%deadcrootc_storage_to_xfer
   gresp_storage_to_xfer         => clm3%g%l%c%p%pcf%gresp_storage_to_xfer
   leafn_storage_to_xfer         => clm3%g%l%c%p%pnf%leafn_storage_to_xfer
   frootn_storage_to_xfer        => clm3%g%l%c%p%pnf%frootn_storage_to_xfer
   livestemn_storage_to_xfer     => clm3%g%l%c%p%pnf%livestemn_storage_to_xfer
   deadstemn_storage_to_xfer     => clm3%g%l%c%p%pnf%deadstemn_storage_to_xfer
   livecrootn_storage_to_xfer    => clm3%g%l%c%p%pnf%livecrootn_storage_to_xfer
   deadcrootn_storage_to_xfer    => clm3%g%l%c%p%pnf%deadcrootn_storage_to_xfer

   ! set time steps
   dt = real( get_step_size() )
   fracday = dt/86400.0

   ! critical daylength from Biome-BGC, v4.1.2
   crit_dayl = 39300.
   ndays_on = 30.
   ndays_off = 15.

   ! transfer parameters
   fstor2tran = 0.5

   ! start pft loop
   do fp = 1,num_soilp
      p = filter_soilp(fp)
      c = pcolumn(p)

      if (season_decid(ivt(p)) == 1.) then

         ! set background litterfall rate, background transfer rate, and
         ! long growing season factor to 0 for seasonal deciduous types
         bglfr(p) = 0.
         bgtr(p) = 0.
         lgsf(p) = 0.

         ! onset gdd sum from Biome-BGC, v4.1.2
         crit_onset_gdd = exp(4.8 + 0.13*(annavg_t2m(p) - SHR_CONST_TKFRZ))

         ! set flag for solstice period (winter->summer = 1, summer->winter = 0)
         if (dayl(p) >= prev_dayl(p)) then
            ws_flag = 1.
         else
            ws_flag = 0.
         end if

         ! update offset_counter and test for the end of the offset period
         if (offset_flag(p) == 1.0) then
            ! decrement counter for offset period
            offset_counter(p) = offset_counter(p) - dt

            ! if this is the end of the offset_period, reset phenology
            ! flags and indices
            if (offset_counter(p) <= 0.0) then
               ! this code block was originally handled by call cn_offset_cleanup(p)
               ! inlined during vectorization

               offset_flag(p) = 0.
               offset_counter(p) = 0.
               dormant_flag(p) = 1.
               days_active(p) = 0.

               ! reset the previous timestep litterfall flux memory
               prev_leafc_to_litter(p) = 0.
               prev_frootc_to_litter(p) = 0.
            end if
         end if

         ! update onset_counter and test for the end of the onset period
         if (onset_flag(p) == 1.0) then
            ! decrement counter for onset period
            onset_counter(p) = onset_counter(p) - dt

            ! if this is the end of the onset period, reset phenology
            ! flags and indices
            if (onset_counter(p) <= 0.0) then
               ! this code block was originally handled by call cn_onset_cleanup(p)
               ! inlined during vectorization

               onset_flag(p) = 0.0
               onset_counter(p) = 0.0
               ! set all transfer growth rates to 0.0
               leafc_xfer_to_leafc(p)   = 0.0
               frootc_xfer_to_frootc(p) = 0.0
               leafn_xfer_to_leafn(p)   = 0.0
               frootn_xfer_to_frootn(p) = 0.0
               if (woody(ivt(p)) == 1.0) then
                  livestemc_xfer_to_livestemc(p)   = 0.0
                  deadstemc_xfer_to_deadstemc(p)   = 0.0
                  livecrootc_xfer_to_livecrootc(p) = 0.0
                  deadcrootc_xfer_to_deadcrootc(p) = 0.0
                  livestemn_xfer_to_livestemn(p)   = 0.0
                  deadstemn_xfer_to_deadstemn(p)   = 0.0
                  livecrootn_xfer_to_livecrootn(p) = 0.0
                  deadcrootn_xfer_to_deadcrootn(p) = 0.0
               end if
               ! set transfer pools to 0.0
               leafc_xfer(p) = 0.0
               leafn_xfer(p) = 0.0
               frootc_xfer(p) = 0.0
               frootn_xfer(p) = 0.0
               if (woody(ivt(p)) == 1.0) then
                  livestemc_xfer(p) = 0.0
                  livestemn_xfer(p) = 0.0
                  deadstemc_xfer(p) = 0.0
                  deadstemn_xfer(p) = 0.0
                  livecrootc_xfer(p) = 0.0
                  livecrootn_xfer(p) = 0.0
                  deadcrootc_xfer(p) = 0.0
                  deadcrootn_xfer(p) = 0.0
               end if
            end if
         end if

         ! test for switching from dormant period to growth period
         if (dormant_flag(p) == 1.0) then

            ! Test to turn on growing degree-day sum, if off.
            ! switch on the growing degree day sum on the winter solstice

            if (onset_gddflag(p) == 0. .and. ws_flag == 1.) then
               onset_gddflag(p) = 1.
               onset_gdd(p) = 0.
            end if

            ! Test to turn off growing degree-day sum, if on.
            ! This test resets the growing degree day sum if it gets past
            ! the summer solstice without reaching the threshold value.
            ! In that case, it will take until the next winter solstice
            ! before the growing degree-day summation starts again.

            if (onset_gddflag(p) == 1. .and. ws_flag == 0.) then
               onset_gddflag(p) = 0.
               onset_gdd(p) = 0.
            end if

            ! if the gdd flag is set, and if the soil is above freezing
            ! then accumulate growing degree days for onset trigger

!!!!        soilt = t_soisno(c,2) ! gkw: hardwared code!
            soilt = max(t_grnd(c),t_soisno(c,1))
            if (onset_gddflag(p) == 1.0 .and. soilt > SHR_CONST_TKFRZ) then
               onset_gdd(p) = onset_gdd(p) + (soilt-SHR_CONST_TKFRZ)*fracday
            end if

            ! set onset_flag if critical growing degree-day sum is exceeded
            if (onset_gdd(p) > crit_onset_gdd) then
               onset_flag(p) = 1.0
               dormant_flag(p) = 0.0
               onset_gddflag(p) = 0.0
               onset_gdd(p) = 0.0
               onset_counter(p) = ndays_on * 86400.0

               ! move all the storage pools into transfer pools,
               ! where they will be transfered to displayed growth over the onset period.
               ! this code was originally handled with call cn_storage_to_xfer(p)
               ! inlined during vectorization

               ! set carbon fluxes for shifting storage pools to transfer pools
               leafc_storage_to_xfer(p)  = fstor2tran * leafc_storage(p)/dt
               frootc_storage_to_xfer(p) = fstor2tran * frootc_storage(p)/dt
               if (woody(ivt(p)) == 1.0) then
                  livestemc_storage_to_xfer(p)  = fstor2tran * livestemc_storage(p)/dt
                  deadstemc_storage_to_xfer(p)  = fstor2tran * deadstemc_storage(p)/dt
                  livecrootc_storage_to_xfer(p) = fstor2tran * livecrootc_storage(p)/dt
                  deadcrootc_storage_to_xfer(p) = fstor2tran * deadcrootc_storage(p)/dt
                  gresp_storage_to_xfer(p)      = fstor2tran * gresp_storage(p)/dt
               end if

               ! set nitrogen fluxes for shifting storage pools to transfer pools
               leafn_storage_to_xfer(p)  = fstor2tran * leafn_storage(p)/dt
               frootn_storage_to_xfer(p) = fstor2tran * frootn_storage(p)/dt
               if (woody(ivt(p)) == 1.0) then
                  livestemn_storage_to_xfer(p)  = fstor2tran * livestemn_storage(p)/dt
                  deadstemn_storage_to_xfer(p)  = fstor2tran * deadstemn_storage(p)/dt
                  livecrootn_storage_to_xfer(p) = fstor2tran * livecrootn_storage(p)/dt
                  deadcrootn_storage_to_xfer(p) = fstor2tran * deadcrootn_storage(p)/dt
               end if
            end if

         ! test for switching from growth period to offset period
         else if (offset_flag(p) == 0.0) then

            ! only begin to test for offset daylength once past the summer sol
            if (ws_flag == 0. .and. dayl(p) < crit_dayl) then
               offset_flag(p) = 1.
               offset_counter(p) = ndays_off * 86400.0
               prev_leafc_to_litter(p) = 0.
               prev_frootc_to_litter(p) = 0.
            end if
         end if

      end if ! end if seasonal deciduous

   end do ! end of pft loop

end subroutine CNSeasonDecidPhenology
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNStressDecidPhenology
!
! !INTERFACE:
subroutine CNStressDecidPhenology (num_soilp, filter_soilp)
!
! !DESCRIPTION:
! This routine handles phenology for vegetation types, such as grasses and
! tropical drought deciduous trees, that respond to cold and drought stress
! signals and that can have multiple growing seasons in a given year.
! This routine allows for the possibility that leaves might persist year-round
! in the absence of a suitable stress trigger, by switching to an essentially
! evergreen habit, but maintaining a deciduous leaf longevity, while waiting
! for the next stress trigger.  This is in contrast to the seasonal deciduous
! algorithm (for temperate deciduous trees) that forces a single growing season
! per year.
!
! !USES:
   use clm_time_manager, only: get_step_size
   use shr_const_mod, only: SHR_CONST_TKFRZ
!
! !ARGUMENTS:
   integer, intent(in) :: num_soilp       ! number of soil pfts in filter
   integer, intent(in) :: filter_soilp(:) ! filter for soil pfts
!
! !CALLED FROM:
! subroutine CNPhenology
!
! !REVISION HISTORY:
! 10/27/03: Created by Peter Thornton
! 01/29/04: Made onset_gdd critical sum a function of temperature, as in
!           seasonal deciduous algorithm.
!
! !LOCAL VARIABLES:
! local pointers to implicit in scalars
!
   integer , pointer :: ivt(:)                ! pft vegetation type
   integer , pointer :: pcolumn(:)            ! pft's column index
   integer , pointer :: pgridcell(:)          ! pft's gridcell index
   real, pointer :: leafc_storage(:)      ! (kgC/m2) leaf C storage
   real, pointer :: frootc_storage(:)     ! (kgC/m2) fine root C storage
   real, pointer :: livestemc_storage(:)  ! (kgC/m2) live stem C storage
   real, pointer :: deadstemc_storage(:)  ! (kgC/m2) dead stem C storage
   real, pointer :: livecrootc_storage(:) ! (kgC/m2) live coarse root C storage
   real, pointer :: deadcrootc_storage(:) ! (kgC/m2) dead coarse root C storage
   real, pointer :: gresp_storage(:)      ! (kgC/m2) growth respiration storage
   real, pointer :: leafn_storage(:)      ! (kgN/m2) leaf N storage
   real, pointer :: frootn_storage(:)     ! (kgN/m2) fine root N storage
   real, pointer :: livestemn_storage(:)  ! (kgN/m2) live stem N storage
   real, pointer :: deadstemn_storage(:)  ! (kgN/m2) dead stem N storage
   real, pointer :: livecrootn_storage(:) ! (kgN/m2) live coarse root N storage
   real, pointer :: deadcrootn_storage(:) ! (kgN/m2) dead coarse root N storage
   real, pointer :: t_soisno(:,:)         ! soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)
   real, pointer :: soilpsi(:,:)          ! soil water potential in each soil layer (MPa)
   real, pointer :: psiwilt(:)            ! root-zone soil water potential at wilting point (MPa)
   real, pointer :: leaf_long(:)          ! leaf longevity (yrs)
   real, pointer :: stress_decid(:)       ! binary flag for stress-deciduous leaf habit (0 or 1)
   real, pointer :: woody(:)              ! binary flag for woody lifeform (1=woody, 0=not woody)
   real, pointer :: t_grnd(:)             ! ground temperature (Kelvin)
!
! local pointers to implicit in/out scalars
!
   real, pointer :: dormant_flag(:)    ! dormancy flag
   real, pointer :: days_active(:)     ! number of days since last dormancy
   real, pointer :: onset_flag(:)      ! onset flag
   real, pointer :: onset_counter(:)   ! onset counter (seconds)
   real, pointer :: onset_gddflag(:)   ! onset freeze flag
   real, pointer :: onset_fdd(:)       ! onset freezing degree days counter
   real, pointer :: onset_gdd(:)       ! onset growing degree days
   real, pointer :: onset_swi(:)       ! onset soil water index
   real, pointer :: offset_flag(:)     ! offset flag
   real, pointer :: offset_counter(:)  ! offset counter (seconds)
   real, pointer :: prev_dayl(:)       ! daylength from previous albedo timestep (seconds)
   real, pointer :: dayl(:)            ! daylength (seconds)
   real, pointer :: offset_fdd(:)      ! offset freezing degree days counter
   real, pointer :: offset_swi(:)      ! offset soil water index
   real, pointer :: annavg_t2m(:)      ! annual average 2m air temperature (K)
   real, pointer :: lgsf(:)            ! long growing season factor [0-1]
   real, pointer :: bglfr(:)           ! background litterfall rate (1/s)
   real, pointer :: bgtr(:)            ! background transfer growth rate (1/s)
   real, pointer :: prev_leafc_to_litter(:)  ! previous timestep leaf C litterfall flux (gC/m2/s)
   real, pointer :: prev_frootc_to_litter(:) ! previous timestep froot C litterfall flux (gC/m2/s)
   real, pointer :: leafc_xfer_to_leafc(:)
   real, pointer :: frootc_xfer_to_frootc(:)
   real, pointer :: livestemc_xfer_to_livestemc(:)
   real, pointer :: deadstemc_xfer_to_deadstemc(:)
   real, pointer :: livecrootc_xfer_to_livecrootc(:)
   real, pointer :: deadcrootc_xfer_to_deadcrootc(:)
   real, pointer :: leafn_xfer_to_leafn(:)
   real, pointer :: frootn_xfer_to_frootn(:)
   real, pointer :: livestemn_xfer_to_livestemn(:)
   real, pointer :: deadstemn_xfer_to_deadstemn(:)
   real, pointer :: livecrootn_xfer_to_livecrootn(:)
   real, pointer :: deadcrootn_xfer_to_deadcrootn(:)
   real, pointer :: leafc_xfer(:)      ! (kgC/m2) leaf C transfer
   real, pointer :: frootc_xfer(:)     ! (kgC/m2) fine root C transfer
   real, pointer :: livestemc_xfer(:)  ! (kgC/m2) live stem C transfer
   real, pointer :: deadstemc_xfer(:)  ! (kgC/m2) dead stem C transfer
   real, pointer :: livecrootc_xfer(:) ! (kgC/m2) live coarse root C transfer
   real, pointer :: deadcrootc_xfer(:) ! (kgC/m2) dead coarse root C transfer
   real, pointer :: leafn_xfer(:)      ! (kgN/m2) leaf N transfer
   real, pointer :: frootn_xfer(:)     ! (kgN/m2) fine root N transfer
   real, pointer :: livestemn_xfer(:)  ! (kgN/m2) live stem N transfer
   real, pointer :: deadstemn_xfer(:)  ! (kgN/m2) dead stem N transfer
   real, pointer :: livecrootn_xfer(:) ! (kgN/m2) live coarse root N transfer
   real, pointer :: deadcrootn_xfer(:) ! (kgN/m2) dead coarse root N transfer
   real, pointer :: leafc_storage_to_xfer(:)
   real, pointer :: frootc_storage_to_xfer(:)
   real, pointer :: livestemc_storage_to_xfer(:)
   real, pointer :: deadstemc_storage_to_xfer(:)
   real, pointer :: livecrootc_storage_to_xfer(:)
   real, pointer :: deadcrootc_storage_to_xfer(:)
   real, pointer :: gresp_storage_to_xfer(:)
   real, pointer :: leafn_storage_to_xfer(:)
   real, pointer :: frootn_storage_to_xfer(:)
   real, pointer :: livestemn_storage_to_xfer(:)
   real, pointer :: deadstemn_storage_to_xfer(:)
   real, pointer :: livecrootn_storage_to_xfer(:)
   real, pointer :: deadcrootn_storage_to_xfer(:)
!
! local pointers to implicit out scalars
!
!
! !OTHER LOCAL VARIABLES:
   integer :: c,p             ! indices
   integer :: fp              ! lake filter pft index
   real:: fracday         ! dtime as a fraction of day
   real:: crit_dayl       ! critical daylength for offset (seconds)
   real:: ws_flag         ! winter-summer solstice flag (0 or 1)
   real:: dt              ! radiation time step delta t (seconds)
   real:: crit_onset_fdd  ! critical number of freezing days
   real:: crit_onset_gdd  ! degree days for onset trigger
   real:: crit_offset_fdd ! critical number of freezing degree days
                              ! to trigger offset
   real:: crit_onset_swi  ! water stress days for offset trigger
   real:: crit_offset_swi ! water stress days for offset trigger
   real:: soilpsi_on      ! water potential for onset trigger (MPa)
   real:: soilpsi_off     ! water potential for offset trigger (MPa)
   real:: ndays_on        ! number of days to complete onset
   real:: ndays_off       ! number of days to complete offset
   real:: soilt           ! temperature of top soil layer
   real:: psi             ! water stress of top soil layer
   real:: fstor2tran      ! fraction of storage to move to transfer
                              ! on each onset
!EOP
!-----------------------------------------------------------------------
   ! Assign local pointers to derived type arrays (in)
    ivt                            => clm3%g%l%c%p%itype
    pcolumn                        => clm3%g%l%c%p%column
    pgridcell                      => clm3%g%l%c%p%gridcell
    leafc_storage                  => clm3%g%l%c%p%pcs%leafc_storage
    frootc_storage                 => clm3%g%l%c%p%pcs%frootc_storage
    livestemc_storage              => clm3%g%l%c%p%pcs%livestemc_storage
    deadstemc_storage              => clm3%g%l%c%p%pcs%deadstemc_storage
    livecrootc_storage             => clm3%g%l%c%p%pcs%livecrootc_storage
    deadcrootc_storage             => clm3%g%l%c%p%pcs%deadcrootc_storage
    gresp_storage                  => clm3%g%l%c%p%pcs%gresp_storage
    leafn_storage                  => clm3%g%l%c%p%pns%leafn_storage
    frootn_storage                 => clm3%g%l%c%p%pns%frootn_storage
    livestemn_storage              => clm3%g%l%c%p%pns%livestemn_storage
    deadstemn_storage              => clm3%g%l%c%p%pns%deadstemn_storage
    livecrootn_storage             => clm3%g%l%c%p%pns%livecrootn_storage
    deadcrootn_storage             => clm3%g%l%c%p%pns%deadcrootn_storage
    soilpsi                        => clm3%g%l%c%cps%soilpsi
    psiwilt                        => clm3%g%l%c%cps%psiwilt
    t_soisno                       => clm3%g%l%c%ces%t_soisno
    t_grnd                         => clm3%g%l%c%ces%t_grnd
    leaf_long                      => pftcon%leaf_long
    woody                          => pftcon%woody
    stress_decid                   => pftcon%stress_decid

   ! Assign local pointers to derived type arrays (out)
    dormant_flag                   => clm3%g%l%c%p%pepv%dormant_flag
    days_active                    => clm3%g%l%c%p%pepv%days_active
    onset_flag                     => clm3%g%l%c%p%pepv%onset_flag
    onset_counter                  => clm3%g%l%c%p%pepv%onset_counter
    onset_gddflag                  => clm3%g%l%c%p%pepv%onset_gddflag
    onset_fdd                      => clm3%g%l%c%p%pepv%onset_fdd
    onset_gdd                      => clm3%g%l%c%p%pepv%onset_gdd
    onset_swi                      => clm3%g%l%c%p%pepv%onset_swi
    offset_flag                    => clm3%g%l%c%p%pepv%offset_flag
    offset_counter                 => clm3%g%l%c%p%pepv%offset_counter
    dayl                           => clm3%g%l%c%p%pepv%dayl
    prev_dayl                      => clm3%g%l%c%p%pepv%prev_dayl
    offset_fdd                     => clm3%g%l%c%p%pepv%offset_fdd
    offset_swi                     => clm3%g%l%c%p%pepv%offset_swi
    annavg_t2m                     => clm3%g%l%c%p%pepv%annavg_t2m
    prev_leafc_to_litter           => clm3%g%l%c%p%pepv%prev_leafc_to_litter
    prev_frootc_to_litter          => clm3%g%l%c%p%pepv%prev_frootc_to_litter
    lgsf                           => clm3%g%l%c%p%pepv%lgsf
    bglfr                          => clm3%g%l%c%p%pepv%bglfr
    bgtr                           => clm3%g%l%c%p%pepv%bgtr
    leafc_xfer_to_leafc            => clm3%g%l%c%p%pcf%leafc_xfer_to_leafc
    frootc_xfer_to_frootc          => clm3%g%l%c%p%pcf%frootc_xfer_to_frootc
    livestemc_xfer_to_livestemc    => clm3%g%l%c%p%pcf%livestemc_xfer_to_livestemc
    deadstemc_xfer_to_deadstemc    => clm3%g%l%c%p%pcf%deadstemc_xfer_to_deadstemc
    livecrootc_xfer_to_livecrootc  => clm3%g%l%c%p%pcf%livecrootc_xfer_to_livecrootc
    deadcrootc_xfer_to_deadcrootc  => clm3%g%l%c%p%pcf%deadcrootc_xfer_to_deadcrootc
    leafn_xfer_to_leafn            => clm3%g%l%c%p%pnf%leafn_xfer_to_leafn
    frootn_xfer_to_frootn          => clm3%g%l%c%p%pnf%frootn_xfer_to_frootn
    livestemn_xfer_to_livestemn    => clm3%g%l%c%p%pnf%livestemn_xfer_to_livestemn
    deadstemn_xfer_to_deadstemn    => clm3%g%l%c%p%pnf%deadstemn_xfer_to_deadstemn
    livecrootn_xfer_to_livecrootn  => clm3%g%l%c%p%pnf%livecrootn_xfer_to_livecrootn
    deadcrootn_xfer_to_deadcrootn  => clm3%g%l%c%p%pnf%deadcrootn_xfer_to_deadcrootn
    leafc_xfer                     => clm3%g%l%c%p%pcs%leafc_xfer
    frootc_xfer                    => clm3%g%l%c%p%pcs%frootc_xfer
    livestemc_xfer                 => clm3%g%l%c%p%pcs%livestemc_xfer
    deadstemc_xfer                 => clm3%g%l%c%p%pcs%deadstemc_xfer
    livecrootc_xfer                => clm3%g%l%c%p%pcs%livecrootc_xfer
    deadcrootc_xfer                => clm3%g%l%c%p%pcs%deadcrootc_xfer
    leafn_xfer                     => clm3%g%l%c%p%pns%leafn_xfer
    frootn_xfer                    => clm3%g%l%c%p%pns%frootn_xfer
    livestemn_xfer                 => clm3%g%l%c%p%pns%livestemn_xfer
    deadstemn_xfer                 => clm3%g%l%c%p%pns%deadstemn_xfer
    livecrootn_xfer                => clm3%g%l%c%p%pns%livecrootn_xfer
    deadcrootn_xfer                => clm3%g%l%c%p%pns%deadcrootn_xfer
    leafc_storage_to_xfer          => clm3%g%l%c%p%pcf%leafc_storage_to_xfer
    frootc_storage_to_xfer         => clm3%g%l%c%p%pcf%frootc_storage_to_xfer
    livestemc_storage_to_xfer      => clm3%g%l%c%p%pcf%livestemc_storage_to_xfer
    deadstemc_storage_to_xfer      => clm3%g%l%c%p%pcf%deadstemc_storage_to_xfer
    livecrootc_storage_to_xfer     => clm3%g%l%c%p%pcf%livecrootc_storage_to_xfer
    deadcrootc_storage_to_xfer     => clm3%g%l%c%p%pcf%deadcrootc_storage_to_xfer
    gresp_storage_to_xfer          => clm3%g%l%c%p%pcf%gresp_storage_to_xfer
    leafn_storage_to_xfer          => clm3%g%l%c%p%pnf%leafn_storage_to_xfer
    frootn_storage_to_xfer         => clm3%g%l%c%p%pnf%frootn_storage_to_xfer
    livestemn_storage_to_xfer      => clm3%g%l%c%p%pnf%livestemn_storage_to_xfer
    deadstemn_storage_to_xfer      => clm3%g%l%c%p%pnf%deadstemn_storage_to_xfer
    livecrootn_storage_to_xfer     => clm3%g%l%c%p%pnf%livecrootn_storage_to_xfer
    deadcrootn_storage_to_xfer     => clm3%g%l%c%p%pnf%deadcrootn_storage_to_xfer

   ! set time steps
   dt = real( get_step_size() )
   fracday = dt/86400.0

   ! set some local parameters - these will be moved into
   ! parameter file after testing

   ! onset parameters
!  crit_onset_fdd = 15.0 ! gkw: default value
   crit_onset_fdd =  7.0 ! gkw: this may prevent "January thaw" growth spurt
   ! critical onset gdd now being calculated as a function of annual
   ! average 2m temp.
   ! crit_onset_gdd = 150.0 ! c3 grass value
   ! crit_onset_gdd = 1000.0   ! c4 grass value
   crit_onset_swi = 15.0
   soilpsi_on = -2.0
   ndays_on = 30.0

   ! offset parameters
   crit_offset_fdd = 15.0 ! gkw: changed from 15 to 30 days for run 68a 20110520
   crit_offset_swi = 15.0
   soilpsi_off = -2.0
   ndays_off = 15.0

   crit_dayl = 39300.

   ! transfer parameters
   fstor2tran = 0.5

   do fp = 1,num_soilp
      p = filter_soilp(fp)
      c = pcolumn(p)

      if (stress_decid(ivt(p)) == 1.) then
!!!!     soilt = t_soisno(c,2) ! gkw: hardwared code!
         soilt = max(t_grnd(c),t_soisno(c,1)) ! use TP1 or weighted TG & TP1
!!!!     psi = soilpsi(c,2)    ! gkw: hardwired code!
         psi = soilpsi(c,1)    ! use root-zone

         soilpsi_on  = psiwilt(c) ! gkw 20110510: use catchment wilting point
         soilpsi_off = psiwilt(c) ! gkw 20110510: use catchment wilting point

         ! set flag for solstice period (winter->summer = 1, summer->winter = 0)
         if (dayl(p) >= prev_dayl(p)) then
            ws_flag = 1.
         else
            ws_flag = 0.
         end if

         ! onset gdd sum from Biome-BGC, v4.1.2
         crit_onset_gdd = exp(4.8 + 0.13*(annavg_t2m(p) - SHR_CONST_TKFRZ))

         crit_offset_fdd = 15.0
         if(ivt(p)==14 .or. ivt(p)==18 .or. ivt(p)==10 .or. ivt(p)==16) crit_offset_fdd = 999.0 ! no T stress trigger
         if(ivt(p)==15 .or. ivt(p)==19 .or. ivt(p)==11 .or. ivt(p)==17) crit_offset_fdd = 999.0 ! no T stress trigger

         ! update offset_counter and test for the end of the offset period
         if (offset_flag(p) == 1.) then
            ! decrement counter for offset period
            offset_counter(p) = offset_counter(p) - dt

            ! if this is the end of the offset_period, reset phenology
            ! flags and indices
            if (offset_counter(p) <= 0.) then
               ! this code block was originally handled by call cn_offset_cleanup(p)
               ! inlined during vectorization
               offset_flag(p) = 0.
               offset_counter(p) = 0.
               dormant_flag(p) = 1.
               days_active(p) = 0.

               ! reset the previous timestep litterfall flux memory
               prev_leafc_to_litter(p) = 0.
               prev_frootc_to_litter(p) = 0.
            end if
         end if

         ! update onset_counter and test for the end of the onset period
         if (onset_flag(p) == 1.0) then
            ! decrement counter for onset period
            onset_counter(p) = onset_counter(p) - dt

            ! if this is the end of the onset period, reset phenology
            ! flags and indices
            if (onset_counter(p) <= 0.0) then
               ! this code block was originally handled by call cn_onset_cleanup(p)
               ! inlined during vectorization
               onset_flag(p) = 0.
               onset_counter(p) = 0.
               ! set all transfer growth rates to 0.0
               leafc_xfer_to_leafc(p)   = 0.
               frootc_xfer_to_frootc(p) = 0.
               leafn_xfer_to_leafn(p)   = 0.
               frootn_xfer_to_frootn(p) = 0.
               if (woody(ivt(p)) == 1.0) then
                  livestemc_xfer_to_livestemc(p)   = 0.
                  deadstemc_xfer_to_deadstemc(p)   = 0.
                  livecrootc_xfer_to_livecrootc(p) = 0.
                  deadcrootc_xfer_to_deadcrootc(p) = 0.
                  livestemn_xfer_to_livestemn(p)   = 0.
                  deadstemn_xfer_to_deadstemn(p)   = 0.
                  livecrootn_xfer_to_livecrootn(p) = 0.
                  deadcrootn_xfer_to_deadcrootn(p) = 0.
               end if
               ! set transfer pools to 0.0
               leafc_xfer(p) = 0.
               leafn_xfer(p) = 0.
               frootc_xfer(p) = 0.
               frootn_xfer(p) = 0.
               if (woody(ivt(p)) == 1.0) then
                  livestemc_xfer(p) = 0.
                  livestemn_xfer(p) = 0.
                  deadstemc_xfer(p) = 0.
                  deadstemn_xfer(p) = 0.
                  livecrootc_xfer(p) = 0.
                  livecrootn_xfer(p) = 0.
                  deadcrootc_xfer(p) = 0.
                  deadcrootn_xfer(p) = 0.
               end if
            end if
         end if

         ! test for switching from dormant period to growth period
         if (dormant_flag(p) == 1.) then

            ! keep track of the number of freezing degree days in this
            ! dormancy period (only if the freeze flag has not previously been set
            ! for this dormancy period

            if (onset_gddflag(p) == 0. .and. soilt < SHR_CONST_TKFRZ) onset_fdd(p) = onset_fdd(p) + fracday

            ! if the number of freezing degree days exceeds a critical value,
            ! then onset will require both wet soils and a critical soil
            ! temperature sum.  If this case is triggered, reset any previously
            ! accumulated value in onset_swi, so that onset now depends on
            ! the accumulated soil water index following the freeze trigger

            if (onset_fdd(p) > crit_onset_fdd) then
                onset_gddflag(p) = 1.
                onset_fdd(p) = 0.
                onset_swi(p) = 0.
            end if

            if (ivt(p)==14 .or. ivt(p)==18 .or. ivt(p)==10 .or. ivt(p)==16) then ! gkw: special case; seasonal deciduous

! after winter solstice, allow check for new growth
              if (onset_gddflag(p) == 0. .and. ws_flag == 1.) then
                 onset_gddflag(p) = 1.
                 onset_gdd(p) = 0.
                 onset_fdd(p) = 999.
                 onset_swi(p) = 0.
              end if

! before winter solstice, prevent growth onset
              if (ws_flag == 0.) then
                if (onset_flag(p) == 1. .or. dormant_flag(p) == 1. .or. onset_gddflag(p) == 1.) then
                  onset_flag(p) = 0.
                  onset_gddflag(p) = 0.
                  onset_gdd(p) = 0.
                  onset_fdd(p) = 999.
                  onset_swi(p) = 0.
                end if
	      endif

            endif

            ! if the freeze flag is set, and if the soil is above freezing
            ! then accumulate growing degree days for onset trigger

            if (onset_gddflag(p) == 1. .and. soilt > SHR_CONST_TKFRZ) then
               onset_gdd(p) = onset_gdd(p) + (soilt-SHR_CONST_TKFRZ)*fracday
            end if

            ! if soils are wet, accumulate soil water index for onset trigger
            if (psi > soilpsi_on) onset_swi(p) = onset_swi(p) + fracday

            ! if critical soil water index is exceeded, set onset_flag, and
            ! then test for soil temperature criteria

            if (onset_swi(p) > crit_onset_swi) then
                onset_flag(p) = 1.

                ! only check soil temperature criteria if freeze flag set since
                ! beginning of last dormancy.  If freeze flag set and growing
                ! degree day sum (since freeze trigger) is lower than critical
                ! value, then override the onset_flag set from soil water.

                if (onset_gddflag(p) == 1. .and. onset_gdd(p) < crit_onset_gdd) onset_flag(p) = 0.
            end if
            
            ! only allow onset if dayl > 6hrs
            if (onset_flag(p) == 1. .and. dayl(p) <= 21600.) then
                onset_flag(p) = 0.
            end if

            ! if this is the beginning of the onset period
            ! then reset the phenology flags and indices

            if (onset_flag(p) == 1.) then
               dormant_flag(p) = 0.
               days_active(p) = 0.
               onset_gddflag(p) = 0.
               onset_fdd(p) = 0.
               onset_gdd(p) = 0.
               onset_swi(p) = 0.
               onset_counter(p) = ndays_on * 86400.

               ! call subroutine to move all the storage pools into transfer pools,
               ! where they will be transfered to displayed growth over the onset period.
               ! this code was originally handled with call cn_storage_to_xfer(p)
               ! inlined during vectorization

               ! set carbon fluxes for shifting storage pools to transfer pools
               leafc_storage_to_xfer(p)  = fstor2tran * leafc_storage(p)/dt
               frootc_storage_to_xfer(p) = fstor2tran * frootc_storage(p)/dt
               if (woody(ivt(p)) == 1.0) then
                  livestemc_storage_to_xfer(p)  = fstor2tran * livestemc_storage(p)/dt
                  deadstemc_storage_to_xfer(p)  = fstor2tran * deadstemc_storage(p)/dt
                  livecrootc_storage_to_xfer(p) = fstor2tran * livecrootc_storage(p)/dt
                  deadcrootc_storage_to_xfer(p) = fstor2tran * deadcrootc_storage(p)/dt
                  gresp_storage_to_xfer(p)      = fstor2tran * gresp_storage(p)/dt
               end if

               ! set nitrogen fluxes for shifting storage pools to transfer pools
               leafn_storage_to_xfer(p)  = fstor2tran * leafn_storage(p)/dt
               frootn_storage_to_xfer(p) = fstor2tran * frootn_storage(p)/dt
               if (woody(ivt(p)) == 1.0) then
                  livestemn_storage_to_xfer(p)  = fstor2tran * livestemn_storage(p)/dt
                  deadstemn_storage_to_xfer(p)  = fstor2tran * deadstemn_storage(p)/dt
                  livecrootn_storage_to_xfer(p) = fstor2tran * livecrootn_storage(p)/dt
                  deadcrootn_storage_to_xfer(p) = fstor2tran * deadcrootn_storage(p)/dt
               end if
            end if

         ! test for switching from growth period to offset period
         else if (offset_flag(p) == 0.) then

            ! if soil water potential lower than critical value, accumulate
            ! as stress in offset soil water index

            if (psi <= soilpsi_off) then
               offset_swi(p) = offset_swi(p) + fracday

               ! if the offset soil water index exceeds critical value, and
               ! if this is not the middle of a previously initiated onset period,
               ! then set flag to start the offset period and reset index variables

               if (offset_swi(p) >= crit_offset_swi .and. onset_flag(p) == 0.) offset_flag(p) = 1.

            ! if soil water potential higher than critical value, reduce the
            ! offset water stress index.  By this mechanism, there must be a
            ! sustained period of water stress to initiate offset.

            else if (psi > soilpsi_on) then
               offset_swi(p) = offset_swi(p) - fracday
               offset_swi(p) = max(offset_swi(p),0.)
            end if

            ! decrease freezing day accumulator for warm soil
            if (offset_fdd(p) > 0. .and. soilt > SHR_CONST_TKFRZ) then
                offset_fdd(p) = offset_fdd(p) - fracday
                offset_fdd(p) = max(0., offset_fdd(p))
            end if

            ! increase freezing day accumulator for cold soil
            if (soilt <= SHR_CONST_TKFRZ) then
               offset_fdd(p) = offset_fdd(p) + fracday

               ! if freezing degree day sum is greater than critical value, initiate offset
               if (offset_fdd(p) > crit_offset_fdd .and. onset_flag(p) == 0.) offset_flag(p) = 1.
            end if
            
            ! force offset if daylength is < 6 hrs
            if (dayl(p) <= 21600.) then
            	offset_flag(p) = 1.
            end if

            ! only begin to test for offset daylength once past the summer sol
            if( ivt(p)==14 .or. ivt(p)==18 .or. ivt(p)==10 .or. ivt(p)==16) then ! gkw: special case
              if (ws_flag == 0. .and. dayl(p) < crit_dayl) then
                 offset_flag(p) = 1.
              end if
            endif

            ! if this is the beginning of the offset period
            ! then reset flags and indices
            if (offset_flag(p) == 1.) then
               offset_fdd(p) = 0.
               offset_swi(p) = 0.
               offset_counter(p) = ndays_off * 86400.
               prev_leafc_to_litter(p) = 0.
               prev_frootc_to_litter(p) = 0.
            end if
         end if

         ! keep track of number of days since last dormancy for control on
         ! fraction of new growth to send to storage for next growing season

         if (dormant_flag(p) == 0.0) then
             days_active(p) = days_active(p) + fracday
         end if

         ! calculate long growing season factor (lgsf)
         ! only begin to calculate a lgsf greater than 0.0 once the number
         ! of days active exceeds 365.
         lgsf(p) = max(min((days_active(p)-365.)/365., 1.),0.)

         ! set background litterfall rate, when not in the phenological offset period
         if (offset_flag(p) == 1.) then
            bglfr(p) = 0.
         else
            ! calculate the background litterfall rate (bglfr)
            ! in units 1/s, based on leaf longevity (yrs) and correction for long growing season

            bglfr(p) = (1./(leaf_long(ivt(p))*365.*86400.))*lgsf(p)
         end if

         ! set background transfer rate when active but not in the phenological onset period
         if (onset_flag(p) == 1.) then
            bgtr(p) = 0.
         else
            ! the background transfer rate is calculated as the rate that would result
            ! in complete turnover of the storage pools in one year at steady state,
            ! once lgsf has reached 1.0 (after 730 days active).

            bgtr(p) = (1./(365.*86400.))*lgsf(p)

            ! set carbon fluxes for shifting storage pools to transfer pools

            leafc_storage_to_xfer(p)  = leafc_storage(p) * bgtr(p)
            frootc_storage_to_xfer(p) = frootc_storage(p) * bgtr(p)
            if (woody(ivt(p)) == 1.0) then
               livestemc_storage_to_xfer(p)  = livestemc_storage(p) * bgtr(p)
               deadstemc_storage_to_xfer(p)  = deadstemc_storage(p) * bgtr(p)
               livecrootc_storage_to_xfer(p) = livecrootc_storage(p) * bgtr(p)
               deadcrootc_storage_to_xfer(p) = deadcrootc_storage(p) * bgtr(p)
               gresp_storage_to_xfer(p)      = gresp_storage(p) * bgtr(p)
            end if

            ! set nitrogen fluxes for shifting storage pools to transfer pools
            leafn_storage_to_xfer(p)  = leafn_storage(p) * bgtr(p)
            frootn_storage_to_xfer(p) = frootn_storage(p) * bgtr(p)
            if (woody(ivt(p)) == 1.0) then
               livestemn_storage_to_xfer(p)  = livestemn_storage(p) * bgtr(p)
               deadstemn_storage_to_xfer(p)  = deadstemn_storage(p) * bgtr(p)
               livecrootn_storage_to_xfer(p) = livecrootn_storage(p) * bgtr(p)
               deadcrootn_storage_to_xfer(p) = deadcrootn_storage(p) * bgtr(p)
            end if
         end if

      end if ! end if stress deciduous

   end do ! end of pft loop

end subroutine CNStressDecidPhenology
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNOnsetGrowth
!
! !INTERFACE:
subroutine CNOnsetGrowth (num_soilp, filter_soilp)
!
! !DESCRIPTION:
! Determines the flux of stored C and N from transfer pools to display
! pools during the phenological onset period.
!
! !USES:
   use clm_time_manager, only: get_step_size
!
! !ARGUMENTS:
   integer, intent(in) :: num_soilp       ! number of soil pfts in filter
   integer, intent(in) :: filter_soilp(:) ! filter for soil pfts
!
! !CALLED FROM:
! subroutine CNPhenology
!
! !REVISION HISTORY:
! 10/27/03: Created by Peter Thornton
!
! !LOCAL VARIABLES:
! local pointers to implicit in scalars
!
   integer , pointer :: ivt(:)             ! pft vegetation type
   real, pointer :: onset_flag(:)      ! onset flag
   real, pointer :: onset_counter(:)   ! onset days counter
   real, pointer :: leafc_xfer(:)      ! (kgC/m2) leaf C transfer
   real, pointer :: frootc_xfer(:)     ! (kgC/m2) fine root C transfer
   real, pointer :: livestemc_xfer(:)  ! (kgC/m2) live stem C transfer
   real, pointer :: deadstemc_xfer(:)  ! (kgC/m2) dead stem C transfer
   real, pointer :: livecrootc_xfer(:) ! (kgC/m2) live coarse root C transfer
   real, pointer :: deadcrootc_xfer(:) ! (kgC/m2) dead coarse root C transfer
   real, pointer :: leafn_xfer(:)      ! (kgN/m2) leaf N transfer
   real, pointer :: frootn_xfer(:)     ! (kgN/m2) fine root N transfer
   real, pointer :: livestemn_xfer(:)  ! (kgN/m2) live stem N transfer
   real, pointer :: deadstemn_xfer(:)  ! (kgN/m2) dead stem N transfer
   real, pointer :: livecrootn_xfer(:) ! (kgN/m2) live coarse root N transfer
   real, pointer :: deadcrootn_xfer(:) ! (kgN/m2) dead coarse root N transfer
   real, pointer :: woody(:)           ! binary flag for woody lifeform (1=woody, 0=not woody)
   real, pointer :: bgtr(:)            ! background transfer growth rate (1/s)
!
! local pointers to implicit in/out scalars
!
   real, pointer :: leafc_xfer_to_leafc(:)
   real, pointer :: frootc_xfer_to_frootc(:)
   real, pointer :: livestemc_xfer_to_livestemc(:)
   real, pointer :: deadstemc_xfer_to_deadstemc(:)
   real, pointer :: livecrootc_xfer_to_livecrootc(:)
   real, pointer :: deadcrootc_xfer_to_deadcrootc(:)
   real, pointer :: leafn_xfer_to_leafn(:)
   real, pointer :: frootn_xfer_to_frootn(:)
   real, pointer :: livestemn_xfer_to_livestemn(:)
   real, pointer :: deadstemn_xfer_to_deadstemn(:)
   real, pointer :: livecrootn_xfer_to_livecrootn(:)
   real, pointer :: deadcrootn_xfer_to_deadcrootn(:)
!
! local pointers to implicit out scalars
!
! !OTHER LOCAL VARIABLES:
   integer :: p            ! indices
   integer :: fp           ! lake filter pft index
   real:: dt           ! radiation time step delta t (seconds)
   real:: t1           ! temporary variable

!EOP
!-----------------------------------------------------------------------
   ! assign local pointers to derived type arrays (in)
    ivt                            => clm3%g%l%c%p%itype
    onset_flag                     => clm3%g%l%c%p%pepv%onset_flag
    onset_counter                  => clm3%g%l%c%p%pepv%onset_counter
    leafc_xfer                     => clm3%g%l%c%p%pcs%leafc_xfer
    frootc_xfer                    => clm3%g%l%c%p%pcs%frootc_xfer
    livestemc_xfer                 => clm3%g%l%c%p%pcs%livestemc_xfer
    deadstemc_xfer                 => clm3%g%l%c%p%pcs%deadstemc_xfer
    livecrootc_xfer                => clm3%g%l%c%p%pcs%livecrootc_xfer
    deadcrootc_xfer                => clm3%g%l%c%p%pcs%deadcrootc_xfer
    leafn_xfer                     => clm3%g%l%c%p%pns%leafn_xfer
    frootn_xfer                    => clm3%g%l%c%p%pns%frootn_xfer
    livestemn_xfer                 => clm3%g%l%c%p%pns%livestemn_xfer
    deadstemn_xfer                 => clm3%g%l%c%p%pns%deadstemn_xfer
    livecrootn_xfer                => clm3%g%l%c%p%pns%livecrootn_xfer
    deadcrootn_xfer                => clm3%g%l%c%p%pns%deadcrootn_xfer
    bgtr                           => clm3%g%l%c%p%pepv%bgtr
    woody                          => pftcon%woody

   ! assign local pointers to derived type arrays (out)
    leafc_xfer_to_leafc            => clm3%g%l%c%p%pcf%leafc_xfer_to_leafc
    frootc_xfer_to_frootc          => clm3%g%l%c%p%pcf%frootc_xfer_to_frootc
    livestemc_xfer_to_livestemc    => clm3%g%l%c%p%pcf%livestemc_xfer_to_livestemc
    deadstemc_xfer_to_deadstemc    => clm3%g%l%c%p%pcf%deadstemc_xfer_to_deadstemc
    livecrootc_xfer_to_livecrootc  => clm3%g%l%c%p%pcf%livecrootc_xfer_to_livecrootc
    deadcrootc_xfer_to_deadcrootc  => clm3%g%l%c%p%pcf%deadcrootc_xfer_to_deadcrootc
    leafn_xfer_to_leafn            => clm3%g%l%c%p%pnf%leafn_xfer_to_leafn
    frootn_xfer_to_frootn          => clm3%g%l%c%p%pnf%frootn_xfer_to_frootn
    livestemn_xfer_to_livestemn    => clm3%g%l%c%p%pnf%livestemn_xfer_to_livestemn
    deadstemn_xfer_to_deadstemn    => clm3%g%l%c%p%pnf%deadstemn_xfer_to_deadstemn
    livecrootn_xfer_to_livecrootn  => clm3%g%l%c%p%pnf%livecrootn_xfer_to_livecrootn
    deadcrootn_xfer_to_deadcrootn  => clm3%g%l%c%p%pnf%deadcrootn_xfer_to_deadcrootn

   ! set time steps
   dt = real( get_step_size() )

   ! pft loop
   do fp = 1,num_soilp
      p = filter_soilp(fp)

      ! only calculate these fluxes during onset period
      if (onset_flag(p) == 1.) then

         ! The transfer rate is a linearly decreasing function of time,
         ! going to zero on the last timestep of the onset period

         if (onset_counter(p) == dt) then
             t1 = 1.0 / dt
         else
             t1 = 2.0 / (onset_counter(p))
         end if
         leafc_xfer_to_leafc(p)   = t1 * leafc_xfer(p)
         frootc_xfer_to_frootc(p) = t1 * frootc_xfer(p)
         leafn_xfer_to_leafn(p)   = t1 * leafn_xfer(p)
         frootn_xfer_to_frootn(p) = t1 * frootn_xfer(p)
         if (woody(ivt(p)) == 1.0) then
             livestemc_xfer_to_livestemc(p)   = t1 * livestemc_xfer(p)
             deadstemc_xfer_to_deadstemc(p)   = t1 * deadstemc_xfer(p)
             livecrootc_xfer_to_livecrootc(p) = t1 * livecrootc_xfer(p)
             deadcrootc_xfer_to_deadcrootc(p) = t1 * deadcrootc_xfer(p)
             livestemn_xfer_to_livestemn(p)   = t1 * livestemn_xfer(p)
             deadstemn_xfer_to_deadstemn(p)   = t1 * deadstemn_xfer(p)
             livecrootn_xfer_to_livecrootn(p) = t1 * livecrootn_xfer(p)
             deadcrootn_xfer_to_deadcrootn(p) = t1 * deadcrootn_xfer(p)
         end if

      end if ! end if onset period

      ! calculate the background rate of transfer growth (used for stress
      ! deciduous algorithm). In this case, all of the mass in the transfer
      ! pools should be moved to displayed growth in each timestep.

      if (bgtr(p) > 0.) then
         leafc_xfer_to_leafc(p)   = leafc_xfer(p) / dt
         frootc_xfer_to_frootc(p) = frootc_xfer(p) / dt
         leafn_xfer_to_leafn(p)   = leafn_xfer(p) / dt
         frootn_xfer_to_frootn(p) = frootn_xfer(p) / dt
         if (woody(ivt(p)) == 1.0) then
             livestemc_xfer_to_livestemc(p)   = livestemc_xfer(p) / dt
             deadstemc_xfer_to_deadstemc(p)   = deadstemc_xfer(p) / dt
             livecrootc_xfer_to_livecrootc(p) = livecrootc_xfer(p) / dt
             deadcrootc_xfer_to_deadcrootc(p) = deadcrootc_xfer(p) / dt
             livestemn_xfer_to_livestemn(p)   = livestemn_xfer(p) / dt
             deadstemn_xfer_to_deadstemn(p)   = deadstemn_xfer(p) / dt
             livecrootn_xfer_to_livecrootn(p) = livecrootn_xfer(p) / dt
             deadcrootn_xfer_to_deadcrootn(p) = deadcrootn_xfer(p) / dt
         end if
      end if ! end if bgtr

   end do ! end pft loop

end subroutine CNOnsetGrowth
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNOffsetLitterfall
!
! !INTERFACE:
subroutine CNOffsetLitterfall (num_soilp, filter_soilp)
!
! !DESCRIPTION:
! Determines the flux of C and N from displayed pools to litter
! pools during the phenological offset period.
!
! !USES:
   use clm_time_manager, only: get_step_size
!
! !ARGUMENTS:
   integer, intent(in) :: num_soilp       ! number of soil pfts in filter
   integer, intent(in) :: filter_soilp(:) ! filter for soil pfts
!
! !CALLED FROM:
! subroutine CNPhenology
!
! !REVISION HISTORY:
! 10/27/03: Created by Peter Thornton
!
! !LOCAL VARIABLES:
! local pointers to implicit in scalars
!
   integer , pointer :: ivt(:)                   ! pft vegetation type
   real, pointer :: offset_flag(:)           ! offset flag
   real, pointer :: offset_counter(:)        ! offset days counter
   real, pointer :: leafc(:)                 ! (kgC/m2) leaf C
   real, pointer :: frootc(:)                ! (kgC/m2) fine root C
   real, pointer :: cpool_to_leafc(:)
   real, pointer :: cpool_to_frootc(:)
   real, pointer :: leafcn(:)                ! leaf C:N (gC/gN)
   real, pointer :: lflitcn(:)               ! leaf litter C:N (gC/gN)
   real, pointer :: frootcn(:)               ! fine root C:N (gC/gN)
!
! local pointers to implicit in/out scalars
!
   real, pointer :: prev_leafc_to_litter(:)  ! previous timestep leaf C litterfall flux (gC/m2/s)
   real, pointer :: prev_frootc_to_litter(:) ! previous timestep froot C litterfall flux (gC/m2/s)
   real, pointer :: leafc_to_litter(:)
   real, pointer :: frootc_to_litter(:)
   real, pointer :: leafn_to_litter(:)
   real, pointer :: leafn_to_retransn(:)
   real, pointer :: frootn_to_litter(:)
!
! local pointers to implicit out scalars
!
!
! !OTHER LOCAL VARIABLES:
   integer :: p, c         ! indices
   integer :: fp           ! lake filter pft index
   real:: dt           ! radiation time step delta t (seconds)
   real:: t1           ! temporary variable

!EOP
!-----------------------------------------------------------------------
   ! assign local pointers to derived type arrays (in)
    ivt                            => clm3%g%l%c%p%itype
    offset_flag                    => clm3%g%l%c%p%pepv%offset_flag
    offset_counter                 => clm3%g%l%c%p%pepv%offset_counter
    leafc                          => clm3%g%l%c%p%pcs%leafc
    frootc                         => clm3%g%l%c%p%pcs%frootc
    cpool_to_leafc                 => clm3%g%l%c%p%pcf%cpool_to_leafc
    cpool_to_frootc                => clm3%g%l%c%p%pcf%cpool_to_frootc
    leafcn                         => pftcon%leafcn
    lflitcn                        => pftcon%lflitcn
    frootcn                        => pftcon%frootcn

   ! assign local pointers to derived type arrays (out)
    prev_leafc_to_litter           => clm3%g%l%c%p%pepv%prev_leafc_to_litter
    prev_frootc_to_litter          => clm3%g%l%c%p%pepv%prev_frootc_to_litter
    leafc_to_litter                => clm3%g%l%c%p%pcf%leafc_to_litter
    frootc_to_litter               => clm3%g%l%c%p%pcf%frootc_to_litter
    leafn_to_litter                => clm3%g%l%c%p%pnf%leafn_to_litter
    leafn_to_retransn              => clm3%g%l%c%p%pnf%leafn_to_retransn
    frootn_to_litter               => clm3%g%l%c%p%pnf%frootn_to_litter

   ! set time steps
   dt = real( get_step_size() )

   ! The litterfall transfer rate starts at 0.0 and increases linearly
   ! over time, with displayed growth going to 0.0 on the last day of litterfall

   do fp = 1,num_soilp
      p = filter_soilp(fp)

      ! only calculate fluxes during offset period
      if (offset_flag(p) == 1.) then

         if (offset_counter(p) == dt) then
             t1 = 1.0 / dt
             leafc_to_litter(p)  = t1 * leafc(p)  + cpool_to_leafc(p)
             frootc_to_litter(p) = t1 * frootc(p) + cpool_to_frootc(p)
         else
             t1 = dt * 2.0 / (offset_counter(p) * offset_counter(p))
             leafc_to_litter(p)  = prev_leafc_to_litter(p)  + t1*(leafc(p)  - prev_leafc_to_litter(p)*offset_counter(p))
             frootc_to_litter(p) = prev_frootc_to_litter(p) + t1*(frootc(p) - prev_frootc_to_litter(p)*offset_counter(p))
         end if

         ! calculate the leaf N litterfall and retranslocation
         leafn_to_litter(p)   = leafc_to_litter(p)  / lflitcn(ivt(p))
         leafn_to_retransn(p) = (leafc_to_litter(p) / leafcn(ivt(p))) - leafn_to_litter(p)

         ! calculate fine root N litterfall (no retranslocation of fine root N)
         frootn_to_litter(p) = frootc_to_litter(p) / frootcn(ivt(p))

         ! save the current litterfall fluxes
         prev_leafc_to_litter(p)  = leafc_to_litter(p)
         prev_frootc_to_litter(p) = frootc_to_litter(p)

      end if ! end if offset period

   end do ! end pft loop

end subroutine CNOffsetLitterfall
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNBackgroundLitterfall
!
! !INTERFACE:
subroutine CNBackgroundLitterfall (num_soilp, filter_soilp)
!
! !DESCRIPTION:
! Determines the flux of C and N from displayed pools to litter
! pools as the result of background litter fall.
!
! !USES:
   use clm_time_manager, only: get_step_size
!
! !ARGUMENTS:
   integer, intent(in) :: num_soilp       ! number of soil pfts in filter
   integer, intent(in) :: filter_soilp(:) ! filter for soil pfts
!
! !CALLED FROM:
! subroutine CNPhenology
!
! !REVISION HISTORY:
! 10/2/03: Created by Peter Thornton
! 10/24/03, Peter Thornton: migrated to vector data structures
!
! !LOCAL VARIABLES:
! local pointers to implicit in scalars
!
   ! pft level
   integer , pointer :: ivt(:)       ! pft vegetation type
   real, pointer :: bglfr(:)     ! background litterfall rate (1/s)
   real, pointer :: leafc(:)     ! (kgC/m2) leaf C
   real, pointer :: frootc(:)    ! (kgC/m2) fine root C
   ! ecophysiological constants
   real, pointer :: leafcn(:)    ! leaf C:N (gC/gN)
   real, pointer :: lflitcn(:)   ! leaf litter C:N (gC/gN)
   real, pointer :: frootcn(:)   ! fine root C:N (gC/gN)
!
! local pointers to implicit in/out scalars
!
   real, pointer :: leafc_to_litter(:)
   real, pointer :: frootc_to_litter(:)
   real, pointer :: leafn_to_litter(:)
   real, pointer :: leafn_to_retransn(:)
   real, pointer :: frootn_to_litter(:)
!
! local pointers to implicit out scalars
!
!
! !OTHER LOCAL VARIABLES:
   integer :: p            ! indices
   integer :: fp           ! lake filter pft index
   real:: dt           ! decomp timestep (seconds)

!EOP
!-----------------------------------------------------------------------
   ! assign local pointers to derived type arrays (in)
    ivt                            => clm3%g%l%c%p%itype
    bglfr                          => clm3%g%l%c%p%pepv%bglfr
    leafc                          => clm3%g%l%c%p%pcs%leafc
    frootc                         => clm3%g%l%c%p%pcs%frootc
    leafcn                         => pftcon%leafcn
    lflitcn                        => pftcon%lflitcn
    frootcn                        => pftcon%frootcn

   ! assign local pointers to derived type arrays (out)
    leafc_to_litter                => clm3%g%l%c%p%pcf%leafc_to_litter
    frootc_to_litter               => clm3%g%l%c%p%pcf%frootc_to_litter
    leafn_to_litter                => clm3%g%l%c%p%pnf%leafn_to_litter
    leafn_to_retransn              => clm3%g%l%c%p%pnf%leafn_to_retransn
    frootn_to_litter               => clm3%g%l%c%p%pnf%frootn_to_litter

   ! set time steps
   dt = real( get_step_size() )

   ! pft loop
   do fp = 1,num_soilp
      p = filter_soilp(fp)

      ! only calculate these fluxes if the background litterfall rate is non-zero
      if (bglfr(p) > 0.) then
         ! units for bglfr are already 1/s
         leafc_to_litter(p)  = bglfr(p) * leafc(p)
         frootc_to_litter(p) = bglfr(p) * frootc(p)

         ! calculate the leaf N litterfall and retranslocation
         leafn_to_litter(p)   = leafc_to_litter(p)  / lflitcn(ivt(p))
         leafn_to_retransn(p) = (leafc_to_litter(p) / leafcn(ivt(p))) - leafn_to_litter(p)

         ! calculate fine root N litterfall (no retranslocation of fine root N)
         frootn_to_litter(p) = frootc_to_litter(p) / frootcn(ivt(p))

      end if

   end do

end subroutine CNBackgroundLitterfall
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNLivewoodTurnover
!
! !INTERFACE:
subroutine CNLivewoodTurnover (num_soilp, filter_soilp)
!
! !DESCRIPTION:
! Determines the flux of C and N from live wood to
! dead wood pools, for stem and coarse root.
!
! !USES:
   use clm_time_manager, only: get_step_size
!
! !ARGUMENTS:
   integer, intent(in) :: num_soilp       ! number of soil pfts in filter
   integer, intent(in) :: filter_soilp(:) ! filter for soil pfts
!
! !CALLED FROM:
! subroutine CNPhenology
!
! !REVISION HISTORY:
! 12/5/03: created by Peter Thornton
!
! !LOCAL VARIABLES:
! local pointers to implicit in scalars
!
   ! pft level
   integer , pointer :: ivt(:)         ! pft vegetation type
   real, pointer :: livestemc(:)   ! (gC/m2) live stem C
   real, pointer :: livecrootc(:)  ! (gC/m2) live coarse root C
   real, pointer :: livestemn(:)   ! (gN/m2) live stem N
   real, pointer :: livecrootn(:)  ! (gN/m2) live coarse root N
   ! ecophysiological constants
   real, pointer :: woody(:)       ! binary flag for woody lifeform (1=woody, 0=not woody)
   real, pointer :: livewdcn(:)    ! live wood (phloem and ray parenchyma) C:N (gC/gN)
   real, pointer :: deadwdcn(:)    ! dead wood (xylem and heartwood) C:N (gC/gN)
!
! local pointers to implicit in/out scalars
!
   real, pointer :: livestemc_to_deadstemc(:)
   real, pointer :: livecrootc_to_deadcrootc(:)
   real, pointer :: livestemn_to_deadstemn(:)
   real, pointer :: livestemn_to_retransn(:)
   real, pointer :: livecrootn_to_deadcrootn(:)
   real, pointer :: livecrootn_to_retransn(:)
!
! local pointers to implicit out scalars
!
!
! !OTHER LOCAL VARIABLES:
   integer :: p            ! indices
   integer :: fp           ! lake filter pft index
   real:: dt           ! decomp timestep (seconds)
   real:: lwtop        ! live wood turnover proportion (annual fraction)
   real:: ctovr        ! temporary variable for carbon turnover
   real:: ntovr        ! temporary variable for nitrogen turnover

!EOP
!-----------------------------------------------------------------------
   ! assign local pointers to derived type arrays (in)
    ivt                            => clm3%g%l%c%p%itype
    livestemc                      => clm3%g%l%c%p%pcs%livestemc
    livecrootc                     => clm3%g%l%c%p%pcs%livecrootc
    livestemn                      => clm3%g%l%c%p%pns%livestemn
    livecrootn                     => clm3%g%l%c%p%pns%livecrootn
    woody                          => pftcon%woody
    livewdcn                       => pftcon%livewdcn
    deadwdcn                       => pftcon%deadwdcn

   ! assign local pointers to derived type arrays (out)
    livestemc_to_deadstemc         => clm3%g%l%c%p%pcf%livestemc_to_deadstemc
    livecrootc_to_deadcrootc       => clm3%g%l%c%p%pcf%livecrootc_to_deadcrootc
    livestemn_to_deadstemn         => clm3%g%l%c%p%pnf%livestemn_to_deadstemn
    livestemn_to_retransn          => clm3%g%l%c%p%pnf%livestemn_to_retransn
    livecrootn_to_deadcrootn       => clm3%g%l%c%p%pnf%livecrootn_to_deadcrootn
    livecrootn_to_retransn         => clm3%g%l%c%p%pnf%livecrootn_to_retransn

   ! set time steps
   dt = real( get_step_size() )

   ! set the global parameter for livewood turnover rate
   ! define as an annual fraction (0.7), and convert to fraction per second
   lwtop = 0.7 / 31536000.0

   ! pft loop
   do fp = 1,num_soilp
      p = filter_soilp(fp)

      ! only calculate these fluxes for woody types
      if (woody(ivt(p)) > 0.) then

         ! live stem to dead stem turnover

         ctovr = livestemc(p) * lwtop
         ntovr = ctovr / livewdcn(ivt(p))
         livestemc_to_deadstemc(p) = ctovr
         livestemn_to_deadstemn(p) = ctovr / deadwdcn(ivt(p))
         livestemn_to_retransn(p)  = ntovr - livestemn_to_deadstemn(p)

         ! live coarse root to dead coarse root turnover

         ctovr = livecrootc(p) * lwtop
         ntovr = ctovr / livewdcn(ivt(p))
         livecrootc_to_deadcrootc(p) = ctovr
         livecrootn_to_deadcrootn(p) = ctovr / deadwdcn(ivt(p))
         livecrootn_to_retransn(p)  = ntovr - livecrootn_to_deadcrootn(p)

      end if

   end do

end subroutine CNLivewoodTurnover
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNLitterToColumn
!
! !INTERFACE:
subroutine CNLitterToColumn (num_soilc, filter_soilc)
!
! !DESCRIPTION:
! called at the end of cn_phenology to gather all pft-level litterfall fluxes
! to the column level and assign them to the three litter pools
!
! !USES:
  use clm_varpar, only : max_pft_per_col
!
! !ARGUMENTS:
  integer, intent(in) :: num_soilc       ! number of soil columns in filter
  integer, intent(in) :: filter_soilc(:) ! filter for soil columns
!
! !CALLED FROM:
! subroutine CNPhenology
!
! !REVISION HISTORY:
! 9/8/03: Created by Peter Thornton
!
! !LOCAL VARIABLES:
! local pointers to implicit in scalars
!
   integer , pointer :: ivt(:)          ! pft vegetation type
   real, pointer :: wtcol(:)        ! weight (relative to column) for this pft (0-1)
   real, pointer :: pwtgcell(:)     ! weight of pft relative to corresponding gridcell
   real, pointer :: leafc_to_litter(:)
   real, pointer :: frootc_to_litter(:)
   real, pointer :: leafn_to_litter(:)
   real, pointer :: frootn_to_litter(:)
   real, pointer :: lf_flab(:)      ! leaf litter labile fraction
   real, pointer :: lf_fcel(:)      ! leaf litter cellulose fraction
   real, pointer :: lf_flig(:)      ! leaf litter lignin fraction
   real, pointer :: fr_flab(:)      ! fine root litter labile fraction
   real, pointer :: fr_fcel(:)      ! fine root litter cellulose fraction
   real, pointer :: fr_flig(:)      ! fine root litter lignin fraction
   integer , pointer :: npfts(:)        ! number of pfts for each column
   integer , pointer :: pfti(:)         ! beginning pft index for each column
!
! local pointers to implicit in/out scalars
!
   real, pointer :: leafc_to_litr1c(:)
   real, pointer :: leafc_to_litr2c(:)
   real, pointer :: leafc_to_litr3c(:)
   real, pointer :: frootc_to_litr1c(:)
   real, pointer :: frootc_to_litr2c(:)
   real, pointer :: frootc_to_litr3c(:)
   real, pointer :: leafn_to_litr1n(:)
   real, pointer :: leafn_to_litr2n(:)
   real, pointer :: leafn_to_litr3n(:)
   real, pointer :: frootn_to_litr1n(:)
   real, pointer :: frootn_to_litr2n(:)
   real, pointer :: frootn_to_litr3n(:)
!
! local pointers to implicit out scalars
!
!
! !OTHER LOCAL VARIABLES:
    integer :: fc,c,pi,p
!EOP
!-----------------------------------------------------------------------
   ! assign local pointers to derived type arrays (in)
    ivt                            => clm3%g%l%c%p%itype
    wtcol                          => clm3%g%l%c%p%wtcol
    pwtgcell                       => clm3%g%l%c%p%wtgcell  
    leafc_to_litter                => clm3%g%l%c%p%pcf%leafc_to_litter
    frootc_to_litter               => clm3%g%l%c%p%pcf%frootc_to_litter
    leafn_to_litter                => clm3%g%l%c%p%pnf%leafn_to_litter
    frootn_to_litter               => clm3%g%l%c%p%pnf%frootn_to_litter
    npfts                          => clm3%g%l%c%npfts
    pfti                           => clm3%g%l%c%pfti
    lf_flab                        => pftcon%lf_flab
    lf_fcel                        => pftcon%lf_fcel
    lf_flig                        => pftcon%lf_flig
    fr_flab                        => pftcon%fr_flab
    fr_fcel                        => pftcon%fr_fcel
    fr_flig                        => pftcon%fr_flig

   ! assign local pointers to derived type arrays (out)
    leafc_to_litr1c                => clm3%g%l%c%ccf%leafc_to_litr1c
    leafc_to_litr2c                => clm3%g%l%c%ccf%leafc_to_litr2c
    leafc_to_litr3c                => clm3%g%l%c%ccf%leafc_to_litr3c
    frootc_to_litr1c               => clm3%g%l%c%ccf%frootc_to_litr1c
    frootc_to_litr2c               => clm3%g%l%c%ccf%frootc_to_litr2c
    frootc_to_litr3c               => clm3%g%l%c%ccf%frootc_to_litr3c
    leafn_to_litr1n                => clm3%g%l%c%cnf%leafn_to_litr1n
    leafn_to_litr2n                => clm3%g%l%c%cnf%leafn_to_litr2n
    leafn_to_litr3n                => clm3%g%l%c%cnf%leafn_to_litr3n
    frootn_to_litr1n               => clm3%g%l%c%cnf%frootn_to_litr1n
    frootn_to_litr2n               => clm3%g%l%c%cnf%frootn_to_litr2n
    frootn_to_litr3n               => clm3%g%l%c%cnf%frootn_to_litr3n

   do pi = 1,max_pft_per_col
      do fc = 1,num_soilc
         c = filter_soilc(fc)

         if ( pi <=  npfts(c) ) then
            p = pfti(c) + pi - 1
            if (pwtgcell(p)>0.) then

               ! leaf litter carbon fluxes
               leafc_to_litr1c(c) = leafc_to_litr1c(c) + leafc_to_litter(p) * lf_flab(ivt(p)) * wtcol(p)
               leafc_to_litr2c(c) = leafc_to_litr2c(c) + leafc_to_litter(p) * lf_fcel(ivt(p)) * wtcol(p)
               leafc_to_litr3c(c) = leafc_to_litr3c(c) + leafc_to_litter(p) * lf_flig(ivt(p)) * wtcol(p)

               ! leaf litter nitrogen fluxes
               leafn_to_litr1n(c) = leafn_to_litr1n(c) + leafn_to_litter(p) * lf_flab(ivt(p)) * wtcol(p)
               leafn_to_litr2n(c) = leafn_to_litr2n(c) + leafn_to_litter(p) * lf_fcel(ivt(p)) * wtcol(p)
               leafn_to_litr3n(c) = leafn_to_litr3n(c) + leafn_to_litter(p) * lf_flig(ivt(p)) * wtcol(p)

               ! fine root litter carbon fluxes
               frootc_to_litr1c(c) = frootc_to_litr1c(c) + frootc_to_litter(p) * fr_flab(ivt(p)) * wtcol(p)
               frootc_to_litr2c(c) = frootc_to_litr2c(c) + frootc_to_litter(p) * fr_fcel(ivt(p)) * wtcol(p)
               frootc_to_litr3c(c) = frootc_to_litr3c(c) + frootc_to_litter(p) * fr_flig(ivt(p)) * wtcol(p)

               ! fine root litter nitrogen fluxes
               frootn_to_litr1n(c) = frootn_to_litr1n(c) + frootn_to_litter(p) * fr_flab(ivt(p)) * wtcol(p)
               frootn_to_litr2n(c) = frootn_to_litr2n(c) + frootn_to_litter(p) * fr_fcel(ivt(p)) * wtcol(p)
               frootn_to_litr3n(c) = frootn_to_litr3n(c) + frootn_to_litter(p) * fr_flig(ivt(p)) * wtcol(p)

            end if
         end if

      end do

   end do

end subroutine CNLitterToColumn
!-----------------------------------------------------------------------

end module CNPhenologyMod
