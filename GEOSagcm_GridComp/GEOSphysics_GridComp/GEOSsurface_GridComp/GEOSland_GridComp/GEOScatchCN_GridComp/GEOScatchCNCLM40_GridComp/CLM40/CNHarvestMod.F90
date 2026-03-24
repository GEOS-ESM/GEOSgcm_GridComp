module CNHarvestMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: CNHarvestMod
!
! !DESCRIPTION:
! Harvest mortality routine for coupled carbon-nitrogen code (CN)
!
! !USES:
    implicit none
    save
    private
! !PUBLIC MEMBER FUNCTIONS:
    public :: CNHarvest
!
! !REVISION HISTORY:
! 3/29/2004: Created by Peter Thornton
!
!EOP

! ! PRIVATE TYPES
! real, parameter :: days_per_year = 365.
! integer , pointer   :: yearspft(:)
! real, pointer   :: wtpft1(:,:)   
! real, pointer   :: wtpft2(:,:)
! real, pointer   :: harvest(:)   
! real, pointer   :: wtcol_old(:)
! integer :: nt1
! integer :: nt2
! integer :: ntimes
! logical :: do_harvest

!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNHarvest
!
! !INTERFACE:
subroutine CNHarvest (num_soilc, filter_soilc, num_soilp, filter_soilp)
!
! !DESCRIPTION:
! Harvest mortality routine for coupled carbon-nitrogen code (CN)
!
! !USES:
   use clmtype
   use pftvarcon, only : noveg, nbrdlf_evr_shrub
!
! !ARGUMENTS:
   implicit none
   integer, intent(in) :: num_soilc       ! number of soil columns in filter
   integer, intent(in) :: filter_soilc(:) ! column filter for soil points
   integer, intent(in) :: num_soilp       ! number of soil pfts in filter
   integer, intent(in) :: filter_soilp(:) ! pft filter for soil points
!
! !CALLED FROM:
! subroutine CNEcosystemDyn
!
! !REVISION HISTORY:
! 3/29/04: Created by Peter Thornton
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arrays
   integer , pointer :: pgridcell(:)   ! pft-level index into gridcell-level quantities
   integer , pointer :: ivt(:)         ! pft vegetation type

   real, pointer :: leafc(:)              ! (gC/m2) leaf C
   real, pointer :: frootc(:)             ! (gC/m2) fine root C
   real, pointer :: livestemc(:)          ! (gC/m2) live stem C
   real, pointer :: deadstemc(:)          ! (gC/m2) dead stem C
   real, pointer :: livecrootc(:)         ! (gC/m2) live coarse root C
   real, pointer :: deadcrootc(:)         ! (gC/m2) dead coarse root C
   real, pointer :: xsmrpool(:)           ! (gC/m2) abstract C pool to meet excess MR demand
   real, pointer :: leafc_storage(:)      ! (gC/m2) leaf C storage
   real, pointer :: frootc_storage(:)     ! (gC/m2) fine root C storage
   real, pointer :: livestemc_storage(:)  ! (gC/m2) live stem C storage
   real, pointer :: deadstemc_storage(:)  ! (gC/m2) dead stem C storage
   real, pointer :: livecrootc_storage(:) ! (gC/m2) live coarse root C storage
   real, pointer :: deadcrootc_storage(:) ! (gC/m2) dead coarse root C storage
   real, pointer :: gresp_storage(:)      ! (gC/m2) growth respiration storage
   real, pointer :: leafc_xfer(:)         ! (gC/m2) leaf C transfer
   real, pointer :: frootc_xfer(:)        ! (gC/m2) fine root C transfer
   real, pointer :: livestemc_xfer(:)     ! (gC/m2) live stem C transfer
   real, pointer :: deadstemc_xfer(:)     ! (gC/m2) dead stem C transfer
   real, pointer :: livecrootc_xfer(:)    ! (gC/m2) live coarse root C transfer
   real, pointer :: deadcrootc_xfer(:)    ! (gC/m2) dead coarse root C transfer
   real, pointer :: gresp_xfer(:)         ! (gC/m2) growth respiration transfer
   real, pointer :: leafn(:)              ! (gN/m2) leaf N
   real, pointer :: frootn(:)             ! (gN/m2) fine root N
   real, pointer :: livestemn(:)          ! (gN/m2) live stem N
   real, pointer :: deadstemn(:)          ! (gN/m2) dead stem N
   real, pointer :: livecrootn(:)         ! (gN/m2) live coarse root N
   real, pointer :: deadcrootn(:)         ! (gN/m2) dead coarse root N
   real, pointer :: retransn(:)           ! (gN/m2) plant pool of retranslocated N
   real, pointer :: leafn_storage(:)      ! (gN/m2) leaf N storage
   real, pointer :: frootn_storage(:)     ! (gN/m2) fine root N storage
   real, pointer :: livestemn_storage(:)  ! (gN/m2) live stem N storage
   real, pointer :: deadstemn_storage(:)  ! (gN/m2) dead stem N storage
   real, pointer :: livecrootn_storage(:) ! (gN/m2) live coarse root N storage
   real, pointer :: deadcrootn_storage(:) ! (gN/m2) dead coarse root N storage
   real, pointer :: leafn_xfer(:)         ! (gN/m2) leaf N transfer
   real, pointer :: frootn_xfer(:)        ! (gN/m2) fine root N transfer
   real, pointer :: livestemn_xfer(:)     ! (gN/m2) live stem N transfer
   real, pointer :: deadstemn_xfer(:)     ! (gN/m2) dead stem N transfer
   real, pointer :: livecrootn_xfer(:)    ! (gN/m2) live coarse root N transfer
   real, pointer :: deadcrootn_xfer(:)    ! (gN/m2) dead coarse root N transfer
!
! local pointers to implicit in/out arrays
!
! local pointers to implicit out arrays
   real, pointer :: hrv_leafc_to_litter(:)
   real, pointer :: hrv_frootc_to_litter(:)
   real, pointer :: hrv_livestemc_to_litter(:)
   real, pointer :: hrv_deadstemc_to_prod10c(:)
   real, pointer :: hrv_deadstemc_to_prod100c(:)
   real, pointer :: hrv_livecrootc_to_litter(:)
   real, pointer :: hrv_deadcrootc_to_litter(:)
   real, pointer :: hrv_xsmrpool_to_atm(:)
   real, pointer :: hrv_leafc_storage_to_litter(:)
   real, pointer :: hrv_frootc_storage_to_litter(:)
   real, pointer :: hrv_livestemc_storage_to_litter(:)
   real, pointer :: hrv_deadstemc_storage_to_litter(:)
   real, pointer :: hrv_livecrootc_storage_to_litter(:)
   real, pointer :: hrv_deadcrootc_storage_to_litter(:)
   real, pointer :: hrv_gresp_storage_to_litter(:)
   real, pointer :: hrv_leafc_xfer_to_litter(:)
   real, pointer :: hrv_frootc_xfer_to_litter(:)
   real, pointer :: hrv_livestemc_xfer_to_litter(:)
   real, pointer :: hrv_deadstemc_xfer_to_litter(:)
   real, pointer :: hrv_livecrootc_xfer_to_litter(:)
   real, pointer :: hrv_deadcrootc_xfer_to_litter(:)
   real, pointer :: hrv_gresp_xfer_to_litter(:)
   real, pointer :: hrv_leafn_to_litter(:)
   real, pointer :: hrv_frootn_to_litter(:)
   real, pointer :: hrv_livestemn_to_litter(:)
   real, pointer :: hrv_deadstemn_to_prod10n(:)
   real, pointer :: hrv_deadstemn_to_prod100n(:)
   real, pointer :: hrv_livecrootn_to_litter(:)
   real, pointer :: hrv_deadcrootn_to_litter(:)
   real, pointer :: hrv_retransn_to_litter(:)
   real, pointer :: hrv_leafn_storage_to_litter(:)
   real, pointer :: hrv_frootn_storage_to_litter(:)
   real, pointer :: hrv_livestemn_storage_to_litter(:)
   real, pointer :: hrv_deadstemn_storage_to_litter(:)
   real, pointer :: hrv_livecrootn_storage_to_litter(:)
   real, pointer :: hrv_deadcrootn_storage_to_litter(:)
   real, pointer :: hrv_leafn_xfer_to_litter(:)
   real, pointer :: hrv_frootn_xfer_to_litter(:)
   real, pointer :: hrv_livestemn_xfer_to_litter(:)
   real, pointer :: hrv_deadstemn_xfer_to_litter(:)
   real, pointer :: hrv_livecrootn_xfer_to_litter(:)
   real, pointer :: hrv_deadcrootn_xfer_to_litter(:)
!
! !OTHER LOCAL VARIABLES:
   integer :: p                         ! pft index
   integer :: g                         ! gridcell index
   integer :: fp                        ! pft filter index
   real:: am                        ! rate for fractional harvest mortality (1/yr)
   real:: m                         ! rate for fractional harvest mortality (1/s)
   real :: pprod10(1:8)   ! proportion of deadstem to 10-yr product pool  (for tree pfts - 1 through 8)
!EOP
!-----------------------------------------------------------------------

   ! assign local pointers to pft-level arrays
   pgridcell                      => clm3%g%l%c%p%gridcell
   
   ivt                            => clm3%g%l%c%p%itype
   leafc                          => clm3%g%l%c%p%pcs%leafc
   frootc                         => clm3%g%l%c%p%pcs%frootc
   livestemc                      => clm3%g%l%c%p%pcs%livestemc
   deadstemc                      => clm3%g%l%c%p%pcs%deadstemc
   livecrootc                     => clm3%g%l%c%p%pcs%livecrootc
   deadcrootc                     => clm3%g%l%c%p%pcs%deadcrootc
   xsmrpool                       => clm3%g%l%c%p%pcs%xsmrpool
   leafc_storage                  => clm3%g%l%c%p%pcs%leafc_storage
   frootc_storage                 => clm3%g%l%c%p%pcs%frootc_storage
   livestemc_storage              => clm3%g%l%c%p%pcs%livestemc_storage
   deadstemc_storage              => clm3%g%l%c%p%pcs%deadstemc_storage
   livecrootc_storage             => clm3%g%l%c%p%pcs%livecrootc_storage
   deadcrootc_storage             => clm3%g%l%c%p%pcs%deadcrootc_storage
   gresp_storage                  => clm3%g%l%c%p%pcs%gresp_storage
   leafc_xfer                     => clm3%g%l%c%p%pcs%leafc_xfer
   frootc_xfer                    => clm3%g%l%c%p%pcs%frootc_xfer
   livestemc_xfer                 => clm3%g%l%c%p%pcs%livestemc_xfer
   deadstemc_xfer                 => clm3%g%l%c%p%pcs%deadstemc_xfer
   livecrootc_xfer                => clm3%g%l%c%p%pcs%livecrootc_xfer
   deadcrootc_xfer                => clm3%g%l%c%p%pcs%deadcrootc_xfer
   gresp_xfer                     => clm3%g%l%c%p%pcs%gresp_xfer
   leafn                          => clm3%g%l%c%p%pns%leafn
   frootn                         => clm3%g%l%c%p%pns%frootn
   livestemn                      => clm3%g%l%c%p%pns%livestemn
   deadstemn                      => clm3%g%l%c%p%pns%deadstemn
   livecrootn                     => clm3%g%l%c%p%pns%livecrootn
   deadcrootn                     => clm3%g%l%c%p%pns%deadcrootn
   retransn                       => clm3%g%l%c%p%pns%retransn
   leafn_storage                  => clm3%g%l%c%p%pns%leafn_storage
   frootn_storage                 => clm3%g%l%c%p%pns%frootn_storage
   livestemn_storage              => clm3%g%l%c%p%pns%livestemn_storage
   deadstemn_storage              => clm3%g%l%c%p%pns%deadstemn_storage
   livecrootn_storage             => clm3%g%l%c%p%pns%livecrootn_storage
   deadcrootn_storage             => clm3%g%l%c%p%pns%deadcrootn_storage
   leafn_xfer                     => clm3%g%l%c%p%pns%leafn_xfer
   frootn_xfer                    => clm3%g%l%c%p%pns%frootn_xfer
   livestemn_xfer                 => clm3%g%l%c%p%pns%livestemn_xfer
   deadstemn_xfer                 => clm3%g%l%c%p%pns%deadstemn_xfer
   livecrootn_xfer                => clm3%g%l%c%p%pns%livecrootn_xfer
   deadcrootn_xfer                => clm3%g%l%c%p%pns%deadcrootn_xfer
   hrv_leafc_to_litter              => clm3%g%l%c%p%pcf%hrv_leafc_to_litter
   hrv_frootc_to_litter             => clm3%g%l%c%p%pcf%hrv_frootc_to_litter
   hrv_livestemc_to_litter          => clm3%g%l%c%p%pcf%hrv_livestemc_to_litter
   hrv_deadstemc_to_prod10c         => clm3%g%l%c%p%pcf%hrv_deadstemc_to_prod10c
   hrv_deadstemc_to_prod100c        => clm3%g%l%c%p%pcf%hrv_deadstemc_to_prod100c
   hrv_livecrootc_to_litter         => clm3%g%l%c%p%pcf%hrv_livecrootc_to_litter
   hrv_deadcrootc_to_litter         => clm3%g%l%c%p%pcf%hrv_deadcrootc_to_litter
   hrv_xsmrpool_to_atm              => clm3%g%l%c%p%pcf%hrv_xsmrpool_to_atm
   hrv_leafc_storage_to_litter      => clm3%g%l%c%p%pcf%hrv_leafc_storage_to_litter
   hrv_frootc_storage_to_litter     => clm3%g%l%c%p%pcf%hrv_frootc_storage_to_litter
   hrv_livestemc_storage_to_litter  => clm3%g%l%c%p%pcf%hrv_livestemc_storage_to_litter
   hrv_deadstemc_storage_to_litter  => clm3%g%l%c%p%pcf%hrv_deadstemc_storage_to_litter
   hrv_livecrootc_storage_to_litter => clm3%g%l%c%p%pcf%hrv_livecrootc_storage_to_litter
   hrv_deadcrootc_storage_to_litter => clm3%g%l%c%p%pcf%hrv_deadcrootc_storage_to_litter
   hrv_gresp_storage_to_litter      => clm3%g%l%c%p%pcf%hrv_gresp_storage_to_litter
   hrv_leafc_xfer_to_litter         => clm3%g%l%c%p%pcf%hrv_leafc_xfer_to_litter
   hrv_frootc_xfer_to_litter        => clm3%g%l%c%p%pcf%hrv_frootc_xfer_to_litter
   hrv_livestemc_xfer_to_litter     => clm3%g%l%c%p%pcf%hrv_livestemc_xfer_to_litter
   hrv_deadstemc_xfer_to_litter     => clm3%g%l%c%p%pcf%hrv_deadstemc_xfer_to_litter
   hrv_livecrootc_xfer_to_litter    => clm3%g%l%c%p%pcf%hrv_livecrootc_xfer_to_litter
   hrv_deadcrootc_xfer_to_litter    => clm3%g%l%c%p%pcf%hrv_deadcrootc_xfer_to_litter
   hrv_gresp_xfer_to_litter         => clm3%g%l%c%p%pcf%hrv_gresp_xfer_to_litter
   hrv_leafn_to_litter              => clm3%g%l%c%p%pnf%hrv_leafn_to_litter
   hrv_frootn_to_litter             => clm3%g%l%c%p%pnf%hrv_frootn_to_litter
   hrv_livestemn_to_litter          => clm3%g%l%c%p%pnf%hrv_livestemn_to_litter
   hrv_deadstemn_to_prod10n         => clm3%g%l%c%p%pnf%hrv_deadstemn_to_prod10n
   hrv_deadstemn_to_prod100n        => clm3%g%l%c%p%pnf%hrv_deadstemn_to_prod100n
   hrv_livecrootn_to_litter         => clm3%g%l%c%p%pnf%hrv_livecrootn_to_litter
   hrv_deadcrootn_to_litter         => clm3%g%l%c%p%pnf%hrv_deadcrootn_to_litter
   hrv_retransn_to_litter           => clm3%g%l%c%p%pnf%hrv_retransn_to_litter
   hrv_leafn_storage_to_litter      => clm3%g%l%c%p%pnf%hrv_leafn_storage_to_litter
   hrv_frootn_storage_to_litter     => clm3%g%l%c%p%pnf%hrv_frootn_storage_to_litter
   hrv_livestemn_storage_to_litter  => clm3%g%l%c%p%pnf%hrv_livestemn_storage_to_litter
   hrv_deadstemn_storage_to_litter  => clm3%g%l%c%p%pnf%hrv_deadstemn_storage_to_litter
   hrv_livecrootn_storage_to_litter => clm3%g%l%c%p%pnf%hrv_livecrootn_storage_to_litter
   hrv_deadcrootn_storage_to_litter => clm3%g%l%c%p%pnf%hrv_deadcrootn_storage_to_litter
   hrv_leafn_xfer_to_litter         => clm3%g%l%c%p%pnf%hrv_leafn_xfer_to_litter
   hrv_frootn_xfer_to_litter        => clm3%g%l%c%p%pnf%hrv_frootn_xfer_to_litter
   hrv_livestemn_xfer_to_litter     => clm3%g%l%c%p%pnf%hrv_livestemn_xfer_to_litter
   hrv_deadstemn_xfer_to_litter     => clm3%g%l%c%p%pnf%hrv_deadstemn_xfer_to_litter
   hrv_livecrootn_xfer_to_litter    => clm3%g%l%c%p%pnf%hrv_livecrootn_xfer_to_litter
   hrv_deadcrootn_xfer_to_litter    => clm3%g%l%c%p%pnf%hrv_deadcrootn_xfer_to_litter

   ! set deadstem proportions to 10-year product pool. 
   ! remainder (1-pprod10) is assumed to go to 100-year product pool
   ! veg type:       1        2        3       4        5       6        7        8      
   pprod10 =    (/0.75, 0.75, 0.75, 1.0, 0.75, 1.0, 0.75, 0.75/)

   ! pft loop
   do fp = 1,num_soilp
      p = filter_soilp(fp)
      g = pgridcell(p)
      
      ! If this is a tree pft, then
      ! get the annual harvest "mortality" rate (am) from harvest array
      ! and convert to rate per second
      if (ivt(p) > noveg .and. ivt(p) < nbrdlf_evr_shrub) then

!        if (do_harvest) then
!           am = harvest(g)
!           m  = am/(365. * 86400.)
!        else
            m = 0. ! gkw force to zero
!        end if   

         ! pft-level harvest carbon fluxes
         ! displayed pools
         hrv_leafc_to_litter(p)               = leafc(p)               * m
         hrv_frootc_to_litter(p)              = frootc(p)              * m
         hrv_livestemc_to_litter(p)           = livestemc(p)           * m
         hrv_deadstemc_to_prod10c(p)          = deadstemc(p)           * m * pprod10(ivt(p))
         hrv_deadstemc_to_prod100c(p)         = deadstemc(p)           * m * (1.0 - pprod10(ivt(p)))
         hrv_livecrootc_to_litter(p)          = livecrootc(p)          * m
         hrv_deadcrootc_to_litter(p)          = deadcrootc(p)          * m
         hrv_xsmrpool_to_atm(p)               = xsmrpool(p)            * m

         ! storage pools
         hrv_leafc_storage_to_litter(p)       = leafc_storage(p)       * m
         hrv_frootc_storage_to_litter(p)      = frootc_storage(p)      * m
         hrv_livestemc_storage_to_litter(p)   = livestemc_storage(p)   * m
         hrv_deadstemc_storage_to_litter(p)   = deadstemc_storage(p)   * m
         hrv_livecrootc_storage_to_litter(p)  = livecrootc_storage(p)  * m
         hrv_deadcrootc_storage_to_litter(p)  = deadcrootc_storage(p)  * m
         hrv_gresp_storage_to_litter(p)       = gresp_storage(p)       * m

         ! transfer pools
         hrv_leafc_xfer_to_litter(p)          = leafc_xfer(p)          * m
         hrv_frootc_xfer_to_litter(p)         = frootc_xfer(p)         * m
         hrv_livestemc_xfer_to_litter(p)      = livestemc_xfer(p)      * m
         hrv_deadstemc_xfer_to_litter(p)      = deadstemc_xfer(p)      * m
         hrv_livecrootc_xfer_to_litter(p)     = livecrootc_xfer(p)     * m
         hrv_deadcrootc_xfer_to_litter(p)     = deadcrootc_xfer(p)     * m
         hrv_gresp_xfer_to_litter(p)          = gresp_xfer(p)          * m

         ! pft-level harvest mortality nitrogen fluxes
         ! displayed pools
         hrv_leafn_to_litter(p)               = leafn(p)               * m
         hrv_frootn_to_litter(p)              = frootn(p)              * m
         hrv_livestemn_to_litter(p)           = livestemn(p)           * m
         hrv_deadstemn_to_prod10n(p)          = deadstemn(p)           * m * pprod10(ivt(p))
         hrv_deadstemn_to_prod100n(p)         = deadstemn(p)           * m * (1.0 - pprod10(ivt(p)))
         hrv_livecrootn_to_litter(p)          = livecrootn(p)          * m
         hrv_deadcrootn_to_litter(p)          = deadcrootn(p)          * m
         hrv_retransn_to_litter(p)            = retransn(p)            * m

         ! storage pools
         hrv_leafn_storage_to_litter(p)       = leafn_storage(p)       * m
         hrv_frootn_storage_to_litter(p)      = frootn_storage(p)      * m
         hrv_livestemn_storage_to_litter(p)   = livestemn_storage(p)   * m
         hrv_deadstemn_storage_to_litter(p)   = deadstemn_storage(p)   * m
         hrv_livecrootn_storage_to_litter(p)  = livecrootn_storage(p)  * m
         hrv_deadcrootn_storage_to_litter(p)  = deadcrootn_storage(p)  * m

         ! transfer pools
         hrv_leafn_xfer_to_litter(p)          = leafn_xfer(p)          * m
         hrv_frootn_xfer_to_litter(p)         = frootn_xfer(p)         * m
         hrv_livestemn_xfer_to_litter(p)      = livestemn_xfer(p)      * m
         hrv_deadstemn_xfer_to_litter(p)      = deadstemn_xfer(p)      * m
         hrv_livecrootn_xfer_to_litter(p)     = livecrootn_xfer(p)     * m
         hrv_deadcrootn_xfer_to_litter(p)     = deadcrootn_xfer(p)     * m
         
      end if  ! end tree block

   end do ! end of pft loop

   ! gather all pft-level litterfall fluxes from harvest to the column
   ! for litter C and N inputs

   call CNHarvestPftToColumn(num_soilc, filter_soilc)

end subroutine CNHarvest
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNHarvestPftToColumn
!
! !INTERFACE:
subroutine CNHarvestPftToColumn (num_soilc, filter_soilc)
!
! !DESCRIPTION:
! called at the end of CNHarvest to gather all pft-level harvest litterfall fluxes
! to the column level and assign them to the three litter pools
!
! !USES:
  use clmtype
  use clm_varpar, only : max_pft_per_col, maxpatch_pft
!
! !ARGUMENTS:
  implicit none
  integer, intent(in) :: num_soilc       ! number of soil columns in filter
  integer, intent(in) :: filter_soilc(:) ! soil column filter
!
! !CALLED FROM:
! subroutine CNphenology
!
! !REVISION HISTORY:
! 9/8/03: Created by Peter Thornton
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in scalars
   integer , pointer :: ivt(:)      ! pft vegetation type
   real, pointer :: wtcol(:)    ! pft weight relative to column (0-1)
   real, pointer :: pwtgcell(:) ! weight of pft relative to corresponding gridcell
   real, pointer :: lf_flab(:)  ! leaf litter labile fraction
   real, pointer :: lf_fcel(:)  ! leaf litter cellulose fraction
   real, pointer :: lf_flig(:)  ! leaf litter lignin fraction
   real, pointer :: fr_flab(:)  ! fine root litter labile fraction
   real, pointer :: fr_fcel(:)  ! fine root litter cellulose fraction
   real, pointer :: fr_flig(:)  ! fine root litter lignin fraction
   integer , pointer :: npfts(:)    ! number of pfts for each column
   integer , pointer :: pfti(:)     ! beginning pft index for each column
   real, pointer :: hrv_leafc_to_litter(:)
   real, pointer :: hrv_frootc_to_litter(:)
   real, pointer :: hrv_livestemc_to_litter(:)
   real, pointer :: phrv_deadstemc_to_prod10c(:)
   real, pointer :: phrv_deadstemc_to_prod100c(:)
   real, pointer :: hrv_livecrootc_to_litter(:)
   real, pointer :: hrv_deadcrootc_to_litter(:)
   real, pointer :: hrv_leafc_storage_to_litter(:)
   real, pointer :: hrv_frootc_storage_to_litter(:)
   real, pointer :: hrv_livestemc_storage_to_litter(:)
   real, pointer :: hrv_deadstemc_storage_to_litter(:)
   real, pointer :: hrv_livecrootc_storage_to_litter(:)
   real, pointer :: hrv_deadcrootc_storage_to_litter(:)
   real, pointer :: hrv_gresp_storage_to_litter(:)
   real, pointer :: hrv_leafc_xfer_to_litter(:)
   real, pointer :: hrv_frootc_xfer_to_litter(:)
   real, pointer :: hrv_livestemc_xfer_to_litter(:)
   real, pointer :: hrv_deadstemc_xfer_to_litter(:)
   real, pointer :: hrv_livecrootc_xfer_to_litter(:)
   real, pointer :: hrv_deadcrootc_xfer_to_litter(:)
   real, pointer :: hrv_gresp_xfer_to_litter(:)
   real, pointer :: hrv_leafn_to_litter(:)
   real, pointer :: hrv_frootn_to_litter(:)
   real, pointer :: hrv_livestemn_to_litter(:)
   real, pointer :: phrv_deadstemn_to_prod10n(:)
   real, pointer :: phrv_deadstemn_to_prod100n(:)
   real, pointer :: hrv_livecrootn_to_litter(:)
   real, pointer :: hrv_deadcrootn_to_litter(:)
   real, pointer :: hrv_retransn_to_litter(:)
   real, pointer :: hrv_leafn_storage_to_litter(:)
   real, pointer :: hrv_frootn_storage_to_litter(:)
   real, pointer :: hrv_livestemn_storage_to_litter(:)
   real, pointer :: hrv_deadstemn_storage_to_litter(:)
   real, pointer :: hrv_livecrootn_storage_to_litter(:)
   real, pointer :: hrv_deadcrootn_storage_to_litter(:)
   real, pointer :: hrv_leafn_xfer_to_litter(:)
   real, pointer :: hrv_frootn_xfer_to_litter(:)
   real, pointer :: hrv_livestemn_xfer_to_litter(:)
   real, pointer :: hrv_deadstemn_xfer_to_litter(:)
   real, pointer :: hrv_livecrootn_xfer_to_litter(:)
   real, pointer :: hrv_deadcrootn_xfer_to_litter(:)
!
! local pointers to implicit in/out arrays
   real, pointer :: hrv_leafc_to_litr1c(:)
   real, pointer :: hrv_leafc_to_litr2c(:)
   real, pointer :: hrv_leafc_to_litr3c(:)
   real, pointer :: hrv_frootc_to_litr1c(:)
   real, pointer :: hrv_frootc_to_litr2c(:)
   real, pointer :: hrv_frootc_to_litr3c(:)
   real, pointer :: hrv_livestemc_to_cwdc(:)
   real, pointer :: chrv_deadstemc_to_prod10c(:)
   real, pointer :: chrv_deadstemc_to_prod100c(:)
   real, pointer :: hrv_livecrootc_to_cwdc(:)
   real, pointer :: hrv_deadcrootc_to_cwdc(:)
   real, pointer :: hrv_leafc_storage_to_litr1c(:)
   real, pointer :: hrv_frootc_storage_to_litr1c(:)
   real, pointer :: hrv_livestemc_storage_to_litr1c(:)
   real, pointer :: hrv_deadstemc_storage_to_litr1c(:)
   real, pointer :: hrv_livecrootc_storage_to_litr1c(:)
   real, pointer :: hrv_deadcrootc_storage_to_litr1c(:)
   real, pointer :: hrv_gresp_storage_to_litr1c(:)
   real, pointer :: hrv_leafc_xfer_to_litr1c(:)
   real, pointer :: hrv_frootc_xfer_to_litr1c(:)
   real, pointer :: hrv_livestemc_xfer_to_litr1c(:)
   real, pointer :: hrv_deadstemc_xfer_to_litr1c(:)
   real, pointer :: hrv_livecrootc_xfer_to_litr1c(:)
   real, pointer :: hrv_deadcrootc_xfer_to_litr1c(:)
   real, pointer :: hrv_gresp_xfer_to_litr1c(:)
   real, pointer :: hrv_leafn_to_litr1n(:)
   real, pointer :: hrv_leafn_to_litr2n(:)
   real, pointer :: hrv_leafn_to_litr3n(:)
   real, pointer :: hrv_frootn_to_litr1n(:)
   real, pointer :: hrv_frootn_to_litr2n(:)
   real, pointer :: hrv_frootn_to_litr3n(:)
   real, pointer :: hrv_livestemn_to_cwdn(:)
   real, pointer :: chrv_deadstemn_to_prod10n(:)
   real, pointer :: chrv_deadstemn_to_prod100n(:)
   real, pointer :: hrv_livecrootn_to_cwdn(:)
   real, pointer :: hrv_deadcrootn_to_cwdn(:)
   real, pointer :: hrv_retransn_to_litr1n(:)
   real, pointer :: hrv_leafn_storage_to_litr1n(:)
   real, pointer :: hrv_frootn_storage_to_litr1n(:)
   real, pointer :: hrv_livestemn_storage_to_litr1n(:)
   real, pointer :: hrv_deadstemn_storage_to_litr1n(:)
   real, pointer :: hrv_livecrootn_storage_to_litr1n(:)
   real, pointer :: hrv_deadcrootn_storage_to_litr1n(:)
   real, pointer :: hrv_leafn_xfer_to_litr1n(:)
   real, pointer :: hrv_frootn_xfer_to_litr1n(:)
   real, pointer :: hrv_livestemn_xfer_to_litr1n(:)
   real, pointer :: hrv_deadstemn_xfer_to_litr1n(:)
   real, pointer :: hrv_livecrootn_xfer_to_litr1n(:)
   real, pointer :: hrv_deadcrootn_xfer_to_litr1n(:)
!
! local pointers to implicit out arrays
!
!
! !OTHER LOCAL VARIABLES:
   integer :: fc,c,pi,p               ! indices
!EOP
!-----------------------------------------------------------------------

   ! assign local pointers
   lf_flab                        => pftcon%lf_flab
   lf_fcel                        => pftcon%lf_fcel
   lf_flig                        => pftcon%lf_flig
   fr_flab                        => pftcon%fr_flab
   fr_fcel                        => pftcon%fr_fcel
   fr_flig                        => pftcon%fr_flig

   ! assign local pointers to column-level arrays
   npfts                          => clm3%g%l%c%npfts
   pfti                           => clm3%g%l%c%pfti
   hrv_leafc_to_litr1c              => clm3%g%l%c%ccf%hrv_leafc_to_litr1c
   hrv_leafc_to_litr2c              => clm3%g%l%c%ccf%hrv_leafc_to_litr2c
   hrv_leafc_to_litr3c              => clm3%g%l%c%ccf%hrv_leafc_to_litr3c
   hrv_frootc_to_litr1c             => clm3%g%l%c%ccf%hrv_frootc_to_litr1c
   hrv_frootc_to_litr2c             => clm3%g%l%c%ccf%hrv_frootc_to_litr2c
   hrv_frootc_to_litr3c             => clm3%g%l%c%ccf%hrv_frootc_to_litr3c
   hrv_livestemc_to_cwdc            => clm3%g%l%c%ccf%hrv_livestemc_to_cwdc
   chrv_deadstemc_to_prod10c        => clm3%g%l%c%ccf%hrv_deadstemc_to_prod10c
   chrv_deadstemc_to_prod100c       => clm3%g%l%c%ccf%hrv_deadstemc_to_prod100c
   hrv_livecrootc_to_cwdc           => clm3%g%l%c%ccf%hrv_livecrootc_to_cwdc
   hrv_deadcrootc_to_cwdc           => clm3%g%l%c%ccf%hrv_deadcrootc_to_cwdc
   hrv_leafc_storage_to_litr1c      => clm3%g%l%c%ccf%hrv_leafc_storage_to_litr1c
   hrv_frootc_storage_to_litr1c     => clm3%g%l%c%ccf%hrv_frootc_storage_to_litr1c
   hrv_livestemc_storage_to_litr1c  => clm3%g%l%c%ccf%hrv_livestemc_storage_to_litr1c
   hrv_deadstemc_storage_to_litr1c  => clm3%g%l%c%ccf%hrv_deadstemc_storage_to_litr1c
   hrv_livecrootc_storage_to_litr1c => clm3%g%l%c%ccf%hrv_livecrootc_storage_to_litr1c
   hrv_deadcrootc_storage_to_litr1c => clm3%g%l%c%ccf%hrv_deadcrootc_storage_to_litr1c
   hrv_gresp_storage_to_litr1c      => clm3%g%l%c%ccf%hrv_gresp_storage_to_litr1c
   hrv_leafc_xfer_to_litr1c         => clm3%g%l%c%ccf%hrv_leafc_xfer_to_litr1c
   hrv_frootc_xfer_to_litr1c        => clm3%g%l%c%ccf%hrv_frootc_xfer_to_litr1c
   hrv_livestemc_xfer_to_litr1c     => clm3%g%l%c%ccf%hrv_livestemc_xfer_to_litr1c
   hrv_deadstemc_xfer_to_litr1c     => clm3%g%l%c%ccf%hrv_deadstemc_xfer_to_litr1c
   hrv_livecrootc_xfer_to_litr1c    => clm3%g%l%c%ccf%hrv_livecrootc_xfer_to_litr1c
   hrv_deadcrootc_xfer_to_litr1c    => clm3%g%l%c%ccf%hrv_deadcrootc_xfer_to_litr1c
   hrv_gresp_xfer_to_litr1c         => clm3%g%l%c%ccf%hrv_gresp_xfer_to_litr1c
   hrv_leafn_to_litr1n              => clm3%g%l%c%cnf%hrv_leafn_to_litr1n
   hrv_leafn_to_litr2n              => clm3%g%l%c%cnf%hrv_leafn_to_litr2n
   hrv_leafn_to_litr3n              => clm3%g%l%c%cnf%hrv_leafn_to_litr3n
   hrv_frootn_to_litr1n             => clm3%g%l%c%cnf%hrv_frootn_to_litr1n
   hrv_frootn_to_litr2n             => clm3%g%l%c%cnf%hrv_frootn_to_litr2n
   hrv_frootn_to_litr3n             => clm3%g%l%c%cnf%hrv_frootn_to_litr3n
   hrv_livestemn_to_cwdn            => clm3%g%l%c%cnf%hrv_livestemn_to_cwdn
   chrv_deadstemn_to_prod10n        => clm3%g%l%c%cnf%hrv_deadstemn_to_prod10n
   chrv_deadstemn_to_prod100n       => clm3%g%l%c%cnf%hrv_deadstemn_to_prod100n
   hrv_livecrootn_to_cwdn           => clm3%g%l%c%cnf%hrv_livecrootn_to_cwdn
   hrv_deadcrootn_to_cwdn           => clm3%g%l%c%cnf%hrv_deadcrootn_to_cwdn
   hrv_retransn_to_litr1n           => clm3%g%l%c%cnf%hrv_retransn_to_litr1n
   hrv_leafn_storage_to_litr1n      => clm3%g%l%c%cnf%hrv_leafn_storage_to_litr1n
   hrv_frootn_storage_to_litr1n     => clm3%g%l%c%cnf%hrv_frootn_storage_to_litr1n
   hrv_livestemn_storage_to_litr1n  => clm3%g%l%c%cnf%hrv_livestemn_storage_to_litr1n
   hrv_deadstemn_storage_to_litr1n  => clm3%g%l%c%cnf%hrv_deadstemn_storage_to_litr1n
   hrv_livecrootn_storage_to_litr1n => clm3%g%l%c%cnf%hrv_livecrootn_storage_to_litr1n
   hrv_deadcrootn_storage_to_litr1n => clm3%g%l%c%cnf%hrv_deadcrootn_storage_to_litr1n
   hrv_leafn_xfer_to_litr1n         => clm3%g%l%c%cnf%hrv_leafn_xfer_to_litr1n
   hrv_frootn_xfer_to_litr1n        => clm3%g%l%c%cnf%hrv_frootn_xfer_to_litr1n
   hrv_livestemn_xfer_to_litr1n     => clm3%g%l%c%cnf%hrv_livestemn_xfer_to_litr1n
   hrv_deadstemn_xfer_to_litr1n     => clm3%g%l%c%cnf%hrv_deadstemn_xfer_to_litr1n
   hrv_livecrootn_xfer_to_litr1n    => clm3%g%l%c%cnf%hrv_livecrootn_xfer_to_litr1n
   hrv_deadcrootn_xfer_to_litr1n    => clm3%g%l%c%cnf%hrv_deadcrootn_xfer_to_litr1n

   ! assign local pointers to pft-level arrays
   ivt                            => clm3%g%l%c%p%itype
   wtcol                          => clm3%g%l%c%p%wtcol
   pwtgcell                       => clm3%g%l%c%p%wtgcell  
   hrv_leafc_to_litter              => clm3%g%l%c%p%pcf%hrv_leafc_to_litter
   hrv_frootc_to_litter             => clm3%g%l%c%p%pcf%hrv_frootc_to_litter
   hrv_livestemc_to_litter          => clm3%g%l%c%p%pcf%hrv_livestemc_to_litter
   phrv_deadstemc_to_prod10c        => clm3%g%l%c%p%pcf%hrv_deadstemc_to_prod10c
   phrv_deadstemc_to_prod100c       => clm3%g%l%c%p%pcf%hrv_deadstemc_to_prod100c
   hrv_livecrootc_to_litter         => clm3%g%l%c%p%pcf%hrv_livecrootc_to_litter
   hrv_deadcrootc_to_litter         => clm3%g%l%c%p%pcf%hrv_deadcrootc_to_litter
   hrv_leafc_storage_to_litter      => clm3%g%l%c%p%pcf%hrv_leafc_storage_to_litter
   hrv_frootc_storage_to_litter     => clm3%g%l%c%p%pcf%hrv_frootc_storage_to_litter
   hrv_livestemc_storage_to_litter  => clm3%g%l%c%p%pcf%hrv_livestemc_storage_to_litter
   hrv_deadstemc_storage_to_litter  => clm3%g%l%c%p%pcf%hrv_deadstemc_storage_to_litter
   hrv_livecrootc_storage_to_litter => clm3%g%l%c%p%pcf%hrv_livecrootc_storage_to_litter
   hrv_deadcrootc_storage_to_litter => clm3%g%l%c%p%pcf%hrv_deadcrootc_storage_to_litter
   hrv_gresp_storage_to_litter      => clm3%g%l%c%p%pcf%hrv_gresp_storage_to_litter
   hrv_leafc_xfer_to_litter         => clm3%g%l%c%p%pcf%hrv_leafc_xfer_to_litter
   hrv_frootc_xfer_to_litter        => clm3%g%l%c%p%pcf%hrv_frootc_xfer_to_litter
   hrv_livestemc_xfer_to_litter     => clm3%g%l%c%p%pcf%hrv_livestemc_xfer_to_litter
   hrv_deadstemc_xfer_to_litter     => clm3%g%l%c%p%pcf%hrv_deadstemc_xfer_to_litter
   hrv_livecrootc_xfer_to_litter    => clm3%g%l%c%p%pcf%hrv_livecrootc_xfer_to_litter
   hrv_deadcrootc_xfer_to_litter    => clm3%g%l%c%p%pcf%hrv_deadcrootc_xfer_to_litter
   hrv_gresp_xfer_to_litter         => clm3%g%l%c%p%pcf%hrv_gresp_xfer_to_litter
   hrv_leafn_to_litter              => clm3%g%l%c%p%pnf%hrv_leafn_to_litter
   hrv_frootn_to_litter             => clm3%g%l%c%p%pnf%hrv_frootn_to_litter
   hrv_livestemn_to_litter          => clm3%g%l%c%p%pnf%hrv_livestemn_to_litter
   phrv_deadstemn_to_prod10n        => clm3%g%l%c%p%pnf%hrv_deadstemn_to_prod10n
   phrv_deadstemn_to_prod100n       => clm3%g%l%c%p%pnf%hrv_deadstemn_to_prod100n
   hrv_livecrootn_to_litter         => clm3%g%l%c%p%pnf%hrv_livecrootn_to_litter
   hrv_deadcrootn_to_litter         => clm3%g%l%c%p%pnf%hrv_deadcrootn_to_litter
   hrv_retransn_to_litter           => clm3%g%l%c%p%pnf%hrv_retransn_to_litter
   hrv_leafn_storage_to_litter      => clm3%g%l%c%p%pnf%hrv_leafn_storage_to_litter
   hrv_frootn_storage_to_litter     => clm3%g%l%c%p%pnf%hrv_frootn_storage_to_litter
   hrv_livestemn_storage_to_litter  => clm3%g%l%c%p%pnf%hrv_livestemn_storage_to_litter
   hrv_deadstemn_storage_to_litter  => clm3%g%l%c%p%pnf%hrv_deadstemn_storage_to_litter
   hrv_livecrootn_storage_to_litter => clm3%g%l%c%p%pnf%hrv_livecrootn_storage_to_litter
   hrv_deadcrootn_storage_to_litter => clm3%g%l%c%p%pnf%hrv_deadcrootn_storage_to_litter
   hrv_leafn_xfer_to_litter         => clm3%g%l%c%p%pnf%hrv_leafn_xfer_to_litter
   hrv_frootn_xfer_to_litter        => clm3%g%l%c%p%pnf%hrv_frootn_xfer_to_litter
   hrv_livestemn_xfer_to_litter     => clm3%g%l%c%p%pnf%hrv_livestemn_xfer_to_litter
   hrv_deadstemn_xfer_to_litter     => clm3%g%l%c%p%pnf%hrv_deadstemn_xfer_to_litter
   hrv_livecrootn_xfer_to_litter    => clm3%g%l%c%p%pnf%hrv_livecrootn_xfer_to_litter
   hrv_deadcrootn_xfer_to_litter    => clm3%g%l%c%p%pnf%hrv_deadcrootn_xfer_to_litter

   do pi = 1,maxpatch_pft
      do fc = 1,num_soilc
         c = filter_soilc(fc)

         if (pi <=  npfts(c)) then
            p = pfti(c) + pi - 1

            if (pwtgcell(p)>0.) then

               ! leaf harvest mortality carbon fluxes
               hrv_leafc_to_litr1c(c) = hrv_leafc_to_litr1c(c) + &
                  hrv_leafc_to_litter(p) * lf_flab(ivt(p)) * wtcol(p)
               hrv_leafc_to_litr2c(c) = hrv_leafc_to_litr2c(c) + &
                  hrv_leafc_to_litter(p) * lf_fcel(ivt(p)) * wtcol(p)
               hrv_leafc_to_litr3c(c) = hrv_leafc_to_litr3c(c) + &
                  hrv_leafc_to_litter(p) * lf_flig(ivt(p)) * wtcol(p)

               ! fine root harvest mortality carbon fluxes
               hrv_frootc_to_litr1c(c) = hrv_frootc_to_litr1c(c) + &
                  hrv_frootc_to_litter(p) * fr_flab(ivt(p)) * wtcol(p)
               hrv_frootc_to_litr2c(c) = hrv_frootc_to_litr2c(c) + &
                  hrv_frootc_to_litter(p) * fr_fcel(ivt(p)) * wtcol(p)
               hrv_frootc_to_litr3c(c) = hrv_frootc_to_litr3c(c) + &
                  hrv_frootc_to_litter(p) * fr_flig(ivt(p)) * wtcol(p)

               ! wood harvest mortality carbon fluxes
               hrv_livestemc_to_cwdc(c)  = hrv_livestemc_to_cwdc(c)  + &
                  hrv_livestemc_to_litter(p)  * wtcol(p)
               chrv_deadstemc_to_prod10c(c)  = chrv_deadstemc_to_prod10c(c)  + &
                  phrv_deadstemc_to_prod10c(p)  * wtcol(p)
               chrv_deadstemc_to_prod100c(c)  = chrv_deadstemc_to_prod100c(c)  + &
                  phrv_deadstemc_to_prod100c(p)  * wtcol(p)
               hrv_livecrootc_to_cwdc(c) = hrv_livecrootc_to_cwdc(c) + &
                  hrv_livecrootc_to_litter(p) * wtcol(p)
               hrv_deadcrootc_to_cwdc(c) = hrv_deadcrootc_to_cwdc(c) + &
                  hrv_deadcrootc_to_litter(p) * wtcol(p)

               ! storage harvest mortality carbon fluxes
               hrv_leafc_storage_to_litr1c(c)      = hrv_leafc_storage_to_litr1c(c)      + &
                  hrv_leafc_storage_to_litter(p)      * wtcol(p)
               hrv_frootc_storage_to_litr1c(c)     = hrv_frootc_storage_to_litr1c(c)     + &
                  hrv_frootc_storage_to_litter(p)     * wtcol(p)
               hrv_livestemc_storage_to_litr1c(c)  = hrv_livestemc_storage_to_litr1c(c)  + &
                  hrv_livestemc_storage_to_litter(p)  * wtcol(p)
               hrv_deadstemc_storage_to_litr1c(c)  = hrv_deadstemc_storage_to_litr1c(c)  + &
                  hrv_deadstemc_storage_to_litter(p)  * wtcol(p)
               hrv_livecrootc_storage_to_litr1c(c) = hrv_livecrootc_storage_to_litr1c(c) + &
                  hrv_livecrootc_storage_to_litter(p) * wtcol(p)
               hrv_deadcrootc_storage_to_litr1c(c) = hrv_deadcrootc_storage_to_litr1c(c) + &
                  hrv_deadcrootc_storage_to_litter(p) * wtcol(p)
               hrv_gresp_storage_to_litr1c(c)      = hrv_gresp_storage_to_litr1c(c)      + &
                  hrv_gresp_storage_to_litter(p)      * wtcol(p)

               ! transfer harvest mortality carbon fluxes
               hrv_leafc_xfer_to_litr1c(c)      = hrv_leafc_xfer_to_litr1c(c)      + &
                  hrv_leafc_xfer_to_litter(p)      * wtcol(p)
               hrv_frootc_xfer_to_litr1c(c)     = hrv_frootc_xfer_to_litr1c(c)     + &
                  hrv_frootc_xfer_to_litter(p)     * wtcol(p)
               hrv_livestemc_xfer_to_litr1c(c)  = hrv_livestemc_xfer_to_litr1c(c)  + &
                  hrv_livestemc_xfer_to_litter(p)  * wtcol(p)
               hrv_deadstemc_xfer_to_litr1c(c)  = hrv_deadstemc_xfer_to_litr1c(c)  + &
                  hrv_deadstemc_xfer_to_litter(p)  * wtcol(p)
               hrv_livecrootc_xfer_to_litr1c(c) = hrv_livecrootc_xfer_to_litr1c(c) + &
                  hrv_livecrootc_xfer_to_litter(p) * wtcol(p)
               hrv_deadcrootc_xfer_to_litr1c(c) = hrv_deadcrootc_xfer_to_litr1c(c) + &
                  hrv_deadcrootc_xfer_to_litter(p) * wtcol(p)
               hrv_gresp_xfer_to_litr1c(c)      = hrv_gresp_xfer_to_litr1c(c)      + &
                  hrv_gresp_xfer_to_litter(p)      * wtcol(p)

               ! leaf harvest mortality nitrogen fluxes
               hrv_leafn_to_litr1n(c) = hrv_leafn_to_litr1n(c) + &
                  hrv_leafn_to_litter(p) * lf_flab(ivt(p)) * wtcol(p)
               hrv_leafn_to_litr2n(c) = hrv_leafn_to_litr2n(c) + &
                  hrv_leafn_to_litter(p) * lf_fcel(ivt(p)) * wtcol(p)
               hrv_leafn_to_litr3n(c) = hrv_leafn_to_litr3n(c) + &
                  hrv_leafn_to_litter(p) * lf_flig(ivt(p)) * wtcol(p)

               ! fine root litter nitrogen fluxes
               hrv_frootn_to_litr1n(c) = hrv_frootn_to_litr1n(c) + &
                  hrv_frootn_to_litter(p) * fr_flab(ivt(p)) * wtcol(p)
               hrv_frootn_to_litr2n(c) = hrv_frootn_to_litr2n(c) + &
                  hrv_frootn_to_litter(p) * fr_fcel(ivt(p)) * wtcol(p)
               hrv_frootn_to_litr3n(c) = hrv_frootn_to_litr3n(c) + &
                  hrv_frootn_to_litter(p) * fr_flig(ivt(p)) * wtcol(p)

               ! wood harvest mortality nitrogen fluxes
               hrv_livestemn_to_cwdn(c)  = hrv_livestemn_to_cwdn(c)  + &
                  hrv_livestemn_to_litter(p)  * wtcol(p)
               chrv_deadstemn_to_prod10n(c)  = chrv_deadstemn_to_prod10n(c)  + &
                  phrv_deadstemn_to_prod10n(p)  * wtcol(p)
               chrv_deadstemn_to_prod100n(c)  = chrv_deadstemn_to_prod100n(c)  + &
                  phrv_deadstemn_to_prod100n(p)  * wtcol(p)
               hrv_livecrootn_to_cwdn(c) = hrv_livecrootn_to_cwdn(c) + &
                  hrv_livecrootn_to_litter(p) * wtcol(p)
               hrv_deadcrootn_to_cwdn(c) = hrv_deadcrootn_to_cwdn(c) + &
                  hrv_deadcrootn_to_litter(p) * wtcol(p)

               ! retranslocated N pool harvest mortality fluxes
               hrv_retransn_to_litr1n(c) = hrv_retransn_to_litr1n(c) + &
                  hrv_retransn_to_litter(p) * wtcol(p)

               ! storage harvest mortality nitrogen fluxes
               hrv_leafn_storage_to_litr1n(c)      = hrv_leafn_storage_to_litr1n(c)      + &
                  hrv_leafn_storage_to_litter(p)      * wtcol(p)
               hrv_frootn_storage_to_litr1n(c)     = hrv_frootn_storage_to_litr1n(c)     + &
                  hrv_frootn_storage_to_litter(p)     * wtcol(p)
               hrv_livestemn_storage_to_litr1n(c)  = hrv_livestemn_storage_to_litr1n(c)  + &
                  hrv_livestemn_storage_to_litter(p)  * wtcol(p)
               hrv_deadstemn_storage_to_litr1n(c)  = hrv_deadstemn_storage_to_litr1n(c)  + &
                  hrv_deadstemn_storage_to_litter(p)  * wtcol(p)
               hrv_livecrootn_storage_to_litr1n(c) = hrv_livecrootn_storage_to_litr1n(c) + &
                  hrv_livecrootn_storage_to_litter(p) * wtcol(p)
               hrv_deadcrootn_storage_to_litr1n(c) = hrv_deadcrootn_storage_to_litr1n(c) + &
                  hrv_deadcrootn_storage_to_litter(p) * wtcol(p)

               ! transfer harvest mortality nitrogen fluxes
               hrv_leafn_xfer_to_litr1n(c)      = hrv_leafn_xfer_to_litr1n(c)      + &
                  hrv_leafn_xfer_to_litter(p)      * wtcol(p)
               hrv_frootn_xfer_to_litr1n(c)     = hrv_frootn_xfer_to_litr1n(c)     + &
                  hrv_frootn_xfer_to_litter(p)     * wtcol(p)
               hrv_livestemn_xfer_to_litr1n(c)  = hrv_livestemn_xfer_to_litr1n(c)  + &
                  hrv_livestemn_xfer_to_litter(p)  * wtcol(p)
               hrv_deadstemn_xfer_to_litr1n(c)  = hrv_deadstemn_xfer_to_litr1n(c)  + &
                  hrv_deadstemn_xfer_to_litter(p)  * wtcol(p)
               hrv_livecrootn_xfer_to_litr1n(c) = hrv_livecrootn_xfer_to_litr1n(c) + &
                  hrv_livecrootn_xfer_to_litter(p) * wtcol(p)
               hrv_deadcrootn_xfer_to_litr1n(c) = hrv_deadcrootn_xfer_to_litr1n(c) + &
                  hrv_deadcrootn_xfer_to_litter(p) * wtcol(p)

            end if
         end if

      end do

   end do

end subroutine CNHarvestPftToColumn
!-----------------------------------------------------------------------

end module CNHarvestMod
