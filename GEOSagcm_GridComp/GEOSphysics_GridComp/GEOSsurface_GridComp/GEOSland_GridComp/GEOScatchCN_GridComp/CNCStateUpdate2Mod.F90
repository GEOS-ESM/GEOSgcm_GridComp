module CNCStateUpdate2Mod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: CStateUpdate2Mod
!
! !DESCRIPTION:
! Module for carbon state variable update, mortality fluxes.
!
! !USES:
    implicit none
    save
    private
! !PUBLIC MEMBER FUNCTIONS:
    public:: CStateUpdate2
    public:: CStateUpdate2h
!
! !REVISION HISTORY:
! 4/23/2004: Created by Peter Thornton
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CStateUpdate2
!
! !INTERFACE:
subroutine CStateUpdate2(num_soilc, filter_soilc, num_soilp, filter_soilp)
!
! !DESCRIPTION:
! On the radiation time step, update all the prognostic carbon state
! variables affected by gap-phase mortality fluxes
!
! !USES:
   use clmtype
   use clm_time_manager, only: get_step_size
!
! !ARGUMENTS:
   implicit none
   integer, intent(in) :: num_soilc       ! number of soil columns in filter
   integer, intent(in) :: filter_soilc(:) ! filter for soil columns
   integer, intent(in) :: num_soilp       ! number of soil pfts in filter
   integer, intent(in) :: filter_soilp(:) ! filter for soil pfts
!
! !CALLED FROM:
! subroutine CNEcosystemDyn
!
! !REVISION HISTORY:
! 3/29/04: Created by Peter Thornton
!
! !LOCAL VARIABLES:
! local pointers to implicit in arrays
   real, pointer :: m_deadcrootc_storage_to_litr1c(:)
   real, pointer :: m_deadcrootc_to_cwdc(:)
   real, pointer :: m_deadcrootc_xfer_to_litr1c(:)
   real, pointer :: m_deadstemc_storage_to_litr1c(:)
   real, pointer :: m_deadstemc_to_cwdc(:)
   real, pointer :: m_deadstemc_xfer_to_litr1c(:)
   real, pointer :: m_frootc_storage_to_litr1c(:)
   real, pointer :: m_frootc_to_litr1c(:)
   real, pointer :: m_frootc_to_litr2c(:)
   real, pointer :: m_frootc_to_litr3c(:)
   real, pointer :: m_frootc_xfer_to_litr1c(:)
   real, pointer :: m_gresp_storage_to_litr1c(:)
   real, pointer :: m_gresp_xfer_to_litr1c(:)
   real, pointer :: m_leafc_storage_to_litr1c(:)
   real, pointer :: m_leafc_to_litr1c(:)
   real, pointer :: m_leafc_to_litr2c(:)
   real, pointer :: m_leafc_to_litr3c(:)
   real, pointer :: m_leafc_xfer_to_litr1c(:)
   real, pointer :: m_livecrootc_storage_to_litr1c(:)
   real, pointer :: m_livecrootc_to_cwdc(:)
   real, pointer :: m_livecrootc_xfer_to_litr1c(:)
   real, pointer :: m_livestemc_storage_to_litr1c(:)
   real, pointer :: m_livestemc_to_cwdc(:)
   real, pointer :: m_livestemc_xfer_to_litr1c(:)
   real, pointer :: m_deadcrootc_storage_to_litter(:)
   real, pointer :: m_deadcrootc_to_litter(:)
   real, pointer :: m_deadcrootc_xfer_to_litter(:)
   real, pointer :: m_deadstemc_storage_to_litter(:)
   real, pointer :: m_deadstemc_to_litter(:)
   real, pointer :: m_deadstemc_xfer_to_litter(:)
   real, pointer :: m_frootc_storage_to_litter(:)
   real, pointer :: m_frootc_to_litter(:)
   real, pointer :: m_frootc_xfer_to_litter(:)
   real, pointer :: m_gresp_storage_to_litter(:)
   real, pointer :: m_gresp_xfer_to_litter(:)
   real, pointer :: m_leafc_storage_to_litter(:)
   real, pointer :: m_leafc_to_litter(:)
   real, pointer :: m_leafc_xfer_to_litter(:)
   real, pointer :: m_livecrootc_storage_to_litter(:)
   real, pointer :: m_livecrootc_to_litter(:)
   real, pointer :: m_livecrootc_xfer_to_litter(:)
   real, pointer :: m_livestemc_storage_to_litter(:)
   real, pointer :: m_livestemc_to_litter(:)
   real, pointer :: m_livestemc_xfer_to_litter(:)
!
! local pointers to implicit in/out arrays
   real, pointer :: cwdc(:)               ! (gC/m2) coarse woody debris C
   real, pointer :: litr1c(:)             ! (gC/m2) litter labile C
   real, pointer :: litr2c(:)             ! (gC/m2) litter cellulose C
   real, pointer :: litr3c(:)             ! (gC/m2) litter lignin C
   real, pointer :: deadcrootc(:)         ! (gC/m2) dead coarse root C
   real, pointer :: deadcrootc_storage(:) ! (gC/m2) dead coarse root C storage
   real, pointer :: deadcrootc_xfer(:)    !(gC/m2) dead coarse root C transfer
   real, pointer :: deadstemc(:)          ! (gC/m2) dead stem C
   real, pointer :: deadstemc_storage(:)  ! (gC/m2) dead stem C storage
   real, pointer :: deadstemc_xfer(:)     ! (gC/m2) dead stem C transfer
   real, pointer :: frootc(:)             ! (gC/m2) fine root C
   real, pointer :: frootc_storage(:)     ! (gC/m2) fine root C storage
   real, pointer :: frootc_xfer(:)        ! (gC/m2) fine root C transfer
   real, pointer :: gresp_storage(:)      ! (gC/m2) growth respiration storage
   real, pointer :: gresp_xfer(:)         ! (gC/m2) growth respiration transfer
   real, pointer :: leafc(:)              ! (gC/m2) leaf C
   real, pointer :: leafc_storage(:)      ! (gC/m2) leaf C storage
   real, pointer :: leafc_xfer(:)         ! (gC/m2) leaf C transfer
   real, pointer :: livecrootc(:)         ! (gC/m2) live coarse root C
   real, pointer :: livecrootc_storage(:) ! (gC/m2) live coarse root C storage
   real, pointer :: livecrootc_xfer(:)    !(gC/m2) live coarse root C transfer
   real, pointer :: livestemc(:)          ! (gC/m2) live stem C
   real, pointer :: livestemc_storage(:)  ! (gC/m2) live stem C storage
   real, pointer :: livestemc_xfer(:)     ! (gC/m2) live stem C transfer
!
!
! local pointers to implicit out arrays
!
!
! !OTHER LOCAL VARIABLES:
   integer :: c,p      ! indices
   integer :: fp,fc    ! lake filter indices
   real:: dt       ! radiation time step (seconds)
!
!EOP
!-----------------------------------------------------------------------
    ! assign local pointers at the column level
    m_deadcrootc_storage_to_litr1c => clm3%g%l%c%ccf%m_deadcrootc_storage_to_litr1c
    m_deadcrootc_to_cwdc           => clm3%g%l%c%ccf%m_deadcrootc_to_cwdc
    m_deadcrootc_xfer_to_litr1c    => clm3%g%l%c%ccf%m_deadcrootc_xfer_to_litr1c
    m_deadstemc_storage_to_litr1c  => clm3%g%l%c%ccf%m_deadstemc_storage_to_litr1c
    m_deadstemc_to_cwdc            => clm3%g%l%c%ccf%m_deadstemc_to_cwdc
    m_deadstemc_xfer_to_litr1c     => clm3%g%l%c%ccf%m_deadstemc_xfer_to_litr1c
    m_frootc_storage_to_litr1c     => clm3%g%l%c%ccf%m_frootc_storage_to_litr1c
    m_frootc_to_litr1c             => clm3%g%l%c%ccf%m_frootc_to_litr1c
    m_frootc_to_litr2c             => clm3%g%l%c%ccf%m_frootc_to_litr2c
    m_frootc_to_litr3c             => clm3%g%l%c%ccf%m_frootc_to_litr3c
    m_frootc_xfer_to_litr1c        => clm3%g%l%c%ccf%m_frootc_xfer_to_litr1c
    m_gresp_storage_to_litr1c      => clm3%g%l%c%ccf%m_gresp_storage_to_litr1c
    m_gresp_xfer_to_litr1c         => clm3%g%l%c%ccf%m_gresp_xfer_to_litr1c
    m_leafc_storage_to_litr1c      => clm3%g%l%c%ccf%m_leafc_storage_to_litr1c
    m_leafc_to_litr1c              => clm3%g%l%c%ccf%m_leafc_to_litr1c
    m_leafc_to_litr2c              => clm3%g%l%c%ccf%m_leafc_to_litr2c
    m_leafc_to_litr3c              => clm3%g%l%c%ccf%m_leafc_to_litr3c
    m_leafc_xfer_to_litr1c         => clm3%g%l%c%ccf%m_leafc_xfer_to_litr1c
    m_livecrootc_storage_to_litr1c => clm3%g%l%c%ccf%m_livecrootc_storage_to_litr1c
    m_livecrootc_to_cwdc           => clm3%g%l%c%ccf%m_livecrootc_to_cwdc
    m_livecrootc_xfer_to_litr1c    => clm3%g%l%c%ccf%m_livecrootc_xfer_to_litr1c
    m_livestemc_storage_to_litr1c  => clm3%g%l%c%ccf%m_livestemc_storage_to_litr1c
    m_livestemc_to_cwdc            => clm3%g%l%c%ccf%m_livestemc_to_cwdc
    m_livestemc_xfer_to_litr1c     => clm3%g%l%c%ccf%m_livestemc_xfer_to_litr1c
    cwdc                           => clm3%g%l%c%ccs%cwdc
    litr1c                         => clm3%g%l%c%ccs%litr1c
    litr2c                         => clm3%g%l%c%ccs%litr2c
    litr3c                         => clm3%g%l%c%ccs%litr3c

    ! assign local pointers at the pft level
    m_deadcrootc_storage_to_litter => clm3%g%l%c%p%pcf%m_deadcrootc_storage_to_litter
    m_deadcrootc_to_litter         => clm3%g%l%c%p%pcf%m_deadcrootc_to_litter
    m_deadcrootc_xfer_to_litter    => clm3%g%l%c%p%pcf%m_deadcrootc_xfer_to_litter
    m_deadstemc_storage_to_litter  => clm3%g%l%c%p%pcf%m_deadstemc_storage_to_litter
    m_deadstemc_to_litter          => clm3%g%l%c%p%pcf%m_deadstemc_to_litter
    m_deadstemc_xfer_to_litter     => clm3%g%l%c%p%pcf%m_deadstemc_xfer_to_litter
    m_frootc_storage_to_litter     => clm3%g%l%c%p%pcf%m_frootc_storage_to_litter
    m_frootc_to_litter             => clm3%g%l%c%p%pcf%m_frootc_to_litter
    m_frootc_xfer_to_litter        => clm3%g%l%c%p%pcf%m_frootc_xfer_to_litter
    m_gresp_storage_to_litter      => clm3%g%l%c%p%pcf%m_gresp_storage_to_litter
    m_gresp_xfer_to_litter         => clm3%g%l%c%p%pcf%m_gresp_xfer_to_litter
    m_leafc_storage_to_litter      => clm3%g%l%c%p%pcf%m_leafc_storage_to_litter
    m_leafc_to_litter              => clm3%g%l%c%p%pcf%m_leafc_to_litter
    m_leafc_xfer_to_litter         => clm3%g%l%c%p%pcf%m_leafc_xfer_to_litter
    m_livecrootc_storage_to_litter => clm3%g%l%c%p%pcf%m_livecrootc_storage_to_litter
    m_livecrootc_to_litter         => clm3%g%l%c%p%pcf%m_livecrootc_to_litter
    m_livecrootc_xfer_to_litter    => clm3%g%l%c%p%pcf%m_livecrootc_xfer_to_litter
    m_livestemc_storage_to_litter  => clm3%g%l%c%p%pcf%m_livestemc_storage_to_litter
    m_livestemc_to_litter          => clm3%g%l%c%p%pcf%m_livestemc_to_litter
    m_livestemc_xfer_to_litter     => clm3%g%l%c%p%pcf%m_livestemc_xfer_to_litter
    deadcrootc                     => clm3%g%l%c%p%pcs%deadcrootc
    deadcrootc_storage             => clm3%g%l%c%p%pcs%deadcrootc_storage
    deadcrootc_xfer                => clm3%g%l%c%p%pcs%deadcrootc_xfer
    deadstemc                      => clm3%g%l%c%p%pcs%deadstemc
    deadstemc_storage              => clm3%g%l%c%p%pcs%deadstemc_storage
    deadstemc_xfer                 => clm3%g%l%c%p%pcs%deadstemc_xfer
    frootc                         => clm3%g%l%c%p%pcs%frootc
    frootc_storage                 => clm3%g%l%c%p%pcs%frootc_storage
    frootc_xfer                    => clm3%g%l%c%p%pcs%frootc_xfer
    gresp_storage                  => clm3%g%l%c%p%pcs%gresp_storage
    gresp_xfer                     => clm3%g%l%c%p%pcs%gresp_xfer
    leafc                          => clm3%g%l%c%p%pcs%leafc
    leafc_storage                  => clm3%g%l%c%p%pcs%leafc_storage
    leafc_xfer                     => clm3%g%l%c%p%pcs%leafc_xfer
    livecrootc                     => clm3%g%l%c%p%pcs%livecrootc
    livecrootc_storage             => clm3%g%l%c%p%pcs%livecrootc_storage
    livecrootc_xfer                => clm3%g%l%c%p%pcs%livecrootc_xfer
    livestemc                      => clm3%g%l%c%p%pcs%livestemc
    livestemc_storage              => clm3%g%l%c%p%pcs%livestemc_storage
    livestemc_xfer                 => clm3%g%l%c%p%pcs%livestemc_xfer

    ! set time steps
    dt = real( get_step_size() )

    ! column loop
    do fc = 1,num_soilc
       c = filter_soilc(fc)

       ! column level carbon fluxes from gap-phase mortality

       ! leaf to litter
       litr1c(c) = litr1c(c) + m_leafc_to_litr1c(c) * dt
       litr2c(c) = litr2c(c) + m_leafc_to_litr2c(c) * dt
       litr3c(c) = litr3c(c) + m_leafc_to_litr3c(c) * dt

       ! fine root to litter
       litr1c(c) = litr1c(c) + m_frootc_to_litr1c(c) * dt
       litr2c(c) = litr2c(c) + m_frootc_to_litr2c(c) * dt
       litr3c(c) = litr3c(c) + m_frootc_to_litr3c(c) * dt

       ! wood to CWD
       cwdc(c) = cwdc(c) + m_livestemc_to_cwdc(c)  * dt
       cwdc(c) = cwdc(c) + m_deadstemc_to_cwdc(c)  * dt
       cwdc(c) = cwdc(c) + m_livecrootc_to_cwdc(c) * dt
       cwdc(c) = cwdc(c) + m_deadcrootc_to_cwdc(c) * dt

       ! storage pools to litter
       litr1c(c) = litr1c(c) + m_leafc_storage_to_litr1c(c)      * dt
       litr1c(c) = litr1c(c) + m_frootc_storage_to_litr1c(c)     * dt
       litr1c(c) = litr1c(c) + m_livestemc_storage_to_litr1c(c)  * dt
       litr1c(c) = litr1c(c) + m_deadstemc_storage_to_litr1c(c)  * dt
       litr1c(c) = litr1c(c) + m_livecrootc_storage_to_litr1c(c) * dt
       litr1c(c) = litr1c(c) + m_deadcrootc_storage_to_litr1c(c) * dt
       litr1c(c) = litr1c(c) + m_gresp_storage_to_litr1c(c)      * dt

       ! transfer pools to litter
       litr1c(c) = litr1c(c) + m_leafc_xfer_to_litr1c(c)      * dt
       litr1c(c) = litr1c(c) + m_frootc_xfer_to_litr1c(c)     * dt
       litr1c(c) = litr1c(c) + m_livestemc_xfer_to_litr1c(c)  * dt
       litr1c(c) = litr1c(c) + m_deadstemc_xfer_to_litr1c(c)  * dt
       litr1c(c) = litr1c(c) + m_livecrootc_xfer_to_litr1c(c) * dt
       litr1c(c) = litr1c(c) + m_deadcrootc_xfer_to_litr1c(c) * dt
       litr1c(c) = litr1c(c) + m_gresp_xfer_to_litr1c(c)      * dt

    end do ! end of columns loop

    ! pft loop
    do fp = 1,num_soilp
       p = filter_soilp(fp)

       ! pft-level carbon fluxes from gap-phase mortality
       ! displayed pools
       leafc(p)               = leafc(p)              - m_leafc_to_litter(p)              * dt
       frootc(p)              = frootc(p)             - m_frootc_to_litter(p)             * dt
       livestemc(p)           = livestemc(p)          - m_livestemc_to_litter(p)          * dt
       deadstemc(p)           = deadstemc(p)          - m_deadstemc_to_litter(p)          * dt
       livecrootc(p)          = livecrootc(p)         - m_livecrootc_to_litter(p)         * dt
       deadcrootc(p)          = deadcrootc(p)         - m_deadcrootc_to_litter(p)         * dt

       ! storage pools
       leafc_storage(p)       = leafc_storage(p)      - m_leafc_storage_to_litter(p)      * dt
       frootc_storage(p)      = frootc_storage(p)     - m_frootc_storage_to_litter(p)     * dt
       livestemc_storage(p)   = livestemc_storage(p)  - m_livestemc_storage_to_litter(p)  * dt
       deadstemc_storage(p)   = deadstemc_storage(p)  - m_deadstemc_storage_to_litter(p)  * dt
       livecrootc_storage(p)  = livecrootc_storage(p) - m_livecrootc_storage_to_litter(p) * dt
       deadcrootc_storage(p)  = deadcrootc_storage(p) - m_deadcrootc_storage_to_litter(p) * dt
       gresp_storage(p)       = gresp_storage(p)      - m_gresp_storage_to_litter(p)      * dt

       ! transfer pools
       leafc_xfer(p)          = leafc_xfer(p)         - m_leafc_xfer_to_litter(p)         * dt
       frootc_xfer(p)         = frootc_xfer(p)        - m_frootc_xfer_to_litter(p)        * dt
       livestemc_xfer(p)      = livestemc_xfer(p)     - m_livestemc_xfer_to_litter(p)     * dt
       deadstemc_xfer(p)      = deadstemc_xfer(p)     - m_deadstemc_xfer_to_litter(p)     * dt
       livecrootc_xfer(p)     = livecrootc_xfer(p)    - m_livecrootc_xfer_to_litter(p)    * dt
       deadcrootc_xfer(p)     = deadcrootc_xfer(p)    - m_deadcrootc_xfer_to_litter(p)    * dt
       gresp_xfer(p)          = gresp_xfer(p)         - m_gresp_xfer_to_litter(p)         * dt

    end do ! end of pft loop

end subroutine CStateUpdate2
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CStateUpdate2h
!
! !INTERFACE:
subroutine CStateUpdate2h(num_soilc, filter_soilc, num_soilp, filter_soilp)
!
! !DESCRIPTION:
! Update all the prognostic carbon state
! variables affected by harvest mortality fluxes
!
! !USES:
   use clmtype
   use clm_time_manager, only: get_step_size
!
! !ARGUMENTS:
   implicit none
   integer, intent(in) :: num_soilc       ! number of soil columns in filter
   integer, intent(in) :: filter_soilc(:) ! filter for soil columns
   integer, intent(in) :: num_soilp       ! number of soil pfts in filter
   integer, intent(in) :: filter_soilp(:) ! filter for soil pfts
!
! !CALLED FROM:
! subroutine CNEcosystemDyn
!
! !REVISION HISTORY:
! 5/20/09: Created by Peter Thornton
!
! !LOCAL VARIABLES:
! local pointers to implicit in arrays
   real, pointer :: hrv_deadcrootc_storage_to_litr1c(:)
   real, pointer :: hrv_deadcrootc_to_cwdc(:)
   real, pointer :: hrv_deadcrootc_xfer_to_litr1c(:)
   real, pointer :: hrv_deadstemc_storage_to_litr1c(:)
   real, pointer :: hrv_deadstemc_xfer_to_litr1c(:)
   real, pointer :: hrv_frootc_storage_to_litr1c(:)
   real, pointer :: hrv_frootc_to_litr1c(:)
   real, pointer :: hrv_frootc_to_litr2c(:)
   real, pointer :: hrv_frootc_to_litr3c(:)
   real, pointer :: hrv_frootc_xfer_to_litr1c(:)
   real, pointer :: hrv_gresp_storage_to_litr1c(:)
   real, pointer :: hrv_gresp_xfer_to_litr1c(:)
   real, pointer :: hrv_leafc_storage_to_litr1c(:)
   real, pointer :: hrv_leafc_to_litr1c(:)
   real, pointer :: hrv_leafc_to_litr2c(:)
   real, pointer :: hrv_leafc_to_litr3c(:)
   real, pointer :: hrv_leafc_xfer_to_litr1c(:)
   real, pointer :: hrv_livecrootc_storage_to_litr1c(:)
   real, pointer :: hrv_livecrootc_to_cwdc(:)
   real, pointer :: hrv_livecrootc_xfer_to_litr1c(:)
   real, pointer :: hrv_livestemc_storage_to_litr1c(:)
   real, pointer :: hrv_livestemc_to_cwdc(:)
   real, pointer :: hrv_livestemc_xfer_to_litr1c(:)
   real, pointer :: hrv_deadcrootc_storage_to_litter(:)
   real, pointer :: hrv_deadcrootc_to_litter(:)
   real, pointer :: hrv_deadcrootc_xfer_to_litter(:)
   real, pointer :: hrv_deadstemc_storage_to_litter(:)
   real, pointer :: hrv_deadstemc_to_prod10c(:)
   real, pointer :: hrv_deadstemc_to_prod100c(:)
   real, pointer :: hrv_deadstemc_xfer_to_litter(:)
   real, pointer :: hrv_frootc_storage_to_litter(:)
   real, pointer :: hrv_frootc_to_litter(:)
   real, pointer :: hrv_frootc_xfer_to_litter(:)
   real, pointer :: hrv_gresp_storage_to_litter(:)
   real, pointer :: hrv_gresp_xfer_to_litter(:)
   real, pointer :: hrv_leafc_storage_to_litter(:)
   real, pointer :: hrv_leafc_to_litter(:)
   real, pointer :: hrv_leafc_xfer_to_litter(:)
   real, pointer :: hrv_livecrootc_storage_to_litter(:)
   real, pointer :: hrv_livecrootc_to_litter(:)
   real, pointer :: hrv_livecrootc_xfer_to_litter(:)
   real, pointer :: hrv_livestemc_storage_to_litter(:)
   real, pointer :: hrv_livestemc_to_litter(:)
   real, pointer :: hrv_livestemc_xfer_to_litter(:)
   real, pointer :: hrv_xsmrpool_to_atm(:)
!
! local pointers to implicit in/out arrays
   real, pointer :: cwdc(:)               ! (gC/m2) coarse woody debris C
   real, pointer :: litr1c(:)             ! (gC/m2) litter labile C
   real, pointer :: litr2c(:)             ! (gC/m2) litter cellulose C
   real, pointer :: litr3c(:)             ! (gC/m2) litter lignin C
   real, pointer :: deadcrootc(:)         ! (gC/m2) dead coarse root C
   real, pointer :: deadcrootc_storage(:) ! (gC/m2) dead coarse root C storage
   real, pointer :: deadcrootc_xfer(:)    ! (gC/m2) dead coarse root C transfer
   real, pointer :: deadstemc(:)          ! (gC/m2) dead stem C
   real, pointer :: deadstemc_storage(:)  ! (gC/m2) dead stem C storage
   real, pointer :: deadstemc_xfer(:)     ! (gC/m2) dead stem C transfer
   real, pointer :: frootc(:)             ! (gC/m2) fine root C
   real, pointer :: frootc_storage(:)     ! (gC/m2) fine root C storage
   real, pointer :: frootc_xfer(:)        ! (gC/m2) fine root C transfer
   real, pointer :: gresp_storage(:)      ! (gC/m2) growth respiration storage
   real, pointer :: gresp_xfer(:)         ! (gC/m2) growth respiration transfer
   real, pointer :: leafc(:)              ! (gC/m2) leaf C
   real, pointer :: leafc_storage(:)      ! (gC/m2) leaf C storage
   real, pointer :: leafc_xfer(:)         ! (gC/m2) leaf C transfer
   real, pointer :: livecrootc(:)         ! (gC/m2) live coarse root C
   real, pointer :: livecrootc_storage(:) ! (gC/m2) live coarse root C storage
   real, pointer :: livecrootc_xfer(:)    ! (gC/m2) live coarse root C transfer
   real, pointer :: livestemc(:)          ! (gC/m2) live stem C
   real, pointer :: livestemc_storage(:)  ! (gC/m2) live stem C storage
   real, pointer :: livestemc_xfer(:)     ! (gC/m2) live stem C transfer
   real, pointer :: xsmrpool(:)           ! (gC/m2) abstract C pool to meet excess MR demand
!
!
! local pointers to implicit out arrays
!
!
! !OTHER LOCAL VARIABLES:
   integer :: c,p      ! indices
   integer :: fp,fc    ! lake filter indices
   real:: dt       ! radiation time step (seconds)
!
!EOP
!-----------------------------------------------------------------------
    ! assign local pointers at the column level
    hrv_deadcrootc_storage_to_litr1c => clm3%g%l%c%ccf%hrv_deadcrootc_storage_to_litr1c
    hrv_deadcrootc_to_cwdc           => clm3%g%l%c%ccf%hrv_deadcrootc_to_cwdc
    hrv_deadcrootc_xfer_to_litr1c    => clm3%g%l%c%ccf%hrv_deadcrootc_xfer_to_litr1c
    hrv_deadstemc_storage_to_litr1c  => clm3%g%l%c%ccf%hrv_deadstemc_storage_to_litr1c
    hrv_deadstemc_xfer_to_litr1c     => clm3%g%l%c%ccf%hrv_deadstemc_xfer_to_litr1c
    hrv_frootc_storage_to_litr1c     => clm3%g%l%c%ccf%hrv_frootc_storage_to_litr1c
    hrv_frootc_to_litr1c             => clm3%g%l%c%ccf%hrv_frootc_to_litr1c
    hrv_frootc_to_litr2c             => clm3%g%l%c%ccf%hrv_frootc_to_litr2c
    hrv_frootc_to_litr3c             => clm3%g%l%c%ccf%hrv_frootc_to_litr3c
    hrv_frootc_xfer_to_litr1c        => clm3%g%l%c%ccf%hrv_frootc_xfer_to_litr1c
    hrv_gresp_storage_to_litr1c      => clm3%g%l%c%ccf%hrv_gresp_storage_to_litr1c
    hrv_gresp_xfer_to_litr1c         => clm3%g%l%c%ccf%hrv_gresp_xfer_to_litr1c
    hrv_leafc_storage_to_litr1c      => clm3%g%l%c%ccf%hrv_leafc_storage_to_litr1c
    hrv_leafc_to_litr1c              => clm3%g%l%c%ccf%hrv_leafc_to_litr1c
    hrv_leafc_to_litr2c              => clm3%g%l%c%ccf%hrv_leafc_to_litr2c
    hrv_leafc_to_litr3c              => clm3%g%l%c%ccf%hrv_leafc_to_litr3c
    hrv_leafc_xfer_to_litr1c         => clm3%g%l%c%ccf%hrv_leafc_xfer_to_litr1c
    hrv_livecrootc_storage_to_litr1c => clm3%g%l%c%ccf%hrv_livecrootc_storage_to_litr1c
    hrv_livecrootc_to_cwdc           => clm3%g%l%c%ccf%hrv_livecrootc_to_cwdc
    hrv_livecrootc_xfer_to_litr1c    => clm3%g%l%c%ccf%hrv_livecrootc_xfer_to_litr1c
    hrv_livestemc_storage_to_litr1c  => clm3%g%l%c%ccf%hrv_livestemc_storage_to_litr1c
    hrv_livestemc_to_cwdc            => clm3%g%l%c%ccf%hrv_livestemc_to_cwdc
    hrv_livestemc_xfer_to_litr1c     => clm3%g%l%c%ccf%hrv_livestemc_xfer_to_litr1c
    cwdc                           => clm3%g%l%c%ccs%cwdc
    litr1c                         => clm3%g%l%c%ccs%litr1c
    litr2c                         => clm3%g%l%c%ccs%litr2c
    litr3c                         => clm3%g%l%c%ccs%litr3c

    ! assign local pointers at the pft level
    hrv_deadcrootc_storage_to_litter => clm3%g%l%c%p%pcf%hrv_deadcrootc_storage_to_litter
    hrv_deadcrootc_to_litter         => clm3%g%l%c%p%pcf%hrv_deadcrootc_to_litter
    hrv_deadcrootc_xfer_to_litter    => clm3%g%l%c%p%pcf%hrv_deadcrootc_xfer_to_litter
    hrv_deadstemc_storage_to_litter  => clm3%g%l%c%p%pcf%hrv_deadstemc_storage_to_litter
    hrv_deadstemc_to_prod10c         => clm3%g%l%c%p%pcf%hrv_deadstemc_to_prod10c
    hrv_deadstemc_to_prod100c        => clm3%g%l%c%p%pcf%hrv_deadstemc_to_prod100c
    hrv_deadstemc_xfer_to_litter     => clm3%g%l%c%p%pcf%hrv_deadstemc_xfer_to_litter
    hrv_frootc_storage_to_litter     => clm3%g%l%c%p%pcf%hrv_frootc_storage_to_litter
    hrv_frootc_to_litter             => clm3%g%l%c%p%pcf%hrv_frootc_to_litter
    hrv_frootc_xfer_to_litter        => clm3%g%l%c%p%pcf%hrv_frootc_xfer_to_litter
    hrv_gresp_storage_to_litter      => clm3%g%l%c%p%pcf%hrv_gresp_storage_to_litter
    hrv_gresp_xfer_to_litter         => clm3%g%l%c%p%pcf%hrv_gresp_xfer_to_litter
    hrv_leafc_storage_to_litter      => clm3%g%l%c%p%pcf%hrv_leafc_storage_to_litter
    hrv_leafc_to_litter              => clm3%g%l%c%p%pcf%hrv_leafc_to_litter
    hrv_leafc_xfer_to_litter         => clm3%g%l%c%p%pcf%hrv_leafc_xfer_to_litter
    hrv_livecrootc_storage_to_litter => clm3%g%l%c%p%pcf%hrv_livecrootc_storage_to_litter
    hrv_livecrootc_to_litter         => clm3%g%l%c%p%pcf%hrv_livecrootc_to_litter
    hrv_livecrootc_xfer_to_litter    => clm3%g%l%c%p%pcf%hrv_livecrootc_xfer_to_litter
    hrv_livestemc_storage_to_litter  => clm3%g%l%c%p%pcf%hrv_livestemc_storage_to_litter
    hrv_livestemc_to_litter          => clm3%g%l%c%p%pcf%hrv_livestemc_to_litter
    hrv_livestemc_xfer_to_litter     => clm3%g%l%c%p%pcf%hrv_livestemc_xfer_to_litter
    hrv_xsmrpool_to_atm              => clm3%g%l%c%p%pcf%hrv_xsmrpool_to_atm
    deadcrootc                     => clm3%g%l%c%p%pcs%deadcrootc
    deadcrootc_storage             => clm3%g%l%c%p%pcs%deadcrootc_storage
    deadcrootc_xfer                => clm3%g%l%c%p%pcs%deadcrootc_xfer
    deadstemc                      => clm3%g%l%c%p%pcs%deadstemc
    deadstemc_storage              => clm3%g%l%c%p%pcs%deadstemc_storage
    deadstemc_xfer                 => clm3%g%l%c%p%pcs%deadstemc_xfer
    frootc                         => clm3%g%l%c%p%pcs%frootc
    frootc_storage                 => clm3%g%l%c%p%pcs%frootc_storage
    frootc_xfer                    => clm3%g%l%c%p%pcs%frootc_xfer
    gresp_storage                  => clm3%g%l%c%p%pcs%gresp_storage
    gresp_xfer                     => clm3%g%l%c%p%pcs%gresp_xfer
    leafc                          => clm3%g%l%c%p%pcs%leafc
    leafc_storage                  => clm3%g%l%c%p%pcs%leafc_storage
    leafc_xfer                     => clm3%g%l%c%p%pcs%leafc_xfer
    livecrootc                     => clm3%g%l%c%p%pcs%livecrootc
    livecrootc_storage             => clm3%g%l%c%p%pcs%livecrootc_storage
    livecrootc_xfer                => clm3%g%l%c%p%pcs%livecrootc_xfer
    livestemc                      => clm3%g%l%c%p%pcs%livestemc
    livestemc_storage              => clm3%g%l%c%p%pcs%livestemc_storage
    livestemc_xfer                 => clm3%g%l%c%p%pcs%livestemc_xfer
    xsmrpool                       => clm3%g%l%c%p%pcs%xsmrpool

    ! set time steps
    dt = real( get_step_size() )

    ! column loop
    do fc = 1,num_soilc
       c = filter_soilc(fc)

       ! column level carbon fluxes from harvest mortality

       ! leaf to litter
       litr1c(c) = litr1c(c) + hrv_leafc_to_litr1c(c) * dt
       litr2c(c) = litr2c(c) + hrv_leafc_to_litr2c(c) * dt
       litr3c(c) = litr3c(c) + hrv_leafc_to_litr3c(c) * dt

       ! fine root to litter
       litr1c(c) = litr1c(c) + hrv_frootc_to_litr1c(c) * dt
       litr2c(c) = litr2c(c) + hrv_frootc_to_litr2c(c) * dt
       litr3c(c) = litr3c(c) + hrv_frootc_to_litr3c(c) * dt

       ! wood to CWD
       cwdc(c) = cwdc(c) + hrv_livestemc_to_cwdc(c)  * dt
       cwdc(c) = cwdc(c) + hrv_livecrootc_to_cwdc(c) * dt
       cwdc(c) = cwdc(c) + hrv_deadcrootc_to_cwdc(c) * dt

       ! wood to product pools - states updated in CNWoodProducts()

       ! storage pools to litter
       litr1c(c) = litr1c(c) + hrv_leafc_storage_to_litr1c(c)      * dt
       litr1c(c) = litr1c(c) + hrv_frootc_storage_to_litr1c(c)     * dt
       litr1c(c) = litr1c(c) + hrv_livestemc_storage_to_litr1c(c)  * dt
       litr1c(c) = litr1c(c) + hrv_deadstemc_storage_to_litr1c(c)  * dt
       litr1c(c) = litr1c(c) + hrv_livecrootc_storage_to_litr1c(c) * dt
       litr1c(c) = litr1c(c) + hrv_deadcrootc_storage_to_litr1c(c) * dt
       litr1c(c) = litr1c(c) + hrv_gresp_storage_to_litr1c(c)      * dt

       ! transfer pools to litter
       litr1c(c) = litr1c(c) + hrv_leafc_xfer_to_litr1c(c)      * dt
       litr1c(c) = litr1c(c) + hrv_frootc_xfer_to_litr1c(c)     * dt
       litr1c(c) = litr1c(c) + hrv_livestemc_xfer_to_litr1c(c)  * dt
       litr1c(c) = litr1c(c) + hrv_deadstemc_xfer_to_litr1c(c)  * dt
       litr1c(c) = litr1c(c) + hrv_livecrootc_xfer_to_litr1c(c) * dt
       litr1c(c) = litr1c(c) + hrv_deadcrootc_xfer_to_litr1c(c) * dt
       litr1c(c) = litr1c(c) + hrv_gresp_xfer_to_litr1c(c)      * dt

    end do ! end of columns loop

    ! pft loop
    do fp = 1,num_soilp
       p = filter_soilp(fp)

       ! pft-level carbon fluxes from harvest mortality
       ! displayed pools
       leafc(p)               = leafc(p)              - hrv_leafc_to_litter(p)              * dt
       frootc(p)              = frootc(p)             - hrv_frootc_to_litter(p)             * dt
       livestemc(p)           = livestemc(p)          - hrv_livestemc_to_litter(p)          * dt
       deadstemc(p)           = deadstemc(p)          - hrv_deadstemc_to_prod10c(p)         * dt
       deadstemc(p)           = deadstemc(p)          - hrv_deadstemc_to_prod100c(p)        * dt
       livecrootc(p)          = livecrootc(p)         - hrv_livecrootc_to_litter(p)         * dt
       deadcrootc(p)          = deadcrootc(p)         - hrv_deadcrootc_to_litter(p)         * dt
       
       ! xsmrpool
       xsmrpool(p)            = xsmrpool(p)           - hrv_xsmrpool_to_atm(p)              * dt

       ! storage pools
       leafc_storage(p)       = leafc_storage(p)      - hrv_leafc_storage_to_litter(p)      * dt
       frootc_storage(p)      = frootc_storage(p)     - hrv_frootc_storage_to_litter(p)     * dt
       livestemc_storage(p)   = livestemc_storage(p)  - hrv_livestemc_storage_to_litter(p)  * dt
       deadstemc_storage(p)   = deadstemc_storage(p)  - hrv_deadstemc_storage_to_litter(p)  * dt
       livecrootc_storage(p)  = livecrootc_storage(p) - hrv_livecrootc_storage_to_litter(p) * dt
       deadcrootc_storage(p)  = deadcrootc_storage(p) - hrv_deadcrootc_storage_to_litter(p) * dt
       gresp_storage(p)       = gresp_storage(p)      - hrv_gresp_storage_to_litter(p)      * dt

       ! transfer pools
       leafc_xfer(p)          = leafc_xfer(p)         - hrv_leafc_xfer_to_litter(p)         * dt
       frootc_xfer(p)         = frootc_xfer(p)        - hrv_frootc_xfer_to_litter(p)        * dt
       livestemc_xfer(p)      = livestemc_xfer(p)     - hrv_livestemc_xfer_to_litter(p)     * dt
       deadstemc_xfer(p)      = deadstemc_xfer(p)     - hrv_deadstemc_xfer_to_litter(p)     * dt
       livecrootc_xfer(p)     = livecrootc_xfer(p)    - hrv_livecrootc_xfer_to_litter(p)    * dt
       deadcrootc_xfer(p)     = deadcrootc_xfer(p)    - hrv_deadcrootc_xfer_to_litter(p)    * dt
       gresp_xfer(p)          = gresp_xfer(p)         - hrv_gresp_xfer_to_litter(p)         * dt

    end do ! end of pft loop

end subroutine CStateUpdate2h
!-----------------------------------------------------------------------

end module CNCStateUpdate2Mod
