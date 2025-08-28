module clmtypeInitMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: clmtypeInitMod
!
! !DESCRIPTION:
! Allocate clmtype components and initialize them to signaling NaN.
!
! !USES:
  use nanMod      , only : nan, bigint
  use clmtype
  use clm_varpar  , only : nlevsno, nlevgrnd, numpft
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: initClmtype
!
! !REVISION HISTORY:
! Created by Peter Thornton and Mariana Vertenstein
! Modified by Colette L. Heald (05/06) for VOC emission factors
! 3/17/08 David Lawrence, changed nlevsoi to nlevgrnd where appropriate
!
! !PRIVATE MEMBER FUNCTIONS:
  private :: init_pft_type
  private :: init_column_type
  private :: init_landunit_type
  private :: init_gridcell_type
  private :: init_pft_ecophys_constants
  private :: init_pft_pstate_type
  private :: init_pft_epv_type
  private :: init_pft_estate_type
  private :: init_pft_cstate_type
  private :: init_pft_nstate_type
  private :: init_pft_cflux_type
  private :: init_pft_nflux_type
  private :: init_column_pstate_type
  private :: init_column_estate_type
  private :: init_column_wstate_type
  private :: init_column_cstate_type
  private :: init_column_nstate_type
  private :: init_column_wflux_type
  private :: init_column_cflux_type
  private :: init_column_nflux_type
!EOP
!----------------------------------------------------

contains

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: initClmtype
!
! !INTERFACE:
  subroutine initClmtype(begg,endg,begl,endl,begc,endc,begp,endp)
!
! !DESCRIPTION:
! Initialize clmtype components to signaling nan
! The following clmtype components should NOT be initialized here
! since they are set in routine clm_map which is called before this
! routine is invoked
!    *%area, *%wt, *%wtlnd, *%wtxy, *%ixy, *%jxy, *%mxy, %snindex
!    *%ifspecial, *%ityplun, *%itype
!    *%pfti, *%pftf, *%pftn
!    *%coli, *%colf, *%coln
!    *%luni, *%lunf, *%lunn
!
! !USES:
!   use decompMod , only : get_proc_bounds, get_proc_global : gkw: ???
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! LOCAL VARAIBLES:
    integer :: begp, endp   ! per-proc beginning and ending pft indices
    integer :: begc, endc   ! per-proc beginning and ending column indices
    integer :: begl, endl   ! per-proc beginning and ending landunit indices
    integer :: begg, endg   ! per-proc gridcell ending gridcell indices
    integer :: numg         ! total number of gridcells across all processors
    integer :: numl         ! total number of landunits across all processors
    integer :: numc         ! total number of columns across all processors
    integer :: nump         ! total number of pfts across all processors
!------------------------------------------------------------------------

    ! Determine necessary indices

! gkw: ???

!   call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)
!   call get_proc_global(numg, numl, numc, nump)

!   begg = 1 ; endg = 1 ! gkw: all this must come from driver or elsewhere
!   begl = 1 ; endl = 1
!   begc = 1 ; endc = 1
!   begp = 1 ; endp = 1
    
!   numg = 1 ; numl = 1 ; numc = 1 ; nump = 1

    call init_pft_type     (begp, endp, clm3%g%l%c%p)
    call init_column_type  (begc, endc, clm3%g%l%c)
    call init_landunit_type(begl, endl, clm3%g%l)
    call init_gridcell_type(begg, endg, clm3%g)

    ! pft ecophysiological constants

    call init_pft_ecophys_constants()

    ! carbon balance structures (pft and column levels)

    call init_carbon_balance_type(begp, endp, clm3%g%l%c%p%pcbal)
    call init_carbon_balance_type(begc, endc, clm3%g%l%c%ccbal)

    ! nitrogen balance structures (pft and column levels)

    call init_nitrogen_balance_type(begp, endp, clm3%g%l%c%p%pnbal)
    call init_nitrogen_balance_type(begc, endc, clm3%g%l%c%cnbal)

    ! pft physical state variables at pft level

    call init_pft_pstate_type(begp, endp, clm3%g%l%c%p%pps)

    ! pft ecophysiological variables (only at the pft level for now)
    call init_pft_epv_type(begp, endp, clm3%g%l%c%p%pepv)

    ! pft energy state variables at the pft level

    call init_pft_estate_type(begp, endp, clm3%g%l%c%p%pes)

    ! pft carbon state variables at the pft level and averaged to the column

    call init_pft_cstate_type(begp, endp, clm3%g%l%c%p%pcs)
    call init_pft_cstate_type(begc, endc, clm3%g%l%c%ccs%pcs_a)

    ! pft nitrogen state variables at the pft level and averaged to the column

    call init_pft_nstate_type(begp, endp, clm3%g%l%c%p%pns)
    call init_pft_nstate_type(begc, endc, clm3%g%l%c%cns%pns_a)

    ! pft carbon flux variables at pft level and averaged to column

    call init_pft_cflux_type(begp, endp, clm3%g%l%c%p%pcf)
    call init_pft_cflux_type(begc, endc, clm3%g%l%c%ccf%pcf_a)

    ! pft nitrogen flux variables at pft level and averaged to column

    call init_pft_nflux_type(begp, endp, clm3%g%l%c%p%pnf)
    call init_pft_nflux_type(begc, endc, clm3%g%l%c%cnf%pnf_a)

    ! column physical state variables at column level

    call init_column_pstate_type(begc, endc, clm3%g%l%c%cps)

    ! column energy state variables at column level

    call init_column_estate_type(begc, endc, clm3%g%l%c%ces)

    ! column water state variables at column level

    call init_column_wstate_type(begc, endc, clm3%g%l%c%cws)

    ! column carbon state variables at column level

    call init_column_cstate_type(begc, endc, clm3%g%l%c%ccs)

    ! column nitrogen state variables at column level

    call init_column_nstate_type(begc, endc, clm3%g%l%c%cns)

    ! column water flux variables at column level

    call init_column_wflux_type(begc, endc, clm3%g%l%c%cwf)

    ! column carbon flux variables at column level

    call init_column_cflux_type(begc, endc, clm3%g%l%c%ccf)

    ! column nitrogen flux variables at column level

    call init_column_nflux_type(begc, endc, clm3%g%l%c%cnf)

  end subroutine initClmtype

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_type
!
! !INTERFACE:
  subroutine init_pft_type (beg, end, p)
!
! !DESCRIPTION:
! Initialize components of pft_type structure
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type(pft_type), intent(inout):: p
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(p%gridcell(beg:end),p%wtgcell(beg:end))
    allocate(p%landunit(beg:end))
    allocate(p%column  (beg:end),p%wtcol  (beg:end))

    allocate(p%itype(beg:end))

  end subroutine init_pft_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_column_type
!
! !INTERFACE:
  subroutine init_column_type (beg, end, c)
!
! !DESCRIPTION:
! Initialize components of column_type structure
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type(column_type), intent(inout):: c
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

   allocate(c%gridcell(beg:end),c%wtgcell(beg:end))
   allocate(c%landunit(beg:end),c%wtlunit(beg:end))

   allocate(c%pfti(beg:end),c%pftf(beg:end),c%npfts(beg:end))

   allocate(c%itype(beg:end))

  end subroutine init_column_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_landunit_type
!
! !INTERFACE:
  subroutine init_landunit_type (beg, end,l)
!
! !DESCRIPTION:
! Initialize components of landunit_type structure
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type(landunit_type), intent(inout):: l
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

   allocate(l%itype(beg:end))
   allocate(l%ifspecial(beg:end))

  end subroutine init_landunit_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_gridcell_type
!
! !INTERFACE:
  subroutine init_gridcell_type (beg, end,g)
!
! !DESCRIPTION:
! Initialize components of gridcell_type structure
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type(gridcell_type), intent(inout):: g
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

   allocate(g%gindex(beg:end))

   allocate(g%forc_ndep(beg:end))

  end subroutine init_gridcell_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_carbon_balance_type
!
! !INTERFACE:
  subroutine init_carbon_balance_type(beg, end, cbal)
!
! !DESCRIPTION:
! Initialize carbon balance variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type(carbon_balance_type), intent(inout):: cbal
!
! !REVISION HISTORY:
! Created by Peter Thornton, 12/11/2003
!
!EOP
!------------------------------------------------------------------------

    allocate(cbal%begcb(beg:end))
    allocate(cbal%endcb(beg:end))
    allocate(cbal%errcb(beg:end))

    cbal%begcb(beg:end) = nan
    cbal%endcb(beg:end) = nan
    cbal%errcb(beg:end) = nan

  end subroutine init_carbon_balance_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_nitrogen_balance_type
!
! !INTERFACE:
  subroutine init_nitrogen_balance_type(beg, end, nbal)
!
! !DESCRIPTION:
! Initialize nitrogen balance variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type(nitrogen_balance_type), intent(inout):: nbal
!
! !REVISION HISTORY:
! Created by Peter Thornton, 12/11/2003
!
!EOP
!------------------------------------------------------------------------

    allocate(nbal%begnb(beg:end))
    allocate(nbal%endnb(beg:end))
    allocate(nbal%errnb(beg:end))

    nbal%begnb(beg:end) = nan
    nbal%endnb(beg:end) = nan
    nbal%errnb(beg:end) = nan

  end subroutine init_nitrogen_balance_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_ecophys_constants
!
! !INTERFACE:
  subroutine init_pft_ecophys_constants()
!
! !DESCRIPTION:
! Initialize pft physical state
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(pftcon%noveg(0:numpft))
    allocate(pftcon%tree(0:numpft))
    allocate(pftcon%fnitr(0:numpft))
    allocate(pftcon%c3psn(0:numpft))
    allocate(pftcon%vcmx25(0:numpft))
    allocate(pftcon%mp(0:numpft))
    allocate(pftcon%qe25(0:numpft))
    allocate(pftcon%z0mr(0:numpft))
    allocate(pftcon%displar(0:numpft))
    allocate(pftcon%slatop(0:numpft))
    allocate(pftcon%dsladlai(0:numpft))
    allocate(pftcon%leafcn(0:numpft))
    allocate(pftcon%flnr(0:numpft))
    allocate(pftcon%woody(0:numpft))
    allocate(pftcon%lflitcn(0:numpft))
    allocate(pftcon%frootcn(0:numpft))
    allocate(pftcon%livewdcn(0:numpft))
    allocate(pftcon%deadwdcn(0:numpft))
    allocate(pftcon%froot_leaf(0:numpft))
    allocate(pftcon%stem_leaf(0:numpft))
    allocate(pftcon%croot_stem(0:numpft))
    allocate(pftcon%flivewd(0:numpft))
    allocate(pftcon%fcur(0:numpft))
    allocate(pftcon%lf_flab(0:numpft))
    allocate(pftcon%lf_fcel(0:numpft))
    allocate(pftcon%lf_flig(0:numpft))
    allocate(pftcon%fr_flab(0:numpft))
    allocate(pftcon%fr_fcel(0:numpft))
    allocate(pftcon%fr_flig(0:numpft))
    allocate(pftcon%dw_fcel(0:numpft))
    allocate(pftcon%dw_flig(0:numpft))
    allocate(pftcon%leaf_long(0:numpft))
    allocate(pftcon%evergreen(0:numpft))
    allocate(pftcon%stress_decid(0:numpft))
    allocate(pftcon%season_decid(0:numpft))
    allocate(pftcon%resist(0:numpft))
    allocate(pftcon%dwood(0:numpft))

    allocate(pftcon%xl(0:numpft))
    allocate(pftcon%rhol(0:numpft))
    allocate(pftcon%rhos(0:numpft))
    allocate(pftcon%taul(0:numpft))
    allocate(pftcon%taus(0:numpft))

    pftcon%noveg(:) = bigint
    pftcon%tree(:) = bigint
    pftcon%fnitr(:) = nan
    pftcon%c3psn(:) = nan
    pftcon%vcmx25(:) = nan
    pftcon%mp(:) = nan
    pftcon%qe25(:) = nan
    pftcon%z0mr(:) = nan
    pftcon%displar(:) = nan
    pftcon%slatop(:) = nan
    pftcon%dsladlai(:) = nan
    pftcon%leafcn(:) = nan
    pftcon%flnr(:) = nan
    pftcon%woody(:) = nan
    pftcon%lflitcn(:) = nan
    pftcon%frootcn(:) = nan
    pftcon%livewdcn(:) = nan
    pftcon%deadwdcn(:) = nan
    pftcon%froot_leaf(:) = nan
    pftcon%stem_leaf(:) = nan
    pftcon%croot_stem(:) = nan
    pftcon%flivewd(:) = nan
    pftcon%fcur(:) = nan
    pftcon%lf_flab(:) = nan
    pftcon%lf_fcel(:) = nan
    pftcon%lf_flig(:) = nan
    pftcon%fr_flab(:) = nan
    pftcon%fr_fcel(:) = nan
    pftcon%fr_flig(:) = nan
    pftcon%dw_fcel(:) = nan
    pftcon%dw_flig(:) = nan
    pftcon%leaf_long(:) = nan
    pftcon%evergreen(:) = nan
    pftcon%stress_decid(:) = nan
    pftcon%season_decid(:) = nan
    pftcon%resist(:) = nan
    pftcon%dwood(:) = nan

    pftcon%xl(:) = nan
    pftcon%rhol(:) = nan
    pftcon%rhos(:) = nan
    pftcon%taul(:) = nan
    pftcon%taus(:) = nan

  end subroutine init_pft_ecophys_constants

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_pstate_type
!
! !INTERFACE:
  subroutine init_pft_pstate_type(beg, end, pps)
!
! !DESCRIPTION:
! Initialize pft physical state
!
! !USES:
    use clm_varcon, only : spval
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_pstate_type), intent(inout):: pps
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(pps%frac_veg_nosno(beg:end))
    allocate(pps%frac_veg_nosno_alb(beg:end))
    allocate(pps%rootfr(beg:end,1:nlevgrnd))
    allocate(pps%laisun(beg:end))
    allocate(pps%laisha(beg:end))
    allocate(pps%fsun(beg:end))
    allocate(pps%tlai(beg:end))
    allocate(pps%tsai(beg:end))
    allocate(pps%elai(beg:end))
    allocate(pps%esai(beg:end))
    allocate(pps%fwet(beg:end))
    allocate(pps%fdry(beg:end))
    allocate(pps%htop(beg:end))
    allocate(pps%hbot(beg:end))
    allocate(pps%forc_hgt_u_pft(beg:end))

    pps%frac_veg_nosno(beg:end) = bigint
    pps%frac_veg_nosno_alb(beg:end) = 0
    pps%rootfr(beg:end,:nlevgrnd) = spval
    pps%laisun(beg:end) = nan
    pps%laisha(beg:end) = nan
    pps%fsun(beg:end) = spval
    pps%tlai(beg:end) = 0.
    pps%tsai(beg:end) = 0.
    pps%elai(beg:end) = 0.
    pps%esai(beg:end) = 0.
    pps%fwet(beg:end) = nan
    pps%fdry(beg:end) = nan
    pps%htop(beg:end) = 0.
    pps%hbot(beg:end) = 0.
    pps%forc_hgt_u_pft(beg:end) = nan

  end subroutine init_pft_pstate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_epv_type
!
! !INTERFACE:
  subroutine init_pft_epv_type(beg, end, pepv)
!
! !DESCRIPTION:
! Initialize pft ecophysiological variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_epv_type), intent(inout):: pepv
!
! !REVISION HISTORY:
! Created by Peter Thornton
!
!EOP
!------------------------------------------------------------------------

    allocate(pepv%dormant_flag(beg:end))
    allocate(pepv%days_active(beg:end))
    allocate(pepv%onset_flag(beg:end))
    allocate(pepv%onset_counter(beg:end))
    allocate(pepv%onset_gddflag(beg:end))
    allocate(pepv%onset_fdd(beg:end))
    allocate(pepv%onset_gdd(beg:end))
    allocate(pepv%onset_swi(beg:end))
    allocate(pepv%offset_flag(beg:end))
    allocate(pepv%offset_counter(beg:end))
    allocate(pepv%offset_fdd(beg:end))
    allocate(pepv%offset_swi(beg:end))
    allocate(pepv%lgsf(beg:end))
    allocate(pepv%bglfr(beg:end))
    allocate(pepv%bgtr(beg:end))
    allocate(pepv%dayl(beg:end))
    allocate(pepv%prev_dayl(beg:end))
    allocate(pepv%annavg_t2m(beg:end))
    allocate(pepv%tempavg_t2m(beg:end))
    allocate(pepv%gpp(beg:end))
    allocate(pepv%availc(beg:end))
    allocate(pepv%xsmrpool_recover(beg:end))
    allocate(pepv%alloc_pnow(beg:end))
    allocate(pepv%c_allometry(beg:end))
    allocate(pepv%n_allometry(beg:end))
    allocate(pepv%plant_ndemand(beg:end))
    allocate(pepv%tempsum_potential_gpp(beg:end))
    allocate(pepv%annsum_potential_gpp(beg:end))
    allocate(pepv%tempmax_retransn(beg:end))
    allocate(pepv%annmax_retransn(beg:end))
    allocate(pepv%avail_retransn(beg:end))
    allocate(pepv%plant_nalloc(beg:end))
    allocate(pepv%plant_calloc(beg:end))
    allocate(pepv%excess_cflux(beg:end))
    allocate(pepv%downreg(beg:end))
    allocate(pepv%prev_leafc_to_litter(beg:end))
    allocate(pepv%prev_frootc_to_litter(beg:end))
    allocate(pepv%tempsum_npp(beg:end))
    allocate(pepv%annsum_npp(beg:end))

    pepv%dormant_flag(beg:end) = nan
    pepv%days_active(beg:end) = nan
    pepv%onset_flag(beg:end) = nan
    pepv%onset_counter(beg:end) = nan
    pepv%onset_gddflag(beg:end) = nan
    pepv%onset_fdd(beg:end) = nan
    pepv%onset_gdd(beg:end) = nan
    pepv%onset_swi(beg:end) = nan
    pepv%offset_flag(beg:end) = nan
    pepv%offset_counter(beg:end) = nan
    pepv%offset_fdd(beg:end) = nan
    pepv%offset_swi(beg:end) = nan
    pepv%lgsf(beg:end) = nan
    pepv%bglfr(beg:end) = nan
    pepv%bgtr(beg:end) = nan
    pepv%dayl(beg:end) = nan
    pepv%prev_dayl(beg:end) = nan
    pepv%annavg_t2m(beg:end) = nan
    pepv%tempavg_t2m(beg:end) = nan
    pepv%gpp(beg:end) = nan
    pepv%availc(beg:end) = nan
    pepv%xsmrpool_recover(beg:end) = nan
    pepv%alloc_pnow(beg:end) = nan
    pepv%c_allometry(beg:end) = nan
    pepv%n_allometry(beg:end) = nan
    pepv%plant_ndemand(beg:end) = nan
    pepv%tempsum_potential_gpp(beg:end) = nan
    pepv%annsum_potential_gpp(beg:end) = nan
    pepv%tempmax_retransn(beg:end) = nan
    pepv%annmax_retransn(beg:end) = nan
    pepv%avail_retransn(beg:end) = nan
    pepv%plant_nalloc(beg:end) = nan
    pepv%plant_calloc(beg:end) = nan
    pepv%excess_cflux(beg:end) = nan
    pepv%downreg(beg:end) = nan
    pepv%prev_leafc_to_litter(beg:end) = nan
    pepv%prev_frootc_to_litter(beg:end) = nan
    pepv%tempsum_npp(beg:end) = nan
    pepv%annsum_npp(beg:end) = nan
    
  end subroutine init_pft_epv_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_estate_type
!
! !INTERFACE:
  subroutine init_pft_estate_type(beg, end, pes)
!
! !DESCRIPTION:
! Initialize pft energy state
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_estate_type), intent(inout):: pes
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(pes%t_ref2m(beg:end))

    pes%t_ref2m(beg:end) = nan

  end subroutine init_pft_estate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_cstate_type
!
! !INTERFACE:
  subroutine init_pft_cstate_type(beg, end, pcs)
!
! !DESCRIPTION:
! Initialize pft carbon state
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_cstate_type), intent(inout):: pcs !pft carbon state
!
! !REVISION HISTORY:
! Created by Peter Thornton
!
!EOP
!------------------------------------------------------------------------

    allocate(pcs%leafc(beg:end))
    allocate(pcs%leafc_storage(beg:end))
    allocate(pcs%leafc_xfer(beg:end))
    allocate(pcs%frootc(beg:end))
    allocate(pcs%frootc_storage(beg:end))
    allocate(pcs%frootc_xfer(beg:end))
    allocate(pcs%livestemc(beg:end))
    allocate(pcs%livestemc_storage(beg:end))
    allocate(pcs%livestemc_xfer(beg:end))
    allocate(pcs%deadstemc(beg:end))
    allocate(pcs%deadstemc_storage(beg:end))
    allocate(pcs%deadstemc_xfer(beg:end))
    allocate(pcs%livecrootc(beg:end))
    allocate(pcs%livecrootc_storage(beg:end))
    allocate(pcs%livecrootc_xfer(beg:end))
    allocate(pcs%deadcrootc(beg:end))
    allocate(pcs%deadcrootc_storage(beg:end))
    allocate(pcs%deadcrootc_xfer(beg:end))
    allocate(pcs%gresp_storage(beg:end))
    allocate(pcs%gresp_xfer(beg:end))
    allocate(pcs%cpool(beg:end))
    allocate(pcs%xsmrpool(beg:end))
    allocate(pcs%pft_ctrunc(beg:end))
    allocate(pcs%dispvegc(beg:end))
    allocate(pcs%storvegc(beg:end))
    allocate(pcs%totvegc(beg:end))
    allocate(pcs%totpftc(beg:end))
    allocate(pcs%leafcmax(beg:end))

    pcs%leafc(beg:end) = nan
    pcs%leafc_storage(beg:end) = nan
    pcs%leafc_xfer(beg:end) = nan
    pcs%frootc(beg:end) = nan
    pcs%frootc_storage(beg:end) = nan
    pcs%frootc_xfer(beg:end) = nan
    pcs%livestemc(beg:end) = nan
    pcs%livestemc_storage(beg:end) = nan
    pcs%livestemc_xfer(beg:end) = nan
    pcs%deadstemc(beg:end) = nan
    pcs%deadstemc_storage(beg:end) = nan
    pcs%deadstemc_xfer(beg:end) = nan
    pcs%livecrootc(beg:end) = nan
    pcs%livecrootc_storage(beg:end) = nan
    pcs%livecrootc_xfer(beg:end) = nan
    pcs%deadcrootc(beg:end) = nan
    pcs%deadcrootc_storage(beg:end) = nan
    pcs%deadcrootc_xfer(beg:end) = nan
    pcs%gresp_storage(beg:end) = nan
    pcs%gresp_xfer(beg:end) = nan
    pcs%cpool(beg:end) = nan
    pcs%xsmrpool(beg:end) = nan
    pcs%pft_ctrunc(beg:end) = nan
    pcs%dispvegc(beg:end) = nan
    pcs%storvegc(beg:end) = nan
    pcs%totvegc(beg:end) = nan
    pcs%totpftc(beg:end) = nan
    pcs%leafcmax(beg:end) = nan

  end subroutine init_pft_cstate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_nstate_type
!
! !INTERFACE:
  subroutine init_pft_nstate_type(beg, end, pns)
!
! !DESCRIPTION:
! Initialize pft nitrogen state
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_nstate_type), intent(inout):: pns !pft nitrogen state
!
! !REVISION HISTORY:
! Created by Peter Thornton
!
!EOP
!------------------------------------------------------------------------

    allocate(pns%leafn(beg:end))
    allocate(pns%leafn_storage(beg:end))
    allocate(pns%leafn_xfer(beg:end))
    allocate(pns%frootn(beg:end))
    allocate(pns%frootn_storage(beg:end))
    allocate(pns%frootn_xfer(beg:end))
    allocate(pns%livestemn(beg:end))
    allocate(pns%livestemn_storage(beg:end))
    allocate(pns%livestemn_xfer(beg:end))
    allocate(pns%deadstemn(beg:end))
    allocate(pns%deadstemn_storage(beg:end))
    allocate(pns%deadstemn_xfer(beg:end))
    allocate(pns%livecrootn(beg:end))
    allocate(pns%livecrootn_storage(beg:end))
    allocate(pns%livecrootn_xfer(beg:end))
    allocate(pns%deadcrootn(beg:end))
    allocate(pns%deadcrootn_storage(beg:end))
    allocate(pns%deadcrootn_xfer(beg:end))
    allocate(pns%retransn(beg:end))
    allocate(pns%npool(beg:end))
    allocate(pns%pft_ntrunc(beg:end))
    allocate(pns%dispvegn(beg:end))
    allocate(pns%storvegn(beg:end))
    allocate(pns%totvegn(beg:end))
    allocate(pns%totpftn(beg:end))

    pns%leafn(beg:end) = nan
    pns%leafn_storage(beg:end) = nan
    pns%leafn_xfer(beg:end) = nan
    pns%frootn(beg:end) = nan
    pns%frootn_storage(beg:end) = nan
    pns%frootn_xfer(beg:end) = nan
    pns%livestemn(beg:end) = nan
    pns%livestemn_storage(beg:end) = nan
    pns%livestemn_xfer(beg:end) = nan
    pns%deadstemn(beg:end) = nan
    pns%deadstemn_storage(beg:end) = nan
    pns%deadstemn_xfer(beg:end) = nan
    pns%livecrootn(beg:end) = nan
    pns%livecrootn_storage(beg:end) = nan
    pns%livecrootn_xfer(beg:end) = nan
    pns%deadcrootn(beg:end) = nan
    pns%deadcrootn_storage(beg:end) = nan
    pns%deadcrootn_xfer(beg:end) = nan
    pns%retransn(beg:end) = nan
    pns%npool(beg:end) = nan
    pns%pft_ntrunc(beg:end) = nan
    pns%dispvegn(beg:end) = nan
    pns%storvegn(beg:end) = nan
    pns%totvegn(beg:end) = nan
    pns%totpftn(beg:end) = nan

  end subroutine init_pft_nstate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_cflux_type
!
! !INTERFACE:
  subroutine init_pft_cflux_type(beg, end, pcf)
!
! !DESCRIPTION:
! Initialize pft carbon flux variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_cflux_type), intent(inout) :: pcf
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(pcf%psnsun(beg:end))
    allocate(pcf%psnsha(beg:end))
    allocate(pcf%fco2(beg:end))

    allocate(pcf%m_leafc_to_litter(beg:end))
    allocate(pcf%m_frootc_to_litter(beg:end))
    allocate(pcf%m_leafc_storage_to_litter(beg:end))
    allocate(pcf%m_frootc_storage_to_litter(beg:end))
    allocate(pcf%m_livestemc_storage_to_litter(beg:end))
    allocate(pcf%m_deadstemc_storage_to_litter(beg:end))
    allocate(pcf%m_livecrootc_storage_to_litter(beg:end))
    allocate(pcf%m_deadcrootc_storage_to_litter(beg:end))
    allocate(pcf%m_leafc_xfer_to_litter(beg:end))
    allocate(pcf%m_frootc_xfer_to_litter(beg:end))
    allocate(pcf%m_livestemc_xfer_to_litter(beg:end))
    allocate(pcf%m_deadstemc_xfer_to_litter(beg:end))
    allocate(pcf%m_livecrootc_xfer_to_litter(beg:end))
    allocate(pcf%m_deadcrootc_xfer_to_litter(beg:end))
    allocate(pcf%m_livestemc_to_litter(beg:end))
    allocate(pcf%m_deadstemc_to_litter(beg:end))
    allocate(pcf%m_livecrootc_to_litter(beg:end))
    allocate(pcf%m_deadcrootc_to_litter(beg:end))
    allocate(pcf%m_gresp_storage_to_litter(beg:end))
    allocate(pcf%m_gresp_xfer_to_litter(beg:end))
    allocate(pcf%hrv_leafc_to_litter(beg:end))             
    allocate(pcf%hrv_leafc_storage_to_litter(beg:end))     
    allocate(pcf%hrv_leafc_xfer_to_litter(beg:end))        
    allocate(pcf%hrv_frootc_to_litter(beg:end))            
    allocate(pcf%hrv_frootc_storage_to_litter(beg:end))    
    allocate(pcf%hrv_frootc_xfer_to_litter(beg:end))       
    allocate(pcf%hrv_livestemc_to_litter(beg:end))         
    allocate(pcf%hrv_livestemc_storage_to_litter(beg:end)) 
    allocate(pcf%hrv_livestemc_xfer_to_litter(beg:end))    
    allocate(pcf%hrv_deadstemc_to_prod10c(beg:end))        
    allocate(pcf%hrv_deadstemc_to_prod100c(beg:end))       
    allocate(pcf%hrv_deadstemc_storage_to_litter(beg:end)) 
    allocate(pcf%hrv_deadstemc_xfer_to_litter(beg:end))    
    allocate(pcf%hrv_livecrootc_to_litter(beg:end))        
    allocate(pcf%hrv_livecrootc_storage_to_litter(beg:end))
    allocate(pcf%hrv_livecrootc_xfer_to_litter(beg:end))   
    allocate(pcf%hrv_deadcrootc_to_litter(beg:end))        
    allocate(pcf%hrv_deadcrootc_storage_to_litter(beg:end))
    allocate(pcf%hrv_deadcrootc_xfer_to_litter(beg:end))   
    allocate(pcf%hrv_gresp_storage_to_litter(beg:end))     
    allocate(pcf%hrv_gresp_xfer_to_litter(beg:end))        
    allocate(pcf%hrv_xsmrpool_to_atm(beg:end))                 
    allocate(pcf%m_leafc_to_fire(beg:end))
    allocate(pcf%m_frootc_to_fire(beg:end))
    allocate(pcf%m_leafc_storage_to_fire(beg:end))
    allocate(pcf%m_frootc_storage_to_fire(beg:end))
    allocate(pcf%m_livestemc_storage_to_fire(beg:end))
    allocate(pcf%m_deadstemc_storage_to_fire(beg:end))
    allocate(pcf%m_livecrootc_storage_to_fire(beg:end))
    allocate(pcf%m_deadcrootc_storage_to_fire(beg:end))
    allocate(pcf%m_leafc_xfer_to_fire(beg:end))
    allocate(pcf%m_frootc_xfer_to_fire(beg:end))
    allocate(pcf%m_livestemc_xfer_to_fire(beg:end))
    allocate(pcf%m_deadstemc_xfer_to_fire(beg:end))
    allocate(pcf%m_livecrootc_xfer_to_fire(beg:end))
    allocate(pcf%m_deadcrootc_xfer_to_fire(beg:end))
    allocate(pcf%m_livestemc_to_fire(beg:end))
    allocate(pcf%m_deadstemc_to_fire(beg:end))
    allocate(pcf%m_deadstemc_to_litter_fire(beg:end))
    allocate(pcf%m_livecrootc_to_fire(beg:end))
    allocate(pcf%m_deadcrootc_to_fire(beg:end))
    allocate(pcf%m_deadcrootc_to_litter_fire(beg:end))
    allocate(pcf%m_gresp_storage_to_fire(beg:end))
    allocate(pcf%m_gresp_xfer_to_fire(beg:end))
    allocate(pcf%leafc_xfer_to_leafc(beg:end))
    allocate(pcf%frootc_xfer_to_frootc(beg:end))
    allocate(pcf%livestemc_xfer_to_livestemc(beg:end))
    allocate(pcf%deadstemc_xfer_to_deadstemc(beg:end))
    allocate(pcf%livecrootc_xfer_to_livecrootc(beg:end))
    allocate(pcf%deadcrootc_xfer_to_deadcrootc(beg:end))
    allocate(pcf%leafc_to_litter(beg:end))
    allocate(pcf%frootc_to_litter(beg:end))
    allocate(pcf%leaf_mr(beg:end))
    allocate(pcf%froot_mr(beg:end))
    allocate(pcf%livestem_mr(beg:end))
    allocate(pcf%livecroot_mr(beg:end))
    allocate(pcf%leaf_curmr(beg:end))
    allocate(pcf%froot_curmr(beg:end))
    allocate(pcf%livestem_curmr(beg:end))
    allocate(pcf%livecroot_curmr(beg:end))
    allocate(pcf%leaf_xsmr(beg:end))
    allocate(pcf%froot_xsmr(beg:end))
    allocate(pcf%livestem_xsmr(beg:end))
    allocate(pcf%livecroot_xsmr(beg:end))
    allocate(pcf%psnsun_to_cpool(beg:end))
    allocate(pcf%psnshade_to_cpool(beg:end))
    allocate(pcf%cpool_to_xsmrpool(beg:end))
    allocate(pcf%cpool_to_leafc(beg:end))
    allocate(pcf%cpool_to_leafc_storage(beg:end))
    allocate(pcf%cpool_to_frootc(beg:end))
    allocate(pcf%cpool_to_frootc_storage(beg:end))
    allocate(pcf%cpool_to_livestemc(beg:end))
    allocate(pcf%cpool_to_livestemc_storage(beg:end))
    allocate(pcf%cpool_to_deadstemc(beg:end))
    allocate(pcf%cpool_to_deadstemc_storage(beg:end))
    allocate(pcf%cpool_to_livecrootc(beg:end))
    allocate(pcf%cpool_to_livecrootc_storage(beg:end))
    allocate(pcf%cpool_to_deadcrootc(beg:end))
    allocate(pcf%cpool_to_deadcrootc_storage(beg:end))
    allocate(pcf%cpool_to_gresp_storage(beg:end))
    allocate(pcf%cpool_leaf_gr(beg:end))
    allocate(pcf%cpool_leaf_storage_gr(beg:end))
    allocate(pcf%transfer_leaf_gr(beg:end))
    allocate(pcf%cpool_froot_gr(beg:end))
    allocate(pcf%cpool_froot_storage_gr(beg:end))
    allocate(pcf%transfer_froot_gr(beg:end))
    allocate(pcf%cpool_livestem_gr(beg:end))
    allocate(pcf%cpool_livestem_storage_gr(beg:end))
    allocate(pcf%transfer_livestem_gr(beg:end))
    allocate(pcf%cpool_deadstem_gr(beg:end))
    allocate(pcf%cpool_deadstem_storage_gr(beg:end))
    allocate(pcf%transfer_deadstem_gr(beg:end))
    allocate(pcf%cpool_livecroot_gr(beg:end))
    allocate(pcf%cpool_livecroot_storage_gr(beg:end))
    allocate(pcf%transfer_livecroot_gr(beg:end))
    allocate(pcf%cpool_deadcroot_gr(beg:end))
    allocate(pcf%cpool_deadcroot_storage_gr(beg:end))
    allocate(pcf%transfer_deadcroot_gr(beg:end))
    allocate(pcf%leafc_storage_to_xfer(beg:end))
    allocate(pcf%frootc_storage_to_xfer(beg:end))
    allocate(pcf%livestemc_storage_to_xfer(beg:end))
    allocate(pcf%deadstemc_storage_to_xfer(beg:end))
    allocate(pcf%livecrootc_storage_to_xfer(beg:end))
    allocate(pcf%deadcrootc_storage_to_xfer(beg:end))
    allocate(pcf%gresp_storage_to_xfer(beg:end))
    allocate(pcf%livestemc_to_deadstemc(beg:end))
    allocate(pcf%livecrootc_to_deadcrootc(beg:end))
    allocate(pcf%gpp(beg:end))
    allocate(pcf%mr(beg:end))
    allocate(pcf%current_gr(beg:end))
    allocate(pcf%transfer_gr(beg:end))
    allocate(pcf%storage_gr(beg:end))
    allocate(pcf%gr(beg:end))
    allocate(pcf%ar(beg:end))
    allocate(pcf%rr(beg:end))
    allocate(pcf%npp(beg:end))
    allocate(pcf%agnpp(beg:end))
    allocate(pcf%bgnpp(beg:end))
    allocate(pcf%litfall(beg:end))
    allocate(pcf%vegfire(beg:end))
    allocate(pcf%wood_harvestc(beg:end))
    allocate(pcf%pft_cinputs(beg:end))
    allocate(pcf%pft_coutputs(beg:end))
    allocate(pcf%pft_fire_closs(beg:end))

    pcf%psnsun(beg:end) = nan
    pcf%psnsha(beg:end) = nan
    pcf%fco2(beg:end) = 0.

    pcf%m_leafc_to_litter(beg:end) = nan
    pcf%m_frootc_to_litter(beg:end) = nan
    pcf%m_leafc_storage_to_litter(beg:end) = nan
    pcf%m_frootc_storage_to_litter(beg:end) = nan
    pcf%m_livestemc_storage_to_litter(beg:end) = nan
    pcf%m_deadstemc_storage_to_litter(beg:end) = nan
    pcf%m_livecrootc_storage_to_litter(beg:end) = nan
    pcf%m_deadcrootc_storage_to_litter(beg:end) = nan
    pcf%m_leafc_xfer_to_litter(beg:end) = nan
    pcf%m_frootc_xfer_to_litter(beg:end) = nan
    pcf%m_livestemc_xfer_to_litter(beg:end) = nan
    pcf%m_deadstemc_xfer_to_litter(beg:end) = nan
    pcf%m_livecrootc_xfer_to_litter(beg:end) = nan
    pcf%m_deadcrootc_xfer_to_litter(beg:end) = nan
    pcf%m_livestemc_to_litter(beg:end) = nan
    pcf%m_deadstemc_to_litter(beg:end) = nan
    pcf%m_livecrootc_to_litter(beg:end) = nan
    pcf%m_deadcrootc_to_litter(beg:end) = nan
    pcf%m_gresp_storage_to_litter(beg:end) = nan
    pcf%m_gresp_xfer_to_litter(beg:end) = nan
    pcf%hrv_leafc_to_litter(beg:end) = nan             
    pcf%hrv_leafc_storage_to_litter(beg:end) = nan     
    pcf%hrv_leafc_xfer_to_litter(beg:end) = nan        
    pcf%hrv_frootc_to_litter(beg:end) = nan            
    pcf%hrv_frootc_storage_to_litter(beg:end) = nan    
    pcf%hrv_frootc_xfer_to_litter(beg:end) = nan       
    pcf%hrv_livestemc_to_litter(beg:end) = nan         
    pcf%hrv_livestemc_storage_to_litter(beg:end) = nan 
    pcf%hrv_livestemc_xfer_to_litter(beg:end) = nan    
    pcf%hrv_deadstemc_to_prod10c(beg:end) = nan        
    pcf%hrv_deadstemc_to_prod100c(beg:end) = nan       
    pcf%hrv_deadstemc_storage_to_litter(beg:end) = nan 
    pcf%hrv_deadstemc_xfer_to_litter(beg:end) = nan    
    pcf%hrv_livecrootc_to_litter(beg:end) = nan        
    pcf%hrv_livecrootc_storage_to_litter(beg:end) = nan
    pcf%hrv_livecrootc_xfer_to_litter(beg:end) = nan   
    pcf%hrv_deadcrootc_to_litter(beg:end) = nan        
    pcf%hrv_deadcrootc_storage_to_litter(beg:end) = nan
    pcf%hrv_deadcrootc_xfer_to_litter(beg:end) = nan   
    pcf%hrv_gresp_storage_to_litter(beg:end) = nan     
    pcf%hrv_gresp_xfer_to_litter(beg:end) = nan        
    pcf%hrv_xsmrpool_to_atm(beg:end) = nan                 
    pcf%m_leafc_to_fire(beg:end) = nan
    pcf%m_frootc_to_fire(beg:end) = nan
    pcf%m_leafc_storage_to_fire(beg:end) = nan
    pcf%m_frootc_storage_to_fire(beg:end) = nan
    pcf%m_livestemc_storage_to_fire(beg:end) = nan
    pcf%m_deadstemc_storage_to_fire(beg:end) = nan
    pcf%m_livecrootc_storage_to_fire(beg:end) = nan
    pcf%m_deadcrootc_storage_to_fire(beg:end) = nan
    pcf%m_leafc_xfer_to_fire(beg:end) = nan
    pcf%m_frootc_xfer_to_fire(beg:end) = nan
    pcf%m_livestemc_xfer_to_fire(beg:end) = nan
    pcf%m_deadstemc_xfer_to_fire(beg:end) = nan
    pcf%m_livecrootc_xfer_to_fire(beg:end) = nan
    pcf%m_deadcrootc_xfer_to_fire(beg:end) = nan
    pcf%m_livestemc_to_fire(beg:end) = nan
    pcf%m_deadstemc_to_fire(beg:end) = nan
    pcf%m_deadstemc_to_litter_fire(beg:end) = nan
    pcf%m_livecrootc_to_fire(beg:end) = nan
    pcf%m_deadcrootc_to_fire(beg:end) = nan
    pcf%m_deadcrootc_to_litter_fire(beg:end) = nan
    pcf%m_gresp_storage_to_fire(beg:end) = nan
    pcf%m_gresp_xfer_to_fire(beg:end) = nan
    pcf%leafc_xfer_to_leafc(beg:end) = nan
    pcf%frootc_xfer_to_frootc(beg:end) = nan
    pcf%livestemc_xfer_to_livestemc(beg:end) = nan
    pcf%deadstemc_xfer_to_deadstemc(beg:end) = nan
    pcf%livecrootc_xfer_to_livecrootc(beg:end) = nan
    pcf%deadcrootc_xfer_to_deadcrootc(beg:end) = nan
    pcf%leafc_to_litter(beg:end) = nan
    pcf%frootc_to_litter(beg:end) = nan
    pcf%leaf_mr(beg:end) = nan
    pcf%froot_mr(beg:end) = nan
    pcf%livestem_mr(beg:end) = nan
    pcf%livecroot_mr(beg:end) = nan
    pcf%leaf_curmr(beg:end) = nan
    pcf%froot_curmr(beg:end) = nan
    pcf%livestem_curmr(beg:end) = nan
    pcf%livecroot_curmr(beg:end) = nan
    pcf%leaf_xsmr(beg:end) = nan
    pcf%froot_xsmr(beg:end) = nan
    pcf%livestem_xsmr(beg:end) = nan
    pcf%livecroot_xsmr(beg:end) = nan
    pcf%psnsun_to_cpool(beg:end) = nan
    pcf%psnshade_to_cpool(beg:end) = nan
    pcf%cpool_to_xsmrpool(beg:end) = nan
    pcf%cpool_to_leafc(beg:end) = nan
    pcf%cpool_to_leafc_storage(beg:end) = nan
    pcf%cpool_to_frootc(beg:end) = nan
    pcf%cpool_to_frootc_storage(beg:end) = nan
    pcf%cpool_to_livestemc(beg:end) = nan
    pcf%cpool_to_livestemc_storage(beg:end) = nan
    pcf%cpool_to_deadstemc(beg:end) = nan
    pcf%cpool_to_deadstemc_storage(beg:end) = nan
    pcf%cpool_to_livecrootc(beg:end) = nan
    pcf%cpool_to_livecrootc_storage(beg:end) = nan
    pcf%cpool_to_deadcrootc(beg:end) = nan
    pcf%cpool_to_deadcrootc_storage(beg:end) = nan
    pcf%cpool_to_gresp_storage(beg:end) = nan
    pcf%cpool_leaf_gr(beg:end) = nan
    pcf%cpool_leaf_storage_gr(beg:end) = nan
    pcf%transfer_leaf_gr(beg:end) = nan
    pcf%cpool_froot_gr(beg:end) = nan
    pcf%cpool_froot_storage_gr(beg:end) = nan
    pcf%transfer_froot_gr(beg:end) = nan
    pcf%cpool_livestem_gr(beg:end) = nan
    pcf%cpool_livestem_storage_gr(beg:end) = nan
    pcf%transfer_livestem_gr(beg:end) = nan
    pcf%cpool_deadstem_gr(beg:end) = nan
    pcf%cpool_deadstem_storage_gr(beg:end) = nan
    pcf%transfer_deadstem_gr(beg:end) = nan
    pcf%cpool_livecroot_gr(beg:end) = nan
    pcf%cpool_livecroot_storage_gr(beg:end) = nan
    pcf%transfer_livecroot_gr(beg:end) = nan
    pcf%cpool_deadcroot_gr(beg:end) = nan
    pcf%cpool_deadcroot_storage_gr(beg:end) = nan
    pcf%transfer_deadcroot_gr(beg:end) = nan
    pcf%leafc_storage_to_xfer(beg:end) = nan
    pcf%frootc_storage_to_xfer(beg:end) = nan
    pcf%livestemc_storage_to_xfer(beg:end) = nan
    pcf%deadstemc_storage_to_xfer(beg:end) = nan
    pcf%livecrootc_storage_to_xfer(beg:end) = nan
    pcf%deadcrootc_storage_to_xfer(beg:end) = nan
    pcf%gresp_storage_to_xfer(beg:end) = nan
    pcf%livestemc_to_deadstemc(beg:end) = nan
    pcf%livecrootc_to_deadcrootc(beg:end) = nan
    pcf%gpp(beg:end) = nan
    pcf%mr(beg:end) = nan
    pcf%current_gr(beg:end) = nan
    pcf%transfer_gr(beg:end) = nan
    pcf%storage_gr(beg:end) = nan
    pcf%gr(beg:end) = nan
    pcf%ar(beg:end) = nan
    pcf%rr(beg:end) = nan
    pcf%npp(beg:end) = nan
    pcf%agnpp(beg:end) = nan
    pcf%bgnpp(beg:end) = nan
    pcf%litfall(beg:end) = nan
    pcf%vegfire(beg:end) = nan
    pcf%wood_harvestc(beg:end) = nan
    pcf%pft_cinputs(beg:end) = nan
    pcf%pft_coutputs(beg:end) = nan
    pcf%pft_fire_closs(beg:end) = nan

  end subroutine init_pft_cflux_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_nflux_type
!
! !INTERFACE:
  subroutine init_pft_nflux_type(beg, end, pnf)
!
! !DESCRIPTION:
! Initialize pft nitrogen flux variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_nflux_type), intent(inout) :: pnf
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(pnf%m_leafn_to_litter(beg:end))
    allocate(pnf%m_frootn_to_litter(beg:end))
    allocate(pnf%m_leafn_storage_to_litter(beg:end))
    allocate(pnf%m_frootn_storage_to_litter(beg:end))
    allocate(pnf%m_livestemn_storage_to_litter(beg:end))
    allocate(pnf%m_deadstemn_storage_to_litter(beg:end))
    allocate(pnf%m_livecrootn_storage_to_litter(beg:end))
    allocate(pnf%m_deadcrootn_storage_to_litter(beg:end))
    allocate(pnf%m_leafn_xfer_to_litter(beg:end))
    allocate(pnf%m_frootn_xfer_to_litter(beg:end))
    allocate(pnf%m_livestemn_xfer_to_litter(beg:end))
    allocate(pnf%m_deadstemn_xfer_to_litter(beg:end))
    allocate(pnf%m_livecrootn_xfer_to_litter(beg:end))
    allocate(pnf%m_deadcrootn_xfer_to_litter(beg:end))
    allocate(pnf%m_livestemn_to_litter(beg:end))
    allocate(pnf%m_deadstemn_to_litter(beg:end))
    allocate(pnf%m_livecrootn_to_litter(beg:end))
    allocate(pnf%m_deadcrootn_to_litter(beg:end))
    allocate(pnf%m_retransn_to_litter(beg:end))
    allocate(pnf%hrv_leafn_to_litter(beg:end))             
    allocate(pnf%hrv_frootn_to_litter(beg:end))            
    allocate(pnf%hrv_leafn_storage_to_litter(beg:end))     
    allocate(pnf%hrv_frootn_storage_to_litter(beg:end))    
    allocate(pnf%hrv_livestemn_storage_to_litter(beg:end)) 
    allocate(pnf%hrv_deadstemn_storage_to_litter(beg:end)) 
    allocate(pnf%hrv_livecrootn_storage_to_litter(beg:end))
    allocate(pnf%hrv_deadcrootn_storage_to_litter(beg:end))
    allocate(pnf%hrv_leafn_xfer_to_litter(beg:end))        
    allocate(pnf%hrv_frootn_xfer_to_litter(beg:end))       
    allocate(pnf%hrv_livestemn_xfer_to_litter(beg:end))    
    allocate(pnf%hrv_deadstemn_xfer_to_litter(beg:end))    
    allocate(pnf%hrv_livecrootn_xfer_to_litter(beg:end))   
    allocate(pnf%hrv_deadcrootn_xfer_to_litter(beg:end))   
    allocate(pnf%hrv_livestemn_to_litter(beg:end))         
    allocate(pnf%hrv_deadstemn_to_prod10n(beg:end))        
    allocate(pnf%hrv_deadstemn_to_prod100n(beg:end))       
    allocate(pnf%hrv_livecrootn_to_litter(beg:end))        
    allocate(pnf%hrv_deadcrootn_to_litter(beg:end))        
    allocate(pnf%hrv_retransn_to_litter(beg:end))              
    allocate(pnf%m_leafn_to_fire(beg:end))
    allocate(pnf%m_frootn_to_fire(beg:end))
    allocate(pnf%m_leafn_storage_to_fire(beg:end))
    allocate(pnf%m_frootn_storage_to_fire(beg:end))
    allocate(pnf%m_livestemn_storage_to_fire(beg:end))
    allocate(pnf%m_deadstemn_storage_to_fire(beg:end))
    allocate(pnf%m_livecrootn_storage_to_fire(beg:end))
    allocate(pnf%m_deadcrootn_storage_to_fire(beg:end))
    allocate(pnf%m_leafn_xfer_to_fire(beg:end))
    allocate(pnf%m_frootn_xfer_to_fire(beg:end))
    allocate(pnf%m_livestemn_xfer_to_fire(beg:end))
    allocate(pnf%m_deadstemn_xfer_to_fire(beg:end))
    allocate(pnf%m_livecrootn_xfer_to_fire(beg:end))
    allocate(pnf%m_deadcrootn_xfer_to_fire(beg:end))
    allocate(pnf%m_livestemn_to_fire(beg:end))
    allocate(pnf%m_deadstemn_to_fire(beg:end))
    allocate(pnf%m_deadstemn_to_litter_fire(beg:end))
    allocate(pnf%m_livecrootn_to_fire(beg:end))
    allocate(pnf%m_deadcrootn_to_fire(beg:end))
    allocate(pnf%m_deadcrootn_to_litter_fire(beg:end))
    allocate(pnf%m_retransn_to_fire(beg:end))
    allocate(pnf%leafn_xfer_to_leafn(beg:end))
    allocate(pnf%frootn_xfer_to_frootn(beg:end))
    allocate(pnf%livestemn_xfer_to_livestemn(beg:end))
    allocate(pnf%deadstemn_xfer_to_deadstemn(beg:end))
    allocate(pnf%livecrootn_xfer_to_livecrootn(beg:end))
    allocate(pnf%deadcrootn_xfer_to_deadcrootn(beg:end))
    allocate(pnf%leafn_to_litter(beg:end))
    allocate(pnf%leafn_to_retransn(beg:end))
    allocate(pnf%frootn_to_litter(beg:end))
    allocate(pnf%retransn_to_npool(beg:end))
    allocate(pnf%sminn_to_npool(beg:end))
    allocate(pnf%npool_to_leafn(beg:end))
    allocate(pnf%npool_to_leafn_storage(beg:end))
    allocate(pnf%npool_to_frootn(beg:end))
    allocate(pnf%npool_to_frootn_storage(beg:end))
    allocate(pnf%npool_to_livestemn(beg:end))
    allocate(pnf%npool_to_livestemn_storage(beg:end))
    allocate(pnf%npool_to_deadstemn(beg:end))
    allocate(pnf%npool_to_deadstemn_storage(beg:end))
    allocate(pnf%npool_to_livecrootn(beg:end))
    allocate(pnf%npool_to_livecrootn_storage(beg:end))
    allocate(pnf%npool_to_deadcrootn(beg:end))
    allocate(pnf%npool_to_deadcrootn_storage(beg:end))
    allocate(pnf%leafn_storage_to_xfer(beg:end))
    allocate(pnf%frootn_storage_to_xfer(beg:end))
    allocate(pnf%livestemn_storage_to_xfer(beg:end))
    allocate(pnf%deadstemn_storage_to_xfer(beg:end))
    allocate(pnf%livecrootn_storage_to_xfer(beg:end))
    allocate(pnf%deadcrootn_storage_to_xfer(beg:end))
    allocate(pnf%livestemn_to_deadstemn(beg:end))
    allocate(pnf%livestemn_to_retransn(beg:end))
    allocate(pnf%livecrootn_to_deadcrootn(beg:end))
    allocate(pnf%livecrootn_to_retransn(beg:end))
    allocate(pnf%ndeploy(beg:end))
    allocate(pnf%pft_ninputs(beg:end))
    allocate(pnf%pft_noutputs(beg:end))
    allocate(pnf%wood_harvestn(beg:end))
    allocate(pnf%pft_fire_nloss(beg:end))

    pnf%m_leafn_to_litter(beg:end) = nan
    pnf%m_frootn_to_litter(beg:end) = nan
    pnf%m_leafn_storage_to_litter(beg:end) = nan
    pnf%m_frootn_storage_to_litter(beg:end) = nan
    pnf%m_livestemn_storage_to_litter(beg:end) = nan
    pnf%m_deadstemn_storage_to_litter(beg:end) = nan
    pnf%m_livecrootn_storage_to_litter(beg:end) = nan
    pnf%m_deadcrootn_storage_to_litter(beg:end) = nan
    pnf%m_leafn_xfer_to_litter(beg:end) = nan
    pnf%m_frootn_xfer_to_litter(beg:end) = nan
    pnf%m_livestemn_xfer_to_litter(beg:end) = nan
    pnf%m_deadstemn_xfer_to_litter(beg:end) = nan
    pnf%m_livecrootn_xfer_to_litter(beg:end) = nan
    pnf%m_deadcrootn_xfer_to_litter(beg:end) = nan
    pnf%m_livestemn_to_litter(beg:end) = nan
    pnf%m_deadstemn_to_litter(beg:end) = nan
    pnf%m_livecrootn_to_litter(beg:end) = nan
    pnf%m_deadcrootn_to_litter(beg:end) = nan
    pnf%m_retransn_to_litter(beg:end) = nan
    pnf%hrv_leafn_to_litter(beg:end) = nan             
    pnf%hrv_frootn_to_litter(beg:end) = nan            
    pnf%hrv_leafn_storage_to_litter(beg:end) = nan     
    pnf%hrv_frootn_storage_to_litter(beg:end) = nan    
    pnf%hrv_livestemn_storage_to_litter(beg:end) = nan 
    pnf%hrv_deadstemn_storage_to_litter(beg:end) = nan 
    pnf%hrv_livecrootn_storage_to_litter(beg:end) = nan
    pnf%hrv_deadcrootn_storage_to_litter(beg:end) = nan
    pnf%hrv_leafn_xfer_to_litter(beg:end) = nan        
    pnf%hrv_frootn_xfer_to_litter(beg:end) = nan       
    pnf%hrv_livestemn_xfer_to_litter(beg:end) = nan    
    pnf%hrv_deadstemn_xfer_to_litter(beg:end) = nan    
    pnf%hrv_livecrootn_xfer_to_litter(beg:end) = nan   
    pnf%hrv_deadcrootn_xfer_to_litter(beg:end) = nan   
    pnf%hrv_livestemn_to_litter(beg:end) = nan         
    pnf%hrv_deadstemn_to_prod10n(beg:end) = nan        
    pnf%hrv_deadstemn_to_prod100n(beg:end) = nan       
    pnf%hrv_livecrootn_to_litter(beg:end) = nan        
    pnf%hrv_deadcrootn_to_litter(beg:end) = nan        
    pnf%hrv_retransn_to_litter(beg:end) = nan           
    pnf%m_leafn_to_fire(beg:end) = nan
    pnf%m_frootn_to_fire(beg:end) = nan
    pnf%m_leafn_storage_to_fire(beg:end) = nan
    pnf%m_frootn_storage_to_fire(beg:end) = nan
    pnf%m_livestemn_storage_to_fire(beg:end) = nan
    pnf%m_deadstemn_storage_to_fire(beg:end) = nan
    pnf%m_livecrootn_storage_to_fire(beg:end) = nan
    pnf%m_deadcrootn_storage_to_fire(beg:end) = nan
    pnf%m_leafn_xfer_to_fire(beg:end) = nan
    pnf%m_frootn_xfer_to_fire(beg:end) = nan
    pnf%m_livestemn_xfer_to_fire(beg:end) = nan
    pnf%m_deadstemn_xfer_to_fire(beg:end) = nan
    pnf%m_livecrootn_xfer_to_fire(beg:end) = nan
    pnf%m_deadcrootn_xfer_to_fire(beg:end) = nan
    pnf%m_livestemn_to_fire(beg:end) = nan
    pnf%m_deadstemn_to_fire(beg:end) = nan
    pnf%m_deadstemn_to_litter_fire(beg:end) = nan
    pnf%m_livecrootn_to_fire(beg:end) = nan
    pnf%m_deadcrootn_to_fire(beg:end) = nan
    pnf%m_deadcrootn_to_litter_fire(beg:end) = nan
    pnf%m_retransn_to_fire(beg:end) = nan
    pnf%leafn_xfer_to_leafn(beg:end) = nan
    pnf%frootn_xfer_to_frootn(beg:end) = nan
    pnf%livestemn_xfer_to_livestemn(beg:end) = nan
    pnf%deadstemn_xfer_to_deadstemn(beg:end) = nan
    pnf%livecrootn_xfer_to_livecrootn(beg:end) = nan
    pnf%deadcrootn_xfer_to_deadcrootn(beg:end) = nan
    pnf%leafn_to_litter(beg:end) = nan
    pnf%leafn_to_retransn(beg:end) = nan
    pnf%frootn_to_litter(beg:end) = nan
    pnf%retransn_to_npool(beg:end) = nan
    pnf%sminn_to_npool(beg:end) = nan
    pnf%npool_to_leafn(beg:end) = nan
    pnf%npool_to_leafn_storage(beg:end) = nan
    pnf%npool_to_frootn(beg:end) = nan
    pnf%npool_to_frootn_storage(beg:end) = nan
    pnf%npool_to_livestemn(beg:end) = nan
    pnf%npool_to_livestemn_storage(beg:end) = nan
    pnf%npool_to_deadstemn(beg:end) = nan
    pnf%npool_to_deadstemn_storage(beg:end) = nan
    pnf%npool_to_livecrootn(beg:end) = nan
    pnf%npool_to_livecrootn_storage(beg:end) = nan
    pnf%npool_to_deadcrootn(beg:end) = nan
    pnf%npool_to_deadcrootn_storage(beg:end) = nan
    pnf%leafn_storage_to_xfer(beg:end) = nan
    pnf%frootn_storage_to_xfer(beg:end) = nan
    pnf%livestemn_storage_to_xfer(beg:end) = nan
    pnf%deadstemn_storage_to_xfer(beg:end) = nan
    pnf%livecrootn_storage_to_xfer(beg:end) = nan
    pnf%deadcrootn_storage_to_xfer(beg:end) = nan
    pnf%livestemn_to_deadstemn(beg:end) = nan
    pnf%livestemn_to_retransn(beg:end) = nan
    pnf%livecrootn_to_deadcrootn(beg:end) = nan
    pnf%livecrootn_to_retransn(beg:end) = nan
    pnf%ndeploy(beg:end) = nan
    pnf%pft_ninputs(beg:end) = nan
    pnf%pft_noutputs(beg:end) = nan
    pnf%wood_harvestn(beg:end) = nan
    pnf%pft_fire_nloss(beg:end) = nan

  end subroutine init_pft_nflux_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_column_pstate_type
!
! !INTERFACE:
  subroutine init_column_pstate_type(beg, end, cps)
!
! !DESCRIPTION:
! Initialize column physical state variables
!
! !USES:
    use clm_varcon, only : spval
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (column_pstate_type), intent(inout):: cps
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(cps%snowdp(beg:end))
    allocate(cps%dz(beg:end,-nlevsno+1:nlevgrnd))
    allocate(cps%wf(beg:end))
    allocate(cps%psisat(beg:end,nlevgrnd))
    allocate(cps%psiwilt(beg:end))
    allocate(cps%soilpsi(beg:end,nlevgrnd))
    allocate(cps%fpi(beg:end))
    allocate(cps%fpg(beg:end))
    allocate(cps%annsum_counter(beg:end))
    allocate(cps%cannsum_npp(beg:end))
    allocate(cps%cannavg_t2m(beg:end))
    allocate(cps%me(beg:end))
    allocate(cps%fire_prob(beg:end))
    allocate(cps%mean_fire_prob(beg:end))
    allocate(cps%fireseasonl(beg:end))
    allocate(cps%farea_burned(beg:end))
    allocate(cps%ann_farea_burned(beg:end))

    cps%snowdp(beg:end) = nan
    cps%dz(beg:end,-nlevsno+1:nlevgrnd) = nan
    cps%wf(beg:end) = nan
    cps%psisat(beg:end,1:nlevgrnd) = nan
    cps%psiwilt(beg:end) = nan
    cps%soilpsi(beg:end,1:nlevgrnd) = spval
    cps%fpi(beg:end) = nan
    cps%fpg(beg:end) = nan
    cps%annsum_counter(beg:end) = nan
    cps%cannsum_npp(beg:end) = nan
    cps%cannavg_t2m(beg:end) = nan
    cps%me(beg:end) = nan
    cps%fire_prob(beg:end) = nan
    cps%mean_fire_prob(beg:end) = nan
    cps%fireseasonl(beg:end) = nan
    cps%farea_burned(beg:end) = nan
    cps%ann_farea_burned(beg:end) = nan
  end subroutine init_column_pstate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_column_estate_type
!
! !INTERFACE:
  subroutine init_column_estate_type(beg, end, ces)
!
! !DESCRIPTION:
! Initialize column energy state variables
!
! !USES:
    use clm_varcon, only : spval
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (column_estate_type), intent(inout):: ces
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------
    allocate(ces%t_grnd(beg:end))
    allocate(ces%t_soisno(beg:end,-nlevsno+1:nlevgrnd))

    ces%t_grnd(beg:end)    = nan
    ces%t_soisno(beg:end,-nlevsno+1:nlevgrnd) = spval

  end subroutine init_column_estate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_column_wstate_type
!
! !INTERFACE:
  subroutine init_column_wstate_type(beg, end, cws)
!
! !DESCRIPTION:
! Initialize column water state variables
!
! !USES:
    use clm_varcon, only : spval
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (column_wstate_type), intent(inout):: cws !column water state
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(cws%h2osoi_liq(beg:end))

    cws%h2osoi_liq(beg:end)= spval

  end subroutine init_column_wstate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_column_cstate_type
!
! !INTERFACE:
  subroutine init_column_cstate_type(beg, end, ccs)
!
! !DESCRIPTION:
! Initialize column carbon state variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (column_cstate_type), intent(inout):: ccs
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(ccs%cwdc(beg:end))
    allocate(ccs%litr1c(beg:end))
    allocate(ccs%litr2c(beg:end))
    allocate(ccs%litr3c(beg:end))
    allocate(ccs%soil1c(beg:end))
    allocate(ccs%soil2c(beg:end))
    allocate(ccs%soil3c(beg:end))
    allocate(ccs%soil4c(beg:end))
    allocate(ccs%seedc(beg:end))
    allocate(ccs%col_ctrunc(beg:end))
    allocate(ccs%prod10c(beg:end))
    allocate(ccs%prod100c(beg:end))
    allocate(ccs%totprodc(beg:end))
    allocate(ccs%totlitc(beg:end))
    allocate(ccs%totsomc(beg:end))
    allocate(ccs%totecosysc(beg:end))
    allocate(ccs%totcolc(beg:end))

!   ccs%soilc(beg:end) = nan
    ccs%cwdc(beg:end) = nan
    ccs%litr1c(beg:end) = nan
    ccs%litr2c(beg:end) = nan
    ccs%litr3c(beg:end) = nan
    ccs%soil1c(beg:end) = nan
    ccs%soil2c(beg:end) = nan
    ccs%soil3c(beg:end) = nan
    ccs%soil4c(beg:end) = nan
    ccs%seedc(beg:end) = nan
    ccs%col_ctrunc(beg:end) = nan
    ccs%prod10c(beg:end) = nan
    ccs%prod100c(beg:end) = nan
    ccs%totprodc(beg:end) = nan
    ccs%totlitc(beg:end) = nan
    ccs%totsomc(beg:end) = nan
    ccs%totecosysc(beg:end) = nan
    ccs%totcolc(beg:end) = nan

  end subroutine init_column_cstate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_column_nstate_type
!
! !INTERFACE:
  subroutine init_column_nstate_type(beg, end, cns)
!
! !DESCRIPTION:
! Initialize column nitrogen state variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (column_nstate_type), intent(inout):: cns
!
! !REVISION HISTORY:
! Created by Peter Thornton
!
!EOP
!------------------------------------------------------------------------

    allocate(cns%cwdn(beg:end))
    allocate(cns%litr1n(beg:end))
    allocate(cns%litr2n(beg:end))
    allocate(cns%litr3n(beg:end))
    allocate(cns%soil1n(beg:end))
    allocate(cns%soil2n(beg:end))
    allocate(cns%soil3n(beg:end))
    allocate(cns%soil4n(beg:end))
    allocate(cns%sminn(beg:end))
    allocate(cns%col_ntrunc(beg:end))
    allocate(cns%seedn(beg:end))
    allocate(cns%prod10n(beg:end))
    allocate(cns%prod100n(beg:end))
    allocate(cns%totprodn(beg:end))
    allocate(cns%totlitn(beg:end))
    allocate(cns%totsomn(beg:end))
    allocate(cns%totecosysn(beg:end))
    allocate(cns%totcoln(beg:end))

    cns%cwdn(beg:end) = nan
    cns%litr1n(beg:end) = nan
    cns%litr2n(beg:end) = nan
    cns%litr3n(beg:end) = nan
    cns%soil1n(beg:end) = nan
    cns%soil2n(beg:end) = nan
    cns%soil3n(beg:end) = nan
    cns%soil4n(beg:end) = nan
    cns%sminn(beg:end) = nan
    cns%col_ntrunc(beg:end) = nan
    cns%seedn(beg:end) = nan
    cns%prod10n(beg:end) = nan
    cns%prod100n(beg:end) = nan
    cns%totprodn(beg:end) = nan
    cns%totlitn(beg:end) = nan
    cns%totsomn(beg:end) = nan
    cns%totecosysn(beg:end) = nan
    cns%totcoln(beg:end) = nan

  end subroutine init_column_nstate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_column_wflux_type
!
! !INTERFACE:
  subroutine init_column_wflux_type(beg, end, cwf)
!
! !DESCRIPTION:
! Initialize column water flux variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (column_wflux_type), intent(inout):: cwf
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(cwf%qflx_drain(beg:end))

    cwf%qflx_drain(beg:end) = nan

  end subroutine init_column_wflux_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_column_cflux_type
!
! !INTERFACE:
  subroutine init_column_cflux_type(beg, end, ccf)
!
! !DESCRIPTION:
! Initialize column carbon flux variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (column_cflux_type), intent(inout):: ccf
!
! !REVISION HISTORY:
! Created by Peter Thornton
!
!EOP
!------------------------------------------------------------------------

    allocate(ccf%m_leafc_to_litr1c(beg:end))
    allocate(ccf%m_leafc_to_litr2c(beg:end))
    allocate(ccf%m_leafc_to_litr3c(beg:end))
    allocate(ccf%m_frootc_to_litr1c(beg:end))
    allocate(ccf%m_frootc_to_litr2c(beg:end))
    allocate(ccf%m_frootc_to_litr3c(beg:end))
    allocate(ccf%m_leafc_storage_to_litr1c(beg:end))
    allocate(ccf%m_frootc_storage_to_litr1c(beg:end))
    allocate(ccf%m_livestemc_storage_to_litr1c(beg:end))
    allocate(ccf%m_deadstemc_storage_to_litr1c(beg:end))
    allocate(ccf%m_livecrootc_storage_to_litr1c(beg:end))
    allocate(ccf%m_deadcrootc_storage_to_litr1c(beg:end))
    allocate(ccf%m_leafc_xfer_to_litr1c(beg:end))
    allocate(ccf%m_frootc_xfer_to_litr1c(beg:end))
    allocate(ccf%m_livestemc_xfer_to_litr1c(beg:end))
    allocate(ccf%m_deadstemc_xfer_to_litr1c(beg:end))
    allocate(ccf%m_livecrootc_xfer_to_litr1c(beg:end))
    allocate(ccf%m_deadcrootc_xfer_to_litr1c(beg:end))
    allocate(ccf%m_livestemc_to_cwdc(beg:end))
    allocate(ccf%m_deadstemc_to_cwdc(beg:end))
    allocate(ccf%m_livecrootc_to_cwdc(beg:end))
    allocate(ccf%m_deadcrootc_to_cwdc(beg:end))
    allocate(ccf%m_gresp_storage_to_litr1c(beg:end))
    allocate(ccf%m_gresp_xfer_to_litr1c(beg:end))
    allocate(ccf%m_deadstemc_to_cwdc_fire(beg:end))
    allocate(ccf%m_deadcrootc_to_cwdc_fire(beg:end))
    allocate(ccf%hrv_leafc_to_litr1c(beg:end))             
    allocate(ccf%hrv_leafc_to_litr2c(beg:end))             
    allocate(ccf%hrv_leafc_to_litr3c(beg:end))             
    allocate(ccf%hrv_frootc_to_litr1c(beg:end))            
    allocate(ccf%hrv_frootc_to_litr2c(beg:end))            
    allocate(ccf%hrv_frootc_to_litr3c(beg:end))            
    allocate(ccf%hrv_livestemc_to_cwdc(beg:end))           
    allocate(ccf%hrv_deadstemc_to_prod10c(beg:end))        
    allocate(ccf%hrv_deadstemc_to_prod100c(beg:end))       
    allocate(ccf%hrv_livecrootc_to_cwdc(beg:end))          
    allocate(ccf%hrv_deadcrootc_to_cwdc(beg:end))          
    allocate(ccf%hrv_leafc_storage_to_litr1c(beg:end))     
    allocate(ccf%hrv_frootc_storage_to_litr1c(beg:end))    
    allocate(ccf%hrv_livestemc_storage_to_litr1c(beg:end)) 
    allocate(ccf%hrv_deadstemc_storage_to_litr1c(beg:end)) 
    allocate(ccf%hrv_livecrootc_storage_to_litr1c(beg:end))
    allocate(ccf%hrv_deadcrootc_storage_to_litr1c(beg:end))
    allocate(ccf%hrv_gresp_storage_to_litr1c(beg:end))     
    allocate(ccf%hrv_leafc_xfer_to_litr1c(beg:end))        
    allocate(ccf%hrv_frootc_xfer_to_litr1c(beg:end))       
    allocate(ccf%hrv_livestemc_xfer_to_litr1c(beg:end))    
    allocate(ccf%hrv_deadstemc_xfer_to_litr1c(beg:end))    
    allocate(ccf%hrv_livecrootc_xfer_to_litr1c(beg:end))   
    allocate(ccf%hrv_deadcrootc_xfer_to_litr1c(beg:end))   
    allocate(ccf%hrv_gresp_xfer_to_litr1c(beg:end))        
    allocate(ccf%m_litr1c_to_fire(beg:end))
    allocate(ccf%m_litr2c_to_fire(beg:end))
    allocate(ccf%m_litr3c_to_fire(beg:end))
    allocate(ccf%m_cwdc_to_fire(beg:end))
    allocate(ccf%leafc_to_litr1c(beg:end))
    allocate(ccf%leafc_to_litr2c(beg:end))
    allocate(ccf%leafc_to_litr3c(beg:end))
    allocate(ccf%frootc_to_litr1c(beg:end))
    allocate(ccf%frootc_to_litr2c(beg:end))
    allocate(ccf%frootc_to_litr3c(beg:end))
    allocate(ccf%cwdc_to_litr2c(beg:end))
    allocate(ccf%cwdc_to_litr3c(beg:end))
    allocate(ccf%litr1_hr(beg:end))
    allocate(ccf%litr1c_to_soil1c(beg:end))
    allocate(ccf%litr2_hr(beg:end))
    allocate(ccf%litr2c_to_soil2c(beg:end))
    allocate(ccf%litr3_hr(beg:end))
    allocate(ccf%litr3c_to_soil3c(beg:end))
    allocate(ccf%soil1_hr(beg:end))
    allocate(ccf%soil1c_to_soil2c(beg:end))
    allocate(ccf%soil2_hr(beg:end))
    allocate(ccf%soil2c_to_soil3c(beg:end))
    allocate(ccf%soil3_hr(beg:end))
    allocate(ccf%soil3c_to_soil4c(beg:end))
    allocate(ccf%soil4_hr(beg:end))
    allocate(ccf%dwt_seedc_to_leaf(beg:end))
    allocate(ccf%dwt_seedc_to_deadstem(beg:end))
    allocate(ccf%dwt_conv_cflux(beg:end))
    allocate(ccf%dwt_prod10c_gain(beg:end))
    allocate(ccf%dwt_prod100c_gain(beg:end))
    allocate(ccf%dwt_frootc_to_litr1c(beg:end))
    allocate(ccf%dwt_frootc_to_litr2c(beg:end))
    allocate(ccf%dwt_frootc_to_litr3c(beg:end))
    allocate(ccf%dwt_livecrootc_to_cwdc(beg:end))
    allocate(ccf%dwt_deadcrootc_to_cwdc(beg:end))
    allocate(ccf%dwt_closs(beg:end))
    allocate(ccf%landuseflux(beg:end))
    allocate(ccf%landuptake(beg:end))
    allocate(ccf%prod10c_loss(beg:end))
    allocate(ccf%prod100c_loss(beg:end))
    allocate(ccf%product_closs(beg:end))
    allocate(ccf%lithr(beg:end))
    allocate(ccf%somhr(beg:end))
    allocate(ccf%hr(beg:end))
    allocate(ccf%sr(beg:end))
    allocate(ccf%er(beg:end))
    allocate(ccf%litfire(beg:end))
    allocate(ccf%somfire(beg:end))
    allocate(ccf%totfire(beg:end))
    allocate(ccf%nep(beg:end))
    allocate(ccf%nbp(beg:end))
    allocate(ccf%nee(beg:end))
    allocate(ccf%col_cinputs(beg:end))
    allocate(ccf%col_coutputs(beg:end))
    allocate(ccf%col_fire_closs(beg:end))

    ccf%m_leafc_to_litr1c(beg:end)                = nan
    ccf%m_leafc_to_litr2c(beg:end)                = nan
    ccf%m_leafc_to_litr3c(beg:end)                = nan
    ccf%m_frootc_to_litr1c(beg:end)               = nan
    ccf%m_frootc_to_litr2c(beg:end)               = nan
    ccf%m_frootc_to_litr3c(beg:end)               = nan
    ccf%m_leafc_storage_to_litr1c(beg:end)        = nan
    ccf%m_frootc_storage_to_litr1c(beg:end)       = nan
    ccf%m_livestemc_storage_to_litr1c(beg:end)    = nan
    ccf%m_deadstemc_storage_to_litr1c(beg:end)    = nan
    ccf%m_livecrootc_storage_to_litr1c(beg:end)   = nan
    ccf%m_deadcrootc_storage_to_litr1c(beg:end)   = nan
    ccf%m_leafc_xfer_to_litr1c(beg:end)           = nan
    ccf%m_frootc_xfer_to_litr1c(beg:end)          = nan
    ccf%m_livestemc_xfer_to_litr1c(beg:end)       = nan
    ccf%m_deadstemc_xfer_to_litr1c(beg:end)       = nan
    ccf%m_livecrootc_xfer_to_litr1c(beg:end)      = nan
    ccf%m_deadcrootc_xfer_to_litr1c(beg:end)      = nan
    ccf%m_livestemc_to_cwdc(beg:end)              = nan
    ccf%m_deadstemc_to_cwdc(beg:end)              = nan
    ccf%m_livecrootc_to_cwdc(beg:end)             = nan
    ccf%m_deadcrootc_to_cwdc(beg:end)             = nan
    ccf%m_gresp_storage_to_litr1c(beg:end)        = nan
    ccf%m_gresp_xfer_to_litr1c(beg:end)           = nan
    ccf%m_deadstemc_to_cwdc_fire(beg:end)         = nan
    ccf%m_deadcrootc_to_cwdc_fire(beg:end)        = nan
    ccf%hrv_leafc_to_litr1c(beg:end)              = nan             
    ccf%hrv_leafc_to_litr2c(beg:end)              = nan             
    ccf%hrv_leafc_to_litr3c(beg:end)              = nan             
    ccf%hrv_frootc_to_litr1c(beg:end)             = nan            
    ccf%hrv_frootc_to_litr2c(beg:end)             = nan            
    ccf%hrv_frootc_to_litr3c(beg:end)             = nan            
    ccf%hrv_livestemc_to_cwdc(beg:end)            = nan           
    ccf%hrv_deadstemc_to_prod10c(beg:end)         = nan        
    ccf%hrv_deadstemc_to_prod100c(beg:end)        = nan       
    ccf%hrv_livecrootc_to_cwdc(beg:end)           = nan          
    ccf%hrv_deadcrootc_to_cwdc(beg:end)           = nan          
    ccf%hrv_leafc_storage_to_litr1c(beg:end)      = nan     
    ccf%hrv_frootc_storage_to_litr1c(beg:end)     = nan    
    ccf%hrv_livestemc_storage_to_litr1c(beg:end)  = nan 
    ccf%hrv_deadstemc_storage_to_litr1c(beg:end)  = nan 
    ccf%hrv_livecrootc_storage_to_litr1c(beg:end) = nan
    ccf%hrv_deadcrootc_storage_to_litr1c(beg:end) = nan
    ccf%hrv_gresp_storage_to_litr1c(beg:end)      = nan     
    ccf%hrv_leafc_xfer_to_litr1c(beg:end)         = nan        
    ccf%hrv_frootc_xfer_to_litr1c(beg:end)        = nan       
    ccf%hrv_livestemc_xfer_to_litr1c(beg:end)     = nan    
    ccf%hrv_deadstemc_xfer_to_litr1c(beg:end)     = nan    
    ccf%hrv_livecrootc_xfer_to_litr1c(beg:end)    = nan   
    ccf%hrv_deadcrootc_xfer_to_litr1c(beg:end)    = nan   
    ccf%hrv_gresp_xfer_to_litr1c(beg:end)         = nan        
    ccf%m_litr1c_to_fire(beg:end)                 = nan
    ccf%m_litr2c_to_fire(beg:end)                 = nan
    ccf%m_litr3c_to_fire(beg:end)                 = nan
    ccf%m_cwdc_to_fire(beg:end)                   = nan
    ccf%leafc_to_litr1c(beg:end)                  = nan
    ccf%leafc_to_litr2c(beg:end)                  = nan
    ccf%leafc_to_litr3c(beg:end)                  = nan
    ccf%frootc_to_litr1c(beg:end)                 = nan
    ccf%frootc_to_litr2c(beg:end)                 = nan
    ccf%frootc_to_litr3c(beg:end)                 = nan
    ccf%cwdc_to_litr2c(beg:end)                   = nan
    ccf%cwdc_to_litr3c(beg:end)                   = nan
    ccf%litr1_hr(beg:end)                         = nan
    ccf%litr1c_to_soil1c(beg:end)                 = nan
    ccf%litr2_hr(beg:end)                         = nan
    ccf%litr2c_to_soil2c(beg:end)                 = nan
    ccf%litr3_hr(beg:end)                         = nan
    ccf%litr3c_to_soil3c(beg:end)                 = nan
    ccf%soil1_hr(beg:end)                         = nan
    ccf%soil1c_to_soil2c(beg:end)                 = nan
    ccf%soil2_hr(beg:end)                         = nan
    ccf%soil2c_to_soil3c(beg:end)                 = nan
    ccf%soil3_hr(beg:end)                         = nan
    ccf%soil3c_to_soil4c(beg:end)                 = nan
    ccf%soil4_hr(beg:end)                         = nan
    ccf%dwt_seedc_to_leaf(beg:end)                = 0. !nan !gkw: disable dwt
    ccf%dwt_seedc_to_deadstem(beg:end)            = 0. !nan
    ccf%dwt_conv_cflux(beg:end)                   = 0. !nan
    ccf%dwt_prod10c_gain(beg:end)                 = 0. !nan
    ccf%dwt_prod100c_gain(beg:end)                = 0. !nan
    ccf%dwt_frootc_to_litr1c(beg:end)             = 0. !nan
    ccf%dwt_frootc_to_litr2c(beg:end)             = 0. !nan
    ccf%dwt_frootc_to_litr3c(beg:end)             = 0. !nan
    ccf%dwt_livecrootc_to_cwdc(beg:end)           = 0. !nan
    ccf%dwt_deadcrootc_to_cwdc(beg:end)           = 0. !nan
    ccf%dwt_closs(beg:end)                        = 0. !nan
    ccf%landuseflux(beg:end)                      = nan
    ccf%landuptake(beg:end)                       = nan
    ccf%prod10c_loss(beg:end)                     = nan
    ccf%prod100c_loss(beg:end)                    = nan
    ccf%product_closs(beg:end)                    = nan
    ccf%lithr(beg:end)                            = nan
    ccf%somhr(beg:end)                            = nan
    ccf%hr(beg:end)                               = nan
    ccf%sr(beg:end)                               = nan
    ccf%er(beg:end)                               = nan
    ccf%litfire(beg:end)                          = nan
    ccf%somfire(beg:end)                          = nan
    ccf%totfire(beg:end)                          = nan
    ccf%nep(beg:end)                              = nan
    ccf%nbp(beg:end)                              = nan
    ccf%nee(beg:end)                              = nan
    ccf%col_cinputs(beg:end)                      = nan
    ccf%col_coutputs(beg:end)                     = nan
    ccf%col_fire_closs(beg:end)                   = nan

  end subroutine init_column_cflux_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_column_nflux_type
!
! !INTERFACE:
  subroutine init_column_nflux_type(beg, end, cnf)
!
! !DESCRIPTION:
! Initialize column nitrogen flux variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (column_nflux_type), intent(inout):: cnf
!
! !REVISION HISTORY:
! Created by Peter Thornton
!
!EOP
!------------------------------------------------------------------------

    allocate(cnf%ndep_to_sminn(beg:end))
    allocate(cnf%nfix_to_sminn(beg:end))
    allocate(cnf%m_leafn_to_litr1n(beg:end))
    allocate(cnf%m_leafn_to_litr2n(beg:end))
    allocate(cnf%m_leafn_to_litr3n(beg:end))
    allocate(cnf%m_frootn_to_litr1n(beg:end))
    allocate(cnf%m_frootn_to_litr2n(beg:end))
    allocate(cnf%m_frootn_to_litr3n(beg:end))
    allocate(cnf%m_leafn_storage_to_litr1n(beg:end))
    allocate(cnf%m_frootn_storage_to_litr1n(beg:end))
    allocate(cnf%m_livestemn_storage_to_litr1n(beg:end))
    allocate(cnf%m_deadstemn_storage_to_litr1n(beg:end))
    allocate(cnf%m_livecrootn_storage_to_litr1n(beg:end))
    allocate(cnf%m_deadcrootn_storage_to_litr1n(beg:end))
    allocate(cnf%m_leafn_xfer_to_litr1n(beg:end))
    allocate(cnf%m_frootn_xfer_to_litr1n(beg:end))
    allocate(cnf%m_livestemn_xfer_to_litr1n(beg:end))
    allocate(cnf%m_deadstemn_xfer_to_litr1n(beg:end))
    allocate(cnf%m_livecrootn_xfer_to_litr1n(beg:end))
    allocate(cnf%m_deadcrootn_xfer_to_litr1n(beg:end))
    allocate(cnf%m_livestemn_to_cwdn(beg:end))
    allocate(cnf%m_deadstemn_to_cwdn(beg:end))
    allocate(cnf%m_livecrootn_to_cwdn(beg:end))
    allocate(cnf%m_deadcrootn_to_cwdn(beg:end))
    allocate(cnf%m_retransn_to_litr1n(beg:end))
    allocate(cnf%hrv_leafn_to_litr1n(beg:end))             
    allocate(cnf%hrv_leafn_to_litr2n(beg:end))             
    allocate(cnf%hrv_leafn_to_litr3n(beg:end))             
    allocate(cnf%hrv_frootn_to_litr1n(beg:end))            
    allocate(cnf%hrv_frootn_to_litr2n(beg:end))            
    allocate(cnf%hrv_frootn_to_litr3n(beg:end))            
    allocate(cnf%hrv_livestemn_to_cwdn(beg:end))           
    allocate(cnf%hrv_deadstemn_to_prod10n(beg:end))        
    allocate(cnf%hrv_deadstemn_to_prod100n(beg:end))       
    allocate(cnf%hrv_livecrootn_to_cwdn(beg:end))          
    allocate(cnf%hrv_deadcrootn_to_cwdn(beg:end))          
    allocate(cnf%hrv_retransn_to_litr1n(beg:end))          
    allocate(cnf%hrv_leafn_storage_to_litr1n(beg:end))     
    allocate(cnf%hrv_frootn_storage_to_litr1n(beg:end))    
    allocate(cnf%hrv_livestemn_storage_to_litr1n(beg:end)) 
    allocate(cnf%hrv_deadstemn_storage_to_litr1n(beg:end)) 
    allocate(cnf%hrv_livecrootn_storage_to_litr1n(beg:end))
    allocate(cnf%hrv_deadcrootn_storage_to_litr1n(beg:end))
    allocate(cnf%hrv_leafn_xfer_to_litr1n(beg:end))        
    allocate(cnf%hrv_frootn_xfer_to_litr1n(beg:end))       
    allocate(cnf%hrv_livestemn_xfer_to_litr1n(beg:end))    
    allocate(cnf%hrv_deadstemn_xfer_to_litr1n(beg:end))    
    allocate(cnf%hrv_livecrootn_xfer_to_litr1n(beg:end))   
    allocate(cnf%hrv_deadcrootn_xfer_to_litr1n(beg:end))   
    allocate(cnf%m_deadstemn_to_cwdn_fire(beg:end))
    allocate(cnf%m_deadcrootn_to_cwdn_fire(beg:end))
    allocate(cnf%m_litr1n_to_fire(beg:end))
    allocate(cnf%m_litr2n_to_fire(beg:end))
    allocate(cnf%m_litr3n_to_fire(beg:end))
    allocate(cnf%m_cwdn_to_fire(beg:end))
    allocate(cnf%leafn_to_litr1n(beg:end))
    allocate(cnf%leafn_to_litr2n(beg:end))
    allocate(cnf%leafn_to_litr3n(beg:end))
    allocate(cnf%frootn_to_litr1n(beg:end))
    allocate(cnf%frootn_to_litr2n(beg:end))
    allocate(cnf%frootn_to_litr3n(beg:end))
    allocate(cnf%cwdn_to_litr2n(beg:end))
    allocate(cnf%cwdn_to_litr3n(beg:end))
    allocate(cnf%litr1n_to_soil1n(beg:end))
    allocate(cnf%sminn_to_soil1n_l1(beg:end))
    allocate(cnf%litr2n_to_soil2n(beg:end))
    allocate(cnf%sminn_to_soil2n_l2(beg:end))
    allocate(cnf%litr3n_to_soil3n(beg:end))
    allocate(cnf%sminn_to_soil3n_l3(beg:end))
    allocate(cnf%soil1n_to_soil2n(beg:end))
    allocate(cnf%sminn_to_soil2n_s1(beg:end))
    allocate(cnf%soil2n_to_soil3n(beg:end))
    allocate(cnf%sminn_to_soil3n_s2(beg:end))
    allocate(cnf%soil3n_to_soil4n(beg:end))
    allocate(cnf%sminn_to_soil4n_s3(beg:end))
    allocate(cnf%soil4n_to_sminn(beg:end))
    allocate(cnf%sminn_to_denit_l1s1(beg:end))
    allocate(cnf%sminn_to_denit_l2s2(beg:end))
    allocate(cnf%sminn_to_denit_l3s3(beg:end))
    allocate(cnf%sminn_to_denit_s1s2(beg:end))
    allocate(cnf%sminn_to_denit_s2s3(beg:end))
    allocate(cnf%sminn_to_denit_s3s4(beg:end))
    allocate(cnf%sminn_to_denit_s4(beg:end))
    allocate(cnf%sminn_to_denit_excess(beg:end))
    allocate(cnf%sminn_leached(beg:end))
    allocate(cnf%dwt_seedn_to_leaf(beg:end))
    allocate(cnf%dwt_seedn_to_deadstem(beg:end))
    allocate(cnf%dwt_conv_nflux(beg:end))
    allocate(cnf%dwt_prod10n_gain(beg:end))
    allocate(cnf%dwt_prod100n_gain(beg:end))
    allocate(cnf%dwt_frootn_to_litr1n(beg:end))
    allocate(cnf%dwt_frootn_to_litr2n(beg:end))
    allocate(cnf%dwt_frootn_to_litr3n(beg:end))
    allocate(cnf%dwt_livecrootn_to_cwdn(beg:end))
    allocate(cnf%dwt_deadcrootn_to_cwdn(beg:end))
    allocate(cnf%dwt_nloss(beg:end))
    allocate(cnf%prod10n_loss(beg:end))
    allocate(cnf%prod100n_loss(beg:end))
    allocate(cnf%product_nloss(beg:end))
    allocate(cnf%potential_immob(beg:end))
    allocate(cnf%actual_immob(beg:end))
    allocate(cnf%sminn_to_plant(beg:end))
    allocate(cnf%supplement_to_sminn(beg:end))
    allocate(cnf%gross_nmin(beg:end))
    allocate(cnf%net_nmin(beg:end))
    allocate(cnf%denit(beg:end))
    allocate(cnf%col_ninputs(beg:end))
    allocate(cnf%col_noutputs(beg:end))
    allocate(cnf%col_fire_nloss(beg:end))

    cnf%ndep_to_sminn(beg:end) = nan
    cnf%nfix_to_sminn(beg:end) = nan
    cnf%m_leafn_to_litr1n(beg:end) = nan
    cnf%m_leafn_to_litr2n(beg:end) = nan
    cnf%m_leafn_to_litr3n(beg:end) = nan
    cnf%m_frootn_to_litr1n(beg:end) = nan
    cnf%m_frootn_to_litr2n(beg:end) = nan
    cnf%m_frootn_to_litr3n(beg:end) = nan
    cnf%m_leafn_storage_to_litr1n(beg:end) = nan
    cnf%m_frootn_storage_to_litr1n(beg:end) = nan
    cnf%m_livestemn_storage_to_litr1n(beg:end) = nan
    cnf%m_deadstemn_storage_to_litr1n(beg:end) = nan
    cnf%m_livecrootn_storage_to_litr1n(beg:end) = nan
    cnf%m_deadcrootn_storage_to_litr1n(beg:end) = nan
    cnf%m_leafn_xfer_to_litr1n(beg:end) = nan
    cnf%m_frootn_xfer_to_litr1n(beg:end) = nan
    cnf%m_livestemn_xfer_to_litr1n(beg:end) = nan
    cnf%m_deadstemn_xfer_to_litr1n(beg:end) = nan
    cnf%m_livecrootn_xfer_to_litr1n(beg:end) = nan
    cnf%m_deadcrootn_xfer_to_litr1n(beg:end) = nan
    cnf%m_livestemn_to_cwdn(beg:end) = nan
    cnf%m_deadstemn_to_cwdn(beg:end) = nan
    cnf%m_livecrootn_to_cwdn(beg:end) = nan
    cnf%m_deadcrootn_to_cwdn(beg:end) = nan
    cnf%m_retransn_to_litr1n(beg:end) = nan
    cnf%hrv_leafn_to_litr1n(beg:end) = nan             
    cnf%hrv_leafn_to_litr2n(beg:end) = nan             
    cnf%hrv_leafn_to_litr3n(beg:end) = nan             
    cnf%hrv_frootn_to_litr1n(beg:end) = nan            
    cnf%hrv_frootn_to_litr2n(beg:end) = nan            
    cnf%hrv_frootn_to_litr3n(beg:end) = nan            
    cnf%hrv_livestemn_to_cwdn(beg:end) = nan           
    cnf%hrv_deadstemn_to_prod10n(beg:end) = nan        
    cnf%hrv_deadstemn_to_prod100n(beg:end) = nan       
    cnf%hrv_livecrootn_to_cwdn(beg:end) = nan          
    cnf%hrv_deadcrootn_to_cwdn(beg:end) = nan          
    cnf%hrv_retransn_to_litr1n(beg:end) = nan          
    cnf%hrv_leafn_storage_to_litr1n(beg:end) = nan     
    cnf%hrv_frootn_storage_to_litr1n(beg:end) = nan    
    cnf%hrv_livestemn_storage_to_litr1n(beg:end) = nan 
    cnf%hrv_deadstemn_storage_to_litr1n(beg:end) = nan 
    cnf%hrv_livecrootn_storage_to_litr1n(beg:end) = nan
    cnf%hrv_deadcrootn_storage_to_litr1n(beg:end) = nan
    cnf%hrv_leafn_xfer_to_litr1n(beg:end) = nan        
    cnf%hrv_frootn_xfer_to_litr1n(beg:end) = nan       
    cnf%hrv_livestemn_xfer_to_litr1n(beg:end) = nan    
    cnf%hrv_deadstemn_xfer_to_litr1n(beg:end) = nan    
    cnf%hrv_livecrootn_xfer_to_litr1n(beg:end) = nan   
    cnf%hrv_deadcrootn_xfer_to_litr1n(beg:end) = nan   
    cnf%m_deadstemn_to_cwdn_fire(beg:end) = nan
    cnf%m_deadcrootn_to_cwdn_fire(beg:end) = nan
    cnf%m_litr1n_to_fire(beg:end) = nan
    cnf%m_litr2n_to_fire(beg:end) = nan
    cnf%m_litr3n_to_fire(beg:end) = nan
    cnf%m_cwdn_to_fire(beg:end) = nan
    cnf%leafn_to_litr1n(beg:end) = nan
    cnf%leafn_to_litr2n(beg:end) = nan
    cnf%leafn_to_litr3n(beg:end) = nan
    cnf%frootn_to_litr1n(beg:end) = nan
    cnf%frootn_to_litr2n(beg:end) = nan
    cnf%frootn_to_litr3n(beg:end) = nan
    cnf%cwdn_to_litr2n(beg:end) = nan
    cnf%cwdn_to_litr3n(beg:end) = nan
    cnf%litr1n_to_soil1n(beg:end) = nan
    cnf%sminn_to_soil1n_l1(beg:end) = nan
    cnf%litr2n_to_soil2n(beg:end) = nan
    cnf%sminn_to_soil2n_l2(beg:end) = nan
    cnf%litr3n_to_soil3n(beg:end) = nan
    cnf%sminn_to_soil3n_l3(beg:end) = nan
    cnf%soil1n_to_soil2n(beg:end) = nan
    cnf%sminn_to_soil2n_s1(beg:end) = nan
    cnf%soil2n_to_soil3n(beg:end) = nan
    cnf%sminn_to_soil3n_s2(beg:end) = nan
    cnf%soil3n_to_soil4n(beg:end) = nan
    cnf%sminn_to_soil4n_s3(beg:end) = nan
    cnf%soil4n_to_sminn(beg:end) = nan
    cnf%sminn_to_denit_l1s1(beg:end) = nan
    cnf%sminn_to_denit_l2s2(beg:end) = nan
    cnf%sminn_to_denit_l3s3(beg:end) = nan
    cnf%sminn_to_denit_s1s2(beg:end) = nan
    cnf%sminn_to_denit_s2s3(beg:end) = nan
    cnf%sminn_to_denit_s3s4(beg:end) = nan
    cnf%sminn_to_denit_s4(beg:end) = nan
    cnf%sminn_to_denit_excess(beg:end) = nan
    cnf%sminn_leached(beg:end) = nan
    cnf%dwt_seedn_to_leaf(beg:end) = 0. !nan
    cnf%dwt_seedn_to_deadstem(beg:end) = 0. !nan
    cnf%dwt_conv_nflux(beg:end) = 0. !nan
    cnf%dwt_prod10n_gain(beg:end) = 0. !nan
    cnf%dwt_prod100n_gain(beg:end) = 0. !nan
    cnf%dwt_frootn_to_litr1n(beg:end) = 0. !nan
    cnf%dwt_frootn_to_litr2n(beg:end) = 0. !nan
    cnf%dwt_frootn_to_litr3n(beg:end) = 0. !nan
    cnf%dwt_livecrootn_to_cwdn(beg:end) = 0. !nan
    cnf%dwt_deadcrootn_to_cwdn(beg:end) = 0. !nan
    cnf%dwt_nloss(beg:end) = 0. !nan
    cnf%prod10n_loss(beg:end) = nan
    cnf%prod100n_loss(beg:end) = nan
    cnf%product_nloss(beg:end) = nan
    cnf%potential_immob(beg:end) = nan
    cnf%actual_immob(beg:end) = nan
    cnf%sminn_to_plant(beg:end) = nan
    cnf%supplement_to_sminn(beg:end) = nan
    cnf%gross_nmin(beg:end) = nan
    cnf%net_nmin(beg:end) = nan
    cnf%denit(beg:end) = nan
    cnf%col_ninputs(beg:end) = nan
    cnf%col_noutputs(beg:end) = nan
    cnf%col_fire_nloss(beg:end) = nan

  end subroutine init_column_nflux_type

end module clmtypeInitMod
