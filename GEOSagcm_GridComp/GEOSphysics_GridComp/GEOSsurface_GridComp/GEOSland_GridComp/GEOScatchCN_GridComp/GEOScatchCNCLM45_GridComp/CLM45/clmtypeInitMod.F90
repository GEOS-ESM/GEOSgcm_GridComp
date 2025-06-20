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
  use MAPL_ConstantsMod, ONLY: r8 => MAPL_R4
! use shr_infnan_mod, only : nan => shr_infnan_nan, assignment(=)
  use nanMod      , only : nan, bigint
  use clmtype
  use clm_varpar  , only : maxpatch_pft, nlevsno, nlevgrnd, &
                           numpft, nlevsoi, nlevdecomp, nlevdecomp_full, &
                           ndecomp_cascade_transitions, ndecomp_pools, nlevcan
  use clm_varctl  , only : use_c13, use_c14

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
!!F. Li and S. Levis (11/06/12)
! !PRIVATE MEMBER FUNCTIONS:
  private :: init_pft_type
  private :: init_column_type
  private :: init_landunit_type
  private :: init_gridcell_type
  private :: init_pft_ecophys_constants
  private :: init_decomp_cascade_constants
#if (defined CNDV)
  private :: init_pft_DGVMecophys_constants
#endif
  private :: init_pft_pstate_type
  private :: init_pft_epv_type
#if (defined CNDV)
  private :: init_pft_pdgvstate_type
#endif
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
#ifdef LCH4
  private :: init_column_ch4_type
#endif
  private :: init_column_nflux_type
#ifdef LCH4
  private :: init_gridcell_ch4_type
#endif
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
!    *%area, *%wtlnd, *%wtxy, *%ixy, *%jxy, *%mxy, %snindex
!    *%ifspecial, *%ityplun, *%itype
!    *%pfti, *%pftf, *%pftn
!    *%coli, *%colf, *%coln
!    *%luni, *%lunf, *%lunn
!
! !USES:
!   use decompMod , only : get_proc_bounds, get_proc_global
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
    character(len=32), parameter :: subname = "initClmtype"
!------------------------------------------------------------------------

    ! Determine necessary indices

!   call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)
!   call get_proc_global(numg, numl, numc, nump)

    call init_pft_type     (begp, endp, pft)
    call init_column_type  (begc, endc, col)
    call init_landunit_type(begl, endl, lun)
    call init_gridcell_type(begg, endg, grc)

    ! pft ecophysiological constants

    call init_pft_ecophys_constants()

    call init_decomp_cascade_constants()

#if (defined CNDV)
    ! pft DGVM-specific ecophysiological constants

    call init_pft_DGVMecophys_constants()
#endif

    ! carbon balance structures (pft and column levels)

    call init_carbon_balance_type(begp, endp, pcbal)
    call init_carbon_balance_type(begc, endc, ccbal)

    ! nitrogen balance structures (pft and column levels)

    call init_nitrogen_balance_type(begp, endp, pnbal)
    call init_nitrogen_balance_type(begc, endc, cnbal)

    ! pft physical state variables at pft level 

    call init_pft_pstate_type(begp, endp, pps)

    ! pft ecophysiological variables (only at the pft level for now)
    call init_pft_epv_type(begp, endp, pepv)

    !pft photosynthesis relevant variables
    call init_pft_psynstate_type(begp, endp, ppsyns)
#if (defined CNDV)
    ! pft DGVM state variables at pft level and averaged to column

    call init_pft_pdgvstate_type(begp, endp, pdgvs)
#endif

    ! pft energy state variables at the pft level and averaged to the column

    call init_pft_estate_type(begp, endp, pes)

    ! pft carbon state variables at the pft level and averaged to the column

    call init_pft_cstate_type(begp, endp, pcs)
    call init_pft_cstate_type(begc, endc, pcs_a)
    
!    if ( use_c13 ) then       
!       call init_pft_cstate_type(begp, endp, pc13s)
!       call init_pft_cstate_type(begc, endc, pc13s_a)
!#ifdef CROP
!       stop 'initClmtype ERROR:: CROP and C13 can NOT be on at the same time'
!#endif
!    endif

!    if ( use_c14 ) then
!       call init_pft_cstate_type(begp, endp, pc14s)
!       call init_pft_cstate_type(begc, endc, pc14s_a)
!#ifdef CROP
!       stop 'initClmtype ERROR:: CROP and C14 can NOT be on at the same time'
!#endif
!    endif

    ! pft nitrogen state variables at the pft level and averaged to the column

    call init_pft_nstate_type(begp, endp, pns)
    call init_pft_nstate_type(begc, endc, pns_a)

    ! pft carbon flux variables at pft level and averaged to column

    call init_pft_cflux_type(begp, endp, pcf)
    call init_pft_cflux_type(begc, endc, pcf_a)
    
!    if ( use_c13 ) then       
!       call init_pft_cflux_type(begp, endp, pc13f)
!       call init_pft_cflux_type(begc, endc, pc13f_a)
!    endif
    
!    if ( use_c14 ) then
!       call init_pft_cflux_type(begp, endp, pc14f)
!       call init_pft_cflux_type(begc, endc, pc14f_a)
!    endif

    ! pft nitrogen flux variables at pft level and averaged to column

    call init_pft_nflux_type(begp, endp, pnf)
    call init_pft_nflux_type(begc, endc, pnf_a)

    ! column physical state variables at column level 

    call init_column_pstate_type(begc, endc, cps)

    ! column energy state variables at column level

    call init_column_estate_type(begc, endc, ces)

    ! column water state variables at column level 

    call init_column_wstate_type(begc, endc, cws)

    ! column carbon state variables at column level

    call init_column_cstate_type(begc, endc, ccs)
    
!    if ( use_c13 ) then       
!       call init_column_cstate_type(begc, endc, cc13s)
!    endif

!    if ( use_c14 ) then       
!       call init_column_cstate_type(begc, endc, cc14s)
!    endif

    ! column nitrogen state variables at column level

    call init_column_nstate_type(begc, endc, cns)

    ! column water flux variables at column level 

    call init_column_wflux_type(begc, endc, cwf)

    ! column carbon flux variables at column level

    call init_column_cflux_type(begc, endc, ccf)
    
!    if ( use_c13 ) then       
!       call init_column_cflux_type(begc, endc, cc13f)
!    endif
    
!    if ( use_c14 ) then       
!       call init_column_cflux_type(begc, endc, cc14f)
!    endif

#if (defined LCH4)
    ! column CH4 flux variables at column level
    call init_column_ch4_type(begc, endc, cch4)
#endif

    ! column nitrogen flux variables at column level

    call init_column_nflux_type(begc, endc, cnf)

#if (defined CNDV)
    ! gridcell DGVM variables

    call init_gridcell_dgvstate_type(begg, endg, gdgvs)
#endif

#if (defined LCH4)
    ! gridcell: ch4 variables

    call init_gridcell_ch4_type(begg, endg, gch4)
#endif

  end subroutine initClmtype

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_type
!
! !INTERFACE:
  subroutine init_pft_type (beg, end, pft)
!
! !DESCRIPTION:
! Initialize components of pft_type structure
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type(pft_type), intent(inout):: pft
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(pft%gridcell(beg:end),pft%wtgcell(beg:end))
    allocate(pft%landunit(beg:end))
    allocate(pft%column  (beg:end),pft%wtcol  (beg:end))

    allocate(pft%itype(beg:end))
    allocate(pft%active(beg:end))  

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

   allocate(col%gridcell(beg:end),col%wtgcell(beg:end))
   allocate(col%landunit(beg:end),col%wtlunit(beg:end))

   allocate(col%pfti(beg:end),col%pftf(beg:end),col%npfts(beg:end))

   allocate(col%itype(beg:end))
   allocate(col%active(beg:end))

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

   allocate(lun%itype(beg:end))
   allocate(lun%ifspecial(beg:end))

  end subroutine init_landunit_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_gridcell_type
!
! !INTERFACE:
  subroutine init_gridcell_type (beg, end,grc)
!
! !DESCRIPTION:
! Initialize components of gridcell_type structure
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type(gridcell_type), intent(inout):: grc
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

   ! Modified by fzeng, 2017
   allocate(grc%gindex(beg:end))
   allocate(grc%forc_ndep(beg:end))   
   allocate(grc%forc_rh  (beg:end))
   allocate(grc%forc_wind(beg:end))
   allocate(grc%forc_t   (beg:end))
   allocate(grc%forc_rain(beg:end))
   allocate(grc%forc_snow(beg:end))   
   allocate(grc%latdeg(beg:end))
   allocate(grc%londeg(beg:end))
   allocate(grc%forc_hdm(beg:end))
   allocate(grc%forc_lnfm(beg:end))

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
    allocate(pftcon%xl(0:numpft))
    allocate(pftcon%rhol(0:numpft))
    allocate(pftcon%rhos(0:numpft))
    allocate(pftcon%taul(0:numpft))
    allocate(pftcon%taus(0:numpft))
    allocate(pftcon%z0mr(0:numpft))
    allocate(pftcon%displar(0:numpft))
    allocate(pftcon%roota_par(0:numpft))
    allocate(pftcon%rootb_par(0:numpft))
    allocate(pftcon%slatop(0:numpft))
    allocate(pftcon%dsladlai(0:numpft))
    allocate(pftcon%leafcn(0:numpft))
    allocate(pftcon%flnr(0:numpft))
    allocate(pftcon%woody(0:numpft))
    allocate(pftcon%lflitcn(0:numpft))
    allocate(pftcon%frootcn(0:numpft))
    allocate(pftcon%livewdcn(0:numpft))
    allocate(pftcon%deadwdcn(0:numpft))
    allocate(pftcon%graincn(0:numpft))
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
    allocate(pftcon%leaf_long(0:numpft))
    allocate(pftcon%evergreen(0:numpft))
    allocate(pftcon%stress_decid(0:numpft))
    allocate(pftcon%season_decid(0:numpft))
    allocate(pftcon%dwood(0:numpft))
    allocate(pftcon%cc_dstem(0:numpft))
    allocate(pftcon%cc_leaf(0:numpft)) 
    allocate(pftcon%cc_lstem(0:numpft))
    allocate(pftcon%cc_other(0:numpft))
    allocate(pftcon%fm_dstem(0:numpft))
    allocate(pftcon%fm_leaf(0:numpft)) 
    allocate(pftcon%fm_lstem(0:numpft))
    allocate(pftcon%fm_other(0:numpft))
    allocate(pftcon%fm_root(0:numpft)) 
    allocate(pftcon%fm_lroot(0:numpft))
    allocate(pftcon%fm_droot(0:numpft))
    allocate(pftcon%rootprof_beta(0:numpft))
    allocate(pftcon%fertnitro(0:numpft))
    allocate(pftcon%fleafcn(0:numpft))
    allocate(pftcon%ffrootcn(0:numpft))
    allocate(pftcon%fstemcn(0:numpft))
    allocate(pftcon%laimx(0:numpft)) 
    allocate(pftcon%ztopmx(0:numpft))

    ! fzeng:    
    allocate(pftcon%declfact(0:numpft))
    allocate(pftcon%bfact(0:numpft))
    allocate(pftcon%aleaff(0:numpft))
    allocate(pftcon%arootf(0:numpft))
    allocate(pftcon%astemf(0:numpft))
    allocate(pftcon%arooti(0:numpft))
    allocate(pftcon%fleafi(0:numpft))
    allocate(pftcon%allconsl(0:numpft))
    allocate(pftcon%allconss(0:numpft))
    allocate(pftcon%grperc(0:numpft))
    allocate(pftcon%grpnow(0:numpft))
    allocate(pftcon%fsr_pft(0:numpft))
    allocate(pftcon%fd_pft(0:numpft))
    allocate(pftcon%mnNHplantdate(0:numpft))
    allocate(pftcon%mxNHplantdate(0:numpft))
    allocate(pftcon%mnSHplantdate(0:numpft))
    allocate(pftcon%mxSHplantdate(0:numpft)) 
    allocate(pftcon%gddmin(0:numpft))
    allocate(pftcon%hybgdd(0:numpft))
    allocate(pftcon%lfemerg(0:numpft))
    allocate(pftcon%grnfill(0:numpft))
    allocate(pftcon%mxmat(0:numpft))
    allocate(pftcon%minplanttemp(0:numpft))
    allocate(pftcon%planttemp(0:numpft))
    allocate(pftcon%mxtmp(0:numpft))
    allocate(pftcon%baset(0:numpft))
    allocate(pftcon%qe25(0:numpft))    

    pftcon%noveg(:) = huge(1)
    pftcon%tree(:) = huge(1)
    pftcon%fnitr(:) = nan
    pftcon%c3psn(:) = nan
    pftcon%xl(:) = nan
    pftcon%rhol(:) = nan
    pftcon%rhos(:) = nan
    pftcon%taul(:) = nan
    pftcon%taus(:) = nan
    pftcon%z0mr(:) = nan
    pftcon%displar(:) = nan
    pftcon%roota_par(:) = nan
    pftcon%rootb_par(:) = nan
    pftcon%slatop(:) = nan
    pftcon%dsladlai(:) = nan
    pftcon%leafcn(:) = nan
    pftcon%flnr(:) = nan
    pftcon%woody(:) = nan
    pftcon%lflitcn(:) = nan
    pftcon%frootcn(:) = nan
    pftcon%livewdcn(:) = nan
    pftcon%deadwdcn(:) = nan
    pftcon%graincn(:) = nan
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
    pftcon%leaf_long(:) = nan
    pftcon%evergreen(:) = nan
    pftcon%stress_decid(:) = nan
    pftcon%season_decid(:) = nan
    pftcon%dwood(:) = nan
    pftcon%cc_dstem(:) = nan
    pftcon%cc_leaf(:) = nan
    pftcon%cc_lstem(:) = nan
    pftcon%cc_other(:) = nan
    pftcon%fm_dstem(:) = nan
    pftcon%fm_leaf(:) = nan
    pftcon%fm_lstem(:) = nan
    pftcon%fm_other(:) = nan
    pftcon%fm_root(:) = nan
    pftcon%fm_lroot(:) = nan
    pftcon%fm_droot(:) = nan
    pftcon%rootprof_beta(:) = nan
    pftcon%fertnitro(:) = nan
    pftcon%fleafcn(:)   = nan
    pftcon%ffrootcn(:)  = nan
    pftcon%fstemcn(:)   = nan
    pftcon%laimx(:)   = nan
    pftcon%ztopmx(:)   = nan
 
    ! fzeng:    
    pftcon%declfact(:) = nan
    pftcon%bfact(:)    = nan
    pftcon%aleaff(:)   = nan
    pftcon%arootf(:)   = nan
    pftcon%astemf(:)   = nan
    pftcon%arooti(:)   = nan
    pftcon%fleafi(:)   = nan
    pftcon%allconsl(:) = nan
    pftcon%allconss(:) = nan
    pftcon%grperc(:)   = nan
    pftcon%grpnow(:)   = nan
    pftcon%fsr_pft(:)  = nan
    pftcon%fd_pft(:)   = nan
    pftcon%mnNHplantdate(:)   = huge(1)
    pftcon%mxNHplantdate(:)   = huge(1)
    pftcon%mnSHplantdate(:)   = huge(1)
    pftcon%mxSHplantdate(:)   = huge(1)
    pftcon%gddmin(:)   = nan
    pftcon%hybgdd(:)   = nan
    pftcon%lfemerg(:)  = nan
    pftcon%grnfill(:)  = nan
    pftcon%mxmat(:)    = nan
    pftcon%minplanttemp(:)    = nan
    pftcon%planttemp(:)= nan
    pftcon%mxtmp(:)    = nan
    pftcon%baset(:)    = nan
    pftcon%qe25(:)     = nan
    
  end subroutine init_pft_ecophys_constants

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_decomp_cascade_constants
!
! !INTERFACE:
  subroutine init_decomp_cascade_constants()
!
! !DESCRIPTION:
! Initialize decomposition cascade state
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Created by Charlie Koven
!
!EOP
!------------------------------------------------------------------------

    !-- properties of each pathway along decomposition cascade 
    allocate(decomp_cascade_con%cascade_step_name(1:ndecomp_cascade_transitions))
    allocate(decomp_cascade_con%cascade_donor_pool(1:ndecomp_cascade_transitions))
    allocate(decomp_cascade_con%cascade_receiver_pool(1:ndecomp_cascade_transitions))
    !-- properties of each decomposing pool
    allocate(decomp_cascade_con%floating_cn_ratio_decomp_pools(0:ndecomp_pools))
    allocate(decomp_cascade_con%decomp_pool_name_restart(0:ndecomp_pools))
    allocate(decomp_cascade_con%decomp_pool_name_history(0:ndecomp_pools))
    allocate(decomp_cascade_con%decomp_pool_name_long(0:ndecomp_pools))
    allocate(decomp_cascade_con%decomp_pool_name_short(0:ndecomp_pools))
    allocate(decomp_cascade_con%is_litter(0:ndecomp_pools))
    allocate(decomp_cascade_con%is_soil(0:ndecomp_pools))
    allocate(decomp_cascade_con%is_cwd(0:ndecomp_pools))
    allocate(decomp_cascade_con%initial_cn_ratio(0:ndecomp_pools))
    allocate(decomp_cascade_con%initial_stock(0:ndecomp_pools))
    allocate(decomp_cascade_con%is_metabolic(0:ndecomp_pools))
    allocate(decomp_cascade_con%is_cellulose(0:ndecomp_pools))
    allocate(decomp_cascade_con%is_lignin(0:ndecomp_pools))
    allocate(decomp_cascade_con%spinup_factor(0:ndecomp_pools))
    !-- properties of each pathway along decomposition cascade 
    decomp_cascade_con%cascade_step_name(1:ndecomp_cascade_transitions) = ''
    decomp_cascade_con%cascade_donor_pool(1:ndecomp_cascade_transitions) = 0
    decomp_cascade_con%cascade_receiver_pool(1:ndecomp_cascade_transitions) = 0
    !-- properties of each decomposing pool
    decomp_cascade_con%floating_cn_ratio_decomp_pools(0:ndecomp_pools) = .false.
    decomp_cascade_con%decomp_pool_name_history(0:ndecomp_pools) = ''
    decomp_cascade_con%decomp_pool_name_restart(0:ndecomp_pools) = ''
    decomp_cascade_con%decomp_pool_name_long(0:ndecomp_pools) = ''
    decomp_cascade_con%decomp_pool_name_short(0:ndecomp_pools) = ''
    decomp_cascade_con%is_litter(0:ndecomp_pools) = .false.
    decomp_cascade_con%is_soil(0:ndecomp_pools) = .false.
    decomp_cascade_con%is_cwd(0:ndecomp_pools) = .false.
    decomp_cascade_con%initial_cn_ratio(0:ndecomp_pools) = nan
    decomp_cascade_con%initial_stock(0:ndecomp_pools) = nan
    decomp_cascade_con%is_metabolic(0:ndecomp_pools) = .false.
    decomp_cascade_con%is_cellulose(0:ndecomp_pools) = .false.
    decomp_cascade_con%is_lignin(0:ndecomp_pools) = .false.
    decomp_cascade_con%spinup_factor(0:ndecomp_pools) = nan

  end subroutine init_decomp_cascade_constants


#if (defined CNDV)
!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_DGVMecophys_constants
!
! !INTERFACE:
  subroutine init_pft_DGVMecophys_constants()
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

    allocate(dgv_pftcon%crownarea_max(0:numpft))
    allocate(dgv_pftcon%tcmin(0:numpft))
    allocate(dgv_pftcon%tcmax(0:numpft))
    allocate(dgv_pftcon%gddmin(0:numpft))
    allocate(dgv_pftcon%twmax(0:numpft))
    allocate(dgv_pftcon%reinickerp(0:numpft))
    allocate(dgv_pftcon%allom1(0:numpft))
    allocate(dgv_pftcon%allom2(0:numpft))
    allocate(dgv_pftcon%allom3(0:numpft))

    dgv_pftcon%crownarea_max(:) = nan
    dgv_pftcon%tcmin(:) = nan
    dgv_pftcon%tcmax(:) = nan
    dgv_pftcon%gddmin(:) = nan
    dgv_pftcon%twmax(:) = nan
    dgv_pftcon%reinickerp(:) = nan
    dgv_pftcon%allom1(:) = nan
    dgv_pftcon%allom2(:) = nan
    dgv_pftcon%allom3(:) = nan

  end subroutine init_pft_DGVMecophys_constants
#endif

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
    use clm_varcon, only : spval,ispval
    use clm_varctl , only : crop_prog
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
    allocate(pps%prec10(beg:end)) !F. Li and S. Levis
    allocate(pps%prec60(beg:end)) !F. Li and S. Levis
    allocate(pps%frac_veg_nosno(beg:end))
    allocate(pps%frac_veg_nosno_alb(beg:end))
    allocate(pps%rootfr(beg:end,1:nlevgrnd))
    allocate(pps%rootr(beg:end,1:nlevgrnd))
    allocate(pps%laisun(beg:end))
    allocate(pps%laisha(beg:end))
    allocate(pps%btran2(beg:end))   ! F. Li and S. Levis
    allocate(pps%tlai(beg:end))
    allocate(pps%tsai(beg:end))
    allocate(pps%elai(beg:end))
    allocate(pps%esai(beg:end))
    allocate(pps%htop(beg:end))
    allocate(pps%hbot(beg:end))
    allocate(pps%burndate(beg:end))   ! F. Li and S. Levis
    if ( crop_prog )then
       allocate(pps%hdidx(beg:end))
       allocate(pps%cumvd(beg:end))
       allocate(pps%htmx(beg:end))
       allocate(pps%vf(beg:end))
       allocate(pps%gddmaturity(beg:end))
       allocate(pps%gdd0(beg:end))
       allocate(pps%gdd8(beg:end))
       allocate(pps%gdd10(beg:end))
       allocate(pps%gdd020(beg:end))
       allocate(pps%gdd820(beg:end))
       allocate(pps%gdd1020(beg:end))
       allocate(pps%gddplant(beg:end))
       allocate(pps%gddtsoi(beg:end))
       allocate(pps%huileaf(beg:end))
       allocate(pps%huigrain(beg:end))
       allocate(pps%aleafi(beg:end))
       allocate(pps%astemi(beg:end))
       allocate(pps%aleaf(beg:end))
       allocate(pps%astem(beg:end))
       allocate(pps%croplive(beg:end))
       allocate(pps%cropplant(beg:end)) !,numpft)) ! make 2-D if using
       allocate(pps%harvdate(beg:end))  !,numpft)) ! crop rotation
       allocate(pps%idop(beg:end))
       allocate(pps%peaklai(beg:end))
    end if
    allocate(pps%forc_hgt_u_pft(beg:end))
    allocate(pps%lfpftd(beg:end))      !F. Li and S. Levis

    ! 4/14/05: PET
    ! Adding isotope code
    
!    if ( use_c13 ) then       
!       allocate(pps%alphapsnsun(beg:end))
!       allocate(pps%alphapsnsha(beg:end))
!    endif
    
#if (defined LCH4)
    ! CH4 code
    allocate(pps%grnd_ch4_cond(beg:end))
    allocate(pps%canopy_cond(beg:end))
#endif
   ! and vertical profiles for calculating fluxes
    allocate(pps%leaf_prof(beg:end,1:nlevdecomp_full))
    allocate(pps%froot_prof(beg:end,1:nlevdecomp_full))
    allocate(pps%croot_prof(beg:end,1:nlevdecomp_full))
    allocate(pps%stem_prof(beg:end,1:nlevdecomp_full))
    pps%prec10(beg:end) = nan   ! F. Li and S. Levis
    pps%prec60(beg:end) = nan   ! F. Li and S. Levis
    pps%frac_veg_nosno(beg:end) = huge(1)
    pps%frac_veg_nosno_alb(beg:end) = 0
    pps%rootfr(beg:end,:nlevgrnd) = spval
    pps%rootr (beg:end,:nlevgrnd) = spval
    pps%laisun(beg:end) = nan
    pps%laisha(beg:end) = nan
    pps%btran2(beg:end) = spval       !F. Li and S. Levis
    pps%tlai(beg:end) = 0._r8
    pps%tsai(beg:end) = 0._r8
    pps%elai(beg:end) = 0._r8
    pps%esai(beg:end) = 0._r8
    pps%htop(beg:end) = 0._r8
    pps%hbot(beg:end) = 0._r8
    pps%burndate(beg:end)    = ispval   ! F. Li and S. Levis
    if ( crop_prog )then
       pps%hdidx(beg:end)       = nan
       pps%cumvd(beg:end)       = nan
       pps%htmx(beg:end)        = 0.0_r8
       pps%vf(beg:end)          = 0.0_r8
       pps%gddmaturity(beg:end) = spval
       pps%gdd0(beg:end)        = spval
       pps%gdd8(beg:end)        = spval
       pps%gdd10(beg:end)       = spval
       pps%gdd020(beg:end)      = spval
       pps%gdd820(beg:end)      = spval
       pps%gdd1020(beg:end)     = spval
       pps%gddplant(beg:end)    = spval
       pps%gddtsoi(beg:end)     = spval
       pps%huileaf(beg:end)     = nan
       pps%huigrain(beg:end)    = nan
       pps%aleafi(beg:end)      = nan
       pps%astemi(beg:end)      = nan
       pps%aleaf(beg:end)       = nan
       pps%astem(beg:end)       = nan
       pps%croplive(beg:end)    = .false.
       pps%cropplant(beg:end)   = .false.
       pps%harvdate(beg:end)    = huge(1)
       pps%idop(beg:end)        = huge(1)
       pps%peaklai(beg:end)     = 0
    end if
    pps%forc_hgt_u_pft(beg:end) = nan

    ! 4/14/05: PET
    ! Adding isotope code    ! EBK Check this!
    !!!pps%cisun(beg:end) = spval
    !!!pps%cisha(beg:end) = spval
    
!    if ( use_c13 ) then       
!       pps%alphapsnsun(beg:end) = spval
!       pps%alphapsnsha(beg:end) = spval
!    endif

#if defined (LCH4)
    ! CH4 code
    pps%grnd_ch4_cond(beg:end) = nan
    pps%canopy_cond(beg:end) = nan
#endif
   ! and vertical profiles for calculating fluxes
    pps%leaf_prof(beg:end,1:nlevdecomp_full) = spval
    pps%froot_prof(beg:end,1:nlevdecomp_full) = spval
    pps%croot_prof(beg:end,1:nlevdecomp_full) = spval
    pps%stem_prof(beg:end,1:nlevdecomp_full) = spval

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
! !USES:
    use clm_varcon, only : spval
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
    allocate(pepv%fert_counter(beg:end))
    allocate(pepv%grain_flag(beg:end))
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
!    if ( use_c13 ) then
!       allocate(pepv%xsmrpool_c13ratio(beg:end))
!    endif
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
#if (defined CNDV)
    allocate(pepv%tempsum_litfall(beg:end))
    allocate(pepv%annsum_litfall(beg:end))
#endif
!    if ( use_c13 ) then
!       allocate(pepv%rc13_canair(beg:end))
!       allocate(pepv%rc13_psnsun(beg:end))
!       allocate(pepv%rc13_psnsha(beg:end))
!    endif
    
!    if ( use_c14 ) then
!       allocate(pepv%rc14_atm(beg:end))
!    endif

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
    pepv%fert_counter(beg:end) = nan
    pepv%grain_flag(beg:end) = nan
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
!    if ( use_c13 ) then
!       pepv%xsmrpool_c13ratio(beg:end) = nan
!    endif
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
#if (defined CNDV)
    pepv%tempsum_litfall(beg:end) = nan
    pepv%annsum_litfall(beg:end) = nan
#endif
!   if ( use_c13 ) then
!      pepv%rc13_canair(beg:end) = spval
!      pepv%rc13_psnsun(beg:end) = spval
!      pepv%rc13_psnsha(beg:end) = spval
!   endif

!   if ( use_c14 ) then
!      pepv%rc14_atm(beg:end) = nan
!      ! pepv%rc14_canair(beg:end) = nan
!      ! pepv%rc14_psnsun(beg:end) = nan
!      ! pepv%rc14_psnsha(beg:end) = nan
!   endif
    
  end subroutine init_pft_epv_type

#if (defined CNDV)
!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_pdgvstate_type
!
! !INTERFACE:
  subroutine init_pft_pdgvstate_type(beg, end, pdgvs)
!
! !DESCRIPTION:
! Initialize pft DGVM state variables
!
! !USES:
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_dgvstate_type), intent(inout):: pdgvs
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(pdgvs%agddtw(beg:end))
    allocate(pdgvs%agdd(beg:end))
    allocate(pdgvs%t_mo(beg:end))
    allocate(pdgvs%t_mo_min(beg:end))
    allocate(pdgvs%prec365(beg:end))
    allocate(pdgvs%present(beg:end))
    allocate(pdgvs%pftmayexist(beg:end))
    allocate(pdgvs%nind(beg:end))
    allocate(pdgvs%lm_ind(beg:end))
    allocate(pdgvs%lai_ind(beg:end))
    allocate(pdgvs%fpcinc(beg:end))
    allocate(pdgvs%fpcgrid(beg:end))
    allocate(pdgvs%fpcgridold(beg:end))
    allocate(pdgvs%crownarea(beg:end))
    allocate(pdgvs%greffic(beg:end))
    allocate(pdgvs%heatstress(beg:end))

    pdgvs%agddtw(beg:end)           = nan
    pdgvs%agdd(beg:end)             = nan
    pdgvs%t_mo(beg:end)             = nan
    pdgvs%t_mo_min(beg:end)         = nan
    pdgvs%prec365(beg:end)          = nan
    pdgvs%present(beg:end)          = .false.
    pdgvs%pftmayexist(beg:end)      = .true.
    pdgvs%nind(beg:end)             = nan
    pdgvs%lm_ind(beg:end)           = nan
    pdgvs%lai_ind(beg:end)          = nan
    pdgvs%fpcinc(beg:end)           = nan
    pdgvs%fpcgrid(beg:end)          = nan
    pdgvs%fpcgridold(beg:end)       = nan
    pdgvs%crownarea(beg:end)        = nan
    pdgvs%greffic(beg:end)          = nan
    pdgvs%heatstress(beg:end)       = nan

  end subroutine init_pft_pdgvstate_type
#endif

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_psynstate_type
!
! !INTERFACE:
  subroutine init_pft_psynstate_type(beg, end, ppsyns)
!
! !DESCRIPTION:
! Initialize pft energy state
!
! !USES:
    use clm_varcon, only : spval 
    use clm_varctl, only : crop_prog
! !AGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_psynstate_type), intent(inout):: ppsyns
!
! !REVISION HISTORY:
! Created by Jinyun Tang
!
!EOP
!-----------------------------------------------------------------------

   allocate(ppsyns%c3flag(beg:end))
   allocate(ppsyns%ac(beg:end,1:nlevcan))
   allocate(ppsyns%aj(beg:end,1:nlevcan))
   allocate(ppsyns%ap(beg:end,1:nlevcan))
   allocate(ppsyns%ag(beg:end,1:nlevcan))
   allocate(ppsyns%an(beg:end,1:nlevcan))
   allocate(ppsyns%vcmax_z(beg:end,1:nlevcan))
   allocate(ppsyns%cp(beg:end))
   allocate(ppsyns%kc(beg:end))
   allocate(ppsyns%ko(beg:end))
   allocate(ppsyns%qe(beg:end))
   allocate(ppsyns%tpu_z(beg:end,1:nlevcan))
   allocate(ppsyns%kp_z(beg:end,1:nlevcan))   
   allocate(ppsyns%theta_cj(beg:end))
   allocate(ppsyns%bbb(beg:end))
   allocate(ppsyns%mbb(beg:end))
   allocate(ppsyns%gb_mol(beg:end))
   allocate(ppsyns%gs_mol(beg:end,1:nlevcan))

   ppsyns%c3flag = .false.
   ppsyns%ac(beg:end,1:nlevcan) = nan
   ppsyns%aj(beg:end,1:nlevcan) = nan
   ppsyns%ap(beg:end,1:nlevcan) = nan
   ppsyns%ag(beg:end,1:nlevcan) = nan
   ppsyns%an(beg:end,1:nlevcan) = nan
   ppsyns%vcmax_z(beg:end,1:nlevcan) = nan
   ppsyns%cp(beg:end) = nan
   ppsyns%kc(beg:end) = nan
   ppsyns%ko(beg:end) = nan
   ppsyns%qe(beg:end) = nan
   ppsyns%tpu_z(beg:end,1:nlevcan) = nan
   ppsyns%kp_z(beg:end,1:nlevcan) = nan
   ppsyns%theta_cj(beg:end) = nan
   ppsyns%bbb(beg:end) = nan
   ppsyns%mbb(beg:end) = nan
   ppsyns%gb_mol(beg:end) = nan
   ppsyns%gs_mol(beg:end,1:nlevcan) = nan

  end subroutine init_pft_psynstate_type


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
! !USES:
    use clm_varcon, only : spval
    use clm_varctl, only : crop_prog
! !AGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_estate_type), intent(inout):: pes
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------

    allocate(pes%t_ref2m(beg:end))
    allocate(pes%t_ref2m_min(beg:end))
    allocate(pes%t_ref2m_max(beg:end))
    allocate(pes%t10(beg:end))
    if ( crop_prog )then
       allocate(pes%a10tmin(beg:end))
       allocate(pes%a5tmin(beg:end))
    end if

    pes%t_ref2m(beg:end) = nan
    pes%t_ref2m_min(beg:end) = nan
    pes%t_ref2m_max(beg:end) = nan
    pes%t10(beg:end)                = spval
    if ( crop_prog )then
       pes%a10tmin(beg:end)     = spval
       pes%a5tmin(beg:end)      = spval
    end if

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
! !USES:
    use clm_varctl, only : crop_prog
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
    if ( crop_prog )then
       allocate(pcs%grainc(beg:end))
       allocate(pcs%grainc_storage(beg:end))
       allocate(pcs%grainc_xfer(beg:end))
    end if
!#ifdef CN
    allocate(pcs%woodc(beg:end))
!#endif

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
    if ( crop_prog )then
       pcs%grainc(beg:end)         = nan
       pcs%grainc_storage(beg:end) = nan
       pcs%grainc_xfer(beg:end)    = nan
    end if
!#ifdef CN
    pcs%woodc(beg:end) = nan
!#endif

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
! !USES:
    use clm_varctl, only : crop_prog
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

    if ( crop_prog )then
       allocate(pns%grainn(beg:end))
       allocate(pns%grainn_storage(beg:end))
       allocate(pns%grainn_xfer(beg:end))
    end if
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

    if ( crop_prog )then
       pns%grainn(beg:end)         = nan
       pns%grainn_storage(beg:end) = nan
       pns%grainn_xfer(beg:end)    = nan
    end if
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
! !USES:
    use clm_varcon, only : spval
    use clm_varctl , only : crop_prog
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
    allocate(pcf%psnsun_z(beg:end,1:nlevcan))
    allocate(pcf%psnsha_z(beg:end,1:nlevcan))
    allocate(pcf%cisun_z(beg:end,1:nlevcan))
    allocate(pcf%cisha_z(beg:end,1:nlevcan))
    allocate(pcf%lmrsun(beg:end))
    allocate(pcf%lmrsha(beg:end))
    allocate(pcf%lmrsun_z(beg:end,1:nlevcan))
    allocate(pcf%lmrsha_z(beg:end,1:nlevcan))
    allocate(pcf%fpsn(beg:end))
    allocate(pcf%fco2(beg:end))
    allocate(pcf%psnsun_wc(beg:end))
    allocate(pcf%psnsha_wc(beg:end))
    allocate(pcf%fpsn_wc(beg:end))
    allocate(pcf%psnsun_wj(beg:end))
    allocate(pcf%psnsha_wj(beg:end))
    allocate(pcf%fpsn_wj(beg:end))
    allocate(pcf%psnsun_wp(beg:end))
    allocate(pcf%psnsha_wp(beg:end))
    allocate(pcf%fpsn_wp(beg:end))

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
             
    ! fire related variables changed by F. Li and S. Levis           
    allocate(pcf%m_leafc_to_fire(beg:end))
    allocate(pcf%m_leafc_storage_to_fire(beg:end))
    allocate(pcf%m_leafc_xfer_to_fire(beg:end))
    allocate(pcf%m_livestemc_to_fire(beg:end))
    allocate(pcf%m_livestemc_storage_to_fire(beg:end))
    allocate(pcf%m_livestemc_xfer_to_fire(beg:end))
    allocate(pcf%m_deadstemc_to_fire(beg:end))
    allocate(pcf%m_deadstemc_storage_to_fire(beg:end))
    allocate(pcf%m_deadstemc_xfer_to_fire(beg:end))
    allocate(pcf%m_frootc_to_fire(beg:end))
    allocate(pcf%m_frootc_storage_to_fire(beg:end))
    allocate(pcf%m_frootc_xfer_to_fire(beg:end))
    allocate(pcf%m_livecrootc_to_fire(beg:end))
    allocate(pcf%m_livecrootc_storage_to_fire(beg:end))
    allocate(pcf%m_livecrootc_xfer_to_fire(beg:end))
    allocate(pcf%m_deadcrootc_to_fire(beg:end))
    allocate(pcf%m_deadcrootc_storage_to_fire(beg:end))
    allocate(pcf%m_deadcrootc_xfer_to_fire(beg:end))
    allocate(pcf%m_gresp_storage_to_fire(beg:end))
    allocate(pcf%m_gresp_xfer_to_fire(beg:end))
    allocate(pcf%m_leafc_to_litter_fire(beg:end))
    allocate(pcf%m_leafc_storage_to_litter_fire(beg:end))
    allocate(pcf%m_leafc_xfer_to_litter_fire(beg:end))
    allocate(pcf%m_livestemc_to_litter_fire(beg:end))
    allocate(pcf%m_livestemc_storage_to_litter_fire(beg:end))
    allocate(pcf%m_livestemc_xfer_to_litter_fire(beg:end))
    allocate(pcf%m_livestemc_to_deadstemc_fire(beg:end))
    allocate(pcf%m_deadstemc_to_litter_fire(beg:end))
    allocate(pcf%m_deadstemc_storage_to_litter_fire(beg:end))
    allocate(pcf%m_deadstemc_xfer_to_litter_fire(beg:end))
    allocate(pcf%m_frootc_to_litter_fire(beg:end))
    allocate(pcf%m_frootc_storage_to_litter_fire(beg:end))
    allocate(pcf%m_frootc_xfer_to_litter_fire(beg:end))
    allocate(pcf%m_livecrootc_to_litter_fire(beg:end))
    allocate(pcf%m_livecrootc_storage_to_litter_fire(beg:end))
    allocate(pcf%m_livecrootc_xfer_to_litter_fire(beg:end))
    allocate(pcf%m_livecrootc_to_deadcrootc_fire(beg:end))
    allocate(pcf%m_deadcrootc_to_litter_fire(beg:end))
    allocate(pcf%m_deadcrootc_storage_to_litter_fire(beg:end))
    allocate(pcf%m_deadcrootc_xfer_to_litter_fire(beg:end))
    allocate(pcf%m_gresp_storage_to_litter_fire(beg:end))
    allocate(pcf%m_gresp_xfer_to_litter_fire(beg:end))


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
    allocate(pcf%grain_mr(beg:end))
    allocate(pcf%leaf_curmr(beg:end))
    allocate(pcf%froot_curmr(beg:end))
    allocate(pcf%livestem_curmr(beg:end))
    allocate(pcf%livecroot_curmr(beg:end))
    allocate(pcf%grain_curmr(beg:end))
    allocate(pcf%leaf_xsmr(beg:end))
    allocate(pcf%froot_xsmr(beg:end))
    allocate(pcf%livestem_xsmr(beg:end))
    allocate(pcf%livecroot_xsmr(beg:end))
    allocate(pcf%grain_xsmr(beg:end))
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
    if ( crop_prog )then
       allocate(pcf%xsmrpool_to_atm(beg:end))
       allocate(pcf%grainc_xfer_to_grainc(beg:end))
       allocate(pcf%livestemc_to_litter(beg:end))
       allocate(pcf%grainc_to_food(beg:end))
       allocate(pcf%cpool_to_grainc(beg:end))
       allocate(pcf%cpool_to_grainc_storage(beg:end))
       allocate(pcf%cpool_grain_gr(beg:end))
       allocate(pcf%cpool_grain_storage_gr(beg:end))
       allocate(pcf%transfer_grain_gr(beg:end))
       allocate(pcf%grainc_storage_to_xfer(beg:end))
    end if
!#ifdef CN
    allocate(pcf%frootc_alloc(beg:end))
    allocate(pcf%frootc_loss(beg:end))
    allocate(pcf%leafc_alloc(beg:end))
    allocate(pcf%leafc_loss(beg:end))
    allocate(pcf%woodc_alloc(beg:end))
    allocate(pcf%woodc_loss(beg:end))
!#endif
#ifdef LCH4
    allocate(pcf%tempavg_agnpp(beg:end))
    allocate(pcf%tempavg_bgnpp(beg:end))
    allocate(pcf%annavg_agnpp(beg:end))
    allocate(pcf%annavg_bgnpp(beg:end))
#endif

    pcf%psnsun(beg:end) = nan
    pcf%psnsha(beg:end) = nan
    pcf%psnsun_z(beg:end,:nlevcan) = nan
    pcf%psnsha_z(beg:end,:nlevcan) = nan
    pcf%cisun_z(beg:end,:nlevcan) = nan
    pcf%cisha_z(beg:end,:nlevcan) = nan
    pcf%lmrsun(beg:end) = nan
    pcf%lmrsha(beg:end) = nan
    pcf%lmrsun_z(beg:end,:nlevcan) = nan
    pcf%lmrsha_z(beg:end,:nlevcan) = nan
    pcf%fpsn(beg:end) = spval
    pcf%fco2(beg:end) = 0._r8
    pcf%psnsun_wc(beg:end) = nan
    pcf%psnsha_wc(beg:end) = nan
    pcf%fpsn_wc(beg:end) = nan
    pcf%psnsun_wj(beg:end) = nan
    pcf%psnsha_wj(beg:end) = nan
    pcf%fpsn_wj(beg:end) = nan
    pcf%psnsun_wp(beg:end) = nan
    pcf%psnsha_wp(beg:end) = nan
    pcf%fpsn_wp(beg:end) = nan

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
     
    ! fire variable changed by F. Li and S. Levis
    pcf%m_leafc_to_fire(beg:end) = nan
    pcf%m_leafc_storage_to_fire(beg:end) = nan
    pcf%m_leafc_xfer_to_fire(beg:end) = nan
    pcf%m_livestemc_to_fire(beg:end) = nan
    pcf%m_livestemc_storage_to_fire(beg:end) = nan
    pcf%m_livestemc_xfer_to_fire(beg:end) = nan
    pcf%m_deadstemc_to_fire(beg:end) = nan
    pcf%m_deadstemc_storage_to_fire(beg:end) = nan
    pcf%m_deadstemc_xfer_to_fire(beg:end) = nan
    pcf%m_frootc_to_fire(beg:end) = nan
    pcf%m_frootc_storage_to_fire(beg:end) = nan
    pcf%m_frootc_xfer_to_fire(beg:end) = nan
    pcf%m_livecrootc_to_fire(beg:end) = nan
    pcf%m_livecrootc_storage_to_fire(beg:end) = nan
    pcf%m_livecrootc_xfer_to_fire(beg:end) = nan
    pcf%m_deadcrootc_to_fire(beg:end) = nan
    pcf%m_deadcrootc_storage_to_fire(beg:end) = nan
    pcf%m_deadcrootc_xfer_to_fire(beg:end) = nan
    pcf%m_gresp_storage_to_fire(beg:end) = nan
    pcf%m_gresp_xfer_to_fire(beg:end) = nan
    
    pcf%m_leafc_to_litter_fire(beg:end) = nan
    pcf%m_leafc_storage_to_litter_fire(beg:end) = nan
    pcf%m_leafc_xfer_to_litter_fire(beg:end) = nan
    pcf%m_livestemc_to_litter_fire(beg:end) = nan
    pcf%m_livestemc_storage_to_litter_fire(beg:end) = nan
    pcf%m_livestemc_xfer_to_litter_fire(beg:end) = nan
    pcf%m_livestemc_to_deadstemc_fire(beg:end) = nan
    pcf%m_deadstemc_to_litter_fire(beg:end) = nan
    pcf%m_deadstemc_storage_to_litter_fire(beg:end) = nan
    pcf%m_deadstemc_xfer_to_litter_fire(beg:end) = nan
    pcf%m_frootc_to_litter_fire(beg:end) = nan
    pcf%m_frootc_storage_to_litter_fire(beg:end) = nan
    pcf%m_frootc_xfer_to_litter_fire(beg:end) = nan
    pcf%m_livecrootc_to_litter_fire(beg:end) = nan
    pcf%m_livecrootc_storage_to_litter_fire(beg:end) = nan
    pcf%m_livecrootc_xfer_to_litter_fire(beg:end) = nan
    pcf%m_livecrootc_to_deadcrootc_fire(beg:end) = nan
    pcf%m_deadcrootc_to_litter_fire(beg:end) = nan
    pcf%m_deadcrootc_storage_to_litter_fire(beg:end) = nan
    pcf%m_deadcrootc_xfer_to_litter_fire(beg:end) = nan
    pcf%m_gresp_storage_to_litter_fire(beg:end) = nan
    pcf%m_gresp_xfer_to_litter_fire(beg:end) = nan


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
    pcf%grain_mr(beg:end) = nan
    pcf%leaf_curmr(beg:end) = nan
    pcf%froot_curmr(beg:end) = nan
    pcf%livestem_curmr(beg:end) = nan
    pcf%livecroot_curmr(beg:end) = nan
    pcf%grain_curmr(beg:end) = nan
    pcf%leaf_xsmr(beg:end) = nan
    pcf%froot_xsmr(beg:end) = nan
    pcf%livestem_xsmr(beg:end) = nan
    pcf%livecroot_xsmr(beg:end) = nan
    pcf%grain_xsmr(beg:end) = nan
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
    if ( crop_prog )then
       pcf%xsmrpool_to_atm(beg:end)         = nan
       pcf%grainc_xfer_to_grainc(beg:end)   = nan
       pcf%livestemc_to_litter(beg:end)     = nan
       pcf%grainc_to_food(beg:end)          = nan
       pcf%cpool_to_grainc(beg:end)         = nan
       pcf%cpool_to_grainc_storage(beg:end) = nan
       pcf%cpool_grain_gr(beg:end)          = nan
       pcf%cpool_grain_storage_gr(beg:end)  = nan
       pcf%transfer_grain_gr(beg:end)       = nan
       pcf%grainc_storage_to_xfer(beg:end)  = nan
    end if
!#if (defined CN)
    pcf%frootc_alloc(beg:end) = nan
    pcf%frootc_loss(beg:end) = nan
    pcf%leafc_alloc(beg:end) = nan
    pcf%leafc_loss(beg:end) = nan
    pcf%woodc_alloc(beg:end) = nan
    pcf%woodc_loss(beg:end) = nan
!#endif
#if (defined LCH4)
    pcf%tempavg_agnpp(beg:end) = spval ! For back-compatibility
    pcf%tempavg_bgnpp(beg:end) = spval ! For back-compatibility
    pcf%annavg_agnpp(beg:end) = spval ! To detect first year
    pcf%annavg_bgnpp(beg:end) = spval ! To detect first year
#endif


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
! !USES:
    use clm_varctl , only : crop_prog
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
           
    ! fire variables changed by F. Li and S. Levis            
   allocate(pnf%m_leafn_to_fire(beg:end))
    allocate(pnf%m_leafn_storage_to_fire(beg:end))
    allocate(pnf%m_leafn_xfer_to_fire(beg:end))
    allocate(pnf%m_livestemn_to_fire(beg:end))
    allocate(pnf%m_livestemn_storage_to_fire(beg:end))
    allocate(pnf%m_livestemn_xfer_to_fire(beg:end))
    allocate(pnf%m_deadstemn_to_fire(beg:end))
    allocate(pnf%m_deadstemn_storage_to_fire(beg:end))
    allocate(pnf%m_deadstemn_xfer_to_fire(beg:end))
    allocate(pnf%m_frootn_to_fire(beg:end))
    allocate(pnf%m_frootn_storage_to_fire(beg:end))
    allocate(pnf%m_frootn_xfer_to_fire(beg:end))
    allocate(pnf%m_livecrootn_to_fire(beg:end))
    allocate(pnf%m_livecrootn_storage_to_fire(beg:end))
    allocate(pnf%m_livecrootn_xfer_to_fire(beg:end))
    allocate(pnf%m_deadcrootn_to_fire(beg:end))
    allocate(pnf%m_deadcrootn_storage_to_fire(beg:end))
    allocate(pnf%m_deadcrootn_xfer_to_fire(beg:end))
    allocate(pnf%m_retransn_to_fire(beg:end))

    allocate(pnf%m_leafn_to_litter_fire(beg:end))
    allocate(pnf%m_leafn_storage_to_litter_fire(beg:end))
    allocate(pnf%m_leafn_xfer_to_litter_fire(beg:end))
    allocate(pnf%m_livestemn_to_litter_fire(beg:end))
    allocate(pnf%m_livestemn_storage_to_litter_fire(beg:end))
    allocate(pnf%m_livestemn_xfer_to_litter_fire(beg:end))
    allocate(pnf%m_livestemn_to_deadstemn_fire(beg:end))
    allocate(pnf%m_deadstemn_to_litter_fire(beg:end))
    allocate(pnf%m_deadstemn_storage_to_litter_fire(beg:end))
    allocate(pnf%m_deadstemn_xfer_to_litter_fire(beg:end))
    allocate(pnf%m_frootn_to_litter_fire(beg:end))
    allocate(pnf%m_frootn_storage_to_litter_fire(beg:end))
    allocate(pnf%m_frootn_xfer_to_litter_fire(beg:end))
    allocate(pnf%m_livecrootn_to_litter_fire(beg:end))
    allocate(pnf%m_livecrootn_storage_to_litter_fire(beg:end))
    allocate(pnf%m_livecrootn_xfer_to_litter_fire(beg:end))
    allocate(pnf%m_livecrootn_to_deadcrootn_fire(beg:end))
    allocate(pnf%m_deadcrootn_to_litter_fire(beg:end))
    allocate(pnf%m_deadcrootn_storage_to_litter_fire(beg:end))
    allocate(pnf%m_deadcrootn_xfer_to_litter_fire(beg:end))
    allocate(pnf%m_retransn_to_litter_fire(beg:end))



    allocate(pnf%leafn_xfer_to_leafn(beg:end))
    allocate(pnf%frootn_xfer_to_frootn(beg:end))
    allocate(pnf%livestemn_xfer_to_livestemn(beg:end))
    allocate(pnf%deadstemn_xfer_to_deadstemn(beg:end))
    allocate(pnf%livecrootn_xfer_to_livecrootn(beg:end))
    allocate(pnf%deadcrootn_xfer_to_deadcrootn(beg:end))
    allocate(pnf%leafn_to_litter(beg:end))
    allocate(pnf%leafn_to_retransn(beg:end))
    allocate(pnf%frootn_to_retransn(beg:end))
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
    if ( crop_prog )then
       allocate(pnf%grainn_xfer_to_grainn(beg:end))
       allocate(pnf%livestemn_to_litter(beg:end))
       allocate(pnf%grainn_to_food(beg:end))
       allocate(pnf%npool_to_grainn(beg:end))
       allocate(pnf%npool_to_grainn_storage(beg:end))
       allocate(pnf%grainn_storage_to_xfer(beg:end))
       allocate(pnf%fert(beg:end))
       allocate(pnf%soyfixn(beg:end))
    end if

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
        
    ! fire varibles changed by F. Li and S. Levis         
   pnf%m_leafn_to_fire(beg:end) = nan
    pnf%m_leafn_storage_to_fire(beg:end) = nan
    pnf%m_leafn_xfer_to_fire(beg:end) = nan
    pnf%m_livestemn_to_fire(beg:end) = nan
    pnf%m_livestemn_storage_to_fire(beg:end) = nan
    pnf%m_livestemn_xfer_to_fire(beg:end) = nan
    pnf%m_deadstemn_to_fire(beg:end) = nan
    pnf%m_deadstemn_storage_to_fire(beg:end) = nan
    pnf%m_deadstemn_xfer_to_fire(beg:end) = nan
    pnf%m_frootn_to_fire(beg:end) = nan
    pnf%m_frootn_storage_to_fire(beg:end) = nan
    pnf%m_frootn_xfer_to_fire(beg:end) = nan
    pnf%m_livestemn_to_fire(beg:end) = nan
    pnf%m_livecrootn_storage_to_fire(beg:end) = nan
    pnf%m_livecrootn_xfer_to_fire(beg:end) = nan
    pnf%m_deadcrootn_to_fire(beg:end) = nan
    pnf%m_deadcrootn_storage_to_fire(beg:end) = nan
    pnf%m_deadcrootn_xfer_to_fire(beg:end) = nan
    pnf%m_retransn_to_fire(beg:end) = nan

    pnf%m_leafn_to_litter_fire(beg:end) = nan
    pnf%m_leafn_storage_to_litter_fire(beg:end) = nan
    pnf%m_leafn_xfer_to_litter_fire(beg:end) = nan
    pnf%m_livestemn_to_litter_fire(beg:end) = nan
    pnf%m_livestemn_storage_to_litter_fire(beg:end) = nan
    pnf%m_livestemn_xfer_to_litter_fire(beg:end) = nan
    pnf%m_livestemn_to_deadstemn_fire(beg:end) = nan
    pnf%m_deadstemn_to_litter_fire(beg:end) = nan
    pnf%m_deadstemn_storage_to_litter_fire(beg:end) = nan
    pnf%m_deadstemn_xfer_to_litter_fire(beg:end) = nan
    pnf%m_frootn_to_litter_fire(beg:end) = nan
    pnf%m_frootn_storage_to_litter_fire(beg:end) = nan
    pnf%m_frootn_xfer_to_litter_fire(beg:end) = nan
    pnf%m_livecrootn_to_litter_fire(beg:end) = nan
    pnf%m_livecrootn_storage_to_litter_fire(beg:end) = nan
    pnf%m_livecrootn_xfer_to_litter_fire(beg:end) = nan
    pnf%m_livecrootn_to_deadcrootn_fire(beg:end) = nan
    pnf%m_deadcrootn_to_litter_fire(beg:end) = nan
    pnf%m_deadcrootn_storage_to_litter_fire(beg:end) = nan
    pnf%m_deadcrootn_xfer_to_litter_fire(beg:end) = nan
    pnf%m_retransn_to_litter_fire(beg:end) = nan



    pnf%leafn_xfer_to_leafn(beg:end) = nan
    pnf%frootn_xfer_to_frootn(beg:end) = nan
    pnf%livestemn_xfer_to_livestemn(beg:end) = nan
    pnf%deadstemn_xfer_to_deadstemn(beg:end) = nan
    pnf%livecrootn_xfer_to_livecrootn(beg:end) = nan
    pnf%deadcrootn_xfer_to_deadcrootn(beg:end) = nan
    pnf%leafn_to_litter(beg:end) = nan
    pnf%leafn_to_retransn(beg:end) = nan
    pnf%frootn_to_retransn(beg:end) = nan
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
    if ( crop_prog )then
       pnf%grainn_xfer_to_grainn(beg:end)   = nan
       pnf%livestemn_to_litter(beg:end)     = nan
       pnf%grainn_to_food(beg:end)          = nan
       pnf%npool_to_grainn(beg:end)         = nan
       pnf%npool_to_grainn_storage(beg:end) = nan
       pnf%grainn_storage_to_xfer(beg:end)  = nan
       pnf%fert(beg:end)                    = nan
       pnf%soyfixn(beg:end)                 = nan
    end if

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

    allocate(cps%snl(beg:end))      !* cannot be averaged up
   
   !F. Li and S. Levis
    allocate(cps%gdp_lf(beg:end))  
    allocate(cps%peatf_lf(beg:end))  
    allocate(cps%abm_lf(beg:end))  
    allocate(cps%lgdp_col(beg:end))  
    allocate(cps%lgdp1_col(beg:end))
    allocate(cps%lpop_col(beg:end))  

    allocate(cps%bsw(beg:end,nlevgrnd))
    allocate(cps%watsat(beg:end,nlevgrnd))
    allocate(cps%watfc(beg:end,nlevgrnd))
    allocate(cps%sucsat(beg:end,nlevgrnd))
    allocate(cps%snow_depth(beg:end))   
    allocate(cps%snowdp(beg:end))
    allocate(cps%zi(beg:end,-nlevsno+0:nlevgrnd))
    allocate(cps%dz(beg:end,-nlevsno+1:nlevgrnd))
    allocate(cps%z (beg:end,-nlevsno+1:nlevgrnd))
    allocate(cps%wf(beg:end))
    allocate(cps%wf2(beg:end))
    allocate(cps%psisat(beg:end,-nlevsno+1:nlevgrnd))       ! added by fzeng
    allocate(cps%psiwilt(beg:end))                          ! added by fzeng
    allocate(cws%fsat(beg:end))                             ! added by fzeng
    allocate(cps%soilpsi(beg:end,nlevgrnd))
    allocate(cps%bd(beg:end,nlevgrnd))
    allocate(cps%fpi(beg:end))
    allocate(cps%fpi_vr(beg:end,1:nlevdecomp_full))
    allocate(cps%fpg(beg:end))
    allocate(cps%annsum_counter(beg:end))
    allocate(cps%cannsum_npp(beg:end))
    allocate(cps%col_lag_npp(beg:end))
    allocate(cps%cannavg_t2m(beg:end))
    

    ! fire-related variables changed by F. Li and S. Levis
    allocate(cps%nfire(beg:end))
    allocate(cps%farea_burned(beg:end))
    allocate(cps%fsr_col(beg:end))
    allocate(cps%fd_col(beg:end))
    allocate(cps%cropf_col(beg:end))
    allocate(cps%prec10_col(beg:end))
    allocate(cps%prec60_col(beg:end))
    allocate(cps%lfc(beg:end))
    allocate(cps%lfc2(beg:end))
    allocate(cps%trotr1_col(beg:end))
    allocate(cps%trotr2_col(beg:end))
    allocate(cps%dtrotr_col(beg:end))
    allocate(cps%baf_crop(beg:end))
    allocate(cps%baf_peatf(beg:end))
    allocate(cps%fbac(beg:end))
    allocate(cps%fbac1(beg:end))
    allocate(cps%btran_col(beg:end))
    allocate(cps%wtlf(beg:end))
    allocate(cps%lfwt(beg:end))
#ifdef LCH4
! New variable for finundated parameterization
    allocate(cps%zwt0(beg:end))
    allocate(cps%f0(beg:end))
    allocate(cps%p3(beg:end))
! New variable for methane
    allocate(cps%pH(beg:end))
#endif

    allocate(cps%rf_decomp_cascade(beg:end,1:nlevdecomp_full,1:ndecomp_cascade_transitions))
    allocate(cps%pathfrac_decomp_cascade(beg:end,1:nlevdecomp_full,1:ndecomp_cascade_transitions))
    allocate(cps%nfixation_prof(beg:end,1:nlevdecomp_full))
    allocate(cps%ndep_prof(beg:end,1:nlevdecomp_full))
    allocate(cps%alt(beg:end))
    allocate(cps%altmax(beg:end))
    allocate(cps%altmax_lastyear(beg:end))
    allocate(cps%alt_indx(beg:end))
    allocate(cps%altmax_indx(beg:end))
    allocate(cps%altmax_lastyear_indx(beg:end))
    allocate(cps%som_adv_coef(beg:end,1:nlevdecomp_full))
    allocate(cps%som_diffus_coef(beg:end,1:nlevdecomp_full))
  
    cps%snl(beg:end) = huge(1)

    !F. Li and S. Levis
    cps%gdp_lf(beg:end) = nan
    cps%peatf_lf(beg:end) = nan
    cps%abm_lf(beg:end) = 13 
    cps%lgdp_col(beg:end) = nan
    cps%lgdp1_col(beg:end) = nan
    cps%lpop_col(beg:end) = nan

    cps%bsw(beg:end,1:nlevgrnd) = nan
    cps%watsat(beg:end,1:nlevgrnd) = nan
    cps%watfc(beg:end,1:nlevgrnd) = nan
    cps%sucsat(beg:end,1:nlevgrnd) = nan
    cps%snow_depth(beg:end) = nan
    cps%snowdp(beg:end) = nan
    cps%zi(beg:end,-nlevsno+0:nlevgrnd) = nan
    cps%dz(beg:end,-nlevsno+1:nlevgrnd) = nan
    cps%z (beg:end,-nlevsno+1:nlevgrnd) = nan
    cps%wf(beg:end) = nan
    cps%wf2(beg:end) = nan
    cps%psisat(beg:end,-nlevsno+1:nlevgrnd) = nan                ! added by fzeng
    cps%psiwilt(beg:end) = nan                                   ! added by fzeng
    cws%fsat(beg:end) = nan                                      ! added by fzeng
    cps%soilpsi(beg:end,1:nlevgrnd) = spval
    cps%bd(beg:end,1:nlevgrnd) = spval
    cps%fpi(beg:end) = nan
    cps%fpi_vr(beg:end,1:nlevdecomp_full) = nan
    cps%fpg(beg:end) = nan
    cps%annsum_counter(beg:end) = nan
    cps%cannsum_npp(beg:end) = nan
    cps%col_lag_npp(beg:end) = spval
    cps%cannavg_t2m(beg:end) = nan

    ! fire-related varibles changed by F. Li and S. Levis
    cps%nfire(beg:end) = spval
    cps%farea_burned(beg:end) = nan
    cps%btran_col(beg:end) = nan
    cps%wtlf(beg:end) = nan
    cps%lfwt(beg:end) = nan
    cps%fsr_col(beg:end) = nan
    cps%fd_col(beg:end) = nan
    cps%cropf_col(beg:end) = nan
    cps%baf_crop(beg:end) = nan
    cps%baf_peatf(beg:end) = nan
    cps%fbac(beg:end) = nan
    cps%fbac1(beg:end) = nan
    cps%trotr1_col(beg:end) = 0._r8
    cps%trotr2_col(beg:end) = 0._r8
    cps%dtrotr_col(beg:end) = 0._r8
    cps%prec10_col(beg:end) = nan
    cps%prec60_col(beg:end) = nan
    cps%lfc(beg:end) = spval
    cps%lfc2(beg:end) = 0._r8
#ifdef LCH4
    cps%zwt0(beg:end) = nan
    cps%f0(beg:end)   = nan
    cps%p3(beg:end)   = nan
! New variable for methane
    cps%pH(beg:end)   = nan
#endif

    cps%rf_decomp_cascade(beg:end,1:nlevdecomp_full,1:ndecomp_cascade_transitions) = nan
    cps%pathfrac_decomp_cascade(beg:end,1:nlevdecomp_full,1:ndecomp_cascade_transitions) = nan
    cps%nfixation_prof(beg:end,1:nlevdecomp_full) = spval
    cps%ndep_prof(beg:end,1:nlevdecomp_full) = spval
    cps%alt(beg:end) = spval
    cps%altmax(beg:end) = spval
    cps%altmax_lastyear(beg:end) = spval
    cps%alt_indx(beg:end) = huge(1)
    cps%altmax_indx(beg:end) = huge(1)
    cps%altmax_lastyear_indx(beg:end) = huge(1)
    cps%som_adv_coef(beg:end,1:nlevdecomp_full) = spval
    cps%som_diffus_coef(beg:end,1:nlevdecomp_full) = spval

    allocate(cps%frac_h2osfc(beg:end))
    cps%frac_h2osfc(beg:end) = spval

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
    allocate(ces%tsoi17(beg:end))

    ces%t_grnd(beg:end)    = nan
    ces%t_soisno(beg:end,-nlevsno+1:nlevgrnd) = spval
    ces%tsoi17(beg:end) = spval

    allocate(ces%t_h2osfc(beg:end))

    ces%t_h2osfc(beg:end) = spval

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

    allocate(cws%h2osoi_liq(beg:end,-nlevsno+1:nlevgrnd))
!New variable for methane code
#ifdef LCH4
    allocate(cws%finundated(beg:end))
#endif

    cws%h2osoi_liq(beg:end,-nlevsno+1:nlevgrnd)= spval
#ifdef LCH4
    cws%finundated(beg:end) = nan
#endif

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
! !USES:
    use clm_varcon, only : spval
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
    allocate(ccs%col_ctrunc(beg:end))
    allocate(ccs%decomp_cpools_vr(beg:end,1:nlevdecomp_full,1:ndecomp_pools))
    allocate(ccs%decomp_cpools(beg:end,1:ndecomp_pools))
    allocate(ccs%decomp_cpools_1m(beg:end,1:ndecomp_pools))
    allocate(ccs%col_ctrunc_vr(beg:end,1:nlevdecomp_full))
    allocate(ccs%seedc(beg:end))
    allocate(ccs%prod10c(beg:end))
    allocate(ccs%prod100c(beg:end))
    allocate(ccs%totprodc(beg:end))
    allocate(ccs%totlitc(beg:end))
    allocate(ccs%totsomc(beg:end))
    allocate(ccs%totlitc_1m(beg:end))
    allocate(ccs%totsomc_1m(beg:end))
    allocate(ccs%totecosysc(beg:end))
    allocate(ccs%totcolc(beg:end))

    !F. Li and S. Levis
    allocate(ccs%rootc_col(beg:end))
    allocate(ccs%totvegc_col(beg:end))
    allocate(ccs%leafc_col(beg:end))
    allocate(ccs%fuelc(beg:end))
    allocate(ccs%fuelc_crop(beg:end))
    allocate(ccs%cpool_col(beg:end))

    ccs%cwdc(beg:end) = nan
    ccs%decomp_cpools_vr(beg:end,1:nlevdecomp_full,1:ndecomp_pools) = nan
    ccs%decomp_cpools(beg:end,1:ndecomp_pools) = nan
    ccs%decomp_cpools_1m(beg:end,1:ndecomp_pools) = nan
    ccs%col_ctrunc(beg:end) = nan
    ccs%col_ctrunc_vr(beg:end,1:nlevdecomp_full) = nan
    ccs%seedc(beg:end) = nan
    ccs%prod10c(beg:end) = nan
    ccs%prod100c(beg:end) = nan
    ccs%totprodc(beg:end) = nan
    ccs%totlitc(beg:end) = nan
    ccs%totsomc(beg:end) = nan
    ccs%totlitc_1m(beg:end) = nan
    ccs%totsomc_1m(beg:end) = nan
    ccs%totecosysc(beg:end) = nan
    ccs%totcolc(beg:end) = nan

    ccs%rootc_col(beg:end) = nan
    ccs%totvegc_col(beg:end) = nan
    ccs%leafc_col(beg:end) = nan
    ccs%fuelc(beg:end) = spval
    ccs%fuelc_crop(beg:end) = nan
    ccs%cpool_col(beg:end) = nan

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
! !USES:
    use clm_varcon, only : spval
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


    allocate(cns%decomp_npools(beg:end,1:ndecomp_pools))
    allocate(cns%decomp_npools_1m(beg:end,1:ndecomp_pools))
    allocate(cns%decomp_npools_vr(beg:end,1:nlevdecomp_full,1:ndecomp_pools))
    allocate(cns%sminn_vr(beg:end,1:nlevdecomp_full))
    allocate(cns%col_ntrunc_vr(beg:end,1:nlevdecomp_full))
#ifdef NITRIF_DENITRIF
    allocate(cns%smin_no3_vr(beg:end,1:nlevdecomp_full))
    allocate(cns%smin_nh4_vr(beg:end,1:nlevdecomp_full))
    allocate(cns%smin_no3(beg:end))
    allocate(cns%smin_nh4(beg:end))
#endif
    allocate(cns%cwdn(beg:end))
    allocate(cns%sminn(beg:end))
    allocate(cns%col_ntrunc(beg:end))
    allocate(cns%seedn(beg:end))
    allocate(cns%prod10n(beg:end))
    allocate(cns%prod100n(beg:end))
    allocate(cns%totprodn(beg:end))
    allocate(cns%totlitn(beg:end))
    allocate(cns%totsomn(beg:end))
    allocate(cns%totlitn_1m(beg:end))
    allocate(cns%totsomn_1m(beg:end))
    allocate(cns%totecosysn(beg:end))
    allocate(cns%totcoln(beg:end))

    cns%decomp_npools(beg:end,1:ndecomp_pools) = nan
    cns%decomp_npools_1m(beg:end,1:ndecomp_pools) = nan
    cns%decomp_npools_vr(beg:end,1:nlevdecomp_full,1:ndecomp_pools) = nan
    cns%sminn_vr(beg:end,1:nlevdecomp_full) = nan
    cns%col_ntrunc_vr(beg:end,1:nlevdecomp_full) = nan
#ifdef NITRIF_DENITRIF
    cns%smin_no3_vr(beg:end,1:nlevdecomp_full) = nan
    cns%smin_nh4_vr(beg:end,1:nlevdecomp_full) = nan
    cns%smin_no3(beg:end) = nan
    cns%smin_nh4(beg:end) = nan
#endif
    cns%cwdn(beg:end) = nan
    cns%sminn(beg:end) = nan
    cns%col_ntrunc(beg:end) = nan
    cns%seedn(beg:end) = nan
    cns%prod10n(beg:end) = nan
    cns%prod100n(beg:end) = nan
    cns%totprodn(beg:end) = nan
    cns%totlitn(beg:end) = nan
    cns%totsomn(beg:end) = nan
    cns%totlitn_1m(beg:end) = nan
    cns%totsomn_1m(beg:end) = nan
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
! !USES:
    use clm_varcon, only : spval
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

    allocate(cwf%qflx_surf(beg:end))
    allocate(cwf%qflx_drain(beg:end))

    cwf%qflx_surf(beg:end) = nan
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
  use clm_varcon, only : spval
!
! !DESCRIPTION:
! Initialize column carbon flux variables
!
! !USES:
    use clm_varctl , only : crop_prog
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

    allocate(ccf%hrv_deadstemc_to_prod10c(beg:end))        
    allocate(ccf%hrv_deadstemc_to_prod100c(beg:end))       
    allocate(ccf%m_decomp_cpools_to_fire_vr(beg:end,1:nlevdecomp_full,1:ndecomp_pools))
    allocate(ccf%m_decomp_cpools_to_fire(beg:end,1:ndecomp_pools))
    allocate(ccf%decomp_cascade_hr_vr(beg:end,1:nlevdecomp_full,1:ndecomp_cascade_transitions))
    allocate(ccf%decomp_cascade_hr(beg:end,1:ndecomp_cascade_transitions))
    allocate(ccf%decomp_cascade_ctransfer_vr(beg:end,1:nlevdecomp_full,1:ndecomp_cascade_transitions))
    allocate(ccf%decomp_cascade_ctransfer(beg:end,1:ndecomp_cascade_transitions))
    allocate(ccf%decomp_cpools_sourcesink(beg:end,1:nlevdecomp_full,1:ndecomp_pools))
    allocate(ccf%decomp_k(beg:end,1:nlevdecomp_full,1:ndecomp_cascade_transitions))
    allocate(ccf%t_scalar(beg:end,1:nlevdecomp_full))
    allocate(ccf%w_scalar(beg:end,1:nlevdecomp_full))
    allocate(ccf%hr_vr(beg:end,1:nlevdecomp_full))
    allocate(ccf%o_scalar(beg:end,1:nlevdecomp_full))
    allocate(ccf%som_c_leached(beg:end))
    allocate(ccf%decomp_cpools_leached(beg:end,1:ndecomp_pools))
    allocate(ccf%decomp_cpools_transport_tendency(beg:end,1:nlevdecomp_full,1:ndecomp_pools))

    allocate(ccf%phenology_c_to_litr_met_c(beg:end, 1:nlevdecomp_full))
    allocate(ccf%phenology_c_to_litr_cel_c(beg:end, 1:nlevdecomp_full))
    allocate(ccf%phenology_c_to_litr_lig_c(beg:end, 1:nlevdecomp_full))
    allocate(ccf%gap_mortality_c_to_litr_met_c(beg:end, 1:nlevdecomp_full))
    allocate(ccf%gap_mortality_c_to_litr_cel_c(beg:end, 1:nlevdecomp_full))
    allocate(ccf%gap_mortality_c_to_litr_lig_c(beg:end, 1:nlevdecomp_full))
    allocate(ccf%gap_mortality_c_to_cwdc(beg:end, 1:nlevdecomp_full))
    allocate(ccf%fire_mortality_c_to_cwdc(beg:end, 1:nlevdecomp_full))
    allocate(ccf%m_c_to_litr_met_fire(beg:end,1:nlevdecomp_full))
    allocate(ccf%m_c_to_litr_cel_fire(beg:end,1:nlevdecomp_full))
    allocate(ccf%m_c_to_litr_lig_fire(beg:end,1:nlevdecomp_full))
    allocate(ccf%harvest_c_to_litr_met_c(beg:end, 1:nlevdecomp_full))
    allocate(ccf%harvest_c_to_litr_cel_c(beg:end, 1:nlevdecomp_full))
    allocate(ccf%harvest_c_to_litr_lig_c(beg:end, 1:nlevdecomp_full))
    allocate(ccf%harvest_c_to_cwdc(beg:end, 1:nlevdecomp_full))

#ifdef NITRIF_DENITRIF
    allocate(ccf%phr_vr(beg:end,1:nlevdecomp_full))
#endif

!#ifdef CN
   !F. Li and S. Levis
    allocate(ccf%somc_fire(beg:end))
    allocate(ccf%lf_conv_cflux(beg:end))
    allocate(ccf%dwt_seedc_to_leaf(beg:end))
    allocate(ccf%dwt_seedc_to_deadstem(beg:end))
    allocate(ccf%dwt_conv_cflux(beg:end))
    allocate(ccf%dwt_prod10c_gain(beg:end))
    allocate(ccf%dwt_prod100c_gain(beg:end))
    allocate(ccf%dwt_frootc_to_litr_met_c(beg:end,1:nlevdecomp_full))
    allocate(ccf%dwt_frootc_to_litr_cel_c(beg:end,1:nlevdecomp_full))
    allocate(ccf%dwt_frootc_to_litr_lig_c(beg:end,1:nlevdecomp_full))
    allocate(ccf%dwt_livecrootc_to_cwdc(beg:end,1:nlevdecomp_full))
    allocate(ccf%dwt_deadcrootc_to_cwdc(beg:end,1:nlevdecomp_full))
    allocate(ccf%dwt_closs(beg:end))
    allocate(ccf%landuseflux(beg:end))
    allocate(ccf%landuptake(beg:end))
    allocate(ccf%prod10c_loss(beg:end))
    allocate(ccf%prod100c_loss(beg:end))
    allocate(ccf%product_closs(beg:end))
!#endif
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

!#if (defined CN)
    allocate(ccf%cwdc_hr(beg:end))
    allocate(ccf%cwdc_loss(beg:end))
    allocate(ccf%litterc_loss(beg:end))
!#endif

    ccf%m_c_to_litr_met_fire(beg:end,1:nlevdecomp_full)             = nan
    ccf%m_c_to_litr_cel_fire(beg:end,1:nlevdecomp_full)             = nan
    ccf%m_c_to_litr_lig_fire(beg:end,1:nlevdecomp_full)             = nan
    ccf%hrv_deadstemc_to_prod10c(beg:end)         = nan        
    ccf%hrv_deadstemc_to_prod100c(beg:end)        = nan       
    ccf%m_decomp_cpools_to_fire_vr(beg:end,1:nlevdecomp_full,1:ndecomp_pools)                = nan
    ccf%m_decomp_cpools_to_fire(beg:end,1:ndecomp_pools)                                     = nan
    ccf%decomp_cascade_hr_vr(beg:end,1:nlevdecomp_full,1:ndecomp_cascade_transitions)        = nan
    ccf%decomp_cascade_hr(beg:end,1:ndecomp_cascade_transitions)                             = nan
    ccf%decomp_cascade_ctransfer_vr(beg:end,1:nlevdecomp_full,1:ndecomp_cascade_transitions) = nan
    ccf%decomp_cascade_ctransfer(beg:end,1:ndecomp_cascade_transitions)                      = nan
    ccf%decomp_cpools_sourcesink(beg:end,1:nlevdecomp_full,1:ndecomp_pools)                  = nan
    ccf%decomp_k(beg:end,1:nlevdecomp_full,1:ndecomp_cascade_transitions)                    = spval
! Initialize these four below to spval to allow history to not average over inactive points.
    ccf%t_scalar(beg:end,1:nlevdecomp_full)                         = spval
    ccf%w_scalar(beg:end,1:nlevdecomp_full)                         = spval
    ccf%hr_vr(beg:end, 1:nlevdecomp_full)                           = spval
    ccf%o_scalar(beg:end, 1:nlevdecomp_full)                        = spval
    ccf%som_c_leached(beg:end)                                                      = nan 
    ccf%decomp_cpools_leached(beg:end,1:ndecomp_pools)                              = nan
    ccf%decomp_cpools_transport_tendency(beg:end,1:nlevdecomp_full,1:ndecomp_pools) = nan

    ccf%phenology_c_to_litr_met_c(beg:end, 1:nlevdecomp_full)                     = nan
    ccf%phenology_c_to_litr_cel_c(beg:end, 1:nlevdecomp_full)                     = nan
    ccf%phenology_c_to_litr_lig_c(beg:end, 1:nlevdecomp_full)                     = nan
    ccf%gap_mortality_c_to_litr_met_c(beg:end, 1:nlevdecomp_full)                 = nan
    ccf%gap_mortality_c_to_litr_cel_c(beg:end, 1:nlevdecomp_full)                 = nan
    ccf%gap_mortality_c_to_litr_lig_c(beg:end, 1:nlevdecomp_full)                 = nan
    ccf%gap_mortality_c_to_cwdc(beg:end, 1:nlevdecomp_full)                       = nan
    ccf%fire_mortality_c_to_cwdc(beg:end, 1:nlevdecomp_full)                      = nan
    ccf%harvest_c_to_litr_met_c(beg:end, 1:nlevdecomp_full)                       = nan
    ccf%harvest_c_to_litr_cel_c(beg:end, 1:nlevdecomp_full)                       = nan
    ccf%harvest_c_to_litr_lig_c(beg:end, 1:nlevdecomp_full)                       = nan
    ccf%harvest_c_to_cwdc(beg:end, 1:nlevdecomp_full)                             = nan

#ifdef NITRIF_DENITRIF
    ccf%phr_vr(beg:end,1:nlevdecomp_full)                              = nan
#endif
!#if (defined CN)
    !F. Li and S. Levis
    ccf%somc_fire(beg:end)                        = nan
    ccf%lf_conv_cflux(beg:end)                    = nan
    ccf%dwt_seedc_to_leaf(beg:end)                = 0.   !nan, followed what Greg did to disable dwt, fzeng
    ccf%dwt_seedc_to_deadstem(beg:end)            = 0.   !nan, followed what Greg did to disable dwt, fzeng
    ccf%dwt_conv_cflux(beg:end)                   = nan
    ccf%dwt_prod10c_gain(beg:end)                 = nan
    ccf%dwt_prod100c_gain(beg:end)                = nan
    ccf%dwt_frootc_to_litr_met_c(beg:end,1:nlevdecomp_full)             = nan
    ccf%dwt_frootc_to_litr_cel_c(beg:end,1:nlevdecomp_full)             = nan
    ccf%dwt_frootc_to_litr_lig_c(beg:end,1:nlevdecomp_full)             = nan
    ccf%dwt_livecrootc_to_cwdc(beg:end,1:nlevdecomp_full)           = nan
    ccf%dwt_deadcrootc_to_cwdc(beg:end,1:nlevdecomp_full)           = nan
    ccf%dwt_closs(beg:end)                        = nan
    ccf%landuseflux(beg:end)                      = nan
    ccf%landuptake(beg:end)                       = nan
    ccf%prod10c_loss(beg:end)                     = nan
    ccf%prod100c_loss(beg:end)                    = nan
    ccf%product_closs(beg:end)                    = nan
!#endif
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

!#if (defined CN)
    ccf%cwdc_hr(beg:end)                          = nan
    ccf%cwdc_loss(beg:end)                        = nan
    ccf%litterc_loss(beg:end)                     = nan
!#endif

  end subroutine init_column_cflux_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_column_ch4_type
!
! !INTERFACE:
#ifdef LCH4
  subroutine init_column_ch4_type(beg, end, cch4)
!
! !DESCRIPTION:
! Initialize column methane flux variables
!
  use clm_varcon, only : spval
  use clm_varpar, only : ngases
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (column_ch4_type), intent(inout):: cch4
!
! !REVISION HISTORY:
! Created by William J. Riley
!
!EOP
!------------------------------------------------------------------------

    allocate(cch4%ch4_prod_depth_sat(beg:end,1:nlevgrnd))
    allocate(cch4%ch4_prod_depth_unsat(beg:end,1:nlevgrnd))
    allocate(cch4%ch4_prod_depth_lake(beg:end,1:nlevgrnd))
    allocate(cch4%ch4_oxid_depth_sat(beg:end,1:nlevgrnd))
    allocate(cch4%ch4_oxid_depth_unsat(beg:end,1:nlevgrnd))
    allocate(cch4%ch4_oxid_depth_lake(beg:end,1:nlevgrnd))
    allocate(cch4%o2_oxid_depth_sat(beg:end,1:nlevgrnd))
    allocate(cch4%o2_oxid_depth_unsat(beg:end,1:nlevgrnd))
    allocate(cch4%o2_decomp_depth_sat(beg:end,1:nlevgrnd))
    allocate(cch4%o2_decomp_depth_unsat(beg:end,1:nlevgrnd))
    allocate(cch4%o2_aere_depth_sat(beg:end,1:nlevgrnd))
    allocate(cch4%o2_aere_depth_unsat(beg:end,1:nlevgrnd))
    allocate(cch4%co2_decomp_depth_sat(beg:end,1:nlevgrnd))
    allocate(cch4%co2_decomp_depth_unsat(beg:end,1:nlevgrnd))
    allocate(cch4%co2_oxid_depth_sat(beg:end,1:nlevgrnd))
    allocate(cch4%co2_oxid_depth_unsat(beg:end,1:nlevgrnd))
    allocate(cch4%ch4_aere_depth_sat(beg:end,1:nlevgrnd))
    allocate(cch4%ch4_aere_depth_unsat(beg:end,1:nlevgrnd))
    allocate(cch4%ch4_tran_depth_sat(beg:end,1:nlevgrnd))
    allocate(cch4%ch4_tran_depth_unsat(beg:end,1:nlevgrnd))
    allocate(cch4%co2_aere_depth_sat(beg:end,1:nlevgrnd))
    allocate(cch4%co2_aere_depth_unsat(beg:end,1:nlevgrnd))
    allocate(cch4%ch4_surf_aere_sat(beg:end))
    allocate(cch4%ch4_surf_aere_unsat(beg:end))
    allocate(cch4%ch4_ebul_depth_sat(beg:end,1:nlevgrnd))
    allocate(cch4%ch4_ebul_depth_unsat(beg:end,1:nlevgrnd))
    allocate(cch4%ch4_ebul_total_sat(beg:end))
    allocate(cch4%ch4_ebul_total_unsat(beg:end))
    allocate(cch4%ch4_surf_ebul_sat(beg:end))
    allocate(cch4%ch4_surf_ebul_unsat(beg:end))
    allocate(cch4%ch4_surf_ebul_lake(beg:end))
    allocate(cch4%conc_ch4_sat(beg:end,1:nlevgrnd))
    allocate(cch4%conc_ch4_unsat(beg:end,1:nlevgrnd))
    allocate(cch4%conc_ch4_lake(beg:end,1:nlevgrnd))
    allocate(cch4%ch4_surf_diff_sat(beg:end))
    allocate(cch4%ch4_surf_diff_unsat(beg:end))
    allocate(cch4%ch4_surf_diff_lake(beg:end))
    allocate(cch4%conc_o2_sat(beg:end,1:nlevgrnd))
    allocate(cch4%conc_o2_unsat(beg:end,1:nlevgrnd))
    allocate(cch4%conc_o2_lake(beg:end,1:nlevgrnd))
    allocate(cch4%ch4_dfsat_flux(beg:end))
    allocate(cch4%zwt_ch4_unsat(beg:end))
    allocate(cch4%fsat_bef(beg:end))
    allocate(cch4%lake_soilc(beg:end,1:nlevgrnd))
    allocate(cch4%lake_raw(beg:end))
    allocate(cch4%totcolch4(beg:end))
    allocate(cch4%fphr(beg:end,1:nlevgrnd))
    allocate(cch4%annsum_counter(beg:end))
    allocate(cch4%tempavg_somhr(beg:end))
    allocate(cch4%annavg_somhr(beg:end))
    allocate(cch4%tempavg_finrw(beg:end))
    allocate(cch4%annavg_finrw(beg:end))
    allocate(cch4%sif(beg:end))
    allocate(cch4%o2stress_unsat(beg:end,1:nlevgrnd))
    allocate(cch4%o2stress_sat(beg:end,1:nlevgrnd))
    allocate(cch4%ch4stress_unsat(beg:end,1:nlevgrnd))
    allocate(cch4%ch4stress_sat(beg:end,1:nlevgrnd))
    allocate(cch4%qflx_surf_lag(beg:end))
    allocate(cch4%finundated_lag(beg:end))
    allocate(cch4%layer_sat_lag(beg:end,1:nlevgrnd))


    cch4%ch4_prod_depth_sat(beg:end,1:nlevgrnd) = nan
    cch4%ch4_prod_depth_unsat(beg:end,1:nlevgrnd) = nan
    cch4%ch4_prod_depth_lake(beg:end,1:nlevgrnd) = nan
    cch4%ch4_oxid_depth_sat(beg:end,1:nlevgrnd) = nan
    cch4%ch4_oxid_depth_unsat(beg:end,1:nlevgrnd) = nan
    cch4%ch4_oxid_depth_lake(beg:end,1:nlevgrnd) = nan
    cch4%o2_oxid_depth_sat(beg:end,1:nlevgrnd) = nan
    cch4%o2_oxid_depth_unsat(beg:end,1:nlevgrnd) = nan
    cch4%o2_decomp_depth_sat(beg:end,1:nlevgrnd) = nan
    cch4%o2_decomp_depth_unsat(beg:end,1:nlevgrnd) = spval ! To detect first time-step for denitrification code
    cch4%o2_aere_depth_sat(beg:end,1:nlevgrnd) = nan
    cch4%o2_aere_depth_unsat(beg:end,1:nlevgrnd) = nan
    cch4%co2_decomp_depth_sat(beg:end,1:nlevgrnd) = nan
    cch4%co2_decomp_depth_unsat(beg:end,1:nlevgrnd) = nan
    cch4%co2_oxid_depth_sat(beg:end,1:nlevgrnd) = nan
    cch4%co2_oxid_depth_unsat(beg:end,1:nlevgrnd) = nan
    cch4%ch4_aere_depth_sat(beg:end,1:nlevgrnd) = nan
    cch4%ch4_aere_depth_unsat(beg:end,1:nlevgrnd) = nan
    cch4%ch4_tran_depth_sat(beg:end,1:nlevgrnd) = nan
    cch4%ch4_tran_depth_unsat(beg:end,1:nlevgrnd) = nan
    cch4%co2_aere_depth_sat(beg:end,1:nlevgrnd) = nan
    cch4%co2_aere_depth_unsat(beg:end,1:nlevgrnd) = nan
    cch4%ch4_surf_aere_sat(beg:end) = nan
    cch4%ch4_surf_aere_unsat(beg:end) = nan
    cch4%ch4_ebul_depth_sat(beg:end,1:nlevgrnd) = nan
    cch4%ch4_ebul_depth_unsat(beg:end,1:nlevgrnd) = nan
    cch4%ch4_ebul_total_sat(beg:end) = nan
    cch4%ch4_ebul_total_unsat(beg:end) = nan
    cch4%ch4_surf_ebul_sat(beg:end) = nan
    cch4%ch4_surf_ebul_unsat(beg:end) = nan
    cch4%ch4_surf_ebul_lake(beg:end) = nan
    cch4%conc_ch4_sat(beg:end,1:nlevgrnd) = spval ! To detect file input
    cch4%conc_ch4_unsat(beg:end,1:nlevgrnd) = spval ! To detect file input
    cch4%conc_ch4_lake(beg:end,1:nlevgrnd) = nan ! Just a diagnostic, so nan is fine
    cch4%ch4_surf_diff_sat(beg:end) = nan
    cch4%ch4_surf_diff_unsat(beg:end) = nan
    cch4%ch4_surf_diff_lake(beg:end) = nan
    cch4%conc_o2_sat(beg:end,1:nlevgrnd) = spval ! To detect file input
    cch4%conc_o2_unsat(beg:end,1:nlevgrnd) = spval ! To detect file input and detect first time-step for denitrification code
    cch4%conc_o2_lake(beg:end,1:nlevgrnd) = nan ! Just a diagnostic, so nan is fine
    cch4%ch4_dfsat_flux(beg:end) = nan
    cch4%zwt_ch4_unsat(beg:end) = nan
    cch4%fsat_bef(beg:end) = spval ! To detect first time-step
    cch4%lake_soilc(beg:end,1:nlevgrnd) = spval ! To detect file input
    cch4%lake_raw(beg:end) = nan
    cch4%totcolch4(beg:end) = spval ! To detect first time-step
    cch4%fphr(beg:end,1:nlevgrnd) = nan
    cch4%annsum_counter(beg:end) = spval ! To detect first time-step
    cch4%tempavg_somhr(beg:end) = nan
    cch4%annavg_somhr(beg:end) = spval ! To detect first year
    cch4%tempavg_finrw(beg:end) = nan
    cch4%annavg_finrw(beg:end) = spval ! To detect first year
    cch4%sif(beg:end) = nan
    cch4%o2stress_unsat(beg:end,1:nlevgrnd) = spval ! To detect file input
    cch4%o2stress_sat(beg:end,1:nlevgrnd) = spval ! To detect file input
    cch4%ch4stress_unsat(beg:end,1:nlevgrnd) = nan
    cch4%ch4stress_sat(beg:end,1:nlevgrnd) = nan
    cch4%qflx_surf_lag(beg:end) = spval ! To detect file input
    cch4%finundated_lag(beg:end) = spval ! To detect file input
    cch4%layer_sat_lag(beg:end,1:nlevgrnd) = spval ! To detect file input


! def CH4

  end subroutine init_column_ch4_type
#endif

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_column_nflux_type
!
! !INTERFACE:
  subroutine init_column_nflux_type(beg, end, cnf)
!
  use clm_varcon, only : spval
! !DESCRIPTION:
! Initialize column nitrogen flux variables
!
! !USES:
    use clm_varctl , only : crop_prog
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
    allocate(cnf%fert_to_sminn(beg:end))
    allocate(cnf%soyfixn_to_sminn(beg:end))    
    allocate(cnf%hrv_deadstemn_to_prod10n(beg:end))        
    allocate(cnf%hrv_deadstemn_to_prod100n(beg:end))       

    allocate(cnf%m_n_to_litr_met_fire(beg:end,1:nlevdecomp_full))
    allocate(cnf%m_n_to_litr_cel_fire(beg:end,1:nlevdecomp_full))
    allocate(cnf%m_n_to_litr_lig_fire(beg:end,1:nlevdecomp_full))
    allocate(cnf%sminn_to_plant(beg:end))
    allocate(cnf%potential_immob(beg:end))
    allocate(cnf%actual_immob(beg:end))
    allocate(cnf%gross_nmin(beg:end))
    allocate(cnf%net_nmin(beg:end))
    allocate(cnf%denit(beg:end))
    allocate(cnf%supplement_to_sminn(beg:end))
    allocate(cnf%m_decomp_npools_to_fire_vr(beg:end,1:nlevdecomp_full,1:ndecomp_pools))
    allocate(cnf%m_decomp_npools_to_fire(beg:end,1:ndecomp_pools))
    allocate(cnf%decomp_cascade_ntransfer_vr(beg:end,1:nlevdecomp_full,1:ndecomp_cascade_transitions))
    allocate(cnf%decomp_cascade_ntransfer(beg:end,1:ndecomp_cascade_transitions))
    allocate(cnf%decomp_cascade_sminn_flux_vr(beg:end,1:nlevdecomp_full,1:ndecomp_cascade_transitions))
    allocate(cnf%decomp_cascade_sminn_flux(beg:end,1:ndecomp_cascade_transitions))
    allocate(cnf%decomp_npools_sourcesink(beg:end,1:nlevdecomp_full,1:ndecomp_pools))

    allocate(cnf%phenology_n_to_litr_met_n(beg:end, 1:nlevdecomp_full))
    allocate(cnf%phenology_n_to_litr_cel_n(beg:end, 1:nlevdecomp_full))
    allocate(cnf%phenology_n_to_litr_lig_n(beg:end, 1:nlevdecomp_full))
    allocate(cnf%gap_mortality_n_to_litr_met_n(beg:end, 1:nlevdecomp_full))
    allocate(cnf%gap_mortality_n_to_litr_cel_n(beg:end, 1:nlevdecomp_full))
    allocate(cnf%gap_mortality_n_to_litr_lig_n(beg:end, 1:nlevdecomp_full))
    allocate(cnf%gap_mortality_n_to_cwdn(beg:end, 1:nlevdecomp_full))
    allocate(cnf%fire_mortality_n_to_cwdn(beg:end, 1:nlevdecomp_full))
    allocate(cnf%harvest_n_to_litr_met_n(beg:end, 1:nlevdecomp_full))
    allocate(cnf%harvest_n_to_litr_cel_n(beg:end, 1:nlevdecomp_full))
    allocate(cnf%harvest_n_to_litr_lig_n(beg:end, 1:nlevdecomp_full))
    allocate(cnf%harvest_n_to_cwdn(beg:end, 1:nlevdecomp_full))

#ifndef NITRIF_DENITRIF
    allocate(cnf%sminn_to_denit_decomp_cascade_vr(beg:end,1:nlevdecomp_full,1:ndecomp_cascade_transitions))
    allocate(cnf%sminn_to_denit_decomp_cascade(beg:end,1:ndecomp_cascade_transitions))
    allocate(cnf%sminn_to_denit_excess_vr(beg:end,1:nlevdecomp_full))
    allocate(cnf%sminn_to_denit_excess(beg:end))
    allocate(cnf%sminn_leached_vr(beg:end,1:nlevdecomp_full))
    allocate(cnf%sminn_leached(beg:end))
#else
    allocate(cnf%f_nit_vr(beg:end,1:nlevdecomp_full))
    allocate(cnf%f_denit_vr(beg:end,1:nlevdecomp_full))
    allocate(cnf%smin_no3_leached_vr(beg:end,1:nlevdecomp_full))
    allocate(cnf%smin_no3_leached(beg:end))
    allocate(cnf%smin_no3_runoff_vr(beg:end,1:nlevdecomp_full))
    allocate(cnf%smin_no3_runoff(beg:end))
    allocate(cnf%pot_f_nit_vr(beg:end,1:nlevdecomp_full))
    allocate(cnf%pot_f_nit(beg:end))
    allocate(cnf%pot_f_denit_vr(beg:end,1:nlevdecomp_full))
    allocate(cnf%pot_f_denit(beg:end))
    allocate(cnf%actual_immob_no3_vr(beg:end,1:nlevdecomp_full))
    allocate(cnf%actual_immob_nh4_vr(beg:end,1:nlevdecomp_full))
    allocate(cnf%smin_no3_to_plant_vr(beg:end,1:nlevdecomp_full))
    allocate(cnf%smin_nh4_to_plant_vr(beg:end,1:nlevdecomp_full))
    allocate(cnf%f_nit(beg:end))
    allocate(cnf%f_denit(beg:end))
    allocate(cnf%n2_n2o_ratio_denit_vr(beg:end,1:nlevdecomp_full))
    allocate(cnf%f_n2o_denit(beg:end))
    allocate(cnf%f_n2o_denit_vr(beg:end,1:nlevdecomp_full))
    allocate(cnf%f_n2o_nit(beg:end))
    allocate(cnf%f_n2o_nit_vr(beg:end,1:nlevdecomp_full))

    allocate(cnf%smin_no3_massdens_vr(beg:end,1:nlevdecomp_full))
    allocate(cnf%soil_bulkdensity(beg:end,1:nlevdecomp_full))
    allocate(cnf%k_nitr_t_vr(beg:end,1:nlevdecomp_full))
    allocate(cnf%k_nitr_ph_vr(beg:end,1:nlevdecomp_full))
    allocate(cnf%k_nitr_h2o_vr(beg:end,1:nlevdecomp_full)) 
    allocate(cnf%k_nitr_vr(beg:end,1:nlevdecomp_full)) 
    allocate(cnf%wfps_vr(beg:end,1:nlevdecomp_full)) 
    allocate(cnf%fmax_denit_carbonsubstrate_vr(beg:end,1:nlevdecomp_full)) 
    allocate(cnf%fmax_denit_nitrate_vr(beg:end,1:nlevdecomp_full)) 
    allocate(cnf%f_denit_base_vr(beg:end,1:nlevdecomp_full))
    allocate(cnf%diffus(beg:end,1:nlevdecomp_full))
    allocate(cnf%ratio_k1(beg:end,1:nlevdecomp_full))
    allocate(cnf%ratio_no3_co2(beg:end,1:nlevdecomp_full))
    allocate(cnf%soil_co2_prod(beg:end,1:nlevdecomp_full))
    allocate(cnf%fr_WFPS(beg:end,1:nlevdecomp_full))

    allocate(cnf%r_psi(beg:end,1:nlevdecomp_full))
    allocate(cnf%anaerobic_frac(beg:end,1:nlevdecomp_full))
#endif
    allocate(cnf%potential_immob_vr(beg:end,1:nlevdecomp_full))
    allocate(cnf%actual_immob_vr(beg:end,1:nlevdecomp_full))
    allocate(cnf%sminn_to_plant_vr(beg:end,1:nlevdecomp_full))
    allocate(cnf%supplement_to_sminn_vr(beg:end,1:nlevdecomp_full))
    allocate(cnf%gross_nmin_vr(beg:end,1:nlevdecomp_full))
    allocate(cnf%net_nmin_vr(beg:end,1:nlevdecomp_full))
    allocate(cnf%dwt_seedn_to_leaf(beg:end))
    allocate(cnf%dwt_seedn_to_deadstem(beg:end))
    allocate(cnf%dwt_conv_nflux(beg:end))
    allocate(cnf%dwt_prod10n_gain(beg:end))
    allocate(cnf%dwt_prod100n_gain(beg:end))
    allocate(cnf%dwt_frootn_to_litr_met_n(beg:end,1:nlevdecomp_full))
    allocate(cnf%dwt_frootn_to_litr_cel_n(beg:end,1:nlevdecomp_full))
    allocate(cnf%dwt_frootn_to_litr_lig_n(beg:end,1:nlevdecomp_full))
    allocate(cnf%dwt_livecrootn_to_cwdn(beg:end,1:nlevdecomp_full))
    allocate(cnf%dwt_deadcrootn_to_cwdn(beg:end,1:nlevdecomp_full))
    allocate(cnf%dwt_nloss(beg:end))
    allocate(cnf%prod10n_loss(beg:end))
    allocate(cnf%prod100n_loss(beg:end))
    allocate(cnf%product_nloss(beg:end))
    allocate(cnf%col_ninputs(beg:end))
    allocate(cnf%col_noutputs(beg:end))
    allocate(cnf%col_fire_nloss(beg:end))
    allocate(cnf%som_n_leached(beg:end))
    allocate(cnf%decomp_npools_leached(beg:end,1:ndecomp_pools))
    allocate(cnf%decomp_npools_transport_tendency(beg:end,1:nlevdecomp_full,1:ndecomp_pools))
    
    cnf%ndep_to_sminn(beg:end) = nan
    cnf%nfix_to_sminn(beg:end) = nan
    cnf%fert_to_sminn(beg:end) = nan
    cnf%soyfixn_to_sminn(beg:end) = nan
    cnf%hrv_deadstemn_to_prod10n(beg:end) = nan        
    cnf%hrv_deadstemn_to_prod100n(beg:end) = nan       
    cnf%m_n_to_litr_met_fire(beg:end,1:nlevdecomp_full) = nan
    cnf%m_n_to_litr_cel_fire(beg:end,1:nlevdecomp_full) = nan
    cnf%m_n_to_litr_lig_fire(beg:end,1:nlevdecomp_full) = nan
    cnf%sminn_to_plant(beg:end) = nan
    cnf%potential_immob(beg:end) = nan
    cnf%actual_immob(beg:end) = nan
    cnf%gross_nmin(beg:end) = nan
    cnf%net_nmin(beg:end) = nan
    cnf%denit(beg:end) = nan
    cnf%supplement_to_sminn(beg:end) = nan
    cnf%m_decomp_npools_to_fire_vr(beg:end,1:nlevdecomp_full,1:ndecomp_pools)                 = nan
    cnf%m_decomp_npools_to_fire(beg:end,1:ndecomp_pools)                                      = nan
    cnf%decomp_cascade_ntransfer_vr(beg:end,1:nlevdecomp_full,1:ndecomp_cascade_transitions)  = nan
    cnf%decomp_cascade_ntransfer(beg:end,1:ndecomp_cascade_transitions)                       = nan
    cnf%decomp_cascade_sminn_flux_vr(beg:end,1:nlevdecomp_full,1:ndecomp_cascade_transitions) = nan
    cnf%decomp_cascade_sminn_flux(beg:end,1:ndecomp_cascade_transitions)                      = nan
    cnf%decomp_npools_sourcesink(beg:end,1:nlevdecomp_full,1:ndecomp_pools)                   = nan
    
    cnf%phenology_n_to_litr_met_n(beg:end, 1:nlevdecomp_full)                     = nan
    cnf%phenology_n_to_litr_cel_n(beg:end, 1:nlevdecomp_full)                     = nan
    cnf%phenology_n_to_litr_lig_n(beg:end, 1:nlevdecomp_full)                     = nan
    cnf%gap_mortality_n_to_litr_met_n(beg:end, 1:nlevdecomp_full)                 = nan
    cnf%gap_mortality_n_to_litr_cel_n(beg:end, 1:nlevdecomp_full)                 = nan
    cnf%gap_mortality_n_to_litr_lig_n(beg:end, 1:nlevdecomp_full)                 = nan
    cnf%gap_mortality_n_to_cwdn(beg:end, 1:nlevdecomp_full)                       = nan
    cnf%fire_mortality_n_to_cwdn(beg:end, 1:nlevdecomp_full)                      = nan
    cnf%harvest_n_to_litr_met_n(beg:end, 1:nlevdecomp_full)                       = nan
    cnf%harvest_n_to_litr_cel_n(beg:end, 1:nlevdecomp_full)                       = nan
    cnf%harvest_n_to_litr_lig_n(beg:end, 1:nlevdecomp_full)                       = nan
    cnf%harvest_n_to_cwdn(beg:end, 1:nlevdecomp_full)                             = nan

#ifndef NITRIF_DENITRIF
    cnf%sminn_to_denit_decomp_cascade_vr(beg:end,1:nlevdecomp_full,1:ndecomp_cascade_transitions) = nan
    cnf%sminn_to_denit_decomp_cascade(beg:end,1:ndecomp_cascade_transitions) = nan
    cnf%sminn_to_denit_excess_vr(beg:end,1:nlevdecomp_full) = nan
    cnf%sminn_to_denit_excess(beg:end) = nan
    cnf%sminn_leached_vr(beg:end,1:nlevdecomp_full) = nan
    cnf%sminn_leached(beg:end) = nan
#else
    cnf%f_nit_vr(beg:end,1:nlevdecomp_full) = nan
    cnf%f_denit_vr(beg:end,1:nlevdecomp_full) = nan
    cnf%smin_no3_leached_vr(beg:end,1:nlevdecomp_full) = nan
    cnf%smin_no3_leached(beg:end) = nan
    cnf%smin_no3_runoff_vr(beg:end,1:nlevdecomp_full) = nan
    cnf%smin_no3_runoff(beg:end) = nan
    cnf%pot_f_nit_vr(beg:end,1:nlevdecomp_full) = nan
    cnf%pot_f_nit(beg:end) = nan
    cnf%pot_f_denit_vr(beg:end,1:nlevdecomp_full) = nan
    cnf%pot_f_denit(beg:end) = nan
    cnf%actual_immob_no3_vr(beg:end,1:nlevdecomp_full) = nan
    cnf%actual_immob_nh4_vr(beg:end,1:nlevdecomp_full) = nan
    cnf%smin_no3_to_plant_vr(beg:end,1:nlevdecomp_full) = nan
    cnf%smin_nh4_to_plant_vr(beg:end,1:nlevdecomp_full) = nan
    cnf%f_nit(beg:end) = nan
    cnf%f_denit(beg:end) = nan
    cnf%n2_n2o_ratio_denit_vr(beg:end,1:nlevdecomp_full) = nan
    cnf%f_n2o_denit(beg:end) = nan
    cnf%f_n2o_denit_vr(beg:end,1:nlevdecomp_full) = nan
    cnf%f_n2o_nit(beg:end) = nan
    cnf%f_n2o_nit_vr(beg:end,1:nlevdecomp_full) = nan

    cnf%smin_no3_massdens_vr(beg:end,1:nlevdecomp_full) = nan
    cnf%soil_bulkdensity(beg:end,1:nlevdecomp_full) = nan
    cnf%k_nitr_t_vr(beg:end,1:nlevdecomp_full) = nan
    cnf%k_nitr_ph_vr(beg:end,1:nlevdecomp_full) = nan
    cnf%k_nitr_h2o_vr(beg:end,1:nlevdecomp_full) = nan
    cnf%k_nitr_vr(beg:end,1:nlevdecomp_full) = nan 
    cnf%wfps_vr(beg:end,1:nlevdecomp_full) = nan 
    cnf%fmax_denit_carbonsubstrate_vr(beg:end,1:nlevdecomp_full) = nan 
    cnf%fmax_denit_nitrate_vr(beg:end,1:nlevdecomp_full) = nan 
    cnf%f_denit_base_vr(beg:end,1:nlevdecomp_full) = nan
    cnf%diffus(beg:end,1:nlevdecomp_full) = spval
    cnf%ratio_k1(beg:end,1:nlevdecomp_full) = nan
    cnf%ratio_no3_co2(beg:end,1:nlevdecomp_full) = spval
    cnf%soil_co2_prod(beg:end,1:nlevdecomp_full) = nan
    cnf%fr_WFPS(beg:end,1:nlevdecomp_full) = spval

    cnf%r_psi(beg:end,1:nlevdecomp_full) = spval
    cnf%anaerobic_frac(beg:end,1:nlevdecomp_full) = spval
#endif
    cnf%potential_immob_vr(beg:end,1:nlevdecomp_full) = nan
    cnf%actual_immob_vr(beg:end,1:nlevdecomp_full) = nan
    cnf%sminn_to_plant_vr(beg:end,1:nlevdecomp_full) = nan
    cnf%supplement_to_sminn_vr(beg:end,1:nlevdecomp_full) = nan
    cnf%gross_nmin_vr(beg:end,1:nlevdecomp_full) = nan
    cnf%net_nmin_vr(beg:end,1:nlevdecomp_full) = nan
    cnf%dwt_seedn_to_leaf(beg:end) = nan
    cnf%dwt_seedn_to_deadstem(beg:end) = nan
    cnf%dwt_conv_nflux(beg:end) = nan
    cnf%dwt_prod10n_gain(beg:end) = nan
    cnf%dwt_prod100n_gain(beg:end) = nan
    cnf%dwt_frootn_to_litr_met_n(beg:end,1:nlevdecomp_full) = nan
    cnf%dwt_frootn_to_litr_cel_n(beg:end,1:nlevdecomp_full) = nan
    cnf%dwt_frootn_to_litr_lig_n(beg:end,1:nlevdecomp_full) = nan
    cnf%dwt_livecrootn_to_cwdn(beg:end,1:nlevdecomp_full) = nan
    cnf%dwt_deadcrootn_to_cwdn(beg:end,1:nlevdecomp_full) = nan
    cnf%dwt_nloss(beg:end) = nan
    cnf%prod10n_loss(beg:end) = nan
    cnf%prod100n_loss(beg:end) = nan
    cnf%product_nloss(beg:end) = nan
    cnf%col_ninputs(beg:end) = nan
    cnf%col_noutputs(beg:end) = nan
    cnf%col_fire_nloss(beg:end) = nan
    cnf%som_n_leached(beg:end)                                                      = nan 
    cnf%decomp_npools_leached(beg:end,1:ndecomp_pools)                              = nan
    cnf%decomp_npools_transport_tendency(beg:end,1:nlevdecomp_full,1:ndecomp_pools) = nan


  end subroutine init_column_nflux_type




#if (defined CNDV)
!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_gridcell_dgvstate_type
!
! !INTERFACE:
  subroutine init_gridcell_dgvstate_type(beg, end, gps)
!
! !DESCRIPTION:
! Initialize gridcell DGVM variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (gridcell_dgvstate_type), intent(inout):: gps
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(gps%agdd20(beg:end))
    allocate(gps%tmomin20(beg:end))
    allocate(gps%t10min(beg:end))
    gps%agdd20(beg:end) = nan
    gps%tmomin20(beg:end) = nan
    gps%t10min(beg:end) = nan

  end subroutine init_gridcell_dgvstate_type
#endif

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_gridcell_ch4_type
!
! !INTERFACE:
#ifdef LCH4
  subroutine init_gridcell_ch4_type(beg, end, gch4)
!
! !DESCRIPTION:
! Initialize gridcell ch4 variables
!
  use clm_varpar, only: ngases
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (gridcell_ch4_type), intent(inout):: gch4
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------
    allocate(gch4%c_atm(beg:end,1:ngases))
    allocate(gch4%ch4co2f(beg:end))
    allocate(gch4%ch4prodg(beg:end))
    allocate(gch4%nem(beg:end))

    gch4%c_atm(beg:end,1:ngases) = nan
    gch4%ch4co2f(beg:end) = nan
    gch4%ch4prodg(beg:end) = nan
    gch4%nem(beg:end) = nan

  end subroutine init_gridcell_ch4_type
#endif

end module clmtypeInitMod
