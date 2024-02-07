
module CNFireMod
!#ifdef CN

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: CNFireMod
!
! !DESCRIPTION:
! module for fire dynamics 
! created in Nov, 2012  and revised in Apr, 2013 by F. Li and S. Levis
! based on Li et al. (2012a,b; 2013)"
! Fire-related parameters were calibrated or tuned in Apr, 2013 based on the 
! 20th Century transient simulations at f19_g16 with (newfire05_clm45sci15_clm4_0_58) 
! a CLM4.5 version, Qian et al. (2006) atmospheric forcing, and
! climatological lightning data.
!
! !USES:
  use MAPL_ConstantsMod, ONLY: r8 => MAPL_R4
  use shr_const_mod  , only: SHR_CONST_PI,SHR_CONST_TKFRZ
!  use shr_strdata_mod, only: shr_strdata_type, shr_strdata_create, shr_strdata_print, &
!                             shr_strdata_advance
  use subgridAveMod  , only: p2c
  use clm_varctl     , only: iulog, use_nofire
  use clm_varpar     , only: nlevdecomp, ndecomp_pools
  use clm_varcon     , only: dzsoi_decomp
  use pftvarcon      , only: noveg, nbrdlf_evr_trp_tree
!  use catch_types    , only: cn_param_type
!  use spmdMod        , only: masterproc, mpicom, comp_id
!  use fileutils      , only: getavu, relavu
!  use controlMod     , only: NLFilename
!  use decompMod      , only: gsmap_lnd_gdc2glo
!  use domainMod      , only: ldomain
!  use mct_mod
  implicit none
  save
  private
! !PUBLIC TYPES:

!! fzeng: we don't use CNFireInit or CNFireInterp

! !PUBLIC MEMBER FUNCTIONS:
! public :: CNFireInit    ! Initialization of CNFire
! public :: CNFireInterp  ! Interpolate fire data
  public :: CNFireArea    ! Calculate fire area
  public :: CNFireFluxes  ! Calculate fire fluxes

! !PRIVATE MEMBER FUNCTIONS:
! private :: hdm_init    ! position datasets for dynamic human population density
! private :: hdm_interp  ! interpolates between two years of human pop. density file data
! private :: lnfm_init   ! position datasets for Lightning
! private :: lnfm_interp ! interpolates between two years of Lightning file data

! !PRIVATE MEMBER DATA:
  real(r8), pointer     :: forc_lnfm(:)        ! Lightning frequency, unit must be counts/km2/hr, fzeng, 24 Aug 2017
  real(r8), pointer     :: forc_hdm(:)         ! Human population density, unit must be person/km2, fzeng, 24 Aug 2017
  real(r8), parameter   :: secsphr = 3600._r8  ! Seconds in an hour
  real(r8), parameter   :: borealat = 40._r8   ! Latitude for boreal peat fires

!! hdm and lnfm data are read in another way, fzeng 
!  type(shr_strdata_type) :: sdat_hdm    ! Human population density input data stream
!  type(shr_strdata_type) :: sdat_lnfm   ! Lightning input data stream
!
! !REVISION HISTORY:
!
!EOP
!-----------------------------------------------------------------------

contains

!!-----------------------------------------------------------------------
!!BOP
!!
!! !IROUTINE: CNFireInit
!!
!! !INTERFACE:
!subroutine CNFireInit( begg, endg )
!
!! !DESCRIPTION:
!! Initialize CN Fire module
!! Not used in Catchment-CN. hdm and lnfm data are read in another way. fzeng
!!
!! !USES:
!!
!! !ARGUMENTS:
!   implicit none
!   integer, intent(IN) :: begg, endg   ! gridcell index bounds
!
!! !REVISION HISTORY:
!! !LOCAL VARIABLES:
!!EOP
!!-----------------------------------------------------------------------
!  call hdm_init(   begg, endg )
!  call hdm_interp( )
!  call lnfm_init(  begg, endg )
!  call lnfm_interp()

!-----------------------------------------------------------------------
!end subroutine CNFireInit


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNFireInterp
!
! !INTERFACE:
!subroutine CNFireInterp()
!
! !DESCRIPTION:
! Interpolate CN Fire datasets
!
! !USES:
!
! !ARGUMENTS:
!   implicit none
!
! !REVISION HISTORY:
! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------
!   call hdm_interp()
!   call lnfm_interp()

!-----------------------------------------------------------------------
!end subroutine CNFireInterp

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNFireArea
!
! !INTERFACE:
subroutine CNFireArea (num_soilc, filter_soilc, num_soilp, filter_soilp)
!
! !DESCRIPTION:
! Computes column-level burned area in each timestep

! fzeng, 24 Jan 2018
! Change nc4_grass to nc4_grass2 throughout CNFireArea because Catchment-CN has split types for c4 grass
!
! !USES:
   use clmtype
   use clm_time_manager, only: get_step_size, get_days_per_year, get_curr_date, get_nstep
   use clm_varpar      , only: max_pft_per_col, numpft
   use clm_varcon      , only: secspday
!   use shr_infnan_mod  , only: shr_infnan_isnan
   use, intrinsic :: ieee_arithmetic, only: shr_infnan_isnan => ieee_is_nan
   use clm_varctl      , only: fpftdyn
   use pftvarcon       , only: nc4_grass2, nc3crop, ndllf_evr_tmp_tree, &
                               nbrdlf_evr_trp_tree, nbrdlf_dcd_trp_tree,      &
                               nbrdlf_evr_shrub
!
! !ARGUMENTS:
   implicit none
   integer, intent(in) :: num_soilc       ! number of soil columns in filter
   integer, intent(in) :: filter_soilc(:) ! filter for soil columns
   integer, intent(in) :: num_soilp       ! number of soil pfts in filter
   integer, intent(in) :: filter_soilp(:) ! filter for soil pfts
!
! !CALLED FROM:
! subroutine CNEcosystemDyn in module CNEcosystemDynMod.F90
!
! !REVISION HISTORY:
! !LOCAL VARIABLES:
! local pointers to implicit in scalars
!
   ! pft-level
   real(r8), pointer :: prec10(:)       ! 10-day running mean of tot. precipitation
   real(r8), pointer :: prec60(:)       ! 60-day running mean of tot. precipitation
   real(r8), pointer :: lfpftd(:)       ! decrease of pft weight (0-1) on the col. for timestep 
   real(r8), pointer :: wtcol(:)        ! pft weight on the column
   integer , pointer :: ivt(:)          ! vegetation type for this pft
   real(r8), pointer :: fveg(:)         ! vegetation fraction for this pft
   real(r8), pointer :: deadcrootc(:)         ! (gC/m2) dead coarse root C
   real(r8), pointer :: deadcrootc_storage(:) ! (gC/m2) dead coarse root C storage
   real(r8), pointer :: deadcrootc_xfer(:)    ! (gC/m2) dead coarse root C transfer
   real(r8), pointer :: frootc(:)             ! (gC/m2) fine root C
   real(r8), pointer :: frootc_storage(:)     ! (gC/m2) fine root C storage
   real(r8), pointer :: frootc_xfer(:)        ! (gC/m2) fine root C transfer
   real(r8), pointer :: livecrootc(:)         ! (gC/m2) live coarse root C
   real(r8), pointer :: livecrootc_storage(:) ! (gC/m2) live coarse root C storage
   real(r8), pointer :: livecrootc_xfer(:)    ! (gC/m2) live coarse root C transfer
   real(r8), pointer :: totvegc(:)            ! (gC/m2) total vegetation carbon, excluding cpool 
   real(r8), pointer :: btran2(:)             ! root zone soil wetness
   real(r8), pointer :: leafc(:)              ! (gC/m2) leaf C
   real(r8), pointer :: leafc_storage(:)      ! (gC/m2) leaf C storage
   real(r8), pointer :: leafc_xfer(:)         ! (gC/m2) leaf C transfer
   integer , pointer :: burndate(:)           ! burn date for crop
   real(r8), pointer :: cpool(:)              ! (gC/m2) temporary photosynthate C pool


   ! column-level
   real(r8), pointer :: fsat(:)        ! fractional area with water table at surface
   real(r8), pointer :: lfc(:)         ! conversion area frac. of BET+BDT that haven't burned before
   real(r8), pointer :: cwtgcell(:)    ! column's weight relative to corresponding gridcell
   real(r8), pointer :: dtrotr_col(:)  ! annual decreased fraction coverage of BET+BDT on gridcell
   real(r8), pointer :: trotr1_col(:)  ! pft weight of BET on the gridcell (0-1)
   real(r8), pointer :: trotr2_col(:)  ! pft weight of BDT on the gridcell (0-1)
   real(r8), pointer :: prec10_col(:)  ! 10-day running mean of tot. precipitation
   real(r8), pointer :: prec60_col(:)  ! 60-day running mean of tot. precipitation
   integer , pointer :: npfts(:)       ! number of pfts on the column
   integer , pointer :: pfti(:)        ! pft index array
   integer , pointer :: cgridcell(:)   ! gridcell of corresponding column
   real(r8), pointer :: wf(:)          ! soil water as frac. of whc for top 0.05 m   
   real(r8), pointer :: wf2(:)         ! soil water as frac. of whc for top 0.17 m   
   real(r8), pointer :: tsoi17(:)      ! soil T for top 0.17 m   
   real(r8), pointer :: gdp_lf(:)      ! gdp data
   real(r8), pointer :: peatf_lf(:)    ! peatland fraction data
   integer, pointer :: abm_lf(:)       ! prescribed crop fire time
   real(r8), pointer :: totlitc(:)     ! (gC/m2) total lit C (column-level mean)
   real(r8), pointer :: fsr_col(:)     ! fire spread rate at column level (m/s)
   real(r8), pointer :: fd_col(:)      ! fire duration at column level (second) 
   real(r8), pointer :: rootc_col(:)   ! root carbon 
   real(r8), pointer :: baf_crop(:)    ! burned area fraction for cropland
   real(r8), pointer :: baf_peatf(:)   ! burned area fraction for peatland
   real(r8), pointer :: fbac(:)        ! total burned area out of conversion 
   real(r8), pointer :: fbac1(:)       ! burned area out of conversion region due to land use fire
   real(r8), pointer :: cropf_col(:)   ! cropland fraction in veg column
   real(r8), pointer :: btran_col(:)   ! transpiration wetness factor (0 to 1)
   real(r8), pointer :: wtlf(:)        ! fractional coverage of non-crop PFTs
   real(r8), pointer :: lfwt(:)        ! fractional coverage of non-crop and non-bare-soil PFTs
   real(r8), pointer :: totvegc_col(:) ! totvegc at column level
   real(r8), pointer :: leafc_col(:)   ! leaf carbon at column level
   real(r8), pointer :: lgdp_col(:)    ! gdp limitation factor for nfire
   real(r8), pointer :: lgdp1_col(:)   ! gdp limitation factor for baf per fire
   real(r8), pointer :: lpop_col(:)    ! pop limitation factor for baf per fire
   real(r8), pointer :: fuelc(:)       ! fuel avalability factor for Reg.C
   real(r8), pointer :: fuelc_crop(:)  ! fuel avalability factor for Reg.A
   real(r8), pointer  :: decomp_cpools_vr(:,:,:) ! (gC/m3)  vert.-resolved decomposing c pools
   real(r8), pointer :: cpool_col(:)   ! temporary photosynthate C pool at column level

   ! grid-level
   real(r8), pointer :: latdeg(:)      ! latitude (degrees)
   real(r8), pointer :: forc_rain(:)   ! rain rate [mm/s]
   real(r8), pointer :: forc_snow(:)   ! snow rate [mm/s]
   real(r8), pointer :: forc_rh(:)     ! atmospheric relative humidity (%)  
   real(r8), pointer :: forc_t(:)      ! atmospheric temperature (Kelvin)
   real(r8), pointer :: forc_wind(:)   ! atmospheric wind speed (m/s)
   
   ! ecophysiological constants, fzeng
   real(r8), pointer :: fsr_pft(:)     ! Fire spread rate (m/s)
   real(r8), pointer :: fd_pft(:)      ! Fire duration (hr)
   
! local pointers to implicit in/out scalars
!
   ! column-level
   real(r8), pointer :: nfire(:)       ! fire counts (count/km2/second), valid only in Reg. C
   real(r8), pointer :: farea_burned(:)! fractional area burned per second
   logical, pointer :: is_cwd(:)       ! TRUE => pool is a cwd pool
!
! !OTHER LOCAL VARIABLES:
   ! modified according to clm4_6_00, fzeng, 7 May 2019
   real(r8), parameter  :: lfuel=75._r8 !110._r8   ! lower threshold of fuel mass (gC/m2) for ignition, Li et al., BG 2014, fzeng, 7 May 2019
   real(r8), parameter  :: ufuel=1050._r8  ! upper threshold of fuel mass(gC/m2) for ignition 
   real(r8), parameter  :: g0=0.05_r8      ! g(W) when W=0 m/s

   ! a1 parameter for cropland fire in (Li et. al., 2014), but changed from
   ! /timestep to /hr
   real(r8), parameter :: cropfire_a1 = 0.3_r8 !0.31_r8    ! hr-1, see clm4_6_00 and Li et al., BG 2014, fzeng, 7 May 2019

   ! c parameter for peatland fire in Li et. al. (2013)
   ! boreal peat fires (was different in paper),changed from /timestep to /hr
   real(r8), parameter :: boreal_peatfire_c = 4.2e-5_r8 !4.2d-5    ! hr-1, see clm4_6_00, fzeng, 7 May 2019
   ! non-boreal peat fires (was different in paper)
   real(r8), parameter :: non_boreal_peatfire_c = 0.001_r8     ! hr-1, fzeng, 24 Aug 2017

   integer :: g,l,c,p,pi,j,fc,fp,kyr, kmo, kda, mcsec   ! index variables
   real(r8):: dt        ! time step variable (s)
   real(r8):: m         ! top-layer soil moisture (proportion)
   real(r8):: dayspyr   ! days per year
   real(r8) ::cli       !
   real(r8), parameter ::cli_scale = 0.035_r8   !global constant for deforestation fires (/d), see clm4_6_00, fzeng, 7 May 2019
   real(r8) ::cri       !
   real(r8):: fb        ! availability of fuel 
   real(r8):: fhd       ! impact of hd on agricultural fire
   real(r8):: fgdp      ! impact of gdp on agricultural fire
   real(r8):: fire_m    ! combustability of fuel on fire occurrence
   real(r8):: spread_m  ! combustability of fuel on fire spread
   real(r8):: Lb_lf     ! length-to-breadth ratio added by Lifang 
   integer :: i_cwd     ! cwd pool
   real(r8) :: lh       !
   real(r8) :: fs       !
   real(r8) :: ig       !
   real(r8) :: hdmlf    ! human density

   ! local variables for PFT fire modifications
   logical :: pft_4_flag ! flag for presence of PFT 4
   real(r8):: fire_m_tmp
   real, dimension(0:numpft) :: fire_m_fac ! array to parameterize combustibility sensitivity to soil moisture separately for each PFT
!EOP
!-----------------------------------------------------------------------
! declare fire_m_fac array used to compute PFT-dependent combustibility (jkolassa 08/2020)  

data     fire_m_fac    /  SHR_CONST_PI, SHR_CONST_PI,SHR_CONST_PI,SHR_CONST_PI,SHR_CONST_PI,SHR_CONST_PI,SHR_CONST_PI,SHR_CONST_PI,SHR_CONST_PI,SHR_CONST_PI,SHR_CONST_PI,SHR_CONST_PI,SHR_CONST_PI,SHR_CONST_PI,SHR_CONST_PI,SHR_CONST_PI,SHR_CONST_PI,SHR_CONST_PI,SHR_CONST_PI,SHR_CONST_PI/

! only allow fires in very dry conditions in PFTs 4 and 6 (jkolassa 08/2020)
! 02/2024: changed value for fire_m_fac(4,6) from 10.e15 to 10.e5 to avoid floating point exception in fire_m_tmp calculation with GNU on NCCS/SLES15.
fire_m_fac(4) = 10.e5 
fire_m_fac(6) = 10.e5 

! assign local pointers to derived type members (pft-level)
  wtcol              =>pft%wtcol
  ivt                =>pft%itype 
  prec60             => pps%prec60
  prec10             => pps%prec10
  deadcrootc         => pcs%deadcrootc
  deadcrootc_storage => pcs%deadcrootc_storage
  deadcrootc_xfer    => pcs%deadcrootc_xfer
  frootc             => pcs%frootc
  frootc_storage     => pcs%frootc_storage
  frootc_xfer        => pcs%frootc_xfer
  livecrootc         => pcs%livecrootc
  livecrootc_storage => pcs%livecrootc_storage
  livecrootc_xfer    => pcs%livecrootc_xfer
  totvegc            => pcs%totvegc
  btran2             => pps%btran2
  leafc              => pcs%leafc
  leafc_storage      => pcs%leafc_storage
  leafc_xfer         => pcs%leafc_xfer
  lfpftd             => pps%lfpftd
  burndate           => pps%burndate
  cpool              => pcs%cpool


 ! assign local pointers to derived type members (column-level)
  cwtgcell         =>col%wtgcell
  npfts            =>col%npfts
  pfti             =>col%pfti
  wf               => cps%wf
  wf2              => cps%wf2
  tsoi17           => ces%tsoi17
  farea_burned     => cps%farea_burned
  baf_crop         => cps%baf_crop 
  baf_peatf        => cps%baf_peatf 
  fbac             => cps%fbac
  fbac1            => cps%fbac1
  cropf_col        => cps%cropf_col 
  gdp_lf           => cps%gdp_lf
  peatf_lf         => cps%peatf_lf
  abm_lf           => cps%abm_lf
  nfire            => cps%nfire             
  totlitc          => ccs%totlitc
  fsr_col          => cps%fsr_col 
  fd_col           => cps%fd_col   
  rootc_col        => ccs%rootc_col 
  totvegc_col      => ccs%totvegc_col     
  leafc_col        => ccs%leafc_col     
  lgdp_col         => cps%lgdp_col  
  lgdp1_col        => cps%lgdp1_col  
  lpop_col         => cps%lpop_col  
  fuelc            => ccs%fuelc 
  fuelc_crop       => ccs%fuelc_crop  
  btran_col        => cps%btran_col
  wtlf             => cps%wtlf  
  lfwt             => cps%lfwt    
  cgridcell        =>col%gridcell
  trotr1_col       => cps%trotr1_col
  trotr2_col       => cps%trotr2_col
  dtrotr_col       => cps%dtrotr_col
  prec60_col       => cps%prec60_col
  prec10_col       => cps%prec10_col
  lfc              => cps%lfc
  fsat             => cws%fsat
  is_cwd           => decomp_cascade_con%is_cwd
  decomp_cpools_vr => ccs%decomp_cpools_vr
  cpool_col        => ccs%cpool_col 
 
!  write(*,*) 'size PFTi ', size(pfti)
!  write(*,*) 'PFTi vals ', pfti(1:10)

  !assign local pointers to derived type members (grid-level) 
  forc_rh       => grc%forc_rh
  forc_wind     => grc%forc_wind
  forc_t        => grc%forc_t
  forc_rain     => grc%forc_rain
  forc_snow     => grc%forc_snow   
  latdeg        => grc%latdeg
  forc_hdm      => grc%forc_hdm
  forc_lnfm     => grc%forc_lnfm
  
  !assign local pointers to ecophysiological constants, fzeng
  fsr_pft       => pftcon%fsr_pft
  fd_pft        => pftcon%fd_pft
 
  !pft to column average 
  call p2c(num_soilc, filter_soilc, prec10,  prec10_col)
  call p2c(num_soilc, filter_soilc, prec60,  prec60_col)
  call p2c(num_soilc, filter_soilc,totvegc, totvegc_col)
  call p2c(num_soilc, filter_soilc,leafc, leafc_col)
  call p2c(num_soilc, filter_soilc,cpool, cpool_col)
  call get_curr_date (kyr, kmo, kda, mcsec)
  dayspyr = get_days_per_year()
  ! Get model step size
  dt      = real( get_step_size(), r8 )
  !
  ! On first time-step, just set area burned to zero and exit
  !
  if ( get_nstep() == 0 )then
     do fc = 1,num_soilc
        c = filter_soilc(fc)
        farea_burned(c) = 0._r8
        baf_crop(c)     = 0._r8
        baf_peatf(c)    = 0._r8
        fbac(c)         = 0._r8
        fbac1(c)        = 0._r8
        cropf_col(c)    = 0._r8  ! see clm4_6_00, fzeng, 7 May 2019
     end do
     return
  end if
  !
  ! Calculate fraction of crop (cropf_col) and non-crop and non-bare-soil 
  ! vegetation (lfwt) in vegetated column
  !
  do fc = 1,num_soilc
     c = filter_soilc(fc)
     cropf_col(c) = 0._r8 
     lfwt(c)      = 0._r8   
  end do
  do pi = 1,max_pft_per_col
     do fc = 1,num_soilc
        c = filter_soilc(fc)
        if (pi <=  npfts(c)) then
           p = pfti(c) + pi - 1
           ! For crop veg types
           if( ivt(p) > nc4_grass2 )then
               cropf_col(c) = cropf_col(c) + wtcol(p)
           end if
           ! For natural vegetation (non-crop)
           if( ivt(p) >= ndllf_evr_tmp_tree .and. ivt(p) <= nc4_grass2 )then
               lfwt(c) = lfwt(c) + wtcol(p)
           end if
        end if
     end do
  end do
  ! 
  ! Calculate crop fuel   
  !
  do fc = 1,num_soilc
     c = filter_soilc(fc)
     fuelc_crop(c)=0._r8
  end do
  do pi = 1,max_pft_per_col
     do fc = 1,num_soilc
        c = filter_soilc(fc)
        if (pi <=  npfts(c)) then
           p = pfti(c) + pi - 1
           ! For crop PFT's
           if( ivt(p) > nc4_grass2 .and. wtcol(p) > 0._r8 .and. leafc_col(c) > 0._r8 )then
              fuelc_crop(c)=fuelc_crop(c) + (leafc(p) + leafc_storage(p) + &
                            leafc_xfer(p))*wtcol(p)/cropf_col(c)     + &
                            totlitc(c)*leafc(p)/leafc_col(c)*wtcol(p)/cropf_col(c)
           end if
        end if
     end do
  end do          
  !   
  ! Calculate noncrop column variables
  !
  do fc = 1,num_soilc
     c = filter_soilc(fc)
     fsr_col(c)   = 0._r8
     fd_col(c)    = 0._r8
     rootc_col(c) = 0._r8
     lgdp_col(c)  = 0._r8
     lgdp1_col(c) = 0._r8
     lpop_col(c)  = 0._r8
     btran_col(c) = 0._r8
     wtlf(c)      = 0._r8
     trotr1_col(c)= 0._r8 ! see clm4_6_00 and Li et al., BG 2014, fzeng, 7 May 2019
     trotr2_col(c)= 0._r8 ! see clm4_6_00 and Li et al., BG 2014, fzeng, 7 May 2019
     if (fpftdyn /= ' ') then    !true when landuse data is used
        dtrotr_col(c)=0._r8
     end if
  end do
  do pi = 1,max_pft_per_col
     do fc = 1,num_soilc
        c = filter_soilc(fc)
        g = cgridcell(c)
        if (pi <=  npfts(c)) then
           p = pfti(c) + pi - 1
           ! For non-crop -- natural vegetation and bare-soil
           if( ivt(p) .lt. nc3crop .and. cropf_col(c) .lt. 1.0_r8 )then
              if( .not. shr_infnan_isnan(btran2(p)) .and. btran2(p) .le. 1._r8 )then
                 btran_col(c) = btran_col(c)+btran2(p)*wtcol(p)
                 wtlf(c)      = wtlf(c)+wtcol(p)
              end if
              if( ivt(p) == nbrdlf_evr_trp_tree .and. wtcol(p) .gt. 0._r8 )then
!                trotr1_col(c)=trotr1_col(c)+wtcol(p)*cwtgcell(c)  
                 trotr1_col(c)=trotr1_col(c)+wtcol(p) ! wtcol is already the pft weight in the gridcell/tile in Catchment-CN4.5, fzeng, 12 July 2019  
              end if
              if( ivt(p) == nbrdlf_dcd_trp_tree .and. wtcol(p) .gt. 0._r8 )then
!                trotr2_col(c)=trotr2_col(c)+wtcol(p)*cwtgcell(c)
                 trotr2_col(c)=trotr2_col(c)+wtcol(p) ! wtcol is already the pft weight in the gridcell/tile in Catchment-CN4.5, fzeng, 12 July 2019  
              end if
              if ( fpftdyn /= ' ' ) then    !true when landuse data is used
                 if( ivt(p) == nbrdlf_evr_trp_tree .or. ivt(p) == nbrdlf_dcd_trp_tree )then
                    if(lfpftd(p).gt.0._r8)then
                       dtrotr_col(c)=dtrotr_col(c)+lfpftd(p)*cwtgcell(c)
                    end if
                 end if
              end if
              if (wtcol(p) > 0.) then ! Exclude the PFTs that do not exist in the column to avoid adding nan. F Zeng, 23 June 2017 
                 rootc_col(c) = rootc_col(c) + (frootc(p) + frootc_storage(p) + &
                                frootc_xfer(p) + deadcrootc(p) +                &
                                deadcrootc_storage(p) + deadcrootc_xfer(p) +    &
                                livecrootc(p)+livecrootc_storage(p) +           &
                                livecrootc_xfer(p))*wtcol(p)
              end if

              fsr_col(c) = fsr_col(c) + fsr_pft(ivt(p))*wtcol(p)/(1.0_r8-cropf_col(c))
    
              if( lfwt(c) .ne. 0.0_r8 )then    
                 hdmlf=forc_hdm(g)
              
                 ! all these constants are in Li et al. BG (2012a,b;2013)

                 if( hdmlf .gt. 0.1_r8 )then            
                    ! For NOT bare-soil
                    if( ivt(p) .ne. noveg )then
                       ! For shrub and grass (crop already excluded above)
                       if( ivt(p) .ge. nbrdlf_evr_shrub )then      !for shurb and grass
                          lgdp_col(c)  = lgdp_col(c) + (0.1_r8 + 0.9_r8*    &
                                         exp(-1._r8*SHR_CONST_PI* &
                                         (gdp_lf(c)/8._r8)**0.5_r8))*wtcol(p) &
                                         /(1.0_r8 - cropf_col(c))
                          lgdp1_col(c) = lgdp1_col(c) + (0.2_r8 + 0.8_r8*   &
                                         exp(-1._r8*SHR_CONST_PI* &
                                         (gdp_lf(c)/7._r8)))*wtcol(p)/lfwt(c)
                          lpop_col(c)  = lpop_col(c) + (0.2_r8 + 0.8_r8*    &
                                         exp(-1._r8*SHR_CONST_PI* &
                                         (hdmlf/450._r8)**0.5_r8))*wtcol(p)/lfwt(c)
                       else   ! for trees
                          if( gdp_lf(c) .gt. 20._r8 )then
                             lgdp_col(c)  =lgdp_col(c)+0.2_r8*wtcol(p)/(1.0_r8 - cropf_col(c))    ! GDP 0.2
                          else    
                             lgdp_col(c) = lgdp_col(c)+wtcol(p)/(1.0_r8 - cropf_col(c))
                          end if
                          if( gdp_lf(c) .gt. 20._r8 )then   
                             lgdp1_col(c) = lgdp1_col(c)+0.62_r8*wtcol(p)/lfwt(c)
                          else
                             if( gdp_lf(c) .gt. 8._r8 ) then
                                lgdp1_col(c)=lgdp1_col(c)+0.83_r8*wtcol(p)/lfwt(c)
                             else
                                lgdp1_col(c)=lgdp1_col(c)+wtcol(p)/lfwt(c)
                             end if
                          end if
                          lpop_col(c) = lpop_col(c) + (0.4_r8 + 0.6_r8*    &
                                        exp(-1._r8*SHR_CONST_PI* &
                                        (hdmlf/125._r8)))*wtcol(p)/lfwt(c) 
                       end if
                    end if
                 else
                    lgdp_col(c)  = lgdp_col(c)+wtcol(p)/(1.0_r8 - cropf_col(c))
                    lgdp1_col(c) = lgdp1_col(c)+wtcol(p)/lfwt(c)
                    lpop_col(c)  = lpop_col(c)+wtcol(p)/lfwt(c)
                 end if   
              end if
           
              fd_col(c) = fd_col(c) + fd_pft(ivt(p))*wtcol(p)*secsphr/(1.0_r8-cropf_col(c))         
           end if                  
        end if
     end do
  end do

  if (fpftdyn /= ' ') then    !true when landuse data is used
     do fc = 1,num_soilc
        c = filter_soilc(fc)
        if( dtrotr_col(c) .gt. 0._r8 )then
           if( kmo == 1 .and. kda == 1 .and. mcsec == 0)then
              lfc(c) = 0._r8
           end if
           if( kmo == 1 .and. kda == 1 .and. mcsec == dt)then
              lfc(c) = dtrotr_col(c)*dayspyr*secspday/dt
           end if
        else
           lfc(c)=0._r8
        end if
     end do
  end if
  !
  ! calculate burned area fraction in cropland
  !
  do fc = 1,num_soilc
     c = filter_soilc(fc)
     baf_crop(c)=0._r8
  end do

  do fp = 1,num_soilp
     p = filter_soilp(fp)  
     if( kmo == 1 .and. kda == 1 .and. mcsec == 0 )then
         burndate(p) = 10000 ! init. value; actual range [0 365]
     end if
  end do
 
  do pi = 1,max_pft_per_col
     do fc = 1,num_soilc
        c = filter_soilc(fc)
        g= cgridcell(c)
        hdmlf=forc_hdm(g)
        if (pi <=  npfts(c)) then
           p = pfti(c) + pi - 1
           ! For crop
           if( forc_t(g) .ge. SHR_CONST_TKFRZ .and. ivt(p) .gt. nc4_grass2 .and.  &
              kmo == abm_lf(c) .and. forc_rain(g)+forc_snow(g) .eq. 0._r8  .and.  &
              burndate(p) >= 999 .and. wtcol(p) .gt. 0._r8 )then ! catch  crop burn time
              ! calculate human density impact on ag. fire
              fhd = 0.04_r8+0.96_r8*exp(-1._r8*SHR_CONST_PI*(hdmlf/350._r8)**0.5_r8)
              ! calculate impact of GDP on ag. fire
              fgdp = 0.01_r8+0.99_r8*exp(-1._r8*SHR_CONST_PI*(gdp_lf(c)/10._r8))
              ! calculate burned area
              fb   = max(0.0_r8,min(1.0_r8,(fuelc_crop(c)-lfuel)/(ufuel-lfuel)))
              ! crop fire only for generic crop types at this time
              ! managed crops are treated as grasses if crop model is turned on
              baf_crop(c) = baf_crop(c) + cropfire_a1/secsphr*fb*fhd*fgdp*wtcol(p)   ! removed *dt, see clm4_6_00, fzeng, 7 May 2019, so baf_crop is per second instead of per time step now
              if( fb*fhd*fgdp*wtcol(p) .gt. 0._r8)then
                 burndate(p)=kda
              end if
           end if
        end if
     end do
  end do
  !
  ! calculate peatland fire
  !
  do fc = 1, num_soilc
     c = filter_soilc(fc)
     g= cgridcell(c)
     if(grc%latdeg(g).lt.borealat )then
        baf_peatf(c) = non_boreal_peatfire_c/secsphr*max(0._r8, &   ! removed *dt, see clm4_6_00, fzeng, 7 May 2019
                       min(1._r8,(4.0_r8-prec60_col(c)*secspday)/ &
                       4.0_r8))**2*peatf_lf(c)*(1._r8-fsat(c))
     else
        baf_peatf(c) = boreal_peatfire_c/secsphr*exp(-SHR_CONST_PI*(max(wf2(c),0._r8)/0.3_r8))* &  ! removed *dt, see clm4_6_00, fzeng, 7 May 2019
        max(0._r8,min(1._r8,(tsoi17(c)-SHR_CONST_TKFRZ)/10._r8))*peatf_lf(c)* &
        (1._r8-fsat(c))
     end if
  end do
  !
  ! calculate other fires
  !

  ! Set the number of timesteps for e-folding.
  ! When the simulation has run fewer than this number of steps,
  ! re-scale the e-folding time to get a stable early estimate.
  
  ! find which pool is the cwd pool
  i_cwd = 0
  do l = 1, ndecomp_pools
     if ( is_cwd(l) ) then
        i_cwd = l
     endif
  end do
 
  !
  ! begin column loop to calculate fractional area affected by fire
  !
  fuelc = 0.  ! added to avoid floating overflow, fzeng, 31 July 2019
  do fc = 1, num_soilc
     c = filter_soilc(fc)
     g = cgridcell(c)
     hdmlf=forc_hdm(g)
     nfire(c) = 0._r8  ! This is done in CNInitMod.F90 in clm4_6_00, fzeng, 9 May 2019
     if( cropf_col(c) .lt. 1.0 )then
       pft_4_flag = .false.                                                                                                                           
       do pi = 1,max_pft_per_col ! this loop can be avoided by directly targeting PFT 4                                                 
                                 ! kept for now to be able to extend modifications to other PFTs;jkolassa Jun 2020                      
          p = pfti(c) + pi - 1                                                                                                          
          if ((ivt(p)==nbrdlf_evr_trp_tree).and.(wtcol(p).gt.0.)) then                                                                
              pft_4_flag = .true.                                                                                         
           end if                                                                                                                        
        end do  
       ! see clm4_6_00 and Li et al., BG 2014, fzeng, 7 May 2019
!      if (trotr1_col(c)+trotr2_col(c)>0.6_r8) then
       ! Allow trotr2 (i.e. broadleaf deciduous tropical trees in ESA)
       ! in Africa to burn, because this type is classified as woody savanna 
       ! and savanna in MODIS land cover which is more consistent with
       ! CLM4.5CN tree and grass fractions in Africa, fzeng, 12 July 2019; 
       ! additional modification by jkolassa (pft_4_flag=.false.) enables
       ! fires for PFT 4 again in attempt to improve fire carbon emissions over Africa  
       ! changes (pft_4_flag==.false.) to (.not. pft_4_flag). W.J.
       if ((trotr1_col(c)>0.6_r8) .and. ( .not. pft_4_flag )) then
            farea_burned(c)=min(1.0_r8,baf_crop(c)+baf_peatf(c))
         else     
           fuelc(c) = totlitc(c)+totvegc_col(c)-rootc_col(c)-fuelc_crop(c)*cropf_col(c)        
           do j = 1, nlevdecomp  
              fuelc(c) = fuelc(c)+decomp_cpools_vr(c,j,i_cwd) * dzsoi_decomp(j)
           end do
           fuelc(c) = fuelc(c)/(1._r8-cropf_col(c))
           fb       = max(0.0_r8,min(1.0_r8,(fuelc(c)-lfuel)/(ufuel-lfuel)))
           m        = max(0._r8,wf(c))
           fire_m = 0.     ! jkolassa 09/2020: made fire combustability computation PFT-dependent in order to be able to control combustibility sensitivity to soil moisture separately for each PFT
           do pi = 1,max_pft_per_col
              p = pfti(c) + pi - 1
              fire_m_tmp = exp(-fire_m_fac(ivt(p))*(m/0.69_r8)**2)*(1.0_r8 - max(0._r8, &
                      min(1._r8,(forc_rh(g)-30._r8)/(80._r8-30._r8))))*  &
                      min(1._r8,exp(SHR_CONST_PI*(forc_t(g)-SHR_CONST_TKFRZ)/10._r8))
              fire_m = fire_m + (fire_m_tmp*wtcol(p))
           end do
           lh       = 0.0035_r8*6.8_r8*hdmlf**(0.43_r8)/30._r8/24._r8
           fs       = 1._r8-(0.01_r8+0.98_r8*exp(-0.025_r8*hdmlf))
!          ig       = (lh+forc_lnfm(g)/24/(5.16_r8+2.16_r8*cos(3._r8*grc%latdeg(g)))*0.25_r8)*(1._r8-fs)*(1._r8-cropf_col(c))                         ! There is a bug here. The input of "cos" must be in radians!! fzeng, 2 Aug 2019              
           ig       = (lh+forc_lnfm(g)/24/(5.16_r8+2.16_r8*cos(3._r8*min(60.,abs(grc%latdeg(g)))*(SHR_CONST_PI/180._r8)))*0.25_r8)*(1._r8-fs)*(1._r8-cropf_col(c))  ! our forc_lnfm is counts/km2/day, so divided by 24 to get per hr, fzeng, 7 May 2019      
           nfire(c) = ig/secsphr*fb*fire_m*lgdp_col(c) !fire counts/km2/sec
           Lb_lf    = 1._r8+10.0_r8*(1._r8-EXP(-0.06_r8*forc_wind(g)))
           if ( wtlf(c) > 0.0_r8 )then
!             spread_m = (1.0_r8 - max(0._r8,min(1._r8,(btran_col(c)/wtlf(c)-0.3_r8)/ &
!                        (0.7_r8-0.3_r8))))*(1.0-max(0._r8, &
!                        min(1._r8,(forc_rh(g)-30._r8)/(80._r8-30._r8))))
              spread_m = (1.0_r8 - max(0._r8,min(1._r8,(btran_col(c)/wtlf(c)-0.4_r8)/ &   ! To reduce burned area biases, changed the btran lower threshold from 0.3 to 0.4, 
                         (0.8_r8-0.4_r8))))*(1.0-max(0._r8, &                             ! and the btran upper threshold from 0.7 to 0.8, fzeng, 9 Aug 2019
                         min(1._r8,(forc_rh(g)-30._r8)/(80._r8-30._r8))))           
           else
              spread_m = 0.0_r8
           end if
           farea_burned(c) = min(1._r8,(g0*spread_m*fsr_col(c)* &
                          fd_col(c)/1000._r8)**2*lgdp1_col(c)* &
                          lpop_col(c)*nfire(c)*SHR_CONST_PI*Lb_lf+ &
                          baf_crop(c)+baf_peatf(c))  ! fraction (0-1) per sec
       end if   ! (trotr1_col(c)+trotr2_col(c)>0.6_r8)       
       !
       ! if landuse change data is used, calculate deforestation fires and 
       ! add it in the total of burned area fraction
       !
       if (fpftdyn /= ' ') then    !true when landuse change data is used
          if( trotr1_col(c)+trotr2_col(c) > 0.6_r8 )then
             if(( kmo == 1 .and. kda == 1 .and. mcsec == 0) .or. &
                  dtrotr_col(c) <=0._r8 )then
                fbac1(c)        = 0._r8
                farea_burned(c) = baf_crop(c)+baf_peatf(c)
             else
                cri = (4.0_r8*trotr1_col(c)+1.8_r8*trotr2_col(c))/(trotr1_col(c)+trotr2_col(c))
                cli = (max(0._r8,min(1._r8,(cri-prec60_col(c)*secspday)/cri))**0.5)* &
                      (max(0._r8,min(1._r8,(cri-prec10_col(c)*secspday)/cri))**0.5)* &
                      max(0.0005_r8,min(1._r8,19._r8*dtrotr_col(c)*dayspyr*secspday/dt-0.001_r8))* &
                      max(0._r8,min(1._r8,(0.25_r8-(forc_rain(g)+forc_snow(g))*secsphr)/0.25_r8))
                farea_burned(c) = cli*(cli_scale/secspday)+baf_crop(c)+baf_peatf(c)
                ! burned area out of conversion region due to land use fire
                fbac1(c) = max(0._r8,cli*(cli_scale/secspday) - 2.0_r8*lfc(c)/dt)   
             end if
             ! total burned area out of conversion 
             fbac(c) = fbac1(c)+baf_crop(c)+baf_peatf(c) 
          else
             fbac(c) = farea_burned(c)
          end if
       end if  ! (fpftdyn /= ' ')

     else
        farea_burned(c) = min(1._r8,baf_crop(c)+baf_peatf(c))
     end if  ! ( cropf_col(c) .lt. 1.0 )

     if (use_nofire) then
        ! zero out the fire area if NOFIRE flag is on

        farea_burned(c) = 0._r8
        baf_crop(c)     = 0._r8
        baf_peatf(c)    = 0._r8
        fbac(c)         = 0._r8
        fbac1(c)        = 0._r8
        ! with NOFIRE, tree carbon is still removed in landuse change regions by the
        ! landuse code
     end if

  end do  ! end of column loop

end subroutine CNFireArea
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNFireFluxes
!
! !INTERFACE:
subroutine CNFireFluxes (num_soilc, filter_soilc, num_soilp, filter_soilp)
!
! !DESCRIPTION:
! Fire effects routine for coupled carbon-nitrogen code (CN).
! Relies primarily on estimate of fractional area burned in this
! timestep, from CNFireArea().
!
! Total fire carbon emissions (g C/m2 land area/yr) 
!  =avg(COL_FIRE_CLOSS)*seconds_per_year + avg(SOMC_FIRE)*seconds_per_year + 
!   avg(LF_CONV_CFLUX)*seconds_per_year*min(1.0,avg(LFC2)/dt*seconds_per_year)*0.8
! where dt is the time step size (sec),avg means the temporal average in a year
! seconds_per_year is the number of seconds in a year.
!
! !USES:
   use clmtype
   use pftvarcon, only: nc3crop
   use clm_time_manager, only: get_step_size,get_days_per_year,get_curr_date
   use clm_varpar, only : max_pft_per_col
   use clm_varctl  , only: fpftdyn
   use clm_varcon  , only: secspday
!
! !ARGUMENTS:
   implicit none
   integer, intent(in) :: num_soilc       ! number of soil columns in filter
   integer, intent(in) :: filter_soilc(:) ! filter for soil columns
   integer, intent(in) :: num_soilp       ! number of soil pfts in filter
   integer, intent(in) :: filter_soilp(:) ! filter for soil pfts
!
! !CALLED FROM:
! subroutine CNEcosystemDyn()
!
!
! !LOCAL VARIABLES:
! local pointers to implicit in scalars
!
#if (defined CNDV)
   real(r8), pointer :: nind(:)         ! number of individuals (#/m2)
#endif
   real(r8), pointer :: woody(:)        ! woody lifeform (1=woody, 0=not woody) 
   logical , pointer :: pactive(:)      ! true=>do computations on this pft (see reweightMod for details)
   integer , pointer :: ivt(:)          ! pft vegetation type
   real(r8), pointer :: wtcol(:)        ! pft weight relative to column 
   real(r8), pointer :: latdeg(:)       ! latitude (degrees)
   integer , pointer :: cgridcell(:)    ! gridcell of corresponding column
   integer , pointer :: npfts(:)        ! number of pfts for each column
   integer , pointer :: pfti(:)         ! beginning pft index for each column
   integer , pointer :: pcolumn(:)      ! pft's column index
   real(r8), pointer :: farea_burned(:) ! timestep fractional area burned (proportion)
   real(r8), pointer :: fire_mortality_c_to_cwdc(:,:)              ! C fluxes associated with fire mortality to CWD pool (gC/m3/s)
   real(r8), pointer :: decomp_cpools_vr(:,:,:)    ! (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) c pools
   real(r8), pointer :: decomp_npools_vr(:,:,:)    ! (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
   real(r8), pointer :: fire_mortality_n_to_cwdn(:,:)              ! N fluxes associated with fire mortality to CWD pool (gN/m3/s)
   real(r8), pointer :: lfc(:)          ! conversion area frac. of BET+BDT that haven't burned before
   real(r8), pointer :: lfc2(:)         ! conversion area frac. of BET+BDT that burned this timestep
   real(r8), pointer :: fbac1(:)        ! burned area out of conversion region due to land use fire
   real(r8), pointer :: baf_crop(:)     ! baf for cropland
   real(r8), pointer :: baf_peatf(:)    ! baf for peatlabd
   real(r8), pointer :: leafcmax(:)     ! (gC/m2) ann max leaf C
   real(r8), pointer :: fbac(:)         ! total burned area out of conversion 
   real(r8), pointer :: cropf_col(:)    ! cropland fraction in veg column
   
   real(r8), pointer :: dtrotr_col(:)   ! annual decreased fraction coverage of BET+BDT (0-1) on the gridcell
   real(r8), pointer :: trotr1_col(:)   ! pft weight of BET on the gridcell (0-1)
   real(r8), pointer :: trotr2_col(:)   ! pft weight of BDT on the gridcell (0-1)

   real(r8), pointer :: totsomc(:)            ! (gC/m2) total soil organic matter carbon 
   real(r8), pointer :: somc_fire(:)          ! (gC/m2/s)fire carbon emissions due to peat burning
   
   real(r8), pointer :: leafc(:)              ! (gC/m2) leaf C
   real(r8), pointer :: leafc_storage(:)      ! (gC/m2) leaf C storage
   real(r8), pointer :: leafc_xfer(:)         ! (gC/m2) leaf C transfer
   real(r8), pointer :: livestemc(:)          ! (gC/m2) live stem C
   real(r8), pointer :: livestemc_storage(:)  ! (gC/m2) live stem C storage
   real(r8), pointer :: livestemc_xfer(:)     ! (gC/m2) live stem C transfer

   real(r8), pointer :: deadstemc(:)          ! (gC/m2) dead stem C
   real(r8), pointer :: deadstemc_storage(:)  ! (gC/m2) dead stem C storage
   real(r8), pointer :: deadstemc_xfer(:)     ! (gC/m2) dead stem C transfer
   real(r8), pointer :: frootc(:)             ! (gC/m2) fine root C
   real(r8), pointer :: frootc_storage(:)     ! (gC/m2) fine root C storage
   real(r8), pointer :: frootc_xfer(:)        ! (gC/m2) fine root C transfer
   real(r8), pointer :: deadcrootc(:)         ! (gC/m2) dead coarse root C
   real(r8), pointer :: deadcrootc_storage(:) ! (gC/m2) dead coarse root C storage
   real(r8), pointer :: deadcrootc_xfer(:)    ! (gC/m2) dead coarse root C transfer
   real(r8), pointer :: livecrootc(:)         ! (gC/m2) live coarse root C
   real(r8), pointer :: livecrootc_storage(:) ! (gC/m2) live coarse root C storage
   real(r8), pointer :: livecrootc_xfer(:)    ! (gC/m2) live coarse root C transfer
   real(r8), pointer :: gresp_storage(:)      ! (gC/m2) growth respiration storage
   real(r8), pointer :: gresp_xfer(:)         ! (gC/m2) growth respiration transfer

   real(r8), pointer :: leafn(:)              ! (gN/m2) leaf N 
   real(r8), pointer :: leafn_storage(:)      ! (gN/m2) leaf N storage
   real(r8), pointer :: leafn_xfer(:)         ! (gN/m2) leaf N transfer
   real(r8), pointer :: livestemn(:)          ! (gN/m2) live stem N
   real(r8), pointer :: livestemn_storage(:)  ! (gN/m2) live stem N storage
   real(r8), pointer :: livestemn_xfer(:)     ! (gN/m2) live stem N transfer
   real(r8), pointer :: deadstemn(:)          ! (gN/m2) dead stem N
   real(r8), pointer :: deadstemn_storage(:)  ! (gN/m2) dead stem N storage
   real(r8), pointer :: deadstemn_xfer(:)     ! (gN/m2) dead stem N transfer
   real(r8), pointer :: frootn(:)             ! (gN/m2) fine root N
   real(r8), pointer :: frootn_storage(:)     ! (gN/m2) fine root N storage
   real(r8), pointer :: frootn_xfer(:)        ! (gN/m2) fine root N transfer
   real(r8), pointer :: livecrootn(:)         ! (gN/m2) live coarse root N
   real(r8), pointer :: livecrootn_storage(:) ! (gN/m2) live coarse root N storage
   real(r8), pointer :: livecrootn_xfer(:)    ! (gN/m2) live coarse root N transfer
   real(r8), pointer :: deadcrootn(:)         ! (gN/m2) dead coarse root N
   real(r8), pointer :: deadcrootn_storage(:) ! (gN/m2) dead coarse root N storage
   real(r8), pointer :: deadcrootn_xfer(:)    ! (gN/m2) dead coarse root N transfer
   real(r8), pointer :: retransn(:)           ! (gN/m2) plant pool of retranslocated N

   real(r8), pointer :: m_leafc_to_fire(:)             ! (gC/m2/s) fire C emissions from leafc 
   real(r8), pointer :: m_leafc_storage_to_fire(:)     ! (gC/m2/s) fire C emissions from leafc_storage             
   real(r8), pointer :: m_leafc_xfer_to_fire(:)        ! (gC/m2/s) fire C emissions from leafc_xfer
   real(r8), pointer :: m_livestemc_to_fire(:)         ! (gC/m2/s) fire C emissions from livestemc
   real(r8), pointer :: m_livestemc_storage_to_fire(:) ! (gC/m2/s) fire C emissions from livestemc_storage       
   real(r8), pointer :: m_livestemc_xfer_to_fire(:)    ! (gC/m2/s) fire C emissions from livestemc_xfer
   real(r8), pointer :: m_deadstemc_to_fire(:)         ! (gC/m2/s) fire C emissions from deadstemc_xfer
   real(r8), pointer :: m_deadstemc_storage_to_fire(:) ! (gC/m2/s) fire C emissions from deadstemc_storage         
   real(r8), pointer :: m_deadstemc_xfer_to_fire(:)    ! (gC/m2/s) fire C emissions from deadstemc_xfer
   real(r8), pointer :: m_frootc_to_fire(:)            ! (gC/m2/s) fire C emissions from frootc
   real(r8), pointer :: m_frootc_storage_to_fire(:)    ! (gC/m2/s) fire C emissions from frootc_storage
   real(r8), pointer :: m_frootc_xfer_to_fire(:)       ! (gC/m2/s) fire C emissions from frootc_xfer
   real(r8), pointer :: m_livecrootc_to_fire(:)        ! (gC/m2/s) fire C emissions from livecrootc
   real(r8), pointer :: m_livecrootc_storage_to_fire(:)! (gC/m2/s) fire C emissions from livecrootc_storage     
   real(r8), pointer :: m_livecrootc_xfer_to_fire(:)   ! (gC/m2/s) fire C emissions from livecrootc_xfer
   real(r8), pointer :: m_deadcrootc_to_fire(:)        ! (gC/m2/s) fire C emissions from deadcrootc
   real(r8), pointer :: m_deadcrootc_storage_to_fire(:)! (gC/m2/s) fire C emissions from deadcrootc_storage 
   real(r8), pointer :: m_deadcrootc_xfer_to_fire(:)   ! (gC/m2/s) fire C emissions from deadcrootc_xfer
   real(r8), pointer :: m_gresp_storage_to_fire(:)     ! (gC/m2/s) fire C emissions from gresp_storage 
   real(r8), pointer :: m_gresp_xfer_to_fire(:)        ! (gC/m2/s) fire C emissions from gresp_xfer
   real(r8), pointer :: m_decomp_cpools_to_fire_vr(:,:,:) ! (gC/m3/s) vertically-resolved decomposing C fire loss 

   real(r8), pointer :: m_leafn_to_fire(:)             ! (gN/m2/s) fire N emissions from leafn 
   real(r8), pointer :: m_leafn_storage_to_fire(:)     ! (gN/m2/s) fire N emissions from leafn_storage            
   real(r8), pointer :: m_leafn_xfer_to_fire(:)        ! (gN/m2/s) fire N emissions from leafn_xfer     
   real(r8), pointer :: m_livestemn_to_fire(:)         ! (gN/m2/s) fire N emissions from livestemn 
   real(r8), pointer :: m_livestemn_storage_to_fire(:) ! (gN/m2/s) fire N emissions from livestemn_storage      
   real(r8), pointer :: m_livestemn_xfer_to_fire(:)    ! (gN/m2/s) fire N emissions from livestemn_xfer
   real(r8), pointer :: m_deadstemn_to_fire(:)         ! (gN/m2/s) fire N emissions from deadstemn
   real(r8), pointer :: m_deadstemn_storage_to_fire(:) ! (gN/m2/s) fire N emissions from deadstemn_storage         
   real(r8), pointer :: m_deadstemn_xfer_to_fire(:)    ! (gN/m2/s) fire N emissions from deadstemn_xfer
   real(r8), pointer :: m_frootn_to_fire(:)            ! (gN/m2/s) fire N emissions from frootn
   real(r8), pointer :: m_frootn_storage_to_fire(:)    ! (gN/m2/s) fire N emissions from frootn_storage
   real(r8), pointer :: m_frootn_xfer_to_fire(:)       ! (gN/m2/s) fire N emissions from frootn_xfer
   real(r8), pointer :: m_livecrootn_to_fire(:)        ! (gN/m2/s) fire N emissions from m_livecrootn_to_fire
   real(r8), pointer :: m_livecrootn_storage_to_fire(:)! (gN/m2/s) fire N emissions from livecrootn_storage     
   real(r8), pointer :: m_livecrootn_xfer_to_fire(:)   ! (gN/m2/s) fire N emissions from livecrootn_xfer
   real(r8), pointer :: m_deadcrootn_to_fire(:)        ! (gN/m2/s) fire N emissions from deadcrootn
   real(r8), pointer :: m_deadcrootn_storage_to_fire(:)! (gN/m2/s) fire N emissions from deadcrootn_storage  
   real(r8), pointer :: m_deadcrootn_xfer_to_fire(:)   ! (gN/m2/s) fire N emissions from deadcrootn_xfer
   real(r8), pointer :: m_retransn_to_fire(:)          ! (gN/m2/s) fire N emissions from retransn
   real(r8), pointer :: m_decomp_npools_to_fire_vr(:,:,:)  ! vertically-resolved decomposing N fire loss (gN/m3/s)

!  (gC/m2/s) C transfers from various C pools to litter and cwd pools due to fire mortality
   real(r8), pointer :: m_leafc_to_litter_fire(:)   
   real(r8), pointer :: m_leafc_storage_to_litter_fire(:)                
   real(r8), pointer :: m_leafc_xfer_to_litter_fire(:)  
   real(r8), pointer :: m_livestemc_to_litter_fire(:)    
   real(r8), pointer :: m_livestemc_storage_to_litter_fire(:)        
   real(r8), pointer :: m_livestemc_xfer_to_litter_fire(:) 
   real(r8), pointer :: m_livestemc_to_deadstemc_fire(:)    
   real(r8), pointer :: m_deadstemc_to_litter_fire(:) 
   real(r8), pointer :: m_deadstemc_storage_to_litter_fire(:)           
   real(r8), pointer :: m_deadstemc_xfer_to_litter_fire(:) 
   real(r8), pointer :: m_frootc_to_litter_fire(:)        
   real(r8), pointer :: m_frootc_storage_to_litter_fire(:)  
   real(r8), pointer :: m_frootc_xfer_to_litter_fire(:)
   real(r8), pointer :: m_livecrootc_to_litter_fire(:)    
   real(r8), pointer :: m_livecrootc_storage_to_litter_fire(:)      
   real(r8), pointer :: m_livecrootc_xfer_to_litter_fire(:)
   real(r8), pointer :: m_livecrootc_to_deadcrootc_fire(:)    
   real(r8), pointer :: m_deadcrootc_to_litter_fire(:)        
   real(r8), pointer :: m_deadcrootc_storage_to_litter_fire(:)  
   real(r8), pointer :: m_deadcrootc_xfer_to_litter_fire(:)
   real(r8), pointer :: m_gresp_storage_to_litter_fire(:)      
   real(r8), pointer :: m_gresp_xfer_to_litter_fire(:)    
   real(r8), pointer :: m_c_to_litr_met_fire(:,:)
   real(r8), pointer :: m_c_to_litr_cel_fire(:,:)
   real(r8), pointer :: m_c_to_litr_lig_fire(:,:)

!  (gN/m2/s) N transfers from various C pools to litter and cwd pools due to fire mortality   
   real(r8), pointer :: m_leafn_to_litter_fire(:)   
   real(r8), pointer :: m_leafn_storage_to_litter_fire(:)                
   real(r8), pointer :: m_leafn_xfer_to_litter_fire(:)  
   real(r8), pointer :: m_livestemn_to_litter_fire(:)    
   real(r8), pointer :: m_livestemn_storage_to_litter_fire(:)        
   real(r8), pointer :: m_livestemn_xfer_to_litter_fire(:) 
   real(r8), pointer :: m_livestemn_to_deadstemn_fire(:)    
   real(r8), pointer :: m_deadstemn_to_litter_fire(:) 
   real(r8), pointer :: m_deadstemn_storage_to_litter_fire(:)           
   real(r8), pointer :: m_deadstemn_xfer_to_litter_fire(:) 
   real(r8), pointer :: m_frootn_to_litter_fire(:)        
   real(r8), pointer :: m_frootn_storage_to_litter_fire(:)  
   real(r8), pointer :: m_frootn_xfer_to_litter_fire(:)
   real(r8), pointer :: m_livecrootn_to_litter_fire(:)    
   real(r8), pointer :: m_livecrootn_storage_to_litter_fire(:)      
   real(r8), pointer :: m_livecrootn_xfer_to_litter_fire(:)
   real(r8), pointer :: m_livecrootn_to_deadcrootn_fire(:)    
   real(r8), pointer :: m_deadcrootn_to_litter_fire(:)        
   real(r8), pointer :: m_deadcrootn_storage_to_litter_fire(:)  
   real(r8), pointer :: m_deadcrootn_xfer_to_litter_fire(:)
   real(r8), pointer :: m_retransn_to_litter_fire(:)      
   real(r8), pointer :: m_n_to_litr_met_fire(:,:)
   real(r8), pointer :: m_n_to_litr_cel_fire(:,:)
   real(r8), pointer :: m_n_to_litr_lig_fire(:,:)
  
   logical, pointer  :: is_cwd(:)               ! TRUE => pool is a cwd pool
   logical, pointer  :: is_litter(:)            ! TRUE => pool is a litter pool
   real(r8), pointer :: froot_prof(:,:)         ! (1/m) profile of fine roots
   real(r8), pointer :: croot_prof(:,:)         ! (1/m) profile of coarse roots
   real(r8), pointer :: stem_prof(:,:)          ! (1/m) profile of stems
   real(r8), pointer :: leaf_prof(:,:)          ! (1/m) profile of leaves
   
   real(r8), pointer :: cc_dstem(:)             ! Combustion completeness factor for dead stem (0 to 1)
   real(r8), pointer :: cc_leaf(:)              ! Combustion completeness factor for leaf (0 to 1) 
   real(r8), pointer :: cc_lstem(:)             ! Combustion completeness factor for live stem (0 to 1)
   real(r8), pointer :: cc_other(:)             ! Combustion completeness factor for other plant (0 to 1)
   real(r8), pointer :: fm_leaf(:)              ! Fire-related mortality factor for leaf (0 to 1)
   real(r8), pointer :: fm_lstem(:)             ! Fire-related mortality factor for live stem (0 to 1) 
   real(r8), pointer :: fm_other(:)             ! Fire-related mortality factor for other plant (0 to 1)
   real(r8), pointer :: fm_root(:)              ! Fire-related mortality factor for fine roots (0 to 1)
   real(r8), pointer :: fm_lroot(:)             ! Fire-related mortality factor for live roots (0 to 1) 
   real(r8), pointer :: fm_droot(:)             ! Fire-related mortality factor for dead roots (0 to 1)
   real(r8), pointer :: lf_flab(:)              ! leaf litter labile fraction
   real(r8), pointer :: lf_fcel(:)              ! leaf litter cellulose fraction
   real(r8), pointer :: lf_flig(:)              ! leaf litter lignin fraction
   real(r8), pointer :: fr_flab(:)              ! fine root litter labile fraction
   real(r8), pointer :: fr_fcel(:)              ! fine root litter cellulose fraction
   real(r8), pointer :: fr_flig(:)              ! fine root litter lignin fraction
!
! !OTHER LOCAL VARIABLES:
   integer :: g,c,p,j,l,pi,kyr, kmo, kda, mcsec   ! indices
   integer :: fp,fc                ! filter indices
   real(r8):: f                    ! rate for fire effects (1/s)
   real(r8):: dt                   ! time step variable (s)
   real(r8):: dayspyr              ! days per year
!EOP
!-----------------------------------------------------------------------

   ! assign local pointers

#if (defined CNDV)
   nind                           => pdgvs%nind
#endif
   pcolumn                        =>pft%column
   cgridcell                      =>col%gridcell
   farea_burned                   => cps%farea_burned
   woody                          => pftcon%woody
   fire_mortality_c_to_cwdc       => ccf%fire_mortality_c_to_cwdc
   fire_mortality_n_to_cwdn       => cnf%fire_mortality_n_to_cwdn
   
   lfc                            => cps%lfc
   lfc2                           => cps%lfc2
   fbac1                          => cps%fbac1
   fbac                           => cps%fbac
   baf_crop                       => cps%baf_crop
   baf_peatf                      => cps%baf_peatf
   leafcmax                       => pcs%leafcmax
   latdeg                         =>  grc%latdeg
   wtcol                          =>pft%wtcol   
   pfti                           =>col%pfti 
   cropf_col                      => cps%cropf_col
   
   ivt                            =>pft%itype
   npfts                          =>col%npfts
   
   trotr1_col                     => cps%trotr1_col
   trotr2_col                     => cps%trotr2_col
   dtrotr_col                     => cps%dtrotr_col
   
   
   somc_fire                      => ccf%somc_fire
   totsomc                        => ccs%totsomc
   decomp_cpools_vr               => ccs%decomp_cpools_vr
   decomp_npools_vr               => cns%decomp_npools_vr

   leafc                          => pcs%leafc
   leafc_storage                  => pcs%leafc_storage
   leafc_xfer                     => pcs%leafc_xfer
   livestemc                      => pcs%livestemc
   livestemc_storage              => pcs%livestemc_storage
   livestemc_xfer                 => pcs%livestemc_xfer
   deadstemc                      => pcs%deadstemc
   deadstemc_storage              => pcs%deadstemc_storage
   deadstemc_xfer                 => pcs%deadstemc_xfer
   frootc                         => pcs%frootc
   frootc_storage                 => pcs%frootc_storage
   frootc_xfer                    => pcs%frootc_xfer
   livecrootc                     => pcs%livecrootc
   livecrootc_storage             => pcs%livecrootc_storage
   livecrootc_xfer                => pcs%livecrootc_xfer
   deadcrootc                     => pcs%deadcrootc
   deadcrootc_storage             => pcs%deadcrootc_storage
   deadcrootc_xfer                => pcs%deadcrootc_xfer
   gresp_storage                  => pcs%gresp_storage
   gresp_xfer                     => pcs%gresp_xfer
    
   leafn                          => pns%leafn
   leafn_storage                  => pns%leafn_storage
   leafn_xfer                     => pns%leafn_xfer
   livestemn                      => pns%livestemn
   livestemn_storage              => pns%livestemn_storage
   livestemn_xfer                 => pns%livestemn_xfer
   deadstemn                      => pns%deadstemn
   deadstemn_storage              => pns%deadstemn_storage
   deadstemn_xfer                 => pns%deadstemn_xfer
   frootn                         => pns%frootn
   frootn_storage                 => pns%frootn_storage
   frootn_xfer                    => pns%frootn_xfer
   livecrootn                     => pns%livecrootn
   livecrootn_storage             => pns%livecrootn_storage
   livecrootn_xfer                => pns%livecrootn_xfer
   deadcrootn                     => pns%deadcrootn
   deadcrootn_storage             => pns%deadcrootn_storage
   deadcrootn_xfer                => pns%deadcrootn_xfer
   retransn                       => pns%retransn  
   pactive                        => pft%active

   m_leafc_to_fire                => pcf%m_leafc_to_fire
   m_leafc_storage_to_fire        => pcf%m_leafc_storage_to_fire
   m_leafc_xfer_to_fire           => pcf%m_leafc_xfer_to_fire
   m_livestemc_to_fire            => pcf%m_livestemc_to_fire
   m_livestemc_storage_to_fire    => pcf%m_livestemc_storage_to_fire
   m_livestemc_xfer_to_fire       => pcf%m_livestemc_xfer_to_fire
   m_deadstemc_to_fire            => pcf%m_deadstemc_to_fire
   m_deadstemc_storage_to_fire    => pcf%m_deadstemc_storage_to_fire
   m_deadstemc_xfer_to_fire       => pcf%m_deadstemc_xfer_to_fire
   m_frootc_to_fire               => pcf%m_frootc_to_fire
   m_frootc_storage_to_fire       => pcf%m_frootc_storage_to_fire
   m_frootc_xfer_to_fire          => pcf%m_frootc_xfer_to_fire
   m_livecrootc_to_fire           => pcf%m_livecrootc_to_fire
   m_livecrootc_storage_to_fire   => pcf%m_livecrootc_storage_to_fire
   m_livecrootc_xfer_to_fire      => pcf%m_livecrootc_xfer_to_fire
   m_deadcrootc_to_fire           => pcf%m_deadcrootc_to_fire
   m_deadcrootc_storage_to_fire   => pcf%m_deadcrootc_storage_to_fire
   m_deadcrootc_xfer_to_fire      => pcf%m_deadcrootc_xfer_to_fire
   m_gresp_storage_to_fire        => pcf%m_gresp_storage_to_fire
   m_gresp_xfer_to_fire           => pcf%m_gresp_xfer_to_fire

   m_leafn_to_fire                => pnf%m_leafn_to_fire
   m_leafn_storage_to_fire        => pnf%m_leafn_storage_to_fire
   m_leafn_xfer_to_fire           => pnf%m_leafn_xfer_to_fire
   m_livestemn_to_fire            => pnf%m_livestemn_to_fire
   m_livestemn_storage_to_fire    => pnf%m_livestemn_storage_to_fire
   m_livestemn_xfer_to_fire       => pnf%m_livestemn_xfer_to_fire
   m_deadstemn_to_fire            => pnf%m_deadstemn_to_fire
   m_deadstemn_storage_to_fire    => pnf%m_deadstemn_storage_to_fire
   m_deadstemn_xfer_to_fire       => pnf%m_deadstemn_xfer_to_fire
   m_frootn_to_fire               => pnf%m_frootn_to_fire
   m_frootn_storage_to_fire       => pnf%m_frootn_storage_to_fire
   m_frootn_xfer_to_fire          => pnf%m_frootn_xfer_to_fire
   m_livecrootn_to_fire           => pnf%m_livecrootn_to_fire
   m_livecrootn_storage_to_fire   => pnf%m_livecrootn_storage_to_fire
   m_livecrootn_xfer_to_fire      => pnf%m_livecrootn_xfer_to_fire
   m_deadcrootn_to_fire           => pnf%m_deadcrootn_to_fire
   m_deadcrootn_storage_to_fire   => pnf%m_deadcrootn_storage_to_fire
   m_deadcrootn_xfer_to_fire      => pnf%m_deadcrootn_xfer_to_fire
   m_retransn_to_fire             => pnf%m_retransn_to_fire

   m_leafc_to_litter_fire                => pcf%m_leafc_to_litter_fire
   m_leafc_storage_to_litter_fire        => pcf%m_leafc_storage_to_litter_fire
   m_leafc_xfer_to_litter_fire           => pcf%m_leafc_xfer_to_litter_fire
   m_livestemc_to_litter_fire            => pcf%m_livestemc_to_litter_fire
   m_livestemc_storage_to_litter_fire    => pcf%m_livestemc_storage_to_litter_fire
   m_livestemc_xfer_to_litter_fire       => pcf%m_livestemc_xfer_to_litter_fire
   m_livestemc_to_deadstemc_fire         => pcf%m_livestemc_to_deadstemc_fire
   m_deadstemc_to_litter_fire            => pcf%m_deadstemc_to_litter_fire
   m_deadstemc_storage_to_litter_fire    => pcf%m_deadstemc_storage_to_litter_fire
   m_deadstemc_xfer_to_litter_fire       => pcf%m_deadstemc_xfer_to_litter_fire
   m_frootc_to_litter_fire               => pcf%m_frootc_to_litter_fire
   m_frootc_storage_to_litter_fire       => pcf%m_frootc_storage_to_litter_fire
   m_frootc_xfer_to_litter_fire          => pcf%m_frootc_xfer_to_litter_fire
   m_livecrootc_to_litter_fire           => pcf%m_livecrootc_to_litter_fire
   m_livecrootc_storage_to_litter_fire   => pcf%m_livecrootc_storage_to_litter_fire
   m_livecrootc_xfer_to_litter_fire      => pcf%m_livecrootc_xfer_to_litter_fire
   m_livecrootc_to_deadcrootc_fire       => pcf%m_livecrootc_to_deadcrootc_fire
   m_deadcrootc_to_litter_fire           => pcf%m_deadcrootc_to_litter_fire
   m_deadcrootc_storage_to_litter_fire   => pcf%m_deadcrootc_storage_to_litter_fire
   m_deadcrootc_xfer_to_litter_fire      => pcf%m_deadcrootc_xfer_to_litter_fire
   m_gresp_storage_to_litter_fire        => pcf%m_gresp_storage_to_litter_fire
   m_gresp_xfer_to_litter_fire           => pcf%m_gresp_xfer_to_litter_fire
   m_decomp_cpools_to_fire_vr            => ccf%m_decomp_cpools_to_fire_vr
   m_c_to_litr_met_fire                  => ccf%m_c_to_litr_met_fire
   m_c_to_litr_cel_fire                  => ccf%m_c_to_litr_cel_fire
   m_c_to_litr_lig_fire                  => ccf%m_c_to_litr_lig_fire

   m_leafn_to_litter_fire                => pnf%m_leafn_to_litter_fire
   m_leafn_storage_to_litter_fire        => pnf%m_leafn_storage_to_litter_fire
   m_leafn_xfer_to_litter_fire           => pnf%m_leafn_xfer_to_litter_fire
   m_livestemn_to_litter_fire            => pnf%m_livestemn_to_litter_fire
   m_livestemn_storage_to_litter_fire    => pnf%m_livestemn_storage_to_litter_fire
   m_livestemn_xfer_to_litter_fire       => pnf%m_livestemn_xfer_to_litter_fire
   m_livestemn_to_deadstemn_fire         => pnf%m_livestemn_to_deadstemn_fire
   m_deadstemn_to_litter_fire            => pnf%m_deadstemn_to_litter_fire
   m_deadstemn_storage_to_litter_fire    => pnf%m_deadstemn_storage_to_litter_fire
   m_deadstemn_xfer_to_litter_fire       =>pnf%m_deadstemn_xfer_to_litter_fire
   m_frootn_to_litter_fire               => pnf%m_frootn_to_litter_fire
   m_frootn_storage_to_litter_fire       => pnf%m_frootn_storage_to_litter_fire
   m_frootn_xfer_to_litter_fire          => pnf%m_frootn_xfer_to_litter_fire
   m_livecrootn_to_litter_fire           => pnf%m_livecrootn_to_litter_fire
   m_livecrootn_storage_to_litter_fire   => pnf%m_livecrootn_storage_to_litter_fire
   m_livecrootn_xfer_to_litter_fire      => pnf%m_livecrootn_xfer_to_litter_fire
   m_livecrootn_to_deadcrootn_fire       => pnf%m_livecrootn_to_deadcrootn_fire
   m_deadcrootn_to_litter_fire           => pnf%m_deadcrootn_to_litter_fire
   m_deadcrootn_storage_to_litter_fire   => pnf%m_deadcrootn_storage_to_litter_fire
   m_deadcrootn_xfer_to_litter_fire      => pnf%m_deadcrootn_xfer_to_litter_fire
   m_retransn_to_litter_fire             => pnf%m_retransn_to_litter_fire
   m_decomp_npools_to_fire_vr            => cnf%m_decomp_npools_to_fire_vr
   m_n_to_litr_met_fire                  => cnf%m_n_to_litr_met_fire
   m_n_to_litr_cel_fire                  => cnf%m_n_to_litr_cel_fire
   m_n_to_litr_lig_fire                  => cnf%m_n_to_litr_lig_fire
   
   is_cwd                                => decomp_cascade_con%is_cwd
   is_litter                             => decomp_cascade_con%is_litter
   croot_prof                            => pps%croot_prof
   stem_prof                             => pps%stem_prof
   froot_prof                            => pps%froot_prof
   leaf_prof                             => pps%leaf_prof
   
   cc_dstem                              => pftcon%cc_dstem 
   cc_leaf                               => pftcon%cc_leaf  
   cc_lstem                              => pftcon%cc_lstem 
   cc_other                              => pftcon%cc_other 
   fm_leaf                               => pftcon%fm_leaf  
   fm_lstem                              => pftcon%fm_lstem 
   fm_other                              => pftcon%fm_other 
   fm_root                               => pftcon%fm_root  
   fm_lroot                              => pftcon%fm_lroot 
   fm_droot                              => pftcon%fm_droot
   lf_flab                               => pftcon%lf_flab
   lf_fcel                               => pftcon%lf_fcel
   lf_flig                               => pftcon%lf_flig
   fr_flab                               => pftcon%fr_flab
   fr_fcel                               => pftcon%fr_fcel
   fr_flig                               => pftcon%fr_flig

   ! Get model step size
   ! calculate burned area fraction per sec
   dt = real( get_step_size(), r8 )

   dayspyr = get_days_per_year()
   !
   ! pft loop
   !
   do fp = 1,num_soilp
      p = filter_soilp(fp)
      c = pcolumn(p)

      ! For non-crop (bare-soil and natural vegetation)
      ! modified according to clm4_6_00, fzeng, 7 May 2019
      if( ivt(p) .lt. nc3crop .and. cropf_col(c) < 1.0_r8)then
         if (fpftdyn /= ' ') then    !true when landuse data is used
            f = (fbac(c)-baf_crop(c))/(1.0_r8-cropf_col(c))
         else
            f = (farea_burned(c)-baf_crop(c))/(1.0_r8-cropf_col(c))        
         end if
      else
         ! For crops
           if(cropf_col(c) > 0._r8)then
             f = baf_crop(c) /cropf_col(c)
           else
             f = 0._r8
           end if
      end if
      
      ! apply this rate to the pft state variables to get flux rates
      ! biomass burning
      ! carbon fluxes
      m_leafc_to_fire(p)               =  leafc(p)              * f * cc_leaf(ivt(p))
      m_leafc_storage_to_fire(p)       =  leafc_storage(p)      * f * cc_other(ivt(p))
      m_leafc_xfer_to_fire(p)          =  leafc_xfer(p)         * f * cc_other(ivt(p))
      m_livestemc_to_fire(p)           =  livestemc(p)          * f * cc_lstem(ivt(p))
      m_livestemc_storage_to_fire(p)   =  livestemc_storage(p)  * f * cc_other(ivt(p))
      m_livestemc_xfer_to_fire(p)      =  livestemc_xfer(p)     * f * cc_other(ivt(p))
      m_deadstemc_to_fire(p)           =  deadstemc(p)          * f * cc_dstem(ivt(p))
      m_deadstemc_storage_to_fire(p)   =  deadstemc_storage(p)  * f * cc_other(ivt(p))
      m_deadstemc_xfer_to_fire(p)      =  deadstemc_xfer(p)     * f * cc_other(ivt(p))
      m_frootc_to_fire(p)              =  frootc(p)             * f * 0._r8
      m_frootc_storage_to_fire(p)      =  frootc_storage(p)     * f * cc_other(ivt(p)) 
      m_frootc_xfer_to_fire(p)         =  frootc_xfer(p)        * f * cc_other(ivt(p))
      m_livecrootc_to_fire(p)          =  livecrootc(p)         * f * 0._r8
      m_livecrootc_storage_to_fire(p)  =  livecrootc_storage(p) * f * cc_other(ivt(p)) 
      m_livecrootc_xfer_to_fire(p)     =  livecrootc_xfer(p)    * f * cc_other(ivt(p)) 
      m_deadcrootc_to_fire(p)          =  deadcrootc(p)         * f * 0._r8
      m_deadcrootc_storage_to_fire(p)  =  deadcrootc_storage(p) * f*  cc_other(ivt(p)) 
      m_deadcrootc_xfer_to_fire(p)     =  deadcrootc_xfer(p)    * f * cc_other(ivt(p)) 
      m_gresp_storage_to_fire(p)       =  gresp_storage(p)      * f * cc_other(ivt(p))
      m_gresp_xfer_to_fire(p)          =  gresp_xfer(p)         * f * cc_other(ivt(p))


      ! nitrogen fluxes
      m_leafn_to_fire(p)               =  leafn(p)              * f * cc_leaf(ivt(p))
      m_leafn_storage_to_fire(p)       =  leafn_storage(p)      * f * cc_other(ivt(p))
      m_leafn_xfer_to_fire(p)          =  leafn_xfer(p)         * f * cc_other(ivt(p))
      m_livestemn_to_fire(p)           =  livestemn(p)          * f * cc_lstem(ivt(p))
      m_livestemn_storage_to_fire(p)   =  livestemn_storage(p)  * f * cc_other(ivt(p))
      m_livestemn_xfer_to_fire(p)      =  livestemn_xfer(p)     * f * cc_other(ivt(p))
      m_deadstemn_to_fire(p)           =  deadstemn(p)          * f * cc_dstem(ivt(p))
      m_deadstemn_storage_to_fire(p)   =  deadstemn_storage(p)  * f * cc_other(ivt(p))
      m_deadstemn_xfer_to_fire(p)      =  deadstemn_xfer(p)     * f * cc_other(ivt(p))
      m_frootn_to_fire(p)              =  frootn(p)             * f * 0._r8
      m_frootn_storage_to_fire(p)      =  frootn_storage(p)     * f * cc_other(ivt(p))
      m_frootn_xfer_to_fire(p)         =  frootn_xfer(p)        * f * cc_other(ivt(p))
      m_livecrootn_to_fire(p)          =  livecrootn(p)         * f * 0._r8 
      m_livecrootn_storage_to_fire(p)  =  livecrootn_storage(p) * f * cc_other(ivt(p)) 
      m_livecrootn_xfer_to_fire(p)     =  livecrootn_xfer(p)    * f * cc_other(ivt(p))
      m_deadcrootn_to_fire(p)          =  deadcrootn(p)         * f * 0._r8
      m_deadcrootn_xfer_to_fire(p)     =  deadcrootn_xfer(p)    * f * cc_other(ivt(p)) 
      m_deadcrootn_storage_to_fire(p)  =  deadcrootn_storage(p) * f * cc_other(ivt(p))
      m_retransn_to_fire(p)            =  retransn(p)           * f * cc_other(ivt(p))
      
      ! mortality due to fire
      ! carbon bool 
      m_leafc_to_litter_fire(p)                   =  leafc(p) * f * &
                                                     (1._r8 - cc_leaf(ivt(p))) * &
                                                     fm_leaf(ivt(p))
      m_leafc_storage_to_litter_fire(p)           =  leafc_storage(p) * f * &
                                                     (1._r8 - cc_other(ivt(p))) * &
                                                     fm_other(ivt(p))
      m_leafc_xfer_to_litter_fire(p)              =  leafc_xfer(p) * f * &
                                                     (1._r8 - cc_other(ivt(p))) * &
                                                     fm_other(ivt(p))
      m_livestemc_to_litter_fire(p)               =  livestemc(p) * f * &
                                                     (1._r8 - cc_lstem(ivt(p))) * &
                                                     fm_droot(ivt(p))    
      m_livestemc_storage_to_litter_fire(p)       =  livestemc_storage(p) * f * &
                                                     (1._r8 - cc_other(ivt(p))) * &
                                                     fm_other(ivt(p))
      m_livestemc_xfer_to_litter_fire(p)          =  livestemc_xfer(p) * f * &
                                                     (1._r8 - cc_other(ivt(p))) * &
                                                     fm_other(ivt(p)) 
      m_livestemc_to_deadstemc_fire(p)            =  livestemc(p) * f * &
                                                     (1._r8 - cc_lstem(ivt(p))) * &
                                                     (fm_lstem(ivt(p))-fm_droot(ivt(p)))
      m_deadstemc_to_litter_fire(p)               =  deadstemc(p) * f * &
                                                     (1._r8 - cc_dstem(ivt(p))) * &
                                                     fm_droot(ivt(p))    
      m_deadstemc_storage_to_litter_fire(p)       =  deadstemc_storage(p) * f * &
                                                     (1._r8 - cc_other(ivt(p))) * &
                                                     fm_other(ivt(p))
      m_deadstemc_xfer_to_litter_fire(p)          =  deadstemc_xfer(p) * f * &
                                                     (1._r8 - cc_other(ivt(p))) * &
                                                     fm_other(ivt(p))
      m_frootc_to_litter_fire(p)                  =  frootc(p)             * f * &
                                                     fm_root(ivt(p))
      m_frootc_storage_to_litter_fire(p)          =  frootc_storage(p)     * f * &
                                                     fm_other(ivt(p))
      m_frootc_xfer_to_litter_fire(p)             =  frootc_xfer(p)        * f * &
                                                     fm_other(ivt(p))
      m_livecrootc_to_litter_fire(p)              =  livecrootc(p)         * f * &
                                                     fm_droot(ivt(p))
      m_livecrootc_storage_to_litter_fire(p)      =  livecrootc_storage(p) * f * &
                                                     fm_other(ivt(p)) 
      m_livecrootc_xfer_to_litter_fire(p)         =  livecrootc_xfer(p)    * f * &
                                                     fm_other(ivt(p)) 
      m_livecrootc_to_deadcrootc_fire(p)          =  livecrootc(p)         * f * &
                                                     (fm_lroot(ivt(p))-fm_droot(ivt(p)))
      m_deadcrootc_to_litter_fire(p)              =  deadcrootc(p)         * f * &
                                                     fm_droot(ivt(p))
      m_deadcrootc_storage_to_litter_fire(p)      =  deadcrootc_storage(p) * f * &
                                                     fm_other(ivt(p))
      m_deadcrootc_xfer_to_litter_fire(p)         =  deadcrootc_xfer(p)    * f * &
                                                     fm_other(ivt(p))      
      m_gresp_storage_to_litter_fire(p)           =  gresp_storage(p) * f * &
                                                     (1._r8 - cc_other(ivt(p))) * &
                                                     fm_other(ivt(p))  
      m_gresp_xfer_to_litter_fire(p)              =  gresp_xfer(p) * f * &
                                                     (1._r8 - cc_other(ivt(p))) * &
                                                     fm_other(ivt(p)) 
   
       
                   
      ! nitrogen pools    
      m_leafn_to_litter_fire(p)                  =  leafn(p) * f * &
                                                    (1._r8 - cc_leaf(ivt(p))) * &
                                                    fm_leaf(ivt(p))
      m_leafn_storage_to_litter_fire(p)          =  leafn_storage(p) * f * &
                                                    (1._r8 - cc_other(ivt(p))) * &
                                                    fm_other(ivt(p))  
      m_leafn_xfer_to_litter_fire(p)             =  leafn_xfer(p) * f * &
                                                    (1._r8 - cc_other(ivt(p))) * &
                                                    fm_other(ivt(p))
      m_livestemn_to_litter_fire(p)              =  livestemn(p) * f * &
                                                    (1._r8 - cc_lstem(ivt(p))) * &
                                                    fm_droot(ivt(p))
      m_livestemn_storage_to_litter_fire(p)      =  livestemn_storage(p) * f * &
                                                    (1._r8 - cc_other(ivt(p))) * &
                                                    fm_other(ivt(p))   
      m_livestemn_xfer_to_litter_fire(p)         =  livestemn_xfer(p) * f * &
                                                    (1._r8 - cc_other(ivt(p))) * &
                                                    fm_other(ivt(p))
      m_livestemn_to_deadstemn_fire(p)           =  livestemn(p) * f * &
                                                    (1._r8 - cc_lstem(ivt(p))) * &
                                                    (fm_lstem(ivt(p))-fm_droot(ivt(p)))
      m_frootn_to_litter_fire(p)                 =  frootn(p)             * f * &
                                                    fm_root(ivt(p))
      m_frootn_storage_to_litter_fire(p)         =  frootn_storage(p)     * f * &
                                                    fm_other(ivt(p))
      m_frootn_xfer_to_litter_fire(p)            =  frootn_xfer(p)        * f * &
                                                    fm_other(ivt(p))
      m_livecrootn_to_litter_fire(p)             =  livecrootn(p)         * f * &
                                                    fm_droot(ivt(p))
      m_livecrootn_storage_to_litter_fire(p)     =  livecrootn_storage(p) * f * &
                                                    fm_other(ivt(p))
      m_livecrootn_xfer_to_litter_fire(p)        =  livecrootn_xfer(p)    * f * &
                                                    fm_other(ivt(p)) 
      m_livecrootn_to_deadcrootn_fire(p)         =  livecrootn(p)         * f * &
                                                    (fm_lroot(ivt(p))-fm_droot(ivt(p)))
      m_deadcrootn_to_litter_fire(p)             =  deadcrootn(p)         * f * &
                                                    fm_droot(ivt(p))
      m_deadcrootn_storage_to_litter_fire(p)     =  deadcrootn_storage(p) * f * &
                                                    fm_other(ivt(p))
      m_deadcrootn_xfer_to_litter_fire(p)        =  deadcrootn_xfer(p)    * f * &
                                                    fm_other(ivt(p))
      m_retransn_to_litter_fire(p)               =  retransn(p)           * f * &
                                                    (1._r8 - cc_other(ivt(p))) * &
                                                    fm_other(ivt(p)) 
      
#if (defined CNDV)
      if ( woody(ivt(p)) == 1._r8 )then
          if ( livestemc(p)+deadstemc(p) > 0._r8 )then
             nind(p) = nind(p)*(1._r8-1._r8*fm_droot(ivt(p))*f) 
          else
             nind(p) = 0._r8
          end if
      end if
      leafcmax(p) = max(leafc(p)-m_leafc_to_fire(p)*dt, leafcmax(p))
      if (ivt(p) == noveg) leafcmax(p) = 0._r8
#endif

   end do  ! end of pfts loop  
   !
   ! fire-affected carbon to litter and cwd
   !
   do j = 1,nlevdecomp
      do pi = 1,max_pft_per_col
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            if (pi <=  npfts(c)) then
               p = pfti(c) + pi - 1
               if ( pactive(p) ) then
                  
                  fire_mortality_c_to_cwdc(c,j) = fire_mortality_c_to_cwdc(c,j) + &
                        m_deadstemc_to_litter_fire(p) * wtcol(p) * stem_prof(p,j)
                  fire_mortality_c_to_cwdc(c,j) = fire_mortality_c_to_cwdc(c,j) + &
                        m_deadcrootc_to_litter_fire(p) * wtcol(p) * croot_prof(p,j)
                  fire_mortality_n_to_cwdn(c,j) = fire_mortality_n_to_cwdn(c,j) + &
                        m_deadstemn_to_litter_fire(p) * wtcol(p) * stem_prof(p,j)
                  fire_mortality_n_to_cwdn(c,j) = fire_mortality_n_to_cwdn(c,j) + &
                        m_deadcrootn_to_litter_fire(p) * wtcol(p) * croot_prof(p,j)
                  

                  fire_mortality_c_to_cwdc(c,j) = fire_mortality_c_to_cwdc(c,j) + &
                        m_livestemc_to_litter_fire(p) * wtcol(p) * stem_prof(p,j)
                  fire_mortality_c_to_cwdc(c,j) = fire_mortality_c_to_cwdc(c,j) + &
                        m_livecrootc_to_litter_fire(p) * wtcol(p) * croot_prof(p,j)
                  fire_mortality_n_to_cwdn(c,j) = fire_mortality_n_to_cwdn(c,j) + &
                        m_livestemn_to_litter_fire(p) * wtcol(p) * stem_prof(p,j)
                  fire_mortality_n_to_cwdn(c,j) = fire_mortality_n_to_cwdn(c,j) + &
                        m_livecrootn_to_litter_fire(p) * wtcol(p) * croot_prof(p,j)
               
                  
                  m_c_to_litr_met_fire(c,j)=m_c_to_litr_met_fire(c,j) + &
                          ((m_leafc_to_litter_fire(p)*lf_flab(ivt(p)) &
                          +m_leafc_storage_to_litter_fire(p) + &
                           m_leafc_xfer_to_litter_fire(p) + &
                           m_gresp_storage_to_litter_fire(p) &
                          +m_gresp_xfer_to_litter_fire(p))*leaf_prof(p,j) + &
                           (m_frootc_to_litter_fire(p)*fr_flab(ivt(p)) &
                          +m_frootc_storage_to_litter_fire(p) + &
                           m_frootc_xfer_to_litter_fire(p))*froot_prof(p,j) &
                          +(m_livestemc_storage_to_litter_fire(p) + &
                           m_livestemc_xfer_to_litter_fire(p) &
                          +m_deadstemc_storage_to_litter_fire(p) + &
                           m_deadstemc_xfer_to_litter_fire(p))* stem_prof(p,j)&
                          +(m_livecrootc_storage_to_litter_fire(p) + &
                           m_livecrootc_xfer_to_litter_fire(p) &
                          +m_deadcrootc_storage_to_litter_fire(p) + &
                           m_deadcrootc_xfer_to_litter_fire(p))* croot_prof(p,j))* wtcol(p)    
                  m_c_to_litr_cel_fire(c,j)=m_c_to_litr_cel_fire(c,j) + &
                           (m_leafc_to_litter_fire(p)*lf_fcel(ivt(p))*leaf_prof(p,j) + &
                          m_frootc_to_litter_fire(p)*fr_fcel(ivt(p))*froot_prof(p,j))* wtcol(p) 
                  m_c_to_litr_lig_fire(c,j)=m_c_to_litr_lig_fire(c,j) + &
                           (m_leafc_to_litter_fire(p)*lf_flig(ivt(p))*leaf_prof(p,j) + &
                          m_frootc_to_litter_fire(p)*fr_flig(ivt(p))*froot_prof(p,j))* wtcol(p)  
                 
                  m_n_to_litr_met_fire(c,j)=m_n_to_litr_met_fire(c,j) + &
                           ((m_leafn_to_litter_fire(p)*lf_flab(ivt(p)) &
                          +m_leafn_storage_to_litter_fire(p) + &
                           m_leafn_xfer_to_litter_fire(p)+m_retransn_to_litter_fire(p)) &
                          *leaf_prof(p,j) +(m_frootn_to_litter_fire(p)*fr_flab(ivt(p)) &
                          +m_frootn_storage_to_litter_fire(p) + &
                           m_frootn_xfer_to_litter_fire(p))*froot_prof(p,j) &
                          +(m_livestemn_storage_to_litter_fire(p) + &
                           m_livestemn_xfer_to_litter_fire(p) &
                          +m_deadstemn_storage_to_litter_fire(p) + &
                           m_deadstemn_xfer_to_litter_fire(p))* stem_prof(p,j)&
                          +(m_livecrootn_storage_to_litter_fire(p) + &
                           m_livecrootn_xfer_to_litter_fire(p) &
                          +m_deadcrootn_storage_to_litter_fire(p) + &
                           m_deadcrootn_xfer_to_litter_fire(p))* croot_prof(p,j))* wtcol(p)    
                  m_n_to_litr_cel_fire(c,j)=m_n_to_litr_cel_fire(c,j) + &
                           (m_leafn_to_litter_fire(p)*lf_fcel(ivt(p))*leaf_prof(p,j) + &
                          m_frootn_to_litter_fire(p)*fr_fcel(ivt(p))*froot_prof(p,j))* wtcol(p) 
                  m_n_to_litr_lig_fire(c,j)=m_n_to_litr_lig_fire(c,j) + &
                           (m_leafn_to_litter_fire(p)*lf_flig(ivt(p))*leaf_prof(p,j) + &
                          m_frootn_to_litter_fire(p)*fr_flig(ivt(p))*froot_prof(p,j))* wtcol(p) 
               end if
            end if
         end do
      end do
   end do
   !
   ! vertically-resolved decomposing C/N fire loss   
   ! column loop, modified according to clm4_6_00, fzeng, 7 May 2019
   !
   do fc = 1,num_soilc
      c = filter_soilc(fc)

      f = farea_burned(c)

      ! apply this rate to the column state variables to get flux rates
      ! change CC for litter from 0.4_r8 to 0.5_r8 and CC for CWD from 0.2_r8
      ! to 0.25_r8 according to Li et al.(2014), fzeng, 7 May 2019

      do j = 1, nlevdecomp
         ! carbon fluxes
         do l = 1, ndecomp_pools
            if ( is_litter(l) ) then
               m_decomp_cpools_to_fire_vr(c,j,l) = decomp_cpools_vr(c,j,l) * f * 0.5_r8
            end if
            if ( is_cwd(l) ) then
               m_decomp_cpools_to_fire_vr(c,j,l) = decomp_cpools_vr(c,j,l) * &
                                                   (f-baf_crop(c)) * 0.25_r8
            end if
         end do
         
         ! nitrogen fluxes
         do l = 1, ndecomp_pools
            if ( is_litter(l) ) then
               m_decomp_npools_to_fire_vr(c,j,l) = decomp_npools_vr(c,j,l) * f * 0.5_r8
            end if
            if ( is_cwd(l) ) then
               m_decomp_npools_to_fire_vr(c,j,l) = decomp_npools_vr(c,j,l) * &
                                                   (f-baf_crop(c)) * 0.25_r8
            end if
         end do

      end do
   end do  ! end of column loop

   !
   ! carbon loss due to deforestation fires
   ! modified according to clm4_6_00, fzeng, 7 May 2019
   !
   if (fpftdyn /= ' ') then    !true when landuse data is used
      call get_curr_date (kyr, kmo, kda, mcsec)
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         lfc2(c)=0._r8
         if( .not. (kmo == 1 .and. kda == 1 .and. mcsec == 0) )then
            if( trotr1_col(c)+trotr2_col(c) > 0.6_r8 .and. dtrotr_col(c) > 0._r8 .and. &
                lfc(c) > 0._r8 .and. fbac1(c) == 0._r8) then
               lfc2(c) = max(0._r8,min(lfc(c),(farea_burned(c)-baf_crop(c) - &
                         baf_peatf(c))/2.0*dt))/(dtrotr_col(c)*dayspyr*secspday/dt)/dt
               lfc(c)  = lfc(c)-max(0._r8,min(lfc(c),(farea_burned(c)-baf_crop(c) - &
                         baf_peatf(c))*dt/2.0_r8))
            end if
         end if
      end do
   end if
   !
   ! Carbon loss due to peat fires
   !
   ! somc_fire is not connected to clm45 soil carbon pool, ie does not decrease
   ! soil carbon b/c clm4 soil carbon was very low in peatland areas
   ! Fang Li has not checked clm45 soil carbon in peatland areas
   !
   ! modified according to clm4_6_00, fzeng, 7 May 2019
   
   do fc = 1,num_soilc
      c = filter_soilc(fc)
      g = cgridcell(c)
      if( grc%latdeg(g) .lt. borealat)then
         somc_fire(c)= totsomc(c)*baf_peatf(c)*6.0_r8/33.9_r8
      else
         somc_fire(c)= baf_peatf(c)*2.2e3_r8
      end if
   end do

   ! Fang Li has not added aerosol and trace gas emissions due to fire, yet
   ! They will be added here in proportion to the carbon emission
   ! Emission factors differ for various fire types

end subroutine CNFireFluxes

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: hdm_init
!
! !INTERFACE:
!subroutine hdm_init( begg, endg )
!
! !DESCRIPTION:
!
! Initialize data stream information for population density.

! Not used in Catchment-CN. fzeng
!
! !USES:
!  use clm_varctl       , only : inst_name
!  use clm_time_manager , only : get_calendar
!  use ncdio_pio        , only : pio_subsystem
!  use shr_pio_mod      , only : shr_pio_getiotype
!  use clm_nlUtilsMod   , only : find_nlgroup_name
!  use ndepStreamMod    , only : clm_domain_mct
!  use histFileMod      , only : hist_addfld1d
!
! !ARGUMENTS:
!   implicit none
!   integer, intent(IN) :: begg, endg   ! gridcell index bounds
!
! !LOCAL VARIABLES:
!  integer            :: stream_year_first_popdens   ! first year in pop. dens. stream to use
!  integer            :: stream_year_last_popdens    ! last year in pop. dens. stream to use
!  integer            :: model_year_align_popdens    ! align stream_year_first_hdm with 
!  integer            :: nu_nml                      ! unit for namelist file
!  integer            :: nml_error                   ! namelist i/o error flag
!  type(mct_ggrid)    :: dom_clm                     ! domain information 
!  character(len=CL)  :: stream_fldFileName_popdens  ! population density streams filename
!  character(len=CL)  :: popdensmapalgo = 'bilinear' ! mapping alogrithm for population density
!  character(*), parameter :: subName = "('hdmdyn_init')"
!  character(*), parameter :: F00 = "('(hdmdyn_init) ',4a)"
!-----------------------------------------------------------------------
!  namelist /popd_streams/          &
!       stream_year_first_popdens,  &
!       stream_year_last_popdens,   &
!       model_year_align_popdens,   &
!       popdensmapalgo,             &
!       stream_fldFileName_popdens
!EOP
!-----------------------------------------------------------------------

!  ! Allocate pop dens forcing data
!  allocate( forc_hdm(begg:endg) )

!  ! Default values for namelist
!   stream_year_first_popdens  = 1       ! first year in stream to use
!   stream_year_last_popdens   = 1       ! last  year in stream to use
!   model_year_align_popdens   = 1       ! align stream_year_first_popdens with this model year
!   stream_fldFileName_popdens = ' '

!  ! Read popd_streams namelist
!  if (masterproc) then
!     nu_nml = getavu()
!     open( nu_nml, file=trim(NLFilename), status='old', iostat=nml_error )
!     call find_nlgroup_name(nu_nml, 'popd_streams', status=nml_error)
!     if (nml_error == 0) then
!        read(nu_nml, nml=popd_streams,iostat=nml_error)
!        if (nml_error /= 0) then
!            call endrun(subname // ':: ERROR reading popd_streams namelist')
!           stop 'ERROR reading popd_streams namelist'
!        end if
!     end if
!     close(nu_nml)
!     call relavu( nu_nml )
!  endif

!  call shr_mpi_bcast(stream_year_first_popdens, mpicom)
!  call shr_mpi_bcast(stream_year_last_popdens, mpicom)
!  call shr_mpi_bcast(model_year_align_popdens, mpicom)
!  call shr_mpi_bcast(stream_fldFileName_popdens, mpicom)

!  if (masterproc) then
!     write(iulog,*) ' '
!     write(iulog,*) 'popdens_streams settings:'
!     write(iulog,*) '  stream_year_first_popdens  = ',stream_year_first_popdens  
!     write(iulog,*) '  stream_year_last_popdens   = ',stream_year_last_popdens   
!     write(iulog,*) '  model_year_align_popdens   = ',model_year_align_popdens   
!     write(iulog,*) '  stream_fldFileName_popdens = ',stream_fldFileName_popdens
!     write(iulog,*) ' '
!  endif

!  call clm_domain_mct (dom_clm)

!  call shr_strdata_create(sdat_hdm,name="clmhdm",     &
!       pio_subsystem=pio_subsystem,                   & 
!       pio_iotype=shr_pio_getiotype(inst_name),       &
!       mpicom=mpicom, compid=comp_id,                 &
!       gsmap=gsmap_lnd_gdc2glo, ggrid=dom_clm,        &
!       nxg=ldomain%ni, nyg=ldomain%nj,                &
!       yearFirst=stream_year_first_popdens,           &
!       yearLast=stream_year_last_popdens,             &
!       yearAlign=model_year_align_popdens,            &
!       offset=0,                                      &
!       domFilePath='',                                &
!       domFileName=trim(stream_fldFileName_popdens),  &
!       domTvarName='time',                            &
!       domXvarName='lon' ,                            &
!       domYvarName='lat' ,                            &  
!       domAreaName='area',                            &
!       domMaskName='mask',                            &
!       filePath='',                                   &
!       filename=(/trim(stream_fldFileName_popdens)/), &
!       fldListFile='hdm',                             &
!       fldListModel='hdm',                            &
!       fillalgo='none',                               &
!       mapalgo=popdensmapalgo,                        &
!       calendar=get_calendar(),                       &
!       tintalgo='nearest',                            &
!       taxmode='extend'                           )

!  if (masterproc) then
!     call shr_strdata_print(sdat_hdm,'population density data')
!  endif

!  ! Add history fields
!  call hist_addfld1d (fname='HDM', units='counts/km^2',      &
!        avgflag='A', long_name='human population density',   &
!        ptr_lnd=forc_hdm, default='inactive')

!end subroutine hdm_init
  
!================================================================

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: hdm_interp
!
! !INTERFACE:
!subroutine hdm_interp( )
!
! !DESCRIPTION:
!
! Interpolate data stream information for population density.

! Not used in Catchment-CN. fzeng
!
! !USES:
!  use decompMod       , only : get_proc_bounds
!  use clm_time_manager, only : get_curr_date
!
! !ARGUMENTS:
!  implicit none
!
! !LOCAL VARIABLES:
!  integer :: g, ig, begg, endg
!  integer :: year    ! year (0, ...) for nstep+1
!  integer :: mon     ! month (1, ..., 12) for nstep+1
!  integer :: day     ! day of month (1, ..., 31) for nstep+1
!  integer :: sec     ! seconds into current date for nstep+1
!  integer :: mcdate  ! Current model date (yyyymmdd)
!EOP
!-----------------------------------------------------------------------

!  call get_curr_date(year, mon, day, sec)
!  mcdate = year*10000 + mon*100 + day

!  call shr_strdata_advance(sdat_hdm, mcdate, sec, mpicom, 'hdmdyn')

!  call get_proc_bounds(begg, endg)
!  ig = 0
!  do g = begg,endg
!     ig = ig+1
!     forc_hdm(g) = sdat_hdm%avs(1)%rAttr(1,ig)
!  end do
!  
!end subroutine hdm_interp

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: lnfm_init
!
! !INTERFACE:
!subroutine lnfm_init( begg, endg )
!
! !DESCRIPTION:
!
! Initialize data stream information for Lightning.

! Not used in Catchment-CN. fzeng
!
! !USES:
!  use clm_varctl       , only : inst_name
!  use clm_time_manager , only : get_calendar
!  use ncdio_pio        , only : pio_subsystem
!  use shr_pio_mod      , only : shr_pio_getiotype
!  use clm_nlUtilsMod   , only : find_nlgroup_name
!  use ndepStreamMod    , only : clm_domain_mct
!  use histFileMod      , only : hist_addfld1d
!
! !ARGUMENTS:
!  implicit none
!  integer, intent(IN) :: begg, endg   ! gridcell index bounds
!
! !LOCAL VARIABLES:
!  integer            :: stream_year_first_lightng  ! first year in Lightning stream to use
!  integer            :: stream_year_last_lightng   ! last year in Lightning stream to use
!  integer            :: model_year_align_lightng   ! align stream_year_first_lnfm with 
!  integer            :: nu_nml                     ! unit for namelist file
!  integer            :: nml_error                  ! namelist i/o error flag
!  type(mct_ggrid)    :: dom_clm                    ! domain information 
!  character(len=CL)  :: stream_fldFileName_lightng ! lightning stream filename to read
!  character(len=CL)  :: lightngmapalgo = 'bilinear'! Mapping alogrithm
!  character(*), parameter :: subName = "('lnfmdyn_init')"
!  character(*), parameter :: F00 = "('(lnfmdyn_init) ',4a)"
!-----------------------------------------------------------------------
!  namelist /light_streams/         &
!       stream_year_first_lightng,  &
!       stream_year_last_lightng,   &
!       model_year_align_lightng,   &
!       lightngmapalgo,             &
!       stream_fldFileName_lightng
!EOP
!-----------------------------------------------------------------------
!  ! Allocate lightning forcing data
!  allocate( forc_lnfm(begg:endg) )

!  ! Default values for namelist
!   stream_year_first_lightng  = 1      ! first year in stream to use
!   stream_year_last_lightng   = 1      ! last  year in stream to use
!   model_year_align_lightng   = 1      ! align stream_year_first_lnfm with this model year
!   stream_fldFileName_lightng = ' '

!  ! Read light_streams namelist
!  if (masterproc) then
!     nu_nml = getavu()
!     open( nu_nml, file=trim(NLFilename), status='old', iostat=nml_error )
!     call find_nlgroup_name(nu_nml, 'light_streams', status=nml_error)
!     if (nml_error == 0) then
!        read(nu_nml, nml=light_streams,iostat=nml_error)
!        if (nml_error /= 0) then
!            call endrun(subname // ':: ERROR reading light_streams namelist')
!           stop 'ERROR reading light_streams namelist'
!        end if
!     end if
!     close(nu_nml)
!     call relavu( nu_nml )
!  endif

!  call shr_mpi_bcast(stream_year_first_lightng, mpicom)
!  call shr_mpi_bcast(stream_year_last_lightng, mpicom)
!  call shr_mpi_bcast(model_year_align_lightng, mpicom)
!  call shr_mpi_bcast(stream_fldFileName_lightng, mpicom)

!  if (masterproc) then
!     write(iulog,*) ' '
!     write(iulog,*) 'light_stream settings:'
!     write(iulog,*) '  stream_year_first_lightng  = ',stream_year_first_lightng  
!     write(iulog,*) '  stream_year_last_lightng   = ',stream_year_last_lightng   
!     write(iulog,*) '  model_year_align_lightng   = ',model_year_align_lightng   
!     write(iulog,*) '  stream_fldFileName_lightng = ',stream_fldFileName_lightng
!     write(iulog,*) ' '
!  endif
!     write(iulog,*) '  stream_year_last_lightng   = ',stream_year_last_lightng   
!     write(iulog,*) '  stream_year_last_lightng   = ',stream_year_last_lightng   
!     write(iulog,*) '  model_year_align_lightng   = ',model_year_align_lightng   
!     write(iulog,*) '  stream_fldFileName_lightng = ',stream_fldFileName_lightng
!     write(iulog,*) ' '
!  endif

!  call clm_domain_mct (dom_clm)

!  call shr_strdata_create(sdat_lnfm,name="clmlnfm",  &
!       pio_subsystem=pio_subsystem,                  & 
!       pio_iotype=shr_pio_getiotype(inst_name),      &
!       mpicom=mpicom, compid=comp_id,                &
!       gsmap=gsmap_lnd_gdc2glo, ggrid=dom_clm,       &
!       nxg=ldomain%ni, nyg=ldomain%nj,               &
!       yearFirst=stream_year_first_lightng,          &
!       yearLast=stream_year_last_lightng,            &
!       yearAlign=model_year_align_lightng,           &
!       offset=0,                                     &
!       domFilePath='',                               &
!       domFileName=trim(stream_fldFileName_lightng), &
!       domTvarName='time',                           &
!       domXvarName='lon' ,                           &
!       domYvarName='lat' ,                           &  
!       domAreaName='area',                           &
!       domMaskName='mask',                           &
!       filePath='',                                  &
!       filename=(/trim(stream_fldFileName_lightng)/),&
!       fldListFile='lnfm',                           &
!       fldListModel='lnfm',                          &
!       fillalgo='none',                              &
!       mapalgo=lightngmapalgo,                       &
!       calendar=get_calendar(),                      &
!       taxmode='cycle'                            )

!  if (masterproc) then
!     call shr_strdata_print(sdat_lnfm,'Lightning data')
!  endif

!  ! Add history fields
!  call hist_addfld1d (fname='LNFM', units='counts/km^2/hr',  &
!        avgflag='A', long_name='Lightning frequency',        &
!        ptr_lnd=forc_lnfm, default='inactive')

!end subroutine lnfm_init
  
!================================================================

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: lnfm_interp
!
! !INTERFACE:
!subroutine lnfm_interp( )
!
! !DESCRIPTION:
!
! Interpolate data stream information for Lightning.

! Not used in Catchment-CN. fzeng
!
! !USES:
!  use decompMod       , only : get_proc_bounds
!  use clm_time_manager, only : get_curr_date
!
! !ARGUMENTS:
!  implicit none
!
! !LOCAL VARIABLES:
!  integer :: g, ig, begg, endg
!  integer :: year    ! year (0, ...) for nstep+1
!  integer :: mon     ! month (1, ..., 12) for nstep+1
!  integer :: day     ! day of month (1, ..., 31) for nstep+1
!  integer :: sec     ! seconds into current date for nstep+1
!  integer :: mcdate  ! Current model date (yyyymmdd)
!EOP
!-----------------------------------------------------------------------

!  call get_curr_date(year, mon, day, sec)
!  mcdate = year*10000 + mon*100 + day

!  call shr_strdata_advance(sdat_lnfm, mcdate, sec, mpicom, 'lnfmdyn')

!  call get_proc_bounds(begg, endg)
!  ig = 0
!  do g = begg,endg
!     ig = ig+1
!     forc_lnfm(g) = sdat_lnfm%avs(1)%rAttr(1,ig)
!  end do
!  
!end subroutine lnfm_interp

!-----------------------------------------------------------------------
!#endif

end module CNFireMod
