module CNVegStateType

  use shr_kind_mod     , only : r8 => shr_kind_r8
  use nanMod           , only : nan
  use clm_varpar       , only : nlevsno, nlevgrnd, nlevlak, nlevsoi, &
                                num_zon, num_veg, var_col, var_pft, numpft
  use clm_varcon       , only : spval, ispval
  use decompMod        , only : bounds_type
  use AnnualFluxDribbler, only : annual_flux_dribbler_type, annual_flux_dribbler_patch


  ! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:

  type, public :: cnveg_state_type

     integer  , pointer :: burndate_patch              (:)     ! patch crop burn date
     type(annual_flux_dribbler_type) :: dwt_dribbler_patch     ! object to convert instantaneous dwt values into values that are smoothed / dribbled throughout the year
     real(r8) , pointer :: dwt_smoothed_patch          (:)     ! change in patch weight (-1 to 1) on the gridcell in this time step; changes in first time step of year are smoothed (dribbled) over the whole year

     ! Prognostic crop model
     !
     ! TODO(wjs, 2016-02-22) Most / all of these crop-specific state variables should be
     ! moved to CropType
     real(r8) , pointer :: hdidx_patch                 (:)     ! patch cold hardening index?
     real(r8) , pointer :: cumvd_patch                 (:)     ! patch cumulative vernalization d?ependence?
     real(r8) , pointer :: gddmaturity_patch           (:)     ! patch growing degree days (gdd) needed to harvest (ddays)
     real(r8) , pointer :: huileaf_patch               (:)     ! patch heat unit index needed from planting to leaf emergence
     real(r8) , pointer :: huigrain_patch              (:)     ! patch heat unit index needed to reach vegetative maturity
     real(r8) , pointer :: aleafi_patch                (:)     ! patch saved leaf allocation coefficient from phase 2
     real(r8) , pointer :: astemi_patch                (:)     ! patch saved stem allocation coefficient from phase 2
     real(r8) , pointer :: aleaf_patch                 (:)     ! patch leaf allocation coefficient
     real(r8) , pointer :: astem_patch                 (:)     ! patch stem allocation coefficient
     real(r8) , pointer :: htmx_patch                  (:)     ! patch max hgt attained by a crop during yr (m)
     integer  , pointer :: peaklai_patch               (:)     ! patch 1: max allowed lai; 0: not at max

     integer  , pointer :: idop_patch                  (:)     ! patch date of planting

     real(r8) , pointer :: gdp_lf_col                  (:)     ! col global real gdp data (k US$/capita)
     real(r8) , pointer :: peatf_lf_col                (:)     ! col global peatland fraction data (0-1)
     integer  , pointer :: abm_lf_col                  (:)     ! col global peak month of crop fire emissions 

     real(r8) , pointer :: lgdp_col                    (:)     ! col gdp limitation factor for fire occurrence (0-1)
     real(r8) , pointer :: lgdp1_col                   (:)     ! col gdp limitation factor for fire spreading (0-1)
     real(r8) , pointer :: lpop_col                    (:)     ! col pop limitation factor for fire spreading (0-1)

     real(r8) , pointer :: tempavg_t2m_patch           (:)     ! patch temporary average 2m air temperature (K)
     real(r8) , pointer :: annavg_t2m_patch            (:)     ! patch annual average 2m air temperature (K)
     real(r8) , pointer :: annavg_t2m_col              (:)     ! col annual average of 2m air temperature, averaged from patch-level (K)
     real(r8) , pointer :: annsum_counter_col          (:)     ! col seconds since last annual accumulator turnover

     ! Fire
     real(r8) , pointer :: nfire_col                   (:)     ! col fire counts (count/km2/sec), valid only in Reg. C
     real(r8) , pointer :: fsr_col                     (:)     ! col fire spread rate at column level (m/s)
     real(r8) , pointer :: fd_col                      (:)     ! col fire duration at column level (hr)
     real(r8) , pointer :: lfc_col                     (:)     ! col conversion area fraction of BET and BDT that haven't burned before (/timestep)
     real(r8) , pointer :: lfc2_col                    (:)     ! col conversion area fraction of BET and BDT that burned (/sec)
     real(r8) , pointer :: dtrotr_col                  (:)     ! col annual decreased fraction coverage of BET on the gridcell (0-1)
     real(r8) , pointer :: trotr1_col                  (:)     ! col patch weight of BET on the column (0-1)
     real(r8) , pointer :: trotr2_col                  (:)     ! col patch weight of BDT on the column (0-1)
     real(r8) , pointer :: cropf_col                   (:)     ! col crop fraction in veg column (0-1)
     real(r8) , pointer :: baf_crop_col                (:)     ! col baf for cropland(/sec)
     real(r8) , pointer :: baf_peatf_col               (:)     ! col baf for peatland (/sec)
     real(r8) , pointer :: fbac_col                    (:)     ! col total burned area out of conversion (/sec)
     real(r8) , pointer :: fbac1_col                   (:)     ! col burned area out of conversion region due to land use fire (/sec)
     real(r8) , pointer :: wtlf_col                    (:)     ! col fractional coverage of non-crop Patches (0-1)
     real(r8) , pointer :: lfwt_col                    (:)     ! col fractional coverage of non-crop and non-bare-soil Patches (0-1)
     real(r8) , pointer :: farea_burned_col            (:)     ! col fractional area burned (/sec) 

     real(r8), pointer :: dormant_flag_patch           (:)     ! patch dormancy flag
     real(r8), pointer :: days_active_patch            (:)     ! patch number of days since last dormancy
     real(r8), pointer :: onset_flag_patch             (:)     ! patch onset flag
     real(r8), pointer :: onset_counter_patch          (:)     ! patch onset days counter
     real(r8), pointer :: onset_gddflag_patch          (:)     ! patch onset flag for growing degree day sum
     real(r8), pointer :: onset_fdd_patch              (:)     ! patch onset freezing degree days counter
     real(r8), pointer :: onset_gdd_patch              (:)     ! patch onset growing degree days
     real(r8), pointer :: onset_swi_patch              (:)     ! patch onset soil water index
     real(r8), pointer :: offset_flag_patch            (:)     ! patch offset flag
     real(r8), pointer :: offset_counter_patch         (:)     ! patch offset days counter
     real(r8), pointer :: offset_fdd_patch             (:)     ! patch offset freezing degree days counter
     real(r8), pointer :: offset_swi_patch             (:)     ! patch offset soil water index
     real(r8), pointer :: grain_flag_patch             (:)     ! patch 1: grain fill stage; 0: not
     real(r8), pointer :: lgsf_patch                   (:)     ! patch long growing season factor [0-1]
     real(r8), pointer :: bglfr_patch                  (:)     ! patch background litterfall rate (1/s)
     real(r8), pointer :: bgtr_patch                   (:)     ! patch background transfer growth rate (1/s)
     real(r8), pointer :: c_allometry_patch            (:)     ! patch C allocation index (DIM)
     real(r8), pointer :: n_allometry_patch            (:)     ! patch N allocation index (DIM)

     real(r8), pointer :: tempsum_potential_gpp_patch  (:)     ! patch temporary annual sum of potential GPP
     real(r8), pointer :: annsum_potential_gpp_patch   (:)     ! patch annual sum of potential GPP
     real(r8), pointer :: tempmax_retransn_patch       (:)     ! patch temporary annual max of retranslocated N pool (gN/m2)
     real(r8), pointer :: annmax_retransn_patch        (:)     ! patch annual max of retranslocated N pool (gN/m2)
     real(r8), pointer :: downreg_patch                (:)     ! patch fractional reduction in GPP due to N limitation (DIM)
     real(r8), pointer :: leafcn_offset_patch          (:)     ! patch leaf C:N used by FUN
     real(r8), pointer :: plantCN_patch                (:)     ! patch plant C:N used by FUN

    contains

     procedure, public :: Init

   end type cnveg_state_type

   type(cnveg_state_type), public, target, save :: cnveg_state_inst

contains

!-----------------------------------------------------
!----------------------------------------------
  subroutine Init(this, bounds, nch, ityp, fveg, cncol, cnpft)

! !DESCRIPTION:
! Initialize CTSM vegetation states
! jk Apr 2021: type is allocated and initialized to NaN;
! if data arrays from restart file are passed (cncol and cnpft), the type is then initialized with these values
!
! !ARGUMENTS:
    implicit none

  ! INPUT
    type(bounds_type),                            intent(in) :: bounds
    integer,                                      intent(in) :: nch ! number of tiles
    integer, dimension(nch,NUM_VEG,NUM_ZON),      intent(in) :: ityp ! PFT index
    real, dimension(nch,NUM_VEG,NUM_ZON),         intent(in) :: fveg    ! PFT fraction
    real, dimension(nch,NUM_ZON,VAR_COL),         intent(in) :: cncol ! gkw: column CN restart
    real, dimension(nch,NUM_ZON,NUM_VEG,VAR_PFT), intent(in) :: cnpft ! gkw: PFT CN restart
    class(cnveg_state_type)                                  :: this

    ! LOCAL
    integer  :: begp, endp
    integer  :: begc, endc
    integer  :: np, nc, nz, p, nv, n
    !--------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc

    allocate(this%burndate_patch      (begp:endp))                   ; this%burndate_patch      (:)   = ispval
    allocate(this%dwt_smoothed_patch  (begp:endp))                   ; this%dwt_smoothed_patch  (:)   = nan

    allocate(this%hdidx_patch         (begp:endp))                   ; this%hdidx_patch         (:)   = nan
    allocate(this%cumvd_patch         (begp:endp))                   ; this%cumvd_patch         (:)   = nan
    allocate(this%gddmaturity_patch   (begp:endp))                   ; this%gddmaturity_patch   (:)   = spval
    allocate(this%huileaf_patch       (begp:endp))                   ; this%huileaf_patch       (:)   = nan
    allocate(this%huigrain_patch      (begp:endp))                   ; this%huigrain_patch      (:)   = 0.0_r8
    allocate(this%aleafi_patch        (begp:endp))                   ; this%aleafi_patch        (:)   = nan
    allocate(this%astemi_patch        (begp:endp))                   ; this%astemi_patch        (:)   = nan
    allocate(this%aleaf_patch         (begp:endp))                   ; this%aleaf_patch         (:)   = nan
    allocate(this%astem_patch         (begp:endp))                   ; this%astem_patch         (:)   = nan
    allocate(this%htmx_patch          (begp:endp))                   ; this%htmx_patch          (:)   = 0.0_r8
    allocate(this%peaklai_patch       (begp:endp))                   ; this%peaklai_patch       (:)   = 0

    allocate(this%idop_patch          (begp:endp))                   ; this%idop_patch          (:)   = huge(1)

    allocate(this%gdp_lf_col          (begc:endc))                   ;
    allocate(this%peatf_lf_col        (begc:endc))                   ;
    allocate(this%abm_lf_col          (begc:endc))                   ;

    allocate(this%lgdp_col            (begc:endc))                   ;
    allocate(this%lgdp1_col           (begc:endc))                   ;
    allocate(this%lpop_col            (begc:endc))                   ;

    allocate(this%tempavg_t2m_patch   (begp:endp))                   ; this%tempavg_t2m_patch   (:)   = nan
    allocate(this%annsum_counter_col  (begc:endc))                   ; this%annsum_counter_col  (:)   = nan
    allocate(this%annavg_t2m_col      (begc:endc))                   ; this%annavg_t2m_col      (:)   = nan
    allocate(this%annavg_t2m_patch    (begp:endp))                   ; this%annavg_t2m_patch    (:)   = spval

    allocate(this%nfire_col           (begc:endc))                   ; this%nfire_col           (:)   = spval
    allocate(this%fsr_col             (begc:endc))                   ; this%fsr_col             (:)   = nan
    allocate(this%fd_col              (begc:endc))                   ; this%fd_col              (:)   = nan
    allocate(this%lfc_col             (begc:endc))                   ; this%lfc_col             (:)   = spval
    allocate(this%lfc2_col            (begc:endc))                   ; this%lfc2_col            (:)   = 0._r8
    allocate(this%dtrotr_col          (begc:endc))                   ; this%dtrotr_col          (:)   = 0._r8
    allocate(this%trotr1_col          (begc:endc))                   ; this%trotr1_col          (:)   = 0._r8
    allocate(this%trotr2_col          (begc:endc))                   ; this%trotr2_col          (:)   = 0._r8
    allocate(this%cropf_col           (begc:endc))                   ; this%cropf_col           (:)   = nan
    allocate(this%baf_crop_col        (begc:endc))                   ; this%baf_crop_col        (:)   = nan
    allocate(this%baf_peatf_col       (begc:endc))                   ; this%baf_peatf_col       (:)   = nan
    allocate(this%fbac_col            (begc:endc))                   ; this%fbac_col            (:)   = nan
    allocate(this%fbac1_col           (begc:endc))                   ; this%fbac1_col           (:)   = nan
    allocate(this%wtlf_col            (begc:endc))                   ; this%wtlf_col            (:)   = nan
    allocate(this%lfwt_col            (begc:endc))                   ; this%lfwt_col            (:)   = nan
    allocate(this%farea_burned_col    (begc:endc))                   ; this%farea_burned_col    (:)   = nan

    allocate(this%dormant_flag_patch          (begp:endp)) ;    this%dormant_flag_patch          (:) = nan
    allocate(this%days_active_patch           (begp:endp)) ;    this%days_active_patch           (:) = nan
    allocate(this%onset_flag_patch            (begp:endp)) ;    this%onset_flag_patch            (:) = nan
    allocate(this%onset_counter_patch         (begp:endp)) ;    this%onset_counter_patch         (:) = nan
    allocate(this%onset_gddflag_patch         (begp:endp)) ;    this%onset_gddflag_patch         (:) = nan
    allocate(this%onset_fdd_patch             (begp:endp)) ;    this%onset_fdd_patch             (:) = nan
    allocate(this%onset_gdd_patch             (begp:endp)) ;    this%onset_gdd_patch             (:) = nan
    allocate(this%onset_swi_patch             (begp:endp)) ;    this%onset_swi_patch             (:) = nan
    allocate(this%offset_flag_patch           (begp:endp)) ;    this%offset_flag_patch           (:) = nan
    allocate(this%offset_counter_patch        (begp:endp)) ;    this%offset_counter_patch        (:) = nan
    allocate(this%offset_fdd_patch            (begp:endp)) ;    this%offset_fdd_patch            (:) = nan
    allocate(this%offset_swi_patch            (begp:endp)) ;    this%offset_swi_patch            (:) = nan
    allocate(this%grain_flag_patch            (begp:endp)) ;    this%grain_flag_patch            (:) = nan
    allocate(this%lgsf_patch                  (begp:endp)) ;    this%lgsf_patch                  (:) = nan
    allocate(this%bglfr_patch                 (begp:endp)) ;    this%bglfr_patch                 (:) = nan
    allocate(this%bgtr_patch                  (begp:endp)) ;    this%bgtr_patch                  (:) = nan
    allocate(this%c_allometry_patch           (begp:endp)) ;    this%c_allometry_patch           (:) = nan
    allocate(this%n_allometry_patch           (begp:endp)) ;    this%n_allometry_patch           (:) = nan
    allocate(this%tempsum_potential_gpp_patch (begp:endp)) ;    this%tempsum_potential_gpp_patch (:) = nan
    allocate(this%annsum_potential_gpp_patch  (begp:endp)) ;    this%annsum_potential_gpp_patch  (:) = nan
    allocate(this%tempmax_retransn_patch      (begp:endp)) ;    this%tempmax_retransn_patch      (:) = nan
    allocate(this%annmax_retransn_patch       (begp:endp)) ;    this%annmax_retransn_patch       (:) = nan
    allocate(this%downreg_patch               (begp:endp)) ;    this%downreg_patch               (:) = nan
    allocate(this%leafcn_offset_patch         (begp:endp)) ;    this%leafcn_offset_patch         (:) = nan
    allocate(this%plantCN_patch               (begp:endp)) ;    this%plantCN_patch               (:) = nan

 ! initialize variables from restart file or set to cold start value
 n = 0
 np = 0
    do nc = 1,nch           ! catchment tile loop
       do nz = 1,num_zon    ! CN zone loop
          n = n + 1
          this%annsum_counter_col (n) = cncol(nc,nz, 31)
          this%annavg_t2m_col     (n) = cncol(nc,nz, 32) 
          this%farea_burned_col   (n) = cncol(nc,nz, 34)

          do p = 0,numpft  ! PFT index loop
             np = np + 1
             do nv = 1,num_veg ! defined veg loop
                if(ityp(nc,nv,nz)==p .and. fveg(nc,nv,nz)>1.e-4) then
                  
                  this%annavg_t2m_patch            (np) = cnpft(nc,nz,nv,  24)
                  this%annmax_retransn_patch       (np) = cnpft(nc,nz,nv, 25)
                  this%annsum_potential_gpp_patch  (np) = cnpft(nc,nz,nv, 27)
                  this%days_active_patch           (np) = cnpft(nc,nz,nv, 29)
                  this%dormant_flag_patch          (np) = cnpft(nc,nz,nv, 30)
                  this%offset_counter_patch        (np) = cnpft(nc,nz,nv, 31)
                  this%offset_fdd_patch            (np) = cnpft(nc,nz,nv, 32)
                  this%offset_flag_patch           (np) = cnpft(nc,nz,nv, 33)
                  this%offset_swi_patch            (np) = cnpft(nc,nz,nv, 34)
                  this%onset_counter_patch         (np) = cnpft(nc,nz,nv, 35)
                  this%onset_fdd_patch             (np) = cnpft(nc,nz,nv, 36)
                  this%onset_flag_patch            (np) = cnpft(nc,nz,nv, 37)
                  this%onset_gdd_patch             (np) = cnpft(nc,nz,nv, 38)
                  this%onset_gddflag_patch         (np) = cnpft(nc,nz,nv, 39)
                  this%onset_swi_patch             (np) = cnpft(nc,nz,nv, 40)
                  this%tempavg_t2m_patch           (np) = cnpft(nc,nz,nv, 43)
                  this%tempmax_retransn_patch      (np) = cnpft(nc,nz,nv, 44)
                  this%tempsum_potential_gpp_patch (np) = cnpft(nc,nz,nv, 46)

                 end if
            end do !nv
       end do ! p
     end do ! nz
  end do ! nc   

  end subroutine Init


end module CNVegStateType
