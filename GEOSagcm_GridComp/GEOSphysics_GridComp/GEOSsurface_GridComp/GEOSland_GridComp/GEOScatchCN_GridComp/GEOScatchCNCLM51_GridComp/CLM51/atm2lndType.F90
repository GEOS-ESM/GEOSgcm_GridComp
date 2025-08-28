module atm2lndType

  use shr_kind_mod     , only : r8 => shr_kind_r8
  use clm_varpar       , only : numrad
  use clm_varctl       , only : use_fates, use_luna
  use nanMod           , only : nan
  use decompMod        , only : bounds_type
  
  ! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:

     !
  type, public :: atm2lnd_type

     ! atm->lnd not downscaled                                                                                          
     real(r8), pointer :: forc_u_grc                    (:)   => null() ! atm wind speed, east direction (m/s)          
     real(r8), pointer :: forc_v_grc                    (:)   => null() ! atm wind speed, north direction (m/s)         
     real(r8), pointer :: forc_wind_grc                 (:)   => null() ! atmospheric wind speed                        
     real(r8), pointer :: forc_hgt_grc                  (:)   => null() ! atmospheric reference height (m)              
     real(r8), pointer :: forc_topo_grc                 (:)   => null() ! atmospheric surface height (m)                
     real(r8), pointer :: forc_hgt_u_grc                (:)   => null() ! obs height of wind [m] (new)                  
     real(r8), pointer :: forc_hgt_t_grc                (:)   => null() ! obs height of temperature [m] (new)           
     real(r8), pointer :: forc_hgt_q_grc                (:)   => null() ! obs height of humidity [m] (new)              
     real(r8), pointer :: forc_vp_grc                   (:)   => null() ! atmospheric vapor pressure (Pa)               
     real(r8), pointer :: forc_pco2_grc                 (:)   => null() ! CO2 partial pressure (Pa)                     
     real(r8), pointer :: forc_pco2_240_patch           (:)   => null() ! 10-day mean CO2 partial pressure (Pa)         
     real(r8), pointer :: forc_solad_grc                (:,:) => null() ! direct beam radiation (numrad) (vis=forc_sols , nir=forc_soll )
     real(r8), pointer :: forc_solai_grc                (:,:) => null() ! diffuse radiation (numrad) (vis=forc_solsd, nir=forc_solld)  
     real(r8), pointer :: forc_solar_grc                (:)   => null() ! incident solar radiation                      
     real(r8), pointer :: forc_ndep_grc                 (:)   => null() ! nitrogen deposition rate (gN/m2/s)            
     real(r8), pointer :: forc_pc13o2_grc               (:)   => null() ! C13O2 partial pressure (Pa)                   
     real(r8), pointer :: forc_po2_grc                  (:)   => null() ! O2 partial pressure (Pa)                      
     real(r8), pointer :: forc_po2_240_patch            (:)   => null() ! 10-day mean O2 partial pressure (Pa)          
     real(r8), pointer :: forc_aer_grc                  (:,:) => null() ! aerosol deposition array                      
     real(r8), pointer :: forc_pch4_grc                 (:)   => null() ! CH4 partial pressure (Pa)                     
                                                                                                                        
     real(r8), pointer :: forc_t_not_downscaled_grc     (:)   => null() ! not downscaled atm temperature (Kelvin)       
     real(r8), pointer :: forc_th_not_downscaled_grc    (:)   => null() ! not downscaled atm potential temperature (Kelvin)            
     real(r8), pointer :: forc_pbot_not_downscaled_grc  (:)   => null() ! not downscaled atm pressure (Pa)              
     real(r8), pointer :: forc_pbot240_downscaled_patch (:)   => null() ! 10-day mean downscaled atm pressure (Pa)                     
     real(r8), pointer :: forc_rho_not_downscaled_grc   (:)   => null() ! not downscaled atm density (kg/m**3)                         
     real(r8), pointer :: forc_lwrad_not_downscaled_grc (:)   => null() ! not downscaled atm downwrd IR longwave radiation (W/m**2)

     ! atm->lnd downscaled                                                                                              
     real(r8), pointer :: forc_t_downscaled_col         (:)   => null() ! downscaled atm temperature (Kelvin)           
     real(r8), pointer :: forc_th_downscaled_col        (:)   => null() ! downscaled atm potential temperature (Kelvin) 
     real(r8), pointer :: forc_pbot_downscaled_col      (:)   => null() ! downscaled atm pressure (Pa)                  
     real(r8), pointer :: forc_rho_downscaled_col       (:)   => null() ! downscaled atm density (kg/m**3)              
     real(r8), pointer :: forc_lwrad_downscaled_col     (:)   => null() ! downscaled atm downwrd IR longwave radiation (W/m**2)        
                                                                                                                        
                                                                                                                        
     ! time averaged quantities                                                                                         
     real(r8) , pointer :: fsd24_patch                  (:)   => null() ! patch 24hr average of direct beam radiation   
     real(r8) , pointer :: fsd240_patch                 (:)   => null() ! patch 240hr average of direct beam radiation  
     real(r8) , pointer :: fsi24_patch                  (:)   => null() ! patch 24hr average of diffuse beam radiation  
     real(r8) , pointer :: fsi240_patch                 (:)   => null() ! patch 240hr average of diffuse beam radiation 
     real(r8) , pointer :: wind24_patch                 (:)   => null() ! patch 24-hour running mean of wind            
     real(r8) , pointer :: t_mo_patch                   (:)   => null() ! patch 30-day average temperature (Kelvin)     
     real(r8) , pointer :: t_mo_min_patch               (:)   => null() ! patch annual min of t_mo (Kelvin) 

   contains

     procedure, public :: Init

end type atm2lnd_type                                                                                       
type(atm2lnd_type), public, target, save :: atm2lnd_inst 

contains

!---------------------------------------------------------------------
  subroutine Init(this, bounds)

  ! !DESCRIPTION:
! Initialize CTSM atmosphere2land (forcing type) needed for calling CTSM routines                                 
! jk Apr 2021: type is allocated and initialized to NaN; values are assigned from Catchment states before calls to CLM subroutines are made
! this type is only used to be able to pass Catchment states and fluxes to CLM subroutines in the format they expect         
!                                                                                                                       
! !ARGUMENTS:                                                                                                           
    implicit none                                                                                                       
    type(bounds_type), intent(in) :: bounds
    class(atm2lnd_type)           :: this     

    ! LOCAL:
    real(r8) :: ival  = 0.0_r8  ! initial value
    integer  :: begg, endg
    integer  :: begc, endc
    integer  :: begp, endp
    !------------------------------------------------------------------------

    begg = bounds%begg; endg= bounds%endg
    begc = bounds%begc; endc= bounds%endc
    begp = bounds%begp; endp= bounds%endp


        ! atm->lnd
    allocate(this%forc_u_grc                    (begg:endg))        ; this%forc_u_grc                    (:)   = ival
    allocate(this%forc_v_grc                    (begg:endg))        ; this%forc_v_grc                    (:)   = ival   
    allocate(this%forc_wind_grc                 (begg:endg))        ; this%forc_wind_grc                 (:)   = ival   
    allocate(this%forc_hgt_grc                  (begg:endg))        ; this%forc_hgt_grc                  (:)   = ival   
    allocate(this%forc_topo_grc                 (begg:endg))        ; this%forc_topo_grc                 (:)   = ival   
    allocate(this%forc_hgt_u_grc                (begg:endg))        ; this%forc_hgt_u_grc                (:)   = ival   
    allocate(this%forc_hgt_t_grc                (begg:endg))        ; this%forc_hgt_t_grc                (:)   = ival   
    allocate(this%forc_hgt_q_grc                (begg:endg))        ; this%forc_hgt_q_grc                (:)   = ival   
    allocate(this%forc_vp_grc                   (begg:endg))        ; this%forc_vp_grc                   (:)   = ival   
    allocate(this%forc_pco2_grc                 (begg:endg))        ; this%forc_pco2_grc                 (:)   = ival   
    allocate(this%forc_solad_grc                (begg:endg,numrad)) ; this%forc_solad_grc                (:,:) = ival   
    allocate(this%forc_solai_grc                (begg:endg,numrad)) ; this%forc_solai_grc                (:,:) = ival   
    allocate(this%forc_solar_grc                (begg:endg))        ; this%forc_solar_grc                (:)   = ival   
    allocate(this%forc_ndep_grc                 (begg:endg))        ; this%forc_ndep_grc                 (:)   = ival   
    allocate(this%forc_pc13o2_grc               (begg:endg))        ; this%forc_pc13o2_grc               (:)   = ival   
    allocate(this%forc_po2_grc                  (begg:endg))        ; this%forc_po2_grc                  (:)   = ival   
    allocate(this%forc_aer_grc                  (begg:endg,14))     ; this%forc_aer_grc                  (:,:) = ival   
    allocate(this%forc_pch4_grc                 (begg:endg))        ; this%forc_pch4_grc                 (:)   = ival   
    if(use_luna)then                                                                                                    
     allocate(this%forc_pco2_240_patch          (begp:endp))        ; this%forc_pco2_240_patch           (:)   = ival   
     allocate(this%forc_po2_240_patch           (begp:endp))        ; this%forc_po2_240_patch            (:)   = ival   
     allocate(this%forc_pbot240_downscaled_patch(begp:endp))        ; this%forc_pbot240_downscaled_patch (:)   = ival   
    end if
    ! atm->lnd not downscaled                                                                                           
    allocate(this%forc_t_not_downscaled_grc     (begg:endg))        ; this%forc_t_not_downscaled_grc     (:)   = ival   
    allocate(this%forc_pbot_not_downscaled_grc  (begg:endg))        ; this%forc_pbot_not_downscaled_grc  (:)   = ival   
    allocate(this%forc_th_not_downscaled_grc    (begg:endg))        ; this%forc_th_not_downscaled_grc    (:)   = ival   
    allocate(this%forc_rho_not_downscaled_grc   (begg:endg))        ; this%forc_rho_not_downscaled_grc   (:)   = ival   
    allocate(this%forc_lwrad_not_downscaled_grc (begg:endg))        ; this%forc_lwrad_not_downscaled_grc (:)   = ival   
                                                                                                                        
    ! atm->lnd downscaled                                                                                               
    allocate(this%forc_t_downscaled_col         (begc:endc))        ; this%forc_t_downscaled_col         (:)   = ival   
    allocate(this%forc_pbot_downscaled_col      (begc:endc))        ; this%forc_pbot_downscaled_col      (:)   = ival   
    allocate(this%forc_th_downscaled_col        (begc:endc))        ; this%forc_th_downscaled_col        (:)   = ival   
    allocate(this%forc_rho_downscaled_col       (begc:endc))        ; this%forc_rho_downscaled_col       (:)   = ival   
    allocate(this%forc_lwrad_downscaled_col     (begc:endc))        ; this%forc_lwrad_downscaled_col     (:)   = ival   
                                                                                                                        
    allocate(this%fsd24_patch                   (begp:endp))        ; this%fsd24_patch                   (:)   = nan    
    allocate(this%fsd240_patch                  (begp:endp))        ; this%fsd240_patch                  (:)   = nan    
    allocate(this%fsi24_patch                   (begp:endp))        ; this%fsi24_patch                   (:)   = nan    
    allocate(this%fsi240_patch                  (begp:endp))        ; this%fsi240_patch                  (:)   = nan    
    if (use_fates) then                                                                                                 
       allocate(this%wind24_patch               (begp:endp))        ; this%wind24_patch                  (:)   = nan    
    end if                                                                                                              
    allocate(this%t_mo_patch                    (begp:endp))        ; this%t_mo_patch               (:)   = nan         
    allocate(this%t_mo_min_patch                (begp:endp))        ; this%t_mo_min_patch           (:)   = nan !

  end subroutine Init 

end module atm2lndType
