
!   $Id$
module gw_drag_ncar

!---------------------------------------------------------------------------------
! Purpose:
!
! Module to compute the forcing due to parameterized gravity waves. Both an 
! orographic and an internal source spectrum are considered.
!
! Author: Byron Boville
!         In-Sun Song
!
!---------------------------------------------------------------------------------

  use MAPL_Constants, only: MAPL_RGAS, MAPL_GRAV

  use gw_rdg, only     : gw_rdg_ifc
  use gw_oro, only     : gw_oro_ifc
  use gw_convect, only : BeresSourceDesc, gw_beres_ifc
  use gw_common, only  : GWBand,gw_prof

  private                          ! Make default type private to the module
!
! PUBLIC: interfaces
!
  public gw_intr_ncar                   ! interface to actual parameterization

!
! PRIVATE: Rest of the data and interfaces are private to this module
!
  real, parameter :: KWVB    = 6.28e-5        ! effective horizontal wave number for background
  real, parameter :: KWVBEQ  = 6.28e-5/7.     ! effective horizontal wave number for background
  real, parameter :: KWVO    = 6.28e-5        ! effective horizontal wave number for orographic
  real, parameter :: FRACLDV = 0.0            ! fraction of stress deposited in low level region

  real, parameter :: MXASYM  = 0.1            ! max asymmetry between tau(c) and tau(-c)
  real, parameter :: MXRANGE = 0.001          ! max range of tau for all c
  real, parameter :: N2MIN   = 1.e-8          ! min value of bouyancy frequency
  real, parameter :: FCRIT2  = 0.5            ! critical froude number
  real, parameter :: OROHMIN = 10.            ! min surface displacment height for orographic waves
  real, parameter :: OROVMIN = 2.0            ! min wind speed for orographic waves
  real, parameter :: TAUBGND = 6.4            ! background source strength (/TAUSCAL)
  real, parameter :: TAUMIN  = 1.e-10         ! minimum (nonzero) stress
  real, parameter :: TAUSCAL = 0.001          ! scale factor for background stress source
  real, parameter :: UMCFAC  = 0.5            ! factor to limit tendency to prevent reversing u-c
  real, parameter :: UBMC2MN = 0.01           ! min (u-c)**2
  real, parameter :: ZLDVCON = 10.            ! constant for determining zldv from tau0

  real, parameter :: ROG     = MAPL_RGAS/MAPL_GRAV
  real, parameter :: OROKO2  = 0.5 * KWVO     ! 1/2 * horizontal wavenumber
  real, parameter :: PI_GWD  = 4.0*atan(1.0)  ! This is *not* MAPL_PI
contains

!===============================================================================

subroutine gw_intr_ncar(pcols,      pver,         dt,         nrdg,                &    
          beres_dc_desc, beres_band,   oro_band, rdg_band,                           &
          pint_dev,      t_dev,         u_dev,        v_dev,                         &
          ht_dc_dev,     dtdtm_dev,     speed_dev,                                   &
          sgh_dev,       mxdis_dev,     hwdth_dev,    clngt_dev,  angll_dev,         &
          anixy_dev,     gbxar_dev,     kwvrdg_dev,   effrdg_dev, pref_dev,          & 
          pmid_dev,      pdel_dev,      rpdel_dev,    lnpint_dev, zm_dev,  rlat_dev, &
          phis_dev,                                                                  &
          bkg_tau, bkg_tau_cnv, bkg_tau_dry, bkg_tau_mst, &
          dudt_gwd_dev,  dvdt_gwd_dev,  dtdt_gwd_dev,                                &
          dudt_org_dev,  dvdt_org_dev,  dtdt_org_dev,                                &
          tauorox_dev,   tauoroy_dev,   &
          taubkgx_dev,   taubkgy_dev,   &
          taugwx, taugwy, fegw, fepgw, &
          taugwx_east, taugwx_west, taugwy_east, taugwy_west, &
          fegw_east, fegw_west, fepgw_east, fepgw_west, &
          effgworo,      effgwbkg,      alpha, rc            )

  !-----------------------------------------------------------------------
  ! NCAR Gravity Wave Drag Parameterization Interface
  !
  ! This routine orchestrates the complete gravity wave drag calculation
  ! from multiple sources:
  !   1. Deep convection and frontal sources (Beres et al. 2004)
  !   2. Orographic sources (ridge-based or isotropic)
  !
  ! The parameterization computes:
  !   - Wind tendencies (dudt, dvdt) from wave drag
  !   - Temperature tendencies (dtdt) from wave dissipation
  !   - Diagnostic momentum and energy fluxes
  !   - Surface stress from each source
  !
  ! References:
  ! Beres et al. (2004): "A method of specifying the gravity wave spectrum
  ! above convection based on latent heating properties and background wind"
  ! J. Atmos. Sci., Vol 61, No. 3, pp. 324-337.
  !-----------------------------------------------------------------------

  !------------------------------Arguments--------------------------------

  !------ Dimensions ------
  integer, intent(in) :: pcols                    ! Number of atmospheric columns
  integer, intent(in) :: pver                     ! Number of vertical layers
  integer, intent(in) :: nrdg                     ! Number of ridges per grid box
  real, intent(in) :: dt                          ! Time step (s)

  !------ Configuration ------
  type(GWBand), intent(inout) :: beres_band      ! Beres convective band descriptor
  type(GWBand), intent(inout) :: oro_band        ! Isotropic orographic band descriptor
  type(GWBand), intent(inout) :: rdg_band        ! Ridge-based orographic band descriptor
  type(BeresSourceDesc), intent(inout) :: beres_dc_desc ! Beres source configuration
  real, intent(in) :: effgwbkg                   ! Background GWD efficiency factor
  real, intent(in) :: effgworo                   ! Orographic GWD efficiency factor
  real, intent(in) :: alpha(:)                   ! Vertical damping profile

  !------ Pressure Coordinates ------
  real, intent(in) :: pint_dev(pcols,pver+1)     ! Pressure at layer interfaces (Pa)
  real, intent(in) :: pmid_dev(pcols,pver)       ! Pressure at layer midpoints (Pa)
  real, intent(in) :: pdel_dev(pcols,pver)       ! Pressure layer thickness (Pa)
  real, intent(in) :: rpdel_dev(pcols,pver)      ! Inverse pressure thickness (Pa^-1)
  real, intent(in) :: lnpint_dev(pcols,pver+1)   ! Log of interface pressures
  real, intent(in) :: pref_dev(pver+1)           ! Reference pressure at interfaces (Pa)

  !------ Thermodynamic Fields ------
  real, intent(in) :: t_dev(pcols,pver)          ! Temperature (K)
  real, intent(in) :: u_dev(pcols,pver)          ! Zonal wind (m/s)
  real, intent(in) :: v_dev(pcols,pver)          ! Meridional wind (m/s)

  !------ Altitude Coordinates ------
  real, intent(in) :: zm_dev(pcols,pver)         ! Height at layer midpoints (m)
  real, intent(in) :: phis_dev(pcols)            ! Surface geopotential (m^2/s^2)

  !------ Stability Parameters ------
  real, intent(in) :: rlat_dev(pcols)            ! Latitude (radians)

  !------ Convective Heating ------
  real, intent(in) :: ht_dc_dev(pcols,pver)      ! Deep convection heating rate (K/s)
  real, intent(in) :: dtdtm_dev(pcols,pver)      ! Microphysics heating rate (K/s)
  real, intent(in) :: speed_dev(pcols)           ! Max wind in stable surface layer (m/s)

  !------ Orographic Parameters ------
  real, intent(in) :: sgh_dev(pcols)             ! Standard deviation of orography (m)
  real, intent(in) :: mxdis_dev(pcols,nrdg)      ! Ridge/obstacle height (m)
  real, intent(in) :: hwdth_dev(pcols,nrdg)      ! Ridge width (m)
  real, intent(in) :: clngt_dev(pcols,nrdg)      ! Ridge crest length (m)
  real, intent(in) :: angll_dev(pcols,nrdg)      ! Ridge orientation (degrees)
  real, intent(in) :: anixy_dev(pcols,nrdg)      ! Ridge anisotropy parameter
  real, intent(in) :: gbxar_dev(pcols)           ! Grid box area (m^2)
  real, intent(in) :: kwvrdg_dev(pcols,nrdg)     ! Horizontal wavenumber (m^-1)
  real, intent(in) :: effrdg_dev(pcols,nrdg)     ! Ridge efficiency factor

  !------ Output: Background Stress Diagnostics ------
  real, intent(out) :: bkg_tau(pcols)            ! Total background stress (Pa)
  real, intent(out) :: bkg_tau_cnv(pcols)        ! Convective source stress (Pa)
  real, intent(out) :: bkg_tau_dry(pcols)        ! Dry/katabatic source stress (Pa)
  real, intent(out) :: bkg_tau_mst(pcols)        ! Moist source stress (Pa)

  !------ Output: Wind and Temperature Tendencies ------
  real, intent(out) :: dudt_gwd_dev(pcols,pver)  ! Total zonal wind tendency (m/s^2)
  real, intent(out) :: dvdt_gwd_dev(pcols,pver)  ! Total meridional wind tendency (m/s^2)
  real, intent(out) :: dtdt_gwd_dev(pcols,pver)  ! Total temperature tendency (K/s)

  real, intent(out) :: dudt_org_dev(pcols,pver)  ! Orographic zonal wind tendency (m/s^2)
  real, intent(out) :: dvdt_org_dev(pcols,pver)  ! Orographic meridional wind tendency (m/s^2)
  real, intent(out) :: dtdt_org_dev(pcols,pver)  ! Orographic temperature tendency (K/s)

  !------ Output: Surface Stress ------
  real, intent(out) :: tauorox_dev(pcols)        ! Zonal orographic surface stress (Pa)
  real, intent(out) :: tauoroy_dev(pcols)        ! Meridional orographic surface stress (Pa)
  real, intent(out) :: taubkgx_dev(pcols)        ! Zonal background surface stress (Pa)
  real, intent(out) :: taubkgy_dev(pcols)        ! Meridional background surface stress (Pa)

  !------ Output: Momentum and Energy Flux Diagnostics ------
  ! Total fluxes (all phase speeds combined)
  real, intent(out) :: taugwx(pcols,pver)        ! Zonal momentum flux (Pa)
  real, intent(out) :: taugwy(pcols,pver)        ! Meridional momentum flux (Pa)
  real, intent(out) :: fegw(pcols,pver)          ! Kinetic energy flux (W/m^2)
  real, intent(out) :: fepgw(pcols,pver)         ! Phase-speed weighted energy flux (W/m^2)

  ! Eastward-propagating waves (c > 0)
  real, intent(out) :: taugwx_east(pcols,pver)   ! Zonal momentum flux (Pa)
  real, intent(out) :: taugwy_east(pcols,pver)   ! Meridional momentum flux (Pa)
  real, intent(out) :: fegw_east(pcols,pver)     ! Kinetic energy flux (W/m^2)
  real, intent(out) :: fepgw_east(pcols,pver)    ! Phase-speed weighted energy flux (W/m^2)

  ! Westward-propagating waves (c < 0)
  real, intent(out) :: taugwx_west(pcols,pver)   ! Zonal momentum flux (Pa)
  real, intent(out) :: taugwy_west(pcols,pver)   ! Meridional momentum flux (Pa)
  real, intent(out) :: fegw_west(pcols,pver)     ! Kinetic energy flux (W/m^2)
  real, intent(out) :: fepgw_west(pcols,pver)    ! Phase-speed weighted energy flux (W/m^2)

  !------ Return Code ------
  integer, optional, intent(out) :: RC           ! Return code (0 = success)

  !---------------------------Local Storage-------------------------------

  !------ Dimensions ------
  integer :: pverp, pcnst
  integer :: i, k

  !------ Wind and Temperature Tendencies ------
  real :: utgw(pcols,pver)                       ! Zonal wind tendency (m/s^2)
  real :: vtgw(pcols,pver)                       ! Meridional wind tendency (m/s^2)
  real :: ttgw(pcols,pver)                       ! Temperature tendency (K/s)

  !------ Atmospheric Profiles ------
  real :: z(pcols)                               ! Surface elevation (m)
  real :: zi(pcols,pver+1)                       ! Interface heights above ground (m)
  real :: ni(pcols,pver+1)                       ! Interface Brunt-Vaisala frequency (s^-1)
  real :: nm(pcols,pver)                         ! Midpoint Brunt-Vaisala frequency (s^-1)
  real :: rhoi(pcols,pver+1)                     ! Interface density (kg/m^3)
  real :: kvtt(pcols,pver+1)                     ! Molecular thermal diffusivity (m^2/s)

  !------ Energy and Orographic Parameters ------
  real :: flx_heat(pcols)                        ! Energy change (J/m^2)
  real :: rdg_cd_llb                             ! Drag coefficient for low-level flow
  logical :: trpd_leewv                          ! Flag for trapped lee waves

  !-----------------------------------------------------------------------
  ! Initialize dimensions and output arrays
  !-----------------------------------------------------------------------
  pverp = pver + 1
  pcnst = 1

  ! Initialize wind and temperature tendencies
  dudt_gwd_dev(:,:) = 0.0
  dvdt_gwd_dev(:,:) = 0.0
  dtdt_gwd_dev(:,:) = 0.0

  ! Initialize momentum and energy flux diagnostics
  taugwx(:,:) = 0.0
  taugwy(:,:) = 0.0
  fegw(:,:) = 0.0
  fepgw(:,:) = 0.0
  taugwx_east(:,:) = 0.0
  taugwx_west(:,:) = 0.0
  taugwy_east(:,:) = 0.0
  taugwy_west(:,:) = 0.0
  fegw_east(:,:) = 0.0
  fegw_west(:,:) = 0.0
  fepgw_east(:,:) = 0.0
  fepgw_west(:,:) = 0.0

  ! Initialize surface stress diagnostics
  tauorox_dev(:) = 0.0
  tauoroy_dev(:) = 0.0
  taubkgx_dev(:) = 0.0
  taubkgy_dev(:) = 0.0

  ! Initialize working arrays
  kvtt(:,:) = 0.0
  flx_heat(:) = 0.0

  !-----------------------------------------------------------------------
  ! Compute atmospheric profiles (Brunt-Vaisala frequency, density)
  !-----------------------------------------------------------------------
  call gw_prof(pcols, pver, pint_dev, pmid_dev, t_dev, rhoi, nm, ni)

  ! Convert surface geopotential to elevation
  z = phis_dev / MAPL_GRAV

  ! Compute interface heights above ground
  zi(:,pver+1) = 0.0
  do k = 2, pver
     zi(:,k) = 0.5 * (zm_dev(:,k-1) + zm_dev(:,k))
  end do
  zi(:,1) = zi(:,2) + 0.5 * (zm_dev(:,1) - zm_dev(:,2))

  !-----------------------------------------------------------------------
  ! STEP 1: Deep Convection and Frontal Background GWD (Beres Scheme)
  !-----------------------------------------------------------------------
  if (beres_dc_desc%active .and. effgwbkg > 0.0) then

     call gw_beres_ifc(beres_band, &
          pcols, pver, dt, effgwbkg, &
          u_dev, v_dev, t_dev, &
          pref_dev, pint_dev, &
          pdel_dev, rpdel_dev, lnpint_dev, &
          zm_dev, zi, &
          nm, ni, rhoi, kvtt, &
          ht_dc_dev, beres_dc_desc, alpha, &
          bkg_tau, bkg_tau_cnv, bkg_tau_dry, bkg_tau_mst, &
          utgw, vtgw, ttgw, flx_heat, dtdtm_dev, speed_dev, &
          taugwx, taugwy, fegw, fepgw, &
          taugwx_east, taugwx_west, taugwy_east, taugwy_west, &
          fegw_east, fegw_west, fepgw_east, fepgw_west, &
          taubkgx_dev, taubkgy_dev)

     ! Accumulate background GWD tendencies
     dudt_gwd_dev = dudt_gwd_dev + utgw
     dvdt_gwd_dev = dvdt_gwd_dev + vtgw
     dtdt_gwd_dev = dtdt_gwd_dev + ttgw

  endif

  !-----------------------------------------------------------------------
  ! STEP 2: Orographic GWD (Ridge-Based or Isotropic)
  !-----------------------------------------------------------------------
  if (effgworo > 0.0) then

     if (nrdg > 0) then
        !--------------------------------------------------------------------
        ! Ridge-based orographic scheme (anisotropic orography)
        !--------------------------------------------------------------------
        trpd_leewv = .FALSE.
        rdg_cd_llb = 1.0

        call gw_rdg_ifc(rdg_band, &
             pcols, pver, pverp, pcnst, nrdg, dt, &
             u_dev, v_dev, t_dev, &
             pint_dev, pmid_dev, &
             pdel_dev, rpdel_dev, &
             lnpint_dev, zm_dev, zi, z, &
             ni, nm, rhoi, &
             kvtt, &
             kwvrdg_dev, effrdg_dev, &
             hwdth_dev, clngt_dev, gbxar_dev, &
             mxdis_dev, angll_dev, anixy_dev, &
             rdg_cd_llb, trpd_leewv, alpha, &
             utgw, vtgw, ttgw, flx_heat, &
             tauorox_dev, tauoroy_dev)

     else
        !--------------------------------------------------------------------
        ! Isotropic orographic scheme (standard deviation of orography)
        !--------------------------------------------------------------------
        call gw_oro_ifc(oro_band, &
             pcols, pver, dt, effgworo, &
             u_dev, v_dev, t_dev, &
             pint_dev, pmid_dev, &
             pdel_dev, rpdel_dev, lnpint_dev, &
             zm_dev, zi, &
             nm, ni, rhoi, kvtt, &
             sgh_dev, rlat_dev, alpha, &
             utgw, vtgw, ttgw)

     endif

     ! Save orographic tendencies separately
     dudt_org_dev = utgw
     dvdt_org_dev = vtgw
     dtdt_org_dev = ttgw

     ! Accumulate orographic GWD tendencies into total
     dudt_gwd_dev = dudt_gwd_dev + dudt_org_dev
     dvdt_gwd_dev = dvdt_gwd_dev + dvdt_org_dev
     dtdt_gwd_dev = dtdt_gwd_dev + dtdt_org_dev

  endif

  !-----------------------------------------------------------------------
  ! End of gravity wave drag calculation
  !-----------------------------------------------------------------------

end subroutine gw_intr_ncar

end module gw_drag_ncar

