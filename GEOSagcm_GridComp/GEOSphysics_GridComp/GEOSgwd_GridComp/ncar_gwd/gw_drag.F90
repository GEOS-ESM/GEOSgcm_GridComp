
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

  use MAPL_ConstantsMod, only: MAPL_P00,  MAPL_CP, MAPL_GRAV, &
                               MAPL_RGAS, MAPL_VIREPS

  use gw_oro, only     : gw_oro_ifc
  use gw_convect, only : BeresSourceDesc, gw_beres_ifc
  use gw_common, only  : GWBand,gw_prof

  implicit none

  !save
  private                          ! Make default type private to the module
!
! PUBLIC: interfaces
!
  public gw_intr_ncar                   ! interface to actual parameterization

!
! PRIVATE: Rest of the data and interfaces are private to this module
!
integer,parameter :: r8 = selected_real_kind(12)

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
  real, parameter :: TNDMAX  = 500. / 86400.  ! maximum wind tendency
  real, parameter :: UMCFAC  = 0.5            ! factor to limit tendency to prevent reversing u-c
  real, parameter :: UBMC2MN = 0.01           ! min (u-c)**2
  real, parameter :: ZLDVCON = 10.            ! constant for determining zldv from tau0

  real, parameter :: ROG     = MAPL_RGAS/MAPL_GRAV
  real, parameter :: OROKO2  = 0.5 * KWVO     ! 1/2 * horizontal wavenumber
  real, parameter :: PI_GWD  = 4.0*atan(1.0)  ! This is *not* MAPL_PI
contains

!===============================================================================


  subroutine gw_intr_ncar(pcols,      pver,         dt,         pgwv,              &
          beres_desc,   beres_band,   oro_band,                                    &
          pint_dev,     t_dev,        u_dev,        v_dev,      ht_dpc_dev,        &
          sgh_dev,      pref_dev,                                                  & 
          pmid_dev,     pdel_dev,     rpdel_dev,    lnpint_dev, zm_dev,  rlat_dev, &
          dudt_gwd_dev, dvdt_gwd_dev, dtdt_gwd_dev,                                &
          dudt_org_dev, dvdt_org_dev, dtdt_org_dev,                                &
          taugwdx_dev,  taugwdy_dev,  tauox_dev,    tauoy_dev,  feo_dev,           &
          taubkgx_dev,  taubkgy_dev,  taubx_dev,    tauby_dev,  feb_dev,           &
          fepo_dev,     fepb_dev,     utbsrc_dev,   vtbsrc_dev, ttbsrc_dev,        &
          bgstressmax,  effgworo,     effgwbkg,     rc            )

!-----------------------------------------------------------------------
! Interface for multiple gravity wave drag parameterization.
!-----------------------------------------------------------------------

!------------------------------Arguments--------------------------------
    integer, intent(in   ) :: pcols                    ! number of columns
    integer, intent(in   ) :: pver                     ! number of vertical layers
    real,    intent(in   ) :: dt                       ! time step
    integer, intent(in   ) :: pgwv                     ! number of waves allowed                (Default = 4, 0 nullifies)
    type(GWBand),          intent(inout) :: oro_band   ! Band descriptor
    type(GWBand),          intent(inout) :: beres_band ! Band descriptor
    type(BeresSourceDesc), intent(inout) :: beres_desc ! Table descriptor for Beres scheme
    real,    intent(in   ) :: bgstressmax              ! Max of equatorial profile of BG stress factor
    real,    intent(in   ) :: effgwbkg                 ! tendency efficiency for background gwd (Default = 0.125)
    real,    intent(in   ) :: effgworo                 ! tendency efficiency for orographic gwd (Default = 0.125)
    real,    intent(in   ) :: pint_dev(pcols,pver+1)   ! pressure at the layer edges
    real,    intent(in   ) :: t_dev(pcols,pver)        ! temperature at layers
    real,    intent(in   ) :: u_dev(pcols,pver)        ! zonal wind at layers
    real,    intent(in   ) :: v_dev(pcols,pver)        ! meridional wind at layers
    real,    intent(in   ) :: ht_dpc_dev(pcols,pver)   ! moist heating in layers
    real,    intent(in   ) :: sgh_dev(pcols)           ! standard deviation of orography
    real,    intent(in   ) :: pref_dev(pver+1)         ! reference pressure at the layeredges
    real,    intent(in   ) :: pmid_dev(pcols,pver)     ! pressure at the layers
    real,    intent(in   ) :: pdel_dev(pcols,pver)     ! pressure thickness at the layers
    real,    intent(in   ) :: rpdel_dev(pcols,pver)    ! 1.0 / pdel
    real,    intent(in   ) :: lnpint_dev(pcols,pver+1) ! log(pint)
    real,    intent(in   ) :: zm_dev(pcols,pver)       ! height above surface at layers
    real,    intent(in   ) :: rlat_dev(pcols)          ! latitude in radian
    
    real,    intent(  out) :: dudt_gwd_dev(pcols,pver) ! zonal wind tendency at layer 
    real,    intent(  out) :: dvdt_gwd_dev(pcols,pver) ! meridional wind tendency at layer 
    real,    intent(  out) :: dtdt_gwd_dev(pcols,pver) ! temperature tendency at layer
    real,    intent(  out) :: dudt_org_dev(pcols,pver) ! zonal wind tendency at layer due to orography GWD
    real,    intent(  out) :: dvdt_org_dev(pcols,pver) ! meridional wind tendency at layer  due to orography GWD
    real,    intent(  out) :: dtdt_org_dev(pcols,pver) ! temperature tendency at layer  due to orography GWD
    real,    intent(  out) :: taugwdx_dev(pcols)       ! zonal      gravity wave surface    stress
    real,    intent(  out) :: taugwdy_dev(pcols)       ! meridional gravity wave surface    stress
    real,    intent(  out) :: tauox_dev(pcols,pver+1)  ! zonal      orographic gravity wave stress
    real,    intent(  out) :: tauoy_dev(pcols,pver+1)  ! meridional orographic gravity wave stress
    real,    intent(  out) :: feo_dev  (pcols,pver+1)  ! energy flux of orographic gravity waves
    real,    intent(  out) :: fepo_dev (pcols,pver+1)  ! pseudoenergy flux of orographic gravity waves
    real,    intent(  out) :: taubkgx_dev(pcols)       ! zonal      gravity wave background stress
    real,    intent(  out) :: taubkgy_dev(pcols)       ! meridional gravity wave background stress
    real,    intent(  out) :: taubx_dev(pcols,pver+1)  ! zonal      background gravity wave stress
    real,    intent(  out) :: tauby_dev(pcols,pver+1)  ! meridional background gravity wave stress
    real,    intent(  out) :: feb_dev  (pcols,pver+1)  ! energy flux of background gravity waves
    real,    intent(  out) :: fepb_dev (pcols,pver+1)  ! pseudoenergy flux of background gravity waves
    real,    intent(  out) :: utbsrc_dev(pcols,pver)   ! dU/dt below background launch level
    real,    intent(  out) :: vtbsrc_dev(pcols,pver)   ! dV/dt below background launch level
    real,    intent(  out) :: ttbsrc_dev(pcols,pver)   ! dT/dt below background launch level

    integer, optional, intent(out) :: RC               ! return code


#ifndef GPU_MAXLEVS
#define GPU_MAXLEVS pver
#endif

#ifndef MAXPGWV
#define MAXPGWV pgwv
#endif

#define CAMGWCODE

#ifdef CAMGWCODE
#define IFCDIMS 1:pver+1
#else
#define IFCDIMS 0:pver
#endif

!---------------------------Local storage-------------------------------

    integer :: i,k,kc                   ! loop indexes
    integer :: kbotoro                  ! launch-level index for orographic
    integer :: kbotbg                   ! launch-level index for background
    integer :: ktopbg, ktoporo          ! top interface of gwd region
    integer :: kldv                     ! top interface of low level stress divergence region
    integer :: kldvmn                   ! min value of kldv
    integer :: ksrc                     ! index of top interface of source region
    integer :: ksrcmn                   ! min value of ksrc

    integer :: src_level(pcols)
    integer :: tend_level(pcols)


    real(r8)    :: ttgw(pcols,pver)            ! temperature tendency
    real(r8)    :: utgw(pcols,pver)            ! zonal wind tendency
    real(r8)    :: vtgw(pcols,pver)            ! meridional wind tendency

    real(r8)    :: zi(pcols,pver+1)                   ! interface heights above ground
    real(r8)    :: ni(pcols,pver+1)                   ! interface Brunt-Vaisalla frequency
    real(r8)    :: nm(pcols,pver)                     ! midpoint Brunt-Vaisalla frequency
    !!!real    :: rdpldv                              ! 1/dp across low level divergence region
    real(r8)    :: rhoi(pcols,pver+1)                 ! interface density
    real(r8)    :: tau0x(pcols)                       ! c=0 sfc. stress (zonal)
    real(r8)    :: tau0y(pcols)                       ! c=0 sfc. stress (meridional)
    real(r8)    :: ti(pcols,pver+1)                   ! interface temperature
    real(r8)    :: ubi(pcols,pver+1)                  ! projection of wind at interfaces
    real(r8)    :: ubm(pcols,pver)                    ! projection of wind at midpoints
    real(r8)    :: xv(pcols)                          ! unit vectors of source wind (x)
    real(r8)    :: yv(pcols)                          ! unit vectors of source wind (y)
    real(r8)    :: kvtt(pcols,pver+1) ! Molecular thermal diffusivity.

    real(r8)    :: utosrc(pcols,pver)
    real(r8)    :: vtosrc(pcols,pver)
    real(r8)    :: ttosrc(pcols,pver)

    real(r8)    :: maxq0(pcols),hdepth(pcols)

    real(r8)    :: flx_heat(pcols)


    real(r8), allocatable    :: c  (:,:)    ! wave phase speeds
    real(r8), allocatable    :: tau(:,:,:)  ! wave Reynolds stress
    real(r8), allocatable    :: gwut(:,:,:) ! wind speed tendency from each wave


    !!real(r8) ::  pint_dev_r8(pcols,pver+1) , pmid_dev_r8(pcols,pver) , t_dev_r8(pcols,pver) 

    real(r8)  :: bgstressmax_ff              ! Max of equatorial profile of BG stress factor
    real(r8)  :: effgwbkg_ff                 ! tendency efficiency for background gwd (Default = 0.125)
    real(r8)  :: effgworo_ff                 ! tendency efficiency for orographic gwd (Default = 0.125)
    real(r8)  :: pint_dev_ff(pcols,pver+1)   ! pressure at the layer edges
    real(r8)  :: t_dev_ff(pcols,pver)        ! temperature at layers
    real(r8)  :: u_dev_ff(pcols,pver)        ! zonal wind at layers
    real(r8)  :: v_dev_ff(pcols,pver)        ! meridional wind at layers
    real(r8)  :: ht_dpc_dev_ff(pcols,pver)   ! moist heating in layers
    real(r8)  :: sgh_dev_ff(pcols)           ! standard deviation of orography
    real(r8)  :: pref_dev_ff(pver+1)         ! reference pressure at the layeredges
    real(r8)  :: pmid_dev_ff(pcols,pver)     ! pressure at the layers
    real(r8)  :: pdel_dev_ff(pcols,pver)     ! pressure thickness at the layers
    real(r8)  :: rpdel_dev_ff(pcols,pver)    ! 1.0 / pdel
    real(r8)  :: lnpint_dev_ff(pcols,pver+1) ! log(pint)
    real(r8)  :: zm_dev_ff(pcols,pver)       ! height above surface at layers
    real(r8)  :: rlat_dev_ff(pcols)          ! latitude in radian

    real(r8)  :: t_gwt_dc_ff(pcols,pver)     ! temperature tendency at layers from deep conv GW
    real(r8)  :: u_gwt_dc_ff(pcols,pver)     ! zonal wind tendency at layers  "
    real(r8)  :: v_gwt_dc_ff(pcols,pver)     ! meridional tendency wind at layers  "

    real(r8)  :: t_gwt_org_ff(pcols,pver)     ! temperature tendency at layers from orographic GW
    real(r8)  :: u_gwt_org_ff(pcols,pver)     ! zonal wind tendency at layers   "
    real(r8)  :: v_gwt_org_ff(pcols,pver)     ! meridional tendency wind at layers  "

    real(r8)  :: t_gwt_ff(pcols,pver)        ! temperature tendency at layers
    real(r8)  :: u_gwt_ff(pcols,pver)        ! zonal wind tendency at layers
    real(r8)  :: v_gwt_ff(pcols,pver)        ! meridional tendency wind at layers


    real(r8)  :: effgw_dp, dt_ff
!-----------------------------------------------------------------------------

! Initialize accumulated tendencies
! and other things ...
  t_gwt_ff(:,:) = 0._r8
  u_gwt_ff(:,:) = 0._r8
  v_gwt_ff(:,:) = 0._r8

  kvtt(:,:)  = 0._r8

! Calling CAM6 GW codes 

 dt_ff         =  dt
 rlat_dev_ff   =  rlat_dev

 pref_dev_ff   =  pref_dev 
 pint_dev_ff   =  pint_dev 
 pmid_dev_ff   =  pmid_dev
 pdel_dev_ff   =  pdel_dev
 rpdel_dev_ff  =  rpdel_dev
 lnpint_dev_ff =  lnpint_dev

! Met profiles u,v,t,zm
!-----------------------
 t_dev_ff      =  t_dev
 u_dev_ff      =  u_dev
 v_dev_ff      =  v_dev
 zm_dev_ff     =  zm_dev

! Heating 
!----------
ht_dpc_dev_ff  =  ht_dpc_dev

! SGH
!----------
sgh_dev_ff  =  sgh_dev


call gw_prof (pcols , pver, pint_dev_ff , pmid_dev_ff , t_dev_ff , rhoi, nm, ni )


   ! Allocate wavenumber fields.
   allocate(tau(pcols,-beres_band%ngwv:beres_band%ngwv,pver+1))
   !!allocate(gwut(ncol,pver,-band%ngwv:band%ngwv))
   allocate(c(pcols,-beres_band%ngwv:beres_band%ngwv))


    zi(:,pver+1) = 0.0
    do k=2,pver 
       zi(:,k)  =  0.5 * ( zm_dev_ff(:,k-1) + zm_dev_ff(:,k) )  
    end do
    zi(:,1) = zi(:,2) + 0.5*( zm_dev_ff(:,1) - zm_dev_ff(:,2)  )



!get rid of lchnk
    effgw_dp = 0.5_r8
    call gw_beres_ifc( beres_band, &
       pcols, pver, dt_ff , effgw_dp,  &
       u_dev_ff , v_dev_ff, t_dev_ff, &
       pref_dev_ff, pint_dev_ff, & 
       pdel_dev_ff , rpdel_dev_ff, lnpint_dev_ff, &
       zm_dev_ff, zi, &
       nm, ni, rhoi, kvtt,  &
       ht_dpc_dev_ff,beres_desc,rlat_dev_ff, &
       u_gwt_dc_ff, v_gwt_dc_ff, t_gwt_dc_ff, &
       flx_heat)

       u_gwt_ff = u_gwt_ff + u_gwt_dc_ff
       v_gwt_ff = v_gwt_ff + v_gwt_dc_ff
       t_gwt_ff = t_gwt_ff + t_gwt_dc_ff


     effgworo_ff=0.125
     call gw_oro_ifc( oro_band, &
       pcols, pver, dt_ff , effgworo_ff,  &
       u_dev_ff , v_dev_ff, t_dev_ff, &
       pint_dev_ff, pmid_dev_ff, & 
       pdel_dev_ff , rpdel_dev_ff, lnpint_dev_ff, &
       zm_dev_ff, zi, &
       nm, ni, rhoi, kvtt,  &
       sgh_dev_ff   ,rlat_dev_ff, &
       u_gwt_org_ff, v_gwt_org_ff, t_gwt_org_ff, &
       flx_heat)


       u_gwt_ff = u_gwt_ff + u_gwt_org_ff
       v_gwt_ff = v_gwt_ff + v_gwt_org_ff
       t_gwt_ff = t_gwt_ff + t_gwt_org_ff


     dudt_gwd_dev(1:pcols,1:pver) = REAL( u_gwt_ff(1:pcols,1:pver))  !zonal wind tendency at layer 
     dvdt_gwd_dev(1:pcols,1:pver) = REAL( v_gwt_ff(1:pcols,1:pver))  !meridional wind tendency at layer 
     dtdt_gwd_dev(1:pcols,1:pver) = REAL( t_gwt_ff(1:pcols,1:pver))  !temperature tendency at layer
     dudt_org_dev(1:pcols,1:pver) = REAL( u_gwt_org_ff(1:pcols,1:pver))  !zonal wind tendency at layer due to orography GWD
     dvdt_org_dev(1:pcols,1:pver) = REAL( v_gwt_org_ff(1:pcols,1:pver))  !meridional wind tendency at layer  due to orography GWD
     dtdt_org_dev(1:pcols,1:pver) = REAL( t_gwt_org_ff(1:pcols,1:pver))  !temperature tendency at layer  due to orography GWD



     taugwdx_dev(1:pcols)         = 0.0  !zonal      gravity wave surface    stress
     taugwdy_dev(1:pcols)         = 0.0  !meridional gravity wave surface    stress
     tauox_dev(1:pcols,1:pver+1)  = 0.0  !zonal      orographic gravity wave stress
     tauoy_dev(1:pcols,1:pver+1)  = 0.0  !meridional orographic gravity wave stress
     feo_dev  (1:pcols,1:pver+1)  = 0.0  !energy flux of orographic gravity waves
     fepo_dev (1:pcols,1:pver+1)  = 0.0  !pseudoenergy flux of orographic gravity waves
     taubkgx_dev(1:pcols)         = 0.0  !zonal      gravity wave background stress
     taubkgy_dev(1:pcols)         = 0.0  !meridional gravity wave background stress
     taubx_dev(1:pcols,1:pver+1)  = 0.0  !zonal      background gravity wave stress
     tauby_dev(1:pcols,1:pver+1)  = 0.0  !meridional background gravity wave stress
     feb_dev  (1:pcols,1:pver+1)  = 0.0  !energy flux of background gravity waves
     fepb_dev (1:pcols,1:pver+1)  = 0.0  !pseudoenergy flux of background gravity waves
     utbsrc_dev(1:pcols,1:pver)   = 0.0  !dU/dt below background launch level
     vtbsrc_dev(1:pcols,1:pver)   = 0.0  !dV/dt below background launch level
     ttbsrc_dev(1:pcols,1:pver)   = 0.0  !dT/dt below background launch level


    return
  end subroutine gw_intr_ncar

!================================================================================
  subroutine gw_prof_geos (i, k, pcols, pver, u, v, t, pm, pi, rhoi, ni, ti, nm)
!-----------------------------------------------------------------------
! Compute profiles of background state quantities for the multiple
! gravity wave drag parameterization.
! 
! The parameterization is assumed to operate only where water vapor 
! concentrations are negligible in determining the density.
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------
    integer, intent(in)  :: i                  ! current atmospheric column
    integer, intent(in)  :: k                  ! current atmospheric layer
    integer, intent(in)  :: pcols              ! number of atmospheric columns
    integer, intent(in)  :: pver               ! number of vertical layers

    real,    intent(in)  :: u(pcols,pver)      ! midpoint zonal wind
    real,    intent(in)  :: v(pcols,pver)      ! midpoint meridional wind
    real,    intent(in)  :: t(pcols,pver)      ! midpoint temperatures
    real,    intent(in)  :: pm(pcols,pver)     ! midpoint pressures
    real,    intent(in)  :: pi(pcols,0:pver)   ! interface pressures

    real,    intent(out) :: rhoi(0:pver)       ! interface density
    real,    intent(out) :: ni(0:pver)         ! interface Brunt-Vaisalla frequency
    real,    intent(out) :: ti(0:pver)         ! interface temperature
    real,    intent(out) :: nm(pver)           ! midpoint Brunt-Vaisalla frequency

!---------------------------Local storage-------------------------------

    real(r8)    :: dtdp
    real(r8)    :: n2                              ! Brunt-Vaisalla frequency squared

!-----------------------------------------------------------------------------
! Determine the interface densities and Brunt-Vaisala frequencies.
!-----------------------------------------------------------------------------

! The top interface values are calculated assuming an isothermal atmosphere 
! above the top level.
    if (k == 0) then
       ti(k)   = t(i,k+1)
       rhoi(k) = pi(i,k) / (MAPL_RGAS*ti(k))
       ni(k)   = sqrt (MAPL_GRAV*MAPL_GRAV / (MAPL_CP*ti(k)))

! Interior points use centered differences
    else if (k > 0 .and. k < pver) then
       ti(k)   = 0.5 * (t(i,k) + t(i,k+1))
       rhoi(k) = pi(i,k) / (MAPL_RGAS*ti(k))
       dtdp    = (t(i,k+1)-t(i,k)) / (pm(i,k+1)-pm(i,k))
       n2      = MAPL_GRAV*MAPL_GRAV/ti(k) * (1./MAPL_CP - rhoi(k)*dtdp)
       ni(k)   = sqrt (max (N2MIN, n2))

! Bottom interface uses bottom level temperature, density; next interface
! B-V frequency.
    else if (k == pver) then
       ti(k)   = t(i,k)
       rhoi(k) = pi(i,k) / (MAPL_RGAS*ti(k))
       ni(k)   = ni(k-1)
    end if

!-----------------------------------------------------------------------------
! Determine the midpoint Brunt-Vaisala frequencies.
!-----------------------------------------------------------------------------
    if (k > 0) then
      nm(k) = 0.5 * (ni(k-1) + ni(k))
    end if

    return
  end subroutine gw_prof_geos

!================================================================================

  subroutine gw_oro_geos (i, pcols, pver, pgwv, &
       u, v, t, sgh, pm, pi, dpm, zm, nm,  &
       kldv, kldvmn, ksrc, ksrcmn, rdpldv, &
       tau, ubi, ubm, xv, yv, kbot, rlat)
!-----------------------------------------------------------------------
! Orographic source for multiple gravity wave drag parameterization.
! 
! The stress is returned for a single wave with c=0, over orography.
! For points where the orographic variance is small (including ocean),
! the returned stress is zero. 
!------------------------------Arguments--------------------------------
    integer, intent(in)  :: i                    ! number of current column
    integer, intent(in)  :: pcols                ! number of atmospheric columns
    integer, intent(in)  :: pver                 ! number of atmospheric columns
    integer, intent(in)  :: pgwv                 ! number of waves allowed

    real,    intent(in)  :: u(pcols,pver)        ! midpoint zonal wind
    real,    intent(in)  :: v(pcols,pver)        ! midpoint meridional wind
    real,    intent(in)  :: t(pcols,pver)        ! midpoint temperatures
    real,    intent(in)  :: sgh(pcols)           ! standard deviation of orography
    real,    intent(in)  :: pm(pcols,pver)       ! midpoint pressures
    real,    intent(in)  :: pi(pcols,0:pver)     ! interface pressures
    real,    intent(in)  :: dpm(pcols,pver)      ! midpoint delta p (pi(k)-pi(k-1))
    real,    intent(in)  :: zm(pcols,pver)       ! midpoint heights
    real,    intent(in)  :: nm(pver)             ! midpoint Brunt-Vaisalla frequency

    integer, intent(out) :: kldv                 ! top interface of low level stress div region
    integer, intent(out) :: kldvmn               ! min value of kldv
    integer, intent(out) :: ksrc                 ! index of top interface of source region
    integer, intent(out) :: ksrcmn               ! min value of ksrc

    real,    intent(out) :: rdpldv               ! 1/dp across low level divergence region
    real,    intent(out) :: tau(-pgwv:pgwv,0:pver)! wave Reynolds stress
    real,    intent(out) :: ubi(0:pver)          ! projection of wind at interfaces
    real,    intent(out) :: ubm(pver)            ! projection of wind at midpoints
    real,    intent(out) :: xv                   ! unit vectors of source wind (x)
    real,    intent(out) :: yv                   ! unit vectors of source wind (y)
    integer, intent(inout) :: kbot
    real,    intent(in)    :: rlat(pcols)

!---------------------------Local storage-------------------------------
    integer :: k                                 ! loop indexes

    real    :: ubsrc                             ! Source-layer basic-state wind
    real    :: hdsp                              ! surface streamline displacment height (2*sgh)
    real    :: sghmax                            ! max orographic sdv to use
    real    :: tauoro                            ! c=0 stress from orography
!    real    :: zldv                              ! top of the low level stress divergence region
    real    :: nsrc                              ! b-f frequency averaged over source region
    real    :: psrc                              ! interface pressure at top of source region
    real    :: rsrc                              ! density averaged over source region
    real    :: usrc                              ! u wind averaged over source region
    real    :: vsrc                              ! v wind averaged over source region

! Begins

!---------------------------------------------------------------------------
! Average the basic state variables for the wave source over the depth of
! the orographic standard deviation. Here we assume that the apropiate
! values of wind, stability, etc. for determining the wave source are 
! averages over the depth of the atmosphere pentrated by the typical mountain.
! Reduces to the bottom midpoint values when sgh=0, such as over ocean.
! 
! Also determine the depth of the low level stress divergence region, as
! the max of the boundary layer depth and the source region depth. This
! can be done here if the stress magnitude does not determine the depth,
! otherwise it must be done below.
!---------------------------------------------------------------------------

    ksrc = pver-1
    kldv = pver-1
    psrc = pi(i,pver-1)
    rsrc = pm(i,pver)/(MAPL_RGAS*t(i,pver)) * dpm(i,pver)
    usrc = u(i,pver) * dpm(i,pver)
    vsrc = v(i,pver) * dpm(i,pver)
    nsrc = nm(pver)* dpm(i,pver)
    hdsp = 2.0 * sgh(i)

    do k = pver-1, pver/2, -1
       if (hdsp > sqrt(zm(i,k)*zm(i,k+1))) then
          ksrc = k-1
          kldv = k-1
          psrc = pi(i,k-1)
          rsrc = rsrc + pm(i,k) / (MAPL_RGAS*t(i,k))* dpm(i,k)
          usrc = usrc + u(i,k) * dpm(i,k)
          vsrc = vsrc + v(i,k) * dpm(i,k)
          nsrc = nsrc + nm(k)* dpm(i,k)
       end if
    end do

    rsrc = rsrc / (pi(i,pver) - psrc)
    usrc = usrc / (pi(i,pver) - psrc)
    vsrc = vsrc / (pi(i,pver) - psrc)
    nsrc = nsrc / (pi(i,pver) - psrc)

    if ( usrc == 0. .and. vsrc == 0. ) then
       ubsrc = sqrt(UBMC2MN)
       xv = 1.
       yv = 0.
    else
       ubsrc = sqrt(usrc**2+vsrc**2)
       xv = usrc/ubsrc
       yv = vsrc/ubsrc
    end if

! Project the local wind at midpoints onto the source wind.
    do k = 1, pver
       ubm(k) = u(i,k) * xv + v(i,k) * yv
    end do

! Compute the interface wind projection by averaging the midpoint winds.
! Use the top level wind at the top interface.
    ubi(0) = ubm(1)
    do k = 1, pver
       ubi(k) = ubm(k)
    end do

! Determine the orographic c=0 source term following McFarlane (1987).
! Set the source top interface index to pver, if the orographic term is zero.
    if ((ubsrc .gt. OROVMIN) .and. (hdsp .gt. OROHMIN)) then
       sghmax = FCRIT2 * (ubsrc / nsrc)**2
       tauoro = OROKO2 * min(hdsp**2, sghmax) * rsrc * nsrc * ubsrc
    else
       tauoro = 0.
       ksrc   = pver
       kldv   = pver
    end if

! tauoro is nontrivial when ubsrc is positive. However, if ubi(ksrc) is negative
! [note that the sign of ubsrc is irrelevant to the sign of ubi(ksrc)], orographic
! GWs can propagative upward passing through the negative basic-state wind.
! The following is to prevent this physically unjustified simulation.
! In addition, even if ubi(ksrc) > 0, if ubm(ksrc) < 0 .and. ubi(ksrc-1) < 0,
! negative wave stress leads to the acceleration of the negative ubm(ksrc).
! This result is physically inconsistent. However, if ubm(ksrc) > 0 .and.
! ubi(ksrc-1) < 0, GWs are filtered in physically consistent way, and
! decelerate the positive ubm(ksrc). Therefore, GWs are also assumed not to be
! launched when ubm(ksrc) < 0.
    if (ubi(kbot) < 0. .or. ubm(kbot) < 0.) then
       tauoro = 0.
       ksrc   = pver
       kldv   = pver
    end if

! Sets kbot equal to ksrc
    kbot = ksrc

! Set the phase speeds and wave numbers in the direction of the source wind.
! Set the source stress magnitude (positive only, note that the sign of the 
! stress is the same as (c-u).
    tau(0,kbot) = tauoro
! +
! Find the top interface of the low level stress divergence region according
! to the maximum depth of three criterion.
! 1. source region depth
! 2. planetary boundary layer depth
! 3. 10 * (u_*) / N where u_* is defined from the gravity wave stresss
! = sqrt(tau/rho) using source region values
! -
!      if (kbot .lt. pver) then
!         kldv = kbot
!      else
!         zldv    = max (pblh(i), sgh(i)
!         zldv    = max (zdlv(i), ZLDVCON * sqrt(tau(0,k)/rsrc) / nsrc)
!         kldv = pver-1
!         do k = pver-1, pver/2, -1
!            if (zldv .gt. sqrt(zm(i,k)*zm(i,k+1))) kldv = k-1
!         end do
!      end if

! Determine the min value of kldv and ksrc for limiting later loops
! and the pressure at the top interface of the low level stress divergence
! region.

    ksrcmn = pver
    kldvmn = pver

    ksrcmn = min(ksrcmn, ksrc)
    kldvmn = min(kldvmn, kldv)
    
    if (kldv .ne. pver) then
       rdpldv = 1. / (pi(i,kldv) - pi(i,pver))
    end if
! kldvmn is always pver because FRACLDV == 0.
    if (FRACLDV .le. 0.) kldvmn = pver

    return
  end subroutine gw_oro_geos

!===============================================================================
  subroutine gw_bgnd (i, pcols, pver, cw,       &
       u, v, t, pm, pi, dpm, rdpm, piln, rlat,  &
       kldv, kldvmn, ksrc, ksrcmn, rdpldv, tau, &
       ubi, ubm, xv, yv, ngwv, kbot, bgstressmax)
!-----------------------------------------------------------------------
! Driver for multiple gravity wave drag parameterization.
! 
! The parameterization is assumed to operate only where water vapor 
! concentrations are negligible in determining the density.
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------
    integer, intent(in)  :: i                     ! number of current column
    integer, intent(in)  :: pcols                 ! number of atmospheric columns
    integer, intent(in)  :: pver                  ! number of atmospheric columns

    integer, intent(in)  :: kbot                  ! index of bottom (source) interface
    integer, intent(in)  :: ngwv                  ! number of gravity waves to use
    real,    intent(in)  :: cw(-ngwv:ngwv)        ! wave weights
    real,    intent(in)  :: u(pcols,pver)         ! midpoint zonal wind
    real,    intent(in)  :: v(pcols,pver)         ! midpoint meridional wind
    real,    intent(in)  :: t(pcols,pver)         ! midpoint temperatures
    real,    intent(in)  :: pm(pcols,pver)        ! midpoint pressures
    real,    intent(in)  :: pi(pcols,0:pver)      ! interface pressures
    real,    intent(in)  :: dpm(pcols,pver)       ! midpoint delta p (pi(k)-pi(k-1))
    real,    intent(in)  :: rdpm(pcols,pver)      ! 1. / (pi(k)-pi(k-1))
    real,    intent(in)  :: piln(pcols,0:pver)    ! ln(interface pressures)
    real,    intent(in)  :: rlat(pcols)           ! latitude in radians for columns

    integer, intent(out) :: kldv                  ! top interface of low level stress divergence region
    integer, intent(out) :: kldvmn                ! min value of kldv
    integer, intent(out) :: ksrc                  ! index of top interface of source region
    integer, intent(out) :: ksrcmn                ! min value of ksrc

    real,    intent(in)  :: rdpldv                ! 1/dp across low level divergence region
    real,    intent(out) :: tau(-ngwv:ngwv,0:pver)! wave Reynolds stress
    real,    intent(out) :: ubi(0:pver)           ! projection of wind at interfaces
    real,    intent(out) :: ubm(pver)             ! projection of wind at midpoints
    real,    intent(out) :: xv                    ! unit vectors of source wind (x)
    real,    intent(out) :: yv                    ! unit vectors of source wind (y)

    real,    intent(in)  :: bgstressmax           ! Max of equatorial profile of BG stress factor
!---------------------------Local storage-------------------------------
    integer :: k,l                                ! loop indexes

    real    :: tauback                            ! background stress at c=0
    real    :: usrc                               ! u wind averaged over source region
    real    :: vsrc                               ! v wind averaged over source region
    real    :: ubsrc          
    real    :: al0                                ! Used in lat dependence of GW spec. 
    real    :: dlat0                              ! Used in lat dependence of GW spec.
    real    :: latdeg           
    real    :: flat_gw                            ! The actual lat dependence of GW spec.

!---------------------------------------------------------------------------
! Determine the source layer wind and unit vectors, then project winds.
!---------------------------------------------------------------------------

! Just use the source level interface values for the source
! wind speed and direction (unit vector).

    ksrc  = kbot
    kldv  = kbot
    usrc  = 0.5*(u(i,kbot+1)+u(i,kbot))
    vsrc  = 0.5*(v(i,kbot+1)+v(i,kbot))
    ubsrc = max(sqrt (usrc**2 + vsrc**2), sqrt (UBMC2MN))
    if (usrc == 0. .and. vsrc == 0.) then
       xv = 1.0
       yv = 0.0
    else
       xv = usrc / ubsrc
       yv = vsrc / ubsrc
    end if

! Project the local wind at midpoints onto the source wind.
    do k = 1, pver
       ubm(k) = u(i,k) * xv + v(i,k) * yv
    end do

! Compute the bottom interface wind projection using the midpoint winds.
    ubi(0) = ubm(1)
    do k = 1, pver
       ubi(k) = ubm(k)
    end do

!-----------------------------------------------------------------------
! Gravity wave sources
!-----------------------------------------------------------------------

! Determine the background stress at c=0
    tauback = TAUBGND * TAUSCAL

! Include dependence on latitude:

    latdeg = rlat(i)*180./PI_GWD
!
    if (-15.3 < latdeg .and. latdeg < 15.3) then
!!AMM  flat_gw = 1.2*dexp(-dble((abs(latdeg)-3.)/8.0)**2) 
!!AMM  if (flat_gw < 1.2 .and. abs(latdeg) <= 3.) flat_gw = 1.2
!!AMM  flat_gw = 2.5*dexp(-dble((abs(latdeg)-3.)/8.0)**2) 
!!AMM  if (flat_gw < 2.5 .and. abs(latdeg) <= 3.) flat_gw = 2.5
       flat_gw = bgstressmax*dexp(-dble((abs(latdeg)-3.)/8.0)**2) 
       if (flat_gw < bgstressmax .and. abs(latdeg) <= 3.) flat_gw = bgstressmax
    else if (latdeg > -31. .and. latdeg <= -15.3) then
       flat_gw =  0.10
    else if (latdeg <  31. .and. latdeg >=  15.3) then
       flat_gw =  0.10
    else if (latdeg > -60. .and. latdeg <= -31.) then
       flat_gw =  0.50*dexp(-dble((abs(latdeg)-60.)/23.)**2)
    else if (latdeg <  60. .and. latdeg >=  31.) then
       flat_gw =  0.50*dexp(-dble((abs(latdeg)-60.)/23.)**2)
    else if (latdeg <= -60.) then
       flat_gw =  0.50*dexp(-dble((abs(latdeg)-60.)/70.)**2)
    else if (latdeg >=  60.) then
       flat_gw =  0.50*dexp(-dble((abs(latdeg)-60.)/70.)**2)
    end if
    tauback=tauback*flat_gw

! Set the phase speeds and wave numbers in the direction of the source wind.
! Set the source stress magnitude (positive only, note that the sign of the 
! stress is the same as (c-u).

    do l = 1, ngwv
       tau( l,kbot) = tauback * cw(l)
       tau(-l,kbot) = tau( l,kbot)
    end do
    tau(0,kbot) = tauback

! Determine the min value of kldv and ksrc for limiting later loops
! and the pressure at the top interface of the low level stress divergence
! region.

    ksrcmn = pver
    kldvmn = pver

    return
  end subroutine gw_bgnd

!===============================================================================
  subroutine gw_drag_prof_geos (i,     pcols, pver,      &
             pgwv,  ngwv,  kbot,  ktop,  c,     u,       &
             v,     t,     pi,    dpm,   rdpm,  piln,    &
             rlat,  rhoi,  ni,    ti,    nm,    dt,      &
             alpha, dback, kldv,  kldvmn,ksrc,  ksrcmn,  &
             rdpldv,tau,   ubi,   ubm,   xv,    yv,      & 
             ut,    vt,    tt,    taugwx,taugwy, fegw,   & 
             fepgw, dusrc, dvsrc, dtsrc,                 &
             tau0x, tau0y, effgw )
!-----------------------------------------------------------------------
! Solve for the drag profile from the multiple gravity wave drag
! parameterization.
! 1. scan up from the wave source to determine the stress profile
! 2. scan down the stress profile to determine the tendencies
!     => apply bounds to the tendency
!          a. from wkb solution
!          b. from computational stability constraints
!     => adjust stress on interface below to reflect actual bounded tendency
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------
    integer, intent(in) :: i                     ! number of current column
    integer, intent(in) :: pcols                 ! number of atmospheric columns
    integer, intent(in) :: pver                  ! number of atmospheric columns
    integer, intent(in) :: kbot                  ! index of bottom (source) interface
    integer, intent(in) :: ktop                  ! index of top interface of gwd region
    integer, intent(in) :: pgwv                  ! number of gravity waves possible
    integer, intent(in) :: ngwv                  ! number of gravity waves to use
    integer, intent(in) :: kldv                  ! top interface of low level stress  divergence region
    integer, intent(in) :: kldvmn                ! min value of kldv
    integer, intent(in) :: ksrc                  ! index of top interface of source region
    integer, intent(in) :: ksrcmn                ! min value of ksrc

    real,    intent(in) :: c(-pgwv:pgwv)         ! wave phase speeds
    real,    intent(in) :: u(pcols,pver)         ! midpoint zonal wind
    real,    intent(in) :: v(pcols,pver)         ! midpoint meridional wind
    real,    intent(in) :: t(pcols,pver)         ! midpoint temperatures
    real,    intent(in) :: pi(pcols,0:pver)      ! interface pressures
    real,    intent(in) :: dpm(pcols,pver)       ! midpoint delta p (pi(k)-pi(k-1))
    real,    intent(in) :: rdpm(pcols,pver)      ! 1. / (pi(k)-pi(k-1))
    real,    intent(in) :: piln(pcols,0:pver)    ! ln(interface pressures)
    real,    intent(in) :: rlat(pcols)
    real,    intent(in) :: rhoi(0:pver)          ! interface density
    real,    intent(in) :: ni(0:pver)            ! interface Brunt-Vaisalla frequency
    real,    intent(in) :: ti(0:pver)            ! interface temperature
    real,    intent(in) :: nm(pver)              ! midpoint Brunt-Vaisalla frequency
    real,    intent(in) :: dt                    ! time step
    real,    intent(in) :: alpha(0:pver)         ! newtonian cooling coefficients
    real,    intent(in) :: dback(0:pver)         ! newtonian cooling coefficients
    real,    intent(in) :: rdpldv                ! 1/dp across low level divergence region
    real,    intent(in) :: ubi(0:pver)           ! projection of wind at interfaces
    real,    intent(in) :: ubm(pver)             ! projection of wind at midpoints
    real,    intent(in) :: xv                    ! unit vectors of source wind (x)
    real,    intent(in) :: yv                    ! unit vectors of source wind (y)
    real,    intent(in) :: effgw                 ! tendency efficiency for gwd 

    real,    intent(inout) :: tau(-pgwv:pgwv,0:pver)! wave Reynolds stress

    real,    intent(out) :: ut(pver)             ! zonal wind tendency
    real,    intent(out) :: vt(pver)             ! meridional wind tendency
    real,    intent(out) :: tt(pver)             ! temperature tendency
    real,    intent(out) :: taugwx(pcols,0:pver) ! Total zonal GW momentum flux
    real,    intent(out) :: taugwy(pcols,0:pver) ! Total meridional GW momentum flux
    real,    intent(out) :: fegw (pcols,0:pver)  ! Total GW energy flux
    real,    intent(out) :: fepgw(pcols,0:pver)  ! Total GW pseudo energy flux
    real,    intent(out) :: dusrc(pver)          ! Total U tendency below launch level
    real,    intent(out) :: dvsrc(pver)          ! Total V tendency below launch level
    real,    intent(out) :: dtsrc(pver)          ! Total V tendency below launch level
    real,    intent(out) :: tau0x                ! c=0 sfc. stress (zonal)
    real,    intent(out) :: tau0y                ! c=0 sfc. stress (meridional)
!---------------------------Local storage-------------------------------
    integer :: k,l                               ! loop indexes

    real    :: d !MATMAT Is this used?           ! "total" diffusivity 
    real    :: dsat !MATMAT Is this used?        ! saturation diffusivity
    real    :: dscal                             ! fraction of dsat to use
    real    :: mi                                ! imaginary part of vertical wavenumber
    real    :: taudmp                            ! stress after damping
    !real    :: taumax(pcols)                     ! max(tau) for any l
    real    :: tausat                            ! saturation stress
    real    :: ubmc                              ! (ub-c)
    real    :: ubmc2                             ! (ub-c)**2
    real    :: ubt                               ! ubar tendency
    real    :: tbt                               ! tbar tendency
    real    :: ubtl                              ! ubar tendency from wave l
    real    :: ubtlsat                           ! saturation tendency

    real    :: pm                                ! layer pressure
    real    :: rhom                              ! layer density
    real    :: zlb                               ! launch level height
    real    :: cmu                               ! c-u
    real    :: dzm, hscal, tautmp
    real    :: utl
    real    :: ttl
    real    :: fpmx                              ! zonal pseudomomentum flux spectrum
    real    :: fpmy                              ! meridional pseudomomentum flux spectrum
    real    :: fe                                ! energy flux (p'w') spectrum
    real    :: fpe                               ! pseudoenergy flux (p'w'+U rho u'w') spectrum
    real    :: fpml
    real    :: fpmt
    real    :: fpel
    real    :: fpet
    real    :: dusrcl
    real    :: dvsrcl
    real    :: dtsrcl

    real    :: zi
    real    :: effkwvmap
    real    :: zfac
    real    :: uhtmax
    real    :: utfac


! Initialize gravity wave drag tendencies to zero

    do k=1,pver
       ut(k)    = 0.
       vt(k)    = 0.
       tt(k)    = 0.
       dusrc(k) = 0.
       dvsrc(k) = 0.
       dtsrc(k) = 0.
    end do

! Initialize total momentum and energy fluxes to zero

    do k=0,pver
       taugwx(i,k) = 0. 
       taugwy(i,k) = 0. 
       fegw  (i,k) = 0. 
       fepgw (i,k) = 0.
    end do

! Initialize surface wave stress at c = 0 to zero

    tau0x = 0.
    tau0y = 0.

!---------------------------------------------------------------------------
! Compute the stress profiles and diffusivities
!---------------------------------------------------------------------------

! Determine the absolute value of the saturation stress and the diffusivity
! for each wave.
! Define critical levels where the sign of (u-c) changes between interfaces.

! Loop from bottom to top to get stress profiles
    do l = -ngwv, ngwv
       do k = pver-1, ktop, -1
          if (k <= kbot-1) then
             d = dback(k)
             ubmc = ubi(k) - c(l)

             if ( ngwv > 0 ) then
                effkwvmap = FCRIT2*KWVB*  &
                            (0.5*(1.0+tanh( (rlat(i)*180./PI_GWD-20.0)/6.0)) +  &
                             0.5*(1.0+tanh(-(rlat(i)*180./PI_GWD+20.0)/6.0)))
                if (-15.0 < rlat(i)*180./PI_GWD .and. rlat(i)*180./PI_GWD < 15.0) then
                   effkwvmap = FCRIT2*KWVBEQ
                end if
             else
                if (pi(i,k) < 1000.0) then
                   zfac = (pi(i,k)/1000.0)**3
                else
                   zfac = 1.0
                end if
                effkwvmap = FCRIT2*KWVO*zfac
             end if

             tausat = abs (effkwvmap * rhoi(k) * ubmc**3 / (2.*ni(k)) )
             if (tausat .le. TAUMIN) tausat = 0.0
             if (ubmc * (ubi(k+1) - c(l)) .le. 0.0) tausat = 0.0
             if (k .eq. ktop) tausat = 0.
!
             if (k == ktop+1) tausat = tausat*0.02
             if (k == ktop+2) tausat = tausat*0.05
             if (k == ktop+3) tausat = tausat*0.10
             if (k == ktop+4) tausat = tausat*0.20
             if (k == ktop+5) tausat = tausat*0.50
!
             tau(l,k) = min (tau(l,k+1), tausat)
             dsat = (ubmc / ni(k))**2 * &
                  (effkwvmap * ubmc**2 / (2. * ROG * ti(k) * ni(k)) - alpha(k))
             if ( tau(l,k+1) .ge. tausat ) then
                d = dsat
             else
                d = 0.
             end if
          end if
       end do

! The orographic stress term does not change across the source region

!         if (ngwv == 0 .and. k .ge. ksrcmn) then
!            if (k .ge. ksrc) then
!               tau(0,k) = tau(0,kbot)
!            end if
!         end if

! Require that the orographic term decrease linearly (with pressure) 
! within the low level stress divergence region. This supersedes the above
! requirment of constant stress within the source region.
! Note that k ge kldvmn cannot occur without an orographic source term, since
! kldvmn=pver then and k<=pver-1

!         if (ngwv == 0 .and. k .ge. kldvmn) then
!            if (k .ge. kldv) then
!               tau(0,k) = min (tau(0,k), tau(0,kbot)  * &
!                    (1. - FRACLDV * (pi(i,k)-pi(i,pver)) * rdpldv))
!            end if
!         end if

    end do

!---------------------------------------------------------------------------
! Compute the tendencies from the stress divergence.
!---------------------------------------------------------------------------

! Accumulate the mean wind tendency over wavenumber.

! Loop over levels from top to bottom
    do k = ktop+1, pver

       ubt = 0.0
       tbt = 0.0
       do l = -ngwv, ngwv
          if (k <= kbot) then

! Determine the wind tendency including excess stress carried down from above.
             ubtl = MAPL_GRAV * (tau(l,k)-tau(l,k-1)) * rdpm(i,k)

! Calculate the sign of wind tendency
             utl = sign(ubtl, c(l)-ubi(k))

! Accumulate the mean wind tendency over wavenumber.
             ubt = ubt + utl

! Calculate irreversible temperature tendency associated with gravity wave breaking.
             ttl = (c(l)-ubm(k))*utl/MAPL_CP

! Adding frictional heating associated with the GW momentum forcing
             tbt = tbt + ttl
          end if
       end do

! Project the mean wind tendency onto the components and scale by "efficiency".
       if (k <= kbot) then
          ut(k) = ubt * xv * effgw
          vt(k) = ubt * yv * effgw
          tt(k) = tbt      * effgw
       end if
    end do

!-----------------------------------------------------------------------
! Calculates wind and temperature tendencies below launch level for
! energy and momentum conservation (does nothing for orographic GWs).
!-----------------------------------------------------------------------

! Calculate launch level height
    zlb = 0.

    do k = ktop+1, pver
       if (k >= kbot+1) then

! Define layer pressure and density
          pm   = (pi(i,k-1)+pi(i,k))*0.5
          rhom = pm/(MAPL_RGAS*t(i,k))

          zlb  = zlb + dpm(i,k)/MAPL_GRAV/rhom
       end if
    end do

   !-----------------------------------------------------------------------
   ! Calculates energy and momentum flux profiles
   !-----------------------------------------------------------------------

    do l = -ngwv, ngwv
       do k = ktop, pver
          if ( k <= kbot ) then
             cmu  = c(l)-ubi(k)
             fpmx =      sign(1.0,cmu)*tau(l,k)*xv*effgw
             fpmy =      sign(1.0,cmu)*tau(l,k)*yv*effgw
             fe   =  cmu*sign(1.0,cmu)*tau(l,k)*effgw
             fpe  = c(l)*sign(1.0,cmu)*tau(l,k)*effgw

             if (k == kbot) then
                fpml = fpmx*xv+fpmy*yv
                fpel = fpe
             end if

             if (k == ktop) then
                fpmt = fpmx*xv+fpmy*yv
                fpet = fpe
             end if

             ! Record outputs for GW fluxes
             taugwx(i,k) = taugwx(i,k) + fpmx
             taugwy(i,k) = taugwy(i,k) + fpmy
             fegw  (i,k) = fegw  (i,k) + fe
             fepgw (i,k) = fepgw (i,k) + fpe
          end if
       end do

       do k = ktop+1, pver
          if (k >= kbot+1) then

! Define layer pressure and density
             pm   = (pi(i,k-1)+pi(i,k))*0.5
             rhom = pm/(MAPL_RGAS*t(i,k))

             dusrcl = - (fpml-fpmt)/(rhom*zlb)*xv
             dvsrcl = - (fpml-fpmt)/(rhom*zlb)*yv
             dtsrcl = -((fpel-fpet)-ubm(k)*(fpml-fpmt))/  &
                              (rhom*zlb*MAPL_CP)

! Add sub-source wind and temperature tendencies
             dusrc(k) = dusrc(k) + dusrcl
             dvsrc(k) = dvsrc(k) + dvsrcl
             dtsrc(k) = dtsrc(k) + dtsrcl
          end if
       end do
    end do

! For orographic waves, sub-source tendencies are set equal to zero.
    do k = 1, pver
       if (ngwv == 0) then
          dusrc(k) = 0.0
          dvsrc(k) = 0.0
          dtsrc(k) = 0.0
       end if
    end do

!-----------------------------------------------------------------------
! Adjust efficiency factor to prevent unrealistically strong forcing
!-----------------------------------------------------------------------

    uhtmax = 0.0 
    utfac  = 1.0

    do k = 1, pver
       uhtmax = max(sqrt(ut(k)**2 + vt(k)**2), uhtmax)
    end do

    if (uhtmax > TNDMAX) utfac = TNDMAX/uhtmax

    do k = 1, pver
       ut   (k) = ut   (k)*utfac
       vt   (k) = vt   (k)*utfac
       tt   (k) = tt   (k)*utfac
       dusrc(k) = dusrc(k)*utfac
       dvsrc(k) = dvsrc(k)*utfac
       dtsrc(k) = dtsrc(k)*utfac
    end do

!-----------------------------------------------------------------------
! Project the c=0 stress (scaled) in the direction of the source wind
! for recording on the output file.
!-----------------------------------------------------------------------
    tau0x = tau(0,kbot) * xv * effgw*utfac
    tau0y = tau(0,kbot) * yv * effgw*utfac

    return
  end subroutine gw_drag_prof_geos

end module gw_drag_ncar

