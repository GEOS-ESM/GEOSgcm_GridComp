
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

  use MAPL_ConstantsMod, only: MAPL_RGAS, MAPL_GRAV

  use gw_rdg, only     : gw_rdg_ifc
  use gw_oro, only     : gw_oro_ifc
  use gw_convect, only : BeresSourceDesc, gw_beres_ifc
  use gw_common, only  : GWBand,gw_prof

  !save
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


  subroutine gw_intr_ncar(pcols,      pver,         dt,         nrdg,              &    
          beres_dc_desc, beres_sc_desc, beres_band,   oro_band,                      &
          pint_dev,      t_dev,         u_dev,        v_dev,                         &
          ht_dc_dev,     ht_sc_dev,     dqcdt_dev,                                   &
          sgh_dev,       mxdis_dev,     hwdth_dev,    clngt_dev,  angll_dev,         &
          anixy_dev,     gbxar_dev,     kwvrdg_dev,   effrdg_dev, pref_dev,          & 
          pmid_dev,      pdel_dev,      rpdel_dev,    lnpint_dev, zm_dev,  rlat_dev, &
          phis_dev,                                                                  &
          dudt_gwd_dev,  dvdt_gwd_dev,  dtdt_gwd_dev,                                &
          dudt_org_dev,  dvdt_org_dev,  dtdt_org_dev,                                &
          taugwdx_dev,   taugwdy_dev,   &
          taubkgx_dev,   taubkgy_dev,   &
          effgworo,      effgwbkg,      rc            )

!-----------------------------------------------------------------------
! Interface for multiple gravity wave drag parameterization.
!-----------------------------------------------------------------------

!------------------------------Arguments--------------------------------
    integer, intent(in   ) :: pcols                    ! number of columns
    integer, intent(in   ) :: pver                     ! number of vertical layers
!++jtb 01/25/21
    integer, intent(in   ) :: nrdg                     ! number of Ridges per gridbox
!--jtb
    real,    intent(in   ) :: dt                       ! time step
    type(GWBand),          intent(inout) :: oro_band   ! Band descriptor
    type(GWBand),          intent(inout) :: beres_band ! Band descriptor
    type(BeresSourceDesc), intent(inout) :: beres_dc_desc ! Table descriptor for DeepCu Beres scheme
    type(BeresSourceDesc), intent(inout) :: beres_sc_desc ! Table descriptor for ShallowCu Beres scheme
    real,    intent(in   ) :: effgwbkg                 ! tendency efficiency for background gwd (Default = 0.125)
    real,    intent(in   ) :: effgworo                 ! tendency efficiency for orographic gwd (Default = 0.125)
    real,    intent(in   ) :: pint_dev(pcols,pver+1)   ! pressure at the layer edges
    real,    intent(in   ) :: t_dev(pcols,pver)        ! temperature at layers
    real,    intent(in   ) :: u_dev(pcols,pver)        ! zonal wind at layers
    real,    intent(in   ) :: v_dev(pcols,pver)        ! meridional wind at layers
    real,    intent(in   ) :: ht_dc_dev(pcols,pver)    ! DeepCu heating in layers
    real,    intent(in   ) :: ht_sc_dev(pcols,pver)    ! ShallowCu heating in layers
    real,    intent(in   ) :: dqcdt_dev(pcols,pver)    ! Condensate tendencies due to large-scale
    real,    intent(in   ) :: sgh_dev(pcols)           ! standard deviation of orography
!++jtb 01/25/21 New topo vars
    real,    intent(in   ) :: mxdis_dev(pcols,nrdg)     ! obstacle/ridge height 
    real,    intent(in   ) :: hwdth_dev(pcols,nrdg)     ! obstacle width
    real,    intent(in   ) :: clngt_dev(pcols,nrdg)     ! obstacle along-crest length
    real,    intent(in   ) :: angll_dev(pcols,nrdg)     ! obstacle orientation
    real,    intent(in   ) :: anixy_dev(pcols,nrdg)     ! obstacle ansitropy param 
    real,    intent(in   ) :: gbxar_dev(pcols)          ! duplicate grid box area
    real,    intent(in   ) :: kwvrdg_dev(pcols,nrdg)    ! horizontal wavenumber
    real,    intent(in   ) :: effrdg_dev(pcols,nrdg)    ! efficiency of ridge scheme
!!--jtb
    real,    intent(in   ) :: pref_dev(pver+1)         ! reference pressure at the layeredges
    real,    intent(in   ) :: pmid_dev(pcols,pver)     ! pressure at the layers
    real,    intent(in   ) :: pdel_dev(pcols,pver)     ! pressure thickness at the layers
    real,    intent(in   ) :: rpdel_dev(pcols,pver)    ! 1.0 / pdel
    real,    intent(in   ) :: lnpint_dev(pcols,pver+1) ! log(pint)
    real,    intent(in   ) :: zm_dev(pcols,pver)       ! height above surface at layers
    real,    intent(in   ) :: rlat_dev(pcols)          ! latitude in radian
    real,    intent(in   ) :: phis_dev(pcols)          ! surface geopotential
 
    real,    intent(  out) :: dudt_gwd_dev(pcols,pver) ! zonal wind tendency at layer 
    real,    intent(  out) :: dvdt_gwd_dev(pcols,pver) ! meridional wind tendency at layer 
    real,    intent(  out) :: dtdt_gwd_dev(pcols,pver) ! temperature tendency at layer
    real,    intent(  out) :: dudt_org_dev(pcols,pver) ! zonal wind tendency at layer due to orography GWD
    real,    intent(  out) :: dvdt_org_dev(pcols,pver) ! meridional wind tendency at layer  due to orography GWD
    real,    intent(  out) :: dtdt_org_dev(pcols,pver) ! temperature tendency at layer  due to orography GWD
    real,    intent(  out) :: taugwdx_dev(pcols)       ! zonal      gravity wave surface    stress
    real,    intent(  out) :: taugwdy_dev(pcols)       ! meridional gravity wave surface    stress
    real,    intent(  out) :: taubkgx_dev(pcols)       ! zonal      gravity wave background stress
    real,    intent(  out) :: taubkgy_dev(pcols)       ! meridional gravity wave background stress

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

    real :: ttgw(pcols,pver)            ! temperature tendency
    real :: utgw(pcols,pver)            ! zonal wind tendency
    real :: vtgw(pcols,pver)            ! meridional wind tendency

    real :: z(pcols)                           ! Surface elevation
    real :: zi(pcols,pver+1)                   ! interface heights above ground
    real :: ni(pcols,pver+1)                   ! interface Brunt-Vaisalla frequency
    real :: nm(pcols,pver)                     ! midpoint Brunt-Vaisalla frequency
    !!!real    :: rdpldv                              ! 1/dp across low level divergence region
    real :: rhoi(pcols,pver+1)                 ! interface density
    real :: ti(pcols,pver+1)                   ! interface temperature
    real :: ubi(pcols,pver+1)                  ! projection of wind at interfaces
    real :: ubm(pcols,pver)                    ! projection of wind at midpoints
    real :: xv(pcols)                          ! unit vectors of source wind (x)
    real :: yv(pcols)                          ! unit vectors of source wind (y)
    real :: kvtt(pcols,pver+1)                 ! Molecular thermal diffusivity.

    real :: maxq0(pcols),hdepth(pcols)
    real :: flx_heat(pcols)

    real         :: rdg_cd_llb
    integer      :: pverp, pcnst
    logical      :: trpd_leewv

!-----------------------------------------------------------------------------

! Misc dimensions needed
    pverp=pver+1
    pcnst=1

! Initialize accumulated tendencies
! and other things ...
    dtdt_gwd_dev(:,:) = 0.
    dudt_gwd_dev(:,:) = 0.
    dvdt_gwd_dev(:,:) = 0.
    kvtt(:,:)  = 0.
    flx_heat(:) = 0.0

    call gw_prof (pcols , pver, pint_dev , pmid_dev , t_dev , rhoi, nm, ni )

    z = phis_dev/MAPL_GRAV

    zi(:,pver+1) = 0.0
    do k=2,pver 
       zi(:,k)  =  0.5 * ( zm_dev(:,k-1) + zm_dev(:,k) )  
    end do
    zi(:,1) = zi(:,2) + 0.5*( zm_dev(:,1) - zm_dev(:,2)  )

   ! DeepCu BKG
    if (beres_dc_desc%active .and. effgwbkg>0.0) then
    call gw_beres_ifc( beres_band, &
       pcols, pver, dt , effgwbkg,  &
       u_dev , v_dev, t_dev, &
       pref_dev, pint_dev, & 
       pdel_dev , rpdel_dev, lnpint_dev, &
       zm_dev, zi, &
       nm, ni, rhoi, kvtt,  &
       dqcdt_dev, &
       ht_dc_dev,beres_dc_desc,rlat_dev, &
       utgw, vtgw, ttgw, flx_heat)
       dudt_gwd_dev = dudt_gwd_dev + utgw
       dvdt_gwd_dev = dvdt_gwd_dev + vtgw
       dtdt_gwd_dev = dtdt_gwd_dev + ttgw
    endif
   ! ShallowCu BKG
    if (beres_sc_desc%active .and. effgwbkg>0.0) then
    call gw_beres_ifc( beres_band, &
       pcols, pver, dt , effgwbkg,  &
       u_dev , v_dev, t_dev, &
       pref_dev, pint_dev, &
       pdel_dev , rpdel_dev, lnpint_dev, &
       zm_dev, zi, &
       nm, ni, rhoi, kvtt,  &
       dqcdt_dev, &
       ht_sc_dev,beres_sc_desc,rlat_dev, &
       utgw, vtgw, ttgw, flx_heat)
       dudt_gwd_dev = dudt_gwd_dev + utgw
       dvdt_gwd_dev = dvdt_gwd_dev + vtgw
       dtdt_gwd_dev = dtdt_gwd_dev + ttgw
     endif
    ! Orographic
     if (effgworo > 0.0) then
     if (nrdg > 0) then
       trpd_leewv    = .FALSE.
       rdg_cd_llb    = 1.0
       call gw_rdg_ifc( &
         pcols, pver, pverp, pcnst, nrdg, dt, &
         u_dev , v_dev, t_dev, &
         pint_dev, pmid_dev, &
         pdel_dev, rpdel_dev, &
         lnpint_dev, zm_dev, zi, z, &
         ni, nm, rhoi, &
         kvtt, &
         kwvrdg_dev, effrdg_dev, &
         hwdth_dev, clngt_dev, gbxar_dev, &
         mxdis_dev, angll_dev, anixy_dev, &
         rdg_cd_llb, trpd_leewv, &
         utgw, vtgw, ttgw, flx_heat)
     else
       call gw_oro_ifc( oro_band, &
         pcols, pver, dt , effgworo,  &
         u_dev , v_dev, t_dev, &
         pint_dev, pmid_dev, &
         pdel_dev , rpdel_dev, lnpint_dev, &
         zm_dev, zi, &
         nm, ni, rhoi, kvtt,  &
         sgh_dev   ,rlat_dev, &
         utgw, vtgw, ttgw)
     endif
     dudt_org_dev = utgw
     dvdt_org_dev = vtgw
     dtdt_org_dev = ttgw
     dudt_gwd_dev = dudt_gwd_dev + dudt_org_dev
     dvdt_gwd_dev = dvdt_gwd_dev + dvdt_org_dev
     dtdt_gwd_dev = dtdt_gwd_dev + dtdt_org_dev
     endif

     taugwdx_dev(1:pcols)         = 0.0  !zonal      gravity wave surface    stress
     taugwdy_dev(1:pcols)         = 0.0  !meridional gravity wave surface    stress
     taubkgx_dev(1:pcols)         = 0.0  !zonal      gravity wave background stress
     taubkgy_dev(1:pcols)         = 0.0  !meridional gravity wave background stress

    return
  end subroutine gw_intr_ncar

end module gw_drag_ncar

