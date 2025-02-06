module shoc

! Implementation of the Simplified High Order Closure (SHOC) scheme
! of Bogenschutz and Krueger (2013), J. Adv. Model. Earth Syst, 5, 195-211,
! doi: 10.1002/jame.200118. (further referred to as BK13)
! in a single column form suitable for use in a GCM physics package.
! Alex Belochitski, heavily based on the code of Peter Bogenschutz.
! S Moorthi - optimization, cleanup, improve and customize for gsm
!
! N Arnold - implemented in GEOS
!

 use MAPL_ConstantsMod, only: ggr    => MAPL_GRAV,   &
                              cp     => MAPL_CP,     &
                              rgas   => MAPL_RGAS,   &
                              rv     => MAPL_RVAP,   &
                              lcond  => MAPL_ALHL,   &
                              lfus   => MAPL_ALHS,   &
                              pi     => MAPL_PI,     &
                              MAPL_H2OMW, MAPL_AIRMW

 use MAPL,              only: MAPL_UNDEF

 use MAPL_SatVaporMod,  only: MAPL_EQsat

 implicit none

 type SHOCPARAMS_TYPE
    integer :: LENOPT
    integer :: BUOYOPT
    real    :: PRNUM
    real    :: LAMBDA
    real    :: TSCALE
    real    :: CKVAL
    real    :: CEFAC
    real    :: CESFAC
    real    :: LENFAC1
    real    :: LENFAC2
    real    :: LENFAC3
 endtype SHOCPARAMS_TYPE
 type (SHOCPARAMS_TYPE) :: shocparams

 private

 public run_shoc, update_moments, shocparams

 contains

 subroutine run_shoc( nx, ny, nzm, nz, dtn,                      &  ! in
                 prsl_inv, phii_inv, phil_inv,                   &  ! in
                 u_inv, v_inv,                                   &  ! in
                 omega_inv, tabs_inv, qwv_inv,                   &  ! in
                 qi_inv, qc_inv, qpi_inv,                        &  ! in
                 qpl_inv, cld_sgs_inv, wthv_sec_inv,             &  ! in
                 wthv_mf_inv, tke_mf, zpbl,                      &  ! in
                 tke_inv, tkh_inv,                               &  ! inout
                 tkm_inv, isotropy_inv,                          &  ! out
                 tkesbdiss_inv, tkesbbuoy_inv,                   &  ! out
                 tkesbshear_inv,                                 &  ! out
                 smixt_inv, lmix_out, smixt1_inv,                &  ! out
                 smixt2_inv,smixt3_inv,                          &  ! out
                 bruntmst_inv, ri_inv, prnum_inv,                &  ! out
!                 bruntmst_inv, bruntdry_inv, bruntedg_inv, ri_inv, prnum_inv,  &  ! out
                 shocparams )


  real, parameter :: lsub = lcond+lfus,         &
                     fac_cond = lcond/cp,       &
                     fac_fus = lfus/cp,         &
                     fac_sub = lsub/cp,         &
                     ggri = 1.0/ggr,            &
                     kapa = rgas/cp,            &
                     gocp = ggr/cp,             &
                     rog = rgas*ggri,           &
                     tkef1=0.5,                 &
                     tkef2=1.0-tkef1,           &
                     tkhmax=1000.0,             &
                     epsv=MAPL_H2OMW/MAPL_AIRMW

  integer, intent(in) :: nx    ! Number of points in the physics window in the x
  integer, intent(in) :: ny    ! and y directions

  integer, intent(in) :: nzm   ! Number of vertical layers
  integer, intent(in) :: nz    ! Number of layer interfaces  (= nzm + 1)
  type(shocparams_type), intent(in) :: shocparams
  real, intent(in   ) :: dtn          ! Physics time step, s
  real, intent(in   ) :: prsl_inv (nx,ny,nzm) ! mean layer presure
  real, intent(in   ) :: phii_inv (nx,ny,nz ) ! interface geopotential height
  real, intent(in   ) :: phil_inv (nx,ny,nzm) ! layer geopotential height
  real, intent(in   ) :: u_inv    (nx,ny,nzm) ! u-wind, m/s
  real, intent(in   ) :: v_inv    (nx,ny,nzm) ! v-wind, m/s
  real, intent(in   ) :: omega_inv(nx,ny,nzm) ! omega, Pa/s
  real, intent(in   ) :: wthv_sec_inv(nx,ny,nzm) ! Buoyancy flux, K*m/s
  real, intent(in   ) :: wthv_mf_inv(nx,ny,nzm) ! Buoyancy flux, K*m/s

  real, intent(in   ) :: tabs_inv   (nx,ny,nzm) ! temperature, K
  real, intent(in   ) :: qwv_inv    (nx,ny,nzm) ! water vapor mixing ratio, kg/kg
  real, intent(in   ) :: qc_inv     (nx,ny,nzm) ! cloud water mixing ratio, kg/kg
  real, intent(in   ) :: qi_inv     (nx,ny,nzm) ! cloud ice   mixing ratio, kg/kg
  real, intent(in   ) :: qpl_inv    (nx,ny,nzm) ! rain mixing ratio, kg/kg
  real, intent(in   ) :: qpi_inv    (nx,ny,nzm) ! snow mixing ratio, kg/kg
  real, intent(in   ) :: cld_sgs_inv(nx,ny,nzm) ! sgs cloud fraction
!  real, intent(in   ) :: dtdtrad    (nx,ny,nzm) ! radiative cooling tendency
  real, intent(in   ) :: tke_mf     (nx,ny,nz)  ! MF vertical velocity on edges, m/s
  real, intent(in   ) :: zpbl       (nx,ny)     ! PBLH diagnosed in TurbGridComp
  real, intent(inout) :: tke_inv    (nx,ny,nzm) ! turbulent kinetic energy. m**2/s**2
  real, intent(inout) :: tkh_inv    (nx,ny,nzm) ! eddy scalar diffusivity
  real, intent(  out) :: tkm_inv    (nx,ny,nzm) ! eddy momentum diffusivity
  real, intent(  out) :: isotropy_inv(nx,ny,nzm) ! return to isotropy timescale

  real, dimension(:,:,:), pointer :: tkesbdiss_inv  ! dissipation
  real, dimension(:,:,:), pointer :: tkesbbuoy_inv  ! buoyancy production
  real, dimension(:,:,:), pointer :: tkesbshear_inv ! shear production

  real, dimension(:,:,:), pointer :: smixt_inv    ! dissipation length scale
  real, dimension(:,:),   pointer :: lmix_out     ! mixed layer depth
  real, dimension(:,:,:), pointer :: smixt1_inv   ! length scale, term 1
  real, dimension(:,:,:), pointer :: smixt2_inv   ! length scale, term 2
  real, dimension(:,:,:), pointer :: smixt3_inv   ! length scale, term 3
  real, dimension(:,:,:), pointer :: bruntmst_inv ! moist Brunt vaisala frequency
!  real, dimension(:,:,:), pointer :: bruntdry_inv ! Brunt vaisala frequency on edges
!  real, dimension(:,:,:), pointer :: bruntedg_inv ! Brunt vaisala frequency on edges
  real, dimension(:,:,:), pointer :: ri_inv
  real, dimension(:,:,:), pointer :: prnum_inv

! SHOC tunable parameters
  real :: lambda
  real, parameter :: min_tke = 1e-4  ! Minumum TKE value, m**2/s**2
  real, parameter :: max_tke = 10.    ! Maximum TKE value, m**2/s**2

! Maximum turbulent eddy length scale, m
  real, parameter :: max_eddy_length_scale  = 2000.

! Maximum "return-to-isotropy" time scale, s
  real, parameter :: max_eddy_dissipation_time_scale = 2000.

! Constants for the TKE dissipation term based on Deardorff (1980)
  real, parameter :: pt19=0.19, pt51=0.51
  real, parameter :: Cs  = 0.15
  real :: Ck, Ce, Ces
  real :: tscale

! Number iterations for TKE solution
  integer, parameter :: nitr=6

! Local variables

  real zl      (nx,ny,nzm)  ! height of the pressure levels above surface, m
  real zi      (nx,ny,nz)   ! height of the interface levels, m
  real adzl    (nx,ny,nzm)  ! layer thickness i.e. zi(k+1)-zi(k) - defined at levels
  real adzi    (nx,ny,nz)   ! level thickness i.e. zl(k)-zl(k-1) - defined at interface
  real ri      (nx,ny,nz)   ! Local Richardson number
  real hl      (nx,ny,nzm)  ! liquid/ice water static energy , K
  real qv      (nx,ny,nzm)  ! water vapor, kg/kg
  real qcl     (nx,ny,nzm)  ! liquid water  (condensate), kg/kg
  real qci     (nx,ny,nzm)  ! ice water  (condensate), kg/kg
  real w       (nx,ny,nzm)  ! z-wind, m/s
  real bet     (nx,ny,nzm)  ! ggr/tv0
  real gamaz   (nx,ny,nzm)  ! ggr/cp*z
  real prsl    (nx,ny,nzm)  ! pressure, Pa
  real u       (nx,ny,nzm)  ! zonal velocity, m/s
  real v       (nx,ny,nzm)  ! meridional velocity, m/s
  real omega   (nx,ny,nzm)  ! pressure velocity, Pa/s
  real tabs    (nx,ny,nzm)  ! absolute temperature, K
  real qwv     (nx,ny,nzm)  ! specific humidity, kg/kg
  real qpl     (nx,ny,nzm)
  real qpi     (nx,ny,nzm)
  real cld_sgs (nx,ny,nzm)  ! cloud fraction
  real tke     (nx,ny,nzm)  ! turbulent kinetic energy, m2/s2
  real tkh     (nx,ny,nzm)
  real prnum   (nx,ny,nz)   ! Prandtl number
  real wthv_sec(nx,ny,nzm)  ! Total buoyancy flux
  real wthv_mf(nx,ny,nzm)   ! Buoyancy flux diagnosed from MF
  real tkesbdiss(nx,ny,nzm) ! TKE tendency from dissipation
  real tkesbbuoy(nx,ny,nzm)  ! TKE tendency from buoyancy
  real tkesbshear(nx,ny,nzm) ! TKE tendency from shear

! Eddy length formulation
  real smixt      (nx,ny,nzm)  ! Turbulent length scale, m
  real smixt1     (nx,ny,nzm)  ! Turbulent length scale, m
  real smixt2     (nx,ny,nzm)  ! Turbulent length scale, m
  real smixt3     (nx,ny,nzm)  ! Turbulent length scale, m
  real isotropy    (nx,ny,nzm) ! "Return-to-isotropy" eddy dissipation time scale, s
  real brunt       (nx,ny,nzm) ! Moist Brunt-Vaisalla frequency, s^-1

  real, dimension(nx,ny,nzm) :: total_water, brunt2, def2, thv, l_par
  real, dimension(nx,ny,nz)  :: brunt_edge !, brunt_dry

  real, dimension(nx,ny)     :: l_inf, l_mix, zcb, lts!, l_par!, denom, numer, cldarr

  real lstarn,    depth,    omn,      betdz, betdze,    bbb,        &
       term,      qsatt,    dqsat,    thedz,    conv_var,   &
       tkes,      pval,     pkap,     thlsec,   qwsec,      &
       qwthlsec,  wqwsec,   wthlsec,  dum,      sm,         &
       prespot,   wrk,      wrk1,     wrk2,     wrk3,       &
       tkeavg,    dtqw,     dtqi !,     l_par

  integer i,j,k,km1,ku,kd,ka,kb,kinv,strt,fnsh,cnvl
!  integer, dimension(nx,ny) :: cldbasek

  real, parameter :: bruntmin = 1e-7
  real, parameter :: vonk = 0.4

! set parameter values
! real, parameter :: Ce  = Ck**3/(0.7*Cs**4)
! real, parameter :: Ces = Ce/0.7*3.0
  lambda   = shocparams%LAMBDA              ! used in return-to-isotropy timescale
  Ck       = shocparams%CKVAL               ! Coeff in the eddy diffusivity - TKE relationship, see Eq. 7 in BK13
  Ce       = shocparams%CEFAC*Ck**3/Cs**4   ! diss ~= Ce * sqrt(tke)
  Ces      = shocparams%CESFAC*Ce           ! Ce surface factor
  tscale   = shocparams%TSCALE              ! time scale set based off of similarity results of BK13, s

! Map GEOS variables to those of SHOC
  do k=1,nz
    do j=1,ny
      do i=1,nx
        kinv = nz-k+1
        zi(i,j,kinv)    = phii_inv(i,j,k)-phii_inv(i,j,nz)
      enddo
    enddo
  enddo

  do k=1,nzm
    do j=1,ny
      do i=1,nx
        kinv = nzm-k+1
        zl(i,j,kinv)       = phil_inv(i,j,k)-phii_inv(i,j,nz)
        tkh(i,j,kinv)      = tkh_inv(i,j,k)
        prsl(i,j,kinv)     = prsl_inv(i,j,k)
        u(i,j,kinv)        = u_inv(i,j,k)
        v(i,j,kinv)        = v_inv(i,j,k)
        omega(i,j,kinv)    = omega_inv(i,j,k)
        tabs(i,j,kinv)     = tabs_inv(i,j,k)
        qwv(i,j,kinv)      = qwv_inv(i,j,k)
        qcl(i,j,kinv)      = qc_inv(i,j,k)
        qci(i,j,kinv)      = qi_inv(i,j,k)
        cld_sgs(i,j,kinv)  = cld_sgs_inv(i,j,k)
        tke(i,j,kinv)      = tke_inv(i,j,k)
        wthv_sec(i,j,kinv) = wthv_sec_inv(i,j,k)
        wthv_mf(i,j,kinv)  = wthv_mf_inv(i,j,k)
      enddo
    enddo
  enddo

  do k=1,nzm
    do j=1,ny
      do i=1,nx
        wrk            = 1.0 / prsl(i,j,k)
        qv(i,j,k)      = max(qwv(i,j,k), 0.0)
        thv(i,j,k)     = tabs(i,j,k) * (1.0+epsv*qv(i,j,k)-qcl(i,j,k)-qci(i,j,k))
        w(i,j,k)       = - rog * omega(i,j,k) * thv(i,j,k) * wrk
        qpl(i,j,k)     = 0.0  ! comment or remove when using with prognostic rain/snow
        qpi(i,j,k)     = 0.0  ! comment or remove when using with prognostic rain/snow
        total_water(i,j,k) = qcl(i,j,k) + qci(i,j,k) + qv(i,j,k)
        prespot        = (100000.0*wrk) ** kapa        ! Exner function
        bet(i,j,k)     = ggr/(tabs(i,j,k)*prespot)     ! Moorthi
        thv(i,j,k)     = thv(i,j,k)*prespot            ! Moorthi
!      
! Lapse rate * height = reference temperature
        gamaz(i,j,k) = gocp * zl(i,j,k)

! Liquid/ice water static energy - ! Note the the units are degrees K
        hl(i,j,k) = tabs(i,j,k) + gamaz(i,j,k) - fac_cond*(qcl(i,j,k)+qpl(i,j,k)) &
                                               - fac_fus *(qci(i,j,k)+qpi(i,j,k))
      enddo
    enddo
  enddo


! Define vertical grid increments for later use in the vertical differentiation

  do k=2,nzm
    km1 = k - 1
    do j=1,ny
      do i=1,nx
        adzi(i,j,k)   = (zl(i,j,k) - zl(i,j,km1))
        adzl(i,j,km1) = (zi(i,j,k) - zi(i,j,km1))   ! from 1 to nzm-1
      enddo
    enddo
  enddo
  do j=1,ny
    do i=1,nx
      adzi(i,j,1)   = (zl(i,j,1)-zi(i,j,1)) !  unused in the code
      adzi(i,j,nz)  = zi(i,j,nz)-zl(i,j,nzm) ! at the top - probably unused
      adzl(i,j,nzm) = adzi(i,j,nzm)
    enddo
  enddo

  call tke_shoc()        ! Integrate prognostic TKE equation forward in time

!=== Assign exports and flip vertical ===!

  tkh_inv(:,:,1:nzm)       = tkh(:,:,nzm:1:-1)
  tkm_inv(:,:,1:nzm)       = min(tkhmax,tkh(:,:,nzm:1:-1)*prnum(:,:,nzm:1:-1))
  isotropy_inv(:,:,1:nzm)  = isotropy(:,:,nzm:1:-1)
  tke_inv(:,:,1:nzm)       = tke(:,:,nzm:1:-1)

  ! Below exports are optional

  if (associated(tkesbdiss_inv))  tkesbdiss_inv(:,:,1:nzm)  = tkesbdiss(:,:,nzm:1:-1)
  if (associated(tkesbbuoy_inv))  tkesbbuoy_inv(:,:,1:nzm)  = tkesbbuoy(:,:,nzm:1:-1)
  if (associated(tkesbshear_inv)) tkesbshear_inv(:,:,1:nzm) = tkesbshear(:,:,nzm:1:-1)

  if (associated(lmix_out))     lmix_out(:,:)           = l_mix
  if (associated(smixt_inv))    smixt_inv(:,:,1:nzm)    = smixt(:,:,nzm:1:-1)
  if (associated(smixt1_inv))   smixt1_inv(:,:,1:nzm)   = smixt1(:,:,nzm:1:-1)
  if (associated(smixt2_inv))   smixt2_inv(:,:,1:nzm)   = smixt2(:,:,nzm:1:-1)
  if (associated(smixt3_inv))   smixt3_inv(:,:,1:nzm)   = smixt3(:,:,nzm:1:-1)

!  if (associated(bruntedg_inv)) bruntedg_inv(:,:,0:nzm) = brunt_edge(:,:,nz:1:-1)
!  if (associated(bruntdry_inv)) bruntdry_inv(:,:,0:nzm) = brunt_dry(:,:,nz:1:-1)
  if (associated(bruntmst_inv)) bruntmst_inv(:,:,1:nzm) = brunt(:,:,nzm:1:-1)
  if (associated(prnum_inv))    prnum_inv(:,:,0:nz-1)     = prnum(:,:,nz:1:-1)
  if (associated(ri_inv))       ri_inv(:,:,0:nz-1)        = ri(:,:,nz:1:-1)

!========================================!


contains

  subroutine tke_shoc()

! This subroutine solves the TKE equation,
! Heavily based on SAM's tke_full.f90 by Marat Khairoutdinov

    real grd,betdz,betdze,Cek,Cee,lstarn, bbb, omn, omp,qsatt,dqsat, smix,         &
         buoy_sgs,a_prod_sh,a_prod_bu,a_diss,a_prod_bu_debug, buoy_sgs_debug, &
         wrk, wrk1, wtke, wtk2, rdtn, tke_env
    integer i,j,k,ku,kd,itr
    real, dimension(nx,ny,nzm) :: tscale1

    rdtn = 1.0 / dtn

    call tke_shear_prod(def2)   ! Calculate shear production of TKE

    call calc_numbers()     ! returns RI and PRNUM

    do k=1,nzm
      do j=1,ny
        do i=1,nx
          tke(i,j,k)        = max(min_tke,tke(i,j,k))
          tkesbdiss(i,j,k)  = 0.
          tkesbshear(i,j,k) = 0.
          tkesbbuoy(i,j,k)  = 0.
        enddo
      enddo
    enddo

    call eddy_length()   ! Find turbulent mixing length, brunt-vaisala freq

    do k=1,nzm
      ku = k+1
      kd = k

      Cek = Ce/0.7

      if(k == 1) then
        ku = 2
        kd = 2
        Cek = Ces
      elseif(k == nzm) then
        ku = k
        kd = k
        Cek = Ces
      end if

      do j=1,ny
        do i=1,nx
          grd = adzl(i,j,k)             !  adzl(k) = zi(k+1)-zi(k)

! TKE boyancy production term. wthv_sec (buoyancy flux) is calculated in Moist GridComp.

          if (shocparams%BUOYOPT==2) then
            a_prod_bu = (ggr / thv(i,j,k)) * wthv_sec(i,j,k)
          else
            wrk  = 0.5 * (tkh(i,j,ku)*brunt_edge(i,j,ku)+tkh(i,j,kd)*brunt_edge(i,j,kd))
            a_prod_bu = -1.*wrk + (ggr / thv(i,j,k))*wthv_mf(i,j,k)
          end if

          buoy_sgs = 0.5*(brunt_edge(i,j,ku)+brunt_edge(i,j,kd))
!          buoy_sgs = - a_prod_bu / (wrk + 0.0001)   ! tkh is eddy thermal diffussivity

!Compute $c_k$ (variable Cee) for the TKE dissipation term following Eq. 11 in Deardorff (1980)
          if (buoy_sgs <= 0.0) then
            smix = grd
          else
            smix = min(grd,max(0.1*grd, 0.76*sqrt(tke(i,j,k)/(buoy_sgs+1.e-10))))
          end if

          Cee = Cek* (pt19 + pt51*smix/grd)

          wrk   = 0.25 * (tkh(i,j,ku)+tkh(i,j,kd)) * (prnum(i,j,ku) + prnum(i,j,kd))
          if (nx.eq.1) then
            a_prod_sh = min(min(tkhmax,wrk)*def2(i,j,k),0.0001)    ! TKE shear production term
          else
            a_prod_sh = min(min(tkhmax,wrk)*def2(i,j,k),0.01)    ! TKE shear production term
          end if

! Semi-implicitly integrate TKE equation forward in time
          wtke = tke(i,j,k)
          wtk2 = wtke
          wrk  = (dtn*Cee)/smixt(i,j,k)
          wrk1 = wtke + dtn*(a_prod_sh+a_prod_bu)

          wrk2 = min_tke+0.5*(tke_mf(i,j,nz-k+1)+tke_mf(i,j,nz-k))
          do itr=1,nitr                        ! iterate for implicit solution
            wtke   = min(max(wrk2, wtke), max_tke)
            a_diss = wrk*sqrt(wtke)            ! Coefficient in the TKE dissipation term
!            if (a_diss.ne.-1.) then
              wtke   = wrk1 / (1.+a_diss)
!            else
!              wtke   = wrk1 / (1.01+a_diss)
!            end if
            wtke   = tkef1*wtke + tkef2*wtk2   ! tkef1+tkef2 = 1.0
            wtk2   = wtke

          enddo

          tke(i,j,k) = min(max(wrk2, wtke), max_tke)

          tscale1(i,j,k) = (dtn+dtn) / a_diss        ! See Eq 8 in BK13 (note typo, flipped num/denom)

          a_diss     = (a_diss/dtn)*tke(i,j,k)    ! TKE dissipation term, epsilon

! TKE budget terms
          tkesbdiss(i,j,k)       = -a_diss
          tkesbshear(i,j,k)      = a_prod_sh
          tkesbbuoy(i,j,k)       = a_prod_bu

        end do ! i loop
      end do   ! j loop
    end do     ! k


    ! Below averages TKE from adjacent levels and subtracts TKE diagnosed from EDMF updrafts
    ! to estimate an environmental TKE on edge. Isotropy from adjacent levels is similarly
    ! averaged, and product of edge isotropy and environmental TKE provides diffusivity.
!    wrk = ck
    ! tkh defined 1:nzm (corresponding to edges 1:nzm), isotropy defined 1:nz, brunt_edge defined 1:nz, tscale1 defined 1:nzm 
    do k=2,nzm
      do j=1,ny
        do i=1,nx
          ! Calculate "return-to-isotropy" eddy dissipation time scale, see Eq. 8 in BK13
          if (brunt_edge(i,j,k) <= 1e-5 .or. zl(i,j,k).lt.0.5*zpbl(i,j)) then
            isotropy(i,j,k) = max(30.,min(max_eddy_dissipation_time_scale,0.5*(tscale1(i,j,k)+tscale1(i,j,k-1))))
          else
            wrk = 0.5*(tscale1(i,j,k)+tscale1(i,j,k-1))
            isotropy(i,j,k) = max(30.,min(max_eddy_dissipation_time_scale,wrk/(1.0+lambda*brunt_edge(i,j,k)*wrk*wrk)))
!            isotropy(i,j,k) = max(30.,min(max_eddy_dissipation_time_scale,wrk/(1.0+lambda*0.5*(brunt(i,j,k)+brunt(i,j,k-1))*wrk*wrk)))
          endif
          if (tke(i,j,k).lt.2e-4) isotropy(i,j,k) = 30.

          wrk1 = ck / prnum(i,j,k)

!          tkh(i,j,k) = 0.5*( smixt(i,j,k)*sqrt(tke(i,j,k)) &      ! alternate form
!                           + smixt(i,j,k-1)*sqrt(tke(i,j,k-1)) )

!          tke_env = max(min_tke,0.5*(tke(i,j,k)+tke(i,j,k-1))-0.*tke_mf(i,j,nz-k+1))
          tkh(i,j,k) = wrk1*isotropy(i,j,k)*0.5*(tke(i,j,k)+tke(i,j,k-1)) !   &
!                       *(tke_env)             ! remove MF TKE
          tkh(i,j,k) = min(tkh(i,j,k),tkhmax)
        end do ! i
      end do ! j
    end do ! k
    isotropy(:,:,1) = isotropy(:,:,2)

! add radiation-driven entrainment (simplified)
!  do j=1,ny
!    do i=1,nx
!      do k=2,nzm-1
!        if (zl(i,j,k).gt.4000.) exit
!        if (cld_sgs(i,j,k).gt.0.1 .and. cld_sgs(i,j,k+1).lt.0.1) then

!          dbuoy = max( 0.1, (thv(i,j,k)/ggr)*brunt(i,j,k) )
!          krad = -1.*adzl(i,j,k)*shocparams%kradfac*minval(dtdtrad(i,j,nzm-k-1:nzm-k+1))/dbuoy
!          tkh(i,j,k+1) = tkh(i,j,k+1) + krad  !add to KH at interface between cloud>0.1 and <0.1
!          print *,'krad=',krad,'  ztop=',zi(i,j,k+1)

!          exit
!        end if
!      end do
!    end do
!  end do

  end subroutine tke_shoc

  subroutine calc_numbers()
     ! Defines Richardson number and Prandtl number on edges
    real, dimension(nx,ny,nzm-1) :: DU

     DU  = (U(:,:,1:nzm-1) - U(:,:,2:nzm))**2 + &    ! shear on edges
           (V(:,:,1:nzm-1) - V(:,:,2:nzm))**2
     DU  = MIN( MAX( SQRT(DU) / adzi(:,:,1:nzm-1), 0.005 ), 0.025 )

     RI = 0.0
     RI(:,:,2:nz-1) = ggr*( (THV(:,:,2:nzm) - THV(:,:,1:nzm-1)) / adzi(:,:,1:nzm-1) ) &
                      / ( 0.5*( THV(:,:,1:nzm-1)+THV(:,:,2:nzm) ) * (DU**2) )

     if (SHOCPARAMS%PRNUM.lt.0.) then
!        where (RI.le.0. .or. tke_mf(:,:,nz:1:-1).gt.1e-6)
        where (RI.le.0. .or. tke_mf(:,:,nz:1:-1).gt.1e-4)
          PRNUM = -1.*SHOCPARAMS%PRNUM
        elsewhere
          ! He et al 2019
!             tmp3de = RI*(1.+6.*RI)
!          PRNUM = (0.9+4.*tmp3de*SQRT(1.-SHOCPARAMS%PRNUM*8.*tmp3de/3.))/(1.+4.*tmp3de)
        ! Han and Bretherton 2019
          PRNUM = -1.*SHOCPARAMS%PRNUM+2.1*MIN(10.,RI) ! limit RI to avoid instability
        end where
     else
        PRNUM = SHOCPARAMS%PRNUM
     end if

  end subroutine calc_numbers


  subroutine tke_shear_prod(def2)

! Calculate TKE shear production term

    real, intent(out):: def2(nx,ny,nzm)

    real    rdzw_up, rdzw_dn, wrku(2), wrkv(2), wrkw(2)
    real    txd(nx,ny)
    integer i,j,k,kb,kc

    do k=1,nzm
      do j=1,ny
        do i=1,nx
          def2(i,j,k) = 0.0
        enddo
      enddo
    enddo

    do k=1,nzm

      kb = k-1
      kc = k+1

      if (k == 1) then

        do j=1,ny
          do i=1,nx
            rdzw_up     = 1./adzi(i,j,kc)
            wrku(1)     = (u(i,j,kc)-u(i,j,k))*rdzw_up
            wrkv(1)     = (v(i,j,kc)-v(i,j,k))*rdzw_up
!           wrkw(1)     = (w(i,j,kc)-w(i,j,k))*rdzw_up
            def2(i,j,1) = wrku(1)*wrku(1) + wrkv(1)*wrkv(1) !+ 2*wrkw(1) * wrkw(1)
            txd(i,j)    = rdzw_up
          enddo
        enddo

      elseif (k < nzm ) then
        do j=1,ny
          do i=1,nx
            rdzw_up     = 1./adzi(i,j,kc)
            rdzw_dn     = txd(i,j)
            wrku(1)     = (u(i,j,kc)-u(i,j,k))*rdzw_up
            wrku(2)     = (u(i,j,k)-u(i,j,kb))*rdzw_dn
            wrkv(1)     = (v(i,j,kc)-v(i,j,k))*rdzw_up
            wrkv(2)     = (v(i,j,k)-v(i,j,kb))*rdzw_dn
!           wrkw(1)     = (w(i,j,kc)-w(i,j,k))*rdzw_up
!           wrkw(2)     = (w(i,j,k)-w(i,j,kb))*rdzw_dn

            def2(i,j,k) = 0.5 * (wrku(1)*wrku(1) + wrku(2)*wrku(2)     &
                               + wrkv(1)*wrkv(1) + wrkv(2)*wrkv(2))  ! &
!                              + wrkw(1)*wrkw(1) + wrkw(2)*wrkw(2)
            txd(i,j)    = rdzw_up
          enddo
        enddo
      else
        do j=1,ny
          do i=1,nx
            rdzw_dn     = txd(i,j)
            wrku(2)     = (u(i,j,k)-u(i,j,kb))*rdzw_dn
            wrkv(2)     = (v(i,j,k)-v(i,j,kb))*rdzw_dn
!           wrkw(2)     = (w(i,j,k)-w(i,j,kb))*rdzw_dn
            def2(i,j,k) = wrku(2)*wrku(2) + wrkv(2)*wrkv(2)  !+ 2*wrkw(2) * wrkw(2)
          enddo
        enddo
      endif

    end do     ! k  loop


  end subroutine tke_shear_prod

  subroutine eddy_length()

! This subroutine computes the turbulent length scale

! Local variables
    real    wrk, wrk1, wrk2, wrk3!, zdecay
    integer i, j, k, kk, kl, ku, kb, kc!, ktop

!    do j=1,ny
!      do i=1,nx
!        cldarr(i,j) = 0.0
!        numer(i,j)  = 0.0
!        denom(i,j)  = 0.0
!      enddo
!    enddo


    do k=1,nzm
      kb = k-1
      kc = k+1

      do j=1,ny
        do i=1,nx

          if (k == 1) then
            kb = 1
            kc = 2
            thedz = adzi(i,j,kc)
          elseif (k == nzm) then
            kb = nzm-1
            kc = nzm
            thedz = adzi(i,j,k)
          else
            thedz = (adzi(i,j,kc)+adzi(i,j,k)) !  = (z(k+1)-z(k-1))
          endif
          betdz = bet(i,j,k) / thedz


!          brunt_edge(i,j,k) = (2.*ggr/(thv(i,j,k)+thv(i,j,kb)))*(thv(i,j,k)-thv(i,j,kb))/adzi(i,j,k)  !g/thv/dz *(thv-thv)
!          brunt_dry(i,j,k) = (thv(i,j,kc)-thv(i,j,kb))*betdz

! Reinitialize the mixing length related arrays to zero
          smixt(i,j,k)    = 1.0   ! shoc_mod module variable smixt
          brunt(i,j,k)    = 0.0

!Eq. 11 in BK13 (Eq. 4.13 in Pete's dissertation)
!Outside of cloud, integrate from the surface to the cloud base
!Should the 'if' below check if the cloud liquid < a small constant instead?

!          if (qcl(i,j,k)+qci(i,j,k) <= 1e-6 .and. cldarr(i,j).eq.0.0) then
!          if (qcl(i,j,k)+qci(i,j,k) <= 0.) then
!            tkes       = sqrt(tke(i,j,k)) * adzl(i,j,k)
!            numer(i,j) = numer(i,j) + tkes*zl(i,j,k) ! Numerator in Eq. 11 in BK13
!            denom(i,j) = denom(i,j) + tkes           ! Denominator in Eq. 11 in BK13

!          else
!            cldarr(i,j) = 1.0   ! Take note of columns containing cloud.
!          endif
        enddo
      enddo
    enddo
!    brunt_edge(:,:,nz) = brunt_edge(:,:,nzm)

! Calculate the measure of PBL depth,  Eq. 11 in BK13 (Is this really PBL depth?)
!    cldbasek(:,:) = 1
!    do j=1,ny
!      do i=1,nx

!         do k=1,nzm
!           if (zl(i,j,k).gt.3000. .or. cld_sgs(i,j,k).gt.0.01) exit
!           tkes       = sqrt(tke(i,j,k)) * adzl(i,j,k)
!           numer(i,j) = numer(i,j) + tkes*zl(i,j,k) ! Numerator in Eq. 11 in BK13
!           denom(i,j) = denom(i,j) + tkes           ! Denominator in Eq. 11 in BK13
!         end do

!         if (denom(i,j) >  0.0 .and. numer(i,j) > 0.0) then
!!           l_inf(i,j) = max(min(0.1 * (numer(i,j)/denom(i,j)),300.),10.)
!           l_inf(i,j) = max(min( numer(i,j)/denom(i,j), 1000. ),100.)
!         else
!           l_inf(i,j) = 100.
!         endif

        ! Identify mixed layer top as level where THV exceeds THV(3) + 0.4 K
        ! Interpolate for final height based on gradient
        ! Ignore single isolated levels
!        kk = 4
!        do while (thv(i,j,3)+0.4 .gt. thv(i,j,kk) .or. thv(i,j,3)+0.4 .gt. thv(i,j,kk+1))
!          kk = kk+1
!        end do
!        dum = (thv(i,j,kk-1)-thv(i,j,kk-2))

!        if (abs(dum) .gt. 1e-3) then
!          l_mix(i,j) = max(zl(i,j,kk-1)+0.*(thv(i,j,3)+0.4-thv(i,j,kk-1))*(zl(i,j,kk-1)-zl(i,j,kk-2))/dum,100.)
!        else
!          l_mix(i,j) = max(zl(i,j,kk-1),100.)
!        end if

!        do while ((zl(i,j,cldbasek(i,j)).lt.300.) .or. (cld_sgs(i,j,cldbasek(i,j)).lt.0.001 .and. cldbasek(i,j).lt.nzm))
!          cldbasek(i,j) = cldbasek(i,j) + 1
!        end do

!        kk = 1
!        do while (zl(i,j,kk) .lt. 3000. .or. kk.eq.nzm)
!          kk = kk + 1
!        end do
!        lts(i,j) = thv(i,j,kk) - thv(i,j,1)

!       Alternate cloud base calculation
!        tep = tabs(i,j,1)
!        qsp = MAPL_EQsat(tabs(i,j,1),prsl(i,j,1),dtqw)
!        kk = 1
!        do while (qsp .gt. total_water(i,j,1) .and. zl(i,j,kk).lt.1500.)
!          kk = kk+1
!          tep = tep - ggr*( zl(i,j,kk)-zl(i,j,kk-1) )/cp
!          qsp = MAPL_EQsat(tep,prsl(i,j,kk),dtqw)
!        end do
!        zcb(i,j) = max(200.,zl(i,j,kk-1))  !kk-1 is highest level *before* condensation
!        if (nx.eq.1) print *,'zcb=',zcb(i,j)
!      enddo
!    enddo


!Calculate length scale outside of cloud, Eq. 10 in BK13 (Eq. 4.12 in Pete's dissertation)
    do k=1,nzm

      kb = k-1
      kc = k+1

      do j=1,ny
        do i=1,nx

!  vars module variable bet (=ggr/tv0) ; grid module variable  adzi

          if (k == 1) then
            kb = 1
            kc = 2
            thedz = adzi(i,j,kc)
          elseif (k == nzm) then
            kb = nzm-1
            kc = nzm
            thedz = adzi(i,j,k)
            betdze = 0.5*(bet(i,j,k)-bet(i,j,kb)) / adzi(i,j,k)
          else
            thedz = (adzi(i,j,kc)+adzi(i,j,k)) !  = (z(k+1)-z(k-1))
            betdze = 0.5*(bet(i,j,k)-bet(i,j,kb)) / adzi(i,j,k)
          endif
          betdz = bet(i,j,k) / thedz


! Compute local Brunt-Vaisalla frequency

          wrk = qcl(i,j,k) + qci(i,j,k)

! Find the in-cloud Brunt-Vaisalla frequency

! ideally we should use fQi or ice_fraction() from MoistGC here
             omn = qcl(i,j,k) / (wrk+1.e-20) ! Ratio of liquid water to total water

! Latent heat of phase transformation based on relative water phase content
! fac_cond = lcond/cp, fac_fus = lfus/cp

             lstarn = fac_cond + (1.-omn)*fac_fus

! Saturation mixing ratio over water/ice wrt temp  based on relative water phase content
!             qsatt =     omn  * qsatw(tabs(i,j,k),prsl(i,j,k))                &
!                   + (1.-omn) * qsati(tabs(i,j,k),prsl(i,j,k))
             qsatt =     omn  * MAPL_EQsat(tabs(i,j,k),prsl(i,j,k),dtqw)     &
                   + (1.-omn) * MAPL_EQsat(tabs(i,j,k),prsl(i,j,k),dtqi,OverIce=.TRUE.)

! Derivative of saturation mixing ratio over water/ice wrt temp. based on relative water phase content
!             dqsat =     omn  * dtqsatw(tabs(i,j,k),prsl(i,j,k))             &
!                   + (1.-omn) * dtqsati(tabs(i,j,k),prsl(i,j,k))
             dqsat =  omn * dtqw + (1.-omn) * dtqi

! liquid/ice moist static energy static energy divided by cp?
             bbb = (1. + epsv*qsatt-wrk-qpl(i,j,k)-qpi(i,j,k)                &
                 + 1.61*tabs(i,j,k)*dqsat) / (1.+lstarn*dqsat)

! Calculate Brunt-Vaisalla frequency using centered differences in the vertical

             brunt(i,j,k) = cld_sgs(i,j,k)*betdz*(bbb*(hl(i,j,kc)-hl(i,j,kb))               &
                          + (bbb*lstarn - (1.+lstarn*dqsat)*tabs(i,j,k))     &
                          * (total_water(i,j,kc)-total_water(i,j,kb))        &
                          + (bbb*fac_cond - (1.+fac_cond*dqsat)*tabs(i,j,k))*(qpl(i,j,kc)-qpl(i,j,kb))  &
                          + (bbb*fac_sub  - (1.+fac_sub*dqsat)*tabs(i,j,k))*(qpi(i,j,kc)-qpi(i,j,kb)) )

             bbb = 0.5*(bbb + (1. + epsv*qsatt-wrk-qpl(i,j,k-1)-qpi(i,j,k-1)                &
                 + 1.61*tabs(i,j,k-1)*dqsat) / (1.+lstarn*dqsat) )
             if (k.gt.1) then
             brunt_edge(i,j,k) = 0.5*(cld_sgs(i,j,k)+cld_sgs(i,j,k-1))*betdz*(bbb*(hl(i,j,k)-hl(i,j,k-1))               &
                          + (bbb*lstarn - (1.+lstarn*dqsat)*tabs(i,j,k))     &
                          * (total_water(i,j,k)-total_water(i,j,k-1))        &
                          + (bbb*fac_cond - (1.+fac_cond*dqsat)*tabs(i,j,k))*(qpl(i,j,k)-qpl(i,j,k-1))  &
                          + (bbb*fac_sub  - (1.+fac_sub*dqsat)*tabs(i,j,k))*(qpi(i,j,k)-qpi(i,j,k-1)) )
             end if

! Find outside-of-cloud Brunt-Vaisalla frequency
! Only unsaturated air, rain and snow contribute to virt. pot. temp.
! liquid/ice moist static energy divided by cp?

             bbb = 1. + epsv*qv(i,j,k) - qpl(i,j,k) - qpi(i,j,k)
             brunt(i,j,k) = brunt(i,j,k) + (1.-cld_sgs(i,j,k))*betdz*( bbb*(hl(i,j,kc)-hl(i,j,kb))                        &
                          + epsv*tabs(i,j,k)*(total_water(i,j,kc)-total_water(i,j,kb)) &
                          + (bbb*fac_cond-tabs(i,j,k))*(qpl(i,j,kc)-qpl(i,j,kb))       &
                          + (bbb*fac_sub -tabs(i,j,k))*(qpi(i,j,kc)-qpi(i,j,kb)) )

             bbb = 0.5*(bbb + 1. + epsv*qv(i,j,k-1) - qpl(i,j,k-1) - qpi(i,j,k-1))
             if (k.gt.1) then
             brunt_edge(i,j,k) = brunt_edge(i,j,k) + (1.-0.5*(cld_sgs(i,j,k)+cld_sgs(i,j,k-1)))*betdz*( bbb*(hl(i,j,k)-hl(i,j,k-1))                        &
                          + epsv*tabs(i,j,k)*(total_water(i,j,k)-total_water(i,j,k-1)) &
                          + (bbb*fac_cond-tabs(i,j,k))*(qpl(i,j,k)-qpl(i,j,k-1))       &
                          + (bbb*fac_sub -tabs(i,j,k))*(qpi(i,j,k)-qpi(i,j,k-1)) )
             end if
          end do

        end do
      end do

! Reduction of mixing length in the stable regions (where B.-V. freq. > 0) is required.
! Here we find regions of Brunt-Vaisalla freq. > 0 for later use.

      brunt_edge(:,:,1) = brunt_edge(:,:,2)
      brunt_edge(:,:,nz) = brunt_edge(:,:,nzm)
      do i=1,nx
        do j=1,ny
          do k=1,nzm
            brunt2(i,j,k) = (1.-cld_sgs(i,j,k))*0.5*(brunt_edge(i,j,k-1)+brunt_edge(i,j,k)) + cld_sgs(i,j,k)*min(brunt_edge(i,j,k-1),brunt_edge(i,j,k))
            if (brunt2(i,j,k) < 1e-5 .or. zl(i,j,k).lt.0.5*zpbl(i,j)) then
              brunt2(i,j,k) = 1e-10
            endif
          end do
        end do
      end do
      
      brunt2(:,:,1) = brunt2(:,:,2)
      brunt2(:,:,nzm) = brunt2(:,:,nzm-1)

!      do j=1,ny
!        do i=1,nx
!          kk = 1
!          l_mix(i,j) = 0.
!          do while (l_mix(i,j).lt.1e-5 .and. kk.lt.nzm)
!            l_mix(i,j) = l_mix(i,j) + brunt(i,j,kk)*adzl(i,j,kk)
!            kk = kk+1
!          end do
!          l_mix(i,j) = zl(i,j,kk)
!        end do
!      end do

!         brunt_dry = max( bruntmin, brunt_dry )

      do k=1,nzm
        do j=1,ny
          do i=1,nx

            tkes = sqrt(tke(i,j,k))

! Calculate turbulent length scale in the boundary layer.
! See Eq. 10 in BK13 (Eq. 4.12 in Pete's dissertation)


          !----------------------------------
          ! calculate parcel mixing length
          !----------------------------------
              kk = k
              wrk = thv(i,j,k)+0.2 !max(0.001,3.*wthv_sec(i,j,k)) / max(0.05,sqrt(0.667)*tkes) ! upward T perturbation
              do while (wrk .gt. thv(i,j,kk+1) .and. kk.lt.nzm)
                kk = kk+1
              end do
              l_par(i,j,k) = zl(i,j,kk) + max(0.,(wrk-thv(i,j,kk))* &
                             (zl(i,j,kk+1)-zl(i,j,kk)) / (thv(i,j,kk+1)-thv(i,j,kk)))
              kk = k
              wrk = thv(i,j,k)-0.2 !max(0.001,3.*wthv_sec(i,j,k)) / max(0.05,sqrt(0.667)*tkes) ! downward T perturbation
              do while (wrk .lt. thv(i,j,kk-1) .and. kk .gt. 1)
                kk = kk-1
              end do
              l_par(i,j,k) = l_par(i,j,k) - zl(i,j,kk) + max(0.,(thv(i,j,kk)-wrk)* &
                             (zl(i,j,kk)-zl(i,j,kk-1))/(thv(i,j,kk)-thv(i,j,kk-1)))
              l_par(i,j,k) = max(min(l_par(i,j,k),1500.),25.)


          !----------------------------------
          ! calculate 'TKE' mixing length
          !----------------------------------
!                 l_mix(i,j) = 0.
!                 kinv = k
!                 do while (tke(i,j,kinv).gt.0.02) 
!                    l_mix(i,j) = l_mix(i,j) + adzl(i,j,kinv)
!                    kinv = kinv+1
!                 end do
!                 kinv = k-1
!                 do while (tke(i,j,kinv).gt.0.02) 
!                    l_mix(i,j) = l_mix(i,j) + adzl(i,j,kinv) 
!                    kinv = kinv-1
!                 end do
!                 l_mix(i,j) = 0.1*l_mix(i,j)


            if ( shocparams%LENOPT .lt. 4 ) then  ! SHOC-MF length scale

                 ! Surface length scale
                 smixt1(i,j,k) = vonk*zl(i,j,k)*shocparams%LENFAC1
!                 smixt1(i,j,k) = sqrt(400.*tkes*vonk*zl(i,j,k))*shocparams%LENFAC1  ! original SHOC, includes TKE

                 ! Turbulent length scale
                 smixt2(i,j,k) = sqrt(l_par(i,j,k)*400.*tkes)*shocparams%LENFAC2

                 ! Stability length scale, reduced influence below 300m
                 smixt3(i,j,k) = max(0.05,tkes)*shocparams%LENFAC3/(sqrt(brunt2(i,j,k))) 

                 !=== Combine component length scales ===
                 if (shocparams%LENOPT .eq. 1) then  ! JPL blending approach (w/SHOC length scales)
                      wrk1 = sqrt(3./(1./smixt2(i,j,k)**2+1./smixt3(i,j,k)**2))
                      if (zl(i,j,k).lt.300.) then
                         smixt(i,j,k) = wrk1 + (smixt1(i,j,k)-wrk1)*exp(-(zl(i,j,k)/60.))
                      else
                         smixt(i,j,k) = wrk1
                      end if
                 else if (shocparams%LENOPT .eq. 2) then  ! Geometric average
                    smixt(i,j,k) = min(max_eddy_length_scale, 3./(1./smixt1(i,j,k)+1./smixt2(i,j,k)+1./smixt3(i,j,k)) )
                 else if (shocparams%LENOPT .eq. 3) then  ! SHOC classic approach
                    smixt(i,j,k) = min(max_eddy_length_scale, SQRT(3.)/SQRT(1./smixt1(i,j,k)**2+1./smixt2(i,j,k)**2+1./smixt3(i,j,k)**2) )
                 end if
           else if (shocparams%LENOPT .eq. 4) then  ! JPL Length scale (Suselj et al 2012)
              wrk2 = 1.0/(400.*tkes)
              wrk3 = sqrt(brunt2(i,j,k))/(0.7*tkes)
              wrk1 = 1.0/(wrk2+wrk3)
              smixt(i,j,k) = 3.3*shocparams%LENFAC1*(wrk1 + (vonk*zl(i,j,k)-wrk1)*exp(-zl(i,j,k)/(0.1*zpbl(i,j))))
              smixt1(i,j,k) = 3.3*shocparams%LENFAC1/wrk2
              smixt2(i,j,k) = 3.3*shocparams%LENFAC1/wrk3
              smixt3(i,j,k) = 3.3*shocparams%LENFAC1*vonk*zl(i,j,k)
           end if

           ! Enforce minimum and maximum length scales
           wrk = 20. !0.5*min(100.,adzl(i,j,k))     ! Minimum 0.1 of local dz (up to 200 m)
           if (zl(i,j,k) .lt. 2000.) then
              smixt(i,j,k) = max(wrk, smixt(i,j,k))
           else if (zl(i,j,k).gt.zpbl(i,j)) then ! if above 2 km and dry CBL top, cap length scale
              smixt(i,j,k) = max(wrk, min(200.,smixt(i,j,k)))
           end if
        end do
      end do
    end do



  end subroutine eddy_length



!  subroutine buoyancy_single_gaussian()

! Compute SGS buoyancy flux using analytic SINGLE-gaussian PDF
! for moisture and liquid water static energy.

! Local variables

!    integer i,j,k,ku,kd
!    real a,b,c,alpha,beta,cc,s,stds

!    DO k=1,nzm

!      kd = k
!      ku = k + 1
!      if (k == nzm) ku = k

!      DO j=1,ny
!        DO i=1,nx

!          pval = prsl(i,j,k)
!          pkap = (pval/100000.0) ** kapa

!          wqwsec   = 0.5 * (wqw_sec(i,j,kd)   + wqw_sec(i,j,ku))
!          wthlsec  = 0.5 * (wthl_sec(i,j,kd)  + wthl_sec(i,j,ku))

          ! following Bechtold et al 1995
!          b = tabs(i,j,k) - fac_cond*qcl(i,j,k)  ! Bechtold eqn 4
!          c = MAPL_EQsat(b,prsl(i,j,k),dtqw)
!          s = total_water(i,j,k) - MAPL_EQsat(tabs(i,j,k),prsl(i,j,k))

!          a = 1. / (1.+fac_cond*dtqw)            ! Bechtold eqn 6
!          b = a * (1./pkap) * dtqw
!          c = a*(total_water(i,j,k) - c)

!          alpha = 0.61 * (tabs(i,j,k)/pkap)      ! Bechtold eqn 12
!          beta = (1./pkap) * fac_cond - 1.61*(tabs(i,j,k)/pkap)

          ! Bechtold eqn
!          stds = alpha**2 *qwsec + beta**2 * thlsec &
!                 - 2.*alpha*beta*qwthlsec

          ! Cloud fraction, assuming single Gaussian qt distribution
!          cc = 0.5 + 0.5*erf(s/(sqrt2*stds))

!          wthv_sec(i,j,k) = wthlsec*( 1.+0.61*total_water(i,j,k) &
!                                     - beta*b*cc )               &
!                           + wqwsec*( alpha + beta*a*cc )

!        ENDDO
!      ENDDO
!    ENDDO

!  end subroutine buoyancy_single_gaussian

 end subroutine run_shoc


 subroutine update_moments( IM, JM, LM, & ! in
                             DT,       &  ! in
                             SH,       &  ! in
                             EVAP,     &  ! in
                             ZL,       &  ! in
                             ZLE,      &  ! in
                             KH,       &  ! in
                             BRUNT,    &  ! in
                             TKE,      &  ! in
                             ISOTROPY, &  ! in
                             QT,       &  ! in
                             HL,       &  ! in
                             MFFRC,    &  ! in
                             MFQT3,    &  ! in
                             MFHL3,    &  ! in
                             MFW2,     &  ! in
                             MFW3,     &  ! in
                             MFWQT,    &  ! in
                             MFWHL,    &  ! in
                             MFHLQT,   &  ! in
                             WQT_DC,   &  ! in
                             PDF_A,    &  ! inout
                             qt2,      &  ! inout
                             qt3,      &  ! inout
                             hl2,      &  ! out
                             hl3,      &  ! out
                             w2,       &  ! out
                             w3,       &  ! out
                             w3can,    &  ! out
                             wqt,      &  ! out
                             whl,      &  ! out
                             hlqt,     &  ! out
                             qt2diag,  &  ! out
                             hl2diag,  &  ! out
                             hlqtdiag, &  ! out
                           doprogqt2,  &  ! tuning parameters
                           hl2tune,    &
                           qt2tune,    &
                           hlqt2tune,  &
                           skew_tgen,  &
                           skew_tdis,  &
                           docanuto )


    integer, intent(in   ) :: IM, JM, LM       ! dimensions
    real,    intent(in   ) :: DT               ! timestep [s]
    real,    intent(in   ) :: SH   (IM,JM)     ! surface sensible heat flux
    real,    intent(in   ) :: EVAP (IM,JM)     ! surface evaporation
    real,    intent(in   ) :: ZL   (IM,JM,LM)  ! heights [m]
    real,    intent(in   ) :: ZLE  (IM,JM,0:LM)  ! edge heights [m]
    real,    intent(in   ) :: KH   (IM,JM,0:LM)  ! diffusivity
    real,    intent(in   ) :: BRUNT(IM,JM,LM)  ! Brunt-Vaisala frequency
    real,    intent(in   ) :: TKE  (IM,JM,LM)  ! turbulent kinetic energy
    real,    intent(in   ) :: ISOTROPY(IM,JM,0:LM)  ! isotropy timescale
    real,    intent(in   ) :: QT   (IM,JM,LM)  ! total water
    real,    intent(in   ) :: HL   (IM,JM,LM)  ! liquid water static energy
    real,    intent(in   ) :: MFFRC(IM,JM,LM)  ! mass flux area fraction
    real,    intent(in   ) :: MFQT3(IM,JM,LM)  !
    real,    intent(in   ) :: MFHL3(IM,JM,LM)  !
    real,    intent(in   ) :: MFW2 (IM,JM,LM)  !
    real,    intent(in   ) :: MFW3 (IM,JM,LM)  !
    real,    intent(in   ) :: MFWQT(IM,JM,0:LM)  !
    real,    intent(in   ) :: MFWHL(IM,JM,0:LM)  !
    real,    intent(in   ) :: MFHLQT(IM,JM,LM) !
    real,    intent(in   ) :: WQT_DC(IM,JM,0:LM)  !
    real,    intent(inout) :: PDF_A(IM,JM,LM)  ! first plume area fraction
    real,    intent(inout) :: qt2  (IM,JM,LM)  ! total water variance
    real,    intent(inout) :: qt3  (IM,JM,LM)  ! third moment of total water
    real,    intent(  out) :: hl2  (IM,JM,LM)  ! liquid water static energy variance
    real,    intent(  out) :: hl3  (IM,JM,LM)  ! third moment static energy
    real,    intent(  out) :: w2   (IM,JM,LM)  ! vertical velocity variance
    real,    intent(  out) :: w3   (IM,JM,LM)  ! third moment vertical velocity
    real,    intent(  out) :: w3can(IM,JM,LM)  ! third moment vertical velocity
    real,    intent(  out) :: wqt  (IM,JM,LM)  ! vertical flux of total water
    real,    intent(  out) :: whl  (IM,JM,LM)  ! vertical flux of liquid water static energy
    real,    intent(  out) :: hlqt (IM,JM,LM)  ! total water, static energy covariance
    real,    intent(  out) :: qt2diag(IM,JM,LM)
    real,    intent(  out) :: hl2diag(IM,JM,LM)
    real,    intent(  out) :: hlqtdiag(IM,JM,LM)

    real,    intent(in   ) :: HL2TUNE,     &   ! tuning parameters
                              HLQT2TUNE,   &
                              QT2TUNE,     &
                              SKEW_TGEN,   &
                              SKEW_TDIS

    integer, intent(in   ) :: DOPROGQT2,   &   ! prognostic QT2 switch
                              DOCANUTO

    real, parameter :: HL2MIN = 0.0005
    real, parameter :: HL2MAX = 2.0

    ! Local variables
    integer :: k, kd, ku
    real, dimension(IM,JM) :: wrk1, wrk2, wrk3
    real, dimension(IM,JM) :: sm, onemmf
    real, dimension(IM,JM,0:LM) :: qt2_edge, &
                                   qt2_edge_nomf, &
                                   hl2_edge, &
                                   hl2_edge_nomf, &
                                   wqt_edge, &
                                   whl_edge, &
                                   hlqt_edge,&
                                   qtgrad
    real, dimension(IM,JM,LM) :: adzl, bet, whl_can
!======= Canuto variables
    integer i, j, kb, kc, km1
    real bet2,   f0,     f1,  f2,    f3,   f4,  f5,  iso, isosqr,             &
         omega0,  omega1, omega2, X0,  Y0,    X1,   Y1,  AA0, AA1, buoy_sgs2, &
         thedz,   thedz2, cond,   wrk, avew, wrk1b, wrk2b, wrk3b, dum
! See Eq. 7 in C01 (B.7 in Pete's dissertation)
    real, parameter :: c=7.0, a0=0.52/(c*c*(c-2.)), a1=0.87/(c*c),      &
                       a2=0.5/c, a3=0.6/(c*(c-2.)), a4=2.4/(3.*c+5.),   &
                       a5=0.6/(c*(3.*c+5))
!========

    bet = 9.806/300.

    qt2diag = 0.
    hl2diag = 0.
    hlqtdiag = 0.

    ! define resolved gradients on edges
    do k=1,LM-1
        wrk1 = 1.0 / (ZL(:,:,k)-ZL(:,:,k+1))
        wrk3 = KH(:,:,k) * wrk1

        sm   = 0.5*ISOTROPY(:,:,k)*wrk1*wrk3 !Tau*Kh/dz^2

        ! SGS vertical flux liquid/ice water static energy. Eq 1 in BK13
        wrk1            = HL(:,:,k) - HL(:,:,k+1)
        whl_edge(:,:,k) = - wrk3 * wrk1

        ! SGS vertical flux of total water. Eq 2 in BK13
        wrk2            = QT(:,:,k) - QT(:,:,k+1)
        wqt_edge(:,:,k) = - wrk3 * wrk2

        ! Second moment of liquid/ice water static energy. Eq 4 in BK13
        hl2_edge_nomf(:,:,k) = HL2TUNE * sm * wrk1 * wrk1
        hl2_edge(:,:,k) = HL2TUNE * 0.5*ISOTROPY(:,:,k) * &
                          (wrk3*wrk1-MFWHL(:,:,k)) * wrk1/(ZL(:,:,k)-ZL(:,:,k+1))

        ! Second moment of total water mixing ratio.  Eq 3 in BK13
        qtgrad(:,:,k)   = wrk2 / (ZL(:,:,k)-ZL(:,:,k+1))
        qt2_edge(:,:,k) = (KH(:,:,k)*qtgrad(:,:,k)-MFWQT(:,:,k)-0.*WQT_DC(:,:,k))*qtgrad(:,:,k)  ! gradient production
        qt2_edge_nomf(:,:,k) = (KH(:,:,k)*qtgrad(:,:,k))*qtgrad(:,:,k)  ! gradient production

        ! Covariance of total water mixing ratio and liquid/ice water static energy.  Eq 5 in BK13
        hlqt_edge(:,:,k) = HLQT2TUNE * sm * wrk1 * wrk2
    end do

    ! set lower boundary conditions
    whl_edge(:,:,LM)  = SH(:,:)/cp
    wqt_edge(:,:,LM)  = EVAP(:,:)
    hl2_edge(:,:,LM)  = hl2_edge(:,:,LM-1)
    hl2_edge_nomf(:,:,LM)  = hl2_edge_nomf(:,:,LM-1)
    qt2_edge(:,:,LM)  = qt2_edge(:,:,LM-1)
    qt2_edge_nomf(:,:,LM)  = qt2_edge_nomf(:,:,LM-1)
    hlqt_edge(:,:,LM) = hlqt_edge(:,:,LM-1)
    qtgrad(:,:,LM)    = qtgrad(:,:,LM-1)
    qtgrad(:,:,0)     = qtgrad(:,:,1)




    do k=1,LM
        kd = k-1
        ku = k
        if (k==1) kd = k

        if (DOCANUTO/=0) then
          w2(:,:,k) = 0.667*TKE(:,:,k)

          hl2(:,:,k) = 0.5*( hl2_edge(:,:,kd) + hl2_edge(:,:,ku) )

          wrk1 = 0.5*(qt2_edge(:,:,kd)+qt2_edge(:,:,ku))              ! averaging gradient production term
          if (DOPROGQT2 /= 0) then
            qt2(:,:,k) = (qt2(:,:,k)+DT*wrk1) / (1. + DT*QT2TUNE*1.5e-4)
          else
            qt2(:,:,k) = QT2TUNE*ISOTROPY(:,:,k)*wrk1
          end if

          hlqt(:,:,k) = 0.5*( hlqt_edge(:,:,kd) + hlqt_edge(:,:,ku) )

          wqt(:,:,k)  = 0.5*( wqt_edge(:,:,kd) + wqt_edge(:,:,ku) )
        else
          onemmf = 1.0 - MFFRC(:,:,k)

          w2(:,:,k) = onemmf*0.667*TKE(:,:,k) !+ MFW2(:,:,k)

!          hl2(:,:,k) = onemmf*0.5*( hl2_edge(:,:,kd) + hl2_edge(:,:,ku) )  !+ MFHL2(:,:,k)
          hl2(:,:,k) = 0.5*( hl2_edge(:,:,kd) + hl2_edge(:,:,ku) )
          hl2diag(:,:,k) = 0.5*( hl2_edge_nomf(:,:,kd) + hl2_edge_nomf(:,:,ku) )

          wrk1 = 0.5*(qt2_edge(:,:,kd)+qt2_edge(:,:,ku))              ! averaging gradient production term
          if (DOPROGQT2 /= 0) then
            wrk3 = QT2TUNE*1.5e-4 ! dissipation
            qt2(:,:,k) = (qt2(:,:,k)+DT*wrk1) / (1. + DT*wrk3)
            qt2diag(:,:,k) = QT2TUNE*ISOTROPY(:,:,k)*0.5*(qt2_edge_nomf(:,:,kd)+qt2_edge_nomf(:,:,ku))
          else
!            qt2(:,:,k) = QT2TUNE*ISOTROPY(:,:,k)*wrk1 + MFQT2(:,:,k)
            qt2(:,:,k) = QT2TUNE*ISOTROPY(:,:,k)*wrk1
            qt2diag(:,:,k) = 1.0*ISOTROPY(:,:,k)*0.5*(qt2_edge_nomf(:,:,kd)+qt2_edge_nomf(:,:,ku))
          end if

          hlqt(:,:,k) = onemmf*0.5*( hlqt_edge(:,:,kd) + hlqt_edge(:,:,ku) ) + MFHLQT(:,:,k)
          hlqtdiag(:,:,k) = 0.5*( hlqt_edge(:,:,kd) + hlqt_edge(:,:,ku) )

          wqt(:,:,k)  = onemmf*0.5*( wqt_edge(:,:,kd) + wqt_edge(:,:,ku) ) + MFWQT(:,:,k)

        end if

        whl(:,:,k)  = onemmf*0.5*( whl_edge(:,:,kd) + whl_edge(:,:,ku) ) + MFWHL(:,:,k)
        whl_can(:,:,k) = onemmf*0.5*( whl_edge(:,:,kd) + whl_edge(:,:,ku) + mfwhl(:,:,kd) + mfwhl(:,:,ku))

        ! Restrict QT variance, 2-25% of total water.
        qt2(:,:,k) = max(min(qt2(:,:,k),(0.25*QT(:,:,k))**2),(0.02*QT(:,:,k))**2)
        qt2diag(:,:,k) = max(min(qt2diag(:,:,k),(0.25*QT(:,:,k))**2),(0.02*QT(:,:,k))**2)

        hl2(:,:,k) = max(min(hl2(:,:,k),HL2MAX),HL2MIN)
        hl2diag(:,:,k) = max(min(hl2diag(:,:,k),HL2MAX),HL2MIN)

        ! Ensure realizibility
        hlqt(:,:,k) = sign( min( abs(hlqt(:,:,k)), sqrt(hl2(:,:,k)*qt2(:,:,k)) ), hlqt(:,:,k) )
        hlqtdiag(:,:,k) = sign( min( abs(hlqtdiag(:,:,k)), sqrt(hl2diag(:,:,k)*qt2diag(:,:,k)) ), hlqtdiag(:,:,k) )

    end do

    ! Update PDF_A
    if (SKEW_TDIS.gt.0.) then
!      pdf_a = (pdf_a+mffrc+2.*0.5*(cnv_mfc(:,:,1:LM)+cnv_mfc(:,:,0:LM-1)))/(1.+DT/AFRC_TSCALE)
      pdf_a = (pdf_a+mffrc*DT/SKEW_TGEN)/(1.+DT/SKEW_TDIS)
    else
      pdf_a = pdf_a/(1.-DT/SKEW_TDIS)
    end if
    where (mffrc.gt.pdf_a)
      pdf_a = mffrc
    end where
    pdf_a = min(0.5,max(0.,pdf_a))

  if (DOCANUTO==0) then
    qt3 = ( qt3 + max(MFQT3,0.)*DT/SKEW_TGEN ) / ( 1. + DT/SKEW_TDIS )
    hl3 = MFHL3
    w3  = MFW3
  else

! pre-define adzl,
  do k=2,LM
    km1 = k - 1
    do j=1,JM
      do i=1,IM
        adzl(i,j,km1) = (ZLE(i,j,k) - ZLE(i,j,km1))  ! level thickness
      enddo
    end do
   end do

    do k=2,LM

      kb = k-1
      kc = k+1

      do j=1,JM
        do i=1,IM

          if(k == 1) then
            kb = 1
            kc = 2
            thedz  = adzl(i,j,kc)
            thedz2 = thedz
          elseif(k == LM) then
            kb = LM-1
            kc = LM
            thedz  = adzl(i,j,k)
            thedz2 = thedz
          else
            thedz  = adzl(i,j,k)
            thedz2 = adzl(i,j,k)+adzl(i,j,kb)
          endif

!          brunt = (bet(i,j,k)/thedz)*(thv(i,j,kc)-thv(i,j,kb))

          thedz     = 1. / thedz
          thedz2    = 1. / thedz2

          iso       = 0.5*(isotropy(i,j,k)+isotropy(i,j,kb))
          isosqr    = iso*iso ! Two-level average of "return-to-isotropy" time scale squared
          buoy_sgs2 = isosqr*0.5*(brunt(i,j,k)+brunt(i,j,kb))
          bet2      = 0.5*(bet(i,j,k)+bet(i,j,kb))  !Two-level average of BV frequency squared
! Compute functions f0-f5, see Eq, 8 in C01 (B.8 in Pete's dissertation)

          avew = 0.5*(0.667*TKE(i,j,k)+0.667*TKE(i,j,kb))
          if (abs(avew).ge.1e10) avew = sign(1e10,avew)

          cond = 1.2*sqrt(max(1.0e-20,2.*avew*avew*avew))
          wrk1b = bet2*iso
          wrk2b = thedz2*wrk1b*wrk1b*iso
          wrk3b = hl2diag(i,j,kc) - hl2diag(i,j,kb)

          f0   = wrk2b * wrk1b * whl_can(i,j,k) * wrk3b

          wrk  = whl_can(i,j,kc) - whl_can(i,j,kb)

          f1   = wrk2b * (wrk*whl_can(i,j,k) + 0.5*avew*wrk3b)

          wrk1b = bet2*isosqr
          f2   = thedz*wrk1b*whl_can(i,j,k)*0.667*(TKE(i,j,k)-TKE(i,j,kb))     &
               + (thedz2+thedz2)*bet(i,j,k)*isosqr*avew*wrk

          f3   = thedz2*wrk1b*wrk*avew + thedz*bet2*isosqr*(whl_can(i,j,k)*(tke(i,j,k)-tke(i,j,kb)))

          wrk1b = thedz*iso*avew
          f4   = wrk1b*(0.667*TKE(i,j,k)-0.667*TKE(i,j,kb) + tke(i,j,k)-tke(i,j,kb))

          f5   = wrk1b*0.667*(TKE(i,j,k)-TKE(i,j,kb))

! Compute the "omega" terms, see Eq. 6 in C01 (B.6 in Pete's dissertation)
          dum = 1.-a5*buoy_sgs2
          if (abs(dum).le.1e-20) dum = sign(1e-20,dum)

          omega0 = a4 / dum
          omega1 = omega0 / (c+c)
          omega2 = omega1*f3+(5./4.)*omega0*f4

! Compute the X0, Y0, X1, Y1 terms,  see Eq. 5 a-b in C01  (B.5 in Pete's dissertation)
          dum = 1.-(a1+a3)*buoy_sgs2
          if (abs(dum).le.1e-20) dum = sign(1e-20,dum)

          wrk1b = 1.0 / dum
          dum = 1.-a3*buoy_sgs2
          if (abs(dum).le.1e-20) dum = sign(1e-20,dum)

          wrk2b = 1.0 / dum
          X0   = wrk1b * (a2*buoy_sgs2*(1.-a3*buoy_sgs2))
          Y0   = wrk2b * (2.*a2*buoy_sgs2*X0)
          X1   = wrk1b * (a0*f0+a1*f1+a2*(1.-a3*buoy_sgs2)*f2)
          Y1   = wrk2b * (2.*a2*(buoy_sgs2*X1+(a0/a1)*f0+f1))

! Compute the A0, A1 terms,  see Eq. 5d in C01 (B.5 in Pete's dissertation)
          AA0 = omega0*X0 + omega1*Y0
          AA1 = omega0*X1 + omega1*Y1 + omega2

! Finally, we have the third moment of w, see Eq. 4c in C01 (B.4 in Pete's dissertation)
! cond is an estimate of third moment from second oment - If the third moment is larger
! than the estimate - limit w3.

           dum = c-1.2*X0+AA0
           if (abs(dum).le.1e-20) dum = sign(1e-20,dum)

           w3can(i,j,k) = max(-cond, min(cond, (AA1-1.2*X1-1.5*f5)/dum))
! Implemetation of the C01 approach in this subroutine is nearly complete
! (the missing part are Eqs. 5c and 5e which are very simple)
! therefore it's easy to diagnose other third order moments obtained in C01 using this code.

           end do
      end do
    end do
    do j=1,JM
      do i=1,IM
        w3can(i,j,LM) = w3can(i,j,LM-1)
      enddo
    enddo
!    w3 = w3can

!!   skew_w = w3 / w2**1.5
!   qt3 = 1.2*w3*(qt2/w2)**1.5
!   hl3 = w3 * (hl2 / w2)**1.5

  end if ! DOCANUTO conditional

 end subroutine update_moments

end module shoc
