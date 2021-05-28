#define EDMF_DIAG 1

module edmf_mod

use MAPL_ConstantsMod, only: mapl_grav, mapl_cp, mapl_alhl, mapl_p00, mapl_vireps, mapl_alhs, mapl_kappa
!use MAPL_Mod, only:          mapl_undef

implicit none

real, parameter ::     &
     WSTARmin = 1.e-3, &
     zpblmin  = 100.,  &
     onethird = 1./3., &
     au0      = 0.2,  &
     delta    = 2.E-3

contains

!
! edmf
!
subroutine run_edmf(IM, JM, LM, numup, iras, jras, kbotp, &                         ! in
                    discrete_type, implicit_flag, stochastic_flag, &                ! in
                    plume_type, test_flag, debug_flag, &                            ! in
                    th00, dt, zl, zle, ple, rho, rhoe, exf, &                       ! in
                    u, v, thl, qt, qv, ql, qi, thv, &                               ! in         
                    ui, vi, thli, qti, qvi, qli, qii, thvi, &                       ! in
                    ustar, sh, evap, ice_ramp, &                                    ! in
                    pwmin, pwmax, AlphaW, AlphaQT, AlphaTH, c_kh_mf, &              ! in
                    ET, L0, ENT0, EDfac, EntWFac, Wa, Wb, &                         ! in
                    zpbl, &                                                         ! inout
                    edmfdrya, edmfmoista, &                                         ! out
                    edmfdryw, edmfmoistw, &                                         ! out
                    edmfdryqt, edmfmoistqt, &                                       ! out
                    edmfdrythl, edmfmoistthl, &                                     ! out
                    edmfdryu, edmfmoistu,  &                                        ! out
                    edmfdryv, edmfmoistv,  &                                        ! out
                    edmfmoistqc, &                                                  ! out
                    tke_mf, &                                                       ! out (diagnostics)
                    ae, awu, awv, aw, aws, awqv, awql, awqi, Kh_mf, Kh_t, Kh_q, &   ! out (for solver)
                    whl_mf, wqt_mf, wthv_mf, &                                      ! out (for MYNN-EDMF inconsistent partitioning)
                    buoyf, mfw2, mfw3, mfqt3, mfwqt, mfqt2, mfhl2, mfqthl, mfwhl, & ! out (for SHOC)
                    au_full, hlu_full, qtu_full, acu_full, Tu_full, qlu_full, &     ! out (for MOIST)
#ifdef EDMF_DIAG
                    qt_plume1,qt_plume2,qt_plume3,qt_plume4,qt_plume5, &
                    qt_plume6,qt_plume7,qt_plume8,qt_plume9,qt_plume10, &
                    thl_plume1,thl_plume2,thl_plume3,thl_plume4,thl_plume5, &
                    thl_plume6,thl_plume7,thl_plume8,thl_plume9,thl_plume10, &
#endif
                    au, wu, Mu, E, D, wdet)                                         ! out (for MYNN-EDMF consistent partitioning)
  
  ! Inputs
  integer, intent(in)                     :: IM, JM, LM, numup, discrete_type, implicit_flag, &
                                             stochastic_flag, plume_type, ET, kbotp, test_flag, debug_flag
  integer, dimension(IM,JM), intent(in)   :: iras, jras
  real, dimension(IM,JM,LM), intent(in)   :: u, v, thl, qt, thv, qv, ql, qi, zl, exf, rho
  real, dimension(IM,JM,0:LM), intent(in) :: zle, ple, rhoe, ui, vi, thli, qti, qvi, qli, qii, thvi
  real, dimension(IM,JM), intent(in)      :: ustar, sh, evap, L0
  real, dimension(IM,JM), intent(inout)   :: zpbl
  real, intent(in)                        :: th00, ice_ramp, dt, EntWFac, ENT0, Wa, Wb, pwmin, pwmax, &
                                             AlphaW, AlphaQT, AlphaTH, EDfac, c_kh_mf
 
  ! Outputs
  real, dimension(IM,JM,0:LM), intent(out) :: edmfdrya, edmfmoista, edmfdryw, edmfmoistw, &
                                              edmfdryqt, edmfmoistqt, edmfdrythl, edmfmoistthl, &
                                              edmfdryu, edmfmoistu, edmfdryv, edmfmoistv, &
                                              edmfmoistqc, tke_mf, &
                                              ae, aw, aws, awqv, awql, awqi, awu, awv, &
#ifdef EDMF_DIAG
                                              qt_plume1,qt_plume2,qt_plume3,qt_plume4,qt_plume5, &
                                              qt_plume6,qt_plume7,qt_plume8,qt_plume9,qt_plume10, &
                                              thl_plume1,thl_plume2,thl_plume3,thl_plume4,thl_plume5, &
                                              thl_plume6,thl_plume7,thl_plume8,thl_plume9,thl_plume10, &
#endif
                                              whl_mf, wqt_mf, wthv_mf, au, Mu, wu, Kh_mf, Kh_t, Kh_q

  real, dimension(IM,JM,LM), intent(out) :: buoyf, mfw2, mfw3, mfqt3, mfqt2, mfwqt, &
                                            mfhl2, mfqthl, mfwhl, E, D, wdet, &
                                            au_full, hlu_full, qtu_full, acu_full, &
                                            Tu_full, qlu_full

  real, dimension(numup,IM,JM) :: upm, upw, upthl, upqt, upql, upqi, upa, upu, upv, upthv

  integer :: i, j, k, km1, kp1, iup, kbot

  real :: wthv, wstar, qstar, thstar, sigmaw, sigmaqt, sigmath, z0, wmin, wmax, wlv, wtv, wp, &
          B, QTn, THLn, THVn, QCn, Un, Vn, Wn, Wn2, Mn, EntEXP, EntEXPU, EntW, wf, &
          stmp, QTsrfF, THVsrfF, mft, mfthvt, dzle, idzle, ifac, test, mft_work, mfthvt_work, &
          goth00, thlu_full, work, exfh, dsdz, dqdz, wu0, dthv0, dqt0, thvu0, qtu0, thv_high, thv_low, &
          thvmin, thvmax

  ! Temporary (too slow; need to figure out how random number generator works)
  integer, dimension(numup,IM,JM,LM)  :: enti
  real, dimension(numup,IM,JM,LM)     :: entf, ent

  real, dimension(IM,JM,LM) :: thlu, qtu, qlu, edmfmoistql

  real, dimension(IM,JM,0:LM) :: aw2, ahl2, aqt2, aw3, aqt3, aqthl

  integer, dimension(2)  :: seedmf, the_seed

  logical :: stop_cond

  real, dimension(numup) :: upa_out
  logical, dimension(IM,JM)       :: work_flag
  logical, dimension(numup,IM,JM) :: active_updraft

  ! Thermal plume stuff 
  real                      :: E_work, D_work
  real, dimension(IM,JM,LM) :: f_thermal, ent_sl

  kbot = LM - kbotp

  goth00 = mapl_grav/th00

  !
  ! Initialize arrays
  !

  ! Outputs that are diagnostic updraft statistics
  edmfdrya     = 0.
  edmfmoista   = 0.
  edmfdryw     = 0.
  edmfmoistw   = 0.
  edmfdryqt    = 0.
  edmfmoistqt  = 0.
  edmfdrythl   = 0.
  edmfmoistthl = 0.
  edmfdryu     = 0.
  edmfmoistu   = 0.
  edmfdryv     = 0.
  edmfmoistv   = 0.
  edmfmoistqc  = 0.
  edmfmoistql  = 0.

  ! Outputs needed for solver in TURBULENCE
  ae   = EDfac 
  aw   = 0.
  awu  = 0.
  awv  = 0.
  aws  = 0.
  awqv = 0.
  awql = 0.
  awqi = 0.

  ! Outputs need for MYNN-ED
  whl_mf  = 0.
  wqt_mf  = 0.
  wthv_mf = 0.

  ! 
  tke_mf = 0.

  ! Intermediate quantities for SHOC cloud scheme in MOIST
  aw2   = 0.
  ahl2  = 0.
  aqt2  = 0.
  aqthl = 0.
  aw3   = 0.
  aqt3  = 0.

  ! Outputs needed for SHOC in TURBULENCE
  buoyf = 0.

  ! Outputs needed for SHOC cloud scheme in MOIST
  mfw2   = 0.
  mfw3   = 0.
  mfqt3  = 0.
  mfqt2  = 0.
  mfwqt  = 0.
  mfhl2  = 0.
  mfqthl = 0.
  mfwhl  = 0.

  ! 
  au   = 0.
  wu   = 0.
  Mu   = 0.
  E    = 0.
  D    = 0.
  wdet = 0.

#ifdef EDMF_DIAG
  qt_plume1  = qti
  qt_plume2  = qti
  qt_plume3  = qti
  qt_plume4  = qti
  qt_plume5  = qti
  qt_plume6  = qti
  qt_plume7  = qti
  qt_plume8  = qti
  qt_plume9  = qti
  qt_plume10 = qti

  thl_plume1  = thvi
  thl_plume2  = thvi
  thl_plume3  = thvi
  thl_plume4  = thvi
  thl_plume5  = thvi
  thl_plume6  = thvi
  thl_plume7  = thvi
  thl_plume8  = thvi
  thl_plume9  = thvi
  thl_plume10 = thvi
#endif

  !
  !
  !

  ! Loop across surface tiles
  do j = 1,JM
  do i = 1,IM
     seedmf(1) = 1000000*( 100*thl(i,j,1) - INT(100*thl(i,j,1)) )
     seedmf(2) = 1000000*( 100*thl(i,j,2) - INT(100*thl(i,j,2)) )

     the_seed(1) = seedmf(1)*iras(i,j) + seedmf(2)*jras(i,j)
     the_seed(2) = seedmf(1)*jras(i,j) + seedmf(2)*iras(i,j)
     the_seed(1) = the_seed(1)*seedmf(1)/( seedmf(2) + 10 )
     the_seed(2) = the_seed(2)*seedmf(1)/( seedmf(2) + 10 )
     if(the_seed(1) == 0) the_seed(1) =  5
     if(the_seed(2) == 0) the_seed(2) = -5

     do k = 1,LM
        km1 = k - 1
        do iup = 1,numup
           entf(iup,i,j,k) = ( zle(i,j,km1) - zle(i,j,k) )/L0(i,j)
        end do
     end do
  end do
  end do

  ! Get surface layer organized entrainment
  if ( plume_type == 1 ) then
     call A_star_closure(IM, JM, LM, th00, zle, zl, ple, ice_ramp, & ! in
                         rho, rhoe, thl, qt, thv, debug_flag, &      ! in
                         E, f_thermal, ent_sl)                       ! out
  end if

  !
  ! Initialize updrafts
  !
  do j = 1,JM
  do i = 1,IM
     ent = 0.

     zpbl(i,j) = max( zpbl(i,j), zpblmin )

     wthv = sh(i,j)/mapl_cp + mapl_vireps*th00*evap(i,j)

     if ( wthv > 0. .and. thv(i,j,LM-1) < thv(i,j,LM) ) then
        wstar = max( wstarmin, (goth00*wthv*zpbl(i,j))**onethird )

        qstar  = evap(i,j)/wstar
        thstar = wthv/wstar
        
        sigmaW  = AlphaW*wstar
        sigmaQT = AlphaQT*max( qstar, 0. )
        sigmaTH = AlphaTH*max( thstar, 0. )

        wmin = sigmaW*pwmin ! vertical velocity of least energetic updraft
        wmax = sigmaW*pwmax ! "        "        "  most  "         "
        
        work = 1./( sqrt(2.)*sigmaW )

        if ( plume_type == 1 ) then
           wu0   = E(i,j,LM)*( zle(i,j,LM-1) - zle(i,j,LM) )/( au0*rhoe(i,j,LM-1) )
           dthv0 = th00*wu0**2./( 2.*mapl_grav*( zle(i,j,LM-1) - zle(i,j,LM) ) )
           dqt0  = dthv0*evap(i,j)/wthv

           thvu0 = thv(i,j,LM) + dthv0
           qtu0  = qt(i,j,LM)  + dqt0

           ! Temporary
           sigmaTH = 0.75*dthv0
           sigmaQT = 1.3*dqt0

           thvmin = thv(i,j,LM)
           thvmax = thvu0 + 3.*sigmaTH

           work = 1./( sqrt(2.)*sigmaTH )
        end if

        QTsrfF  = 0.
        THVsrfF = 0.
        do iup = 1,numup
           active_updraft(iup,i,j) = .true.
           
           wlv = wmin + ( wmax - wmin )*real(iup-1)/real(numup)
           wtv = wmin + ( wmax - wmin )*real(iup)/real(numup)

           if ( plume_type == 0 ) then
              upw(iup,i,j)   = 0.5*( wlv + wtv ) 
              upa(iup,i,j)   = 0.5*erf(wtv*work) - 0.5*erf(wlv*work)       
              upu(iup,i,j)   = ui(i,j,kbot)
              upv(iup,i,j)   = vi(i,j,kbot)
              upqt(iup,i,j)  = qti(i,j,kbot)  + 0.32*upw(iup,i,j)*sigmaQT/sigmaW ! 0.32~=0.58*0.55 (Stull 1988)
              upthv(iup,i,j) = thvi(i,j,kbot) + 0.58*upw(iup,i,j)*sigmaTH/sigmaW ! Stull 1988
           else
              thv_low  = thvmin + ( thvmax - thvmin )*real(iup-1)/real(numup)
              thv_high = thvmin + ( thvmax - thvmin )*real(iup)/real(numup)

              upthv(iup,i,j) = 0.5*( thv_low + thv_high )
              upa(iup,i,j)   = au0*( 0.5*erf(( thv_high - thvu0 )*work) - 0.5*erf( ( thv_low - thvu0 )*work) )
              upw(iup,i,j)   = wu0
              upu(iup,i,j)   = u(i,j,LM)
              upv(iup,i,j)   = v(i,j,LM)
              upqt(iup,i,j)  = max( 0., qtu0 + 0.4*( upthv(iup,i,j) - thvu0 )*sigmaQT/sigmaTH )
!              upqt(iup,i,j)  = qtu0
!              upa(iup,i,j)   = au0/real(numup)
!              upthv(iup,i,j) = thv(i,j,LM) + dthv0
           end if

           upm(iup,i,j) = rhoe(i,j,kbot)*upa(iup,i,j)*upw(iup,i,j)

           ! For stability make sure that the surface mass-fluxes are not more than their 
           ! values computed from the surface scheme
           QTsrfF  = QTsrfF  + upa(iup,i,j)*upw(iup,i,j)*( upqt(iup,i,j)  - qti(i,j,kbot) )
           THVsrfF = THVsrfF + upa(iup,i,j)*upw(iup,i,j)*( upthv(iup,i,j) - thvi(i,j,kbot) )
        end do

        ! Change surface THV so that the fluxes from the mass flux equal prescribed values
        if ( kbot == 0 .and. THVsrfF > wthv ) then
           upthv(:,i,j) = ( upthv(:,i,j) - thvi(i,j,kbot) )*wthv/THVsrfF + thvi(i,j,kbot)
        end if

        ! Change surface QT so that the fluxes from the mass flux equal prescribed values 
        ! We do not need to worry about the negative values as they should not exist
        if ( kbot == 0. .and. QTsrfF > evap(i,j) .and. evap(i,j) > 0. )  then
           upqt(:,i,j) = ( upqt(:,i,j) - qti(i,j,kbot) )*evap(i,j)/QTsrfF + qti(i,j,kbot)
        end if

        ! Compute condensation and thl, ql, qi
        do iup = 1,numup
           if ( test_flag == 0 ) then
              call condensation_edmfA(upthv(iup,i,j), upqt(iup,i,j), ple(i,j,kbot), &
                                      upthl(iup,i,j), upql(iup,i,j), upqi(iup,i,j), ice_ramp)
           else
              upthl(iup,i,j) = upthv(iup,i,j)/( 1. + mapl_vireps*upqt(iup,i,j) )
              upql(iup,i,j)  = 0.
              upqi(iup,i,j)  = 0.
           end if
        end do
     else
        active_updraft(:,i,j) = .false.
     end if    ! wthv > 0.
  end do
  end do

  if ( debug_flag /= 0 ) then
     write(*,'(I4,12F7.2)') LM-1, sum(100.*upa(:,1,1)), 0., 100.*upa(:,1,1)
  end if

  !
  ! Integrate updraft equations and save updraft statistics
  !
  do k = kbot,2,-1
     km1 = k - 1
     kp1 = k + 1
     
     do j = 1,JM
     do i = 1,IM
        dzle  = zle(i,j,km1) - zle(i,j,k)
        idzle = 1./dzle
        
        exfh = (ple(i,j,k)/mapl_p00)**mapl_kappa

        mfthvt = 0.
        mft    = 0.

        do iup = 1,numup
           if ( active_updraft(iup,i,j) ) then
              ! Initialize
              
              !
              ! Sample updraft statistics
              !

              ! Quantities required for solver
              aw(i,j,k)    = aw(i,j,k)    + upa(iup,i,j)*upw(iup,i,j)
              aw2(i,j,k)   = aw2(i,j,k)   + upa(iup,i,j)*upw(iup,i,j)**2.
              ahl2(i,j,k)  = ahl2(i,j,k)  + upa(iup,i,j)*( exfh*( upthl(iup,i,j) - thli(i,j,k) ) )**2.
              aqt2(i,j,k)  = aqt2(i,j,k)  + upa(iup,i,j)*( upqt(iup,i,j) - qti(i,j,k) )**2.
              aqthl(i,j,k) = aqthl(i,j,k) + upa(iup,i,j)*exfh*( upthl(iup,i,j) - thli(i,j,k) )&
                                                                    *( upqt(iup,i,j) - qti(i,j,k) )
              aw3(i,j,k)   = aw3(i,j,k)   + upa(iup,i,j)*upw(iup,i,j)**3.
              aqt3(i,j,k)  = aqt3(i,j,k)  + upa(iup,i,j)*( upqt(iup,i,j) - qti(i,j,k) )**3.

              if ( implicit_flag == 0 ) then
                 stmp =   mapl_cp*exfh*( upthl(iup,i,j) - thli(i,j,k) ) &
                        + mapl_alhs*( upqi(iup,i,j) - qii(i,j,k) )             &
                        + mapl_alhl*( upql(iup,i,j) - qli(i,j,k) ) 

                 awu(i,j,k)  = awu(i,j,k)  + upa(iup,i,j)*upw(iup,i,j)*( upu(iup,i,j)  - ui(i,j,k) )
                 awv(i,j,k)  = awv(i,j,k)  + upa(iup,i,j)*upw(iup,i,j)*( upv(iup,i,j)  - vi(i,j,k) )
                 aws(i,j,k)  = aws(i,j,k)  + upa(iup,i,j)*upw(iup,i,j)*stmp
                 awqv(i,j,k) = awqv(i,j,k) + upa(iup,i,j)*upw(iup,i,j)&
                                             *( upqt(iup,i,j) - upql(iup,i,j) - upqi(iup,i,j) - qvi(i,j,k))
                 awql(i,j,k) = awql(i,j,k) + upa(iup,i,j)*upw(iup,i,j)*( upql(iup,i,j) - qli(i,j,k) )
                 awqi(i,j,k) = awqi(i,j,k) + upa(iup,i,j)*upw(iup,i,j)*( upqi(iup,i,j) - qii(i,j,k) )
              else
                 stmp = mapl_cp*exfh*upthl(iup,i,j) + mapl_grav*zle(i,j,k) + mapl_alhl*upql(iup,i,j) + mapl_alhs*upqi(iup,i,j) 
                        
                 awu(i,j,k)  = awu(i,j,k)  + upa(iup,i,j)*upw(iup,i,j)*upu(iup,i,j)
                 awv(i,j,k)  = awv(i,j,k)  + upa(iup,i,j)*upw(iup,i,j)*upv(iup,i,j)
                 aws(i,j,k)  = aws(i,j,k)  + upa(iup,i,j)*upw(iup,i,j)*stmp
                 awqv(i,j,k) = awqv(i,j,k) + upa(iup,i,j)*upw(iup,i,j)&
                                             *( upqt(iup,i,j) - upql(iup,i,j) - upqi(iup,i,j) )
                 awql(i,j,k) = awql(i,j,k) + upa(iup,i,j)*upw(iup,i,j)*upql(iup,i,j)
                 awqi(i,j,k) = awqi(i,j,k) + upa(iup,i,j)*upw(iup,i,j)*upqi(iup,i,j)
              end if

              ! Quantities required for MYNN
              whl_mf(i,j,k)  = whl_mf(i,j,k)  + upa(iup,i,j)*upw(iup,i,j)*exfh*( upthl(iup,i,j) - thli(i,j,k) )
              wqt_mf(i,j,k)  = wqt_mf(i,j,k)  + upa(iup,i,j)*upw(iup,i,j)*( upqt(iup,i,j) - qti(i,j,k) )
              wthv_mf(i,j,k) = wthv_mf(i,j,k) + upa(iup,i,j)*upw(iup,i,j)*( upthv(iup,i,j) - thvi(i,j,k) )

              ! Sample dry and moist updraft statistics at half levels
              if ( upql(iup,i,j) > 0. .or. upqi(iup,i,j) > 0. ) then
                 edmfmoista(i,j,k)   = edmfmoista(i,j,k)   + upa(iup,i,j)
                 edmfmoistw(i,j,k)   = edmfmoistw(i,j,k)   + upa(iup,i,j)*upw(iup,i,j)
                 edmfmoistqt(i,j,k)  = edmfmoistqt(i,j,k)  + upa(iup,i,j)*upqt(iup,i,j)
                 edmfmoistthl(i,j,k) = edmfmoistthl(i,j,k) + upa(iup,i,j)*upthl(iup,i,j)
                 edmfmoistu(i,j,k)   = edmfmoistu(i,j,k)   + upa(iup,i,j)*upu(iup,i,j)
                 edmfmoistv(i,j,k)   = edmfmoistv(i,j,k)   + upa(iup,i,j)*upv(iup,i,j)
                 edmfmoistqc(i,j,k)  = edmfmoistqc(i,j,k)  + upa(iup,i,j)*( upql(iup,i,j) + upqi(iup,i,j) )
                 edmfmoistql(i,j,k)  = edmfmoistql(i,j,k)  + upa(iup,i,j)*upql(iup,i,j)
              else
                 edmfdrya(i,j,k)   = edmfdrya(i,j,k)   + upa(iup,i,j)
                 edmfdryw(i,j,k)   = edmfdryw(i,j,k)   + upa(iup,i,j)*upw(iup,i,j)
                 edmfdryqt(i,j,k)  = edmfdryqt(i,j,k)  + upa(iup,i,j)*upqt(iup,i,j)
                 edmfdrythl(i,j,k) = edmfdrythl(i,j,k) + upa(iup,i,j)*upthl(iup,i,j)
                 edmfdryu(i,j,k)   = edmfdryu(i,j,k)   + upa(iup,i,j)*upu(iup,i,j)
                 edmfdryv(i,j,k)   = edmfdryv(i,j,k)   + upa(iup,i,j)*upv(iup,i,j)
              end if

#ifdef EDMF_DIAG
              if ( iup == 1 ) then
                 qt_plume1(i,j,k)  = upqt(iup,i,j)
                 thl_plume1(i,j,k) = upthv(iup,i,j)
              elseif ( iup == 2 ) then
                 qt_plume2(i,j,k)  = upqt(iup,i,j)
                 thl_plume2(i,j,k) = upthv(iup,i,j)
              elseif ( iup == 3 ) then
                 qt_plume3(i,j,k)  = upqt(iup,i,j)
                 thl_plume3(i,j,k) = upthv(iup,i,j)
              elseif ( iup == 4 ) then
                 qt_plume4(i,j,k)  = upqt(iup,i,j)
                 thl_plume4(i,j,k) = upthv(iup,i,j)
              elseif ( iup == 5 ) then
                 qt_plume5(i,j,k)  = upqt(iup,i,j)
                 thl_plume5(i,j,k) = upthv(iup,i,j)
              elseif ( iup == 6 ) then
                 qt_plume6(i,j,k)  = upqt(iup,i,j)
                 thl_plume6(i,j,k) = upthv(iup,i,j)
              elseif ( iup == 7 ) then
                 qt_plume7(i,j,k)  = upqt(iup,i,j)
                 thl_plume7(i,j,k) = upthv(iup,i,j)
              elseif ( iup == 8 ) then
                 qt_plume8(i,j,k)  = upqt(iup,i,j)
                 thl_plume8(i,j,k) = upthv(iup,i,j)
              elseif ( iup == 9 ) then
                 qt_plume9(i,j,k)  = upqt(iup,i,j)
                 thl_plume9(i,j,k) = upthv(iup,i,j)
              elseif ( iup == 10 ) then
                 qt_plume10(i,j,k)  = upqt(iup,i,j)
                 thl_plume10(i,j,k) = upthv(iup,i,j)
              end if
#endif


              !
              !
              !
              
              ! Compute fractional entrainment rate
              if ( L0(i,j) > 0. ) then
                 if ( stochastic_flag /= 0 ) then
                    ent(iup,i,j,k) = entf(iup,i,j,k)*ent0*idzle
                 else
                    call Poisson(1, 1, 1, 1, entf(iup,i,j,k), enti(iup,i,j,k), the_seed)
                    ent(iup,i,j,k) = real(enti(iup,i,j,k))*ent0*idzle
                 end if
              else
                 ! negative L0 means 0 entrainment  
                 ent(:,i,j,:) = 0. ! check
              end if

              ! Integrate plume budgets acroos one vertical step
              if ( plume_type == 0 ) then ! JPL entraining plume model
                 ! sample entrainment
                 E(i,j,k) = E(i,j,k) + rhoe(i,j,k)*upa(iup,i,j)*upw(iup,i,j)*ent(iup,i,j,k)
                 
                 !
                 EntExp  = exp(-ent(iup,i,j,k)*dzle)
                 EntExpU = exp(-ent(iup,i,j,k)*dzle*EntWFac)
                 
                 ! Thermodynamic variables in updraft
                 QTn  = qt(i,j,k)*( 1. - EntExp )  + upqt(iup,i,j)*EntExp
                 THLn = thl(i,j,k)*( 1. - EntExp ) + upthl(iup,i,j)*EntExp
                 Un   = u(i,j,k)*( 1. - EntExpU )  + upu(iup,i,j)*EntExpU
                 Vn   = v(i,j,k)*( 1. - EntExpU )  + upv(iup,i,j)*EntExpU

                 ! Condensation
                 if ( test_flag == 0 ) then
                    call condensation_edmf(QTn, THLn, ple(i,j,km1), THVn, QCn, wf, ice_ramp)
                 else
                    QCn  = 0.
                    THVn = THLn*( 1. + mapl_vireps*QTn )
                 end if

                 ! Buoyancy
                 B = goth00*( 0.5*( THVn + upthv(iup,i,j) ) - thv(i,j,k) )

                 ! Vertical velocity
                 WP = Wb*ent(iup,i,j,k)
                 if (WP == 0.) then
                    Wn2 = upw(iup,i,j)**2. + 2.*Wa*B*dzle
                 else
                    EntW = exp(-2.*WP*dzle)
                    Wn2  = EntW*upw(iup,i,j)**2. + Wa*B/WP*( 1. - EntW )
                 end if

              elseif ( plume_type == 1 ) then ! thermal plume
                 ! Integrate vertical velocity budget
                 B = goth00*( upthv(iup,i,j) - thv(i,j,k) )

                 Wn2 = f_thermal(i,j,k)**2.*upw(iup,i,j)**2. + 2.*dzle*B                 

                 if ( Wn2 > 0. ) then
                    if ( B > 0. ) then
                       Mn     = upm(iup,i,j)/f_thermal(i,j,k)
                       E_work = max( 0., ( upm(iup,i,j)/f_thermal(i,j,k) - upm(iup,i,j) )/dzle )
                       D_work = 0.
                    else
                       Mn     = rhoe(i,j,k-1)*upa(iup,i,j)*sqrt(Wn2)
                       E_work = max( 0., ( upm(iup,i,j)/f_thermal(i,j,k) - upm(iup,i,j) )/dzle )
                       D_work = max( 0., E_work - ( Mn - upm(iup,i,j) )/dzle )
                    end if

                    ! Integrate scalar budgets
                    THLn = ( upm(iup,i,j)*upthl(iup,i,j) + dzle*E_work*thl(i,j,k) )/( Mn + dzle*D_work )
                    QTn  = ( upm(iup,i,j)*upqt(iup,i,j) + dzle*E_work*qt(i,j,k) )/( Mn + dzle*D_work )
                    Un   = ( upm(iup,i,j)*upu(iup,i,j) + dzle*E_work*u(i,j,k) )/( Mn + dzle*D_work )
                    Vn   = ( upm(iup,i,j)*upu(iup,i,j) + dzle*E_work*v(i,j,k) )/( Mn + dzle*D_work )

                    ! Condensation
                    if ( test_flag == 0 ) then
                       call condensation_edmf(QTn, THLn, ple(i,j,km1), THVn, QCn, wf, ice_ramp)
                    else
                       QCn  = 0.
                       THVn = THLn*( 1. + mapl_vireps*QTn )
                    end if

                    if ( debug_flag /= 0 .and. iup == 1 ) then
!                       write(*,*) k, B, upw(iup,i,j), sqrt(Wn2)
                    end if
                 end if
              end if
              
              ! Intermediate quantities 
              mfthvt_work = upa(iup,i,j)*upw(iup,i,j)*upthv(iup,i,j)
              mft_work    = upa(iup,i,j)*upw(iup,i,j)

              ! Determine stopping condition
              stop_cond = Wn2 > 0.
!              if ( plume_type == 0 ) then
!                 stop_cond = Wn2 > 0.
!              else
!                 stop_cond = B > 0.
!              end if

              ! Sample detrainment velocity
              wdet(i,j,k) = wdet(i,j,k) + upa(iup,i,j)*upw(iup,i,j)

              ! Update plume
              if ( stop_cond ) then
                 upw(iup,i,j)   = sqrt(Wn2)
                 upthv(iup,i,j) = THVn
                 upthl(iup,i,j) = THLn
                 upqt(iup,i,j)  = QTn
                 upql(iup,i,j)  = QCn*wf
                 upqi(iup,i,j)  = QCn*( 1. - wf )
                 upu(iup,i,j)   = Un
                 upv(iup,i,j)   = Vn

                 if ( plume_type == 1 ) then
                    upm(iup,i,j) = Mn
                    upa(iup,i,j) = upm(iup,i,j)/( rhoe(i,j,km1)*upw(iup,i,j) )
                 end if

                 mfthvt = mfthvt + 0.5*( upa(iup,i,j)*upw(iup,i,j)*upthv(iup,i,j) + mfthvt_work )
                 mft    = mft    + 0.5*(                upa(iup,i,j)*upw(iup,i,j) + mft_work )

                 upa_out(iup) = upa(iup,i,j)
              else
                 mfthvt = mfthvt + 0.5*mfthvt_work
                 mft    = mft    + 0.5*mft_work

                 active_updraft(iup,i,j) = .false.
                 upa_out(iup)            = 0.
              end if
           end if ! active_updraft(iup,i,j)
        end do ! iup = 1,numup

        ! Outputs needed for SHOC in TURBULENCE
        buoyf(i,j,k) = buoyf(i,j,k) + exf(i,j,k)*( mfthvt - mft*thv(i,j,k) )

        ! Average dry updraft samples
        if ( edmfdrya(i,j,k) > 0. ) then
           edmfdryw(i,j,k)   = edmfdryw(i,j,k)/edmfdrya(i,j,k)
           edmfdryqt(i,j,k)  = edmfdryqt(i,j,k)/edmfdrya(i,j,k)
           edmfdrythl(i,j,k) = edmfdrythl(i,j,k)/edmfdrya(i,j,k)
           edmfdryu(i,j,k)   = edmfdryu(i,j,k)/edmfdrya(i,j,k)
           edmfdryv(i,j,k)   = edmfdryv(i,j,k)/edmfdrya(i,j,k)
        else
           edmfdryw(i,j,k)   = 0.
           edmfdryqt(i,j,k)  = 0.
           edmfdrythl(i,j,k) = thli(i,j,k)
           edmfdryu(i,j,k)   = ui(i,j,k)
           edmfdryv(i,j,k)   = vi(i,j,k)
        end if

        ! Average moist updraft samples
        if ( edmfmoista(i,j,k) > 0. ) then
           edmfmoistw(i,j,k)   = edmfmoistw(i,j,k)/edmfmoista(i,j,k)
           edmfmoistqt(i,j,k)  = edmfmoistqt(i,j,k)/edmfmoista(i,j,k)
           edmfmoistthl(i,j,k) = edmfmoistthl(i,j,k)/edmfmoista(i,j,k)
           edmfmoistu(i,j,k)   = edmfmoistu(i,j,k)/edmfmoista(i,j,k)
           edmfmoistv(i,j,k)   = edmfmoistv(i,j,k)/edmfmoista(i,j,k)
           edmfmoistqc(i,j,k)  = edmfmoistqc(i,j,k)/edmfmoista(i,j,k)
           edmfmoistql(i,j,k)  = edmfmoistql(i,j,k)/edmfmoista(i,j,k)
        else
           edmfmoistw(i,j,k)   = 0.
           edmfmoistqt(i,j,k)  = 0.
           edmfmoistthl(i,j,k) = thli(i,j,k)
           edmfmoistu(i,j,k)   = ui(i,j,k)
           edmfmoistv(i,j,k)   = vi(i,j,k)
           edmfmoistqc(i,j,k)  = 0.
        end if

        ! Average both dry and moist updraft samples together
        au(i,j,k) = edmfdrya(i,j,k) + edmfmoista(i,j,k)
        if ( au(i,j,k) > 0. ) then
           wdet(i,j,k) = wdet(i,j,k)/au(i,j,k)

           if ( edmfmoista(i,j,k) > 0. ) then
              wu(i,j,k)   = (  edmfdrya(i,j,k)*edmfdryw(i,j,k)  &
                             + edmfmoista(i,j,k)*edmfmoistw(i,j,k) )/au(i,j,k)
              thlu(i,j,k) = (  edmfdrya(i,j,k)*edmfdrythl(i,j,k)  &
                             + edmfmoista(i,j,k)*edmfmoistthl(i,j,k) )/au(i,j,k)
              qtu(i,j,k)  = (  edmfdrya(i,j,k)*edmfdryqt(i,j,k)  &
                             + edmfmoista(i,j,k)*edmfmoistqt(i,j,k) )/au(i,j,k)
              qlu(i,j,k)  = edmfmoistql(i,j,k)
           else
              wu(i,j,k)   = edmfdryw(i,j,k)
              thlu(i,j,k) = edmfdrythl(i,j,k)
              qtu(i,j,k)  = edmfdryqt(i,j,k)
              qlu(i,j,k)  = 0.
           end if

           if ( debug_flag /= 0 ) then
              write(*,'(I4,12F7.2)') k-1, sum(100.*upa_out(:)), 100.*edmfmoista(1,1,k-1), 100.*upa_out(:)
!              write(*,*) k-1, 100*upa(1,1,1), upw(1,1,1)
           end if
        else
           wu(i,j,k)   = 0.
           thlu(i,j,k) = thli(i,j,k)
           qtu(i,j,k)  = qti(i,j,k)
           qlu(i,j,k)  = 0.
        end if
        Mu(i,j,k)     = rhoe(i,j,k)*au(i,j,k)*wu(i,j,k)
        tke_mf(i,j,k) = au(i,j,k)/( 1. - au(i,j,k) )*wu(i,j,k)**2.
        ae(i,j,k)     = ( 1. - edmfdrya(i,j,k) - edmfmoista(i,j,k) )*EDfac 
     end do ! i = 1,IM
     end do ! j = 1,JM
  end do ! k = LM,2,-1

  ! Interpolate updraft properties to full levels
  do k = 1,LM
     km1 = k - 1

     do j = 1,JM
     do i = 1,IM
        if ( k >= 2 .and. k <= LM-1 ) then
           ! Compute detrainment rate as residual
           D(i,j,k) = E(i,j,k) - ( Mu(i,j,km1) - Mu(i,j,k) )/( zle(i,j,km1) - zle(i,j,k) )

           if ( au(i,j,k) > 0. .and. au(i,j,km1) > 0. ) then
              au_full(i,j,k)  = 0.5*( au(i,j,k)   + au(i,j,km1) )
              thlu_full       = 0.5*( thlu(i,j,k) + thlu(i,j,km1) )
              qtu_full(i,j,k) = 0.5*( qtu(i,j,k)  + qtu(i,j,km1) )
              
              if ( edmfmoista(i,j,k) > 0. .and. edmfmoista(i,j,km1) > 0. ) then
                 acu_full(i,j,k) = min( au_full(i,j,k), 0.5*( edmfmoista(i,j,k) + edmfmoista(i,j,km1) ) )
                 qlu_full(i,j,k) = 0.5*( qlu(i,j,k) + qlu(i,j,km1) )
!              else if (      edmfmoista(i,j,k+1) == 0. &
!                       .and. edmfmoista(i,j,k) > 0. .and. edmfmoista(i,j,km1) == 0. ) then
!                 acu_full(i,j,k) = min( au_full(i,j,k), edmfmoista(i,j,k) )
!                 qlu_full(i,j,k) = qlu(i,j,k)
              else
                 acu_full(i,j,k) = 0.
                 qlu_full(i,j,k) = 0.
              end if
           else
              au_full(i,j,k)  = 0.
              thlu_full       = thl(i,j,k)
              qtu_full(i,j,k) = qt(i,j,k)
                 
              acu_full(i,j,k) = 0.
              qlu_full(i,j,k) = 0.
           end if
        elseif ( k == LM ) then
           ! Compute entrainment rate near surface as residual
           E(i,j,LM) = Mu(i,j,LM-1)/zle(i,j,LM-1)

           if ( au(i,j,LM-1) > 0. ) then
              au_full(i,j,LM)  = au(i,j,LM-1)
              thlu_full        = thlu(i,j,LM-1)
              qtu_full(i,j,LM) = qtu(i,j,LM-1)
                 
              if ( edmfmoista(i,j,LM-1) > 0. ) then
                 acu_full(i,j,LM) = min( au_full(i,j,LM), edmfmoista(i,j,LM-1) )
                 qlu_full(i,j,LM) = qlu(i,j,LM-1)
              else
                 acu_full(i,j,LM) = 0.
                 qlu_full(i,j,LM) = 0.
              end if
           else
              au_full(i,j,LM)  = 0.
              thlu_full        = thl(i,j,LM)
              qtu_full(i,j,LM) = qt(i,j,LM)
              
              acu_full(i,j,LM) = 0.
              qlu_full(i,j,LM) = 0.
           end if
        elseif ( k == 1 ) then
           au_full(i,j,1)  = 0.
           thlu_full       = thl(i,j,1)
           qtu_full(i,j,1) = qt(i,j,1)

           acu_full(i,j,1) = 0.
           qlu_full(i,j,1) = 0.
        end if

        hlu_full(i,j,k) = exf(i,j,k)*thlu_full + (mapl_grav/mapl_cp)*zl(i,j,k)
        Tu_full(i,j,k)  = exf(i,j,k)*thlu_full + (mapl_alhl/mapl_cp)*qlu_full(i,j,k)

        ! Outputs needed for SHOC
        if ( k >= 2 ) then
           mfw2(i,j,k)   = 0.5*( aw2(i,j,k)    + aw2(i,j,km1) )
           mfw3(i,j,k)   = 0.5*( aw3(i,j,k)    + aw3(i,j,km1) )
           mfhl2(i,j,k)  = 0.5*( ahl2(i,j,k)   + ahl2(i,j,km1) )
           mfqt2(i,j,k)  = 0.5*( aqt2(i,j,k)   + aqt2(i,j,km1) )
           mfqt3(i,j,k)  = 0.5*( aqt3(i,j,k)   + aqt3(i,j,km1) )
           mfwqt(i,j,k)  = 0.5*( wqt_mf(i,j,k) + wqt_mf(i,j,km1) )
           mfqthl(i,j,k) = 0.5*( aqthl(i,j,k)  + aqthl(i,j,km1) )
           mfwhl(i,j,k)  = 0.5*( whl_mf(i,j,k) + whl_mf(i,j,km1) )
        elseif ( k == 1 ) then
           mfw2(i,j,1)   = aw2(i,j,2)
           mfw3(i,j,1)   = aw3(i,j,2)
           mfhl2(i,j,1)  = ahl2(i,j,2)
           mfqt2(i,j,1)  = aqt2(i,j,2)
           mfqt3(i,j,1)  = aqt3(i,j,2)
           mfwqt(i,j,1)  = wqt_mf(i,j,2)
           mfqthl(i,j,1) = aqthl(i,j,2)
           mfwhl(i,j,1)  = whl_mf(i,j,2)
        end if
     end do

     end do
  end do

  ! Test
  Kh_mf(:,:,:)      = 0.
  Kh_mf(:,:,1:LM-1) = c_kh_mf*( zl(:,:,1:LM-1) - zl(:,:,2:LM) )*(Mu(:,:,1:LM-1)/rhoe(:,:,1:LM-1))

  ! Compute zpbl of thermal plume
  if ( plume_type == 1 ) then
     work_flag(:,:) = .true.
     do k = LM-1,1,-1
        
        do j = 1,JM
        do i = 1,IM
           if ( work_flag(i,j) ) then
              zpbl(i,j) = zle(i,j,k)
              if ( edmfdrya(i,j,k) == 0. ) then
                 work_flag(i,j) = .false.
              end if
           end if
        end do
        end do
     end do
  end if

  ! Text output
  KH_t(:,:,:) = 0.
  KH_q(:,:,:) = 0.
  do i = 1,IM
  do j = 1,JM
     do k = 1,LM-1
         dsdz = (  ( mapl_cp*exf(i,j,k)*thl(i,j,k)     + mapl_grav*zl(i,j,k)   + mapl_alhs*qi(i,j,k)   + mapl_alhl*ql(i,j,k) ) &
                 - ( mapl_cp*exf(i,j,k+1)*thl(i,j,k+1) + mapl_grav*zl(i,j,k+1) + mapl_alhs*qi(i,j,k+1) + mapl_alhl*ql(i,j,k+1) ) )&
                /( zl(i,j,k) - zl(i,j,k+1) )
         work = -aws(i,j,k)/dsdz
         if ( work > 0. ) then
!            KH_t(i,j,k) = min( 1., 0.01*work )
!            KH_t(i,j,k) = 0.05*work
            aws(i,j,k)  = aws(i,j,k) + KH_t(i,j,k)*dsdz
         end if

         dqdz = ( qv(i,j,k) - qv(i,j,k+1) )/( zl(i,j,k) - zl(i,j,k+1) )
         work = -awqv(i,j,k)/dqdz
         if ( work > 0. ) then
!            KH_q(i,j,k) = min( 1., 0.01*work )
!            KH_q(i,j,k) = 0.05*work
            awqv(i,j,k) = awqv(i,j,k) + KH_q(i,j,k)*dqdz
         end if
      end do
   end do
   end do

end subroutine run_edmf

!
! A_star_closure
!
subroutine A_star_closure(IM, JM, LM, th00, zle, zl, ple, ice_ramp, & ! in
                          rho, rhoe, thl, qt, thv, debug_flag, &      ! in
                          A, f, ent_sl)                               ! out

  integer, intent(in)                      :: IM, JM, LM, debug_flag
  real, intent(in)                         :: ice_ramp, th00
  real, dimension(IM,JM,LM), intent(in)    :: zl, thl, qt, thv, rho
  real, dimension(IM,JM,0:LM), intent(in)  :: zle, rhoe, ple
  real, dimension(IM,JM,LM), intent(out)   :: A, f, ent_sl

  integer                     :: i, j, k, km1, iter
  real                        :: dthvdz, dz, B, wf, qcu, Mu_next, wu2_next, &
                                 thlu_next, qtu_next, thvu_next, goth00
  integer, dimension(IM,JM)   :: izsl
  real, dimension(IM,JM)      :: A_star_sum, A_star2_int, zi, Mu0, Mu, thlu, qtu, thvu, wu2, wu0
  logical, dimension(IM,JM)   :: conv_flag, sl_flag, test_flag

  goth00 = mapl_grav/th00

  A(:,:,:)        = 0.
  A_star_sum(:,:) = 0.
  Mu0(:,:)        = 0.

  ! Determine if column is convective (surface-driven) and compute A_star
  ! (before normalization) at first model level
  do j = 1,JM
  do i = 1,IM
     conv_flag(i,j) = thv(i,j,LM-1) < thv(i,j,LM)
     sl_flag(i,j)   = conv_flag(i,j)

     if ( .not. conv_flag(i,j) ) then
        izsl(i,j) = -1
     end if
  end do
  end do

  ! Compute A_star (before normalization) and find top of surface layer
  do k = LM,1,-1
     km1 = k - 1

     do j = 1,JM
     do i = 1,IM
        if ( sl_flag(i,j) ) then
           dthvdz = ( thv(i,j,km1) - thv(i,j,k) )/( zl(i,j,km1) - zl(i,j,k) )

           if ( dthvdz < 0. ) then
              A(i,j,k)        = -sqrt(zle(i,j,km1))*dthvdz
              A_star_sum(i,j) = A_star_sum(i,j) + A(i,j,k)*( zle(i,j,km1) - zle(i,j,k) )
           else
              izsl(i,j) = k

              sl_flag(i,j) = .false.
           end if
        end if
     end do
     end do
  end do

  ! Normalize A_star and initialize test plume
  do j = 1,JM
  do i = 1,IM
     if ( conv_flag(i,j) ) then
        wu0(i,j)         = 0.
        A_star2_int(i,j) = 0.
        do k = LM,izsl(i,j)+1,-1         ! Note: this loop needs to be refactored for column-major order
           dz = zle(i,j,k-1) - zle(i,j,k)

           A(i,j,k)         = A(i,j,k)/A_star_sum(i,j)
           A_star2_int(i,j) = A_star2_int(i,j) + A(i,j,k)**2.*dz/rho(i,j,k)
        end do
     end if
  end do
  end do

  ! Find magnitude and height of maximum plume velocity and compute Mu0 accordingly
  wu0(i,j) = 0.
  do iter = 1,1

     do j = 1,JM
     do i = 1,IM
        ! Initialize plume
        Mu(i,j)   = A(i,j,LM)*( zle(i,j,LM-1) - zle(i,j,LM) )
        wu2(i,j)  = wu0(i,j)**2.
        thlu(i,j) = thl(i,j,LM)
        qtu(i,j)  = qt(i,j,LM)
        thvu(i,j) = thv(i,j,LM)

        test_flag(i,j) = conv_flag(i,j)
     end do
     end do

     f(:,:,:) = 1.
     do k = LM-1,1,-1
        km1 = k - 1

        do j = 1,JM
        do i = 1,IM
           if ( test_flag(i,j) ) then
              dz = zle(i,j,km1) - zle(i,j,k)
              
              Mu_next = Mu(i,j) + A(i,j,k)*dz

              f(i,j,k)      = Mu(i,j)/Mu_next
              ent_sl(i,j,k) = -log(f(i,j,k))/dz
              
              B = goth00*( thvu(i,j) - thv(i,j,k) )
              
              thlu_next = f(i,j,k)*thlu(i,j) + ( 1. - f(i,j,k) )*thl(i,j,k)
              qtu_next  = f(i,j,k)*qtu(i,j)  + ( 1. - f(i,j,k) )*qt(i,j,k)
              wu2_next  = f(i,j,k)**2.*wu2(i,j) + 2.*dz*B

              if ( debug_flag /= 0 ) then
!                 write(*,'(I4,2F7.2)') k, B, sqrt(wu2(i,j)), sqrt(wu2_next)
!                 write(*,*) k, B, sqrt(wu2(i,j)), sqrt(wu2_next)
              end if
              
              ! call condensation_edmf(qtu_next, thlu_next, ple(i,j,km1), thvu_next, qcu, wf, ice_ramp)           
              thvu_next = thlu_next*( 1. + mapl_vireps*qtu_next )

              if ( B <= 0. ) then
                 zi(i,j)        = zle(i,j,k)
                 Mu0(i,j)       = sqrt( wu2(i,j) )/( 2.*zi(i,j)*A_star2_int(i,j) )
                 wu0(i,j)       = Mu0(i,j)*A(i,j,LM)*( zle(i,j,LM-1) - zle(i,j,LM) )/( au0*rhoe(i,j,LM-1) )
                 test_flag(i,j) = .false.

                 if ( debug_flag /= 0 ) then
!                    write(*,*) '*', iter, wu0(i,j) 
                 end if
              else
                 Mu(i,j)   = Mu_next
                 wu2(i,j)  = wu2_next
                 thlu(i,j) = thlu_next
                 qtu(i,j)  = qtu_next
                 thvu(i,j) = thvu_next
              end if
           end if
        end do
        end do
     end do
  end do

  do k = 1,LM
     do j = 1,JM
     do i = 1,IM
        A(i,j,k) = Mu0(i,j)*A(i,j,k)
     end do
     end do
  end do

end subroutine A_star_closure

!
! condensation_edmf
!
subroutine condensation_edmf(QT, THL, P, THV, QC, wf, ice_ramp)
!
! zero or one condensation for edmf: calculates THV and QC
!
use GEOS_UtilsMod, only : GEOS_Qsat

real,intent(in) :: QT,THL,P
real,intent(in) :: ice_ramp
real,intent(out):: THV,QC,wf

integer :: niter, i
real    :: diff, exn, t, qs, qcold
 
! max number of iterations
niter=50
! minimum difference
diff=2.e-5

EXN=(P/mapl_p00)**mapl_kappa
QC=0. 
T=EXN*THL

do i=1,NITER
  T=EXN*THL+get_alhl(T,ice_ramp)/mapl_cp*QC
! qsat, p is in pascal
  QS=geos_qsat(T,P,pascals=.true.,ramp=ice_ramp)
  QCOLD=QC
  QC=max(0.5*QC+0.5*(QT-QS),0.)
if (abs(QC-QCOLD)<Diff) exit
enddo

T=EXN*THL+get_alhl(T,ice_ramp)/mapl_cp*QC
QS=geos_qsat(T,P,pascals=.true.,ramp=ice_ramp)
QC=max(QT-QS,0.)
THV=(THL+get_alhl(T,ice_ramp)/mapl_cp*QC/EXN)*(1.+(mapl_vireps)*(QT-QC)-QC)
wf=water_f(T,ice_ramp)

end subroutine condensation_edmf


subroutine condensation_edmfA(THV,QT,P,THL,QL,QI,ice_ramp)
!
! zero or one condensation for edmf: calculates QL,QI from THV and QT
!

use GEOS_UtilsMod, only : GEOS_Qsat

real,intent(in) :: THV,QT,P
real,intent(in) :: ice_ramp
real,intent(out):: THL,QL,QI


integer :: niter,i
real :: diff,exn,t,qs,qcold,wf,qc
 
! max number of iterations
niter=50
! minimum difference
diff=2.e-5

EXN=(P/mapl_p00)**mapl_kappa
QC=0. 

T=EXN*THL

do i=1,NITER
   T=EXN*THV/(1.+mapl_vireps*(QT-QC)-QC)
   QS=geos_qsat(T,P,pascals=.true.,ramp=ice_ramp)
   QCOLD=QC
   QC=max(0.5*QC+0.5*(QT-QS),0.)
if (abs(QC-QCOLD)<Diff) exit
enddo

 THL=(T-get_alhl(T,ice_ramp)/mapl_cp*QC)/EXN
 wf=water_f(T,ice_ramp)
 QL=QC*wf
 QI=QC*(1.-wf) 

end subroutine condensation_edmfA


function  get_alhl3(T,IM,JM,LM,iceramp)

real,dimension(IM,JM,LM) ::  T,get_alhl3
real :: iceramp
integer :: IM,JM,LM
integer :: I,J,L

do i=1,im
  do j=1,jm
    do l=1,lm
        get_alhl3(i,j,l)=get_alhl(T(i,j,l),iceramp)
    enddo
  enddo
enddo

end function get_alhl3

function get_alhl(T,iceramp)
   real :: T,get_alhl,iceramp,wf
    wf=water_f(T,iceramp)
    get_alhl=wf*mapl_alhl+(1.-wf)*mapl_alhs
end function get_alhl


function water_f(T,iceramp)
!
! computes water fraction
!
real ::T,iceramp,water_f,Tw
real :: Tmax,Tmin

  Tmax=0.
  Tmin=Tmax-abs(iceramp)
  Tw=T-273.16

! water fraction
  IF (Tw>Tmax) THEN
    water_f=1.
  ELSE IF (Tw<Tmin) THEN
    water_f=0.
  ELSE
    water_f=(Tw-Tmin)/(Tmax-Tmin);
  END IF

end function water_f


subroutine Poisson(istart,iend,jstart,jend,mu,POI,seed)

integer, intent(in) :: istart,iend,jstart,jend 
real,dimension(istart:iend,jstart:jend),intent(in) :: MU
integer, dimension(istart:iend,jstart:jend), intent(out) :: POI
integer :: sed_len
integer,dimension(2),  intent(in) :: seed

integer :: seed_len
integer :: i,j,idum,p
integer,allocatable :: theseed(:)

call random_seed(SIZE=seed_len)
allocate(theseed(seed_len))

theseed(1:2)=seed
! Gfortran uses longer seeds, so fill the rest with zero
if (seed_len > 2) theseed(3:) = seed(2)
 
 
call random_seed(put=theseed)


do i=istart,iend
 do j=jstart,jend
    poi(i,j)=poidev(mu(i,j),idum)
enddo
 enddo

end subroutine Poisson

subroutine Poisson_new(mu, poi, seed)

  real, intent(in)                  :: mu
  integer, dimension(2), intent(in) :: seed
  
  integer, intent(out) :: poi
  
  integer                            :: sed_len, seed_len, idum, p 
  integer, dimension(:), allocatable :: theseed

  call random_seed(SIZE=seed_len)
  allocate(theseed(seed_len))

  theseed(1:2) = seed
  ! Gfortran uses longer seeds, so fill the rest with zero
  if (seed_len > 2) theseed(3:) = seed(2)
  
  call random_seed(put=theseed)

  poi = poidev(mu,idum)

end subroutine Poisson_new



      FUNCTION poidev(xm,idum)
      INTEGER idum
      REAL poidev,xm,PI
      PARAMETER (PI=3.141592654)
!CU    USES gammln,ran1
      REAL alxm,em,g,oldm,sq,t,y
      SAVE alxm,g,oldm,sq
      DATA oldm /-1./
      if (xm.lt.12.)then
        if (xm.ne.oldm) then
          oldm=xm
          g=exp(-xm)
        endif
        em=-1
        t=1.
2       em=em+1.
        t=t*ran1(idum)
        if (t.gt.g) goto 2
      else
        if (xm.ne.oldm) then
          oldm=xm
          sq=sqrt(2.*xm)
          alxm=log(xm)
          g=xm*alxm-gammln(xm+1.)
        endif
1       y=tan(PI*ran1(idum))
        em=sq*y+xm
        if (em.lt.0.) goto 1
        em=int(em)
        t=0.9*(1.+y**2)*exp(em*alxm-gammln(em+1.)-g)
        if (ran1(idum).gt.t) goto 1
      endif
      poidev=em
      return
      END FUNCTION poidev

      FUNCTION gammln(xx)
      REAL gammln,xx
      INTEGER j
      DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0, &
     24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2, &
     -.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do 11 j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
11    continue
      gammln=tmp+log(stp*ser/x)
      return
      END FUNCTION gammln

      function ran1(idum)
      real ran1
      integer idum

      call random_number(ran1)
      end function ran1


end module edmf_mod

