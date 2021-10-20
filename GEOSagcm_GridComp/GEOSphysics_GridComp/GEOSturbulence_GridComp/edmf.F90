!#define EDMF_DIAG 1

module edmf_mod

use MAPL_ConstantsMod, only: mapl_grav, mapl_cp, mapl_alhl, mapl_p00, mapl_vireps, mapl_alhs, mapl_kappa, mapl_rgas
!use MAPL_Mod, only:          mapl_undef

implicit none

real, parameter ::     &
     WSTARmin = 1.e-3, &
     zpblmin  = 100.,  &
     onethird = 1./3., &
     r        = 2.

contains

!
! edmf
!
subroutine run_edmf(IM, JM, LM, nup, mup, iras, jras, &                                 ! in
                    discrete_type, implicit_flag, stochastic, plume_type, debug_flag, & ! in
                    dt, zl, zle, ple, rho, rhoe, exf, &                                 ! in
                    u, v, thl, qt, qv, ql, qi, thv, &                                   ! in         
                    ustar, sh, evap, frland, zpbl_in, ice_ramp, &                       ! in
                    pwmin, pwmax, AlphaW, AlphaQT, AlphaTH, &                           ! in
                    ET, L0, ENT0, EDfac, EntWFac, Wa, Wb, &                             ! in
                    au0, cth1, cth2, rho_qb, &                                          ! in
                    edmfdrya, edmfmoista, &                                             ! out
                    edmfdryw, edmfmoistw, &                                             ! out
                    edmfdryqt, edmfmoistqt, &                                           ! out
                    edmfdrythl, edmfmoistthl, &                                         ! out
                    edmfdryu, edmfmoistu,  &                                            ! out
                    edmfdryv, edmfmoistv,  &                                            ! out
                    edmfmoistqc, &                                                      ! out
                    ae, awu, awv, aw, aws, awqv, awql, awqi, &                          ! out (for solver)
                    whl_mf, wqt_mf, wthv_mf, &                                          ! out (for inconsistent partitioning)
                    buoyf, mfw2, mfw3, mfqt3, mfwqt, mfqt2, mfhl2, mfqthl, mfwhl, &     ! out (for SHOC)
#ifdef EDMF_DIAG
                    w_plume1,w_plume2,w_plume3,w_plume4,w_plume5, &
                    w_plume6,w_plume7,w_plume8,w_plume9,w_plume10, &
                    qt_plume1,qt_plume2,qt_plume3,qt_plume4,qt_plume5, &
                    qt_plume6,qt_plume7,qt_plume8,qt_plume9,qt_plume10, &
                    thl_plume1,thl_plume2,thl_plume3,thl_plume4,thl_plume5, &
                    thl_plume6,thl_plume7,thl_plume8,thl_plume9,thl_plume10, &
#endif
                    wu0, dthv0)                                                              ! out
  
  ! Inputs
  integer, intent(in)                     :: IM, JM, LM, nup, mup, discrete_type, implicit_flag, plume_type, debug_flag, ET
  real, intent(in)                        :: stochastic
  integer, dimension(IM,JM), intent(in)   :: iras, jras
  real, dimension(IM,JM,LM), intent(in)   :: u, v, thl, qt, thv, qv, ql, qi, zl, exf, rho
  real, dimension(IM,JM,0:LM), intent(in) :: zle, ple, rhoe
  real, dimension(IM,JM), intent(in)      :: ustar, sh, evap, frland, zpbl_in, L0
  real, intent(in)                        :: ice_ramp, dt, EntWFac, ENT0, Wa, Wb, au0, cth1, cth2, rho_qb, pwmin, pwmax, &
                                             AlphaW, AlphaQT, AlphaTH, EDfac
 
  ! Outputs
  real, dimension(IM,JM,0:LM), intent(out) :: edmfdrya, edmfmoista, edmfdryw, edmfmoistw, &
                                              edmfdryqt, edmfmoistqt, edmfdrythl, edmfmoistthl, &
                                              edmfdryu, edmfmoistu, edmfdryv, edmfmoistv, &
                                              edmfmoistqc, &
                                              ae, aw, aws, awqv, awql, awqi, awu, awv, &
                                              whl_mf, wqt_mf, wthv_mf
#ifdef EDMF_DIAG
                                              w_plume1,w_plume2,w_plume3,w_plume4,w_plume5, &
                                              w_plume6,w_plume7,w_plume8,w_plume9,w_plume10, &
                                              qt_plume1,qt_plume2,qt_plume3,qt_plume4,qt_plume5, &
                                              qt_plume6,qt_plume7,qt_plume8,qt_plume9,qt_plume10, &
                                              thl_plume1,thl_plume2,thl_plume3,thl_plume4,thl_plume5, &
                                              thl_plume6,thl_plume7,thl_plume8,thl_plume9,thl_plume10, &
#endif
  real, dimension(IM,JM,LM), intent(out) :: buoyf, mfw2, mfw3, mfqt3, mfqt2, mfwqt, mfhl2, mfqthl, mfwhl                                            
  real, dimension(IM,JM), intent(out)    :: wu0, dthv0

  ! Local variables
  integer :: i, j, k, km1, kp1, iup, jup, iup1, iup2, kbot, numup
  real    :: wthv, qstar, thstar, sigmaw, sigmaqt, sigmath, sigmaQT_cond, z0, wp, &
             B, QTn, THLn, THVn, QCn, Un, Vn, Wn, Wn2, Mn, Mn_test, au_test, EntEXP, EntEXPU, EntW, wf, &
             stmp, QTsrfF, THVsrfF, mft, mfthvt, dzle, idzle, ifac, test, mft_work, mfthvt_work, &
             thlu_full, work1, work2, exfh, dsdz, dqdz, wmin, wmax, wlv, wtv, &
             thv_high, thv_low, thvmin, thvmax, qt_high, qt_low, qtmin, qtmax, upa_mean, upqt_mean, upthv_mean, &
             thvu0, qtu0, foo, eps, eps_org, eps_turb, pdf_fac1, pdf_fac2, pdf_fac3, pdf_fac4, Mdet, E_org, E_turb, D_turb, &
             udet_test, vdet_test, wdet_test, E_plume, dp, An

  logical, dimension(IM,JM) :: work_flag
  real, dimension(IM,JM)    :: zpbl, zi, Phi, au_fac, wstar, dqt0, factor

  real, dimension(IM,JM,LM) :: thlu, qtu, qlu, edmfmoistql, A, E, D, f_thermal, &
                               au_full, hlu_full, qtu_full, acu_full, Tu_full, qlu_full

  real, dimension(IM,JM,0:LM) :: aw2, ahl2, aqt2, aw3, aqt3, aqthl, &
                                 ui, vi, thli, qti, qvi, qli, qii, thvi, &
                                 au, Mu, uu, vu, wu

  integer, dimension(2) :: seedmf, the_seed

  real, dimension(mup*nup,IM,JM) :: upm, upw, upthl, upqt, upql, upqi, upa, upu, upv, upthv, D_frac

  integer, dimension(mup*nup,IM,JM,LM) :: enti
  real, dimension(mup*nup,IM,JM,LM)    :: entf, ent

  real, dimension(mup*nup) :: upa_out, upac_out

  logical, dimension(mup*nup,IM,JM) :: active_updraft

  ! Interpolate EDMF profiles to half levels 
  call interp_edmf(IM, JM, LM, &                           ! in
                   discrete_type, implicit_flag, &         ! in
                   zle, zl, &                              ! in
                   u, v, thl, qt, qv, ql, qi, thv, &       ! in
                   ui, vi, thli, qti, qvi, qli, qii, thvi) ! out

  numup = mup*nup

  if ( plume_type == 0 ) then
     kbot = LM
  else
     kbot = LM - 1
  end if

  !
  ! Initialize arrays
  !

  wstar(:,:) = 0.

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
  uu   = ui
  vu   = vi
  wu   = 0.
  Mu   = 0.
  E    = 0.
  D    = 0.

#ifdef EDMF_DIAG
  w_plume1  = 0.
  w_plume2  = 0.
  w_plume3  = 0.
  w_plume4  = 0.
  w_plume5  = 0.
  w_plume6  = 0.
  w_plume7  = 0.
  w_plume8  = 0.
  w_plume9  = 0.
  w_plume10 = 0.

  qt_plume1  = 0.
  qt_plume2  = 0.
  qt_plume3  = 0.
  qt_plume4  = 0.
  qt_plume5  = 0.
  qt_plume6  = 0.
  qt_plume7  = 0.
  qt_plume8  = 0.
  qt_plume9  = 0.
  qt_plume10 = 0.

  thl_plume1  = 0.
  thl_plume2  = 0.
  thl_plume3  = 0.
  thl_plume4  = 0.
  thl_plume5  = 0.
  thl_plume6  = 0.
  thl_plume7  = 0.
  thl_plume8  = 0.
  thl_plume9  = 0.
  thl_plume10 = 0.
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
     call A_star_closure(IM, JM, LM, &                       ! in
                         wb, au0, debug_flag, &              ! in
                         zle, zl, rho, rhoe, thl, qt, thv, & ! in
                         A, f_thermal, zi, Phi)              ! out
     E = A
  end if

  !
  ! Initialize updrafts
  !
  do j = 1,JM
  do i = 1,IM
     ent = 0.

     if ( plume_type == 0 ) then
        zpbl(i,j) = max( zpbl_in(i,j), zpblmin )
     end if

     ! Compute surface THV flux
     wthv = ( sh(i,j)/mapl_cp + mapl_vireps*thv(i,j,LM)*evap(i,j) )/rhoe(i,j,LM)

     ! Test for surface-driven convection
!     if ( wthv > 0. .and. thv(i,j,LM-1) < thv(i,j,LM) ) then
     if ( wthv > 0. ) then
        ! Initialize plumes at surface
        if ( plume_type == 0 ) then ! JPL scheme
           wstar(i,j) = max( wstarmin, ((mapl_grav/thv(i,j,LM))*wthv*zpbl(i,j))**onethird )

           qstar  = evap(i,j)/(rhoe(i,j,LM)*wstar(i,j))
           thstar = wthv/wstar(i,j)
        
           sigmaW  = AlphaW*wstar(i,j)
           sigmaQT = AlphaQT*max( qstar, 0. )
           sigmaTH = AlphaTH*max( thstar, 0. )

           wmin = sigmaW*pwmin ! vertical velocity of least energetic updraft
           wmax = sigmaW*pwmax ! "        "        "  most  "         "
        
!           pdf_fac1 = 1./( sqrt(2.)*sigmaW )

           QTsrfF  = 0.
           THVsrfF = 0.
           do iup = 1,numup
              active_updraft(iup,i,j) = .true.         
              
!              wlv = wmin + ( wmax - wmin )*real(iup-1)/real(numup)
!              wtv = wmin + ( wmax - wmin )*real(iup)/real(numup) 
              wlv = wmin + ( wmax - wmin )/(real(numup))*( real(iup) - 1. )
              wtv = wmin + ( wmax - wmin )/(real(numup))*real(iup)

              upw(iup,i,j)   = 0.5*( wlv + wtv ) 
              upa(iup,i,j)   = 0.5*( erf(wtv/(sqrt(2.)*sigmaW)) - erf(wlv/(sqrt(2.)*sigmaW)) )
!              upa(iup,i,j)   = 0.5*( erf(wtv*pdf_fac1) - erf(wlv*pdf_fac1) )
              upu(iup,i,j)   = ui(i,j,kbot)
              upv(iup,i,j)   = vi(i,j,kbot)
              upqt(iup,i,j)  = qti(i,j,kbot)  + 0.32*upw(iup,i,j)*sigmaQT/sigmaW ! 0.32~=0.58*0.55 (Stull 1988)
              upthv(iup,i,j) = thvi(i,j,kbot) + 0.58*upw(iup,i,j)*sigmaTH/sigmaW ! Stull 1988

              ! For stability make sure that the surface mass-fluxes are not more than their 
              ! values computed from the surface scheme
              if ( kbot == LM ) then
                 QTsrfF  = QTsrfF  + upa(iup,i,j)*upw(iup,i,j)*( upqt(iup,i,j)  - qti(i,j,kbot) )
                 THVsrfF = THVsrfF + upa(iup,i,j)*upw(iup,i,j)*( upthv(iup,i,j) - thvi(i,j,kbot) )
              end if
           end do ! iup = 1,num

           ! Change surface THV and QT so that the fluxes from the mass flux equal prescribed values
           if ( kbot == LM ) then
              if ( THVsrfF > wthv ) then
                 upthv(:,i,j) = ( upthv(:,i,j) - thvi(i,j,kbot) )*wthv/THVsrfF + thvi(i,j,kbot)
              end if

              ! We do not need to worry about the negative values as they should not exist
              if ( QTsrfF > evap(i,j)/rhoe(i,j,LM) .and. evap(i,j) > 0. )  then
                 upqt(:,i,j) = ( upqt(:,i,j) - qti(i,j,kbot) )*evap(i,j)/(rhoe(i,j,LM)*QTsrfF) + qti(i,j,kbot)
              end if
           end if
        elseif ( plume_type == 1 ) then ! Thermal plume
           work1  = A(i,j,LM)*( zle(i,j,LM-1) - zle(i,j,LM) )/(  A(i,j,LM)*( zle(i,j,LM-1) - zle(i,j,LM) ) &
                                                              + A(i,j,LM-1)*( zle(i,j,LM-2) - zle(i,j,LM-1) ) )
           work2 = rhoe(i,j,LM-2)/rhoe(i,j,LM-1) 

           au_fac(i,j) = ( work2/work1 )**2. - work1**2.

           ! Get mean vertical velocity of updraft distribution 
           wu0(i,j) = A(i,j,LM)*( zle(i,j,LM-1) - zle(i,j,LM) )/( au0*rhoe(i,j,LM-1) )

           ! Get CBL scales
           wstar(i,j) = ( (mapl_grav/thv(i,j,k))*wthv*zi(i,j) )**onethird
           thstar     = wthv/wstar(i,j) 

           !
           ! Specify updraft THV distribution
           !

           dthv0(i,j) = cth1*wu0(i,j)*thstar/wstar(i,j)

           thvu0   = thvi(i,j,LM-1) + dthv0(i,j)
           sigmaTH = cth2*dthv0(i,j)

           thvmin = thvi(i,j,LM-1)
           thvmax = thvu0 + 3.*sigmaTH

           pdf_fac1 = 1./( sqrt(2.)*sigmaTH )

           pdf_fac2 = au0/( 0.5*( erf(( thvmax - thvu0 )*pdf_fac1) - erf(( thvmin - thvu0 )*pdf_fac1) ) )

           !
           ! Specify updraft QT distribution
           !

           dqt0(i,j) = dthv0(i,j)*evap(i,j)/(rhoe(i,j,LM-1)*wthv) ! keep Bowen ratio same as surface

           qtu0 = qti(i,j,LM-1) + dqt0(i,j)

           sigmaQT = sigmaTH/rho_qb * evap(i,j)/(rhoe(i,j,LM-1)*wthv) 

           sigmaQT_cond = sqrt( 1. - rho_qb**2. )*sigmaQT

           !
           ! Initialize plumes
           !

           upa_out(:) = 0.

           do jup = 1,nup
              iup1 = ( jup - 1 )*mup + 1
              iup2 = jup*mup

              thv_low  = thvmin + ( thvmax - thvmin )*real(jup-1)/real(nup)
              thv_high = thvmin + ( thvmax - thvmin )*real(jup)/real(nup)

              upa_mean   = 0.5*pdf_fac2*( erf(( thv_high - thvu0 )*pdf_fac1) - erf(( thv_low - thvu0 )*pdf_fac1) )                 
              upthv_mean = 0.5*( thv_low + thv_high )

              upqt_mean = qti(i,j,LM-1) + ( upthv_mean - thvi(i,j,LM-1) )*evap(i,j)/(rhoe(i,j,LM)*wthv)
!              upqt_mean = qtu0 + sigmaQT/sigmaTH*rho_qb*( upthv_mean - thvu0 )
              
              qtmin = upqt_mean - 3.*sigmaQT_cond
              qtmax = upqt_mean + 3.*sigmaQT_cond

              pdf_fac3 = 1./( sqrt(2.)*sigmaQT_cond )

              pdf_fac4 = upa_mean/( 0.5*( erf(( qtmax - upqt_mean )*pdf_fac3) - erf(( qtmin - upqt_mean )*pdf_fac3) ) )

              upa_out(jup) = upa_mean

              do iup = iup1,iup2
                 active_updraft(iup,i,j) = .true.         

                 upthv(iup,i,j) = upthv_mean
                 upu(iup,i,j)   = ui(i,j,LM-1)
                 upv(iup,i,j)   = vi(i,j,LM-1)                

                 upw(iup,i,j) = wu0(i,j)
!                 if ( au_fac(i,j) > 0. ) then
!                    upw(iup,i,j) = sqrt( 2.*wb*(mapl_grav/thv(i,j,LM-1))*( zle(i,j,LM-2) - zle(i,j,LM-1) )*( upthv(iup,i,j) - thvi(i,j,LM-1) )/au_fac(i,j) )
!                 else
!                    upw(iup,i,j) = wu0(i,j)
!                 end if

                 qt_low  = qtmin + ( qtmax - qtmin )*real(iup-iup1)/real(mup)
                 qt_high = qtmin + ( qtmax - qtmin )*real(iup-iup1+1)/real(mup)
                 
                 upa(iup,i,j)  = 0.5*pdf_fac4*( erf(( qt_high - upqt_mean )*pdf_fac3) - erf(( qt_low - upqt_mean )*pdf_fac3) )
                 upqt(iup,i,j) = 0.5*( qt_low + qt_high )

                 upm(iup,i,j) = rhoe(i,j,LM-1)*upa(iup,i,j)*upw(iup,i,j)
              end do
           end do

           ! Adjust vertical velocity so that mass flux is correct at half level nearest surface
!           upw(:,i,j) = upw(:,i,j) * sum(upa(:,i,j))*wu0(i,j)/sum(upa(:,i,j)*upw(:,i,j))
!           upm(:,i,j) = rhoe(:,i,j)*upa(:,i,j)*upw(:,i,j)
        end if

        ! Compute condensation and thl, ql, qi
        do iup = 1,numup
           call condensation_edmfA(upthv(iup,i,j), upqt(iup,i,j), ple(i,j,kbot), &
                                   upthl(iup,i,j), upql(iup,i,j), upqi(iup,i,j), ice_ramp)
        end do
     else
        active_updraft(:,i,j) = .false.
     end if ! wthv > 0. .and. thv(i,j,LM-1) < thv(i,j,LM)
  end do ! i = 1,IM
  end do ! j = 1,JM

  if ( debug_flag /= 0 ) then
     write(*,'(I4,12F7.2)') LM-1, sum(100.*upa_out(:)), 0., 100.*upa_out(1:10)
  end if

  ! Check that mass flux does not exceed 2x of layer mass at any level  
  ! If it does, rescale updraft area.
  ! See discussion in Beljaars et al 2018 [ECMWF Tech Memo]
!  factor(:,:) = 1.0
!  do k = 1,LM
!     km1 = k - 1

!     do j = 1,JM
!     do i = 1,IM
!        mf = sum( rhoe(i,j,k)*upa(i,j,k)*upw(i,j,k) )
!        dp = ple(i,j,km1) - ple(i,j,k)
!        if ( mf > dp/(mapl_grav*dt) ) then
!           factor(i,j) = min( factor(i,j), dp)/(mf*mapl_grav*dt) )
!        end if
!     end do
!     end do
!  end do
!  upa = factor*upa

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
        
        exfh = ( ple(i,j,k)/mapl_p00 )**mapl_kappa

        mfthvt = 0.
        mft    = 0.

        Mdet      = 0.
        udet_test = 0.
        vdet_test = 0.
        wdet_test = 0.

        upa_out(:)  = 0.
        upac_out(:) = 0.

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
              if ( mod(iup,mup) == (mup+1)/2 ) then
                 jup = iup/mup + 1
              else
                 jup = -1
              end if

              if ( jup == 1 ) then
                 w_plume1(i,j,k)   = upw(iup,i,j)
                 qt_plume1(i,j,k)  = upqt(iup,i,j) - qti(i,j,k)
                 thl_plume1(i,j,k) = upthv(iup,i,j) - thvi(i,j,k)
!                 thl_plume1(i,j,k) = ( upthl(iup,i,j) - thli(i,j,k) )*exfh
              elseif ( jup == 2 ) then
                 w_plume2(i,j,k)   = upw(iup,i,j)
                 qt_plume2(i,j,k)  = upqt(iup,i,j) - qti(i,j,k)
                 thl_plume2(i,j,k) = upthv(iup,i,j) - thvi(i,j,k)
!                 thl_plume2(i,j,k) = ( upthl(iup,i,j) - thli(i,j,k) )*exfh
              elseif ( jup == 3 ) then
                 w_plume3(i,j,k)   = upw(iup,i,j)
                 qt_plume3(i,j,k)  = upqt(iup,i,j) - qti(i,j,k)
                 thl_plume3(i,j,k) = upthv(iup,i,j) - thvi(i,j,k)
!                 thl_plume3(i,j,k) = ( upthl(iup,i,j) - thli(i,j,k) )*exfh
              elseif ( jup == 4 ) then
                 w_plume4(i,j,k)   = upw(iup,i,j)
                 qt_plume4(i,j,k)  = upqt(iup,i,j) - qti(i,j,k)
                 thl_plume4(i,j,k) = upthv(iup,i,j) - thvi(i,j,k)
!                 thl_plume4(i,j,k) = ( upthl(iup,i,j) - thli(i,j,k) )*exfh
              elseif ( jup == 5 ) then
                 w_plume5(i,j,k)   = upw(iup,i,j)
                 qt_plume5(i,j,k)  = upqt(iup,i,j) - qti(i,j,k)
                 thl_plume5(i,j,k) = upthv(iup,i,j) - thvi(i,j,k)
!                 thl_plume5(i,j,k) = ( upthl(iup,i,j) - thli(i,j,k) )*exfh
              elseif ( jup == 6 ) then
                 w_plume6(i,j,k)   = upw(iup,i,j)
                 qt_plume6(i,j,k)  = upqt(iup,i,j) - qti(i,j,k)
                 thl_plume6(i,j,k) = upthv(iup,i,j) - thvi(i,j,k)
!                 thl_plume6(i,j,k) = ( upthl(iup,i,j) - thli(i,j,k) )*exfh
              elseif ( jup == 7 ) then
                 w_plume7(i,j,k)  = upw(iup,i,j)
                 qt_plume7(i,j,k)  = upqt(iup,i,j) - qti(i,j,k)
                 thl_plume7(i,j,k) = upthv(iup,i,j) - thvi(i,j,k)
!                 thl_plume7(i,j,k) = ( upthl(iup,i,j) - thli(i,j,k) )*exfh
              elseif ( jup == 8 ) then
                 w_plume8(i,j,k)   = upw(iup,i,j)
                 qt_plume8(i,j,k)  = upqt(iup,i,j) - qti(i,j,k)
                 thl_plume8(i,j,k) = upthv(iup,i,j) - thvi(i,j,k)
!                 thl_plume8(i,j,k) = ( upthl(iup,i,j) - thli(i,j,k) )*exfh
              elseif ( jup == 9 ) then
                 w_plume9(i,j,k)   = upw(iup,i,j)
                 qt_plume9(i,j,k)  = upqt(iup,i,j) - qti(i,j,k)
                 thl_plume9(i,j,k) = upthv(iup,i,j) - thvi(i,j,k)
!                 thl_plume9(i,j,k) = ( upthl(iup,i,j) - thli(i,j,k) )*exfh
              elseif ( jup == 10 ) then
                 w_plume10(i,j,k)   = upw(iup,i,j)
                 qt_plume10(i,j,k)  = upqt(iup,i,j) - qti(i,j,k)
                 thl_plume10(i,j,k) = upthv(iup,i,j) - thvi(i,j,k)
!                 thl_plume10(i,j,k) = ( upthl(iup,i,j) - thli(i,j,k) )*exfh
              end if
#endif

              !
              !
              !
              
              ! Compute fractional entrainment rate
              if ( L0(i,j) > 0. ) then
                 call Poisson(1, 1, 1, 1, entf(iup,i,j,k), enti(iup,i,j,k), the_seed)

                 ent(iup,i,j,k) = ( 1. - stochastic )*ent0/L0(i,j) + stochastic*real( enti(iup,i,j,k) )*ent0/( zle(i,j,km1) - zle(i,j,k) )
              else
                 ! negative L0 means 0 entrainment  
                 ent(:,i,j,:) = 0. ! check
              end if

              ! Double entrainment over land
              ent(iup,i,j,k) = ( 1. + frland(i,j) )*ent(iup,i,j,k)

              !
!              if ( k <= LM-1 .and. k >= 2 .and. thv(i,j,k) < thv(i,j,k+1) .and. thv(i,j,k) > thv(i,j,k-1) ) then
!                 ent(iup,i,j,k) = ent(iup,i,j,k) + 5.*ent0/L0(i,j)
!              end if

              ! Integrate plume budgets across one vertical step
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
                 call condensation_edmf(QTn, THLn, ple(i,j,km1), THVn, QCn, wf, ice_ramp)

                 ! Buoyancy
!                 B = (mapl_grav/thv(i,j,k))*( 0.5*( THVn + upthv(iup,i,j) ) - thv(i,j,k) )
                 B = mapl_grav*( 0.5*( THVn + upthv(iup,i,j) )/thv(i,j,k) - 1. )

                 ! Vertical velocity
                 WP = Wb*ent(iup,i,j,k)
                 if (WP == 0.) then
                    Wn2 = upw(iup,i,j)**2. + 2.*Wa*B*dzle
                 else
                    EntW = exp(-2.*WP*dzle)
                    Wn2  = EntW*upw(iup,i,j)**2. + Wa*B/WP*( 1. - EntW )
                 end if

              elseif ( plume_type == 1 ) then ! thermal plume
                 ! Compute entrainment rate 
                 eps_org  = ( 1./f_thermal(i,j,k) - 1. )/dzle
                 eps_turb = max( 0., 1./L0(i,j) - eps_org )

                 eps = eps_org + eps_turb

                 ! Integrate scalar budgets
                 THLn = ( upthl(iup,i,j) + dzle*eps*thl(i,j,k) )/( 1. + dzle*eps )
                 QTn  = ( upqt(iup,i,j) + dzle*eps*qt(i,j,k) )/( 1. + dzle*eps )
                 Un   = ( upu(iup,i,j) + dzle*eps*u(i,j,k) )/( 1. + dzle*eps )
                 Vn   = ( upv(iup,i,j) + dzle*eps*v(i,j,k) )/( 1. + dzle*eps )

                 ! Condensation
                 call condensation_edmf(QTn, THLn, ple(i,j,km1), THVn, QCn, wf, ice_ramp)

                 ! Compute buoyancy
                 B = wb*(mapl_grav/thv(i,j,k))*( 0.5*( THVn + upthv(iup,i,j) ) - thv(i,j,k) )

                 ! Integrate vertical velocity budget
                 Wn2 = upw(iup,i,j)**2./( 1./f_thermal(i,j,k) + dzle*eps_turb  )**2. + 2.*dzle*B                 

                 if ( Wn2 > 0. ) then
                    E_org  = upm(iup,i,j)*eps_org
                    E_turb = upm(iup,i,j)*eps_turb

                    D_turb = E_turb

                    E(i,j,k) = E(i,j,k) + E_turb

                    Mn_test = upm(iup,i,j) + dzle*( E_org + E_turb )

                    if ( B > 0. ) then
                       Mn = max( 0., Mn_test - dzle*D_turb )
                       if ( Mn == 0. ) then
                          Wn2 = 0.
                       end if
                    else
                       au_test = max( 0., Mn_test - dzle*D_turb )/( rhoe(i,j,km1)*sqrt(Wn2) )

                       if ( au_test > upa(iup,i,j) ) then
                          Mn = rhoe(i,j,km1)*upa(iup,i,j)*sqrt(Wn2)

                          D_turb = D_turb + rhoe(i,j,km1)*( au_test - upa(iup,i,j) )*sqrt(Wn2)/dzle
                       else
                          Mn = max( 0., Mn_test - dzle*D_turb )
                       end if
                    end if
                    
                    Mdet      = Mdet      + dzle*D_turb
                    udet_test = udet_test + dzle*D_turb*Un
                    vdet_test = vdet_test + dzle*D_turb*Vn
                    wdet_test = wdet_test + dzle*D_turb*sqrt(Wn2)
                 end if
              end if
              
              ! Intermediate quantities 
              if ( plume_type == 1 ) then
                 An = Mn/( rhoe(i,j,km1)*upw(iup,i,j) )
              else
                 An = upa(iup,i,j)
              end if

              ! Update plume
              if ( Wn2 > 0. ) then
                 mfthvt = 0.5*( upa(iup,i,j)*upw(iup,i,j)*upthv(iup,i,j) + An*sqrt(Wn2)*THVn )
                 mft    = 0.5*( upa(iup,i,j)*upw(iup,i,j) + An*sqrt(Wn2) )

                 buoyf(i,j,k) = buoyf(i,j,k) + ( mfthvt - mft*thv(i,j,k) )*exf(i,j,k)

                 ! Update plume properties
                 upm(iup,i,j)   = Mn
                 upa(iup,i,j)   = An
                 upw(iup,i,j)   = sqrt(Wn2)
                 upthv(iup,i,j) = THVn
                 upthl(iup,i,j) = THLn
                 upqt(iup,i,j)  = QTn
                 upql(iup,i,j)  = QCn*wf
                 upqi(iup,i,j)  = QCn*( 1. - wf )
                 upu(iup,i,j)   = Un
                 upv(iup,i,j)   = Vn

!                 if ( plume_type == 1 ) then
!                    upm(iup,i,j) = Mn
!                    upa(iup,i,j) = upm(iup,i,j)/( rhoe(i,j,km1)*upw(iup,i,j) )
!                 end if

                 !
!                 mfthvt = mfthvt + 0.5*( upa(iup,i,j)*upw(iup,i,j)*upthv(iup,i,j) + mfthvt_work )
!                 mft    = mft    + 0.5*(                upa(iup,i,j)*upw(iup,i,j) + mft_work )

                 ! Save area fractions for text output
                 jup = iup/mup + 1
                 upa_out(jup) = upa_out(jup) + upa(iup,i,j)
                 if ( upql(iup,i,j) > 0. ) then
                    upac_out(jup) = upac_out(jup) + upa(iup,i,j)
                 end if
              else
                 mfthvt = 0.5*( upa(iup,i,j)*upw(iup,i,j)*upthv(iup,i,j) )
                 mft    = 0.5*( upa(iup,i,j)*upw(iup,i,j) )

                 buoyf(i,j,k) = buoyf(i,j,k) + ( mfthvt - mft*thv(i,j,k) )*exf(i,j,k)

                 active_updraft(iup,i,j) = .false.
              end if
           end if ! active_updraft(iup,i,j)
        end do ! iup = 1,numup

        ! Average detrainment velocity samples
!!$        if ( Mdet > 0. ) then
!!$           udet(i,j,k) = udet_test/Mdet
!!$           vdet(i,j,k) = vdet_test/Mdet
!!$           wdet(i,j,k) = wdet_test/Mdet
!!$        else
!!$           udet(i,j,k) = ui(i,j,k)
!!$           vdet(i,j,k) = vi(i,j,k)
!!$           wdet(i,j,k) = 0.
!!$        end if

        ! Outputs needed for SHOC in TURBULENCE
!        buoyf(i,j,k) = buoyf(i,j,k) + exf(i,j,k)*( mfthvt - mft*thv(i,j,k) )

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
           if ( edmfmoista(i,j,k) > 0. ) then
              uu(i,j,k)   = (  edmfdrya(i,j,k)*edmfdryu(i,j,k)  &
                             + edmfmoista(i,j,k)*edmfmoistu(i,j,k) )/au(i,j,k)
              vu(i,j,k)   = (  edmfdrya(i,j,k)*edmfdryv(i,j,k)  &
                             + edmfmoista(i,j,k)*edmfmoistv(i,j,k) )/au(i,j,k)
              wu(i,j,k)   = (  edmfdrya(i,j,k)*edmfdryw(i,j,k)  &
                             + edmfmoista(i,j,k)*edmfmoistw(i,j,k) )/au(i,j,k)
              thlu(i,j,k) = (  edmfdrya(i,j,k)*edmfdrythl(i,j,k)  &
                             + edmfmoista(i,j,k)*edmfmoistthl(i,j,k) )/au(i,j,k)
              qtu(i,j,k)  = (  edmfdrya(i,j,k)*edmfdryqt(i,j,k)  &
                             + edmfmoista(i,j,k)*edmfmoistqt(i,j,k) )/au(i,j,k)
              qlu(i,j,k)  = edmfmoistql(i,j,k)
           else
              uu(i,j,k)   = edmfdryu(i,j,k)
              vu(i,j,k)   = edmfdryv(i,j,k)
              wu(i,j,k)   = edmfdryw(i,j,k)
              thlu(i,j,k) = edmfdrythl(i,j,k)
              qtu(i,j,k)  = edmfdryqt(i,j,k)
              qlu(i,j,k)  = 0.
           end if

           if ( debug_flag /= 0 ) then
              write(*,'(I4,12F7.2)') k-1, sum(100.*upa_out(:)), sum(100.*upac_out(:)), 100.*upa_out(1:10)
           end if
        else
           wu(i,j,k)   = 0.
           thlu(i,j,k) = thli(i,j,k)
           qtu(i,j,k)  = qti(i,j,k)
           qlu(i,j,k)  = 0.
        end if
        Mu(i,j,k) = rhoe(i,j,k)*au(i,j,k)*wu(i,j,k)

        ae(i,j,k)     = ( 1. - edmfdrya(i,j,k) - edmfmoista(i,j,k) )*EDfac 
     end do ! i = 1,IM
     end do ! j = 1,JM
  end do ! k = kbot,2,-1

  ! Interpolate updraft properties to full levels
  do k = 1,LM
     km1 = k - 1

     do j = 1,JM
     do i = 1,IM
        if ( k >= 2 .and. k <= LM-1 ) then
           ! Compute detrainment rate as residual
           D(i,j,k) = E(i,j,k) - ( Mu(i,j,km1) - Mu(i,j,k) )/( zle(i,j,km1) - zle(i,j,k) )

           au_full(i,j,k)  = 0.5*( au(i,j,k)   + au(i,j,km1) )
           thlu_full       = 0.5*( thlu(i,j,k) + thlu(i,j,km1) )
           qtu_full(i,j,k) = 0.5*( qtu(i,j,k)  + qtu(i,j,km1) )

           acu_full(i,j,k) = min( au_full(i,j,k), 0.5*( edmfmoista(i,j,k) + edmfmoista(i,j,km1) ) )
           qlu_full(i,j,k) = 0.5*( qlu(i,j,k) + qlu(i,j,km1) )
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

end subroutine run_edmf

!
! interp_edmf
!
subroutine interp_edmf(IM, JM, LM, &                           ! in
                       discrete_type, implicit_flag, &         ! in
                       zle, zl, &                              ! in
                       u, v, thl, qt, qv, ql, qi, thv, &       ! in
                       ui, vi, thli, qti, qvi, qli, qii, thvi) ! out
  
  implicit none

  integer, intent(in)         :: IM, JM, LM, discrete_type, implicit_flag
  real, dimension(IM,JM,LM)   :: zl, u, v, thl, qt, qv, ql, qi, thv
  real, dimension(IM,JM,0:LM) :: zle, ui, vi, thli, qti, qvi, qli, qii, thvi

  integer :: i, j, k, kp1
  real    :: ifac

  do k = 1,LM-1
     kp1 = k + 1
           
     do j = 1,JM
     do i = 1,IM
        if ( discrete_type == 0 ) then ! centered
           ! This is temporary until I fix the interplation for implicit mass flux discretization
           if ( implicit_flag == 0 ) then
              ifac = ( zle(i,j,k) - zl(i,j,kp1) )/( zl(i,j,k) - zl(i,j,kp1) )
                    
              ui(i,j,k)   = u(i,j,kp1)   + ifac*( u(i,j,k)   - u(i,j,kp1) )
              vi(i,j,k)   = v(i,j,kp1)   + ifac*( v(i,j,k)   - v(i,j,kp1) )
              thli(i,j,k) = thl(i,j,kp1) + ifac*( thl(i,j,k) - thl(i,j,kp1) )
              qti(i,j,k)  = qt(i,j,kp1)  + ifac*( qt(i,j,k)  - qt(i,j,kp1) )
              qvi(i,j,k)  = qv(i,j,kp1)  + ifac*( qv(i,j,k)  - qv(i,j,kp1) )
              qli(i,j,k)  = ql(i,j,kp1)  + ifac*( ql(i,j,k)  - ql(i,j,kp1) )
              qii(i,j,k)  = qi(i,j,kp1)  + ifac*( qi(i,j,k)  - qi(i,j,kp1) )
              thvi(i,j,k) = thv(i,j,kp1) + ifac*( thv(i,j,k) - thv(i,j,kp1) )
           else
              ui(i,j,k)   = 0.5*( u(i,j,kp1)   + u(i,j,k) )
              vi(i,j,k)   = 0.5*( v(i,j,kp1)   + v(i,j,k) )
              thli(i,j,k) = 0.5*( thl(i,j,kp1) + thl(i,j,k) )
              qti(i,j,k)  = 0.5*( qt(i,j,kp1)  + qt(i,j,k) )
              qvi(i,j,k)  = 0.5*( qv(i,j,kp1)  + qv(i,j,k) )
              qli(i,j,k)  = 0.5*( ql(i,j,kp1)  + ql(i,j,k) )
              qii(i,j,k)  = 0.5*( qi(i,j,kp1)  + qi(i,j,k) )
              thvi(i,j,k) = 0.5*( thv(i,j,kp1) + thv(i,j,k) )
           end if
        elseif ( discrete_type == 1 ) then ! upwind
           ui(i,j,k)   = u(i,j,k)
           vi(i,j,k)   = v(i,j,k)
           thli(i,j,k) = thl(i,j,k)
           qti(i,j,k)  = qt(i,j,k)
           qvi(i,j,k)  = qv(i,j,k)
           qli(i,j,k)  = ql(i,j,k)
           qii(i,j,k)  = qi(i,j,k)
           thvi(i,j,k) = thv(i,j,k)
        elseif ( discrete_type == 2 ) then ! upwind (Carpeneter et al. 1990)
           ui(i,j,k)   = interp_carpenter1990_up(IM, JM, LM, i, j, k, zle, zl, u)
           vi(i,j,k)   = interp_carpenter1990_up(IM, JM, LM, i, j, k, zle, zl, v)
           thli(i,j,k) = interp_carpenter1990_up(IM, JM, LM, i, j, k, zle, zl, thl)
           qti(i,j,k)  = interp_carpenter1990_up(IM, JM, LM, i, j, k, zle, zl, qt)
           qvi(i,j,k)  = interp_carpenter1990_up(IM, JM, LM, i, j, k, zle, zl, qv)
           qli(i,j,k)  = interp_carpenter1990_up(IM, JM, LM, i, j, k, zle, zl, ql)
           qii(i,j,k)  = interp_carpenter1990_up(IM, JM, LM, i, j, k, zle, zl, qi)
           thvi(i,j,k) = interp_carpenter1990_up(IM, JM, LM, i, j, k, zle, zl, thv)
        end if
     end do
     end do
  end do

  do j = 1,JM
  do i = 1,IM
     ui(i,j,0)   = u(i,j,1)
     vi(i,j,0)   = v(i,j,1)
     thli(i,j,0) = thl(i,j,1)
     qti(i,j,0)  = qt(i,j,1)
     qvi(i,j,0)  = qv(i,j,1)
     qli(i,j,0)  = ql(i,j,1)
     qii(i,j,0)  = qi(i,j,1)
     thvi(i,j,0) = thv(i,j,1)

     ui(i,j,LM)   = u(i,j,LM)
     vi(i,j,LM)   = v(i,j,LM)
     thli(i,j,LM) = thl(i,j,LM)
     qti(i,j,LM)  = qt(i,j,LM)
     qvi(i,j,LM)  = qv(i,j,LM)
     qli(i,j,LM)  = ql(i,j,LM)
     qii(i,j,LM)  = qi(i,j,LM)
     thvi(i,j,LM) = thv(i,j,LM)
  end do
  end do
end subroutine interp_edmf

!
!
!
real function interp_carpenter1990_up(IM, JM, LM, i, j, k, zle, zl, var)

  implicit none

  integer, intent(in)         :: IM, JM, LM, i, j, k
  real, dimension(IM,JM,LM)   :: zl, var
  real, dimension(IM,JM,0:LM) :: zle 
        
  real :: s_up, s_down, s
  
  if ( k /= 1 ) then 
     s_up   = ( var(i,j,k-1) - var(i,j,k) )/( zl(i,j,k-1) - zl(i,j,k) )
     s_down = ( var(i,j,k) - var(i,j,k+1) )/( zl(i,j,k) - zl(i,j,k+1) )
           
     if ( sign(s_up,s_down) /= s_up ) then
        s = 0.
     else
        s = sign( min(abs(s_up),abs(s_down)), s_up )
     end if
           
     interp_carpenter1990_up = var(i,j,k) + s*( zle(i,j,k) - zl(i,j,k) )
  else
     interp_carpenter1990_up = var(i,j,k)
  end if

end function interp_carpenter1990_up

!
!
!
real function interp_carpenter1990_down(IM, JM, LM, i, j, k, zle, zl, var)

  implicit none
  
  integer, intent(in)         :: IM, JM, LM, i, j, k
  real, dimension(IM,JM,LM)   :: zl, var
  real, dimension(IM,JM,0:LM) :: zle 
        
  real :: s_up, s_down, s

  if ( k /= LM - 1 ) then
     s_up   = ( var(i,j,k) - var(i,j,k+1) )/( zl(i,j,k) - zl(i,j,k+1) )
     s_down = ( var(i,j,k+1) - var(i,j,k+2) )/( zl(i,j,k+1) - zl(i,j,k+2) )
        
     if ( sign(s_up,s_down) /= s_up ) then
        s = 0.
     else
        s = sign( min(abs(s_up),abs(s_down)), s_up )
     end if
           
     interp_carpenter1990_down = var(i,j,k+1) + s*( zle(i,j,k) - zl(i,j,k+1) )
  else
     interp_carpenter1990_down = var(i,j,k+1)
  end if

end function interp_carpenter1990_down

!
! A_star_closure
!
subroutine A_star_closure(IM, JM, LM, &                       ! in
                          wb, au0, debug_flag, &              ! in
                          zle, zl, rho, rhoe, thl, qt, thv, & ! in
                          A, f, zi, Mu0)                      ! out

  integer, intent(in)                      :: IM, JM, LM, debug_flag
  real, intent(in)                         :: wb, au0
  real, dimension(IM,JM,LM), intent(in)    :: zl, thl, qt, thv, rho
  real, dimension(IM,JM,0:LM), intent(in)  :: zle, rhoe
  real, dimension(IM,JM), intent(out)      :: zi, Mu0
  real, dimension(IM,JM,LM), intent(out)   :: A, f

  integer                     :: i, j, k, km1, iter
  real                        :: dthvdz, dzle, B, Mu_next, wu2_next, &
                                 thlu_next, qtu_next, thvu_next, E
  integer, dimension(IM,JM)   :: izsl
  real, dimension(IM,JM)      :: A_star_sum, A_star2_int, Mu, thlu, qtu, thvu, wu2, wu0
  logical, dimension(IM,JM)   :: conv_flag, sl_flag, test_flag

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
              A(i,j,k)        = -sqrt(zl(i,j,k))*dthvdz
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
        wu0(i,j)         = 0.01
        A_star2_int(i,j) = 0.
        do k = LM,izsl(i,j)+1,-1         ! Note: this loop needs to be refactored for column-major order
           dzle = zle(i,j,k-1) - zle(i,j,k)

           A(i,j,k)         = A(i,j,k)/A_star_sum(i,j)
           A_star2_int(i,j) = A_star2_int(i,j) + A(i,j,k)**2.*dzle/rho(i,j,k)
        end do
     end if
  end do
  end do

  ! Find magnitude and height of maximum plume velocity and compute Mu0 accordingly
  wu0(:,:) = 0.01
  do iter = 1,1

     do j = 1,JM
     do i = 1,IM
        ! Initialize plume
        Mu(i,j)   = A(i,j,LM)*( zle(i,j,LM-1) - zle(i,j,LM) )
        wu2(i,j)  = wu0(i,j)**2.
        thlu(i,j) = thl(i,j,LM-1)
        qtu(i,j)  = qt(i,j,LM-1)
        thvu(i,j) = thv(i,j,LM-1)

        test_flag(i,j) = conv_flag(i,j)
     end do
     end do

     f(:,:,:) = 1.
     do k = LM-1,1,-1
        km1 = k - 1

        do j = 1,JM
        do i = 1,IM
           if ( test_flag(i,j) ) then
              dzle = zle(i,j,km1) - zle(i,j,k)
              
              f(i,j,k) = Mu(i,j)/( Mu(i,j) + A(i,j,k)*dzle )

              E = Mu(i,j)*( 1./f(i,j,k) - 1. )/dzle

              Mu_next   =   Mu(i,j)           + dzle*E
              thlu_next = ( Mu(i,j)*thlu(i,j) + dzle*E*thl(i,j,k) )/Mu_next
              qtu_next  = ( Mu(i,j)*qtu(i,j)  + dzle*E*qt(i,j,k) )/Mu_next
              
              thvu_next = thlu_next*( 1. + mapl_vireps*qtu_next )
              
              B = wb*(mapl_grav/thv(i,j,k))*( 0.5*( thvu(i,j) + thvu_next ) - thv(i,j,k) )

              wu2_next = f(i,j,k)**2.*wu2(i,j) + 2.*dzle*B

              if ( B < 0. ) then
                 zi(i,j)        = zle(i,j,k)
                 Mu0(i,j)       = sqrt( wu2(i,j) )/( 2.*zi(i,j)*A_star2_int(i,j) )
                 wu0(i,j)       = Mu0(i,j)*A(i,j,LM)*( zle(i,j,LM-1) - zle(i,j,LM) )/( au0*rhoe(i,j,LM-1) )
                 test_flag(i,j) = .false.
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

